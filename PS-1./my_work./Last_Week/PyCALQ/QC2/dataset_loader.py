#Loads energy levels from the HDF5 dataset and prepares them for the fitting pipeline.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import h5py as h5
import numpy as np
from general.data_reader import LQCD_DATA_READER
from profiler import profiler

class DataLoader:
    def __init__(self, file_path, L, use_ref=True):
        self.file_path = file_path
        self.L = L
        self.use_ref = use_ref
        self.reader = LQCD_DATA_READER(file_path, L)
        self.data = self.reader.load_data()
        top_keys = list(self.data.keys())
        channel_keys = [k for k in top_keys if k.startswith('iso')]
        if channel_keys:
            self.has_channel_layer = True
            self.channel = channel_keys[0]
            print(f"Using channel layer-'{self.channel}'")
        else:
            self.has_channel_layer = False
            self.channel = None
            print("Using direct structure (no channel layer)")

    def _get_energy_data(self, psq, irrep, key):
        if self.has_channel_layer:
            if self.use_ref:
                return self.reader.ref_energy_data(psq, irrep, key)
            else:
                return self.reader.energy_data(psq, irrep, key)
        else:
            try:
                return self.data[psq][irrep][key][:]
            except KeyError:
                print(f"Could not find {psq}/{irrep}/{key}")
                return None

    def _get_free_levels(self, psq, irrep, level):
        if self.has_channel_layer:
            return self.reader.free_levels(psq, irrep, level)
        else:
            try:
                return self.data[psq][irrep].attrs['free_levels'][level]
            except (KeyError, AttributeError):
                return None

    def scan_levels(self):
        levels = []
        psq_keys = [key for key in self.data.keys() if key.startswith('PSQ')]
        for psq in psq_keys:
            if self.has_channel_layer:
                group = self.data[self.channel][psq]
            else:
                group = self.data[psq]
            irreps = list(group.keys())
            for irrep in irreps:
                if irrep.startswith('__') or irrep.startswith('.'):
                    continue
                irrep_group = group[irrep]
                ecm_keys = [key for key in irrep_group.keys() if key.startswith('ecm_')]
                def extract_level(k):
                    num_str = k.replace('ecm_', '').replace('_ref', '')
                    try:
                        return int(num_str)
                    except ValueError:
                        return -1
                ecm_keys.sort(key=extract_level)
                for key in ecm_keys:
                    level_int = extract_level(key)
                    if level_int < 0:
                        continue
                    levels.append({'psq': psq, 'irrep': irrep, 'level': level_int, 'key': key})
        for idx, item in enumerate(levels):
            item['index'] = idx
        if not levels:
            print("Dataset appears to be empty or has an unexpected layout")
        return levels

    def print_levels(self, levels):
        if not levels:
            print("No levels were found")
            return
        print(f"{'Index':<6} {'PSQ':<8} {'Irrep':<10} {'Level':<6} {'Key':<20}")
        print("-"*60)
        for item in levels:
            print(f"{item['index']:<6} {item['psq']:<8} {item['irrep']:<10} {item['level']:<6} {item['key']:<20}")

    def print_levels_sorted_by_energy(self, levels, e_min=None, e_max=None):
        level_data = []
        for item in levels:
            psq, irrep, key = item['psq'], item['irrep'], item['key']
            if not key.endswith('_ref'):
                continue
            arr = self._get_energy_data(psq, irrep, key)
            if arr is None or len(arr)==0:
                continue
            mean = arr[0]
            level_data.append({'index':item['index'], 'mean':mean, 'psq':psq, 'irrep':irrep,
                               'level':item['level'], 'key':key})
        if e_min is not None:
            level_data = [d for d in level_data if d['mean'] >= e_min]
        if e_max is not None:
            level_data = [d for d in level_data if d['mean'] <= e_max]
        level_data.sort(key=lambda x: x['mean'])
        print(f"{'Index':<8} {'Mean Energy':<14} {'PSQ':<8} {'Irrep':<10} {'Level':<6} {'Key':<20}")
        print("-"*80)
        for d in level_data:
            print(f"{d['index']:<8} {d['mean']:14.6f} {d['psq']:<8} {d['irrep']:<10} {d['level']:<6} {d['key']:<20}")
        print(f"\nTotal levels shown: {len(level_data)}")

    def _compute_free_cm_energy(self, free_levels, psq, m1, m2):
        if free_levels is None:
            return None
        import re
        k_values = []
        for fl in free_levels:
            match = re.search(r'\((\d+)\)', fl)
            if match:
                k = float(match.group(1)) * 2.0 * np.pi / self.L
                k_values.append(k)
            else:
                k_values.append(0.0)
        if len(k_values) < 2:
            return None
        E_lab_free = np.sqrt(m1**2 + k_values[0]**2) + np.sqrt(m2**2 + k_values[1]**2)
        psq_num = int(psq.replace("PSQ", ""))
        P_mag = np.sqrt(psq_num) * 2.0 * np.pi / self.L
        E_cm_free = np.sqrt(max(E_lab_free**2 - P_mag**2, 0))
        return float(E_cm_free)

    @profiler.decorator('Dataset Loading')
    def build_dataset(self, indices, levels_scan, m1=1.0, m2=5.862544):
        level_map = {item["index"]: item for item in levels_scan}
        selected = []
        for idx in indices:
            level_info = level_map.get(idx)
            if level_info is None:
                raise ValueError(f"Index {idx} not found in scan list.")
            selected.append(level_info)

        means = []
        bootstrap_rows = []
        metadata = []
        free_energies = []

        for info in selected:
            psq = info['psq']; irrep = info['irrep']; key = info['key']
            arr = self._get_energy_data(psq, irrep, key)
            if arr is None or len(arr)==0:
                raise ValueError(f"Level {psq} {irrep} {key} not found.")
            mean = arr[0]; boots = arr[1:]
            means.append(mean); bootstrap_rows.append(boots)
            metadata.append({'psq': psq, 'irrep': irrep, 'level': info['level']})
            free_levels = self._get_free_levels(psq, irrep, info['level'])
            E_free = self._compute_free_cm_energy(free_levels, psq, m1, m2)
            free_energies.append(E_free)

        lengths = [len(b) for b in bootstrap_rows]
        if len(set(lengths)) != 1:
            raise ValueError(f"Inconsistent bootstrap lengths: {lengths}")

        bootstrap_matrix = np.array(bootstrap_rows)
        mean_vector = np.array(means)
        with profiler.context('Covariance'):
            covariance_matrix = np.cov(bootstrap_matrix, rowvar=True)

        dataset = DataSet(metadata, mean_vector, bootstrap_matrix, covariance_matrix,
                          free_energies=free_energies)
        self._print_diagnostics(dataset)
        return dataset

    def _print_diagnostics(self, dataset):
        print("\nDataset Diagnostics")
        print("Selected levels:")
        for i, meta in enumerate(dataset.metadata):
            print(f"  {i}: {meta['psq']} {meta['irrep']} level {meta['level']}")
            if dataset.free_energies[i] is not None:
                print(f"      free CM energy = {dataset.free_energies[i]:.6f}")
        print(f"Bootstrap matrix shape- {dataset.bootstrap.shape}")
        print(f"Covariance matrix shape- {dataset.covariance.shape}")
        try:
            cond = np.linalg.cond(dataset.covariance)
            print(f"Condition number- {cond:.3e}")
        except:
            pass

class DataSet:
    def __init__(self, metadata, means, bootstrap, covariance, free_energies=None):
        self.metadata = metadata
        self.means = means
        self.bootstrap = bootstrap
        self.covariance = covariance
        self.free_energies = free_energies if free_energies is not None else [None]*len(means)
        self.n_levels = len(means)
        self.n_bootstrap = bootstrap.shape[1] if bootstrap.ndim == 2 else 0
