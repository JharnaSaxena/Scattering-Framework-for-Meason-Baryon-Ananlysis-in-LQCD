"""
plot_figure8.py

Recreate a lattice-QCD-style finite-volume spectrum plot:
- green circles with error bars = data points (can be several per category)
- gray horizontal bands = non-interacting levels
- dashed horizontal lines = physical thresholds, labeled on the right
- categorical x-axis with irrep-style labels
- y-axis ticks at 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25

Fully automated reads data from HDF5 and generates the plot.
Matches the style of paper 2407.13471.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
from collections import defaultdict
from dataset_loader import DataLoader
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 10


def get_noninteracting_energies(L, mass_pion=1.0, psq_max=3, n_shells=4):
    #Compute non-interacting two-particle energies for all momentum frames.

    two_pi_over_L = 2.0 * np.pi / L
    
    psq_to_mom = {
        0: np.array([0, 0, 0]),
        1: np.array([0, 0, 1]),
        2: np.array([1, 1, 0]),
        3: np.array([1, 1, 1]),
    }
    
    ni_energies = {}
    
    for psq in range(psq_max + 1):
        if psq not in psq_to_mom:
            continue
        P = psq_to_mom[psq]
        energies = []
        
        for nx in range(-n_shells, n_shells + 1):
            for ny in range(-n_shells, n_shells + 1):
                for nz in range(-n_shells, n_shells + 1):
                    n = np.array([nx, ny, nz])
                    n1 = n
                    n2 = P - n
                    
                    E1 = np.sqrt(mass_pion**2 + np.dot(n1, n1) * two_pi_over_L**2)
                    E2 = np.sqrt(mass_pion**2 + np.dot(n2, n2) * two_pi_over_L**2)
                    Etot = E1 + E2
                    Psq = np.dot(P, P) * two_pi_over_L**2
                    Ecm = np.sqrt(max(Etot**2 - Psq, 0)) / mass_pion
                    
                    if 1.9 < Ecm < 10.0:
                        Ecm_rounded = round(Ecm, 6)
                        if Ecm_rounded not in energies:
                            energies.append(Ecm_rounded)
        
        ni_energies[psq] = sorted(energies)
    
    return ni_energies


def get_threshold_energies(L, mass_pion=1.0):
    #Compute threshold energies for different channels for πΣ, ππΛ, and K̄N threshold
    two_pi_over_L = 2.0 * np.pi / L
    # ππ threshold at rest
    pi_pi_threshold = 2.0
    # Moving ππ thresholds
    moving_thresholds = {}
    for psq in [1, 2, 3]:
        P_mag = np.sqrt(psq) * two_pi_over_L
        E_pi = np.sqrt(mass_pion**2 + (P_mag/2)**2)
        Etot = 2 * E_pi
        Psq = psq * two_pi_over_L**2
        Ecm = np.sqrt(max(Etot**2 - Psq, 0)) / mass_pion
        moving_thresholds[psq] = Ecm
    return {
        'pi_pi': pi_pi_threshold,
        'moving': moving_thresholds
    }
def plot_figure8_from_hdf5(hdf5_path, L, use_ref=True, energy_window=(6.5, 8.5), save_path=None):
    # Load data
    loader = DataLoader(hdf5_path, L, use_ref=use_ref)
    levels = loader.scan_levels()
    # Select levels in the energy window 
    selected_levels = []
    for item in levels:
        psq = item['psq']
        irrep = item['irrep']
        key = item['key']
        if not key.endswith('_ref'):
            continue
        
        arr = loader._get_energy_data(psq, irrep, key)
        if arr is None or len(arr) == 0:
            continue
        mean = arr[0]
        if energy_window[0] <= mean <= energy_window[1]:
            selected_levels.append({
                'index': item['index'],
                'psq': psq,
                'irrep': irrep,
                'level': item['level'],
                'mean': mean,
                'boots': arr[1:],
                'key': key
            })
    
    if not selected_levels:
        print(f"Warning: No levels found in energy window {energy_window}")
        print("Try adjusting the energy window.")
        return None, None
    
    # Group by (PSQ, irrep) for plotting
    grouped = defaultdict(list)
    for item in selected_levels:
        key = (item['psq'], item['irrep'])
        grouped[key].append(item)
    
    for key in grouped:
        grouped[key] = sorted(grouped[key], key=lambda x: x['mean'])
    
    # Determine x-axis categories
    psq_order = ['PSQ0', 'PSQ1', 'PSQ2', 'PSQ3']
    irrep_priority = {
        'PSQ0': ['G1g', 'G1u', 'Hu'],
        'PSQ1': ['G1', 'G2'],
        'PSQ2': ['G'],
        'PSQ3': ['F1', 'F2', 'G'],
    }
    
    categories = []
    for psq in psq_order:
        for irrep in irrep_priority.get(psq, []):
            key = (psq, irrep)
            if key in grouped:
                label = f"{irrep}({psq[-1]})"
                categories.append(label)
    
    if not categories:
        print("Warning: No categories found for plotting.")
        return None, None
    
    # Build data structures for plotting
    points = {cat: [] for cat in categories}
    bands = {cat: [] for cat in categories}
    level_numbers = {cat: [] for cat in categories}
    ni_energies = get_noninteracting_energies(L)
    thresholds = get_threshold_energies(L)
    for cat in categories:
        psq_num = cat.split('(')[1].strip(')')
        psq = f"PSQ{psq_num}"
        irrep = cat.split('(')[0]
        key = (psq, irrep)
        if key in grouped:
            for item in grouped[key]:
                mean = item['mean']
                err = np.std(item['boots']) if len(item['boots']) > 1 else 0.01
                points[cat].append((mean, err))
                level_numbers[cat].append(item['level'])
        
        psq_int = int(psq_num)
        for ni_e in ni_energies.get(psq_int, [])[:3]:
            if energy_window[0] - 0.1 <= ni_e <= energy_window[1] + 0.1:
                bands[cat].append((ni_e, 0.015))
    
    # Threshold labels
    # These are the physical thresholds from the paper
    threshold_labels = {
        r'$\pi\Sigma$': 6.85,      # πΣ threshold
        r'$\bar{K}N$': 7.20,       # K̄N threshold
        r'$\pi\pi\Lambda$': 7.55,  # ππΛ threshold
    }
    if energy_window[0] - 0.2 <= thresholds['pi_pi'] <= energy_window[1] + 0.2:
        threshold_labels[r'$\pi\pi$'] = thresholds['pi_pi']
    threshold_labels = dict(sorted(threshold_labels.items(), key=lambda x: x[1]))
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(10, 6.5))
    
    n_cat = len(categories)
    x_positions = np.arange(n_cat)
    band_width = 0.5
    for i, cat in enumerate(categories):
        x0 = x_positions[i]
        for center, halfwidth in bands.get(cat, []):
            ax.add_patch(
                plt.Rectangle(
                    (x0 - band_width / 2, center - halfwidth),
                    band_width, 2 * halfwidth,
                    facecolor="0.85", edgecolor="none", zorder=1,
                    alpha=0.7,
                )
            )
    
    xmin, xmax = -0.6, n_cat - 0.4
    threshold_colors = {
        r'$\pi\Sigma$': '#1a1a1a',
        r'$\bar{K}N$': '#2c3e50',
        r'$\pi\pi\Lambda$': '#34495e',
        r'$\pi\pi$': '#4a4a4a',
    }
    
    for label, y in threshold_labels.items():
        color = threshold_colors.get(label, 'black')
        ax.hlines(y, xmin, xmax, colors=color, linestyles="dashed", 
                 linewidth=1.2, zorder=2, alpha=0.8)
        ax.text(
            xmax + 0.05, y, label, 
            va="center", ha="left", fontsize=11,
            color=color
        )
    for i, cat in enumerate(categories):
        x0 = x_positions[i]
        pts = points.get(cat, [])
        n_pts = len(pts)
        
        if n_pts == 0:
            continue
        
        jitter = np.linspace(-0.12, 0.12, n_pts) if n_pts > 1 else [0]
        
        for (val, err), dx in zip(pts, jitter):
            # Error bar
            ax.errorbar(
                x0 + dx, val, yerr=err,
                fmt="o", color="forestgreen", markersize=8,
                elinewidth=1.2, capsize=4, zorder=3,
                markeredgecolor='black', markeredgewidth=0.7,
            )
    for i, cat in enumerate(categories):
        x0 = x_positions[i]
        pts = points.get(cat, [])
        levels = level_numbers.get(cat, [])
        
        n_pts = len(pts)
        jitter = np.linspace(-0.12, 0.12, n_pts) if n_pts > 1 else [0]
        
        for (val, err), dx, lvl in zip(pts, jitter, levels):
            ax.annotate(
                f'{lvl}',
                xy=(x0 + dx, val),
                xytext=(0, 10),
                textcoords='offset points',
                ha='center',
                fontsize=8,
                color='black',
                fontweight='bold'
            )
    
    ax.set_xlim(xmin, xmax + 0.9)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(categories, color="steelblue", fontsize=11)
    y_ticks = [6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{t:.2f}' for t in y_ticks], fontsize=10)
    ax.set_ylabel(r"$E_{\rm cm}/m_\pi$", fontsize=14)
    ax.set_ylim(6.5, 8.4)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.tick_params(direction="in", top=True, right=True, length=5, width=0.8)
    ax.grid(axis='y', alpha=0.15, linestyle='-', linewidth=0.5)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='forestgreen', marker='o', linestyle='None',
               markersize=8, label='Lattice data', markeredgecolor='black'),
        plt.Rectangle((0, 0), 1, 1, facecolor='0.85', edgecolor='none', 
                      label='NI levels', alpha=0.7),
        Line2D([0], [0], color='#2c3e50', linestyle='dashed', linewidth=1.2,
               label='Thresholds'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10, 
             framealpha=0.9, fancybox=True)
    
    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path) if os.path.dirname(save_path) else '.', exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    plt.show()
    return fig, ax

if __name__ == "__main__":
    import os
    
    hdf5_path = os.path.expanduser("~/Desktop/Last_Week/my_work/DataSet.hdf5")
    L = 48
    
    print("Generating Figure 8 with all layers...")
    
    fig, ax = plot_figure8_from_hdf5(
        hdf5_path=hdf5_path,
        L=L,
        use_ref=True,
        energy_window=(6.5, 8.5),
        save_path=os.path.expanduser("~/Desktop/Last_Week/PyCALQ/QC2/figure8_final.pdf")
    )
    
    print(f"Y-axis ticks: 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25")
    print("Thresholds: πΣ, K̄N, ππΛ")
    print("X-axis labels: G1g(0), G1u(0), G1(1), G2(1), G(2), F1(3), F2(3), G(3)")
