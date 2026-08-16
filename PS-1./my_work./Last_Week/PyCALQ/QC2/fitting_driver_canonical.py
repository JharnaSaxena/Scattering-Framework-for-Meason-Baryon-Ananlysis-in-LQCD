#!/usr/bin/env python3
"""
Main fitting driver energy‑based fit with root finding
Now uses covariance estimation via stats.parameter_covariance.
"""
import numpy as np
from typing import List, Tuple, Optional, Callable
import warnings
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from profiler import profiler
from root_finder import RootFinder
from ere import ERE
from stats import (
    bootstrap_covariance, chi2, parameter_covariance,
    parameter_errors, reduced_chi2, standardized_residuals,
    correlation_matrix, aic, bic, FitResult
)
from tools.kinematics import KinematicVars, compute_kinematics as kinematics_compute

try:
    import scipy.optimize as opt
except ImportError:
    raise RuntimeError("scipy.optimize is required.")

try:
    multiprocessing.set_start_method('fork', force=True)
except RuntimeError:
    pass


def _predict_single_level(
    irrep: str,
    d: Tuple[int, int, int],
    e_obs: float,
    params: np.ndarray,
    L: float,
    m1: float,
    m2: float,
    free_energy: Optional[float],
    level_index: int,
    verbose: bool = False,
) -> float:
    if verbose:
        print(f"[DEBUG] Predicting level {level_index}: irrep={irrep}, d={d}, e_obs={e_obs:.6f}")

    try:
        from morningstar_bmatrix import SingleChannelBMatrix as MorningstarBMatrix
        from fitting_driver_canonical import PhysicsModule
        from root_finder import RootFinder
        from ere import ERE
    except Exception as e:
        warnings.warn(f"Import failed in _predict_single_level for level {level_index}: {e}")
        return e_obs

    bmatrix_impl = MorningstarBMatrix()
    physics = PhysicsModule(L, m1, m2, bmatrix_impl=bmatrix_impl)
    root_finder = RootFinder(
        benchmark=True,
        root_tolerance=1e-4,
        continuity_tol=0.3,
        local_scan_points=50,
        verbose=False,
        debug=False,
    )
    ere = ERE(params)
    omega = physics.build_omega(irrep, d, e_obs, ere)

    root = root_finder.find_root_near_guess(
        omega,
        x_guess=e_obs,
        prev_root=None,
        reference_energy=free_energy,
        level_index=level_index,
    )
    if verbose:
        print(f"[DEBUG] Level {level_index} root found: {root:.6f}")
    return root


class PhysicsModule:
    def __init__(self, L: float, m1: float, m2: float, bmatrix_impl=None):
        self.L = float(L)
        self.m1 = float(m1)
        self.m2 = float(m2)
        self.bmatrix = bmatrix_impl
        self._kin_cache = {}
        if self.bmatrix is None:
            warnings.warn("No B-matrix using dummy B=0.", stacklevel=2)
            self._setup_dummy_bmatrix()

    def _setup_dummy_bmatrix(self):
        class DummyBMatrix:
            def compute(self, irrep, kin):
                return 0.0
        self.bmatrix = DummyBMatrix()

    @profiler.decorator('Kinematics')
    def compute_kinematics(self, E_cm: float, d: Tuple[int, int, int]) -> KinematicVars:
        key = (round(E_cm, 12), d)
        if key in self._kin_cache:
            return self._kin_cache[key]
        kin = kinematics_compute(E_cm, d, self.m1, self.m2, self.L, Mref=1.0)
        self._kin_cache[key] = kin
        return kin

    def compute_bmatrix(self, irrep: str, kin: KinematicVars) -> float:
        return float(self.bmatrix.compute(irrep, kin))

    def compute_kinv(self, kin: KinematicVars, ere: ERE) -> float:
        return ere.compute_kinv(kin)

    def build_omega(self, irrep: str, d: Tuple[int, int, int],
                     E_cm: float, ere: ERE) -> Callable[[float], float]:
        def omega(E: float) -> float:
            with profiler.context('Omega'):
                if E <= 0.0:
                    return 1e12
                kin = self.compute_kinematics(E, d)
                B = self.compute_bmatrix(irrep, kin)
                Kinv = self.compute_kinv(kin, ere)
                return Kinv - B
        return omega


class LuscherFitter:
    def __init__(
        self,
        observed_mean: np.ndarray,
        bootstrap_samples: np.ndarray,
        irrep_list: List[str],
        d_list: List[Tuple[int, int, int]],
        L: float,
        m1: float,
        m2: float,
        physics: Optional[PhysicsModule] = None,
        root_finder: Optional[RootFinder] = None,
        cov_matrix: Optional[np.ndarray] = None,
        free_energies: Optional[np.ndarray] = None,
        debug_objective: bool = False,
        max_workers: Optional[int] = None,
        verbose_predict: bool = False,
    ):
        self.observed_mean = np.asarray(observed_mean, dtype=float)
        self.bootstrap_samples = np.asarray(bootstrap_samples)
        self.irrep_list = irrep_list
        self.d_list = d_list
        self.L = float(L)
        self.m1 = float(m1)
        self.m2 = float(m2)
        self.n_data = len(observed_mean)
        self.physics = physics or PhysicsModule(L, m1, m2)
        self.root_finder = root_finder or RootFinder(
            benchmark=True, root_tolerance=1e-4,
            continuity_tol=0.3, local_scan_points=50, verbose=False
        )
        self.cov_matrix = cov_matrix if cov_matrix is not None else bootstrap_covariance(bootstrap_samples)
        self._n_evaluations = 0
        self._result = None
        self.free_energies = np.asarray(free_energies) if free_energies is not None else None
        self.debug_objective = debug_objective
        self._last_roots = [None] * self.n_data
        self.threshold = self.m1 + self.m2
        self.verbose_predict = verbose_predict
        
        # Set number of workers
        if max_workers is None:
            import os
            self.max_workers = min(len(self.irrep_list), os.cpu_count() or 4)
        else:
            self.max_workers = max_workers
        
        # For serial fallback
        self._use_parallel = self.max_workers > 1

    def _get_param_labels(self, params):
        return ERE(params)._get_coeff_labels()

    @profiler.decorator('Predict Energies')
    def predict_energies(self, params: np.ndarray) -> np.ndarray:
        """
        Predict energy levels in parallel using ProcessPoolExecutor.
        Each level is independent, so we can parallelize across levels.
        """
        # Build argument list for each level
        args_list = []
        for idx, (irrep, d, e_obs) in enumerate(zip(self.irrep_list, self.d_list, self.observed_mean)):
            free_energy = self.free_energies[idx] if self.free_energies is not None else None
            args_list.append((irrep, d, e_obs, params, self.L, self.m1, self.m2, free_energy, idx, self.verbose_predict))
        if len(args_list) <= 1 or self.max_workers <= 1:
            return np.array([_predict_single_level(*args) for args in args_list])

        # Use ProcessPoolExecutor for parallel prediction
        try:
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [executor.submit(_predict_single_level, *args) for args in args_list]
                results = [future.result() for future in as_completed(futures)]
                idx_future = [(args_list[i][8], futures[i]) for i in range(len(futures))]
                results = [0] * len(futures)
                for idx, future in idx_future:
                    results[idx] = future.result()
                return np.array(results)
        except Exception as e:
            warnings.warn(f"Parallel prediction failed ({e}); falling back to serial")
            return np.array([_predict_single_level(*args) for args in args_list])

    def objective(self, params):
        self._n_evaluations += 1
        profiler.increment_counter('chi2_evaluations')
        try:
            with profiler.context('Chi2'):
                predicted = self.predict_energies(params)
                chi2_val = chi2(self.observed_mean, predicted, self.cov_matrix)
                if self.debug_objective:
                    print(f"\n[OBJECTIVE DEBUG] params={params}")
                    print(f"  predicted={predicted}")
                    print(f"  residuals={self.observed_mean - predicted}")
                    print(f"  chi2={chi2_val:.6f}")
                return float(chi2_val)
        except Exception as e:
            warnings.warn(f"Objective failed: {e}")
            return 1e12

    def objective_with_prediction(self, params):
        predicted = self.predict_energies(params)
        chi2_val = chi2(self.observed_mean, predicted, self.cov_matrix)
        return chi2_val, predicted

    def vij(self, params: np.ndarray, epsilon: float = 1e-5) -> np.ndarray:
        """
        Compute parameter covariance using the Fisher information matrix.
        Uses stats.parameter_covariance with PDG-style χ² scaling.
        """
        params = np.asarray(params, dtype=float)
        n_params = len(params)
        n_data = self.n_data
        jac = np.zeros((n_data, n_params))
        for i in range(n_params):
            h = epsilon
            params_plus = params.copy(); params_plus[i] += h
            params_minus = params.copy(); params_minus[i] -= h
            pred_plus = self.predict_energies(params_plus)
            pred_minus = self.predict_energies(params_minus)
            jac[:, i] = (pred_plus - pred_minus) / (2.0 * h)
        predicted = self.predict_energies(params)
        chi2_val = chi2(self.observed_mean, predicted, self.cov_matrix)
        cov_par = parameter_covariance(
            jacobian=jac,
            data_cov=self.cov_matrix,
            chi2_value=chi2_val,
            n_data=n_data,
            n_params=n_params,
        )
        return cov_par

    def bootstrap_parameter_errors(self, params, n_refit=50, verbose=True):
        bootstrap_samples = self.bootstrap_samples  # shape (n_data, n_bootstrap)
        n_data, n_boot = bootstrap_samples.shape
        n_use = min(n_refit, n_boot)

        if verbose:
            print(f"Refitting {n_use} bootstrap replicas (this may take a while)...")
            print(f"Using {self.max_workers} workers for parallel predictions within each fit.")

        param_list = []
        for i in range(n_use):
            replica_mean = bootstrap_samples[:, i]

            # Create a temporary fitter with this replica's data
            temp_fitter = LuscherFitter(
                observed_mean=replica_mean,
                bootstrap_samples=bootstrap_samples,
                irrep_list=self.irrep_list,
                d_list=self.d_list,
                L=self.L,
                m1=self.m1,
                m2=self.m2,
                physics=self.physics,
                cov_matrix=self.cov_matrix,
                free_energies=self.free_energies,
                max_workers=self.max_workers,  # Reuse parallel settings
                debug_objective=False,
                verbose_predict=False,
            )
            # Fit with reduced iterations for speed
            result = temp_fitter.fit(
                initial_guess=params,
                bounds=[(-10, 10), (-10, 10)],
                verbose=False,
                method='nelder-mead',
                maxiter=1000,
                compute_vij=False,
            )
            param_list.append(result.params)
            if verbose and (i+1) % 10 == 0:
                print(f"  Completed {i+1} bootstrap fits")

        if len(param_list) == 0:
            return None, None, None

        param_array = np.array(param_list)
        errors = np.std(param_array, axis=0)
        cov_par = np.cov(param_array, rowvar=False)
        corr_par = np.corrcoef(param_array, rowvar=False)
        return errors, cov_par, corr_par

    def fit(self, initial_guess, bounds=None, verbose=True, param_labels=None,
            method='nelder-mead', maxiter=5000, compute_vij=True):
        if param_labels is None:
            param_labels = self._get_param_labels(initial_guess)
        if verbose:
            print(f"\n  Optimizer: {method.upper()}")
            print(f"  Initial guess: {initial_guess}")
            print(f"  Parallel workers: {self.max_workers}")
            if bounds:
                print(f"  Limits: {bounds}")

        self._n_evaluations = 0
        self._last_roots = [None] * self.n_data

        bounds_scipy = bounds if bounds is not None else None

        result_opt = opt.minimize(
            self.objective,
            initial_guess,
            method='Nelder-Mead',
            bounds=bounds_scipy,
            options={'maxiter': maxiter, 'xatol': 1e-6, 'fatol': 1e-6},
        )
        if verbose:
            print(f"\n  FULL OPTIMIZER RESULT (Nelder-Mead):")
            print(f"    success: {result_opt.success}")
            print(f"    message: {result_opt.message}")
            print(f"    nfev: {result_opt.nfev}")
            print(f"    nit: {result_opt.nit}")
            print(f"    x: {result_opt.x}")
            print(f"    fun: {result_opt.fun}")

        best_params = result_opt.x
        best_chi2 = result_opt.fun
        success = result_opt.success
        message = result_opt.message
        n_iter = result_opt.nfev

        final_chi2, predicted = self.objective_with_prediction(best_params)
        n_params = len(best_params)
        dof = self.n_data - n_params

        residuals = self.observed_mean - predicted
        pulls = standardized_residuals(self.observed_mean, predicted, self.cov_matrix)
        red_chi2 = reduced_chi2(final_chi2, self.n_data, n_params)
        aic_val = aic(final_chi2, n_params, self.n_data)
        bic_val = bic(final_chi2, n_params, self.n_data)

        if compute_vij:
            try:
                cov_par = self.vij(best_params)
                errors = np.sqrt(np.abs(np.diag(cov_par)))
                corr_par = correlation_matrix(cov_par)
            except Exception as e:
                if verbose:
                    print(f"  Warning: vij failed: {e}")
                cov_par = np.eye(n_params) * np.nan
                errors = np.full(n_params, np.nan)
                corr_par = np.eye(n_params) * np.nan
        else:
            cov_par = np.eye(n_params) * np.nan
            errors = np.full(n_params, np.nan)
            corr_par = np.eye(n_params) * np.nan

        self._result = FitResult(
            params=best_params,
            chi2=final_chi2,
            ndof=dof,
            reduced_chi2=red_chi2,
            cov_params=cov_par,
            errors=errors,
            corr_params=corr_par,
            pulls=pulls,
            residuals=residuals,
            predicted=predicted,
            aic=aic_val,
            bic=bic_val,
            success=success,
            message=message,
            n_iter=n_iter,
            n_evaluations=self._n_evaluations,
            param_labels=param_labels,
        )

        if verbose:
            self._print_results()
        return self._result

    def _print_results(self):
        r = self._result
        print("FIT RESULTS")
        print(f"  Success: {r.success}")
        print(f"  χ² = {r.chi2:.6f}, ndof = {r.ndof}, χ²/ndof = {r.reduced_chi2:.6f}")
        print(f"  AIC = {r.aic:.6f}, BIC = {r.bic:.6f}")
        print(f"  Evaluations: {r.n_evaluations}")
        print("\n  Parameters:")
        for label, val, err in zip(r.param_labels, r.params, r.errors):
            err_str = f"{err:12.6f}" if np.isfinite(err) else "      nan"
            print(f"    {label:>10} = {val:>12.6f} +/- {err_str}")
        print("\n  Pulls:")
        for irrep, pull in zip(self.irrep_list, r.pulls):
            print(f"    {irrep:>10}: {pull:>8.4f}")
        print("=" * 70)

    def save_results(self, filename):
        if self._result is None:
            print("No results to save")
            return
        data = {
            'fitter': 'LuscherFitter',
            'fitting_method': 'energy_based',
            'max_workers': self.max_workers,
            'parameters': {'L': self.L, 'm1': self.m1, 'm2': self.m2},
            'result': self._result.to_dict(),
            'input_data': {
                'observed_mean': self.observed_mean.tolist(),
                'irreps': self.irrep_list,
                'd_vectors': [list(d) for d in self.d_list],
                'cov_matrix': self.cov_matrix.tolist(),
            },
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Results saved to: {filename}")

if __name__ == "__main__":
    print("Fitting Driver (with unified covariance)")
