#!/usr/bin/env python3
"""
run_fit_from_dataset.py - Fits the ERE parametrization to lattice data.
Uses parallel predictions (ProcessPoolExecutor) for speed.
No bootstrap error estimation – just the fit.
"""
import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dataset_loader import DataLoader
from pipeline_adapter import PSQ_TO_D, full_irrep_label
from morningstar_bmatrix import SingleChannelBMatrix as MorningstarBMatrix
from fitting_driver_canonical import LuscherFitter, PhysicsModule
from profiler import profiler
# CONFIGURATION
HDF5_PATH = os.path.expanduser("~/Desktop/my_work./Last_Week/my_work/DataSet.hdf5")
L = 64.0
m1 = 1.0                     # pion mass (mπ/mπ = 1)
m2 = 5.862544007347314       # Sigma mass (mΣ/mπ = 0.3830/0.06533)

SELECTED_INDICES = [11, 35, 67, 121]   # G1u(0), G1(1), G(2), G(3)
INITIAL_GUESS = np.array([0.047, 0.65])
BOUNDS = [(-10.0, 10.0), (-10.0, 10.0)]
# LOAD DATA
loader = DataLoader(HDF5_PATH, L, use_ref=True)
levels = loader.scan_levels()
dataset = loader.build_dataset(SELECTED_INDICES, levels, m1=m1, m2=m2)
observed_mean = dataset.means
bootstrap_samples = dataset.bootstrap
irrep_list = [full_irrep_label(meta['psq'], meta['irrep']) for meta in dataset.metadata]
d_list = [PSQ_TO_D[meta['psq']] for meta in dataset.metadata]
total_cov_matrix = dataset.covariance
# PHYSICS MODULE
bmatrix_impl = MorningstarBMatrix()
physics = PhysicsModule(L, m1, m2, bmatrix_impl=bmatrix_impl)
# FITTER – PARALLEL PREDICTIONS (max_workers = number of levels)
fitter = LuscherFitter(
    observed_mean=observed_mean,
    bootstrap_samples=bootstrap_samples,
    irrep_list=irrep_list,
    d_list=d_list,
    L=L,
    m1=m1,
    m2=m2,
    physics=physics,
    cov_matrix=total_cov_matrix,
    free_energies=dataset.free_energies,
    max_workers=None,   # auto = 4 workers for 4 levels
    debug_objective=False,
    verbose_predict=False,
)

print("RUNNING OPTIMIZER (parallel predictions)")
print(f"  Levels: {len(irrep_list)}")
print(f"  Workers: {fitter.max_workers}")
start_time = time.perf_counter()
result = fitter.fit(
    initial_guess=INITIAL_GUESS,
    bounds=BOUNDS,
    verbose=True,
    method='nelder-mead',
    maxiter=5000,
    compute_vij=True,
)
fit_time = time.perf_counter() - start_time
print(f"\nOptimization finished in {fit_time:.2f} seconds")
fitter.save_results("fit_results.json")
print("PROFILING REPORT")
profiler.report()
print("FIT SUMMARY")
print(f"  a = {result.params[0]:.6f} +/- {result.errors[0]:.6f}")
print(f"  b = {result.params[1]:.6f} +/- {result.errors[1]:.6f}")
print(f"  χ² = {result.chi2:.6f}, ndof = {result.ndof}")
print(f"  χ²/ndof = {result.reduced_chi2:.6f}")
print(f"  Parallel workers: {fitter.max_workers}")
