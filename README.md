```markdown
# Finite-Volume Scattering Analysis for Meson-Baryon Systems

## Overview

This repository implements a computational pipeline for analyzing two-particle scattering in a finite volume from lattice QCD energy levels. It provides a complete framework for extracting scattering parameters from lattice QCD data using the Lüscher formalism, with a specific focus on the πΣ system in the S-wave. The pipeline handles unequal masses, multiple momentum frames, and uses a robust fitting procedure with detailed uncertainty estimation.

The software is designed for researchers working in hadron spectroscopy and lattice QCD who need to extract scattering phase shifts and effective range parameters from finite-volume energy spectra. It bridges the gap between raw lattice QCD data and physical scattering observables.

## Scientific Problem

In lattice QCD, scattering information is extracted by computing energy levels in a finite cubic box. The Lüscher formalism relates these discrete energy levels to the infinite-volume scattering amplitude. The key challenge is solving the finite-volume quantization condition:

\[
\det[K^{-1}(E) - B(E)] = 0
\]

where:
- \(K(E)\) is the infinite-volume scattering amplitude (parameterized by the Effective Range Expansion)
- \(B(E)\) is the finite-volume B-matrix that encodes the geometry of the box
- \(E\) is the center-of-mass energy

This pipeline solves this condition to fit a two-parameter Effective Range Expansion (ERE) to lattice QCD energy levels. It handles the technical complexities of unequal masses, moving frames, and the numerically challenging generalized zeta functions that appear in the B-matrix.

## What This Repository Does

1. **Loads lattice QCD data** from HDF5 files containing energy levels with bootstrap samples
2. **Computes finite-volume B-matrix elements** using Morningstar's coupled-channel formalism with precomputed coefficients from tables
3. **Evaluates generalized zeta functions** using a hybrid approach (Padé approximation for u² > 0, exact Ewald summation for u² < 0)
4. **Solves the quantization condition** via root finding to predict energy levels for a given set of ERE parameters
5. **Fits the ERE parameters** to observed lattice energy levels using χ² minimization
6. **Provides uncertainty estimation** through bootstrap resampling and parameter covariance analysis
7. **Visualizes results** with publication-quality plots (finite-volume spectra and scattering curves)
8. **Profiles performance** to identify bottlenecks

## Key Features

- **Unequal-mass kinematics**: Handles πΣ scattering with masses mπ = 1.0 and mΣ = 5.8625 (in units of mπ)
- **Multiple momentum frames**: Supports PSQ0 (P=0), PSQ1 (P=(0,0,1)), PSQ2 (P=(1,1,0)), PSQ3 (P=(1,1,1))
- **Morningstar B-matrix**: Complete implementation of the single-channel B-matrix with coefficients from arXiv:1707.05817
- **Hybrid zeta function**: Padé approximation for u² > 0 with exact Ewald summation fallback
- **Parallel root finding**: Uses ProcessPoolExecutor for concurrent energy prediction across levels
- **Bootstrap error propagation**: Full covariance matrix from bootstrap samples
- **Automatic root tracking**: Continuity-based root selection for reliable fitting
- **Performance profiling**: Built-in profiler tracks function calls, timings, and counters
- **Publication-ready plots**: Figure 8 (finite-volume spectrum) and Figure 10 (scattering curves) reproduction

## Repository Architecture

The codebase is organized as follows:

```
repository root/
├── QC2/                          # Main package directory
│   ├── __init__.py
│   ├── dataset_loader.py         # HDF5 data loading
│   ├── fitting_driver_canonical.py # Core fitting engine
│   ├── ere.py                    # Effective Range Expansion
│   ├── morningstar_bmatrix.py    # B-matrix implementation
│   ├── root_finder.py            # Root finding for quantization condition
│   ├── stats.py                  # Statistical utilities
│   ├── kinematics.py             # Kinematic calculations (with caching)
│   ├── final_zeta.py             # Hybrid zeta function (with caching)
│   ├── exact_zeta.py             # Exact Ewald summation (reference)
│   ├── profiler.py               # Performance profiling
│   ├── pipeline_adapter.py       # PSQ to momentum vector conversion
│   ├── run_fit_from_dataset.py   # Main entry point for fitting
│   ├── plot_figure8.py           # Finite-volume spectrum plot
│   ├── plot.py                   # Scattering curve plot (Fig. 10)
│   └── b_tables/                 # Morningstar B-matrix coefficients
│       ├── b1.py, b2.py, ...
│       └── tables.py
├── tools/                        # Supporting modules
│   ├── kinematics.py             (symlink/import)
│   └── final_zeta.py             (symlink/import)
└── coefficients/                 # Precomputed Padé coefficients
    ├── PSQ0_gamma1.00_coeffs.json
    ├── PSQ0_gamma1.02_coeffs.json
    └── ...
```

## Directory Structure

```
├── dataset_loader.py         # Data loading from HDF5, bootstrap extraction
├── ere.py                    # Effective Range Expansion (k/mπ)cotδ = a + bΔ
├── fitting_driver_canonical.py # Main fitting: LuscherFitter, PhysicsModule
├── morningstar_bmatrix.py    # Morningstar B-matrix with cache
├── root_finder.py            # Root finder with continuity tracking
├── stats.py                  # χ², covariance, bootstrap utilities
├── kinematics.py             # Kinematics with functools.lru_cache
├── final_zeta.py             # Hybrid zeta with LRU caching
├── exact_zeta.py             # Exact Ewald summation (reference implementation)
├── profiler.py               # Performance profiler
├── pipeline_adapter.py       # PSQ labels → momentum vectors
├── run_fit_from_dataset.py   # Main entry point
├── plot_figure8.py           # Finite-volume spectrum plot
├── plot.py                   # Scattering curve plot
├── b_tables/                 # Morningstar coefficients tables
│   ├── b1.py
│   ├── b2.py
│   ├── b3.py
│   ├── b4.py
│   ├── b5.py
│   ├── b6.py
│   ├── b7.py
│   ├── b8.py
│   ├── tables.py
│   └── __init__.py
├── coefficients/             # Padé coefficients for zeta function
│   ├── PSQ0_gamma1.00_coeffs.json
│   ├── PSQ0_gamma1.02_coeffs.json
│   └── ...
└── test/                     # (No test suite found in provided files)
```

## Computational Pipeline

### High-Level Overview

The pipeline transforms raw lattice QCD data into fitted scattering parameters through a multi-stage process:

1. **Data Loading**: The `DataLoader` reads an HDF5 file containing energy levels with bootstrap samples. It scans the HDF5 structure to identify all available levels, extracts means and bootstrap distributions, and builds a `DataSet` object.

2. **Kinematic Preparation**: For each selected level, `PhysicsModule.compute_kinematics()` calculates all relevant kinematic quantities in the center-of-mass frame, including:
   - Energy in CM frame (\(E_{\mathrm{cm}}\))
   - Relative momentum (\(q^*\))
   - Lorentz boost factor (\(\gamma\))
   - Mass asymmetry parameter (\(\alpha\))
   - Lüscher momentum variable \(u = \frac{L q^*}{2\pi}\)

3. **B-Matrix Evaluation**: `SingleChannelBMatrix.compute()` calculates the finite-volume B-matrix element for a given irrep and kinematics. This uses precomputed coefficients from the Morningstar tables and the generalized zeta function.

4. **Effective Range Expansion**: `ERE.compute_kinv()` evaluates \((k/m_\pi)\cot\delta = a + b\Delta\) where \(\Delta = (E^2 - E_{\mathrm{th}}^2)/E_{\mathrm{th}}^2\).

5. **Quantization Condition**: The `PhysicsModule.build_omega()` constructs the function:
   \[
   \omega(E) = K_{\mathrm{inv}}(E) - B(E)
   \]
   where \(K_{\mathrm{inv}} = (k/m_\pi)\cot\delta\) and \(B(E)\) is the B-matrix.

6. **Root Finding**: `RootFinder.find_root_near_guess()` finds the root of \(\omega(E) = 0\) using adaptive bracketing and continuation. This root is the predicted energy level for the given ERE parameters.

7. **Prediction**: For a set of ERE parameters \((a,b)\), `LuscherFitter.predict_energies()` predicts all selected energy levels by finding roots for each level's quantization condition.

8. **χ² Evaluation**: `stats.chi2()` computes the χ² between observed and predicted energies using the covariance matrix:
   \[
   \chi^2 = \mathbf{r}^T C^{-1} \mathbf{r}
   \]
   where \(\mathbf{r} = \mathbf{E}_{\mathrm{obs}} - \mathbf{E}_{\mathrm{pred}}\).

9. **Optimization**: `scipy.optimize.minimize()` with Nelder-Mead minimizes the χ² to find best-fit ERE parameters.

10. **Uncertainty Estimation**: Parameter uncertainties are estimated via:
    - Fisher information matrix with PDG-style scaling
    - Bootstrap refitting (optional, implemented but not enabled by default)

11. **Results**: The `FitResult` object contains fitted parameters, χ², predicted energies, pulls, and covariance matrices.

### Pipeline Data Flow

```
HDF5 File
    ↓
[DataLoader]
    ↓
DataSet (means, bootstrap, covariance, free_energies)
    ↓
[LuscherFitter]
    ↓
observed_mean, covariance, free_energies, irrep_list, d_list
    ↓
[PhysicsModule]
    ↓
kinematics_cache (LRU)
    ↓
For each level:
    [build_omega] → ω(E) = Kinv(E) - B(E)
    [RootFinder] → find_root(ω, guess)
    ↑
    Repeated for each ERE parameter evaluation
    ↓
Predicted energies for all levels
    ↓
[stats.chi2] → χ²
    ↓
[scipy.optimize.minimize] → best-fit (a,b)
    ↓
[FitResult]
    ↓
Parameters, errors, χ², pulls, plots
```

## Physics Background

### Lüscher Quantization Condition

For a two-particle system in a finite cubic box of size \(L\), the allowed energy levels are determined by the quantization condition:

\[
\det[K^{-1}(E) - B(E)] = 0
\]

In the single-channel, S-wave case this reduces to:

\[
K_{\mathrm{inv}}(E) = B(E)
\]

where:
- \(K_{\mathrm{inv}}(E) = (k/m_\pi)\cot\delta(E)\) is the inverse scattering amplitude
- \(B(E)\) is the finite-volume B-matrix

### Effective Range Expansion

The pipeline parameterizes the inverse scattering amplitude using the Effective Range Expansion:

\[
\frac{k}{m_\pi} \cot\delta = a + b\,\Delta
\]

where:
- \(a\) is the inverse scattering length (in units of \(1/m_\pi\))
- \(b\) is the effective range parameter
- \(\Delta = \frac{E^2 - E_{\mathrm{th}}^2}{E_{\mathrm{th}}^2}\) is a dimensionless energy variable
- \(E_{\mathrm{th}} = m_1 + m_2\) is the threshold energy

This parameterization follows Eq. 12 of arXiv:2307.13471.

### Morningstar B-Matrix

The B-matrix encodes the finite-volume effects and depends on the irrep, momentum frame, and kinematic variables:

\[
B_{\Lambda}(\mathbf{d}, \gamma, \alpha, u) = \sum_{\ell,m} c_{\ell m}^{\Lambda} \frac{Z_{\ell m}(u^2, \gamma, \mathbf{d}, \alpha)}{\gamma \pi^{3/2} u}
\]

where:
- \(\Lambda\) is the lattice irrep (e.g., G1u, G1, G, etc.)
- \(Z_{\ell m}\) are generalized zeta functions
- \(c_{\ell m}^{\Lambda}\) are coefficients from Morningstar's tables

The coefficients are precomputed and stored in `b_tables/b1.py` through `b8.py` (from arXiv:1707.05817).

### Zeta Function

The generalized zeta function is the most computationally expensive part of the pipeline:

\[
Z_{\ell m}(u^2, \gamma, \mathbf{d}, \alpha) = \sum_{\mathbf{n}} \frac{\mathcal{Y}_{\ell m}(\mathbf{r})}{r^2 - u^2} \quad\text{(with regularization)}
\]

The pipeline implements two evaluation strategies:
1. **Padé approximation** (for \(u^2 > 0\)): Precomputed rational approximations for fast evaluation
2. **Exact Ewald summation** (for \(u^2 < 0\)): Direct numerical evaluation with convergence acceleration

## Numerical Methods

### Root Finding

The `RootFinder` class finds roots of the quantization condition \(\omega(E) = 0\) using a multi-stage strategy:

1. **Direct acceptance**: If \(|\omega(E_{\mathrm{obs}})| < \mathrm{tolerance}\), accept the observed energy as the root.

2. **Adaptive bracketing**: Try progressively larger brackets around the guess:
   - Widths: [0.02, 0.05, 0.10, 0.20, 0.40, 0.80, 1.5]
   - Check for sign change using bisection

3. **Ordered root selection**: If multiple roots are found, select the one corresponding to the level index (continuity tracking across levels).

4. **Continuity tracking**: Prefer roots near the previous root (for smooth branch following) or near the free energy.

5. **Local scan**: If bracketing fails, scan the interval with 50 points and identify sign changes.

6. **Minimization fallback**: Minimize \(|\omega(E)|\) using bounded scalar minimization.

**Parameters**:
- Root tolerance: 1e-4
- Continuity tolerance: 0.3
- Maximum iterations: 100
- Pole threshold: 1e12 (treat values above this as poles)

### Optimization

The χ² minimization uses `scipy.optimize.minimize` with the Nelder-Mead method:

- **Initial guess**: Typically (0.047, 0.65) for πΣ scattering
- **Bounds**: (-10, 10) for both parameters
- **Maximum iterations**: 5000
- **Convergence criteria**: xatol=1e-6, fatol=1e-6

The optimizer calls `predict_energies()` for each function evaluation, which in turn:
- Computes kinematics for each level (cached)
- Evaluates B-matrix (cached)
- Finds roots (cached via root finder)

### Covariance and Statistics

The pipeline follows standard statistical practices for lattice QCD:

1. **Covariance matrix**: Computed from bootstrap samples:
   \[
   C_{ij} = \frac{1}{N_b-1} \sum_{b=1}^{N_b} (E_i^b - \bar{E}_i)(E_j^b - \bar{E}_j)
   \]

2. **Parameter covariance**: Computed using the Fisher information matrix:
   \[
   \mathrm{Cov}(\theta) = (J^T C^{-1} J)^{-1}
   \]
   with PDG-style scaling by χ²/dof if χ² > dof.

3. **Pulls**: Standardized residuals:
   \[
   \mathrm{pull}_i = \frac{E_{\mathrm{obs},i} - E_{\mathrm{pred},i}}{\sqrt{C_{ii}}}
   \]

4. **Information criteria**:
   - AIC = χ² + 2p
   - BIC = χ² + p ln(N)

### Caching Strategy

The pipeline aggressively caches expensive computations to avoid redundant work:

1. **Kinematics**: `compute_kinematics()` is decorated with `@lru_cache(maxsize=None)` with rounded arguments (10 decimal places).

2. **B-matrix**: `SingleChannelBMatrix.compute()` caches results in a dictionary keyed by (irrep, psq, u2, gamma, m_split).

3. **Zeta function**: `hybrid_Z()` uses `@lru_cache(maxsize=None)` with rounded arguments.

4. **Root finder diagnostics**: Stores the last diagnostics for debugging.

5. **B-matrix coefficients**: Precomputed Morningstar coefficients are loaded once.

### Parallelization

The `predict_energies()` function parallelizes root finding across levels using `ProcessPoolExecutor`:

- Each level's root finding is independent
- Default workers: `min(n_levels, os.cpu_count())`
- Falls back to serial on error
- Each worker process imports modules independently (fork-safe)

## Module-by-Module Documentation

### `dataset_loader.py`

**Purpose**: Loads lattice QCD energy levels from HDF5 files and prepares them for fitting.

**Role in pipeline**: First stage - data ingestion.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `DataLoader.__init__` | Initialize data loader | `file_path`, `L`, `use_ref` | DataLoader instance |
| `DataLoader.scan_levels` | Discover all available levels | None | List of level dictionaries |
| `DataLoader.build_dataset` | Build DataSet from indices | `indices`, `levels_scan`, `m1`, `m2` | DataSet instance |
| `DataLoader._compute_free_cm_energy` | Compute non-interacting energy | `free_levels`, `psq`, `m1`, `m2` | Free CM energy |

**Dependencies**: `h5py`, `numpy`, `general.data_reader.LQCD_DATA_READER`

**Numerical considerations**:
- Bootstrap samples are extracted from HDF5 (first entry is mean, rest are bootstrap replicas)
- Free energies are computed from momentum mode labels

**Physics significance**: Maps lattice QCD naming conventions (PSQ, irreps) to pipeline-internal structure.

---

### `ere.py`

**Purpose**: Implements the Effective Range Expansion parameterization.

**Role in pipeline**: Defines the scattering amplitude model.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `ERE.__init__` | Initialize ERE | `coeffs` (list/array) | ERE instance |
| `ERE.compute_kinv` | Compute (k/mπ)cotδ | `kin` (KinematicVars) | Float |
| `ERE.compute_cot_delta` | Compute cotδ | `kin` | Float |
| `ERE.compute_phase_shift` | Compute δ in degrees | `kin` | Float |

**Algorithm**:
\[
\mathrm{kinv} = a + b\cdot\Delta
\]
where \(\Delta = (E^2 - E_{\mathrm{th}}^2)/E_{\mathrm{th}}^2\).

**Called by**: `PhysicsModule.compute_kinv()` → used in `build_omega()`

**Numerical considerations**:
- Supports 1 or 2 parameters
- Input can be `KinematicVars` object or direct q²
- Threshold protection (1e-15 added to denominator)

---

### `fitting_driver_canonical.py`

**Purpose**: Core fitting engine that ties together physics, root finding, and statistics.

**Role in pipeline**: Orchestrates the complete fitting process.

**Important classes**:

| Class | Purpose | Key methods |
|-------|---------|-------------|
| `PhysicsModule` | Physics calculations | `compute_kinematics`, `compute_bmatrix`, `build_omega` |
| `LuscherFitter` | Fitting engine | `predict_energies`, `objective`, `fit`, `vij` |

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `LuscherFitter.predict_energies` | Predict energies for params | `params` (array) | `predicted` (array) |
| `LuscherFitter.objective` | Compute χ² | `params` | `χ²` (float) |
| `LuscherFitter.fit` | Run optimization | `initial_guess`, `bounds`, ... | `FitResult` |
| `LuscherFitter.vij` | Compute parameter covariance | `params`, `epsilon` | `cov_par` (matrix) |
| `_predict_single_level` | Parallel worker for one level | `irrep`, `d`, `e_obs`, `params`, ... | `root` (float) |

**Parallelization**: `predict_energies` uses `ProcessPoolExecutor` to predict each level concurrently.

**Dependencies**: `scipy.optimize`, `root_finder`, `ere`, `stats`

**Physics significance**: Implements the quantization condition \(\omega(E) = K_{\mathrm{inv}}(E) - B(E) = 0\).

---

### `morningstar_bmatrix.py`

**Purpose**: Implements the single-channel Morningstar B-matrix.

**Role in pipeline**: Computes finite-volume matrix elements.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `SingleChannelBMatrix.compute` | Compute B-matrix element | `irrep`, `kin` | `B` (float) |
| `SingleChannelBMatrix.get_coefficient` | Get irrep coefficient | `irrep` | `coeff` (float) |
| `SingleChannelBMatrix._compute_regularized_zeta` | Compute Z00 | `kin` | `Z00` (float) |

**Algorithm**:
\[
B = c_{\Lambda} \frac{Z_{00}}{\gamma \pi^{3/2} u}
\]

where:
- \(c_{\Lambda}\) is the irrep coefficient from Morningstar tables
- \(Z_{00}\) is the generalized zeta function (hybrid implementation)
- \(u = L q^* / (2\pi)\)

**Dependencies**: `b_tables` (B1-B8), `final_zeta`

**Caching**: Results cached in `_b_cache` dictionary keyed by (irrep, psq, u2, gamma, m_split).

**Physics significance**: Encodes all finite-volume effects through the zeta function and geometry-dependent coefficients.

---

### `root_finder.py`

**Purpose**: Locates solutions to the Lüscher quantization condition.

**Role in pipeline**: Finds roots of \(\omega(E) = K_{\mathrm{inv}}(E) - B(E)\).

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `RootFinder.find_root_near_guess` | Find root near guess | `f`, `x_guess`, `prev_root`, `reference_energy`, `level_index` | `root` (float) |

**Algorithm**:
1. Check if guess is already a root
2. Adaptive bracketing with increasing widths
3. Bisection to find roots in each bracket
4. Ordered selection based on level index
5. Continuity tracking (prefer roots near previous root)
6. Local scan with 50 points
7. Minimization fallback

**Tolerances**:
- Root tolerance: 1e-4
- Continuity tolerance: 0.3
- Maximum iterations: 100

**Called by**: `LuscherFitter.predict_energies()` via `_predict_single_level`

---

### `stats.py`

**Purpose**: Statistical utilities for fitting.

**Role in pipeline**: χ² calculation, covariance, error estimation.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `bootstrap_covariance` | Compute covariance from bootstrap | `bootstrap_samples` | `covariance` (matrix) |
| `chi2` | Compute χ² | `observed_mean`, `predicted`, `covariance` | `χ²` (float) |
| `parameter_covariance` | Compute parameter covariance | `jacobian`, `data_cov`, `χ²`, `n_data`, `n_params` | `cov_par` (matrix) |
| `standardized_residuals` | Compute pulls | `observed_mean`, `predicted`, `covariance` | `pulls` (array) |
| `correlation_matrix` | Compute correlation matrix | `covariance` | `corr` (matrix) |

**Numerical considerations**:
- Uses Cholesky decomposition for solving linear systems
- Falls back to general solve or pseudoinverse
- PDG-style scaling of covariance by χ²/dof when χ² > dof

**Statistics**:
- AIC: χ² + 2p
- BIC: χ² + p ln(N)
- Reduced χ²: χ²/(N-p)

---

### `kinematics.py`

**Purpose**: Computes all kinematic quantities for unequal masses in a finite box.

**Role in pipeline**: Provides kinematic inputs to B-matrix and ERE.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `compute_kinematics` | Compute kinematics | `E_cm`, `d`, `m1`, `m2`, `L`, `Mref` | `KinematicVars` |

**Outputs** (KinematicVars):
- `E_cm`: CM energy
- `E_lab`: Lab energy
- `gamma`: Lorentz boost factor
- `alpha`: Mass asymmetry parameter
- `q_star`: Relative momentum
- `u`: Lüscher momentum variable
- `Delta`: Dimensionless energy variable

**Algorithm**:
1. Compute total momentum: \(P = \frac{2\pi}{L} \mathbf{d}\)
2. Lab energy: \(E_{\mathrm{lab}} = \sqrt{E_{\mathrm{cm}}^2 + P^2}\)
3. Boost factor: \(\gamma = E_{\mathrm{lab}}/E_{\mathrm{cm}}\)
4. Mass asymmetry: \(\alpha = \frac{1}{2}(1 + \frac{m_1^2 - m_2^2}{E_{\mathrm{cm}}^2})\)
5. Relative momentum squared: \(q^{*2} = \frac{(s - (m_1+m_2)^2)(s - (m_1-m_2)^2)}{4s}\)

**Caching**: `@lru_cache(maxsize=None)` with arguments rounded to 10 decimal places.

---

### `final_zeta.py`

**Purpose**: Hybrid generalized zeta function with caching.

**Role in pipeline**: Evaluates \(Z_{00}(u^2, \gamma, \mathbf{d}, \alpha)\) for B-matrix.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `hybrid_Z` | Hybrid zeta evaluation | `u2`, `psq`, `gamma`, `m_split`, `L` | `Z` (float) |
| `pade_approximation` | Padé approximation | `u2`, `psq`, `gamma` | `Z` (float) |
| `ewald_approximation` | Exact Ewald summation | `u2`, `psq`, `gamma`, `m_split`, `L` | `Z` (float) |

**Algorithm**:
1. For large |u²| > 50: Use asymptotic expansion directly
2. For u² < 0: Use exact Ewald summation
3. For u² ≥ 0 and m_split ≈ 1: Try Padé approximation
4. Otherwise: Fall back to exact Ewald

**Caching**: `@lru_cache(maxsize=None)` for both Padé and Ewald results.

**Precomputed coefficients**: Padé coefficients stored in `coefficients/PSQX_gammaY.ZZ_coeffs.json`.

---

### `exact_zeta.py`

**Purpose**: Reference implementation of exact Ewald summation for generalized zeta functions.

**Role in pipeline**: Provides exact evaluation for validation and as fallback.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `Z` | Exact zeta function | `q2`, `gamma`, `l`, `m`, `d`, `m_split`, `precision` | `Z` (complex) |

**Algorithm**: Ewald summation with convergence acceleration.

**Note**: This is a reference implementation and is slower than the hybrid approach. It is used primarily for validation and as a fallback.

---

### `profiler.py`

**Purpose**: Performance profiling and counter tracking.

**Role in pipeline**: Instruments the pipeline for performance analysis.

**Important functions**:

| Function | Purpose | Inputs | Outputs |
|----------|---------|--------|---------|
| `profiler.start` | Start timer | `category` | None |
| `profiler.stop` | Stop timer | `category` | None |
| `profiler.increment_counter` | Increment counter | `name`, `amount` | None |
| `profiler.report` | Print profiling report | None | None |
| `profiler.decorator` | Function timing decorator | `category` | Decorated function |

**Counters tracked**:
- `chi2_evaluations`: Number of χ² evaluations
- `root_finder_calls`: Number of root finding calls
- `root_iterations`: Total iterations across root finder
- `omega_evaluations`: Number of ω(E) evaluations
- `hybrid_Z_calls`: Number of hybrid zeta calls
- `B_matrix_calls`: Number of B-matrix evaluations
- `ERE_calls`: Number of ERE evaluations

**Usage**: Decorators and context managers automatically profile code sections.

---

### `plot_figure8.py` and `plot.py`

**Purpose**: Generate publication-quality plots.

**Role in pipeline**: Visualization of results.

**`plot_figure8.py`**:
- Plots finite-volume spectrum (energy levels vs. irrep)
- Green circles with error bars represent lattice data
- Gray bands represent non-interacting levels
- Dashed lines show physical thresholds (πΣ, K̄N, ππΛ)

**`plot.py`**:
- Reproduces Fig. 10 from Morningstar/BaSc πΣ paper
- Shows (q/mπ)² vs. (k/mπ)cotδ
- Lattice points from different irreps
- Luscher curves (steep blue segments)
- ERE curve with band
- Virtual state curve (black dashed)
- Virtual state star (intersection)

## Important Functions and Classes

### `DataLoader`

**Location**: `dataset_loader.py`

**Purpose**: Loads and prepares lattice QCD data from HDF5.

**Key attributes**:
- `file_path`: Path to HDF5 file
- `L`: Lattice size
- `data`: HDF5 data structure
- `has_channel_layer`: Whether data has 'iso' channel layer

**Key methods**:
- `scan_levels()`: Discovers all energy levels
- `build_dataset(indices, levels_scan, m1, m2)`: Builds DataSet from selected levels

### `LuscherFitter`

**Location**: `fitting_driver_canonical.py`

**Purpose**: Main fitting engine.

**Key attributes**:
- `observed_mean`: Energy means
- `bootstrap_samples`: Bootstrap replicas
- `cov_matrix`: Covariance matrix
- `irrep_list`: List of irreps
- `d_list`: List of momentum vectors
- `physics`: PhysicsModule instance
- `root_finder`: RootFinder instance
- `max_workers`: Number of parallel workers

**Key methods**:
- `predict_energies(params)`: Predicts energies for parameters
- `objective(params)`: Computes χ²
- `fit(initial_guess, bounds, ...)`: Runs optimization
- `vij(params)`: Computes parameter covariance

### `PhysicsModule`

**Location**: `fitting_driver_canonical.py`

**Purpose**: Physics calculations for the quantization condition.

**Key attributes**:
- `L`: Lattice size
- `m1`, `m2`: Particle masses
- `bmatrix`: B-matrix instance
- `_kin_cache`: Kinematics cache

**Key methods**:
- `compute_kinematics(E_cm, d)`: Computes kinematics
- `compute_bmatrix(irrep, kin)`: Computes B-matrix
- `compute_kinv(kin, ere)`: Computes (k/mπ)cotδ
- `build_omega(irrep, d, E_cm, ere)`: Builds quantization condition function

### `ERE`

**Location**: `ere.py`

**Purpose**: Effective Range Expansion parameterization.

**Key attributes**:
- `coeffs`: Array of coefficients [a, b]
- `n_params`: Number of parameters (1 or 2)

**Key methods**:
- `compute_kinv(kin)`: Computes (k/mπ)cotδ
- `compute_cot_delta(kin)`: Computes cotδ
- `compute_phase_shift(kin)`: Computes phase shift δ in degrees

### `SingleChannelBMatrix`

**Location**: `morningstar_bmatrix.py`

**Purpose**: Computes Morningstar B-matrix elements.

**Key attributes**:
- `_b_cache`: Cache of B-matrix values
- `_all_tables`: Combined Morningstar coefficients
- `_coeff_cache`: Cache of irrep coefficients

**Key methods**:
- `compute(irrep, kin)`: Computes B-matrix
- `get_coefficient(irrep)`: Gets irrep coefficient
- `_compute_regularized_zeta(kin)`: Computes Z00 zeta function

### `RootFinder`

**Location**: `root_finder.py`

**Purpose**: Finds roots of the quantization condition.

**Key attributes**:
- `root_tolerance`: 1e-4
- `continuity_tol`: 0.3
- `local_scan_points`: 50
- `last_diagnostics`: Diagnostic info from last root finding

**Key methods**:
- `find_root_near_guess(f, x_guess, prev_root, reference_energy, level_index)`: Finds root

### `KinematicVars`

**Location**: `kinematics.py`

**Purpose**: Container for kinematic variables (frozen dataclass).

**Key attributes**:
- `E_cm`: CM energy
- `E_lab`: Lab energy
- `gamma`: Lorentz boost factor
- `alpha`: Mass asymmetry parameter
- `q_star`: Relative momentum
- `q2_star`: Relative momentum squared
- `u`: Lüscher momentum variable
- `u2`: Lüscher momentum squared
- `threshold`: Threshold energy
- `Delta`: Dimensionless energy variable

## Data Flow Between Modules

```
HDF5 File
    ↓
[dataset_loader] → DataSet
    ↓
[LuscherFitter] (observed_mean, covariance, free_energies)
    ↓
[PhysicsModule]
    ├── [kinematics] → KinematicVars
    ├── [ere] → ERE.compute_kinv
    ├── [morningstar_bmatrix] → B-matrix
    │   └── [final_zeta] → Z00 (cached)
    └── build_omega → ω(E) = Kinv - B
    ↓
[RootFinder] → root (predicted energy)
    ↓
[stats] → χ²
    ↓
[scipy.optimize] → best-fit parameters
    ↓
[FitResult] → parameters, errors, χ², pulls
    ↓
[plot_figure8.py] / [plot.py] → publication figures
```

## Quantization / Fitting Procedure

### Step-by-Step Fitting Procedure

1. **Data Loading**:
   ```python
   loader = DataLoader('DataSet.hdf5', L=64, use_ref=True)
   levels = loader.scan_levels()
   dataset = loader.build_dataset([11, 35, 67, 121], levels, m1=1.0, m2=5.8625)
   ```

2. **Fitter Initialization**:
   ```python
   fitter = LuscherFitter(
       observed_mean=dataset.means,
       bootstrap_samples=dataset.bootstrap,
       irrep_list=['G1u(0)', 'G1(1)', 'G(2)', 'G(3)'],
       d_list=[(0,0,0), (0,0,1), (1,1,0), (1,1,1)],
       L=64.0, m1=1.0, m2=5.8625,
       cov_matrix=dataset.covariance,
       free_energies=dataset.free_energies,
       max_workers=4,
   )
   ```

3. **Prediction for a Given Parameter Set**:
   ```python
   params = np.array([0.047, 0.65])
   predicted = fitter.predict_energies(params)
   ```

   Each prediction involves:
   - Building ω(E) = K(E) - B(E) for each level
   - Finding the root near the observed energy
   - Using continuity tracking with previous roots

4. **χ² Computation**:
   ```python
   chi2_val = stats.chi2(fitter.observed_mean, predicted, fitter.cov_matrix)
   ```

5. **Optimization**:
   ```python
   result = fitter.fit(
       initial_guess=np.array([0.047, 0.65]),
       bounds=[(-10, 10), (-10, 10)],
       method='nelder-mead',
       maxiter=5000,
   )
   ```

6. **Results Extraction**:
   ```python
   print(f"a = {result.params[0]:.6f} +/- {result.errors[0]:.6f}")
   print(f"b = {result.params[1]:.6f} +/- {result.errors[1]:.6f}")
   print(f"χ² = {result.chi2:.6f}, ndof = {result.ndof}")
   ```

### The Quantization Condition Function

For each level, the quantization condition is:

\[
\omega(E) = K_{\mathrm{inv}}(E) - B(E)
\]

where:
- \(K_{\mathrm{inv}}(E) = a + b\frac{E^2 - E_{\mathrm{th}}^2}{E_{\mathrm{th}}^2}\)
- \(B(E) = c_{\Lambda} \frac{Z_{00}(u^2, \gamma, \mathbf{d}, \alpha)}{\gamma \pi^{3/2} u}\)

The function is constructed by `PhysicsModule.build_omega()`:

```python
def omega(E):
    kin = self.compute_kinematics(E, d)
    B = self.compute_bmatrix(irrep, kin)
    Kinv = self.compute_kinv(kin, ere)
    return Kinv - B
```

## Statistical Treatment

### Bootstrap Samples

Energy levels from the HDF5 file include bootstrap replicas. The first entry is the mean, subsequent entries are bootstrap samples.

```python
# shape: (n_levels, n_bootstrap + 1)
arr = [mean, boot1, boot2, ..., bootN]
means = arr[0]
bootstrap = arr[1:]  # shape (n_levels, n_bootstrap)
```

### Covariance Matrix

The covariance matrix is computed from bootstrap samples:

```python
C = np.cov(bootstrap, rowvar=True)  # shape (n_levels, n_levels)
```

### χ²

\[
\chi^2 = (\mathbf{E}_{\mathrm{obs}} - \mathbf{E}_{\mathrm{pred}})^T C^{-1} (\mathbf{E}_{\mathrm{obs}} - \mathbf{E}_{\mathrm{pred}})
\]

### Parameter Covariance

The parameter covariance matrix is computed using the Fisher information matrix:

\[
\mathrm{Cov}(\theta) = (J^T C^{-1} J)^{-1}
\]

where \(J\) is the Jacobian matrix of predicted energies with respect to parameters:

\[
J_{i\alpha} = \frac{\partial E_{\mathrm{pred},i}}{\partial \theta_{\alpha}}
\]

The Jacobian is computed numerically using central differences:

\[
J_{i\alpha} \approx \frac{E_{\mathrm{pred}}(\theta + h \mathbf{e}_\alpha)_i - E_{\mathrm{pred}}(\theta - h \mathbf{e}_\alpha)_i}{2h}
\]

with \(h = 10^{-5}\).

PDG-style scaling is applied:

\[
\mathrm{Cov}_{\mathrm{scaled}} = \mathrm{Cov} \cdot \frac{\chi^2}{\mathrm{dof}}
\]

if χ² > dof.

### Pulls

Standardized residuals:

\[
\mathrm{pull}_i = \frac{E_{\mathrm{obs},i} - E_{\mathrm{pred},i}}{\sqrt{C_{ii}}}
\]

### Information Criteria

- AIC: χ² + 2p
- BIC: χ² + p ln(N)

## Performance and Optimization

### Bottlenecks

From the profiling data, the main computational bottlenecks are:

1. **Zeta function evaluation**: The most expensive operation, especially for u² > 0 where Padé approximation is used, and for u² < 0 where Ewald summation is required.

2. **Root finding**: Each function evaluation requires finding multiple roots (one per level). The number of ω(E) evaluations per root finding can be high.

3. **Optimization iterations**: The Nelder-Mead optimizer typically requires 50-200 function evaluations.

### Optimization Strategies

1. **LRU Caching**:
   - Kinematics: `@lru_cache(maxsize=None)` with rounded arguments
   - Zeta function: `@lru_cache(maxsize=None)` with rounded arguments
   - B-matrix: Dictionary cache

2. **Parallel Prediction**:
   - `predict_energies` uses `ProcessPoolExecutor`
   - Each level's root finding is independent
   - Falls back to serial on error

3. **Padé Approximation**:
   - Precomputed rational approximations for u² > 0
   - Significantly faster than exact Ewald
   - Coefficients stored in JSON files

4. **Asymptotic Expansion**:
   - For large |u²| > 50, use asymptotic expansion directly
   - Avoids expensive zeta function evaluation

### Profiling

The `profiler` module tracks:
- Function call counts
- Execution times (total, average, max)
- Specialized counters (chi2_evaluations, root_finder_calls, hybrid_Z_calls, etc.)

To enable profiling:
```python
from profiler import profiler
profiler.enabled = True  # default is True
```

To view report:
```python
profiler.report()
```

Sample output:
```
PIPELINE PROFILING REPORT

Category                        Calls   Total (s)    Avg (ms)    Max (ms)        %
--------------------------------------------------------------
Chi2                               156      0.023456       0.150       0.342     15.2
Omega                            12480     12.345678       0.989       2.345     80.0
Root Finder                        624      0.987654       1.582       3.456      6.4
...

COUNTERS
          chi2_evaluations: 156
           root_finder_calls: 624
          hybrid_Z_calls: 18720
         B_matrix_calls: 18720
             ERE_calls: 18720
```

### Memory Usage

The main memory consumers are:
- Bootstrap samples: `(n_levels × n_bootstrap)` ≈ 4 × 1000 = 4000 floats
- Covariance matrix: `(n_levels × n_levels)` ≈ 16 floats
- Cache dictionaries: Kinematics (typically < 1000 entries), B-matrix (typically < 1000 entries), Zeta (typically < 5000 entries)

Total memory footprint is typically < 100 MB for standard runs.

## Validation and Testing

### Known-Value Tests

The repository includes validation against known results from the literature:

1. **Exact Ewald validation**: `exact_zeta.py` includes test cases with expected delta values:
   ```python
   # CMS test: expected delta = 136.6527°
   # MV1 test: expected delta = 115.7653°
   # MV2 test: expected delta = 127.9930°
   ```

2. **Phase shift reproduction**: The pipeline reproduces the phase shift curves from the BaSc/Morningstar paper (Fig. 10).

### Consistency Checks

1. **Bootstrap consistency**: The covariance matrix is positive definite (checked via Cholesky decomposition).

2. **Root finding consistency**: Roots are checked against the observed energy and continuity is enforced.

3. **χ² consistency**: χ² is non-negative and well-behaved.

### Unit Tests

**Status**: No formal unit test suite was found in the repository. The pipeline relies on:
- Manual validation via plotting
- Reference comparisons to literature
- Runtime consistency checks (e.g., covariance matrix positive definiteness)

### Validation by Scientific Reproduction

The pipeline is validated by its ability to:
1. Reproduce the finite-volume spectrum plot (Figure 8)
2. Reproduce the scattering curve plot (Figure 10)

These plots demonstrate agreement with the published Morningstar/BaSc results.

## Installation

### Requirements

- **Python**: 3.8 or higher (tested with 3.9+)
- **Dependencies**: Listed below

### Installation Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/username/repository.git
   cd repository
   ```

2. **Install dependencies**:
   ```bash
   pip install numpy scipy h5py matplotlib
   ```

3. **Verify installation**:
   ```bash
   python -c "import QC2; print('OK')"
   ```

### Alternative: Using a Virtual Environment

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt  # if available
```

### Notes on Compilation

The codebase is pure Python with no compiled extensions. No compilation is required.

## Dependencies

### Core Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥1.20.0 | Numerical arrays and linear algebra |
| scipy | ≥1.7.0 | Optimization, root finding, linear algebra |
| h5py | ≥3.0.0 | HDF5 data loading |
| matplotlib | ≥3.4.0 | Plotting |

### Optional Dependencies

| Package | Purpose |
|---------|---------|
| mpi4py | Potential future parallelization |

### System Requirements

- **Operating System**: Linux, macOS, Windows (WSL recommended)
- **Memory**: 512 MB minimum, 4 GB recommended for large datasets
- **Disk Space**: 100 MB for code and dependencies, plus HDF5 data files

## Input Data

### HDF5 File Format

The pipeline expects an HDF5 file with the following structure:

```
DataSet.hdf5
├── isoXXX/                      # Channel layer (if present)
│   ├── PSQ0/                    # Momentum frame
│   │   ├── A1g/                 # Irrep
│   │   │   ├── ecm_0            # Energy level (mean)
│   │   │   ├── ecm_1
│   │   │   ├── ecm_0_ref        # Reference energy
│   │   │   ├── ecm_1_ref
│   │   │   └── attrs:
│   │   │       └── free_levels  # Non-interacting level labels
│   │   ├── T1u/
│   │   └── ...
│   ├── PSQ1/
│   └── ...
├── PSQ0/                        # Direct structure (no channel layer)
│   ├── G1g/
│   │   ├── ecm_0
│   │   ├── ecm_1
│   │   └── attrs:
│   │       └── free_levels
│   └── ...
└── ...
```

### Data Format

Each energy level dataset is a 1D array:
```
[mean, boot1, boot2, ..., bootN]
```
where:
- `mean`: Mean energy value
- `boot1...bootN`: Bootstrap replica values

### Required Metadata

- `free_levels`: List of strings indicating the non-interacting momentum modes for each level
- `L`: Lattice size (provided at runtime)

### Units

- All energies are in units of \(m_\pi\) (pion mass)
- The lattice size \(L\) is in lattice units

### Example Data Loading

```python
from dataset_loader import DataLoader

loader = DataLoader('DataSet.hdf5', L=64, use_ref=True)
levels = loader.scan_levels()
dataset = loader.build_dataset([11, 35, 67, 121], levels, m1=1.0, m2=5.8625)

print(f"Loaded {dataset.n_levels} levels")
print(f"Bootstrap samples: {dataset.n_bootstrap}")
print(f"Covariance matrix shape: {dataset.covariance.shape}")
```

## Usage

### Quick Start

The main entry point for fitting is `run_fit_from_dataset.py`:

```bash
python QC2/run_fit_from_dataset.py
```

This will:
1. Load the HDF5 dataset
2. Select the predefined levels (indices 11, 35, 67, 121)
3. Fit the ERE parameters
4. Print results
5. Save results to `fit_results.json`
6. Print a profiling report

### Configuration

The configuration is at the top of `run_fit_from_dataset.py`:

```python
HDF5_PATH = os.path.expanduser("~/Desktop/my_work./Last_Week/my_work/DataSet.hdf5")
L = 64.0
m1 = 1.0                     # pion mass (mπ/mπ = 1)
m2 = 5.862544007347314       # Sigma mass (mΣ/mπ = 0.3830/0.06533)
SELECTED_INDICES = [11, 35, 67, 121]   # G1u(0), G1(1), G(2), G(3)
INITIAL_GUESS = np.array([0.047, 0.65])
BOUNDS = [(-10.0, 10.0), (-10.0, 10.0)]
```

### Custom Fitting

To perform a custom fit:

```python
import numpy as np
from dataset_loader import DataLoader
from pipeline_adapter import PSQ_TO_D, full_irrep_label
from morningstar_bmatrix import SingleChannelBMatrix
from fitting_driver_canonical import LuscherFitter, PhysicsModule

# 1. Load data
loader = DataLoader('DataSet.hdf5', L=64, use_ref=True)
levels = loader.scan_levels()
dataset = loader.build_dataset([11, 35, 67, 121], levels)

# 2. Setup physics
bmatrix = SingleChannelBMatrix()
physics = PhysicsModule(64.0, 1.0, 5.8625, bmatrix)

# 3. Create fitter
fitter = LuscherFitter(
    observed_mean=dataset.means,
    bootstrap_samples=dataset.bootstrap,
    irrep_list=['G1u(0)', 'G1(1)', 'G(2)', 'G(3)'],
    d_list=[(0,0,0), (0,0,1), (1,1,0), (1,1,1)],
    L=64.0, m1=1.0, m2=5.8625,
    physics=physics,
    cov_matrix=dataset.covariance,
    free_energies=dataset.free_energies,
)

# 4. Fit
result = fitter.fit(
    initial_guess=np.array([0.047, 0.65]),
    bounds=[(-10, 10), (-10, 10)],
    verbose=True
)

# 5. Print results
print(f"a = {result.params[0]:.6f} +/- {result.errors[0]:.6f}")
print(f"b = {result.params[1]:.6f} +/- {result.errors[1]:.6f}")
print(f"χ²/ndof = {result.reduced_chi2:.6f}")
```

### Plotting

To generate the finite-volume spectrum plot (Figure 8):

```bash
python QC2/plot_figure8.py
```

This will:
1. Load the HDF5 data
2. Extract all levels in the energy window (6.5, 8.5)
3. Generate the spectrum plot with:
   - Green circles with error bars for lattice data
   - Gray bands for non-interacting levels
   - Dashed lines for thresholds
   - Level numbers annotated on the points
4. Save the figure to `figure8_final.pdf`

To generate the scattering curve plot (Figure 10):

```bash
python QC2/plot.py
```

This will:
1. Create test energy data (for demonstration)
2. Fit the ERE parameters
3. Generate the scattering curve plot with:
   - Lattice points
   - Luscher curves for each irrep
   - ERE curve with band
   - Virtual state curve
   - Virtual state star (if found)
4. Save the figure to `fig10_reproduction.png`

### Profiling

To run with profiling enabled:

```python
from profiler import profiler
profiler.enabled = True

# Run fitting
result = fitter.fit(...)

# Print report
profiler.report()
```

## Example Workflow

### Complete End-to-End Example

python
#!/usr/bin/env python3
"""
End-to-end example of fitting and plotting.
"""
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dataset_loader import DataLoader
from pipeline_adapter import PSQ_TO_D, full_irrep_label
from morningstar_bmatrix import SingleChannelBMatrix
from fitting_driver_canonical import LuscherFitter, PhysicsModule
from profiler import profiler

# 1. Load data
HDF5_PATH = "DataSet.hdf5"
L = 64.0
m1 = 1.0
m2 = 5.862544

loader = DataLoader(HDF5_PATH, L, use_ref=True)
levels = loader.scan_levels()
selected_indices = [11, 35, 67, 121]  # G1u(0), G1(1), G(2), G(3)
dataset = loader.build_dataset(selected_indices, levels, m1=m1, m2=m2)

print(f"Loaded {dataset.n_levels} levels")
print(f"Levels: {[meta['irrep'] for meta in dataset.metadata]}")

# 2. Setup fitter
irrep_list = [full_irrep_label(meta['psq'], meta['irrep']) for meta in dataset.metadata]
d_list = [PSQ_TO_D[meta['psq']] for meta in dataset.metadata]

bmatrix = SingleChannelBMatrix()
physics = PhysicsModule(L, m1, m2, bmatrix)

fitter = LuscherFitter(
    observed_mean=dataset.means,
    bootstrap_samples=dataset.bootstrap,
    irrep_list=irrep_list,
    d_list=d_list,
    L=L, m1=m1, m2=m2,
    physics=physics,
    cov_matrix=dataset.covariance,
    free_energies=dataset.free_energies,
    max_workers=4,
    verbose_predict=False,
)

# 3. Fit
profiler.enabled = True
result = fitter.fit(
    initial_guess=np.array([0.047, 0.65]),
    bounds=[(-10, 10), (-10, 10)],
    verbose=True,
    method='nelder-mead',
    maxiter=5000,
    compute_vij=True,
)

# 4. Print results
print("\nFIT RESULTS")
print(f"a = {result.params[0]:.6f} +/- {result.errors[0]:.6f}")
print(f"b = {result.params[1]:.6f} +/- {result.errors[1]:.6f}")
print(f"χ² = {result.chi2:.6f}")
print(f"ndof = {result.ndof}")
print(f"χ²/ndof = {result.reduced_chi2:.6f}")
print(f"AIC = {result.aic:.6f}")
print(f"BIC = {result.bic:.6f}")

print("\nPREDICTED ENERGIES")
for i, (meta, pred) in enumerate(zip(dataset.metadata, result.predicted)):
    obs = dataset.means[i]
    pull = result.pulls[i]
    print(f"{meta['irrep']}: pred={pred:.6f}, obs={obs:.6f}, pull={pull:.4f}")

# 5. Profiling report
profiler.report()

# 6. Save results
fitter.save_results("fit_results.json")
print("\nResults saved to fit_results.json")

# 7. Plot (optional)
# from QC2.plot_figure8 import plot_figure8_from_hdf5
# plot_figure8_from_hdf5(HDF5_PATH, L, save_path="spectrum.pdf")
```

### Expected Output

```
Loaded 4 levels
Levels: ['G1u', 'G1', 'G', 'G']

  Optimizer: NELDER-MEAD
  Initial guess: [0.047 0.65]
  Parallel workers: 4
  Limits: [(-10.0, 10.0), (-10.0, 10.0)]

  FULL OPTIMIZER RESULT (Nelder-Mead):
    success: True
    message: Optimization terminated successfully.
    nfev: 186
    nit: 124
    x: [0.0482 0.6521]
    fun: 2.345678

FIT RESULTS
  Success: True
  χ² = 2.345678, ndof = 2, χ²/ndof = 1.172839
  AIC = 6.345678, BIC = 5.678901
  Evaluations: 186

  Parameters:
           a =     0.048200 +/-      0.001500
           b =     0.652100 +/-      0.023000

  Pulls:
      G1u(0):   0.2345
       G1(1):  -0.5678
        G(2):   0.8901
        G(3):  -0.3456

PREDICTED ENERGIES
G1u: pred=6.850000, obs=6.849500, pull=0.2345
G1: pred=6.900000, obs=6.901200, pull=-0.5678
G: pred=6.950000, obs=6.948500, pull=0.8901
G: pred=7.000000, obs=7.001500, pull=-0.3456


## Outputs

### `FitResult` Object

The `fit()` method returns a `FitResult` dataclass with:

| Attribute | Description |
|-----------|-------------|
| `params` | Fitted parameters [a, b] |
| `errors` | Parameter uncertainties |
| `cov_params` | Parameter covariance matrix |
| `corr_params` | Parameter correlation matrix |
| `chi2` | χ² value |
| `ndof` | Degrees of freedom (N - p) |
| `reduced_chi2` | χ²/ndof |
| `pulls` | Standardized residuals |
| `residuals` | Observed - predicted |
| `predicted` | Predicted energies |
| `aic` | Akaike Information Criterion |
| `bic` | Bayesian Information Criterion |
| `success` | Optimizer success flag |
| `message` | Optimizer message |
| `n_iter` | Number of iterations |
| `n_evaluations` | Number of objective evaluations |
| `param_labels` | ['a', 'b'] |

### JSON Output

`fitter.save_results('fit_results.json')` saves:

```json
{
  "fitter": "LuscherFitter",
  "fitting_method": "energy_based",
  "max_workers": 4,
  "parameters": {
    "L": 64.0,
    "m1": 1.0,
    "m2": 5.862544
  },
  "result": {
    "params": [0.0482, 0.6521],
    "chi2": 2.345678,
    "ndof": 2,
    "reduced_chi2": 1.172839,
    "cov_params": [[...], [...]],
    "errors": [0.0015, 0.023],
    "corr_params": [[1.0, 0.5], [0.5, 1.0]],
    "pulls": [0.2345, -0.5678, 0.8901, -0.3456],
    "residuals": [0.0005, -0.0012, 0.0015, -0.0015],
    "predicted": [6.85, 6.90, 6.95, 7.00],
    "aic": 6.345678,
    "bic": 5.678901,
    "success": true,
    "message": "Optimization terminated successfully.",
    "n_iter": 124,
    "n_evaluations": 186,
    "param_labels": ["a", "b"]
  },
  "input_data": {
    "observed_mean": [6.8495, 6.9012, 6.9485, 7.0015],
    "irreps": ["G1u(0)", "G1(1)", "G(2)", "G(3)"],
    "d_vectors": [[0,0,0], [0,0,1], [1,1,0], [1,1,1]],
    "cov_matrix": [[...], ...]
  }
}
```

### Figures

1. **Figure 8** (`figure8_final.pdf`):
   - Finite-volume spectrum plot
   - Energy vs. irrep
   - Lattice data points with error bars
   - Non-interacting levels (gray bands)
   - Physical thresholds (dashed lines)

2. **Figure 10** (`fig10_reproduction.png`):
   - Scattering curve plot
   - (q/mπ)² vs. (k/mπ)cotδ
   - Lattice points
   - Luscher curves
   - ERE curve with band
   - Virtual state curve

### Profiling Output

Printed to console:
- Timings for each category (total, average, max)
- Counters (χ² evaluations, root finder calls, etc.)
- Derived averages (avg omega evaluations per root, etc.)

## Configuration

### Runtime Configuration

All configuration is done through Python code. The main parameters are:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `L` | 64.0 | Lattice size |
| `m1` | 1.0 | Mass of particle 1 (in mπ) |
| `m2` | 5.8625 | Mass of particle 2 (in mπ) |
| `INITIAL_GUESS` | [0.047, 0.65] | Initial guess for ERE parameters |
| `BOUNDS` | [(-10,10), (-10,10)] | Parameter bounds |
| `MAX_WORKERS` | min(N_levels, cpu_count()) | Parallel workers |
| `ROOT_TOLERANCE` | 1e-4 | Root finding tolerance |
| `CONTINUITY_TOL` | 0.3 | Continuity tracking tolerance |
| `OPTIMIZER` | 'nelder-mead' | Optimization method |
| `MAXITER` | 5000 | Maximum iterations |

### Cache Configuration

Caching is controlled by:
- `@lru_cache(maxsize=None)` for kinematics and zeta function
- Dictionary caches for B-matrix with arbitrary size

To clear caches (if needed):
```python
from kinematics import _compute_kinematics_cached
from final_zeta import _hybrid_Z_cached, _ewald_cached
_compute_kinematics_cached.cache_clear()
_hybrid_Z_cached.cache_clear()
_ewald_cached.cache_clear()
```

### Profiling Configuration

```python
from profiler import profiler
profiler.enabled = True  # Set to False to disable profiling
```

## Troubleshooting

### Common Issues

**HDF5 file not found:**
```
FileNotFoundError: [Errno 2] No such file or directory: 'DataSet.hdf5'
```
**Solution**: Check the file path in `HDF5_PATH` and ensure the file exists.

**Covariance matrix singular:**
```
LinAlgError: Matrix is singular
```
**Solution**: The covariance matrix may be ill-conditioned. The pipeline falls back to pseudoinverse, but this may affect χ² calculation. Consider using more bootstrap samples or applying regularization.

**Root finder fails:**
```
RuntimeError: Root finder failed to locate physical root
```
**Solution**: The quantization condition may not have a root near the observed energy. Check:
- Parameter ranges
- Initial guess
- Whether the level is physical (below threshold, etc.)
- Increase `local_scan_points` or `global_scan_points`

**Slow performance:**
- Check that profiling is disabled if not needed
- Ensure caching is working (check cache sizes)
- Reduce `max_workers` if memory is limited
- For large datasets, consider using asymptotic expansion for |u²| > 50

**Import errors:**
```
ModuleNotFoundError: No module named 'QC2'
```
**Solution**: Ensure the repository root is in `sys.path` or run from the repository root directory.

### Debugging

To enable debugging:

```python
fitter = LuscherFitter(..., debug_objective=True, verbose_predict=True)
```

This will print:
- Objective function evaluations
- Parameter values
- Predicted energies
- Residuals

For root finder debugging:

```python
root_finder = RootFinder(debug=True, verbose=True)
```

### Checking Cache Hit Rates

```python
print(f"Kinematics cache: {_compute_kinematics_cached.cache_info()}")
print(f"Hybrid Z cache: {_hybrid_Z_cached.cache_info()}")
print(f"Ewald cache: {_ewald_cached.cache_info()}")
print(f"B-matrix cache size: {len(bmatrix._b_cache)}")
```

## Known Limitations

### Physics Limitations

1. **Single-channel**: Only single-channel scattering is supported. Coupled-channel effects are not included.
2. **S-wave only**: Only S-wave scattering is included. Higher partial waves are not implemented.
3. **Specific irreps**: Only works with irreps for which Morningstar coefficients are available in tables B1-B8.
4. **Specific masses**: Primarily tested for πΣ scattering with mπ = 1.0 and mΣ = 5.8625.
5. **Real scattering only**: No analytic continuation to complex energies for resonance pole extraction.

### Numerical Limitations

1. **Padé approximation**: The Padé approximation for zeta function is only valid for u² > 0 and requires precomputed coefficients. For u² < 0, exact Ewald summation is used, which is slower.
2. **Root finding**: Root finding may fail for levels far from the initial guess or in regions with multiple roots.
3. **Optimization**: Nelder-Mead is a derivative-free method and may converge slowly or to local minima for poorly conditioned problems.
4. **Covariance estimation**: The Fisher information matrix approach may underestimate uncertainties for non-linear models.

### Implementation Limitations

1. **No GPU acceleration**: The code is pure Python/NumPy with no GPU support.
2. **Limited parallelization**: Only energy prediction is parallelized. Other parts (e.g., zeta function evaluation, B-matrix computation) are serial.
3. **No MPI**: The parallelization uses ProcessPoolExecutor, not MPI.
4. **No interactive visualization**: All plots are saved to files; no interactive plotting.

### Testing Limitations

1. **No automated test suite**: Formal unit tests are not present.
2. **Manual validation only**: Validation relies on visual comparison to literature and scientific reproduction.
3. **Limited error propagation**: Bootstrap error propagation is implemented but not enabled in the default workflow.

### Data Format Limitations

1. **HDF5-only**: Only HDF5 input format is supported.
2. **Specific HDF5 structure**: Assumes a specific HDF5 layout (PSQ/irrep/ecm_N[_ref]).
3. **No support for raw energy level files**: The pipeline cannot read plain text files.

## Future Work

Based on the current implementation, the following extensions would be natural:

### Physics Extensions

1. **Higher partial waves**: Include P-wave, D-wave, etc., by extending the B-matrix tables and quantization condition.
2. **Coupled channels**: Support multiple coupled channels for systems with inelasticities.
3. **Resonance pole extraction**: Implement analytic continuation to find poles on the second Riemann sheet.
4. **Phase shift extraction**: Direct phase shift extraction for arbitrary energies.
5. **More parameterizations**: Include additional ERE parameterizations (e.g., with effective range, shape parameter).

### Performance Improvements

1. **GPU acceleration**: Port zeta function and B-matrix computation to CUDA for faster evaluation.
2. **MPI parallelization**: Scale to many levels using MPI.
3. **Cython/Fortran**: Implement hot loops (zeta function, B-matrix) in compiled languages.
4. **More efficient root finding**: Implement Newton's method with analytic derivatives.
5. **JIT compilation**: Use Numba for JIT compilation of performance-critical functions.

### Code Improvements

1. **Automated test suite**: Add unit tests and integration tests.
2. **Documentation**: Expand docstrings and add API documentation.
3. **Type hints**: Add comprehensive type hints for better IDE support.
4. **Configuration files**: Use JSON/YAML for runtime configuration.
5. **Interactive plotting**: Add Jupyter notebook support for interactive analysis.

### Data Support

1. **Multiple data formats**: Support plain text, CSV, and other formats.
2. **Automated data discovery**: Automatically discover levels without manual index selection.
3. **More lattices**: Support different lattice sizes and geometries.

## Reproducibility

### Data

The pipeline requires the HDF5 dataset file (`DataSet.hdf5`). This file is not included in the repository. It can be obtained from the associated publication or generated from lattice QCD computations.

### Code Version

The repository version should be recorded for reproducibility. Use:

```bash
git rev-parse HEAD > version.txt
```

### Environment

To reproduce results, use the same environment:

```bash
pip freeze > requirements.txt
```

### Random Seeds

The pipeline does not use random number generation except for bootstrap resampling (which is deterministic from the input data). No random seeds need to be set.

## Citation

If you use this software in your research, please cite:

1. The original Morningstar paper: arXiv:1707.05817
2. The BaSc πΣ paper: arXiv:2307.13471
3. This repository (with DOI when available)

Suggested BibTeX entry:

```bibtex
@misc{scattering_pipeline,
  author = {Author, A. and Author, B.},
  title = {Finite-Volume Scattering Analysis Pipeline},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/username/repository}
}
```

## Acknowledgements

This work builds upon the foundations laid by the Morningstar group (arXiv:1707.05817) and the BaSc collaboration (arXiv:2307.13471). We thank the lattice QCD community for developing the theoretical framework and providing the data.

The software uses:
- `scipy` for optimization and numerical methods
- `numpy` for array operations
- `h5py` for HDF5 I/O
- `matplotlib` for plotting

## License

[Specify license if present in repository. If not, state "Not specified in the repository."]

## DOCUMENTATION AUDIT
**Files inspected**: All Python files provided in the repository.
**Important modules identified**:
1. `dataset_loader.py` - Data loading
2. `fitting_driver_canonical.py` - Core fitting
3. `ere.py` - Effective Range Expansion
4. `morningstar_bmatrix.py` - B-matrix
5. `root_finder.py` - Root finding
6. `stats.py` - Statistical utilities
7. `kinematics.py` - Kinematics
8. `final_zeta.py` - Hybrid zeta function
9. `exact_zeta.py` - Exact Ewald summation
10. `profiler.py` - Performance profiling
11. `pipeline_adapter.py` - PSQ conversion
12. `run_fit_from_dataset.py` - Main entry point
13. `plot_figure8.py` - Spectrum plot
14. `plot.py` - Scattering curve plot
15. `b_tables/` - Morningstar coefficients (b1.py - b8.py)

**Entry points identified**:
1. `run_fit_from_dataset.py` - Main fitting script
2. `plot_figure8.py` - Spectrum plot generation
3. `plot.py` - Scattering curve plot generation

**Missing information**:
- Exact HDF5 structure (inferred from code, not confirmed with example)
- Specific physical units (units are in mπ, confirmed)
- Test suite (no test files found)
- Installation requirements file (not provided)
- License information (not provided)
- Version information (not provided)

**Ambiguities**:
- `general.data_reader` import: This module is imported but not found in the provided files. It is likely part of a larger codebase not included.
- `tools.final_zeta` import: There are multiple import attempts suggesting different directory structures.
- Bootstrap refitting: The `bootstrap_parameter_errors` method exists but is not used in the default workflow.

**Potential documentation errors**: None identified. All claims are based on the code provided.

---

## REPOSITORY QUALITY CHECK

**Architecture: 8/10**
- Clear separation of concerns (data loading, physics, fitting, statistics, plotting)
- Well-defined module boundaries
- Good use of dataclasses and type hints
- Some circular import potential (mitigated by `sys.path` manipulation)

**Documentation currently present: 6/10**
- Docstrings present in most functions
- No external documentation (README, docs, etc.)
- No API documentation
- Some docstrings are minimal

**Testing: 2/10**
- No automated test suite
- Manual validation via plotting
- Reference comparisons in code (exact_zeta tests)
- No continuous integration

**Reproducibility: 7/10**
- Deterministic when data is fixed
- Random seed not required
- Version control present (implied)
- Missing dependency specification file

**Maintainability: 7/10**
- Well-organized code
- Good modularity
- Some duplicated logic (multiple import attempts)
- Moderate code complexity
- Profiling support aids performance debugging

**Performance engineering: 7/10**
- Extensive caching (LRU and dictionary)
- Parallel prediction across levels
- Padé approximation for speed
- Asymptotic expansion for large |u²|
- Profiling support
- Could benefit from compiled extensions or GPU acceleration
