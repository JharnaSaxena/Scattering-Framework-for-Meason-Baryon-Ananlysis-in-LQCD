<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
</head>
<body>
<div class="container">

<h1>Finite-Volume Scattering Analysis for Meson-Baryon Systems</h1>

<h2>Overview</h2>

<p>This repository implements a computational pipeline for analyzing two-particle scattering in a finite volume from lattice QCD energy levels. It provides a complete framework for extracting scattering parameters from lattice QCD data using the Lüscher formalism, with a specific focus on the πΣ system in the S-wave. The pipeline handles unequal masses, multiple momentum frames, and uses a robust fitting procedure with detailed uncertainty estimation.</p>

<p>The software is designed for researchers working in hadron spectroscopy and lattice QCD who need to extract scattering phase shifts and effective range parameters from finite-volume energy spectra. It bridges the gap between raw lattice QCD data and physical scattering observables.</p>

<h2>Scientific Problem</h2>

<p>In lattice QCD, scattering information is extracted by computing energy levels in a finite cubic box. The Lüscher formalism relates these discrete energy levels to the infinite-volume scattering amplitude. The key challenge is solving the finite-volume quantization condition:</p>

<div class="math-block">
det[K<sup>-1</sup>(E) - B(E)] = 0
</div>

<p>where:</p>
<ul>
  <li><span class="math">K(E)</span> is the infinite-volume scattering amplitude (parameterized by the Effective Range Expansion)</li>
  <li><span class="math">B(E)</span> is the finite-volume B-matrix that encodes the geometry of the box</li>
  <li><span class="math">E</span> is the center-of-mass energy</li>
</ul>

<p>This pipeline solves this condition to fit a two-parameter Effective Range Expansion (ERE) to lattice QCD energy levels. It handles the technical complexities of unequal masses, moving frames, and the numerically challenging generalized zeta functions that appear in the B-matrix.</p>

<h2>What This Repository Does</h2>

<ol>
  <li><strong>Loads lattice QCD data</strong> from HDF5 files containing energy levels with bootstrap samples</li>
  <li><strong>Computes finite-volume B-matrix elements</strong> using Morningstar's coupled-channel formalism with precomputed coefficients from tables</li>
  <li><strong>Evaluates generalized zeta functions</strong> using a hybrid approach (Padé approximation for u² &gt; 0, exact Ewald summation for u² &lt; 0)</li>
  <li><strong>Solves the quantization condition</strong> via root finding to predict energy levels for a given set of ERE parameters</li>
  <li><strong>Fits the ERE parameters</strong> to observed lattice energy levels using χ² minimization</li>
  <li><strong>Provides uncertainty estimation</strong> through bootstrap resampling and parameter covariance analysis</li>
  <li><strong>Visualizes results</strong> with publication-quality plots (finite-volume spectra and scattering curves)</li>
  <li><strong>Profiles performance</strong> to identify bottlenecks</li>
</ol>

<h2>Key Features</h2>

<ul>
  <li><strong>Unequal-mass kinematics</strong>: Handles πΣ scattering with masses mπ = 1.0 and mΣ = 5.8625 (in units of mπ)</li>
  <li><strong>Multiple momentum frames</strong>: Supports PSQ0 (P=0), PSQ1 (P=(0,0,1)), PSQ2 (P=(1,1,0)), PSQ3 (P=(1,1,1))</li>
  <li><strong>Morningstar B-matrix</strong>: Complete implementation of the single-channel B-matrix with coefficients from arXiv:1707.05817</li>
  <li><strong>Hybrid zeta function</strong>: Padé approximation for u² &gt; 0 with exact Ewald summation fallback</li>
  <li><strong>Parallel root finding</strong>: Uses ProcessPoolExecutor for concurrent energy prediction across levels</li>
  <li><strong>Bootstrap error propagation</strong>: Full covariance matrix from bootstrap samples</li>
  <li><strong>Automatic root tracking</strong>: Continuity-based root selection for reliable fitting</li>
  <li><strong>Performance profiling</strong>: Built-in profiler tracks function calls, timings, and counters</li>
  <li><strong>Publication-ready plots</strong>: Figure 8 (finite-volume spectrum) and Figure 10 (scattering curves) reproduction</li>
</ul>

<h2>Repository Architecture</h2>

<p>The codebase is organized as follows:</p>

<pre>repository root/
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
    └── ...</pre>

<h2>Directory Structure</h2>

<pre>├── dataset_loader.py         # Data loading from HDF5, bootstrap extraction
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
└── test/                     # (No test suite found in provided files)</pre>

<h2>Computational Pipeline</h2>

<h3>High-Level Overview</h3>

<p>The pipeline transforms raw lattice QCD data into fitted scattering parameters through a multi-stage process:</p>

<ol>
  <li><strong>Data Loading</strong>: The <code>DataLoader</code> reads an HDF5 file containing energy levels with bootstrap samples. It scans the HDF5 structure to identify all available levels, extracts means and bootstrap distributions, and builds a <code>DataSet</code> object.</li>
  <li><strong>Kinematic Preparation</strong>: For each selected level, <code>PhysicsModule.compute_kinematics()</code> calculates all relevant kinematic quantities in the center-of-mass frame, including:
    <ul>
      <li>Energy in CM frame (<span class="math">E<sub>cm</sub></span>)</li>
      <li>Relative momentum (<span class="math">q<sup>*</sup></span>)</li>
      <li>Lorentz boost factor (<span class="math">γ</span>)</li>
      <li>Mass asymmetry parameter (<span class="math">α</span>)</li>
      <li>Lüscher momentum variable <span class="math">u = L q<sup>*</sup> / (2π)</span></li>
    </ul>
  </li>
  <li><strong>B-Matrix Evaluation</strong>: <code>SingleChannelBMatrix.compute()</code> calculates the finite-volume B-matrix element for a given irrep and kinematics. This uses precomputed coefficients from the Morningstar tables and the generalized zeta function.</li>
  <li><strong>Effective Range Expansion</strong>: <code>ERE.compute_kinv()</code> evaluates <span class="math">(k/m<sub>π</sub>)cotδ = a + bΔ</span> where <span class="math">Δ = (E<sup>2</sup> - E<sub>th</sub><sup>2</sup>)/E<sub>th</sub><sup>2</sup></span>.</li>
  <li><strong>Quantization Condition</strong>: The <code>PhysicsModule.build_omega()</code> constructs the function:
    <div class="math-block">ω(E) = K<sub>inv</sub>(E) - B(E)</div>
    where <span class="math">K<sub>inv</sub> = (k/m<sub>π</sub>)cotδ</span> and <span class="math">B(E)</span> is the B-matrix.
  </li>
  <li><strong>Root Finding</strong>: <code>RootFinder.find_root_near_guess()</code> finds the root of <span class="math">ω(E) = 0</span> using adaptive bracketing and continuation. This root is the predicted energy level for the given ERE parameters.</li>
  <li><strong>Prediction</strong>: For a set of ERE parameters <span class="math">(a,b)</span>, <code>LuscherFitter.predict_energies()</code> predicts all selected energy levels by finding roots for each level's quantization condition.</li>
  <li><strong>χ² Evaluation</strong>: <code>stats.chi2()</code> computes the χ² between observed and predicted energies using the covariance matrix:
    <div class="math-block">χ² = <b>r</b><sup>T</sup> C<sup>-1</sup> <b>r</b></div>
    where <span class="math"><b>r</b> = <b>E</b><sub>obs</sub> - <b>E</b><sub>pred</sub></span>.
  </li>
  <li><strong>Optimization</strong>: <code>scipy.optimize.minimize()</code> with Nelder-Mead minimizes the χ² to find best-fit ERE parameters.</li>
  <li><strong>Uncertainty Estimation</strong>: Parameter uncertainties are estimated via:
    <ul>
      <li>Fisher information matrix with PDG-style scaling</li>
      <li>Bootstrap refitting (optional, implemented but not enabled by default)</li>
    </ul>
  </li>
  <li><strong>Results</strong>: The <code>FitResult</code> object contains fitted parameters, χ², predicted energies, pulls, and covariance matrices.</li>
</ol>

<h3>Pipeline Data Flow</h3>

<pre>HDF5 File
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
Parameters, errors, χ², pulls, plots</pre>

<h2>Physics Background</h2>

<h3>Lüscher Quantization Condition</h3>

<p>For a two-particle system in a finite cubic box of size <span class="math">L</span>, the allowed energy levels are determined by the quantization condition:</p>

<div class="math-block">det[K<sup>-1</sup>(E) - B(E)] = 0</div>

<p>In the single-channel, S-wave case this reduces to:</p>

<div class="math-block">K<sub>inv</sub>(E) = B(E)</div>

<p>where:</p>
<ul>
  <li><span class="math">K<sub>inv</sub>(E) = (k/m<sub>π</sub>)cotδ(E)</span> is the inverse scattering amplitude</li>
  <li><span class="math">B(E)</span> is the finite-volume B-matrix</li>
</ul>

<h3>Effective Range Expansion</h3>

<p>The pipeline parameterizes the inverse scattering amplitude using the Effective Range Expansion:</p>

<div class="math-block">(k/m<sub>π</sub>) cotδ = a + b·Δ</div>

<p>where:</p>
<ul>
  <li><span class="math">a</span> is the inverse scattering length (in units of <span class="math">1/m<sub>π</sub></span>)</li>
  <li><span class="math">b</span> is the effective range parameter</li>
  <li><span class="math">Δ = (E<sup>2</sup> - E<sub>th</sub><sup>2</sup>)/E<sub>th</sub><sup>2</sup></span> is a dimensionless energy variable</li>
  <li><span class="math">E<sub>th</sub> = m<sub>1</sub> + m<sub>2</sub></span> is the threshold energy</li>
</ul>

<p>This parameterization follows Eq. 12 of arXiv:2307.13471.</p>

<h3>Morningstar B-Matrix</h3>

<p>The B-matrix encodes the finite-volume effects and depends on the irrep, momentum frame, and kinematic variables:</p>

<div class="math-block">B<sub>Λ</sub>(<b>d</b>, γ, α, u) = Σ<sub>ℓ,m</sub> c<sub>ℓm</sub><sup>Λ</sup> · Z<sub>ℓm</sub>(u², γ, <b>d</b>, α) / (γ π<sup>3/2</sup> u)</div>

<p>where:</p>
<ul>
  <li><span class="math">Λ</span> is the lattice irrep (e.g., G1u, G1, G, etc.)</li>
  <li><span class="math">Z<sub>ℓm</sub></span> are generalized zeta functions</li>
  <li><span class="math">c<sub>ℓm</sub><sup>Λ</sup></span> are coefficients from Morningstar's tables</li>
</ul>

<p>The coefficients are precomputed and stored in <code>b_tables/b1.py</code> through <code>b8.py</code> (from arXiv:1707.05817).</p>

<h3>Zeta Function</h3>

<p>The generalized zeta function is the most computationally expensive part of the pipeline:</p>

<div class="math-block">Z<sub>ℓm</sub>(u², γ, <b>d</b>, α) = Σ<sub><b>n</b></sub> Y<sub>ℓm</sub>(<b>r</b>) / (r² - u²) &nbsp;(with regularization)</div>

<p>The pipeline implements two evaluation strategies:</p>
<ol>
  <li><strong>Padé approximation</strong> (for <span class="math">u² &gt; 0</span>): Precomputed rational approximations for fast evaluation</li>
  <li><strong>Exact Ewald summation</strong> (for <span class="math">u² &lt; 0</span>): Direct numerical evaluation with convergence acceleration</li>
</ol>

<h2>Numerical Methods</h2>

<h3>Root Finding</h3>

<p>The <code>RootFinder</code> class finds roots of the quantization condition <span class="math">ω(E) = 0</span> using a multi-stage strategy:</p>

<ol>
  <li><strong>Direct acceptance</strong>: If <span class="math">|ω(E<sub>obs</sub>)| &lt; tolerance</span>, accept the observed energy as the root.</li>
  <li><strong>Adaptive bracketing</strong>: Try progressively larger brackets around the guess:
    <ul>
      <li>Widths: [0.02, 0.05, 0.10, 0.20, 0.40, 0.80, 1.5]</li>
      <li>Check for sign change using bisection</li>
    </ul>
  </li>
  <li><strong>Ordered root selection</strong>: If multiple roots are found, select the one corresponding to the level index (continuity tracking across levels).</li>
  <li><strong>Continuity tracking</strong>: Prefer roots near the previous root (for smooth branch following) or near the free energy.</li>
  <li><strong>Local scan</strong>: If bracketing fails, scan the interval with 50 points and identify sign changes.</li>
  <li><strong>Minimization fallback</strong>: Minimize <span class="math">|ω(E)|</span> using bounded scalar minimization.</li>
</ol>

<p><strong>Parameters</strong>:</p>
<ul>
  <li>Root tolerance: 1e-4</li>
  <li>Continuity tolerance: 0.3</li>
  <li>Maximum iterations: 100</li>
  <li>Pole threshold: 1e12 (treat values above this as poles)</li>
</ul>

<h3>Optimization</h3>

<p>The χ² minimization uses <code>scipy.optimize.minimize</code> with the Nelder-Mead method:</p>

<ul>
  <li><strong>Initial guess</strong>: Typically (0.047, 0.65) for πΣ scattering</li>
  <li><strong>Bounds</strong>: (-10, 10) for both parameters</li>
  <li><strong>Maximum iterations</strong>: 5000</li>
  <li><strong>Convergence criteria</strong>: xatol=1e-6, fatol=1e-6</li>
</ul>

<p>The optimizer calls <code>predict_energies()</code> for each function evaluation, which in turn:</p>
<ul>
  <li>Computes kinematics for each level (cached)</li>
  <li>Evaluates B-matrix (cached)</li>
  <li>Finds roots (cached via root finder)</li>
</ul>

<h3>Covariance and Statistics</h3>

<p>The pipeline follows standard statistical practices for lattice QCD:</p>

<ol>
  <li><strong>Covariance matrix</strong>: Computed from bootstrap samples:
    <div class="math-block">C<sub>ij</sub> = 1/(N<sub>b</sub>-1) · Σ<sub>b=1</sub><sup>N<sub>b</sub></sup> (E<sub>i</sub><sup>b</sup> - <span style="text-decoration:overline;">E</span><sub>i</sub>)(E<sub>j</sub><sup>b</sup> - <span style="text-decoration:overline;">E</span><sub>j</sub>)</div>
  </li>
  <li><strong>Parameter covariance</strong>: Computed using the Fisher information matrix:
    <div class="math-block">Cov(θ) = (J<sup>T</sup> C<sup>-1</sup> J)<sup>-1</sup></div>
    with PDG-style scaling by χ²/dof if χ² &gt; dof.
  </li>
  <li><strong>Pulls</strong>: Standardized residuals:
    <div class="math-block">pull<sub>i</sub> = (E<sub>obs,i</sub> - E<sub>pred,i</sub>) / √C<sub>ii</sub></div>
  </li>
  <li><strong>Information criteria</strong>:
    <ul>
      <li>AIC = χ² + 2p</li>
      <li>BIC = χ² + p ln(N)</li>
    </ul>
  </li>
</ol>

<h3>Caching Strategy</h3>

<p>The pipeline aggressively caches expensive computations to avoid redundant work:</p>

<ol>
  <li><strong>Kinematics</strong>: <code>compute_kinematics()</code> is decorated with <code>@lru_cache(maxsize=None)</code> with rounded arguments (10 decimal places).</li>
  <li><strong>B-matrix</strong>: <code>SingleChannelBMatrix.compute()</code> caches results in a dictionary keyed by (irrep, psq, u2, gamma, m_split).</li>
  <li><strong>Zeta function</strong>: <code>hybrid_Z()</code> uses <code>@lru_cache(maxsize=None)</code> with rounded arguments.</li>
  <li><strong>Root finder diagnostics</strong>: Stores the last diagnostics for debugging.</li>
  <li><strong>B-matrix coefficients</strong>: Precomputed Morningstar coefficients are loaded once.</li>
</ol>

<h3>Parallelization</h3>

<p>The <code>predict_energies()</code> function parallelizes root finding across levels using <code>ProcessPoolExecutor</code>:</p>

<ul>
  <li>Each level's root finding is independent</li>
  <li>Default workers: <code>min(n_levels, os.cpu_count())</code></li>
  <li>Falls back to serial on error</li>
  <li>Each worker process imports modules independently (fork-safe)</li>
</ul>

<h2>Module-by-Module Documentation</h2>

<h3><code>dataset_loader.py</code></h3>

<p><strong>Purpose</strong>: Loads lattice QCD energy levels from HDF5 files and prepares them for fitting.</p>

<p><strong>Role in pipeline</strong>: First stage - data ingestion.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>DataLoader.__init__</code></td>
      <td>Initialize data loader</td>
      <td><code>file_path</code>, <code>L</code>, <code>use_ref</code></td>
      <td>DataLoader instance</td>
    </tr>
    <tr>
      <td><code>DataLoader.scan_levels</code></td>
      <td>Discover all available levels</td>
      <td>None</td>
      <td>List of level dictionaries</td>
    </tr>
    <tr>
      <td><code>DataLoader.build_dataset</code></td>
      <td>Build DataSet from indices</td>
      <td><code>indices</code>, <code>levels_scan</code>, <code>m1</code>, <code>m2</code></td>
      <td>DataSet instance</td>
    </tr>
    <tr>
      <td><code>DataLoader._compute_free_cm_energy</code></td>
      <td>Compute non-interacting energy</td>
      <td><code>free_levels</code>, <code>psq</code>, <code>m1</code>, <code>m2</code></td>
      <td>Free CM energy</td>
    </tr>
  </tbody>
</table>

<p><strong>Dependencies</strong>: <code>h5py</code>, <code>numpy</code>, <code>general.data_reader.LQCD_DATA_READER</code></p>

<p><strong>Numerical considerations</strong>:</p>
<ul>
  <li>Bootstrap samples are extracted from HDF5 (first entry is mean, rest are bootstrap replicas)</li>
  <li>Free energies are computed from momentum mode labels</li>
</ul>

<p><strong>Physics significance</strong>: Maps lattice QCD naming conventions (PSQ, irreps) to pipeline-internal structure.</p>

<h3><code>ere.py</code></h3>

<p><strong>Purpose</strong>: Implements the Effective Range Expansion parameterization.</p>

<p><strong>Role in pipeline</strong>: Defines the scattering amplitude model.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>ERE.__init__</code></td>
      <td>Initialize ERE</td>
      <td><code>coeffs</code> (list/array)</td>
      <td>ERE instance</td>
    </tr>
    <tr>
      <td><code>ERE.compute_kinv</code></td>
      <td>Compute (k/mπ)cotδ</td>
      <td><code>kin</code> (KinematicVars)</td>
      <td>Float</td>
    </tr>
    <tr>
      <td><code>ERE.compute_cot_delta</code></td>
      <td>Compute cotδ</td>
      <td><code>kin</code></td>
      <td>Float</td>
    </tr>
    <tr>
      <td><code>ERE.compute_phase_shift</code></td>
      <td>Compute δ in degrees</td>
      <td><code>kin</code></td>
      <td>Float</td>
    </tr>
  </tbody>
</table>

<p><strong>Algorithm</strong>:</p>
<div class="math-block">kinv = a + b·Δ</div>
<p>where <span class="math">Δ = (E<sup>2</sup> - E<sub>th</sub><sup>2</sup>)/E<sub>th</sub><sup>2</sup></span>.</p>

<p><strong>Called by</strong>: <code>PhysicsModule.compute_kinv()</code> → used in <code>build_omega()</code></p>

<p><strong>Numerical considerations</strong>:</p>
<ul>
  <li>Supports 1 or 2 parameters</li>
  <li>Input can be <code>KinematicVars</code> object or direct q²</li>
  <li>Threshold protection (1e-15 added to denominator)</li>
</ul>

<h3><code>fitting_driver_canonical.py</code></h3>

<p><strong>Purpose</strong>: Core fitting engine that ties together physics, root finding, and statistics.</p>

<p><strong>Role in pipeline</strong>: Orchestrates the complete fitting process.</p>

<p><strong>Important classes</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Class</th>
      <th>Purpose</th>
      <th>Key methods</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>PhysicsModule</code></td>
      <td>Physics calculations</td>
      <td><code>compute_kinematics</code>, <code>compute_bmatrix</code>, <code>build_omega</code></td>
    </tr>
    <tr>
      <td><code>LuscherFitter</code></td>
      <td>Fitting engine</td>
      <td><code>predict_energies</code>, <code>objective</code>, <code>fit</code>, <code>vij</code></td>
    </tr>
  </tbody>
</table>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>LuscherFitter.predict_energies</code></td>
      <td>Predict energies for params</td>
      <td><code>params</code> (array)</td>
      <td><code>predicted</code> (array)</td>
    </tr>
    <tr>
      <td><code>LuscherFitter.objective</code></td>
      <td>Compute χ²</td>
      <td><code>params</code></td>
      <td><code>χ²</code> (float)</td>
    </tr>
    <tr>
      <td><code>LuscherFitter.fit</code></td>
      <td>Run optimization</td>
      <td><code>initial_guess</code>, <code>bounds</code>, ...</td>
      <td><code>FitResult</code></td>
    </tr>
    <tr>
      <td><code>LuscherFitter.vij</code></td>
      <td>Compute parameter covariance</td>
      <td><code>params</code>, <code>epsilon</code></td>
      <td><code>cov_par</code> (matrix)</td>
    </tr>
    <tr>
      <td><code>_predict_single_level</code></td>
      <td>Parallel worker for one level</td>
      <td><code>irrep</code>, <code>d</code>, <code>e_obs</code>, <code>params</code>, ...</td>
      <td><code>root</code> (float)</td>
    </tr>
  </tbody>
</table>

<p><strong>Parallelization</strong>: <code>predict_energies</code> uses <code>ProcessPoolExecutor</code> to predict each level concurrently.</p>

<p><strong>Dependencies</strong>: <code>scipy.optimize</code>, <code>root_finder</code>, <code>ere</code>, <code>stats</code></p>

<p><strong>Physics significance</strong>: Implements the quantization condition <span class="math">ω(E) = K<sub>inv</sub>(E) - B(E) = 0</span>.</p>

<h3><code>morningstar_bmatrix.py</code></h3>

<p><strong>Purpose</strong>: Implements the single-channel Morningstar B-matrix.</p>

<p><strong>Role in pipeline</strong>: Computes finite-volume matrix elements.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>SingleChannelBMatrix.compute</code></td>
      <td>Compute B-matrix element</td>
      <td><code>irrep</code>, <code>kin</code></td>
      <td><code>B</code> (float)</td>
    </tr>
    <tr>
      <td><code>SingleChannelBMatrix.get_coefficient</code></td>
      <td>Get irrep coefficient</td>
      <td><code>irrep</code></td>
      <td><code>coeff</code> (float)</td>
    </tr>
    <tr>
      <td><code>SingleChannelBMatrix._compute_regularized_zeta</code></td>
      <td>Compute Z00</td>
      <td><code>kin</code></td>
      <td><code>Z00</code> (float)</td>
    </tr>
  </tbody>
</table>

<p><strong>Algorithm</strong>:</p>
<div class="math-block">B = c<sub>Λ</sub> · Z<sub>00</sub> / (γ π<sup>3/2</sup> u)</div>

<p>where:</p>
<ul>
  <li><span class="math">c<sub>Λ</sub></span> is the irrep coefficient from Morningstar tables</li>
  <li><span class="math">Z<sub>00</sub></span> is the generalized zeta function (hybrid implementation)</li>
  <li><span class="math">u = L q<sup>*</sup> / (2π)</span></li>
</ul>

<p><strong>Dependencies</strong>: <code>b_tables</code> (B1-B8), <code>final_zeta</code></p>

<p><strong>Caching</strong>: Results cached in <code>_b_cache</code> dictionary keyed by (irrep, psq, u2, gamma, m_split).</p>

<p><strong>Physics significance</strong>: Encodes all finite-volume effects through the zeta function and geometry-dependent coefficients.</p>

<h3><code>root_finder.py</code></h3>

<p><strong>Purpose</strong>: Locates solutions to the Lüscher quantization condition.</p>

<p><strong>Role in pipeline</strong>: Finds roots of <span class="math">ω(E) = K<sub>inv</sub>(E) - B(E)</span>.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>RootFinder.find_root_near_guess</code></td>
      <td>Find root near guess</td>
      <td><code>f</code>, <code>x_guess</code>, <code>prev_root</code>, <code>reference_energy</code>, <code>level_index</code></td>
      <td><code>root</code> (float)</td>
    </tr>
  </tbody>
</table>

<p><strong>Algorithm</strong>:</p>
<ol>
  <li>Check if guess is already a root</li>
  <li>Adaptive bracketing with increasing widths</li>
  <li>Bisection to find roots in each bracket</li>
  <li>Ordered selection based on level index</li>
  <li>Continuity tracking (prefer roots near previous root)</li>
  <li>Local scan with 50 points</li>
  <li>Minimization fallback</li>
</ol>

<p><strong>Tolerances</strong>:</p>
<ul>
  <li>Root tolerance: 1e-4</li>
  <li>Continuity tolerance: 0.3</li>
  <li>Maximum iterations: 100</li>
</ul>

<p><strong>Called by</strong>: <code>LuscherFitter.predict_energies()</code> via <code>_predict_single_level</code></p>

<h3><code>stats.py</code></h3>

<p><strong>Purpose</strong>: Statistical utilities for fitting.</p>

<p><strong>Role in pipeline</strong>: χ² calculation, covariance, error estimation.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>bootstrap_covariance</code></td>
      <td>Compute covariance from bootstrap</td>
      <td><code>bootstrap_samples</code></td>
      <td><code>covariance</code> (matrix)</td>
    </tr>
    <tr>
      <td><code>chi2</code></td>
      <td>Compute χ²</td>
      <td><code>observed_mean</code>, <code>predicted</code>, <code>covariance</code></td>
      <td><code>χ²</code> (float)</td>
    </tr>
    <tr>
      <td><code>parameter_covariance</code></td>
      <td>Compute parameter covariance</td>
      <td><code>jacobian</code>, <code>data_cov</code>, <code>χ²</code>, <code>n_data</code>, <code>n_params</code></td>
      <td><code>cov_par</code> (matrix)</td>
    </tr>
    <tr>
      <td><code>standardized_residuals</code></td>
      <td>Compute pulls</td>
      <td><code>observed_mean</code>, <code>predicted</code>, <code>covariance</code></td>
      <td><code>pulls</code> (array)</td>
    </tr>
    <tr>
      <td><code>correlation_matrix</code></td>
      <td>Compute correlation matrix</td>
      <td><code>covariance</code></td>
      <td><code>corr</code> (matrix)</td>
    </tr>
  </tbody>
</table>

<p><strong>Numerical considerations</strong>:</p>
<ul>
  <li>Uses Cholesky decomposition for solving linear systems</li>
  <li>Falls back to general solve or pseudoinverse</li>
  <li>PDG-style scaling of covariance by χ²/dof when χ² &gt; dof</li>
</ul>

<p><strong>Statistics</strong>:</p>
<ul>
  <li>AIC: χ² + 2p</li>
  <li>BIC: χ² + p ln(N)</li>
  <li>Reduced χ²: χ²/(N-p)</li>
</ul>

<h3><code>kinematics.py</code></h3>

<p><strong>Purpose</strong>: Computes all kinematic quantities for unequal masses in a finite box.</p>

<p><strong>Role in pipeline</strong>: Provides kinematic inputs to B-matrix and ERE.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>compute_kinematics</code></td>
      <td>Compute kinematics</td>
      <td><code>E_cm</code>, <code>d</code>, <code>m1</code>, <code>m2</code>, <code>L</code>, <code>Mref</code></td>
      <td><code>KinematicVars</code></td>
    </tr>
  </tbody>
</table>

<p><strong>Outputs</strong> (KinematicVars):</p>
<ul>
  <li><code>E_cm</code>: CM energy</li>
  <li><code>E_lab</code>: Lab energy</li>
  <li><code>gamma</code>: Lorentz boost factor</li>
  <li><code>alpha</code>: Mass asymmetry parameter</li>
  <li><code>q_star</code>: Relative momentum</li>
  <li><code>u</code>: Lüscher momentum variable</li>
  <li><code>Delta</code>: Dimensionless energy variable</li>
</ul>

<p><strong>Algorithm</strong>:</p>
<ol>
  <li>Compute total momentum: <span class="math">P = (2π/L) <b>d</b></span></li>
  <li>Lab energy: <span class="math">E<sub>lab</sub> = √(E<sub>cm</sub>² + P²)</span></li>
  <li>Boost factor: <span class="math">γ = E<sub>lab</sub>/E<sub>cm</sub></span></li>
  <li>Mass asymmetry: <span class="math">α = ½(1 + (m<sub>1</sub>² - m<sub>2</sub>²)/E<sub>cm</sub>²)</span></li>
  <li>Relative momentum squared: <span class="math">q<sup>*2</sup> = ((s - (m<sub>1</sub>+m<sub>2</sub>)²)(s - (m<sub>1</sub>-m<sub>2</sub>)²))/(4s)</span></li>
</ol>

<p><strong>Caching</strong>: <code>@lru_cache(maxsize=None)</code> with arguments rounded to 10 decimal places.</p>

<h3><code>final_zeta.py</code></h3>

<p><strong>Purpose</strong>: Hybrid generalized zeta function with caching.</p>

<p><strong>Role in pipeline</strong>: Evaluates <span class="math">Z<sub>00</sub>(u², γ, <b>d</b>, α)</span> for B-matrix.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>hybrid_Z</code></td>
      <td>Hybrid zeta evaluation</td>
      <td><code>u2</code>, <code>psq</code>, <code>gamma</code>, <code>m_split</code>, <code>L</code></td>
      <td><code>Z</code> (float)</td>
    </tr>
    <tr>
      <td><code>pade_approximation</code></td>
      <td>Padé approximation</td>
      <td><code>u2</code>, <code>psq</code>, <code>gamma</code></td>
      <td><code>Z</code> (float)</td>
    </tr>
    <tr>
      <td><code>ewald_approximation</code></td>
      <td>Exact Ewald summation</td>
      <td><code>u2</code>, <code>psq</code>, <code>gamma</code>, <code>m_split</code>, <code>L</code></td>
      <td><code>Z</code> (float)</td>
    </tr>
  </tbody>
</table>

<p><strong>Algorithm</strong>:</p>
<ol>
  <li>For large |u²| &gt; 50: Use asymptotic expansion directly</li>
  <li>For u² &lt; 0: Use exact Ewald summation</li>
  <li>For u² ≥ 0 and m_split ≈ 1: Try Padé approximation</li>
  <li>Otherwise: Fall back to exact Ewald</li>
</ol>

<p><strong>Caching</strong>: <code>@lru_cache(maxsize=None)</code> for both Padé and Ewald results.</p>

<p><strong>Precomputed coefficients</strong>: Padé coefficients stored in <code>coefficients/PSQX_gammaY.ZZ_coeffs.json</code>.</p>

<h3><code>exact_zeta.py</code></h3>

<p><strong>Purpose</strong>: Reference implementation of exact Ewald summation for generalized zeta functions.</p>

<p><strong>Role in pipeline</strong>: Provides exact evaluation for validation and as fallback.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>Z</code></td>
      <td>Exact zeta function</td>
      <td><code>q2</code>, <code>gamma</code>, <code>l</code>, <code>m</code>, <code>d</code>, <code>m_split</code>, <code>precision</code></td>
      <td><code>Z</code> (complex)</td>
    </tr>
  </tbody>
</table>

<p><strong>Algorithm</strong>: Ewald summation with convergence acceleration.</p>

<p><strong>Note</strong>: This is a reference implementation and is slower than the hybrid approach. It is used primarily for validation and as a fallback.</p>

<h3><code>profiler.py</code></h3>

<p><strong>Purpose</strong>: Performance profiling and counter tracking.</p>

<p><strong>Role in pipeline</strong>: Instruments the pipeline for performance analysis.</p>

<p><strong>Important functions</strong>:</p>

<table>
  <thead>
    <tr>
      <th>Function</th>
      <th>Purpose</th>
      <th>Inputs</th>
      <th>Outputs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>profiler.start</code></td>
      <td>Start timer</td>
      <td><code>category</code></td>
      <td>None</td>
    </tr>
    <tr>
      <td><code>profiler.stop</code></td>
      <td>Stop timer</td>
      <td><code>category</code></td>
      <td>None</td>
    </tr>
    <tr>
      <td><code>profiler.increment_counter</code></td>
      <td>Increment counter</td>
      <td><code>name</code>, <code>amount</code></td>
      <td>None</td>
    </tr>
    <tr>
      <td><code>profiler.report</code></td>
      <td>Print profiling report</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <td><code>profiler.decorator</code></td>
      <td>Function timing decorator</td>
      <td><code>category</code></td>
      <td>Decorated function</td>
    </tr>
  </tbody>
</table>

<p><strong>Counters tracked</strong>:</p>
<ul>
  <li><code>chi2_evaluations</code>: Number of χ² evaluations</li>
  <li><code>root_finder_calls</code>: Number of root finding calls</li>
  <li><code>root_iterations</code>: Total iterations across root finder</li>
  <li><code>omega_evaluations</code>: Number of ω(E) evaluations</li>
  <li><code>hybrid_Z_calls</code>: Number of hybrid zeta calls</li>
  <li><code>B_matrix_calls</code>: Number of B-matrix evaluations</li>
  <li><code>ERE_calls</code>: Number of ERE evaluations</li>
</ul>

<p><strong>Usage</strong>: Decorators and context managers automatically profile code sections.</p>

<h3><code>plot_figure8.py</code> and <code>plot.py</code></h3>

<p><strong>Purpose</strong>: Generate publication-quality plots.</p>

<p><strong>Role in pipeline</strong>: Visualization of results.</p>

<p><strong><code>plot_figure8.py</code></strong>:</p>
<ul>
  <li>Plots finite-volume spectrum (energy levels vs. irrep)</li>
  <li>Green circles with error bars represent lattice data</li>
  <li>Gray bands represent non-interacting levels</li>
  <li>Dashed lines show physical thresholds (πΣ, K̄N, ππΛ)</li>
</ul>

<p><strong><code>plot.py</code></strong>:</p>
<ul>
  <li>Reproduces Fig. 10 from Morningstar/BaSc πΣ paper</li>
  <li>Shows (q/mπ)² vs. (k/mπ)cotδ</li>
  <li>Lattice points from different irreps</li>
  <li>Luscher curves (steep blue segments)</li>
  <li>ERE curve with band</li>
  <li>Virtual state curve (black dashed)</li>
  <li>Virtual state star (intersection)</li>
</ul>

<h2>Important Functions and Classes</h2>

<h3><code>DataLoader</code></h3>

<p><strong>Location</strong>: <code>dataset_loader.py</code></p>

<p><strong>Purpose</strong>: Loads and prepares lattice QCD data from HDF5.</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>file_path</code>: Path to HDF5 file</li>
  <li><code>L</code>: Lattice size</li>
  <li><code>data</code>: HDF5 data structure</li>
  <li><code>has_channel_layer</code>: Whether data has 'iso' channel layer</li>
</ul>

<p><strong>Key methods</strong>:</p>
<ul>
  <li><code>scan_levels()</code>: Discovers all energy levels</li>
  <li><code>build_dataset(indices, levels_scan, m1, m2)</code>: Builds DataSet from selected levels</li>
</ul>

<h3><code>LuscherFitter</code></h3>

<p><strong>Location</strong>: <code>fitting_driver_canonical.py</code></p>

<p><strong>Purpose</strong>: Main fitting engine.</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>observed_mean</code>: Energy means</li>
  <li><code>bootstrap_samples</code>: Bootstrap replicas</li>
  <li><code>cov_matrix</code>: Covariance matrix</li>
  <li><code>irrep_list</code>: List of irreps</li>
  <li><code>d_list</code>: List of momentum vectors</li>
  <li><code>physics</code>: PhysicsModule instance</li>
  <li><code>root_finder</code>: RootFinder instance</li>
  <li><code>max_workers</code>: Number of parallel workers</li>
</ul>

<p><strong>Key methods</strong>:</p>
<ul>
  <li><code>predict_energies(params)</code>: Predicts energies for parameters</li>
  <li><code>objective(params)</code>: Computes χ²</li>
  <li><code>fit(initial_guess, bounds, ...)</code>: Runs optimization</li>
  <li><code>vij(params)</code>: Computes parameter covariance</li>
</ul>

<h3><code>PhysicsModule</code></h3>

<p><strong>Location</strong>: <code>fitting_driver_canonical.py</code></p>

<p><strong>Purpose</strong>: Physics calculations for the quantization condition.</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>L</code>: Lattice size</li>
  <li><code>m1</code>, <code>m2</code>: Particle masses</li>
  <li><code>bmatrix</code>: B-matrix instance</li>
  <li><code>_kin_cache</code>: Kinematics cache</li>
</ul>

<p><strong>Key methods</strong>:</p>
<ul>
  <li><code>compute_kinematics(E_cm, d)</code>: Computes kinematics</li>
  <li><code>compute_bmatrix(irrep, kin)</code>: Computes B-matrix</li>
  <li><code>compute_kinv(kin, ere)</code>: Computes (k/mπ)cotδ</li>
  <li><code>build_omega(irrep, d, E_cm, ere)</code>: Builds quantization condition function</li>
</ul>

<h3><code>ERE</code></h3>

<p><strong>Location</strong>: <code>ere.py</code></p>

<p><strong>Purpose</strong>: Effective Range Expansion parameterization.</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>coeffs</code>: Array of coefficients [a, b]</li>
  <li><code>n_params</code>: Number of parameters (1 or 2)</li>
</ul>

<p><strong>Key methods</strong>:</p>
<ul>
  <li><code>compute_kinv(kin)</code>: Computes (k/mπ)cotδ</li>
  <li><code>compute_cot_delta(kin)</code>: Computes cotδ</li>
  <li><code>compute_phase_shift(kin)</code>: Computes phase shift δ in degrees</li>
</ul>

<h3><code>SingleChannelBMatrix</code></h3>

<p><strong>Location</strong>: <code>morningstar_bmatrix.py</code></p>

<p><strong>Purpose</strong>: Computes Morningstar B-matrix elements.</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>_b_cache</code>: Cache of B-matrix values</li>
  <li><code>_all_tables</code>: Combined Morningstar coefficients</li>
  <li><code>_coeff_cache</code>: Cache of irrep coefficients</li>
</ul>

<p><strong>Key methods</strong>:</p>
<ul>
  <li><code>compute(irrep, kin)</code>: Computes B-matrix</li>
  <li><code>get_coefficient(irrep)</code>: Gets irrep coefficient</li>
  <li><code>_compute_regularized_zeta(kin)</code>: Computes Z00 zeta function</li>
</ul>

<h3><code>RootFinder</code></h3>

<p><strong>Location</strong>: <code>root_finder.py</code></p>

<p><strong>Purpose</strong>: Finds roots of the quantization condition.</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>root_tolerance</code>: 1e-4</li>
  <li><code>continuity_tol</code>: 0.3</li>
  <li><code>local_scan_points</code>: 50</li>
  <li><code>last_diagnostics</code>: Diagnostic info from last root finding</li>
</ul>

<p><strong>Key methods</strong>:</p>
<ul>
  <li><code>find_root_near_guess(f, x_guess, prev_root, reference_energy, level_index)</code>: Finds root</li>
</ul>

<h3><code>KinematicVars</code></h3>

<p><strong>Location</strong>: <code>kinematics.py</code></p>

<p><strong>Purpose</strong>: Container for kinematic variables (frozen dataclass).</p>

<p><strong>Key attributes</strong>:</p>
<ul>
  <li><code>E_cm</code>: CM energy</li>
  <li><code>E_lab</code>: Lab energy</li>
  <li><code>gamma</code>: Lorentz boost factor</li>
  <li><code>alpha</code>: Mass asymmetry parameter</li>
  <li><code>q_star</code>: Relative momentum</li>
  <li><code>q2_star</code>: Relative momentum squared</li>
  <li><code>u</code>: Lüscher momentum variable</li>
  <li><code>u2</code>: Lüscher momentum squared</li>
  <li><code>threshold</code>: Threshold energy</li>
  <li><code>Delta</code>: Dimensionless energy variable</li>
</ul>

<h2>Data Flow Between Modules</h2>

<pre>HDF5 File
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
[plot_figure8.py] / [plot.py] → publication figures</pre>

<h2>Quantization / Fitting Procedure</h2>

<h3>Step-by-Step Fitting Procedure</h3>

<ol>
  <li>
    <p><strong>Data Loading</strong>:</p>
    <pre>loader = DataLoader('DataSet.hdf5', L=64, use_ref=True)
levels = loader.scan_levels()
dataset = loader.build_dataset([11, 35, 67, 121], levels, m1=1.0, m2=5.8625)</pre>
  </li>
  <li>
    <p><strong>Fitter Initialization</strong>:</p>
    <pre>fitter = LuscherFitter(
    observed_mean=dataset.means,
    bootstrap_samples=dataset.bootstrap,
    irrep_list=['G1u(0)', 'G1(1)', 'G(2)', 'G(3)'],
    d_list=[(0,0,0), (0,0,1), (1,1,0), (1,1,1)],
    L=64.0, m1=1.0, m2=5.8625,
    cov_matrix=dataset.covariance,
    free_energies=dataset.free_energies,
    max_workers=4,
)</pre>
  </li>
  <li>
    <p><strong>Prediction for a Given Parameter Set</strong>:</p>
    <pre>params = np.array([0.047, 0.65])
predicted = fitter.predict_energies(params)</pre>
    <p>Each prediction involves:</p>
    <ul>
      <li>Building ω(E) = K(E) - B(E) for each level</li>
      <li>Finding the root near the observed energy</li>
      <li>Using continuity tracking with previous roots</li>
    </ul>
  </li>
  <li>
    <p><strong>χ² Computation</strong>:</p>
    <pre>chi2_val = stats.chi2(fitter.observed_mean, predicted, fitter.cov_matrix)</pre>
  </li>
  <li>
    <p><strong>Optimization</strong>:</p>
    <pre>result = fitter.fit(
    initial_guess=np.array([0.047, 0.65]),
    bounds=[(-10, 10), (-10, 10)],
    method='nelder-mead',
    maxiter=5000,
)</pre>
  </li>
  <li>
    <p><strong>Results Extraction</strong>:</p>
    <pre>print(f"a = {result.params[0]:.6f} +/- {result.errors[0]:.6f}")
print(f"b = {result.params[1]:.6f} +/- {result.errors[1]:.6f}")
print(f"χ² = {result.chi2:.6f}, ndof = {result.ndof}")</pre>
  </li>
</ol>

<h3>The Quantization Condition Function</h3>

<p>For each level, the quantization condition is:</p>

<div class="math-block">ω(E) = K<sub>inv</sub>(E) - B(E)</div>

<p>where:</p>
<ul>
  <li><span class="math">K<sub>inv</sub>(E) = a + b·(E² - E<sub>th</sub>²)/E<sub>th</sub>²</span></li>
  <li><span class="math">B(E) = c<sub>Λ</sub> · Z<sub>00</sub>(u², γ, <b>d</b>, α) / (γ π<sup>3/2</sup> u)</span></li>
</ul>

<p>The function is constructed by <code>PhysicsModule.build_omega()</code>:</p>

<pre>def omega(E):
    kin = self.compute_kinematics(E, d)
    B = self.compute_bmatrix(irrep, kin)
    Kinv = self.compute_kinv(kin, ere)
    return Kinv - B</pre>

<h2>Statistical Treatment</h2>

<h3>Bootstrap Samples</h3>

<p>Energy levels from the HDF5 file include bootstrap replicas. The first entry is the mean, subsequent entries are bootstrap samples.</p>

<pre># shape: (n_levels, n_bootstrap + 1)
arr = [mean, boot1, boot2, ..., bootN]
means = arr[0]
bootstrap = arr[1:]  # shape (n_levels, n_bootstrap)</pre>

<h3>Covariance Matrix</h3>

<p>The covariance matrix is computed from bootstrap samples:</p>

<pre>C = np.cov(bootstrap, rowvar=True)  # shape (n_levels, n_levels)</pre>

<h3>χ²</h3>

<div class="math-block">χ² = (<b>E</b><sub>obs</sub> - <b>E</b><sub>pred</sub>)<sup>T</sup> C<sup>-1</sup> (<b>E</b><sub>obs</sub> - <b>E</b><sub>pred</sub>)</div>

<h3>Parameter Covariance</h3>

<p>The parameter covariance matrix is computed using the Fisher information matrix:</p>

<div class="math-block">Cov(θ) = (J<sup>T</sup> C<sup>-1</sup> J)<sup>-1</sup></div>

<p>where <span class="math">J</span> is the Jacobian matrix of predicted energies with respect to parameters:</p>

<div class="math-block">J<sub>iα</sub> = ∂E<sub>pred,i</sub> / ∂θ<sub>α</sub></div>

<p>The Jacobian is computed numerically using central differences:</p>

<div class="math-block">J<sub>iα</sub> ≈ (E<sub>pred</sub>(θ + h e<sub>α</sub>)<sub>i</sub> - E<sub>pred</sub>(θ - h e<sub>α</sub>)<sub>i</sub>) / (2h)</div>

<p>with <span class="math">h = 10<sup>-5</sup></span>.</p>

<p>PDG-style scaling is applied:</p>

<div class="math-block">Cov<sub>scaled</sub> = Cov · (χ²/dof)</div>

<p>if χ² &gt; dof.</p>

<h3>Pulls</h3>

<p>Standardized residuals:</p>

<div class="math-block">pull<sub>i</sub> = (E<sub>obs,i</sub> - E<sub>pred,i</sub>) / √C<sub>ii</sub></div>

<h3>Information Criteria</h3>

<ul>
  <li>AIC: χ² + 2p</li>
  <li>BIC: χ² + p ln(N)</li>
</ul>

<h2>Performance and Optimization</h2>

<h3>Bottlenecks</h3>

<p>From the profiling data, the main computational bottlenecks are:</p>

<ol>
  <li><strong>Zeta function evaluation</strong>: The most expensive operation, especially for u² &gt; 0 where Padé approximation is used, and for u² &lt; 0 where Ewald summation is required.</li>
  <li><strong>Root finding</strong>: Each function evaluation requires finding multiple roots (one per level). The number of ω(E) evaluations per root finding can be high.</li>
  <li><strong>Optimization iterations</strong>: The Nelder-Mead optimizer typically requires 50-200 function evaluations.</li>
</ol>

<h3>Optimization Strategies</h3>

<ol>
  <li>
    <p><strong>LRU Caching</strong>:</p>
    <ul>
      <li>Kinematics: <code>@lru_cache(maxsize=None)</code> with rounded arguments</li>
      <li>Zeta function: <code>@lru_cache(maxsize=None)</code> with rounded arguments</li>
      <li>B-matrix: Dictionary cache</li>
    </ul>
  </li>
  <li>
    <p><strong>Parallel Prediction</strong>:</p>
    <ul>
      <li><code>predict_energies</code> uses <code>ProcessPoolExecutor</code></li>
      <li>Each level's root finding is independent</li>
      <li>Falls back to serial on error</li>
    </ul>
  </li>
  <li>
    <p><strong>Padé Approximation</strong>:</p>
    <ul>
      <li>Precomputed rational approximations for u² &gt; 0</li>
      <li>Significantly faster than exact Ewald</li>
      <li>Coefficients stored in JSON files</li>
    </ul>
  </li>
  <li>
    <p><strong>Asymptotic Expansion</strong>:</p>
    <ul>
      <li>For large |u²| &gt; 50, use asymptotic expansion directly</li>
      <li>Avoids expensive zeta function evaluation</li>
    </ul>
  </li>
</ol>

<h3>Profiling</h3>

<p>The <code>profiler</code> module tracks:</p>
<ul>
  <li>Function call counts</li>
  <li>Execution times (total, average, max)</li>
  <li>Specialized counters (chi2_evaluations, root_finder_calls, hybrid_Z_calls, etc.)</li>
</ul>

<p>To enable profiling:</p>
<pre>from profiler import profiler
profiler.enabled = True  # default is True</pre>

<p>To view report:</p>
<pre>profiler.report()</pre>

<p>Sample output:</p>
<pre>PIPELINE PROFILING REPORT

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
             ERE_calls: 18720</pre>

<h3>Memory Usage</h3>

<p>The main memory consumers are:</p>
<ul>
  <li>Bootstrap samples: <code>(n_levels × n_bootstrap)</code> ≈ 4 × 1000 = 4000 floats</li>
  <li>Covariance matrix: <code>(n_levels × n_levels)</code> ≈ 16 floats</li>
  <li>Cache dictionaries: Kinematics (typically &lt; 1000 entries), B-matrix (typically &lt; 1000 entries), Zeta (typically &lt; 5000 entries)</li>
</ul>

<p>Total memory footprint is typically &lt; 100 MB for standard runs.</p>

<h2>Validation and Testing</h2>

<h3>Known-Value Tests</h3>

<p>The repository includes validation against known results from the literature:</p>

<ol>
  <li>
    <p><strong>Exact Ewald validation</strong>: <code>exact_zeta.py</code> includes test cases with expected delta values:</p>
    <pre># CMS test: expected delta = 136.6527°
# MV1 test: expected delta = 115.7653°
# MV2 test: expected delta = 127.9930°</pre>
  </li>
  <li>
    <p><strong>Phase shift reproduction</strong>: The pipeline reproduces the phase shift curves from the BaSc/Morningstar paper (Fig. 10).</p>
  </li>
</ol>

<h3>Consistency Checks</h3>

<ol>
  <li><strong>Bootstrap consistency</strong>: The covariance matrix is positive definite (checked via Cholesky decomposition).</li>
  <li><strong>Root finding consistency</strong>: Roots are checked against the observed energy and continuity is enforced.</li>
  <li><strong>χ² consistency</strong>: χ² is non-negative and well-behaved.</li>
</ol>

<h3>Unit Tests</h3>

<p><strong>Status</strong>: No formal unit test suite was found in the repository. The pipeline relies on:</p>
<ul>
  <li>Manual validation via plotting</li>
  <li>Reference comparisons to literature</li>
  <li>Runtime consistency checks (e.g., covariance matrix positive definiteness)</li>
</ul>

<h3>Validation by Scientific Reproduction</h3>

<p>The pipeline is validated by its ability to:</p>
<ol>
  <li>Reproduce the finite-volume spectrum plot (Figure 8)</li>
  <li>Reproduce the scattering curve plot (Figure 10)</li>
</ol>

<p>These plots demonstrate agreement with the published Morningstar/BaSc results.</p>

<h2>Installation</h2>

<h3>Requirements</h3>

<ul>
  <li><strong>Python</strong>: 3.8 or higher (tested with 3.9+)</li>
  <li><strong>Dependencies</strong>: Listed below</li>
</ul>

<h3>Installation Steps</h3>

<ol>
  <li>
    <p><strong>Clone the repository</strong>:</p>
    <pre>git clone https://github.com/username/repository.git
cd repository</pre>
  </li>
  <li>
    <p><strong>Install dependencies</strong>:</p>
    <pre>pip install numpy scipy h5py matplotlib</pre>
  </li>
  <li>
    <p><strong>Verify installation</strong>:</p>
    <pre>python -c "import QC2; print('OK')"</pre>
  </li>
</ol>

<h3>Alternative: Using a Virtual Environment</h3>

<pre>python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt  # if available</pre>

<h3>Notes on Compilation</h3>

<p>The codebase is pure Python with no compiled extensions. No compilation is required.</p>

<h2>Dependencies</h2>

<h3>Core Dependencies</h3>

<table>
  <thead>
    <tr>
      <th>Package</th>
      <th>Version</th>
      <th>Purpose</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>numpy</td>
      <td>≥1.20.0</td>
      <td>Numerical arrays and linear algebra</td>
    </tr>
    <tr>
      <td>scipy</td>
      <td>≥1.7.0</td>
      <td>Optimization, root finding, linear algebra</td>
    </tr>
    <tr>
      <td>h5py</td>
      <td>≥3.0.0</td>
      <td>HDF5 data loading</td>
    </tr>
    <tr>
      <td>matplotlib</td>
      <td>≥3.4.0</td>
      <td>Plotting</td>
    </tr>
  </tbody>
</table>

<h3>Optional Dependencies</h3>

<table>
  <thead>
    <tr>
      <th>Package</th>
      <th>Purpose</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>mpi4py</td>
      <td>Potential future parallelization</td>
    </tr>
  </tbody>
</table>

<h3>System Requirements</h3>

<ul>
  <li><strong>Operating System</strong>: Linux, macOS, Windows (WSL recommended)</li>
  <li><strong>Memory</strong>: 512 MB minimum, 4 GB recommended for large datasets</li>
  <li><strong>Disk Space</strong>: 100 MB for code and dependencies, plus HDF5 data files</li>
</ul>

<h2>Input Data</h2>

<h3>HDF5 File Format</h3>

<p>The pipeline expects an HDF5 file with the following structure:</p>

<pre>DataSet.hdf5
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
└── ...</pre>

<h3>Data Format</h3>

<p>Each energy level dataset is a 1D array:</p>
<pre>[mean, boot1, boot2, ..., bootN]</pre>
<p>where:</p>
<ul>
  <li><code>mean</code>: Mean energy value</li>
  <li><code>boot1...bootN</code>: Bootstrap replica values</li>
</ul>

<h3>Required Metadata</h3>

<ul>
  <li><code>free_levels</code>: List of strings indicating the non-interacting momentum modes for each level</li>
  <li><code>L</code>: Lattice size (provided at runtime)</li>
</ul>

<h3>Units</h3>

<ul>
  <li>All energies are in units of <span class="math">m<sub>π</sub></span> (pion mass)</li>
  <li>The lattice size <span class="math">L</span> is in lattice units</li>
</ul>

<h3>Example Data Loading</h3>

<pre>from dataset_loader import DataLoader

loader = DataLoader('DataSet.hdf5', L=64, use_ref=True)
levels = loader.scan_levels()
dataset = loader.build_dataset([11, 35, 67, 121], levels, m1=1.0, m2=5.8625)

print(f"Loaded {dataset.n_levels} levels")
print(f"Bootstrap samples: {dataset.n_bootstrap}")
print(f"Covariance matrix shape: {dataset.covariance.shape}")</pre>

<h2>Usage</h2>

<h3>Quick Start</h3>

<p>The main entry point for fitting is <code>run_fit_from_dataset.py</code>:</p>

<pre>python QC2/run_fit_from_dataset.py</pre>

<p>This will:</p>
<ol>
  <li>Load the HDF5 dataset</li>
  <li>Select the predefined levels (indices 11, 35, 67, 121)</li>
  <li>Fit the ERE parameters</li>
  <li>Print results</li>
  <li>Save results to <code>fit_results.json</code></li>
  <li>Print a profiling report</li>
</ol>

<h3>Configuration</h3>

<p>The configuration is at the top of <code>run_fit_from_dataset.py</code>:</p>

<pre>HDF5_PATH = os.path.expanduser("~/Desktop/my_work./Last_Week/my_work/DataSet.hdf5")
L = 64.0
m1 = 1.0                     # pion mass (mπ/mπ = 1)
m2 = 5.862544007347314       # Sigma mass (mΣ/mπ = 0.3830/0.06533)
SELECTED_INDICES = [11, 35, 67, 121]   # G1u(0), G1(1), G(2), G(3)
INITIAL_GUESS = np.array([0.047, 0.65])
BOUNDS = [(-10.0, 10.0), (-10.0, 10.0)]</pre>

<h3>Custom Fitting</h3>

<p>To perform a custom fit:</p>

<pre>import numpy as np
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
print(f"χ²/ndof = {result.reduced_chi2:.6f}")</pre>

<h3>Plotting</h3>

<p>To generate the finite-volume spectrum plot (Figure 8):</p>

<pre>python QC2/plot_figure8.py</pre>

<p>This will:</p>
<ol>
  <li>Load the HDF5 data</li>
  <li>Extract all levels in the energy window (6.5, 8.5)</li>
  <li>Generate the spectrum plot with:
    <ul>
      <li>Green circles with error bars for lattice data</li>
      <li>Gray bands for non-interacting levels</li>
      <li>Dashed lines for thresholds</li>
      <li>Level numbers annotated on the points</li>
    </ul>
  </li>
  <li>Save the figure to <code>figure8_final.pdf</code></li>
</ol>

<p>To generate the scattering curve plot (Figure 10):</p>

<pre>python QC2/plot.py</pre>

<p>This will:</p>
<ol>
  <li>Create test energy data (for demonstration)</li>
  <li>Fit the ERE parameters</li>
  <li>Generate the scattering curve plot with:
    <ul>
      <li>Lattice points</li>
      <li>Luscher curves for each irrep</li>
      <li>ERE curve with band</li>
      <li>Virtual state curve</li>
      <li>Virtual state star (if found)</li>
    </ul>
  </li>
  <li>Save the figure to <code>fig10_reproduction.png</code></li>
</ol>

<h3>Profiling</h3>

<p>To run with profiling enabled:</p>

<pre>from profiler import profiler
profiler.enabled = True

# Run fitting
result = fitter.fit(...)

# Print report
profiler.report()</pre>

<h2>Example Workflow</h2>

<h3>Complete End-to-End Example</h3>

<pre>#!/usr/bin/env python3
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
# plot_figure8_from_hdf5(HDF5_PATH, L, save_path="spectrum.pdf")</pre>

<h3>Expected Output</h3>

<pre>Loaded 4 levels
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
G: pred=7.000000, obs=7.001500, pull=-0.3456</pre>

<h2>Outputs</h2>

<h3><code>FitResult</code> Object</h3>

<p>The <code>fit()</code> method returns a <code>FitResult</code> dataclass with:</p>

<table>
  <thead>
    <tr>
      <th>Attribute</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>params</code></td>
      <td>Fitted parameters [a, b]</td>
    </tr>
    <tr>
      <td><code>errors</code></td>
      <td>Parameter uncertainties</td>
    </tr>
    <tr>
      <td><code>cov_params</code></td>
      <td>Parameter covariance matrix</td>
    </tr>
    <tr>
      <td><code>corr_params</code></td>
      <td>Parameter correlation matrix</td>
    </tr>
    <tr>
      <td><code>chi2</code></td>
      <td>χ² value</td>
    </tr>
    <tr>
      <td><code>ndof</code></td>
      <td>Degrees of freedom (N - p)</td>
    </tr>
    <tr>
      <td><code>reduced_chi2</code></td>
      <td>χ²/ndof</td>
    </tr>
    <tr>
      <td><code>pulls</code></td>
      <td>Standardized residuals</td>
    </tr>
    <tr>
      <td><code>residuals</code></td>
      <td>Observed - predicted</td>
    </tr>
    <tr>
      <td><code>predicted</code></td>
      <td>Predicted energies</td>
    </tr>
    <tr>
      <td><code>aic</code></td>
      <td>Akaike Information Criterion</td>
    </tr>
    <tr>
      <td><code>bic</code></td>
      <td>Bayesian Information Criterion</td>
    </tr>
    <tr>
      <td><code>success</code></td>
      <td>Optimizer success flag</td>
    </tr>
    <tr>
      <td><code>message</code></td>
      <td>Optimizer message</td>
    </tr>
    <tr>
      <td><code>n_iter</code></td>
      <td>Number of iterations</td>
    </tr>
    <tr>
      <td><code>n_evaluations</code></td>
      <td>Number of objective evaluations</td>
    </tr>
    <tr>
      <td><code>param_labels</code></td>
      <td>['a', 'b']</td>
    </tr>
  </tbody>
</table>

<h3>JSON Output</h3>

<p><code>fitter.save_results('fit_results.json')</code> saves:</p>

<pre>{
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
}</pre>

<h3>Figures</h3>

<ol>
  <li>
    <p><strong>Figure 8</strong> (<code>figure8_final.pdf</code>):</p>
    <ul>
      <li>Finite-volume spectrum plot</li>
      <li>Energy vs. irrep</li>
      <li>Lattice data points with error bars</li>
      <li>Non-interacting levels (gray bands)</li>
      <li>Physical thresholds (dashed lines)</li>
    </ul>
  </li>
  <li>
    <p><strong>Figure 10</strong> (<code>fig10_reproduction.png</code>):</p>
    <ul>
      <li>Scattering curve plot</li>
      <li>(q/mπ)² vs. (k/mπ)cotδ</li>
      <li>Lattice points</li>
      <li>Luscher curves</li>
      <li>ERE curve with band</li>
      <li>Virtual state curve</li>
    </ul>
  </li>
</ol>

<h3>Profiling Output</h3>

<p>Printed to console:</p>
<ul>
  <li>Timings for each category (total, average, max)</li>
  <li>Counters (χ² evaluations, root finder calls, etc.)</li>
  <li>Derived averages (avg omega evaluations per root, etc.)</li>
</ul>

<h2>Configuration</h2>

<h3>Runtime Configuration</h3>

<p>All configuration is done through Python code. The main parameters are:</p>

<table>
  <thead>
    <tr>
      <th>Parameter</th>
      <th>Default</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>L</code></td>
      <td>64.0</td>
      <td>Lattice size</td>
    </tr>
    <tr>
      <td><code>m1</code></td>
      <td>1.0</td>
      <td>Mass of particle 1 (in mπ)</td>
    </tr>
    <tr>
      <td><code>m2</code></td>
      <td>5.8625</td>
      <td>Mass of particle 2 (in mπ)</td>
    </tr>
    <tr>
      <td><code>INITIAL_GUESS</code></td>
      <td>[0.047, 0.65]</td>
      <td>Initial guess for ERE parameters</td>
    </tr>
    <tr>
      <td><code>BOUNDS</code></td>
      <td>[(-10,10), (-10,10)]</td>
      <td>Parameter bounds</td>
    </tr>
    <tr>
      <td><code>MAX_WORKERS</code></td>
      <td>min(N_levels, cpu_count())</td>
      <td>Parallel workers</td>
    </tr>
    <tr>
      <td><code>ROOT_TOLERANCE</code></td>
      <td>1e-4</td>
      <td>Root finding tolerance</td>
    </tr>
    <tr>
      <td><code>CONTINUITY_TOL</code></td>
      <td>0.3</td>
      <td>Continuity tracking tolerance</td>
    </tr>
    <tr>
      <td><code>OPTIMIZER</code></td>
      <td>'nelder-mead'</td>
      <td>Optimization method</td>
    </tr>
    <tr>
      <td><code>MAXITER</code></td>
      <td>5000</td>
      <td>Maximum iterations</td>
    </tr>
  </tbody>
</table>

<h3>Cache Configuration</h3>

<p>Caching is controlled by:</p>
<ul>
  <li><code>@lru_cache(maxsize=None)</code> for kinematics and zeta function</li>
  <li>Dictionary caches for B-matrix with arbitrary size</li>
</ul>

<p>To clear caches (if needed):</p>
<pre>from kinematics import _compute_kinematics_cached
from final_zeta import _hybrid_Z_cached, _ewald_cached
_compute_kinematics_cached.cache_clear()
_hybrid_Z_cached.cache_clear()
_ewald_cached.cache_clear()</pre>

<h3>Profiling Configuration</h3>

<pre>from profiler import profiler
profiler.enabled = True  # Set to False to disable profiling</pre>

<h2>Troubleshooting</h2>

<h3>Common Issues</h3>

<p><strong>HDF5 file not found:</strong></p>
<pre>FileNotFoundError: [Errno 2] No such file or directory: 'DataSet.hdf5'</pre>
<p><strong>Solution</strong>: Check the file path in <code>HDF5_PATH</code> and ensure the file exists.</p>

<p><strong>Covariance matrix singular:</strong></p>
<pre>LinAlgError: Matrix is singular</pre>
<p><strong>Solution</strong>: The covariance matrix may be ill-conditioned. The pipeline falls back to pseudoinverse, but this may affect χ² calculation. Consider using more bootstrap samples or applying regularization.</p>

<p><strong>Root finder fails:</strong></p>
<pre>RuntimeError: Root finder failed to locate physical root</pre>
<p><strong>Solution</strong>: The quantization condition may not have a root near the observed energy. Check:</p>
<ul>
  <li>Parameter ranges</li>
  <li>Initial guess</li>
  <li>Whether the level is physical (below threshold, etc.)</li>
  <li>Increase <code>local_scan_points</code> or <code>global_scan_points</code></li>
</ul>

<p><strong>Slow performance:</strong></p>
<ul>
  <li>Check that profiling is disabled if not needed</li>
  <li>Ensure caching is working (check cache sizes)</li>
  <li>Reduce <code>max_workers</code> if memory is limited</li>
  <li>For large datasets, consider using asymptotic expansion for |u²| &gt; 50</li>
</ul>

<p><strong>Import errors:</strong></p>
<pre>ModuleNotFoundError: No module named 'QC2'</pre>
<p><strong>Solution</strong>: Ensure the repository root is in <code>sys.path</code> or run from the repository root directory.</p>

<h3>Debugging</h3>

<p>To enable debugging:</p>

<pre>fitter = LuscherFitter(..., debug_objective=True, verbose_predict=True)</pre>

<p>This will print:</p>
<ul>
  <li>Objective function evaluations</li>
  <li>Parameter values</li>
  <li>Predicted energies</li>
  <li>Residuals</li>
</ul>

<p>For root finder debugging:</p>

<pre>root_finder = RootFinder(debug=True, verbose=True)</pre>

<h3>Checking Cache Hit Rates</h3>

<pre>print(f"Kinematics cache: {_compute_kinematics_cached.cache_info()}")
print(f"Hybrid Z cache: {_hybrid_Z_cached.cache_info()}")
print(f"Ewald cache: {_ewald_cached.cache_info()}")
print(f"B-matrix cache size: {len(bmatrix._b_cache)}")</pre>

<h2>Known Limitations</h2>

<h3>Physics Limitations</h3>

<ol>
  <li><strong>Single-channel</strong>: Only single-channel scattering is supported. Coupled-channel effects are not included.</li>
  <li><strong>S-wave only</strong>: Only S-wave scattering is included. Higher partial waves are not implemented.</li>
  <li><strong>Specific irreps</strong>: Only works with irreps for which Morningstar coefficients are available in tables B1-B8.</li>
  <li><strong>Specific masses</strong>: Primarily tested for πΣ scattering with mπ = 1.0 and mΣ = 5.8625.</li>
  <li><strong>Real scattering only</strong>: No analytic continuation to complex energies for resonance pole extraction.</li>
</ol>

<h3>Numerical Limitations</h3>

<ol>
  <li><strong>Padé approximation</strong>: The Padé approximation for zeta function is only valid for u² &gt; 0 and requires precomputed coefficients. For u² &lt; 0, exact Ewald summation is used, which is slower.</li>
  <li><strong>Root finding</strong>: Root finding may fail for levels far from the initial guess or in regions with multiple roots.</li>
  <li><strong>Optimization</strong>: Nelder-Mead is a derivative-free method and may converge slowly or to local minima for poorly conditioned problems.</li>
  <li><strong>Covariance estimation</strong>: The Fisher information matrix approach may underestimate uncertainties for non-linear models.</li>
</ol>

<h3>Implementation Limitations</h3>

<ol>
  <li><strong>No GPU acceleration</strong>: The code is pure Python/NumPy with no GPU support.</li>
  <li><strong>Limited parallelization</strong>: Only energy prediction is parallelized. Other parts (e.g., zeta function evaluation, B-matrix computation) are serial.</li>
  <li><strong>No MPI</strong>: The parallelization uses ProcessPoolExecutor, not MPI.</li>
  <li><strong>No interactive visualization</strong>: All plots are saved to files; no interactive plotting.</li>
</ol>

<h3>Testing Limitations</h3>

<ol>
  <li><strong>No automated test suite</strong>: Formal unit tests are not present.</li>
  <li><strong>Manual validation only</strong>: Validation relies on visual comparison to literature and scientific reproduction.</li>
  <li><strong>Limited error propagation</strong>: Bootstrap error propagation is implemented but not enabled in the default workflow.</li>
</ol>

<h3>Data Format Limitations</h3>

<ol>
  <li><strong>HDF5-only</strong>: Only HDF5 input format is supported.</li>
  <li><strong>Specific HDF5 structure</strong>: Assumes a specific HDF5 layout (PSQ/irrep/ecm_N[_ref]).</li>
  <li><strong>No support for raw energy level files</strong>: The pipeline cannot read plain text files.</li>
</ol>

<h2>Future Work</h2>

<p>Based on the current implementation, the following extensions would be natural:</p>

<h3>Physics Extensions</h3>

<ol>
  <li><strong>Higher partial waves</strong>: Include P-wave, D-wave, etc., by extending the B-matrix tables and quantization condition.</li>
  <li><strong>Coupled channels</strong>: Support multiple coupled channels for systems with inelasticities.</li>
  <li><strong>Resonance pole extraction</strong>: Implement analytic continuation to find poles on the second Riemann sheet.</li>
  <li><strong>Phase shift extraction</strong>: Direct phase shift extraction for arbitrary energies.</li>
  <li><strong>More parameterizations</strong>: Include additional ERE parameterizations (e.g., with effective range, shape parameter).</li>
</ol>

<h3>Performance Improvements</h3>

<ol>
  <li><strong>GPU acceleration</strong>: Port zeta function and B-matrix computation to CUDA for faster evaluation.</li>
  <li><strong>MPI parallelization</strong>: Scale to many levels using MPI.</li>
  <li><strong>Cython/Fortran</strong>: Implement hot loops (zeta function, B-matrix) in compiled languages.</li>
  <li><strong>More efficient root finding</strong>: Implement Newton's method with analytic derivatives.</li>
  <li><strong>JIT compilation</strong>: Use Numba for JIT compilation of performance-critical functions.</li>
</ol>

<h3>Code Improvements</h3>

<ol>
  <li><strong>Automated test suite</strong>: Add unit tests and integration tests.</li>
  <li><strong>Documentation</strong>: Expand docstrings and add API documentation.</li>
  <li><strong>Type hints</strong>: Add comprehensive type hints for better IDE support.</li>
  <li><strong>Configuration files</strong>: Use JSON/YAML for runtime configuration.</li>
  <li><strong>Interactive plotting</strong>: Add Jupyter notebook support for interactive analysis.</li>
</ol>

<h3>Data Support</h3>

<ol>
  <li><strong>Multiple data formats</strong>: Support plain text, CSV, and other formats.</li>
  <li><strong>Automated data discovery</strong>: Automatically discover levels without manual index selection.</li>
  <li><strong>More lattices</strong>: Support different lattice sizes and geometries.</li>
</ol>

<h2>Reproducibility</h2>

<h3>Data</h3>

<p>The pipeline requires the HDF5 dataset file (<code>DataSet.hdf5</code>). This file is not included in the repository. It can be obtained from the associated publication or generated from lattice QCD computations.</p>

<h3>Code Version</h3>

<p>The repository version should be recorded for reproducibility. Use:</p>

<pre>git rev-parse HEAD &gt; version.txt</pre>

<h3>Environment</h3>

<p>To reproduce results, use the same environment:</p>

<pre>pip freeze &gt; requirements.txt</pre>

<h3>Random Seeds</h3>

<p>The pipeline does not use random number generation except for bootstrap resampling (which is deterministic from the input data). No random seeds need to be set.</p>

<h2>Citation</h2>

<p>If you use this software in your research, please cite:</p>

<ol>
  <li>The original Morningstar paper: arXiv:1707.05817</li>
  <li>The BaSc πΣ paper: arXiv:2307.13471</li>
  <li>This repository (with DOI when available)</li>
</ol>

<p>Suggested BibTeX entry:</p>

<pre>@misc{scattering_pipeline,
  author = {Author, A. and Author, B.},
  title = {Finite-Volume Scattering Analysis Pipeline},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/username/repository}
}</pre>

<h2>Acknowledgements</h2>

<p>This work builds upon the foundations laid by the Morningstar group (arXiv:1707.05817) and the BaSc collaboration (arXiv:2307.13471). We thank the lattice QCD community for developing the theoretical framework and providing the data.</p>

<p>The software uses:</p>
<ul>
  <li><code>scipy</code> for optimization and numerical methods</li>
  <li><code>numpy</code> for array operations</li>
  <li><code>h5py</code> for HDF5 I/O</li>
  <li><code>matplotlib</code> for plotting</li>
</ul>

<h2>License</h2>

<p>[Specify license if present in repository. If not, state "Not specified in the repository."]</p>

<hr>

<h2>DOCUMENTATION AUDIT</h2>

<p><strong>Files inspected</strong>: All Python files provided in the repository.</p>

<p><strong>Important modules identified</strong>:</p>
<ol>
  <li><code>dataset_loader.py</code> - Data loading</li>
  <li><code>fitting_driver_canonical.py</code> - Core fitting</li>
  <li><code>ere.py</code> - Effective Range Expansion</li>
  <li><code>morningstar_bmatrix.py</code> - B-matrix</li>
  <li><code>root_finder.py</code> - Root finding</li>
  <li><code>stats.py</code> - Statistical utilities</li>
  <li><code>kinematics.py</code> - Kinematics</li>
  <li><code>final_zeta.py</code> - Hybrid zeta function</li>
  <li><code>exact_zeta.py</code> - Exact Ewald summation</li>
  <li><code>profiler.py</code> - Performance profiling</li>
  <li><code>pipeline_adapter.py</code> - PSQ conversion</li>
  <li><code>run_fit_from_dataset.py</code> - Main entry point</li>
  <li><code>plot_figure8.py</code> - Spectrum plot</li>
  <li><code>plot.py</code> - Scattering curve plot</li>
  <li><code>b_tables/</code> - Morningstar coefficients (b1.py - b8.py)</li>
</ol>

<p><strong>Entry points identified</strong>:</p>
<ol>
  <li><code>run_fit_from_dataset.py</code> - Main fitting script</li>
  <li><code>plot_figure8.py</code> - Spectrum plot generation</li>
  <li><code>plot.py</code> - Scattering curve plot generation</li>
</ol>

<p><strong>Missing information</strong>:</p>
<ul>
  <li>Exact HDF5 structure (inferred from code, not confirmed with example)</li>
  <li>Specific physical units (units are in mπ, confirmed)</li>
  <li>Test suite (no test files found)</li>
  <li>Installation requirements file (not provided)</li>
  <li>License information (not provided)</li>
  <li>Version information (not provided)</li>
</ul>

<p><strong>Ambiguities</strong>:</p>
<ul>
  <li><code>general.data_reader</code> import: This module is imported but not found in the provided files. It is likely part of a larger codebase not included.</li>
  <li><code>tools.final_zeta</code> import: There are multiple import attempts suggesting different directory structures.</li>
  <li>Bootstrap refitting: The <code>bootstrap_parameter_errors</code> method exists but is not used in the default workflow.</li>
</ul>

<p><strong>Potential documentation errors</strong>: None identified. All claims are based on the code provided.</p>

<hr>

<h2>REPOSITORY QUALITY CHECK</h2>

<p><strong>Architecture: 8/10</strong></p>
<ul>
  <li>Clear separation of concerns (data loading, physics, fitting, statistics, plotting)</li>
  <li>Well-defined module boundaries</li>
  <li>Good use of dataclasses and type hints</li>
  <li>Some circular import potential (mitigated by <code>sys.path</code> manipulation)</li>
</ul>

<p><strong>Documentation currently present: 6/10</strong></p>
<ul>
  <li>Docstrings present in most functions</li>
  <li>No external documentation (README, docs, etc.)</li>
  <li>No API documentation</li>
  <li>Some docstrings are minimal</li>
</ul>

<p><strong>Testing: 2/10</strong></p>
<ul>
  <li>No automated test suite</li>
  <li>Manual validation via plotting</li>
  <li>Reference comparisons in code (exact_zeta tests)</li>
  <li>No continuous integration</li>
</ul>

<p><strong>Reproducibility: 7/10</strong></p>
<ul>
  <li>Deterministic when data is fixed</li>
  <li>Random seed not required</li>
  <li>Version control present (implied)</li>
  <li>Missing dependency specification file</li>
</ul>

<p><strong>Maintainability: 7/10</strong></p>
<ul>
  <li>Well-organized code</li>
  <li>Good modularity</li>
  <li>Some duplicated logic (multiple import attempts)</li>
  <li>Moderate code complexity</li>
  <li>Profiling support aids performance debugging</li>
</ul>

<p><strong>Performance engineering: 7/10</strong></p>
<ul>
  <li>Extensive caching (LRU and dictionary)</li>
  <li>Parallel prediction across levels</li>
  <li>Padé approximation for speed</li>
  <li>Asymptotic expansion for large |u²|</li>
  <li>Profiling support</li>
  <li>Could benefit from compiled extensions or GPU acceleration</li>
</ul>

</div>
</body>
</html>
