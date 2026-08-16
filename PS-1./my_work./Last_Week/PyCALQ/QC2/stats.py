"""
Statistical utilities used throughout the fitting pipeline
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve, LinAlgError
from typing import Optional, Tuple, List, Dict, Union
from dataclasses import dataclass
from profiler import profiler

def _solve_spd(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    try:
        c, low = cho_factor(matrix, lower=True)
        return cho_solve((c, low), rhs)
    except LinAlgError:
        pass
    try:
        return np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(matrix) @ rhs

@profiler.decorator('Covariance')
def bootstrap_covariance(bootstrap_samples: np.ndarray) -> np.ndarray:
    samples = np.asarray(bootstrap_samples)
    if samples.ndim == 1:
        return np.array([[np.var(samples[1:], ddof=1)]])
    replicas = samples[:, 1:]
    return np.cov(replicas)

@profiler.decorator('Chi2')
def chi2(observed_mean: np.ndarray, predicted: np.ndarray, covariance: np.ndarray) -> float:
    observed_mean = np.asarray(observed_mean, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if len(observed_mean) != len(predicted):
        raise ValueError(f"Length mismatch: {len(observed_mean)} vs {len(predicted)}")
    residual = observed_mean - predicted
    x = _solve_spd(covariance, residual)
    return float(residual.T @ x)

@profiler.decorator('Parameter Covariance')
def parameter_covariance(jacobian: np.ndarray, data_cov: np.ndarray,
                         chi2_value: Optional[float] = None,
                         n_data: Optional[int] = None,
                         n_params: Optional[int] = None) -> np.ndarray:
    jacobian = np.asarray(jacobian, dtype=float)
    data_cov = np.asarray(data_cov, dtype=float)
    X = _solve_spd(data_cov, jacobian)
    fisher = jacobian.T @ X
    if chi2_value is not None and n_data is not None and n_params is not None:
        dof = n_data - n_params
        if dof > 0:
            reduced = chi2_value / dof
            if reduced > 1.0:
                fisher = fisher / reduced
    try:
        return np.linalg.inv(fisher)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(fisher)

def parameter_errors(cov_par: np.ndarray) -> np.ndarray:
    return np.sqrt(np.abs(np.diag(cov_par)))

def reduced_chi2(chi2_value: float, n_data: int, n_params: int) -> float:
    dof = n_data - n_params
    return chi2_value / dof if dof > 0 else float('inf')

def standardized_residuals(observed_mean: np.ndarray, predicted: np.ndarray,
                            covariance: np.ndarray) -> np.ndarray:
    residual = np.asarray(observed_mean) - np.asarray(predicted)
    diag = np.diag(np.asarray(covariance))
    return residual / np.sqrt(diag)

def correlation_matrix(covariance: np.ndarray) -> np.ndarray:
    cov = np.asarray(covariance, dtype=float)
    diag = np.sqrt(np.diag(cov))
    corr = cov / np.outer(diag, diag)
    return np.clip(corr, -1.0, 1.0)

def pull_distribution(observed_mean: np.ndarray, predicted: np.ndarray,
                       covariance: np.ndarray) -> Dict[str, float]:
    pulls = standardized_residuals(observed_mean, predicted, covariance)
    valid = pulls[np.isfinite(pulls)]
    if len(valid) == 0:
        return {'mean': float('nan'), 'std': float('nan'), 'max': float('nan'),
                'min': float('nan'), 'rms': float('nan'), 'n_valid': 0}
    return {
        'mean': float(np.mean(valid)),
        'std': float(np.std(valid)),
        'max': float(np.max(valid)),
        'min': float(np.min(valid)),
        'rms': float(np.sqrt(np.mean(valid ** 2))),
        'n_valid': len(valid),
    }

def aic(chi2_value: float, n_params: int, n_data: int) -> float:
    return chi2_value + 2 * n_params

def aicc(chi2_value: float, n_params: int, n_data: int) -> float:
    aic_val = chi2_value + 2 * n_params
    dof = n_data - n_params - 1
    if dof > 0:
        aic_val += (2 * n_params * (n_params + 1)) / dof
    return aic_val

def bic(chi2_value: float, n_params: int, n_data: int) -> float:
    return chi2_value + n_params * np.log(n_data)

@dataclass
class FitResult:
    params: np.ndarray
    chi2: float
    ndof: int
    reduced_chi2: float
    cov_params: np.ndarray
    errors: np.ndarray
    corr_params: np.ndarray
    pulls: np.ndarray
    residuals: np.ndarray
    predicted: np.ndarray
    aic: float
    bic: float
    success: bool
    message: str
    n_iter: int
    n_evaluations: int
    param_labels: List[str] = None

    def to_dict(self) -> Dict:
        return {
            'params': self.params.tolist(),
            'chi2': float(self.chi2),
            'ndof': self.ndof,
            'reduced_chi2': float(self.reduced_chi2),
            'cov_params': self.cov_params.tolist(),
            'errors': self.errors.tolist(),
            'corr_params': self.corr_params.tolist(),
            'pulls': self.pulls.tolist(),
            'residuals': self.residuals.tolist(),
            'predicted': self.predicted.tolist(),
            'aic': float(self.aic),
            'bic': float(self.bic),
            'success': self.success,
            'message': self.message,
            'n_iter': self.n_iter,
            'n_evaluations': self.n_evaluations,
            'param_labels': self.param_labels,
        }

if __name__ == "__main__":
    print("stats.py")
