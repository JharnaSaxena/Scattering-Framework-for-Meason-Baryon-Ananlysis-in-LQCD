"""
Effective Range Expansion implementation,
This module follows Eq.12 of arXiv:2307.13471 for the single channel πΣ scattering analysis.
The ERE is parameterized as
(k/mπ) cotdelta = (E/mπ) [a + bdelta]
where
delta = (E^2- threshold^2) / threshold^2
"""
import numpy as np
from typing import List, Union, Protocol
from dataclasses import dataclass
from profiler import profiler

class KinematicVarsLike(Protocol):
    q_star_squared: float
    E_cm: float
    threshold: float

@dataclass
class ERECoefficients:
    coeffs: np.ndarray
    labels: List[str]
    def __post_init__(self):
        self.coeffs = np.asarray(self.coeffs, dtype=float)
    @property
    def n_params(self): return len(self.coeffs)
    def __getitem__(self, idx): return self.coeffs[idx]
    def __len__(self): return len(self.coeffs)

class ERE:
    def __init__(self, coeffs: Union[List[float], np.ndarray]):
        self.coeffs = np.asarray(coeffs, dtype=float)
        self.n_params = len(self.coeffs)
        if self.n_params == 0:
            raise ValueError("ERE requires at least one coefficient")
        if self.n_params > 2:
            raise ValueError("ERE currently supports up to 2 parameters (a, b)")
        if np.any(~np.isfinite(self.coeffs)):
            raise ValueError("ERE coefficients contain NaN or Inf")

    @profiler.decorator('ERE')
    def compute_kinv(self, kin: Union[KinematicVarsLike, float]) -> float:
        profiler.increment_counter('ERE_calls')
        if isinstance(kin, (float, int, np.floating, np.integer)):
            q2 = kin
            kinv = self.coeffs[0] + (self.coeffs[1] * q2 if self.n_params > 1 else 0.0)
            return float(kinv)
        if hasattr(kin, "E_cm") and hasattr(kin, "threshold"):
            E = kin.E_cm
            th = kin.threshold
            # Dimensionless energy variable defined in Eq.12
            Delta = (E*E - th*th) / (th*th + 1e-15)
            kinv = E * (self.coeffs[0] + self.coeffs[1] * Delta)
            return float(kinv)
        else:
            q2 = self._extract_q_star_squared(kin)
            kinv = self.coeffs[0] + (self.coeffs[1] * q2 if self.n_params > 1 else 0.0)
            return float(kinv)

    def _extract_q_star_squared(self, kin):
        if isinstance(kin, (float, int, np.floating, np.integer)):
            return float(kin)
        if hasattr(kin, "q_star_squared"):
            return float(kin.q_star_squared)
        if hasattr(kin, "q2_star"):
            return float(kin.q2_star)
        raise TypeError("Cannot extract q² from kin")

    def compute_cot_delta(self, kin):
        q2 = self._extract_q_star_squared(kin)
        kinv = self.compute_kinv(kin)
        if q2 >= 0:
            q_star = np.sqrt(q2)
            return kinv / (q_star + 1e-15)
        else:
            kappa = np.sqrt(-q2)
            return kinv / (kappa + 1e-15)

    def compute_phase_shift(self, kin):
        q2 = self._extract_q_star_squared(kin)
        if q2 < 0: return 0.0
        q_star = np.sqrt(q2)
        kinv = self.compute_kinv(kin)
        if abs(kinv) < 1e-15:
            return 90.0
        delta = np.arctan2(q_star, kinv)
        return np.degrees(delta)

    def get_coeffs(self) -> ERECoefficients:
        labels = ['a'] if self.n_params >= 1 else []
        if self.n_params > 1: labels.append('b')
        return ERECoefficients(self.coeffs.copy(), labels)

    def _get_coeff_labels(self):
        labels = ['a']
        if self.n_params > 1: labels.append('b')
        return labels

    def __call__(self, kin): return self.compute_kinv(kin)

    def __repr__(self):
        labels = self._get_coeff_labels()
        terms = [f"{labels[i]} = {self.coeffs[i]:.6f}" for i in range(len(self.coeffs))]
        return f"ERE({', '.join(terms)})"

def standard_ere(a0=None, r0=None): return ERE([0.0, 0.0])
def constant_ere(c0=0.0): return ERE([c0])
def ere_from_coeffs(coeffs): return ERE(coeffs)

if __name__ == "__main__":
    print("ere.py")
