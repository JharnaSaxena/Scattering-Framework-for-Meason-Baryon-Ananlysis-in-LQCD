#!/usr/bin/env python3
"""
Kinematics module – dimensionless finite volume scattering.
All quantities are dimensionless with respect to Mref, L is in lattice units.
Uses functools.lru_cache for caching with rounded arguments.
"""
import numpy as np
import math
from typing import Tuple, Optional, Dict, Any
from dataclasses import dataclass
import warnings
import functools

@functools.lru_cache(maxsize=None)
def _compute_kinematics_cached(
    E_cm_rounded: float,
    d: Tuple[int, int, int],
    m1_rounded: float,
    m2_rounded: float,
    L_rounded: float,
    Mref_rounded: float,
) -> 'KinematicVars':
    # Reconstruct rounded values (they are the same as inputs)
    E_cm = E_cm_rounded
    m1 = m1_rounded
    m2 = m2_rounded
    L = L_rounded
    Mref = Mref_rounded
    d_arr = np.array(d, dtype=float)
    # Total momentum (dimensionless)
    P = d_arr * (2.0 * np.pi / L)
    P2 = float(np.sum(P * P))
    # Lab frame energy (dimensionless)
    E_lab = np.sqrt(max(0, E_cm * E_cm + P2))
    # Lorentz boost factor
    gamma = E_lab / (E_cm + 1e-15)
    # Threshold (dimensionless)
    threshold = m1 + m2
    # Mass asymmetry parameter (unequal masses)
    if E_cm > 0:
        alpha = 0.5 * (1.0 + (m1 * m1 - m2 * m2) / (E_cm * E_cm + 1e-15))
    else:
        alpha = 0.5
    # Kallen function (dimensionless)
    s = E_cm * E_cm
    m1_2 = m1 * m1
    m2_2 = m2 * m2
    lam = s * s + m1_2 * m1_2 + m2_2 * m2_2 - 2 * s * m1_2 - 2 * s * m2_2 - 2 * m1_2 * m2_2
    # CM momentum (dimensionless)
    q2_star = lam / (4.0 * s + 1e-15)
    if q2_star < 0:
        q_star = -np.sqrt(max(0, -q2_star))
        q2_star = -np.abs(q2_star)
    else:
        q_star = np.sqrt(q2_star)
    # Dimensionless momentum u
    u = L * q_star / (2.0 * np.pi)
    u2 = u * u
    # BaSc Delta variable
    Delta = (E_cm * E_cm - threshold * threshold) / (threshold * threshold + 1e-15)
    # Build the KinematicVars object (frozen dataclass)
    kin = KinematicVars(
        E_cm=E_cm,
        m1=m1,
        m2=m2,
        L=L,
        Mref=Mref,
        d=d,
        P=P,
        P2=P2,
        E_lab=E_lab,
        gamma=gamma,
        alpha=alpha,
        threshold=threshold,
        Delta=Delta,
        q_star=q_star,
        q2_star=q2_star,
        u=u,
        u2=u2,
    )
    return kin
def compute_kinematics(
    E_cm: float,
    d: Tuple[int, int, int],
    m1: float,
    m2: float,
    L: float,
    Mref: float = 1.0,
) -> 'KinematicVars':
    """
    Compute all dimensionless kinematic quantities, with caching.

    Argiments:
    E_cm: Center-of-mass energy / Mref
    d: Momentum vector (integer multiples of 2π/L)
    m1: Mass of particle 1 / Mref
    m2: Mass of particle 2 / Mref
    L: Lattice size (lattice units)
    Mref: Reference mass (physical units)

    Returns: KinematicVars object with all quantities dimensionless.
    """
    key = (round(E_cm, 10), d, round(m1, 10), round(m2, 10), round(L, 10), round(Mref, 10))
    return _compute_kinematics_cached(*key)


# Helper functions (unchanged they now call the cached version)
def compute_kinematics_from_physical(
    E_cm_phys: float,
    d: Tuple[int, int, int],
    m1_phys: float,
    m2_phys: float,
    L: float,
    Mref: float = 1.0,
) -> 'KinematicVars':
    E_cm = E_cm_phys / Mref
    m1 = m1_phys / Mref
    m2 = m2_phys / Mref
    return compute_kinematics(E_cm, d, m1, m2, L, Mref)

def q2_cm(E_cm: float, m1: float, m2: float) -> float:
    #alias for q2_star for compatibility
    kin = compute_kinematics(E_cm, (0, 0, 0), m1, m2, L=48.0)
    return kin.q2_star

def q_cm(E_cm: float, m1: float, m2: float) -> float:
    #alias for q_star for compatibility
    kin = compute_kinematics(E_cm, (0, 0, 0), m1, m2, L=48.0)
    return kin.q_star

def alpha(E_cm: float, m1: float, m2: float) -> float:
    #compute alpha for unequal masses
    if E_cm <= 0:
        return 0.5
    return 0.5 * (1.0 + (m1 * m1 - m2 * m2) / (E_cm * E_cm + 1e-15))

def gamma(E_lab: float, E_cm: float) -> float:
    #compute Lorentz boost factor
    return E_lab / (E_cm + 1e-15)
# KinematicVars dataclass (unchanged, but must be defined before use)
@dataclass(frozen=True)
class KinematicVars:
    E_cm: float
    m1: float
    m2: float
    L: float
    Mref: float
    d: Tuple[int, int, int]
    P: np.ndarray
    P2: float
    E_lab: float
    gamma: float
    alpha: float
    threshold: float
    Delta: float
    q_star: float
    q2_star: float
    u: float
    u2: float

    def __post_init__(self):
        if self.E_cm <= 0:
            raise ValueError(f"E_cm must be positive, got {self.E_cm}")
        if self.L <= 0:
            raise ValueError(f"L must be positive, got {self.L}")
        if self.gamma < 1.0 - 1e-8:
            raise ValueError(
                f"gamma = {self.gamma} << 1.0. E_lab must be >= E_cm. "
                f"Check E_cm and P²."
            )

    @property
    def q_star_squared(self):
        return self.q2_star

    def is_below_threshold(self) -> bool:
        return self.E_cm < self.threshold

    def is_at_threshold(self) -> bool:
        return abs(self.E_cm - self.threshold) < 1e-10

    def get_psq_label(self) -> str:
        d2 = int(np.sum(np.array(self.d) ** 2))
        return f"PSQ{d2}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            'E_cm': float(self.E_cm),
            'm1': float(self.m1),
            'm2': float(self.m2),
            'L': float(self.L),
            'Mref': float(self.Mref),
            'd': tuple(self.d),
            'P': self.P.tolist() if hasattr(self.P, 'tolist') else list(self.P),
            'P2': float(self.P2),
            'E_lab': float(self.E_lab),
            'gamma': float(self.gamma),
            'alpha': float(self.alpha),
            'threshold': float(self.threshold),
            'Delta': float(self.Delta),
            'q_star': float(self.q_star),
            'q2_star': float(self.q2_star),
            'u': float(self.u),
            'u2': float(self.u2),
            'psq': self.get_psq_label(),
        }
# Self test
if __name__ == "__main__":
    print("KINEMATICS MODULE with functools.lru_cache")
    
