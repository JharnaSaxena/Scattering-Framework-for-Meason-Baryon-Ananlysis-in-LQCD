"""
HYBRID ZETA FUNCTION with LRU caching for both Exact Ewald and Hybrid Zeta.
u^2 < 0 uses exact Ewald (if available) or asymptotic for large |u²|
u^2 >= 0 uses Pade if m_split=1, else exact Ewald or asymptotic
"""

import numpy as np
import json
import os
import warnings
import sys
from functools import lru_cache
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from .exact_zeta import Z as exact_Z
    print("DEBUG: exact_zeta imported successfully from tools/")
except ImportError:
    try:
        from exact_zeta import Z as exact_Z
        print("DEBUG: exact_zeta imported from current directory")
    except ImportError:
        warnings.warn("exact_zeta.py not found. Will use asymptotic for all u².")
        exact_Z = None

try:
    from profiler import profiler
except ImportError:
    class DummyProfiler:
        @staticmethod
        def decorator(name):
            def wrapper(func):
                return func
            return wrapper
    profiler = DummyProfiler()

PSQ_D = {'PSQ0': (0,0,0), 'PSQ1': (0,0,1), 'PSQ2': (1,1,0), 'PSQ3': (1,1,1), 'PSQ4': (0,0,2)}

COEFF_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'coefficients')

def load_coeffs_for_gamma(psq, gamma):
    gamma_str = f"gamma{gamma:.2f}".replace('.', '_')
    path = os.path.join(COEFF_DIR, f'{psq}_{gamma_str}_coeffs.json')
    if not os.path.exists(path):
        return None
    with open(path, 'r') as f:
        data = json.load(f)
    return data.get('intervals', None)

def pade_func(x, coeffs, order):
    n = order
    m = order
    a = coeffs[:n+1]
    b = [1.0] + coeffs[n+1:]
    P = sum(a[i] * x**i for i in range(n+1))
    Q = sum(b[i] * x**i for i in range(m+1))
    return P / Q

@profiler.decorator('Padé')
def eval_pade_interval(u2, interval):
    a = interval['a']; b = interval['b']
    if u2 < a or u2 > b:
        return None
    mid = interval['mid']; half = interval['half']
    x = (u2 - mid) / half
    coeffs = interval['coeffs']; order = interval['order']
    val = pade_func(x, coeffs, order)
    residues = interval.get('residues', {'left': 0.0, 'right': 0.0})
    if abs(u2 - a) > 1e-10:
        val += residues['left'] / (u2 - a)
    if abs(u2 - b) > 1e-10:
        val += residues['right'] / (u2 - b)
    return val

@profiler.decorator('Padé')
def pade_approximation(u2, psq, gamma):
    intervals = load_coeffs_for_gamma(psq, gamma)
    if not intervals:
        return None
    for interval in intervals:
        if interval['a'] <= u2 <= interval['b']:
            val = eval_pade_interval(u2, interval)
            if val is not None:
                return val
    return None

# ----- CACHED EXACT EWALD -----
@lru_cache(maxsize=None)
def _ewald_cached(u2, psq, gamma, m_split, L):
    """
    Cached version of exact Ewald summation.
    All arguments are rounded in the caller to avoid floating-point noise.
    """
    if exact_Z is None:
        return asymptotic_Z00(u2, psq, gamma, m_split, L)
    d = np.array(PSQ_D[psq], dtype=float)
    try:
        # Use exact Ewald sum
        z = exact_Z(u2, gamma=gamma, l=0, m=0, d=d, m_split=m_split, precision=1e-8)
        return float(np.real(z))
    except Exception as e:
        warnings.warn(f"Ewald failed: {e}; using asymptotic")
        return asymptotic_Z00(u2, psq, gamma, m_split, L)

@profiler.decorator('Exact Ewald')
def ewald_approximation(u2, psq, gamma, m_split=1.0, L=64.0):
    """
    Wrapper for exact Ewald with caching.
    Rounds arguments to 10 decimal places to maximize cache hits.
    """
    # Round to stable keys
    key = (round(u2, 10), psq, round(gamma, 10), round(m_split, 10), round(L, 10))
    return _ewald_cached(*key)

def asymptotic_Z00(u2, psq, gamma, m_split, L=64.0):
    """Asymptotic expansion for large |u²|"""
    d = np.array(PSQ_D[psq], dtype=float)
    d2 = np.sum(d**2)
    term1 = 1.0 / u2
    term2 = (m_split**2 / 4.0 - 1.0 / gamma**2 + (4.0 * np.pi**2 * d2) / (3.0 * gamma**2 * L**2)) / (2.0 * u2**2)
    return (term1 + term2) / (gamma * np.pi**1.5)

@lru_cache(maxsize=None)
def _hybrid_Z_cached(u2, psq, gamma, m_split, L):
    """
    Cached hybrid zeta function core logic, all arguments are already rounded before calling.
    """
    # For very large |u²|, use asymptotic directly (faster and stable)
    if abs(u2) > 50.0:
        return asymptotic_Z00(u2, psq, gamma, m_split, L)

    if u2 < 0:
        return ewald_approximation(u2, psq, gamma, m_split, L)
    else:
        if abs(m_split - 1.0) < 1e-12:
            val = pade_approximation(u2, psq, gamma)
            if val is not None and np.isfinite(val):
                return val
        return ewald_approximation(u2, psq, gamma, m_split, L)

@profiler.decorator('Hybrid Zeta')
def hybrid_Z(u2, psq, gamma=1.0, m_split=1.0, L=64.0):
    """
    Hybrid zeta function with caching rounds arguments to 10 decimal places to maximize cache hits.
    """
    # Round to stable keys
    key = (round(u2, 10), psq, round(gamma, 10), round(m_split, 10), round(L, 10))
    return _hybrid_Z_cached(*key)

def hybrid_Z_from_kinematics(kin):
    psq = kin.get_psq_label()
    return hybrid_Z(kin.u2, psq, gamma=kin.gamma, m_split=2*kin.alpha, L=kin.L)

if __name__ == "__main__":
    print("HYBRID ZETA with LRU caching for both Ewald and Hybrid Zeta")
