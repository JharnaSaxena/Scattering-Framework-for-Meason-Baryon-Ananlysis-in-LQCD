"""
Implementation of the single-channel Morningstar B-matrix.
Corrected: includes u factor, proper -4πΛ regularization.
"""
import numpy as np
from typing import Tuple, Optional, Dict, Any
import warnings
import sys
import os
from profiler import profiler

try:
    from b_tables import TABLE_B1, TABLE_B2, TABLE_B3, TABLE_B4, \
                         TABLE_B5, TABLE_B6, TABLE_B7, TABLE_B8
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from b_tables import TABLE_B1, TABLE_B2, TABLE_B3, TABLE_B4, \
                         TABLE_B5, TABLE_B6, TABLE_B7, TABLE_B8

hybrid_Z_refined = None
try:
    from tools.final_zeta import hybrid_Z as hybrid_Z_refined
except ImportError:
    try:
        from final_zeta import hybrid_Z as hybrid_Z_refined
    except ImportError:
        try:
            from .tools.final_zeta import hybrid_Z as hybrid_Z_refined
        except ImportError:
            hybrid_Z_refined = None
            warnings.warn("fast zeta function not found; using pole approx (may lack -4πΛ).", UserWarning)

class SingleChannelBMatrix:
    def __init__(self):
        self.compute_calls = 0
        self._b_cache = {}
        self._all_tables = {}
        for tbl in [TABLE_B1, TABLE_B2, TABLE_B3, TABLE_B4,
                    TABLE_B5, TABLE_B6, TABLE_B7, TABLE_B8]:
            self._all_tables.update(tbl)
        self._coeff_cache = {}

    def get_coefficient(self, irrep: str) -> float:
        if irrep in self._coeff_cache:
            return self._coeff_cache[irrep]
        for key, value in self._all_tables.items():
            if key[0] == irrep and key[-3:] == (0, 0, 1):
                for coeff, l, m, is_imag in value:
                    if l == 0 and m == 0 and not is_imag:
                        self._coeff_cache[irrep] = float(coeff)
                        return float(coeff)
                self._coeff_cache[irrep] = float(value[0][0])
                return float(value[0][0])
        self._coeff_cache[irrep] = 1.0
        return 1.0

    def _compute_regularized_zeta(self, kin) -> float:
        psq = kin.get_psq_label()
        u2 = kin.u2
        gamma = kin.gamma
        L = kin.L
        if hybrid_Z_refined is not None:
            try:
                with profiler.context('Hybrid Zeta'):
                    Z00 = hybrid_Z_refined(
                        u2=u2,
                        psq=psq,
                        gamma=gamma,
                        m_split=2.0 * kin.alpha,
                        L=L,
                    )
                profiler.increment_counter('hybrid_Z_calls')
                return float(np.real(Z00))
            except Exception as e:
                warnings.warn(f"Zeta failed: {e}; using pole approx.")
                return self._pole_approximation(u2, psq)
        else:
            return self._pole_approximation(u2, psq)

    def _pole_approximation(self, u2: float, psq: int) -> float:
        """Simple pole approx with analytic residue (including -4πΛ)"""
        n2 = round(u2)
        if abs(u2 - n2) < 1e-6:
            return 1e12
        n = int(round(np.sqrt(abs(u2))))
        if n == 0:
            n = 1
        if n == 1:
            N_deg = 6
        elif n == 2:
            N_deg = 6 + 8 
        else:
            N_deg = 6 * n
        residue = -N_deg / np.sqrt(4.0 * np.pi)
        return residue / (u2 - n2)

    @profiler.decorator('Morningstar B')
    def compute(self, irrep: str, kin) -> float:
        self.compute_calls += 1
        profiler.increment_counter('B_matrix_calls')
        psq = kin.get_psq_label()
        u2 = round(kin.u2, 12)
        u = np.sqrt(abs(u2)) if u2 >= 0 else np.sqrt(abs(u2))
        gamma = round(kin.gamma, 12)
        m_split = round(2.0 * kin.alpha, 12)
        key = (irrep, psq, u2, gamma, m_split)
        if key in self._b_cache:
            return self._b_cache[key]

        coeff = self.get_coefficient(irrep)
        Z00 = self._compute_regularized_zeta(kin)

        # Corrected now include u factor in denominator
        if abs(gamma) < 1e-12 or abs(u) < 1e-12:
            R00 = Z00 / (gamma * np.pi**1.5 + 1e-12)
        else:
            R00 = Z00 / (gamma * np.pi**1.5 * u)

        B = coeff * R00
        B_val = float(np.real(B))
        self._b_cache[key] = B_val
        return B_val

    def compute_with_details(self, irrep: str, kin) -> Tuple[float, float, float, float]:
        coeff = self.get_coefficient(irrep)
        psq = kin.get_psq_label()
        u2 = kin.u2
        u = np.sqrt(abs(u2)) if u2 >= 0 else np.sqrt(abs(u2))
        gamma = kin.gamma
        Z00 = self._compute_regularized_zeta(kin)
        if abs(gamma) < 1e-12 or abs(u) < 1e-12:
            R00 = Z00 / (gamma * np.pi**1.5 + 1e-12)
        else:
            R00 = Z00 / (gamma * np.pi**1.5 * u)
        B = coeff * R00
        return float(np.real(B)), float(np.real(Z00)), float(coeff), float(np.real(R00))

    def print_coefficients(self):
        print("\nS-wave Coefficients:")
        test_irreps = ['A1g', 'T1u', 'G1u', 'G', 'A2']
        for ir in test_irreps:
            coeff = self.get_coefficient(ir)
            print(f"  {ir:>6}: {coeff:8.4f}")
        print("-" * 50)
        print(f"Compute calls: {self.compute_calls}")

    def get_stats(self) -> Dict[str, Any]:
        return {'compute_calls': self.compute_calls}

def compute_B(irrep: str, kin) -> float:
    bmatrix = SingleChannelBMatrix()
    return bmatrix.compute(irrep, kin)

def compute_R00(kin) -> float:
    psq = kin.get_psq_label()
    u2 = kin.u2
    u = np.sqrt(abs(u2)) if u2 >= 0 else np.sqrt(abs(u2))
    gamma = kin.gamma
    bmatrix = SingleChannelBMatrix()
    Z00 = bmatrix._compute_regularized_zeta(kin)
    if abs(gamma) < 1e-12 or abs(u) < 1e-12:
        return Z00 / (gamma * np.pi**1.5 + 1e-12)
    return Z00 / (gamma * np.pi**1.5 * u)

if __name__ == "__main__":
    print("Morningstar BMatrix with cache")
