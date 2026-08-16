"""
Utilities for locating solutions to the Lüscher quantization condition

The algorithm first searches for a nearby root using adaptive bracketing, If this is unsuccessful, progressively more robust
fallback strategies are used to recover a physically reasonable solution
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from typing import Callable, Optional, List
from profiler import profiler

class RootFinder:
    def __init__(self, xtol=1e-12, rtol=1e-12, maxiter=100,
                 pole_threshold=1e12, root_tolerance=1e-4,
                 local_scan_points=100,
                 global_scan_min=0.1,
                 global_scan_max=50.0, global_scan_points=500,
                 verbose=False, benchmark=False,
                 debug=False, continuity_tol=0.3, observed_tol=0.5):
        self.xtol = xtol
        self.rtol = rtol
        self.maxiter = maxiter
        self.pole_threshold = pole_threshold
        self.root_tolerance = root_tolerance
        self.local_scan_points = local_scan_points
        self.global_scan_min = global_scan_min
        self.global_scan_max = global_scan_max
        self.global_scan_points = global_scan_points
        self.verbose = verbose
        self.benchmark = benchmark
        self.debug = debug
        self.continuity_tol = continuity_tol
        self.observed_tol = observed_tol
        self.last_diagnostics = {}
        self.selection_reason = None

    @profiler.decorator('Root Finder')
    def find_root_near_guess(
        self,
        f: Callable[[float], float],
        x_guess: float,
        prev_root: Optional[float] = None,
        exclude_points: Optional[List[float]] = None,
        tolerance: float = 1e-4,
        reference_energy: Optional[float] = None,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        level_index: Optional[int] = None,
    ) -> float:
        profiler.increment_counter('root_finder_calls')
        diag = {'method': None, 'bracket': None, 'sign_changes': 0,
                'roots_found': [], 'omega_at_root': None, 'converged': False,
                'selection_reason': None}
        self.selection_reason = None

        guess = x_guess
        if guess <= 0:
            guess = max(0.1, x_guess)

        cache = {}
        def f_cached(x):
            if x not in cache:
                cache[x] = f(x)
            return cache[x]

        f_guess = f_cached(guess)
        if abs(f_guess) < self.root_tolerance:
            diag['method'] = 'guess_accept'
            diag['converged'] = True
            diag['omega_at_root'] = f_guess
            self.last_diagnostics = diag
            return guess

        # Determine bounds
        if lower_bound is None:
            lower_bound = 0.1
        if upper_bound is None:
            upper_bound = reference_energy if reference_energy is not None else guess + 0.5
        if upper_bound <= lower_bound:
            upper_bound = lower_bound + 0.5
        # Adaptive bracketing
        widths = [0.02, 0.05, 0.10, 0.20, 0.40, 0.80, 1.5]
        adaptive_brackets = []
        for w in widths:
            a = max(guess - w, 1e-6)
            b = guess + w
            if lower_bound is not None:
                a = max(a, lower_bound)
            if upper_bound is not None:
                b = min(b, upper_bound)
            if b <= a:
                continue
            fa = f_cached(a)
            fb = f_cached(b)
            if not (np.isfinite(fa) and np.isfinite(fb)):
                continue
            if abs(fa) > self.pole_threshold or abs(fb) > self.pole_threshold:
                continue
            if fa * fb < 0:
                adaptive_brackets.append((a, b))

        all_roots = []
        for a, b in adaptive_brackets:
            try:
                root = brentq(f_cached, a, b, xtol=self.xtol, rtol=self.rtol, maxiter=self.maxiter)
                if root is not None and abs(f_cached(root)) < self.root_tolerance:
                    all_roots.append(root)
            except Exception:
                continue
        if all_roots:
            unique = []
            for r in sorted(all_roots):
                if not any(abs(r - u) < 1e-6 for u in unique):
                    unique.append(r)
            all_roots = unique
        # Ordered root selection
        if level_index is not None and all_roots:
            sorted_roots = sorted(all_roots)
            if len(sorted_roots) > level_index:
                chosen = sorted_roots[level_index]
                reason = f'ordered_root_{level_index}'
                diag['method'] = 'ordered'
                diag['converged'] = True
                diag['omega_at_root'] = f_cached(chosen)
                diag['selection_reason'] = reason
                self.selection_reason = reason
                self.last_diagnostics = diag
                if self.debug:
                    print(f"[RootFinder DEBUG] Ordered root selected- {chosen:.6f} (index {level_index})")
                return chosen
            else:
                if self.debug:
                    print(f"[RootFinder DEBUG] Not enough roots found ({len(sorted_roots)}) for level_index {level_index}")
        # Check if any root is acceptable (near prev_root or reference)
        chosen = None
        reason = None
        if all_roots:
            f_vals = [f_cached(r) for r in all_roots]
            # Priority: near prev_root
            if prev_root is not None:
                candidates = [(r, fr) for r, fr in zip(all_roots, f_vals) if abs(r - prev_root) < self.continuity_tol]
                if candidates:
                    best = min(candidates, key=lambda x: abs(x[1]))
                    chosen = best[0]
                    reason = 'near_prev'
            # If no prev_root, or none near it, use reference_energy
            if chosen is None and reference_energy is not None:
                candidates = [(r, fr) for r, fr in zip(all_roots, f_vals) if abs(r - reference_energy) < 0.5]
                if candidates:
                    best = min(candidates, key=lambda x: abs(x[1]))
                    chosen = best[0]
                    reason = 'near_reference'
            # Fallback closest to x_guess
            if chosen is None:
                candidates = [(r, fr) for r, fr in zip(all_roots, f_vals)]
                best = min(candidates, key=lambda x: (abs(x[0] - x_guess), abs(x[1])))
                chosen = best[0]
                reason = 'closest_to_observed'

        # If no acceptable root, run local scan (fallback)
        if chosen is None:
            points = np.linspace(lower_bound, upper_bound, self.local_scan_points)
            f_vals = [f_cached(p) for p in points]
            finite_mask = np.isfinite(f_vals)
            roots = []
            for i in range(len(points)-1):
                if finite_mask[i] and finite_mask[i+1]:
                    if f_vals[i] * f_vals[i+1] < 0:
                        try:
                            root = brentq(f_cached, points[i], points[i+1],
                                          xtol=self.xtol, rtol=self.rtol, maxiter=self.maxiter)
                            if root is not None and abs(f_cached(root)) < self.root_tolerance:
                                roots.append(root)
                        except Exception:
                            continue
            if roots:
                unique_roots = []
                for r in sorted(roots):
                    if not any(abs(r - u) < 1e-6 for u in unique_roots):
                        unique_roots.append(r)
                roots = unique_roots

                # Ordered root selection for scan
                if level_index is not None and roots:
                    sorted_roots = sorted(roots)
                    if len(sorted_roots) > level_index:
                        chosen = sorted_roots[level_index]
                        reason = f'ordered_scan_{level_index}'
                        diag['method'] = 'local_scan_ordered'
                        diag['converged'] = True
                        diag['omega_at_root'] = f_cached(chosen)
                        diag['selection_reason'] = reason
                        self.selection_reason = reason
                        self.last_diagnostics = diag
                        if self.debug:
                            print(f"[RootFinder DEBUG] Ordered scan root: {chosen:.6f} (index {level_index})")
                        return chosen
                # Fallback to closest
                if prev_root is not None:
                    root = min(roots, key=lambda r: abs(r - prev_root))
                    reason = 'near_prev_scan'
                elif reference_energy is not None:
                    root = min(roots, key=lambda r: abs(r - reference_energy))
                    reason = 'near_reference_scan'
                else:
                    root = min(roots, key=lambda r: abs(r - x_guess))
                    reason = 'closest_scan'
                chosen = root
                diag['method'] = 'local_scan'
            else:
                # Ultimate fallback: minimize |f|
                def abs_f(x):
                    return abs(f_cached(x))
                try:
                    res = minimize_scalar(abs_f, bounds=(lower_bound, upper_bound), method='bounded')
                    if res.fun < self.root_tolerance:
                        chosen = res.x
                        reason = 'min_abs'
                        diag['method'] = 'min_abs'
                    else:
                        # No root found – raise an error instead of returning guess
                        raise RuntimeError(
                            f"Root finder failed to locate physical root for f(x) near "
                            f"x_guess={guess:.6f}. Check bracket or function behaviour."
                        )
                except Exception as e:
                    # Re-raise any unexpected error (including our own RuntimeError)
                    raise RuntimeError(
                        f"Root finder failed: {e}\n"
                        f"x_guess={guess:.6f}, lower_bound={lower_bound}, upper_bound={upper_bound}"
                    ) from e

        if self.debug:
            print(f"\n[RootFinder DEBUG] x_guess={x_guess:.6f}, prev_root={prev_root if prev_root is None else prev_root:.6f}")
            print(f"[RootFinder DEBUG] All roots found: {[f'{r:.6f}' for r in all_roots]}")
            print(f"[RootFinder DEBUG] Chosen: {chosen:.6f}, reason: {reason}")

        diag['converged'] = True
        diag['omega_at_root'] = f_cached(chosen)
        diag['selection_reason'] = reason
        self.selection_reason = reason
        self.last_diagnostics = diag
        return chosen

    def _scan_and_solve(self, f, center, n_points, exclude, tol):
        pass

    def _global_scan(self, f, guess, exclude, tol):
        pass

    @staticmethod
    def _not_excluded(root, exclude, tol):
        if exclude is None: return True
        return not any(abs(root - x) < tol for x in exclude)

def find_root_near_guess(f, x_guess, search_width=0.1, exclude_points=None,
                         tolerance=1e-4, reference_energy=None, level_index=None):
    finder = RootFinder()
    return finder.find_root_near_guess(f, x_guess, exclude_points=exclude_points,
                                       tolerance=tolerance, reference_energy=reference_energy,
                                       level_index=level_index)

if __name__ == "__main__":
    def f(x): return (x-1)*(x-10)
    fr = RootFinder(debug=True)
    root = fr.find_root_near_guess(f, 5.0, reference_energy=10.0, level_index=0)
    print("Root closest to 10:", root)
