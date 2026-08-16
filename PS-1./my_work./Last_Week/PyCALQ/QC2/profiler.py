"""
Simple profiling utility for the lattice scattering pipeline
The profiler records execution times for different parts of the code and keeps track of useful counters such as the number of χ^2 evaluations, root-finder calls and B-matrix evaluations.
Profiling can be enabled or disabled using the global flag below
"""
import time
import functools
from collections import defaultdict
ENABLE_PROFILING = True

class Profiler:
    #central profiler that records timing and counters

    def __init__(self):
        self.enabled = ENABLE_PROFILING
        self.stats = defaultdict(lambda: {
            'calls': 0,
            'total': 0.0,
            'max': 0.0,
            'times': [],
        })
        self.counters = {
            'optimizer_iterations': 0,
            'chi2_evaluations': 0,
            'omega_evaluations': 0,
            'root_finder_calls': 0,
            'root_iterations': 0,
            'root_solver_attempts': 0,
            'hybrid_Z_calls': 0,
            'B_matrix_calls': 0,
            'ERE_calls': 0,
            'pade_evaluations': 0,
            'exact_ewald_evaluations': 0,
        }
        self.current_timers = {}

    def start(self, category):
        #Start a timer for the given category
        if not self.enabled:
            return
        self.current_timers[category] = time.perf_counter()

    def stop(self, category):
        #Stop the timer and record the elapsed time
        if not self.enabled:
            return
        if category in self.current_timers:
            elapsed = time.perf_counter() - self.current_timers[category]
            del self.current_timers[category]
            self.record(category, elapsed)

    def record(self, category, elapsed):
        #Record a single timing measurement
        if not self.enabled:
            return
        stat = self.stats[category]
        stat['calls'] += 1
        stat['total'] += elapsed
        if elapsed > stat['max']:
            stat['max'] = elapsed
        stat['times'].append(elapsed)

    def increment_counter(self, name, amount=1):
        #Increment a named counter
        if not self.enabled:
            return
        if name in self.counters:
            self.counters[name] += amount

    def context(self, category):
        #Return a context manager for timing a block
        class Context:
            def __enter__(self_ctx):
                self.start(category)
            def __exit__(self_ctx, *args):
                self.stop(category)
        return Context()

    def decorator(self, category):
        #Return a decorator that times the decorated function
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                if not self.enabled:
                    return func(*args, **kwargs)
                self.start(category)
                try:
                    return func(*args, **kwargs)
                finally:
                    self.stop(category)
            return wrapper
        return decorator

    def report(self):
        #Print a formatted profiling report
        if not self.enabled:
            print("\nProfiling is disabled. Set ENABLE_PROFILING = True to collect data.")
            return

        print("PIPELINE PROFILING REPORT")
        total_time = sum(stat['total'] for stat in self.stats.values())
        rows = []
        for cat, stat in sorted(self.stats.items()):
            calls = stat['calls']
            total = stat['total']
            avg_ms = (total / calls * 1000) if calls else 0.0
            max_ms = stat['max'] * 1000
            pct = (total / total_time * 100) if total_time > 0 else 0.0
            rows.append((cat, calls, total, avg_ms, max_ms, pct))
        print(f"\n{'Category':<30} {'Calls':>8} {'Total (s)':>12} "
              f"{'Avg (ms)':>12} {'Max (ms)':>12} {'%':>8}")
        for row in rows:
            print(f"{row[0]:<30} {row[1]:>8} {row[2]:>12.6f} {row[3]:>12.3f} "
                  f"{row[4]:>12.3f} {row[5]:>8.1f}")

        print(f"{'TOTAL':<30} {'':>8} {total_time:>12.6f} {'':>12} {'':>12} {'':>8}")
        print("COUNTERS")
        for name, value in self.counters.items():
            print(f"{name:>30}: {value}")

        # Derive some averages
        if self.counters['root_finder_calls'] > 0:
            avg_omega_per_root = self.counters['omega_evaluations'] / self.counters['root_finder_calls']
            avg_root_iter = self.counters['root_iterations'] / self.counters['root_finder_calls']
            print(f"{'avg omega evaluations per root':>30}: {avg_omega_per_root:.2f}")
            print(f"{'avg root iterations per root':>30}: {avg_root_iter:.2f}")

        if self.counters['chi2_evaluations'] > 0:
            avg_hybrid_per_chi2 = self.counters['hybrid_Z_calls'] / self.counters['chi2_evaluations']
            avg_B_per_chi2 = self.counters['B_matrix_calls'] / self.counters['chi2_evaluations']
            avg_ERE_per_chi2 = self.counters['ERE_calls'] / self.counters['chi2_evaluations']
            print(f"{'avg hybrid_Z calls per χ² eval':>30}: {avg_hybrid_per_chi2:.2f}")
            print(f"{'avg B-matrix calls per χ² eval':>30}: {avg_B_per_chi2:.2f}")
            print(f"{'avg ERE calls per χ² eval':>30}: {avg_ERE_per_chi2:.2f}")
        if 'Omega' in self.stats:
            omega_stat = self.stats['Omega']
            omega_calls = omega_stat['calls']
            omega_total = omega_stat['total']
            if omega_calls > 0:
                avg_omega_time_ms = (omega_total / omega_calls) * 1000
                print(f"{'avg time per Omega (ms)':>30}: {avg_omega_time_ms:.3f}")

# Global profiler instance
profiler = Profiler()
