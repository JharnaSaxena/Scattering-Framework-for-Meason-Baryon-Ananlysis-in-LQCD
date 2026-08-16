#!/usr/bin/env python3
"""
plot_fig10.py

Reproduces the structure of Fig. 10 (Morningstar / BaSc pi-Sigma paper):

  1. Lattice points:      (q/m_pi)^2  vs  (k/m_pi) cot(delta)   [from the 4 irreps]
  2. Luscher curves:      for each irrep, sweep E_cm near the lattice level and
                           plot B(E) (the finite-volume quantization curve) --
                           these are the steep near-vertical blue segments.
  3. ERE curve (+ band):  (k/m_pi) cot(delta) = p0 + p1*Delta(E_cm), evaluated
                           on a dense E_cm grid and converted to q^2 -- the
                           smooth dashed curve across the whole range.
  4. Virtual-state curve: below threshold, k cot(delta) = i k  =>  y = -sqrt(-x).
                           Purely kinematic, not fitted -- the black dashed line.
  5. Virtual-state star:  intersection of the ERE curve and the black curve.

Usage python3 QC2/plot_fig10.py

Run from the parent of QC2/ (same convention as run_generalized.py).
"""

import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from QC2.generalized_fit import GeneralizedFitModel
from QC2.morningstar_bmatrix import MorningstarBMatrix, IRREP_CONFIG
from QC2.unequal_kinematics import compute_kinematics
# Config reuse the same test energies as run_generalized.py
def create_test_energies():
    test_data = {
        "G1u(0)": {"E_cm": 0.637, "sigma": 0.001},
        "G1(1)":  {"E_cm": 0.645, "sigma": 0.002},
        "G(2)":   {"E_cm": 0.652, "sigma": 0.003},
        "G(3)":   {"E_cm": 0.660, "sigma": 0.004},
    }
    with open("test_energies.json", "w") as f:
        json.dump(test_data, f, indent=2)
    return "test_energies.json"


def q2_over_mpi2(E_cm, m_pi, m_sigma):
    """(q_cm/m_pi)^2 from E_cm, signed (negative below threshold)."""
    kin_like = compute_kinematics(E_cm, (0, 0, 0), m_pi, m_sigma, L=1.0)
    # q2_cm from compute_kinematics is already signed correctly
    return kin_like.q2_cm / m_pi**2


def ere_curve(E_cm, p0, p1, threshold):
    """(k/m_pi) cot(delta) = p0 + p1 * Delta(E_cm)"""
    Delta = (E_cm**2 - threshold**2) / threshold**2
    return p0 + p1 * Delta


def virtual_state_curve(x):
    """
    Kinematic curve for a virtual state below threshold
    Below threshold k = i*kappa, so k*cot(delta) = -i*k = kappa = sqrt(-k^2)
    i.e. y = +sqrt(-x)
    """
    x = np.asarray(x, dtype=float)
    y = np.full_like(x, np.nan)
    mask = x < 0
    y[mask] = np.sqrt(-x[mask])
    return y


def main():
    # Build model, load energies, fit (or reuse existing p0,p1)
    energy_file = create_test_energies()
    bmatrix = MorningstarBMatrix()

    model = GeneralizedFitModel(
        L=48.0, mpi=0.13957, msigma=0.49761,
        channel="Pion-Sigma",
        energy_file=energy_file,
        bmatrix=bmatrix,
    )

    results = model.fit(initial_guess=(0.047, 0.65))
    model.print_results()

    p0, p1 = results["p0"], results["p1"]
    m_pi, m_sigma = model.mpi, model.msigma
    threshold = model.threshold
    os.remove(energy_file)
    # Lattice points
    xs_pts, ys_pts = model.predict_points(p0, p1)
    # ERE curve + band 
    q2_grid = np.linspace(-0.25, 0.10, 400)               # (q/m_pi)^2 grid
    Ecm_grid = m_pi * np.sqrt(q2_grid + (threshold / m_pi) ** 2)
    ere_y = ere_curve(Ecm_grid, p0, p1, threshold)

    # crude error band from parameter covariance if available, else fixed width
    if getattr(model, "covariance", None) is not None:
        # propagate p0,p1 uncertainty (diagonal approx) if you have sigma_p0/p1
        band = 0.03 * np.ones_like(ere_y)
    else:
        band = 0.03 * np.ones_like(ere_y)

    # Per-irrep Luscher curves 
    from QC2.unequal_kinematics import compute_Ecm_from_k2, compute_E_lab_from_E_cm
    luscher_curves = {}
    for irrep in model.levels:
        d = IRREP_CONFIG[irrep]["d"]
        E_lattice = model.energies[irrep][0]
        # find the lattice point's own q^2 so we can center the sweep on it
        kin_lat = compute_kinematics(
            compute_E_lab_from_E_cm(E_lattice, d, model.L), d, m_pi, m_sigma, model.L
        )
        q2_lat = kin_lat.q2_cm / m_pi ** 2
        # sweep a SMALL, fixed window in q^2 NOT in E_cm (q^2 ~ E^2 near
        # threshold, so a tiny E window already covers a huge q^2 range here)
        q2_sweep = np.linspace(q2_lat - 0.03, q2_lat + 0.03, 150)

        xs_curve, ys_curve = [], []
        for q2 in q2_sweep:
            E_cm = compute_Ecm_from_k2(q2, m_pi, m_sigma)
            E_lab = compute_E_lab_from_E_cm(E_cm, d, model.L)
            kin = compute_kinematics(E_lab, d, m_pi, m_sigma, model.L)
            try:
                B = bmatrix.compute(
                    irrep=irrep, u=kin.u, gamma=kin.gamma,
                    d=tuple(int(v) for v in kin.d), E_cm=kin.E_cm,
                )
            except Exception:
                continue
            if not np.isfinite(B) or abs(B) > 5:
                xs_curve.append(np.nan)
                ys_curve.append(np.nan)
                continue
            xs_curve.append(kin.q2_cm / m_pi ** 2)
            ys_curve.append(B)

        luscher_curves[irrep] = (np.array(xs_curve), np.array(ys_curve))
    black_y = virtual_state_curve(q2_grid)

    def diff(q2):
        Ecm = m_pi * np.sqrt(q2 + (threshold / m_pi) ** 2)
        return ere_curve(Ecm, p0, p1, threshold) - np.sqrt(-q2)

    star_x, star_y = None, None
    neg_mask = q2_grid < 0
    d_vals = diff(q2_grid[neg_mask])
    sign_changes = np.where(np.diff(np.sign(d_vals)) != 0)[0]
    if len(sign_changes):
        i = sign_changes[0]
        x_lo, x_hi = q2_grid[neg_mask][i], q2_grid[neg_mask][i + 1]
        try:
            star_x = brentq(diff, x_lo, x_hi)
            star_y = np.sqrt(-star_x)
        except Exception:
            pass
    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.fill_between(q2_grid, ere_y - band, ere_y + band,
                     color="tab:blue", alpha=0.2)
    ax.plot(q2_grid, ere_y, "b--", lw=1.5,
            label=r"Single channel: $(k/m_\pi)\cot\delta = p_0+p_1\Delta$")
    ax.plot(q2_grid, black_y, "k--", lw=1.2, label="Free particle / virtual state")
    markers = {"G1u(0)": "o", "G1(1)": "s", "G(2)": "D", "G(3)": "v"}
    for irrep, (xc, yc) in luscher_curves.items():
        ax.plot(xc, yc, "-", color="tab:blue", lw=1.2)
    for irrep, x, y in zip(model.levels, xs_pts, ys_pts):
        ax.plot(x, y, marker=markers.get(irrep, "o"), color="tab:blue",
                 mfc="white", mec="tab:blue", ms=8, mew=1.5,
                 label=irrep, linestyle="None")
    if star_x is not None:
        ax.plot(star_x, star_y, marker="*", color="black", ms=16,
                 mfc="white", mec="black", label="Virtual state")
    ax.axhline(0, color="gray", lw=0.8, zorder=0)
    ax.axvline(0, color="gray", lw=0.8, zorder=0)
    ax.set_xlabel(r"$(k_{\pi\Sigma}/m_\pi)^2$")
    ax.set_ylabel(r"$\frac{k_{\pi\Sigma}}{m_\pi}\cot\delta_{\pi\Sigma}$")
    ax.set_xlim(-0.25, 0.10)
    ax.set_ylim(0.0, 0.65)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig10_reproduction.png")
    fig.savefig(out_path, dpi=200)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
