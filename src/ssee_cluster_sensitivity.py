"""
SSEE-V3.6 — Cluster sensitivity analysis
Varies dM_obs by ±50% and checks how chi2_r changes.
Tests robustness of the cluster fit to observational uncertainty.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

phi  = (1 + 5**0.5) / 2
pi   = np.pi
beta = (pi + phi) / 2
KAL0 = beta + pi          # 5.5214
FNU  = 0.053              # neutrino fraction used in Paper 2

CLUSTERS = [
    {"name": "Coma",   "M_ig": 1.8, "dM_obs": 1.0, "M_obs": 9.8 },
    {"name": "A2029",  "M_ig": 2.2, "dM_obs": 1.2, "M_obs": 12.0},
    {"name": "A478",   "M_ig": 1.5, "dM_obs": 1.0, "M_obs": 8.0 },
    {"name": "Bullet", "M_ig": 1.2, "dM_obs": 1.0, "M_obs": 6.5 },
]

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "results", "figures")
os.makedirs(FIG_DIR, exist_ok=True)


def chi2_r(clusters, KAL=KAL0, fnu=FNU, scale=1.0):
    """chi2_r with error bars scaled by `scale`."""
    vals = [((c["M_ig"] * KAL * (1 + fnu) - c["M_obs"]) /
             (c["dM_obs"] * scale))**2
            for c in clusters]
    return sum(vals) / len(vals)


def sensitivity_scan():
    scales = np.linspace(0.5, 1.5, 201)   # 50% tighter to 50% looser
    chi2_vals = [chi2_r(CLUSTERS, scale=s) for s in scales]

    chi2_nom = chi2_r(CLUSTERS, scale=1.0)
    chi2_50p = chi2_r(CLUSTERS, scale=1.5)   # 50% larger errors
    chi2_50m = chi2_r(CLUSTERS, scale=0.5)   # 50% smaller errors

    print("=" * 50)
    print("SSEE Cluster Sensitivity Analysis")
    print("=" * 50)
    print(f"\nNominal    (scale=1.0):  χ²_r = {chi2_nom:.3f}")
    print(f"Errors ×1.5 (+50%):      χ²_r = {chi2_50p:.3f}")
    print(f"Errors ×0.5 (−50%):      χ²_r = {chi2_50m:.3f}")

    # Individual cluster residuals
    print("\nIndividual cluster residuals (nominal):")
    print(f"  {'Cluster':<10} {'M_SSEE':>8} {'M_obs':>8} {'σ':>6} {'residual':>10}")
    for c in CLUSTERS:
        M_pred = c["M_ig"] * KAL0 * (1 + FNU)
        res    = (M_pred - c["M_obs"]) / c["dM_obs"]
        print(f"  {c['name']:<10} {M_pred:>8.2f} {c['M_obs']:>8.1f} "
              f"{c['dM_obs']:>6.1f} {res:>10.3f}σ")

    # Scale at which chi2_r crosses 1.0 and 2.0
    for threshold, label in [(1.0, "χ²_r = 1.0"), (2.0, "χ²_r = 2.0")]:
        idx = np.argmin(np.abs(np.array(chi2_vals) - threshold))
        print(f"\n  {label} crossed at scale = {scales[idx]:.2f}  "
              f"(errors {(scales[idx]-1)*100:+.0f}%)")

    return scales, chi2_vals, chi2_nom, chi2_50p, chi2_50m


def plot_sensitivity(scales, chi2_vals, chi2_nom):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(scales, chi2_vals, color="tab:blue", lw=2)
    ax.axhline(chi2_nom, color="gray", lw=1, ls="--",
               label=f"Nominal χ²_r = {chi2_nom:.3f}")
    ax.axhline(1.0, color="green", lw=1, ls=":", label="χ²_r = 1.0 (ideal)")
    ax.axhline(2.0, color="red",   lw=1, ls=":", label="χ²_r = 2.0 (marginal)")
    ax.axvline(1.0, color="gray",  lw=0.8, ls="-", alpha=0.5)
    ax.fill_between(scales, chi2_vals, 2.0,
                    where=np.array(chi2_vals) < 2.0,
                    alpha=0.12, color="tab:blue", label="χ²_r < 2 (acceptable)")
    ax.set_xlabel("Error bar scale factor", fontsize=13)
    ax.set_ylabel(r"$\chi^2_r$ (clusters)", fontsize=13)
    ax.set_title("SSEE Cluster Fit: Sensitivity to Observational Uncertainties",
                 fontsize=12)
    ax.legend(fontsize=10)
    ax.set_xlim(0.5, 1.5)
    ax.grid(True, alpha=0.3)
    out = os.path.join(FIG_DIR, "fig_cluster_sensitivity.pdf")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"\n  Figura guardada: {out}")


def plot_individual_residuals():
    """Bar chart of per-cluster residuals."""
    fig, ax = plt.subplots(figsize=(7, 4))
    names  = [c["name"] for c in CLUSTERS]
    M_pred = [c["M_ig"] * KAL0 * (1 + FNU) for c in CLUSTERS]
    M_obs  = [c["M_obs"] for c in CLUSTERS]
    dM     = [c["dM_obs"] for c in CLUSTERS]
    res    = [(p - o) / d for p, o, d in zip(M_pred, M_obs, dM)]
    colors = ["tab:green" if abs(r) < 1 else
              "tab:orange" if abs(r) < 2 else "tab:red"
              for r in res]
    x = np.arange(len(names))
    ax.bar(x, res, color=colors, alpha=0.8, edgecolor="k", linewidth=0.7)
    ax.axhline(0,  color="k",    lw=0.8)
    ax.axhline(+1, color="gray", lw=0.7, ls="--")
    ax.axhline(-1, color="gray", lw=0.7, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=12)
    ax.set_ylabel(r"Residual $(M_{\rm SSEE} - M_{\rm obs})/\sigma$", fontsize=11)
    ax.set_title("SSEE Cluster Residuals (nominal uncertainties)", fontsize=12)
    ax.grid(True, axis="y", alpha=0.3)
    out = os.path.join(FIG_DIR, "fig_cluster_residuals.pdf")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Figura guardada: {out}")


if __name__ == "__main__":
    scales, chi2_vals, chi2_nom, chi2_50p, chi2_50m = sensitivity_scan()
    plot_sensitivity(scales, chi2_vals, chi2_nom)
    plot_individual_residuals()
    print("\n✓ Sensibilidad de cúmulos completada.")
