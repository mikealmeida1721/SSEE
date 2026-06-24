#!/usr/bin/env python3
"""
SSEE-V3.6 — B1: Full CMB Posterior Analysis via Cobaya + plik TTTEEE (full)
=============================================================================
Runs full MCMC posteriors for SSEE and ΛCDM against the Planck 2018
plik TTTEEE (full, plik_rd12_HM_v22b) + lowl TT + lowl EE + lensing.native.

All ~20 plik nuisance parameters (dust, SZ, CIB, calibration) are included
automatically with their recommended Gaussian/uniform priors as defined by
the Planck Collaboration (Aghanim et al. 2020, A&A 641 A5).

SSEE model parameters (k=2, with --k2):
  logA     — log(10^10 As) amplitude   (genuinely sampled)
  tau      — reionization optical depth (constrained by lowl EE / SimAll)
  [default without --k2: H0 also floated as a k=3 validation chain that
   recovers the anchor, confirming the k=2 count is sound]
  Full plik TTTEEE + lowl TT + lowl EE + lensing.native via clipy (no clik needed).

SSEE fixed (algebraically — ω_m-directo reframe):
  H0    = 3(φ+π)²  = 67.962   (SH0ES–f_screen inversion, Paper 9; derived)
  w0    = -Tr/Mv   ≈ -0.8399
  wa    = -P_sc/Kv ≈ -0.6700
  ns    = 1 − (1/φ)^7 ≈ 0.96556
  ombh2 = (π−φ)/(3Ω²)        = 0.022423  (ω_b directo)
  omch2 = KAL₀·ombh2·ns      = 0.11951   (ω_c forward, NOT derived from H0)

ΛCDM free parameters (k=6):
  H0, ombh2, omch2, ns, logA, tau

Nuisance parameters (~20 plik + 1 lensing A_planck) are shared between
both models and are NOT counted in k for BIC purposes.

Usage:
  python ssee_paper3_b1_mcmc.py [--mode ssee|lcdm|both] [--fast]

  --fast: R−1 < 0.05  (validation/debug)
  default: R−1 < 0.02  (production, publishable)

Outputs:
  results/chains/ssee_cmb.*   — SSEE chain files
  results/chains/lcdm_cmb.*   — ΛCDM chain files
  results/figures/fig_b1_h0_posterior.pdf
  results/figures/fig_b1_corner_ssee.pdf
  results/figures/fig_b1_corner_lcdm.pdf
  results/tables/b1_mcmc_constraints.txt
"""

import numpy as np
import os
import sys
import argparse
import time
import json
import logging

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGES_PATH = os.environ.get("COBAYA_PACKAGES_PATH", os.path.expanduser("~/cobaya_packages"))
CHAINS_DIR = os.path.join(ROOT, "results", "chains")
FIGS_DIR   = os.path.join(ROOT, "results", "figures")
TABS_DIR   = os.path.join(ROOT, "results", "tables")
for d in [CHAINS_DIR, FIGS_DIR, TABS_DIR]:
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# SSEE algebraic constants
# ---------------------------------------------------------------------------
phi  = (1 + 5**0.5) / 2                # golden ratio
pi   = np.pi
Omega_ssee = pi + phi                  # 4.7596
beta       = (pi + phi) / 2            # 2.3798
P_sc       = Omega_ssee + phi          # 6.3776  (PYROS)
Kv         = phi + pi + Omega_ssee     # 9.5192  (KRYSTOS)
Tr         = 3 * (phi + beta)          # 11.9935 (TRIAL)
Mv         = phi + pi + Kv             # 14.2788 (MIKAEL_V)
BIAL       = (phi + pi) / 2
AURA       = phi + BIAL
MIRA       = AURA / 2                  # 1.998924

w0_ssee    = -Tr / Mv                  # −0.83989
wa_ssee    = -P_sc / Kv                # −0.66990
OmDE_ssee  = Tr / Mv                   # 0.83989
Omm_dyn    = 1.0 - OmDE_ssee          # 0.16011

# Reframe ω_m-DIRECTO (OP-8 cerrado, 2026-06-18): NO hay factor materia.
# ω_b y ω_c son densidades físicas FIJAS algebraicamente; Ω_m,CMB = ω_m/h² es DERIVADO.
ns_ssee    = 1.0 - (1.0 / phi) ** 7   # 0.96556 — algebraic spectral index
KAL0       = beta + pi                 # 5.5214 — Structural Viscosity
ombh2_ssee = (pi - phi) / (3.0 * Omega_ssee**2)   # 0.022423 — ω_b directo
omch2_ssee = KAL0 * ombh2_ssee * ns_ssee          # 0.11951 — ω_c forward (KAL₀·ω_b·n_s)
mnu_ssee   = 0.0690                    # Σm_ν canónico (R₂·0.960318 eV); ω_ν = Σm_ν/93.14
# Ω_m,CMB derivado: (ω_b+ω_c+ω_ν)/h² → 0.30889 @ H_alg=67.962 (era 0.31993 vía MIRA, retirado)
Omm_cmb    = (ombh2_ssee + omch2_ssee + mnu_ssee/93.14) / (67.962/100.0)**2  # ≈0.30889 (diagnóstico)

# ΛCDM Planck 2018 best-fit (TT+TE+EE+lowE, Table 2, arXiv:1807.06209)
H0_lcdm    = 67.36
ombh2_lcdm = 0.02237
omch2_lcdm = 0.1200
ns_lcdm    = 0.9649
logA_lcdm  = 3.044
tau_lcdm   = 0.0544

# ---------------------------------------------------------------------------
# Shared likelihoods + theory block
# ---------------------------------------------------------------------------
# Full plik (plik_rd12_HM_v22b_TTTEEE.clik) — the gold-standard Planck 2018
# high-ℓ likelihood used in all Planck parameter papers.  The ~20 nuisance
# parameters (dust TT/TE/EE, SZ, CIB, calibration) are defined automatically
# by Cobaya from the likelihood defaults; we do NOT override them here.
# clik no es necesario: Cobaya carga el plik FULL vía `clipy` (reimplementación
# pura-python de clik, verificada got=expected −1172.47). Usamos el likelihood
# completo off-diagonal — el mismo estándar de los papers de parámetros Planck.
# lowl EE constriñe τ con datos reales (NO un prior prestado).
LIKELIHOODS = {
    "planck_2018_highl_plik.TTTEEE": None,   # full plik TTTEEE (off-diagonal), ~6413 bins
    "planck_2018_lowl.TT": None,             # commander lowl TT, 28 bins
    "planck_2018_lowl.EE": None,             # SimAll lowl EE (constriñe τ), 28 bins
    "planck_2018_lensing.native": None,      # MV lensing, 9 bins
}

THEORY_CAMB = {
    "camb": {
        "extra_args": {
            "dark_energy_model": "ppf",
            "halofit_version": "mead",
            "WantTensors": False,
            "lens_potential_accuracy": 2,
        }
    }
}

# Effective data points for BIC (full plik suite):
#   plik full TTTEEE high-ℓ  ≈ 6413 multipole bins (TT+TE+EE)
#   lowl TT  (ℓ=2–29)          = 28 bins
#   lowl EE  (ℓ=2–29)          = 28 bins
#   lensing.native              = 9 bandpower bins
N_DATA = 6413 + 28 + 28 + 9  # = 6478

# ---------------------------------------------------------------------------
# Build Cobaya info dicts
# ---------------------------------------------------------------------------

def _ssee_info(output_prefix, Rminus1_stop=0.02, burn_in=300, h0_fixed=False):
    """Cobaya MCMC info for SSEE CMB analysis.

    h0_fixed=False  → k=3 validation chain: H0 floated on [64,72] (recovers anchor).
    h0_fixed=True   → k=2 model chain: H0 fixed at 3(φ+π)²=67.962, only {logA,τ} free.
    """
    H0_block = (
        67.962 if h0_fixed
        else {
            "prior": {"min": 64.0, "max": 72.0},
            "ref":  {"dist": "norm", "loc": 67.96, "scale": 0.4},   # ancla reframe H_alg=67.962
            "proposal": 0.3,
            "latex": r"H_0",
        }
    )
    return {
        "packages_path": PACKAGES_PATH,
        "likelihood": LIKELIHOODS,
        "theory": THEORY_CAMB,
        # τ se constriñe con lowl EE real (SimAll) — sin prior prestado.
        "params": {
            # --- Geometry anchor: derived (k=2) or floated for validation (k=3) ---
            "H0": H0_block,
            "logA": {
                "prior": {"min": 2.80, "max": 3.25},
                "ref":  {"dist": "norm", "loc": 3.044, "scale": 0.015},
                "proposal": 0.015,
                "drop": True,            # not passed to CAMB directly
                "latex": r"\log(10^{10}A_s)",
            },
            "tau": {
                "prior": {"min": 0.02, "max": 0.12},
                "ref":  {"dist": "norm", "loc": 0.054, "scale": 0.008},
                "proposal": 0.006,
                "latex": r"\tau_\mathrm{reio}",
            },
            # A_planck and all other plik nuisance parameters (~20 total)
            # are defined automatically by the full plik likelihood defaults.
            # --- Fixed algebraic ---
            "ombh2":  ombh2_ssee,        # 0.02242 — ω_b directo (fijo)
            "omch2":  omch2_ssee,        # 0.11951 — ω_c forward KAL₀·ω_b·n_s (fijo, NO derivado de H0)
            "ns":     ns_ssee,
            "w":      w0_ssee,
            "wa":     wa_ssee,
            "mnu":    mnu_ssee,          # 0.0690 — Σm_ν canónico
            "omk":    0.0,
            # --- Derived from sampled ---
            "As": {
                "value": lambda logA: 1e-10 * np.exp(logA),
                "derived": False,
                "latex": r"A_s",
            },
            # Ω_m,CMB = ω_m/h² es ahora un OUTPUT derivado (reframe), no un input.
            # --- Derived outputs (sigma8 only — Omega_m not a CAMB output param) ---
            "sigma8": {"derived": True, "latex": r"\sigma_8"},
        },
        "sampler": {
            "mcmc": {
                "Rminus1_stop": Rminus1_stop,
                "Rminus1_cl_stop": 0.2,
                "max_tries": 40 * 7,      # 40 × n_params
                "burn_in": burn_in,
                "oversample_power": 0.4,
            }
        },
        "output": output_prefix,
        "debug": False,
        "resume": True,
    }


def _lcdm_info(output_prefix, Rminus1_stop=0.02, burn_in=500):
    """Cobaya MCMC info for ΛCDM CMB analysis."""
    return {
        "packages_path": PACKAGES_PATH,
        "likelihood": LIKELIHOODS,
        "theory": THEORY_CAMB,
        "params": {
            "H0": {
                "prior": {"min": 60.0, "max": 80.0},
                "ref":  {"dist": "norm", "loc": H0_lcdm, "scale": 0.6},
                "proposal": 0.5,
                "latex": r"H_0",
            },
            "ombh2": {
                "prior": {"min": 0.019, "max": 0.026},
                "ref":  {"dist": "norm", "loc": ombh2_lcdm, "scale": 0.0002},
                "proposal": 0.0002,
                "latex": r"\Omega_b h^2",
            },
            "omch2": {
                "prior": {"min": 0.09, "max": 0.15},
                "ref":  {"dist": "norm", "loc": omch2_lcdm, "scale": 0.002},
                "proposal": 0.002,
                "latex": r"\Omega_c h^2",
            },
            "ns": {
                "prior": {"min": 0.88, "max": 1.06},
                "ref":  {"dist": "norm", "loc": ns_lcdm, "scale": 0.005},
                "proposal": 0.004,
                "latex": r"n_s",
            },
            "logA": {
                "prior": {"min": 2.80, "max": 3.25},
                "ref":  {"dist": "norm", "loc": logA_lcdm, "scale": 0.015},
                "proposal": 0.015,
                "drop": True,
                "latex": r"\log(10^{10}A_s)",
            },
            "tau": {
                "prior": {"min": 0.02, "max": 0.12},
                "ref":  {"dist": "norm", "loc": tau_lcdm, "scale": 0.008},
                "proposal": 0.006,
                "latex": r"\tau_\mathrm{reio}",
            },
            # A_planck and all other plik nuisance parameters defined by likelihood.
            # Derived:
            "As": {
                "value": lambda logA: 1e-10 * np.exp(logA),
                "derived": False,
            },
            "w":  -1.0,
            "wa":  0.0,
            "mnu": 0.06,
            "omk": 0.0,
            "sigma8": {"derived": True, "latex": r"\sigma_8"},
        },
        "sampler": {
            "mcmc": {
                "Rminus1_stop": Rminus1_stop,
                "Rminus1_cl_stop": 0.2,
                "max_tries": 40 * 10,
                "burn_in": burn_in,
                "oversample_power": 0.4,
            }
        },
        "output": output_prefix,
        "debug": False,
        "resume": True,
    }


# ---------------------------------------------------------------------------
# Run MCMC
# ---------------------------------------------------------------------------

def run_mcmc(label, info_fn, output_prefix, args):
    """Run or resume Cobaya MCMC chain."""
    from cobaya.run import run

    logging.getLogger("cobaya").setLevel(logging.WARNING)

    Rminus1 = 0.05 if args.fast else 0.02
    burn_in  = 100 if args.fast else (300 if label == "SSEE" else 500)

    if label == "SSEE":
        info = info_fn(output_prefix, Rminus1_stop=Rminus1, burn_in=burn_in,
                       h0_fixed=getattr(args, "k2", False))
    else:
        info = info_fn(output_prefix, Rminus1_stop=Rminus1, burn_in=burn_in)

    print(f"\n{'='*65}")
    print(f"  Running MCMC: {label}")
    print(f"  Output prefix: {output_prefix}")
    print(f"  Convergence target R−1 < {Rminus1}")
    print(f"{'='*65}")
    t0 = time.time()

    updated_info, sampler = run(info)

    elapsed = time.time() - t0
    print(f"\n  {label} MCMC finished in {elapsed/60:.1f} min")
    return updated_info, sampler


# ---------------------------------------------------------------------------
# Analysis: load chains and compute statistics
# ---------------------------------------------------------------------------

def analyse_chains(ssee_prefix, lcdm_prefix):
    """Load chains, compute constraints, ΔBIC, write table."""
    try:
        import getdist
        from getdist import loadMCSamples, plots
    except ImportError:
        print("  [WARNING] getdist not installed — skipping corner plots.")
        print("  Install with: pip install getdist")
        return None, None

    print("\n" + "="*65)
    print("  POST-PROCESSING CHAINS")
    print("="*65)

    results = {}

    for label, prefix in [("SSEE", ssee_prefix), ("LCDM", lcdm_prefix)]:
        chain_file = prefix + ".1.txt"
        if not os.path.exists(chain_file):
            print(f"  Chain not found: {chain_file}")
            continue

        samples = loadMCSamples(prefix, settings={"ignore_rows": 0.3})
        stats = samples.getMargeStats()

        params_of_interest = ["H0", "logA", "tau", "sigma8"]
        if label == "LCDM":
            params_of_interest += ["ombh2", "omch2", "ns"]

        print(f"\n  [{label}] Posterior constraints (68% CI):")
        row = {"label": label}
        for p in params_of_interest:
            try:
                m  = samples.mean(p)
                lo = samples.confidence(p, 0.16, upper=False)
                hi = samples.confidence(p, 0.84, upper=True)
                print(f"    {p:15s} = {m:.4f}  [{lo:.4f}, {hi:.4f}]")
                row[p] = (m, lo, hi)
            except Exception:
                pass

        # best-fit chi2
        try:
            bf = samples.getBestFit()
            chi2 = -2.0 * bf.logLike
            row["chi2_bf"] = chi2
            print(f"    chi2_bestfit      = {chi2:.3f}")
        except Exception:
            row["chi2_bf"] = None

        results[label] = row

    # Compute ΔBIC
    if "SSEE" in results and "LCDM" in results:
        chi2_ssee = results["SSEE"].get("chi2_bf")
        chi2_lcdm = results["LCDM"].get("chi2_bf")
        if chi2_ssee is not None and chi2_lcdm is not None:
            delta_chi2 = chi2_ssee - chi2_lcdm
            k_ssee = 3  # H0, logA, tau (A_planck shared)
            k_lcdm = 6  # H0, ombh2, omch2, ns, logA, tau
            delta_bic = delta_chi2 + (k_ssee - k_lcdm) * np.log(N_DATA)
            print(f"\n  Δχ²(SSEE−ΛCDM) = {delta_chi2:+.3f}")
            print(f"  ΔBIC (k_SSEE={k_ssee}, k_ΛCDM={k_lcdm}) = {delta_bic:+.3f}")
            # Conservative k=4
            delta_bic_cons = delta_chi2 + (4 - k_lcdm) * np.log(N_DATA)
            print(f"  ΔBIC conservative (k_SSEE=4) = {delta_bic_cons:+.3f}")
            results["delta_chi2"] = delta_chi2
            results["delta_bic"]  = delta_bic
            results["delta_bic_conservative"] = delta_bic_cons

    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(ssee_prefix, lcdm_prefix, results):
    """Generate publication-quality figures for Paper 3."""
    try:
        from getdist import loadMCSamples, plots
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError:
        print("  [WARNING] getdist/matplotlib unavailable — skipping figures.")
        return

    # Load samples
    ssee_samples = None
    lcdm_samples = None
    try:
        ssee_samples = loadMCSamples(ssee_prefix, settings={"ignore_rows": 0.3})
    except Exception as e:
        print(f"  [WARNING] Could not load SSEE chains: {e}")
    try:
        lcdm_samples = loadMCSamples(lcdm_prefix, settings={"ignore_rows": 0.3})
    except Exception as e:
        print(f"  [WARNING] Could not load ΛCDM chains: {e}")

    # ---- Figure 1: H0 posterior comparison ----
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = {"SSEE": "#2166ac", "LCDM": "#d6604d"}

    for label, samples, color in [("SSEE", ssee_samples, colors["SSEE"]),
                                   ("ΛCDM", lcdm_samples, colors["LCDM"])]:
        if samples is None:
            continue
        h0_chain = samples.getParams().H0
        weights   = samples.weights
        # Weighted histogram
        counts, bins = np.histogram(h0_chain, bins=60, weights=weights, density=True)
        centres = 0.5 * (bins[:-1] + bins[1:])
        ax.plot(centres, counts / counts.max(), lw=2, color=color, label=label)
        ax.fill_between(centres, 0, counts / counts.max(), alpha=0.15, color=color)

    # Planck 2018 prior (for reference)
    h0_arr = np.linspace(63, 72, 300)
    planck_prior = np.exp(-0.5 * ((h0_arr - 67.36) / 0.54) ** 2)
    ax.plot(h0_arr, planck_prior / planck_prior.max(), "k--", lw=1.5,
            alpha=0.6, label="Planck 2018 TT+TE+EE+lowE")

    # SSEE algebraic prediction
    ax.axvline(67.96, color="#1a9641", lw=1.5, ls=":", label=r"$H_0^{\rm alg}=67.96$")

    ax.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]", fontsize=13)
    ax.set_ylabel("Normalised posterior", fontsize=12)
    ax.set_xlim(63.5, 71.5)
    ax.set_ylim(0, 1.15)
    ax.legend(fontsize=10, framealpha=0.85)
    ax.set_title(r"$H_0$ posterior: SSEE vs $\Lambda$CDM (Planck 2018 plik full TTTEEE+lensing)", fontsize=11)
    plt.tight_layout()
    out = os.path.join(FIGS_DIR, "fig_b1_h0_posterior.pdf")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

    # ---- Figure 2: SSEE corner plot ----
    if ssee_samples is not None:
        g = plots.get_subplot_plotter(width_inch=6)
        g.settings.axes_fontsize = 9
        g.settings.lab_fontsize  = 11
        g.triangle_plot(
            [ssee_samples],
            params=["H0", "logA", "tau", "sigma8"],
            filled=True,
            colors=[colors["SSEE"]],
            legend_labels=["SSEE"],
        )
        out = os.path.join(FIGS_DIR, "fig_b1_corner_ssee.pdf")
        g.export(out)
        print(f"  Saved: {out}")

    # ---- Figure 3: ΛCDM corner plot ----
    if lcdm_samples is not None:
        g = plots.get_subplot_plotter(width_inch=8)
        g.settings.axes_fontsize = 8
        g.settings.lab_fontsize  = 10
        g.triangle_plot(
            [lcdm_samples],
            params=["H0", "ombh2", "omch2", "ns", "logA", "tau"],
            filled=True,
            colors=[colors["LCDM"]],
            legend_labels=[r"$\Lambda$CDM"],
        )
        out = os.path.join(FIGS_DIR, "fig_b1_corner_lcdm.pdf")
        g.export(out)
        print(f"  Saved: {out}")

    # ---- Figure 4: H0 vs sigma8 for SSEE (derived constraints) ----
    if ssee_samples is not None:
        try:
            g = plots.get_subplot_plotter(width_inch=5)
            g.plot_2d([ssee_samples], "H0", "sigma8",
                      filled=True, colors=[colors["SSEE"]])
            g.add_legend(["SSEE"], colored_text=True)
            ax2 = g.get_axes()
            ax2.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]")
            ax2.set_ylabel(r"$\sigma_8$")
            out = os.path.join(FIGS_DIR, "fig_b1_h0_sigma8_ssee.pdf")
            g.export(out)
            print(f"  Saved: {out}")
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Write constraints table
# ---------------------------------------------------------------------------

def write_table(results):
    """Write human-readable and machine-readable constraints."""
    lines = [
        "SSEE-V3.6 — B1 Full CMB Posterior Constraints",
        "Cobaya MCMC + Planck 2018 plik TTTEEE (full) + lowl TT + lowl EE + lensing",
        "=" * 70,
        "",
    ]

    for label in ["SSEE", "LCDM"]:
        if label not in results:
            continue
        row = results[label]
        lines.append(f"[{label}]")
        for p, v in row.items():
            if p in ("label",):
                continue
            if isinstance(v, tuple):
                m, lo, hi = v
                lines.append(f"  {p:15s}  {m:.5f}  [{lo:.5f}, {hi:.5f}]")
            elif isinstance(v, float):
                lines.append(f"  {'chi2_bestfit':15s}  {v:.3f}")
        lines.append("")

    if "delta_bic" in results:
        lines += [
            "BIC COMPARISON (from MCMC best-fits)",
            "-" * 40,
            f"  N_data          = {N_DATA}",
            f"  Δχ²(SSEE−ΛCDM) = {results['delta_chi2']:+.3f}",
            f"  ΔBIC (k=3 vs 6) = {results['delta_bic']:+.3f}",
            f"  ΔBIC (k=4 vs 6) = {results['delta_bic_conservative']:+.3f}",
            "",
            "Jeffreys scale: ΔBIC < −10 → decisive evidence for SSEE",
        ]

    out = os.path.join(TABS_DIR, "b1_mcmc_constraints.txt")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n  Table saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="SSEE B1 Full CMB MCMC")
    parser.add_argument("--mode", choices=["ssee", "lcdm", "both", "analyse"],
                        default="both",
                        help="Which model to run (default: both)")
    parser.add_argument("--fast", action="store_true",
                        help="Fast convergence (R−1<0.05) for debugging")
    parser.add_argument("--k2", action="store_true",
                        help="k=2 model run: fix H0 at the algebraic anchor 67.962 "
                             "(only {logA,τ} sampled). Default off = k=3 validation chain.")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip figure generation")
    args = parser.parse_args()

    ssee_prefix = os.path.join(CHAINS_DIR, "ssee_cmb_k2" if args.k2 else "ssee_cmb")
    lcdm_prefix = os.path.join(CHAINS_DIR, "lcdm_cmb")

    print("\n" + "="*65)
    print("  SSEE-V3.6 — B1 Full CMB Posterior Analysis")
    print(f"  Mode: {args.mode} | Fast: {args.fast}")
    print("="*65)
    print(f"\n  SSEE algebraic predictions:")
    print(f"    w₀    = {w0_ssee:.6f}  (−Tr/Mv)")
    print(f"    wₐ    = {wa_ssee:.6f}  (−P_sc/Kv)")
    print(f"    Ωm_cmb= {Omm_cmb:.6f}  (ω_m/h² derivado, ω_m-directo)")
    print(f"    ns    = {ns_ssee:.6f}  (1 − φ⁻⁷)")
    print(f"    N_data= {N_DATA}  (full plik TTTEEE + lowl TT+EE + lensing, via clipy)")

    # ---- Run MCMC chains ----
    if args.mode in ("ssee", "both"):
        run_mcmc("SSEE", _ssee_info, ssee_prefix, args)

    if args.mode in ("lcdm", "both"):
        run_mcmc("ΛCDM", _lcdm_info, lcdm_prefix, args)

    # ---- Analyse and plot ----
    if args.mode != "ssee" or os.path.exists(lcdm_prefix + ".1.txt"):
        results = analyse_chains(ssee_prefix, lcdm_prefix)
    else:
        results = analyse_chains(ssee_prefix, lcdm_prefix)

    if results and not args.no_plots:
        make_figures(ssee_prefix, lcdm_prefix, results)

    if results:
        write_table(results)

    print("\n" + "="*65)
    print("  B1 COMPLETE")
    print("="*65)


if __name__ == "__main__":
    main()
