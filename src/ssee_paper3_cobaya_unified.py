"""
SSEE-V3.6 — CMB Unified Likelihood Evaluation via Cobaya
Scans H0 to find the minimum chi2_eff against the full Planck 2018 plik_lite TTTEEE + lowT + lowE.
Calculates Delta BIC correctly with k=2 for SSEE vs k=6 for LambdaCDM.
"""

import numpy as np
import os
import time
from scipy.optimize import minimize_scalar

PACKAGES_PATH = os.environ.get("COBAYA_PACKAGES_PATH", os.path.expanduser("~/cobaya_packages"))

# ---------------------------------------------------------------------------
# SSEE constants (algebraically fixed)
# ---------------------------------------------------------------------------
from ssee_core import (
    PHI as phi, PI as pi, OMEGA as Omega, BETA as beta, KAL0,
    P_SC as P_sc, K_V as Kv, T_R as Tr, M_V as Mv,
    W0 as w0_ssee, WA as wa_ssee, OMEGA_DE as OmDE,
    OMEGA_M_DYN as Omm_dyn, AURA, MIRA, OMEGA_M_CMB_MIRA as Omm_cmb,
    N_S as ns_ssee,
)

ombh2_ssee = 0.02237                # Planck 2018 prior (no algebraico)
# ns_ssee = 1 - phi^-7 = 0.96556 — importado de ssee_core
As_ssee    = np.exp(3.044) * 1e-10
tau_ssee   = 0.054

# ΛCDM Planck 2018 best-fit (TT+TE+EE+lowE, Table 2)
H0_lcdm    = 67.36
ombh2_lcdm = 0.02237
omch2_lcdm = 0.1200
ns_lcdm    = 0.9649
As_lcdm    = np.exp(3.044) * 1e-10
tau_lcdm   = 0.0544

def build_cobaya_info(H0, ombh2, omch2, w0, wa, As, ns, tau):
    """Build a Cobaya model info dict for fixed-parameter likelihood evaluation."""
    return {
        "packages_path": PACKAGES_PATH,
        "likelihood": {
            "planck_2018_highl_plik.TTTEEE_lite": None,
            "planck_2018_lowl.TT": None,
            "planck_2018_lowl.EE": None,
        },
        "theory": {
            "camb": {
                "extra_args": {
                    "dark_energy_model": "ppf",
                    "halofit_version": "mead",
                    "WantTensors": False,
                    "lens_potential_accuracy": 2,
                },
            }
        },
        "params": {
            "H0":     H0,
            "ombh2":  ombh2,
            "omch2":  omch2,
            "As":     As,
            "ns":     ns,
            "tau":    tau,
            "mnu":    0.06,
            "omk":    0.0,
            "w":      w0,
            "wa":     wa,
            "A_planck": 1.0,
        },
        "debug": False,
    }

def evaluate_model(H0, ombh2, omch2, w0, wa, As, ns, tau, quiet=False):
    from cobaya.model import get_model
    import logging
    
    # Suppress cobaya output during scan
    logging.getLogger('cobaya').setLevel(logging.ERROR)
    
    info = build_cobaya_info(H0, ombh2, omch2, w0, wa, As, ns, tau)
    model = get_model(info)

    loglikes, derived = model.loglikes({})
    total_loglike = sum(loglikes)
    chi2_eff = -2 * total_loglike

    if not quiet:
        like_names = list(info["likelihood"].keys())
        like_dict = dict(zip(like_names, loglikes))
        print(f"\n  Log-likelihoods at minimum:")
        for name, val in like_dict.items():
            short = name.split(".")[-1]
            print(f"    {short:30s}  logL = {val:10.3f}  chi2_eff = {-2*val:.3f}")
        print(f"  {'TOTAL':30s}  logL = {total_loglike:10.3f}  chi2_eff = {chi2_eff:.3f}")

    return chi2_eff

def get_ssee_chi2(H0):
    omch2 = Omm_cmb * (H0 / 100)**2 - ombh2_ssee
    print(f"Evaluating SSEE at H0 = {H0:.3f} (omch2 = {omch2:.5f})... ", end='', flush=True)
    t0 = time.time()
    chi2 = evaluate_model(H0, ombh2_ssee, omch2, w0_ssee, wa_ssee, As_ssee, ns_ssee, tau_ssee, quiet=True)
    print(f"chi2_eff = {chi2:.3f} ({time.time()-t0:.1f}s)")
    return chi2

def main():
    print("\n" + "="*65)
    print("  SSEE-V3.6 — CMB Unified Likelihood Optimization")
    print("  Minimizing chi2_eff over H0 via Cobaya (plik_lite TTTEEE + lowl)")
    print("="*65)

    print("\n1. Evaluating ΛCDM baseline...")
    chi2_lcdm = evaluate_model(H0_lcdm, ombh2_lcdm, omch2_lcdm, -1.0, 0.0, As_lcdm, ns_lcdm, tau_lcdm, quiet=True)
    print(f"   ΛCDM chi2_eff = {chi2_lcdm:.3f}")

    print("\n2. Scanning H0 for SSEE...")
    # SciPy minimize_scalar with bounded method
    res = minimize_scalar(get_ssee_chi2, bounds=(66.5, 67.5), method='bounded', options={'xatol': 0.05})
    
    H0_opt = res.x
    chi2_ssee = res.fun
    
    print("\n3. Re-evaluating SSEE at optimal H0 for detailed breakdown...")
    omch2_opt = Omm_cmb * (H0_opt / 100)**2 - ombh2_ssee
    evaluate_model(H0_opt, ombh2_ssee, omch2_opt, w0_ssee, wa_ssee, As_ssee, ns_ssee, tau_ssee, quiet=False)

    N_data  = 6413 + 28 + 28   # 6469
    k_ssee  = 2                # H0, ombh2
    k_lcdm  = 6                # H0, ombh2, omch2, ns, As, tau

    delta_chi2 = chi2_ssee - chi2_lcdm
    delta_bic  = delta_chi2 + (k_ssee - k_lcdm) * np.log(N_data)

    print(f"\n{'='*65}")
    print(f"  RESULTADO FINAL UNIFICADO")
    print(f"{'='*65}")
    print(f"  H0 óptimo SSEE  = {H0_opt:.3f}")
    print(f"  chi2_eff SSEE   = {chi2_ssee:.3f}")
    print(f"  chi2_eff ΛCDM   = {chi2_lcdm:.3f}")
    print(f"  Δchi2           = {delta_chi2:+.3f}  (SSEE − ΛCDM)")
    print(f"  ΔBIC            = {delta_bic:+.3f}  (k_SSEE={k_ssee} vs k_ΛCDM={k_lcdm})")
    print()

    out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "planck_cobaya_unified.txt")
    with open(out_path, "w") as f:
        f.write("SSEE-V3.6 — Planck 2018 Unified Evaluation (Cobaya)\n")
        f.write("="*65 + "\n")
        f.write(f"Optimal SSEE H0 = {H0_opt:.3f}\n")
        f.write(f"chi2_eff SSEE   = {chi2_ssee:.4f}\n")
        f.write(f"chi2_eff ΛCDM   = {chi2_lcdm:.4f}\n")
        f.write(f"delta_chi2      = {delta_chi2:+.4f}\n")
        f.write(f"delta_BIC       = {delta_bic:+.4f}\n")
        f.write(f"k_SSEE = {k_ssee}, k_ΛCDM = {k_lcdm}\n")
    
    print(f"  Guardado en {out_path}")

if __name__ == "__main__":
    main()
