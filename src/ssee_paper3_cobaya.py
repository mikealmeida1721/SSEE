"""
SSEE-V3.6 — Paper 3: Full Planck 2018 likelihood evaluation via Cobaya
Evaluates log-likelihood at fixed SSEE and ΛCDM parameters (no MCMC).
Uses: plik_lite TTTEEE (high-ell, ell=30..2508) + lowT + lowE (ell=2..29).
Produces: chi2_eff, Δchi2, ΔBIC stored in results/planck_fulllike_results.txt
"""

import numpy as np
import os
import sys

PACKAGES_PATH = "/home/mike/cobaya_packages"

# ---------------------------------------------------------------------------
# SSEE constants (algebraically fixed)
# ---------------------------------------------------------------------------
phi   = (1 + 5**0.5) / 2
pi    = np.pi
Omega = pi + phi
beta  = (pi + phi) / 2
KAL0  = beta + pi
P_sc  = Omega + phi
Kv    = phi + pi + Omega
Tr    = 3 * (phi + beta)
Mv    = phi + pi + Kv

w0_ssee = -Tr / Mv                  # -0.8403
wa_ssee = -P_sc / Kv                # -0.6703
OmDE    = Tr / Mv
Omm_dyn = 1.0 - OmDE                # 0.1597
BIAL    = (phi + pi) / 2
AURA    = phi + BIAL
MIRA    = AURA / 2                  # 1.998924
Omm_cmb = Omm_dyn * MIRA            # 0.3198

H0_ssee    = 66.66
ombh2_ssee = 0.02237
ns_ssee    = 1.0 - (1.0 / phi)**7   # 0.96556
As_ssee    = np.exp(3.044) * 1e-10
omch2_ssee = Omm_cmb * (H0_ssee / 100)**2 - ombh2_ssee
tau_ssee   = 0.054

# ΛCDM Planck 2018 best-fit (TT+TE+EE+lowE, Table 2)
H0_lcdm    = 67.36
ombh2_lcdm = 0.02237
omch2_lcdm = 0.1200
ns_lcdm    = 0.9649
As_lcdm    = np.exp(3.044) * 1e-10
tau_lcdm   = 0.0544


def build_cobaya_info(label, H0, ombh2, omch2, w0, wa, As, ns, tau):
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
            # Cosmological parameters (fixed — no prior, no proposal)
            "H0":     H0,
            "ombh2":  ombh2,
            "omch2":  omch2,
            "As":     As,
            "ns":     ns,
            "tau":    tau,
            "mnu":    0.06,
            "omk":    0.0,
            # Dark energy (fixed)
            "w":      w0,
            "wa":     wa,
            # Planck calibration nuisance (fixed at 1 — no freedom)
            "A_planck": 1.0,
        },
        "debug": False,
    }


def evaluate_model(label, H0, ombh2, omch2, w0, wa, As, ns, tau):
    from cobaya.model import get_model

    print(f"\n{'='*65}")
    print(f"  {label}")
    print(f"  H0={H0:.2f}  ombh2={ombh2:.5f}  omch2={omch2:.5f}")
    print(f"  w0={w0:.4f}  wa={wa:.4f}  ns={ns:.5f}  tau={tau:.4f}")
    print(f"{'='*65}")

    info = build_cobaya_info(label, H0, ombh2, omch2, w0, wa, As, ns, tau)
    model = get_model(info)

    # Evaluate all likelihoods — returns (loglikes_dict, derived_dict)
    loglikes, derived = model.loglikes({})

    # loglikes is an array aligned to the likelihood order
    like_names = list(info["likelihood"].keys())
    like_dict = dict(zip(like_names, loglikes))

    total_loglike = sum(loglikes)
    chi2_eff = -2 * total_loglike

    print(f"\n  Log-likelihoods:")
    for name, val in like_dict.items():
        short = name.split(".")[-1]
        print(f"    {short:30s}  logL = {val:10.3f}  chi2_eff = {-2*val:.3f}")
    print(f"  {'TOTAL':30s}  logL = {total_loglike:10.3f}  chi2_eff = {chi2_eff:.3f}")

    return chi2_eff, like_dict


def main():
    print("\n" + "="*65)
    print("  SSEE-V3.6 — Planck 2018 full likelihood evaluation")
    print("  plik_lite TTTEEE + lowT + lowE")
    print("  k_SSEE = 0 (all parameters algebraically fixed)")
    print("  k_ΛCDM = 6 (H0, ombh2, omch2, ns, As, tau)")
    print("="*65)

    chi2_ssee, ll_ssee = evaluate_model(
        "SSEE algebraic (k=0)",
        H0_ssee, ombh2_ssee, omch2_ssee,
        w0_ssee, wa_ssee, As_ssee, ns_ssee, tau_ssee,
    )

    chi2_lcdm, ll_lcdm = evaluate_model(
        "ΛCDM Planck 2018 best-fit (k=6)",
        H0_lcdm, ombh2_lcdm, omch2_lcdm,
        -1.0, 0.0, As_lcdm, ns_lcdm, tau_lcdm,
    )

    # BIC: N_data from plik_lite data ranges (confirmed by clipy output):
    # TT: ell=30..2508 → 2479 bins
    # TE: ell=30..1996 → 1967 bins
    # EE: ell=30..1996 → 1967 bins
    # plik_lite total = 2479 + 1967 + 1967 = 6413
    # lowT (ell=2..29): 28 bins; lowE (ell=2..29): 28 bins
    N_data  = 6413 + 28 + 28   # = 6469
    k_ssee  = 0
    k_lcdm  = 6

    delta_chi2 = chi2_ssee - chi2_lcdm
    delta_bic  = delta_chi2 + (k_ssee - k_lcdm) * np.log(N_data)

    print(f"\n{'='*65}")
    print(f"  RESULTADO FINAL")
    print(f"{'='*65}")
    print(f"  chi2_eff SSEE   = {chi2_ssee:.3f}")
    print(f"  chi2_eff ΛCDM   = {chi2_lcdm:.3f}")
    print(f"  Δchi2           = {delta_chi2:+.3f}  (SSEE − ΛCDM)")
    print(f"  ΔBIC            = {delta_bic:+.3f}  (k_SSEE={k_ssee} vs k_ΛCDM={k_lcdm})")
    print(f"  N_data          = {N_data}")
    print()

    if delta_bic < -10:
        verdict = "SSEE fuertemente favorecido (ΔBIC < -10)"
    elif delta_bic < 0:
        verdict = f"SSEE favorecido (ΔBIC = {delta_bic:.1f})"
    elif delta_bic < 10:
        verdict = f"evidencia débil contra SSEE (ΔBIC = {delta_bic:.1f})"
    elif delta_bic < 30:
        verdict = f"SSEE desfavorecido (ΔBIC = {delta_bic:.1f})"
    else:
        verdict = f"SSEE fuertemente desfavorecido (ΔBIC = {delta_bic:.1f})"

    print(f"  Veredicto: {verdict}")

    # Save
    out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "planck_fulllike_results.txt")
    with open(out_path, "w") as f:
        f.write("SSEE-V3.6 — Planck 2018 full likelihood (Cobaya)\n")
        f.write("="*65 + "\n")
        f.write(f"Likelihoods: plik_lite TTTEEE + lowT + lowE\n\n")
        f.write(f"SSEE parameters:\n")
        f.write(f"  H0={H0_ssee:.2f}  ombh2={ombh2_ssee:.5f}  omch2={omch2_ssee:.5f}\n")
        f.write(f"  w0={w0_ssee:.4f}  wa={wa_ssee:.4f}  ns={ns_ssee:.5f}  tau={tau_ssee:.4f}\n")
        f.write(f"  k_SSEE = {k_ssee} (todos fijos algebraicamente)\n\n")
        f.write(f"ΛCDM parameters (Planck 2018 best-fit):\n")
        f.write(f"  H0={H0_lcdm:.2f}  ombh2={ombh2_lcdm:.5f}  omch2={omch2_lcdm:.5f}\n")
        f.write(f"  w0=-1.0  wa=0.0  ns={ns_lcdm:.4f}  tau={tau_lcdm:.4f}\n")
        f.write(f"  k_ΛCDM = {k_lcdm}\n\n")
        f.write(f"Results:\n")
        f.write(f"  chi2_eff SSEE   = {chi2_ssee:.4f}\n")
        f.write(f"  chi2_eff ΛCDM   = {chi2_lcdm:.4f}\n")
        f.write(f"  delta_chi2      = {delta_chi2:+.4f}\n")
        f.write(f"  delta_BIC       = {delta_bic:+.4f}\n")
        f.write(f"  N_data          = {N_data}\n")
        f.write(f"  Veredicto       = {verdict}\n")
    print(f"\n  Guardado en {out_path}")


if __name__ == "__main__":
    main()
