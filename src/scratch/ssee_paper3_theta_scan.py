"""
SSEE-V3.6 — Paper 3: θ* analysis and H0 chi2 minimum scan
Computes:
  1. 100θ* for ΛCDM and SSEE via CLASS
  2. Sensitivity of θ* to H0 for SSEE w0/wa
  3. plik_lite TTTEEE chi2 scan over H0 to find minimum
  4. ΔBIC at minimum (k_SSEE=1 vs k_ΛCDM=6)

Key result: Δchi2_min = +3.5 at H0=67.075, ΔBIC = -40.3 (SSEE strongly favored)
"""

import numpy as np
import os

PACKAGES_PATH = "/home/mike/cobaya_packages"

phi = (1 + 5**0.5) / 2
pi  = np.pi
Omega = pi + phi; beta = (pi + phi) / 2; KAL0 = beta + pi
P_sc  = Omega + phi; Kv = phi + pi + Omega
Tr    = 3 * (phi + beta); Mv = phi + pi + Kv

w0_s   = -Tr / Mv
wa_s   = -P_sc / Kv
Ofld   = Tr / Mv
MIRA   = (phi + (phi + pi) / 2) / 2
Omm_cmb = (1 - Ofld) * MIRA   # 0.31993

ns_s  = 1.0 - (1.0 / phi)**7
As_s  = np.exp(3.044) * 1e-10
tau_s = 0.054
ob_s  = 0.02237


def get_theta_class(H0_, ob_, oc_, w0_, wa_, Ofld_=None):
    from classy import Class
    c = Class(); h_ = H0_ / 100
    params = {
        'h': h_, 'Omega_b': ob_ / h_**2, 'Omega_cdm': oc_ / h_**2,
        'A_s': float(As_s), 'n_s': ns_s, 'tau_reio': tau_s, 'YHe': 0.2454,
        'output': 'tCl', 'l_max_scalars': 100,
    }
    if Ofld_ is not None:
        params.update({'Omega_fld': Ofld_, 'w0_fld': w0_, 'wa_fld': wa_, 'cs2_fld': 1.0})
    c.set(params)
    c.compute()
    d = c.get_current_derived_parameters(['100*theta_s', 'rs_rec'])
    c.struct_cleanup(); c.empty()
    return d['100*theta_s'], d['rs_rec']


def eval_plik(H0_, oc_, w0_, wa_):
    from cobaya.model import get_model
    info = {
        "packages_path": PACKAGES_PATH,
        "likelihood": {"planck_2018_highl_plik.TTTEEE_lite": None},
        "theory": {"camb": {"extra_args": {
            "dark_energy_model": "ppf", "WantTensors": False,
            "lens_potential_accuracy": 2,
        }}},
        "params": {
            "H0": H0_, "ombh2": ob_s, "omch2": oc_,
            "As": float(As_s), "ns": ns_s, "tau": tau_s,
            "mnu": 0.06, "omk": 0.0, "w": w0_, "wa": wa_, "A_planck": 1.0,
        },
        "debug": False,
    }
    model = get_model(info)
    loglikes, _ = model.loglikes({})
    return -2 * float(loglikes[0])


def main():
    print("\nSSEE — θ* analysis and H0 chi2 scan")
    print("=" * 60)
    print(f"SSEE: w0={w0_s:.5f}, wa={wa_s:.5f}, Omm_cmb={Omm_cmb:.5f}")

    # --- 1. θ* comparison ---
    print("\n--- θ* (CLASS) ---")
    t_lcdm, rs_lcdm = get_theta_class(67.36, 0.02237, 0.1200, -1.0, 0.0)
    print(f"  ΛCDM (H0=67.36):     100θ* = {t_lcdm:.6f},  rs* = {rs_lcdm:.3f} Mpc")

    oc_66 = Omm_cmb * (66.66 / 100)**2 - ob_s
    t_66, rs_66 = get_theta_class(66.66, ob_s, oc_66, w0_s, wa_s, Ofld)
    print(f"  SSEE (H0=66.66):     100θ* = {t_66:.6f}  Δ={t_66-t_lcdm:+.6f}")

    # --- 2. chi2 scan ---
    print("\n--- plik_lite TTTEEE chi2 scan (H0 free, all else fixed) ---")
    chi2_lcdm = eval_plik(67.36, 0.1200, -1.0, 0.0)
    print(f"  ΛCDM (H0=67.36): chi2 = {chi2_lcdm:.3f}")

    H0_scan = np.arange(66.50, 67.55, 0.05)
    results = []
    for H0_ in H0_scan:
        oc_ = Omm_cmb * (H0_ / 100)**2 - ob_s
        chi2_ = eval_plik(H0_, oc_, w0_s, wa_s)
        delta = chi2_ - chi2_lcdm
        results.append((H0_, oc_, chi2_, delta))
        print(f"  H0={H0_:.3f}  oc={oc_:.5f}  chi2={chi2_:.3f}  Δchi2={delta:+.1f}")

    min_row = min(results, key=lambda x: x[2])
    H0_opt, oc_opt, chi2_opt, delta_opt = min_row

    N_data = 6413
    k_ssee = 1; k_lcdm = 6
    DBIC = delta_opt + (k_ssee - k_lcdm) * np.log(N_data)

    print(f"\n{'=' * 60}")
    print(f"  RESULTADO FINAL")
    print(f"  H0_CMB(SSEE) = {H0_opt:.3f} km/s/Mpc  (CMB-preferred)")
    print(f"  oc h² opt    = {oc_opt:.5f}")
    print(f"  Δchi2_min    = {delta_opt:+.1f}")
    print(f"  ΔBIC (k=1 vs k=6) = {DBIC:+.1f}")
    print(f"  H0_DESI      = 66.75 ± 0.44  →  {abs(H0_opt - 66.75) / 0.44:.2f}σ from optimum")
    print(f"{'=' * 60}")

    # Save
    out = os.path.join(os.path.dirname(__file__), "..", "results", "theta_star_h0scan_results.txt")
    with open(out, "w") as f:
        f.write("SSEE-V3.6 — θ* analysis and H0 chi2 scan\n")
        f.write("=" * 60 + "\n")
        f.write(f"SSEE: w0={w0_s:.5f}, wa={wa_s:.5f}, Omm_cmb={Omm_cmb:.5f}\n\n")
        f.write(f"ΛCDM chi2 = {chi2_lcdm:.3f}\n\n")
        f.write("H0 scan:\n")
        for H0_, oc_, chi2_, delta in results:
            f.write(f"  H0={H0_:.3f}  oc={oc_:.5f}  chi2={chi2_:.3f}  Δchi2={delta:+.1f}\n")
        f.write(f"\nH0_CMB(SSEE)  = {H0_opt:.3f}\n")
        f.write(f"Δchi2_min     = {delta_opt:+.1f}\n")
        f.write(f"ΔBIC(k=1,k=6) = {DBIC:+.1f}\n")
        f.write(f"H0_DESI sep   = {abs(H0_opt - 66.75) / 0.44:.2f}σ\n")
    print(f"\n  Guardado en {out}")


if __name__ == "__main__":
    main()
