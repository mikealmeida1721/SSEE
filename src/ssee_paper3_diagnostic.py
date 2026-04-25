"""
SSEE-V3.6 — Diagnóstico de desajuste CMB alto ell
Evalúa plik_lite TTTEEE en variantes de parámetros para aislar
cuál es el responsable del Δchi2≈145 vs ΛCDM.

Variantes (un parámetro cambiado a la vez hacia ΛCDM):
  A) SSEE completo          (baseline)
  B) w0/wa → ΛCDM           ¿es la energía oscura?
  C) H0 → ΛCDM              ¿es la distancia angular?
  D) omch2 → ΛCDM           ¿es la densidad de materia oscura?
  E) ns → ΛCDM              ¿es el índice espectral?
  F) tau → ΛCDM             ¿es la reionización?
  G) w0/wa + H0 → ΛCDM      ¿combinado?
  H) ΛCDM completo          (referencia)
"""

import numpy as np
import os

PACKAGES_PATH = "/home/mike/cobaya_packages"

# --- SSEE algebraico ---
phi   = (1 + 5**0.5) / 2
pi    = np.pi
Omega = pi + phi; beta = (pi + phi) / 2; KAL0 = beta + pi
P_sc  = Omega + phi; Kv = phi + pi + Omega
Tr    = 3 * (phi + beta); Mv = phi + pi + Kv
w0_s  = -Tr / Mv         # -0.8399
wa_s  = -P_sc / Kv       # -0.6700
MIRA  = (phi + (phi + pi) / 2) / 2   # 1.998924
Omm_c = (1 - Tr / Mv) * MIRA        # 0.3198
H0_s  = 66.66
ob_s  = 0.02237
oc_s  = Omm_c * (H0_s / 100)**2 - ob_s   # 0.11979
ns_s  = 1.0 - (1.0 / phi)**7        # 0.96556
As_s  = np.exp(3.044) * 1e-10
tau_s = 0.054

# --- ΛCDM best-fit ---
H0_l  = 67.36
ob_l  = 0.02237
oc_l  = 0.1200
ns_l  = 0.9649
As_l  = np.exp(3.044) * 1e-10
tau_l = 0.0544
w0_l  = -1.0
wa_l  = 0.0

VARIANTS = [
    # (label,  H0,   ombh2, omch2, w0,   wa,   ns,   tau,   As)
    ("A SSEE completo",          H0_s, ob_s, oc_s, w0_s, wa_s, ns_s, tau_s, As_s),
    ("B w0/wa → ΛCDM",           H0_s, ob_s, oc_s, w0_l, wa_l, ns_s, tau_s, As_s),
    ("C H0 → ΛCDM",              H0_l, ob_s, oc_s, w0_s, wa_s, ns_s, tau_s, As_s),
    ("D omch2 → ΛCDM",           H0_s, ob_s, oc_l, w0_s, wa_s, ns_s, tau_s, As_s),
    ("E ns → ΛCDM",              H0_s, ob_s, oc_s, w0_s, wa_s, ns_l, tau_s, As_s),
    ("F tau → ΛCDM",             H0_s, ob_s, oc_s, w0_s, wa_s, ns_s, tau_l, As_s),
    ("G w0/wa + H0 → ΛCDM",      H0_l, ob_s, oc_s, w0_l, wa_l, ns_s, tau_s, As_s),
    ("H ΛCDM completo",          H0_l, ob_l, oc_l, w0_l, wa_l, ns_l, tau_l, As_l),
]


def eval_plik(H0, ombh2, omch2, w0, wa, ns, tau, As):
    from cobaya.model import get_model
    info = {
        "packages_path": PACKAGES_PATH,
        "likelihood": {"planck_2018_highl_plik.TTTEEE_lite": None},
        "theory": {"camb": {"extra_args": {
            "dark_energy_model": "ppf",
            "WantTensors": False,
            "lens_potential_accuracy": 2,
        }}},
        "params": {
            "H0": H0, "ombh2": ombh2, "omch2": omch2,
            "As": float(As), "ns": ns, "tau": tau,
            "mnu": 0.06, "omk": 0.0,
            "w": w0, "wa": wa, "A_planck": 1.0,
        },
        "debug": False,
    }
    model = get_model(info)
    loglikes, _ = model.loglikes({})
    return -2 * float(loglikes[0])


def main():
    print("\nSSEE — Diagnóstico de desajuste plik_lite TTTEEE")
    print("="*65)
    print(f"{'Variante':<30}  {'chi2_eff':>10}  {'Δchi2 vs ΛCDM':>14}")
    print("-"*65)

    results = []
    for row in VARIANTS:
        label = row[0]
        H0, ob, oc, w0, wa, ns, tau, As = row[1:]
        print(f"  Evaluando {label}...", end="", flush=True)
        chi2 = eval_plik(H0, ob, oc, w0, wa, ns, tau, As)
        results.append((label, chi2))
        print(f"\r  {label:<30}  {chi2:>10.3f}", flush=True)

    chi2_lcdm = results[-1][1]
    print("="*65)
    for label, chi2 in results:
        delta = chi2 - chi2_lcdm
        bar = "▓" * int(abs(delta) / 5) if abs(delta) > 2 else "≈0"
        sign = "+" if delta >= 0 else ""
        print(f"  {label:<30}  {chi2:>10.3f}  {sign}{delta:>+8.1f}  {bar}")
    print("="*65)

    # Save
    out = os.path.join(os.path.dirname(__file__), "..",
                       "results", "planck_diagnostic_results.txt")
    with open(out, "w") as f:
        f.write("SSEE — Diagnóstico desajuste plik_lite TTTEEE\n")
        f.write("="*65 + "\n")
        f.write(f"{'Variante':<30}  {'chi2_eff':>10}  {'Δchi2 vs ΛCDM':>14}\n")
        f.write("-"*65 + "\n")
        for label, chi2 in results:
            delta = chi2 - chi2_lcdm
            f.write(f"  {label:<30}  {chi2:>10.3f}  {delta:>+14.1f}\n")
    print(f"\n  Guardado en {out}")


if __name__ == "__main__":
    main()
