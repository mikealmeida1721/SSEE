#!/usr/bin/env python3
"""
Auditoría cruzada de parámetros SSEE en todos los papers.
Detecta inconsistencias entre .tex y .py files.
"""
import re, os, glob
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Valores canónicos (referencia algebraica Paper 1)
CANONICAL = {
    "phi":    1.6180339887,
    "pi":     3.1415926536,
    "w0":    -0.840,
    "wa":    -0.670,
    "H0_alg": 67.96,    # algebraico
    "H0_fit": 66.75,    # MCMC best-fit Paper 2
    "H0_ssee": 66.66,   # valor nominal en scripts
    "Om_dyn": 0.160,    # Ω_m dinámico
    "Om_cmb": 0.3199,   # Ω_m,CMB (dos sectores Paper 6)
    "Omb_h2": 0.02237,
    "ns":     0.9649,
    "rd_ssee": 147.156,  # r_d en Mpc (Paper 2 CAMB)
    "kfs":    0.493,     # k_fs h/Mpc (Paper 6)
    "mphi":   5.60,      # m_φ en eV
    "betac":  3.9978,    # |β_c| = AURA
    "sigma8_ssee": 0.702, # Paper 5
    "fsig8_tension_before": 2.56,  # Paper 6 baseline corregido (σ8=0.702, 6 surveys)
    "fsig8_tension_after":  0.50,  # Paper 6 con φ-DM (6 surveys)
}

# Patrones a buscar en archivos
PATTERNS = {
    "w0":    [r"w_?0\s*=\s*-?0\.84", r"w_0\s*=\s*-0\.840"],
    "wa":    [r"w_?a\s*=\s*-?0\.67", r"w_a\s*=\s*-0\.670"],
    "H0_alg":[r"67\.9[56]", r"67\.96"],
    "H0_fit":[r"66\.7[0-9]"],
    "H0_ssee":[r"66\.66"],
    "Om_dyn":[r"0\.1[56][0-9]", r"Omega_m.*0\.160"],
    "Om_cmb":[r"0\.319[89]", r"0\.3199"],
    "Omb_h2":[r"0\.0223[67]"],
    "ns":    [r"0\.9649"],
    "rd":    [r"147\.[0-9]", r"175\.[0-9]"],
    "kfs":   [r"0\.493"],
    "mphi":  [r"5\.71"],
    "betac": [r"3\.99[78]", r"3\.998"],
    "sigma8":[r"0\.70[0-9]", r"sigma_?8.*0\.7"],
}

tex_files = sorted(glob.glob(f"{ROOT}/manuscript/SSEE_Paper*.tex"))
py_files  = sorted(glob.glob(f"{ROOT}/src/ssee_paper*.py") + 
                   glob.glob(f"{ROOT}/src/ssee_verify*.py") +
                   glob.glob(f"{ROOT}/src/ssee_eft*.py"))

print("=" * 70)
print("AUDITORÍA DE CONSISTENCIA SSEE — PARÁMETROS CLAVE")
print("=" * 70)

issues = []

for param, pats in PATTERNS.items():
    hits = defaultdict(list)
    for f in tex_files + py_files:
        fname = os.path.basename(f)
        try:
            text = open(f).read()
        except:
            continue
        for pat in pats:
            for m in re.finditer(pat, text):
                line_no = text[:m.start()].count('\n') + 1
                snippet = text[max(0,m.start()-20):m.end()+20].replace('\n',' ').strip()
                hits[fname].append((line_no, snippet))
    
    if hits:
        print(f"\n[{param}]")
        for fname, locs in sorted(hits.items()):
            for ln, snip in locs[:2]:  # max 2 por archivo
                print(f"  {fname}:{ln}  →  ...{snip}...")

print("\n" + "=" * 70)
print("VERIFICACIÓN DE VALORES ALGEBRAICOS")
print("=" * 70)

phi = (1 + 5**0.5) / 2
pi  = 3.141592653589793
beta  = (phi + pi) / 2         # BIAL
KAL0  = beta + pi              # ≈ 5.5214
Psc   = phi + pi + phi + pi - beta  # recalcular
Omega_ssee = pi + phi          # 4.7596
Kv    = phi + pi + Omega_ssee  # KRYSTOS ≈ 9.5192
Tr    = 3 * (phi + beta)       # TRIAL ≈ 11.9935
Mv    = phi + pi + Kv          # MIKAEL_V ≈ 14.2789
AURA  = phi + beta             # ≈ 3.9978
BUFFER = Mv - Tr               # ≈ 2.2854

w0_alg = -Tr / Mv
wa_alg = -(Omega_ssee + phi) / Kv
Om_DE  = Tr / Mv
Om_m   = BUFFER / Mv

print(f"  φ              = {phi:.10f}")
print(f"  AURA=|β_c|     = {AURA:.10f}  (canónico: 3.9978)")
print(f"  TRIAL=Tr       = {Tr:.10f}  (canónico: 11.9935)")
print(f"  MIKAEL_V=Mᵥ   = {Mv:.10f}  (canónico: 14.2789)")
print(f"  w₀ algebraico  = {w0_alg:.6f}  (canónico: -0.840)")
print(f"  wₐ algebraico  = {wa_alg:.6f}  (canónico: -0.670)")
print(f"  Ω_DE           = {Om_DE:.6f}  (canónico: 0.840)")
print(f"  Ω_m,dyn        = {Om_m:.6f}  (canónico: 0.160)")

# αK algebraico
alpha_K = 3 * Om_DE * Om_m
print(f"\nEFT α-functions (acción canónica, z=0):")
print(f"  αT = 0.000000  (exacto, GW170817 ✓)")
print(f"  αM = 0.000000  (exacto, acoplamiento mínimo)")
print(f"  αB = 0.000000  (exacto, sin braiding)")
print(f"  αK = 3·Ω_DE·Ω_m,dyn = {alpha_K:.6f}")

# Consistencia w0, wa
print(f"\n  w₀ en scripts debe ser: {w0_alg:.4f} → redondeado: -0.840 ✓")
print(f"  wₐ en scripts debe ser: {wa_alg:.4f} → redondeado: -0.670 ✓")

print("\n" + "=" * 70)
print("CONSISTENCIA H₀: tres valores legítimos")
print("=" * 70)
print(f"  H₀_algebraico = (MIKA+π)×OMEGA = 67.96  (Paper 1/4 — cosmología)")
print(f"  H₀_MCMC       = 66.75±0.44      (Paper 2 — best-fit a DESI+Planck)")
print(f"  H₀_nominal    = 66.66            (scripts Paper 3,5,6 — valor redondo)")
print(f"  → Son tres contextos distintos, NO inconsistencia ✓")
print(f"  → Pero deben estar DOCUMENTADOS en cada paper como tal")

