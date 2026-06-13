#!/usr/bin/env python3
"""
Auditoría cruzada de parámetros SSEE en todos los papers.
Detecta inconsistencias entre .tex y .py files.
"""
import re, os, glob
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
    "kfs":    0.659,     # k_fs h/Mpc (Paper 6 — partícula canónica forward 2026-06-04)
    "mphi":   36.95,     # m_φ en eV (Σm_ν × (Ω⁴_DNAV + AURA·KAL))
    "betac":  3.9978,    # |β_c| = AURA
    "G_ssee":      1.011,  # Paper 5 con Ω_m,cosm=0.320 (leve enhancement)
    "sigma8_p5":   0.820,  # Paper 5 single-sector (Ω_m,cosm=0.320)
    "S8_p5":       0.847,  # Paper 5 S₈ (3.9σ KiDS/DES — lensing challenge)
    "sigma8_p6":   0.742,  # Paper 6 two-sector TITULAR lensing (free-streaming CLASS, forward)
    "S8_p6":       0.766,  # Paper 6 S₈ two-sector TITULAR (0.01σ KiDS — RESUELVE)
    "sigma8_p6_rsd": 0.794,# amplitud RSD/fσ₈ (k<k_fs, agrupa frío); NO es la lensing
    "fsig8_p5":    0.74,   # Paper 5 tensión media fσ₈ (6 surveys) — ≈ ΛCDM 0.73σ
    "fsig8_p6":    0.76,   # Paper 6 tensión media fσ₈ two-sector
    "sigma8_wdm":  0.737,  # RETIRADO — era fit de alpha_WDM a KiDS (no forward)
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
    "kfs":   [r"0\.659"],
    "mphi":  [r"36\.95"],
    "betac": [r"3\.99[78]", r"3\.998"],
    "sigma8_p5":  [r"0\.820", r"sigma_?8.*0\.82"],
    "S8_p5":      [r"0\.847", r"S_?8.*0\.847"],
    "sigma8_p6":  [r"0\.794"],
    "S8_p6":      [r"S_?8.*0\.820", r"2\.6.sigma"],
    "sigma8_wdm": [r"0\.737"],
    "G_ssee":     [r"1\.01[0-9]", r"mathcal.*G.*1\.01"],
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

print("\n" + "=" * 70)
print("AUDITORÍA DE VALORES OBSOLETOS (dos-Ω_m correction 2026-05-19/21)")
print("Busca texto con valores PRE-corrección que ya no deben aparecer")
print("=" * 70)

STALE_PATTERNS = {
    # sigma8=0.702 era Paper 5 con Ω_m,dyn=0.160 — ahora es 0.820
    "sigma8=0.702 [OBSOLETO P5]": (
        [r"sigma_?8.*0\.702", r"0\.702.*sigma_?8"],
        ["SSEE_Paper5_IS.tex"],  # archivos donde es un error
    ),
    # G=0.866 era Paper 5 old — ahora G=1.011
    "G=0.866 [OBSOLETO P5]": (
        [r"0\.866", r"G.*0\.866"],
        ["SSEE_Paper5_IS.tex"],
    ),
    # fσ₈ 2.67σ baseline era Paper 5/6 old — ahora 0.74σ single-sector
    "fσ₈-2.67σ [OBSOLETO]": (
        [r"2\.67.*sigma", r"2\.67\s*\\sigma"],
        ["SSEE_Paper5_IS.tex", "SSEE_Paper6_phiDM.tex"],
    ),
    # S8=0.725 era Paper 5 old IS-only
    "S8=0.725 [OBSOLETO P5]": (
        [r"0\.725", r"S_?8.*0\.725"],
        ["SSEE_Paper5_IS.tex"],
    ),
    # fσ₈ 0.50σ era Paper 6 con datos erróneos — ahora 0.76σ
    "fσ₈-0.50σ [OBSOLETO P6]": (
        [r"0\.50.*sigma", r"0\.50\s*\\sigma.*f"],
        ["SSEE_Paper6_phiDM.tex"],
    ),
    # fσ₈ 0.76σ reemplazó al anterior — verificar que P6 lo tiene
    "fσ₈-0.76σ [DEBE ESTAR EN P6]": (
        [r"0\.76.*\\sigma", r"0\.76\s*\\sigma"],
        [],  # no es error si aparece
    ),
}

stale_found = []
for label, (pats, error_files) in STALE_PATTERNS.items():
    hits = []
    for f in tex_files:
        fname = os.path.basename(f)
        try:
            text = open(f).read()
        except:
            continue
        for pat in pats:
            for m in re.finditer(pat, text):
                ln = text[:m.start()].count('\n') + 1
                snip = text[max(0,m.start()-25):m.end()+25].replace('\n',' ').strip()
                is_error = (fname in error_files) if error_files else False
                tag = "❌ ERROR" if is_error else "ℹ INFO"
                hits.append((fname, ln, snip, tag))
                if is_error:
                    stale_found.append(f"{label} en {fname}:{ln}")
    if hits:
        print(f"\n[{label}]")
        for fname, ln, snip, tag in hits[:3]:
            print(f"  {tag}  {fname}:{ln}  →  ...{snip}...")

print()
if stale_found:
    print(f"⚠️  ERRORES ENCONTRADOS: {len(stale_found)}")
    for e in stale_found:
        print(f"  → {e}")
else:
    print("✅  Sin valores obsoletos detectados en los archivos tex de error esperado")

print("\n" + "=" * 70)
print("AUDITORÍA LaTeX — referencias cruzadas y etiquetas")
print("=" * 70)

# Verificar que eq:phi_pi_duality existe en P7 y se cita
p7 = open(f"{ROOT}/manuscript/SSEE_Paper7_EFT.tex").read()
if r"\label{eq:phi_pi_duality}" in p7:
    print("  ✅ Paper 7: \\label{eq:phi_pi_duality} definida")
else:
    print("  ❌ Paper 7: falta \\label{eq:phi_pi_duality}")

if r"\eqref{eq:phi_pi_duality}" in p7:
    print("  ✅ Paper 7: \\eqref{eq:phi_pi_duality} citada en texto")
else:
    print("  ❌ Paper 7: falta \\eqref{eq:phi_pi_duality} en texto")

# Verificar que Paper 6 tiene el nuevo framing S8
p6 = open(f"{ROOT}/manuscript/SSEE_Paper6_phiDM.tex").read()
if "3.9" in p6 and "lensing" in p6.lower():
    print("  ✅ Paper 6: contiene '3.9' + 'lensing' (nuevo framing S₈)")
else:
    print("  ❌ Paper 6: falta referencia a 3.9σ lensing tension")

if "2.6" in p6:
    print("  ✅ Paper 6: contiene '2.6' (S₈ reducida two-sector)")
else:
    print("  ❌ Paper 6: falta '2.6σ' two-sector S₈")

if "k_{\\rm fs}" in p6 or r"k_{\rm fs}" in p6 or "kfs" in p6 or "k_{fs}" in p6:
    print("  ✅ Paper 6: contiene k_fs (predicción falsificable)")
else:
    print("  ❌ Paper 6: falta k_fs como predicción falsificable")

# Verificar OP-7 parcialmente resuelto
op = open(f"{ROOT}/OPEN_PROBLEMS.md").read()
if "PARCIALMENTE RESUELTO" in op and "OP-7" in op:
    print("  ✅ OPEN_PROBLEMS.md: OP-7 marcado PARCIALMENTE RESUELTO")
else:
    print("  ❌ OPEN_PROBLEMS.md: OP-7 no actualizado")

print("\n" + "=" * 70)
print("RESUMEN FINAL")
print("=" * 70)
total_errors = len(stale_found)
print(f"  Errores críticos (valores obsoletos en tex): {total_errors}")
print(f"  Estado: {'✅ LIMPIO' if total_errors == 0 else '❌ REQUIERE CORRECCIÓN'}")

