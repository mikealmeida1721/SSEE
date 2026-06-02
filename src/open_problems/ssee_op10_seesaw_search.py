"""
OP-10 Fase 2d — Búsqueda see-saw: el "punto medio" entre escalas SSEE
========================================================================

INTUICIÓN (Mike, 2026-05-24): "para unificar dos extremos, ¿no se necesitaría
un punto medio?" — la falla de las 4 familias single-field es precisamente que
intentan saltar directo entre dos escalas separadas por 10³-10⁴, sin puente.

MECANISMO SEE-SAW: m_intermediate² = m_light · m_heavy
Si los extremos son escalas SSEE definidas, el "punto medio" es derivable.

ESCALAS PRESENTES EN SSEE:
  M_pl     = 2.435×10²⁷ eV       (Planck reducida)
  M_UV     = 8.81 meV            (cutoff Paper 10)
  V₀^(1/4) = 2.87 meV            (escala DE)
  H₀       = 1.45×10⁻³³ eV       (Hubble hoy)
  Σm_ν     = 0.0824 eV           (neutrino sum)
  m_φ_target = 5.60 eV           (φ-DM, OP-9)

PREGUNTA: ¿alguna combinación see-saw de las escalas conocidas da 5.60 eV?
"""
import numpy as np
from itertools import combinations

PHI    = (1 + np.sqrt(5)) / 2
PI     = np.pi
KAL    = (PHI + 3*PI) / 2

# Escalas físicas SSEE (todas en eV)
scales = {
    "M_pl":           2.435e27,
    "M_UV":           8.81e-3,
    "V0_quarter":     (0.840 * 8.10e-11)**0.25,   # ≈ 2.87 meV
    "H0":             1.45e-33,
    "Sigma_mnu":      0.0824,
    "m_e":            5.11e5,                      # electron mass (BSM ref)
    "T_CMB":          2.35e-4,                     # 2.725 K en eV
    "Lambda_QCD":     0.220,                       # 220 MeV
    "M_W":            8.04e10,                     # 80.4 GeV
}

TARGET = 5.60  # eV

print("=" * 76)
print("OP-10 FASE 2d — Búsqueda see-saw para m_φ = 5.60 eV")
print("=" * 76)
print("\nEscalas SSEE consideradas (eV):")
for name, val in scales.items():
    print(f"  {name:<14} = {val:.4e}")

# ============================================================
# TEST 1 — see-saw clásico: m_φ² = m₁ × m₂
# ============================================================
print("\n" + "=" * 76)
print("TEST 1 — m_φ² = m_light · m_heavy (see-saw clásico)")
print("=" * 76)
print(f"\nTarget: m_φ² = {TARGET**2:.4f} eV²")
print(f"\n{'par':<35} {'producto [eV²]':<18} {'m_φ predicho':<14} {'ratio'}")
print("-" * 80)

results_t1 = []
for (n1, v1), (n2, v2) in combinations(scales.items(), 2):
    product = v1 * v2
    m_pred = np.sqrt(abs(product))
    ratio = m_pred / TARGET
    if 0.1 < ratio < 10:
        results_t1.append((f"{n1} · {n2}", product, m_pred, ratio))
        print(f"  {n1+' · '+n2:<33} {product:<18.3e} {m_pred:<14.3e} {ratio:<.4f}")

# ============================================================
# TEST 2 — see-saw modificado: m_φ² = m₁ · m₂ × factor algebraico
# ============================================================
print("\n" + "=" * 76)
print("TEST 2 — m_φ² = m_light · m_heavy × C (factor algebraico)")
print("=" * 76)

# Factores algebraicos SSEE
factors = {
    "1":          1,
    "φ":          PHI,
    "π":          PI,
    "KAL":        KAL,
    "φ·π":        PHI*PI,
    "Ω":          PHI + PI,
    "1/φ":        1/PHI,
    "1/π":        1/PI,
    "1/KAL":      1/KAL,
    "4π":         4*PI,
    "(2π)²":      (2*PI)**2,
    "(φ+π)²":     (PHI+PI)**2,
    "φ⁴":         PHI**4,
}

print(f"\n{'par + factor':<45} {'m_φ predicho':<14} {'ratio'}")
print("-" * 80)

best_candidates = []
for (n1, v1), (n2, v2) in combinations(scales.items(), 2):
    product = v1 * v2
    if product <= 0: continue
    base = np.sqrt(product)
    for fname, fval in factors.items():
        m_pred = base * fval
        ratio = m_pred / TARGET
        if 0.95 < ratio < 1.05:  # 5% tolerance
            best_candidates.append((f"√({n1}·{n2})·{fname}", m_pred, ratio))
            print(f"  √({n1}·{n2})·{fname:<25} {m_pred:<14.4e} {ratio:.4f}")
        # También probar dividir por factor
        m_pred2 = base / fval
        ratio2 = m_pred2 / TARGET
        if 0.95 < ratio2 < 1.05 and fname != "1":
            best_candidates.append((f"√({n1}·{n2})/{fname}", m_pred2, ratio2))
            print(f"  √({n1}·{n2})/{fname:<25} {m_pred2:<14.4e} {ratio2:.4f}")

# ============================================================
# TEST 3 — see-saw triple: m_φ³ = m₁ · m₂ · m₃
# ============================================================
print("\n" + "=" * 76)
print("TEST 3 — see-saw triple m_φ³ = m₁·m₂·m₃ (raíz cúbica)")
print("=" * 76)

print(f"\n{'triple':<55} {'m_φ predicho':<14} {'ratio'}")
print("-" * 80)
for combo in combinations(scales.items(), 3):
    names = [c[0] for c in combo]
    vals  = [c[1] for c in combo]
    product = np.prod(vals)
    m_pred = product**(1/3)
    ratio = m_pred / TARGET
    if 0.5 < ratio < 2.0:
        print(f"  ∛({'·'.join(names)}):<25 {m_pred:<14.4e} {ratio:.4f}")

# ============================================================
# TEST 4 — Análisis dimensional Paper 6: m_φ = Σm_ν · H₀ con H₀ "no físico"
# ============================================================
print("\n" + "=" * 76)
print("TEST 4 — ¿Qué unidad de referencia rescata m_φ = Σm_ν · H₀?")
print("=" * 76)

# Si la fórmula de P6 m_φ = Σm_ν × H₀ funciona numéricamente como
#     5.60 eV = 0.0824 eV × 67.96
# entonces el "67.96" debería ser una razón adimensional. ¿Cuál?
#
# H₀_alg = 67.96 km/s/Mpc. Para que sea adimensional, dividir por una escala.

# H₀ en eV: 1.45e-33. Razón H₀_eV / X = 67.96 → X = H₀_eV/67.96 = 2.13e-35 eV
H0_eV = 1.45e-33
unit_ref = H0_eV / 67.96
print(f"\nSi m_φ = Σm_ν × (H₀/X) con X tal que (H₀/X) = 67.96:")
print(f"  X = H₀/67.96     = {unit_ref:.3e} eV")
print(f"  X·M_pl           = {unit_ref * 2.435e27:.3e}")
print(f"  X/H₀             = {unit_ref/H0_eV:.6f}  (= 1/67.96 = {1/67.96:.6f})")

# Equivalentemente: m_φ = Σm_ν · h adimensional (h = 0.6796)... pero 0.6796 ≠ 67.96.
# El factor 100 viene de la conversión km/s/Mpc → unidades naturales.
# Más probable: m_φ = Σm_ν × H₀_alg [km/s/Mpc] / (1 km/s/Mpc) en sistema CGS-like

print(f"""
Diagnóstico: la fórmula m_φ = Σm_ν · H₀^alg (Paper 6) trata H₀ como número
67.96, no como energía 1.45e-33 eV. Esto es físicamente inconsistente
(viola análisis dimensional).

La fórmula tiene UN sentido si se reescribe como:
  m_φ = Σm_ν · [H₀^alg / (km·s⁻¹·Mpc⁻¹)]
      = 0.0824 eV · 67.96  [adimensional]
      = 5.60 eV

Pero entonces el "67.96" es un número específico de la convención de unidades
cosmológicas, NO una constante algebraica universal. Esta es la raíz del
problema dimensional de OP-9: m_φ depende de la convención.

INSIGHT: si m_φ no es realmente "5.60 eV físicos" sino "Σm_ν × H₀/H_unit",
entonces el target real para Fase 2 puede ser diferente.

Tres posibles reinterpretaciones físicas:
  (a) m_φ_real = 5.60 eV — el "número adimensional" 67.96 viene de algo
      profundo (e.g., M_v = 3(φ+π) ≈ 14.28; M_v² ≈ 204; ratio M_v²/π ≈ 65)
  (b) m_φ_real = 1.20e-34 eV — interpretación literal SI unidades
  (c) m_φ_real = Σm_ν (unidades eV directas, sin H₀) — entonces m_φ ≈ 0.08 eV
""")

print("=" * 76)
print("RESUMEN")
print("=" * 76)

if len(best_candidates) > 0:
    print(f"\n  ✓ {len(best_candidates)} candidatos see-saw dentro de ±5% del target:")
    for name, val, ratio in best_candidates[:8]:
        print(f"    {name} = {val:.4e}  (ratio = {ratio:.4f})")
else:
    print("\n  ✗ Ningún see-saw simple (de 2 escalas × factor algebraico) alcanza")
    print("    m_φ = 5.60 eV dentro de 5%.")
    print()
    print("  Esto sugiere que m_φ NO sale de un see-saw entre las escalas")
    print("  fenomenológicas obvias. Posibilidades:")
    print("    1) Se necesita una escala adicional (¿escala EW? ¿escala QCD?)")
    print("    2) m_φ = 5.60 eV es interpretación incorrecta dimensionalmente")
    print("    3) La estructura no es see-saw clásica sino algo más exótico")
