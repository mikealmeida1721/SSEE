"""
OP-10 Fase 2d — Búsqueda see-saw: el "punto medio" entre escalas SSEE
========================================================================

⚠ SUPERADO (2026-06-06): la masa canónica φ-DM ya NO se busca por see-saw.
Tiene forward-derivation limpia con cero fiteo (Paper 6 partícula canónica):
    m_φ = Σm_ν · MULT,   MULT = Ω⁴ + AURA·KAL = 535.26,   Σm_ν = 0.06903 eV
    ⇒ m_φ = 36.95 eV  (Lagrangiano escalar libre — cierra OP-9)
Este script se conserva como exploración histórica; el target se re-apuntó al
valor canónico 36.95 eV para que sus conclusiones no queden contra el viejo
5.60 eV (obsoleto). NO es una búsqueda viva ni una predicción.

INTUICIÓN ORIGINAL (Mike, 2026-05-24): "para unificar dos extremos, ¿no se
necesitaría un punto medio?" — la falla de las 4 familias single-field es que
intentan saltar directo entre dos escalas separadas por 10³-10⁴, sin puente.

MECANISMO SEE-SAW: m_intermediate² = m_light · m_heavy
Si los extremos son escalas SSEE definidas, el "punto medio" es derivable.

ESCALAS PRESENTES EN SSEE:
  M_pl     = 2.435×10²⁷ eV       (Planck reducida)
  M_UV     = 9.62 meV            (cutoff Paper 10, ρ_crit×h_MIRA²)
  V₀^(1/4) = 2.35 meV            (escala DE, ρ_crit×h_MIRA²)
  H₀       = 1.45×10⁻³³ eV       (Hubble hoy)
  Σm_ν     = 0.06903 eV          (neutrino sum, R2·0.960318)
  m_φ      = 36.95 eV            (φ-DM canónico, forward-derivation)

PREGUNTA (histórica): ¿alguna combinación see-saw de las escalas conocidas
reproduce el m_φ canónico? (Hoy moot: ya hay derivación directa.)
"""
import numpy as np
from itertools import combinations

PHI    = (1 + np.sqrt(5)) / 2
PI     = np.pi
KAL    = (PHI + 3*PI) / 2
OMEGA  = PHI + PI
AURA   = (3*PHI + PI) / 2
TRIAL  = 3 * (PHI + (PHI + PI) / 2)

# ρ_crit canónico SSEE: 8.10e-11 eV⁴ es el coeficiente pelado (h=1); la densidad
# física se escala por h_MIRA²=0.67037² (Paper 3 Cobaya).
H_MIRA       = 67.037
h_MIRA       = H_MIRA / 100.0
RHO_CRIT_eV4 = 8.10e-11 * h_MIRA**2
M_UV_SSEE    = (5 * PHI**8 * RHO_CRIT_eV4)**0.25    # ≈ 9.62 meV (Paper 10)

# m_φ canónico por forward-derivation (Paper 6 partícula canónica), cero fiteo
SIGMA_MNU = (OMEGA / (KAL * TRIAL)) * 0.960318     # = 0.06903 eV
MULT      = OMEGA**4 + AURA * KAL                  # = 535.26
M_PHI_CAN = SIGMA_MNU * MULT                       # = 36.95 eV

# Escalas físicas SSEE (todas en eV)
scales = {
    "M_pl":           2.435e27,
    "M_UV":           M_UV_SSEE,
    "V0_quarter":     (0.840 * RHO_CRIT_eV4)**0.25,   # ≈ 2.35 meV
    "H0":             1.45e-33,
    "Sigma_mnu":      SIGMA_MNU,
    "m_e":            5.11e5,                      # electron mass (BSM ref)
    "T_CMB":          2.35e-4,                     # 2.725 K en eV
    "Lambda_QCD":     0.220,                       # 220 MeV
    "M_W":            8.04e10,                     # 80.4 GeV
}

TARGET = M_PHI_CAN  # eV (canónico 36.95)

print("=" * 76)
print(f"OP-10 FASE 2d — Búsqueda see-saw para m_φ canónico = {TARGET:.2f} eV")
print("(SUPERADA: m_φ ya tiene forward-derivation Σm_ν·MULT, cero fiteo)")
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
# TEST 4 — Origen canónico: forward-derivation (NO see-saw, NO H₀)
# ============================================================
print("\n" + "=" * 76)
print("TEST 4 — Origen canónico de m_φ: forward-derivation Σm_ν · MULT")
print("=" * 76)

# La masa φ-DM canónica (Paper 6 partícula canónica) NO sale de un see-saw ni
# de la vieja fórmula Σm_ν·H₀ (obsoleta). Sale de un producto algebraico limpio
# con cero parámetros de ajuste:
#     m_φ = Σm_ν · MULT
#     Σm_ν = (Ω/(KAL·TRIAL))·0.960318 = R2·0.960318
#     MULT = Ω⁴ + AURA·KAL
print(f"""
Derivación canónica (cero fiteo, Lagrangiano escalar libre — cierra OP-9):
  Σm_ν    = (Ω/(KAL·TRIAL))·0.960318 = {SIGMA_MNU:.5f} eV
  MULT    = Ω⁴ + AURA·KAL            = {MULT:.4f}
  m_φ     = Σm_ν · MULT              = {M_PHI_CAN:.4f} eV   ✓ canónico

Esta es una predicción algebraica dimensionalmente consistente: una masa (eV)
multiplicada por un número adimensional (MULT). NO depende de la convención de
unidades de H₀. El see-saw de abajo es exploración histórica — el origen de la
masa ya está fijado por la forward-derivation, no por estos productos.
""")

print("=" * 76)
print("RESUMEN")
print("=" * 76)

if len(best_candidates) > 0:
    print(f"\n  ✓ {len(best_candidates)} candidatos see-saw dentro de ±5% del target:")
    for name, val, ratio in best_candidates[:8]:
        print(f"    {name} = {val:.4e}  (ratio = {ratio:.4f})")
else:
    print(f"\n  ✗ Ningún see-saw simple (de 2 escalas × factor algebraico) alcanza")
    print(f"    m_φ = {TARGET:.2f} eV dentro de 5%.")
    print()
    print("  Resultado esperado y SIN consecuencia: el origen de m_φ NO es see-saw.")
    print("  La masa canónica ya está fijada por forward-derivation limpia")
    print(f"  (m_φ = Σm_ν·MULT = {M_PHI_CAN:.2f} eV, cero fiteo, Lagrangiano escalar")
    print("  libre que cierra OP-9). Este script queda como nota histórica.")
