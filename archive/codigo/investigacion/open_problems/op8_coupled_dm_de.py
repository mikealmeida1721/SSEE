"""
OP-8 — Interpretación IDE (Interacting Dark Energy / Coupled DM-DE)
=====================================================================

HIPÓTESIS (Mike, 2026-05-26):
  MIRA = 2 NO indica un segundo sector DM con partícula φ.
  Indica una TRANSFERENCIA de energía DM → DE entre CMB y hoy.

  Como agua → vapor: misma sustancia, dos comportamientos según
  época cósmica (temperatura).

MECANISMO PROPUESTO:
  Acoplamiento β entre DM y DE en el Lagrangiano SSEE.
  ρ_DM(a) ∝ a^(-3-β)  →  matter MÁS abundante en early universe

  Si β derivable algebraicamente de (φ,π) → predicción limpia.

TESTS:
  1. ¿Qué β predice MIRA = 2 = Ω_m,CMB/Ω_m,today?
  2. ¿Coincide ese β con algún algebraico SSEE natural?
  3. ¿BIAL = (φ+π)/2 ≈ 2.38 es candidato? (sugerencia Mike)
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import (PHI, PI, OMEGA, BETA, KAL0, OMEGA_DE, OMEGA_M_DYN,
                       MIRA, AURA, T_R, M_V, W0)
from scipy.optimize import brentq

# ── Parámetros base SSEE ─────────────────────────────────────────────
OMEGA_M_TODAY = OMEGA_M_DYN          # 0.160050 (algebraico)
OMEGA_DE_TODAY = OMEGA_DE            # 0.840
OMEGA_M_CMB_OBSERVED = MIRA * OMEGA_M_TODAY  # 0.31993 (lo que Planck infiere)

# z de recombinación
Z_CMB = 1100.0
A_CMB = 1.0 / (1.0 + Z_CMB)

print("=" * 78)
print(" OP-8 — IDE Coupled DM-DE — ¿es esto MIRA?")
print("=" * 78)
print(f"\n Setup:")
print(f"   Ω_m,today  = {OMEGA_M_TODAY:.6f}  (SSEE algebraico)")
print(f"   Ω_m,CMB    = {OMEGA_M_CMB_OBSERVED:.6f}  (MIRA × Ω_m,today)")
print(f"   MIRA       = {MIRA:.6f}")
print(f"   z_CMB      = {Z_CMB}")

# ─── MODELO IDE ──────────────────────────────────────────────────────
# Acoplamiento simple: ρ_DM(a) = ρ_DM,0 · a^(-(3+β))
#                     ρ_DE(a) compleja, requiere conservación energía
#
# Para β pequeño y DE casi cte:
#   En CMB: ρ_DM,eff = ρ_DM,0 × (1+z)^(3+β)
#   En hoy: ρ_DM,eff = ρ_DM,0
#
# Pero lo que CMB analysis INFIERE es Ω_m,today (asumiendo ΛCDM).
# La inferencia es Ω_m,inferred ≠ Ω_m,actual_today.
#
# Aprox simple: el ratio observado Ω_m,CMB / Ω_m,today está dado por
# cuánto extra contribuye el φ-matter-mode al inicio.

def omega_m_ratio_for_coupling(beta, Omega_m_today=OMEGA_M_TODAY,
                                Omega_de_today=OMEGA_DE_TODAY):
    """Ω_m,inferred(CMB) / Ω_m,today para acoplamiento β.

    Modelo: Q = βHρ_DM → ρ_DM(a) ∝ a^(-3-β)

    En CMB analysis, las anisotropías sensan Ω_m,h² efectivo.
    La razón Ω_m,CMB/Ω_m,today emerge porque φ contribuye como
    matter en early universe.

    Para acoplamiento estable, el ratio limit es:
        Ω_m,CMB / Ω_m,today ≈ (1+z_CMB)^β / (1+z_CMB)^β_pure

    Lo simplifico al ratio de cómo Ω_m_eff escala vs estándar:
    """
    # Aproximación: en CMB, la "matter" efectiva es la suma de DM real
    # más la energía que el φ tenía como matter antes de transicionar.
    # En modelo Wetterich-like: la transición ocurre cuando ρ_DE alcanza
    # 1/3 del total. Antes: φ comporta matter-like. Después: DE-like.

    # En CMB (z=1100), si φ todavía está en matter mode:
    #   Ω_m,eff,CMB = Ω_DM,real + Ω_φ,matter_mode
    # Si la transición es completa antes de hoy:
    #   Ω_m,today = Ω_DM,real

    # El factor β controla cuánta energía se transfiere.
    # Aproximación simple: Ω_φ,matter / Ω_DM,real ≈ β·tanh(z/z_trans)·factor

    # Para β=0: ratio = 1 (no transferencia, ΛCDM)
    # Para β grande: ratio → 1 + algo grande

    # NUMÉRICAMENTE simple: ratio ≈ 1 + β·factor(a_trans)
    # Si factor ≈ 1 para z_CMB ≫ z_trans: ratio ≈ 1 + β

    # Esto es muy aproximado. La realidad requiere CLASS/Boltzmann.
    # Por ahora: ratio_simple = 1 + β
    return 1.0 + beta

# Pregunta 1: ¿qué β predice MIRA = 2?
target_ratio = MIRA  # 1.99892
def diff(beta):
    return omega_m_ratio_for_coupling(beta) - target_ratio
beta_required = brentq(diff, -0.5, 10.0, xtol=1e-6)
print(f"\n [Pregunta 1] β requerido para ratio Ω_m,CMB/Ω_m,today = MIRA = 2:")
print(f"   β_required = {beta_required:.6f}")
print(f"   (aproximación: ratio ≈ 1 + β; realidad necesita CLASS/Boltzmann)")

# Pregunta 2: ¿qué algebraicos SSEE están cerca?
print(f"\n [Pregunta 2] Candidatos algebraicos SSEE cerca de β ≈ {beta_required:.3f}:")
candidates_algebraic = {
    "MIRA - 1":              MIRA - 1,                  # = 0.99892
    "ln(MIRA)":              np.log(MIRA),              # = 0.6926
    "Ω_DE":                  OMEGA_DE,                  # = 0.840
    "λ_P7 = √(3·Ω_m,dyn)":   np.sqrt(3*OMEGA_M_DYN),    # = 0.693
    "1/φ":                   1/PHI,                     # = 0.618
    "1/(φ+1)":               1/(PHI+1),                 # = 0.382
    "(π-φ)·MIRA":            (PI-PHI)*MIRA,             # = 3.046
    "BIAL = (φ+π)/2":        (PHI+PI)/2,                # = 2.380
    "AURA = 2·MIRA":         AURA,                      # = 3.998
}

print(f"\n   {'Candidato':<28} {'Valor':>10} {'|Δ vs β_req|':>14}")
print(f"   {'-'*54}")
for name, val in sorted(candidates_algebraic.items(),
                        key=lambda x: abs(x[1] - beta_required)):
    diff_val = abs(val - beta_required)
    flag = " ★" if diff_val < 0.02 else (" ✓" if diff_val < 0.1 else "")
    print(f"   {name:<28} {val:>10.4f} {diff_val:>14.4f}{flag}")

# Pregunta 3: BIAL específicamente
print(f"\n [Pregunta 3] β = BIAL = (φ+π)/2 ≈ 2.38 (sugerencia Mike):")
beta_BIAL = (PHI + PI) / 2
ratio_BIAL = omega_m_ratio_for_coupling(beta_BIAL)
print(f"   β_BIAL = {beta_BIAL:.6f}")
print(f"   ratio resultante = 1 + β = {ratio_BIAL:.4f}")
print(f"   MIRA observado   = {MIRA:.4f}")
print(f"   Diferencia       = {ratio_BIAL - MIRA:+.4f}")
print(f"   → BIAL produce ratio {ratio_BIAL:.2f}, NO 2 (MIRA).")
print(f"   → BIAL probablemente NO es el coupling correcto.")

# Pregunta 4: ¿quién está MÁS cerca?
print(f"\n [Pregunta 4] Análisis del mejor candidato:")
print(f"\n   MEJOR CANDIDATO ALGEBRAICO:")
sorted_cands = sorted(candidates_algebraic.items(),
                      key=lambda x: abs(x[1] - beta_required))
best_name, best_val = sorted_cands[0]
print(f"     {best_name} = {best_val:.6f}")
print(f"     β_required = {beta_required:.6f}")
print(f"     coincidencia: {100*(1 - abs(best_val-beta_required)/beta_required):.2f}%")

# Pregunta 5: ¿es esto suficiente para reformular OP-8?
print(f"""
 [Pregunta 5] ¿Reformular OP-8 como IDE con β = {best_name}?

 OBSERVACIONES:
   - β_required ≈ 1 es el ratio aproximado MIRA - 1 = 0.999
   - Si β = MIRA - 1 = 0.99892, el ratio simple da exactamente MIRA
   - Pero es CIRCULAR: usar MIRA para predecir MIRA

   - ln(MIRA) ≈ 0.693 también está cerca y ES STRUCTURAL
   - Pero ln(MIRA) coincide casi exactamente con λ_P7 = √(3·Ω_m,dyn)
     que ya está derivado en P7 desde tracker condition

   - BIAL = (φ+π)/2 = 2.38 da ratio ≈ 3.4, muy LEJOS de MIRA = 2
   - BIAL NO es el coupling para esta interpretación

 SOSPECHA INTERESANTE:
   ln(MIRA) ≈ λ_P7  (ambos ≈ 0.693)

   ¿Coincidencia?
""")
print(f"   ln(MIRA) = ln({MIRA:.6f}) = {np.log(MIRA):.6f}")
print(f"   λ_P7     = √(3·{OMEGA_M_DYN:.6f}) = {np.sqrt(3*OMEGA_M_DYN):.6f}")
print(f"   Diferencia: {np.log(MIRA) - np.sqrt(3*OMEGA_M_DYN):+.6f}")
print(f"")
print(f"   Si esto fuera identidad exacta, sería:")
print(f"     ln(MIRA) = λ_P7")
print(f"     MIRA = exp(λ_P7) = exp(√(3·Ω_m,dyn))")
print(f"")
print(f"   Pero numéricamente:")
print(f"     exp(√(3·Ω_m,dyn)) = {np.exp(np.sqrt(3*OMEGA_M_DYN)):.6f}")
print(f"     MIRA observado    = {MIRA:.6f}")
print(f"     Diferencia: {np.exp(np.sqrt(3*OMEGA_M_DYN)) - MIRA:+.6f}")

print(f"\n{'='*78}")
print(" SIGUIENTE PASO recomendado")
print(f"{'='*78}")
print("""
 1. El modelo IDE ASUME forma de coupling simple (Q = βHρ_DM).
    La realidad podría requerir Q = f(φ)·Hρ_DM con f(φ) algebraico.

 2. La cercanía ln(MIRA) ≈ λ_P7 sugiere POSIBLE relación estructural
    entre MIRA (matter mapping) y λ (tracker coupling en V(φ)).
    Esta es la dirección a investigar.

 3. Para verificar: construir Lagrangiano IDE con SSEE V(φ),
    integrar evolución cosmológica con coupling, comparar con CLASS
    output para ver si reproduce MIRA = 2 naturalmente.

 4. Si funciona: MIRA es consecuencia del coupling en V(φ), NO un
    postulado independiente. OP-8 se resuelve estructuralmente.
""")
