"""
OP-8 Fase 2 — λ_P7 como coupling IDE
=====================================

HIPÓTESIS (tras descartar BIAL en op8_coupled_dm_de.py):
  La cuasi-identidad ln(MIRA) ≈ λ_P7 = √(3·Ω_m,dyn) sugiere que MIRA
  NO es un postulado independiente sino una CONSECUENCIA del tracker
  coupling λ en V(φ) (Paper 7).

  Si exacto:  MIRA = exp(λ_P7) = exp(√(3·Ω_m,dyn))

  Implicación: OP-8 se resuelve estructuralmente; MIRA deja de ser
  postulado auxiliar (M) y pasa a corolario de Paper 7.

TESTS:
  1. Verificar identidad numérica ln(MIRA) ?= √(3·Ω_m,dyn)
  2. ¿Cuál es la diferencia relativa? ¿es 0 a precisión machine
     o es coincidencia numérica?
  3. Si NO es identidad exacta: ¿hay corrección algebraica natural
     que la cierre?
  4. ¿Qué predice modelo IDE con β = λ_P7 para ratio Ω_m,CMB/Ω_m,today?
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import (PHI, PI, OMEGA, BETA, KAL0, OMEGA_DE, OMEGA_M_DYN,
                       MIRA, AURA)

print("=" * 78)
print(" OP-8 Fase 2 — λ_P7 = √(3·Ω_m,dyn) como coupling IDE")
print("=" * 78)

# ─── Test 1: identidad numérica ─────────────────────────────────────
lambda_P7 = np.sqrt(3 * OMEGA_M_DYN)
ln_MIRA = np.log(MIRA)
exp_lambda = np.exp(lambda_P7)

print(f"\n [Test 1] Cuasi-identidad ln(MIRA) ?= λ_P7")
print(f"   MIRA              = {MIRA:.10f}")
print(f"   ln(MIRA)          = {ln_MIRA:.10f}")
print(f"   λ_P7 = √(3Ω_m,dyn) = {lambda_P7:.10f}")
print(f"   Δ absoluta        = {ln_MIRA - lambda_P7:+.10f}")
print(f"   Δ relativa        = {(ln_MIRA - lambda_P7)/lambda_P7*100:+.6f} %")

print(f"\n   exp(λ_P7)         = {exp_lambda:.10f}")
print(f"   MIRA observado    = {MIRA:.10f}")
print(f"   Δ                 = {exp_lambda - MIRA:+.10f}")
print(f"   Δ relativa        = {(exp_lambda - MIRA)/MIRA*100:+.6f} %")

# ─── Test 2: ¿identidad exacta o coincidencia? ──────────────────────
print(f"\n [Test 2] ¿Identidad algebraica exacta?")
print(f"   Si MIRA = (3φ+π)/4 y Ω_m,dyn = M_v/(3·algo), debe haber")
print(f"   una relación trigonométrica/algebraica exacta.")
print(f"")
print(f"   Ω_m,dyn algebraico = {OMEGA_M_DYN:.10f}")
print(f"   3·Ω_m,dyn          = {3*OMEGA_M_DYN:.10f}")
print(f"   (ln MIRA)²         = {ln_MIRA**2:.10f}")
print(f"   Δ                  = {3*OMEGA_M_DYN - ln_MIRA**2:+.10f}")
print(f"")
print(f"   → Si 3·Ω_m,dyn = (ln MIRA)² exacto: identidad estructural")
print(f"   → Si no: cuasi-identidad numérica (Type-P)")

# ─── Test 3: ¿qué necesita corrección? ─────────────────────────────
# Reverso: si ASUMIMOS MIRA = exp(√(3Ω_m,dyn)), ¿qué Ω_m,dyn requiere?
Omega_m_required = ln_MIRA**2 / 3.0
print(f"\n [Test 3] Si MIRA = exp(√(3Ω_m)) EXACTO, ¿qué Ω_m necesitamos?")
print(f"   Ω_m,required = (ln MIRA)²/3 = {Omega_m_required:.10f}")
print(f"   Ω_m,dyn SSEE = {OMEGA_M_DYN:.10f}")
print(f"   Δ            = {Omega_m_required - OMEGA_M_DYN:+.10f}")
print(f"   Δ relativa   = {(Omega_m_required - OMEGA_M_DYN)/OMEGA_M_DYN*100:+.4f} %")

# ─── Test 4: predicción IDE con β = λ_P7 ───────────────────────────
print(f"\n [Test 4] Modelo IDE con β = λ_P7:")
print(f"   En aprox 'ratio = 1 + β':")
print(f"     ratio = 1 + λ_P7 = {1 + lambda_P7:.6f}")
print(f"     MIRA observado  = {MIRA:.6f}")
print(f"     → β = λ_P7 NO da ratio = MIRA bajo aprox lineal")
print(f"")
print(f"   En aprox exponencial 'ratio = exp(β·factor)':")
print(f"     Si factor=1: ratio = exp(λ_P7) = {exp_lambda:.6f}")
print(f"     MIRA observado  = {MIRA:.6f}")
print(f"     Δ relativa      = {(exp_lambda-MIRA)/MIRA*100:+.4f} %")
print(f"     → ESTO SÍ MATCHEA dentro de 0.05%")
print(f"")
print(f"   FORMA FUNCIONAL natural:")
print(f"     ρ_DM(a)/ρ_DM,0 = exp(λ_P7 · |ln a|) = a^(-λ_P7) ... no")
print(f"     Pero ρ_m,CMB/ρ_m,today = (1+z_CMB)^δ con δ pequeño da")
print(f"     ratio enorme. La aprox exponencial requiere mecanismo")
print(f"     que NO sea power-law en a.")

# ─── Test 5: ¿hay corrección algebraica natural? ────────────────────
print(f"\n [Test 5] Δ relativa 0.03% — ¿corrección natural?")
delta = MIRA - exp_lambda
print(f"   Δ = MIRA - exp(λ_P7) = {delta:+.8f}")
print(f"   Δ/MIRA = {delta/MIRA:+.6e}")
print(f"")
candidates_correction = {
    "1/(2π·M_v)": 1/(2*PI*14.2788),
    "1/(φ⁸)":     1/PHI**8,
    "(φ-1)/φ⁸":  (PHI-1)/PHI**8,
    "1/Ω⁴":       1/OMEGA**4,
    "1/(KAL₀·M_v)": 1/(KAL0*14.2788),
}
print(f"   {'Candidato':<20} {'Valor':>14} {'|Δ - cand|':>14}")
print(f"   {'-'*50}")
for name, val in sorted(candidates_correction.items(),
                        key=lambda x: abs(x[1] - delta)):
    print(f"   {name:<20} {val:>14.8f} {abs(val-delta):>14.8f}")

# ─── Conclusión ─────────────────────────────────────────────────────
print(f"\n{'='*78}")
print(" CONCLUSIONES")
print(f"{'='*78}")
rel_diff_pct = abs((exp_lambda - MIRA)/MIRA*100)
print(f"""
 1. La cuasi-identidad MIRA ≈ exp(√(3·Ω_m,dyn)) se cumple a {rel_diff_pct:.4f}%.

 2. Esto es MUY cercano pero NO es identidad algebraica exacta a
    precisión machine. Es del nivel "Type-P coincidence" (similar a
    H_alg = 3(φ+π)² coincidencia numérica con H_MIRA físico).

 3. Si fuera estructural, requeriría: (ln MIRA)² = 3·Ω_m,dyn exacto.
    Numéricamente: {ln_MIRA**2:.6f} vs {3*OMEGA_M_DYN:.6f}
    Δ = {3*OMEGA_M_DYN - ln_MIRA**2:+.6f} (no cero).

 4. INTERPRETACIÓN HONESTA:
    - MIRA y λ_P7 ambos derivan de las mismas constantes (φ,π) y
      Ω_m,dyn (que sale de M_v = φ+π+K_v).
    - Por tanto NO son independientes — comparten origen estructural.
    - La cuasi-identidad refleja esta dependencia común, no una
      ley física nueva.

 5. PARA OP-8 (mecanismo MIRA):
    - λ_P7 NO es el coupling IDE que reproduce MIRA via aprox lineal.
    - La aprox exponencial sí coincide pero requiere forma funcional
      que NO emerge naturalmente de Q = β·H·ρ_DM estándar.
    - OP-8 sigue ABIERTO. La conexión MIRA-V(φ) es sugestiva pero
      requiere Boltzmann completo con IDE para verificar.

 6. PRÓXIMO PASO recomendado:
    - Verificar simbólicamente si MIRA = (3φ+π)/4 y Ω_m,dyn = T_r/M_v
      tienen relación trigonométrica que cierre la cuasi-identidad.
    - Si no: catalogar como Type-P coincidence más, no avanzar OP-8
      por esta ruta sin Boltzmann.
""")
