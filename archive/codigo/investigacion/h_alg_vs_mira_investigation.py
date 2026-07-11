"""
H_alg vs H_MIRA — ¿Coincidencia numerológica o relación estructural?
=====================================================================

Pregunta planteada por Mike (2026-05-24 PM):

  H_alg  = 3(φ+π)² = 67.962 km/s/Mpc   (Type-P algebraico)
  H_MIRA =          67.068 km/s/Mpc    (Planck plik_lite + SSEE bg)
  ────────────────────────────────────
  Δ = 0.894 km/s/Mpc   (ratio = 1.01334, exceso = 1.34%)

¿Es esto:
  (A) Coincidencia numerológica afortunada — H_alg es solo numerología
      que casualmente cae a 1.3% del valor físico
  (B) Una relación estructural — el 1.34% corresponde a alguna combinación
      algebraica del modelo (corrección de loop, screening, factor de
      conversión)
  (C) Una predicción todavía por entender — el "exceso" del 1.34% es un
      observable real del modelo

Esta investigación barre sistemáticamente:
  1. Cuán precisa es la coincidencia H_alg ≈ H_phys (compararla con
     coincidencias triviales)
  2. Si el cociente H_alg/H_MIRA = 1.0133 coincide con alguna razón
     algebraica SSEE conocida (1+f_screen, 1+ε, 1+MIRA-2, etc.)
  3. Si la diferencia Δ corresponde a alguna escala física del modelo

NOTA HISTÓRICA (2026-06-11): investigación archivada con el H_MIRA de la época
(67.068 km/s/Mpc; canónico actual 67.037). La pregunta quedó CERRADA en el
ledger como V-L2-06: H_alg=3(φ+π)² como tasa física de background está
refutado (θ*=7.80σ, r_d=5.35σ) → la cercanía H_alg≈H_MIRA es coincidencia
Type-P. NO re-correr con valores nuevos — registro de procedencia.
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import (PHI, PI, OMEGA, BETA, KAL0, OMEGA_DE, OMEGA_M_DYN,
                       MIRA, AURA, H0_ALG, OMEGA_M_CMB_MIRA)

H_ALG  = 67.962
H_MIRA = 67.068
RATIO  = H_ALG / H_MIRA          # 1.01334
EXCESS = RATIO - 1.0              # 0.01334
DELTA  = H_ALG - H_MIRA           # 0.894 km/s/Mpc

# Otras escalas SSEE clave
F_SCREEN_IR  = (PI - PHI) / OMEGA**2     # 0.06725  (Paper 9)
F_SCREEN_UV  = 0.06952                    # Paper 10
ALPHA_K      = 0.40330                    # Paper 7
ALPHA_K_FULL = 0.41691                    # Paper 10
EPSILON      = 5.67e-4                    # X²/M⁴
MIRA_MINUS_2 = MIRA - 2.0                 # -0.00108
MIRA_DEFECT  = 2.0 - MIRA                 # +0.00108

print("=" * 76)
print(" H_alg vs H_MIRA — análisis estructural")
print("=" * 76)
print(f"\n H_alg  = 3(φ+π)² = {H_ALG:.4f} km/s/Mpc")
print(f" H_MIRA = Planck+SSEE bg = {H_MIRA:.4f} km/s/Mpc")
print(f" Δ      = {DELTA:.4f} km/s/Mpc")
print(f" Ratio  = {RATIO:.6f}")
print(f" Excess = {EXCESS:.6f}  (1.334%)")

# ─────────────────────────────────────────────────────────────────────────
# Test 1: ¿qué tan buena es la coincidencia comparada con coincidencias
# triviales aleatorias?
# ─────────────────────────────────────────────────────────────────────────
print("\n" + "─" * 76)
print(" Test 1: ¿Qué tan buena es la coincidencia?")
print("─" * 76)

# Diccionario de combinaciones algebraicas SSEE simples que dan números O(50-80)
candidates = {
    "3(φ+π)²"             : 3*(PHI+PI)**2,
    "3(π+φ)²"             : 3*(PI+PHI)**2,           # = same
    "4(π+φ)²"             : 4*(PI+PHI)**2,
    "(φ+π)·5φ"            : (PHI+PI)*5*PHI,
    "π·KAL₀·2.5"          : PI*KAL0*2.5,
    "φ²·KAL₀·5"           : PHI**2*KAL0*5,
    "2π·KAL₀·2"           : 2*PI*KAL0*2,
    "φ⁵·5/2"              : PHI**5*5/2,
    "Σ_Sov/π·15"          : (PHI+PI+(PHI+PI+OMEGA))/PI*15,  # M_v/π*15
    "OMEGA²·3"            : OMEGA**2*3,                # = H_alg!
    "AURA·KAL₀·3"         : AURA*KAL0*3,
    "π²·7"                : PI**2*7,
}
print(f"\n {'Combinación':<25} {'Valor':>12} {'Δ vs H_MIRA':>12} {'σ':>8}")
print(" " + "─"*68)
for name, val in sorted(candidates.items(), key=lambda x: abs(x[1]-H_MIRA)):
    delta = val - H_MIRA
    sigma = abs(delta)/0.54  # σ_H_MIRA ≈ 0.54 (Planck-tipo)
    flag = " ★" if name == "3(φ+π)²" else ""
    print(f" {name:<25} {val:12.4f} {delta:+12.4f} {sigma:8.2f}{flag}")

# ─────────────────────────────────────────────────────────────────────────
# Test 2: ¿el cociente 1.01334 coincide con alguna razón SSEE conocida?
# ─────────────────────────────────────────────────────────────────────────
print("\n" + "─" * 76)
print(" Test 2: ¿1.01334 es algo conocido?")
print("─" * 76)

ratio_candidates = {
    "1 + ε (UV correction)"             : 1 + EPSILON,
    "1 + f_screen^IR"                   : 1 + F_SCREEN_IR,
    "1 + f_screen^UV"                   : 1 + F_SCREEN_UV,
    "1 + (MIRA-2)"                      : 1 + MIRA_MINUS_2,
    "1 + (2-MIRA)"                      : 1 + MIRA_DEFECT,
    "1 + Ω_b·h²"                        : 1 + 0.02242,
    "1 + (π-φ)/Ω³"                      : 1 + (PI-PHI)/OMEGA**3,
    "1 + (π-φ)/(φ+π)³"                  : 1 + (PI-PHI)/(PHI+PI)**3,
    "1 + (Ω_m,CMB - 0.3)·0.13"          : 1 + (OMEGA_M_CMB_MIRA-0.3)*0.13,
    "(MIRA/2)·1.01334"                  : MIRA/2 * 1.01334,
    "1 + 1/(KAL₀·15)"                   : 1 + 1/(KAL0*15),
    "1 + 1/π⁴"                          : 1 + 1/PI**4,
    "1 + 1/(φ⁶·π)"                      : 1 + 1/(PHI**6*PI),
    "1 + φ/(π·KAL₀·22)"                 : 1 + PHI/(PI*KAL0*22),
    "1 + 4/(φ+π+KAL₀+OMEGA)²"           : 1 + 4/(PHI+PI+KAL0+OMEGA)**2,
    "1 + (αK·αK^UV)/12"                 : 1 + (ALPHA_K*ALPHA_K_FULL)/12,
    "(1-f_screen)/(1-(αK_full/3MIRA))"  : (1-F_SCREEN_IR)/(1-ALPHA_K_FULL/(3*MIRA)),
}

print(f"\n  {'Combinación':<40} {'Valor':>10} {'|Δ vs 1.01334|':>16}")
print("  " + "─"*68)
for name, val in sorted(ratio_candidates.items(), key=lambda x: abs(x[1]-RATIO)):
    delta = abs(val - RATIO)
    rel = delta / EXCESS * 100  # % of excess
    flag = " ★" if delta < 0.001 else (" ✓" if delta < 0.003 else "")
    print(f"  {name:<40} {val:10.5f} {delta:16.5f}{flag}")

# ─────────────────────────────────────────────────────────────────────────
# Test 3: ¿la DIFERENCIA Δ = 0.894 km/s/Mpc tiene sentido físico?
# ─────────────────────────────────────────────────────────────────────────
print("\n" + "─" * 76)
print(" Test 3: ¿Δ = 0.894 km/s/Mpc tiene sentido físico?")
print("─" * 76)

# Comparar con escalas físicas del modelo
diff_candidates = {
    "H_alg · f_screen^IR"               : H_ALG * F_SCREEN_IR,
    "H_alg · (1-(MIRA-2))"              : H_ALG * MIRA_DEFECT,
    "H_alg · ε (UV)"                    : H_ALG * EPSILON,
    "H_alg · Ω_b·h²/0.5"                : H_ALG * 0.02242/0.5,
    "H_MIRA · (αK_full - αK_IR)/αK_IR"  : H_MIRA * (ALPHA_K_FULL-ALPHA_K)/ALPHA_K,
    "H_MIRA · 0.01334 (= excess)"       : H_MIRA * EXCESS,
}

print(f"\n  {'Combinación':<40} {'Valor':>10} {'|Δ vs 0.894|':>14}")
print("  " + "─"*64)
for name, val in sorted(diff_candidates.items(), key=lambda x: abs(x[1]-DELTA)):
    delta = abs(val - DELTA)
    flag = " ★" if delta < 0.05 else (" ✓" if delta < 0.15 else "")
    print(f"  {name:<40} {val:10.4f} {delta:14.4f}{flag}")

# ─────────────────────────────────────────────────────────────────────────
# Diagnóstico
# ─────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 76)
print(" Diagnóstico estructural")
print("=" * 76)

best_ratio = min(ratio_candidates.items(), key=lambda x: abs(x[1]-RATIO))
best_diff  = min(diff_candidates.items(), key=lambda x: abs(x[1]-DELTA))
print(f"\n• Mejor match al ratio 1.01334:")
print(f"    {best_ratio[0]} = {best_ratio[1]:.5f}")
print(f"    desviación: {abs(best_ratio[1]-RATIO):.5f} = {100*abs(best_ratio[1]-RATIO)/EXCESS:.1f}% del exceso")

print(f"\n• Mejor match a la diferencia 0.894 km/s/Mpc:")
print(f"    {best_diff[0]} = {best_diff[1]:.4f}")
print(f"    desviación: {abs(best_diff[1]-DELTA):.4f} km/s/Mpc")

# Cuán probable es esto por azar
import math
# Si H_alg fuera "número aleatorio en rango cosmológico [60, 80]", prob de caer
# a Δ<1.5% de H_MIRA por azar
prob_random = (2 * 1.5/100 * H_MIRA) / 20   # ventana ±1.5% en rango [60,80]
print(f"\n• Probabilidad de coincidencia aleatoria (H ∈ [60,80] km/s/Mpc):")
print(f"    P(|H_alg - H_phys|/H < 1.5%) ≈ {prob_random*100:.1f}%")
print(f"    → la coincidencia NO es trivial pero tampoco extraordinaria")

print("\n• Interpretación:")
print("  - Las tres formas algebraicas que dan ~H_phys (3(φ+π)², OMEGA²·3, etc.)")
print("    son TODAS la misma identidad — H_alg = 3·OMEGA² es única.")
print("  - El exceso 1.34% no coincide limpiamente con ninguna razón estructural")
print("    conocida del modelo (f_screen, ε, MIRA-2, etc. están todos en escala")
print("    diferente).")
print("  - Sugiere que H_alg es una COINCIDENCIA NUMEROLÓGICA AFORTUNADA (Type P),")
print("    no una predicción estructural. El valor físico es H_MIRA.")
print("  - El modelo no DEBE H_alg como predicción; lo retiene como signal de")
print("    que la combinación 3(φ+π)² está cerca de la escala correcta — pista")
print("    sin teorema.")
print("\n" + "=" * 76)
