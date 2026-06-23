#!/usr/bin/env python3
"""
OP-12 — Rama TÉRMICA: ¿qué época de desacople exige T_φ/T_ν = 0.5385?
====================================================================
Contexto (Mike, 2026-06-20). La partícula canónica φ-DM (m_φ = SOLAR²·KRYSTOS·Σm_ν
= 41.02 eV) entra en CLASS como especie TÉRMICA (ncdm) con temperatura T_φ. Pero esa
T_φ hoy se DESPEJA de la fórmula del relic usando m_φ y Ω_φDM — es bookkeeping, no
una temperatura derivada. Este script pregunta la versión NO circular:

  Si φ es un relic térmico, ¿en qué ÉPOCA tuvo que desacoplarse para quedar a
  T_φ/T_ν = 0.5385 hoy?  Respuesta vía dilución de entropía:
        (T_φ/T_ν)³ = g_*s(ν-dec=10.75) / g_*s(T_dec de φ)

Si esa época es física y natural, la rama térmica es viable y vale calcular el
freeze-out (Paso 2/3). Si cae en un lugar absurdo, la rama térmica muere aquí.

INSIGHT de Mike: la VENTANA (portal) por la que φ se desacopla NO es nueva — es la
VISCOSIDAD que ya está en el modelo (KAL₀, ζ̃=KAL₀/3, Paper 5 IS). Disipación =
interacción con un baño; un escalar libre no tiene viscosidad. SOLAR=BIAL+KAL hereda
ese rol radiativo-disipativo → señala que el baño es la RADIACIÓN.

CAVEAT honesto: esto fija la ÉPOCA requerida (robusto, solo entropía). NO calcula el
freeze-out (Paso 2/3, necesita estructura del portal). Y SOLAR≈7.9 NO es un
acoplamiento perturbativo → el portal es casi seguro suprimido por escala (operador
dim>4); el desconocido real es la escala Λ del portal.
"""
import numpy as np
PHI=1.618033988749895; PI=np.pi; OMEGA=PHI+PI

# --- T_φ/T_ν que exige el modelo (relic; m_φ=41.02, Ω_φDM h²=0.06877) ---
RATIO = 0.5385
GSTAR_NU = 10.75                       # g_*s a desacople de neutrinos
gstar_phi = GSTAR_NU/RATIO**3          # g_*s exigido al desacople de φ

# --- tabla g_*s(T) estándar SM (entropía); Husdal 2016 / Kolb-Turner; T en GeV ---
T  = np.array([1e3, 1e2, 10.0, 1.0, 0.30, 0.20, 0.15, 0.10, 0.020, 0.001])
gs = np.array([106.75,96.0,86.0,75.0,73.0, 68.0, 45.0, 20.0, 10.75, 3.91])
order = np.argsort(gs)
T_dec = np.interp(gstar_phi, gs[order], T[order])

print("="*78)
print("  OP-12 — rama TÉRMICA: época de desacople exigida por T_φ/T_ν=0.5385")
print("="*78)
print(f"  T_φ/T_ν = {RATIO}  ->  g_*s(desacople de φ) = {gstar_phi:.2f}")
print("-"*78)
print("  T (MeV)   g_*s")
for t,g in zip(T,gs):
    mark = "  <== desacople de φ" if abs(g-gstar_phi)<4 else ""
    print(f"  {t*1000:9.1f}  {g:6.2f}{mark}")
print("-"*78)
print(f"  => T_dec ~ {T_dec*1000:.0f} MeV  (escala QCD; Λ_QCD ~ 200 MeV)")
print("-"*78)
print("  VEREDICTO:")
print("  · g_*s~69 cae en la transición QCD (~200 MeV): época FÍSICA, no absurda.")
print("  · La rama térmica es VIABLE → predicción falsable (cuándo se desenchufa φ).")
print(f"  · Coincidencia a VIGILAR (no afirmar): g_*s={gstar_phi:.1f} vs H_alg=3Ω²={3*OMEGA**2:.2f}")
print("  · Pendiente Paso 2/3: portal = viscosidad KAL (Paper 5); freeze-out forward.")
print("    Obstáculo honesto: SOLAR≈7.9 no es coupling perturbativo → portal")
print("    suprimido por escala; el desconocido real es la escala Λ del portal.")
print("="*78)
