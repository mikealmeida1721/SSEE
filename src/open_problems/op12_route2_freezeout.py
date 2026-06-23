#!/usr/bin/env python3
"""
OP-12 Ruta 2 — ¿qué acoplamiento de portal exige el freeze-out a T_dec≈217 MeV?
================================================================================
Bridge 2 (viscosidad↔cross-section), parte calculable: la condición de desacople
Γ = n_target·σ·v = H(T_dec) FIJA la sección eficaz σ requerida, y de ahí el
acoplamiento (dim-4: g) o la escala del portal (dim-5: Λ).

Pregunta clave (correctness burden): ¿es ese acoplamiento = SOLAR (g²v de la masa),
como esperaba Mike? ¿O el portal es un acoplamiento DISTINTO y débil?

Caveat honesto: prefactores (4π, g_target, forma del operador) son orden-de-magnitud;
el RESULTADO ROBUSTO es la ESCALA del acoplamiento, no su 3ª cifra.
"""
import numpy as np
PHI=1.618033988749895; PI=np.pi
BIAL=(PHI+PI)/2; KAL=BIAL+PI; SOLAR=BIAL+KAL   # 7.9012

MPL = 1.2209e19          # GeV, masa de Planck (full)
T   = 0.217              # GeV, T_dec (escala QCD)
gstar = 69.0             # g_* en el desacople
GEV2_TO_MB = 0.389       # 1 GeV^-2 = 0.389 mb

# --- H en el desacople ---
H = 1.66*np.sqrt(gstar)*T**2/MPL          # GeV

# --- densidad del baño (objetivos con que choca phi): n ~ 0.12 g_* T^3 ---
n_target = 0.1218*gstar*T**3              # GeV^3 (suma de especies relativistas)

# --- sigma requerida para Γ=H (v≈1, relativista) ---
sigma_req = H/n_target                     # GeV^-2

# --- traducir a acoplamiento ---
# dim-4 (renormalizable): sigma ~ g^4/(8π T^2)
g_dim4 = (sigma_req*8*PI*T**2)**0.25
# dim-5 (suprimido por escala): sigma ~ T^2/(4π Λ^4)
Lambda_dim5 = (T**2/(4*PI*sigma_req))**0.25  # GeV

print("="*78)
print("  OP-12 Ruta 2 — acoplamiento de portal exigido por freeze-out @ 217 MeV")
print("="*78)
print(f"  T_dec={T*1000:.0f} MeV   g_*={gstar:.0f}")
print(f"  H(T_dec)        = {H:.3e} GeV")
print(f"  n_target        = {n_target:.3e} GeV^3")
print(f"  sigma requerida = {sigma_req:.3e} GeV^-2  = {sigma_req*GEV2_TO_MB:.2e} mb")
print("-"*78)
print("  Acoplamiento implicado (orden de magnitud, prefactores ~O(1)):")
print(f"    · portal dim-4 (renormalizable): g ~ {g_dim4:.2e}")
print(f"    · portal dim-5 (escala):         Λ ~ {Lambda_dim5:.3e} GeV  (~{Lambda_dim5/1000:.1f} TeV)")
print("-"*78)
print(f"  COMPARAR con SOLAR = φ+2π = {SOLAR:.4f}  (acoplamiento de la MASA, g²v)")
print("  VEREDICTO (honesto):")
print(f"    · El portal NO es SOLAR. g_portal~{g_dim4:.0e} es DÉBIL (vs SOLAR~7.9).")
print("      => La masa (SOLAR) y el portal son acoplamientos DISTINTOS. La esperanza")
print("         'el mismo SOLAR hace masa Y portal' NO cierra numéricamente.")
print("    · PERO acoplamiento DÉBIL <=> viscosidad ALTA (Bridge 2, signo correcto):")
print("      un portal débil es CONSISTENTE con que KAL sea una viscosidad grande.")
print("    · La escala dim-5 Λ~pocos TeV es falsable (huella colisionador/indirecta).")
print("="*78)
