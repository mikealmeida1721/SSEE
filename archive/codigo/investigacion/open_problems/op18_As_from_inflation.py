#!/usr/bin/env python3
"""
OP-18 — ¿Sale la amplitud primordial A_s de (φ,π)? — cálculo honesto
=====================================================================
Tras el reframe ω_m-directo, el sector CMB de SSEE es k=2={A_s,τ}. OP-18 ataca A_s.

α-attractor (SSEE, ya fijado): α=φ⁴/3, N_*=2φ⁷, n_s=1−φ⁻⁷, r=φ⁻¹⁰ (OP-2).
Slow-roll en el plateau: ε = 3α/(4N_*²) = φ⁻¹⁰/16  (⇒ r=16ε=φ⁻¹⁰ ✓).
Amplitud escalar:  A_s = V/(24π²ε M_Pl⁴) = (2φ¹⁰/3π²)·(V₀/M_Pl⁴).

⇒ El PREFACTOR (2φ¹⁰/3π²) es puro (φ,π). La NORMALIZACIÓN absoluta depende de la
escala de inflación V₀ — y AHÍ está la pregunta real de OP-18.
"""
import numpy as np
PHI=1.618033988749895; PI=np.pi
MPL_RED=2.435e18  # GeV, masa de Planck reducida

# --- ingredientes SSEE (todos φ,π) ---
alpha   = PHI**4/3
N_star  = 2*PHI**7
ns      = 1 - PHI**-7
r       = PHI**-10
eps     = PHI**-10/16
prefac  = 2*PHI**10/(3*PI**2)        # A_s = prefac · V0/Mpl^4

# --- A_s medido (Planck 2018) ---
As_obs  = np.exp(3.044)*1e-10        # ln(1e10 As)=3.044

# --- dado A_s, derivar la escala de inflación V0 ---
V0_over_Mpl4 = As_obs/prefac
V0_quart_GeV = (V0_over_Mpl4)**0.25 * MPL_RED   # V0^(1/4) en GeV

print("="*76)
print("  OP-18 — ¿A_s desde (φ,π)?  cálculo honesto (α-attractor SSEE)")
print("="*76)
print(f"  α=φ⁴/3={alpha:.4f}   N_*=2φ⁷={N_star:.3f}   n_s=1−φ⁻⁷={ns:.5f}   r=φ⁻¹⁰={r:.3e}")
print(f"  ε=φ⁻¹⁰/16={eps:.3e}   prefactor 2φ¹⁰/3π² = {prefac:.5f}  (PURO φ,π)")
print("-"*76)
print(f"  A_s(Planck) = {As_obs:.3e}")
print(f"  ⇒ V₀/M_Pl⁴   = A_s·3π²/2φ¹⁰ = {V0_over_Mpl4:.3e}")
print(f"  ⇒ V₀^(1/4)   = {V0_quart_GeV:.3e} GeV  (~10^16 GeV, escala GUT/inflación típica)")
print("-"*76)

# --- ¿la jerarquía V0/Mpl^4 sale limpia de (φ,π)? (test negativo honesto) ---
print("  ¿La jerarquía V₀/M_Pl⁴ ≈ {:.3e} sale de (φ,π)?".format(V0_over_Mpl4))
target=V0_over_Mpl4
best=[]
for base,bn in [(PHI,'φ'),(PI,'π'),(PHI+PI,'Ω')]:
    for n in range(-80,0):
        val=base**n
        if val>0 and 0.3<val/target<3.0:
            best.append((abs(np.log(val/target)),f"{bn}^{n}",val))
best.sort()
if best:
    for d,lab,val in best[:3]:
        print(f"    {lab:8s} = {val:.3e}   (factor {val/target:.2f} del objetivo)")
else:
    print("    ninguna potencia simple φ,π,Ω cae dentro de factor 3 del objetivo")
print("-"*76)
print("  VEREDICTO (honesto):")
print("  · SSEE fija la FORMA del espectro primordial por completo: n_s, r, α, N_*,")
print("    y el prefactor 2φ¹⁰/3π² de A_s — todo (φ,π), sin perillas.")
print("  · La NORMALIZACIÓN A_s requiere la escala de inflación V₀^(1/4)~10^16 GeV.")
print("    Esa jerarquía NO sale limpia de potencias (φ,π) — es el residuo real de OP-18,")
print("    análogo al problema de jerarquía (UV).")
print("  · CORRECCIÓN a una nota previa: la escala Λ_SSEE=M=9.68 meV (Paper 10) es IR")
print("    (energía oscura); NO aplica a V₀ (UV, inflación) — están separadas ~10^27.")
print("    Anclar A_s a M=9.68 meV es INCORRECTO. OP-18 sigue abierto de verdad.")
print("  · Lo defendible HOY: A_s es 1 input de escala; lo que SSEE aporta es que la")
print("    FORMA (n_s,r) ya está fijada, y dado A_s la escala de inflación cae en el")
print("    rango GUT esperado (check de consistencia, no derivación).")
print("="*76)
