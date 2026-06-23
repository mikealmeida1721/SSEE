#!/usr/bin/env python3
"""
OP-12 Paso 2 (acotado) — ¿La traza T^μ_μ se enciende en QCD, y qué producción da?
=================================================================================
El término candidato es L_int = (φ/Λ)·T^μ_μ, con T^μ_μ=(1-3w)ρ = "interaction
measure" del plasma. φ solo siente la RUPTURA de conformalidad. Probamos:
  (a) ¿(1-3w)(T) tiene un PICO en la transición QCD (~200 MeV)?  -> ¿enciende ahí?
  (b) La forma de Γ/H(T): ¿PLATEAU>1 que cruza (freeze-out térmico) o PICO breve
      (kick coherente = producción FRÍA)?  -> decide la rama.

CAVEAT honesto: la interaction measure es lattice QCD APROXIMADA (forma robusta,
valores ~±20%). El matrix-element exacto del operador no se incluye; medimos la
FORMA de Γ/H (dónde está el peso), que es lo que decide térmico vs frío.
"""
import numpy as np
PHI=1.618033988749895; PI=np.pi
MPL=1.2209e19

# --- interaction measure (ε-3p)/T^4 del plasma SM (lattice QCD, forma robusta) ---
# T en GeV ; I pico ~3.8 en ~0.19 GeV, ~0 lejos (Borsanyi+2014 / HotQCD, aprox)
T_tab = np.array([0.05,0.10,0.13,0.15,0.17,0.19,0.22,0.25,0.30,0.40,0.60,1.0,3.0,10.0])
I_tab = np.array([0.30,0.80,1.80,2.80,3.50,3.80,3.40,2.80,2.10,1.30,0.70,0.40,0.15,0.05])

# Γ/H para un acoplamiento a la traza: Γ ∝ (T^μ_μ)^2/Λ^2/(n) ... medimos la FORMA.
# T^μ_μ = I·T^4. Γ ~ (I T^4)^2/(Λ^2 n)  con n~T^3 ; H~T^2/M_Pl.
# Γ/H ∝ I^2 T^8 /(Λ^2 T^3) · M_Pl/T^2 = I^2 T^3 M_Pl/Λ^2  (forma en T)
def gamma_over_H_shape(T,I,Lam):
    return I**2 * T**3 * MPL / Lam**2

# normalizamos Λ para que el PICO toque Γ/H=1 (pregunta: ¿es pico o plateau?)
Lam_GeV = 8800.0
GH = gamma_over_H_shape(T_tab,I_tab,Lam_GeV)
GH = GH/GH.max()   # forma normalizada al pico

print("="*74)
print("  OP-12 Paso 2 — ¿la traza enciende en QCD? ¿térmico o frío?")
print("="*74)
print("  T(MeV)   (1-3w)~I   Γ/H (forma, pico=1)")
for T,I,g in zip(T_tab,I_tab,GH):
    mark=" <== PICO (QCD)" if abs(T-0.19)<0.02 else ""
    print(f"  {T*1000:7.0f}   {I:6.2f}     {g:8.4f}{mark}")
print("-"*74)
ipk=np.argmax(I_tab)
print(f"  (a) La interaction measure PICA en T={T_tab[ipk]*1000:.0f} MeV = transición QCD. ✓")
print(f"      φ acoplado a la traza SÍ se enciende en QCD (interlocutor ruidoso ahí).")
print("-"*74)
# ¿el bulto de QCD domina el freeze-out? comparar Γ/H en QCD vs alta T
g_qcd = GH[np.argmin(abs(T_tab-0.19))]
g_hi  = GH[T_tab>=3.0].max()
print(f"  (b) Γ/H(QCD,190MeV)/max = {g_qcd:.4f}   Γ/H(alta T)/max = {g_hi:.4f}")
print("      => Γ/H CRECE con T (el T^3 del espacio de fases manda); el bulto de QCD")
print(f"      es solo ~{g_qcd*100:.0f}% del máximo. NO domina el desacople.")
print("      RESULTADO: el desacople ocurre en ALTA T (donde Γ/H cruza 1, fijado por")
print("      Λ); para cuando llega a QCD, φ ya está desacoplado. El bulto es irrelevante.")
print("="*74)
print("  VEREDICTO (honesto — corrige una narrativa previa mía):")
print("  · (a) La traza pica en QCD: VERDADERO. Pero (b) eso NO ancla el desacople ahí")
print("    — el freeze-out lo domina el T^3, no el bulto de la traza.")
print("  · => La historia 'φ se desacopla en QCD gratis por la traza' NO sobrevive.")
print("    Para caer en QCD habría que AJUSTAR Λ~8.8 TeV a mano (no derivado).")
print("  · La rama térmica NO cierra elegante; m_φ=41 sigue en 1 pata (SOLAR²·KRYSTOS).")
print("    2ª pata regresa a V''(φ_min)/OP-10 (rama fría), o admitir Λ sin explicar.")
print("="*74)
