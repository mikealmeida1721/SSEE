"""
OP-9 — Producción gravitacional de χ (materia oscura φ-DM): ¿cuánta sale?
=========================================================================
Mike (2026-05-31): m_φ ya NO es 5.6 eV (numerología). Vía portal Yukawa con
y=1 (natural, 0 params) sale m_φ = 46.29 meV. Paper 6 dice que χ se puebla por
PRODUCCIÓN GRAVITACIONAL (Parker 1968; de Sitter fluctuations) — sin ángulo de
mezcla, sin φ_i a mano. La abundancia depende SOLO de:
    (a) m_φ            ← 46.29 meV (Yukawa, fijo)
    (b) V(φ) inflación ← α-attractor con r = φ⁻¹⁰ (SSEE, fijo)
Ambos YA están fijos en SSEE. Entonces NO hay parámetro libre: el Ω_χ que salga
es una PREDICCIÓN, no un ajuste. La pregunta honesta: ¿da Ω_χ h² ≈ 0.12?

Mecanismo físico (escalar ligero, m_φ << H_inf):
  Durante inflación, el espacio en expansión "sacude el vacío" y el campo χ
  acumula fluctuaciones de de Sitter. Como m_φ << H_inf, NO alcanza equilibrio
  (eso tomaría ~10⁴⁷ e-folds); crece linealmente:
        ⟨χ²⟩ ≈ (H_inf/2π)² · N_efolds
  χ queda congelado hasta que H baja a ~m_φ, ahí oscila y se diluye como
  materia (a⁻³). De ahí sale la densidad reliquia hoy.
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import PHI, PI

# ── Constantes físicas ───────────────────────────────────────────────────
M_PL   = 2.435e18          # Planck REDUCIDA en GeV
A_S    = 2.1e-9            # amplitud perturbaciones escalares (Planck, medido)
N_EF   = 60.0             # e-folds de inflación (estándar)
G_STAR = 106.75          # grados de libertad relativistas (SM completo, T~TeV)
S0     = 2891.0          # densidad de entropía hoy [cm⁻³]
RHO_C_OVER_H2 = 1.054e-5 # ρ_crit/h² [GeV/cm³]

m_phi_eV = 46.2877        # Yukawa y=1
m_phi    = m_phi_eV * 1e-9  # → GeV

r_tensor = PHI**-10        # SSEE: razón tensor-escalar fija

print("=" * 74)
print(" OP-9 — Producción gravitacional de χ: ¿cuánto Ω_χ sale de m_φ + V?")
print("=" * 74)
print(f"  m_φ (Yukawa y=1)     = {m_phi_eV:.4f} meV = {m_phi:.4e} GeV")
print(f"  r = φ⁻¹⁰ (SSEE fijo) = {r_tensor:.6f}")

# ── Paso 1: escala de inflación H_inf desde r y A_s ──────────────────────
#   espectro tensorial  A_t = (2/π²)(H_inf/M_Pl)²   y   r = A_t/A_s
#   ⟹ H_inf = M_Pl · √( (π²/2)·r·A_s )
H_inf = M_PL * np.sqrt((PI**2/2) * r_tensor * A_S)
print("\n── Paso 1: escala de inflación H_inf (de r y A_s) ──")
print(f"  H_inf = M_Pl·√((π²/2)·r·A_s) = {H_inf:.4e} GeV")
print(f"  (m_φ/H_inf = {m_phi/H_inf:.2e} → χ es ULTRA-ligero vs inflación) ✓ régimen ligero")

# ── Paso 2: amplitud χ_i acumulada por fluctuaciones de de Sitter ─────────
chi_rms = (H_inf/(2*PI)) * np.sqrt(N_EF)
print("\n── Paso 2: amplitud del campo χ_i (acumulación de Sitter, N=60) ──")
print(f"  ⟨χ²⟩ = (H_inf/2π)²·N  →  χ_i = {chi_rms:.4e} GeV")

# ── Paso 3: temperatura de inicio de oscilación (cuando 3H = m_φ) ─────────
#   H_osc = m_φ/3 ;  H² = (π²/90)g_* T⁴/M_Pl²
T_osc = (10.0 * m_phi**2 * M_PL**2 / (PI**2 * G_STAR))**0.25
print("\n── Paso 3: T cuando χ empieza a oscilar (H = m_φ/3) ──")
print(f"  T_osc = {T_osc:.4e} GeV  ({T_osc/1e3:.2f} TeV)")

# ── Paso 4: densidad reliquia hoy (entropía conservada ρ/s) ──────────────
rho_osc = 0.5 * m_phi**2 * chi_rms**2                      # GeV⁴
s_osc   = (2*PI**2/45) * G_STAR * T_osc**3                 # GeV³
rho_over_s = rho_osc / s_osc                               # GeV
Omega_h2 = rho_over_s * S0 / RHO_C_OVER_H2
print("\n── Paso 4: densidad reliquia hoy ──")
print(f"  ρ_osc = ½m_φ²χ_i²        = {rho_osc:.4e} GeV⁴")
print(f"  s_osc = (2π²/45)g_*T³    = {s_osc:.4e} GeV³")
print(f"  ρ_χ/s (conservado)       = {rho_over_s:.4e} GeV")
print(f"  ➤ Ω_χ h² (PREDICHO)      = {Omega_h2:.4e}")
print(f"     objetivo (Planck DM)  = 0.12")
print(f"     ratio predicho/objetivo = {Omega_h2/0.12:.3e}×")

# ── Paso 5: diagnóstico — ¿qué tendría que cambiar para dar 0.12? ─────────
#   Ω h² ∝ H_inf² · √m_φ   (derivado del scaling)
print("\n── Paso 5: diagnóstico de escala (Ω h² ∝ H_inf²·√m_φ) ──")
factor = Omega_h2 / 0.12
H_inf_needed = H_inf / np.sqrt(factor)
r_needed = (H_inf_needed/M_PL)**2 / ((PI**2/2)*A_S)
m_needed = m_phi / factor**2
print(f"  Para Ω h²=0.12 con m_φ=46.29 meV fijo:")
print(f"     H_inf debería ser {H_inf_needed:.3e} GeV  (÷{np.sqrt(factor):.1f})")
print(f"     ⟹ r debería ser {r_needed:.3e}   (SSEE da r=φ⁻¹⁰={r_tensor:.3e})")
print(f"  Para Ω h²=0.12 con r=φ⁻¹⁰ fijo:")
print(f"     m_φ debería ser {m_needed*1e9:.3e} eV  (ultraligero, no 46 meV)")

print("\n" + "=" * 74)
print(" LECTURA")
print("=" * 74)
print(f"""
  Predicción SIN parámetros libres (m_φ=46.29 meV + r=φ⁻¹⁰):
     Ω_χ h² = {Omega_h2:.2e}   vs   0.12 objetivo.

  Si |Ω−0.12| es chico → producción gravitacional CIERRA la DM sin ajustar nada.
  Si Ω >> 0.12 (sobreproduce) → la combinación (46 meV ligero + inflación de alta
     escala r=φ⁻¹⁰) genera demasiada DM por fluctuaciones de de Sitter.
     Eso NO es muerte: es el clásico choque "escalar ligero + inflación alta".
     Señala una pieza faltante (dilución por entropía, o m_φ/H_inf distinto),
     NO una rendición.
""")
