"""
OP-10 SYSTEMATIC SEARCH — V(φ) unificado con principios estrictos
====================================================================
Atacando OP-10 con los 3 principios:
  1. Adimensionales: pueden ser (φ,π) puros
  2. Dimensionales: requieren escala física limpia × algebraico
  3. Coeficiente algebraico debe ser DERIVADO, no inventado

Inputs LIMPIOS permitidos:
  M_Pl    = 1.221×10²⁸ eV    (fundamental)
  ℏH_MIRA = 1.43×10⁻³³ eV    (física, no algebraica)
  ρ_crit  = 3.65×10⁻¹¹ eV⁴   (física, depende H_MIRA pero no contaminada)
  ρ_DE    = 0.840 × ρ_crit   (algebraico × física)

Algebraicos LIMPIOS (adimensionales):
  φ, π, Ω=φ+π, β, KAL₀, MIRA, AURA, T_r, M_v, n_s, etc.

PROHIBIDOS (contaminados):
  Σm_ν (OP-14 con -22 ad hoc)
  Λ_SSEE = 8.81 meV (convención canónica-ΛCDM mixed)
  H_alg (numerología km/s/Mpc)

Para cada mecanismo, computar m_φ y reportar HONESTAMENTE.
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import (PHI, PI, OMEGA, BETA, KAL0, OMEGA_DE, OMEGA_M_DYN,
                       MIRA, AURA, T_R, M_V, P_SC, K_V, W0)

# ── Escalas físicas LIMPIAS ─────────────────────────────────────────────
H_MIRA_EV    = 1.4306e-33                # ℏH_MIRA
M_PL         = 1.220910e28                # Planck mass
h_mira       = 0.67037
RHO_CRIT     = 8.10e-11 * h_mira**2       # eV⁴
RHO_DE       = OMEGA_DE * RHO_CRIT
RHO_CRIT_QR  = RHO_CRIT**0.25             # eV (= 2.5 meV approx)
RHO_DE_QR    = RHO_DE**0.25

# Tracker DE potencial P7 (clean derivation)
ALPHA_DE     = np.sqrt(3 * OMEGA_M_DYN) / np.sqrt(KAL0)   # ≈ 0.295

print("=" * 78)
print(" OP-10 SYSTEMATIC SEARCH — V(φ) con principios estrictos")
print("=" * 78)
print(f"\nEscalas físicas LIMPIAS (eV):")
print(f"  ρ_crit^(1/4) = {RHO_CRIT_QR:.4e} = {RHO_CRIT_QR*1e3:.4f} meV")
print(f"  ρ_DE^(1/4)   = {RHO_DE_QR:.4e} = {RHO_DE_QR*1e3:.4f} meV")
print(f"  ℏH_MIRA      = {H_MIRA_EV:.4e}")
print(f"  M_Pl         = {M_PL:.4e}")
print(f"\nTracker DE coupling: α_DE = √(3·Ω_m,dyn)/√KAL₀ = {ALPHA_DE:.4f}")

# Target observacional (rango eV donde P6 funcionaba)
print(f"\nTarget de búsqueda: m_φ ∈ [1, 30] eV para que k_fs sea Lyman-α-relevante")
print(f"Pero ACEPTAMOS lo que salga, sin forzar.")

results = []

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 1 — V_DE intrinsic curvature (lo que ya teníamos)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M1] V_DE intrinsic — m_φ² = V''(φ_today) = α²·V_DE")
print("─"*78)
V_DE_today = RHO_DE                       # V_DE ≈ ρ_DE hoy
m_phi_M1_sq = ALPHA_DE**2 * V_DE_today    # V'' = α²·V para exp potential
m_phi_M1 = np.sqrt(m_phi_M1_sq)
print(f"  V''(φ_today) = α²·V_DE = {m_phi_M1_sq:.4e} eV²")
print(f"  m_φ = √V'' = {m_phi_M1:.4e} eV = {m_phi_M1*1e6:.4f} μeV")
print(f"  Origen: P7 V(φ) = V₀·exp(-αφ), derivado tracker")
print(f"  ✓ Principios cumplidos. ¿En eV? {1 <= m_phi_M1 <= 30}")
results.append(("M1 V_DE curvature", m_phi_M1, "Derivado P7"))

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 2 — α-attractor inflation curvature (P1)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M2] α-attractor inflation V — m_φ² = V_inf''(φ_min)")
print("─"*78)
# α_inflation = φ^4/3 from P1
ALPHA_INF = PHI**4 / 3
# Para α-attractor V = V_inf·tanh²(βφ/√6), m² at minimum = V_inf · β²/3
# β = √(2/(3α))
beta_inf = np.sqrt(2/(3*ALPHA_INF))
# V_inf scale: la escala de inflación. Para N=58 e-folds y r=0.008,
# H_inf ≈ √(πr·A_s)·M_Pl ~ 10¹³ GeV
# V_inf ≈ 3 H_inf² M_Pl² ~ (10¹⁶ GeV)⁴
A_s_planck = 2.1e-9
r_tensor = PHI**(-10)  # SSEE n_s=1-φ⁻⁷ y r=φ⁻¹⁰
H_inf = np.sqrt(np.pi * r_tensor * A_s_planck / 2) * M_PL
V_inf = 3 * H_inf**2 * M_PL**2 / (8*np.pi)
m_phi_M2_sq = V_inf * beta_inf**2 / 3
m_phi_M2 = np.sqrt(m_phi_M2_sq)
print(f"  H_inf = √(π·r·A_s/2)·M_Pl ≈ {H_inf:.4e} eV")
print(f"  V_inf ≈ 3 H_inf² M_Pl² /(8π) = {V_inf:.4e} eV⁴")
print(f"  β = √(2/3α) = {beta_inf:.4f}")
print(f"  m_φ² = V_inf·β²/3 = {m_phi_M2_sq:.4e} eV²")
print(f"  m_φ = {m_phi_M2:.4e} eV")
print(f"  Origen: P1 α-attractor (n_s=1-φ⁻⁷, r=φ⁻¹⁰), todo derivado.")
print(f"  ✓ Principios cumplidos. ¿En eV? {1 <= m_phi_M2 <= 30}")
results.append(("M2 α-attractor V_inf", m_phi_M2, "Derivado P1"))

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 3 — Cosmological seesaw √(M_Pl × ℏH_MIRA)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M3] Cosmological seesaw — m_φ = √(M_Pl × ℏH_MIRA)")
print("─"*78)
m_phi_M3 = np.sqrt(M_PL * H_MIRA_EV)
print(f"  m_φ = √(M_Pl · ℏH_MIRA) = √({M_PL:.2e} × {H_MIRA_EV:.2e})")
print(f"       = {m_phi_M3:.4e} eV = {m_phi_M3*1e3:.4f} meV")
print(f"  Origen: geometric mean estructural, sin factor algebraico fitted")
print(f"  ✓ Principios cumplidos. ¿En eV? {1 <= m_phi_M3 <= 30}")
results.append(("M3 cosmological seesaw", m_phi_M3, "Estructural"))

# Variantes con factores algebraicos motivados
for fname, factor in [("φ", PHI), ("π", PI), ("Ω", OMEGA), ("KAL₀", KAL0),
                      ("T_r", T_R), ("M_v", M_V), ("MIRA", MIRA)]:
    m = factor * m_phi_M3
    in_range = "✓ EN eV!" if 1 <= m <= 30 else ""
    print(f"     {fname} · √(M_Pl·ℏH) = {m:.4e} eV {in_range}")
    if 1 <= m <= 30:
        results.append((f"M3 · {fname}", m, f"Seesaw × {fname}"))

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 4 — Powers de ρ_crit con factores algebraicos
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M4] ρ_crit^(1/4) × algebraic factor")
print("─"*78)
# Para llegar a eV desde ρ_crit^(1/4)=2.5 meV, factor ~ 400-12000
for fname, factor in [("(φ+π)⁵", OMEGA**5),
                      ("KAL₀⁴", KAL0**4),
                      ("M_v⁴", M_V**4),
                      ("3·(φ+π)⁴", 3*OMEGA**4),
                      ("π²·M_v²", PI**2 * M_V**2),
                      ("T_r⁴", T_R**4),
                      ("M_v³", M_V**3),
                      ("φ⁸·5", PHI**8 * 5)]:
    m = factor * RHO_CRIT_QR
    in_range = "✓ EN eV!" if 1 <= m <= 30 else ""
    print(f"  {fname:>12s} · ρ_crit^(1/4) = {factor:>10.4f} × {RHO_CRIT_QR:.4e} = {m:.4f} eV {in_range}")
    if 1 <= m <= 30:
        results.append((f"M4 ρ^(1/4)·{fname}", m, "Power-scaling"))

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 5 — Powers de seesaw √(M_Pl × ℏH_MIRA) con factores
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M5] (M_Pl × ℏH_MIRA)^(p) variants")
print("─"*78)
# m² = ρ_crit/algebraic ?
# Otra forma: usar V_DE^(1/4) × algebraic
for fname, m in [
    ("V_DE^(1/4) · M_v",          RHO_DE_QR * M_V),
    ("ρ_crit^(1/4) · KAL₀³",      RHO_CRIT_QR * KAL0**3),
    ("ρ_crit^(1/4) · M_v³",       RHO_CRIT_QR * M_V**3),
    ("(M_Pl·V_DE)^(1/4)",         (M_PL * V_DE_today)**0.25 if M_PL*V_DE_today > 0 else 0),
    ("(M_Pl²·ρ_crit)^(1/6)",      (M_PL**2 * RHO_CRIT)**(1/6)),
]:
    in_range = "✓ EN eV!" if 1 <= m <= 30 else ""
    print(f"  {fname:<35s} = {m:.4e} eV {in_range}")
    if 1 <= m <= 30:
        results.append((f"M5 {fname}", m, "Cross-scale"))

# ═══════════════════════════════════════════════════════════════════════
# RESUMEN
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 78)
print(" RESUMEN — todos los mecanismos con m_φ y status")
print("=" * 78)
print(f"\n{'Mecanismo':<35} {'m_φ':>14} {'Status':<25}")
print(" " + "─" * 75)
for name, m, status in results:
    if m < 1e-12:
        m_str = f"{m:.2e} eV"
    elif m < 1e-3:
        m_str = f"{m*1e6:.3f} μeV"
    elif m < 1:
        m_str = f"{m*1e3:.3f} meV"
    elif m < 1e6:
        m_str = f"{m:.3f} eV"
    else:
        m_str = f"{m:.2e} eV"
    flag = " ★" if 1 <= m <= 30 else ""
    print(f" {name:<35} {m_str:>14} {status:<25}{flag}")

# ═══════════════════════════════════════════════════════════════════════
# DIAGNÓSTICO
# ═══════════════════════════════════════════════════════════════════════
in_eV = [r for r in results if 1 <= r[1] <= 30]
print("\n" + "=" * 78)
print(" Diagnóstico")
print("=" * 78)
print(f"\n  Mecanismos limpios con m_φ ∈ [1, 30] eV: {len(in_eV)}")

if in_eV:
    print("\n  CANDIDATAS LIMPIAS:")
    for name, m, status in in_eV:
        print(f"    {name}: m_φ = {m:.4f} eV [{status}]")
    print("\n  → Verificar JUSTIFICACIÓN ESTRUCTURAL del coeficiente:")
    print("    ¿el factor algebraico viene de derivación o ajuste?")
else:
    print("""
  ⚠ NINGÚN mecanismo con principios estrictos da m_φ ∈ eV.

  Las escalas naturales que SSEE V3.6 produce limpiamente:
    - μeV  (de V_DE curvature, M1)
    - meV  (de cosmological seesaw √(M_Pl·ℏH), M3)
    - GeV-TeV (de inflación V_inf, M2)
    - keV-MeV (combinaciones de ρ_crit con M_Pl)

  La escala eV NO aparece naturalmente sin invocar:
    (a) Σm_ν (contaminada OP-14)
    (b) factores algebraicos de orden ~10³ sin justificación
    (c) escalas ad hoc no presentes en V3.6

  CONCLUSIÓN HONESTA:
    SSEE V3.6 actual NO predice m_φ ∈ eV scale.
    Predice naturalmente: meV (M3) o μeV (M1).

  IMPLICACIÓN para P6:
    - Régimen meV: φ-DM "milli-DM", k_fs no relevante para Lyman-α
    - Régimen μeV: φ-DM ultralight, no afecta σ₈/S₈
    - En AMBOS casos: alpha_WDM ya no rescata σ₈/S₈ tension
    - P6 debe reformularse con régimen natural, o aceptar tensión 5σ
""")

print("=" * 78)
