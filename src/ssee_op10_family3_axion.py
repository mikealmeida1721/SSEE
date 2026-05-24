"""
OP-10 Fase 2b — Familia 3: axion-like + exponencial DE
========================================================

Potencial:
    V(φ) = V₀ e^(-λφ) + M_UV⁴ [1 - cos(φ/f)]

Donde:
    V₀     = 0.840 ρ_crit          (bloqueado Paper 7)
    λ      = √(3·Ω_m,dyn) ≈ 0.693  (bloqueado Paper 7)
    M_UV   ≈ 8.81 meV              (bloqueado Paper 10)
    f      = ÚNICO free param

ANÁLISIS:

(1) Mínimo local del axion: cos(φ_min/f) = 1, i.e. φ_min = 2πn·f
    En el mínimo: V_min = V₀ e^(-λ·2πn·f)
    Curvatura: V''(φ_min) = λ² V₀ e^(-λφ_min) + M_UV⁴/f²

(2) Para axion estándar, el término M_UV⁴/f² domina la curvatura:
    m_φ² ≈ M_UV⁴/f² (axion-like mass)

(3) Para m_φ = 5.60 eV:
    f = M_UV²/m_φ ≈ (8.81 meV)² / 5.60 eV ≈ 13.9 μeV

(4) PROBLEMA CRÍTICO — incompatibilidad DM-oscilante vs DE-dinámica:
    Campo se "engancha" en mínimo del coseno cuando m_eff² > H²
    m_eff²/H² ≈ 31 eV² / (10⁻⁶⁶ eV²) = 10⁶⁷
    → engancha INSTANTÁNEAMENTE → DEJA DE RODAR → NO HAY DE DINÁMICA
    → w_a = 0, w₀ = -1 (LambdaCDM efectivo)
    → NO reproduce w₀ = -0.840, w_a = -0.670
"""
import numpy as np

# Constantes
PHI          = (1 + np.sqrt(5)) / 2
PI           = np.pi
KAL          = (PHI + 3*PI) / 2
M_V          = 3*(PHI + PI)

RHO_CRIT_eV4 = 8.10e-11
V0           = 0.840 * RHO_CRIT_eV4
LAMBDA_LOCK  = np.sqrt(3 * 0.160)
M_UV         = 8.81e-3            # eV (valor canónico Paper 10)
M4_UV        = M_UV**4

H0_eV        = 1.45e-33
M_PHI_TARGET = 5.60               # eV
W0_TARGET    = -0.840

print("=" * 76)
print("OP-10 FASE 2b — Familia 3: axion-like + exponencial")
print("=" * 76)
print(f"""
Constantes bloqueadas:
  V₀     = {V0:.3e} eV⁴
  λ      = {LAMBDA_LOCK:.4f}
  M_UV   = {M_UV*1e3:.2f} meV
  M_UV⁴  = {M4_UV:.3e} eV⁴
""")

# ============================================================
# Paso 1 — Calcular f desde m_φ = √(M_UV⁴/f²)
# ============================================================
f_required = M_UV**2 / M_PHI_TARGET
print(f"[1] CONSTRAINT 1 (axion mass):")
print(f"    m_φ² ≈ M_UV⁴/f²  →  f = M_UV²/m_φ = {f_required:.3e} eV  ≈ {f_required*1e6:.2f} μeV")

# ¿f tiene forma algebraica?
print(f"\n[1.1] ¿f algebraico?")
candidates = {
    "M_UV/(2π)":      M_UV/(2*PI),
    "M_UV/(2π·KAL)":  M_UV/(2*PI*KAL),
    "M_UV/(M_v·π)":   M_UV/(M_V*PI),
    "M_UV/M_v²":      M_UV/M_V**2,
    "M_UV²/m_φ":      M_UV**2/M_PHI_TARGET,
    "M_UV/φ⁵":        M_UV/PHI**5,
    "M_UV/(KAL·π²)":  M_UV/(KAL*PI**2),
}
print(f"    {'forma':<22} {'valor [eV]':<14} {'ratio vs f_req':<12}")
for name, val in candidates.items():
    print(f"    {name:<22} {val:<14.3e} {val/f_required:<12.4f}")

# ============================================================
# Paso 2 — Test de "enganche": ¿el campo se engancha en mínimo?
# ============================================================
print(f"\n[2] TEST DE ENGANCHE (¿axion oscila o rueda con DE?)")
print(f"    Campo se engancha cuando m_eff² > H²")
m_eff2 = M4_UV / f_required**2
H02 = H0_eV**2
print(f"    m_eff²              = M_UV⁴/f²       = {m_eff2:.3e} eV²")
print(f"    H₀²                 = {H02:.3e} eV²")
print(f"    Ratio m_eff²/H²     = {m_eff2/H02:.3e}")
print(f"    → m_eff/H           ≈ {np.sqrt(m_eff2/H02):.3e}")

if m_eff2 > 100 * H02:
    print(f"\n    ⚠ ENGANCHE INSTANTÁNEO: m_eff/H ≈ 10^{int(np.log10(np.sqrt(m_eff2/H02)))}")
    print(f"      El axion oscila MUY rápido respecto a H — se comporta como DM coherente.")
    print(f"      Pero entonces el rolling DE viene SOLO del término V₀ e^(-λφ).")

# ============================================================
# Paso 3 — Si axion enganchado, ¿qué DE queda?
# ============================================================
print(f"\n[3] DE EFECTIVA cuando axion engancha:")
print(f"    V_DE residual = V₀ e^(-λφ_now)")
print(f"    Pero ahora φ_now está FIJO en el mínimo del coseno más cercano.")
print(f"    → V_DE = constante  →  w_DE = -1 exacto (no -0.840)")
print(f"    → w_a = 0           (no -0.670)")
print(f"    → DE es ΛCDM efectivo, NO la dinámica SSEE")

# ============================================================
# Paso 4 — ¿Hay rescate? Tres sub-casos
# ============================================================
print(f"\n[4] ¿RESCATES POSIBLES?")

# Sub-caso A: f tan grande que el axion no oscila
f_no_eng = M_UV**2 / H0_eV  # f donde m_eff = H₀
print(f"""
    A) f tan grande que m_eff < H₀ (axion frozen, no DM):
       Requiere f > {f_no_eng:.3e} eV ≈ {f_no_eng/M_UV:.3e} × M_UV
       Con ese f, m_φ = M_UV²/f = {M_UV**2/f_no_eng:.3e} eV
       → m_φ ≈ H₀, contradice φ-DM con m_φ = 5.60 eV
       ❌ Descartado.""")

# Sub-caso B: usar el axion solo para DM y otro mecanismo para DE
print(f"""
    B) Axion oscilante para DM + V₀ e^(-λφ_DE) en OTRO campo:
       Esto es two-field model (no Familia 3 single-field).
       Equivalente a Paper 6 two-sector actual.
       OP-9 NO se resuelve: m_φ axion-like sigue siendo "input" vía f.
       ❌ No avanza sobre estado actual.""")

# Sub-caso C: que axion sea DE oscilante (controvertido)
print(f"""
    C) Axion como DE oscilante:
       Requeriría m_φ ≲ H₀ ≈ 10⁻³³ eV (period mayor que edad universo).
       Pero queremos m_φ = 5.60 eV, NO m_φ < H₀.
       ❌ Descartado.""")

# ============================================================
# Veredicto
# ============================================================
print("=" * 76)
print("VEREDICTO FAMILIA 3")
print("=" * 76)
print(f"""
  ✗ DESCARTADA por mismo patrón estructural que Familia 2.

  CAUSA: dos escalas separadas (V₀^(1/4) ≈ {V0**0.25:.2e} eV vs m_φ = 5.60 eV)
  no pueden coexistir en un solo campo:

  - Si axion engancha (m_eff >> H): φ se fija → V_DE = const → w=-1 (no -0.84)
  - Si axion no engancha (m_eff < H): m_φ ≈ H₀ (no 5.60 eV)

  El problema es topológico: el potencial periódico [1-cos] tiene mínimos
  discretos donde el campo se queda atrapado. Para que SSEE-DE funcione
  (rolling exponencial λ≈0.693), el campo debe rodar libremente — pero
  el término axion lo atrapa en CUALQUIER f > μeV.

  ESTADO TRAS FASES 2 + 2b:
    F1 ❌ (contamina plateau)
    F2 ❌ (jerarquía DE/DM, 33 órdenes en f)
    F3 ❌ (enganche topológico vs DE dinámica)
    F4 ❌ (predicción off por 10³³)

  TODAS LAS FAMILIAS SINGLE-FIELD CANÓNICAS AGOTADAS.

  PRÓXIMAS RUTAS (más exóticas, requieren más infraestructura):
    Fase 2c — K-essence no-canónica: K(X)=X/KAL + X²/M⁴ (Paper 10 lock).
              El factor K_X = 1/KAL + 2X/M⁴ renormaliza la "masa efectiva"
              del campo, posiblemente rescatando la jerarquía.
    Fase 2d — Two-field see-saw: m_φ² = m_h·m_l donde m_h, m_l vienen
              de dos escalas SSEE (Λ_UV y Λ_DE). Genera jerarquía
              naturalmente.
    Fase 2e — Conformal/disformal coupling: g_μν → A(φ)g + B(φ)∂φ∂φ
              genera masas efectivas dinámicas.

  RECOMENDACIÓN: Fase 2c (K-essence) primero porque Paper 10 YA tiene K(X)
  bloqueada. Es la ruta más "natural" en términos de la estructura SSEE
  ya construida.
""")
