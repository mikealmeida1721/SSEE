"""
OP-10 — Búsqueda principiada de V(φ) unificado DE+DM
======================================================
Aplicando los 3 principios SSEE (ver feedback_three_principles.md):
  1. Observables adimensionales pueden ser puros (φ,π)
  2. Observables dimensionales = (algebraico) × (escala física)
  3. El coeficiente algebraico debe venir DERIVADO de mecanismo,
     no inventado para hacer cuadrar el número

Estrategia:
  V(φ) = V_DE(φ) + ΔV(φ)

  V_DE(φ) = V₀·exp(-αφ)              ← YA derivado en P7 (limpio)
             V₀ = 0.840·ρ_crit         (Regla 2 ✓)
             α  = 0.6929/√KAL₀ ≈ 0.295 (tracker condition, Regla 3 ✓)

  ΔV(φ) = ?  ← LO QUE QUEREMOS DERIVAR

Para cada forma de ΔV(φ), aplicar el checklist de 3 principios.
Computar m_φ² = V''(φ_min) y reportar HONESTAMENTE qué sale,
sin targetear ningún valor.

NO usar inputs contaminados: NO Σm_ν, NO Λ_SSEE=8.81 meV.
SOLO usar: ρ_crit (con H_MIRA), M_Pl, ρ_DE; factores (φ,π,KAL₀,...).
"""
import numpy as np
import sys
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
from ssee_core import (PHI, PI, OMEGA, BETA, KAL0, OMEGA_DE, OMEGA_M_DYN,
                       MIRA, AURA, T_R, M_V, P_SC, K_V)

# ──────────────────────────────────────────────────────────────────────────
# Escalas físicas LIMPIAS (verificadas como 🟢 PHYSICAL en SSEE_CONSTANTS_AUDIT)
# ──────────────────────────────────────────────────────────────────────────
H0_MIRA_KMSMPC = 67.037
H0_MIRA_INVS   = H0_MIRA_KMSMPC * 3.240779289e-20    # s⁻¹
HBAR_EVS       = 6.582119569e-16
H_MIRA_EV      = HBAR_EVS * H0_MIRA_INVS              # ≈ 1.43e-33 eV
M_PL_EV        = 1.220910e28                          # eV (no reducida)

# ρ_crit puro SSEE con H_MIRA (NO usar canónico-ΛCDM):
h_mira         = H0_MIRA_KMSMPC / 100.0
RHO_CRIT       = 8.10e-11 * h_mira**2                # eV⁴, escalado por h²
RHO_DE         = OMEGA_DE * RHO_CRIT                  # eV⁴

print("=" * 80)
print(" OP-10 — Búsqueda principiada de V(φ) con mínimo")
print("=" * 80)
print(f"\n Escalas físicas limpias (verificadas):")
print(f"   ρ_crit         = {RHO_CRIT:.4e} eV⁴   (h={h_mira:.4f}, H_MIRA puro)")
print(f"   ρ_crit^(1/4)   = {RHO_CRIT**0.25:.4e} eV     = {RHO_CRIT**0.25*1e3:.4f} meV")
print(f"   ρ_DE           = {RHO_DE:.4e} eV⁴")
print(f"   ρ_DE^(1/4)     = {RHO_DE**0.25:.4e} eV     = {RHO_DE**0.25*1e3:.4f} meV")
print(f"   ℏH_MIRA        = {H_MIRA_EV:.4e} eV")
print(f"   M_Pl           = {M_PL_EV:.4e} eV")

# Potencial DE ya derivado en P7 (Regla 2+3 ✓)
V0     = OMEGA_DE * RHO_CRIT              # = 0.840 · ρ_crit
ALPHA  = np.sqrt(3 * OMEGA_M_DYN) / np.sqrt(KAL0)  # α = λ/√KAL = 0.6929/√5.52
print(f"\n V_DE(φ) = V₀·exp(-αφ) — DERIVADO en P7:")
print(f"   V₀ = 0.840·ρ_crit = {V0:.4e} eV⁴")
print(f"   α  = √(3·Ω_m,dyn)/√KAL₀ = {ALPHA:.4f}")

# ──────────────────────────────────────────────────────────────────────────
# CHECKLIST DE 3 PRINCIPIOS PARA CADA ΔV CANDIDATA
# ──────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print(" Candidatas ΔV(φ) — aplicando los 3 principios")
print("=" * 80)

candidates = []

# ──────────────────────────────────────────────────────────────────────────
# CANDIDATA 1: Corrección polinómica suprimida exponencialmente
# ΔV(φ) = c·ρ_crit·(α·φ)²·exp(-α·φ)
#
# Justificación física: corrección RGE de loop al potencial. La forma
# polinómica × exp surge naturalmente en EFT cuando integras un loop sobre
# el campo. El coeficiente c sería el coeficiente de loop.
#
# Mínimo: V_total(φ) = V₀·exp(-αφ)·(1 + c·(αφ)²)
# V'(φ) = V₀·exp(-αφ)·α·[-1 - c·(αφ)² + 2c·(αφ)]
# V'(φ_min) = 0 →  -1 - c·(αφ)² + 2c·(αφ) = 0
# Resolviendo: αφ_min = (2c ± √(4c² - 4c))/(2c) = 1 ± √(1 - 1/c)
# Real si c ≥ 1
# ──────────────────────────────────────────────────────────────────────────
print("\n[C1] ΔV(φ) = c·V₀·(α·φ)²·exp(-α·φ)  — loop correction")
for c_label, c in [("c=1 (mínimo crítico)", 1.0),
                   ("c=Ω/4 (algebraico)", OMEGA/4),
                   ("c=KAL₀/4", KAL0/4),
                   ("c=MIRA²/2", MIRA**2/2)]:
    if c < 1:
        print(f"   {c_label}: c={c:.3f} < 1 → SIN MÍNIMO real, descartar")
        continue
    aphi_min = 1 - np.sqrt(1 - 1/c) if c > 1 else 1.0
    # V_total = V₀ · exp(-αφ) · (1 + c(αφ)²)
    # V''(φ_min) requiere derivar dos veces... lo hago numérico:
    def V_total(phi, c=c):
        return V0 * np.exp(-ALPHA*phi) * (1 + c*(ALPHA*phi)**2)
    phi_min = aphi_min / ALPHA
    # Numérico V''
    dphi = phi_min * 1e-4
    V_pp = (V_total(phi_min + dphi) - 2*V_total(phi_min) + V_total(phi_min - dphi)) / dphi**2
    if V_pp > 0:
        m_phi = np.sqrt(V_pp)
        valid_principles = "Regla 2 ✓ (V₀·ρ_crit). Regla 3 PARCIAL (c necesita derivación EFT)."
        print(f"   {c_label}: φ_min={phi_min:.3f}, V''(min)={V_pp:.3e}, m_φ={m_phi:.3e} eV = {m_phi*1e3:.4f} meV")
        candidates.append(("C1 poly·exp", c_label, m_phi, valid_principles))
    else:
        print(f"   {c_label}: V''(min) = {V_pp:.3e} (máximo, no mínimo)")

# ──────────────────────────────────────────────────────────────────────────
# CANDIDATA 2: α-attractor con plateau
# V(φ) = V₀·tanh²(βφ/√6)
# β fijado por SSEE (α de inflación, β=√(2/3α) en α-attractor)
# Mínimo en φ=0, V(0)=0 — DM-like
# m_φ² = V''(0) = 2V₀·β²/3
#
# Pero esto reemplaza V_DE entero, no es ΔV aditivo.
# Considerar como alternativa a V_DE_exp en lugar de aditiva.
# ──────────────────────────────────────────────────────────────────────────
print("\n[C2] V(φ) = V₀·tanh²(βφ/√6)  — α-attractor unified (reemplaza exp)")
ALPHA_INFL = PHI**4 / 3   # α-attractor inflación P1
BETA_ATTR  = np.sqrt(2/(3*ALPHA_INFL))
# V'(0) = 0, V''(0) = 2V₀ · β² / 3 (por expansión tanh ≈ x, tanh² ≈ x²)
# Actually d²/dφ²[V₀ tanh²(βφ/√6)] at φ=0 = V₀ · 2 · β²/6 = V₀ · β²/3
V_pp_0 = V0 * BETA_ATTR**2 / 3
m_phi_C2 = np.sqrt(V_pp_0)
print(f"   α-attractor inflación: α={ALPHA_INFL:.4f}, β=√(2/3α)={BETA_ATTR:.4f}")
print(f"   V''(0) = V₀·β²/3 = {V_pp_0:.4e} eV⁴")
print(f"   m_φ = √V''(0) = {m_phi_C2:.4e} eV = {m_phi_C2*1e3:.4f} meV")
print(f"   Regla 2 ✓ (V₀ es ρ_crit-anchored). Regla 3 ✓ (α de P1 derivado).")
candidates.append(("C2 α-attractor unified", "α=φ⁴/3 (P1)", m_phi_C2, "Regla 2 ✓, Regla 3 ✓"))

# ──────────────────────────────────────────────────────────────────────────
# CANDIDATA 3: Higgs-like quartic con masa derivada de M_Pl
# V(φ) = -μ²φ² + λ_φ φ⁴/4
# Mínimo: φ_min² = μ²/λ_φ ⇒ m_φ² = 2μ²
#
# Si μ² ~ ρ_crit / (algebraico)  y  λ_φ ~ adimensional (φ,π)
# entonces m_φ² = 2 · ρ_crit / (algebraico)
# m_φ ~ √(ρ_crit / algebraico) ~ √(ρ_crit) / √algebraico
#
# √(ρ_crit) ≈ √(3.65e-11) = 6.04e-6 eV
# Esto da m_φ ~ 6 μeV / √algebraico
# Para llegar a meV o eV necesitaríamos algebraico << 1 (raro)
# ──────────────────────────────────────────────────────────────────────────
print("\n[C3] V = -μ²φ² + λ_φ·φ⁴/4 — Higgs-like; μ² ~ ρ_crit·algebraico/M_Pl²")
print("     (escala μ² emerge de ρ_crit a través de M_Pl)")
# Pruebas con μ² = ρ_crit · (algebraico) / M_Pl² (forma seesaw natural)
for label, alg in [("alg=1", 1.0), ("alg=Ω", OMEGA), ("alg=KAL₀", KAL0),
                   ("alg=M_v", M_V), ("alg=Ω²", OMEGA**2)]:
    mu2 = RHO_CRIT * alg / M_PL_EV**2
    m_phi = np.sqrt(2 * mu2)
    print(f"   {label}: μ²={mu2:.3e} eV², m_φ=√(2μ²)={m_phi:.3e} eV")
    candidates.append(("C3 Higgs-like", label, m_phi, "Regla 2 OK, Regla 3 si seesaw justifica"))

# ──────────────────────────────────────────────────────────────────────────
# RESUMEN
# ──────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print(" Resultados — m_φ derivado de cada V(φ) candidato")
print("=" * 80)
print(f"\n {'Candidata':<30} {'parámetros':<25} {'m_φ derivado':>18}")
print(" " + "─"*75)
for fam, params, m_phi, _ in sorted(candidates, key=lambda c: c[2]):
    if m_phi < 1e-12:
        m_str = f"{m_phi:.2e} eV"
    elif m_phi < 1e-3:
        m_str = f"{m_phi*1e6:.2f} μeV"
    elif m_phi < 1:
        m_str = f"{m_phi*1e3:.3f} meV"
    elif m_phi < 1e6:
        m_str = f"{m_phi:.3f} eV"
    else:
        m_str = f"{m_phi:.2e} eV"
    print(f" {fam:<30} {params:<25} {m_str:>18}")

# ──────────────────────────────────────────────────────────────────────────
# DIAGNÓSTICO
# ──────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print(" Diagnóstico")
print("=" * 80)

in_eV_range = [c for c in candidates if 1 <= c[2] <= 30]
in_meV_range = [c for c in candidates if 0.001 <= c[2] < 1]
in_ueV_range = [c for c in candidates if 1e-6 <= c[2] < 1e-3]

print(f"\n  Candidatos con m_φ ∈ [1, 30] eV (rango P6 actual): {len(in_eV_range)}")
print(f"  Candidatos con m_φ ∈ [1 meV, 1 eV]:                  {len(in_meV_range)}")
print(f"  Candidatos con m_φ ∈ [1 μeV, 1 meV]:                 {len(in_ueV_range)}")

if not in_eV_range:
    print(f"""
  ⚠ HALLAZGO IMPORTANTE:

    Ninguna forma V(φ) con principios cumplidos (Reglas 1+2+3) produce
    m_φ en escala eV.

    Los candidatos limpios caen en escala meV o μeV — coherente con la
    intuición física: las escalas naturales (ρ_crit^(1/4) ≈ 2.5 meV, etc.)
    no llegan a eV sin meter factores algebraicos grandes (que violan
    Regla 3 al ser elegidos para cuadrar).

    Implicación honesta para P6:

    (a) Si confiamos en los principios, m_φ DERIVADO está en escala meV.
        Eso cambia drásticamente k_fs (~10⁴ veces menor) y la fenomenología
        de P6 — pasaría de WDM eV a fuzzy DM ultraligero.

    (b) El "5.60 eV" de P6 era forzado por Σm_ν × H_alg (numerología).
        Para mantener m_φ ~ eV físicamente derivado se requiere una escala
        intermedia que NO existe limpiamente en SSEE actual.

    (c) Esto sugiere que SSEE actual no predice φ-DM en escala eV. O bien:
        - El sector φ-DM en P6 está en escala equivocada (debe rederivarse)
        - O necesitamos una escala física nueva (no presente en V3.6)
""")
print("=" * 80)
