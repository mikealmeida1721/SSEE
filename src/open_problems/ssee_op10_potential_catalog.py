"""
OP-10 Fase 1 — Catálogo de potenciales V(φ) candidatos
=========================================================

OBJETIVO: encontrar V(φ) que satisfaga simultáneamente:
  (R1) Plateau a φ grande → reproduce w_0 = -0.840, w_a = -0.670
  (R2) Mínimo con V''(φ_min) = m_φ² → reproduce m_φ = 36.95 eV  (resuelve OP-9)
  (R3) α-attractor compatible (inflación) → preserva n_s = 1−φ⁻⁷
  (R4) Coeficientes algebraicos en φ, π, KAL, ...                (Postulado D)

CONSTRAINTS DUROS (ya bloqueados por papers previos, NO negociables):
  - α      = φ⁴/3          ≈ 2.2847    (Paper 1)
  - λ      = √(3 Ω_m,dyn)  ≈ 0.6929    (Paper 7)
  - V₀     = 0.840 ρ_crit              (Paper 7)
  - M⁴     = 5 φ⁸ ρ_crit               (Paper 10)
  - m_φ    = 36.95 eV (canónico, forward-derivation Σm_ν·MULT; target, OP-9)

ESPACIO DE BÚSQUEDA: dado que P7 fija V_DE(φ) = V₀ e^(-λφ), la búsqueda real es
    V(φ) = V₀ e^(-λφ) + ΔV(φ)
donde ΔV debe:
  (a) → 0 a φ grande (no contaminar DE plateau)
  (b) crear mínimo a algún φ_min con V''(φ_min) = m_φ²
  (c) ser algebraicamente limpio
"""
import numpy as np

# ============================================================
# CONSTANTES ALGEBRAICAS Y FÍSICAS
# ============================================================
PHI   = (1 + np.sqrt(5)) / 2
PI    = np.pi
OMEGA = PI + PHI
BETA  = (PI + PHI) / 2
KAL   = (PHI + 3*PI) / 2
P_SC  = OMEGA + PHI
K_V   = PHI + PI + OMEGA
T_R   = 3 * (PHI + BETA)
M_V   = PHI + PI + K_V

# Bloqueados por papers previos
ALPHA_LOCK     = PHI**4 / 3        # ≈ 2.2847   (Paper 1)
OM_M_DYN       = 0.160              # Paper 4/7
LAMBDA_LOCK    = np.sqrt(3 * OM_M_DYN)  # ≈ 0.6929  (Paper 7)
OM_DE          = 0.840
RHO_CRIT_eV4   = 8.10e-11           # ρ_crit en eV⁴ (h=0.6796)
V0_LOCK        = OM_DE * RHO_CRIT_eV4
M4_LOCK        = 5 * PHI**8 * RHO_CRIT_eV4
M_UV           = M4_LOCK**0.25      # ≈ 8.81 meV

# Target — m_φ canónico por forward-derivation (Paper 6 partícula canónica)
#   Σm_ν = (Ω/(KAL·T_R))·0.960318 = 0.06903 eV ; MULT = Ω⁴ + AURA·KAL = 535.28
AURA           = (3*PHI + PI) / 2
M_PHI_TARGET   = (OMEGA/(KAL*T_R))*0.960318 * (OMEGA**4 + AURA*KAL)  # = 36.95 eV
M_PHI2_TARGET  = M_PHI_TARGET**2    # eV²

print("=" * 76)
print("OP-10 FASE 1 — Catálogo de potenciales V(φ) candidatos")
print("=" * 76)
print(f"\nConstantes bloqueadas (papers previos):")
print(f"  α = φ⁴/3             = {ALPHA_LOCK:.6f}         (Paper 1)")
print(f"  λ = √(3 Ω_m,dyn)     = {LAMBDA_LOCK:.6f}        (Paper 7)")
print(f"  V₀ = 0.840 ρ_crit    = {V0_LOCK:.3e} eV⁴   (Paper 7)")
print(f"  M⁴ = 5 φ⁸ ρ_crit     = {M4_LOCK:.3e} eV⁴   (Paper 10)")
print(f"  M  = (M⁴)^(1/4)      = {M_UV*1e3:.3f} meV")
print(f"\nTarget a reproducir (OP-9):")
print(f"  m_φ                  = {M_PHI_TARGET} eV")
print(f"  m_φ²                 = {M_PHI2_TARGET:.4f} eV²")

# ============================================================
# CATÁLOGO DE 4 FAMILIAS DE ΔV(φ)
# ============================================================
print("\n" + "=" * 76)
print("CATÁLOGO — 4 familias de ΔV(φ) que se SUMAN al V_DE bloqueado")
print("=" * 76)

# Para análisis cualitativo (mínimo, curvatura), trabajamos en unidades naturales
# donde φ es adimensional (φ/M_pl). El mínimo aparece donde V'(φ) = 0.

# ----------------------------------------------------------------
# FAMILIA 1: Exponencial decreciente + harmónico
#   ΔV₁(φ) = (m²/2) (φ - φ₀)²
# Suma:    V(φ) = V₀ e^(-λφ) + (m²/2)(φ-φ₀)²
#
# Mínimo: V'(φ_min) = 0 → -λ V₀ e^(-λφ) + m²(φ-φ₀) = 0
#         (transcendental, requiere solver numérico)
# Curvatura: V''(φ_min) = λ² V₀ e^(-λφ_min) + m²
#
# Free params: m, φ₀
# Algebraic candidates for (m,φ₀): m=H₀, φ₀=φ (golden ratio en unidades nat),
#   m²=H₀·M_UV,  φ₀=1/KAL, etc.
#
# Pros:    simple, mínimo bien definido si m² > 0
# Contras: ΔV NO → 0 a φ grande (crece cuadráticamente) → contamina DE
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# FAMILIA 2: α-attractor tanh² + binding sech²
#   ΔV₂(φ) = (m²/2) sech²(φ/f)
# Suma:    V(φ) = V₀ tanh²(φ/f) + (m²/2) sech²(φ/f)
# (Aquí cambiamos V_DE de exp→tanh² para preservar α-attractor)
#
# Mínimo: en φ=0 (tanh²(0)=0, sech²(0)=1) → V_min = m²/2
# V''(0) = 2V₀/f² - m²/f² = (2V₀ - m²)/f²
# Para curvatura POSITIVA → m² < 2V₀ → mínimo
# Para curvatura NEGATIVA → m² > 2V₀ → φ=0 es máximo, mínimo en otro lugar
#
# Free params: m, f
# Algebraic candidates for f: f = √(2/3α), f = 1/φ⁵, f = φ²
#
# Pros:    tanh² es α-attractor universalidad (R3 ✓), tiene mínimo natural
# Contras: cambia la forma de V_DE respecto al Paper 7 — requiere recomputar w₀
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# FAMILIA 3: Axion-like + exponencial
#   ΔV₃(φ) = M_UV⁴ [1 - cos(φ/f)]
# Suma:    V(φ) = V₀ e^(-λφ) + M_UV⁴ [1 - cos(φ/f)]
#
# Mínimo: cos(φ_min/f) = 1 (cancelando el cos), pero modulado por exp
# Si f << 1/λ: el coseno oscila rápido sobre el exp suave → muchos mínimos
# Curvatura en mínimo: V'' ≈ M_UV⁴/f² (dominado por axion término)
#
# Free params: f (axion decay constant)
# Algebraic candidates for f: f = 1/M_UV (m_φ ≈ M_UV/f), f = φ/π
#
# Pros:    M_UV ya bloqueado → menos free params, mínimos múltiples (vacuum landscape)
# Contras: φ-DM oscilante (no estático) — requiere análisis condensado
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# FAMILIA 4: α-attractor puro Starobinsky (sin término extra)
#   V₄(φ) = V₀ [1 - e^(-√(2/3α) φ)]²
# Mínimo: en φ=0 (V=0)
# V''(0) = V₀ · 4 · (2/3α) = (8/3α) V₀
#
# Free params: NINGUNO (V₀ y α ya bloqueados!)
# Predicción única: m_φ² = (8/3α) V₀
#
# Esta es la opción más restrictiva — si NO da 36.95 eV, está descartada.
# Si SÍ da 36.95 eV, sería el "milagro" (cero parámetros nuevos).
# ----------------------------------------------------------------

# ============================================================
# EVALUACIÓN PRELIMINAR DE FAMILIA 4 (la más restrictiva)
# ============================================================
print("\n--- Familia 4 (Starobinsky α-attractor puro, cero free params) ---")
print("    V(φ) = V₀ [1 - exp(-√(2/3α) φ)]²")
print("    Mínimo en φ=0, V_min=0")
mphi2_F4 = (8 / (3 * ALPHA_LOCK)) * V0_LOCK  # eV⁴ (V tiene dim de E⁴, m² de E²)
# Pero V'' en mínimo da unidades de energía⁴/φ². Si φ adimensional, V'' tiene
# unidades de eV⁴. La masa física m_φ requiere normalización canónica del campo.
# Si φ es canónicamente normalizado (dimensión de eV), V'' tiene eV² directo.
# Aquí asumo normalización canónica de Starobinsky donde φ = M_pl · canonical_scalar
# Entonces m² = V''(0) en unidades de (V₀/M_pl²).
M_PL_eV = 2.435e27  # eV (Planck mass reducida)
m2_F4_physical = (8 / (3 * ALPHA_LOCK)) * V0_LOCK / M_PL_eV**2  # eV²
m_F4_physical = np.sqrt(m2_F4_physical)  # eV
print(f"    V''(0) = (8/3α) V₀         = {mphi2_F4:.3e} eV⁴")
print(f"    m²_φ   = V''(0) / M_pl²    = {m2_F4_physical:.3e} eV²")
print(f"    m_φ    = √m²_φ             = {m_F4_physical:.3e} eV")
print(f"    Target = {M_PHI_TARGET} eV → ratio m_F4/target = {m_F4_physical/M_PHI_TARGET:.3e}")

if m_F4_physical / M_PHI_TARGET < 1e-10 or m_F4_physical / M_PHI_TARGET > 1e10:
    print("    ⚠ DESCARTADA: m_φ predicho ≠ 36.95 eV por orden(es) de magnitud.")
    print("    Familia 4 (cero parámetros) no funciona — necesitamos free param.")
else:
    print("    ✓ Familia 4 candidata viable — investigar normalización canónica.")

# ============================================================
# TABLA RESUMEN DE LAS 4 FAMILIAS
# ============================================================
print("\n" + "=" * 76)
print("TABLA RESUMEN — 4 familias candidatas")
print("=" * 76)
print(f"""
| # | Forma de ΔV(φ)                | Free params  | Pros                    | Contras                       |
|---|------------------------------|--------------|-------------------------|-------------------------------|
| 1 | (m²/2)(φ-φ₀)²                | m, φ₀ (2)    | Mínimo trivial          | NO → 0 a φ grande → cont. DE  |
| 2 | (m²/2) sech²(φ/f)            | m, f (2)     | α-attractor ✓ R3        | Cambia V_DE forma (recomp w₀) |
| 3 | M_UV⁴ [1−cos(φ/f)]           | f (1)        | M_UV ya bloqueado       | φ oscilante, no estático      |
| 4 | NINGÚN ΔV (Starobinsky puro) | 0            | Cero params (milagro)   | Predice m_φ ≠ 36.95 eV ❌      |
""")

# ============================================================
# DECISIÓN PARA FASE 2
# ============================================================
print("=" * 76)
print("DECISIÓN PARA FASE 2 (próxima sesión)")
print("=" * 76)
print(f"""
  FAMILIA 1 — DESCARTAR en Fase 2 (contamina plateau DE a φ grande).
              Salvo que (m²/2)(φ-φ₀)² se suprima por algún factor exp(-...).

  FAMILIA 2 — PRIORITARIO en Fase 2. tanh² satisface α-attractor (P1, OP-2).
              Ataque: testar 3 valores algebraicos de f: √(2/3α), 1/φ⁵, φ².
              Para cada f, calcular w₀(t) y comparar con −0.840.

  FAMILIA 3 — INVESTIGAR en Fase 2 si tiempo permite. φ-DM oscilante NO es lo
              que P6 asume (P6 trata φ-DM como background coherente). Probable
              descarte tras análisis de coherencia.

  FAMILIA 4 — DESCARTADA ya (predicción m_φ no coincide con 36.95 eV).
              Excepto si la normalización canónica que asumí está mal.

  PRÓXIMO SCRIPT: ssee_op10_family2_dynamics.py
    - Resolver Klein-Gordon en FRW para V(φ) = V₀ tanh²(φ/f) + (m²/2)sech²(φ/f)
    - Barrido en f ∈ {{√(2/3α), 1/φ⁵, φ²}} y m ∈ {{H₀, √(H₀·M_UV), m_φ}}
    - Calcular w₀(z=0), w_a, m_φ² = V''(φ_min)
    - Test combinado: ¿algún (f, m) cumple las 3 restricciones?

  CRITERIO DE ÉXITO FASE 2: encontrar (f, m) ALGEBRAICOS tal que
    |w₀ + 0.840| < 0.01 AND |w_a + 0.670| < 0.05 AND |m_φ − 36.95| < 0.5 eV
""")
