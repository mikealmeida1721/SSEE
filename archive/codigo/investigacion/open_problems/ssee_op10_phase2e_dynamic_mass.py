"""
OP-10 Fase 2e — Masa dinámica por densidad ambiente (chameleon-like)
=====================================================================

IDEA (versión sobreviviente de "Osiris"/propuesta Mike): un SOLO campo φ cuya
masa efectiva depende de la densidad de materia ambiente,

        m_eff²(z) = ρ_m(z) / M_*² ,

de modo que en el universo temprano (denso) el campo es pesado y oscila
(= materia oscura) y al diluirse el universo se libera y rueda (= energía
oscura). Conversión DM→DE con la flecha temporal correcta.

PROTOCOLO FORWARD (cero fiteo):
  M_* = M = 9.62 meV — el cutoff UV YA bloqueado por Paper 10 (M⁴ = 5φ⁸ρ_crit).
  Ningún parámetro nuevo. Si falla, descarte limpio.

TRES REQUISITOS (declarados ANTES de calcular):
  R1 — DM viable: la energía de oscilación debe diluir como a⁻³ (CDM).
  R2 — Masa hoy: m_eff(z=0) debe ser compatible con m_φ = 36.95 eV
       (k_fs = 0.659 h/Mpc de Paper 6 asume masa constante).
  R3 — DE hoy: el campo debe poder rodar (m_eff ≲ H) para dar w₀ = −0.840.
"""
import numpy as np

# ============================================================
# CONSTANTES (todas bloqueadas, NO negociables)
# ============================================================
PHI   = (1 + np.sqrt(5)) / 2
PI    = np.pi
OMEGA = PI + PHI
KAL   = (PHI + 3*PI) / 2
T_R   = 3*(PHI + (PI+PHI)/2)
AURA  = (3*PHI + PI) / 2

H_MIRA       = 67.037                        # km/s/Mpc (canónico 2026-06-09)
h_MIRA       = H_MIRA / 100.0
RHO_CRIT_eV4 = 8.10e-11 * h_MIRA**2          # ρ_crit en eV⁴
OM_M_CMB     = 0.3199                        # Ω_m,CMB (sector Poisson)
H0_eV        = 1.45e-33 * h_MIRA / 0.6736    # H₀ en eV (escala con h)

M_UV         = (5 * PHI**8 * RHO_CRIT_eV4)**0.25   # = 9.62 meV (Paper 10)
M_PHI_TARGET = 36.95                         # eV (canónico Paper 6)

# z de igualdad materia-radiación (referencia física, no parámetro)
OM_R_H2 = 4.15e-5                            # Ω_r h² (fotones+neutrinos rel.)
Z_EQ    = OM_M_CMB * h_MIRA**2 / OM_R_H2 - 1

def rho_m(z):
    """Densidad de materia en eV⁴."""
    return OM_M_CMB * RHO_CRIT_eV4 * (1+z)**3

def m_eff(z, Mstar):
    """Masa efectiva inducida por densidad: m² = ρ_m/M_*²."""
    return np.sqrt(rho_m(z)) / Mstar

print("=" * 76)
print("OP-10 FASE 2e — masa dinámica m_eff²(z) = ρ_m(z)/M²  (M = M_UV, forward)")
print("=" * 76)
print(f"""
Escalas bloqueadas:
  M_* = M_UV          = {M_UV*1e3:.2f} meV   (Paper 10, cero fiteo)
  ρ_crit              = {RHO_CRIT_eV4:.3e} eV⁴
  m_φ target          = {M_PHI_TARGET:.2f} eV  (Paper 6 canónico)
  H₀                  = {H0_eV:.3e} eV
  z_eq (mat-rad)      = {Z_EQ:.0f}
""")

# ============================================================
# PARTE 1 — Historia de m_eff(z) y cruce con 36.95 eV
# ============================================================
print("=" * 76)
print("PARTE 1 — historia de la masa efectiva")
print("=" * 76)
m0 = m_eff(0, M_UV)
# cruce: m_eff(z_x) = M_PHI_TARGET  →  (1+z_x)^1.5 = target/m0
z_cross = (M_PHI_TARGET / m0)**(2/3) - 1
print(f"""
  m_eff(z=0)     = {m0*1e3:.3f} meV
  m_eff(z=1100)  = {m_eff(1100, M_UV):.2f} eV   (recombinación)
  m_eff(z_eq)    = {m_eff(Z_EQ, M_UV):.2f} eV   (igualdad, z={Z_EQ:.0f})

  Cruce m_eff = {M_PHI_TARGET} eV  →  z_x = {z_cross:.0f}
  Comparación: z_x/z_eq = {(1+z_cross)/(1+Z_EQ):.2f}
  → el cambio de fase cae en la época de igualdad materia-radiación
    SIN ajustar nada. (Dato sugestivo — pero sugerencia ≠ viabilidad.)
""")

# ============================================================
# PARTE 2 — los tres requisitos
# ============================================================
print("=" * 76)
print("PARTE 2 — test de los tres requisitos (declarados antes de calcular)")
print("=" * 76)

# --- R1: dilución CDM ---------------------------------------
# Campo oscilando en pozo cuadrático con masa dependiente del tiempo:
# invariante adiabático → número de cuantos n ∝ a⁻³; energía ρ = m(t)·n.
# Con m ∝ √ρ_m ∝ (1+z)^1.5 = a^(-3/2):
#     ρ_osc ∝ a^(-3/2) · a^(-3) = a^(-9/2)
# w_eff de dilución a^(-3(1+w)):  -3(1+w) = -4.5 → w = +0.5
dil_exp = 4.5
w_eff   = dil_exp/3 - 1
print(f"""
[R1] DILUCIÓN COMO CDM (necesita ρ ∝ a⁻³, w=0):
     ρ_osc = m(t)·n,  n ∝ a⁻³ (invariante adiabático)
     m ∝ √ρ_m ∝ a^(-3/2)   ⇒   ρ_osc ∝ a^(-{dil_exp})
     w_eff = +{w_eff:.1f}   →  diluye MÁS RÁPIDO que la radiación (a⁻⁴)
     ❌ FALLA — no es CDM: la 'materia oscura' se evaporaría antes de
        formar galaxias. NOTA: este fallo es INDEPENDIENTE de M_* (estructural:
        cualquier m ∝ ρ_m^(1/2) lo produce).""")

# --- R2: masa hoy vs k_fs -----------------------------------
gap_R2 = M_PHI_TARGET / m0
print(f"""
[R2] MASA HOY vs PAPER 6 (k_fs = 0.659 h/Mpc asume m = 36.95 eV constante):
     m_eff(0) = {m0*1e3:.3f} meV  vs  36.95 eV  →  factor {gap_R2:.2e}
     ❌ FALLA — la masa de hoy queda 10^{np.log10(gap_R2):.1f} por debajo del canónico.""")

# --- R3: ¿puede rodar como DE hoy? --------------------------
ratio_R3 = m0 / H0_eV
print(f"""
[R3] RODAR COMO DE HOY (necesita m_eff ≲ H₀ para w₀ = −0.840):
     m_eff(0)/H₀ = {ratio_R3:.2e}
     ❌ FALLA — el campo sigue 10^{np.log10(ratio_R3):.0f} veces más pesado que H₀:
        queda CLAVADO en su mínimo también hoy → w = −1 exacto, no −0.840.
        (Mismo síndrome de enganche de F1–F3.)""")

# ============================================================
# PARTE 3 — ¿algún M_* lo rescata? (diagnóstico, no fiteo)
# ============================================================
print("=" * 76)
print("PARTE 3 — ¿existe M_* que rescate? (diagnóstico)")
print("=" * 76)
Mstar_R2 = np.sqrt(rho_m(0)) / M_PHI_TARGET   # M_* que daría m_eff(0)=36.95 eV
print(f"""
  Para arreglar R2 haría falta M_* = √ρ_m(0)/m_φ = {Mstar_R2*1e9:.1f} neV
  — escala NUEVA, no bloqueada por ningún paper → sería fiteo. PROHIBIDO.
  Y aunque se permitiera: R1 es ESTRUCTURAL (no depende de M_*) y R3
  empeoraría (masa mayor hoy = más clavado). Ningún M_* rescata.
""")

# ============================================================
# VEREDICTO
# ============================================================
print("=" * 76)
print("VEREDICTO FASE 2e (versión tracking continuo)")
print("=" * 76)
print(f"""
  ✗ DESCARTADA la masa que SIGUE la densidad (m_eff² = ρ_m/M_*²):
      R1 ❌ estructural — ρ_osc ∝ a^(-4.5): no es CDM para ningún M_*
      R2 ❌ masa hoy 10⁵ por debajo del canónico (con M_* = M_UV)
      R3 ❌ campo clavado hoy (m/H ~ 10²⁹): w = −1, no −0.840

  INFORMACIÓN DEL DESCARTE (lo que este fallo enseña):
  1. La masa NO puede 'trackear' la densidad de forma continua. R1 lo
     prohíbe estructuralmente: mientras la masa cambie durante las
     oscilaciones, la energía no diluye como materia.
  2. ⇒ La transición DM↔DE tiene que ser un ESCALÓN (tipo transición de
     fase): masa ≈ constante = 36.95 eV desde z ≳ z_eq hacia atrás...
     no — hacia DELANTE en el tiempo, y el sector DE separado del pozo.
     Es decir: lo que se necesita es que el potencial CAMBIE DE FORMA
     una sola vez cerca de z_eq, no que respire con ρ_m.
  3. El escenario escalón + oscilaciones con masa constante ES el
     misalignment frío ya estudiado (cierre Osiris 2026-05-30): viable
     como DM, pero la amplitud inicial φ_i = 1.17×10⁻⁶ M_Pl no sale
     algebraica de (φ,π) → empata ΛCDM en conteo de parámetros.

  ESTADO OP-10 TRAS FASE 2e:
     F1 ❌  F2 ❌  F3 ❌  F4 ❌  Opción A ❌  2c ❌ (analítico)  2e-tracking ❌
     Patrón consolidado: ni pozo estático NI masa que trackea ρ_m.
     Lo único no excluido: transición de fase única cerca de z_eq cuyo
     disparador esté bloqueado por (φ,π) — mecanismo aún sin candidato.
""")
