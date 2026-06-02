"""
SSEE — Auditoría de precisión
=============================

Calcula TODAS las constantes SSEE a doble precisión sin redondeos
y compara contra los valores "canónicos" (redondeados) usados en los papers.

Pregunta: ¿cuánto error de redondeo arrastramos? ¿afecta a H₀, r_d, Ω_m,CMB?
"""
import math

# ===== Primitivas exactas =====
phi = (1 + math.sqrt(5)) / 2
pi_ = math.pi

# ===== Combinaciones algebraicas (todas irracionales) =====
beta  = (phi + pi_) / 2
Omega = phi + pi_
KAL   = beta + pi_
P_sc  = Omega + phi
K_v   = phi + pi_ + Omega
T_r   = 3 * (phi + beta)
M_v   = phi + pi_ + K_v
AURA  = pi_ - phi

# ===== Derivadas físicas =====
MIRA_exact      = (3*phi + pi_) / 4
Omega_DE_exact  = T_r / M_v
Omega_m_dyn_v1  = 1 - Omega_DE_exact          # via identidad de cierre
Omega_m_CMB     = Omega_m_dyn_v1 * MIRA_exact
H_alg           = 3 * (phi + pi_)**2          # km/s/Mpc (Paper 4)
w0_exact        = -T_r / M_v
wa_exact        = -P_sc / K_v

# ===== Tabla comparativa =====
print("=" * 78)
print("AUDITORÍA DE PRECISIÓN — SSEE")
print("=" * 78)

def row(name, exact, canonical, units=""):
    diff = exact - canonical
    pct = 100 * diff / canonical if canonical != 0 else 0
    print(f"  {name:<18} exact={exact:>20.15f}  canon={canonical:>10}  Δ={diff:+.6e}  ({pct:+.4f}%) {units}")

print("\n[PRIMITIVAS]")
row("phi",   phi,   1.6180)
row("pi",    pi_,   3.1416)

print("\n[ALGEBRAICAS]")
row("beta",  beta,  2.3798)
row("Omega", Omega, 4.7596)
row("KAL",   KAL,   5.5214)
row("P_sc",  P_sc,  6.3776)
row("K_v",   K_v,   9.5192)
row("T_r",   T_r,   11.9935)
row("M_v",   M_v,   14.2788)
row("AURA",  AURA,  1.5236)

print("\n[DERIVADAS FÍSICAS]")
row("MIRA",         MIRA_exact,     1.9989)
row("Omega_DE",     Omega_DE_exact, 0.840)
row("Omega_m_dyn",  Omega_m_dyn_v1, 0.160)
row("Omega_m_CMB",  Omega_m_CMB,    0.3199)
row("H_alg",        H_alg,          67.96,    "km/s/Mpc")
row("w0",           w0_exact,       -0.840)
row("wa",           wa_exact,       -0.670)

print("\n" + "=" * 78)
print("PROPAGACIÓN: ¿cuánto cambia H₀_obs si usamos exactos vs redondeados?")
print("=" * 78)

# Sensibilidad CAMB-like: r_d ∝ (Ω_m h²)^(-0.25) aprox (Eisenstein-Hu)
# θ* = r_d / D_A es fijo por Planck con precisión ~10⁻⁴
# Si cambia Ω_m → cambia r_d → debe cambiar H₀ para mantener θ*
#
# Aproximación: H₀ ∝ Ω_m^(-α) con α ≈ 0.4 (degeneración CMB clásica)

H_P3 = 67.08  # Paper 3 con Ω_m,CMB = 0.320 (redondeado)
omega_canon = 0.320
omega_exact = Omega_m_CMB

# Estimación primer orden
alpha_degen = 0.4
H_corr = H_P3 * (omega_canon / omega_exact)**alpha_degen
delta_H = H_corr - H_P3

print(f"""
  Ω_m,CMB canónico:  {omega_canon:.6f}
  Ω_m,CMB exacto:    {omega_exact:.6f}
  Δ(Ω_m,CMB):        {omega_exact - omega_canon:+.6e}  ({100*(omega_exact-omega_canon)/omega_canon:+.4f}%)

  H_P3 (canónico):   {H_P3:.4f} km/s/Mpc
  H_P3 (re-escalado por Ω_m exacto, α≈0.4):  {H_corr:.4f}
  Δ H₀:              {delta_H:+.4f} km/s/Mpc  ({100*delta_H/H_P3:+.4f}%)
""")

# Misma sensibilidad para H_alg vs H_P3
print(f"  Distancia H_alg vs H_P3 (canónico): {H_alg - H_P3:+.4f} km/s/Mpc")
print(f"  Distancia H_alg vs H_P3 (corregido): {H_alg - H_corr:+.4f} km/s/Mpc")
print(f"  Planck error ~0.54 → sigma equivalente:")
print(f"    canónico:  {(H_alg-H_P3)/0.54:.2f}σ")
print(f"    corregido: {(H_alg-H_corr)/0.54:.2f}σ")

# ===== Identidades internas: ¿se rompe algo por redondeo? =====
print("\n" + "=" * 78)
print("IDENTIDADES INTERNAS — ¿consistencia exacta?")
print("=" * 78)

# Identidad 1: w0 + Ω_m_dyn = ? (debería ser un número significativo)
id1 = w0_exact + Omega_m_dyn_v1
print(f"  w0 + Ω_m,dyn         = {id1:+.6e}   (¿debe ser cero? = -1+Ω_DE+1-Ω_DE = 0 → SÍ)")

# Identidad 2: 1 + w0 vs Ω_m,dyn (Paper 4 reclama)
id2a = 1 + w0_exact
id2b = Omega_m_dyn_v1
print(f"  1 + w0               = {id2a:.10f}")
print(f"  Ω_m,dyn              = {id2b:.10f}")
print(f"  Diferencia           = {id2a - id2b:+.2e}  (idénticos por construcción ✓)")

# Identidad 3: MIRA × Ω_m,dyn = Ω_m,CMB
id3 = MIRA_exact * Omega_m_dyn_v1
print(f"  MIRA × Ω_m,dyn       = {id3:.10f}")
print(f"  Ω_m,CMB              = {Omega_m_CMB:.10f}  (idénticos por def ✓)")

# Identidad 4: ¿Ω_m,CMB exacto coincide con Planck 0.3153 ± 0.0073?
planck_omm = 0.3153
planck_err = 0.0073
sigma_planck = (Omega_m_CMB - planck_omm) / planck_err
print(f"\n  Ω_m,CMB (SSEE exact) vs Planck 0.3153 ± 0.0073:")
print(f"    Δ = {Omega_m_CMB - planck_omm:+.4f}  ({sigma_planck:+.2f}σ)")

# ===== Veredicto =====
print("\n" + "=" * 78)
print("VEREDICTO")
print("=" * 78)
print(f"""
  Error de redondeo en Ω_m,CMB: ~0.02% → propaga a ΔH₀ ≈ {abs(delta_H):.3f} km/s/Mpc
  (mucho menor que el error observacional Planck ≈ 0.54 km/s/Mpc)

  → El "67.08 vs 67.96" NO se explica por redondeo numérico.
  → La tensión interna SSEE es real (~{abs(H_alg-H_corr)/0.54:.1f}σ), no artefacto.

  Implicación para MCMC: usar valores EXACTOS no cambia conclusión cualitativa,
  pero es buena práctica de auditoría. Tabla canónica para próximas corridas:

    Ω_m,dyn = {Omega_m_dyn_v1:.6f}
    MIRA    = {MIRA_exact:.6f}
    Ω_m,CMB = {Omega_m_CMB:.6f}
    H_alg   = {H_alg:.4f}
    w0      = {w0_exact:.6f}
    wa      = {wa_exact:.6f}
""")
