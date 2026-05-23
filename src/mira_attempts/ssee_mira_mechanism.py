"""
ssee_mira_mechanism.py
======================
Test del mecanismo candidato para MIRA — campaña de verificación 2026-05-22.

Hipótesis (Ruta B / "retención"): el "0.320" que el CMB exige no es una
sustancia. Es el campo φ reteniendo energía en el sector materia mientras el
universo es denso — un acoplamiento conformal de quintaesencia (el L_int de
Paper 7, re-apuntado de materia oscura a materia REAL, como exige Paper 1).

Física estándar (quintaesencia acoplada, Amendola 2000):
  ρ̇_m + 3Hρ_m = -β_c φ̇ ρ_m        →  ρ_m = ρ_m,0 a⁻³ · exp[-β_c(φ-φ_0)]
  KG: (P_X+2X P_XX)φ̈ + 3H P_X φ̇ + V_φ = β_c ρ_m

  →  R(a) ≡ ρ_m / (ρ_m,0 a⁻³) = exp[-β_c(φ(a)-φ_0)]   (factor de peso extra)
  →  Γ (tasa de retención) = |β_c φ̇| ;  desacople en Γ = H

Predicciones a verificar (NADA se ajusta — β_c = -AURA fijo):
  R(hoy)      = 1            (normalización)
  R(temprano) = MIRA ≈ 1.999  ⟺  AURA·Δφ_total = ln(MIRA) ≈ 0.693
  z* (desacople) ∈ (2, 1100)  para no chocar con DESI ni con el CMB

Método: integración BACKWARD del fondo acoplado desde hoy (a=1) hasta z≈1100.
El estado de hoy se fija con w_φ(hoy) = w₀ = -0.840 (valor algebraico SSEE).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from ssee_core import (KAL0, AURA, MIRA, W0 as w0,
                       OMEGA_DE as Om_DE, OMEGA_M_DYN as Om_m_dyn)

# ── Parámetros EFT (idénticos a ssee_eft_verification.py) ────────────────────
lam       = np.sqrt(3.0 * Om_m_dyn)      # λ canónico
alpha_pot = lam / np.sqrt(KAL0)          # pendiente del potencial
M4        = 1.0                          # M⁴ = ρ_crit
beta_c    = -AURA                        # acoplamiento conformal (Paper 7)

# ── Unidades H₀ = 1  →  ρ_tot(hoy) = 3 ───────────────────────────────────────
rho_r0   = 3.0 * 9.1e-5                  # radiación hoy
rho_m0   = 3.0 * Om_m_dyn                # materia real acoplada hoy
rho_phi0 = 3.0 - rho_m0 - rho_r0         # φ hoy (universo plano)
V0       = 3.0 * Om_DE                   # amplitud del potencial

def V(phi):   return V0 * np.exp(-alpha_pot * phi)
def dV(phi):  return -alpha_pot * V(phi)
def P_X(X):   return 1.0 / KAL0 + 2.0 * X / M4
P_XX = 2.0 / M4
def rho_phi_f(X, phi): return X / KAL0 + 3.0 * X * X / M4 + V(phi)
def p_phi_f(X, phi):   return X / KAL0 + X * X / M4 - V(phi)

# ── Estado de hoy: imponemos w_φ(hoy) = w₀ = -0.840 ──────────────────────────
p_phi0  = w0 * rho_phi0
V_today = 0.5 * (rho_phi0 - p_phi0)              # ρ_φ - p_φ = 2V
rhs_X   = 0.5 * (rho_phi0 + p_phi0)              # ρ_φ + p_φ = 2(X/KAL0 + 2X²/M4)
# 2/M4 · X² + 1/KAL0 · X - rhs_X = 0
X0   = (-1.0 / KAL0 + np.sqrt((1.0 / KAL0) ** 2 + 4.0 * (2.0 / M4) * rhs_X)) / (2.0 * (2.0 / M4))
phi0 = -np.log(V_today / V0) / alpha_pot
u0   = np.sqrt(2.0 * X0)                          # dφ/dlna hoy (H₀=1); rueda a φ mayor

print("=" * 68)
print("  SSEE — Test del mecanismo MIRA (Ruta B: retención conformal)")
print("=" * 68)
print(f"  β_c        = -AURA = {beta_c:.5f}   (fijo, Paper 7)")
print(f"  α_pot      = {alpha_pot:.5f}    V₀ = {V0:.5f}    M⁴ = {M4:.3f}")
print(f"  Estado hoy : φ₀ = {phi0:.5f}   u₀ = dφ/dlna = {u0:.5f}")
print(f"               X₀ = {X0:.5f}   w_φ(hoy) = {w0:.4f} (impuesto)")
print(f"  Objetivo   : R_temprano = MIRA = {MIRA:.5f}")
print(f"               AURA·Δφ_total = ln(MIRA) = {np.log(MIRA):.5f}")
print(f"               ⟹ Δφ_total requerido = {np.log(MIRA)/AURA:+.5f}")
print("-" * 68)

# ── Friedmann con K(X) cuadrática: 3H² = X/KAL₀+3X²/M⁴+V+ρ_m+ρ_r, X=½H²u² ──
# → 0.75(u⁴/M⁴)·H⁴ + (0.5u²/KAL₀−3)·H² + S = 0   (cuadrática en H²)
# Raíz física = la menor (conectada con Friedmann estándar cuando u→0).
breakdown_N = []

def H2_of(phi, u, rho_m, a):
    S = V(phi) + rho_m + rho_r0 * a ** (-4)
    A = 0.75 * u ** 4 / M4
    B = 0.5 * u ** 2 / KAL0 - 3.0
    if A < 1e-25:
        return -S / B                      # límite lineal (B<0, S>0)
    disc = B * B - 4.0 * A * S
    if disc < 0.0:
        breakdown_N.append(a)              # k-essence sin solución real
        return -B / (2.0 * A)              # vértice (marginal)
    return (-B - np.sqrt(disc)) / (2.0 * A)

# ── Sistema acoplado Friedmann + Klein-Gordon, variable N = ln(a) ────────────
def rhs(N, Y):
    phi, u, rho_m = Y
    a = np.exp(N)
    rho_r = rho_r0 * a ** (-4)

    H2 = H2_of(phi, u, rho_m, a)
    H  = np.sqrt(H2)
    X  = 0.5 * H2 * u * u

    pX  = P_X(X)
    eff = pX + 2.0 * X * P_XX
    rho_tot = 3.0 * H2
    p_tot   = p_phi_f(X, phi) + rho_r / 3.0

    # KG:  eff·H²·u' = β_c ρ_m - V_φ - 3H²P_X u + eff·(ρ+p)·u/2
    uprime = (beta_c * rho_m - dV(phi) - 3.0 * H2 * pX * u
              + eff * (rho_tot + p_tot) * u / 2.0) / (eff * H2)
    drho_m = -3.0 * rho_m - beta_c * u * rho_m
    return [u, uprime, drho_m]

# ── Integración hacia atrás: hoy (N=0) → z≈1100 (N≈-7) ───────────────────────
N_min = np.log(1.0 / (1.0 + 1100.0))
N_eval = np.linspace(0.0, N_min, 4000)
sol = solve_ivp(rhs, [0.0, N_min], [phi0, u0, rho_m0], t_eval=N_eval,
                rtol=1e-10, atol=1e-13, method='DOP853')

if not sol.success:
    print(f"  ERROR de integración: {sol.message}")
    # seguir con lo integrado
a_arr   = np.exp(sol.t)
z_arr   = 1.0 / a_arr - 1.0
phi_arr = sol.y[0]
u_arr   = sol.y[1]
rhom_arr = sol.y[2]
n_ok    = len(sol.t)
print(f"  Integrados {n_ok} puntos; z_max alcanzado = {z_arr.max():.1f}")

# ── Observables ──────────────────────────────────────────────────────────────
R_arr      = rhom_arr * a_arr ** 3 / rho_m0          # factor de peso
R_check    = np.exp(-beta_c * (phi_arr - phi0))       # = exp[AURA·Δφ], cross-check
Gamma_H    = np.abs(beta_c * u_arr)                   # Γ/H = |β_c·u|

dphi_total = phi_arr[-1] - phi0
print(f"  Cross-check R: máx|R - exp(AURA·Δφ)| = {np.max(np.abs(R_arr-R_check)):.2e}")

# z* : desacople donde Γ/H = 1
z_star = None
for i in range(1, n_ok):
    if (Gamma_H[i] - 1.0) * (Gamma_H[i-1] - 1.0) < 0:
        f = (1.0 - Gamma_H[i-1]) / (Gamma_H[i] - Gamma_H[i-1])
        z_star = z_arr[i-1] + f * (z_arr[i] - z_arr[i-1])
        break

def solve_H2(phi, u, rho_m, a):
    return H2_of(phi, u, rho_m, a)

print("-" * 68)
print(f"  {'z':>8} {'R(z)':>10} {'Om_m,ef':>9} {'Gam/H':>9} {'w_phi':>9} {'Om_phi':>8}")
for zt in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1100.0]:
    j = int(np.argmin(np.abs(z_arr - zt)))
    a = a_arr[j]
    H2 = solve_H2(phi_arr[j], u_arr[j], rhom_arr[j], a)
    X  = 0.5 * H2 * u_arr[j] ** 2
    rphi = rho_phi_f(X, phi_arr[j])
    wphi = p_phi_f(X, phi_arr[j]) / rphi
    Om_ef = R_arr[j] * Om_m_dyn
    Om_phi = rphi / (3.0 * H2)
    print(f"  {z_arr[j]:8.2f} {R_arr[j]:10.4f} {Om_ef:9.4f} {Gamma_H[j]:9.4f}"
          f" {wphi:9.4f} {Om_phi:8.4f}")

print("-" * 68)
print(f"  R(hoy)            = {R_arr[0]:.5f}   (debe ser 1)")
print(f"  R(z≈{z_arr[-1]:.0f})       = {R_arr[-1]:.5f}   (objetivo MIRA = {MIRA:.4f})")
print(f"  Ω_m,efectiva temprana = {R_arr[-1]*Om_m_dyn:.5f}   (objetivo ≈ 0.320)")
print(f"  Δφ_total          = {dphi_total:+.5f}   (objetivo {np.log(MIRA)/(-beta_c)*np.sign(1):+.5f})")
print(f"  AURA·Δφ_total     = {-beta_c*dphi_total:+.5f}   (objetivo ln MIRA = {np.log(MIRA):.5f})")
if z_star is not None:
    print(f"  z* (desacople Γ=H)= {z_star:.2f}")
else:
    print(f"  z* (desacople Γ=H)= NO se cruza Γ=H  (Γ/H rango "
          f"[{Gamma_H.min():.3f}, {Gamma_H.max():.3f}])")
if breakdown_N:
    print(f"  ⚠ k-essence sin solución real de H en {len(breakdown_N)} puntos")

# ── Veredicto ────────────────────────────────────────────────────────────────
print("-" * 68)
mag_ratio = abs(-beta_c * dphi_total) / np.log(MIRA)
ok_mach   = abs(R_arr[0] - 1.0) < 1e-6
ok_early  = abs(R_arr[-1] - MIRA) / MIRA < 0.5
print("  VEREDICTO")
print(f"    Maquinaria  R(hoy)=1, cross-check exacto ........ {'OK' if ok_mach else 'FALLA'}")
print(f"    R temprano ≈ MIRA .............................. {'OK' if ok_early else 'FALLA'}"
      f"   (R={R_arr[-1]:.4f} vs {MIRA:.4f})")
print(f"    Magnitud: |AURA·Δφ| / ln(MIRA) ................. ×{mag_ratio:.1f}"
      f"   {'(esperado ×1)' }")
print(f"    Signo: R temprano {'>1 (carga)' if R_arr[-1] > 1 else '<1 (DRENA materia)'}"
      f" .............. {'OK' if R_arr[-1] > 1 else 'FALLA'}")
early_coupled = Gamma_H[-1] > 1.0
print(f"    Timing: retención activa temprano (Γ/H>1) ...... "
      f"{'OK' if early_coupled else 'FALLA'}  (Γ/H(z_max)={Gamma_H[-1]:.3f})")
print(f"\n  → El acoplamiento conformal con β_c=AURA NO reproduce MIRA:")
print(f"    excursión ~×{mag_ratio:.0f} excesiva, signo invertido (drena en vez de")
print(f"    cargar), y timing invertido (campo thawing, no freezing).")
print("=" * 68)
