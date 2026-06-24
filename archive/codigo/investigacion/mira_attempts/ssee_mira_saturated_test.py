"""
ssee_mira_saturated_test.py — test dinámico: ¿la α saturada produce MIRA?

Hipótesis bajo test (veta-2, V-L3-saturacion del Registro):
  α_sat = √(3/(2·KAL₀)) = √(3/(φ+3π)) ≈ 0.5212
  (derivado de la cota Amendola Ω_φ(φ-MDE) ≤ 1, renormalizado al campo SSEE)

Pregunta: con β_c = −α_sat (en lugar de −AURA=−3.998), integrando el fondo
acoplado Friedmann+KG hacia atrás desde hoy hasta z=1100, ¿qué factor de
retención R(z=1100) sale?

  Si R(z=1100) ≈ MIRA ≈ 1.999 → veta-2 deriva MIRA. Resultado positivo.
  Si R(z=1100) ≈ 1.0–1.3      → la cota saturada es muy débil. Negativo.
  Si R(z=1100) ≈ 0 / diverge   → otra falla independiente.

Método y estado base: idénticos a ssee_mira_mechanism.py — solo cambia β_c.
"""
import numpy as np
from scipy.integrate import solve_ivp

from ssee_core import (KAL0, MIRA, W0 as w0,
                       OMEGA_DE as Om_DE, OMEGA_M_DYN as Om_m_dyn)

PHI = (1 + np.sqrt(5))/2
PI  = np.pi

# ── α saturada (la nueva veta-2) ──────────────────────────────────────────────
alpha_sat = np.sqrt(3.0 / (PHI + 3.0*PI))     # = √(3/(2·KAL₀))
# La saturación fija |α|; el SIGNO se determina por la dirección del
# transporte de energía (retención requiere β_c>0 con esta convención).
beta_c    = +alpha_sat

print("="*72)
print("Test dinámico de MIRA con α saturada (veta-2)")
print("="*72)
print(f"\n  α_sat (saturación φ-MDE)        = {alpha_sat:.8f}")
print(f"  √(3/(2·KAL₀)) [check algebraico] = {np.sqrt(3/(2*KAL0)):.8f}")
print(f"  β_c = −α_sat                     = {beta_c:.8f}")
print(f"  AURA (test previo, fallido)      = -3.99794")
print(f"  Cociente β_c_sat / β_c_AURA      = {beta_c/(-3.99794):.4f}  "
      f"({abs(beta_c)/3.99794*100:.1f} % de la AURA original)")
print()
print(f"  Objetivo (MIRA)                  = {MIRA:.5f}")
print(f"  ln(MIRA)                         = {np.log(MIRA):.5f}")
print(f"  Δφ necesario (β_c·Δφ = ln MIRA)  = {np.log(MIRA)/alpha_sat:.4f} M_pl")
print()

# ── Parámetros EFT (idénticos a ssee_mira_mechanism.py) ──────────────────────
lam       = np.sqrt(3.0 * Om_m_dyn)
alpha_pot = lam / np.sqrt(KAL0)
# Régimen lineal puro: K(X) = X/KAL₀ (apagar el término X²/M⁴).
# Es el régimen donde la derivación de α_sat aplica (σ→0).

rho_r0   = 3.0 * 9.1e-5
rho_m0   = 3.0 * Om_m_dyn
rho_phi0 = 3.0 - rho_m0 - rho_r0
V0       = 3.0 * Om_DE

def V(phi):   return V0 * np.exp(-alpha_pot * phi)
def dV(phi):  return -alpha_pot * V(phi)
# Lineal: P_X = 1/KAL₀ constante; P_XX = 0
P_X_const = 1.0 / KAL0

def rho_phi_f(X, phi): return X / KAL0 + V(phi)
def p_phi_f(X, phi):   return X / KAL0 - V(phi)

# Estado de hoy (lineal: X = KAL₀·(ρ_φ+p_φ)/2)
p_phi0  = w0 * rho_phi0
V_today = 0.5 * (rho_phi0 - p_phi0)
rhs_X   = 0.5 * (rho_phi0 + p_phi0)
X0   = KAL0 * rhs_X                              # lineal exacto
phi0 = -np.log(V_today / V0) / alpha_pot
u0   = np.sqrt(2.0 * X0)

print(f"  Estado hoy: φ_0={phi0:.5f}, u_0=dφ/dlna|_hoy={u0:.5f}, X_0={X0:.5f}")
print()

# ── Friedmann lineal: H²(state) en cerrado ──────────────────────────────────
def H2_of(u, phi, rho_m, rho_r):
    """Lineal: 3H² = X/KAL₀ + V + ρ_m + ρ_r con X = H²·u²/2
    → H²·(3 − u²/(2·KAL₀)) = V + ρ_m + ρ_r"""
    denom = 3.0 - u*u/(2.0*KAL0)
    if denom <= 0:
        return 0.0, False
    H2 = (V(phi) + rho_m + rho_r) / denom
    return max(H2, 0.0), True

# ── ODE backward en N = ln(a), desde N=0 (hoy) hacia N=-ln(1101) ────────────
def rhs(N, y):
    phi, u, rho_m = y
    a = np.exp(N)
    rho_r = rho_r0 * a**(-4)
    H2, _ = H2_of(u, phi, rho_m, rho_r)
    H2 = max(H2, 1e-30)
    # KG lineal: u' = (-3·u/KAL₀ - V_φ/H² + β_c·ρ_m/H²) / (1/KAL₀)
    rhs_u = -3.0*u + KAL0 * (-dV(phi)/H2 + beta_c*rho_m/H2)
    # ρ_m: dρ_m/dN = -3ρ_m - β_c·u·ρ_m
    rhs_rho_m = -3.0*rho_m - beta_c*u*rho_m
    return [u, rhs_u, rhs_rho_m]

y0 = [phi0, u0, rho_m0]
N_span = (0.0, -np.log(1101.0))    # de hoy a z=1100, BACKWARD
N_eval = np.linspace(N_span[0], N_span[1], 4000)

sol = solve_ivp(rhs, N_span, y0, t_eval=N_eval, method='Radau',
                rtol=1e-10, atol=1e-13, max_step=0.01)

if not sol.success:
    print(f"  Integración INCOMPLETA: terminó en N={sol.t[-1]:.3f} "
          f"(z={np.exp(-sol.t[-1])-1:.2f})")
else:
    print(f"  Integración completa hasta z=1100 ✓")

# Resultados a varios z
print()
print(f"{'z':>10} {'R(z)':>12} {'φ(z)':>12} {'Δφ':>12} {'β_c·Δφ':>12}")
print("-"*60)
for z_target in [0, 0.5, 1, 2, 5, 10, 100, 500, 1100]:
    N_target = -np.log(1+z_target)
    if N_target < sol.t[-1]:
        continue
    i = np.argmin(np.abs(sol.t - N_target))
    phi_i, u_i, rho_m_i = sol.y[:, i]
    a_i = np.exp(sol.t[i])
    R = rho_m_i / (rho_m0 * a_i**(-3))
    Dphi = phi_i - phi0
    print(f"{z_target:10.2f} {R:12.6f} {phi_i:12.6f} {Dphi:+12.6f} "
          f"{beta_c*Dphi:+12.6f}")

# Veredicto
print()
print("="*72)
print("Veredicto")
print("="*72)
i_z1100 = np.argmin(np.abs(sol.t + np.log(1101.0)))
phi_end, u_end, rho_m_end = sol.y[:, i_z1100]
a_end = np.exp(sol.t[i_z1100])
R_end = rho_m_end / (rho_m0 * a_end**(-3))
ratio_to_mira = R_end / MIRA
print(f"\n  R(z=1100)         = {R_end:.5f}")
print(f"  MIRA (objetivo)   = {MIRA:.5f}")
print(f"  R / MIRA          = {ratio_to_mira:.4f}")
if abs(ratio_to_mira - 1) < 0.1:
    print("\n  >>> RESULTADO POSITIVO: α saturada deriva MIRA dentro del 10%.")
elif R_end > 1.0:
    pct = (R_end - 1.0)/(MIRA - 1.0)*100
    print(f"\n  >>> RESULTADO PARCIAL: R>1 (retención correcta), pero solo "
          f"{pct:.1f}% del camino a MIRA.")
elif R_end < 1.0 and R_end > 0:
    print("\n  >>> RESULTADO NEGATIVO: la materia se drena (signo invertido o "
          "saturación insuficiente).")
else:
    print("\n  >>> RESULTADO PATOLÓGICO: R≤0 o inestable.")
