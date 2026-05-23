"""
ssee_mira_derivative_test.py — mecanismo #6: acoplamiento derivativo.

L_int = (X / M²) · L_DM  con X = (1/2)φ̇²

A diferencia de conformal exp(α·φ/M_pl):
  * No depende del VALOR del campo, sino de su VELOCIDAD al cuadrado.
  * Se enciende cuando φ rueda rápido (era de materia / pasado).
  * Se apaga cuando φ rueda despacio (hoy, w₀ = -0.840 → X chico).

Background coupled equations (Sami-Tsujikawa 2004, derivative coupling
schematic — se usa la forma más simple consistente con eq de Friedmann):

  ρ̇_m + 3H ρ_m + Q_X = 0     con Q_X ∝ (Ẋ/M²) · ρ_m
  KG: (P_X + 2X·P_XX) φ̈ + 3H P_X φ̇ + V_φ = -Q_X/φ̇  (back-reaction)

Simplificación esquemática (efecto leading order): el coupling efectivo
es β_eff(X) = (X/M²)^(1/2) · √2/M_pl ≈ φ̇/M, evaluado para X = ½φ̇².

Equivalente: β_eff = φ̇/M  (en unidades naturales).

Hoy: φ̇₀² = 2X₀, X₀ chico → β_eff hoy ≈ 0.
Era de materia: X grande → β_eff grande automático.

El test honesto: integrar el background y mirar:
  (a) R(z=1100) = exp[∫ β_eff·u dN] vs MIRA
  (b) ¿la integración llega a z=1100? (problema de conectividad)
"""
import numpy as np
from scipy.integrate import solve_ivp

from ssee_core import (KAL0, MIRA, W0 as w0,
                       OMEGA_DE as Om_DE, OMEGA_M_DYN as Om_m_dyn)

PHI = (1 + np.sqrt(5))/2
PI  = np.pi

# Escala M del coupling derivativo.
# SSEE Paper 10 fija M = 8.81 meV = Λ_SSEE.
# En unidades del código (H₀=1, ρ_crit=3), M/H₀ es muchísimo (M >> H₀).
# Pero (M_pl/M)² = (M_pl/Λ_SSEE)² es lo que aparece en el coupling
# adimensional: dim(X/M²) = M⁴/M² = M². Para que (X/M²)·L_DM tenga
# unidades correctas necesitamos otra M en el denominador.
# Convención: L_int = (X / M⁴) · ρ_m  (todo en unidades naturales)
#
# En unidades del código: X tiene unidades de ρ (potencia⁴),
# así que X/M⁴ es adimensional si M⁴ tiene unidades de ρ.
# Tomamos M⁴_code = M⁴/ρ_crit donde M es el valor físico.

# M = 8.81 meV. ρ_crit = (2.5 meV)⁴ aprox. M⁴/ρ_crit = (8.81/2.5)⁴ ≈ 154
# Por tanto M⁴_code ≈ 154 · 3 = 462 (en unidades código ρ_crit_code = 3)
# Esto significa: X/M⁴ ≈ X_code / 462 — coupling chico hoy.

# Para barrer y entender, probamos varios valores de M⁴_code:
M4_values = [1e-2, 1.0, 10.0, 100.0, 462.0, 1e4]   # diferentes escalas

lam       = np.sqrt(3.0 * Om_m_dyn)
alpha_pot = lam / np.sqrt(KAL0)

rho_r0   = 3.0 * 9.1e-5
rho_m0_initial = 3.0 * Om_m_dyn
rho_phi0 = 3.0 - rho_m0_initial - rho_r0
V0       = 3.0 * Om_DE

p_phi0  = w0 * rho_phi0
V_today = 0.5 * (rho_phi0 - p_phi0)
rhs_X   = 0.5 * (rho_phi0 + p_phi0)
X0      = KAL0 * rhs_X
phi0    = -np.log(V_today / V0) / alpha_pot
u0      = np.sqrt(2.0 * X0)


def V(phi):  return V0 * np.exp(-alpha_pot * phi)
def dV(phi): return -alpha_pot * V(phi)


def H2_of(u, phi, rho_m, rho_r):
    """Friedmann lineal puro (sin X²/M⁴ k-essence — es DERIVATIVE coupling,
    no UV completion). 3H² = X/KAL₀ + V + ρ_m + ρ_r, con X = H²·u²/2."""
    denom = 3.0 - u*u/(2.0*KAL0)
    if denom <= 0: return 0.0
    return max((V(phi) + rho_m + rho_r) / denom, 0.0)


def run(M4):
    """Acoplamiento derivativo schematic.

    Coupling efectivo: β_X(X) = u/√(M4)   (φ̇/M en unidades código)
    ρ_m equation: dρ_m/dN = -3ρ_m - β_X·u·ρ_m
    KG back-reaction: rhs_u incluye +β_X·ρ_m/H² · (correction factor)
    """
    def rhs(N, y):
        phi, u, rho_m = y
        a = np.exp(N)
        rho_r = rho_r0 * a**(-4)
        H2 = max(H2_of(u, phi, rho_m, rho_r), 1e-30)
        X = 0.5 * u*u * H2
        # β_eff derivativo: proporcional a u (φ̇/M ~ u·H/M, pero en
        # unidades código H~1, M⁴~M4): tomamos β_eff = u / sqrt(M4)
        beta_eff = u / np.sqrt(M4)
        rhs_u = -3.0*u + KAL0*(-dV(phi)/H2 + beta_eff*rho_m/H2)
        rhs_rho_m = -3.0*rho_m - beta_eff*u*rho_m
        return [u, rhs_u, rhs_rho_m]
    try:
        sol = solve_ivp(rhs, (0.0, -np.log(1101.0)),
                        [phi0, u0, rho_m0_initial],
                        t_eval=np.linspace(0, -np.log(1101.0), 4000),
                        method='Radau', rtol=1e-10, atol=1e-13, max_step=0.01)
        if not sol.success:
            return None, np.exp(-sol.t[-1])-1, None
        # R en varios z
        Rs = []
        for z_target in [0.5, 2, 10, 100, 500, 1100]:
            N_target = -np.log(1+z_target)
            if N_target < sol.t[-1]:
                Rs.append(None); continue
            i = np.argmin(np.abs(sol.t - N_target))
            a_i = np.exp(sol.t[i])
            R = sol.y[2, i] / (rho_m0_initial * a_i**(-3))
            Rs.append(R)
        z_reached = np.exp(-sol.t[-1])-1
        return Rs, z_reached, sol.y[:, -1]
    except Exception as e:
        return None, None, str(e)


print("="*78)
print("Mecanismo #6 — acoplamiento derivativo L_int = (X/M⁴)·ρ_m")
print("="*78)
print(f"\n  MIRA objetivo = {MIRA:.4f}")
print(f"  β_eff hoy ≈ u₀/√M4 = {u0:.3f}/√M4")
print()
print(f"{'M⁴_code':>10} {'R(0.5)':>9} {'R(2)':>9} {'R(10)':>9} "
      f"{'R(100)':>9} {'R(500)':>9} {'R(1100)':>10} {'z_end':>8}")
print("-"*78)

for M4 in M4_values:
    Rs, z_end, _ = run(M4)
    if Rs is None:
        z_str = f"{z_end:.2f}" if z_end else "?"
        print(f"{M4:10.2e}  (falla; integración hasta z={z_str})")
        continue
    line = f"{M4:10.2e}"
    for R in Rs:
        line += f" {R:9.4f}" if R is not None else "      n/d"
    line += f" {z_end:8.1f}"
    print(line)

print()
print("Interpretación:")
print("  R(1100)/MIRA cercano a 1 → mecanismo #6 funciona")
print("  R(1100) lejos de MIRA → no resuelve, sigue #7 o lectura B")
