"""
Diagnóstico: probar ambos signos de α_sat y ambos signos de u₀,
recordar resultados de R(z=1100), reportar matriz 2×2.
"""
import numpy as np
from scipy.integrate import solve_ivp
from ssee_core import (KAL0, MIRA, W0 as w0,
                       OMEGA_DE as Om_DE, OMEGA_M_DYN as Om_m_dyn)

PHI = (1 + np.sqrt(5))/2
PI  = np.pi

alpha_sat = np.sqrt(3.0 / (PHI + 3.0*PI))

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
u0_abs  = np.sqrt(2.0 * X0)


def V(phi):  return V0 * np.exp(-alpha_pot * phi)
def dV(phi): return -alpha_pot * V(phi)


def H2_of(u, phi, rho_m, rho_r):
    denom = 3.0 - u*u/(2.0*KAL0)
    if denom <= 0: return 0.0
    return max((V(phi) + rho_m + rho_r) / denom, 0.0)


def run(beta_c, u0_sign):
    u0 = u0_sign * u0_abs
    def rhs(N, y):
        phi, u, rho_m = y
        a = np.exp(N)
        rho_r = rho_r0 * a**(-4)
        H2 = max(H2_of(u, phi, rho_m, rho_r), 1e-30)
        rhs_u = -3.0*u + KAL0*(-dV(phi)/H2 + beta_c*rho_m/H2)
        rhs_rho_m = -3.0*rho_m - beta_c*u*rho_m
        return [u, rhs_u, rhs_rho_m]
    try:
        sol = solve_ivp(rhs, (0.0, -np.log(1101.0)), [phi0, u0, rho_m0_initial],
                        t_eval=np.linspace(0, -np.log(1101.0), 4000),
                        method='Radau', rtol=1e-10, atol=1e-13, max_step=0.01)
        if not sol.success:
            return None, sol.t[-1]
        N_end = sol.t[-1]
        z_end = np.exp(-N_end) - 1
        # R at final point
        a_end = np.exp(N_end)
        R_end = sol.y[2, -1] / (rho_m0_initial * a_end**(-3))
        return R_end, z_end
    except Exception:
        return None, None


print(f"Test diagnóstico: α_sat = {alpha_sat:.6f}")
print(f"u₀ |magnitud|        = {u0_abs:.6f}")
print(f"MIRA objetivo        = {MIRA:.6f}")
print()
print(f"{'β_c':>10} {'u₀ signo':>10} {'R(z_end)':>12} {'z_end alcanzado':>16}")
print("-"*55)

for beta_c, sign in [(+alpha_sat, +1), (+alpha_sat, -1),
                      (-alpha_sat, +1), (-alpha_sat, -1)]:
    R, z = run(beta_c, sign)
    if R is None:
        print(f"{beta_c:+10.4f} {sign:+10d}  (INTEGRACIÓN FALLÓ)")
    else:
        achieved = "z=1100 ✓" if z > 1000 else f"z={z:.2f} (corte)"
        print(f"{beta_c:+10.4f} {sign:+10d} {R:12.5f} {achieved:>16}")

print()
print("Interpretación:")
print("  R > 1 = retención (materia gravita más que a⁻³)")
print("  R < 1 = drenaje (materia gravita menos)")
print("  z_end < 1100 = la integración rompió antes de matter era")
