"""
ssee_mira_phimde_forward_test.py — mecanismo #7: integración FORWARD desde φ-MDE.

Lógica:
  Hasta ahora todos los tests (β_c=±AURA, β_c=±α_sat, derivativo M físico)
  integran BACKWARD desde hoy (w₀=−0.840, Ω_m=0.16) hacia z=1100 y rompen
  antes de llegar a la era de materia. Eso sugiere que NO existe una
  trayectoria suave que conecte ambos extremos.

  #7 invierte la pregunta: arrancar en el φ-MDE (donde el atractor existe
  por construcción Amendola) y dejar que la dinámica corra HACIA hoy.
  Si la trayectoria llega naturalmente a w_φ≈−0.840 con Ω_m≈0.32, MIRA
  emerge como cociente entre los dos regímenes (lectura B). Si no llega,
  las dos islas dinámicas son inalcanzables entre sí — lectura B se
  refuerza.

  Punto inicial: φ-MDE en variables Copeland (x,y,σ→0):
    x_φMDE = −α√6 · KAL₀ / 3   (raíz lineal, σ→0, ver Step D2)
    y_φMDE = 0   (matter dominated → V≈0)
    Ω_φ(φMDE) = 2·KAL₀·x²/3  → con α=α_sat, Ω_φ=1 (saturado)

  Pero Ω_φ=1 no deja lugar a ρ_m. Hay que arrancar levemente off-attractor:
  tomar Ω_m_init grande (0.99), Ω_φ chico (0.01) en z≈1100, w_φ libre.
  El α elegido determina si el attractor jala hacia DE-fase o se queda
  pegado a φ-MDE.

α valores a barrer:
  • α_sat = √(3/(φ+3π)) ≈ 0.5212   (veta-2)
  • α_pot = √(3Ω_m)/√KAL₀ ≈ 0.295   (Copeland exponencial puro)
  • α grande: 1.0, 1.5                (testear sensibilidad)
"""
import numpy as np
from scipy.integrate import solve_ivp

from ssee_core import KAL0, MIRA, W0 as w0, OMEGA_DE as Om_DE, OMEGA_M_DYN as Om_m_dyn

PHI = (1 + np.sqrt(5))/2
PI  = np.pi

alpha_sat = np.sqrt(3.0 / (PHI + 3.0*PI))
alpha_pot_canonical = np.sqrt(3.0 * Om_m_dyn) / np.sqrt(KAL0)

rho_r0   = 3.0 * 9.1e-5
V0_today = 3.0 * Om_DE


def run_forward(alpha, N_init=-7.0, Om_phi_init=0.01):
    """Integra forward desde N=N_init hacia N=0 (hoy).

    Estado inicial: matter dominated con scalar field encendido al 1%.
    """
    a_init = np.exp(N_init)
    rho_total_init = 3.0  # H=1 normalization in code units
    rho_r_init = rho_r0 * a_init**(-4)
    # Trick: dejamos que ρ_total a este N sea consistente con matter
    # dominated. ρ_m = (1−Ω_φ−Ω_r)·ρ_total a este N, escala como a⁻³.
    # ρ_total_init lo NORMALIZO a 3 (es decir, "H=1 a N_init", no a hoy).
    # Al integrar forward H cambia, el código solo necesita estado.
    Om_r_init = rho_r_init / rho_total_init
    Om_r_init = min(Om_r_init, 0.01)
    Om_m_init = 1.0 - Om_phi_init - Om_r_init

    # φ-MDE: x = −α√6·KAL₀/3, y≈0
    # u = √6·x → u_init = -2·α·KAL₀
    # Mejor parametrizar directo: u_init pequeño negativo (rueda hacia +φ)
    rho_m_init = Om_m_init * rho_total_init
    rho_phi_init = Om_phi_init * rho_total_init
    # potencia exponencial V = V0·exp(−α·φ); inicialmente V chico
    V_init = 0.5 * rho_phi_init  # mitad cinética mitad potencial
    X_init = 0.5 * rho_phi_init  # tentativo
    phi_init = -np.log(V_init / V0_today) / alpha if V_init > 0 else 0.0
    u_init = -np.sqrt(2.0 * KAL0 * X_init)  # u<0: φ rueda hacia +φ (estándar)

    def V(phi):  return V0_today * np.exp(-alpha * phi)
    def dV(phi): return -alpha * V(phi)

    def H2_of(u, phi, rho_m, rho_r):
        denom = 3.0 - u*u/(2.0*KAL0)
        if denom <= 0: return 1e-30
        return max((V(phi) + rho_m + rho_r) / denom, 1e-30)

    def rhs(N, y):
        phi, u, rho_m = y
        a = np.exp(N)
        rho_r = rho_r0 * a**(-4) * (a_init**4 / 1.0)  # recalibrado
        # NB: tratamos rho_r como término residual; lo dominante es matter→DE
        rho_r = rho_r0 * a**(-4)
        H2 = H2_of(u, phi, rho_m, rho_r)
        # coupling conformal con saturación
        beta_c = alpha
        rhs_u = -3.0*u + KAL0*(-dV(phi)/H2 + beta_c*rho_m/H2)
        rhs_rho_m = -3.0*rho_m - beta_c*u*rho_m
        return [u, rhs_u, rhs_rho_m]

    try:
        sol = solve_ivp(rhs, (N_init, 0.0), [phi_init, u_init, rho_m_init],
                        t_eval=np.linspace(N_init, 0.0, 4000),
                        method='Radau', rtol=1e-9, atol=1e-12, max_step=0.01)
        if not sol.success:
            return None, sol.t[-1], "ODE falló"
        # Estado final
        phi_f, u_f, rho_m_f = sol.y[:, -1]
        a_f = 1.0
        rho_r_f = rho_r0
        H2_f = H2_of(u_f, phi_f, rho_m_f, rho_r_f)
        X_f = 0.5 * u_f*u_f * H2_f
        rho_phi_f = X_f/KAL0 + V(phi_f)
        p_phi_f = X_f/KAL0 - V(phi_f)
        w_phi_f = p_phi_f / rho_phi_f if abs(rho_phi_f) > 0 else 0
        rho_total_f = rho_m_f + rho_phi_f + rho_r_f
        Om_m_final = rho_m_f / rho_total_f
        Om_phi_final = rho_phi_f / rho_total_f
        return {'w_phi': w_phi_f, 'Om_m': Om_m_final, 'Om_phi': Om_phi_final,
                'N_end': sol.t[-1]}, sol.t[-1], None
    except Exception as e:
        return None, None, str(e)


print("="*78)
print("Mecanismo #7 — integración FORWARD desde φ-MDE hacia hoy")
print("="*78)
print(f"\n  α_sat (veta-2)           = {alpha_sat:.5f}")
print(f"  α_pot canonical          = {alpha_pot_canonical:.5f}")
print(f"  Objetivo hoy: w₀={w0:.3f}, Ω_m={Om_m_dyn:.3f} (dyn) o {3*Om_m_dyn:.3f} (CMB)")
print()
print(f"{'α':>8} {'w_φ(0)':>10} {'Ω_m(0)':>10} {'Ω_φ(0)':>10} {'N_end':>8} {'estado':>20}")
print("-"*78)

for alpha in [alpha_pot_canonical, alpha_sat, 0.8, 1.0, 1.5]:
    result, N_end, err = run_forward(alpha)
    if result is None:
        z_end = np.exp(-N_end)-1 if N_end else 0
        print(f"{alpha:8.4f}  FALLA (z_corte={z_end:.2f}) — {err[:30] if err else '?'}")
        continue
    arrival_ok = abs(result['w_phi']-w0)<0.1 and abs(result['Om_m']-3*Om_m_dyn)<0.05
    tag = "≈ hoy ✓" if arrival_ok else "no llega"
    print(f"{alpha:8.4f} {result['w_phi']:10.4f} {result['Om_m']:10.4f} "
          f"{result['Om_phi']:10.4f} {result['N_end']:8.2f} {tag:>20}")

print()
print("Interpretación:")
print("  Si algún α aterriza con w_φ≈−0.84 y Ω_m≈0.32 → lectura B viable")
print("  Si ninguno llega → islas dinámicas desconectadas; MIRA no es engranaje")
