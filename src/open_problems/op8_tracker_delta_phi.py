"""
OP-8 Ruta A — Δφ del tracker P7 entre CMB y today
==================================================

OBJETIVO:
  Integrar ecuaciones Klein-Gordon + Friedmann en SSEE V(φ) tracker
  (Paper 7), desde recombinación (z=1100) hasta hoy (z=0).

  Computar Δφ = φ(today) - φ(CMB).

  Comparar con predicción IDE:
      Δφ_IDE = ln(MIRA) / ⟨β⟩
  donde ⟨β⟩ es coupling promedio del Lagrangiano IDE candidato.

  Test: si Δφ_tracker ≈ Δφ_IDE → Lagrangiano consistente.
        si discrepa → interpretación IDE descartada.

ECUACIONES (unidades M_Pl=1, H₀=1):
  3 H² = ρ_m + ρ_φ
  ρ_m = 3·Ω_m,today / a³
  ρ_φ = ½ H² (dφ/dN)² + V(φ)
  V(φ) = V₀ · exp(-α·φ)   con α = √(3(M_v−T_r)/(M_v·KAL))

  Klein-Gordon en e-folds N=ln(a):
    φ'' + (3 + H'/H)·φ' - α·V/H² = 0

  H'/H = -½ [3·ρ_m/(3H²) + (3+w_φ)·ρ_φ/(3H²)]
       (en este script usamos forma directa: differentiar 3H²)

PARÁMETROS SSEE (algebraicos):
  Ω_m,today = T_r/M_v · ... = 0.160050  (Ω_m_dyn)
  Ω_DE,today = 0.840
  α = 0.2948
  λ (canonical) = √(3·Ω_m,dyn) = 0.6929
"""
import numpy as np
import sys
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
from ssee_core import PHI, PI, KAL0, OMEGA_M_DYN, OMEGA_DE, MIRA, AURA, BETA
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ─── Parámetros SSEE ────────────────────────────────────────────────
T_r = 3 * (PHI + BETA)            # = 11.9935
M_v = PHI + PI + (PHI + PI + 4.7596)  # recompute desde algebraico
# Mejor: cargar K_v desde ssee_core
from ssee_core import M_V, T_R
T_r = T_R
M_v = M_V

alpha_P7 = np.sqrt(3*(M_v - T_r) / (M_v * KAL0))
lambda_P7 = np.sqrt(3 * OMEGA_M_DYN)

print("=" * 78)
print(" OP-8 Ruta A — Tracker P7 Δφ integration")
print("=" * 78)
print(f"\n Parámetros tracker:")
print(f"   T_r            = {T_r:.6f}")
print(f"   M_v            = {M_v:.6f}")
print(f"   KAL₀           = {KAL0:.6f}")
print(f"   α (φ-variable) = {alpha_P7:.6f}")
print(f"   λ (χ canonical)= {lambda_P7:.6f}")
print(f"   Ω_m,today      = {OMEGA_M_DYN:.6f}")
print(f"   Ω_DE,today     = {OMEGA_DE:.6f}")

# ─── Ecuaciones de movimiento ───────────────────────────────────────
# Variables de estado: y = [φ, φ']  donde ' = d/dN, N = ln(a)
# Trabajo en unidades M_Pl = 1, normalizamos H₀ = 1.
#
# ρ_m,0 = 3·Ω_m,today (en unidades 3H₀²M_Pl² = ρ_crit_today)
# V₀ se ajustará por shooting para que Ω_φ(N=0) = Ω_DE

Omega_m_today = OMEGA_M_DYN

def H_squared(N, phi, dphi_dN, V0):
    """H²/H₀² dado N, φ, φ', V₀."""
    a = np.exp(N)
    rho_m = 3.0 * Omega_m_today / a**3      # ρ_m / ρ_crit_today
    V = V0 * np.exp(-alpha_P7 * phi)
    # 3 H² = ρ_m + ρ_φ
    # ρ_φ = ½ H² φ'² + V
    # 3 H² = ρ_m + ½ H² φ'² + V
    # H² (3 - ½ φ'²) = ρ_m + V
    H2 = (rho_m + V) / (3.0 - 0.5 * dphi_dN**2)
    return H2

def eom(N, y, V0):
    """Klein-Gordon + Friedmann."""
    phi, dphi = y
    a = np.exp(N)
    rho_m = 3.0 * Omega_m_today / a**3
    V = V0 * np.exp(-alpha_P7 * phi)
    dV_dphi = -alpha_P7 * V

    H2 = (rho_m + V) / (3.0 - 0.5 * dphi**2)

    # d(3H²)/dN = -3 ρ_m  (matter only) - 3(1+w_φ) ρ_φ
    # Usamos forma directa: H'/H = (dH²/dN) / (2 H²)
    # KG: φ'' + (3 + H'/H) φ' + V_φ/H² = 0

    # Calculamos H'/H usando ecuación de continuidad:
    # dρ_m/dN = -3 ρ_m
    # dρ_φ/dN = -3(ρ_φ + p_φ) = -3 H² φ'²
    # d(3H²)/dN = dρ_m/dN + dρ_φ/dN = -3ρ_m - 3 H² φ'²
    # => H'/H = (-3 ρ_m - 3 H² φ'²) / (6 H²) = -½(ρ_m/H² + φ'²)
    Hprime_over_H = -0.5 * (rho_m / H2 + dphi**2)

    ddphi = -(3.0 + Hprime_over_H) * dphi - dV_dphi / H2

    return [dphi, ddphi]

def Omega_phi_today(V0, N_CMB=-np.log(1101), phi_CMB=0.0, dphi_CMB=0.0):
    """Integra desde CMB hasta hoy y devuelve Ω_φ(N=0)."""
    sol = solve_ivp(eom, [N_CMB, 0.0], [phi_CMB, dphi_CMB], args=(V0,),
                    method='Radau', rtol=1e-10, atol=1e-12, max_step=0.01,
                    dense_output=True)
    if not sol.success:
        return None
    phi_f = sol.y[0, -1]
    dphi_f = sol.y[1, -1]
    V_f = V0 * np.exp(-alpha_P7 * phi_f)
    H2_f = (3.0 * Omega_m_today + V_f) / (3.0 - 0.5 * dphi_f**2)
    rho_phi = 0.5 * H2_f * dphi_f**2 + V_f
    return rho_phi / (3.0 * H2_f), sol

# ─── Shooting para fijar V₀ tal que Ω_φ(0) = 0.840 ──────────────────
print(f"\n Shooting V₀ para Ω_φ(today) = {OMEGA_DE:.4f}")

def shoot(V0):
    result = Omega_phi_today(V0)
    if result is None:
        return 1e6
    Om_phi, _ = result
    return Om_phi - OMEGA_DE

# bracket
V0_lo, V0_hi = 0.1, 5.0
try:
    V0_solution = brentq(shoot, V0_lo, V0_hi, xtol=1e-8)
    print(f"   V₀ encontrado = {V0_solution:.6f}")
    print(f"   (compare con 3·Ω_DE = {3*OMEGA_DE:.4f} si campo fuera estático)")
except Exception as e:
    print(f"   ERROR en shooting: {e}")
    sys.exit(1)

# ─── Resultado: Δφ ──────────────────────────────────────────────────
Om_phi_check, sol = Omega_phi_today(V0_solution)
print(f"\n Verificación: Ω_φ(today) = {Om_phi_check:.6f}")
print(f"               (objetivo:   {OMEGA_DE:.6f})")

N_CMB = -np.log(1101)
phi_CMB = sol.y[0, 0]
phi_today = sol.y[0, -1]
delta_phi = phi_today - phi_CMB

print(f"\n{'='*78}")
print(" RESULTADO PRINCIPAL")
print(f"{'='*78}")
print(f"   φ(z=1100)    = {phi_CMB:.6f}   [M_Pl]")
print(f"   φ(z=0)       = {phi_today:.6f}   [M_Pl]")
print(f"   Δφ           = {delta_phi:.6f}   [M_Pl]")

# ─── Comparación con predicción IDE ─────────────────────────────────
print(f"\n{'='*78}")
print(" COMPARACIÓN con predicción IDE candidata")
print(f"{'='*78}")

# Predicción IDE: Δφ · ⟨β⟩ = ln(MIRA)
ln_MIRA = np.log(MIRA)
print(f"\n Predicción IDE: Δφ · ⟨β⟩ = ln(MIRA) = {ln_MIRA:.6f}")
print(f"")
print(f" Para distintos candidatos de ⟨β⟩, Δφ requerido:")
print(f"")
candidates = {
    "BIAL":       BETA,
    "AURA":       AURA,
    "MIRA":       MIRA,
    "AURA/2=MIRA": AURA/2,
    "(AURA+BIAL)/2": (AURA+BETA)/2,
    "λ_P7":       lambda_P7,
    "α_P7":       alpha_P7,
}
print(f"   {'⟨β⟩ candidato':<22} {'valor':>10} {'Δφ_IDE_req':>14} {'Δφ_tracker':>14} {'match':>8}")
print(f"   {'-'*74}")
for name, beta_val in candidates.items():
    dphi_required = ln_MIRA / beta_val
    rel_match = abs(dphi_required - delta_phi) / max(abs(delta_phi), 1e-10) * 100
    flag = "★" if rel_match < 5 else ("✓" if rel_match < 20 else "")
    print(f"   {name:<22} {beta_val:>10.4f} {dphi_required:>14.6f} {delta_phi:>14.6f} {rel_match:>7.1f}% {flag}")

# ─── Conclusión ─────────────────────────────────────────────────────
print(f"\n{'='*78}")
print(" INTERPRETACIÓN")
print(f"{'='*78}")
print(f"""
 Δφ_tracker = {delta_phi:.4f} M_Pl entre z=1100 y z=0.

 Esto fija una constraint sobre el coupling IDE candidato:
     ⟨β⟩_required = ln(MIRA) / Δφ = {ln_MIRA/delta_phi:.4f}

 Si ⟨β⟩_required coincide con un algebraico SSEE natural (en la tabla
 de arriba con 'match' < 5%), el Lagrangiano IDE es consistente y
 OP-8 tiene candidato físico real.

 Si NO coincide: la interpretación IDE-via-tracker no funciona con
 P7 estándar. Hay que (a) modificar el potencial, (b) cambiar mecanismo,
 (c) cerrar OP-8 dejando MIRA como postulado auxiliar.
""")
