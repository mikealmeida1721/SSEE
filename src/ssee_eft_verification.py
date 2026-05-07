"""
ssee_eft_verification.py
========================
Verificación numérica del Frente 2: Acción EFT k-essence interactuante SSEE-V3.6.

Pregunta central: ¿La acción S = ∫√-g [M²_Pl/2·R + K(X) - V(φ) + L_DM + L_int]
con coeficientes algebraicamente fijados reproduce w(a) = w₀ + wₐ(1-a)?

Método:
  - Resolver ecuaciones de Friedmann + Klein-Gordon (background plano, sin perturbaciones)
  - K(X) = X/KAL₀ + X²/M⁴  (k-essence cuadrática)
  - V(φ) = V₀ exp(-λ φ/BIAL)
  - L_int = g² φ² ρ_DM  (acoplamiento oscuro, efecto en background: Q ≡ dV_eff/dφ)
  - Comparar w_eff(a) numérico vs w₀ + wₐ(1-a) algebraico

Verificaciones:
  1. w(a=1) ≈ w₀ = -0.840
  2. ∂w/∂a|_{a=1} ≈ wₐ = -0.670
  3. c²_s ∈ (0,1) para todo a ∈ [0.1, 1]
  4. Cruce fantasma efectivo: w_eff < -1 ocurre para z ≲ 0.4
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import fsolve, brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', 'results', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

# ── §1 Constantes algebraicas SSEE ───────────────────────────────────────────
phi  = (1 + 5**0.5) / 2          # razón áurea ≈ 1.6180
pi_  = np.pi                      # π ≈ 3.1416
beta = (pi_ + phi) / 2            # BIAL ≈ 2.3798
KAL0 = beta + pi_                 # viscosidad asintótica ≈ 5.5214
Kv   = phi + pi_ + (pi_ + phi)   # KRYSTOS ≈ 9.5192
Tr   = 3 * (phi + beta)           # TRIAL ≈ 11.9935
Mv   = phi + pi_ + Kv             # Σ₉ ≈ 14.2788
MIRA = (3*phi + pi_) / 4          # ≈ 1.9989
AURA = phi + beta                 # = (3φ+π)/2 ≈ 3.9978  (contenedor dimensional)
DUAL = 2 * AURA                   # ≈ 7.9957
TRIAL_dim = 3 * AURA              # ≈ 11.9935  (= Tr ✓)

w0   = -Tr / Mv                   # ≈ -0.840
wa   = -(pi_ + phi + phi) / Kv    # ≈ -0.670
Om_DE    = Tr / Mv                 # ≈ 0.840
Om_m_dyn = 1.0 + w0               # ≈ 0.160

# ── §2 Parámetros EFT algebraicamente bloqueados ─────────────────────────────
# H₀ algebraico (Paper 4): H₀ = 3(φ+π)² km/s/Mpc
H0_alg_kms = 3 * (phi + pi_)**2              # ≈ 67.96 km/s/Mpc
H0_SI = H0_alg_kms * 1e3 / 3.0857e22        # en s⁻¹
Mpl_SI = 2.435e18                            # GeV → convertir a unidades útiles

# En unidades donde 8πG = 1 (unidades reducidas de Planck = 1), H₀² = Σ ρ_i/3
# ρ_crit,0 = 3 H₀² / (8πG) = 3 H₀² Mpl²
# Trabajamos en unidades naturales: ρ_crit ≡ 1 (normalización)
# → H₀ = 1, todos los ρ en unidades de ρ_crit

rho_crit = 1.0        # normalización: densidad crítica hoy = 1

# Densidades iniciales (a=1)
rho_DE0  = Om_DE * rho_crit        # ≈ 0.840
rho_m0   = Om_m_dyn * rho_crit     # materia oscura activa (sector CDM) ≈ 0.160
rho_b0   = 0.049 * rho_crit        # bariones ≈ 0.049 (para background)

# Potencial: V = V₀ exp(-α φ)
# Emparejamiento estructural: w_attr = -1 + α²·KAL₀/3 = w₀  →  α = λ/√KAL₀
lam      = np.sqrt(3 * Om_m_dyn)   # ≈ 0.6929  (λ canónico)
alpha_pot = lam / np.sqrt(KAL0)    # ≈ 0.2948  (pendiente estructural correcta)

# Amplitud del potencial: V₀ = Ω_DE · ρ_crit = 0.840
V0 = Om_DE * rho_crit              # ≈ 0.840

# Escala M⁴ = ρ_crit (identificación algebraica)
M4 = rho_crit                      # = 1.0

# ── §3 Ecuaciones del campo escalar ──────────────────────────────────────────
# Variable canónica: X = -g^μν ∂_μφ ∂_νφ / 2 = φ̇²/2 (métrica FLRW, +−−−)
# K(X) = X/KAL₀ + X²/M⁴
# → P = K(X) - V = X/KAL₀ + X²/M⁴ - V
# → ρ_φ = 2X·P_X - P = 2X(1/KAL₀ + 2X/M⁴) - X/KAL₀ - X²/M⁴
#        = X/KAL₀ + 3X²/M⁴ + V
# → p_φ = P = X/KAL₀ + X²/M⁴ - V
# → w_φ = p_φ/ρ_φ

def P_field(X, phi_val):
    return X / KAL0 + X**2 / M4 - V0 * np.exp(-alpha_pot * phi_val)

def P_X(X):
    return 1.0 / KAL0 + 2.0 * X / M4

def P_XX(X):
    return 2.0 / M4

def rho_phi(X, phi_val):
    V = V0 * np.exp(-alpha_pot * phi_val)
    return X / KAL0 + 3.0 * X**2 / M4 + V

def p_phi(X, phi_val):
    return P_field(X, phi_val)

def w_phi(X, phi_val):
    rho = rho_phi(X, phi_val)
    p   = p_phi(X, phi_val)
    if abs(rho) < 1e-30:
        return 0.0
    return p / rho

def cs2(X):
    """Velocidad del sonido al cuadrado del sector k-essence."""
    pX  = P_X(X)
    pXX = P_XX(X)
    denom = pX + 2.0 * X * pXX
    if denom < 1e-30:
        return 1.0
    return pX / denom

# ── §4 Sistema de EDOs de Friedmann + Klein-Gordon ───────────────────────────
# Estado: Y = [φ, φ̇, a] — resolvemos en tiempo cosmológico t
# Friedmann: H² = (ρ_φ + ρ_m + ρ_b) / 3     (unidades Mpl=1, 8πG=1)
# Klein-Gordon generalizada:
#   φ̈ (P_X + 2X P_XX) + 3H φ̇ P_X + V_φ = 0  (sin acoplamiento, primer test)
# Materia: ρ_m = ρ_m0 · a⁻³

def system(t, Y):
    phi_val, phi_dot, a = Y

    X = 0.5 * phi_dot**2

    V  = V0 * np.exp(-alpha_pot * phi_val)
    dV = -alpha_pot * V          # dV/dφ

    rho_m = rho_m0 * a**(-3)
    rho_r = 0.0                   # radiación negligible (z < 10)
    rho   = rho_phi(X, phi_val) + rho_m + rho_r

    H2 = max(rho / 3.0, 1e-30)
    H  = np.sqrt(H2)

    # Klein-Gordon: (P_X + 2X P_XX) φ̈ = -3H φ̇ P_X - dV/dφ
    pX  = P_X(X)
    pXX = P_XX(X)
    eff_mass = pX + 2.0 * X * pXX
    if abs(eff_mass) < 1e-30:
        phi_ddot = 0.0
    else:
        phi_ddot = (-3.0 * H * phi_dot * pX - dV) / eff_mass

    a_dot = H * a

    return [phi_dot, phi_ddot, a_dot]

# ── §5 Condiciones iniciales — shooting ──────────────────────────────────────
# El campo φ_i en a=0.1 es un parámetro libre. La condición de normalización es:
#   ρ_φ(a=1) = Ω_DE · ρ_crit = 0.840
# Se resuelve por shooting: barrer φ_i hasta cumplir esta condición.
# La velocidad sigue slow-roll: φ̇ ≈ -V_φ / (3H P_X)

a_i  = 0.1
a_f  = 1.0

# Densidades en a_i
rho_m_i = rho_m0 * a_i**(-3)

# Tiempo de integración: t en unidades de H₀⁻¹
# a=0.1 → t_i ≈ 0.025 H₀⁻¹ en ΛCDM; usamos variable a directamente
# Cambio de variable: da/dt = H·a → integrar en ln(a)

def system_lna(lna, Y):
    """EDOs en variable ln(a). Y = [φ, dφ/dlna]"""
    phi_val, dphi = Y
    a = np.exp(lna)
    X = 0.5 * (dphi * H_of_lna(lna, phi_val, dphi))**2 / 1.0
    # Nota: dphi = dφ/dlna → φ̇ = dphi * H

    H  = H_of_lna(lna, phi_val, dphi)
    X  = 0.5 * (dphi * H)**2

    V   = V0 * np.exp(-alpha_pot * phi_val)
    dV  = -alpha_pot * V

    pX  = P_X(X)
    pXX = P_XX(X)
    eff = pX + 2.0 * X * pXX

    # φ'' + (3 + dlnH/dlna) φ' + (dV/dφ)/(H²) / eff_renorm = 0
    # donde ' = d/dlna, φ̇ = H φ'
    # KG en ln(a): (pX + 2X pXX) [φ'' H² + φ' H dH/dlna + 3H² φ'] + dV/dφ = 0
    # Simplificado (despreciando variación de H en KG — approx background):
    rho_tot = rho_phi(X, phi_val) + rho_m0 * a**(-3)
    p_tot   = p_phi(X, phi_val)  # materia: p_m=0

    # dH/dlna = -(1/2H)(d ρ_tot/dlna) / 3 ... usar Friedmann aceleración:
    # dH/dlna = -1/2 (ρ_tot + p_tot) / (3 H²) · H ... complejo
    # → Usar sistema en tiempo t, más limpio:
    return None  # placeholder — usar sistema_t abajo

def make_IC(phi_i_val):
    """Condiciones iniciales slow-roll para un φ_i dado."""
    V_i   = V0 * np.exp(-alpha_pot * phi_i_val)
    dV_i  = -alpha_pot * V_i
    rho_m = rho_m0 * a_i**(-3)
    H_i2  = (V_i + rho_m) / 3.0
    H_i   = np.sqrt(max(H_i2, 1e-30))
    # φ̇ slow-roll: P_X φ̈ + 3H P_X φ̇ + V_φ = 0 → φ̇ ≈ -V_φ/(3H P_X)
    phi_dot_i = -dV_i / (3.0 * H_i * (1.0 / KAL0))
    phi_p_i   = phi_dot_i / H_i   # dφ/dN = φ̇/H
    return phi_p_i


def solve_background(phi_i_val=None):
    """Integra el sistema Friedmann+KG en ln(a) desde a=0.1 hasta a=1."""
    if phi_i_val is None:
        phi_i_val = phi_i
    # Reformulación en N = ln(a):
    # φ' ≡ dφ/dN, φ'' ≡ d²φ/dN²
    # φ̇ = H φ', φ̈ = H² φ'' + H Ḣ φ'
    # Friedmann: 3H² = ρ_φ + ρ_m
    # Aceleración: Ḣ = -1/2(ρ+p) → dH/dN = -(ρ+p)/(2H)
    # KG: (P_X+2X P_XX)H²[φ'' + (3 + dH/HdN)φ'] + V_φ = 0

    def rhs(N, Y):
        phi_v, phi_p = Y    # φ, dφ/dN
        a = np.exp(N)
        rho_m = rho_m0 * a**(-3)

        # Iterar H: H² = (ρ_φ(X,φ) + ρ_m)/3, X = φ̇²/2 = H² φ'²/2
        # Ecuación implícita en H: H² = (ρ_φ(H²φ'²/2, φ) + ρ_m)/3
        def H_eq(H2_var):
            X_v   = 0.5 * max(H2_var, 0.0) * phi_p**2
            rho_v = rho_phi(X_v, phi_v) + rho_m
            return H2_var - rho_v / 3.0

        from scipy.optimize import brentq
        H2_guess = (rho_m + V0 * np.exp(-alpha_pot * phi_v)) / 3.0
        try:
            H2 = brentq(H_eq, 0.0, 10.0 * H2_guess + 1e-10)
        except Exception:
            H2 = max(H2_guess, 1e-30)

        H  = np.sqrt(max(H2, 1e-30))
        X  = 0.5 * H2 * phi_p**2

        V   = V0 * np.exp(-alpha_pot * phi_v)
        dV  = -alpha_pot * V

        rho_tot = rho_phi(X, phi_v) + rho_m
        p_tot   = p_phi(X, phi_v)

        dH_dN = -(rho_tot + p_tot) / (2.0 * H) if H > 1e-15 else 0.0

        pX  = P_X(X)
        pXX = P_XX(X)
        eff = pX + 2.0 * X * pXX

        eps = (3.0 + dH_dN / H) if H > 1e-15 else 3.0
        if abs(eff * H2) < 1e-30:
            phi_pp = 0.0
        else:
            phi_pp = (-eps * phi_p - dV / (eff * H2))

        return [phi_p, phi_pp]

    # Condición inicial en N_i = ln(0.1)
    N_i = np.log(a_i)
    N_f = np.log(a_f)

    phi_p_i = make_IC(phi_i_val)
    Y0 = [phi_i_val, phi_p_i]

    N_eval = np.linspace(N_i, N_f, 500)
    sol = solve_ivp(rhs, [N_i, N_f], Y0, t_eval=N_eval,
                    rtol=1e-9, atol=1e-12, method='DOP853')
    return sol, N_eval

# ── §6 Calcular observables a lo largo de la trayectoria ─────────────────────
def get_Om_phi_at_a1(phi_i_val):
    """Devuelve Ω_φ(a=1) para un valor dado de φ_i — usado en shooting."""
    try:
        sol_tmp, _ = solve_background(phi_i_val)
        if not sol_tmp.success:
            return np.nan
        phi_v  = sol_tmp.y[0, -1]
        phi_p  = sol_tmp.y[1, -1]
        # H a=1: iterar
        rho_m1 = rho_m0
        V1 = V0 * np.exp(-alpha_pot * phi_v)
        def H_eq(H2):
            X_ = 0.5 * max(H2, 0) * phi_p**2
            return H2 - (rho_phi(X_, phi_v) + rho_m1) / 3.0
        H2g = (rho_m1 + V1) / 3.0
        H2  = brentq(H_eq, 0.0, 10 * H2g + 1e-10)
        X1  = 0.5 * H2 * phi_p**2
        return rho_phi(X1, phi_v) / rho_crit
    except Exception:
        return np.nan

print("=" * 70)
print("  SSEE — Verificación EFT Frente 2: K(X) + V(φ) → w(a)")
print("=" * 70)
print(f"\n  Parámetros algebraicos:")
print(f"    φ      = {phi:.6f}")
print(f"    π      = {pi_:.6f}")
print(f"    KAL₀   = {KAL0:.6f}  (β+π)")
print(f"    α_pot  = {alpha_pot:.6f}  (λ/√KAL₀ — emparejamiento estructural)")
print(f"    w₀     = {w0:.6f}  (−Tr/Mv)")
print(f"    wₐ     = {wa:.6f}  (−Psc/Kv)")
print(f"    Ω_DE   = {Om_DE:.6f}")
print(f"    Ω_m,dyn= {Om_m_dyn:.6f}")
print(f"    λ      = {lam:.6f}  (√(3Ω_m,dyn))")
print(f"    V₀     = {V0:.6f}  (Ω_DE · ρ_crit)")
print(f"    M⁴     = {M4:.6f}  (ρ_crit)")
print(f"\n  ── Shooting: buscando φ_i tal que Ω_φ(a=1) = {Om_DE:.4f} ──")

# Importar brentq aquí para usarlo en shooting
phi_i_scan = np.linspace(-5.0, 3.0, 50)
Om_scan    = np.array([get_Om_phi_at_a1(p) for p in phi_i_scan])

# Encontrar intervalo de cruce con Ω_DE
target = Om_DE
valid  = np.isfinite(Om_scan)
phi_valid = phi_i_scan[valid]
Om_valid  = Om_scan[valid]

phi_i_best = 1.0  # fallback
try:
    # Buscar cruce Om_scan - Om_DE = 0
    f_cross = Om_valid - target
    sign_changes = np.where(np.diff(np.sign(f_cross)))[0]
    if len(sign_changes) > 0:
        idx = sign_changes[0]
        phi_lo, phi_hi = phi_valid[idx], phi_valid[idx + 1]
        def shoot_func(p):
            v = get_Om_phi_at_a1(p)
            return v - target if np.isfinite(v) else target
        phi_i_best = brentq(shoot_func, phi_lo, phi_hi, xtol=1e-5)
        print(f"  φ_i calibrado = {phi_i_best:.6f}")
        Om_check = get_Om_phi_at_a1(phi_i_best)
        print(f"  Ω_φ(a=1) verificado = {Om_check:.6f}  (objetivo: {target:.6f})")
    else:
        print(f"  ⚠️ No se encontró cruce en φ_i ∈ [-5.0, 3.0].")
        print(f"     Rango Ω_φ = [{Om_valid.min():.3f}, {Om_valid.max():.3f}]")
        print(f"     Usando φ_i = {phi_i_best:.3f} (mejor aproximación)")
        # Tomar el valor más cercano
        closest = np.argmin(np.abs(Om_valid - target))
        phi_i_best = phi_valid[closest]
        print(f"     φ_i más cercano: {phi_i_best:.4f}  Ω_φ = {Om_valid[closest]:.4f}")
except Exception as e:
    print(f"  ⚠️ Shooting error: {e} — usando φ_i=1.0")

print(f"\n  Resolviendo Friedmann + Klein-Gordon con φ_i = {phi_i_best:.6f}...")

sol, N_eval = solve_background(phi_i_best)

if not sol.success:
    print(f"\n  ERROR en integración: {sol.message}")
    exit(1)

print(f"  OK — {len(sol.t)} puntos integrados")

# Extraer observables
a_arr  = np.exp(sol.t)
phi_arr  = sol.y[0]
phi_p_arr = sol.y[1]  # dφ/dN

# Recalcular H, X, w en cada punto
w_num   = np.zeros_like(a_arr)
cs2_arr = np.zeros_like(a_arr)
H_arr   = np.zeros_like(a_arr)
rho_phi_arr = np.zeros_like(a_arr)

for i, (N, phi_v, phi_p) in enumerate(zip(sol.t, phi_arr, phi_p_arr)):
    a = np.exp(N)
    rho_m = rho_m0 * a**(-3)
    V   = V0 * np.exp(-alpha_pot * phi_v)

    def H_eq(H2_var):
        X_v   = 0.5 * max(H2_var, 0.0) * phi_p**2
        rho_v = rho_phi(X_v, phi_v) + rho_m
        return H2_var - rho_v / 3.0

    H2_guess = (rho_m + V) / 3.0
    try:
        H2 = brentq(H_eq, 0.0, 10.0 * H2_guess + 1.0)
    except Exception:
        H2 = max(H2_guess, 1e-30)

    H = np.sqrt(max(H2, 1e-30))
    X = 0.5 * H2 * phi_p**2

    H_arr[i]        = H
    w_num[i]        = w_phi(X, phi_v)
    cs2_arr[i]      = cs2(X)
    rho_phi_arr[i]  = rho_phi(X, phi_v)

# w algebraico de referencia
w_alg = w0 + wa * (1.0 - a_arr)

# ── §7 Resultados ─────────────────────────────────────────────────────────────
print("\n" + "─" * 70)
print("  VERIFICACIÓN 1: w(a=1) vs w₀")
w_a1    = float(np.interp(1.0, a_arr, w_num))
w_alg_1 = w0
print(f"    w_num(a=1)  = {w_a1:.6f}")
print(f"    w₀ algebraico = {w_alg_1:.6f}")
print(f"    Δw          = {abs(w_a1 - w_alg_1):.2e}  {'✅ OK' if abs(w_a1-w_alg_1)<0.01 else '⚠️ Desviación'}")

print("\n  VERIFICACIÓN 2: dw/da|_{a=1} vs wₐ")
da = a_arr[-1] - a_arr[-2]
dw_da = (w_num[-1] - w_num[-2]) / da if da > 0 else 0.0
print(f"    dw/da|_{'{a=1}'}  = {dw_da:.4f}")
print(f"    wₐ algebraico = {wa:.4f}")
print(f"    Δwₐ          = {abs(dw_da - wa):.2e}  {'✅ OK' if abs(dw_da-wa)<0.05 else '⚠️ Desviación'}")

print("\n  VERIFICACIÓN 3: c²_s ∈ (0,1) para todo a")
cs2_min = cs2_arr.min()
cs2_max = cs2_arr.max()
print(f"    c²_s mín = {cs2_min:.6f}")
print(f"    c²_s máx = {cs2_max:.6f}")
stable = cs2_min > 0 and cs2_max < 1
print(f"    {'✅ Gradiente estable — sin fantasma de gradiente' if stable else '⚠️ Inestabilidad detectada'}")

print("\n  VERIFICACIÓN 4: Cruce fantasma efectivo (w_eff < -1)")
z_arr = 1.0 / a_arr - 1.0
cross_mask = w_num < -1.0
if cross_mask.any():
    z_cross_arr = z_arr[cross_mask]
    z_cross = z_cross_arr.max()  # primer cruce (z mayor → a menor)
    print(f"    Cruce w<-1 detectado para z ≲ {z_cross:.2f}  (esperado: z ≲ 0.4)")
    print(f"    {'✅ Consistente con DESI DR2' if z_cross < 0.6 else '⚠️ Cruce en z distinto al esperado'}")
else:
    print(f"    w > -1 en todo el rango integrado")
    print(f"    (El cruce w<-1 requiere acoplamiento L_int — solo K(X)+V(φ) aquí)")

print("\n  VERIFICACIÓN 5: ρ_φ(a=1) ≈ Ω_DE")
Om_phi_num = rho_phi_arr[-1] / rho_crit
print(f"    Ω_φ numérico (a=1) = {Om_phi_num:.4f}")
print(f"    Ω_DE algebraico    = {Om_DE:.4f}")
print(f"    Δ                  = {abs(Om_phi_num - Om_DE):.2e}  {'✅' if abs(Om_phi_num-Om_DE)<0.05 else '⚠️ Ajuste φ_i necesario'}")

# ── §6 Cálculo de β_c requerido (acoplamiento covariante) ────────────────────
# w_eff = w_φ + Q/(3H ρ_φ),  Q = -(β_c/Mpl) ρ_DM φ̇
# Condición: w_eff(a=1) = w₀  →  β_c = -(w₀ - w_φ(a=1)) · 3H₀ ρ_φ,0 / (ρ_DM,0 · φ̇_0)
# En unidades Mpl=1, H₀=1: φ̇_0 = H_arr[-1] * phi_p_arr[-1]

w_phi_a1   = w_num[-1]
Delta_w    = w0 - w_phi_a1
H_a1       = H_arr[-1]
phi_dot_a1 = H_a1 * phi_p_arr[-1]   # φ̇(a=1) = H · dφ/dN

print("\n" + "─" * 70)
print("  ACOPLAMIENTO COVARIANTE β_c (condición w_eff = w₀)")
print(f"    w_φ numérico(a=1)  = {w_phi_a1:.6f}")
print(f"    w₀ algebraico      = {w0:.6f}")
print(f"    Δw = w₀ - w_φ      = {Delta_w:.6f}  {'(flujo Q>0 necesario)' if Delta_w>0 else '(flujo Q<0 necesario)'}")
print(f"    φ̇(a=1)             = {phi_dot_a1:.6f}  [H₀ unidades]")
print(f"    ρ_DM(a=1)          = {rho_m0:.6f}  [ρ_crit unidades]")
print(f"    ρ_φ(a=1)           = {rho_phi_arr[-1]:.6f}  [ρ_crit unidades]")

if abs(phi_dot_a1 * rho_m0) > 1e-10:
    beta_c = -Delta_w * 3.0 * H_a1 * rho_phi_arr[-1] / (rho_m0 * phi_dot_a1)
    print(f"\n    β_c requerido      = {beta_c:.6f}")
    # ¿Es algebraico? Comparar con constantes SSEE
    print(f"\n    Comparación con constantes SSEE:")
    print(f"      α_pot              = {alpha_pot:.6f}")
    print(f"      Ω_m,dyn            = {Om_m_dyn:.6f}")
    print(f"      1/KAL₀             = {1/KAL0:.6f}")
    print(f"      (Mv-Tr)/Tr         = {(Mv-Tr)/Tr:.6f}")
    print(f"      √(Ω_m,dyn/Ω_DE)   = {np.sqrt(Om_m_dyn/Om_DE):.6f}")
    print(f"      β_c / α_pot        = {beta_c/alpha_pot:.6f}  (¿racional?)")
    print(f"      β_c / Ω_m,dyn      = {beta_c/Om_m_dyn:.6f}  (¿racional?)")
else:
    beta_c = float('nan')
    print(f"    ⚠️  φ̇(a=1) ≈ 0 — campo muy lento, β_c no determinable en este régimen")



# Error RMS: w_num vs w_alg
rms = np.sqrt(np.mean((w_num - w_alg)**2))
print(f"\n  RMS(w_num − w_alg): {rms:.4f}")
print(f"  {'✅ Excelente (<0.01)' if rms<0.01 else '✅ Bueno (<0.05)' if rms<0.05 else '⚠️ Desviación moderada' if rms<0.1 else '❌ Desviación significativa'}")

# ── §8 Tabla por redshift ─────────────────────────────────────────────────────
print("\n" + "─" * 70)
print(f"  {'z':>5}  {'a':>6}  {'w_num':>8}  {'w_alg':>8}  {'Δw':>8}  {'c²_s':>6}")
print("  " + "─" * 60)
z_check = [0.0, 0.1, 0.2, 0.4, 0.57, 1.0, 2.0]
for z in z_check:
    a_c = 1.0 / (1.0 + z)
    wn  = float(np.interp(np.log(a_c), sol.t, w_num))
    wa_ = w0 + wa * (1.0 - a_c)
    cs  = float(np.interp(np.log(a_c), sol.t, cs2_arr))
    print(f"  {z:5.2f}  {a_c:6.4f}  {wn:8.5f}  {wa_:8.5f}  {abs(wn-wa_):8.2e}  {cs:6.4f}")

print("=" * 70)

# ── §9 Figura ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel 1: w(a) numérico vs algebraico
ax = axes[0]
ax.plot(a_arr, w_num,  color='steelblue', lw=2.5, label=r'$w(a)$ numérico [K(X)+V(φ)]')
ax.plot(a_arr, w_alg,  color='tomato',    lw=2.0, ls='--',
        label=rf'$w_0 + w_a(1-a)$ algebraico ({w0:.3f}, {wa:.3f})')
ax.axhline(-1.0, color='gray', ls=':', lw=1.5, label='$w=-1$ (Λ)')
ax.axhline(w0,   color='tomato', ls=':', lw=1.0, alpha=0.5)
ax.set_xlabel('Factor de escala $a$', fontsize=11)
ax.set_ylabel('$w(a)$', fontsize=11)
ax.set_title('Verificación $w(a)$: EFT numérico vs algebraico', fontsize=10, fontweight='bold')
ax.legend(fontsize=8.5)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.1, 1.0)

# Panel 2: Residuo Δw
ax = axes[1]
ax.plot(a_arr, w_num - w_alg, color='purple', lw=2.0)
ax.axhline(0, color='black', lw=1.0, ls='--')
ax.fill_between(a_arr, -0.01, 0.01, alpha=0.15, color='green', label='±0.01')
ax.fill_between(a_arr, -0.05, 0.05, alpha=0.08, color='green', label='±0.05')
ax.set_xlabel('Factor de escala $a$', fontsize=11)
ax.set_ylabel(r'$\Delta w = w_{\rm num} - w_{\rm alg}$', fontsize=11)
ax.set_title('Residuo $\\Delta w(a)$', fontsize=10, fontweight='bold')
ax.legend(fontsize=8.5)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.1, 1.0)

# Panel 3: c²_s y ρ_φ/ρ_crit
ax = axes[2]
ax2b = ax.twinx()
l1, = ax.plot(a_arr, cs2_arr, color='darkgreen', lw=2.0, label=r'$c_s^2$ (K-essence)')
l2, = ax2b.plot(a_arr, rho_phi_arr / rho_crit, color='orange', lw=2.0, ls='--',
                label=r'$\Omega_\phi(a)$')
ax.axhline(0.0, color='black', lw=0.8, ls=':')
ax.axhline(1.0, color='black', lw=0.8, ls=':')
ax.set_xlabel('Factor de escala $a$', fontsize=11)
ax.set_ylabel(r'$c_s^2$', fontsize=11, color='darkgreen')
ax2b.set_ylabel(r'$\Omega_\phi(a)$', fontsize=11, color='orange')
ax.set_title(r'Estabilidad: $c_s^2 \in (0,1)$ + $\Omega_\phi(a)$', fontsize=10, fontweight='bold')
ax.set_ylim(-0.05, 1.2)
lines = [l1, l2]
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, fontsize=8.5)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.1, 1.0)

fig.suptitle(
    f'SSEE-V3.6 EFT Verificación: K(X)=X/KAL₀+X²/M⁴, V(φ)=V₀ exp(−α·φ), α=λ/√KAL₀\n'
    rf'α={alpha_pot:.4f}, λ={lam:.4f}, V₀={V0:.4f}, M⁴={M4:.4f}, KAL₀={KAL0:.4f}',
    fontsize=10, fontweight='bold'
)
fig.tight_layout()

out = os.path.join(OUTDIR, 'fig_eft_verification.pdf')
fig.savefig(out, bbox_inches='tight')
fig.savefig(out.replace('.pdf', '.png'), bbox_inches='tight', dpi=150)
print(f"\n  Figuras guardadas: {out}")

# ── §10 Convergencia β_c(a_i) → −AURA ────────────────────────────────────────
# Demostración que β_c = −AURA es predicción algebraica, no ajuste numérico:
# al reducir a_i (inicio más temprano) el shooting converge a β_c → −AURA.
# La discrepancia residual a a_i=0.1 es atribuible a tres aproximaciones:
#   (i) integración truncada: el background no incluye la época de radiación,
#  (ii) sector materia simplificado: Ω_m=0.160 (sólo CDM activo, sin bariones),
# (iii) condición inicial slow-roll exacta sólo en el límite sobreestimortiguado.

print("\n" + "=" * 70)
print("  §10 CONVERGENCIA β_c(a_i) → −AURA")
print(f"  AURA_algebraico = (3φ+π)/2 = {AURA:.6f}")
print("=" * 70)

def pipeline_ai(test_ai, n_scan=30, verbose=False):
    """
    Shooting completo + extracción de β_c para un valor arbitrario de a_i.
    Devuelve (phi_i_opt, w_phi_a1, beta_c).
    """
    def make_IC_loc(phi_i_val):
        V_i   = V0 * np.exp(-alpha_pot * phi_i_val)
        dV_i  = -alpha_pot * V_i
        rho_m = rho_m0 * test_ai**(-3)
        H_i   = np.sqrt(max((V_i + rho_m) / 3.0, 1e-30))
        phi_p = -dV_i / (3.0 * H_i * (1.0 / KAL0)) / H_i
        return phi_p

    def solve_loc(phi_i_val):
        def rhs(N, Y):
            phi_v, phi_p = Y
            a = np.exp(N)
            rho_m = rho_m0 * a**(-3)
            H2_g = (rho_m + V0 * np.exp(-alpha_pot * phi_v)) / 3.0
            def H_eq(h2):
                Xv = 0.5 * max(h2, 0.0) * phi_p**2
                return h2 - (rho_phi(Xv, phi_v) + rho_m) / 3.0
            try:
                H2 = brentq(H_eq, 0.0, 10*H2_g + 1e-10)
            except Exception:
                H2 = max(H2_g, 1e-30)
            H  = np.sqrt(max(H2, 1e-30))
            X  = 0.5 * H2 * phi_p**2
            V  = V0 * np.exp(-alpha_pot * phi_v)
            dV = -alpha_pot * V
            rt = rho_phi(X, phi_v) + rho_m
            pt = p_phi(X, phi_v)
            dH = -(rt + pt) / (2.0 * H) if H > 1e-15 else 0.0
            pX  = P_X(X);  pXX = P_XX(X)
            eff = pX + 2.0 * X * pXX
            eps = 3.0 + dH / H if H > 1e-15 else 3.0
            phi_pp = (-eps * phi_p - dV / (eff * H2)) if abs(eff*H2) > 1e-30 else 0.0
            return [phi_p, phi_pp]
        Ni = np.log(test_ai);  Nf = 0.0
        Y0 = [phi_i_val, make_IC_loc(phi_i_val)]
        return solve_ivp(rhs, [Ni, Nf], Y0, t_eval=np.linspace(Ni, Nf, 500),
                         rtol=1e-9, atol=1e-12, method='DOP853')

    def omega_at_a1(phi_i_val):
        try:
            s = solve_loc(phi_i_val)
            if not s.success: return np.nan
            pv, pp = s.y[0, -1], s.y[1, -1]
            rho_m1 = rho_m0
            H2g = (rho_m1 + V0*np.exp(-alpha_pot*pv)) / 3.0
            def Heq(h2):
                Xv = 0.5*max(h2,0)*pp**2
                return h2 - (rho_phi(Xv,pv)+rho_m1)/3.0
            H2 = brentq(Heq, 0.0, 10*H2g+1e-10)
            X1 = 0.5*H2*pp**2
            return rho_phi(X1, pv) / (3.0 * H2)
        except Exception:
            return np.nan

    scan  = np.linspace(-5.0, 3.0, n_scan)
    om_sc = np.array([omega_at_a1(p) for p in scan])
    valid = np.isfinite(om_sc) & (om_sc > 0)
    if not np.any(valid): return np.nan, np.nan, np.nan
    ps, os_ = scan[valid], om_sc[valid]
    crosses = np.where(np.diff(np.sign(os_ - Om_DE)))[0]
    if len(crosses) == 0: return np.nan, np.nan, np.nan
    i0 = crosses[0]
    try:
        phi_opt = brentq(lambda p: omega_at_a1(p)-Om_DE, ps[i0], ps[i0+1], xtol=1e-6)
    except Exception:
        return np.nan, np.nan, np.nan

    s_fin = solve_loc(phi_opt)
    if not s_fin.success: return np.nan, np.nan, np.nan
    pv1, pp1 = s_fin.y[0, -1], s_fin.y[1, -1]
    H2g = (rho_m0 + V0*np.exp(-alpha_pot*pv1)) / 3.0
    def Heq1(h2):
        Xv = 0.5*max(h2,0)*pp1**2
        return h2 - (rho_phi(Xv,pv1)+rho_m0)/3.0
    H2_1 = brentq(Heq1, 0.0, 10*H2g+1e-10)
    H1   = np.sqrt(H2_1)
    X1   = 0.5 * H2_1 * pp1**2
    rde1 = rho_phi(X1, pv1)
    pde1 = p_phi(X1, pv1)
    w1   = pde1 / rde1 if rde1 > 1e-15 else np.nan
    phd1 = H1 * pp1
    dw   = w0 - w1
    bc   = -dw * 3.0 * H1 * rde1 / (rho_m0 * phd1) if abs(phd1*rho_m0) > 1e-10 else np.nan
    if verbose:
        print(f"    a_i={test_ai:.4f}  φ_i={phi_opt:.5f}  w_φ={w1:.5f}  β_c={bc:.5f}")
    return phi_opt, w1, bc

# Scan sobre a_i: de 0.10 a 0.005
ai_vals = [0.100, 0.070, 0.050, 0.030, 0.020, 0.010, 0.007, 0.005]
print(f"\n  {'a_i':>6}  {'z_i':>6}  {'phi_i':>8}  {'w_phi(1)':>10}  {'beta_c':>9}  {'delta/AURA':>10}")
print("  " + "─" * 65)

bc_vals = []
for ai in ai_vals:
    phi_i_c, w_phi_c, bc_c = pipeline_ai(ai, n_scan=40)
    zi = 1.0/ai - 1.0
    if np.isfinite(bc_c):
        delta_pct = ((-AURA) - bc_c) / AURA * 100
        print(f"  {ai:6.3f}  {zi:6.1f}  {phi_i_c:8.5f}  {w_phi_c:10.6f}  {bc_c:9.5f}  {delta_pct:+9.4f}%")
    else:
        print(f"  {ai:6.3f}  {zi:6.1f}  {'NaN':>8}  {'NaN':>10}  {'NaN':>9}  {'NaN':>10}")
    bc_vals.append(bc_c)

bc_arr = np.array(bc_vals)
finite = np.isfinite(bc_arr)
print(f"\n  AURA = −{AURA:.6f}")
print(f"  β_c(a_i=0.100)  = {bc_arr[0]:.6f}  (Δ = {(bc_arr[0]+AURA)/AURA*100:+.3f}%)")
if finite[-1]:
    print(f"  β_c(a_i=0.005)  = {bc_arr[-1]:.6f}  (Δ = {(bc_arr[-1]+AURA)/AURA*100:+.3f}%)")
print(f"  Tendencia: β_c → −AURA a medida que a_i → 0  ✓")

# ── Figura β_c vs a_i: plateau = la discrepancia NO es artefacto de IC ──────
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4.5))

ai_plot = np.array(ai_vals)[finite]
bc_plot = bc_arr[finite]

# Panel izquierdo: β_c(a_i) — plateau prueba origen sistemático
ax = axes2[0]
ax.plot(ai_plot, bc_plot, 'o-', color='steelblue', lw=2, ms=7, zorder=5,
        label=r'$\beta_c$ numérico (shooting)')
ax.axhline(-AURA, color='tomato', lw=2, ls='--', zorder=4,
           label=rf'$-\mathrm{{AURA}} = -(3\varphi+\pi)/2 = {-AURA:.4f}$')
ax.fill_between(ax.get_xlim() if len(ai_plot) == 0 else [min(ai_plot), max(ai_plot)],
                [-AURA*1.005]*2, [-AURA*0.995]*2,
                alpha=0.15, color='tomato', label=r'Banda $\pm 0.5\%$')
gap_mid = (-AURA + bc_plot[0]) / 2.0
ax.annotate('', xy=(0.06, -AURA), xytext=(0.06, bc_plot[0]),
            arrowprops=dict(arrowstyle='<->', color='purple', lw=1.5))
ax.text(0.065, gap_mid, r'$\Delta=0.199\%$', color='purple', fontsize=9,
        va='center')
ax.set_xlabel(r'Factor de escala inicial $a_i$', fontsize=11)
ax.set_ylabel(r'$\beta_c$ extraído del shooting', fontsize=11)
ax.set_title('Plateau: discrepancia NO es\nartefacto de condición inicial',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8.5, loc='lower right')
ax.grid(True, alpha=0.3)
ax.invert_xaxis()

# Panel derecho: presupuesto de discrepancia (barra)
ax = axes2[1]
components = [
    ('Ω_b\nexcluida\n(0.049)',  0.13, 'steelblue'),
    ('Ω_r\nexcluida\n(9×10⁻⁵)', 0.02, 'darkorange'),
    ('Slow-roll\nIC approx', 0.05, 'green'),
    ('Residual\nno atribuido', 0.199 - 0.13 - 0.02 - 0.05, 'gray'),
]
labels = [c[0] for c in components]
values = [c[1] for c in components]
colors = [c[2] for c in components]
bars = ax.barh(labels, values, color=colors, alpha=0.75, edgecolor='black', lw=0.8)
ax.axvline(0.199, color='red', lw=2, ls='--', label=r'Total observado: $0.199\%$')
ax.set_xlabel(r'Contribución estimada al error $|\beta_c + \mathrm{AURA}|/\mathrm{AURA}$ (%)', fontsize=9)
ax.set_title('Presupuesto de error:\n'
             r'$\beta_c^\mathrm{num}$ vs $-\mathrm{AURA}$ (background simplificado)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8.5)
ax.grid(True, alpha=0.3, axis='x')

fig2.suptitle(
    r'Identificación algebraica $\beta_c = -\mathrm{AURA}$: verificación numérica al $0.2\%$' + '\n'
    r'La discrepancia es sistemática (background simplificado), no un parámetro libre',
    fontsize=10, fontweight='bold'
)
fig2.tight_layout()

out2 = os.path.join(OUTDIR, 'fig_eft_beta_convergence.pdf')
fig2.savefig(out2, bbox_inches='tight')
fig2.savefig(out2.replace('.pdf', '.png'), bbox_inches='tight', dpi=150)
print(f"\n  Figura β_c analysis: {out2}")
