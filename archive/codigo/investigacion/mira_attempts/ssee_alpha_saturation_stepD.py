"""
ssee_alpha_saturation_stepD.py — extender el sistema a K(X) = X/KAL + X²/M⁴.

Variables autónomas (Tsujikawa-Sami style para k-essence):
  x = φ̇/(√6 H M_pl)        (cinética adimensional)
  y = √V/(√3 H M_pl)        (potencial adimensional)
  σ = H² M_pl² / M⁴          (ratio Hubble/k-essence — DINÁMICA)

Derivadas (matter + scalar, sin radiación, conformal coupling):

  Ω_φ = x²/KAL + 9σx⁴ + y²
  Ω_m = 1 − Ω_φ
  w_eff = x²/KAL + 3σx⁴ − y²     (total, materia tiene w=0)
  Ḣ/H² = −(3/2)(1 + w_eff)

  K_X = 1/KAL + 6σx²
  K_X + 2X·K_XX = 1/KAL + 18σx²

  Ecuación KG (acoplada conformal):
    φ̈ = [−(α/M_pl)ρ_m − 3H·K_X·φ̇ − V'] / (K_X + 2X·K_XX)

  En variables adimensionales:
    x' = (√6/2)·[−α·Ω_m − √6·x·K_X + λy²] / (K_X + 18σx²·(2/(2KAL_NA)))
         −x·(Ḣ/H²)

  Específicamente:
    x' = (√6/2)·[−α·Ω_m − √6·x·(1/KAL+6σx²) + λy²]/(1/KAL+18σx²)
         − x·Ḣ/H²

    y' = −(√6/2)λxy + (3/2)y·(1 + x²/KAL + 3σx⁴ − y²)

    σ' = 2σ·(Ḣ/H²) = −3σ(1 + x²/KAL + 3σx⁴ − y²)

Sanity check: con KAL=1, σ=0 → recupera canónico (verificado al final).

Nota: σ EVOLUCIONA. No es parámetro fijo. Por eso buscamos no "puntos fijos
del 3D" (que serían trivialmente σ=0 o de Sitter), sino la frontera de
existencia/estabilidad del φ-MDE tomando σ como condición instantánea.

Esto se hace en dos pasos:
  D1. Para σ fijado a un valor (matter era ≈ 10⁷, hoy ≈ 10⁻³), buscar el
      "instantáneo" FP del subsystema (x,y) con σ congelado.
  D2. Variar σ y α, mapear la frontera de existencia del φ-MDE → α_sat(σ).
"""

import numpy as np
from scipy.optimize import fsolve
from itertools import product

SQ6 = np.sqrt(6)

# SSEE constants
PHI = (1 + np.sqrt(5))/2
PI  = np.pi
BETA = (PHI + PI)/2
KAL0 = BETA + PI                  # 5.5214
AURA = BETA + PHI                 # 3.998
MIRA = AURA/2                     # 1.999


def system_3D(state, alpha, lam, KAL):
    """Lado derecho del sistema 3D (x,y,σ). σ congelado: σ' devuelto pero no usado."""
    x, y, sig = state
    Om_phi = x**2/KAL + 9*sig*x**4 + y**2
    Om_m = 1.0 - Om_phi
    w_eff = x**2/KAL + 3*sig*x**4 - y**2
    Hdot_H2 = -1.5*(1 + w_eff)
    K_X = 1.0/KAL + 6*sig*x**2
    denom = 1.0/KAL + 18*sig*x**2
    # dx/dN
    num = -alpha*Om_m - SQ6*x*K_X + lam*y**2
    dxdN = (SQ6/2)*num/denom - x*Hdot_H2
    # dy/dN
    dydN = -(SQ6/2)*lam*x*y + 1.5*y*(1 + x**2/KAL + 3*sig*x**4 - y**2)
    # dσ/dN
    dsigdN = 2*sig*Hdot_H2
    return np.array([dxdN, dydN, dsigdN])


def subsystem_xy(state_xy, alpha, lam, KAL, sig_fixed):
    """Para encontrar FP de (x,y) con σ congelado."""
    full = np.array([state_xy[0], state_xy[1], sig_fixed])
    rhs = system_3D(full, alpha, lam, KAL)
    return rhs[:2]   # solo dx, dy


def find_FPs_xy(alpha, lam, KAL, sig_fixed, tol=1e-9):
    fps = []
    grid_x = np.linspace(-1.0, 1.0, 11)
    grid_y = np.linspace(0.0, 1.0, 11)
    for x0, y0 in product(grid_x, grid_y):
        try:
            sol = fsolve(subsystem_xy, [x0, y0],
                          args=(alpha, lam, KAL, sig_fixed),
                          full_output=False, xtol=1e-12)
            res = subsystem_xy(sol, alpha, lam, KAL, sig_fixed)
            if np.max(np.abs(res)) > tol:
                continue
            if sol[1] < -1e-6:
                continue
            sol[1] = abs(sol[1])
            is_new = True
            for fp in fps:
                if np.linalg.norm(sol - fp) < 1e-4:
                    is_new = False; break
            if is_new:
                fps.append(sol.copy())
        except Exception:
            pass
    return fps


def classify_xy(fp, alpha, lam, KAL, sig):
    x, y = fp
    Om_phi = x**2/KAL + 9*sig*x**4 + y**2
    Om_m = 1.0 - Om_phi
    if Om_phi > 1e-3:
        w_phi = (x**2/KAL + 3*sig*x**4 - y**2)/Om_phi
    else:
        w_phi = None
    if Om_phi > 0.95: kind = "scalar-dom"
    elif Om_phi > 0.5: kind = "DE-attractor"
    elif Om_phi > 0.01: kind = "φ-MDE/scaling"
    else: kind = "matter-dom"
    return dict(x=x, y=y, Om_phi=Om_phi, Om_m=Om_m, w_phi=w_phi, kind=kind)


# ────────────────────────────────────────────────────────────────
# Sanity check: σ=0, KAL=1 debe reproducir canónico (α_sat = √3/2)
# ────────────────────────────────────────────────────────────────
print("="*78)
print("Sanity check: σ=0, KAL=1 → canónico → α_sat = √(3/2)")
print("="*78)
for alpha in [0.5, 1.0, 1.2, 1.225, 1.25, 1.3, 1.5]:
    fps = find_FPs_xy(alpha, lam=1.0, KAL=1.0, sig_fixed=0.0)
    phi_mde = [classify_xy(fp, alpha, 1.0, 1.0, 0.0) for fp in fps]
    phi_mde = [c for c in phi_mde if c['kind'] == 'φ-MDE/scaling']
    if phi_mde:
        c = phi_mde[0]
        print(f"  α={alpha:+.4f}: φ-MDE x={c['x']:+.4f} y={c['y']:+.4f} "
              f"Ω_φ={c['Om_phi']:.4f}")
    else:
        print(f"  α={alpha:+.4f}: φ-MDE NO existe")

print()
print("Esperado: φ-MDE existe hasta α≈1.225=√(3/2), luego desaparece.")
print()

# ────────────────────────────────────────────────────────────────
# Caso SSEE: KAL=KAL₀, escanear σ y α
# ────────────────────────────────────────────────────────────────
print("="*78)
print(f"SSEE: KAL = KAL₀ = {KAL0:.4f}, escanear σ y α")
print("="*78)

# Para cada σ, encontrar el último α donde existe el φ-MDE
def alpha_max_phimde(sig, KAL, lam=1.0, alpha_grid=None):
    if alpha_grid is None:
        alpha_grid = np.linspace(0.01, 3.0, 200)
    last = None
    for a in alpha_grid:
        fps = find_FPs_xy(a, lam, KAL, sig)
        has = any(classify_xy(fp, a, lam, KAL, sig)['kind']=='φ-MDE/scaling'
                   for fp in fps)
        if has:
            last = a
    return last

print()
print(f"{'σ':>10} {'régimen':>20} {'α_sat(σ)':>12}")
print("-"*50)

# σ values: from σ≈0 (canonical limit) through σ_today (~10⁻³) to σ_eq (~10⁷)
sig_values = [0.0, 1e-6, 1e-3, 1.0, 1e3, 1e6, 1e9]
results = {}
for sig in sig_values:
    a_max = alpha_max_phimde(sig, KAL0)
    # Identificar régimen
    # X/KAL domina cuando 1/KAL ~ 6σx² ; con x~0.5: σ ~ 1/(24·KAL) ~ 0.0075
    if sig < 1e-3:
        regime = "X/KAL dominante"
    elif sig > 100:
        regime = "X²/M⁴ dominante"
    else:
        regime = "transición"
    results[sig] = a_max
    print(f"{sig:10.2e} {regime:>20} {a_max if a_max else 'no encontrado':>12}")

# ────────────────────────────────────────────────────────────────
# Paso E: comparar α_sat con expresiones algebraicas φ,π
# ────────────────────────────────────────────────────────────────
print()
print("="*78)
print("PASO E — comparar α_sat(σ) con expresiones algebraicas SSEE")
print("="*78)

# Expresiones candidatas en φ,π
candidates = {
    "φ/π":               PHI/PI,
    "π/φ²":              PI/PHI**2,
    "1/φ":               1/PHI,
    "2/φ":               2/PHI,
    "√(3/2)":            np.sqrt(1.5),
    "√(3/(φ+3π))":       np.sqrt(3/(PHI+3*PI)),
    "√(3/(2·KAL₀))":     np.sqrt(3/(2*KAL0)),
    "√(MIRA−1)":         np.sqrt(MIRA-1),
    "MIRA−1":            MIRA-1,
    "1/(2β)":            1/(2*BETA),
    "1/β":               1/BETA,
    "π−φ−1":             PI-PHI-1,
    "(π−φ)/2":           (PI-PHI)/2,
    "AURA/π²":           AURA/PI**2,
    "1/φ²·π":            PI/PHI**2,
}

print()
for sig, a_sat in results.items():
    if a_sat is None: continue
    print(f"\nσ = {sig:.2e}, α_sat ≈ {a_sat:.4f}")
    print(f"{'Expresión':>20} {'Valor':>10} {'Desvío':>10}")
    sorted_c = sorted(candidates.items(), key=lambda kv: abs(kv[1]-a_sat))
    for name, val in sorted_c[:5]:
        dev = (val - a_sat)/a_sat * 100
        marker = "  ← <0.1%" if abs(dev) < 0.1 else ""
        print(f"{name:>20} {val:>10.5f} {dev:>+9.3f}%{marker}")
