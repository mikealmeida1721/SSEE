"""
ssee_alpha_saturation_stepD2.py — versión limpia.

El φ-MDE tiene y=0 (purely kinetic during matter era). Lo busco directamente:

Para K(X) = X/KAL + X²/M⁴ con conformal coupling α y σ=H²M_pl²/M⁴:

  dx/dN |_{y=0} = (√6/2)·[−α·Ω_m − √6·x·K_X]/(1/KAL + 18σx²) − x·Ḣ/H²
  con:
    Ω_φ = x²/KAL + 9σx⁴
    Ω_m = 1 − Ω_φ
    K_X = 1/KAL + 6σx²
    Ḣ/H² = −(3/2)(1 + x²/KAL + 3σx⁴)

Saturación α_sat(σ, KAL): el α al cual el φ-MDE alcanza Ω_φ = 1 (deja de ser
era de materia y se vuelve scalar-dominated).

Para σ=0, KAL=1: x_φMDE = −α√6/3, Ω_φ = 2α²/3 = 1 → α_sat = √(3/2) ≈ 1.2247.
"""

import numpy as np
from scipy.optimize import fsolve, brentq

SQ6 = np.sqrt(6)

# SSEE constants
PHI = (1 + np.sqrt(5))/2
PI  = np.pi
BETA = (PHI + PI)/2
KAL0 = BETA + PI                  # ≈ 5.5214
AURA = BETA + PHI                 # ≈ 3.998
MIRA = AURA/2


def dxdN_y0(x, alpha, KAL, sig):
    """dx/dN evaluado en y=0 (φ-MDE)."""
    Om_phi = x**2/KAL + 9*sig*x**4
    Om_m   = 1.0 - Om_phi
    w_eff  = x**2/KAL + 3*sig*x**4
    Hdot   = -1.5*(1 + w_eff)
    K_X    = 1.0/KAL + 6*sig*x**2
    denom  = 1.0/KAL + 18*sig*x**2
    num    = -alpha*Om_m - SQ6*x*K_X
    return (SQ6/2)*num/denom - x*Hdot


def x_phimde(alpha, KAL, sig):
    """Solve for x_φMDE."""
    # Initial guess from canonical: x = -α√6/3
    x_seed = -alpha*SQ6/3
    try:
        sol = fsolve(lambda x: dxdN_y0(x[0], alpha, KAL, sig),
                     [x_seed], xtol=1e-14, full_output=False)
        x = sol[0]
        # Verify residual
        if abs(dxdN_y0(x, alpha, KAL, sig)) > 1e-9:
            return None
        return x
    except Exception:
        return None


def Omega_phi_phimde(alpha, KAL, sig):
    x = x_phimde(alpha, KAL, sig)
    if x is None:
        return None
    return x**2/KAL + 9*sig*x**4


def alpha_sat(KAL, sig, search_lo=0.001, search_hi=10.0):
    """Find α where Ω_φ at φ-MDE = 1."""
    def gap(alpha):
        Om = Omega_phi_phimde(alpha, KAL, sig)
        if Om is None:
            return 100.0
        return Om - 1.0
    # Find sign change first
    try:
        return brentq(gap, search_lo, search_hi, xtol=1e-10)
    except ValueError:
        # No sign change; bracketing failed
        return None


# ────────────────────────────────────────────────────────────────
# Sanity check 1: canónico σ=0, KAL=1 → α_sat = √(3/2)
# ────────────────────────────────────────────────────────────────
print("="*78)
print("Sanity check 1: σ=0, KAL=1 → α_sat = √(3/2)?")
print("="*78)
a_sat_can = alpha_sat(1.0, 0.0)
print(f"\n  α_sat(canónico) = {a_sat_can:.8f}")
print(f"  √(3/2)          = {np.sqrt(1.5):.8f}")
print(f"  diferencia     = {(a_sat_can-np.sqrt(1.5))/np.sqrt(1.5)*100:+.4e}%")
print()

# ────────────────────────────────────────────────────────────────
# Sanity check 2: detalle del φ-MDE para varios α en canónico
# ────────────────────────────────────────────────────────────────
print("="*78)
print("Sanity check 2: φ-MDE en canónico (verifica fórmula x=-α√6/3)")
print("="*78)
print(f"\n{'α':>8} {'x_φMDE':>10} {'x_teo(-α√6/3)':>15} {'Ω_φ':>9} {'2α²/3':>10}")
for a in [0.3, 0.5, 0.8, 1.0, 1.2]:
    x_num = x_phimde(a, 1.0, 0.0)
    x_teo = -a*SQ6/3
    Om = a**2 * 2/3
    Om_num = Omega_phi_phimde(a, 1.0, 0.0)
    print(f"{a:8.3f} {x_num:10.5f} {x_teo:15.5f} {Om_num:9.5f} {Om:10.5f}")
print()

# ────────────────────────────────────────────────────────────────
# SSEE: escanear σ para KAL = KAL₀
# ────────────────────────────────────────────────────────────────
print("="*78)
print(f"SSEE: α_sat(σ) con KAL = KAL₀ = {KAL0:.6f}")
print("="*78)
print()
print(f"{'σ':>12} {'régimen':>22} {'α_sat':>12} {'x_φMDE':>10} {'Ω_φ_check':>10}")
print("-"*68)

sig_values = [0.0, 1e-9, 1e-6, 1e-3, 1.0, 1e3, 1e6, 1e9]
results = {}
for sig in sig_values:
    a_sat = alpha_sat(KAL0, sig)
    if a_sat is None:
        print(f"{sig:12.2e}  {'(no encontrado)':>22}")
        continue
    x_at_sat = x_phimde(a_sat, KAL0, sig)
    Om_check = x_at_sat**2/KAL0 + 9*sig*x_at_sat**4
    # Régimen
    # 1/KAL vs 18σx² at x ~ -0.5:
    lin = 1/KAL0; cuad = 18*sig*0.25
    if cuad < 0.1*lin: regime = "X/KAL dominante"
    elif cuad > 10*lin: regime = "X²/M⁴ dominante"
    else: regime = "transición"
    results[sig] = a_sat
    print(f"{sig:12.2e}  {regime:>22} {a_sat:12.6f} {x_at_sat:10.5f} {Om_check:10.6f}")

# ────────────────────────────────────────────────────────────────
# Comparación con expresiones algebraicas SSEE
# ────────────────────────────────────────────────────────────────
print()
print("="*78)
print("Comparación con expresiones algebraicas φ,π")
print("="*78)

# Candidatos
candidates = {
    "φ/π":               PHI/PI,
    "π/φ²":              PI/PHI**2,
    "1/φ":               1/PHI,
    "2/φ":               2/PHI,
    "√(3/2)":            np.sqrt(1.5),
    "√(3/(2·KAL₀))":     np.sqrt(3/(2*KAL0)),
    "√(KAL₀)":           np.sqrt(KAL0),
    "1/√KAL₀":           1/np.sqrt(KAL0),
    "MIRA−1":            MIRA-1,
    "√(MIRA)":           np.sqrt(MIRA),
    "1/(2β)":            1/(2*BETA),
    "1/β":               1/BETA,
    "β":                 BETA,
    "AURA/π²":           AURA/PI**2,
    "π−φ−1":             PI-PHI-1,
    "(π−φ)/2":           (PI-PHI)/2,
    "1/φ³":              1/PHI**3,
    "π/(φ+π)":           PI/(PHI+PI),
    "φ/(φ+π)":           PHI/(PHI+PI),
}

print()
for sig, a_sat in results.items():
    print(f"\nσ = {sig:.2e}, α_sat = {a_sat:.6f}")
    print(f"{'Expresión':>20} {'Valor':>10} {'Desvío%':>10}")
    print("-"*44)
    sorted_c = sorted(candidates.items(), key=lambda kv: abs(kv[1]-a_sat))
    for name, val in sorted_c[:5]:
        dev = (val - a_sat)/a_sat * 100
        marker = ""
        if abs(dev) < 0.01: marker = "  ← EXACTO"
        elif abs(dev) < 0.1: marker = "  ← <0.1%"
        elif abs(dev) < 1: marker = "  ← <1%"
        print(f"{name:>20} {val:>10.6f} {dev:>+9.4f}%{marker}")
