"""
ssee_alpha_saturation_stepB.py — Pasos B+C numéricos.

Construido sobre las ecuaciones derivadas en el Paso A.
Estrategia: numérico, no simbólico (más robusto a errores de signo).

Sistema autónomo canónico (K=X, potencial exponencial V = V_0 exp(-λφ/M_pl),
acoplamiento conformal Q = (α/M_pl) ρ_m φ̇, sin radiación — apropiado para
el atractor tardío):

  dx/dN = -3x + (√6/2)λ y² - (√6/2)α(1-x²-y²) + (3/2)x(1+x²-y²)
  dy/dN = -(√6/2)λ x y + (3/2)y(1+x²-y²)

Variables:
  x = φ̇/(√6 H M_pl)       (cinética)
  y = √V/(√3 H M_pl)       (potencial)
  Ω_φ = x² + y²
  Ω_m = 1 - x² - y²

Los puntos fijos del atractor de DE de tiempo tardío (FP1) son:
  x* = √6/(2(λ+α))
  y*² = 1 + x*² - (√6/3)λ x*

(deducido analíticamente en el cuaderno de Paso B, anexado abajo)

Paso B: encontrar FP1 numéricamente, calcular jacobiano, autovalores.
Paso C: escanear α y encontrar la frontera donde max(Re(autovalor)) cruza 0.
"""

import numpy as np
from scipy.optimize import fsolve

SQ6 = np.sqrt(6)

# ────────────────────────────────────────────────────────────────
# Sistema autónomo
# ────────────────────────────────────────────────────────────────
def system(state, alpha, lam):
    """Lado derecho del sistema autónomo (dx/dN, dy/dN)."""
    x, y = state
    Om = 1.0 - x**2 - y**2
    dxdN = (-3*x
            + (SQ6/2) * lam * y**2
            - (SQ6/2) * alpha * Om
            + (3/2) * x * (1 + x**2 - y**2))
    dydN = (-(SQ6/2) * lam * x * y
            + (3/2) * y * (1 + x**2 - y**2))
    return np.array([dxdN, dydN])


def find_FP1(alpha, lam):
    """Atractor tardío DE: x* = √6/(2(λ+α)) según derivación analítica."""
    denom = lam + alpha
    if abs(denom) < 1e-10:
        return None
    x0 = SQ6 / (2*denom)
    y0sq = 1 + x0**2 - (SQ6/3)*lam*x0
    if y0sq <= 0:
        return None
    y0 = np.sqrt(y0sq)
    try:
        sol, info, ier, msg = fsolve(system, [x0, y0], args=(alpha, lam),
                                      full_output=True)
        if ier != 1:
            return None
        # Verificar que es fixed point (residuo pequeño)
        res = system(sol, alpha, lam)
        if np.max(np.abs(res)) > 1e-8:
            return None
        return sol
    except Exception:
        return None


def jacobian(state, alpha, lam, eps=1e-7):
    """Jacobiano numérico por diferencias finitas centradas."""
    J = np.zeros((2,2))
    for i in range(2):
        sp_ = state.copy(); sp_[i] += eps
        sm_ = state.copy(); sm_[i] -= eps
        J[:,i] = (system(sp_, alpha, lam) - system(sm_, alpha, lam)) / (2*eps)
    return J


def analyze(alpha, lam):
    FP = find_FP1(alpha, lam)
    if FP is None:
        return None
    x, y = FP
    Om_phi = x**2 + y**2
    Om_m = 1.0 - Om_phi
    w_phi = (x**2 - y**2) / Om_phi if Om_phi > 0 else None
    J = jacobian(np.array(FP), alpha, lam)
    eig = np.linalg.eigvals(J)
    max_re = max(e.real for e in eig)
    return dict(x=x, y=y, Om_phi=Om_phi, Om_m=Om_m, w_phi=w_phi,
                eig=eig, max_re=max_re, stable=(max_re < 0))


# ────────────────────────────────────────────────────────────────
# Paso B — escaneo
# ────────────────────────────────────────────────────────────────
print("="*78)
print("PASO B — Atractor DE de tiempo tardío + estabilidad lineal")
print("="*78)
print()
print(f"{'λ':>5} {'α':>8} {'x*':>8} {'y*':>8} {'Ω_φ*':>7} {'Ω_m*':>7} "
      f"{'w_φ*':>7} {'maxRe':>9} {'stable':>6}")
print("-"*78)

for lam in [1.0, 1.5, 2.0, 2.5]:
    print(f"--- λ = {lam} ---")
    for alpha in np.linspace(-0.5, 3.0, 36):
        r = analyze(alpha, lam)
        if r is None:
            continue
        if not (0 < r['Om_phi'] < 1):
            continue
        st = "✓" if r['stable'] else "✗"
        print(f"{lam:5.2f} {alpha:+8.3f} {r['x']:8.4f} {r['y']:8.4f} "
              f"{r['Om_phi']:7.4f} {r['Om_m']:7.4f} "
              f"{r['w_phi']:+7.4f} {r['max_re']:+9.4f} {st:>6}")

# ────────────────────────────────────────────────────────────────
# Paso C — encontrar α de saturación (frontera de estabilidad)
# ────────────────────────────────────────────────────────────────
print()
print("="*78)
print("PASO C — Frontera de estabilidad α_sat(λ)")
print("="*78)
print()
print("Buscando α donde max Re(autovalor) cruza 0 (transición estable→inestable)")
print()
print(f"{'λ':>6} {'α_sat':>10} {'Ω_φ*':>9} {'w_φ*':>9}")
print("-"*40)

def find_alpha_sat(lam):
    """Bisección: encontrar α donde max_re cruza 0."""
    # Buscar rangos donde cambia de signo
    alphas = np.linspace(-2, 5, 200)
    prev = None
    crossings = []
    for a in alphas:
        r = analyze(a, lam)
        if r is None or not (0 < r['Om_phi'] < 1):
            prev = None
            continue
        if prev is not None and np.sign(r['max_re']) != np.sign(prev['max_re']):
            # Bisección refinada
            a_lo, a_hi = prev['a'], a
            for _ in range(50):
                a_mid = (a_lo + a_hi) / 2
                r_mid = analyze(a_mid, lam)
                if r_mid is None:
                    break
                if np.sign(r_mid['max_re']) == np.sign(prev['max_re']):
                    a_lo = a_mid; prev = {'a': a_mid, 'max_re': r_mid['max_re']}
                else:
                    a_hi = a_mid
            crossings.append((a_mid, r_mid))
        prev = {'a': a, 'max_re': r['max_re']}
    return crossings

for lam in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    crossings = find_alpha_sat(lam)
    if not crossings:
        print(f"{lam:6.2f}  {'(sin frontera dentro del rango escaneado)':>30}")
        continue
    for a_sat, r in crossings:
        print(f"{lam:6.2f} {a_sat:+10.5f} {r['Om_phi']:9.4f} {r['w_phi']:+9.4f}")

# ────────────────────────────────────────────────────────────────
# Sanity check vs Amendola 2000
# ────────────────────────────────────────────────────────────────
print()
print("="*78)
print("Sanity check vs literatura (Amendola 2000)")
print("="*78)
print()
print("""
Amendola 2000 (arXiv:astro-ph/9908023): para coupled quintessence con V
exponencial y β=α (Q = (α/M_pl) ρ_m φ̇), la cota de estabilidad del atractor
de scaling es:

      α² < (3 - λ²) · (algo)   ó    relaciones tipo α² + λα ≷ 3

dependiendo del régimen. El cruce concreto α_sat(λ) calculado arriba debe
coincidir cualitativamente con la literatura.

Si nuestro α_sat sale ≈ ±1.2247 = √(3/2) para λ→0, recuperamos el límite
canónico clásico.
""")
