"""
ssee_alpha_saturation_stepC.py — busca TODOS los puntos fijos, no solo FP1.

Corrección del intento anterior: FP1 (atractor tardío DE) es siempre estable,
así que no da frontera. La frontera de Amendola α²<3/2 vive en otro punto
fijo (φ-MDE = era de materia con escalar activo). Lo busco aquí.

Estrategia:
  1. Gridear (x,y) inicial → fsolve → coleccionar puntos fijos únicos.
  2. Para cada FP: calcular Ω_φ, Ω_m, w_φ, jacobiano, autovalores.
  3. Clasificar: late-time DE (Ω_φ~0.7), φ-MDE (Ω_φ pequeño), scalar-dom (Ω_φ=1).
  4. Para cada tipo: escanear α y encontrar cruce de estabilidad.
"""

import numpy as np
from scipy.optimize import fsolve
from itertools import product

SQ6 = np.sqrt(6)


def system(state, alpha, lam):
    x, y = state
    Om = 1.0 - x**2 - y**2
    dxdN = (-3*x + (SQ6/2)*lam*y**2 - (SQ6/2)*alpha*Om
            + (3/2)*x*(1 + x**2 - y**2))
    dydN = (-(SQ6/2)*lam*x*y + (3/2)*y*(1 + x**2 - y**2))
    return np.array([dxdN, dydN])


def jacobian(state, alpha, lam, eps=1e-7):
    J = np.zeros((2,2))
    for i in range(2):
        sp_ = state.copy(); sp_[i] += eps
        sm_ = state.copy(); sm_[i] -= eps
        J[:,i] = (system(sp_, alpha, lam) - system(sm_, alpha, lam))/(2*eps)
    return J


def find_all_FPs(alpha, lam, tol=1e-9):
    """Buscar todos los puntos fijos por grid de condiciones iniciales."""
    fps = []
    grid_x = np.linspace(-1.0, 1.0, 11)
    grid_y = np.linspace(0.0, 1.0, 11)
    for x0, y0 in product(grid_x, grid_y):
        try:
            sol = fsolve(system, [x0, y0], args=(alpha, lam),
                          full_output=False, xtol=1e-12)
            res = system(sol, alpha, lam)
            if np.max(np.abs(res)) > tol:
                continue
            # Filtrar y reales (y² puede ser <0 en raíz espuria — rechazar y<0)
            if sol[1] < -1e-6:
                continue
            sol[1] = abs(sol[1])  # y = √V definido positivo
            # Filtrar duplicados
            is_new = True
            for fp in fps:
                if np.linalg.norm(sol - fp) < 1e-4:
                    is_new = False; break
            if is_new:
                fps.append(sol.copy())
        except Exception:
            pass
    return fps


def classify(fp, alpha, lam):
    x, y = fp
    Om_phi = x**2 + y**2
    Om_m = 1.0 - Om_phi
    if Om_phi > 0:
        w_phi = (x**2 - y**2) / Om_phi
    else:
        w_phi = None
    J = jacobian(fp, alpha, lam)
    eigs = np.linalg.eigvals(J)
    max_re = max(e.real for e in eigs)
    stable = max_re < -1e-8
    # Etiqueta de tipo
    if Om_phi > 0.95:
        kind = "scalar-dom"
    elif Om_phi > 0.5:
        kind = "DE-attractor"
    elif Om_phi > 0.01:
        kind = "φ-MDE/scaling"
    else:
        kind = "matter-dom"
    return dict(x=x, y=y, Om_phi=Om_phi, Om_m=Om_m, w_phi=w_phi,
                eigs=eigs, max_re=max_re, stable=stable, kind=kind)


# ────────────────────────────────────────────────────────────────
# Escaneo: para cada (α,λ) catalogar TODOS los FPs
# ────────────────────────────────────────────────────────────────
print("="*84)
print("PASO C corregido — buscar TODOS los puntos fijos por (α, λ)")
print("="*84)

# Para entender el bestiario, fijar λ=1 y escanear α
print("\n--- λ = 1.0, escaneando α ---\n")
print(f"{'α':>7} {'kind':>15} {'x':>8} {'y':>8} {'Ω_φ':>7} {'w_φ':>9} "
      f"{'maxRe':>9} {'stable':>6}")
print("-"*84)

for alpha in [0.0, 0.3, 0.5, 0.7, 1.0, 1.225, 1.5, 2.0]:
    fps = find_all_FPs(alpha, 1.0)
    for fp in fps:
        c = classify(fp, alpha, 1.0)
        w_str = f"{c['w_phi']:+9.4f}" if c['w_phi'] is not None else "    —    "
        st = "✓" if c['stable'] else "✗"
        print(f"{alpha:+7.3f} {c['kind']:>15} {c['x']:+8.4f} {c['y']:+8.4f} "
              f"{c['Om_phi']:7.4f} {w_str} {c['max_re']:+9.4f} {st:>6}")
    print()

# ────────────────────────────────────────────────────────────────
# Ahora: para el φ-MDE específicamente, escanear α y buscar transición
# ────────────────────────────────────────────────────────────────
print("="*84)
print("φ-MDE: buscar cruce de estabilidad en función de α (λ fijo)")
print("="*84)

for lam in [0.5, 1.0, 1.5, 2.0]:
    print(f"\n--- λ = {lam} ---")
    print(f"{'α':>8} {'x':>9} {'y':>9} {'Ω_φ':>8} {'w_eff':>9} {'maxRe':>10} stable")
    prev_stable = None
    for alpha in np.linspace(0.01, 2.5, 50):
        fps = find_all_FPs(alpha, lam)
        for fp in fps:
            c = classify(fp, alpha, lam)
            if c['kind'] != 'φ-MDE/scaling':
                continue
            # w_eff total (incluye materia con peso 0)
            w_eff = c['w_phi']*c['Om_phi'] if c['w_phi'] is not None else 0.0
            st = "✓" if c['stable'] else "✗"
            sign_change = ""
            if prev_stable is not None and prev_stable != c['stable']:
                sign_change = "  ← TRANSICIÓN"
            prev_stable = c['stable']
            print(f"{alpha:+8.4f} {c['x']:+9.4f} {c['y']:+9.4f} "
                  f"{c['Om_phi']:8.4f} {w_eff:+9.4f} {c['max_re']:+10.4f} "
                  f"{st:>6}{sign_change}")

# ────────────────────────────────────────────────────────────────
# Existencia del φ-MDE: ¿cuándo deja de existir como FP físico?
# ────────────────────────────────────────────────────────────────
print()
print("="*84)
print("Existencia del φ-MDE: α al cual el FP φ-MDE deja de existir físicamente")
print("="*84)

for lam in [0.5, 1.0, 1.5, 2.0, 2.5]:
    last_alpha_phimde = None
    for alpha in np.linspace(0.01, 3.0, 300):
        fps = find_all_FPs(alpha, lam)
        has_phimde = any(classify(fp, alpha, lam)['kind'] == 'φ-MDE/scaling'
                          for fp in fps)
        if has_phimde:
            last_alpha_phimde = alpha
    print(f"λ={lam:.2f} → último α con φ-MDE físico ≈ {last_alpha_phimde:+.4f}"
          if last_alpha_phimde is not None
          else f"λ={lam:.2f} → φ-MDE no encontrado")
