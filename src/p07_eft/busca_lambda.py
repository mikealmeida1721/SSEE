#!/usr/bin/env python3
"""
PREGUNTA: existe una pendiente del potencial
que, INTEGRADA de verdad (no con la formula
asintotica del atractor), aterrice en
w(a=1) = w0 = -0.839950 SIN acoplamiento?

Si existe, beta_c sobra: la brecha Delta_w
no era fisica, era la distancia entre una
formula asintotica y una corrida.

INCOGNITAS (2):  phi(a_i),  ALPHA
CONDICIONES (2) en a=1:
   Omega_DE = 0.691119
   w_phi    = -0.839950

CONTROL: con ALPHA fijado al valor del
atractor debe reproducir w = -0.927316,
que es lo que ya da fondo_disparo.py.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S

KAL = S.KAL0
W0 = S.W0
OM_M = S.OMEGA_M_TOTAL
OM_PHI = 1.0 - OM_M
LAM_ATR = np.sqrt(3.0 * (1.0 + W0))
ALPHA_ATR = np.sqrt(3.0) * LAM_ATR / np.sqrt(KAL)
A_I = 1.0e-3
V0 = 1.0
PHIP_I = 1e-8


def _campo(phi, phip, a, alpha):
    rho_m = OM_M * a ** -3
    V = V0 * np.exp(-alpha * phi)
    den = 1.0 - phip ** 2 / (2.0 * KAL)
    H2 = (V + rho_m) / den
    X = H2 * phip ** 2 / 2.0
    return rho_m, V, H2, X


def rhs(N, Y, alpha):
    phi, phip = Y
    a = np.exp(N)
    rho_m, V, H2, X = _campo(phi, phip, a, alpha)
    if H2 <= 0:
        return [0.0, 0.0]
    p_phi = X / KAL - V
    rho_tot = X / KAL + V + rho_m
    HpH = -1.5 * (1.0 + p_phi / rho_tot)
    phipp = (-(3.0 + HpH) * phip
             + alpha * V * KAL / H2)
    return [phip, phipp]


def corre(phi_i, alpha, a_fin=1.05, n=1200):
    sol = solve_ivp(rhs, [np.log(A_I), np.log(a_fin)],
                    [phi_i, PHIP_I], args=(alpha,),
                    t_eval=np.linspace(np.log(A_I),
                                       np.log(a_fin), n),
                    rtol=1e-11, atol=1e-14,
                    method='Radau')
    if not sol.success:
        return None
    a = np.exp(sol.t)
    phi, phip = sol.y
    rho_m, V, H2, X = _campo(phi, phip, a, alpha)
    rho_phi = X / KAL + V
    p_phi = X / KAL - V
    tot = rho_phi + rho_m
    return dict(a=a, w=p_phi / rho_phi,
                Om_DE=rho_phi / tot)


def _en_a1(r, k):
    return float(np.interp(1.0, r['a'], r[k]))


def objetivo(u):
    phi_i, alpha = u
    r = corre(phi_i, alpha)
    if r is None:
        return [1e3, 1e3]
    return [_en_a1(r, 'Om_DE') - OM_PHI,
            _en_a1(r, 'w') - W0]


def solo_phi(phi_i, alpha):
    r = corre(float(np.ravel(phi_i)[0]), alpha)
    if r is None:
        return 1e3
    return _en_a1(r, 'Om_DE') - OM_PHI


if __name__ == '__main__':
    print('=' * 62)
    print('  BUSCA LAMBDA — el potencial que llega')
    print('  a w0 SIN acoplamiento')
    print('=' * 62)
    print('  objetivo w0     = %.6f' % W0)
    print('  alpha atractor  = %.6f' % ALPHA_ATR)
    print('  lambda atractor = %.6f' % LAM_ATR)
    print()

    # CONTROL: alpha del atractor
    p = fsolve(solo_phi, 0.5, args=(ALPHA_ATR,),
               xtol=1e-12)[0]
    rc = corre(p, ALPHA_ATR)
    wc = _en_a1(rc, 'w')
    print('  CONTROL (alpha del atractor):')
    print('    phi_i    = %.8f' % p)
    print('    Om_DE    = %.9f' % _en_a1(rc, 'Om_DE'))
    print('    w(a=1)   = %.6f' % wc)
    print('    esperado = -0.927316   dif %.2e'
          % abs(wc + 0.927316))
    print()

    sol = None
    for a0 in (0.51, 0.8, 1.2, 2.0, 3.0, 5.0, 0.3):
        for p0 in (0.5, 0.0, 1.0, 2.0):
            u, info, ier, m = fsolve(
                objetivo, [p0, a0],
                full_output=True, xtol=1e-13)
            if ier == 1 and max(abs(
                    np.array(objetivo(u)))) < 1e-9:
                sol = u
                break
        if sol is not None:
            break

    if sol is None:
        print('  NO HAY SOLUCION en el rango')
        sys.exit(0)

    phi_i, alpha = sol
    lam = alpha * np.sqrt(KAL) / np.sqrt(3.0)
    r = corre(phi_i, alpha)
    print('  SOLUCION SIN ACOPLAMIENTO:')
    print('    phi_i   = %.8f' % phi_i)
    print('    alpha   = %.8f' % alpha)
    print('    lambda  = %.8f' % lam)
    print('    Om_DE   = %.9f  (obj %.9f)'
          % (_en_a1(r, 'Om_DE'), OM_PHI))
    print('    w(a=1)  = %.9f  (obj %.9f)'
          % (_en_a1(r, 'w'), W0))
    print()
    print('  COMPARACION de lambda:')
    print('    del atractor  = %.8f' % LAM_ATR)
    print('    de la corrida = %.8f' % lam)
    print('    razon         = %.6f'
          % (lam / LAM_ATR))
