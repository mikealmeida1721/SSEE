#!/usr/bin/env python3
"""
w(a) y Omega_DE(a) del FONDO REAL de SSEE -- no de la recta CPL.

POR QUE. El numero a* = 1.52103 (la epoca donde Omega_DE = T_r/M_v = 0.839950)
se obtuvo con CPL, w(a) = w0 + wa(1-a). CPL es una RECTA ajustada con datos de
DESI entre a = 0.30 y a = 0.77 -- todo pasado. Extrapolarla a a = 1.52 es
usarla al doble de distancia del dato mas cercano. Si a* depende de la recta,
la "lente temporal" se apoya en un artefacto; si sobrevive al campo real, es
una propiedad del modelo.

COMO. Se integra el campo escalar de SSEE -- K(X) = X/KAL, V = V0 exp(-alpha
phi) -- HACIA ADELANTE desde hoy. Hacia adelante el potencial exponencial es
un ATRACTOR, asi que la integracion es estable; hacia atras es repulsivo y por
eso reventó el intento anterior (E(z=2.5) salia 2x de lo debido).

CONDICIONES DE HOY (a=1), ninguna ajustada:
    Omega_m  = omega_m/h^2 = 0.308881      (ingredientes)
    Omega_DE = 1 - Omega_m = 0.691119      (planitud)
    w_phi    = w0 = -0.839950              (T_r/M_v, lo que ajusta DESI)
de donde salen la energia cinetica K = rho_phi(1+w)/2 y la potencial
V = rho_phi(1-w)/2, y con ellas phi'(1) y V0.

Unidades: rho_crit,0 = 1, H0 = 1  =>  H^2 = rho_tot. N = ln a.
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S

KAL = S.KAL0
W0, WA = S.W0, S.WA
OM_M = S.OMEGA_M_TOTAL
OM_PHI = 1.0 - OM_M
X_KAL_0 = OM_PHI * (1.0 + W0) / 2.0          # densidad cinetica hoy
V0 = OM_PHI * (1.0 - W0) / 2.0               # densidad potencial hoy
PHIP0 = np.sqrt(2.0 * KAL * X_KAL_0)         # H0 = 1
LAM = np.sqrt(3.0 * (1.0 + W0))              # tracker: w = -1 + lam^2/3
ALPHA = LAM / np.sqrt(KAL)


def rhs(N, Y):
    phi, phip = Y
    a = np.exp(N)
    rho_m = OM_M * a ** -3
    V = V0 * np.exp(-ALPHA * phi)
    den = 1.0 - phip ** 2 / (2.0 * KAL)
    H2 = (V + rho_m) / den
    X = H2 * phip ** 2 / 2.0
    p_phi = X / KAL - V
    rho_tot = X / KAL + V + rho_m
    HpH = -1.5 * (1.0 + p_phi / rho_tot)
    phipp = -(3.0 + HpH) * phip + ALPHA * V * KAL / H2
    return [phip, phipp]


def integra(a_max=3.0, n=600):
    sol = solve_ivp(rhs, [0.0, np.log(a_max)], [0.0, PHIP0],
                    t_eval=np.linspace(0.0, np.log(a_max), n),
                    rtol=1e-11, atol=1e-14, method='Radau')
    a = np.exp(sol.t)
    phi, phip = sol.y
    rho_m = OM_M * a ** -3
    V = V0 * np.exp(-ALPHA * phi)
    H2 = (V + rho_m) / (1.0 - phip ** 2 / (2.0 * KAL))
    X = H2 * phip ** 2 / 2.0
    rho_phi = X / KAL + V
    p_phi = X / KAL - V
    return dict(a=a, H=np.sqrt(H2), rho_phi=rho_phi, rho_m=rho_m,
                w=p_phi / rho_phi, Om_DE=rho_phi / (rho_phi + rho_m),
                Om_m=rho_m / (rho_phi + rho_m))


def w_cpl(a):
    return W0 + WA * (1.0 - a)


if __name__ == '__main__':
    r = integra()
    a, w, ODE = r['a'], r['w'], r['Om_DE']
    print('=' * 92)
    print('  FONDO REAL (campo integrado) vs la RECTA CPL')
    print('=' * 92)
    print(f'  alpha = lam/sqrt(KAL) = {ALPHA:.6f}   V0 = {V0:.6f}   '
          f"phi'(1) = {PHIP0:.6f}")
    print(f'  CONTROL en a=1:  w = {w[0]:+.6f}  (debe ser {W0:+.6f}, '
          f'dif {abs(w[0]-W0):.2e})')
    print(f'                   Om_DE = {ODE[0]:.6f}  (debe ser {OM_PHI:.6f})')
    print()
    print(f'  {"a":>7} {"w REAL":>11} {"w CPL":>11} {"dif":>9} '
          f'{"Om_DE REAL":>12} {"Om_m REAL":>11}')
    for at in (1.0, 1.2, 1.5, 1.52103, 2.0, 2.5, 3.0):
        i = int(np.argmin(abs(a - at)))
        print(f'  {a[i]:>7.4f} {w[i]:>+11.6f} {w_cpl(a[i]):>+11.6f} '
              f'{w[i]-w_cpl(a[i]):>+9.4f} {ODE[i]:>12.6f} {r["Om_m"][i]:>11.6f}')
    print()
    if ODE.max() > S.OMEGA_DE > ODE.min():
        f = lambda x: np.interp(x, a, ODE) - S.OMEGA_DE
        a_star = brentq(f, a[0], a[-1])
        print(f'  Omega_DE = T_r/M_v = {S.OMEGA_DE:.6f} en:')
        print(f'     a* = {a_star:.5f}   con el FONDO REAL')
        print(f'     a* = 1.52103          con la recta CPL')
        print(f'     desplazamiento: {100*(a_star-1.52103)/1.52103:+.2f}%')
        print(f'     w REAL ahi = {np.interp(a_star, a, w):+.6f}   '
              f'1+w = {1+np.interp(a_star, a, w):.6f}')
    else:
        print(f'  Omega_DE NUNCA alcanza {S.OMEGA_DE:.6f} en a <= 3 '
              f'(max = {ODE.max():.6f})')
