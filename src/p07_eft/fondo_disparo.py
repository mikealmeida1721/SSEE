#!/usr/bin/env python3
"""
Fondo de SSEE integrado con DISPARO desde epoca temprana. Pasado Y futuro.

POR QUE ASI Y NO COMO ANTES. El intento previo (fondo_real_futuro.py) fijaba las
condiciones HOY y de ahi integraba. Hacia adelante funcionaba; hacia atras daba
resultados imposibles (Omega_DE creciendo al retroceder, w positivo). Eso no es
una propiedad del universo: es que el potencial exponencial es un ATRACTOR hacia
adelante, asi que retroceder amplifica cualquier error. La consecuencia practica
es grave -- el pasado es justo donde ESTAN los datos (DESI mide a = 0.30..0.77),
asi que un metodo que solo mira al futuro no se puede contrastar con nada.

El metodo correcto, que es el que usan CLASS y CAMB: arrancar TEMPRANO y disparar
las condiciones iniciales hasta aterrizar en los valores de hoy. Toda la
trayectoria va con el atractor, asi que pasado y futuro salen igual de estables.

INCOGNITAS (2):   phi(a_i)  y  phi'(a_i)
CONDICIONES (2):  en a = 1,   Omega_DE = 1 - omega_m/h0^2 = 0.691119
                              w_phi    = w0 = -0.839950

Nada mas se ajusta: omega_m viene de los ingredientes, H0 del ancla, y el
potencial V = V0 exp(-alpha phi) tiene alpha fijado por la relacion de tracker
alpha = sqrt(3)*sqrt(3(1+w0))/sqrt(KAL). V0 se lleva junto con phi_i (degenerados:
solo importa V0*exp(-alpha*phi_i)), asi que se fija V0 = 1 y se dispara phi_i.

Unidades: rho_crit,0 = 1, H0 = 1  =>  H^2 = rho_tot.  N = ln a.
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S

KAL = S.KAL0
W0 = S.W0
OM_M = S.OMEGA_M_TOTAL                 # 0.308881, de los ingredientes
OM_PHI = 1.0 - OM_M                    # 0.691119, por planitud
LAM = np.sqrt(3.0 * (1.0 + W0))        # lambda del atractor, en unidades de M_Pl
# FIX 2026-08-10 (falta un sqrt(3)). El exponente del potencial es lambda*phi_c/M_Pl,
# y en estas unidades (H^2 = rho, o sea 3 M_Pl^2 = 1) se tiene M_Pl = 1/sqrt(3).
# El campo canonico es phi_c = phi/sqrt(KAL), asi que ALPHA = sqrt(3)*LAM/sqrt(KAL).
# CONTROL: sin el sqrt(3) la trayectoria asintotaba a w = -0.94665; con el,
# asintota a -0.839965 = w0 (dif 1.5e-5). w0 es el atractor, no un ajuste.
ALPHA = np.sqrt(3.0) * LAM / np.sqrt(KAL)
A_I = 1.0e-3                           # arranque: z = 999
V0 = 1.0                               # degenerado con phi_i; se fija


def _campo(phi, phip, a):
    """Densidades del campo y de la materia en un punto."""
    rho_m = OM_M * a ** -3
    V = V0 * np.exp(-ALPHA * phi)
    den = 1.0 - phip ** 2 / (2.0 * KAL)
    H2 = (V + rho_m) / den
    X = H2 * phip ** 2 / 2.0
    return rho_m, V, H2, X


def rhs(N, Y):
    phi, phip = Y
    a = np.exp(N)
    rho_m, V, H2, X = _campo(phi, phip, a)
    if H2 <= 0:
        return [0.0, 0.0]
    p_phi = X / KAL - V
    rho_tot = X / KAL + V + rho_m
    HpH = -1.5 * (1.0 + p_phi / rho_tot)
    phipp = -(3.0 + HpH) * phip + ALPHA * V * KAL / H2
    return [phip, phipp]


def corre(phi_i, phip_i, a_fin=3.0, n=3000):
    sol = solve_ivp(rhs, [np.log(A_I), np.log(a_fin)], [phi_i, phip_i],
                    t_eval=np.linspace(np.log(A_I), np.log(a_fin), n),
                    rtol=1e-11, atol=1e-14, method='Radau')
    if not sol.success:
        return None
    a = np.exp(sol.t)
    phi, phip = sol.y
    rho_m, V, H2, X = _campo(phi, phip, a)
    rho_phi = X / KAL + V
    p_phi = X / KAL - V
    return dict(a=a, H=np.sqrt(H2), rho_phi=rho_phi, rho_m=rho_m,
                w=p_phi / rho_phi, Om_DE=rho_phi / (rho_phi + rho_m),
                K=X / KAL, V=V)


# ORIGEN de los numeros (R65, 2026-09-19)
# ORIGEN-VALOR: 0.9987 — fila de la rejilla de integracion mas cercana a a=1 (muestreo, no fisica)
# ORIGEN-VALOR: 0.971200 — w_phi de Paper 7 = -0.971202 (manuscript/SSEE_Unified_Journal.tex:716), aqui truncado a 4 cifras
# ORIGEN-VALOR: 0.9712 — el mismo w_phi = -0.971202 a 4 decimales
# ORIGEN-VALOR: 1.52103 — a* con la recta CPL: Omega_DE(a) = T_r/M_v, sale 1.521028 (src/p07_eft/fondo_real_futuro.py)
# ORIGEN-VALOR: 1.3769 — a* con el FONDO REAL = 1.37680 (src/p07_eft/fondo_real_futuro.py); aqui solo es punto de muestreo de la tabla (se toma la fila mas cercana)
def _en_a1(r):
    """Interpolado a a = 1 EXACTO. Leer la fila mas cercana de la rejilla
    (a = 0.9987) metia 1.5e-3 de error puramente de muestreo y rompia la
    normalizacion rho_tot(a=1) = 1."""
    return (float(np.interp(1.0, r['a'], r['Om_DE'])),
            float(np.interp(1.0, r['a'], r['w'])))


def objetivo(u):
    """UNA sola condicion: Omega_DE(a=1) = 0.691119.

    NO se impone w(a=1) = w0. Razon medida: w0 = -0.839950 coincide a precision
    de maquina con el ATRACTOR del potencial, -1 + lambda^2/3. Es donde el campo
    TERMINA, no donde esta hoy. Exigirle ese w hoy, con Omega_phi = 0.691 (aun no
    domina), es pedirle un estado inconsistente -- y por eso el disparo con dos
    condiciones no converge. Asi que w(a=1) sale como PREDICCION.

    El campo arranca congelado por friccion de Hubble (phi'_i ~ 0), asi que la
    unica incognita real es phi_i.
    """
    r = corre(u[0], 1e-8, a_fin=1.05, n=600)
    if r is None:
        return 1e3
    ode, w = _en_a1(r)
    return ode - OM_PHI


if __name__ == '__main__':
    print('=' * 88)
    print('  FONDO SSEE con DISPARO desde a_i = 1e-3 (z = 999)')
    print('=' * 88)
    print(f'  alpha = {ALPHA:.6f}   objetivo en a=1: Om_DE = {OM_PHI:.6f}, '
          f'w = {W0:.6f}')
    from scipy.optimize import brentq
    sol = None
    try:
        phi_i = brentq(lambda p: objetivo([p]), -40.0, 60.0, xtol=1e-13)
        sol = [phi_i, 1e-8]
    except Exception as e:
        print('  brentq:', e)
    if sol is None:
        print('  EL DISPARO NO CONVERGE — no se reporta ningun numero.')
        sys.exit(1)
    print(f"  disparo: phi(a_i) = {sol[0]:.8f}   phi'(a_i) = {sol[1]:.1e} (congelado)")
    r = corre(sol[0], sol[1])
    ode1, w1 = _en_a1(r)
    print(f'  CONTROL en a=1:  Om_DE = {ode1:.9f} (obj {OM_PHI:.9f}, '
          f'dif {abs(ode1-OM_PHI):.2e})')
    print(f'  PREDICCION w(a=1) = {w1:.6f}   (NO impuesto)')
    print(f'     w0 = atractor    = {W0:.6f}   dif = {w1-W0:+.6f}')
    print(f'     Paper 7 numerico = -0.971200   dif = {w1-(-0.9712):+.6f}')
    print()
    a, w, ode, H = r['a'], r['w'], r['Om_DE'], r['H']
    print(f'  {"a":>7} {"z":>9} {"H [km/s/Mpc]":>13} {"Om_m":>10} {"Om_DE":>10} '
          f'{"suma":>10} {"w(a)":>10} {"1+w":>9}')
    for at in (0.001, 0.01, 0.1, 0.30, 0.50, 0.77, 1.0, 1.3769, 1.52103, 2.0, 3.0):
        i = int(np.argmin(abs(a - at)))
        print(f'  {a[i]:>7.4f} {1/a[i]-1:>+9.3f} {S.H0_GLOBAL*H[i]:>13.4f} '
              f'{1-ode[i]:>10.6f} {ode[i]:>10.6f} {1.0:>10.7f} '
              f'{w[i]:>+10.6f} {1+w[i]:>9.6f}')
    np.savez('/home/mike/Proyectos/SSEE/results/logs/fondo_disparo.npz',
             a=a, w=w, Om_DE=ode, H=H, K=r['K'], V=r['V'],
             rho_phi=r['rho_phi'], rho_m=r['rho_m'])
    print('\n  -> results/logs/fondo_disparo.npz')
