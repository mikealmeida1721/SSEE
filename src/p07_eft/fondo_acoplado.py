#!/usr/bin/env python3
"""
Fondo de SSEE ACOPLADO, con disparo desde epoca temprana. Pasado y futuro.

QUE ANADE respecto de fondo_disparo.py: el acoplamiento conformal beta_c entra
DENTRO de las ecuaciones -- en la Klein-Gordon y en la conservacion de la
materia oscura -- y se resuelve autoconsistentemente, en vez de inferirse al
final del Delta_w que haria falta (que es lo que hace el script viejo de
Paper 7 y es solo orden dominante).

DOS CORRECCIONES DE PRECISION respecto del intento anterior:
  1. Las condiciones en a=1 se evaluan por INTERPOLACION a a=1 exacto, no
     leyendo la fila mas cercana de la rejilla (a=0.9987 daba Omega_DE=0.6888
     en vez de 0.691119: 2.3e-3 de error puramente de muestreo).
  2. Los bariones van aparte y SIN acoplar; solo la materia oscura se acopla,
     que es lo que dice el Lagrangiano (L_int = [C(phi)-1] L_DM).

INCOGNITAS (3):  phi(a_i),  rho_c(a_i),  beta_c
CONDICIONES (3) en a=1, todas de los ingredientes o del algebra, ninguna ajustada:
     Omega_DE = 1 - omega_m/h0^2 = 0.691119
     rho_c    = omega_c/h0^2     = 0.258753
     w_eff    = w0               = -0.839950     <- w EFECTIVO, no el del campo

w_eff = w_phi + Q/(3 H rho_phi)  con  Q = -beta_c * rho_DM * phidot.
Esa es la cantidad que ajusta DESI; el campo desnudo tiene otro w (la corrida
sin acoplar predice -0.9761, y Paper 7 da -0.9712 por otro camino).

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
h0 = S.H0_GLOBAL / 100.0
OM_C = S.OMEGA_C_H2 / h0 ** 2                 # 0.258753 materia oscura
OM_BN = (S.OMEGA_B_H2 + S.OMEGA_NU_H2) / h0 ** 2   # bariones + neutrinos
OM_M = S.OMEGA_M_TOTAL
OM_PHI = 1.0 - OM_M
# FIX 2026-08-10: faltaba un sqrt(3). Mismo motivo que en fondo_disparo.py — el
# exponente va como lambda*phi_c/M_Pl y aqui M_Pl = 1/sqrt(3) (porque H^2 = rho).
ALPHA = np.sqrt(3.0) * np.sqrt(3.0 * (1.0 + W0)) / np.sqrt(KAL)
A_I = 1.0e-3
V0 = 1.0
AURA = (3 * S.PHI + S.PI) / 2


def _estado(phi, phip, rho_c, a):
    rho_bn = OM_BN * a ** -3
    V = V0 * np.exp(-ALPHA * phi)
    den = 1.0 - phip ** 2 / (2.0 * KAL)
    H2 = (V + rho_c + rho_bn) / den
    X = H2 * phip ** 2 / 2.0
    return rho_bn, V, H2, X


def rhs(N, Y, bc):
    phi, phip, rho_c = Y
    a = np.exp(N)
    rho_bn, V, H2, X = _estado(phi, phip, rho_c, a)
    if H2 <= 0:
        return [0.0, 0.0, 0.0]
    p_phi = X / KAL - V
    rho_tot = X / KAL + V + rho_c + rho_bn
    HpH = -1.5 * (1.0 + p_phi / rho_tot)
    phipp = (-(3.0 + HpH) * phip + ALPHA * V * KAL / H2
             + bc * rho_c * KAL / H2)
    rho_cp = -3.0 * rho_c + bc * rho_c * phip
    return [phip, phipp, rho_cp]


def corre(phi_i, rho_c_i, bc, a_fin=3.0, n=3000):
    sol = solve_ivp(rhs, [np.log(A_I), np.log(a_fin)], [phi_i, 1e-8, rho_c_i],
                    args=(bc,), t_eval=np.linspace(np.log(A_I), np.log(a_fin), n),
                    rtol=1e-11, atol=1e-14, method='Radau')
    if not sol.success:
        return None
    a = np.exp(sol.t)
    phi, phip, rho_c = sol.y
    rho_bn, V, H2, X = _estado(phi, phip, rho_c, a)
    H = np.sqrt(H2)
    rho_phi = X / KAL + V
    p_phi = X / KAL - V
    w_phi = p_phi / rho_phi
    Q = -bc * rho_c * (H * phip)                 # transferencia de energia
    w_eff = w_phi + Q / (3.0 * H * rho_phi)
    tot = rho_phi + rho_c + rho_bn
    return dict(a=a, H=H, rho_phi=rho_phi, rho_c=rho_c, rho_bn=rho_bn,
                w_phi=w_phi, w_eff=w_eff, Om_DE=rho_phi / tot,
                Om_m=(rho_c + rho_bn) / tot, K=X / KAL, V=V, Q=Q)


def _en_a1(r, clave):
    """Interpolado a a = 1 EXACTO, no la fila mas cercana de la rejilla."""
    return float(np.interp(1.0, r['a'], r[clave]))


def objetivo(u):
    phi_i, l_rho, bc = u
    r = corre(phi_i, np.exp(l_rho), bc, a_fin=1.05, n=900)
    if r is None:
        return [1e3, 1e3, 1e3]
    return [_en_a1(r, 'Om_DE') - OM_PHI,
            _en_a1(r, 'rho_c') - OM_C,
            _en_a1(r, 'w_eff') - W0]


if __name__ == '__main__':
    print('=' * 94)
    print('  FONDO SSEE ACOPLADO — beta_c autoconsistente, disparo desde z=999')
    print('=' * 94)
    print(f'  objetivos en a=1 (interpolados exactos):  Om_DE={OM_PHI:.6f}  '
          f'rho_c={OM_C:.6f}  w_eff={W0:.6f}')
    l0 = np.log(OM_C * A_I ** -3)
    sol = None
    for p0 in (1.158, 0.0, 3.0, 8.0, 20.0):
        for b0 in (0.0, -0.5, -1.0, -2.0, -3.9978, -6.0):
            u, info, ier, msg = fsolve(objetivo, [p0, l0, b0],
                                       full_output=True, xtol=1e-12)
            if ier == 1 and max(abs(np.array(objetivo(u)))) < 1e-8:
                sol = u
                break
        if sol is not None:
            break
    if sol is None:
        print('  EL DISPARO NO CONVERGE — no se reporta ningun numero.')
        sys.exit(1)
    bc = sol[2]
    r = corre(sol[0], np.exp(sol[1]), bc)
    print(f"  disparo: phi_i={sol[0]:.8f}  rho_c_i={np.exp(sol[1]):.6e}  "
          f"beta_c={bc:+.6f}")
    print(f"  CONTROL a=1:  Om_DE={_en_a1(r,'Om_DE'):.9f}  "
          f"rho_c={_en_a1(r,'rho_c'):.9f}  w_eff={_en_a1(r,'w_eff'):.9f}")
    print()
    print(f'  beta_c AUTOCONSISTENTE = {bc:+.6f}')
    print(f'     -AURA               = {-AURA:+.6f}   '
          f'brecha {100*(abs(bc)-AURA)/AURA:+.2f}%')
    print(f'     Paper 7 (orden dom.)= -3.989910')
    print()
    MIRA_bc = abs(bc) / 2.0
    aK1 = 3 * _en_a1(r, 'Om_DE') * (1 + _en_a1(r, 'w_eff'))
    print(f'  alpha_K(a=1) = 3*Om_DE*(1+w_eff) = {aK1:.6f}   '
          f'(canonico papers: 0.403302)')
    print(f'  MIRA = |beta_c|/2 = {MIRA_bc:.6f}   (canonico: {S.MIRA:.6f})')
    f1 = aK1 / (3 * MIRA_bc)
    print(f'  f_screen(a=1) = {f1:.6f}   ->  H_local = '
          f'{S.H0_GLOBAL/(1-f1):.4f}   (SH0ES 73.04+-1.04: '
          f'{abs(S.H0_GLOBAL/(1-f1)-73.04)/1.04:.2f} sigma)')
    print()
    a = r['a']
    print(f'  {"a":>7} {"z":>9} {"H":>11} {"Om_m":>10} {"Om_DE":>10} {"suma":>10} '
          f'{"w_phi":>10} {"w_eff":>10} {"1+w_eff":>9}')
    for at in (0.001, 0.1, 0.30, 0.50, 0.77, 1.0, 1.5, 2.0, 3.0):
        i = int(np.argmin(abs(a - at)))
        print(f'  {a[i]:>7.4f} {1/a[i]-1:>+9.3f} {S.H0_GLOBAL*r["H"][i]:>11.4f} '
              f'{r["Om_m"][i]:>10.6f} {r["Om_DE"][i]:>10.6f} '
              f'{r["Om_m"][i]+r["Om_DE"][i]:>10.7f} {r["w_phi"][i]:>+10.6f} '
              f'{r["w_eff"][i]:>+10.6f} {1+r["w_eff"][i]:>9.6f}')
    np.savez('/home/mike/Proyectos/SSEE/results/logs/fondo_acoplado.npz',
             a=a, w_phi=r['w_phi'], w_eff=r['w_eff'], Om_DE=r['Om_DE'],
             Om_m=r['Om_m'], H=r['H'], K=r['K'], V=r['V'], beta_c=bc)
    print('\n  -> results/logs/fondo_acoplado.npz')
