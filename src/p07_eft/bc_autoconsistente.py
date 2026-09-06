#!/usr/bin/env python3
"""
beta_c AUTOCONSISTENTE: integrar el fondo CON el acoplamiento realimentado.

POR QUE. `ssee_eft_verification.py` integra el fondo DESACOPLADO (beta_c no
entra ni en la Klein-Gordon ni en la conservacion de la materia oscura) y
luego INFIERE beta_c al final, del Delta_w que haria falta. Eso es orden
dominante. Aqui beta_c entra en las ecuaciones y se resuelve el sistema
completo, que es lo unico que decide si -3.99 es el numero o un artefacto.

Y QUE SE VARIA. El otro ingrediente en disputa es rho_DM,0. El script viejo
usa 0.160050, que es 1+w0 -- un numero de la ECUACION DE ESTADO. Ahora bien,
en SSEE Omega_DE = T_r/M_v = 0.839950 y w0 = -T_r/M_v, o sea Omega_DE = -w0
por construccion, y por tanto 1 - Omega_DE = 1 + w0 = 0.160050 es una
IDENTIDAD. Asi que ese 0.160 tiene dos lecturas y hay que separarlas:

  (A) es 1 - Omega_DE, o sea la materia que hace PLANO al universo dado que
      la energia oscura vale 0.840. Bajo esta lectura SI es una densidad.
  (B) es la densidad de la materia oscura que se acopla al escalar, que
      fisicamente es Omega_c = omega_c/h^2 = 0.2588.

No pueden ser las dos: 0.840 + 0.2588 + Omega_b ya pasa de 1. El programa
corre las combinaciones COHERENTES y reporta beta_c en cada una, para que la
eleccion se haga con los numeros a la vista y no por inercia.

Unidades: Mpl = 1, rho_crit,0 = 1, H0 = 1. N = ln(a).
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S

KAL0 = S.KAL0
W0 = S.W0
OM_DE_ALG = S.OMEGA_DE                       # T_r/M_v = 0.839950
OM_C_FIS = S.OMEGA_C_H2 / (S.H0_GLOBAL / 100.0) ** 2   # omega_c/h^2 = 0.2588
OM_B_FIS = S.OMEGA_B_H2 / (S.H0_GLOBAL / 100.0) ** 2
A_I = 1e-3
AURA = (3 * S.PHI + S.PI) / 2


def integra(phi_i, rho_dm_i, bc, om_de, rho_b0, a_i=A_I):
    """Fondo COMPLETO con el acoplamiento dentro. Devuelve el estado en a=1."""
    V0 = om_de                                # V(0) normalizado a Omega_DE
    lam = np.sqrt(3.0 * (1.0 + W0))           # tracker: usa la EoS (legitimo)
    alpha = lam / np.sqrt(KAL0)

    def rhs(N, Y):
        phi, phip, rho_dm = Y
        a = np.exp(N)
        rho_b = rho_b0 * a ** -3
        V = V0 * np.exp(-alpha * phi)
        # H^2 = (rho_phi + rho_dm + rho_b)/3  con rho_phi = X/KAL + V,
        # X = H^2 phip^2/2  ->  despeje lineal en H^2:
        num = V + rho_dm + rho_b
        den = 3.0 - phip ** 2 / (2.0 * KAL0)
        H2 = num / den
        if H2 <= 0 or den <= 0:
            return [0.0, 0.0, 0.0]
        X = H2 * phip ** 2 / 2.0
        rho_phi = X / KAL0 + V
        rho_tot = rho_phi + rho_dm + rho_b
        p_phi = X / KAL0 - V
        # H'/H = -(3/2)(1 + p_tot/rho_tot), con la materia sin presion
        HpH = -1.5 * (1.0 + p_phi / rho_tot)
        V_phi = -alpha * V
        K_X = 1.0 / KAL0
        phipp = (-(3.0 + HpH) * phip - V_phi / (H2 * K_X)
                 + bc * rho_dm / (H2 * K_X))
        rho_dmp = -3.0 * rho_dm + bc * rho_dm * phip
        return [phip, phipp, rho_dmp]

    sol = solve_ivp(rhs, [np.log(a_i), 0.0], [phi_i, 1e-6, rho_dm_i],
                    rtol=1e-10, atol=1e-14, dense_output=True, method='Radau')
    if not sol.success:
        return None
    phi, phip, rho_dm = sol.y[:, -1]
    V = V0 * np.exp(-alpha * phi)
    rho_b = rho_b0
    num = V + rho_dm + rho_b
    den = 3.0 - phip ** 2 / (2.0 * KAL0)
    H2 = num / den
    X = H2 * phip ** 2 / 2.0
    rho_phi = X / KAL0 + V
    p_phi = X / KAL0 - V
    rho_tot = rho_phi + rho_dm + rho_b
    H = np.sqrt(H2)
    phidot = H * phip
    Q = -bc * rho_dm * phidot                  # transferencia de energia
    w_phi = p_phi / rho_phi
    w_eff = w_phi + Q / (3.0 * H * rho_phi)
    return dict(Om_phi=rho_phi / rho_tot, rho_dm=rho_dm, w_phi=w_phi,
                w_eff=w_eff, H=H, phidot=phidot, rho_phi=rho_phi)


def resuelve(om_de, rho_dm_0, rho_b0, etiqueta):
    """3 incognitas (phi_i, rho_dm_i, bc), 3 condiciones en a=1:
         Omega_phi = om_de ,  rho_dm = rho_dm_0 ,  w_eff = w0."""
    def F(u):
        phi_i, ldm_i, bc = u
        r = integra(phi_i, np.exp(ldm_i), bc, om_de, rho_b0)
        if r is None:
            return [1e3, 1e3, 1e3]
        return [r['Om_phi'] - om_de, r['rho_dm'] - rho_dm_0, r['w_eff'] - W0]

    best = None
    for p0 in (0.5, 1.0, 2.0, 0.1):
        for b0 in (-4.0, -2.5, -1.0):
            u, info, ier, msg = fsolve(F, [p0, np.log(rho_dm_0 * A_I ** -3), b0],
                                       full_output=True)
            if ier == 1 and max(abs(np.array(F(u)))) < 1e-7:
                best = u
                break
        if best is not None:
            break
    if best is None:
        print(f'  {etiqueta:52s}  NO CONVERGE')
        return None
    bc = best[2]
    r = integra(best[0], np.exp(best[1]), bc, om_de, rho_b0)
    M = abs(bc) / 2.0
    aK = 3 * AURA * (S.PI - S.PHI) / (2 * (S.PI + S.PHI) ** 2)
    f = aK / (3 * M)
    HL = S.H0_GLOBAL / (1 - f)
    print(f'  {etiqueta:52s}  bc={bc:+9.5f}  vs -AURA {100*(abs(bc)-AURA)/AURA:+7.2f}%'
          f'   MIRA={M:7.5f}  f_scr={f:.5f}  H_loc={HL:7.3f}'
          f'  ({abs(HL-73.04)/1.04:.2f}sig)')
    return bc


if __name__ == '__main__':
    print('=' * 118)
    print('  beta_c AUTOCONSISTENTE — el acoplamiento DENTRO de las ecuaciones')
    print('=' * 118)
    print(f'  -AURA = {-AURA:.6f}   Omega_DE algebraico = {OM_DE_ALG:.6f}   '
          f'Omega_c fisico = {OM_C_FIS:.6f}   Omega_b = {OM_B_FIS:.6f}')
    print(f'  1 + w0 = {1+W0:.6f}  =  1 - Omega_DE  (IDENTIDAD en SSEE, no coincidencia)')
    print()
    print('  Combinaciones COHERENTES (Omega suman 1 salvo donde se indica):')
    resuelve(OM_DE_ALG, 1 + W0, 0.0,
             'A) Om_DE=0.840, rho_DM=0.160 (=1-Om_DE), sin bariones')
    resuelve(OM_DE_ALG, 1 + W0 - OM_B_FIS, OM_B_FIS,
             'B) Om_DE=0.840, rho_DM=0.160-Om_b, con bariones')
    resuelve(1 - S.OMEGA_M_TOTAL, OM_C_FIS, OM_B_FIS,
             'C) Om_DE=1-Om_m=0.691 (geometria real), rho_DM=Om_c')
    resuelve(OM_DE_ALG, OM_C_FIS, OM_B_FIS,
             'D) Om_DE=0.840 y rho_DM=Om_c  (NO cierra: suma 1.149)')
