#!/usr/bin/env python3
"""
¿Cuanta materia CREE ver BAO cuando la materia esta ACOPLADA al campo?

PREMISA DISUELTA 2026-09-05 — LEER ANTES DE CORRER ESTE SCRIPT.
El script nacio para reconciliar "dos numeros de materia que no cuadran".
NO HAY DOS NUMEROS DE MATERIA. Solo hay uno: Omega_m = omega_m/h^2 = 0.308881.
El 0.160 NO es una densidad de materia: es 1+w0 = s_M, la saturacion
complementaria, un numero de la ECUACION DE ESTADO. Medido en la corrida del
fondo (results/logs/fondo_acoplado.npz): lo que satura en 0.8399 es w, no una
densidad; Omega_DE vale 0.691119 hoy y tiende a 1 en el futuro, pasando de
largo por 0.839950 sin detenerse. Ver ssee_core.py y las reglas R52/R52b.
Asi que no habia contradiccion que reconciliar: se comparaba una densidad
medida con un numero de la ecuacion de estado. Es el MISMO error de categoria
que definia la particula phi-DM retirada el 2026-08-01.
Lo que SI quedo abierto sobre beta_c esta en OP-23 (OPEN_PROBLEMS.md), y es
otra cosa: dentro del acoplamiento conformal no existe beta_c que reproduzca
w0 sin pasarse del limite de energia oscura temprana.
Se conserva el script por trazabilidad; su pregunta ya no esta viva.

LA PREGUNTA ORIGINAL (historica, premisa falsa). SSEE tendria dos numeros de
materia que no cuadran: los ingredientes dan Omega_m = omega_m/h^2 = 0.308881,
y la geometria vieja usaba 0.160 (= 1 - T_r/M_v). [PREMISA FALSA: 0.160 es
1+w0, ecuacion de estado, no densidad -- ver arriba, R52 y OP-23.]
La hipotesis, historica y ya no viva: con acoplamiento beta_c != 0 la
materia oscura NO diluye como a^-3 -- el campo le bombea energia --, asi que el
Omega_m que un ajuste ESTANDAR infiere de las distancias no tiene por que ser
el Omega_m real. Si el acoplamiento desplaza lo suficiente, 0.160 y 0.308881
dejan de contradecirse: uno es lo que hay, el otro lo que BAO cree ver.

COMO SE MIDE. Sin elegir nada:
  1. Se integra el fondo COMPLETO con el acoplamiento vivo, partiendo de las
     densidades REALES de hoy (Omega_b, Omega_c de los ingredientes).
  2. De ahi sale E(z) exacto del modelo acoplado.
  3. Se ajusta a ese E(z) la parametrizacion ESTANDAR que usa el analisis BAO,
     E^2 = Om_eff (1+z)^3 + (1-Om_eff) f_CPL(z; w0,wa), dejando SOLO Om_eff
     libre. Eso es literalmente "que Omega_m infiere BAO".
  4. Se compara Om_eff contra el Omega_m real y contra 0.160.

CONTROL. La misma corrida con beta_c = 0 debe devolver Om_eff = Omega_m real.
Si el control falla, el metodo esta mal y el resultado no vale.

Unidades: Mpl = 1, rho_crit,0 = 1, H0 = 1, N = ln(a).
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S

KAL = S.KAL0
W0, WA = S.W0, S.WA
AURA = (3 * S.PHI + S.PI) / 2
h = S.H0_GLOBAL / 100.0
OM_B = S.OMEGA_B_H2 / h ** 2
OM_NU = S.OMEGA_NU_H2 / h ** 2
OM_C = S.OMEGA_C_H2 / h ** 2
OM_M = S.OMEGA_M_TOTAL                    # 0.308881, los ingredientes
OM_PHI = 1.0 - OM_M                       # 0.691119, por planitud
LAM = np.sqrt(3.0 * (1.0 + W0))
ALPHA = LAM / np.sqrt(KAL)


def fondo_acoplado(bc, z_max=2.5, n=400):
    """Integra hacia atras desde hoy. Devuelve (z, E(z)) del modelo acoplado.

    Hoy:  rho_phi = OM_PHI  y  w_phi = W0  fijan X0 y V0;
          rho_c = OM_C (acoplada),  rho_b + rho_nu = OM_B + OM_NU (libres).
    """
    X0_over_KAL = OM_PHI * (1.0 + W0) / 2.0     # = K, densidad cinetica
    V0 = OM_PHI * (1.0 - W0) / 2.0              # = V, densidad potencial
    phip0 = np.sqrt(2.0 * KAL * X0_over_KAL)    # X = phip^2/2 (H0=1)

    def rhs(N, Y):
        phi, phip, rho_c = Y
        a = np.exp(N)
        rho_bn = (OM_B + OM_NU) * a ** -3
        V = V0 * np.exp(-ALPHA * phi)
        den = 1.0 - phip ** 2 / (2.0 * KAL)
        H2 = (V + rho_c + rho_bn) / den
        if H2 <= 0 or den <= 0:
            return [0.0, 0.0, 0.0]
        X = H2 * phip ** 2 / 2.0
        rho_phi = X / KAL + V
        p_phi = X / KAL - V
        rho_tot = rho_phi + rho_c + rho_bn
        HpH = -1.5 * (1.0 + p_phi / rho_tot)   # independiente de unidades
        phipp = (-(3.0 + HpH) * phip + ALPHA * V / (H2 / KAL)
                 + bc * rho_c / (H2 / KAL))
        rho_cp = -3.0 * rho_c + bc * rho_c * phip
        return [phip, phipp, rho_cp]

    N_min = np.log(1.0 / (1.0 + z_max))
    sol = solve_ivp(rhs, [0.0, N_min], [0.0, phip0, OM_C],
                    t_eval=np.linspace(0.0, N_min, n),
                    rtol=1e-10, atol=1e-13, method='Radau')
    if not sol.success:
        return None
    a = np.exp(sol.t)
    phi, phip, rho_c = sol.y
    rho_bn = (OM_B + OM_NU) * a ** -3
    V = V0 * np.exp(-ALPHA * phi)
    den = 1.0 - phip ** 2 / (2.0 * KAL)
    H2 = (V + rho_c + rho_bn) / den
    return 1.0 / a - 1.0, np.sqrt(H2), rho_c


def f_cpl(z):
    a = 1.0 / (1.0 + z)
    return (1 + z) ** (3 * (1 + W0 + WA)) * np.exp(-3 * WA * (1 - a))


def ajusta_om_eff(z, E):
    """El Omega_m que un analisis BAO ESTANDAR inferiria de este E(z)."""
    def chi2(om):
        Em = np.sqrt(om * (1 + z) ** 3 + (1 - om) * f_cpl(z))
        return float(np.sum((Em - E) ** 2))
    r = minimize_scalar(chi2, bounds=(0.01, 0.99), method='bounded',
                        options=dict(xatol=1e-10))
    return float(r.x), float(r.fun)


if __name__ == '__main__':
    print('=' * 96)
    print('  ¿Que Omega_m infiere BAO cuando la materia esta ACOPLADA?')
    print('=' * 96)
    print(f'  Ingredientes:  Omega_b={OM_B:.6f}  Omega_c={OM_C:.6f}  '
          f'Omega_nu={OM_NU:.6f}  ->  Omega_m REAL = {OM_M:.6f}')
    print(f'  Por planitud:  Omega_phi = {OM_PHI:.6f}       w0 = {W0:.6f}')
    print(f'  El numero en disputa: 1 - T_r/M_v = {1-S.OMEGA_DE:.6f}')
    print()
    print(f'  {"beta_c":>10}  {"Om_eff (lo que BAO cree ver)":>30}  '
          f'{"desvio vs real":>16}  {"residuo":>10}')
    print('  ' + '-' * 74)
    for bc in (0.0, -1.0, -2.0, -AURA, -3.9899, -6.0, -10.0):
        r = fondo_acoplado(bc)
        if r is None:
            print(f'  {bc:>10.4f}  {"NO INTEGRA":>30}')
            continue
        z, E, rho_c = r
        om, res = ajusta_om_eff(z, E)
        etq = ' <- CONTROL: debe dar el real' if bc == 0.0 else ''
        print(f'  {bc:>10.4f}  {om:>30.6f}  {100*(om-OM_M)/OM_M:>15.2f}%  '
              f'{res:>10.2e}{etq}')
    print()
    print(f'  Para que BAO viera 0.160050 haria falta un desvio de '
          f'{100*(0.160050-OM_M)/OM_M:.1f}%')
