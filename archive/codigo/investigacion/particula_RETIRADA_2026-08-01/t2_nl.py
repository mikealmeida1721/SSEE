#!/usr/bin/env python3
"""
T2_NL(k,z) = P_phiDM_NOLINEAL / P_cold_NOLINEAL, ambos de CLASS+halofit.

Por que: multiplicar T2 LINEAL por un P(k) NO LINEAL sobre-suprime a k alto,
porque el 52% de materia fria regenera potencia por colapso no lineal y la
supresion lineal no lo sabe. Este cociente si lo sabe (dentro de lo que halofit
puede saber). Es la cota realista; el T2 lineal es la cota pesimista.

CAVEAT: halofit no esta calibrado para una componente caliente al 48%; CLASS lo
avisa. HMcode ni converge. Sin N-body esto es el techo de lo que se puede decir.
"""
import numpy as np, sys, time
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S
from classy import Class

h = S.H0_ALG / 100
M_PHI, T_PHI = 40.70, 0.386633
OM_DYN, OM_PHI = 0.160050228655221, 0.14889
OM_B = S.OMEGA_B_H2 / h**2
OM_NU_ACT = (S.SUM_MNU_EV / 93.14) / h**2
OM_CDM = OM_DYN - OM_B - OM_NU_ACT

ZG = np.array([0.0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
KH = np.logspace(-4, np.log10(19.5), 300)

BASE = {'h': h, 'omega_b': S.OMEGA_B_H2, 'n_s': S.N_S,
        'A_s': 1e-10 * np.exp(3.040704), 'tau_reio': 0.0544, 'T_cmb': 2.7255,
        'N_ur': 2.9610, 'Omega_Lambda': 0.0,
        'fluid_equation_of_state': 'CLP', 'w0_fld': S.W0, 'wa_fld': S.WA,
        'cs2_fld': 1.0, 'output': 'mPk', 'P_k_max_h/Mpc': 25.0,
        'z_max_pk': 3.2, 'non_linear': 'halofit'}


def pk_grid(extra):
    p = dict(BASE); p.update(extra)
    c = Class(); c.set(p); c.compute()
    lin = np.array([[c.pk_lin(k * h, z) * h**3 for k in KH] for z in ZG])
    nl = np.array([[c.pk(k * h, z) * h**3 for k in KH] for z in ZG])
    s8 = c.sigma8()
    c.struct_cleanup(); c.empty()
    return lin, nl, s8


t0 = time.time()
Lp, Np, s8p = pk_grid({'Omega_cdm': OM_CDM, 'N_ncdm': 2,
                       'm_ncdm': f'0.06849,{M_PHI}',
                       'T_ncdm': f'0.71611,{T_PHI}', 'deg_ncdm': '1.0,1.0'})
Lc, Nc, s8c = pk_grid({'Omega_cdm': OM_CDM + OM_PHI, 'N_ncdm': 1,
                       'm_ncdm': '0.06849', 'T_ncdm': '0.71611'})

T2_LIN, T2_NL = Lp / Lc, Np / Nc
print(f'  sigma8 lineal: phi={s8p:.6f} cold={s8c:.6f}  '
      f'cociente={s8p/s8c:.6f}\n')
print(f'  {"k":>8} {"T2 LINEAL":>11} {"T2 NO LIN":>11} {"NL/LIN":>9}   (z=0)')
for kt in (0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 19.0):
    i = np.argmin(abs(KH - kt))
    print(f'  {KH[i]:8.3f} {T2_LIN[0,i]:11.5f} {T2_NL[0,i]:11.5f} '
          f'{T2_NL[0,i]/T2_LIN[0,i]:9.2f}x')
i50l = np.argmin(abs(T2_LIN[0] - 0.5)); i50n = np.argmin(abs(T2_NL[0] - 0.5))
print(f'\n  k al 50% de supresion:  lineal {KH[i50l]:.3f}   '
      f'no lineal {KH[i50n]:.3f} h/Mpc')
np.savez('t2_nl.npz', kh=KH, zg=ZG, T2=T2_NL, T2_lin=T2_LIN,
         s8_phi=s8p, s8_cold=s8c)
print(f'\n  guardado en t2_nl.npz   ({time.time()-t0:.1f}s)')
