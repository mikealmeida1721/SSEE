#!/usr/bin/env python3
"""
T^2(k,z) = P_phiDM(k,z) / P_cold(k,z)  desde Boltzmann, en rejilla de z.

Sin aproximar la supresion como independiente de z: se mide donde KiDS mira.
Mismo fondo en ambas corridas (mismo Omega_m, mismo N_ur) para que el cociente
aisle SOLO el free-streaming del sector phi.
"""
import numpy as np, sys, time
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S
from classy import Class

h = S.H0_ALG / 100
M_PHI, T_PHI = 40.70, 0.386633          # calibrado: Omega_phiDM = 0.14889 exacto
OM_DYN, OM_PHI = 0.160050228655221, 0.14889
OM_B = S.OMEGA_B_H2 / h**2
OM_NU_ACT = (S.SUM_MNU_EV / 93.14) / h**2
OM_CDM = OM_DYN - OM_B - OM_NU_ACT
LOGA = 3.040704
N_UR = 2.9610

ZG = np.array([0.0, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
KH = np.logspace(-4, np.log10(19.5), 300)

BASE = {'h': h, 'omega_b': S.OMEGA_B_H2, 'n_s': S.N_S,
        'A_s': 1e-10 * np.exp(LOGA), 'tau_reio': 0.0544, 'T_cmb': 2.7255,
        'N_ur': N_UR, 'Omega_Lambda': 0.0,
        'fluid_equation_of_state': 'CLP', 'w0_fld': S.W0, 'wa_fld': S.WA,
        'cs2_fld': 1.0, 'output': 'mPk', 'P_k_max_h/Mpc': 25.0,
        'z_max_pk': 3.2}


def pk_grid(extra):
    p = dict(BASE); p.update(extra)
    c = Class(); c.set(p); c.compute()
    out = np.array([[c.pk_lin(k * h, z) * h**3 for k in KH] for z in ZG])
    s8, Om = c.sigma8(), c.Omega_m()
    c.struct_cleanup(); c.empty()
    return out, s8, Om


t0 = time.time()
P_phi, s8_phi, Om_phi = pk_grid(
    {'Omega_cdm': OM_CDM, 'N_ncdm': 2, 'm_ncdm': f'0.06849,{M_PHI}',
     'T_ncdm': f'0.71611,{T_PHI}', 'deg_ncdm': '1.0,1.0'})
P_cold, s8_cold, Om_cold = pk_grid(
    {'Omega_cdm': OM_CDM + OM_PHI, 'N_ncdm': 1, 'm_ncdm': '0.06849',
     'T_ncdm': '0.71611'})

print(f'  phi : Omega_m={Om_phi:.6f}  sigma8={s8_phi:.6f}')
print(f'  cold: Omega_m={Om_cold:.6f}  sigma8={s8_cold:.6f}')
print(f'  cociente sigma8 = {s8_phi/s8_cold:.6f}   [canonico P6: 0.747/0.8334 '
      f'= {0.747/0.8334:.6f}]\n')

T2 = P_phi / P_cold
print(f'  T^2(k,z):   k \\ z ->' + ''.join(f'{z:9.2f}' for z in ZG))
for kt in (0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0):
    i = np.argmin(abs(KH - kt))
    print(f'    k={KH[i]:7.3f}      ' + ''.join(f'{T2[j,i]:9.5f}'
                                                for j in range(len(ZG))))
print(f'\n  dependencia en z (max variacion relativa sobre 0<z<3, k>0.05):')
sel = KH > 0.05
rng = (T2[:, sel].max(axis=0) - T2[:, sel].min(axis=0)) / T2[0, sel]
print(f'    max = {rng.max():.4f}  en k = {KH[sel][np.argmax(rng)]:.3f} h/Mpc')

np.savez('t2_grid.npz', kh=KH, zg=ZG, T2=T2, s8_phi=s8_phi, s8_cold=s8_cold)
print(f'\n  guardado en t2_grid.npz   ({time.time()-t0:.1f}s)')
