#!/usr/bin/env python3
"""
phi-DM en CLASS con la particion CORRECTA y T_ncdm calibrada.

Particion del modelo (Paper 6):
  Omega_m,dyn = 0.160          <- sector dinamico TOTAL (bariones + CDM frio)
  Omega_phiDM = 0.14889        <- sector phi, solo k < k_fs
  suma        = 0.308881       = Omega_m,CMB  (= omega_m/h^2)
=> Omega_cdm = 0.160 - Omega_b - Omega_nu_activos
"""
import numpy as np, sys, time
from scipy.optimize import brentq
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S
from classy import Class

h = S.H0_ALG / 100
M_PHI = 40.70
OM_DYN, OM_PHI = 0.160050228655221, 0.14889
OM_B = S.OMEGA_B_H2 / h**2
OM_NU_ACT = (S.SUM_MNU_EV / 93.14) / h**2
OM_CDM = OM_DYN - OM_B - OM_NU_ACT
LOGA = 3.040704

print(f'  Omega_b        = {OM_B:.6f}')
print(f'  Omega_nu(act)  = {OM_NU_ACT:.6f}')
print(f'  Omega_cdm      = {OM_DYN:.6f} - {OM_B:.6f} - {OM_NU_ACT:.6f} '
      f'= {OM_CDM:.6f}')
print(f'  Omega_phiDM    = {OM_PHI:.6f}   (via ncdm, m={M_PHI} eV)')
print(f'  suma esperada  = {OM_DYN + OM_PHI:.6f}   [Omega_m,CMB = 0.308881]\n')


def corre(T_ncdm, non_linear=None, want_pk=False):
    p = {'h': h, 'omega_b': S.OMEGA_B_H2, 'n_s': S.N_S,
         'A_s': 1e-10 * np.exp(LOGA), 'tau_reio': 0.0544, 'T_cmb': 2.7255,
         'Omega_cdm': OM_CDM, 'N_ur': 2.9610, 'N_ncdm': 2,
         'm_ncdm': f'0.06849,{M_PHI}', 'T_ncdm': f'0.71611,{T_ncdm}',
         'deg_ncdm': '1.0,1.0', 'Omega_Lambda': 0.0,
         'fluid_equation_of_state': 'CLP', 'w0_fld': S.W0, 'wa_fld': S.WA,
         'cs2_fld': 1.0, 'output': 'mPk', 'P_k_max_h/Mpc': 20.0,
         'z_max_pk': 2.0}
    if non_linear:
        p['non_linear'] = non_linear
    c = Class(); c.set(p); c.compute()
    r = dict(Om=c.Omega_m(), Onu=c.Omega_nu, s8=c.sigma8())
    if want_pk:
        kh = np.logspace(-3, np.log10(19.), 500)
        r['kh'] = kh
        r['pk'] = np.array([c.pk(k * h, 0.0) * h**3 for k in kh])
    c.struct_cleanup(); c.empty()
    return r


# calibrar T_ncdm para que Omega_ncdm(phi) = OM_PHI exacto
def resid(T):
    return corre(T)['Onu'] - OM_NU_ACT - OM_PHI


t0 = time.time()
T_cal = brentq(resid, 0.30, 0.48, xtol=1e-6)
r = corre(T_cal, want_pk=True)
print(f'  T_ncdm calibrado = {T_cal:.6f}   (mi estimacion analitica: 0.385363)')
print(f'  Omega_ncdm total = {r["Onu"]:.6f}  ->  phi = {r["Onu"]-OM_NU_ACT:.6f}'
      f'   [objetivo {OM_PHI}]')
print(f'  Omega_m (CLASS)  = {r["Om"]:.6f}   [objetivo 0.308881]')
print(f'  sigma8 LINEAL    = {r["s8"]:.6f}')
print(f'  calibracion en {time.time()-t0:.1f}s\n')

# control: mismo fondo pero SIN free-streaming (phi como CDM frio)
p = {'h': h, 'omega_b': S.OMEGA_B_H2, 'n_s': S.N_S,
     'A_s': 1e-10 * np.exp(LOGA), 'tau_reio': 0.0544, 'T_cmb': 2.7255,
     'Omega_cdm': OM_CDM + OM_PHI, 'N_ur': 2.9610, 'N_ncdm': 1,
     'm_ncdm': '0.06849', 'T_ncdm': '0.71611', 'Omega_Lambda': 0.0,
     'fluid_equation_of_state': 'CLP', 'w0_fld': S.W0, 'wa_fld': S.WA,
     'cs2_fld': 1.0, 'output': 'mPk', 'P_k_max_h/Mpc': 20.0, 'z_max_pk': 2.0}
c = Class(); c.set(p); c.compute()
kh = r['kh']
pk_cold = np.array([c.pk(k * h, 0.0) * h**3 for k in kh])
s8_cold = c.sigma8(); Om_cold = c.Omega_m()
c.struct_cleanup(); c.empty()

print(f'  CONTROL (phi como CDM frio, mismo Omega_m):')
print(f'    Omega_m = {Om_cold:.6f}   sigma8 = {s8_cold:.6f}')
print(f'\n  supresion LINEAL por free-streaming:  sigma8 {s8_cold:.4f} -> '
      f'{r["s8"]:.4f}   ({100*(r["s8"]/s8_cold-1):+.2f}%)')

sup = r['pk'] / pk_cold
print(f'\n  T^2(k) = P_phiDM/P_cold  (supresion de Boltzmann, NO ajustada):')
for kt in (0.05, 0.1, 0.2, 0.35, 0.5, 0.754, 1.0, 2.0, 5.0, 10.0):
    i = np.argmin(abs(kh - kt))
    print(f'    k = {kh[i]:6.3f} h/Mpc :  {sup[i]:.6f}')
i50 = np.argmin(abs(sup - 0.5))
print(f'\n  k donde la supresion llega al 50%: {kh[i50]:.4f} h/Mpc'
      f'   [k_fs analitico del Paper 6 = 0.754]')
np.savez('class_phidm.npz', kh=kh, pk_phi=r['pk'], pk_cold=pk_cold,
         T_ncdm=T_cal, s8_phi=r['s8'], s8_cold=s8_cold, Om=r['Om'])
print(f'\n  guardado en class_phidm.npz')
