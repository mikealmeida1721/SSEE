#!/usr/bin/env python3
"""
COTA INFERIOR A m_phi DESDE LA CIZALLA.

Logica (de Mike): SSEE de un sector con A_s libre ya ajusta KiDS (S8=0.757 vs
0.759). Luego un segundo sector CALIENTE solo puede estropear el ajuste: cambia
la FORMA de P(k), y la forma no la absorbe A_s. Por tanto, si existe un segundo
sector, su borde de free-streaming debe caer FUERA del rango que la cizalla mide.

Fijada la abundancia Omega_phiDM = 0.14889, la masa determina la temperatura
(mas masa -> menos densidad numerica -> mas frio) y por tanto donde muerde el
borde. Este script barre la masa y mide:
   - k_50 : donde T^2 = 0.5      (el borde)
   - dS8  : supresion en sigma8   (con A_s FIJO, para ver el efecto puro)
   - T^2 en k = 0.1 .. 3 h/Mpc    (el rango que domina la cizalla de KiDS)

El criterio de "no estropea": T^2 > 0.99 en todo k < 3 h/Mpc, es decir menos del
1% de distorsion de forma donde estan los datos.
"""
import numpy as np, sys
from scipy.optimize import brentq
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S
from classy import Class

h = S.H0_ALG / 100
OM_DYN, OM_PHI = 0.160050228655221, 0.14889
OM_B = S.OMEGA_B_H2 / h**2
OM_NU_ACT = (S.SUM_MNU_EV / 93.14) / h**2
OM_CDM = OM_DYN - OM_B - OM_NU_ACT
SUM_MNU = S.SUM_MNU_EV

BASE = {'h': h, 'omega_b': S.OMEGA_B_H2, 'n_s': S.N_S,
        'A_s': 1e-10 * np.exp(3.040704), 'tau_reio': 0.0544, 'T_cmb': 2.7255,
        'N_ur': 2.9610, 'Omega_Lambda': 0.0, 'Omega_cdm': OM_CDM,
        'fluid_equation_of_state': 'CLP', 'w0_fld': S.W0, 'wa_fld': S.WA,
        'cs2_fld': 1.0, 'output': 'mPk', 'P_k_max_h/Mpc': 60.0, 'z_max_pk': 1.0}
KH = np.logspace(-3, np.log10(50.), 300)


def corre(m, T, pk=False):
    p = dict(BASE); p.update(
        {'N_ncdm': 2, 'm_ncdm': f'{SUM_MNU},{m}', 'T_ncdm': f'0.71611,{T}',
         'deg_ncdm': '1.0,1.0', 'ncdm_fluid_approximation': 3,
         'l_max_ncdm': 17})
    c = Class(); c.set(p); c.compute()
    r = dict(Onu=c.Omega_nu, s8=c.sigma8(), Om=c.Omega_m())
    if pk:
        r['pk'] = np.array([c.pk_lin(k * h, 0.0) * h**3 for k in KH])
    c.struct_cleanup(); c.empty()
    return r


p = dict(BASE); p['Omega_cdm'] = OM_CDM + OM_PHI
p.update({'N_ncdm': 1, 'm_ncdm': f'{SUM_MNU}', 'T_ncdm': '0.71611'})
c = Class(); c.set(p); c.compute()
S8C, OMC = c.sigma8(), c.Omega_m()
PKC = np.array([c.pk_lin(k * h, 0.0) * h**3 for k in KH])
c.struct_cleanup(); c.empty()
print(f'  control frio: sigma8 = {S8C:.6f}  Omega_m = {OMC:.6f}\n')

KTEST = (0.1, 0.5, 1.0, 3.0, 10.0)
print(f'  {"m_phi(eV)":>10} {"T_ncdm":>9} {"k_50":>8} {"ds8":>8}  ' +
      ' '.join(f'T2(k={k})'.rjust(10) for k in KTEST))
MASAS = [40.70, 100., 300., 1000., 3000., 10000., 30000.]
out = []
for m in MASAS:
    try:
        T = brentq(lambda t: corre(m, t)['Onu'] - OM_NU_ACT - OM_PHI,
                   0.02, 0.60, xtol=1e-7)
        r = corre(m, T, pk=True)
    except Exception as e:
        print(f'  {m:10.1f}  FALLO  {str(e).splitlines()[-1][:60]}')
        continue
    T2 = r['pk'] / PKC
    k50 = KH[np.argmin(abs(T2 - 0.5))]
    ds8 = r['s8'] / S8C - 1
    vals = [T2[np.argmin(abs(KH - k))] for k in KTEST]
    print(f'  {m:10.1f} {T:9.6f} {k50:8.3f} {100*ds8:+7.2f}%  ' +
          ' '.join(f'{v:10.5f}' for v in vals))
    out.append((m, T, k50, ds8, vals))

print(f'\n  CRITERIO "no estropea la cizalla": T^2 > 0.99 para todo k < 3 h/Mpc')
for m, T, k50, ds8, vals in out:
    ok = min(vals[:4]) > 0.99
    print(f'    m = {m:8.1f} eV : min T^2(k<3) = {min(vals[:4]):.5f}  '
          f'{"PASA" if ok else "no pasa"}   (ds8 = {100*ds8:+.2f}%)')
np.savez('cota_masa.npz', masas=[o[0] for o in out], k50=[o[2] for o in out],
         ds8=[o[3] for o in out])
