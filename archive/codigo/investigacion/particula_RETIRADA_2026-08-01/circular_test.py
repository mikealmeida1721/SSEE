#!/usr/bin/env python3
"""
TEST DE CIRCULARIDAD: se eligio m_phi para cerrar S8?

Si la respuesta fuera SI, entonces S8 tiene que DISCRIMINAR entre los candidatos
del multiplicador. Los tres que sobrevivieron a la gramatica de linaje (OP-9,
1/3) son:
    SOLAR^2 * OMEGA     = 297   -> m = 20.34 eV
    SOLAR^2 * PYROS     = 398   -> m = 27.27 eV
    SOLAR^2 * KRYSTOS_V = 594   -> m = 40.70 eV   <- el adoptado
mas los dos retirados historicamente: 535.28 (excluido por forma g^2v) y 615.33.

Omega_phiDM = 0.14889 se mantiene FIJO en todos (es la diferencia omega_m/h^2 -
0.160, no una perilla): para cada masa se recalibra T_ncdm para que la abundancia
salga igual. Asi lo unico que cambia es DONDE muerde el free-streaming.

Si los tres dan S8 parecido -> S8 no pudo haber seleccionado la masa.
Si difieren mucho -> la objecion de Mike es correcta y hay que declararlo.
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
        'cs2_fld': 1.0, 'output': 'mPk', 'P_k_max_h/Mpc': 25.0, 'z_max_pk': 1.0}
KH = np.logspace(-3, np.log10(19.), 250)


def corre(m, T, pk=False):
    p = dict(BASE); p.update(
        {'N_ncdm': 2, 'm_ncdm': f'{SUM_MNU},{m}', 'T_ncdm': f'0.71611,{T}',
         'deg_ncdm': '1.0,1.0'})
    c = Class(); c.set(p); c.compute()
    r = dict(Onu=c.Onu if False else c.Omega_nu, s8=c.sigma8(), Om=c.Omega_m())
    if pk:
        r['pk'] = np.array([c.pk_lin(k * h, 0.0) * h**3 for k in KH])
    c.struct_cleanup(); c.empty()
    return r


# control frio (mismo Omega_m, sin free-streaming)
p = dict(BASE); p['Omega_cdm'] = OM_CDM + OM_PHI
p.update({'N_ncdm': 1, 'm_ncdm': f'{SUM_MNU}', 'T_ncdm': '0.71611'})
c = Class(); c.set(p); c.compute()
S8_COLD, S8C_OM = c.sigma8(), c.Omega_m()
PK_COLD = np.array([c.pk_lin(k * h, 0.0) * h**3 for k in KH])
c.struct_cleanup(); c.empty()
print(f'  CONTROL frio: Omega_m={S8C_OM:.6f}  sigma8={S8_COLD:.6f}  '
      f'S8={S8_COLD*np.sqrt(S8C_OM/0.3):.6f}\n')

CAND = [('SOLAR^2*OMEGA    ', 297.03), ('SOLAR^2*PYROS    ', 398.14),
        ('SOLAR^2*KRYSTOS_V', 594.28), ('[retirado] suma  ', 535.28),
        ('[retirado] PVM   ', 615.33)]
print(f'  {"multiplicador":18} {"mult":>8} {"m_phi(eV)":>10} {"T_ncdm":>9} '
      f'{"sigma8":>9} {"S8":>8} {"k50":>7}  {"n vs KiDS":>10}')
res = []
for nm, mult in CAND:
    m = SUM_MNU * mult
    T = brentq(lambda t: corre(m, t)['Onu'] - OM_NU_ACT - OM_PHI, 0.15, 0.60,
               xtol=1e-6)
    r = corre(m, T, pk=True)
    S8 = r['s8'] * np.sqrt(r['Om'] / 0.3)
    T2 = r['pk'] / PK_COLD
    k50 = KH[np.argmin(abs(T2 - 0.5))]
    n = (S8 - 0.759) / 0.024
    print(f'  {nm:18} {mult:8.2f} {m:10.3f} {T:9.6f} {r["s8"]:9.6f} '
          f'{S8:8.5f} {k50:7.3f}  {n:+10.2f}')
    res.append((nm, mult, m, S8, k50))

S8s = np.array([r[3] for r in res[:3]])
print(f'\n  RANGO de S8 entre los 3 candidatos de la gramatica: '
      f'{S8s.min():.5f} a {S8s.max():.5f}   (spread = {S8s.ptp():.5f})')
print(f'  incertidumbre de KiDS-1000                        : +-0.024')
print(f'  spread / sigma_KiDS = {S8s.ptp()/0.024:.3f}')
print(f'\n  => si spread/sigma << 1, S8 NO puede haber seleccionado la masa.')
