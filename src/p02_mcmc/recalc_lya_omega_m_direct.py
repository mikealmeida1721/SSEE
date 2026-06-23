#!/usr/bin/env python3
"""
Recálculo Ly-α D_H/r_d @ z=2.33 bajo el reframe ω_m-directo (reemplaza el
'Scenario B — MIRA-mapped' del Apéndice B.2 de Paper 2, que usaba Ω_m=0.320).

Física del reframe: la EXPANSIÓN de fondo E(z) cuenta TODA la materia no
relativista. A z=2.33 la φ-DM (m=41 eV, no-rel desde z~4.5e5) cuenta como
materia completa → E(z) usa Ω_m,CMB=0.30889 (CDM 0.160 + φ-DM 0.149).
El free-streaming k_fs solo suprime el CLUSTERING a k>k_fs, NO el fondo.
r_d sale directo de ω_m (CAMB) = 147.17 Mpc (no del mapeo ×|w_0|).
"""
import numpy as np

C_KM_S = 299792.458
PHI = (1+np.sqrt(5))/2; PI = np.pi
w0, wa = -0.840, -0.670
H0 = 3*(PHI+PI)**2            # 67.962 (ancla algebraica ω_m-directo)
z = 2.33
DESI_obs, DESI_err = 8.52, 0.17   # DESI DR2 Ly-α D_H/r_d

def E_of_z(z, Om, Ode):
    fDE = (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*z/(1+z))
    return np.sqrt(Om*(1+z)**3 + Ode*fDE)

def DH_over_rd(Om, Ode, rd):
    E = E_of_z(z, Om, Ode)
    DH = C_KM_S/(H0*E)
    return E, DH, DH/rd

print("="*72)
print("  Ly-α D_H/r_d @ z=2.33 — reframe ω_m-directo vs MIRA (retirado)")
print(f"  H0={H0:.3f}  w0={w0}  wa={wa}  DESI obs={DESI_obs}±{DESI_err}")
print("="*72)
scenarios = [
    ("A  dinámico crudo (0.160, r_d=175.6)", 0.160, 0.840, 175.6),
    ("A' dinámico, r_d ω_m-direct (0.160, r_d=147.17)", 0.160, 0.840, 147.17),
    ("B_viejo MIRA (0.320, r_d=147.2)  [RETIRADO]", 0.320, 0.680, 147.2),
    ("B_nuevo ω_m-direct (0.30889, r_d=147.17)", 0.30889, 0.69111, 147.17),
]
for name, Om, Ode, rd in scenarios:
    E, DH, ratio = DH_over_rd(Om, Ode, rd)
    sig = (ratio-DESI_obs)/DESI_err
    print(f"  {name}")
    print(f"     E(2.33)={E:.3f}  D_H={DH:.1f} Mpc  D_H/r_d={ratio:.2f}  ({sig:+.2f}σ)")
print("="*72)
print("  Nota: el Apéndice B.2 dice 'E=2.206' para Scenario B — es ERRATA;")
print("  con Ω_m=0.320 el E real es ~3.50 (consistente con su propio D_H=1260).")
