#!/usr/bin/env python3
"""CONTROL NEGATIVO: reproducir like = -130.157350 de la cadena oficial KiDS-1000 xi_pm.
Punto de maxima posterior, Blind C (maxpost_multinest_start_C.txt)."""
import numpy as np, time, sys
import kids_shear as K

# --- punto de maxima posterior de la cadena oficial ---
MP = dict(omch2=0.079140, ombh2=0.019120, h0=0.666430, ns=0.926990,
          halo_A=2.816460, A_IA=0.387120, delta_c=-0.000001,
          sigma8_target=0.887050, S8_ref=0.764230, Om_ref=0.222680)
DZ = np.array([0.002320, 0.008922, -0.018278, -0.014158, 0.006854])  # delta_z_out
LIKE_REF = -130.157350

D = K.load_data(); m = K.scale_mask(D)
print(f'puntos usados: {int(m.sum())} de {len(D["d"])}', flush=True)

# --- calibrar A_s para el sigma8 del punto de maxima posterior ---
As = 2.1e-9
for it in range(6):
    res, p, kh, zpk, pk, gr = K.run_camb(As=As, omch2=MP['omch2'], ombh2=MP['ombh2'],
                                         h0=MP['h0'], ns=MP['ns'], halo_A=MP['halo_A'])
    s8 = res.get_sigma8_0()
    print(f'  iter {it}: As={As:.6e}  sigma8={s8:.6f}', flush=True)
    if abs(s8 - MP['sigma8_target']) < 1e-5:
        break
    As *= (MP['sigma8_target'] / s8) ** 2

Om = (p.omch2 + p.ombh2 + p.omnuh2) / p.h**2
print(f'Om calculado = {Om:.6f}   (cadena: {MP["Om_ref"]:.6f})', flush=True)
print(f'S8 calculado = {s8*np.sqrt(Om/0.3):.6f}   (cadena: {MP["S8_ref"]:.6f})', flush=True)

t = time.time()
ells, Cl, idx = K.cl_shear(D, res, p, kh, zpk, pk, gr, MP['A_IA'], DZ)
print(f'C_ell listo en {time.time()-t:.1f}s', flush=True)
th = K.theory_vector(D, ells, Cl, idx, delta_c=MP['delta_c'])
c2 = K.chi2(D, th, m)

# like de CosmoSIS = -0.5*chi2 - 0.5*ln|2 pi C|
Csub = D['C'][np.ix_(m, m)]
sign, logdet = np.linalg.slogdet(2 * np.pi * Csub)
like_full = -0.5 * c2 - 0.5 * logdet
print()
print(f'chi2            = {c2:.4f}   (ndata={int(m.sum())})')
print(f'-0.5*chi2       = {-0.5*c2:.4f}')
print(f'-0.5*ln|2piC|   = {-0.5*logdet:.4f}')
print(f'like completo   = {like_full:.4f}')
print(f'REFERENCIA      = {LIKE_REF:.4f}')
print(f'   -> chi2 implicito por la referencia = {-2*(LIKE_REF + 0.5*logdet):.4f}')
np.save('control_theory.npy', th)
