#!/usr/bin/env python3
"""
Barra de error del lado KiDS bajo SSEE: perfil de chi2 sobre A_s.

Por cada logA se re-optimizan los dos nuisance que degeneran con la amplitud
(halo_A y A_IA); dz y delta_c quedan en el optimo de a1 (su efecto sobre la
amplitud global es de segundo orden). Delta chi2 = 1 da sigma(S8) del lado KiDS.

Homologo del sigma(S8) publicado por KiDS bajo LCDM (+-0.024), para poder
expresar la tension de SSEE en la receta de la literatura sin suponer barras.
"""
import numpy as np, sys, json, time
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-1000/-home-mike-Proyectos-SSEE/'
                   '238cf748-10c3-466e-9be3-e8e03009063e/scratchpad')
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import kids_shear as K
import ssee_core as S

OMB, OMC, HH, NS = S.OMEGA_B_H2, S.OMEGA_C_H2, S.H0_ALG / 100.0, S.N_S
MNU, W0, WA = S.SUM_MNU_EV, S.W0, S.WA

# optimo de a1 (un sector, A_s libre)
LOGA0 = 2.8670623345914974
HALO0, AIA0 = 2.5943533750053858, 0.5836089120795278
DZ0 = np.array([0.0005553766492810671, 3.21174623276652e-05,
                -0.021349310066011436, -0.013879874363113387,
                0.006600669165789613])
DC0 = -2.0122818051196982e-07

SOM_COV = np.loadtxt('/mnt/datos/SSEE_data/kids1000/kcap/runs/3x2pt/'
                     'data_iterated_cov/cosmology/multinest_blindC_EE_nE_w/'
                     'data/KiDS/SOM_cov_multiplied.asc')
DZ_MEAN = np.array([0.000, -0.002, -0.013, -0.011, 0.006])
SOM_INV = np.linalg.inv(SOM_COV)
HALO_RANGE, AIA_RANGE = (2.0, 3.13), (-6.0, 6.0)

D = K.load_data(); MASK = K.scale_mask(D)
CINV = np.linalg.inv(D['C'][np.ix_(MASK, MASK)])
_C = {}


def chi2_at(logA, halo_A, A_IA):
    if not (HALO_RANGE[0] <= halo_A <= HALO_RANGE[1]): return 1e10, 0.0
    if not (AIA_RANGE[0] <= A_IA <= AIA_RANGE[1]): return 1e10, 0.0
    As = 1e-10 * np.exp(logA)
    key = (round(As, 18), round(halo_A, 6))
    if key not in _C:
        if len(_C) > 300: _C.clear()
        _C[key] = K.run_camb(omch2=OMC, ombh2=OMB, h0=HH, ns=NS, As=As,
                             mnu=MNU, w=W0, wa=WA, halo_A=halo_A)
    res, p, kh, zpk, pk, gr = _C[key]
    try:
        ells, Cl, idx = K.cl_shear(D, res, p, kh, zpk, pk, gr, A_IA, DZ0)
        th = K.theory_vector(D, ells, Cl, idx, delta_c=DC0)
    except Exception:
        return 1e10, 0.0
    dv = (th - D['d'])[MASK]
    return float(dv @ CINV @ dv), float(res.get_sigma8_0())


def profile(logA):
    f = lambda x: chi2_at(logA, x[0], x[1])[0]
    r = minimize(f, np.array([HALO0, AIA0]), method='Nelder-Mead',
                 options=dict(maxiter=200, xatol=1e-3, fatol=1e-2))
    c2, s8 = chi2_at(logA, r.x[0], r.x[1])
    return c2, s8, r.x[0], r.x[1]


if __name__ == '__main__':
    OM = S.OMEGA_M_CMB
    grid = LOGA0 + np.array([-0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.09])
    out = []
    t0 = time.time()
    print(f'# perfil chi2(logA) alrededor de logA0={LOGA0:.6f}; Om={OM:.6f}',
          flush=True)
    for lg in grid:
        c2, s8, hA, aI = profile(float(lg))
        s8_ = s8 * np.sqrt(OM / 0.3)
        out.append(dict(logA=float(lg), chi2=c2, sigma8=s8, S8=s8_,
                        halo_A=hA, A_IA=aI))
        print(f'  logA={lg:.5f}  chi2={c2:9.4f}  sigma8={s8:.5f}  '
              f'S8={s8_:.5f}  halo_A={hA:.4f}  A_IA={aI:.4f}', flush=True)
    json.dump(out, open('s8_barra_kids.json', 'w'), indent=1)

    # parabola en (S8, chi2) -> sigma donde Delta chi2 = 1
    a = np.array([o['S8'] for o in out]); b = np.array([o['chi2'] for o in out])
    ok = b < 1e9
    co = np.polyfit(a[ok], b[ok], 2)
    S8_min = -co[1] / (2 * co[0])
    sig = 1.0 / np.sqrt(co[0])          # Delta chi2 = 1  ->  dS8 = 1/sqrt(A)
    print(f'\n  S8 minimo (perfil) = {S8_min:.5f}')
    print(f'  sigma(S8) lado KiDS bajo SSEE = {sig:.5f}')
    print(f'  [KiDS publica +-0.024 bajo LCDM]')
    S8_cmb, sig_cmb = 0.82639, 0.00564
    tot = np.sqrt(sig**2 + sig_cmb**2)
    print(f'\n  tension receta-literatura SSEE = '
          f'({S8_cmb:.5f}-{S8_min:.5f})/{tot:.5f} = '
          f'{(S8_cmb-S8_min)/tot:.2f} sigma')
    print(f'  minutos: {(time.time()-t0)/60:.1f}')
