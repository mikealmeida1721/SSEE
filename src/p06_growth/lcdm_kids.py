#!/usr/bin/env python3
"""
CONTROL SIMETRICO: LCDM ajustado al MISMO CMB que SSEE, evaluado contra los
MISMOS 225 puntos de xi_pm de KiDS-1000, con el MISMO codigo y los MISMOS priors.

Cosmologia de LCDM tomada de la cadena results/chains/lcdm_cmb.[1-4].txt
(likelihood identico al de ssee_cmb_k2: planck_2018_highl_plik.TTTEEE +
lowl.TT + lowl.EE + lensing.native; w=-1, wa=0, mnu=0.06 fijos).

Ninguno de los dos modelos ve la cizalla antes de predecir.

  (c1) LCDM-Planck, A_s FIJO al de su cadena CMB  -> 0 libres cosmologicos
       homologo exacto de b1_un_sector_As_CMB (SSEE)
  (c2) LCDM-Planck, A_s libre                     -> 1 libre cosmologico
       homologo exacto de a1_un_sector_As_libre (SSEE)
"""
import numpy as np, sys, time, json
from scipy.optimize import minimize
sys.path.insert(0, '/tmp/claude-1000/-home-mike-Proyectos-SSEE/'
                   '238cf748-10c3-466e-9be3-e8e03009063e/scratchpad')
import kids_shear as K

# ---------- cosmologia LCDM-Planck: media ponderada de las 4 cadenas ----------
OMB = 0.02238396
OMC = 0.11995244
HH = 0.6737330790
NS = 0.96530732
MNU = 0.06          # fijado en el yaml de la cadena
W0, WA = -1.0, 0.0  # fijados en el yaml de la cadena
LOGA_LCDM = 3.04479445
AS_LCDM = 1e-10 * np.exp(LOGA_LCDM)

# ---------- priors de KiDS, IDENTICOS a los de ssee_kids.py ----------
SOM_COV = np.loadtxt('/mnt/datos/SSEE_data/kids1000/kcap/runs/3x2pt/'
                     'data_iterated_cov/cosmology/multinest_blindC_EE_nE_w/'
                     'data/KiDS/SOM_cov_multiplied.asc')
DZ_MEAN = np.array([0.000, -0.002, -0.013, -0.011, 0.006])
SOM_INV = np.linalg.inv(SOM_COV)
DELTA_C_SIG = 2.3e-4
HALO_A_RANGE = (2.0, 3.13)
A_IA_RANGE = (-6.0, 6.0)

D = K.load_data()
MASK = K.scale_mask(D)
CINV = np.linalg.inv(D['C'][np.ix_(MASK, MASK)])

_CACHE = {}


def _camb_cached(As, halo_A):
    key = (round(As, 16), round(halo_A, 6))
    if key not in _CACHE:
        if len(_CACHE) > 400:
            _CACHE.clear()
        _CACHE[key] = K.run_camb(omch2=OMC, ombh2=OMB, h0=HH, ns=NS, As=As,
                                 mnu=MNU, w=W0, wa=WA, halo_A=halo_A)
    return _CACHE[key]


def model(As, halo_A, A_IA, dz, delta_c):
    res, p, kh, zpk, pk, gr = _camb_cached(As, halo_A)
    ells, Cl, idx = K.cl_shear(D, res, p, kh, zpk, pk, gr, A_IA, dz)
    th = K.theory_vector(D, ells, Cl, idx, delta_c=delta_c)
    return th, res.get_sigma8_0(), (p.omch2 + p.ombh2 + p.omnuh2) / p.h**2


def neglogpost(x, As_fixed):
    i = 0
    if As_fixed is None:
        As = np.exp(x[i]) * 1e-10; i += 1
    else:
        As = As_fixed
    halo_A, A_IA = x[i], x[i + 1]; i += 2
    dz = x[i:i + 5]; i += 5
    delta_c = x[i]
    if not (HALO_A_RANGE[0] <= halo_A <= HALO_A_RANGE[1]): return 1e10
    if not (A_IA_RANGE[0] <= A_IA <= A_IA_RANGE[1]): return 1e10
    try:
        th, _, _ = model(As, halo_A, A_IA, dz, delta_c)
    except Exception:
        return 1e10
    dv = (th - D['d'])[MASK]
    chi2 = float(dv @ CINV @ dv)
    r = dz - DZ_MEAN
    prior = float(r @ SOM_INV @ r) + (delta_c / DELTA_C_SIG) ** 2
    return 0.5 * (chi2 + prior)


def run(label, As_fixed):
    x0 = [] if As_fixed is not None else [LOGA_LCDM]
    x0 += [2.6, 1.0] + list(DZ_MEAN) + [0.0]
    t = time.time()
    r = minimize(neglogpost, np.array(x0), args=(As_fixed,),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-4, fatol=1e-3, disp=False))
    x = r.x; i = 0
    As = (np.exp(x[0]) * 1e-10) if As_fixed is None else As_fixed
    if As_fixed is None: i = 1
    halo_A, A_IA = x[i], x[i + 1]; dz = x[i + 2:i + 7]; delta_c = x[i + 7]
    th, s8, Om = model(As, halo_A, A_IA, dz, delta_c)
    dv = (th - D['d'])[MASK]
    chi2 = float(dv @ CINV @ dv)
    out = dict(label=label, chi2=chi2, ndata=int(MASK.sum()), sigma8=float(s8),
               S8=float(s8 * np.sqrt(Om / 0.3)), Om=float(Om), As=float(As),
               logA=float(np.log(As * 1e10)), halo_A=float(halo_A),
               A_IA=float(A_IA), dz=[float(v) for v in dz],
               delta_c=float(delta_c), nfree=len(x0), secs=time.time() - t)
    print(json.dumps(out), flush=True)
    return out


if __name__ == '__main__':
    print(f'# cosmologia LCDM-Planck fija: omb={OMB:.8f} omc={OMC:.8f} '
          f'h={HH:.8f} ns={NS:.8f} mnu={MNU} w0={W0} wa={WA}', flush=True)
    print(f'# A_s CMB LCDM = {AS_LCDM:.6e}  (logA={LOGA_LCDM})', flush=True)
    R = [run('c1_lcdm_planck_As_CMB', AS_LCDM)]
    json.dump(R, open('lcdm_kids_results.json', 'w'), indent=1)
    R.append(run('c2_lcdm_planck_As_libre', None))
    json.dump(R, open('lcdm_kids_results.json', 'w'), indent=1)
