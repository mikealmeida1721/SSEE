#!/usr/bin/env python3
"""
SSEE vs KiDS-1000 cizarra cosmica — evaluacion DIRECTA contra los datos crudos.

Cosmologia de fondo FIJA por algebra (phi, pi); solo se marginalizan los
parametros de molestia de la MEDIDA (astrofisica + instrumento), con los
priors publicados por KiDS.

Cuatro corridas:
  (a1) un sector,  A_s libre
  (a2) dos sectores, A_s libre
  (b1) un sector,  A_s FIJO al de la cadena CMB de SSEE  -> 0 libres cosmologicos
  (b2) dos sectores, A_s FIJO                             -> 0 libres cosmologicos

Pipeline validado: control negativo contra la cadena oficial KiDS-1000 xi_pm
reproduce like = -130.157 con desviacion 0.476 (chi2 261.27 vs 260.32, 0.4%).
"""
import numpy as np, sys, time, json
from scipy.optimize import minimize
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import kids_shear as K
import ssee_core as S

# ---------- cosmologia SSEE, todo fijo, leido de ssee_core ----------
OMB = S.OMEGA_B_H2          # 0.02241775678  = (pi-phi)/(3 Omega^2)
OMC = S.OMEGA_C_H2          # 0.11951440843  = KAL0 * omega_b * n_s
HH = S.H0_ALG / 100.0       # 0.6796213732   = 3(phi+pi)^2
NS = S.N_S                  # 0.9655581463   = 1 - phi^-7
MNU = S.SUM_MNU_EV          # 0.06849 eV
W0 = S.W0                   # -0.839949771
WA = S.WA                   # -0.669974886

# ---------- A_s de la cadena CMB de SSEE (k=2, 4 cadenas, 111099 muestras) ----------
LOGA_SSEE = 3.040704
AS_SSEE = 1e-10 * np.exp(LOGA_SSEE)

# ---------- segundo sector (forward, sin libres) ----------
M_PHI = 40.70               # eV   = Sigma m_nu * (SOLAR^2 * KRYSTOS_V)
OMEGA_C_DYN = 0.160050228655221
OMEGA_PHI = None            # se calcula abajo desde Omega_m total
ALPHA_WDM = 1.117           # Mpc/h  (calibracion CLASS del proyecto @ m_phi=40.70)
NU_VIEL = 1.12

# ---------- priors de KiDS sobre los desplazamientos de redshift ----------
SOM_COV = np.loadtxt('/mnt/datos/SSEE_data/kids1000/kcap/runs/3x2pt/'
                     'data_iterated_cov/cosmology/multinest_blindC_EE_nE_w/'
                     'data/KiDS/SOM_cov_multiplied.asc')
DZ_MEAN = np.array([0.000, -0.002, -0.013, -0.011, 0.006])   # priors.ini
SOM_INV = np.linalg.inv(SOM_COV)
DELTA_C_SIG = 2.3e-4
HALO_A_RANGE = (2.0, 3.13)
A_IA_RANGE = (-6.0, 6.0)

D = K.load_data()
MASK = K.scale_mask(D)
CINV = np.linalg.inv(D['C'][np.ix_(MASK, MASK)])


def suppression(kh):
    """Fraccion de amplitud que sobrevive al free-streaming, dos sectores.
       T_phi de Viel+2005; la materia total NO cambia, solo su reparto por escala."""
    Om_tot = (OMB + OMC + MNU / 93.14) / HH**2
    Om_phi = Om_tot - OMEGA_C_DYN
    T = (1.0 + (ALPHA_WDM * kh) ** (2 * NU_VIEL)) ** (-5.0 / NU_VIEL)
    return ((OMEGA_C_DYN + Om_phi * T) / Om_tot) ** 2


_CACHE = {}


def _camb_cached(As, halo_A, two_sector):
    """CAMB solo depende de A_s y del feedback bariónico; los otros 7 nuisance no."""
    key = (round(As, 16), round(halo_A, 6), bool(two_sector))
    if key not in _CACHE:
        res, p, kh, zpk, pk, gr = K.run_camb(omch2=OMC, ombh2=OMB, h0=HH, ns=NS,
                                             As=As, mnu=MNU, w=W0, wa=WA,
                                             halo_A=halo_A)
        if two_sector:
            pk = pk * suppression(kh)[None, :]
        if len(_CACHE) > 400:
            _CACHE.clear()
        _CACHE[key] = (res, p, kh, zpk, pk, gr)
    return _CACHE[key]


def model(As, halo_A, A_IA, dz, delta_c, two_sector):
    res, p, kh, zpk, pk, gr = _camb_cached(As, halo_A, two_sector)
    ells, Cl, idx = K.cl_shear(D, res, p, kh, zpk, pk, gr, A_IA, dz)
    th = K.theory_vector(D, ells, Cl, idx, delta_c=delta_c)
    return th, res.get_sigma8_0(), (p.omch2 + p.ombh2 + p.omnuh2) / p.h**2


def neglogpost(x, two_sector, As_fixed):
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
        th, _, _ = model(As, halo_A, A_IA, dz, delta_c, two_sector)
    except Exception:
        return 1e10
    dv = (th - D['d'])[MASK]
    chi2 = float(dv @ CINV @ dv)
    r = dz - DZ_MEAN
    prior = float(r @ SOM_INV @ r) + (delta_c / DELTA_C_SIG) ** 2
    return 0.5 * (chi2 + prior)


def run(label, two_sector, As_fixed):
    x0 = [] if As_fixed is not None else [3.0407]
    x0 += [2.6, 1.0] + list(DZ_MEAN) + [0.0]
    t = time.time()
    r = minimize(neglogpost, np.array(x0), args=(two_sector, As_fixed),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-4, fatol=1e-3, disp=False))
    x = r.x; i = 0
    As = (np.exp(x[0]) * 1e-10) if As_fixed is None else As_fixed
    if As_fixed is None: i = 1
    halo_A, A_IA = x[i], x[i + 1]; dz = x[i + 2:i + 7]; delta_c = x[i + 7]
    th, s8, Om = model(As, halo_A, A_IA, dz, delta_c, two_sector)
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
    print(f'# cosmologia SSEE fija: omb={OMB:.8f} omc={OMC:.8f} h={HH:.8f} '
          f'ns={NS:.8f} mnu={MNU} w0={W0:.6f} wa={WA:.6f}', flush=True)
    print(f'# A_s CMB SSEE = {AS_SSEE:.6e}  (logA={LOGA_SSEE})', flush=True)
    R = []
    # primero las de A_s FIJO: cache activa (rapidas) y son la prueba dura
    R.append(run('b1_un_sector_As_CMB',     False, AS_SSEE))
    R.append(run('b2_dos_sectores_As_CMB',   True, AS_SSEE))
    json.dump(R, open('ssee_kids_results.json', 'w'), indent=1)
    R.append(run('a1_un_sector_As_libre',   False, None))
    R.append(run('a2_dos_sectores_As_libre', True, None))
    json.dump(R, open('ssee_kids_results.json', 'w'), indent=1)
