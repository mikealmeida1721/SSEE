#!/usr/bin/env python3
"""
Likelihood de Cobaya para KiDS-1000 xi_pm — MCMC completo, un sector, con los
20 y pico parametros de molestia MARGINALIZADOS (no minimizados).

Envuelve la fisica ya validada (kids_shear.py: control negativo contra la
cadena oficial, chi2 261.27 vs 260.32 oficial, 0.4%). Este modulo NO reescribe
esa fisica: solo la conecta a un muestreador de verdad.

Dos configuraciones, definidas por si el fondo es libre o fijo:
  SSEE : fondo FIJO por algebra (omega_b, omega_c, h, ns, w0, wa de ssee_core).
         Libres: logA, halo_A, A_IA, 5 dz, delta_c  -> 9 parametros.
  LCDM : fondo LIBRE (omch2, ombh2, h0, ns; w=-1, wa=0, mnu=0.06 fiducial).
         Libres: esos 4 + logA, halo_A, A_IA, 5 dz, delta_c -> 13 parametros.

Priors de nuisance = los oficiales KiDS-1000 (values.ini / priors.ini):
  halo_A         ~ U(2.0, 3.13)
  A_IA           ~ U(-6.0, 6.0)
  dz_i           ~ Gaussianas SOM_cov (medias del SOM KV450, tabla 3)
  delta_c        ~ Gaussiana(0, 2.3e-4)
Priors de fondo LCDM = los oficiales (rangos amplios, planos):
  omch2 in [0.051, 0.255]  ombh2 in [0.019, 0.026]  h0 in [0.64, 0.82]
  ns in [0.84, 1.1]

NO hay segundo sector, NO hay particula: retirados el 2026-07-30 (el 0.160 no
es una densidad, la resta que definia Omega_phiDM no significaba nada, y una
particula solo falsable no es una particula real).
"""
import sys
import time

import numpy as np

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src/p06_growth')
import kids_shear as K
import ssee_core as S

np.seterr(all='ignore')

D = K.load_data()
MASK = K.scale_mask(D)
CINV = np.linalg.inv(D['C'][np.ix_(MASK, MASK)])

DZ_MEAN = np.array([0.000, -0.002, -0.013, -0.011, 0.006])
SOM_COV = np.loadtxt('/mnt/datos/SSEE_data/kids1000/kcap/runs/3x2pt/'
                     'data_iterated_cov/cosmology/multinest_blindC_EE_nE_w/'
                     'data/KiDS/SOM_cov_multiplied.asc')
SOM_INV = np.linalg.inv(SOM_COV)
DELTA_C_SIG = 2.3e-4

# fondo SSEE, fijo, leido del nucleo algebraico
SSEE_BG = dict(ombh2=S.OMEGA_B_H2, omch2=S.OMEGA_C_H2, h0=S.H0_ALG / 100.0,
              ns=S.N_S, mnu=S.SUM_MNU_EV, w0=S.W0, wa=S.WA)

_CACHE = {}


def _camb_cached(bg_key, As, halo_A):
    key = (bg_key, round(As, 15), round(halo_A, 6))
    if key not in _CACHE:
        bg = dict(SSEE_BG) if bg_key == 'SSEE' else bg_key_to_dict(bg_key)
        r, p, kh, zpk, pk, gr = K.run_camb(
            omch2=bg['omch2'], ombh2=bg['ombh2'], h0=bg['h0'], ns=bg['ns'],
            As=As, mnu=bg['mnu'], w=bg['w0'], wa=bg['wa'], halo_A=halo_A)
        if len(_CACHE) > 300:
            _CACHE.clear()
        _CACHE[key] = (r, p, kh, zpk, pk, gr)
    return _CACHE[key]


def bg_key_to_dict(bg_key):
    # bg_key para LCDM: tupla (ombh2, omch2, h0, ns) redondeada -> dict completo
    ombh2, omch2, h0, ns = bg_key
    return dict(ombh2=ombh2, omch2=omch2, h0=h0, ns=ns, mnu=0.06, w0=-1.0, wa=0.0)


def loglike(bg_key, logA, halo_A, A_IA, dz1, dz2, dz3, dz4, dz5, delta_c):
    dz = np.array([dz1, dz2, dz3, dz4, dz5])
    As = np.exp(logA) * 1e-10
    try:
        r, p, kh, zpk, pk, gr = _camb_cached(bg_key, As, halo_A)
        ells, Cl, idx = K.cl_shear(D, r, p, kh, zpk, pk, gr, A_IA, dz)
        th = K.theory_vector(D, ells, Cl, idx, delta_c=delta_c)
    except Exception:
        return -1e10
    dv = (th - D['d'])[MASK]
    chi2 = float(dv @ CINV @ dv)
    rr = dz - DZ_MEAN
    prior_dz = float(rr @ SOM_INV @ rr)
    prior_dc = (delta_c / DELTA_C_SIG) ** 2
    return -0.5 * (chi2 + prior_dz + prior_dc)


def loglike_ssee(logA, halo_A, A_IA, dz1, dz2, dz3, dz4, dz5, delta_c):
    return loglike('SSEE', logA, halo_A, A_IA, dz1, dz2, dz3, dz4, dz5, delta_c)


def loglike_lcdm(ombh2, omch2, h0, ns, logA, halo_A, A_IA,
                 dz1, dz2, dz3, dz4, dz5, delta_c):
    bg_key = (round(ombh2, 6), round(omch2, 6), round(h0, 5), round(ns, 5))
    return loglike(bg_key, logA, halo_A, A_IA, dz1, dz2, dz3, dz4, dz5, delta_c)


# ---------------- configuracion Cobaya ----------------
NUISANCE_PARAMS = dict(
    halo_A=dict(prior=dict(min=2.0, max=3.13), ref=2.6, proposal=0.05,
               latex=r'A_\mathrm{baryon}'),
    A_IA=dict(prior=dict(min=-6.0, max=6.0), ref=1.0, proposal=0.3,
              latex=r'A_\mathrm{IA}'),
    dz1=dict(prior=dict(dist='norm', loc=0.000, scale=1.0), ref=0.0,
             proposal=0.01, latex=r'\delta z_1'),
    dz2=dict(prior=dict(dist='norm', loc=-0.002, scale=1.0), ref=-0.002,
             proposal=0.01, latex=r'\delta z_2'),
    dz3=dict(prior=dict(dist='norm', loc=-0.013, scale=1.0), ref=-0.013,
             proposal=0.01, latex=r'\delta z_3'),
    dz4=dict(prior=dict(dist='norm', loc=-0.011, scale=1.0), ref=-0.011,
             proposal=0.01, latex=r'\delta z_4'),
    dz5=dict(prior=dict(dist='norm', loc=0.006, scale=1.0), ref=0.006,
             proposal=0.01, latex=r'\delta z_5'),
    delta_c=dict(prior=dict(min=-1.2e-3, max=1.2e-3), ref=0.0,
                proposal=1e-4, latex=r'\delta_c'),
)
# NOTA CRITICA: el prior REAL (SOM_INV para dz, gaussiano para delta_c) se
# aplica UNA sola vez, DENTRO de loglike (prior_dz, prior_dc). Los priors
# declarados aqui para Cobaya son deliberadamente ANCHOS/PLANOS sobre el rango
# fisico (dz: gaussiana de escala 1.0, ~100x mas ancha que el shift real de
# ~0.01; delta_c: plano en el rango oficial de values.ini, +-1.2e-3, ~5x el
# sigma real de 2.3e-4) para que Cobaya sepa donde muestrear sin volver a
# contar la restriccion.
#
# Version anterior (BUG encontrado antes de lanzar la corrida real): delta_c
# se declaraba aqui con prior=norm(0, DELTA_C_SIG) -- la MISMA escala fisica
# que el termino ya sumado en loglike. Eso contaba el mismo prior gaussiano
# DOS VECES y estrechaba el posterior de delta_c de mas. Detectado en revision
# antes de lanzar, no en una cadena a medio correr.


def info_ssee(chains_dir):
    p = dict(logA=dict(prior=dict(min=1.5, max=4.5), ref=3.04, proposal=0.05,
                       latex='\\log(10^{10}A_s)'))
    p.update(NUISANCE_PARAMS)
    return dict(
        likelihood={'p06_growth.cobaya_kids.loglike_ssee': {
            'external': loglike_ssee, 'input_params': list(p.keys())}},
        params=p,
        sampler={'mcmc': {'Rminus1_stop': 0.03, 'max_tries': 10000}},
        output=chains_dir + '/ssee', force=True, resume=False)


def info_lcdm(chains_dir):
    p = dict(
        ombh2=dict(prior=dict(min=0.019, max=0.026), ref=0.02237,
                   proposal=0.0005, latex=r'\Omega_b h^2'),
        omch2=dict(prior=dict(min=0.051, max=0.255), ref=0.12,
                  proposal=0.005, latex=r'\Omega_c h^2'),
        h0=dict(prior=dict(min=0.64, max=0.82), ref=0.6736, proposal=0.02,
               latex='h'),
        ns=dict(prior=dict(min=0.84, max=1.1), ref=0.9649, proposal=0.02,
               latex='n_s'),
        logA=dict(prior=dict(min=1.5, max=4.5), ref=3.04, proposal=0.05,
                  latex='\\log(10^{10}A_s)'))
    p.update(NUISANCE_PARAMS)
    return dict(
        likelihood={'p06_growth.cobaya_kids.loglike_lcdm': {
            'external': loglike_lcdm, 'input_params': list(p.keys())}},
        params=p,
        sampler={'mcmc': {'Rminus1_stop': 0.03, 'max_tries': 10000}},
        output=chains_dir + '/lcdm', force=True, resume=False)


if __name__ == '__main__':
    from cobaya.run import run
    model_name = sys.argv[1] if len(sys.argv) > 1 else 'ssee'
    chains_dir = sys.argv[2] if len(sys.argv) > 2 else \
        '/mnt/datos/SSEE_data/chains_p6/kids'
    info = info_ssee(chains_dir) if model_name == 'ssee' else info_lcdm(chains_dir)
    t0 = time.time()
    updated_info, sampler = run(info)
    print(f'\nTERMINADO en {(time.time()-t0)/3600:.2f} h', flush=True)
