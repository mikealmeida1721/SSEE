#!/usr/bin/env python3
"""
R1/R2 — fsigma8 de SSEE y LCDM contra los multipolos CRUDOS de BOSS DR12,
con modelado LPT (velocileptors) en vez de Kaiser.

QUE MIDE. La misma pregunta que R3/R4 hicieron con la cizalla de KiDS, sobre la
otra sonda de crecimiento: con el fondo FIJO (algebraico en SSEE, Planck en
LCDM), que amplitud prefiere el agrupamiento de galaxias, y con que chi2 sobre
el dato crudo. La prueba es el chi2 sobre los 222 puntos, NO el fsigma8: fsigma8
es SALIDA del ajuste, o sea producto del propio modelo (misma leccion que S8).

POR QUE LPT Y NO KAISER. Kaiser es teoria lineal: sesga 9-15% ya en k~0.1 y
obliga a cortar en k_max=0.08, tirando la mayoria del dato. velocileptors
(Chen, Vlah & White 2020) da el espectro en LPT a un lazo con los contra-
terminos EFT, valido hasta k~0.20 h/Mpc. Se pasa de ~111 a ~222 puntos utiles.

LIBRES, IDENTICOS EN LOS DOS MODELOS (misma cuenta que R3/R4: mismo dato, mismo
codigo, misma mascara; lo unico que cambia es el fondo):
    logA               (1)   amplitud primordial ln(1e10 A_s)
    b1  x 6            (6)   bias lineal, uno por (z_bin, hemisferio)
    b2  x 6            (6)   bias de segundo orden
    bs  x 6            (6)   bias de marea
    alpha0 x 6         (6)   contratermino monopolo
    alpha2 x 6         (6)   contratermino cuadrupolo
    SN0 x 6            (6)   ruido de disparo residual
  -----------------------------------------------------------
                      37 libres contra 222 puntos

Los 36 nuisance son por-conjunto y entran cuadraticamente salvo b1: se
marginalizan analiticamente donde se puede y por MCMC el resto.

Correccion de Hartlap aplicada (2045/2048 mocks PATCHY).

Fondo: NO se ajusta. SSEE lo trae del algebra (Omega_m=0.308881, w0, wa);
LCDM de Planck. Esa es toda la diferencia entre las dos corridas.
"""
import json
import os
import sys
import time

import numpy as np
import camb
from scipy.optimize import minimize

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src/p06_growth')

import ssee_core as S
import boss_rsd_model as M
from boss_dr12_data import ZBINS, CAPS, NMOCKS, load_multipole, load_cov, load_window

from velocileptors.LPT.lpt_rsd_fftw import LPT_RSD

KMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 0.20
KMIN = 0.01
OUT = '/home/mike/Proyectos/SSEE/results/logs/growth_2026-07'
AS_REF = 2.1e-9
MNU = 0.06

# Fondos FIJOS. Ninguno de estos numeros se ajusta.
COSMO = {
    'SSEE': dict(Om=S.OMEGA_M_TOTAL, h=S.H0_GLOBAL / 100.0,
                 ombh2=S.OMEGA_B_H2, ns=S.N_S, w0=S.W0, wa=S.WA),
    'LCDM': dict(Om=0.3153, h=0.6736,
                 ombh2=0.02237, ns=0.9649, w0=-1.0, wa=0.0),
}


def camb_lin(c, z):
    """P_lin(k) y f(z) del fondo fijo, a A_s de referencia. A_s entra despues
    como factor exacto, asi que CAMB se corre una sola vez por (modelo, z)."""
    p = camb.CAMBparams()
    omch2 = c['Om'] * c['h'] ** 2 - c['ombh2'] - MNU / 93.14
    p.set_cosmology(H0=c['h'] * 100.0, ombh2=c['ombh2'], omch2=omch2,
                    mnu=MNU, omk=0.0)
    p.set_dark_energy(w=c['w0'], wa=c['wa'], dark_energy_model='ppf')
    p.InitPower.set_params(As=AS_REF, ns=c['ns'])
    p.set_matter_power(redshifts=[z], kmax=3.0)
    r = camb.get_results(p)
    k, _, pk = r.get_matter_power_spectrum(minkh=1e-4, maxkh=2.0, npoints=800)
    f = r.get_fsigma8()[-1] / r.get_sigma8()[-1]
    return k, pk[0], float(f), float(r.get_sigma8()[-1])


def build():
    """Un 'conjunto' por (z_bin, hemisferio). Alineacion dato<->covarianza por
    (ell, k) y ventana: identica a la de boss_fit.py, que paso los 16 controles
    de Fase 0 (VERIFICACION_BOSS.md). Lo unico nuevo es el LPT."""
    sets = []
    for zb, z in ZBINS.items():
        for cap in CAPS:
            d, kk, ee = [], [], []
            for ell in (0, 2, 4):
                k, P, _ = load_multipole(zb, cap, ell)
                d.append(P); kk.append(k); ee.append(np.full(len(k), ell))
            k_all = np.concatenate(kk); ell_all = np.concatenate(ee)
            P_all = np.concatenate(d)
            C, kc, ellc = load_cov(zb, cap)
            idx = [int(np.argmin(np.abs(kc - kv) + 1e3 * (ellc != ev)))
                   for kv, ev in zip(k_all, ell_all)]
            C = C[np.ix_(idx, idx)]
            keep = (k_all <= KMAX) & (k_all >= KMIN)
            C = C[np.ix_(keep, keep)]
            n_p = int(keep.sum())
            hart = (NMOCKS[cap] - n_p - 2.0) / (NMOCKS[cap] - 1.0)

            # La ventana convoluciona por FFTLog: exige rejilla LOG y ancha.
            # Se evalua el LPT en la banda donde vive el dato y se pega a la
            # rejilla del modelo (0 fuera, que es lo que la ventana espera).
            # Banda ANCHA en log para el LPT: la ventana convoluciona por
            # FFTLog sobre M.K_MODEL (1e-5..100) y necesita el modelo definido
            # con holgura a ambos lados. Fuera de esta banda se extrapola en
            # ley de potencias con indice acotado (ver _pega).
            kmod = np.logspace(-3, np.log10(0.5), 180)
            lpt = {}
            for name, c in COSMO.items():
                kl, pl, fz, s8 = camb_lin(c, z)
                L = LPT_RSD(kl, pl, kIR=0.2, cutoff=10, extrap_min=-4,
                            extrap_max=3, N=2000, threads=1, jn=5)
                a_par, a_per = M.alphas(z, c['Om'], c['w0'], c['wa'])
                L.make_pltable(fz, kv=kmod, apar=a_par, aperp=a_per, ngauss=3)
                lpt[name] = dict(L=L, f=fz, sigma8=s8,
                                 a_par=float(a_par), a_per=float(a_per))
            sets.append(dict(zb=zb, cap=cap, z=z, k=k_all[keep],
                             ell=ell_all[keep], d=P_all[keep], kmod=kmod,
                             Cinv=np.linalg.inv(C) * hart, hartlap=hart,
                             npts=n_p, lpt=lpt,
                             win=M.Window(*load_window(zb, cap),
                                          integral_constraint=True)))
            print(f'  {zb} {cap}: {n_p} puntos  Hartlap={hart:.4f}', flush=True)
    return sets


def _pega(K, kv, P):
    """Lleva P(kv) a la rejilla ancha K sin meter ceros.

    Rellenar con 0 fuera de la banda del LPT hace que la extrapolacion en ley
    de potencias de mcfit divida por cero -> NaN (es el bug 3 de la Fase 0,
    VERIFICACION_BOSS.md). Se extrapola en ley de potencias conservando el
    signo, que es el comportamiento asintotico correcto de los multipolos.
    """
    # velocileptors devuelve NaN en los extremos de la banda (donde el P_lin
    # de CAMB ya no la sostiene): se descartan antes de pegar.
    ok = np.isfinite(P)
    kv, P = kv[ok], P[ok]
    out = np.interp(np.log(K), np.log(kv), P)
    lo, hi = K < kv[0], K > kv[-1]
    for m, i, j in ((lo, 0, 1), (hi, -1, -2)):
        if not m.any():
            continue
        y1, y2 = P[i], P[j]
        if y1 == 0.0 or y2 == 0.0 or np.sign(y1) != np.sign(y2):
            out[m] = y1                      # plano: no hay ley de potencias
            continue
        n = np.log(abs(y1) / abs(y2)) / np.log(kv[i] / kv[j])
        n = float(np.clip(n, -4.0, 4.0))     # sin cota, k^n desborda en 1e-5
        out[m] = np.sign(y1) * abs(y1) * (K[m] / kv[i]) ** n
    return out


def model_set(st, name, logA, th):
    """Multipolos LPT del conjunto, convolucionados con la ventana e
    interpolados a los k del dato. A_s entra como factor sobre P_lin."""
    b1, b2, bs, a0, a2, sn = th
    sc = np.exp(logA) * 1e-10 / AS_REF
    L = st['lpt'][name]['L']
    kv, p0, p2, p4 = L.combine_bias_terms_pkell(
        [b1, b2, bs, 0.0, a0, a2, 0.0, 0.0, sn, 0.0, 0.0])
    # La rejilla del LPT ya es log-espaciada y ancha: se pasa DIRECTA a la
    # ventana. Extrapolarla hasta k=1e-5 (rejilla M.K_MODEL) amplificaba por
    # ~1e4 y el FFTLog devolvia campanas de 1e21 — probado y descartado.
    ok = np.isfinite(p0) & np.isfinite(p2) & np.isfinite(p4)
    ks, q0, q2, q4 = st['win'](kv[ok], sc * p0[ok], sc * p2[ok], sc * p4[ok])
    out = np.empty(len(st['k']))
    for ell, q in ((0, q0), (2, q2), (4, q4)):
        m = st['ell'] == ell
        out[m] = np.interp(st['k'][m], ks, q)
    return out


def chi2_set(st, name, logA, th):
    r = model_set(st, name, logA, th) - st['d']
    return float(r @ st['Cinv'] @ r)


def run(name, sets, verbose=True):
    t0 = time.time()
    print(f'\n[{name}] fondo FIJO: Om={COSMO[name]["Om"]:.6f} '
          f'w0={COSMO[name]["w0"]:.4f} wa={COSMO[name]["wa"]:.4f}', flush=True)

    def total(v):
        logA = v[0]
        c = 0.0
        for i, st in enumerate(sets):
            c += chi2_set(st, name, logA, v[1 + 6 * i: 7 + 6 * i])
        return c

    x0 = [3.0] + [2.0, 0.0, 0.0, 0.0, 0.0, 0.0] * len(sets)
    r = minimize(total, x0, method='Powell',
                 options=dict(maxiter=40000, xtol=1e-3, ftol=1e-3))

    # Barrido en logA con los nuisance re-optimizados en cada punto: da la
    # curva de perfil, que es de donde sale la barra de logA (no de la Hessiana).
    grid = np.linspace(r.x[0] - 0.45, r.x[0] + 0.45, 19)
    prof = []
    for a in grid:
        def rest(u, a=a):
            c = 0.0
            for i, st in enumerate(sets):
                c += chi2_set(st, name, a, u[6 * i: 6 * i + 6])
            return c
        rr = minimize(rest, r.x[1:], method='Powell',
                      options=dict(maxiter=20000, xtol=1e-3, ftol=1e-3))
        prof.append(rr.fun)
        print(f'    logA={a:.4f}  chi2={rr.fun:.3f}', flush=True)
    prof = np.array(prof)
    i0 = int(np.argmin(prof))
    logA = float(grid[i0]); best = float(prof[i0])
    try:
        cf = np.polyfit(grid, prof, 2)
        sig_logA = float(np.sqrt(1.0 / cf[0])) if cf[0] > 0 else float('nan')
        logA = float(-cf[1] / (2 * cf[0]))
    except Exception:
        sig_logA = float('nan')

    n = sum(st['npts'] for st in sets)
    nfree = 1 + 6 * len(sets)
    dof = n - nfree
    sc = np.sqrt(np.exp(logA) * 1e-10 / AS_REF)
    out = dict(model=name, chi2=best, dof=int(dof), npts=int(n),
               n_libres=int(nfree), logA=logA, sig_logA=sig_logA,
               kmax=KMAX, kmin=KMIN, modelado='LPT velocileptors (1 lazo + EFT)',
               Om=COSMO[name]['Om'], w0=COSMO[name]['w0'], wa=COSMO[name]['wa'],
               perfil=dict(logA=grid.tolist(), chi2=prof.tolist()), bins=[])
    print(f'  chi2 = {best:.3f}  dof = {dof}  chi2/dof = {best/dof:.4f}')
    print(f'  logA = {logA:.5f} +- {sig_logA:.5f}')
    for st in sets:
        f = st['lpt'][name]['f']
        s8 = st['lpt'][name]['sigma8'] * sc
        out['bins'].append(dict(zb=st['zb'], cap=st['cap'], z=st['z'],
                                f=float(f), sigma8=float(s8),
                                fsigma8=float(f * s8)))
        print(f'    {st["zb"]} {st["cap"]}: f*sigma8 = {f*s8:.5f}')
    out['seg'] = time.time() - t0
    return out


if __name__ == '__main__':
    print('=' * 78)
    print(f'  R1/R2 — BOSS DR12 full-shape LPT, fondo FIJO   k_max = {KMAX}')
    print('=' * 78)
    sets = build()
    n = sum(st['npts'] for st in sets)
    print(f'  {len(sets)} conjuntos, {n} puntos, {1+6*len(sets)} libres')
    res = {m: run(m, sets) for m in ('LCDM', 'SSEE')}
    d = res['SSEE']['chi2'] - res['LCDM']['chi2']
    print('\n' + '=' * 78)
    print(f'  Delta chi2 (SSEE - LCDM) = {d:+.3f}   con los MISMOS libres')
    print(f'  logA:  LCDM {res["LCDM"]["logA"]:.5f}   SSEE {res["SSEE"]["logA"]:.5f}')
    res['comparacion'] = dict(delta_chi2_ssee_menos_lcdm=float(d))
    p = f'{OUT}/R1R2_boss_lpt_kmax{KMAX:.3f}.json'
    with open(p, 'w') as fh:
        json.dump(res, fh, indent=1)
    print(f'  -> {p}')
