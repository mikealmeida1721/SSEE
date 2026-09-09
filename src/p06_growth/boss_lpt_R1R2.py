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

FUENTE: results/logs/growth_2026-07/R1R2_boss_lpt_kmax0.200.json
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
# FIX 2026-09-08 — SSEE LLEVABA LA MASA DE NEUTRINO DE LCDM. Aqui habia un
# `MNU = 0.06` suelto que se usaba para LOS DOS modelos. 0.06 eV es el fiducial
# de Planck; la de SSEE es Sum m_nu = 0.06849 eV, derivada de la clausura
# nu = 93.14 y presente en omega_m = 0.142668. Es el mismo patron del `tau`
# prestado que se corrigio hoy: un modelo evaluado con el ingrediente del otro.
# MEDIDO antes de arreglarlo, mismo fondo, CAMB, z=0.51:
#     sigma8(0)      0.817411 -> 0.815158   -0.276%
#     fsigma8(0.51)  0.474540 -> 0.473366   -0.247%  =  0.058 sigma de la barra
# Pequeno, pero es un sesgo con signo, no ruido. Ahora cada modelo lleva la
# suya y la de SSEE se LEE del nucleo, no se re-teclea (R66).
MNU = {'SSEE': S.SUM_MNU_EV, 'LCDM': 0.06}

# Fondos FIJOS. Ninguno de estos numeros se ajusta.
COSMO = {
    'SSEE': dict(Om=S.OMEGA_M_TOTAL, h=S.H0_GLOBAL / 100.0,
                 ombh2=S.OMEGA_B_H2, ns=S.N_S, w0=S.W0, wa=S.WA,
                 mnu=MNU['SSEE']),
    'LCDM': dict(Om=0.3153, h=0.6736,
                 ombh2=0.02237, ns=0.9649, w0=-1.0, wa=0.0,
                 mnu=MNU['LCDM']),
}


def camb_lin(c, z):
    """P_lin(k) y f(z) del fondo fijo, a A_s de referencia. A_s entra despues
    como factor exacto, asi que CAMB se corre una sola vez por (modelo, z)."""
    p = camb.CAMBparams()
    omch2 = c['Om'] * c['h'] ** 2 - c['ombh2'] - c['mnu'] / 93.14
    p.set_cosmology(H0=c['h'] * 100.0, ombh2=c['ombh2'], omch2=omch2,
                    mnu=c['mnu'], omk=0.0)
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
    # Powell no lleva cotas y se sale del dominio en que el LPT esta definido:
    # con contraterminos absurdos el espectro sale todo no-finito, la mascara
    # queda vacia y el FFTLog de la ventana revienta (Nin=0). No es un fallo
    # fisico sino del optimizador: se devuelve None y chi2_set lo traduce a una
    # penalizacion grande, que es como se le pone pared a un minimizador sin
    # cotas. Umbral 10: por debajo de eso la interpolacion tampoco tendria
    # soporte. Corrida del 2026-09-07: sin esto, aborta a los ~3 min.
    if int(ok.sum()) < 10:
        return None
    ks, q0, q2, q4 = st['win'](kv[ok], sc * p0[ok], sc * p2[ok], sc * p4[ok])
    out = np.empty(len(st['k']))
    for ell, q in ((0, q0), (2, q2), (4, q4)):
        m = st['ell'] == ell
        out[m] = np.interp(st['k'][m], ks, q)
    return out


CHI2_PARED = 1.0e12          # penalizacion fuera del dominio del LPT

# ---------------------------------------------------------------------------
# Marginalizacion analitica de los parametros LINEALES (a0, a2, sn)
# ---------------------------------------------------------------------------
# POR QUE. De los 37 libres, 18 (a0, a2, sn por cada uno de los 6 conjuntos)
# entran LINEALES en el modelo, y son justo los que impedian converger. R-1
# por familia tras 45 min de cadena LCDM, 4 cadenas, mitad de burn-in:
#       a0 1.955   a2 1.502   sn 0.919  |  b1 0.524  b2 0.579  bs 0.618
#       logA 0.218
# Los tres peores son los tres lineales. Muestrearlos no aporta: su integral
# tiene forma cerrada.
#
# CONDICION PREVIA. La linealidad solo vale desde que la ventana lleva relleno
# de indice fijo (boss_rsd_model.Window, NOTA CRITICA). Con el relleno viejo el
# error de superposicion era de 0.2 a 1.2; ahora es 1e-12.
#
# LOS PRIORS SE MUDAN AQUI. En Cobaya estos 18 ya no existen, asi que su prior
# gaussiano N(0, sigma) entra en la matriz LAMBDA de abajo. Sigue habiendo UNA
# sola puerta --- ahora es esta --- que es la leccion del doble conteo de
# delta_c en cobaya_kids.py.
SIGMA_LIN = np.array([50.0, 50.0, 5.0e3])   # a0, a2, sn  (= priors de Cobaya)
_LAMBDA = np.diag(1.0 / SIGMA_LIN ** 2)
_TPL = {}


def templates(st, name):
    """Los tres modelos-plantilla d(modelo)/d(a0,a2,sn), a logA de referencia.

    No dependen de b1/b2/bs (verificado: cambian 1e-12 al pasar b1 de 1 a 2),
    asi que se calculan UNA vez por conjunto y se escalan por exp(logA).

    PERO SI dependen del FONDO. La clave llevaba solo (zb, cap, name), que es
    suficiente mientras el fondo no cambie —el caso de la cadena, donde solo
    varia logA—, y silenciosamente incorrecta en cuanto se barre un parametro
    del fondo: al reconstruir los conjuntos con otro w_c se reutilizaban las
    plantillas del PRIMER w_c calculado. Medido (2026-09-08): los mismos tres
    puntos daban 68.420/80.434/72.461 en un script y 76.828/78.760/81.600 en
    otro, segun donde hubiera arrancado cada barrido.
    sigma8 a A_s de referencia es firma del fondo (depende de w_c, h, w_b,
    n_s, w0, wa), asi que entra en la clave.
    """
    cl = (st['zb'], st['cap'], name,
          round(st['lpt'][name]['sigma8'], 12))
    if cl not in _TPL:
        z = np.zeros(6)
        m0 = model_set(st, name, 0.0, z)
        T = np.empty((len(m0), 3))
        for j in range(3):
            u = np.zeros(6)
            u[3 + j] = 1.0
            T[:, j] = model_set(st, name, 0.0, u) - m0
        _TPL[cl] = T
    return _TPL[cl]


def chi2_marg_set(st, name, logA, th3):
    """chi2 efectivo del conjunto con (a0,a2,sn) integrados analiticamente.

    th3 = (b1, b2, bs). Devuelve  chi2(lambda_hat) + ln det F, que es
    -2 ln de la integral gaussiana salvo una constante que no depende de
    ningun parametro (y por tanto no afecta a ninguna comparacion).
    """
    m0 = model_set(st, name, logA, np.array([th3[0], th3[1], th3[2],
                                             0.0, 0.0, 0.0]))
    if m0 is None or not np.all(np.isfinite(m0)):
        return CHI2_PARED, None
    # templates() ya se evaluo a logA=0, o sea con sc(0); el factor que
    # falta es solo la RAZON sc(logA)/sc(0) = exp(logA). Meter sc(logA)
    # entero metia un 1e-10/AS_REF = 0.0476 de mas y disparaba lambda_hat
    # (sn salia -5.5e4 con prior de escala 5e3).
    T = templates(st, name) * np.exp(logA)
    Ci = st['Cinv']
    r0 = m0 - st['d']
    CT = Ci @ T
    F = T.T @ CT + _LAMBDA
    v = -(r0 @ CT)                      # = T^T Cinv (d - m0)
    try:
        lam = np.linalg.solve(F, v)
    except np.linalg.LinAlgError:
        return CHI2_PARED, None
    c = float(r0 @ Ci @ r0 - v @ lam)
    sgn, ld = np.linalg.slogdet(F)
    if sgn <= 0 or not np.isfinite(c):
        return CHI2_PARED, None
    return c + float(ld), lam



def chi2_set(st, name, logA, th):
    m = model_set(st, name, logA, th)
    if m is None or not np.all(np.isfinite(m)):
        return CHI2_PARED
    r = m - st['d']
    c = float(r @ st['Cinv'] @ r)
    return c if np.isfinite(c) else CHI2_PARED


def run(name, sets, verbose=True):
    t0 = time.time()
    print(f'\n[{name}] fondo FIJO: Om={COSMO[name]["Om"]:.6f} '
          f'w0={COSMO[name]["w0"]:.4f} wa={COSMO[name]["wa"]:.4f}', flush=True)

    # POR QUE NO SE MINIMIZAN LOS 37 A LA VEZ (corregido 2026-09-07).
    # La primera version lanzaba Powell sobre los 37 libres sin cotas y se
    # quedaba en un minimo falso: logA=1.56 con chi2=998 sobre 222 puntos,
    # cuando logA debe rondar 3.0 y el chi2 unos 230. La culpa es de la
    # degeneracion A_s-b1 (a nivel lineal solo b1^2*A_s esta restringido):
    # el optimizador baja la amplitud y sube el bias sin que el chi2 lo frene,
    # y de camino se sale al dominio donde el LPT no esta definido.
    #
    # A logA FIJO el problema FACTORIZA de forma exacta: cada conjunto
    # (z_bin, hemisferio) tiene sus propios 6 nuisance y ninguno comparte
    # parametro con otro. O sea seis problemas de 6 dimensiones, no uno de 36.
    # Eso es lo que se hace aqui, con cotas y con arranque en caliente desde el
    # punto anterior de la rejilla. El perfil en logA es entonces la curva
    # correcta, y de ella sale la barra (no de la Hessiana).
    # ESCALADO. El segundo fallo, y el que de verdad rompia el ajuste: los seis
    # libres viven en escalas que van de b1~1 a SN0~1e4. Un minimizador de
    # busqueda directa da pasos del mismo tamano en todas las direcciones, asi
    # que o no mueve SN0 o destroza b1. Medido el 2026-09-07 sobre z1 NGC:
    #   sin escalar   chi2 = 118.8 / 4564.2 / 10416.8  en logA = 2.9 / 3.0 / 3.1
    #   escalado      chi2 =  40.6 /   40.0 /    39.9   (37 puntos, 31 dof)
    # Cada libre se lleva a O(1) dividiendo por su escala tipica.
    ESC = np.array([1.0, 1.0, 1.0, 50.0, 50.0, 1000.0])
    COTAS = [(0.5, 5.0),      # b1  LAGRANGIANO (b_Euler = b1 + 1)
             (-10.0, 10.0),   # b2  segundo orden
             (-10.0, 10.0),   # bs  marea
             (-4.0, 4.0),     # a0  contratermino monopolo   (x50)
             (-4.0, 4.0),     # a2  contratermino cuadrupolo (x50)
             (-10.0, 10.0)]   # SN0 ruido de disparo         (x1000)
    P0 = np.array([2.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    N_REINICIOS = 4           # Nelder-Mead reiniciado desde su propia salida

    def _baja(st, a, u0):
        """Nelder-Mead reiniciado desde su propia salida, en unidades escaladas."""
        u = np.asarray(u0, float); best = np.inf
        for _ in range(N_REINICIOS):
            rr = minimize(lambda v: chi2_set(st, name, a, v * ESC), u,
                          method='Nelder-Mead', bounds=COTAS,
                          options=dict(maxiter=4000, xatol=1e-4, fatol=1e-4))
            u, best = rr.x, float(rr.fun)
        return best, u

    def ajusta_sets(a, arranque):
        """Nuisance optimos de cada conjunto a logA fijo. Devuelve (chi2, th).

        DOS ARRANQUES por conjunto y se queda el mejor. El caliente (la solucion
        del punto anterior de la rejilla) es barato pero ARRASTRA: la version del
        2026-09-07 barria desde logA=2.60, que es el peor extremo, y llevaba esa
        solucion mala hasta el final — chi2 ~1050-1900 en toda la rejilla cuando
        un solo conjunto, arrancado en frio, ya daba 40 sobre 37 puntos. El frio
        (P0) no arrastra pero puede caer peor donde la superficie es dificil.
        Tomando el minimo de los dos no hay dependencia del camino."""
        c, ths = 0.0, []
        for st, u0 in zip(sets, arranque):
            f_cal, u_cal = _baja(st, a, u0)
            f_fri, u_fri = _baja(st, a, P0)
            if f_fri < f_cal:
                f_cal, u_cal = f_fri, u_fri
            c += f_cal; ths.append(u_cal)
        return c, ths

    # Rejilla fija y ancha en logA: no depende de ningun ajuste previo, asi que
    # no hay forma de que un minimo falso arrastre la ventana del perfil.
    grid = np.linspace(2.60, 3.40, 17)
    prof, arranque = [], [P0.copy() for _ in sets]
    for a in grid:
        c, arranque = ajusta_sets(a, arranque)
        prof.append(c)
        print(f'    logA={a:.4f}  chi2={c:.3f}', flush=True)
    prof = np.array(prof)
    i0 = int(np.argmin(prof))
    logA = float(grid[i0]); best = float(prof[i0])
    if 0 < i0 < len(grid) - 1:                 # parabola de los tres del minimo
        cf = np.polyfit(grid[i0-1:i0+2], prof[i0-1:i0+2], 2)
        sig_logA = float(np.sqrt(1.0 / cf[0])) if cf[0] > 0 else float('nan')
        logA = float(-cf[1] / (2 * cf[0]))
        best = float(np.polyval(cf, logA))
    else:
        sig_logA = float('nan')
        print('    AVISO: el minimo cae en el BORDE de la rejilla de logA')

    n = sum(st['npts'] for st in sets)
    nfree = 1 + 6 * len(sets)
    dof = n - nfree
    sc = np.sqrt(np.exp(logA) * 1e-10 / AS_REF)
    out = dict(model=name, chi2=best, dof=int(dof), npts=int(n),
               n_libres=int(nfree), logA=logA, sig_logA=sig_logA,
               kmax=KMAX, kmin=KMIN, modelado='LPT velocileptors (1 lazo + EFT)',
               Om=COSMO[name]['Om'], w0=COSMO[name]['w0'], wa=COSMO[name]['wa'],
               perfil=dict(logA=grid.tolist(), chi2=prof.tolist()),
               nuisance=[dict(zb=st['zb'], cap=st['cap'],
                              b1_lagrangiano=round(float(u[0] * ESC[0]), 4),
                              b1_euleriano=round(float(u[0] * ESC[0] + 1.0), 4),
                              b2=round(float(u[1] * ESC[1]), 4),
                              bs=round(float(u[2] * ESC[2]), 4),
                              alpha0=round(float(u[3] * ESC[3]), 3),
                              alpha2=round(float(u[4] * ESC[4]), 3),
                              SN0=round(float(u[5] * ESC[5]), 1))
                         for st, u in zip(sets, arranque)], bins=[])
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
