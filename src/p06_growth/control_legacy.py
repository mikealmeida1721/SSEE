#!/usr/bin/env python3
"""CONTROL NEGATIVO KiDS-Legacy (R53) — antes de creerle nada al pipeline.

Reproducir el chi^2 del punto de MAXIMA VEROSIMILITUD de la cadena oficial
`output_nautilus_xipm_Fiducial.txt` (Wright et al. 2025, A&A 703, A158):

    max loglike = -203.823490  ->  chi^2_ref = 407.647   (357 puntos, 20 libres)

Convencion verificada: CosmoSIS reporta loglike = -0.5*chi^2, SIN el termino
-0.5*ln|2 pi C|. Se comprobo contra KiDS-1000, donde -2*(-130.157350) = 260.31
casa con el chi^2 = 260.32 publicado (el `like_full` que calculaba
`control_kids.py` daba 5819.7 y era decorativo: nunca fue lo que casaba).

DIFERENCIAS DECLARADAS de este evaluador contra el pipeline oficial — son el
motivo por el que NO se espera coincidencia exacta, y se cuantifican aqui:
  (1) No lineal: RESUELTA. CAMB 1.6.5 trae `mead2020_feedback` con
      HMCode_logT_AGN, que es el mismo modelo que el oficial sirve por
      emulador. Se corre en el log_T_AGN OFICIAL, sin traduccion libre: ese
      es el titular. El barrido de HMcode-2015 (halo_A) se conserva abajo
      solo como control del otro lado (R53): dice cuanto de la diferencia
      venia del modelo no lineal y cuanto del resto.
  (2) IA: oficial = NLA-M (amplitud A, pendiente beta, 6 masas medias con
      priors gaussianos correlacionados). Aqui se EVALUA ese modelo en el
      punto oficial y se pasa su salida como amplitud POR BIN:
          A_eff,i = A * f_r,i * (M_i/M_piv)^beta
      Ya no es "NLA plano vs NLA-M": es el mismo kernel con las 6 amplitudes.
  (3) Limber extendido + binning simple en theta; el oficial usa
      `bin_xi` con pesos npairs medidos.
"""
import json
import sys
import time

import numpy as np

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src/p06_growth')
import kids_shear as K                                          # noqa: E402

CHI2_REF = 407.6470                # -2 * max loglike de la cadena oficial
LOGLIKE_REF = -203.823490

# --- punto de maxima verosimilitud de la cadena oficial (xipm Fiducial) ---
MP = dict(omch2=0.095019, ombh2=0.0191914, h0=0.769441, ns=1.01013,
          logA=3.7431053611421596, log_t_agn=8.21185,
          S8_ref=0.84347, sigma8_ref=1.0518473627303115,
          Om_ref=0.193999174888498)
DZ = np.array([0.031198, 0.025115, -0.008245, -0.007295, 0.024564, 0.026880])

# --- NLA-M evaluado en ese mismo punto -> amplitud efectiva por bin ---
IA_A = 5.78122
IA_BETA = 0.476046
LOG10_M_PIV = 13.5                                  # values.ini
F_R = np.array([0.158, 0.198, 0.206, 0.258, 0.207, 0.026])       # values.ini
LOG10_M = np.array([11.4971, 12.1786, 12.4623, 12.6081, 12.7464, 12.8492])
A_EFF = IA_A * F_R * 10.0 ** (IA_BETA * (LOG10_M - LOG10_M_PIV))

# HMcode-2015 no tiene log_T_AGN: se barre su unico parametro y se reporta
# la curva entera, no solo el minimo (la traduccion es el resultado, no un
# detalle de implementacion).
HALO_A = [1.2, 1.5, 1.8, 2.0, 2.4, 2.8, 3.13, 3.5]
LOGT_AGN = [7.3, 7.8, 8.0, MP['log_t_agn'], 8.3]


def main():
    cfg = K.set_dataset('legacy')
    D = K.load_data()
    m = K.scale_mask(D)
    As = np.exp(MP['logA']) / 1e10
    print(f'dataset   : {cfg["DATASET"]}  bins={cfg["NZBINS"]}  '
          f'cortes xi+{cfg["KEEP_XIP"]} xi-{cfg["KEEP_XIM"]}', flush=True)
    print(f'puntos    : {int(m.sum())} de {len(D["d"])}', flush=True)
    print(f'A_IA/bin  : {np.round(A_EFF, 4).tolist()}', flush=True)
    print(f'As        : {As:.6e}  (logA={MP["logA"]:.6f})', flush=True)
    print(f'referencia: chi2 = {CHI2_REF:.4f}', flush=True)
    print()

    def evalua(**kw):
        t0 = time.time()
        res, p, kh, zpk, pk, gr = K.run_camb(
            As=As, omch2=MP['omch2'], ombh2=MP['ombh2'], h0=MP['h0'],
            ns=MP['ns'], **kw)
        s8 = float(res.get_sigma8_0())
        Om = float((p.omch2 + p.ombh2 + p.omnuh2) / p.h ** 2)
        ells, Cl, idx = K.cl_shear(D, res, p, kh, zpk, pk, gr, A_EFF, DZ)
        c2 = K.chi2(D, K.theory_vector(D, ells, Cl, idx, delta_c=0.0), m)
        return dict(kw, chi2=c2, sigma8=s8, Om=Om, S8=s8 * np.sqrt(Om / 0.3),
                    desvio_pct=100.0 * (c2 - CHI2_REF) / CHI2_REF,
                    segundos=round(time.time() - t0, 1))

    print('--- TITULAR: HMCode-2020, el MISMO modelo no lineal del oficial ---')
    hm20 = []
    for lT in LOGT_AGN:
        f = evalua(logT_AGN=lT)
        hm20.append(f)
        marca = '  <-- log_T_AGN OFICIAL' if abs(lT - MP['log_t_agn']) < 1e-9 else ''
        print(f'  logT_AGN={lT:7.5f}  chi2={f["chi2"]:9.3f}  '
              f'({f["desvio_pct"]:+6.2f}% vs ref)  sigma8={f["sigma8"]:.5f}'
              f'{marca}', flush=True)
    oficial = [f for f in hm20 if abs(f['logT_AGN'] - MP['log_t_agn']) < 1e-9][0]

    print()
    print('--- CONTROL DEL OTRO LADO (R53): HMcode-2015, que hay que traducir ---')
    salida = []
    for hA in HALO_A:
        f = evalua(halo_A=hA)
        salida.append(f)
        print(f'  halo_A={hA:4.2f}  chi2={f["chi2"]:9.3f}  '
              f'({f["desvio_pct"]:+6.2f}% vs ref)  sigma8={f["sigma8"]:.5f}',
              flush=True)

    mejor = min(salida, key=lambda f: abs(f['chi2'] - CHI2_REF))
    print()
    print(f'VEREDICTO  chi2 = {oficial["chi2"]:.3f} contra {CHI2_REF:.3f} '
          f'-> {oficial["desvio_pct"]:+.2f}%, con TODOS los ingredientes en su '
          f'valor oficial y ningun libre de traduccion.')
    print(f'  control HMcode-2015: mas cerca en halo_A={mejor["halo_A"]} '
          f'({mejor["desvio_pct"]:+.2f}%) -- ahi la traduccion SI es libre, '
          f'por eso no es el titular.')
    print('sigma8/Om/S8 de referencia de la cadena: '
          f'{MP["sigma8_ref"]:.5f} / {MP["Om_ref"]:.5f} / {MP["S8_ref"]:.5f}')

    log = dict(
        prueba='control negativo KiDS-Legacy xipm (R53)',
        fecha='2026-09-19',
        fuente_datos=K.DATA,
        cadena_oficial=('KiDS_Legacy_cosmic_shear_data_release/'
                        'chains_and_config_files/xipm/'
                        'output_nautilus_xipm_Fiducial.txt'),
        loglike_ref=LOGLIKE_REF, chi2_ref=CHI2_REF,
        convencion='CosmoSIS loglike = -0.5*chi2, sin -0.5*ln|2piC|',
        n_puntos=int(m.sum()), n_libres_oficial=20,
        punto=MP, dz=DZ.tolist(), A_IA_por_bin=A_EFF.tolist(),
        titular_hmcode2020=oficial, barrido_logT_AGN=hm20,
        barrido_halo_A_control=salida, mas_cerca_halo_A=mejor,
        diferencias_declaradas=[
            'no lineal: RESUELTA -- CAMB mead2020_feedback(logT_AGN) = el oficial',
            'IA: se pasa la salida del NLA-M oficial como amplitud por bin',
            'binning en theta sin pesos npairs medidos',
        ])
    ruta = ('/home/mike/Proyectos/SSEE/results/logs/growth_2026-07/'
            'control_legacy.json')
    with open(ruta, 'w') as f:
        json.dump(log, f, indent=2)
    print(f'log -> {ruta}')


if __name__ == '__main__':
    main()
