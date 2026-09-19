#!/usr/bin/env python3
"""
Análisis de la cadena R4 (LCDM, control metodológico) — sigma8/Omega_m/S8.

Por qué existe: el chain de Cobaya para LCDM NO lleva sigma8 como columna
derivada (a diferencia del truco sqrt(As) que usó R3 para SSEE, aquí el
fondo es LIBRE — ombh2/omch2/h0/ns varían por muestra, así que sigma8 hay
que recomputarlo por muestra vía CAMB). Con ~47000 filas post burn-in a
~0.77s/llamada (CAMB lineal, sin halofit) eso son ~10h — se submuestrea
sistemáticamente (cada fila post-burn-in ya es un punto ACEPTADO único con
su peso; tomar 1 de cada K no sesga la media, solo sube la varianza Monte
Carlo, que con miles de puntos queda pequeña).

Mismo criterio de burn-in que R3 (30% por cadena) para comparar manzanas
con manzanas.
"""
# ORIGEN-VALOR: 0.025570 — R-1 de las medias de R4, results/logs/R4_lcdm_resume_20260805.log
import json
import time
import numpy as np
import camb

CHAINS_DIR = '/mnt/datos/SSEE_data/chains_p6/kids'
BURN_IN_FRAC = 0.30
THIN_EVERY = 8          # 1 de cada 8 filas post burn-in
MNU = 0.06              # fiducial de esta corrida (fijo en el yaml)


def sigma8_de(ombh2, omch2, h0, ns, As):
    p = camb.CAMBparams()
    p.set_cosmology(H0=h0 * 100.0, ombh2=ombh2, omch2=omch2, mnu=MNU, omk=0.0)
    p.set_dark_energy(w=-1.0, wa=0.0, dark_energy_model='ppf')
    p.InitPower.set_params(As=As, ns=ns)
    p.set_matter_power(redshifts=[0.0], kmax=2.0, nonlinear=False)
    r = camb.get_results(p)
    return float(r.get_sigma8()[-1])


def carga_cadena(path):
    with open(path) as f:
        header = f.readline().lstrip('#').split()
    data = np.loadtxt(path)
    cols = {name: data[:, i] for i, name in enumerate(header)}
    return cols


def main():
    t0 = time.time()
    todas_w, todas_om, todas_s8 = [], [], []
    n_por_cadena = []
    # chi2_min por cadena. Se guarda porque Paper 6 lo PUBLICA (fila LCDM de la
    # tabla S8 y el Delta chi2 de la Eq. dchi2) y la primera version de este
    # script no lo escribia: el numero vivia solo en la cadena, y el log no lo
    # respaldaba. Recuperado a mano el 2026-09-07; desde aqui sale solo.
    chi2_min_cad = []

    for i in range(1, 5):
        path = f'{CHAINS_DIR}/lcdm.{i}.txt'
        cols = carga_cadena(path)
        n = len(cols['weight'])
        i0 = int(n * BURN_IN_FRAC)
        idx = np.arange(i0, n, THIN_EVERY)
        n_por_cadena.append(len(idx))
        chi2_min_cad.append(float(cols['chi2'][i0:].min()))
        print(f'cadena {i}: {n} filas, burn-in {i0}, {len(idx)} submuestreadas')

        for j in idx:
            ombh2 = cols['ombh2'][j]
            omch2 = cols['omch2'][j]
            h0 = cols['h0'][j]
            ns = cols['ns'][j]
            logA = cols['logA'][j]
            As = np.exp(logA) * 1e-10
            w = cols['weight'][j]
            try:
                s8 = sigma8_de(ombh2, omch2, h0, ns, As)
            except Exception as e:
                print('  fallo CAMB en fila', j, ':', e)
                continue
            omnuh2 = MNU / 93.14
            Om = (ombh2 + omch2 + omnuh2) / h0**2
            todas_w.append(w)
            todas_om.append(Om)
            todas_s8.append(s8)  # sigma8; S8 = sigma8*sqrt(Om/0.3) se arma abajo

    todas_w = np.array(todas_w)
    todas_om = np.array(todas_om)
    sigma8_arr = np.array(todas_s8)
    S8_arr = sigma8_arr * np.sqrt(todas_om / 0.3)

    def wmean(x, w):
        return float(np.sum(w * x) / np.sum(w))

    def wstd(x, w):
        m = wmean(x, w)
        return float(np.sqrt(np.sum(w * (x - m) ** 2) / np.sum(w)))

    resultado = {
        'corrida': 'R4 — LCDM control metodologico, KiDS-1000',
        'muestreador': 'Cobaya mcmc, 4 cadenas MPI, reanudada tras corte de luz',
        'n_parametros_libres': 13,
        'burn_in_frac': BURN_IN_FRAC,
        'thin_every': THIN_EVERY,
        'n_muestras_usadas': int(len(todas_w)),
        'n_filas_por_cadena': n_por_cadena,
        'convergencia': {
            'Rminus1_medias': 0.025570,
            'Rminus1_colas': 0.115482,
            'converged': True,
        },
        'dof': 212,
        'chi2_min': round(min(chi2_min_cad), 5),
        'chi2_min_por_dof': min(chi2_min_cad) / 212.0,
        'chi2_min_por_cadena': [round(c, 3) for c in chi2_min_cad],
        'Omega_m': {'media': wmean(todas_om, todas_w), 'sigma': wstd(todas_om, todas_w)},
        'sigma8': {'media': wmean(sigma8_arr, todas_w), 'sigma': wstd(sigma8_arr, todas_w)},
        'S8': {'media': wmean(S8_arr, todas_w), 'sigma': wstd(S8_arr, todas_w)},
        'comparacion_KiDS': {
            'S8_publicado': 0.759,
            'err_publicado': 0.024,
        },
        'tiempo_analisis_seg': time.time() - t0,
    }
    S8m, S8s = resultado['S8']['media'], resultado['S8']['sigma']
    tension = abs(S8m - 0.759) / np.sqrt(S8s**2 + 0.024**2)
    resultado['comparacion_KiDS']['tension_sigma'] = float(tension)

    out = '/home/mike/Proyectos/SSEE/results/logs/growth_2026-07/R4_lcdm_kids_S8.json'
    with open(out, 'w') as f:
        json.dump(resultado, f, indent=2, ensure_ascii=False)

    print(json.dumps(resultado, indent=2, ensure_ascii=False))
    print(f'\nGuardado en {out}')


if __name__ == '__main__':
    main()
