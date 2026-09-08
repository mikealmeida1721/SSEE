"""Perfil de chi2(w_c) contra Planck crudo — la mitad que faltaba (cola #2).

QUE CONTESTA. El Registro cita w_c = 0.119534 +/- 0.000248 del CMB. Ese error
es la barra que hace que la identidad KAL0*w_b*n_s quede a 0.08 sigma. Pero
hasta ahora nadie habia comprobado que la barra la ponga EL DATO y no la
parametrizacion: si el ajuste estuviera sobre-restringido por como esta
escrito el modelo, el +-0.000248 seria artefacto y el 0.08 sigma no
significaria nada.

COMO. Se barre w_c y en cada punto se minimiza sobre logA y tau, con el resto
del fondo fijo. La curvatura de la parabola en el minimo da el error.

CONTROL (R53). El mismo perfil con el fondo de LCDM. Si la barra la pone el
dato, los dos modelos deben salir con anchuras del mismo orden: Planck
restringe w_c con la misma fuerza se le ponga el fondo que se le ponga. Si
SSEE saliera mucho mas estrecho que LCDM, la precision vendria del algebra y
no del dato, y el 0.08 sigma habria que retirarlo.

FIX que hizo falta antes de poder correr esto: el evaluador traia la ecuacion
de estado de SSEE clavada como literal, asi que la fila LCDM habria salido
mal en silencio (ya paso una vez, ver cmb_tau_flotado.json). Ahora w y wa son
argumentos obligatorios. Ver src/p03_cmb/cmb_eval.py y R64.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src/p03_cmb')

import ssee_core as S
from cmb_eval import chi2_y_s8

# fondo de cada modelo: todo fijo salvo el w_c que se barre
FONDOS = {
    'SSEE': dict(ombh2=S.OMEGA_B_H2, H0=S.H0_GLOBAL, ns=S.N_S,
                 w=S.W0, wa=S.WA, centro=0.119514),   # R66-OK
    'LCDM': dict(ombh2=0.02237, H0=67.36, ns=0.9649,
                 w=-1.0, wa=0.0, centro=0.1200),
}
# R66-OK arriba: `centro` NO entra en el calculo, dice DONDE se pone la
# rejilla. Se deja al valor redondeado a proposito, para que la rejilla sea
# la misma que la de la corrida ya hecha; moverlo 4e-07 no cambia la fisica
# pero si dejaria el log sin corresponder con su fuente.
PASO, N = 0.0002, 9          # +-0.0008 alrededor del centro


def mejor(f, wc):
    """chi2 minimo sobre (logA, tau) con este w_c y este fondo."""
    def obj(u):
        if not (1.5 < u[0] < 4.5 and 0.011 < u[1] < 0.16):
            return 1e9
        return chi2_y_s8(dict(ombh2=f['ombh2'], H0=f['H0'], ns=f['ns'],
                              omch2=wc, logA=u[0], tau=u[1]),
                         f['w'], f['wa'])[0]
    r = minimize(obj, [3.044, 0.054], method='Nelder-Mead',
                 options=dict(xatol=1e-5, fatol=1e-4, maxiter=400))
    return float(r.fun)


def perfil(nombre):
    f = FONDOS[nombre]
    rej = f['centro'] + PASO * (np.arange(N) - (N - 1) // 2)
    print(f"\n=== {nombre}  (w={f['w']:.6f}, wa={f['wa']:.6f}) ===",
          flush=True)
    print(' w_c          chi2', flush=True)
    ys = []
    for wc in rej:
        t0 = time.time()
        c = mejor(f, wc)
        ys.append(c)
        print(f'{wc:.6f}  {c:11.4f}   [{time.time()-t0:.0f}s]', flush=True)
    ys = np.asarray(ys)
    cf = np.polyfit(rej, ys, 2)
    wc0 = -cf[1] / (2 * cf[0])
    sg = np.sqrt(1.0 / cf[0]) if cf[0] > 0 else np.nan
    print(f'  parabola: w_c = {wc0:.6f} +- {sg:.6f}', flush=True)
    if nombre == 'SSEE':
        print(f'  identidad KAL0*w_b*n_s = {S.OMEGA_C_H2:.6f}  ->  '
              f'{abs(wc0-S.OMEGA_C_H2)/sg:.2f} sigma', flush=True)
    return dict(wc=rej.tolist(), chi2=ys.tolist(),
                wc0=float(wc0), sigma=float(sg))


if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '1'
    out = {n: perfil(n) for n in ('SSEE', 'LCDM')}
    r = out['SSEE']['sigma'] / out['LCDM']['sigma']
    print(f"\nCONTROL (R53): sigma_SSEE / sigma_LCDM = {r:.3f}", flush=True)
    print("  cerca de 1 => la barra la pone el DATO, no el algebra.",
          flush=True)
    print("  muy < 1     => SSEE sale artificialmente estrecho: el "
          "+-0.000248 seria artefacto y el 0.08 sigma habria que retirarlo.",
          flush=True)
    out['control_razon_sigmas'] = r
    json.dump(out, open('/home/mike/Proyectos/SSEE/results/logs/'
                        'cmb_perfil_wc.json', 'w'), indent=1)
