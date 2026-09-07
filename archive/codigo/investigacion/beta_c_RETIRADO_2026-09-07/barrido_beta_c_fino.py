#!/usr/bin/env python3
"""
Barrido FINO en beta_c entre +0.10 y +0.24: ¿existe un punto que cumpla las
DOS condiciones a la vez?

LA PREGUNTA. El barrido grueso (barrido_beta_c.log, 2026-09-05) dejo el asunto
en un hueco sin explorar:

    beta_c = +0.2351  ->  w_eff = -0.839949 (= w0 EXACTO)  pero DE(z=9) = 4.5%
    beta_c = +0.1000  ->  w_eff = -0.892879 (lejos de w0)  pero DE(z=9) = 1.2%
    limite observacional de energia oscura temprana: ~3%

O sea: el acoplamiento que da el w0 correcto cuesta demasiada energia oscura
temprana, y el que es barato en energia temprana da el w0 equivocado. Entre
+0.10 y +0.24 NO SE HABIA CORRIDO NADA. Este barrido lo cubre.

POR QUE IMPORTA. Si algun punto cumple las dos, no hay problema abierto: beta_c
queda determinado y Paper 7 se escribe con un numero. Si el barrido sale vacio,
ENTONCES si hay un problema abierto, y ademas queda dicho con precision: dentro
del acoplamiento CONFORMAL no existe beta_c que reproduzca w0 sin violar el
limite de energia oscura temprana. Esa es la premisa a atacar despues (Paper 8
ya usa acoplamiento DISFORMAL, Paper 7 conformal).

CONTROL. El barrido incluye +0.1000 y +0.235068, que ya tienen respuesta
conocida del barrido grueso. Si no las reproduce, el barrido esta roto y no hay
que creerle nada de lo demas.
"""
import sys
import time

import numpy as np

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src/p07_eft')
import ssee_core as S
import barrido_beta_c as B

W0 = S.W0                     # -0.839950
LIM_EDE = 0.03   # ⚠️ SIN FUENTE — puesto por Claude, NO citado por
                 # ningun paper ni por CANONICAL_VALUES. Compara la cantidad
                 # equivocada en la epoca equivocada. El criterio VALIDO es
                 # omega_m en recombinacion (ver OP-23). Conservado solo para
                 # reproducir las corridas del 2026-09-05.
LOG = '/home/mike/Proyectos/SSEE/results/logs/barrido_beta_c_fino.log'

# Los dos puntos de control, con su valor conocido del barrido grueso.
CONTROL = {0.100000: (-0.892879, 0.011928),
           0.235068: (-0.839949, 0.045301)}

if __name__ == '__main__':
    f = open(LOG, 'w', buffering=1)

    def p(s=''):
        print(s, flush=True)
        f.write(s + '\n')

    p('=' * 96)
    p('  BARRIDO FINO beta_c en [+0.10, +0.24] — ¿un punto cumple las DOS?')
    p('=' * 96)
    p(f'  objetivo 1:  w_eff(a=1) = w0 = {W0:.6f}')
    p(f'  objetivo 2:  Om_DE(z=9) < {LIM_EDE:.2f}')
    p('')
    p(f'  {"beta_c":>9} {"w_eff(a=1)":>12} {"|w-w0|":>10} '
      f'{"DE(z=9)":>10} {"alpha_K":>10} {"seg":>6}  estado')
    p('  ' + '-' * 88)

    rejilla = sorted(set([round(x, 6) for x in np.arange(0.10, 0.2451, 0.005)]
                         + list(CONTROL)))
    filas, ctrl_ok = [], True
    for bc in rejilla:
        t0 = time.time()
        u = B.resuelve(bc)
        dt = time.time() - t0
        if u is None:
            p(f'  {bc:>9.4f} {"—":>12} {"—":>10} {"—":>10} {"—":>10} '
              f'{dt:>6.1f}  SIN SOLUCION')
            continue
        sol = B._integra(u[0], np.exp(u[1]), bc, 1.02)
        e1, e9 = B._en(sol, 1.0, bc), B._en(sol, 0.1, bc)
        w, ede = e1['w_eff'], (e9['Om_DE'] if e9 else float('nan'))
        aK = 3 * e1['Om_DE'] * (1 + w)
        dw = abs(w - W0)
        # "cumple" = w0 dentro de 0.005 (el ancho con que P7 cita w0) Y DE<3%
        if dw < 0.005 and ede < LIM_EDE:
            est = '*** CUMPLE LAS DOS ***'
        elif ede < LIM_EDE:
            est = f'DE ok, w lejos ({dw:+.4f})'
        else:
            est = f'EXCLUIDO: DE {100*ede:.1f}%'
        p(f'  {bc:>9.4f} {w:>+12.6f} {dw:>10.6f} {ede:>10.6f} '
          f'{aK:>10.6f} {dt:>6.1f}  {est}')
        filas.append((bc, w, dw, ede))

        if bc in CONTROL:
            w_ref, ede_ref = CONTROL[bc]
            ok = abs(w - w_ref) < 2e-5 and abs(ede - ede_ref) < 2e-5
            ctrl_ok &= ok
            p(f'  {"":>9} CONTROL vs barrido grueso: w {w_ref:+.6f} / '
              f'DE {ede_ref:.6f}  ->  {"REPRODUCE" if ok else "NO REPRODUCE"}')

    p('')
    p('  ' + '=' * 88)
    if not ctrl_ok:
        p('  CONTROL FALLIDO: el barrido no reproduce los puntos conocidos.')
        p('  No creerle nada a las filas de arriba.')
    else:
        cumplen = [r for r in filas if r[2] < 0.005 and r[3] < LIM_EDE]
        p(f'  Control: reproduce los {len(CONTROL)} puntos conocidos.')
        if cumplen:
            p(f'  RESULTADO: {len(cumplen)} punto(s) cumplen las dos. '
              'beta_c queda DETERMINADO, no hay problema abierto.')
            for bc, w, dw, ede in cumplen:
                p(f'     beta_c = {bc:.4f}   w_eff = {w:+.6f}   DE = {100*ede:.2f}%')
        else:
            p('  RESULTADO: NINGUN punto cumple las dos.')
            p('  El problema es REAL, y queda dicho con precision: dentro del')
            p('  acoplamiento CONFORMAL no existe beta_c que reproduzca w0 sin')
            p('  pasarse del limite de energia oscura temprana.')
            if filas:
                mej = min(filas, key=lambda r: r[2] if r[3] < LIM_EDE else 9e9)
                if mej[3] < LIM_EDE:
                    p(f'  Lo mas cerca de w0 con DE<3%: beta_c={mej[0]:.4f}, '
                      f'w_eff={mej[1]:+.6f}, a {mej[2]:.4f} de w0 (DE {100*mej[3]:.2f}%)')
    f.close()
