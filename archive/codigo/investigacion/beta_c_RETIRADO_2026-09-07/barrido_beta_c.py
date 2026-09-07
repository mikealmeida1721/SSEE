#!/usr/bin/env python3
"""
Barrido en beta_c: ¿donde deja de haber solucion, y por que?

LA PREGUNTA. El disparo autoconsistente (fondo_acoplado.py) devuelve
beta_c = +0.235068 (valor de 2026-08-10, ya con el fix del sqrt(3) en ALPHA;
el +0.494 de la corrida anterior era previo a ese fix), con signo contrario y
factor ~17 respecto del -AURA = -3.9978 que Paper 7 identifica. Y esa solucion
da 9.3% de energia oscura en z=9, muy por encima del limite observacional
(~3%). Al forzar beta_c hacia -AURA la integracion revienta con NaN. Hay que
distinguir dos cosas MUY distintas:

  (A) mi integrador se rompe  ->  limite numerico, hay que endurecerlo
  (B) no hay solucion fisica  ->  resultado del modelo

Este barrido las separa: para cada beta_c FIJO se dispara (phi_i, rho_c_i) con
las dos condiciones de hoy que vienen de los ingredientes, y se anota si
converge, si revienta, y cuanta energia oscura temprana queda.

COSTE. El disparo anterior usaba t_eval denso dentro del bucle del solver, que
es donde se iba el tiempo. Aqui se integra con dense_output y se evalua SOLO en
a=1 durante el ajuste; la malla completa se pide una vez, al final.
"""
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src/p07_eft')
import ssee_core as S
import fondo_acoplado as F

AURA = (3 * S.PHI + S.PI) / 2
LOG = '/home/mike/Proyectos/SSEE/results/logs/barrido_beta_c.log'


def _integra(phi_i, rho_c_i, bc, a_fin):
    """Integra sin malla densa; devuelve el objeto sol para evaluar donde haga falta."""
    return solve_ivp(F.rhs, [np.log(F.A_I), np.log(a_fin)],
                     [phi_i, 1e-8, rho_c_i], args=(bc,),
                     rtol=1e-9, atol=1e-12, method='Radau', dense_output=True)


def _en(sol, a, bc):
    """Estado del fondo en un a dado, a partir de la solucion densa."""
    phi, phip, rho_c = sol.sol(np.log(a))
    rho_bn, V, H2, X = F._estado(phi, phip, rho_c, a)
    if not np.isfinite(H2) or H2 <= 0:
        return None
    H = np.sqrt(H2)
    rho_phi = X / F.KAL + V
    w_phi = (X / F.KAL - V) / rho_phi
    Q = -bc * rho_c * (H * phip)
    tot = rho_phi + rho_c + rho_bn
    return dict(Om_DE=rho_phi / tot, rho_c=rho_c, H=H,
                w_phi=w_phi, w_eff=w_phi + Q / (3.0 * H * rho_phi))


def objetivo(u, bc):
    try:
        sol = _integra(u[0], np.exp(u[1]), bc, 1.02)
        if not sol.success:
            return [1e3, 1e3]
        e = _en(sol, 1.0, bc)
        if e is None:
            return [1e3, 1e3]
        v = [e['Om_DE'] - F.OM_PHI, e['rho_c'] - F.OM_C]
        return [1e3 if not np.isfinite(x) else x for x in v]
    except Exception:
        return [1e3, 1e3]


def resuelve(bc):
    l0 = np.log(F.OM_C * F.A_I ** -3)
    for p0 in (1.158, -3.543, 0.0, 5.0, -8.0):
        try:
            u, info, ier, msg = fsolve(objetivo, [p0, l0], args=(bc,),
                                       full_output=True, xtol=1e-11,
                                       maxfev=400)
        except Exception:
            continue
        if ier == 1 and max(abs(np.array(objetivo(u, bc)))) < 1e-7:
            return u
    return None


if __name__ == '__main__':
    f = open(LOG, 'w', buffering=1)          # linea a linea, no se pierde nada

    def p(s=''):
        print(s, flush=True)
        f.write(s + '\n')

    p('=' * 104)
    p('  BARRIDO EN beta_c — ¿limite numerico o ausencia de solucion fisica?')
    p('=' * 104)
    p(f'  objetivos en a=1:  Om_DE={F.OM_PHI:.6f}   rho_c={F.OM_C:.6f}')
    p(f'  -AURA = {-AURA:.6f}     Paper 7 (orden dominante) = -3.989910')
    p('')
    p(f'  {"beta_c":>10} {"w_eff(a=1)":>12} {"w_phi(a=1)":>12} '
      f'{"Om_DE(z=9)":>12} {"alpha_K":>10} {"seg":>6}  estado')
    p('  ' + '-' * 98)
    # El disparo autoconsistente cae en +0.235068 (POSITIVO), asi que el
    # barrido tiene que cubrir tambien ese lado o no veria su propia respuesta.
    for bc in (1.0, 0.5, 0.235068, 0.1, 0.0,
               -0.1, -0.2, -0.3, -0.4, -0.494, -0.6, -0.8, -1.0,
               -1.5, -2.0, -2.5, -3.0, -3.5, -3.9978, -5.0):
        t0 = time.time()
        u = resuelve(bc)
        dt = time.time() - t0
        if u is None:
            p(f'  {bc:>10.4f} {"—":>12} {"—":>12} {"—":>12} {"—":>10} '
              f'{dt:>6.1f}  SIN SOLUCION (no converge)')
            continue
        try:
            sol = _integra(u[0], np.exp(u[1]), bc, 1.02)
            e1 = _en(sol, 1.0, bc)
            e9 = _en(sol, 0.1, bc)
            aK = 3 * e1['Om_DE'] * (1 + e1['w_eff'])
            ede = e9['Om_DE'] if e9 else float('nan')
            est = 'OK' if ede < 0.03 else f'EXCLUIDO: DE temprana {100*ede:.1f}%'
            p(f'  {bc:>10.4f} {e1["w_eff"]:>+12.6f} {e1["w_phi"]:>+12.6f} '
              f'{ede:>12.6f} {aK:>10.6f} {dt:>6.1f}  {est}')
        except Exception as ex:
            p(f'  {bc:>10.4f} {"—":>12} {"—":>12} {"—":>12} {"—":>10} '
              f'{dt:>6.1f}  REVIENTA: {type(ex).__name__}')
    p('')
    p('  LECTURA: si converge hasta cierto beta_c y luego falla -> limite numerico.')
    p('           si converge en todo el rango pero todas las de |beta_c| grande')
    p('           dan DE temprana > 3% -> beta_c = -AURA excluido por el CMB.')
    f.close()
