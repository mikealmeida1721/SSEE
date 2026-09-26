#!/usr/bin/env python3
"""
is_growth_gamma_bg.py — gamma_bg del fondo IS, recomputado con la materia REAL.

QUE ARREGLA (hallazgo F1, 2026-09-26)
-------------------------------------
`manuscript/SSEE_EFT_section.tex` eq.(growth_IS) l.439 y el fondo h~^2(a) l.444
ponen Omega_{m,dyn} = 0.160050 como densidad de materia, en DOS ranuras:

    delta'' + [2 + dln h~/dln a] delta' = (3/2) * Omega_{m,dyn} a^-3 / h~^2 * delta
    h~^2(a) = Omega_{m,dyn} a^-3 + rho_DE(a)

Pero 0.160050 es s_m = 1 + w_0: un numero de la ECUACION DE ESTADO, no una
densidad. Es el mismo error de categoria que el bug Omega_m-geometria del
2026-07-09 y que la resta que sostenia la particula. La densidad de materia del
modelo es Omega_m = omega_m/h^2 = 0.308881, la MISMA en fondo y en crecimiento.

De esa ranura mal llenada salen gamma_bg = 0.657 +/- 0.002 y S_8 = 0.837, que
el .tex presenta con numero de ecuacion propio. Este script los re-mide.

QUE MIDE. La ecuacion de crecimiento lineal para CDM sobre el fondo del modelo,
en e-folds (x = ln a):

    delta'' + [2 + dln E/dln x] delta' = (3/2) * Omega_m(a) * delta
    Omega_m(a) = Omega_m a^-3 / E^2(a)

y luego el indice gamma del ajuste f(a) = Omega_m(a)^gamma sobre 0.1 <= a <= 1,
que es la definicion que usa el .tex.

NO es gamma_IS. Este es el indice del fondo SIN realimentacion viscosa en la
ecuacion de delta; el completo lo da Paper 5 (gamma_IS = 0.5504 +/- 0.001).
Aqui solo se corrige QUE DENSIDAD entra, no el formalismo.

CONTROL (R53). Se corre lo mismo tres veces y las tres tienen que dar lo que se
espera, o la medicion no mide nada:
  (a) el fondo BUGGY (Omega_m = 0.160050) debe REPRODUCIR el 0.657 del .tex
      -- si no lo reproduce, el diagnostico del bug es falso;
  (b) el fondo CORRECTO (Omega_m = 0.308881, w0/wa del algebra);
  (c) LCDM (Omega_m = 0.3153, w = -1), cuyo gamma tiene que salir ~0.55, el
      valor de libro. Si (c) falla, el integrador esta roto y (a) y (b) no
      valen nada.

Salida: results/logs/is_growth_gamma_bg.json
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ssee_core import OMEGA_M_TOTAL, S_M, W0, WA          # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(REPO, "results", "logs", "is_growth_gamma_bg.json")

# ORIGEN-VALOR: 0.3153 — Omega_m de Planck 2018 TT,TE,EE+lowE+lensing
#               (arXiv:1807.06209, tabla 2, columna 6)
OM_PLANCK = 0.3153

A_INI, A_FIN = 1.0e-3, 1.0
A_FIT_MIN = 0.1                      # ventana de ajuste del .tex: 0.1 <= a <= 1


def E2(a, om, w0, wa):
    """E^2(a) con energia oscura CPL y plano por construccion."""
    ode = 1.0 - om
    f_de = a ** (-3.0 * (1.0 + w0 + wa)) * math.exp(-3.0 * wa * (1.0 - a))
    return om * a ** -3 + ode * f_de


def dlnE_dlna(a, om, w0, wa, h=1e-5):
    x = math.log(a)
    return (math.log(math.sqrt(E2(math.exp(x + h), om, w0, wa)))
            - math.log(math.sqrt(E2(math.exp(x - h), om, w0, wa)))) / (2 * h)


def crece(om, w0, wa):
    """Integra delta'' + [2 + dlnE/dlna] delta' = (3/2) Om(a) delta."""
    def rhs(x, y):
        a = math.exp(x)
        d, dp = y
        om_a = om * a ** -3 / E2(a, om, w0, wa)
        return [dp, -(2.0 + dlnE_dlna(a, om, w0, wa)) * dp + 1.5 * om_a * d]

    xs = np.log(np.geomspace(A_INI, A_FIN, 4000))
    # En materia dominante delta ~ a  =>  delta = a_ini, delta' = a_ini.
    sol = solve_ivp(rhs, (xs[0], xs[-1]), [A_INI, A_INI], t_eval=xs,
                    rtol=1e-10, atol=1e-14, method="DOP853")
    a = np.exp(sol.t)
    d, dp = sol.y[0], sol.y[1]
    f = dp / d                                     # f = dln delta / dln a
    om_a = om * a ** -3 / np.array([E2(_a, om, w0, wa) for _a in a])
    return a, f, om_a


def gamma_de(om, w0, wa):
    """gamma del ajuste f = Om(a)^gamma en la ventana del .tex."""
    a, f, om_a = crece(om, w0, wa)
    m = (a >= A_FIT_MIN) & (f > 0) & (om_a > 0)
    g = np.log(f[m]) / np.log(om_a[m])
    return float(g.mean()), float(g.std())


def main():
    casos = {
        "buggy_s_m_como_densidad": dict(om=S_M, w0=W0, wa=WA,
                                        nota="Omega_m = s_m = 1+w0 = 0.160050 "
                                             "(la ranura mal llenada del .tex)"),
        "correcto_omega_m_total": dict(om=OMEGA_M_TOTAL, w0=W0, wa=WA,
                                       nota="Omega_m = omega_m/h^2 = 0.308881 "
                                            "(la densidad de verdad)"),
        "control_lcdm": dict(om=OM_PLANCK, w0=-1.0, wa=0.0,
                             nota="Planck 2018; gamma de libro ~0.55"),
    }
    res = {}
    for k, c in casos.items():
        g, s = gamma_de(c["om"], c["w0"], c["wa"])
        res[k] = dict(Omega_m=c["om"], w0=c["w0"], wa=c["wa"],
                      gamma=g, gamma_sd=s, nota=c["nota"])
        print(f"  {k:28s} Om={c['om']:.6f}  gamma = {g:.4f} +/- {s:.4f}")

    # Control (R53): los tres tienen que caer donde se espera.
    ok_buggy = abs(res["buggy_s_m_como_densidad"]["gamma"] - 0.657) < 0.02
    ok_lcdm = abs(res["control_lcdm"]["gamma"] - 0.55) < 0.02
    res["control"] = dict(
        buggy_reproduce_el_0657=bool(ok_buggy),
        lcdm_da_el_055_de_libro=bool(ok_lcdm),
        veredicto=("el integrador es fiable y el 0.657 SI venia de meter s_m "
                   "en la ranura de densidad"
                   if (ok_buggy and ok_lcdm) else
                   "NO concluyente: revisar antes de tocar el .tex"))
    print(f"\n  control: buggy reproduce 0.657 -> {ok_buggy}")
    print(f"  control: LCDM da 0.55 de libro  -> {ok_lcdm}")
    print(f"  {res['control']['veredicto']}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  -> {OUT}")


if __name__ == "__main__":
    main()
