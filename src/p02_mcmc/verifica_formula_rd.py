#!/usr/bin/env python3
"""¿De donde sale la formula de r_d de Paper 2 (ec. eq:rd)?

POR QUE EXISTE (2026-09-19, FP-7). Paper 2 usa
    r_d = 147.27 (w_m/0.1432)^-0.255 (w_b/0.02237)^-0.134 Mpc   (pivote Planck 2018)
y la citaba como «the Eisenstein & Hu 1998 fitting formula». EH98 no es una ley
de potencias. Esto mide que es en realidad: se derivan los exponentes LOCALES
de EH98 (formula completa, ecs. 2-6 de EH98) y de CAMB en el pivote, y se
compara con la formula del paper y con el ajuste global de Aubourg et al. 2015
(arXiv:1411.1074), que usa DESI.

Se corre solo:  python3 src/p02_mcmc/verifica_formula_rd.py
"""
import json
import pathlib

import camb
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
SALIDA = REPO / "results" / "logs" / "formula_rd_exponentes.log"
WM0, WB0 = 0.1432, 0.02237     # pivote de la ec. eq:rd (Planck 2018)
PASO = 0.02                    # derivada centrada a +-2 %


def eh98(wm, wb, T=2.7255):
    """r_s(z_d) de Eisenstein & Hu 1998, ecs. 2-6 (ApJ 496, 605)."""
    th = T / 2.7
    zeq = 2.5e4 * wm * th**-4
    keq = 7.46e-2 * wm * th**-2
    b1 = 0.313 * wm**-0.419 * (1 + 0.607 * wm**0.674)
    b2 = 0.238 * wm**0.223
    zd = 1291 * wm**0.251 / (1 + 0.659 * wm**0.828) * (1 + b1 * wb**b2)
    R = lambda z: 31.5 * wb * th**-4 * (1e3 / z)
    Rd, Req = R(zd), R(zeq)
    return (2 / (3 * keq) * np.sqrt(6 / Req)
            * np.log((np.sqrt(1 + Rd) + np.sqrt(Rd + Req)) / (1 + np.sqrt(Req))))


def camb_rd(wm, wb, h=0.6736):   # h de Planck 2018; r_d no depende de h a w fijos
    p = camb.set_params(H0=100 * h, ombh2=wb, omch2=wm - wb - 0.06 / 93.14,
                        mnu=0.06, As=2.1e-9, ns=0.965)
    return camb.get_background(p).get_derived_params()["rdrag"]


def local(f):
    d = np.log(1 + PASO) - np.log(1 - PASO)
    return dict(
        r_d_pivote=float(f(WM0, WB0)),
        exp_wm=float((np.log(f(WM0 * (1 + PASO), WB0)) - np.log(f(WM0 * (1 - PASO), WB0))) / d),
        exp_wb=float((np.log(f(WM0, WB0 * (1 + PASO))) - np.log(f(WM0, WB0 * (1 - PASO)))) / d))


def main():
    out = {"pivote": {"w_m": WM0, "w_b": WB0},
           "paper2": {"r_d_pivote": 147.27, "exp_wm": -0.255, "exp_wb": -0.134},
           "aubourg2015_desi": {"r_d_pivote": 147.05, "exp_wm": -0.23, "exp_wb": -0.13},
           "eh98_local": local(eh98), "camb_local": local(camb_rd)}
    for k in ("paper2", "eh98_local", "camb_local", "aubourg2015_desi"):
        v = out[k]
        print(f"  {k:18s} r_d={v['r_d_pivote']:7.2f}  exp_wm={v['exp_wm']:+.3f}  exp_wb={v['exp_wb']:+.3f}")
    out["lectura"] = ("los exponentes de Paper 2 son las derivadas LOCALES en el pivote "
                      "(EH98 y CAMB coinciden a ~0.005); la normalizacion es la de "
                      "Boltzmann, no la de EH98 (que da ~2.6 % mas). No es 'la formula "
                      "de EH98' sino una expansion en ley de potencias alrededor de Planck.")
    SALIDA.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"  escrito en {SALIDA.relative_to(REPO)}")


if __name__ == "__main__":
    main()
