#!/usr/bin/env python3
"""Momentos (w0, wa) y rho medidos de las CADENAS OFICIALES de DESI DR2.

POR QUE EXISTE (2026-09-19, R65). `ssee_paper2_analysis.py` lleva escritos a
mano los centros, sigmas y rho de los cuatro contornos w0waCDM de DESI DR2, con
la nota «medido de las cadenas oficiales el 2026-07-08». La medida no dejo log:
el numero no tenia de donde venir. Esto la rehace y la deja escrita.

Fuente: Zenodo 10.5281/zenodo.16644577, cosmology_chains/cobaya/base_w_wa
(copia local en el HDD). Receta, la misma que declara el analisis: se descarta
el 30% inicial de cada cadena (burn-in) y se pesa con la columna `weight`.

Se corre solo:  python3 src/p02_mcmc/momentos_desi_dr2_w0wa.py
"""
import json
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
BASE = pathlib.Path("/mnt/datos/SSEE_data/desi_dr2_official/cosmology_chains/"
                    "cobaya/base_w_wa")
CMB = ("planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_"
       "planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing")
COMBOS = {
    "DESI+CMB (DR2 ec.25)":           f"desi-bao-all_{CMB}",
    "DESI+CMB+Pantheon+ (DR2 ec.26)": f"desi-bao-all_pantheonplus_{CMB}",
    "DESI+CMB+Union3 (DR2 ec.27)":    f"desi-bao-all_union3_{CMB}",
    "DESI+CMB+DESY5 (DR2 ec.28)":     f"desi-bao-all_desy5sn_{CMB}",
}
BURN = 0.30
SALIDA = REPO / "results" / "logs" / "desi_dr2_w0wa_momentos.log"   # JSON dentro; .log para que R35 lo vigile


def momentos(carpeta):
    w, x, y = [], [], []
    for f in sorted(carpeta.glob("chain.[0-9]*.txt")):
        col = open(f).readline().lstrip("#").split()
        d = np.loadtxt(f)
        d = d[int(BURN * len(d)):]
        w.append(d[:, col.index("weight")])
        x.append(d[:, col.index("w")])
        y.append(d[:, col.index("wa")])
    w, x, y = map(np.concatenate, (w, x, y))
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    sx = np.sqrt(np.average((x - mx) ** 2, weights=w))
    sy = np.sqrt(np.average((y - my) ** 2, weights=w))
    rho = np.average((x - mx) * (y - my), weights=w) / (sx * sy)
    return dict(w0_bf=mx, sigma_w0=sx, wa_bf=my, sigma_wa=sy, rho=rho,
                n_filas=int(len(w)))


def main():
    out = {"fuente": "Zenodo 10.5281/zenodo.16644577, cosmology_chains/cobaya/base_w_wa",
           "burn_in": BURN, "pesos": "columna weight", "combos": {}}
    for nombre, sub in COMBOS.items():
        m = momentos(BASE / sub)
        out["combos"][nombre] = m
        # tal como los escribe `ssee_paper2_analysis.py` (4 decimales, rho 3)
        out.setdefault("tal_como_se_usa", {})[nombre] = {
            k: f"{v:.3f}" if k == "rho" else f"{v:.4f}" for k, v in m.items() if k != "n_filas"}
        print(f"  {nombre:32s} w0={m['w0_bf']:+.4f}±{m['sigma_w0']:.4f}  "
              f"wa={m['wa_bf']:+.4f}±{m['sigma_wa']:.4f}  rho={m['rho']:+.3f}")
    SALIDA.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"  escrito en {SALIDA.relative_to(REPO)}")


if __name__ == "__main__":
    main()
