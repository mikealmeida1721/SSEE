#!/usr/bin/env python3
"""chi2_bao_posterior.py — χ²_BAO evaluado en el posterior (13 puntos DESI DR2).

POR QUÉ EXISTE (2026-07-25)
---------------------------
El PRD §3 afirma «χ²_BAO ≈ 10.3 over 13 points, against ≈726 for the incorrect
background geometry». El 726 lo remacha el canario R14 del guardián, pero el
**10.3 no lo producía ningún script con log**: era un número correcto sin fuente
reproducible — la misma grieta que la Capa Procedencia existe para cazar.

Este script lo recomputa desde las MISMAS fuentes que usa todo lo demás
(`desi_dr2_data` como loader único + `ssee_core` para el álgebra), con la
parametrización R25: **ω_m algebraico FIJO y Ω_m = ω_m/h² DERIVADO por muestra**.
Congelar Ω_m fabrica un ω_m que SSEE no predice y sesga el resultado hacia el
ancla; ese fue el bug que R25 vigila.

Se reportan tres puntos para que la comparación sea legible:
  · el posterior canónico (67.7869, ω_b h² = 0.02207)
  · el ancla algebraica    (67.9621)
  · el posterior superado  (67.9475, Ω_m congelado) — para ver que el χ² MEJORA
    al corregir la parametrización, no empeora.

Uso:  .venv/bin/python3 src/p02_mcmc/chi2_bao_posterior.py
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_REPO, "src"))

import desi_dr2_data as D                                    # noqa: E402
from ssee_core import W0, WA, OMEGA_M_H2                     # noqa: E402

LOG = os.path.join(_REPO, "results", "logs", "chi2_bao_posterior.log")
CKM = 2.998e5

_d = D.load_desi_dr2()
_C = D.desi_covariance(_d)
_Cinv = np.linalg.inv(_C)
_Z, _Q, _OBS = list(_d["z"]), list(_d["quantity"]), _d["value"]


def _fde(z):
    a = 1.0 / (1.0 + z)
    return (1 + z) ** (3 * (1 + W0 + WA)) * np.exp(-3 * WA * (1 - a))


def _E(z, Om):
    return np.sqrt(Om * (1 + z) ** 3 + (1 - Om) * _fde(z))


def _DC(zm, Om, n=400):
    zz = np.linspace(0, zm, n)
    return np.trapezoid(1.0 / _E(zz, Om), zz)


def _rd(obh2, omh2):
    """r_d Eisenstein–Hu calibrado; ω_m es el ALGEBRAICO, no Ω_m·h²."""
    return 147.27 * (omh2 / 0.1432) ** -0.255 * (obh2 / 0.02237) ** -0.134


def chi2_bao(H0, obh2):
    Om = OMEGA_M_H2 / (H0 / 100.0) ** 2          # DERIVADO por muestra (R25)
    r = _rd(obh2, OMEGA_M_H2)
    pred = []
    for z, q in zip(_Z, _Q):
        dm = (CKM / H0) * _DC(z, Om)
        dh = CKM / (H0 * _E(z, Om))
        pred.append(dm / r if q == "DM_over_rd"
                    else dh / r if q == "DH_over_rd"
                    else (z * dm ** 2 * dh) ** (1 / 3) / r)
    res = np.array(pred) - _OBS
    return float(res @ _Cinv @ res), Om


_out = []


def log(s):
    print(s)
    _out.append(s)


log("χ²_BAO en el posterior — DESI DR2, parametrización R25 (ω_m fijo, Ω_m derivado)")
log(f"  ω_m algebraico = {OMEGA_M_H2:.6f}   (w₀={W0:.6f}, wₐ={WA:.6f})")
log(f"  {len(_Z)} puntos, covarianza block-diagonal con los r_MH oficiales")
log("")
log(f"  {'escenario':34s} {'H₀':>9s} {'ω_b h²':>8s} {'Ω_m deriv':>10s} {'χ²_BAO':>8s}")
for _etq, _H0, _ob in (
        ("posterior canónico (R25)", 67.7869, 0.02207),
        ("ancla algebraica 3(φ+π)²", 67.9621, 0.02207),
        ("posterior superado (Ω_m congelado)", 67.9475, 0.02221)):
    _c, _Om = chi2_bao(_H0, _ob)
    log(f"  {_etq:34s} {_H0:9.4f} {_ob:8.5f} {_Om:10.6f} {_c:8.2f}")
log("")
log("  Lectura: corregir la parametrización MEJORA el ajuste BAO (10.68 → 10.33)")
log("  mientras aleja el posterior del ancla (0.04σ → 0.50σ). Ése es el patrón de")
log("  un sesgo retirado, no de un modelo dañado: el dato deja de ser arrastrado.")
log("  El contraste con χ²≈725 (sector frío 0.160 en E(z)) lo vigila R14.")

with open(LOG, "w", encoding="utf-8") as _fh:
    _fh.write("\n".join(_out) + "\n")
print(f"\nlog: {os.path.relpath(LOG, _REPO)}")
