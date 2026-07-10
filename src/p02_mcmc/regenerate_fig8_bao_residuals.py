#!/usr/bin/env python3
"""Regenerador standalone de fig8_bao_residuals (Paper 2).

Residuos normalizados (q_pred - q_obs)/σ para los 13 puntos DESI DR2 (Tabla 4,
2503.14738), tres modelos con su MAP propio. Datos: FUENTE ÚNICA
data/raw/desi_dr2_bao.csv (NO hardcodear; drift DR1/DR2 corregido 2026-07-02).
r_d por modelo vía Eisenstein-Hu (eq. 5 del paper). Regenera la figura que estaba
en DR1 (abril) tras el fix DR2.

Uso: python src/p02_mcmc/regenerate_fig8_bao_residuals.py
"""
import os, sys
import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from desi_dr2_data import load_desi_dr2  # noqa: E402

c = 299792.458
d = load_desi_dr2()                      # z, tracer, quantity, value, sigma
OUT = os.path.join(ROOT, "results", "figures")

def rd_EH(om_h2, ob_h2):                 # eq. 5 del paper
    return 147.27 * (om_h2 / 0.1432) ** -0.255 * (ob_h2 / 0.02237) ** -0.134

def f_de(z, w0, wa):
    return (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))

def E(z, Om, w0, wa):
    de = 1.0 if (w0 == -1 and wa == 0) else f_de(z, w0, wa)
    return np.sqrt(Om * (1 + z) ** 3 + (1 - Om) * de)

# ── MAP de las 3 cadenas (paper2_3models, jul-9, DR2) ──
models = {
    "SSEE":  dict(H0=67.62, Om=0.30889, w0=-0.840,  wa=-0.670,  ob=0.02221, c="#c0392b"),
    r"$\Lambda$CDM": dict(H0=68.27, Om=0.3034, w0=-1.0, wa=0.0, ob=0.02233, c="#2c6fbb"),
    "CPL":   dict(H0=67.26, Om=0.3167, w0=-0.8255, wa=-0.5576, ob=0.02237, c="#27ae60"),
}

def predict(z, quantity, m):
    rd = rd_EH(m["Om"] * (m["H0"] / 100) ** 2, m["ob"])
    DH = c / (m["H0"] * E(z, m["Om"], m["w0"], m["wa"]))
    DM = c * quad(lambda zz: 1 / (m["H0"] * E(zz, m["Om"], m["w0"], m["wa"])), 0, z)[0]
    if quantity.startswith("DV"):
        return (z * DM ** 2 * DH) ** (1 / 3) / rd
    return (DM if quantity.startswith("DM") else DH) / rd

# ── figura ──
plt.rcParams.update({"font.size": 12, "savefig.dpi": 300, "savefig.bbox": "tight"})
fig, ax = plt.subplots(figsize=(9, 5))
labels = [f"{zz:.2f} {q.split('_')[0]}" for zz, q in zip(d["z"], d["quantity"])]
x = np.arange(len(labels))
chi2_by_model = {}
for name, m in models.items():
    pulls = [(predict(zz, q, m) - v) / s
             for zz, q, v, s in zip(d["z"], d["quantity"], d["value"], d["sigma"])]
    chi2 = float(np.sum(np.array(pulls) ** 2))
    chi2_by_model[name] = round(chi2, 1)
    ax.plot(x, pulls, "o-", color=m["c"], lw=1.4, ms=6,
            label=fr"{name} ($\chi^2={chi2:.1f}$)")
for lvl in (-2, -1, 1, 2):
    ax.axhline(lvl, color="0.8", lw=0.8, ls="--", zorder=0)
ax.axhline(0, color="0.4", lw=1.0, zorder=0)
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
ax.set_ylabel(r"$(q_{\rm pred}-q_{\rm obs})/\sigma$")
ax.set_title("DESI DR2 BAO residuals (13 points, 2503.14738 Table 4)")
ax.legend(frameon=False); ax.set_ylim(-3.2, 3.2)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"fig8_bao_residuals.{ext}"))
print("fig8_bao_residuals regenerada (DR2). χ² por modelo:", chi2_by_model)
