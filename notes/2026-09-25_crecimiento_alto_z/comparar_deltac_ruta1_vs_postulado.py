#!/usr/bin/env python3
"""
¿Que compra el factor n_s? — Press-Schechter con dos δc, mismo modelo.

Caso A — "Ruta 1" (derivado, forma estandar con el fondo SSEE):
    δc_S(z), δc_L(z) del colapso esferico (spherical_collapse_deltac.py).
    A z=10: 1.68647 vs 1.68646 — indistinguibles.
Caso B — "Postulado n_s" (Paper 4, conjetura):
    δc_S = 1.62839 = 1.68647 × n_s, δc_L = 1.68647.

Todo lo demas identico (fondo SSEE, σ8=0.8153/0.811, D(z) con γ):
la UNICA diferencia entre casos es δc. La brecha B−A cuantifica
la "apuesta observable" de la conjetura: si la Ruta 2 se deriva
algun dia, ESTA curva es la huella que debe producir.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── δc: los dos casos ─────────────────────────────────────────────────────
# Ruta 1 (colapso esferico; z_c -> (SSEE, LCDM)); z>=10: dominado por materia
DC_RUTA1 = {0: (1.67634, 1.67599), 1: (1.68496, 1.68438), 2: (1.68613, 1.68580),
            5: (1.68644, 1.68638), 10: (1.68647, 1.68646),
            12: (1.68647, 1.68646), 15: (1.68647, 1.68646)}
DC_POST = (1.68647 * (1 - ((1 + np.sqrt(5)) / 2) ** (-7)), 1.68647)  # (1.62839, 1.68647)

# ── parametros (canonicos; identicos en ambos casos) ───────────────────────
H0_L, Omm_L, sig8_L, gamma_L = 67.36, 0.3153, 0.811, 0.55
H0_S, Omm_S, sig8_S, gamma_S = 67.9621373234, 0.3088808788, 0.8153, 0.5504
OmDE_S, w0_S, wa_S = 1 - Omm_S, -0.839949771, -0.669974886
h_L, h_S = H0_L / 100, H0_S / 100
alpha = 0.30

def D_gamma(z, gamma, Omm, OmDE, w0, wa):
    def E2(zp):
        if OmDE == 0:
            return Omm * (1 + zp) ** 3 + (1 - Omm)
        return (Omm * (1 + zp) ** 3 + OmDE * (1 + zp) ** (3 * (1 + w0 + wa))
                * np.exp(-3 * wa * zp / (1 + zp)))
    def Om_z(zp):
        return Omm * (1 + zp) ** 3 / E2(zp)
    return np.exp(-quad(lambda zp: Om_z(zp) ** gamma / (1 + zp), 0, z)[0])

def M8(Omm, h):
    return (4 * np.pi / 3) * (Omm * 2.775e11 * h ** 2) * (8 / h) ** 3

M8_L, M8_S = M8(Omm_L, h_L), M8(Omm_S, h_S)
Dz = {}
for z in [10, 12, 15]:
    Dz[z] = (D_gamma(z, gamma_L, Omm_L, 0, -1, 0),
             D_gamma(z, gamma_S, Omm_S, OmDE_S, w0_S, wa_S))

def sigma(M, z, which):
    Dz_L, Dz_S = Dz[z]
    if which == "L":
        return sig8_L * (M / M8_L) ** (-alpha) * Dz_L
    return sig8_S * (M / M8_S) ** (-alpha) * Dz_S

def ps_ratio(M, z, dc_S, dc_L):
    s_s, s_l = sigma(M, z, "S"), sigma(M, z, "L")
    nu_s, nu_l = dc_S / s_s, dc_L / s_l
    return (nu_s / nu_l) * np.exp(-(nu_s ** 2 - nu_l ** 2) / 2)

# ── tablas ────────────────────────────────────────────────────────────────
masses = [3e10, 1e11, 3e11, 1e12, 3e12]
print(f"{'z':>3} {'M [Msol]':>12}  {'Caso A (Ruta 1)':>15}  {'Caso B (n_s)':>13}  {'B/A':>8}")
print("-" * 60)
rows = {}
for z in [10, 12, 15]:
    dcA_S, dcA_L = DC_RUTA1[z]
    dcB_S, dcB_L = DC_POST
    for M in masses:
        rA = ps_ratio(M, z, dcA_S, dcA_L)
        rB = ps_ratio(M, z, dcB_S, dcB_L)
        rows[(z, M)] = (rA, rB)
        print(f"{z:>3} {M:12.2e}  {rA:15.4f}  {rB:13.3f}  {rB/rA:8.2f}")

print(f"\nδc Caso A (z=10): SSEE={DC_RUTA1[10][0]:.5f}, LCDM={DC_RUTA1[10][1]:.5f}")
print(f"δc Caso B       : SSEE={DC_POST[0]:.5f}, LCDM={DC_POST[1]:.5f}")

# ── figura ────────────────────────────────────────────────────────────────
Marr = np.logspace(np.log10(3e10), np.log10(3e12), 120)
fig, ax = plt.subplots(figsize=(7.5, 4.8))
for z, ls, col in [(10, "-", "C0"), (12, "--", "C1"), (15, ":", "C2")]:
    dcA_S, dcA_L = DC_RUTA1[z]
    dcB_S, dcB_L = DC_POST
    rA = [ps_ratio(M, z, dcA_S, dcA_L) for M in Marr]
    rB = [ps_ratio(M, z, dcB_S, dcB_L) for M in Marr]
    ax.plot(Marr, rB, color="crimson", ls=ls, lw=2,
            label=f"Caso B (n_s), z={z}" if z == 10 else None)
    ax.plot(Marr, rA, color="steelblue", ls=ls, lw=2,
            label=f"Caso A (Ruta 1), z={z}" if z == 10 else None)
    if z == 10:
        ax.plot(Marr, rB, color="crimson", ls="-", lw=2)
        ax.plot(Marr, rA, color="steelblue", ls="-", lw=2)
ax.axhline(1, color="k", ls=":", lw=0.8)
ax.set_xscale("log")
ax.set_xlabel(r"Masa del halo $[M_\odot]$")
ax.set_ylabel(r"$n_{\rm SSEE}/n_{\Lambda{\rm CDM}}$ (Press–Schechter)")
ax.set_title("Lo que compra el factor $n_s$: mismo modelo, solo cambia $\\delta_c$")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, which="both")
ax.set_ylim(0.95, 3.2)
plt.tight_layout()
out = "results/figures/fig_deltac_comparacion"
plt.savefig(out + ".pdf", bbox_inches="tight")
plt.savefig(out + ".png", dpi=150, bbox_inches="tight")
print(f"\nFigura → {out}.pdf/.png")
