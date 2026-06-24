"""
SSEE-V3.6 — Fase 1B: Matter Power Spectrum P(k) Comparison
Compares SSEE vs ΛCDM P(k) and identifies the phi-DM free-streaming signature.

⚠️ DEPRECADO / SCRATCH pre-reframe (2026-06-19). Este plot de validación Fase-1B no
   alimenta ningún paper y usa valores retirados (k_fs=0.659 de m_φ=36.95). Canónico:
   k_fs=0.762 h/Mpc, m_φ=41.02 eV (SOLAR²·KRYSTOS). Las figuras φ-DM vigentes salen
   del pipeline forward (ssee_paper6_canonical_particle.py + plot_calibrated_phiDM.py,
   pendiente de re-correr forward). No usar para resultados.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

plt.rcParams.update({
    "figure.facecolor": "#0d1117", "axes.facecolor":   "#161b22",
    "axes.edgecolor":   "#30363d", "axes.labelcolor":  "#e6edf3",
    "xtick.color":      "#8b949e", "ytick.color":      "#8b949e",
    "text.color":       "#e6edf3", "grid.color":       "#21262d",
    "grid.linewidth":   0.6,       "font.family":      "DejaVu Sans",
    "font.size":        11,
})

CLASS_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Load P(k) data ────────────────────────────────────────────────────────────
ssee_pk  = np.loadtxt(os.path.join(CLASS_DIR, "output/ssee_v36__pk.dat"),          comments="#")
lcdm_pk  = np.loadtxt(os.path.join(CLASS_DIR, "output/lcdm_planck2018__pk.dat"),   comments="#")
nomira_pk = np.loadtxt(os.path.join(CLASS_DIR, "output/ssee_v36_nomira__pk.dat"),  comments="#")

k_ssee,  P_ssee  = ssee_pk[:,0],  ssee_pk[:,1]
k_lcdm,  P_lcdm  = lcdm_pk[:,0],  lcdm_pk[:,1]
k_nm,    P_nm    = nomira_pk[:,0], nomira_pk[:,1]

# ── Key scales ────────────────────────────────────────────────────────────────
k_fs   = 0.659   # phi-DM free-streaming scale (Paper 6)
k_eq   = 0.016   # matter-radiation equality
k_bao  = 0.1     # BAO scale

# ── Interpolate ΛCDM onto SSEE k-grid for ratio ──────────────────────────────
P_lcdm_interp = np.interp(k_ssee, k_lcdm, P_lcdm)
ratio          = P_ssee / P_lcdm_interp

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1.2], hspace=0.08)

ax_main = fig.add_subplot(gs[0])
ax_rat  = fig.add_subplot(gs[1], sharex=ax_main)

# Main panel — P(k)
ax_main.loglog(k_lcdm, P_lcdm,  color="#58a6ff", lw=2.0, label=r"$\Lambda$CDM (Planck 2018, CLASS)", alpha=0.85)
ax_main.loglog(k_ssee, P_ssee,  color="#f78166", lw=2.4, label=r"SSEE+MIRA ($\Omega_m=0.320$)", zorder=5)
ax_main.loglog(k_nm,   P_nm,    color="#e3b341", lw=1.4, ls="--", label=r"SSEE-noMIRA ($\Omega_m=0.160$)", alpha=0.8)

# Vertical lines for key scales
ax_main.axvline(k_fs,  color="#3fb950", lw=1.8, ls=":", alpha=0.9, label=rf"$k_{{fs}}=0.659\,h/\mathrm{{Mpc}}$ (φ-DM free-streaming)")
ax_main.axvline(k_bao, color="#d2a8ff", lw=1.2, ls="--", alpha=0.6, label=r"BAO scale $\sim 0.1\,h/\mathrm{Mpc}$")
ax_main.axvline(k_eq,  color="#8b949e", lw=1.0, ls="--", alpha=0.5, label=r"$k_{eq} \sim 0.016\,h/\mathrm{Mpc}$")

# Shaded region: modes suppressed by phi-DM free-streaming
ax_main.axvspan(k_fs, k_ssee[-1], alpha=0.07, color="#3fb950",
                label=r"Suppressed regime $k>k_{fs}$ (φ-DM)")

ax_main.set_ylabel(r"$P(k)$ $[(h^{-1}\mathrm{Mpc})^3]$", fontsize=12)
ax_main.set_xlim(k_ssee[0], k_ssee[-1])
ax_main.grid(True, which="both", alpha=0.2)
ax_main.legend(loc="lower left", fontsize=8.5, framealpha=0.3, ncol=2)
ax_main.set_title(
    "SSEE-V3.6 Matter Power Spectrum — CLASS Boltzmann (Fase 1B)\n"
    r"Identifying the $\varphi$-DM free-streaming signature at $k_{fs}=0.659\,h/\mathrm{Mpc}$",
    fontsize=12, pad=10)

# Ratio panel
ax_rat.semilogx(k_ssee, ratio, color="#f78166", lw=2.0, label="SSEE+MIRA / ΛCDM")
ax_rat.axhline(1.0, color="#58a6ff", lw=1.2, ls="--", alpha=0.7)
ax_rat.axvline(k_fs, color="#3fb950", lw=1.8, ls=":", alpha=0.9)
ax_rat.fill_between(k_ssee, ratio, 1.0,
                    where=(ratio < 1.0), alpha=0.15, color="#f78166",
                    label="SSEE suppression vs ΛCDM")
ax_rat.fill_between(k_ssee, ratio, 1.0,
                    where=(ratio > 1.0), alpha=0.15, color="#3fb950",
                    label="SSEE excess vs ΛCDM")
ax_rat.set_xlabel(r"Wavenumber $k$ $[h\,\mathrm{Mpc}^{-1}]$", fontsize=12)
ax_rat.set_ylabel(r"$P_\mathrm{SSEE}/P_{\Lambda\mathrm{CDM}}$", fontsize=10)
ax_rat.set_ylim(0.0, 2.0)
ax_rat.grid(True, which="both", alpha=0.2)
ax_rat.legend(fontsize=8, framealpha=0.3)

plt.setp(ax_main.get_xticklabels(), visible=False)

out_path = os.path.join(CLASS_DIR, "output/ssee_v36_Pk_comparison.png")
plt.savefig(out_path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print(f"✅  Plot saved: {out_path}")

# ── Numerical report ──────────────────────────────────────────────────────────
print("\n── P(k) Numerical Report ─────────────────────────────────────────")
# Ratio at key scales
for scale, label in [(0.01, "k=0.01 (large scale)"),
                     (0.1,  "k=0.10 (BAO scale)"),
                     (0.659,"k=0.659 (k_fs)"),
                     (1.0,  "k=1.00 (non-linear)"),
                     (5.0,  "k=5.00 (deep non-linear)")]:
    if scale <= k_ssee[-1]:
        r = np.interp(scale, k_ssee, ratio)
        print(f"  SSEE/ΛCDM at {label}: {r:.4f}  ({(r-1)*100:+.1f}%)")

# Suppression at k_fs
supp = np.interp(k_fs, k_ssee, ratio)
print(f"\n  ➤ Suppression at k_fs=0.659: {(1-supp)*100:.1f}% below ΛCDM")
print(f"    This is the expected φ-DM free-streaming signature.")
print("──────────────────────────────────────────────────────────────────")
