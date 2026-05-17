"""
SSEE-V3.6 — Fase 2: Two-Sector phi-DM P(k) con WDM Transfer Function
Aplica la supresión de free-streaming del phi-DM (m=5.60 eV, k_fs=0.493 h/Mpc)
sobre el espectro SSEE+MIRA usando la función de transferencia de Viel et al. 2005
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

# ── Load CLASS outputs ────────────────────────────────────────────────────────
ssee_pk = np.loadtxt(os.path.join(CLASS_DIR, "output/ssee_v36__pk.dat"),         comments="#")
lcdm_pk = np.loadtxt(os.path.join(CLASS_DIR, "output/lcdm_planck2018__pk.dat"),  comments="#")

k_ssee, P_ssee = ssee_pk[:,0], ssee_pk[:,1]
k_lcdm, P_lcdm = lcdm_pk[:,0], lcdm_pk[:,1]

# ── phi-DM Parameters (Paper 6) ───────────────────────────────────────────────
m_phi  = 5.60    # eV
k_fs   = 0.493   # h/Mpc (Dodelson-Widrow derivation)
# Omega split: Omega_CDM = 0.160 / Omega_phiDM = 0.160 / Omega_total = 0.320
f_phiDM = 0.160 / 0.320   # fraction of matter in phi-DM sector = 0.5

# ── WDM Transfer Function (Viel et al. 2005, Eq. 2) ──────────────────────────
# T_wdm(k) = [1 + (alpha * k)^(2*mu)]^(-5/mu)  where mu = 1.12
# alpha calibrated so that T_wdm(k_fs) = 0.5 (half-power point definition)
# From k_fs: alpha = k_fs / (2^(mu/5) - 1)^(1/(2*mu))
mu = 1.12
# Solve: [1 + (alpha*k_fs)^(2mu)]^(-5/mu) = 0.5
# → (alpha*k_fs)^(2mu) = 2^(mu/5) - 1
alpha = k_fs / (2**(mu/5) - 1)**(1/(2*mu))
print(f"WDM transfer function: alpha = {alpha:.4f} h/Mpc (calibrated to k_fs={k_fs})")

def T_wdm(k, alpha=alpha, mu=mu):
    """WDM suppression transfer function (Viel+ 2005)"""
    return (1 + (alpha * k)**(2*mu))**(-5/mu)

def T_wdm_squared(k):
    return T_wdm(k)**2

# ── Two-sector P(k): P_2sector = f_CDM^2 * P_ssee + f_phiDM^2 * T^2 * P_ssee ─
# More precisely: P_2sector(k) = P_ssee(k) * [f_CDM^2 + f_phiDM^2 * T^2(k)]
# But since CDM and phiDM are incoherent (different phases):
# P_2sector(k) ≈ P_ssee(k) * [(1-f_phi) + f_phi * T^2(k)]^2   (coherent approx)
# Or more simply: P_2sector(k) = P_ssee(k) * [1 - f_phi*(1-T^2(k))]  (linear suppression)

# We use the standard approach from mixed dark matter literature:
# P_mix(k) = [f_CDM + f_WDM * T_wdm(k)]^2 * P_CDM(k)
T2 = T_wdm(k_ssee)
P_twosector = P_ssee * (f_phiDM + (1 - f_phiDM) * T2)**2
def sigma8_proxy(k, P):
    """Integral of k^2 P(k) in range 0.05 < k < 5"""
    mask = (k >= 0.05) & (k <= 5.0)
    return np.sqrt(np.trapezoid(P[mask] * k[mask]**2, k[mask]))

P_lcdm_i = np.interp(k_ssee, k_lcdm, P_lcdm)
sig_lcdm   = sigma8_proxy(k_ssee, P_lcdm_i)
sig_ssee   = sigma8_proxy(k_ssee, P_ssee)
sig_two    = sigma8_proxy(k_ssee, P_twosector)

print(f"\n── σ₈ Proxy Results ──────────────────────────────────────────────")
print(f"  ΛCDM (Planck):      σ₈ proxy = {sig_lcdm:.4f}  (normalized reference)")
print(f"  SSEE+MIRA:          σ₈ proxy = {sig_ssee:.4f}  (ratio = {sig_ssee/sig_lcdm:.4f})")
print(f"  SSEE 2-sector:      σ₈ proxy = {sig_two:.4f}   (ratio = {sig_two/sig_lcdm:.4f})")
print(f"\n  Suppression by φ-DM at k_fs={k_fs}: {(1 - T_wdm(k_fs)**2)*100:.1f}% in P(k)")
print(f"  → P(k_fs) ratio vs ΛCDM: {np.interp(k_fs, k_ssee, P_twosector)/np.interp(k_fs, k_lcdm, P_lcdm):.4f}")

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 11))
gs  = gridspec.GridSpec(3, 1, height_ratios=[3, 1.2, 1.0], hspace=0.08)
ax_main = fig.add_subplot(gs[0])
ax_rat  = fig.add_subplot(gs[1], sharex=ax_main)
ax_T    = fig.add_subplot(gs[2], sharex=ax_main)

# Main P(k) panel
ax_main.loglog(k_lcdm, P_lcdm,       color="#58a6ff", lw=2.0, label=r"$\Lambda$CDM (CLASS)", alpha=0.85)
ax_main.loglog(k_ssee, P_ssee,       color="#f78166", lw=1.8, ls="--", label=r"SSEE+MIRA (single sector)", alpha=0.8)
ax_main.loglog(k_ssee, P_twosector,  color="#3fb950", lw=2.4, label=r"SSEE 2-sector ($\varphi$-DM, $k_{fs}=0.493$)", zorder=5)
ax_main.axvline(k_fs, color="#e3b341", lw=1.8, ls=":", alpha=0.9)
ax_main.axvspan(k_fs, k_ssee[-1], alpha=0.06, color="#3fb950")
ax_main.set_ylabel(r"$P(k)$ $[(h^{-1}\mathrm{Mpc})^3]$", fontsize=12)
ax_main.grid(True, which="both", alpha=0.2)
ax_main.legend(loc="lower left", fontsize=9, framealpha=0.3)
ax_main.set_title(
    "SSEE-V3.6: Two-Sector $\\varphi$-DM Matter Power Spectrum (CLASS + WDM Transfer Function)\n"
    r"$m_\varphi = 5.60\,\mathrm{eV}$,  $k_{fs} = 0.493\,h/\mathrm{Mpc}$,  $f_{\varphi\mathrm{DM}} = 0.5$",
    fontsize=12, pad=10)

# Ratio panel
ratio_two  = P_twosector / P_lcdm_i
ratio_ssee = P_ssee / P_lcdm_i
ax_rat.semilogx(k_ssee, ratio_ssee, color="#f78166", lw=1.5, ls="--", alpha=0.8, label="SSEE+MIRA / ΛCDM")
ax_rat.semilogx(k_ssee, ratio_two,  color="#3fb950", lw=2.2, label="SSEE 2-sector / ΛCDM")
ax_rat.axhline(1.0, color="#58a6ff", lw=1.2, ls="--", alpha=0.7)
ax_rat.axvline(k_fs, color="#e3b341", lw=1.5, ls=":", alpha=0.9)
ax_rat.fill_between(k_ssee, ratio_two, 1.0,
                    where=(ratio_two < 1.0), alpha=0.15, color="#3fb950",
                    label=r"φ-DM suppression regime ($k > k_{fs}$)")
ax_rat.set_ylabel(r"$P/P_{\Lambda\mathrm{CDM}}$", fontsize=10)
ax_rat.set_ylim(0, 1.5)
ax_rat.grid(True, which="both", alpha=0.2)
ax_rat.legend(fontsize=8, framealpha=0.3)

# WDM Transfer function panel
k_plot = np.logspace(-3, 0.5, 300)
ax_T.semilogx(k_plot, T_wdm(k_plot)**2, color="#e3b341", lw=2.0,
              label=r"$T^2_{WDM}(k)$ — φ-DM free-streaming suppression")
ax_T.axvline(k_fs, color="#e3b341", lw=1.5, ls=":", alpha=0.9,
             label=rf"$k_{{fs}}=0.493$ (T²=0.5)")
ax_T.axhline(0.5, color="#8b949e", lw=1.0, ls="--", alpha=0.5)
ax_T.set_xlabel(r"Wavenumber $k$ $[h\,\mathrm{Mpc}^{-1}]$", fontsize=12)
ax_T.set_ylabel(r"$T^2_{WDM}(k)$", fontsize=10)
ax_T.set_ylim(-0.05, 1.15)
ax_T.grid(True, which="both", alpha=0.2)
ax_T.legend(fontsize=8, framealpha=0.3)

plt.setp(ax_main.get_xticklabels(), visible=False)
plt.setp(ax_rat.get_xticklabels(),  visible=False)

out_path = os.path.join(CLASS_DIR, "output/ssee_v36_twosector_Pk.png")
plt.savefig(out_path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print(f"\n✅  Plot saved: {out_path}")
