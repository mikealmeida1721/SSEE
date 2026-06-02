#!/usr/bin/env python3
"""
Task 2B: Press-Schechter halo mass function — SSEE vs ΛCDM
δc_SSEE = δc_EdS × n_s = 1.6865 × 0.96556 = 1.6284

Key result: lower δc increases halo abundance exponentially.
Provides quantitative predictions for JWST z>10 galaxy excess.
"""

import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

# ── SSEE algebraic constants ──────────────────────────────────────────────────
PHI    = (1 + np.sqrt(5)) / 2          # golden ratio ≈ 1.6180
n_s    = 1 - PHI**(-7)                 # spectral index  ≈ 0.96556 (SSEE algebraic)
dc_EdS = (3/20) * (12*np.pi)**(2/3)   # EdS collapse threshold ≈ 1.6865
dc_SSEE = dc_EdS * n_s                 # SSEE threshold ≈ 1.6284
dc_LCDM = dc_EdS                       # ΛCDM (EdS limit) ≈ 1.6865

print(f"n_s       = {n_s:.7f}  (1 − φ⁻⁷, Planck: 0.9649)")
print(f"δc(EdS)   = {dc_EdS:.7f}")
print(f"δc(SSEE)  = {dc_SSEE:.7f}")
print(f"δc(ΛCDM)  = {dc_LCDM:.7f}")
print(f"Δδc/δc    = {(dc_LCDM - dc_SSEE)/dc_LCDM*100:.2f}%")

# ── Cosmological parameters ───────────────────────────────────────────────────
# LCDM Parameters (Planck 2018)
H0_L   = 67.36
Omm_L  = 0.3153
OmL_L  = 1 - Omm_L
h_L    = H0_L / 100
sig8_L = 0.811
gamma_L= 0.55

# SSEE Parameters (Paper 1 & Paper 3)
H0_S   = 66.75
Omm_S  = 0.1601
OmDE_S = 0.8399
w0_S   = -0.8399
wa_S   = -0.6699
h_S    = H0_S / 100
sig8_S = 0.792   # IS corrected SSEE sig8 (Paper 3)
gamma_S= 0.657

# ── Linear growth factor D(z) — IS and LCDM integrations ──────────────────────
def D_gamma(z, gamma, Omm, OmDE, w0, wa):
    """Computes growth factor normalized to 1 at z=0 using f = Om(a)^gamma."""
    def E2(zp):
        if OmDE == 0:
            return Omm*(1+zp)**3 + (1-Omm)
        return Omm*(1+zp)**3 + OmDE*(1+zp)**(3*(1+w0+wa))*np.exp(-3*wa*zp/(1+zp))
    
    def Om_z(zp):
        return Omm*(1+zp)**3 / E2(zp)
        
    def integrand(zp):
        return Om_z(zp)**gamma / (1+zp)
        
    val = quad(integrand, 0, z)[0]
    return np.exp(-val)

# Growth factors at key redshifts
z_vals = [0, 5, 6, 7, 8, 9, 10, 12, 15]
Dz_L = {z: D_gamma(z, gamma_L, Omm_L, 0, -1, 0) for z in z_vals}
Dz_S = {z: D_gamma(z, gamma_S, Omm_S, OmDE_S, w0_S, wa_S) for z in z_vals}

print("\nLinear growth factors D(z)/D(0):")
print(f"  z       ΛCDM        SSEE")
for z in [0, 5, 8, 10, 12, 15]:
    print(f"  {z:2d}    {Dz_L[z]:.5f}    {Dz_S[z]:.5f}")

# ── σ_M(M, z=0) via power-law fit ─────────────────────────────────────────────
# We compute M8_L to anchor the power law.
rho_crit0_L = 2.775e11 * h_L**2     # M☉/Mpc³
rho_m0_L    = Omm_L * rho_crit0_L
M8_L        = (4*np.pi/3) * rho_m0_L * (8/h_L)**3  # ≈ 2.8e14 M☉
alpha       = 0.30                 # effective slope d ln σ / d ln M^{-1}

# For a fair comparison at the SAME physical mass M, we anchor both 
# to their respective sigma8 and M8. 
rho_crit0_S = 2.775e11 * h_S**2
rho_m0_S    = Omm_S * rho_crit0_S
M8_S        = (4*np.pi/3) * rho_m0_S * (8/h_S)**3  

def sigma_M_z0_L(M): return sig8_L * (M / M8_L)**(-alpha)
def sigma_M_z0_S(M): return sig8_S * (M / M8_S)**(-alpha)

def sigma_Mz_L(M, z): return sigma_M_z0_L(M) * Dz_L[z]
def sigma_Mz_S(M, z): return sigma_M_z0_S(M) * Dz_S[z]

print(f"\nM_8 (ΛCDM) = {M8_L:.3e} M☉  (σ_8 = {sig8_L})")
print(f"M_8 (SSEE) = {M8_S:.3e} M☉  (σ_8 = {sig8_S})")

# ── Press-Schechter ratio n_SSEE / n_ΛCDM ────────────────────────────────────
def ps_ratio(M, z):
    """n_SSEE(M,z) / n_ΛCDM(M,z) at true SSEE/LCDM variances."""
    sig_s = sigma_Mz_S(M, z)
    sig_l = sigma_Mz_L(M, z)
    nu_s = dc_SSEE / sig_s
    nu_l = dc_LCDM / sig_l
    return (nu_s / nu_l) * np.exp(-(nu_s**2 - nu_l**2) / 2)

# ── Exceedance probability ratio P(δ > δc) ────────────────────────────────────
from scipy.special import erfc

def exceedance_ratio(M, z):
    sig_s = sigma_Mz_S(M, z)
    sig_l = sigma_Mz_L(M, z)
    P_ssee = 0.5 * erfc(dc_SSEE / (np.sqrt(2) * sig_s))
    P_lcdm = 0.5 * erfc(dc_LCDM / (np.sqrt(2) * sig_l))
    return P_ssee / P_lcdm

# ── Mass table at z=10 ────────────────────────────────────────────────────────
masses = [3e10, 1e11, 3e11, 1e12, 3e12]
z_table = 10

print(f"\n{'='*82}")
print(f"Press-Schechter enhancement at z={z_table}")
print(f"{'M [M☉]':>12}  {'σ_M,ΛCDM':>10}  {'σ_M,SSEE':>10}  {'PS ratio':>10}  {'P-ratio':>10}")
print(f"{'-'*82}")
for M in masses:
    sig_l = sigma_Mz_L(M, z_table)
    sig_s = sigma_Mz_S(M, z_table)
    r     = ps_ratio(M, z_table)
    pe    = exceedance_ratio(M, z_table)
    print(f"{M:12.2e}  {sig_l:10.4f}  {sig_s:10.4f}  {r:10.3f}  {pe:10.3f}")
print(f"{'='*82}")

# ── Compare with JWST tension ─────────────────────────────────────────────────
print("\nJWST context (Boylan-Kolchin 2023, Nature Astronomy 7, 728):")
print("  ΛCDM deficit at M>10^10.8 M☉, z~10: ~10–100× in number density")
sig_l_jwst = sigma_Mz_L(10**10.8, z_table)
sig_s_jwst = sigma_Mz_S(10**10.8, z_table)
r_jwst   = ps_ratio(10**10.8, z_table)
print(f"  M = 10^10.8 M☉: σ_M,ΛCDM = {sig_l_jwst:.3f}, σ_M,SSEE = {sig_s_jwst:.3f}")
print(f"  n_SSEE/n_ΛCDM = {r_jwst:.2f}×")

# ── Compact table for Paper 4 (LaTeX-ready) ───────────────────────────────────
print("\nLaTeX table rows:")
tex_masses = [1e11, 3e11, 1e12]
tex_labels = [r"10^{11}", r"3\times10^{11}", r"10^{12}"]
for M, lab in zip(tex_masses, tex_labels):
    sig_l = sigma_Mz_L(M, z_table)
    sig_s = sigma_Mz_S(M, z_table)
    r   = ps_ratio(M, z_table)
    print(f"  ${lab}$ & ${sig_l:.3f}$ & ${sig_s:.3f}$ & ${r:.2f}$ \\\\")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel A: ratio vs sigma_M_LCDM (to show enhancement as a function of the standard variance)
sigma_l_arr = np.linspace(0.15, 1.0, 300)
# To plot against sig_l, we need to map sig_l to sig_s. 
# sig_s / sig_l is independent of M at fixed z.
R_sig = Dz_S[z_table] / Dz_L[z_table] * (sig8_S / sig8_L) * (M8_S / M8_L)**alpha
ratio_arr = []
exceedance_arr = []
for sig_l in sigma_l_arr:
    sig_s = sig_l * R_sig
    nu_s = dc_SSEE / sig_s
    nu_l = dc_LCDM / sig_l
    ratio_arr.append((nu_s / nu_l) * np.exp(-(nu_s**2 - nu_l**2) / 2))
    
    P_ssee = 0.5 * erfc(dc_SSEE / (np.sqrt(2) * sig_s))
    P_lcdm = 0.5 * erfc(dc_LCDM / (np.sqrt(2) * sig_l))
    exceedance_arr.append(P_ssee / P_lcdm)

ax = axes[0]
ax.plot(sigma_l_arr, ratio_arr,     'b-',  lw=2.0, label=r'PS ratio $n_{\rm SSEE}/n_{\Lambda{\rm CDM}}$')
ax.plot(sigma_l_arr, exceedance_arr,'b--', lw=1.5, label=r'Exceedance ratio $P_{\rm SSEE}/P_{\Lambda{\rm CDM}}$')
ax.axhline(1, color='k', ls=':', lw=0.8)
ax.axvspan(0.3, 0.65, alpha=0.12, color='orange', label=r'JWST regime ($z\sim10$)')

# Mark representative points
for M, lab in [(1e11, r'$10^{11}\,M_\odot$'), (3e11, r'$3\times10^{11}\,M_\odot$'), (1e12, r'$10^{12}\,M_\odot$')]:
    sig_l = sigma_Mz_L(M, z_table)
    r   = ps_ratio(M, z_table)
    ax.plot(sig_l, r, 'o', color='darkorange', ms=5, zorder=5)
    ax.annotate(lab, (sig_l, r), xytext=(sig_l+0.02, r+0.05), fontsize=7.5)

ax.set_xlabel(r'$\sigma_{M}^{\Lambda{\rm CDM}}$ (at $z=10$)', fontsize=12)
ax.set_ylabel(r'Halo-count enhancement', fontsize=12)
ax.set_title(r'$\delta_c^{\rm SSEE}=1.6284$ vs $\delta_c^{\Lambda{\rm CDM}}=1.6865$', fontsize=11)
ax.set_ylim(0.9, 5.5)
ax.legend(fontsize=8.5, loc='upper right')
ax.grid(True, alpha=0.3)

# Panel B: ratio vs z at fixed masses
ax2 = axes[1]
colors = ['blue', 'green', 'darkorange']
mass_labels = [r'$M=10^{11}\,M_\odot$', r'$M=3\times10^{11}\,M_\odot$', r'$M=10^{12}\,M_\odot$']
plot_z_vals = [5, 6, 7, 8, 9, 10, 12, 15]

for M, col, lab in zip([1e11, 3e11, 1e12], colors, mass_labels):
    ratios_z = [ps_ratio(M, z) for z in plot_z_vals]
    ax2.plot(plot_z_vals, ratios_z, '-o', color=col, ms=4, label=lab)

ax2.axhline(1, color='k', ls=':', lw=0.8)
ax2.set_xlabel(r'Redshift $z$', fontsize=12)
ax2.set_ylabel(r'$n_{\rm SSEE}/n_{\Lambda{\rm CDM}}$', fontsize=12)
ax2.set_title(r'Enhancement vs redshift', fontsize=11)
ax2.legend(fontsize=8.5)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
outpath = 'results/figures/fig_press_schechter.pdf'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
plt.savefig(outpath.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
print(f"\nFigure → {outpath}")
print("Done.")
