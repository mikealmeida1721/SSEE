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
H0   = 67.36
Omm  = 0.3153   # Planck 2018
OmL  = 1 - Omm
h    = H0 / 100
sig8 = 0.811     # Planck 2018

# ── Linear growth factor D(z) — numerical integration ─────────────────────────
def E(z):
    return np.sqrt(Omm*(1+z)**3 + OmL)

def _integrand_growth(zp):
    return (1 + zp) / E(zp)**3

def growth_factor(z, z_norm=0.0):
    """D(z)/D(z_norm) — normalised to 1 at z=z_norm."""
    def I(z_): return quad(_integrand_growth, z_, 1e4)[0]
    return (E(z) * I(z)) / (E(z_norm) * I(z_norm))

# Growth factors at key redshifts
z_vals = [0, 5, 8, 10, 12, 15]
Dz = {z: growth_factor(z) for z in z_vals}
print("\nLinear growth factors D(z)/D(0):")
for z, d in Dz.items():
    print(f"  D(z={z:2d})/D(0) = {d:.5f}")

# ── σ_M(M, z=0) via power-law fit ─────────────────────────────────────────────
# σ(R) is well-approximated by a power law in M over galaxy scales.
# Calibrated at σ(M_8) = σ_8, M_8 = (4π/3)ρ_m(8/h)³.
rho_crit0 = 2.775e11 * h**2     # M☉/Mpc³
rho_m0    = Omm * rho_crit0
M8        = (4*np.pi/3) * rho_m0 * (8/h)**3  # ≈ 2.8e14 M☉
alpha     = 0.30                 # effective slope d ln σ / d ln M^{-1} in 10^10–10^13 range

def sigma_M_z0(M):
    """RMS matter fluctuation in top-hat sphere of mass M, z=0."""
    return sig8 * (M / M8)**(-alpha)

def sigma_Mz(M, z):
    """σ_M at redshift z."""
    return sigma_M_z0(M) * Dz[z]

print(f"\nM_8 = {M8:.3e} M☉  (σ_8 = {sig8}  at R=8/h Mpc)")

# ── Press-Schechter ratio n_SSEE / n_ΛCDM ────────────────────────────────────
def ps_ratio(sigma_M):
    """
    n_SSEE(M,z) / n_ΛCDM(M,z) at fixed σ_M.
    Derived from n(M) ∝ ν f(ν) where ν = δc/σ_M and f(ν)=exp(−ν²/2).
    """
    nu_s = dc_SSEE / sigma_M
    nu_l = dc_LCDM / sigma_M
    return (nu_s / nu_l) * np.exp(-(nu_s**2 - nu_l**2) / 2)

# ── Exceedance probability ratio P(δ > δc) ────────────────────────────────────
from scipy.special import erfc

def exceedance_ratio(sigma_M):
    """Ratio of collapse probabilities."""
    P_ssee = 0.5 * erfc(dc_SSEE / (np.sqrt(2) * sigma_M))
    P_lcdm = 0.5 * erfc(dc_LCDM / (np.sqrt(2) * sigma_M))
    return P_ssee / P_lcdm

# ── Mass table at z=10 ────────────────────────────────────────────────────────
masses = [3e10, 1e11, 3e11, 1e12, 3e12]
z_table = 10

print(f"\n{'='*68}")
print(f"Press-Schechter enhancement at z={z_table}")
print(f"{'M [M☉]':>12}  {'σ_M(z=10)':>10}  {'PS ratio':>10}  {'P-ratio':>10}")
print(f"{'-'*68}")
for M in masses:
    sig = sigma_Mz(M, z_table)
    r   = ps_ratio(sig)
    pe  = exceedance_ratio(sig)
    print(f"{M:12.2e}  {sig:10.4f}  {r:10.3f}  {pe:10.3f}")
print(f"{'='*68}")

# ── Compare with JWST tension ─────────────────────────────────────────────────
print("\nJWST context (Boylan-Kolchin 2023, Nature Astronomy 7, 728):")
print("  ΛCDM deficit at M>10^10.8 M☉, z~10: ~10–100× in number density")
sig_jwst = sigma_Mz(10**10.8, z_table)
r_jwst   = ps_ratio(sig_jwst)
print(f"  M = 10^10.8 M☉: σ_M = {sig_jwst:.3f}, n_SSEE/n_ΛCDM = {r_jwst:.2f}×")
print(f"  → SSEE partially reduces but does not eliminate the tension")

# ── Compact table for Paper 4 (LaTeX-ready) ───────────────────────────────────
print("\nLaTeX table rows:")
tex_masses = [1e11, 3e11, 1e12]
tex_labels = [r"10^{11}", r"3\times10^{11}", r"10^{12}"]
for M, lab in zip(tex_masses, tex_labels):
    sig = sigma_Mz(M, z_table)
    r   = ps_ratio(sig)
    print(f"  ${lab}$ & ${sig:.3f}$ & ${r:.2f}$ \\\\")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: ratio vs sigma_M
sigma_arr = np.linspace(0.15, 1.0, 300)
ratio_arr = [ps_ratio(s) for s in sigma_arr]
exceedance_arr = [exceedance_ratio(s) for s in sigma_arr]

ax = axes[0]
ax.plot(sigma_arr, ratio_arr,     'b-',  lw=2.0, label=r'PS ratio $n_{\rm SSEE}/n_{\Lambda{\rm CDM}}$')
ax.plot(sigma_arr, exceedance_arr,'b--', lw=1.5, label=r'Exceedance ratio $P_{\rm SSEE}/P_{\Lambda{\rm CDM}}$')
ax.axhline(1, color='k', ls=':', lw=0.8)
ax.axvspan(0.3, 0.65, alpha=0.12, color='orange', label=r'JWST regime ($z\sim10$)')

# Mark representative points
for M, lab in [(1e11, r'$10^{11}\,M_\odot$'), (3e11, r'$3\times10^{11}\,M_\odot$'), (1e12, r'$10^{12}\,M_\odot$')]:
    sig = sigma_Mz(M, z_table)
    r   = ps_ratio(sig)
    ax.plot(sig, r, 'o', color='darkorange', ms=5, zorder=5)
    ax.annotate(lab, (sig, r), xytext=(sig+0.04, r+0.05), fontsize=7.5)

ax.set_xlabel(r'$\sigma_M$ (at $z=10$)', fontsize=12)
ax.set_ylabel(r'Halo-count enhancement', fontsize=12)
ax.set_title(r'$\delta_c^{\rm SSEE}=1.6284$ vs $\delta_c^{\Lambda{\rm CDM}}=1.6865$', fontsize=11)
ax.set_ylim(0.9, 5.5)
ax.legend(fontsize=8.5, loc='upper right')
ax.grid(True, alpha=0.3)

# Panel B: ratio vs z at fixed masses
z_range = np.array([5, 6, 7, 8, 9, 10, 12])
D_range = np.array([growth_factor(z) for z in z_range])

ax2 = axes[1]
colors = ['blue', 'green', 'darkorange']
mass_labels = [r'$M=10^{11}\,M_\odot$', r'$M=3\times10^{11}\,M_\odot$', r'$M=10^{12}\,M_\odot$']
for M, col, lab in zip([1e11, 3e11, 1e12], colors, mass_labels):
    sig0 = sigma_M_z0(M)
    ratios_z = [ps_ratio(sig0 * Dz[z]) for z in z_vals[1:]]
    ax2.plot(z_vals[1:], ratios_z, '-o', color=col, ms=4, label=lab)

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
