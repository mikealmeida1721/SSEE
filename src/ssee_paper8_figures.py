"""
Figures for SSEE Paper 8: Strong gravity, disformal lensing, k-mouflage screening.
Generates:
  - fig_paper8_vainshtein.pdf : r_km vs M_obj (k-mouflage, Brax & Valageas 2014)
  - fig_paper8_lensing_ratio.pdf : theta_E^SSEE / theta_E^GR vs k/k_fs
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')
os.makedirs(OUT, exist_ok=True)

# ── SSEE constants ────────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi   = np.pi
AURA = (3*phi + pi) / 2    # ≈ 3.998
MIRA = AURA / 2             # ≈ 1.999
c    = 3e8                  # m/s
Msun = 1.989e30             # kg
Mpc  = 3.0857e22            # m

# K-mouflage constants (natural units: ℏ=c=1)
# M = 8.81 meV = Λ_SSEE  (Paper 10: M^4 = 5φ^8 ρ_crit)
M_eV   = 8.81e-3            # eV
M_Pl_eV = 2.435e27          # eV  (reduced Planck mass)
Msun_eV = Msun * c**2 / 1.602e-19  # kg→eV
hbar_c_eVm = 1.9733e-7      # ℏc in eV·m  (1 eV^-1 = 1.9733e-7 m)

# ── Figure 1: K-mouflage radius r_km vs object mass ──────────────────────────
# r_km^3 = M_obj / (4π M_Pl M²)   [Brax & Valageas 2014, eq. for static source]
# (natural units, then convert to meters via ℏc)

M_range = np.logspace(-1, 18, 500)  # in Msun
M_obj_eV = M_range * Msun_eV

r_km3_nu = M_obj_eV / (4 * pi * M_Pl_eV * M_eV**2)  # eV^-3
r_km_nu  = r_km3_nu**(1/3)                            # eV^-1
r_km_m   = r_km_nu * hbar_c_eVm                       # meters
r_km_kpc = r_km_m / (Mpc * 1e-3)                      # kpc

fig, ax = plt.subplots(figsize=(7, 5))
ax.loglog(M_range, r_km_kpc, 'k-', lw=2)

# Annotate representative objects with their physical radii for comparison
objects = [
    (1,     r'Sun ($1\,M_\odot$)',           'below'),
    (1e6,   r'M31 core ($10^6\,M_\odot$)',   'above'),
    (1e12,  r'Milky Way ($10^{12}\,M_\odot$)', 'above'),
    (1e15,  r'Cluster ($10^{15}\,M_\odot$)', 'above'),
]
for m, name, pos in objects:
    m_ev = m * Msun_eV
    r_nu = (m_ev / (4*pi * M_Pl_eV * M_eV**2))**(1/3)
    r_kpc = r_nu * hbar_c_eVm / (Mpc * 1e-3)
    ax.plot(m, r_kpc, 'o', color='#d6604d', ms=8, zorder=5)
    dy = 0.6 if pos == 'above' else -0.6
    ax.annotate(name, xy=(m, r_kpc),
                xytext=(m*1.5, r_kpc * 10**dy),
                fontsize=8.5, ha='left', color='#d6604d',
                arrowprops=dict(arrowstyle='->', color='#d6604d', lw=0.8))

# Mark 1 kpc reference line
ax.axhline(1.0, color='#2166ac', lw=1.0, ls='--', alpha=0.7, label=r'$1\,\mathrm{kpc}$')
ax.set_xlabel(r'Object mass $M_{\rm obj}$ [$M_\odot$]', fontsize=11)
ax.set_ylabel(r'K-mouflage radius $r_{\rm km}$ [kpc]', fontsize=11)
ax.set_title(r'SSEE k-mouflage: $r_{\rm km}^3 = M_{\rm obj}\,/\,(4\pi M_{\rm Pl} M^2)$,'
             r'$\quad M=8.81\,\mathrm{meV}$', fontsize=10)
ax.legend(fontsize=9)
ax.grid(which='both', lw=0.4, alpha=0.4)
ax.set_xlim(0.05, 1e18)
fig.tight_layout()
out1 = os.path.join(OUT, 'fig_paper8_vainshtein.pdf')
fig.savefig(out1, bbox_inches='tight')
fig.savefig(out1.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out1}")

# Print numerical check
r_sun_m = (Msun_eV / (4*pi * M_Pl_eV * M_eV**2))**(1/3) * hbar_c_eVm
print(f"  r_km(Sun)       = {r_sun_m:.3e} m  (R_sun = 6.96e8 m)")
print(f"  r_km(MW 1e12)   = {r_sun_m * 1e4:.3e} m  (~1 AU = 1.50e11 m)")
print(f"  r_km(cluster 1e15) = {r_sun_m * 1e5:.3e} m")

# ── Figure 2: Lensing ratio theta_E^SSEE / theta_E^GR vs k/k_fs ─────────────
# From Paper 8: the disformal correction amplifies lensing for k < k_fs.
# Simple model: ratio = 1 + MIRA * T_WDM(k)  (qualitative illustration)
# T_WDM(k) = [1 + (alpha_WDM * k)^{2mu}]^{-5/mu}, mu=1.12, alpha=1.6561 h/Mpc
# k_fs = 0.493 h/Mpc (from Paper 6)

k_fs = 0.493  # h/Mpc
alpha_wdm = 1.6561  # h/Mpc  (Paper 6 calibrated)
mu = 1.12
k = np.logspace(-3, 1, 500)  # h/Mpc

T_sq = (1 + (alpha_wdm * k)**(2*mu))**(-5/mu)
# Lensing ratio: sub-k_fs modes feel extra SSEE contribution
# Qualitative: ratio = 1 + (MIRA-1) * T_sq  (MIRA ~ 1.999, so +1 * T_sq)
ratio = 1 + (MIRA - 1) * T_sq

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.semilogx(k, ratio, 'k-', lw=2)
ax.axvline(k_fs, color='#2166ac', lw=1.2, ls='--', label=fr'$k_{{fs}}={k_fs}\ h/\mathrm{{Mpc}}$')
ax.axhline(1.0, color='gray', lw=0.8, ls=':', alpha=0.6, label='GR baseline')
ax.axhline(MIRA, color='#1a9641', lw=0.8, ls='-.', alpha=0.7,
           label=fr'$\mathcal{{M}}={MIRA:.4f}$ (k→0 limit)')
ax.fill_betweenx([0.9, 2.1], 1e-3, k_fs, alpha=0.07, color='#1a9641',
                 label=r'$k < k_{fs}$: $\phi$-DM active')
ax.set_xlabel(r'Wavenumber $k$ [$h$ Mpc$^{-1}$]', fontsize=11)
ax.set_ylabel(r'$\theta_E^{\rm SSEE}/\theta_E^{\rm GR}$', fontsize=11)
ax.set_title(r'Disformal lensing enhancement vs scale', fontsize=10)
ax.set_ylim(0.9, 2.15)
ax.legend(fontsize=8.5, loc='lower left')
ax.grid(which='both', lw=0.4, alpha=0.4)
fig.tight_layout()
out2 = os.path.join(OUT, 'fig_paper8_lensing_ratio.pdf')
fig.savefig(out2, bbox_inches='tight')
fig.savefig(out2.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out2}")

print(f"\nPaper 8 constants check:")
print(f"  AURA = {AURA:.6f}")
print(f"  MIRA = {MIRA:.6f}")
