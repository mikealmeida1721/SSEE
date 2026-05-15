"""
Figures for SSEE Paper 8: Strong gravity, disformal lensing, Vainshtein.
Generates:
  - fig_paper8_vainshtein.pdf : r_V vs M_obj for solar/stellar/galactic/cluster objects
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
bc   = AURA                 # disformal coupling magnitude
H0   = 67.962e3             # m/s/Mpc  (in SI: 67.96 km/s/Mpc)
Mpc  = 3.0857e22            # m
H0_si = H0 / Mpc            # 1/s ≈ 2.20e-18 s^-1
G    = 6.674e-11            # m^3 kg^-1 s^-2
c    = 3e8                  # m/s
Msun = 1.989e30             # kg

# ── Figure 1: Vainshtein radius r_V vs object mass ───────────────────────────
# r_V = (G M / (3 c^2 H_0^2 / bc^2))^{1/3}
# From Paper 8: r_V^3 = bc^2 G M / (3 c^2 H_0^2)
# (SI units throughout)

M_range = np.logspace(-1, 18, 500)  # in Msun
M_kg = M_range * Msun

r_V_m = (bc**2 * G * M_kg / (3 * c**2 * H0_si**2))**(1/3)
r_V_kpc = r_V_m / (Mpc * 1e-3)   # kpc

fig, ax = plt.subplots(figsize=(7, 5))
ax.loglog(M_range, r_V_kpc, 'k-', lw=2)

# Annotate representative objects
objects = [
    (1,        0.0047,   r'Sun ($1\,M_\odot$)', 'below'),
    (1e6,      0.3,      r'BH ($10^6\,M_\odot$)', 'above'),
    (1e11,     10,       r'Galaxy ($10^{11}\,M_\odot$)', 'above'),
    (1e14,     100,      r'Cluster ($10^{14}\,M_\odot$)', 'above'),
]
for m, r_ref, name, pos in objects:
    r_obj_m = (bc**2 * G * m*Msun / (3*c**2*H0_si**2))**(1/3)
    r_obj_kpc = r_obj_m / (Mpc*1e-3)
    ax.plot(m, r_obj_kpc, 'o', color='#d6604d', ms=8, zorder=5)
    dy = 0.6 if pos == 'above' else -0.6
    ax.annotate(name, xy=(m, r_obj_kpc),
                xytext=(m*1.5, r_obj_kpc * 10**dy),
                fontsize=8.5, ha='left', color='#d6604d',
                arrowprops=dict(arrowstyle='->', color='#d6604d', lw=0.8))

ax.set_xlabel(r'Object mass $M$ [$M_\odot$]', fontsize=11)
ax.set_ylabel(r'Vainshtein radius $r_V$ [kpc]', fontsize=11)
ax.set_title(r'SSEE Vainshtein radius: $r_V = [\beta_c^2 G M\,/\,(3c^2 H_0^2)]^{1/3}$',
             fontsize=10)
ax.grid(which='both', lw=0.4, alpha=0.4)
ax.set_xlim(0.05, 1e18)
fig.tight_layout()
out1 = os.path.join(OUT, 'fig_paper8_vainshtein.pdf')
fig.savefig(out1, bbox_inches='tight')
fig.savefig(out1.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out1}")

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
print(f"  AURA = bc = {bc:.6f}")
print(f"  MIRA    = {MIRA:.6f}")
print(f"  r_V(Sun) = {(bc**2*G*1*Msun/(3*c**2*H0_si**2))**(1/3)/(Mpc*1e-3)*1e3:.4f} pc")
