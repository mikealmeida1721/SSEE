"""
Figures for SSEE Paper 9: Hubble Tension via algebraic screening fraction.
Generates:
  - fig_paper9_h0_tension.pdf  : H0 tension comparison ladder (all measurements)
  - fig_paper9_fscreen_z.pdf   : f_screen(z) vs redshift
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os

OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'figures')
os.makedirs(OUT, exist_ok=True)

# ── SSEE algebraic constants ──────────────────────────────────────────────────
phi = (1 + 5**0.5) / 2
pi  = np.pi
Omega   = phi + pi
AURA    = (3*phi + pi) / 2
MIRA    = AURA / 2
w0      = -AURA / Omega
wa      = -(phi + pi + Omega) / (phi + pi + (phi + pi + Omega))  # ≈ -0.670
H0_alg  = 3 * Omega**2                  # ≈ 67.96 km/s/Mpc  (global background anchor, Postulate D)
H0_global = H0_alg                      # reframe 2026-06-17: H global de fondo = H_alg
H0_MIRA = 67.037                        # km/s/Mpc  (ancla CMB-mapeada VIEJA, superada)
alphaK  = 3*AURA*(pi - phi) / (2*Omega**2)
fscreen = alphaK / (3*MIRA)             # ≈ 0.06725  (MIRA-valor invariante)
# Canonical prediction (reframe): multiplicative screening on the global anchor H_alg
H0_local_phys = H0_global / (1 - fscreen)  # ≈ 72.86  (0.17σ SH0ES) ← titular
# Superseded MIRA-mapped value (old register)
H0_local_old = H0_MIRA / (1 - fscreen)     # ≈ 71.87  (1.12σ SH0ES) ← superado

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: H0 tension ladder
# ─────────────────────────────────────────────────────────────────────────────

# Data: (label, H0, sigma_lo, sigma_hi, color-group)
# group 0 = CMB/early  |  group 1 = distance ladder  |  group 2 = SSEE
measurements = [
    # CMB / early-universe anchors
    ("Planck 2018\n(CMB)",        67.36, 0.54,  0.54,  0),
    ("DESI DR2\n(BAO+Planck)",    67.97, 0.38,  0.38,  0),
    # Distance-ladder anchors
    ("Freedman 2024\n(TRGB+JWST)", 69.96, 1.54, 1.54,  1),
    ("H0LiCOW\n(lensing)",         73.3,  1.7,  1.7,   1),
    ("Masers\n(NGC 4258)",          73.9,  3.0,  3.0,   1),
    ("SH0ES\n(Riess 2022)",         73.04, 1.04, 1.04,  1),
    # SSEE predictions
    ("SSEE IR canonical\n($H_{\\rm alg}$, 0.17$\\sigma$)", H0_local_phys, 0.44, 0.44, 2),
    ("SSEE IR superseded\n($H_{\\rm MIRA}$, 1.12$\\sigma$)", H0_local_old, 0.44, 0.44, 3),
]

colors = {0: '#2166ac', 1: '#d6604d', 2: '#1a9641', 3: '#7fbf7b'}
labels_group = {0: 'CMB / BAO (early)', 1: 'Distance ladder (late)',
                2: 'SSEE prediction (canonical, $H_{\\rm alg}$)', 3: 'SSEE (superseded MIRA)'}

fig, ax = plt.subplots(figsize=(7, 5))

y_pos = list(range(len(measurements)))

for i, (name, h0, slo, shi, grp) in enumerate(measurements):
    c = colors[grp]
    ax.errorbar(h0, i, xerr=[[slo], [shi]], fmt='o', color=c,
                markersize=6, capsize=4, linewidth=1.5, elinewidth=1.5)
    ax.text(h0, i + 0.32, f'{h0:.2f}', ha='center', va='bottom',
            fontsize=7.5, color=c)

ax.set_yticks(y_pos)
ax.set_yticklabels([m[0] for m in measurements], fontsize=8.5)
ax.set_xlabel(r'$H_0$ [km s$^{-1}$ Mpc$^{-1}$]', fontsize=11)
ax.set_title(r'Hubble constant: SSEE predictions vs observations', fontsize=11)
ax.axvline(67.36, color=colors[0], lw=0.8, ls='--', alpha=0.4)
ax.axvline(73.04, color=colors[1], lw=0.8, ls='--', alpha=0.4)
ax.set_xlim(63, 78)
ax.grid(axis='x', lw=0.4, alpha=0.4)

legend_handles = [
    mpatches.Patch(color=colors[0], label=labels_group[0]),
    mpatches.Patch(color=colors[1], label=labels_group[1]),
    mpatches.Patch(color=colors[2], label=labels_group[2]),
    mpatches.Patch(color=colors[3], label=labels_group[3]),
]
ax.legend(handles=legend_handles, loc='lower right', fontsize=8, framealpha=0.9)
ax.invert_yaxis()
fig.tight_layout()
out1 = os.path.join(OUT, 'fig_paper9_h0_tension.pdf')
fig.savefig(out1, bbox_inches='tight')
fig.savefig(out1.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out1}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: f_screen(z) and H0_local(z)
# ─────────────────────────────────────────────────────────────────────────────
# SSEE: w(z) = w0 + wa * z/(1+z) (CPL)
wa_ssee = -(phi + pi + Omega) / (phi + pi + (phi + pi + Omega))
# wa = -P_sc/IGNIS, con P_sc=Omega+phi (PYROS) e IGNIS=pi+P_sc (rama-π).
# El scaffold K_v vale lo mismo (2Ω) pero es otra entidad — guardián R21.
P_sc  = Omega + phi
IGNIS = pi + P_sc
wa_ssee = -P_sc / IGNIS   # ≈ -0.670

z = np.linspace(0, 3, 400)
a = 1 / (1 + z)
w_z = w0 + wa_ssee * z / (1 + z)

# alphaK(z) from Bellini-Sawicki: proportional to (1+w(z)) × rho_DE(z)/rho_tot(z)
# rho_DE / rho_crit = Omega_DE * exp(3 * int_0^z (1+w)/(1+z') dz')
# For CPL: rho_DE(a) = rho_DE0 * a^{-3(1+w0+wa)} * exp(-3 wa (1-a))
Omde0 = AURA / Omega  # ≈ 0.840
Omm0  = 1 - Omde0     # ≈ 0.160

rhoDE_norm = a**(-3*(1+w0+wa_ssee)) * np.exp(-3*wa_ssee*(1-a))
rhotot_norm = Omde0 * rhoDE_norm + Omm0 * a**(-3)
fDE_z = Omde0 * rhoDE_norm / rhotot_norm

# alphaK(z) ≈ alphaK(0) * fDE(z)/fDE(0)  (first-order Bellini-Sawicki scaling)
alphaK_z = alphaK * fDE_z / fDE_z[0]

# f_screen(z) = alphaK(z) / (3*MIRA)
fscreen_z = alphaK_z / (3*MIRA)

# H0_local(z) with multiplicative screening
H0_local_z = H0_alg / (1 - fscreen_z)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

ax1.plot(z, fscreen_z, 'k-', lw=2, label=r'$f_{\rm screen}(z)$')
ax1.axhline(fscreen, color='gray', ls='--', lw=1, alpha=0.7)
ax1.axhline(0, color='k', lw=0.5, alpha=0.3)
ax1.set_ylabel(r'$f_{\rm screen}(z) = \alpha_K(z)\,/\,(3\mathcal{M})$', fontsize=10)
ax1.set_ylim(-0.005, 0.085)
ax1.annotate(fr'$f_{{\rm screen}}(0)={fscreen:.5f}$',
             xy=(0.02, 0.82), xycoords='axes fraction', fontsize=9)
ax1.grid(lw=0.4, alpha=0.4)
ax1.legend(fontsize=9, loc='upper right')

ax2.plot(z, H0_local_z, color='#1a9641', lw=2,
         label=r'$H_0^{\rm local}(z)$ [SSEE mult.]')
ax2.axhline(H0_local_phys, color='#1a9641', ls='--', lw=1, alpha=0.6)
ax2.axhline(H0_alg, color='#2166ac', ls=':', lw=1.2, alpha=0.8,
            label=fr'$H_0^{{\rm alg}}={H0_alg:.2f}$')
ax2.axhline(73.04, color='#d6604d', ls='-.', lw=1.2, alpha=0.8,
            label=r'SH0ES $73.04\pm1.04$')
ax2.fill_between(z, 73.04-1.04, 73.04+1.04, color='#d6604d', alpha=0.10)
ax2.set_ylabel(r'$H_0^{\rm local}$ [km s$^{-1}$ Mpc$^{-1}$]', fontsize=10)
ax2.set_xlabel(r'Redshift $z$', fontsize=11)
ax2.set_ylim(66, 76)
ax2.grid(lw=0.4, alpha=0.4)
ax2.legend(fontsize=8.5, loc='upper right')

fig.suptitle(r'SSEE screening fraction and effective $H_0$ vs redshift', fontsize=11)
fig.tight_layout()
out2 = os.path.join(OUT, 'fig_paper9_fscreen_z.pdf')
fig.savefig(out2, bbox_inches='tight')
fig.savefig(out2.replace('.pdf', '.png'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out2}")

print(f"\nSSEE constants check:")
print(f"  phi     = {phi:.10f}")
print(f"  AURA    = {AURA:.10f}")
print(f"  MIRA    = {MIRA:.10f}")
print(f"  fscreen = {fscreen:.10f}")
print(f"  H0_alg       = {H0_alg:.4f}")
print(f"  H0_MIRA      = {H0_MIRA:.4f}")
print(f"  H0 canonical (H_alg) = {H0_local_phys:.4f}  (0.17 sigma SH0ES)")
print(f"  H0 superseded (MIRA) = {H0_local_old:.4f}  (1.12 sigma SH0ES)")
