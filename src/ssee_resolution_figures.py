"""
Resolution figures — visual closure of claims that previously lived only in text/tables.
All values are canonical (VERIFICATION_LEDGER.md §Valores Canónicos); zero fitting.

Figure A: fig_rd_dual.pdf       — sound horizon r_d: naive dynamic-sector 175.6 Mpc
                                   (Omega_m,dyn=0.160) -> algebraic CMB density 147.17 Mpc
                                   (Omega_m,CMB=0.30889, omega_m direct) vs Planck 147.09±0.26 Mpc.
Figure B: fig_s8_resolution.pdf — S8: single-sector challenge 0.846 -> two-sector
                                   resolution 0.759 (0.01 sigma KiDS) vs data bands.

Outputs: results/figures/
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')
os.makedirs(OUT, exist_ok=True)

# ════════════════════════════════════════════════════════════════════════════
# Figure A — dual r_d
# ════════════════════════════════════════════════════════════════════════════
RD_PURE  = 175.6    # Mpc — SSEE naive dynamic-sector drag scale (Omega_m,dyn = 0.160)
RD_CMB   = 147.17   # Mpc — algebraic CMB density (Omega_m,CMB = 0.30889, omega_m direct), CAMB reframe
RD_PLANCK, RD_PLANCK_ERR = 147.09, 0.26   # Mpc, Planck 2018
RD_LCDM  = 147.3    # Mpc — LCDM Eisenstein-Hu at Omega_m = 0.315

fig, ax = plt.subplots(figsize=(8.0, 4.2))

bars = [
    (r'SSEE naive dynamic' + '\n' + r'($\Omega_{m,\rm dyn}=0.160$)', RD_PURE,  '#d6604d'),
    (r'SSEE algebraic CMB' + '\n' + r'($\Omega_{m,\rm CMB}=0.30889$)',  RD_CMB,  '#1a9641'),
    (r'$\Lambda$CDM (EH98,' + '\n' + r'$\Omega_m=0.315$)',            RD_LCDM,  '#4393c3'),
]
ypos = np.arange(len(bars))[::-1]
for y, (lab, val, col) in zip(ypos, bars):
    ax.barh(y, val, color=col, height=0.55, edgecolor='black', lw=0.6, zorder=3)
    ax.text(val + 1.2, y, f'{val:.1f} Mpc', va='center', fontsize=10.5,
            fontweight='bold', color='#222222')

# Planck band
ax.axvspan(RD_PLANCK - RD_PLANCK_ERR, RD_PLANCK + RD_PLANCK_ERR,
           color='#777777', alpha=0.45, zorder=2)
ax.axvline(RD_PLANCK, color='#444444', lw=1.2, ls='--', zorder=2)
ax.text(RD_PLANCK, 2.62, r'Planck 2018: $147.09\pm0.26$ Mpc',
        ha='center', fontsize=9.5, color='#333333')

# Arrow: naive dynamic -> algebraic CMB density (independent prediction, no MIRA factor)
ax.annotate('', xy=(RD_CMB + 2, 1.42), xytext=(RD_PURE - 2, 2.0),
            arrowprops=dict(arrowstyle='-|>', lw=1.8, color='#222222',
                            connectionstyle='arc3,rad=0.25'))
ax.text(163, 1.86, r'$\omega_m$ direct' + '\n' + r'($\omega_b+\omega_c+\omega_\nu$)',
        fontsize=9.5, ha='center', color='#222222', style='italic')

ax.set_yticks(ypos)
ax.set_yticklabels([b[0] for b in bars], fontsize=10)
ax.set_xlabel(r'Drag-epoch sound horizon $r_d$ [Mpc]', fontsize=11)
ax.set_xlim(140, 185)
ax.set_ylim(-0.55, 2.95)
ax.set_title(r'Sound horizon: naive dynamic-sector vs algebraic CMB-density SSEE value vs Planck',
             fontsize=10.5)
ax.grid(axis='x', lw=0.4, alpha=0.4, zorder=0)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
out_a = os.path.join(OUT, 'fig_rd_dual.pdf')
fig.savefig(out_a, bbox_inches='tight')
print(f"Saved: {out_a}")
plt.close(fig)

# ════════════════════════════════════════════════════════════════════════════
# Figure B — S8: challenge -> resolution
# ════════════════════════════════════════════════════════════════════════════
# Data (mean, err)
S8_PLANCK = (0.832, 0.013)
S8_KIDS   = (0.759, 0.024)   # KiDS-1000 (Asgari+2021)
S8_DES    = (0.776, 0.017)
# Model values (canonical)
S8_SINGLE = (0.846, 0.006)   # single-sector ceiling — "the challenge" (3.5 sigma KiDS)
S8_TWOSEC = 0.759            # two-sector canonical (m_phi=41.02) — 0.01 sigma KiDS (resolution)

fig, ax = plt.subplots(figsize=(8.0, 4.4))

entries = [
    (r'SSEE two-sector $\varphi$-DM' + '\n' + r'(canonical, Paper 6)', S8_TWOSEC, 0.0,   '#1a9641'),
    (r'SSEE single-sector' + '\n' + r'(baseline challenge)',           S8_SINGLE[0], S8_SINGLE[1], '#d6604d'),
    (r'DES-Y3 (3$\times$2pt)',                                          S8_DES[0],  S8_DES[1],  '#888888'),
    (r'KiDS-1000',                                                      S8_KIDS[0], S8_KIDS[1], '#888888'),
    (r'Planck 2018 (CMB)',                                              S8_PLANCK[0], S8_PLANCK[1], '#888888'),
]
ypos = np.arange(len(entries))
for y, (lab, val, err, col) in zip(ypos, entries):
    ax.errorbar(val, y, xerr=err if err > 0 else None, fmt='o', ms=9,
                color=col, ecolor=col, elinewidth=2.2, capsize=5,
                markeredgecolor='black', mew=0.7, zorder=4)

# KiDS band as the lensing reference
ax.axvspan(S8_KIDS[0] - S8_KIDS[1], S8_KIDS[0] + S8_KIDS[1],
           color='#b8d8b8', alpha=0.5, zorder=1)
ax.axvline(S8_KIDS[0], color='#1a9641', lw=1.0, ls='--', zorder=2)

# Annotations
ax.annotate(r'$0.01\sigma$ vs KiDS-1000 — resolved',
            xy=(S8_TWOSEC, 0), xytext=(0.789, -0.32),
            fontsize=10, color='#1a9641', fontweight='bold')
ax.annotate(r'$3.5\sigma$ vs KiDS — the challenge',
            xy=(S8_SINGLE[0], 1), xytext=(0.852, 1.25),
            fontsize=10, color='#d6604d')
ax.annotate('', xy=(S8_TWOSEC + 0.004, 0.18), xytext=(S8_SINGLE[0] - 0.004, 0.85),
            arrowprops=dict(arrowstyle='-|>', lw=1.8, color='#222222',
                            connectionstyle='arc3,rad=-0.3'))
ax.text(0.806, 0.52, 'free-streaming\n' + r'$k_{\rm fs}=0.762\,h/{\rm Mpc}$' + '\n'
        + r'($m_\varphi = 41.02$ eV, forward)',
        fontsize=9, ha='center', color='#222222', style='italic')

ax.set_yticks(ypos)
ax.set_yticklabels([e[0] for e in entries], fontsize=10)
ax.set_xlabel(r'$S_8 \equiv \sigma_8\,(\Omega_m/0.3)^{0.5}$', fontsize=11)
ax.set_xlim(0.735, 0.905)
ax.set_ylim(-0.65, 4.6)
ax.set_title(r'$S_8$: single-sector challenge $\rightarrow$ two-sector resolution (zero fitted parameters)',
             fontsize=10.5)
ax.grid(axis='x', lw=0.4, alpha=0.4, zorder=0)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
out_b = os.path.join(OUT, 'fig_s8_resolution.pdf')
fig.savefig(out_b, bbox_inches='tight')
print(f"Saved: {out_b}")
