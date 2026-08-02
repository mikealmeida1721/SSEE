"""
Resolution figures — visual closure of claims that previously lived only in text/tables.
All values are canonical (VERIFICATION_LEDGER.md §Valores Canónicos); zero fitting.

Figure A: fig_rd_dual.pdf       — sound horizon r_d: SSEE physical value 147.17 Mpc
                                   (total matter Omega_m,CMB=0.308881, omega_m direct) vs
                                   Planck 147.09±0.26 Mpc and LCDM 147.3 Mpc. The 175.6 Mpc
                                   value is flagged as a CATEGORY ERROR: what one gets by
                                   wrongly inserting the cold sector 1+w0=0.160 into the
                                   background geometry (the DR2 bug, chi2_BAO=726).
Figure B: fig_s8_resolution.pdf — S8 con A_s fijado (0.846, techo) vs A_s libre
                                   ajustado al dato CRUDO (0.7555, 0.11 sigma KiDS).

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
RD_BUG   = 175.6    # Mpc — CATEGORY ERROR: cold sector 1+w0=0.160 wrongly in geometry
RD_CMB   = 147.17   # Mpc — SSEE physical value, total matter (Omega_m,CMB=0.308881, omega_m direct), CAMB
RD_PLANCK, RD_PLANCK_ERR = 147.09, 0.26   # Mpc, Planck 2018
RD_LCDM  = 147.3    # Mpc — LCDM Eisenstein-Hu at Omega_m = 0.315

fig, ax = plt.subplots(figsize=(8.0, 4.2))

bars = [
    (r'Category error' + '\n' + r'($1+w_0=0.160$ in geometry)', RD_BUG,  '#d6604d'),
    (r'SSEE physical' + '\n' + r'($\Omega_{m,\rm CMB}=0.308881$)',  RD_CMB,  '#1a9641'),
    (r'$\Lambda$CDM (EH98,' + '\n' + r'$\Omega_m=0.315$)',            RD_LCDM,  '#4393c3'),
]
ypos = np.arange(len(bars))[::-1]
for y, (lab, val, col) in zip(ypos, bars):
    hatch = '//' if col == '#d6604d' else None
    ax.barh(y, val, color=col, height=0.55, edgecolor='black', lw=0.6, zorder=3, hatch=hatch)
    tag = ' ✗' if col == '#d6604d' else (' ✓' if col == '#1a9641' else '')
    ax.text(val + 1.2, y, f'{val:.1f} Mpc{tag}', va='center', fontsize=10.5,
            fontweight='bold', color='#222222')

# Planck band
ax.axvspan(RD_PLANCK - RD_PLANCK_ERR, RD_PLANCK + RD_PLANCK_ERR,
           color='#777777', alpha=0.45, zorder=2)
ax.axvline(RD_PLANCK, color='#444444', lw=1.2, ls='--', zorder=2)
ax.text(RD_PLANCK, 2.62, r'Planck 2018: $147.09\pm0.26$ Mpc',
        ha='center', fontsize=9.5, color='#333333')

# Arrow: the bug -> the physical value (total matter density enters the geometry)
ax.annotate('', xy=(RD_CMB + 2, 1.42), xytext=(RD_BUG - 2, 2.0),
            arrowprops=dict(arrowstyle='-|>', lw=1.8, color='#222222',
                            connectionstyle='arc3,rad=0.25'))
ax.text(163, 1.86, r'total $\omega_m$ in geometry' + '\n' + r'($\omega_b+\omega_c+\omega_\nu$)',
        fontsize=9.5, ha='center', color='#222222', style='italic')

ax.set_yticks(ypos)
ax.set_yticklabels([b[0] for b in bars], fontsize=10)
ax.set_xlabel(r'Drag-epoch sound horizon $r_d$ [Mpc]', fontsize=11)
ax.set_xlim(140, 185)
ax.set_ylim(-0.55, 2.95)
ax.set_title(r'Sound horizon: SSEE physical value matches Planck; $175.6$ Mpc is the $0.160$-in-geometry bug',
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
# Data (mean, err) — OBSERVACIONES. Ojo: el valor de KiDS es 0.759, NO 0.758.
# El 0.758 era la PREDICCION del sector retirado; tenerlo aqui era el bug H1/H2
# (meter la prediccion en el hueco del dato, forzando 0.00 sigma). Corregido
# 2026-08-01 junto con la retraccion. Ver CANONICAL_VALUES.yaml (obs_KiDS_S8).
S8_PLANCK = (0.832, 0.013)
S8_KIDS   = (0.759, 0.024)   # KiDS-1000 (Asgari+2021) — el DATO
S8_DES    = (0.776, 0.017)
# Model values (canonical 2026-08-01)
S8_FIXED_AS = (0.846, 0.006)  # A_s FIJADO a Planck: techo, no prediccion.
                              # El viejo "3.5 sigma desafio" era artefacto de fijarlo.
S8_MCMC     = (0.7555, 0.0192)  # A_s LIBRE, MCMC vs 225 puntos xi_pm CRUDOS de KiDS.
                                # Un solo sector. 0.11 sigma. log R3_ssee_kids_S8.json

fig, ax = plt.subplots(figsize=(8.0, 4.4))

entries = [
    (r'SSEE, $A_s$ free' + '\n' + r'(MCMC vs raw $\xi_\pm$, 1 sector)', S8_MCMC[0], S8_MCMC[1], '#1a9641'),
    (r'SSEE, $A_s$ fixed to Planck' + '\n' + r'(ceiling, not a prediction)', S8_FIXED_AS[0], S8_FIXED_AS[1], '#d6604d'),
    (r'DES-Y3 (3$\times$2pt)',   S8_DES[0],  S8_DES[1],  '#888888'),
    (r'KiDS-1000',               S8_KIDS[0], S8_KIDS[1], '#888888'),
    (r'Planck 2018 (CMB)',       S8_PLANCK[0], S8_PLANCK[1], '#888888'),
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

ax.annotate(r'$0.11\sigma$ vs KiDS-1000 — no tension',
            xy=(S8_MCMC[0], 0), xytext=(0.782, -0.34),
            fontsize=10, color='#1a9641', fontweight='bold')
ax.annotate(r'the old "$3.5\sigma$ challenge"',
            xy=(S8_FIXED_AS[0], 1), xytext=(0.852, 1.28),
            fontsize=10, color='#d6604d')
ax.annotate('', xy=(S8_MCMC[0] + 0.004, 0.18), xytext=(S8_FIXED_AS[0] - 0.004, 0.85),
            arrowprops=dict(arrowstyle='-|>', lw=1.8, color='#222222',
                            connectionstyle='arc3,rad=-0.3'))
ax.text(0.806, 0.52, 'fit $A_s$ to the RAW data\n'
        + r'(single sector, no particle)',
        fontsize=9, ha='center', color='#222222', style='italic')

ax.set_yticks(ypos)
ax.set_yticklabels([e[0] for e in entries], fontsize=10)
ax.set_xlabel(r'$S_8 \equiv \sigma_8\,(\Omega_m/0.3)^{0.5}$', fontsize=11)
ax.set_xlim(0.720, 0.905)
ax.set_ylim(-0.70, 4.6)
ax.set_title(r'$S_8$: fixing $A_s$ imports the Planck--KiDS tension; fitting it removes it',
             fontsize=10.5)
ax.grid(axis='x', lw=0.4, alpha=0.4, zorder=0)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
out_b = os.path.join(OUT, 'fig_s8_resolution.pdf')
fig.savefig(out_b, bbox_inches='tight')
print(f"Saved: {out_b}")
