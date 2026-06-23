"""
README hero figure: canonical SSEE tension summary (horizontal sigma bars).
All values are canonical (VERIFICATION_LEDGER.md §Valores Canónicos).
Output: results/figures/fig_readme_tensions.png
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'figures')
os.makedirs(OUT, exist_ok=True)

# (label, tension_sigma, reference dataset)
entries = [
    (r'$S_8$ (two-sector $\varphi$-DM)',        0.01, 'KiDS-1000'),
    (r'$w_0$–$w_a$ plane',                       0.05, 'DESI DR2 BAO'),
    (r'$n_s = 1-\varphi^{-7}$',                  0.20, 'Planck 2018'),
    (r'$r_d$ (joint posterior)',                 0.32, 'MCMC multi-probe'),
    (r'$\Omega_b h^2 = (\pi-\varphi)/3\Omega^2$', 0.32, 'Planck 2018'),
    (r'$\Omega_{m,\rm CMB} = \omega_m/h^2 = 0.30889$', 0.88, 'Planck 2018'),
    (r'mean $f\sigma_8$ (6 RSD surveys)',        0.76, r'ties $\Lambda$CDM (0.73$\sigma$)'),
    (r'$H_0^{\rm local}$ = 72.86 km/s/Mpc',      0.17, 'SH0ES'),
]
entries = entries[::-1]  # smallest tension on top

labels = [e[0] for e in entries]
sig    = [e[1] for e in entries]
refs   = [e[2] for e in entries]

fig, ax = plt.subplots(figsize=(8.2, 4.6))
colors = ['#1a9641' if s < 1 else '#fdae61' for s in sig]
bars = ax.barh(range(len(sig)), sig, color=colors, height=0.62,
               edgecolor='black', lw=0.5, zorder=3)

for i, (s, r) in enumerate(zip(sig, refs)):
    ax.text(s + 0.04, i, fr'{s:.2f}$\sigma$  ·  {r}', va='center',
            fontsize=9, color='#333333')

ax.axvline(1.0, color='#999999', ls='--', lw=1.0, zorder=2)
ax.axvline(2.0, color='#d6604d', ls='--', lw=1.0, zorder=2)
ax.text(1.02, -0.62, r'1$\sigma$', color='#777777', fontsize=9, ha='left')
ax.text(2.02, -0.62, r'2$\sigma$', color='#d6604d', fontsize=9, ha='left')

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel(r'Tension vs observation [$\sigma$]', fontsize=11)
ax.set_xlim(0, 2.6)
ax.set_title('SSEE — canonical predictions vs data (algebraic background + forward extensions)',
             fontsize=10.5)
ax.set_ylim(-0.9, len(sig)-0.4)
ax.grid(axis='x', lw=0.4, alpha=0.4, zorder=0)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()

out = os.path.join(OUT, 'fig_readme_tensions.png')
fig.savefig(out, dpi=170, bbox_inches='tight')
print(f"Saved: {out}")
