#!/usr/bin/env python3
"""Figura: donde queda el borde y cuanto suprime, contra la masa del sector phi."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M = np.array([40.702, 100.0, 300.0])
K50 = np.array([0.339, 0.901, 2.973])
DS8 = np.array([10.42, 1.56, 0.11])            # supresion de sigma8 en %
SIG_KIDS = 100 * 0.024 / 0.759                 # 3.16% : error de KiDS en S8

# leyes de potencia ajustadas a los 3 puntos medidos
bk = np.polyfit(np.log(M), np.log(K50), 1)
bd = np.polyfit(np.log(M), np.log(DS8), 1)
mm = np.logspace(np.log10(30), np.log10(3000), 200)
k_fit = np.exp(np.polyval(bk, np.log(mm)))
d_fit = np.exp(np.polyval(bd, np.log(mm)))
m_min = np.exp((np.log(SIG_KIDS) - bd[1]) / bd[0])
k_min = np.exp(np.polyval(bk, np.log(m_min)))

fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.0))
C1, C2, CX = '#1f5f8b', '#a8452b', '#4a4a4a'

a = ax[0]
a.loglog(mm, d_fit, color=C1, lw=2, zorder=3)
a.loglog(M, DS8, 'o', color=C1, ms=9, mec='white', mew=1.6, zorder=5,
         label='medido con CLASS')
a.axhspan(0, SIG_KIDS, color='#8fbf6f', alpha=.28, zorder=1)
a.axhline(SIG_KIDS, color='#4f7f2f', lw=1.4, ls='--', zorder=2)
a.text(34, SIG_KIDS * 0.55, f'incertidumbre de KiDS ({SIG_KIDS:.2f}%)\n'
       'debajo de aqui, invisible', fontsize=9, color='#3f6f1f', va='top')
a.axvline(40.702, color=CX, lw=1.4, ls=':', zorder=2)
a.text(42, 0.06, 'SOLAR²·KRYSTOS\n40,70 eV', fontsize=9, color=CX)
a.axvline(m_min, color=C2, lw=1.8, zorder=4)
a.text(m_min * 1.1, 0.06, f'cota inferior\n{m_min:.0f} eV', fontsize=9,
       color=C2, fontweight='bold')
a.set_xlabel('masa del sector φ  [eV]')
a.set_ylabel('supresión de σ₈  [%]')
a.set_title('Cuánto borra estructura', fontsize=11, loc='left')
a.set_ylim(0.05, 40); a.set_xlim(30, 3000)
a.grid(alpha=.22, which='both'); a.legend(fontsize=9, loc='lower left')

a = ax[1]
a.loglog(mm, k_fit, color=C1, lw=2, zorder=3)
a.loglog(M, K50, 'o', color=C1, ms=9, mec='white', mew=1.6, zorder=5)
a.axhspan(0.05, 3.0, color='#c9b8d8', alpha=.30, zorder=1)
a.text(34, 1.3, 'rango que MIDE la cizalla\n(aquí el dato ya dijo: sin supresión)',
       fontsize=9, color='#5a3f7a')
a.axvline(40.702, color=CX, lw=1.4, ls=':', zorder=2)
a.axvline(m_min, color=C2, lw=1.8, zorder=4)
a.plot([m_min], [k_min], 's', color=C2, ms=10, mec='white', mew=1.6, zorder=6)
a.text(m_min * 1.1, k_min * 1.25, f'buscar el borde\nen k ≈ {k_min:.1f} h/Mpc',
       fontsize=9, color=C2, fontweight='bold')
a.set_xlabel('masa del sector φ  [eV]')
a.set_ylabel('k del borde de free-streaming  [h/Mpc]')
a.set_title('Dónde queda el borde', fontsize=11, loc='left')
a.set_xlim(30, 3000); a.set_ylim(0.15, 40)
a.grid(alpha=.22, which='both')

fig.suptitle('Sector φ: la masa fija dónde muerde el borde — y la cizalla pone '
             'una cota inferior', fontsize=12.5, x=0.055, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig('fig_cota_masa.png', dpi=140)
print(f'  k_50   = {np.exp(bk[1]):.5f} * m^{bk[0]:.3f}')
print(f'  ds8(%) = {np.exp(bd[1]):.4f} * m^{bd[0]:.3f}')
print(f'  cota inferior de la masa      = {m_min:.1f} eV')
print(f'  factor vs 40.70 eV            = {m_min/40.702:.2f}x')
print(f'  multiplicador que exigiria    = {m_min/0.06849:.0f}   (canonico 594.28)')
print(f'  k del borde en la cota        = {k_min:.2f} h/Mpc')
