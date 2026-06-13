"""
ssee_candidate1_test.py
=======================
Test del Candidato 1: ¿puede IS viscosity sola (sin WDM) resolver fσ₈?

Descompone la mejora de Paper 6 en sus dos contribuciones:
  (A) Restauración de f(z): Ω_m,growth 0.160→0.320 (dos-sectores, independiente de WDM)
  (B) Supresión σ₈:         WDM T(k) → σ₈_eff 0.702→0.737

Candidato 1 preserva (A), elimina (B), y usa σ₈_IS = 0.702 (Paper 5).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', 'results', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

# ── Constantes SSEE ──────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi_  = np.pi
beta = (pi_ + phi) / 2
Kv   = phi + pi_ + (pi_ + phi)
Tr   = 3 * (phi + beta)
Mv   = phi + pi_ + Kv
KAL0 = beta + pi_
w0   = -Tr / Mv
wa   = -(pi_ + phi + phi) / Kv
MIRA = (3*phi + pi_) / 4

Om_dyn   = 1.0 + w0        # 0.160 — sector CDM dinámico
Om_DE    = 1.0 - Om_dyn    # 0.840
Om_phiDM = (MIRA-1)*Om_dyn # 0.160 — sector φ-DM
Om_both  = Om_dyn + Om_phiDM  # 0.320 — ambos sectores (grandes escalas)

# σ₈ valores clave
sig8_Planck = 0.811   # ΛCDM referencia
sig8_IS     = 0.702   # IS Paper 5 (G = D₁_SSEE/D₁_ΛCDM = 0.866, sin WDM)
sig8_WDM    = 0.737   # Paper 6 (IS + WDM T(k))
sig8_LCDM   = 0.811   # Planck 2018

# ── Ecuaciones de crecimiento ─────────────────────────────────────────────────
def f_DE_raw(a): return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
_n = f_DE_raw(1.0)
def f_DE(a): return f_DE_raw(a) / _n

def make_growth_rhs(Om_m):
    """Factory: crea RHS de la ecuación D'' con Ω_m dado."""
    Om_rest = 1.0 - Om_m  # Ω_DE efectivo (simplificación background)
    def rhs(lna, Y):
        a    = np.exp(lna)
        e2   = max(Om_m * a**(-3) + Om_DE * f_DE(a), 1e-30)
        dlne2= (-3*Om_m*a**(-3)
                + Om_DE*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a) / e2
        h    = 0.5 * dlne2
        Omm_a= Om_m * a**(-3) / e2
        D, Dp= Y
        return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]
    return rhs

def growth_lcdm(lna, Y, Om_lcdm=0.3153):
    a    = np.exp(lna)
    e2   = max(Om_lcdm*a**(-3) + (1-Om_lcdm), 1e-30)
    h    = 0.5*(-3*Om_lcdm*a**(-3))/e2
    Omm_a= Om_lcdm*a**(-3)/e2
    D, Dp= Y
    return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]

lna_pts = np.linspace(np.log(0.01), 0.0, 600)
IC      = [0.01, 0.01]

def integrate(rhs):
    sol = solve_ivp(rhs, [lna_pts[0], 0.0], IC,
                    t_eval=lna_pts, rtol=1e-10, atol=1e-13)
    D  = sol.y[0] / sol.y[0][-1]
    Dp = sol.y[1] / sol.y[0][-1]
    return D, Dp

# Cuatro casos de crecimiento:
D_CDM,   Dp_CDM   = integrate(make_growth_rhs(Om_dyn))   # Ω_m = 0.160 (P5)
D_both,  Dp_both  = integrate(make_growth_rhs(Om_both))  # Ω_m = 0.320 (P6 RSD)
D_lcdm,  Dp_lcdm  = integrate(growth_lcdm)               # ΛCDM

# ── fσ₈ por survey ───────────────────────────────────────────────────────────
Z_RSD     = np.array([0.067, 0.150, 0.320, 0.380, 0.510, 0.570])
surveys   = ['6dFGRS', 'SDSS MGS', 'BOSS z=0.32', 'BOSS z=0.38', 'eBOSS DR16', 'WiggleZ']
fsig8_obs = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
fsig8_err = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])

def get_fsig8(z, D, Dp, sig8):
    lna_t = np.log(1/(1+z))
    f_t   = float(np.interp(lna_t, lna_pts, Dp/D))
    D_t   = float(np.interp(lna_t, lna_pts, D))
    return f_t * sig8 * D_t

# ── Los 4 modelos a comparar ──────────────────────────────────────────────────
models = {
    'P5: IS solo\n(Ω_m=0.160, σ₈=0.702)': {
        'D': D_CDM, 'Dp': Dp_CDM, 'sig8': sig8_IS,
        'color': 'tomato', 'ls': '--',
        'label': rf'P5 IS solo ($\Omega_m=0.160$, $\sigma_8=0.702$)'
    },
    'Cand.1: IS+2sector\n(Ω_m=0.320, σ₈=0.702)': {
        'D': D_both, 'Dp': Dp_both, 'sig8': sig8_IS,
        'color': 'steelblue', 'ls': '-',
        'label': rf'Cand.1: IS+2sector ($\Omega_m=0.320$, $\sigma_8=0.702$, sin WDM)'
    },
    'P6: IS+2sector+WDM\n(Ω_m=0.320, σ₈=0.737)': {
        'D': D_both, 'Dp': Dp_both, 'sig8': sig8_WDM,
        'color': 'darkgreen', 'ls': '-.',
        'label': rf'P6 IS+2sector+WDM ($\Omega_m=0.320$, $\sigma_8=0.737$)'
    },
    'ΛCDM\n(Ω_m=0.315, σ₈=0.811)': {
        'D': D_lcdm, 'Dp': Dp_lcdm, 'sig8': sig8_LCDM,
        'color': 'gray', 'ls': ':',
        'label': r'$\Lambda$CDM ($\Omega_m=0.315$, $\sigma_8=0.811$)'
    },
}

results = {}
for name, cfg in models.items():
    fsig8 = np.array([get_fsig8(z, cfg['D'], cfg['Dp'], cfg['sig8']) for z in Z_RSD])
    tensions = (fsig8_obs - fsig8) / fsig8_err
    results[name] = {'fsig8': fsig8, 'tensions': tensions,
                     'mean_abs_tension': np.mean(np.abs(tensions))}

# ── PRINT ────────────────────────────────────────────────────────────────────
SEP = "=" * 75
sep = "─" * 75

print(SEP)
print("  SSEE — Test Candidato 1: IS sin WDM vs Paper 6")
print("  Pregunta: ¿la restauración Ω_m,growth sola resuelve fσ₈?")
print(SEP)

header = f"  {'Survey':14s} {'z':>5}  {'Obs':>6}  "
for name in models:
    short = name.split('\n')[0]
    header += f"{'fs8':>6} {'σ':>6}  "
print(header)
print(sep)

for i, (z, sv) in enumerate(zip(Z_RSD, surveys)):
    row = f"  {sv:14s} {z:>5.3f}  {fsig8_obs[i]:>6.3f}  "
    for name, r in results.items():
        row += f"{r['fsig8'][i]:>6.3f} {r['tensions'][i]:>+6.2f}σ  "
    print(row)

print(sep)
print(f"\n  Tensión media |σ|:")
for name, r in results.items():
    short = name.split('\n')[0].replace('\n','')
    print(f"    {short:45s}: {r['mean_abs_tension']:.4f}σ")

# ── Diagnóstico central ───────────────────────────────────────────────────────
t_P5   = results['P5: IS solo\n(Ω_m=0.160, σ₈=0.702)']['mean_abs_tension']
t_C1   = results['Cand.1: IS+2sector\n(Ω_m=0.320, σ₈=0.702)']['mean_abs_tension']
t_P6   = results['P6: IS+2sector+WDM\n(Ω_m=0.320, σ₈=0.737)']['mean_abs_tension']
t_LCDM = results['ΛCDM\n(Ω_m=0.315, σ₈=0.811)']['mean_abs_tension']

delta_C1_P6 = t_C1 - t_P6  # cuánto peor es C1 vs P6

print(f"\n{SEP}")
print("  DIAGNÓSTICO")
print(SEP)
print(f"\n  Mejora de IS+2sector sin WDM vs Paper 5 puro:")
print(f"    {t_P5:.3f}σ → {t_C1:.3f}σ  (Δ = {t_P5-t_C1:.3f}σ mejor)")
print(f"\n  Contribución específica del WDM (C1 vs P6):")
print(f"    C1 (sin WDM): {t_C1:.4f}σ")
print(f"    P6 (con WDM): {t_P6:.4f}σ")
print(f"    Diferencia:   {delta_C1_P6:+.4f}σ  ({'+' if delta_C1_P6>0 else ''}{'peor' if delta_C1_P6>0 else 'mejor'} sin WDM)")
print(f"\n  ΛCDM referencia: {t_LCDM:.4f}σ")

if t_C1 < 1.5:
    verdict = f"VIABLE — Candidato 1 resuelve fσ₈ a {t_C1:.2f}σ sin WDM"
elif t_C1 < 2.5:
    verdict = f"PARCIAL — C1 mejora tensión a {t_C1:.2f}σ pero no la resuelve completamente"
else:
    verdict = f"INSUFICIENTE — C1 da {t_C1:.2f}σ, WDM es necesario"

print(f"\n  VEREDICTO: {verdict}")

# Qué porcentaje de la mejora viene de la restauración Ω_m vs WDM
mejora_total = t_P5 - t_P6   # mejora total P6 sobre P5
mejora_Omm   = t_P5 - t_C1   # mejora solo de Ω_m restaurado
mejora_WDM   = t_C1 - t_P6   # mejora adicional del WDM

if mejora_total > 0:
    pct_Omm = mejora_Omm / mejora_total * 100
    pct_WDM = mejora_WDM / mejora_total * 100
    print(f"\n  Descomposición de la mejora P5→P6:")
    print(f"    Restauración Ω_m,growth (0.160→0.320): {mejora_Omm:.3f}σ ({pct_Omm:.1f}% del total)")
    print(f"    WDM σ₈ (0.702→0.737):                  {mejora_WDM:.3f}σ ({pct_WDM:.1f}% del total)")
    print(f"    Mejora total P5→P6:                    {mejora_total:.3f}σ")

print(SEP)

# ── FIGURA ────────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

z_s = np.linspace(0.04, 0.65, 200)
for name, cfg in models.items():
    fs8_s = [get_fsig8(z, cfg['D'], cfg['Dp'], cfg['sig8']) for z in z_s]
    short = cfg['label']
    t_val = results[name]['mean_abs_tension']
    ax1.plot(z_s, fs8_s, color=cfg['color'], ls=cfg['ls'], lw=2.2,
             label=f"{short} | {t_val:.2f}σ")

ax1.errorbar(Z_RSD, fsig8_obs, yerr=fsig8_err, fmt='ko', ms=7,
             capsize=4, zorder=6, label='Observaciones RSD')
ax1.set_xlabel(r'Redshift $z$', fontsize=11)
ax1.set_ylabel(r'$f\sigma_8(z)$', fontsize=11)
ax1.set_title(r'$f\sigma_8$: IS solo vs IS+WDM', fontsize=11, fontweight='bold')
ax1.legend(fontsize=7.5, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0.04, 0.65)

# Panel 2: tensiones por survey
x = np.arange(len(surveys))
width = 0.2
for idx, (name, r) in enumerate(results.items()):
    short = name.split('\n')[0]
    offset = (idx - 1.5) * width
    ax2.bar(x + offset, np.abs(r['tensions']), width, alpha=0.85,
            color=list(models.values())[idx]['color'],
            edgecolor='black', lw=0.8, label=f"{short} ({r['mean_abs_tension']:.2f}σ)")

ax2.axhline(1.0, color='green',  ls=':', lw=1.8, label='1σ')
ax2.axhline(2.0, color='orange', ls='--', lw=1.8, label='2σ')
ax2.set_xticks(x)
ax2.set_xticklabels(surveys, rotation=25, ha='right', fontsize=8)
ax2.set_ylabel(r'$|f\sigma_8|$ tension ($\sigma$)', fontsize=11)
ax2.set_title('Tensión por survey', fontsize=11, fontweight='bold')
ax2.legend(fontsize=7.5)
ax2.grid(True, alpha=0.3, axis='y')

fig.suptitle(
    'Candidato 1: IS+2sector sin WDM — ¿Resuelve fσ₈?\n'
    r'Clave: ¿cuánto aporta $\Omega_{m,growth}$ vs $\sigma_8^{WDM}$?',
    fontsize=11, fontweight='bold'
)
fig.tight_layout()
out = os.path.join(OUTDIR, 'fig_candidate1_fsig8.pdf')
fig.savefig(out, bbox_inches='tight')
fig.savefig(out.replace('.pdf','.png'), bbox_inches='tight', dpi=150)
print(f"\n  Figura guardada: {out}")
