#!/usr/bin/env python3
"""
ssee_paper6_alphaM_scan.py  —  SSEE Paper 6: αB + αM 2D Exploration
====================================================================
Explora el espacio (αB0, αM0) del EFT-DE para encontrar qué combinación
resuelve la tensión fσ8 dentro de los límites observacionales:
  - |αM(z=0)| < 0.02  (lunar laser ranging / equivalence principle)
  - αT = 0            (GW170817)

Física:
  Con αM ≠ 0, la masa de Planck corre: M²_eff(a) = M²_Pl · h(a)
  → G_eff(a) = G_N / h(a)
  Con αM(a) = αM0 · Ω_DE(a)/Ω_DE,0 [tracker]:
    - αM = αM0 hoy (máximo, constrained)
    - αM → 0 en z>>1 (GR en universo temprano)
  Geff/GN combinado (cuasi-estático, αT=0, c_s²=0):
    μ(a) = [1 + (1+αB(a))·Ω_DE(a)/Ω_m(a)] / (1 + αM(a))
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
plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'figure.dpi': 150})

# ── Constantes SSEE ──────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi_  = np.pi
beta = (pi_ + phi) / 2
Kv   = phi + pi_ + (pi_ + phi)
Tr   = 3 * (phi + beta)
Mv   = phi + pi_ + Kv

w0      = -Tr / Mv
wa      = -(pi_ + phi + phi) / Kv
Omm_dyn = 1.0 + w0
OmDE    = 1.0 - Omm_dyn
MIRA    = (3*phi + pi_) / 4

# αB₀ algebraico (derivado en Paper 6a)
aB0_ssee = (MIRA - 1.0) * (pi_ - phi) / (3*phi + pi_) - 1.0

def f_DE_raw(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
_norm = f_DE_raw(1.0)
def f_DE(a):
    return f_DE_raw(a) / _norm

def E2(a):
    return max(Omm_dyn * a**(-3) + OmDE * f_DE(a), 1e-30)

# ── Modelos αB(a) y αM(a) ────────────────────────────────────────────────────
def aB_func(a, aB0):
    return aB0 * f_DE(a)          # tracker DE

def aM_func(a, aM0):
    return aM0 * f_DE(a)          # tracker DE

# ── G_eff/G_N combinado ──────────────────────────────────────────────────────
def G_eff_ratio(a, aB0, aM0):
    OmDE_a = OmDE * f_DE(a) / E2(a)
    Omm_a  = Omm_dyn * a**(-3) / E2(a)
    aB = aB_func(a, aB0)
    aM = aM_func(a, aM0)
    num = 1.0 + (1.0 + aB) * OmDE_a / max(Omm_a, 1e-10)
    den = 1.0 + aM
    return num / max(den, 1e-10)

# ── Ecuación de crecimiento ───────────────────────────────────────────────────
def growth_rhs(lna, Y, aB0, aM0):
    a = np.exp(lna)
    D, Dp = Y
    e2 = E2(a)
    h  = 0.5 * (-3*Omm_dyn*a**(-3) + OmDE*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a) / e2
    Omm_a = Omm_dyn * a**(-3) / e2
    mu    = G_eff_ratio(a, aB0, aM0)
    Dpp   = -(2.0 + h) * Dp + 1.5 * Omm_a * mu * D
    return [Dp, Dpp]

def growth_lcdm_rhs(lna, Y, Om=0.3153):
    a = np.exp(lna)
    D, Dp = Y
    e2 = Om*a**(-3) + (1-Om)
    h  = 0.5*(-3*Om*a**(-3))/e2
    Omm_a = Om*a**(-3)/e2
    return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]

lna_pts = np.linspace(np.log(0.01), 0.0, 500)
sol_lcdm = solve_ivp(growth_lcdm_rhs, [lna_pts[0], 0.0], [0.01, 0.01],
                     t_eval=lna_pts, rtol=1e-10, atol=1e-13)
D_lcdm  = sol_lcdm.y[0] / sol_lcdm.y[0][-1]
Dp_lcdm = sol_lcdm.y[1] / sol_lcdm.y[0][-1]
sig8_lcdm = 0.811

# ── fσ8 observado ────────────────────────────────────────────────────────────
Z_RSD     = np.array([0.067, 0.15, 0.32, 0.38, 0.51, 0.57])
fsig8_obs = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
fsig8_err = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])

def compute_fsig8_tension(aB0, aM0):
    """Devuelve tensión media |σ| de fσ8 para un par (αB0, αM0)."""
    sol = solve_ivp(lambda t, Y: growth_rhs(t, Y, aB0, aM0),
                    [lna_pts[0], 0.0], [0.01, 0.01],
                    t_eval=lna_pts, rtol=1e-9, atol=1e-12)
    if not sol.success:
        return np.nan
    D  = sol.y[0] / sol.y[0][-1]
    Dp = sol.y[1] / sol.y[0][-1]

    # σ8 SSEE: reducido por la baja densidad
    sig8 = sig8_lcdm * np.sqrt(Omm_dyn / 0.3153)

    tensions = []
    for i, z in enumerate(Z_RSD):
        a_t = 1/(1+z)
        lna_t = np.log(a_t)
        f_t = float(np.interp(lna_t, lna_pts, Dp/D))
        D_t = float(np.interp(lna_t, lna_pts, D))
        fs8 = f_t * sig8 * D_t
        tensions.append(abs(fsig8_obs[i] - fs8) / fsig8_err[i])
    return float(np.mean(tensions))

# ── SCAN 2D (αB0, αM0) ───────────────────────────────────────────────────────
print("=" * 65)
print("  SSEE Paper 6 — 2D Scan (αB₀, αM₀)")
print("  Límite observacional: |αM₀| < 0.08 (Planck+GW)")
print("=" * 65)

# Rango observacionalmente motivado
aB_range = np.linspace(-1.0, 0.0, 20)
aM_range = np.linspace(-0.08, 0.08, 20)

tension_grid = np.full((len(aM_range), len(aB_range)), np.nan)
for i, aM0 in enumerate(aM_range):
    for j, aB0 in enumerate(aB_range):
        tension_grid[i, j] = compute_fsig8_tension(aB0, aM0)

best_idx = np.unravel_index(np.nanargmin(tension_grid), tension_grid.shape)
best_aM = aM_range[best_idx[0]]
best_aB = aB_range[best_idx[1]]
best_t  = tension_grid[best_idx]

# Tensión del caso base (sin EFT) y con αB₀ SSEE
t_base       = compute_fsig8_tension(0.0, 0.0)
t_aB_only    = compute_fsig8_tension(aB0_ssee, 0.0)
t_best_combo = best_t

print(f"\n  Sin EFT (Paper 5 original):     {t_base:.2f}σ")
print(f"  Con αB₀ algebraico (Paper 6a):  {t_aB_only:.2f}σ")
print(f"  Mejor combo 2D:                  {t_best_combo:.2f}σ")
print(f"    → αB₀_best = {best_aB:.4f}")
print(f"    → αM₀_best = {best_aM:.4f}")
print(f"\n  αB₀ SSEE algebraico = {aB0_ssee:.4f}")

# ── FIGURAS ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel 1: mapa de tensión 2D
ax = axes[0]
im = ax.contourf(aB_range, aM_range, tension_grid, levels=20, cmap='RdYlGn_r')
plt.colorbar(im, ax=ax, label=r'Mean $|tension|$ ($\sigma$)')
ax.contour(aB_range, aM_range, tension_grid, levels=[1.0, 2.0, 3.0],
           colors='white', linewidths=1, linestyles=['--', '-', ':'])
ax.scatter([best_aB], [best_aM], color='cyan', s=100, zorder=5,
           label=f'Best: {best_t:.2f}σ')
ax.scatter([aB0_ssee], [0.0], color='gold', marker='*', s=200, zorder=5,
           label=rf'SSEE algebraic $\alpha_B={aB0_ssee:.3f}$')
ax.axhline(0.02, color='white', ls=':', lw=1.2, label='Obs. bound |αM|<0.02')
ax.axhline(-0.02, color='white', ls=':', lw=1.2)
ax.set_xlabel(r'$\alpha_{B,0}$')
ax.set_ylabel(r'$\alpha_{M,0}$')
ax.set_title(r'$f\sigma_8$ Tension Map: $(\alpha_B, \alpha_M)$ Scan')
ax.legend(fontsize=7, loc='upper left')
ax.grid(True, alpha=0.2)

# Panel 2: curvas fσ8 comparadas
ax = axes[1]
ax.errorbar(Z_RSD, fsig8_obs, yerr=fsig8_err, fmt='ko', ms=6,
            capsize=4, label='Observed (RSD)', zorder=5)
z_s = np.linspace(0.05, 0.65, 100)

def fsig8_curve(z_arr, aB0, aM0):
    sig8 = sig8_lcdm * np.sqrt(Omm_dyn / 0.3153)
    sol = solve_ivp(lambda t, Y: growth_rhs(t, Y, aB0, aM0),
                    [lna_pts[0], 0.0], [0.01, 0.01],
                    t_eval=lna_pts, rtol=1e-9, atol=1e-12)
    D  = sol.y[0] / sol.y[0][-1]
    Dp = sol.y[1] / sol.y[0][-1]
    out = []
    for z in z_arr:
        f_t = float(np.interp(np.log(1/(1+z)), lna_pts, Dp/D))
        D_t = float(np.interp(np.log(1/(1+z)), lna_pts, D))
        out.append(f_t * sig8 * D_t)
    return np.array(out)

ax.plot(z_s, fsig8_curve(z_s, 0.0, 0.0), 'r-', lw=2,
        label=fr'SSEE base (no EFT), {t_base:.2f}σ')
ax.plot(z_s, fsig8_curve(z_s, aB0_ssee, 0.0), 'b--', lw=2,
        label=fr'SSEE + $\alpha_B={aB0_ssee:.3f}$, {t_aB_only:.2f}σ')
ax.plot(z_s, fsig8_curve(z_s, best_aB, best_aM), 'm-', lw=2,
        label=fr'Best combo ($\alpha_B={best_aB:.2f}, \alpha_M={best_aM:.3f}$), {best_t:.2f}σ')
# ΛCDM
fs8_lcdm_curve = []
for z in z_s:
    f_t = float(np.interp(np.log(1/(1+z)), lna_pts, Dp_lcdm/D_lcdm))
    D_t = float(np.interp(np.log(1/(1+z)), lna_pts, D_lcdm))
    fs8_lcdm_curve.append(f_t * sig8_lcdm * D_t)
ax.plot(z_s, fs8_lcdm_curve, 'g:', lw=2, label=r'$\Lambda$CDM')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$f\sigma_8(z)$')
ax.set_title(r'$f\sigma_8$ Comparison')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

fig.suptitle('SSEE Paper 6: EFT Extension Scan\n'
             r'Exploring $(\alpha_B, \alpha_M)$ with $\alpha_T=0$, $c^2_s=0$',
             fontsize=11)
fig.tight_layout()
out = os.path.join(OUTDIR, 'fig_paper6_alphaM_scan.pdf')
fig.savefig(out, bbox_inches='tight')
print(f"\n  Figura guardada: {out}")
print("=" * 65)
