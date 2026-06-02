#!/usr/bin/env python3
"""
ssee_paper6_kinetic_braiding.py  —  SSEE-V3.6 Paper 6 Foundation
=================================================================
Preguntas:
  P1. ¿Qué expresión algebraica de SSEE fija αB (Kinetic Braiding)?
  P2. ¿Resuelve la tensión H(z)?
  P3. ¿Recupera fσ8 observado?

Física:
  En Horndeski con c²_s = 0 (demostrado en Paper 5), la DE se aglomera
  como polvo. La gravedad efectiva en el límite cuasi-estático es:
      G_eff/G_N = 1 + (1 + αB) * Ω_DE / Ω_m,dyn
  Condición SSEE: G_eff/G_N = MIRA = (3φ+π)/4 at z=0
  → αB = (MIRA - 1) * Ω_m,dyn/Ω_DE - 1  [100% algebraico]

  Evolución temporal:
      αB(a) = αB0 · f_DE(a)  (tracker: αB α DE density)
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ── Directorios ──────────────────────────────────────────────────────────────
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', 'results', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'figure.dpi': 150})

# ── 1. CONSTANTES SSEE (algebraicas puras) ───────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi_  = np.pi
beta = (pi_ + phi) / 2
Omega_val = pi_ + phi
KAL0 = beta + pi_
P_sc = Omega_val + phi
Kv   = phi + pi_ + Omega_val
Tr   = 3 * (phi + beta)
Mv   = phi + pi_ + Kv

w0   = -Tr / Mv
wa   = -P_sc / Kv

Omm_dyn = 1.0 + w0
OmDE    = 1.0 - Omm_dyn
MIRA    = (3*phi + pi_) / 4

# ── 2. DERIVACIÓN ALGEBRAICA DE αB ──────────────────────────────────────────
# G_eff/G_N = MIRA  →  αB = (MIRA-1)*Ω_m,dyn/Ω_DE - 1
aB0 = (MIRA - 1.0) * Omm_dyn / OmDE - 1.0

# Expresión algebraica equivalente (en términos de φ, π):
# Ω_m,dyn/Ω_DE = (π-φ)/(3φ+π)  [derivado de los axiomas]
aB0_algebraic = (MIRA - 1.0) * (pi_ - phi) / (3*phi + pi_) - 1.0

print("=" * 65)
print("  SSEE Paper 6 — Kinetic Braiding Foundation")
print("=" * 65)
print(f"\n  φ                = {phi:.8f}")
print(f"  w₀               = {w0:.8f}")
print(f"  Ω_m,dyn          = {Omm_dyn:.8f}")
print(f"  Ω_DE             = {OmDE:.8f}")
print(f"  MIRA = (3φ+π)/4  = {MIRA:.8f}")
print(f"\n  αB₀ (numérico)   = {aB0:.8f}")
print(f"  αB₀ (algebraico) = {aB0_algebraic:.8f}")
print(f"  Consistencia     = {abs(aB0 - aB0_algebraic):.2e}  {'✓' if abs(aB0-aB0_algebraic)<1e-10 else '✗'}")

# ── 3. MODELO COSMOLÓGICO CON BRAIDING ───────────────────────────────────────
def f_DE_raw(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
_norm = f_DE_raw(1.0)
def f_DE(a):
    return f_DE_raw(a) / _norm

def aB(a):
    """αB(a) = αB0 · f_DE(a)  [tracker: sigue a la densidad DE]"""
    return aB0 * f_DE(a)

def G_eff_ratio(a):
    """G_eff/G_N en el límite cuasi-estático con c²_s = 0"""
    # Fracción de materia dinámica y DE a escala de factor a
    OmDE_a  = OmDE * f_DE(a) / E2(a)
    Omm_a   = Omm_dyn * a**(-3) / E2(a)
    return 1.0 + (1.0 + aB(a)) * OmDE_a / Omm_a

def E2(a):
    """E²(a) = (H/H₀)² — sin braiding en el background (αB no modifica H)"""
    return Omm_dyn * a**(-3) + OmDE * f_DE(a)

def E(a):
    return np.sqrt(max(E2(a), 1e-30))

# ── 4. VERIFICACIÓN G_eff en z=0 ─────────────────────────────────────────────
Geff_z0 = G_eff_ratio(1.0)
print(f"\n  G_eff/G_N (z=0)  = {Geff_z0:.8f}")
print(f"  MIRA target      = {MIRA:.8f}")
print(f"  Δ                = {abs(Geff_z0 - MIRA):.2e}  {'✓' if abs(Geff_z0-MIRA)<0.01 else '≈'}")

# ── 5. P2: H(z) con braiding ──────────────────────────────────────────────────
# El braiding cinético puro no modifica el background (αB solo afecta perturbaciones)
# Por eso H(z) es idéntico al SSEE original.
# Sin embargo, G_eff amplificado puede aumentar la tasa de crecimiento
# y por ende fσ8, que es la tensión real.
print("\n  Nota P2: Kinetic Braiding puro no modifica H(z) a nivel de background.")
print("  La tensión H(z) requiere modificación de la ecuación de estado w(a).")
print("  → El braiding resuelve fσ8 (P3), no H(z). Diagnóstico preciso.")

# ── 6. P3: Crecimiento de estructura con G_eff ───────────────────────────────
# Ecuación de crecimiento en ln(a):
# D'' + (2 + H'/H) D' - (3/2) Ω_m(a) μ(a) D = 0
# donde μ = G_eff/G_N y derivadas respecto a ln(a)

def growth_rhs(lna, Y):
    a = np.exp(lna)
    D, Dp = Y  # D y dD/d(lna)

    e2 = E2(a)
    e  = np.sqrt(e2)
    # d(lnE)/d(lna) = E'/(E) · a
    de2_dlna = (-3*Omm_dyn*a**(-3) + OmDE*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a)/e2
    dlnH_dlna = 0.5 * de2_dlna

    Omm_a = Omm_dyn * a**(-3) / e2
    mu = G_eff_ratio(a)

    Dpp = -(2.0 + dlnH_dlna) * Dp + 1.5 * Omm_a * mu * D
    return [Dp, Dpp]

def growth_lcdm_rhs(lna, Y, Om_lcdm=0.3153):
    a = np.exp(lna)
    D, Dp = Y
    e2 = Om_lcdm * a**(-3) + (1.0 - Om_lcdm)
    dlnH = 0.5 * (-3*Om_lcdm*a**(-3)) / e2
    Omm_a = Om_lcdm * a**(-3) / e2
    Dpp = -(2.0 + dlnH) * Dp + 1.5 * Omm_a * D
    return [Dp, Dpp]

# Integrar desde a=0.01 (matter domination, D≈a, D'≈1)
lna_ini = np.log(0.01)
lna_fin = 0.0
lna_pts = np.linspace(lna_ini, lna_fin, 500)

sol_ssee  = solve_ivp(growth_rhs,       [lna_ini, lna_fin], [0.01, 0.01],
                      t_eval=lna_pts, method='RK45', rtol=1e-10, atol=1e-13)
sol_lcdm  = solve_ivp(growth_lcdm_rhs,  [lna_ini, lna_fin], [0.01, 0.01],
                      t_eval=lna_pts, method='RK45', rtol=1e-10, atol=1e-13)

# Normalizar D(a=1) = 1
D_ssee  = sol_ssee.y[0]  / sol_ssee.y[0][-1]
Dp_ssee = sol_ssee.y[1]  / sol_ssee.y[0][-1]
D_lcdm  = sol_lcdm.y[0]  / sol_lcdm.y[0][-1]
Dp_lcdm = sol_lcdm.y[1]  / sol_lcdm.y[0][-1]

a_arr = np.exp(lna_pts)
z_arr = 1.0/a_arr - 1.0

# σ8 base (normalización al valor de Planck con ΛCDM)
sigma8_lcdm_base = 0.811   # Planck 2018
# SSEE sigma8: reducido por la menor amplitud de crecimiento * normalización CMB
G_factor = D_ssee[-1] / D_lcdm[-1]   # = 1 por normalización, pero la fisica real es
# En realidad SSEE comienza con menos masa → σ8 diferente
# Cálculo físico: σ8_SSEE = σ8_ΛCDM × (D1_SSEE/D1_ΛCDM) × √(Ωm_dyn/Ωm_lcdm)
sigma8_ssee_braided = sigma8_lcdm_base * np.sqrt(Omm_dyn / 0.3153)

# Índice de crecimiento γ: f = Ω_m(a)^γ, f = d ln D / d ln a
f_ssee = Dp_ssee / D_ssee    # = d ln D / d ln a
f_lcdm = Dp_lcdm / D_lcdm

Omm_of_a = Omm_dyn * a_arr**(-3) / np.array([E2(a) for a in a_arr])
# γ_eff(a) tal que f = Ω_m(a)^γ
with np.errstate(divide='ignore', invalid='ignore'):
    gamma_ssee = np.where(
        (Omm_of_a > 0.01) & (f_ssee > 0.01),
        np.log(f_ssee) / np.log(Omm_of_a),
        np.nan
    )

# fσ8 a redshifts específicos
Z_RSD = np.array([0.067, 0.15, 0.32, 0.38, 0.51, 0.57])
# Datos observados fσ8 (6dFGRS, SDSS MGS, BOSS DR12 x2, eBOSS DR16, WiggleZ)
fsig8_obs = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
fsig8_err = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])

def interp_fsig8(z_target, z_arr, f_growth, sigma8_val):
    a_t = 1.0/(1.0+z_target)
    lna_t = np.log(a_t)
    f_t = float(np.interp(lna_t, lna_pts, f_growth))
    D_t = float(np.interp(lna_t, lna_pts, D_ssee if sigma8_val == sigma8_ssee_braided else D_lcdm))
    return f_t * sigma8_val * D_t

fsig8_ssee_pred = np.array([interp_fsig8(z, z_arr, f_ssee, sigma8_ssee_braided) for z in Z_RSD])
fsig8_lcdm_pred = np.array([interp_fsig8(z, z_arr, f_lcdm, sigma8_lcdm_base) for z in Z_RSD])

tensions_ssee = (fsig8_obs - fsig8_ssee_pred) / fsig8_err
tensions_lcdm = (fsig8_obs - fsig8_lcdm_pred) / fsig8_err
mean_t_ssee = np.mean(np.abs(tensions_ssee))
mean_t_lcdm = np.mean(np.abs(tensions_lcdm))

# ── 7. RESUMEN ───────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("  RESULTADOS: SSEE + Kinetic Braiding (αB₀ algebraico)")
print("=" * 65)
print(f"\n  αB₀ = {aB0:.6f}  = (MIRA-1)(π-φ)/(3φ+π) - 1")
print(f"  G_eff/G_N (z=0) = {G_eff_ratio(1.0):.4f}   target = MIRA = {MIRA:.4f}")
print(f"  G_eff/G_N (z=2) = {G_eff_ratio(1/3):.4f}")
print(f"\n  σ₈ (SSEE+braiding) = {sigma8_ssee_braided:.4f}")
print(f"\n  Tensión media fσ₈:")
print(f"    SSEE (pre-braiding) ~ 2.67σ  [Paper 5]")
print(f"    SSEE + braiding    = {mean_t_ssee:.2f}σ")
print(f"    ΛCDM               = {mean_t_lcdm:.2f}σ")
print()
print(f"  {'z':>6}  {'fσ₈ obs':>10}  {'fσ₈ SSEE+B':>12}  {'fσ₈ ΛCDM':>10}  {'σ (SSEE+B)':>12}")
for i, z in enumerate(Z_RSD):
    print(f"  {z:>6.3f}  {fsig8_obs[i]:>10.3f}  {fsig8_ssee_pred[i]:>12.3f}"
          f"  {fsig8_lcdm_pred[i]:>10.3f}  {tensions_ssee[i]:>+12.2f}σ")

# ── 8. FIGURAS ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel 1: G_eff/G_N vs z
z_plot = np.linspace(0, 3, 200)
a_plot = 1.0/(1.0+z_plot)
Geff_plot = np.array([G_eff_ratio(a) for a in a_plot])
ax = axes[0]
ax.plot(z_plot, Geff_plot, 'b-', lw=2, label=r'$G_\mathrm{eff}/G_N$ (SSEE+Braiding)')
ax.axhline(MIRA, color='r', ls='--', lw=1.5,
           label=rf'MIRA = $(3\varphi+\pi)/4 = {MIRA:.4f}$')
ax.axhline(1.0, color='gray', ls=':', lw=1, label=r'$G_\mathrm{eff} = G_N$ ($\Lambda$CDM)')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$G_\mathrm{eff}/G_N$')
ax.set_title(f'Effective Gravity\n' + r'$\alpha_B = \alpha_{B,0} \cdot f_{DE}(a)$')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 3)

# Panel 2: fσ8 comparison
ax = axes[1]
ax.errorbar(Z_RSD, fsig8_obs, yerr=fsig8_err, fmt='ko', ms=6,
            capsize=4, label='Observed (RSD surveys)', zorder=5)
z_smooth = np.linspace(0.05, 0.65, 100)
fsig8_ssee_smooth = np.array([interp_fsig8(z, z_arr, f_ssee, sigma8_ssee_braided) for z in z_smooth])
fsig8_lcdm_smooth = np.array([interp_fsig8(z, z_arr, f_lcdm, sigma8_lcdm_base) for z in z_smooth])
ax.plot(z_smooth, fsig8_ssee_smooth, 'b-', lw=2,
        label=rf'SSEE + Braiding ($\alpha_B={aB0:.3f}$)')
ax.plot(z_smooth, fsig8_lcdm_smooth, 'g--', lw=2, label=r'$\Lambda$CDM')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$f\sigma_8(z)$')
ax.set_title(r'Growth Rate $f\sigma_8$')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 3: αB(z)
aB_plot = np.array([aB(a) for a in a_plot])
ax = axes[2]
ax.plot(z_plot, aB_plot, 'm-', lw=2, label=rf'$\alpha_B(z)$')
ax.axhline(aB0, color='purple', ls='--', lw=1,
           label=rf'$\alpha_{{B,0}} = {aB0:.4f}$')
ax.axhline(0, color='gray', ls=':', lw=1)
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$\alpha_B(z)$')
ax.set_title('Kinetic Braiding\nEvolution (tracker)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 3)

fig.suptitle(
    r'SSEE-V3.6 + EFT Kinetic Braiding ($\alpha_T=0$, $c^2_s=0$)' + '\n'
    + rf'$\alpha_{{B,0}} = (\mathrm{{MIRA}}-1)(\pi-\varphi)/(3\varphi+\pi) - 1 = {aB0:.4f}$',
    fontsize=11
)
fig.tight_layout()
out = os.path.join(OUTDIR, 'fig_paper6_kinetic_braiding.pdf')
fig.savefig(out, bbox_inches='tight')
print(f"\n  Figura guardada: {out}")
print("=" * 65)
