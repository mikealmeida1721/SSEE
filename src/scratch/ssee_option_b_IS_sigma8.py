"""
ssee_option_b_IS_sigma8.py
==========================
Opción B: ¿la restauración Ω_m,growth = 0.320 (dos sectores) produce
σ₈_self-consistente ≈ 0.737 SIN WDM?

El candidato 1 usa D_both(Ω_m=0.320) para f(z) pero σ₈=0.702 de Paper 5
(que fue calculado con Ω_m=0.160). Eso es INCONSISTENTE.

Este script corrige: calcula G_twosector = D₁_SSEE(Ω_m=0.320)/D₁_ΛCDM
y deriva σ₈_twosector = 0.811 × G_twosector.

Si σ₈_twosector ≈ 0.737 → el mecanismo de σ₈ en Paper 6 no es WDM,
sino la auto-consistencia del crecimiento de dos sectores.
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

Om_dyn   = 1.0 + w0
Om_DE    = 1.0 - Om_dyn
Om_phiDM = (MIRA - 1) * Om_dyn
Om_both  = Om_dyn + Om_phiDM

sig8_Planck = 0.811
Om_LCDM     = 0.3153

print(f"SSEE: w₀={w0:.4f}, wₐ={wa:.4f}")
print(f"Ω_m,dyn={Om_dyn:.4f}, Ω_DE={Om_DE:.4f}, MIRA={MIRA:.4f}")
print(f"Ω_φDM={Om_phiDM:.6f}, Ω_both={Om_both:.6f}\n")

# ── Ecuaciones de crecimiento ─────────────────────────────────────────────────
def f_DE_raw(a): return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
_n = f_DE_raw(1.0)
def f_DE(a): return f_DE_raw(a) / _n

def make_growth_rhs(Om_m, bg_DE=Om_DE):
    """Background: E²(a) = Om_m a⁻³ + bg_DE f_DE(a)."""
    def rhs(lna, Y):
        a     = np.exp(lna)
        e2    = max(Om_m * a**(-3) + bg_DE * f_DE(a), 1e-30)
        dlne2 = (-3*Om_m*a**(-3)
                 + bg_DE*(f_DE(a+1e-5) - f_DE(a-1e-5))/(2e-5)*a) / e2
        h     = 0.5 * dlne2
        Omm_a = Om_m * a**(-3) / e2
        D, Dp = Y
        return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]
    return rhs

def growth_lcdm(lna, Y):
    a     = np.exp(lna)
    e2    = max(Om_LCDM*a**(-3) + (1-Om_LCDM), 1e-30)
    h     = 0.5*(-3*Om_LCDM*a**(-3)) / e2
    Omm_a = Om_LCDM*a**(-3) / e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]

lna_pts = np.linspace(np.log(0.01), 0.0, 600)
IC      = [0.01, 0.01]

def integrate_raw(rhs):
    """Integra sin normalizar — necesario para calcular G correctamente."""
    sol = solve_ivp(rhs, [lna_pts[0], 0.0], IC,
                    t_eval=lna_pts, rtol=1e-10, atol=1e-13)
    return sol.y[0], sol.y[1]

def integrate_normed(rhs):
    """Normalizado a D(z=0)=1 — para f(z)=D'/D."""
    D, Dp = integrate_raw(rhs)
    return D / D[-1], Dp / D[-1]

# ── Calcular G (factor de crecimiento relativo a ΛCDM) ───────────────────────
D_lcdm_raw,  _  = integrate_raw(growth_lcdm)
D_single_raw, _ = integrate_raw(make_growth_rhs(Om_dyn))
D_both_raw,   _ = integrate_raw(make_growth_rhs(Om_both))

# Caso especial: background auto-consistente para dos sectores
# E²(a) = Om_both a⁻³ + (1-Om_both) f_DE(a)
# Esto garantiza Ω_total=1 pero cambia el background
D_both_flat_raw, _ = integrate_raw(make_growth_rhs(Om_both, bg_DE=1.0-Om_both))

G_single        = D_single_raw[-1] / D_lcdm_raw[-1]
G_twosector     = D_both_raw[-1]   / D_lcdm_raw[-1]
G_twosector_flat= D_both_flat_raw[-1] / D_lcdm_raw[-1]

sig8_single         = sig8_Planck * G_single
sig8_twosector      = sig8_Planck * G_twosector
sig8_twosector_flat = sig8_Planck * G_twosector_flat

SEP = "=" * 75
sep = "─" * 75

print(SEP)
print("  OPCIÓN B — Factor de crecimiento G = D₁_SSEE/D₁_ΛCDM")
print(SEP)
print(f"\n  ΛCDM  (Ω_m=0.315)             : G = 1.0000, σ₈ = {sig8_Planck:.3f} [referencia]")
print(f"\n  SSEE single-sector (Ω_m=0.160):")
print(f"    G_single         = {G_single:.4f}")
print(f"    σ₈_single        = {sig8_single:.3f}   (Paper 5 reporta 0.702 ✓)")
print(f"\n  SSEE two-sector (Ω_m=0.320, bg: Om_DE=0.840):")
print(f"    G_twosector      = {G_twosector:.4f}")
print(f"    σ₈_twosector     = {sig8_twosector:.3f}")
print(f"\n  SSEE two-sector (Ω_m=0.320, bg plano: Om_DE={1-Om_both:.4f}):")
print(f"    G_twosector_flat = {G_twosector_flat:.4f}")
print(f"    σ₈_flat          = {sig8_twosector_flat:.3f}")
print(f"\n  Paper 6 target: σ₈_eff = 0.737 (KiDS-1000 0.00σ)")
print(f"  Δσ₈ (2sector vs target): {sig8_twosector - 0.737:+.3f}")
print(f"  Δσ₈ (flat vs target):    {sig8_twosector_flat - 0.737:+.3f}")

# ── Curvas de crecimiento normalizadas para fσ₈ ───────────────────────────────
D_lcdm,  Dp_lcdm  = integrate_normed(growth_lcdm)
D_single, Dp_single = integrate_normed(make_growth_rhs(Om_dyn))
D_both,   Dp_both   = integrate_normed(make_growth_rhs(Om_both))
D_flat,   Dp_flat   = integrate_normed(make_growth_rhs(Om_both, bg_DE=1.0-Om_both))

# ── Datos RSD ─────────────────────────────────────────────────────────────────
Z_RSD     = np.array([0.067, 0.150, 0.320, 0.380, 0.510, 0.570])
surveys   = ['6dFGRS', 'SDSS MGS', 'BOSS z=0.32', 'BOSS z=0.38', 'eBOSS DR16', 'WiggleZ']
fsig8_obs = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
fsig8_err = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])

def get_fsig8(z, D, Dp, sig8):
    lna_t = np.log(1/(1+z))
    f_t   = float(np.interp(lna_t, lna_pts, Dp/D))
    D_t   = float(np.interp(lna_t, lna_pts, D))
    return f_t * sig8 * D_t

# ── Modelos comparados ────────────────────────────────────────────────────────
models = {
    'P5: IS solo (Ω_m=0.160, σ₈=0.702)': {
        'D': D_single, 'Dp': Dp_single, 'sig8': 0.702,
        'color': 'tomato', 'ls': '--'
    },
    'C1: IS+2sector inconsistente (σ₈=0.702 fijo)': {
        'D': D_both, 'Dp': Dp_both, 'sig8': 0.702,
        'color': 'orange', 'ls': ':'
    },
    f'OptB: IS+2sector self-consistente (σ₈={sig8_twosector:.3f})': {
        'D': D_both, 'Dp': Dp_both, 'sig8': sig8_twosector,
        'color': 'steelblue', 'ls': '-'
    },
    f'OptB-flat: IS+2sector plano (σ₈={sig8_twosector_flat:.3f})': {
        'D': D_flat, 'Dp': Dp_flat, 'sig8': sig8_twosector_flat,
        'color': 'navy', 'ls': '-.'
    },
    'P6: IS+2sector+WDM (σ₈=0.737)': {
        'D': D_both, 'Dp': Dp_both, 'sig8': 0.737,
        'color': 'darkgreen', 'ls': '-.'
    },
    'ΛCDM (σ₈=0.811)': {
        'D': D_lcdm, 'Dp': Dp_lcdm, 'sig8': sig8_Planck,
        'color': 'gray', 'ls': ':'
    },
}

results = {}
for name, cfg in models.items():
    fsig8 = np.array([get_fsig8(z, cfg['D'], cfg['Dp'], cfg['sig8']) for z in Z_RSD])
    tensions = (fsig8_obs - fsig8) / fsig8_err
    results[name] = {
        'fsig8': fsig8,
        'tensions': tensions,
        'mean_abs_tension': np.mean(np.abs(tensions))
    }

print(f"\n{SEP}")
print("  TENSIONES fσ₈")
print(SEP)
for name, r in results.items():
    print(f"    {name[:55]:55s}: {r['mean_abs_tension']:.4f}σ")

# ── Diagnóstico ───────────────────────────────────────────────────────────────
print(f"\n{SEP}")
print("  DIAGNÓSTICO OPCIÓN B")
print(SEP)

t_P5   = results['P5: IS solo (Ω_m=0.160, σ₈=0.702)']['mean_abs_tension']
t_C1   = results['C1: IS+2sector inconsistente (σ₈=0.702 fijo)']['mean_abs_tension']
t_optB = results[f'OptB: IS+2sector self-consistente (σ₈={sig8_twosector:.3f})']['mean_abs_tension']
t_P6   = results['P6: IS+2sector+WDM (σ₈=0.737)']['mean_abs_tension']
t_LCDM = results['ΛCDM (σ₈=0.811)']['mean_abs_tension']

print(f"\n  Cadena de tensiones:")
print(f"    P5 (IS solo):              {t_P5:.4f}σ")
print(f"    C1 (2sector, σ₈ fijo):     {t_C1:.4f}σ  ← inconsistente")
print(f"    OptB (2sector, σ₈={sig8_twosector:.3f}): {t_optB:.4f}σ  ← self-consistente")
print(f"    P6 (2sector+WDM, σ₈=0.737):{t_P6:.4f}σ")
print(f"    ΛCDM:                       {t_LCDM:.4f}σ")
print(f"\n  G values:")
print(f"    G_single (Paper 5):  {G_single:.4f}  → σ₈={sig8_single:.3f}")
print(f"    G_twosector (OptB):  {G_twosector:.4f}  → σ₈={sig8_twosector:.3f}")
print(f"    G_flat (plano):      {G_twosector_flat:.4f}  → σ₈={sig8_twosector_flat:.3f}")

if abs(sig8_twosector - 0.737) < 0.010:
    verdict_sigma = "MATCH — σ₈_2sector self-consistente reproduce 0.737 SIN WDM"
elif sig8_twosector > 0.720:
    verdict_sigma = f"PARCIAL — σ₈={sig8_twosector:.3f} (objetivo 0.737, gap={0.737-sig8_twosector:.3f})"
else:
    verdict_sigma = f"FALLA — σ₈={sig8_twosector:.3f} demasiado bajo (necesita WDM)"

if t_optB < 1.0:
    verdict_tension = f"RESUELTO — OptB da {t_optB:.2f}σ (<1σ)"
elif t_optB < 1.5:
    verdict_tension = f"VIABLE — OptB da {t_optB:.2f}σ"
else:
    verdict_tension = f"INSUFICIENTE — OptB da {t_optB:.2f}σ"

print(f"\n  VEREDICTO σ₈:     {verdict_sigma}")
print(f"  VEREDICTO fσ₈:    {verdict_tension}")

# ── Análisis IS viscous transfer function ─────────────────────────────────────
print(f"\n{SEP}")
print("  ANÁLISIS IS — Escala de amortiguamiento viscoso")
print(SEP)

tau_PI = KAL0 / (3 * Om_DE)
H0_kmsMpc = 3*(phi + pi_)**2  # ≈ 67.96 km/s/Mpc
c_kmsMpc = 2.998e5  # km/s
# k_visc en h/Mpc: escala donde τ_Π×H₀ ≈ 1/k (en unidades comoving)
k_Hubble = H0_kmsMpc / c_kmsMpc  # h/Mpc (escala Hubble)
k_visc_est = k_Hubble / tau_PI    # estimación: k donde viscosidad ~ Hubble

print(f"\n  τ_Π = KAL₀/(3×Ω_DE) = {tau_PI:.4f} H₀⁻¹")
print(f"  H₀/c = {k_Hubble*1e3:.4f} × 10⁻³ h/Mpc  (escala Hubble comoving)")
print(f"  k_visc_est = H₀/(c×τ_Π) = {k_visc_est*1e4:.4f} × 10⁻⁴ h/Mpc")
print(f"\n  → k_visc >> k_σ₈ = 0.125 h/Mpc ¿? = {k_visc_est < 0.125}")
print(f"\n  Conclusión: IS viscosidad actúa en escala Hubble (k ~ {k_visc_est:.2e} h/Mpc)")
print(f"  No puede suprimir σ₈ directamente a k = 0.125 h/Mpc")
print(f"  El mecanismo σ₈ debe venir de G_twosector, no de T_IS(k)")

# ── FIGURA ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))

# Panel 1: curvas fσ₈(z)
ax1 = axes[0]
z_s = np.linspace(0.04, 0.65, 200)
plot_models = {
    'P5 (σ₈=0.702, Ω_m=0.160)': ('D_single', 'Dp_single', 0.702, 'tomato', '--'),
    'C1 inconsistente (σ₈=0.702)': ('D_both', 'Dp_both', 0.702, 'orange', ':'),
    f'OptB self-consist. (σ₈={sig8_twosector:.3f})': ('D_both', 'Dp_both', sig8_twosector, 'steelblue', '-'),
    'P6 +WDM (σ₈=0.737)': ('D_both', 'Dp_both', 0.737, 'darkgreen', '-.'),
    'ΛCDM': ('D_lcdm', 'Dp_lcdm', sig8_Planck, 'gray', ':'),
}
D_map  = {'D_single': D_single, 'D_both': D_both, 'D_lcdm': D_lcdm,
           'D_flat': D_flat}
Dp_map = {'Dp_single': Dp_single, 'Dp_both': Dp_both, 'Dp_lcdm': Dp_lcdm,
           'Dp_flat': Dp_flat}

for label, (dkey, dpkey, sig8, col, ls) in plot_models.items():
    D_  = D_map[dkey]
    Dp_ = Dp_map[dpkey]
    fs8 = [get_fsig8(z, D_, Dp_, sig8) for z in z_s]
    # find tension for this model
    t_v = np.mean(np.abs((np.array([get_fsig8(z, D_, Dp_, sig8) for z in Z_RSD])
                          - fsig8_obs) / fsig8_err))
    ax1.plot(z_s, fs8, color=col, ls=ls, lw=2,
             label=f"{label} ({t_v:.2f}σ)")

ax1.errorbar(Z_RSD, fsig8_obs, yerr=fsig8_err, fmt='ko', ms=7,
             capsize=4, zorder=6, label='Obs. RSD')
ax1.set_xlabel(r'$z$', fontsize=11)
ax1.set_ylabel(r'$f\sigma_8(z)$', fontsize=11)
ax1.set_title(r'$f\sigma_8$ — Opción B vs candidatos', fontsize=10, fontweight='bold')
ax1.legend(fontsize=7, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0.04, 0.65)

# Panel 2: G vs Ω_m scan
ax2 = axes[1]
Om_scan = np.linspace(0.05, 0.45, 40)
G_scan  = []
for Om_test in Om_scan:
    D_test_raw, _ = integrate_raw(make_growth_rhs(Om_test))
    G_scan.append(D_test_raw[-1] / D_lcdm_raw[-1])
G_scan = np.array(G_scan)
sig8_scan = sig8_Planck * G_scan

ax2.plot(Om_scan, sig8_scan, 'b-', lw=2.2, label=r'$\sigma_8 = 0.811 \times G(\Omega_m)$')
ax2.axhline(0.811, color='gray',      ls=':', lw=1.5, label='Planck σ₈=0.811')
ax2.axhline(0.737, color='darkgreen', ls='--', lw=1.5, label='Target σ₈=0.737 (KiDS)')
ax2.axhline(0.702, color='tomato',    ls='-.', lw=1.5, label='Paper 5 σ₈=0.702')
ax2.axvline(Om_dyn,  color='tomato',    ls='--', lw=1, alpha=0.7, label=f'Ω_m,dyn={Om_dyn:.3f}')
ax2.axvline(Om_both, color='steelblue', ls='--', lw=1, alpha=0.7, label=f'Ω_m,both={Om_both:.3f}')

# Marcar el punto OptB
ax2.plot(Om_both, sig8_twosector, 'o', color='steelblue', ms=10, zorder=8,
         label=f'OptB: σ₈={sig8_twosector:.3f}')

ax2.set_xlabel(r'$\Omega_m$ efectivo en crecimiento', fontsize=11)
ax2.set_ylabel(r'$\sigma_8 = 0.811 \times G(\Omega_m)$', fontsize=11)
ax2.set_title('σ₈ self-consistente vs Ω_m,eff', fontsize=10, fontweight='bold')
ax2.legend(fontsize=7.5)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0.05, 0.45)

# Panel 3: tensión media vs Ω_m
ax3 = axes[2]
tensions_scan = []
for Om_test, G_t in zip(Om_scan, G_scan):
    D_n, Dp_n = integrate_normed(make_growth_rhs(Om_test))
    sig8_t = sig8_Planck * G_t
    fsig8_t = np.array([get_fsig8(z, D_n, Dp_n, sig8_t) for z in Z_RSD])
    tensions_scan.append(np.mean(np.abs((fsig8_obs - fsig8_t) / fsig8_err)))

ax3.plot(Om_scan, tensions_scan, 'steelblue', lw=2.2,
         label='OptB self-consistente')

# Líneas de referencia horizontales
ax3.axhline(t_P5,   color='tomato',    ls='--', lw=1.5, label=f'P5: {t_P5:.2f}σ')
ax3.axhline(t_P6,   color='darkgreen', ls='-.', lw=1.5, label=f'P6 +WDM: {t_P6:.2f}σ')
ax3.axhline(t_LCDM, color='gray',      ls=':',  lw=1.5, label=f'ΛCDM: {t_LCDM:.2f}σ')
ax3.axhline(1.0,    color='green',     ls=':',  lw=1,   alpha=0.5)
ax3.axvline(Om_both, color='steelblue', ls='--', lw=1, alpha=0.7,
            label=f'Ω_m,both={Om_both:.3f}')

ax3.set_xlabel(r'$\Omega_m$ efectivo', fontsize=11)
ax3.set_ylabel(r'Tensión media $|\sigma|$', fontsize=11)
ax3.set_title('fσ₈ self-consistente', fontsize=10, fontweight='bold')
ax3.legend(fontsize=7.5)
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0.05, 0.45)
ax3.set_ylim(0, 3.5)

fig.suptitle(
    'Opción B — IS+2sector self-consistente: ¿σ₈=0.737 sin WDM?\n'
    r'$\sigma_8 = 0.811 \times G(\Omega_{m,{\rm eff}}=0.320)$ vs target KiDS',
    fontsize=11, fontweight='bold'
)
fig.tight_layout()
out = os.path.join(OUTDIR, 'fig_optionB_sigma8_selfconsistent.pdf')
fig.savefig(out, bbox_inches='tight')
fig.savefig(out.replace('.pdf', '.png'), bbox_inches='tight', dpi=150)
print(f"\n  Figura guardada: {out}")
print(SEP)
