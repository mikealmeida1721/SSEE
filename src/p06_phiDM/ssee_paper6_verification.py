#!/usr/bin/env python3
"""
ssee_paper6_verification.py  —  Verificación Completa del Modelo Dos-Sectores
==============================================================================
Física del modelo dos-sectores con partícula φ-DM canónica
(SOLAR²·KRYSTOS, reframe ω_m-directo 2026-06-19, OP-17 adoptado):

  Ω_CDM  = 0.160   → activo a TODAS las escalas
  Ω_φDM  = 0.14889 → Ω_m,CMB − Ω_m,dyn (diferencia, sin factor); activo k < k_fs

  Masa φ-DM por forward-prediction (CERO fiteo):
    m_φ = Σm_ν^active × (SOLAR²·KRYSTOS) = 0.06902 eV × 594.28 = 40.70 eV
    (eV × número puro = eV, dimensionalmente consistente)
    Ref: ssee_paper6_canonical_particle.py

  Amplitud TITULAR (free-streaming CLASS, NO growth-factor × ΛCDM):
    σ₈_eff = 0.747   — leído directamente de la corrida dos-sectores CLASS
                       con la temperatura de relic T_φ propia de la partícula
    S₈_eff = 0.747 × √(0.30889/0.300) = 0.758  → 0.04σ KiDS (RESUELVE S₈)
    k_fs   = 0.754 h/Mpc  (predicción falsable DESI Y3/Euclid)

  El growth-factor dos-sectores G_2s sólo fija la FORMA de f(z),
  no la amplitud (cross-check secundario: 0.811×G_2s).

  Tensión media fσ₈ (6 surveys RSD), two-sector σ₈=0.747 (CON free-streaming):
    0.93σ  (single-sector σ₈=0.8136 da 0.70σ; ΛCDM 0.73σ). El two-sector paga
    ~0.2σ en fσ₈ por bajar σ₈, a cambio de resolver S₈ (0.00σ vs 3.2σ).
    Corregido 2026-06-29: antes se usaba σ₈=0.8136 (single) para el two-sector → 0.70σ espurio.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', '..', 'results', 'figures')  # raíz results/figures (graphicspath del paper)
os.makedirs(OUTDIR, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'figure.dpi': 150})

# ── Constantes SSEE ──────────────────────────────────────────────────────────
phi   = (1 + 5**0.5) / 2
pi_   = np.pi
beta  = (pi_ + phi) / 2
Kv    = phi + pi_ + (pi_ + phi)
Tr    = 3 * (phi + beta)
Mv    = phi + pi_ + Kv
KAL0  = beta + pi_
w0    = -Tr / Mv
wa    = -(pi_ + phi + phi) / Kv
Omm_dyn = 1.0 + w0      # 0.160 (CDM sector)
OmDE    = 1.0 - Omm_dyn
MIRA    = (3*phi + pi_) / 4
AURA    = (3*phi + pi_) / 2
Om_DNAV = pi_ + phi          # 4.7596 (Ω = π+φ)

# ── Masa φ-DM: forward-prediction canónica (CERO fiteo) ───────────────────────
#   Cadena: R₂ = Ω_DNAV/(KAL·TRIAL) → Σm_ν^active = R₂ × 0.960318 eV
#           multiplicador = SOLAR²·KRYSTOS (número PURO; mecanismo g²·v)
#           m_φ = Σm_ν^active × multiplicador  (eV × adimensional = eV ✓)
#   Ref: ssee_paper6_canonical_particle.py (canónico SOLAR²·KRYSTOS 2026-06-19, OP-17)
R2         = Om_DNAV / (KAL0 * Tr)        # 0.071875
mnu_active = R2 * 0.960318                 # 0.06902 eV (0.960318 eV = input SM fijo)
# multiplicador canónico (reframe SOLAR²·KRYSTOS, 2026-06-19, OP-17 adoptado):
# SOLAR = BIAL+KAL = φ+2π (linaje radiativo); KRYSTOS_V = φ+π+Ω (padres, anclado por wₐ); forma g²·v
SOLAR      = beta + KAL0                    # 7.9012 (= φ+2π)
KRYSTOS_V  = phi + pi_ + Om_DNAV            # 9.5192 — padres {φ,π,Ω}, NO 2Ω (colapso)
multiplier = SOLAR**2 * KRYSTOS_V          # 594.28 (puro)
m_phi_eV   = mnu_active * multiplier       # 40.70 eV

# ── Densidad de materia CMB: ω_m-directo (sin factor; OP-8 disuelto) ───────────
H0_GLOBAL = 67.962
h_        = H0_GLOBAL / 100.0
n_s       = 1.0 - phi**-7                  # 0.96556
omega_b   = (pi_ - phi) / (3 * Om_DNAV**2) # 0.02242
omega_c   = KAL0 * omega_b * n_s           # 0.11951 (forward, Paper 1)
omega_nu  = mnu_active / 93.14             # 0.00074
omega_m   = omega_b + omega_c + omega_nu   # 0.14267
Om_m_CMB  = omega_m / h_**2                # 0.30889 (DERIVADO)

# ── Modelo dos-sectores ───────────────────────────────────────────────────────
Om_CDM   = Omm_dyn                  # 0.160
Om_phiDM = Om_m_CMB - Omm_dyn       # 0.14889 (DIFERENCIA, sin factor)
Om_total = Om_m_CMB                 # 0.30889 (ambos sectores = densidad CMB)

sig8_Planck = 0.811

print("=" * 70)
print("  SSEE Paper 6 — Verificación Completa Modelo Dos-Sectores (OptB)")
print("=" * 70)
print(f"\n  Parámetros del modelo:")
print(f"    Ω_CDM  = {Om_CDM:.6f}   (activo siempre)")
print(f"    Ω_φDM  = {Om_phiDM:.6f}   (activo para k < k_fs)")
print(f"    Ω_total= {Om_total:.6f}   (CMB y grandes escalas)")
print(f"    m_φ    = Σm_ν × (SOLAR²·KRYSTOS) = {m_phi_eV:.4f} eV  (escala energía φ-DM)")

# ── Ecuaciones de fondo ───────────────────────────────────────────────────────
def f_DE_raw(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
_norm = f_DE_raw(1.0)
def f_DE(a): return f_DE_raw(a) / _norm

# OmDE debe hacer PLANO cada fondo por separado (E²(a=1)=1). Un OmDE global
# (=1−Omm_dyn=0.840) sólo es plano para el fondo de 0.160; con Om_total=0.30889
# daba E²(hoy)=1.149 → H 7% alto → ahogaba el crecimiento (bug del fσ₈ two-sector
# 0.82σ, artefacto). Corregido 2026-06-22: cada E2 usa su propio Ω_DE=1−Ω_m.
OmDE_total = 1.0 - Om_total        # 0.69111 (fondo de clustering, plano)
OmDE_cdm   = 1.0 - Om_CDM          # 0.84000 (fondo minimal-CDM, plano)

def E2_total(a):
    """H²/H₀² con Ω_m=0.30889 (clustering); PLANO: Ω_DE=1−Ω_m."""
    return max(Om_total * a**(-3) + OmDE_total * f_DE(a), 1e-30)

def E2_CDM(a):
    """H²/H₀² con solo Ω_CDM=0.160 (baseline minimal); PLANO."""
    return max(Om_CDM * a**(-3) + OmDE_cdm * f_DE(a), 1e-30)

# ── Crecimiento dos-sectores (k < k_fs, ambos activos) ───────────────────────
def growth_rhs_twosector(lna, Y):
    """Ambos sectores activos: Ω_m,eff = Ω_total = 0.30889 (ω_m-directo)."""
    a = np.exp(lna)
    e2 = E2_total(a)
    h  = 0.5*(-3*Om_total*a**(-3) + OmDE_total*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a)/e2
    Omm_a = Om_total * a**(-3) / e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]

# ── Crecimiento un sector (Paper 5 baseline, solo CDM) ────────────────────────
def growth_rhs_single(lna, Y):
    """Solo CDM activo: Ω_m,eff = 0.160 (Paper 5)."""
    a = np.exp(lna)
    e2 = E2_CDM(a)
    h  = 0.5*(-3*Om_CDM*a**(-3) + OmDE_cdm*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a)/e2
    Omm_a = Om_CDM * a**(-3) / e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]

# ── ΛCDM referencia ───────────────────────────────────────────────────────────
def growth_lcdm(lna, Y, Om_lcdm=0.3153):
    a = np.exp(lna)
    e2 = max(Om_lcdm*a**(-3) + (1-Om_lcdm), 1e-30)
    h  = 0.5*(-3*Om_lcdm*a**(-3))/e2
    Omm_a = Om_lcdm*a**(-3)/e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm_a*D]

# ── Integración (IC idénticas — ratio D₁ da G directamente) ──────────────────
lna_pts = np.linspace(np.log(0.01), 0.0, 600)
IC = [0.01, 0.01]

sol_twosector = solve_ivp(growth_rhs_twosector, [lna_pts[0], 0.0], IC,
                          t_eval=lna_pts, rtol=1e-10, atol=1e-13)
sol_single    = solve_ivp(growth_rhs_single,    [lna_pts[0], 0.0], IC,
                          t_eval=lna_pts, rtol=1e-10, atol=1e-13)
sol_lcdm      = solve_ivp(growth_lcdm,          [lna_pts[0], 0.0], IC,
                          t_eval=lna_pts, rtol=1e-10, atol=1e-13)

# Factores de crecimiento (ratio de D₁ brutos con IC idénticas)
G_2s    = sol_twosector.y[0][-1] / sol_lcdm.y[0][-1]   # 1.0032 (fondo plano; el 0.979 era el bug fondo no-plano, corregido 2026-06-22)
G_single= sol_single.y[0][-1]    / sol_lcdm.y[0][-1]   # ≈ 0.866 (minimal-CDM Ω_m=0.160; NO es el single-sector canónico)

# ── Amplitud titular: free-streaming CLASS (NO growth-factor × ΛCDM) ──────────
#   σ_eff se lee DIRECTAMENTE de la corrida dos-sectores CLASS con la
#   temperatura de relic T_φ propia de la partícula (m_φ=40.70 eV).
#   Ref manuscript §fσ₈ L563-567 y ssee_paper6_canonical_particle.py.
sig8_eff    = 0.747                   # free-streaming CLASS canónico m_φ=40.70 (TITULAR)
sig8_base   = sig8_Planck * G_single  # 0.702 (minimal-CDM Ω_m=0.160; NO canónico — Paper 5 usa Ω_m=0.30889→σ₈=0.8136, fσ₈ 0.70σ)
sig8_growth = sig8_Planck * G_2s      # 0.8136 (growth-factor Ω_m=0.30889, amplitud RSD canónica)
S8_eff      = sig8_eff * np.sqrt(Om_total / 0.3)   # 0.758 (0.747×√(0.30889/0.3))

print(f"\n  Factores de crecimiento (IC idénticas, ratio D₁):")
print(f"    G_single (minimal-CDM Ω_m=0.160, NO canónico): {G_single:.4f}")
print(f"    G_2s     (Paper 6, Ω_m=0.30889): {G_2s:.4f}  (sólo forma de f(z))")
print(f"\n  σ₈ resultantes:")
print(f"    σ₈_base   (minimal-CDM Ω_m=0.160): {sig8_base:.4f}   (= 0.811 × {G_single:.4f}; NO canónico)")
print(f"    σ₈_eff    (Paper 6 TITULAR):    {sig8_eff:.4f}   (free-streaming CLASS)")
print(f"    σ₈_growth (cross-check 2-sec):  {sig8_growth:.4f}   (= 0.811 × {G_2s:.4f})")
print(f"    S₈_eff    (Paper 6 TITULAR):    {S8_eff:.4f}   (0.04σ KiDS; Ω_total={Om_total:.4f})")

# Normalizar D(z=0) = 1 para cálculo de f(z)
D_2s  = sol_twosector.y[0] / sol_twosector.y[0][-1]
Dp_2s = sol_twosector.y[1] / sol_twosector.y[0][-1]
D_1s  = sol_single.y[0]    / sol_single.y[0][-1]
Dp_1s = sol_single.y[1]    / sol_single.y[0][-1]
D_lcdm= sol_lcdm.y[0]      / sol_lcdm.y[0][-1]
Dp_lcdm= sol_lcdm.y[1]     / sol_lcdm.y[0][-1]

# ── fσ₈ en los 6 surveys RSD (set canónico Paper 5: refs Beutler2012,
#    Howlett2015, Alam2017, Hou2021 — valores idénticos a los citados) ──────────
#
#   AMPLITUD CORRECTA PARA fσ₈ = sig8_eff = 0.747 (CORREGIDO 2026-06-29).
#   El uso previo de sig8_growth=0.8136 era OPTIMISTA e incorrecto: descansaba en
#   suponer que RSD sondea sólo k≲0.1 h/Mpc, por debajo del free-streaming, donde
#   el φ-DM agruparía como frío. VERIFICADO FALSO con CLASS: la supresión muerde
#   DENTRO de la ventana top-hat de σ₈(R=8 Mpc/h) — σ₈ acumulado pierde 3.5% a
#   k=0.2 y 6.9% a k=0.3 (k_half=0.351 cae en medio de la ventana). Además
#   σ₈_cb (cold+baryon) = 0.747 = σ₈_total: no existe una σ₈ alta legítima para RSD.
#   Por convención fσ₈ se normaliza a σ₈(R=8), que para este modelo es 0.747.
#   ⇒ fσ₈ usa 0.747 (igual que S₈). Tensión media sube de 0.6958σ (espuria) a
#   ~0.95σ (honesta), consistente con el MCMC v2 (emulador CLASS, ssee_paper6_mcmc_v2.py).
#   Sigue <1σ y empata aprox. a ΛCDM (~0.73σ): el two-sector baja σ₈ → gana enorme
#   en S₈ (0.00σ vs 3.2σ) y paga poco en fσ₈. El trade-off es real, no un fit.
Z_RSD     = np.array([0.067, 0.150, 0.380, 0.510, 0.610, 1.480])
survey    = ['6dFGRS', 'SDSS MGS', 'BOSS DR12', 'BOSS DR12', 'BOSS DR12', 'eBOSS DR16']
fsig8_obs = np.array([0.423, 0.490, 0.497, 0.458, 0.436, 0.462])
fsig8_err = np.array([0.055, 0.145, 0.045, 0.038, 0.034, 0.045])

def get_fsig8(z, D_arr, Dp_arr, sig8_val):
    lna_t = np.log(1/(1+z))
    f_t = float(np.interp(lna_t, lna_pts, Dp_arr/D_arr))
    D_t = float(np.interp(lna_t, lna_pts, D_arr))
    return f_t * sig8_val * D_t

# Dos escenarios canónicos, MISMA forma de crecimiento (D_2s, Ω_m=0.30889):
#   - two-sector: amplitud σ₈=0.747 (CON free-streaming)  → titular Paper 6
#   - single-sector P5: amplitud σ₈=0.8136 (SIN free-streaming) → baseline Paper 5
# La ÚNICA diferencia es la amplitud — esa es la firma del free-streaming.
fsig8_twosector  = np.array([get_fsig8(z, D_2s,   Dp_2s,   sig8_eff)    for z in Z_RSD])
fsig8_single     = np.array([get_fsig8(z, D_2s,   Dp_2s,   sig8_growth) for z in Z_RSD])
fsig8_lcdm_pred  = np.array([get_fsig8(z, D_lcdm, Dp_lcdm, sig8_Planck) for z in Z_RSD])

t_twosector = (fsig8_obs - fsig8_twosector) / fsig8_err
t_single    = (fsig8_obs - fsig8_single)    / fsig8_err
t_lcdm      = (fsig8_obs - fsig8_lcdm_pred) / fsig8_err

print(f"\n  ── Resultados fσ₈ ───────────────────────────────────────────────")
print(f"\n  {'Survey':12s}  {'z':>5}  {'Obs':>6}  {'2-sector':>9}  {'σ':>7}  "
      f"{'1-sector':>9}  {'σ':>7}  {'ΛCDM':>7}  {'σ':>7}")
for i, z in enumerate(Z_RSD):
    print(f"  {survey[i]:12s}  {z:>5.3f}  {fsig8_obs[i]:>6.3f}  "
          f"{fsig8_twosector[i]:>9.3f}  {t_twosector[i]:>+7.2f}σ  "
          f"{fsig8_single[i]:>9.3f}  {t_single[i]:>+7.2f}σ  "
          f"{fsig8_lcdm_pred[i]:>7.3f}  {t_lcdm[i]:>+7.2f}σ")

mt_2s     = np.mean(np.abs(t_twosector))
mt_single = np.mean(np.abs(t_single))
mt_lcdm   = np.mean(np.abs(t_lcdm))
print(f"\n  Tensión media |σ|:")
print(f"    SSEE two-sector (σ₈=0.747, CON free-streaming):  {mt_2s:.4f}σ  ← titular")
print(f"    SSEE single-sector (σ₈=0.8136, SIN free-stream): {mt_single:.4f}σ  ← baseline P5")
print(f"    ΛCDM:                                            {mt_lcdm:.4f}σ")

# ── Resumen observacional ─────────────────────────────────────────────────────
print(f"\n  ── Resumen Observacional Completo ───────────────────────────────")
kids_s8 = 0.759; kids_s8_err = 0.024   # KiDS-1000 REAL (Asgari+2021); la predicción SSEE es 0.758 → 0.04σ. NO poner 0.758 (era el bug que forzaba 0.00σ)
des_s8  = 0.776; des_s8_err  = 0.017
kids_sig8 = 0.737; kids_sig8_err = 0.020
print(f"  σ₈_eff = {sig8_eff:.4f}  vs KiDS σ₈ {kids_sig8:.3f}±{kids_sig8_err:.3f}  → "
      f"{abs(sig8_eff-kids_sig8)/kids_sig8_err:.2f}σ")
print(f"  S₈_eff = {S8_eff:.4f}   vs KiDS S₈ {kids_s8:.3f}±{kids_s8_err:.3f}  → "
      f"{abs(S8_eff-kids_s8)/kids_s8_err:.2f}σ")
print(f"  S₈_eff = {S8_eff:.4f}   vs DES  S₈ {des_s8:.3f}±{des_s8_err:.3f}  → "
      f"{abs(S8_eff-des_s8)/des_s8_err:.2f}σ")
print(f"  H(z) χ²_r (sin cambio) = 1.861  [igual que SSEE base]")
print(f"  ΔBIC CMB = -32.9  [Paper 3 full plik N=2354, titular; cross-check plik_lite -23.93]")
print(f"  BAO χ²_r ≈ 0.01  [Paper 2, sin cambio]")
print(f"  Cluster χ²_r = 0.122  [Paper 2, sin cambio]")

# ── FIGURA ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Panel 1: fσ₈ comparación
ax = axes[0]
z_smooth = np.linspace(0.05, 0.65, 150)
fs8_2s     = [get_fsig8(z, D_2s, Dp_2s, sig8_eff)    for z in z_smooth]  # two-sector 0.747
fs8_single = [get_fsig8(z, D_2s, Dp_2s, sig8_growth) for z in z_smooth]  # single-sector 0.8136
fs8_lcdm   = [get_fsig8(z, D_lcdm, Dp_lcdm, sig8_Planck) for z in z_smooth]

ax.fill_between(z_smooth, np.array(fs8_2s)*0.92, np.array(fs8_2s)*1.08,
                alpha=0.15, color='blue', label='SSEE 2-sector ±8% envelope')
ax.plot(z_smooth, fs8_2s,   'b-',  lw=2.5,
        label=rf'SSEE two-sector ($\sigma_8=0.747$, free-stream) ($\bar\sigma={mt_2s:.2f}$)')
ax.plot(z_smooth, fs8_single, 'r--', lw=1.8,
        label=rf'SSEE single-sector ($\sigma_8=0.8136$, P5 baseline) ($\bar\sigma={mt_single:.2f}$)')
ax.plot(z_smooth, fs8_lcdm, 'g:',  lw=1.8,
        label=rf'$\Lambda$CDM ($\bar\sigma={mt_lcdm:.2f}$)')
ax.errorbar(Z_RSD, fsig8_obs, yerr=fsig8_err, fmt='ko', ms=7,
            capsize=4, zorder=6, label='RSD Surveys (6dF, SDSS, BOSS, eBOSS)')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$f\sigma_8(z)$')
ax.set_title(r'Growth Rate $f\sigma_8$: Two-Sector Verification')
ax.legend(fontsize=7.5)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.04, 0.65)

# Panel 2: Resumen de tensiones
ax = axes[1]
models   = ['SSEE\ntwo-sector', 'SSEE\nsingle (P5)', 'ΛCDM']
tensions = [mt_2s, mt_single, mt_lcdm]
colors   = ['royalblue', 'tomato', 'seagreen']
bars = ax.bar(models, tensions, color=colors, width=0.5, alpha=0.85,
              edgecolor='black', lw=1.2)
ax.axhline(2.0, color='orange', ls='--', lw=1.5, label='2σ threshold')
ax.axhline(1.0, color='green',  ls=':',  lw=1.5, label='1σ threshold')
for bar, t in zip(bars, tensions):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
            f'{t:.2f}σ', ha='center', va='bottom', fontsize=12, fontweight='bold')
ax.set_ylabel(r'Mean $|f\sigma_8|$ tension ($\sigma$)')
ax.set_title('Tension Summary\n' + r'$f\sigma_8$ vs 6 RSD Surveys')
ax.legend(fontsize=9)
ax.set_ylim(0, max(tensions)*1.25)
ax.grid(True, alpha=0.3, axis='y')

fig.suptitle(
    'SSEE Paper 6: Two-Sector Verification (free-streaming canónico, m_φ=40.70 eV)\n'
    rf'$m_\phi = \Sigma m_\nu^{{act}}\times(\mathrm{{SOLAR}}^2\cdot\mathrm{{KRYSTOS}}) = {m_phi_eV:.2f}$ eV,  '
    rf'$\sigma_8^{{eff}} = {sig8_eff:.4f}$,  '
    rf'$S_8^{{eff}} = {S8_eff:.4f}$ (0.04$\sigma$ KiDS)',
    fontsize=10.5
)
fig.tight_layout()
out = os.path.join(OUTDIR, 'fig_paper6_verification.pdf')
fig.savefig(out, bbox_inches='tight')
print(f"\n  Figura guardada: {out}")

print(f"\n{'='*70}")
print(f"  VEREDICTO FINAL — Modelo SSEE + dos-sectores (OptB)")
print(f"{'='*70}")
print(f"  G_2s:              {G_2s:.4f}  (sólo da la forma de f(z); Ω_m 0.320 vs 0.315)")
print(f"  σ₈_eff:            {sig8_eff:.4f}  (free-streaming CLASS, T_φ propia, sin parámetros libres)")
print(f"  S₈_eff:            {S8_eff:.4f}  (0.04σ KiDS — RESUELVE tensión S₈)")
print(f"  fσ₈ tensión media: {mt_2s:.4f}σ (two-sector, σ₈=0.747)  "
      f"[single-sector {mt_single:.4f}σ; ΛCDM {mt_lcdm:.4f}σ] {'(<1σ OK)' if mt_2s < 1.0 else ''}")
print(f"  σ₈ vs KiDS-1000:   {abs(sig8_eff-kids_sig8)/kids_sig8_err:.4f}σ")
print(f"  S₈ vs KiDS-1000:   {abs(S8_eff-kids_s8)/kids_s8_err:.4f}σ")
print(f"  Parámetros libres nuevos: 0  (m_φ = Σm_ν × (SOLAR²·KRYSTOS), forward-prediction)")
print(f"  Predicción falsable: k_fs = 0.754 h/Mpc (DESI Y3/Euclid 2026–2028)")
print(f"    m_φ = {m_phi_eV:.2f} eV es algebraico; φ-DM no tiene portal con el")
print(f"    Modelo Estándar → observable solo vía k_fs en P(k), no en KATRIN/PTOLEMY")
print(f"{'='*70}")
