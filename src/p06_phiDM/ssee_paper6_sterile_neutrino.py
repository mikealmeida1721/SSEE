#!/usr/bin/env python3
"""
ssee_paper6_sterile_neutrino.py — SSEE Paper 6: φ-DM como Neutrino Estéril
===========================================================================
Vínculo: la masa del neutrino activo (Paper 4) y el φ-DM comparten
el mismo operador τ_Π. ¿Puede un neutrino estéril SSEE resolver σ8?

Método:
  1. Derivar la masa del neutrino estéril desde álgebra SSEE
  2. Calcular k_fs requerido para resolver σ8 (target: σ8_eff ≈ 0.76)
  3. Verificar si la masa SSEE da ese k_fs
  4. Calcular fσ8 resultante y comparar con datos
"""

import numpy as np
from scipy.optimize import brentq
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
KAL0 = beta + pi_

w0   = -Tr / Mv
wa   = -(pi_ + phi + phi) / Kv
Omm_dyn = 1.0 + w0
OmDE    = 1.0 - Omm_dyn
MIRA    = (3*phi + pi_) / 4
tau_Pi  = KAL0 / (3.0 * OmDE)

# Residuo acústico (Paper 4) → masa neutrino activo
R_acoustic   = 4 * KAL0 - 22             # 0.0856
Ob_h2_ssee   = (pi_ - phi) / (3*(phi + pi_)**2)   # 0.02242 = (π−φ)/H₀_SSEE (OP-1)
mnu_active   = R_acoustic * Ob_h2_ssee * 94.07 / tau_Pi  # 0.0824 eV

# φ-DM fracciones
Om_CDM   = Omm_dyn           # 0.160
Om_phiDM = (MIRA - 1)*Omm_dyn   # 0.160
h_ssee   = 0.6706
Om_phiDM_h2 = Om_phiDM * h_ssee**2

# Observacionales
sig8_Planck = 0.811
sig8_KiDS   = 0.737
sig8_err    = 0.020

print("=" * 68)
print("  SSEE Paper 6 — φ-DM como Neutrino Estéril")
print("=" * 68)
print(f"\n  Residuo acústico R = 4·KAL₀−22 = {R_acoustic:.6f}")
print(f"  Σm_ν (activos, Paper 4) = {mnu_active:.6f} eV")
print(f"  τ_Π·H₀ = {tau_Pi:.6f}")

# ── Función de transferencia WDM (Viel et al. 2005) ─────────────────────────
def T_WDM(k_hMpc, k_fs_hMpc, nu=1.12):
    alpha = 1.0 / k_fs_hMpc
    return (1.0 + (alpha * k_hMpc)**(2*nu))**(-5.0/nu)

def sigma8_twosector(k_fs_hMpc):
    """σ8 efectivo con CDM + φ-DM(WDM)."""
    T = T_WDM(0.125, k_fs_hMpc)  # k_σ8 = 1/8 h/Mpc
    return sig8_Planck * (0.5 + 0.5 * T)

# ── P1: k_fs requerido para resolver σ8 ─────────────────────────────────────
target_sig8 = sig8_KiDS  # 0.737 (centro del rango KiDS)
def residual_sig8(k_fs):
    return sigma8_twosector(k_fs) - target_sig8

# Buscar k_fs
k_fs_needed = brentq(residual_sig8, 0.01, 100.0, xtol=1e-6)
print(f"\n  ── P1: k_fs requerido para σ8 = {target_sig8} ──────────────")
print(f"  k_fs needed = {k_fs_needed:.4f} h/Mpc")
print(f"  λ_fs needed = {1/k_fs_needed:.4f} Mpc/h")

# ── P2: Masa del φ-DM para ese k_fs ──────────────────────────────────────────
# Fórmula WDM (Viel+2005, Dodelson-Widrow non-thermal production):
# k_fs [h/Mpc] ≈ 7.68 × (m_s/keV)^0.5 × (Ω_s h²/0.1)^0.5
# Para producción térmica (Bode+2001):
# k_fs [h/Mpc] ≈ 23.7 × (m_s/keV)^1.11 × (Ω_s h²/0.12)^0.15

def k_fs_from_mass_thermal(m_keV):
    return 23.7 * (m_keV)**1.11 * (Om_phiDM_h2/0.12)**0.15

def k_fs_from_mass_DW(m_keV):
    """Dodelson-Widrow (non-thermal)."""
    return 7.68 * (m_keV)**0.5 * (Om_phiDM_h2/0.1)**0.5

# Invertir para encontrar la masa requerida
m_needed_thermal = brentq(lambda m: k_fs_from_mass_thermal(m) - k_fs_needed,
                           1e-5, 100.0, xtol=1e-8)
m_needed_DW      = (k_fs_needed / (7.68 * (Om_phiDM_h2/0.1)**0.5))**2

print(f"\n  Masa φ-DM requerida:")
print(f"    Producción térmica:    m_φ = {m_needed_thermal*1000:.4f} meV  = {m_needed_thermal:.6f} keV")
print(f"    Dodelson-Widrow:       m_φ = {m_needed_DW*1000:.4f} meV  = {m_needed_DW:.6f} keV")

# ── P3: Candidatos algebraicos SSEE para m_φ ─────────────────────────────────
print(f"\n  ── P3: Candidatos algebraicos SSEE para m_φ ─────────────────")
print(f"  {'Expresión':40s}  {'Valor [meV]':>12}  {'Δ thermal [%]':>14}  {'Δ DW [%]':>10}")

# Base: mnu_active = 84 meV
candidates = {
    "Σm_ν":                              mnu_active,
    "Σm_ν × MIRA":                       mnu_active * MIRA,
    "Σm_ν × MIRA²":                      mnu_active * MIRA**2,
    "Σm_ν × KAL₀":                       mnu_active * KAL0,
    "Σm_ν × KAL₀²":                      mnu_active * KAL0**2,
    "Σm_ν × KAL₀ × MIRA":               mnu_active * KAL0 * MIRA,
    "Σm_ν × τ_Π²":                       mnu_active * tau_Pi**2,
    "Σm_ν × (φ+π)²":                     mnu_active * (phi+pi_)**2,
    "Σm_ν × (φ+π)³":                     mnu_active * (phi+pi_)**3,
    "Σm_ν × KAL₀³/MIRA":                mnu_active * KAL0**3 / MIRA,
    "R_acoustic × 94.07 eV":             R_acoustic * 94.07,    # eV
    "R² × 94.07 eV / Ω_b h²":           R_acoustic**2 * 94.07 / Ob_h2_ssee,
}

m_th_meV = m_needed_thermal * 1e6  # meV (keV → meV: ×1e6... wait keV = 1e6 meV)
# Actually: 1 keV = 1e6 meV, 1 eV = 1e3 meV
m_th_eV  = m_needed_thermal * 1e3  # eV
m_DW_eV  = m_needed_DW * 1e3

best_name, best_diff = "", 1e10
for name, val_eV in candidates.items():
    val_meV = val_eV * 1000  # eV → meV
    diff_th = abs(val_eV - m_th_eV) / m_th_eV * 100
    diff_DW = abs(val_eV - m_DW_eV) / m_DW_eV * 100
    star = "★" if diff_th < 5 or diff_DW < 5 else ("~" if diff_th < 15 or diff_DW < 15 else "")
    print(f"  {name:40s}  {val_meV:>12.4f}  {diff_th:>+14.2f}%  {diff_DW:>+10.2f}%  {star}")
    if diff_th < best_diff:
        best_diff = diff_th
        best_name = name

print(f"\n  Mejor candidato (thermal): {best_name}  (Δ = {best_diff:.2f}%)")
print(f"\n  Masa requerida:")
print(f"    Thermal: {m_th_eV*1000:.4f} meV = {m_th_eV:.6f} eV = {m_needed_thermal:.6f} keV")
print(f"    DW:      {m_DW_eV*1000:.4f} meV = {m_DW_eV:.6f} eV = {m_needed_DW:.6f} keV")

# ── P4: Predicciones observacionales del modelo con k_fs requerido ──────────
print(f"\n  ── P4: Predicciones del modelo φ-DM con k_fs = {k_fs_needed:.3f} h/Mpc")
print(f"  σ8 efectivo = {sigma8_twosector(k_fs_needed):.4f}  (target: {sig8_KiDS})")

# Rango de σ8 sin cobertura de incertidumbre de k_fs
sig8_lo = sigma8_twosector(k_fs_needed * 0.8)
sig8_hi = sigma8_twosector(k_fs_needed * 1.2)
tension_s8_lo = abs(sig8_lo - sig8_KiDS) / sig8_err
tension_s8_hi = abs(sig8_hi - sig8_KiDS) / sig8_err
print(f"  Rango σ8 [k_fs±20%] = [{sig8_lo:.4f}, {sig8_hi:.4f}]")
print(f"  Tensión vs KiDS = [{tension_s8_lo:.2f}σ, {tension_s8_hi:.2f}σ]")

# ── FIGURA ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel 1: σ8 vs k_fs
k_scan = np.logspace(-2, 1, 200)
sig8_scan = np.array([sigma8_twosector(k) for k in k_scan])

ax = axes[0]
ax.semilogx(k_scan, sig8_scan, 'b-', lw=2, label=r'$\sigma_8^{eff}(k_{fs})$')
ax.axhline(sig8_KiDS, color='r', ls='--', lw=1.5, label=rf'KiDS $\sigma_8={sig8_KiDS}$')
ax.axhline(sig8_Planck, color='g', ls=':', lw=1.5, label=rf'Planck $\sigma_8={sig8_Planck}$')
ax.fill_between([k_scan[0], k_scan[-1]],
                [sig8_KiDS-sig8_err, sig8_KiDS-sig8_err],
                [sig8_KiDS+sig8_err, sig8_KiDS+sig8_err],
                alpha=0.2, color='red', label=r'KiDS $1\sigma$')
ax.axvline(k_fs_needed, color='purple', ls='-', lw=2,
           label=rf'$k_{{fs}}^{{needed}}={k_fs_needed:.3f}$ h/Mpc')
ax.set_xlabel(r'$k_{fs}$ [h/Mpc]')
ax.set_ylabel(r'$\sigma_8^{eff}$')
ax.set_title('Two-Sector φ-DM:\n' + r'$\sigma_8$ vs Free-streaming Scale')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# Panel 2: masa requerida vs producción
m_scan_eV = np.logspace(-3, 3, 200)
k_fs_thermal = np.array([k_fs_from_mass_thermal(m/1000) for m in m_scan_eV])
k_fs_DW      = np.array([k_fs_from_mass_DW(m/1000) for m in m_scan_eV])

ax = axes[1]
ax.loglog(m_scan_eV, k_fs_thermal, 'b-', lw=2, label='Thermal production')
ax.loglog(m_scan_eV, k_fs_DW, 'r--', lw=2, label='Dodelson-Widrow')
ax.axhline(k_fs_needed, color='purple', ls='-', lw=2,
           label=rf'$k_{{fs}}^{{needed}}={k_fs_needed:.3f}$ h/Mpc')
# Marcar candidatos SSEE
colors_c = ['gold', 'cyan', 'lime', 'orange']
for i, (name, val_eV) in enumerate(list(candidates.items())[:4]):
    val_meV = val_eV * 1000
    ax.axvline(val_meV, color=colors_c[i], ls=':', lw=1.5,
               label=f'{name[:20]}={val_meV:.2f} meV', alpha=0.8)
ax.set_xlabel(r'$m_{\phi DM}$ [meV]')
ax.set_ylabel(r'$k_{fs}$ [h/Mpc]')
ax.set_title(r'φ-DM Mass → Free-streaming Scale')
ax.legend(fontsize=6, loc='lower right')
ax.grid(True, alpha=0.3)

fig.suptitle(
    'SSEE Paper 6: φ-DM Sterile Neutrino\n'
    rf'Σ$m_\nu^{{active}}$ = {mnu_active*1000:.2f} meV (Paper 4)  →  '
    rf'Required: $m_\phi$ = {m_th_eV*1000:.2f} meV (thermal) / {m_DW_eV*1000:.2f} meV (DW)',
    fontsize=10
)
fig.tight_layout()
out = os.path.join(OUTDIR, 'fig_paper6_sterile_neutrino.pdf')
fig.savefig(out, bbox_inches='tight')

print(f"\n  Figura: {out}")

# ── CONCLUSIÓN ───────────────────────────────────────────────────────────────
print(f"\n{'='*68}")
print(f"  VEREDICTO — φ-DM Neutrino Estéril")
print(f"{'='*68}")
print(f"  Σm_ν activos (Paper 4):  {mnu_active*1000:.4f} meV  = {mnu_active:.6f} eV")
print(f"  m_φ requerido (thermal): {m_th_eV*1000:.4f} meV  = {m_th_eV:.6f} eV")
print(f"  m_φ requerido (DW):      {m_DW_eV*1000:.4f} meV  = {m_DW_eV:.6f} eV")
ratio_th = m_th_eV / mnu_active
ratio_DW = m_DW_eV / mnu_active
print(f"\n  Ratio m_φ/Σm_ν (thermal) = {ratio_th:.6f}")
print(f"  Ratio m_φ/Σm_ν (DW)      = {ratio_DW:.6f}")
print(f"\n  Constantes SSEE cercanas al ratio:")
for name, val in {"MIRA": MIRA, "MIRA²": MIRA**2, "KAL₀/MIRA": KAL0/MIRA,
                   "(π-φ)": pi_-phi, "φ": phi, "τ_Π/KAL₀": tau_Pi/KAL0}.items():
    diff = abs(val - ratio_th)/ratio_th*100
    star = "★" if diff < 2 else ("~" if diff < 8 else "")
    print(f"    {name:20s} = {val:.6f}  Δ = {diff:.2f}%  {star}")
print(f"{'='*68}")
