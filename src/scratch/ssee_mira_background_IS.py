#!/usr/bin/env python3
"""
ssee_mira_background_IS.py
Derivación del factor MIRA desde background IS.

Pregunta: ¿La viscosidad bulk IS en el background modifica H(z) de tal manera
que un observador ΛCDM infiere Ω_m,fit ≈ MIRA × Ω_m,dyn = 0.3199?

Método:
  1. Ecuación de Friedmann IS modificada:
       E²(a) = Ω_m,dyn a⁻³ + Ω_DE f_DE(a) + Π̃(a)
  2. Ecuación de transporte IS (background):
       d(Π̃)/d(ln a) = −[Π̃/τ_Π H₀ + ζ̃_bg × Ω_DE f_DE(a) E(a)]
  3. Condición inicial: Π̃(a=1) = 0  (H₀ medido hoy incluye todo)
  4. Integración hacia atrás (a: 1 → 0.01)
  5. Distancias comoving D_M(z) = ∫ da/(a² H/H₀)
  6. Ajuste ΛCDM plano: minimiza χ² en D_M → Ω_m,fit
  7. MIRA_IS = Ω_m,fit / Ω_m,dyn  →  comparar con (3φ+π)/4

Parámetro libre: ζ̃_bg (viscosidad bulk background).
  - Si MIRA emerge con ζ̃_bg = KAL₀/3 (mismo que perturbaciones) → derivación exacta.
  - Si emerge con otro valor SSEE → nuevo resultado algebraico.
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize_scalar, brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', 'results', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'figure.dpi': 150})

# ── 1. CONSTANTES SSEE ────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi_  = np.pi
beta = (pi_ + phi) / 2
KAL0 = beta + pi_
P_sc = (pi_ + phi) + phi
Kv   = phi + pi_ + (pi_ + phi)
Tr   = 3 * (phi + beta)
Mv   = phi + pi_ + Kv

w0   = -Tr / Mv          # ≈ -0.8399
wa   = -P_sc / Kv        # ≈ -0.6699

Omm_dyn  = 1 + w0        # ≈ 0.1601
OmDE     = 1 - Omm_dyn   # ≈ 0.8399
MIRA_alg = (3*phi + pi_) / 4   # ≈ 1.9989
Omm_CMB  = Omm_dyn * MIRA_alg  # ≈ 0.3199

tau_Pi_H0  = KAL0 / (3.0 * OmDE)   # ≈ 2.191
zeta_pert  = KAL0 / 3.0             # perturbación: ζ̃ ≈ 1.841

print("=" * 68)
print("SSEE — Derivación MIRA desde background IS")
print("=" * 68)
print(f"\n  φ              = {phi:.8f}")
print(f"  w₀             = {w0:.8f}")
print(f"  Ω_m,dyn        = {Omm_dyn:.8f}")
print(f"  Ω_m,CMB (MIRA) = {Omm_CMB:.8f}   [target]")
print(f"  MIRA (alg.)    = {MIRA_alg:.8f}")
print(f"  τ_Π H₀         = {tau_Pi_H0:.8f}")
print(f"  ζ̃ (perturb.)   = {zeta_pert:.8f}")

# ── 2. f_DE(a) normalizada ────────────────────────────────────────────────────
def f_DE_raw(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

_fDE_norm = f_DE_raw(1.0)

def f_DE(a):
    return f_DE_raw(a) / _fDE_norm   # f_DE(1) = 1

# ── 3. ODE IS de background ───────────────────────────────────────────────────
def background_IS_rhs(lna, Y, zeta_bg):
    """
    State: Y = [Π̃]   Π̃ ≡ Π/(3H₀²)

    Friedmann: E²(a) = Ω_m,dyn a⁻³ + Ω_DE f_DE(a) + Π̃(a)
    IS transp: dΠ̃/d(ln a) = -[Π̃/τ_Π H₀ + ζ̃_bg × Ω_DE f_DE(a) × E(a)]

    Π̃ > 0 significa presión POSITIVA → más H(z) a z>0.
    Π̃ < 0 sería anti-compresión → menos H(z).
    """
    a   = np.exp(lna)
    Pit = Y[0]
    fde = f_DE(a)

    E2 = Omm_dyn * a**(-3) + OmDE * fde + Pit
    E2 = max(E2, 1e-20)
    E  = np.sqrt(E2)

    dPit_dlna = -(Pit / tau_Pi_H0 + zeta_bg * OmDE * fde * E)
    return [dPit_dlna]

def integrate_IS_background(zeta_bg, lna_ini=0.0, lna_fin=-5.0, n_pts=2000):
    """
    Integra IS de background desde a=1 (lna=0) hacia atrás hasta a=e^{-5}≈0.0067.
    IC: Π̃(1) = 0.
    Devuelve arrays (a, E, Pit).
    """
    sol = solve_ivp(
        lambda x, Y: background_IS_rhs(x, Y, zeta_bg),
        [lna_ini, lna_fin],
        [0.0],
        method='RK45',
        t_eval=np.linspace(lna_ini, lna_fin, n_pts),
        rtol=1e-10, atol=1e-13,
    )
    lna_arr = sol.t
    Pit_arr = sol.y[0]
    a_arr   = np.exp(lna_arr)
    fde_arr = np.array([f_DE(a) for a in a_arr])
    E2_arr  = Omm_dyn * a_arr**(-3) + OmDE * fde_arr + Pit_arr
    E2_arr  = np.maximum(E2_arr, 1e-20)
    E_arr   = np.sqrt(E2_arr)
    return a_arr, E_arr, Pit_arr

# ── 4. Distancias comóviles ───────────────────────────────────────────────────
def comoving_distance(z_target, a_arr, E_arr):
    """D_M(z) = c/H₀ ∫₀ᶻ dz'/E(z')  =  c/H₀ ∫_{a}^{1} da/(a² E(a))"""
    a_target = 1.0 / (1.0 + z_target)
    # Filtrar sólo a < 1 y a >= a_target
    mask = (a_arr <= 1.0) & (a_arr >= a_target)
    if mask.sum() < 2:
        return np.nan
    a_sel = a_arr[mask][::-1]   # orden creciente en a
    E_sel = E_arr[mask][::-1]
    integrand = 1.0 / (a_sel**2 * E_sel)
    return np.trapezoid(integrand, a_sel)   # unidades de c/H₀

# ── 5. Ajuste ΛCDM plano ──────────────────────────────────────────────────────
# ΛCDM plano: E_LCDM²(a) = Ω_m a⁻³ + (1-Ω_m)
# El observador ajusta Ω_m a las distancias D_M(z) del modelo IS.

def E_LCDM(a, Om_fit):
    return np.sqrt(Om_fit * a**(-3) + (1.0 - Om_fit))

def DM_LCDM(z, Om_fit):
    """D_M(z) ΛCDM plano = ∫₀ᶻ dz/E_LCDM"""
    a_arr_L = np.linspace(1.0/(1+z), 1.0, 500)
    integrand = 1.0 / (a_arr_L**2 * np.array([E_LCDM(a, Om_fit) for a in a_arr_L]))
    return np.trapezoid(integrand, a_arr_L)

# Redshifts DESI-like para el ajuste
Z_FIT = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.5, 2.0, 2.5])

def chi2_LCDM_fit(Om_fit, DM_IS_arr):
    """χ² entre D_M(IS) y D_M(ΛCDM plano) en los redshifts Z_FIT."""
    chi2 = 0.0
    for z, dm_IS in zip(Z_FIT, DM_IS_arr):
        dm_lcdm = DM_LCDM(z, Om_fit)
        chi2 += ((dm_IS - dm_lcdm) / dm_IS)**2
    return chi2

def fit_Omm(zeta_bg, verbose=False):
    """
    Para una ζ̃_bg dada:
    1. Integra IS background → H(a)
    2. Calcula D_M(z) a los redshifts DESI
    3. Ajusta ΛCDM plano → Ω_m,fit
    4. Devuelve (Ω_m,fit, MIRA_IS)
    """
    a_arr, E_arr, Pit_arr = integrate_IS_background(zeta_bg)

    DM_IS = np.array([comoving_distance(z, a_arr, E_arr) for z in Z_FIT])
    if np.any(np.isnan(DM_IS)):
        return np.nan, np.nan

    result = minimize_scalar(
        lambda Om: chi2_LCDM_fit(Om, DM_IS),
        bounds=(0.05, 0.70),
        method='bounded',
        options={'xatol': 1e-8}
    )
    Om_fit  = result.x
    MIRA_IS = Om_fit / Omm_dyn

    if verbose:
        print(f"\n  ζ̃_bg = {zeta_bg:.4f}")
        print(f"  Ω_m,fit = {Om_fit:.6f}   MIRA_IS = {MIRA_IS:.6f}"
              f"   Δ = {abs(MIRA_IS - MIRA_alg)/MIRA_alg*100:.2f}%")
        print(f"  Π̃ range: [{Pit_arr.min():.4f}, {Pit_arr.max():.4f}]")

    return Om_fit, MIRA_IS

# ── 6. SCAN ζ̃_bg ─────────────────────────────────────────────────────────────
print("\n" + "─" * 68)
print("Scan ζ̃_bg: efecto sobre Ω_m,fit y MIRA_IS")
print("─" * 68)
print(f"  {'ζ̃_bg':>10}  {'Ω_m,fit':>10}  {'MIRA_IS':>10}  {'Δ MIRA [%]':>12}  nota")

zeta_scan = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.2,
                      1.4, 1.5, 1.6, 1.7, 1.841, 2.0, 2.5, 3.0])

scan_results = []
for zb in zeta_scan:
    Om_fit, MIRA_IS = fit_Omm(zb)
    if np.isnan(Om_fit):
        print(f"  {zb:>10.4f}  {'NaN':>10}  {'NaN':>10}  {'NaN':>12}")
        continue
    delta_pct = (MIRA_IS - MIRA_alg) / MIRA_alg * 100
    nota = ""
    if abs(MIRA_IS - MIRA_alg) < 0.005:
        nota = "★ MIRA EXACTO"
    elif abs(MIRA_IS - MIRA_alg) < 0.02:
        nota = "≈ MIRA (2%)"
    elif abs(MIRA_IS - MIRA_alg) < 0.05:
        nota = "~ MIRA (5%)"
    # Marcar los valores SSEE conocidos
    if abs(zb - KAL0/3) < 0.005:
        nota += "  ← ζ̃ perturb."
    if abs(zb - tau_Pi_H0/3) < 0.005:
        nota += "  ← τ/3"
    scan_results.append((zb, Om_fit, MIRA_IS, delta_pct))
    print(f"  {zb:>10.4f}  {Om_fit:>10.6f}  {MIRA_IS:>10.6f}  {delta_pct:>+12.3f}  {nota}")

# ── 7. Encontrar ζ̃_bg exacto que reproduce MIRA ──────────────────────────────
print("\n" + "─" * 68)
print("Búsqueda exacta del ζ̃_bg que reproduce MIRA_alg")
print("─" * 68)

def MIRA_residual(zb):
    _, MIRA_IS = fit_Omm(zb)
    if np.isnan(MIRA_IS):
        return np.nan
    return MIRA_IS - MIRA_alg

# Buscar dónde MIRA_residual = 0 entre zeta_bg = 0 y 3
try:
    # Encontrar bracket donde cambia signo
    zb_vals = np.linspace(0.01, 3.0, 60)
    res_vals = np.array([MIRA_residual(z) for z in zb_vals])
    # Encontrar cambios de signo
    sign_changes = np.where(np.diff(np.sign(res_vals)))[0]

    if len(sign_changes) > 0:
        idx = sign_changes[0]
        zb_lo, zb_hi = zb_vals[idx], zb_vals[idx+1]
        zb_exact = brentq(MIRA_residual, zb_lo, zb_hi, xtol=1e-6)
        Om_exact, MIRA_exact = fit_Omm(zb_exact)

        print(f"\n  ζ̃_bg exacto que reproduce MIRA = {zb_exact:.6f}")
        print(f"  Ω_m,fit = {Om_exact:.6f}  MIRA_IS = {MIRA_exact:.6f}")
        print(f"  MIRA_alg = {MIRA_alg:.6f}  Δ = {abs(MIRA_exact-MIRA_alg):.2e}")

        # Comparar con constantes SSEE
        candidates = {
            "KAL₀/3 (ζ̃ perturb.)": KAL0/3,
            "τ_Π H₀/3":            tau_Pi_H0/3,
            "Ω_DE":                 OmDE,
            "KAL₀/(2π)":           KAL0/(2*pi_),
            "φ/2":                  phi/2,
            "β":                    beta,
            "KAL₀/Mv":             KAL0/Mv,
            "1/(φ+1)":             1/(phi+1),
            "π/Kv":                pi_/Kv,
        }
        print("\n  Comparación con constantes SSEE:")
        for name, val in candidates.items():
            diff = abs(zb_exact - val)
            star = "★" if diff < 0.01 else ("~" if diff < 0.05 else "")
            print(f"    {name:25s} = {val:.6f}   Δ = {diff:.6f}  {star}")
    else:
        print("\n  No se encontró cruce de signo en [0.01, 3.0].")
        print(f"  Rango de MIRA_IS: [{res_vals.min()+MIRA_alg:.4f}, {res_vals.max()+MIRA_alg:.4f}]")

except Exception as e:
    print(f"  Error en búsqueda: {e}")

# ── 8. Análisis IS con ζ̃_bg = KAL₀/3 (caso principal) ───────────────────────
print("\n" + "─" * 68)
print(f"Caso principal: ζ̃_bg = KAL₀/3 = {KAL0/3:.6f} (mismo que perturbaciones)")
print("─" * 68)

a_arr, E_arr, Pit_arr = integrate_IS_background(KAL0/3)
Om_main, MIRA_main = fit_Omm(KAL0/3, verbose=True)

# Mostrar tabla a redshifts DESI
print(f"\n  {'z':>6}  {'a':>6}  {'E_IS':>8}  {'E_LCDM(Ω_m,dyn)':>18}  {'ΔE/E [%]':>10}")
for z in Z_FIT:
    a = 1/(1+z)
    E_IS = float(np.interp(a, a_arr[::-1], E_arr[::-1]))
    E_std = np.sqrt(Omm_dyn * a**(-3) + OmDE * f_DE(a))
    print(f"  {z:>6.2f}  {a:>6.4f}  {E_IS:>8.5f}  {E_std:>18.5f}  "
          f"{(E_IS-E_std)/E_std*100:>+10.4f}")

# ── 9. SUMMARY ────────────────────────────────────────────────────────────────
print("\n" + "═" * 68)
print("SUMMARY — MIRA desde background IS")
print("═" * 68)
print(f"\n  Ω_m,dyn         = {Omm_dyn:.6f}  (SSEE algebraico)")
print(f"  Ω_m,CMB (MIRA)  = {Omm_CMB:.6f}  (target = MIRA × Ω_m,dyn)")
print(f"  MIRA (alg.)     = {MIRA_alg:.6f}  = (3φ+π)/4")
print(f"\n  Con ζ̃_bg = KAL₀/3 = {KAL0/3:.6f}:")
if not np.isnan(Om_main):
    delta = abs(MIRA_main - MIRA_alg) / MIRA_alg * 100
    print(f"    Ω_m,fit = {Om_main:.6f}   MIRA_IS = {MIRA_main:.6f}   Δ = {delta:.2f}%")
    if delta < 1:
        print(f"    ★ MIRA REPRODUCIDO por background IS (< 1% error)")
    elif delta < 5:
        print(f"    ~ MIRA aproximado (< 5%)")
    else:
        print(f"    ✗ MIRA no reproducido con este ζ̃_bg")

# ── 10. FIGURA ────────────────────────────────────────────────────────────────
if scan_results:
    zb_arr  = np.array([r[0] for r in scan_results])
    Om_arr  = np.array([r[1] for r in scan_results])
    MIRA_arr = np.array([r[2] for r in scan_results])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.plot(zb_arr, MIRA_arr, 'bo-', lw=2, ms=7, label=r'MIRA$_{\rm IS}$ (fitted)')
    ax.axhline(MIRA_alg, color='r', ls='--', lw=2,
               label=rf'MIRA$=(3\varphi+\pi)/4={MIRA_alg:.4f}$')
    ax.axvline(KAL0/3, color='green', ls=':', lw=1.5,
               label=rf'$\tilde{{\zeta}}_{{pert}}=KAL_0/3={KAL0/3:.3f}$')
    ax.axhline(1.0, color='gray', ls=':', lw=1)
    ax.set_xlabel(r'$\tilde{\zeta}_{\rm bg}$ (background IS viscosity)')
    ax.set_ylabel(r'MIRA$_{\rm IS} = \Omega_{m,\rm fit}/\Omega_{m,\rm dyn}$')
    ax.set_title('MIRA from Background IS\n(ΛCDM fit to IS-modified H(z))')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 3)

    ax = axes[1]
    # Evolución Π̃(a) para ζ̃_bg = KAL₀/3
    z_arr_plot = 1/a_arr - 1
    mask_plot = (z_arr_plot >= 0) & (z_arr_plot <= 10)
    ax.plot(z_arr_plot[mask_plot], Pit_arr[mask_plot], 'b-', lw=2,
            label=r'$\tilde{\Pi}(z)$ IS background')
    ax.axhline(0, color='k', lw=0.8, ls=':')
    ax.set_xlabel(r'redshift $z$')
    ax.set_ylabel(r'$\tilde{\Pi}(z) \equiv \Pi/(3H_0^2)$')
    ax.set_title(r'IS Bulk Pressure Background'
                 '\n' rf'($\tilde{{\zeta}}_{{bg}}=KAL_0/3={KAL0/3:.3f}$)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 10)

    fig.suptitle(
        r'SSEE Paper 5b: Background IS $\to$ MIRA Factor'
        '\n'
        rf'$\Omega_{{m,\rm dyn}}={Omm_dyn:.4f}$,  '
        rf'MIRA$_{{alg}}={MIRA_alg:.4f}$,  '
        rf'$\tau_\Pi H_0={tau_Pi_H0:.4f}$',
        fontsize=11
    )
    fig.tight_layout()
    out = os.path.join(OUTDIR, 'fig_paper5b_MIRA_background.pdf')
    fig.savefig(out, bbox_inches='tight')
    print(f"\n  Figura guardada: {out}")

print("\n" + "═" * 68)
print("Análisis background IS completo.")
print("═" * 68)
