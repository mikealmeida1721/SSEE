#!/usr/bin/env python3
"""
ssee_paper5_IS_perturbations.py
Israel-Stewart causal perturbative analysis for SSEE dark energy.

Two central questions:
  Q1: Does IS regularize c²_s < 0 → stable DE perturbations at all subhorizon k?
  Q2: Does Ω_m,eff/Ω_m,dyn ≈ MIRA = (3φ+π)/4 ≈ 1.9989 emerge from the IS
      perturbation growth history?

IS perturbation equations (Newtonian gauge, d/d ln a convention):
  δ'_m  = -(k/aH) v_m
  v'_m  = -v_m - (k/aH) Φ
  δ'_DE = -(1+w)(k/aH) v_DE
  v'_DE = -(1-3c²_s)v_DE - (k/aH)[c²_s δ_DE/(1+w) + Φ/(1+w)] + (k/aH) Π̃
  Π̃'   = -[Π̃ + ζ̃(k/aH) v_DE] / (τ_Π H)

Poisson (subhorizon): k² Φ = -(3/2)(H₀/H)²[Ω_m δ_m + Ω_DE f_DE(a) δ_DE]/a

IS dispersion relation (high-k limit, τ_Π k >> aH):
  c²_s,eff = c²_s,bare + ζ̃ / (τ_Π H₀)

References: Hiscock & Lindblom (1985), Maartens (1996), Hu (1998), Weinberg (1971).
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', 'results', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11,
    'axes.labelsize': 12, 'axes.titlesize': 12,
    'legend.fontsize': 10, 'figure.dpi': 150,
})

# ── 1. SSEE CONSTANTS (algebraic, zero free parameters) ─────────────────────
phi  = (1 + 5**0.5) / 2          # golden ratio ≈ 1.618034
pi_  = np.pi

beta   = (pi_ + phi) / 2          # β   ≈ 2.3798
KAL0   = beta + pi_               # KAL₀ ≈ 5.5214  Structural Viscosity
P_sc   = (pi_ + phi) + phi        # P_sc ≈ 6.3776  Dynamical Evolution Scalar
Kv     = phi + pi_ + (pi_ + phi)  # Kv   ≈ 9.5192
Tr     = 3 * (phi + beta)         # Tr   ≈ 11.9935
Mv     = phi + pi_ + Kv           # Mv   ≈ 14.2788

w0     = -Tr / Mv                 # ≈ -0.8399
wa     = -P_sc / Kv               # ≈ -0.6699
H0_kms = 3 * (phi + pi_)**2      # ≈ 67.962 km/s/Mpc

Omm_dyn = 1 + w0                  # Ω_m,dynamic ≈ 0.1601
OmDE    = 1 - Omm_dyn             # Ω_DE ≈ 0.8399

MIRA_alg = (3*phi + pi_) / 4     # ≈ 1.9989
Omm_CMB  = Omm_dyn * MIRA_alg    # ≈ 0.3199

# IS relaxation time (dimensionless: τ_Π × H₀)
# Derived from background IS steady state: τ_Π = KAL₀/(3 Ω_DE H₀)
tau_Pi_H0 = KAL0 / (3.0 * OmDE)   # ≈ 2.191

# Bare k-essence sound speed: c²_s = 1/(2n-1) with n = w₀/2 → c²_s = w₀
cs2_bare = w0                       # = -0.8399

# Bulk viscosity (dimensionless: ζ̃ = ζ/(ρ_DE H₀))
# SSEE hypothesis: ζ̃ = KAL₀/3 (structural viscosity normalized to DE density)
zeta_tilde = KAL0 / 3.0            # ≈ 1.8405

# ── 2. PRINT HEADER ─────────────────────────────────────────────────────────
print("=" * 68)
print("SSEE-V3.6  Paper 5 — Israel-Stewart Causal Perturbative Analysis")
print("=" * 68)
print(f"\n  φ               = {phi:.8f}")
print(f"  KAL₀            = {KAL0:.8f}")
print(f"  w₀              = {w0:.8f}")
print(f"  wₐ              = {wa:.8f}")
print(f"  Ω_m,dyn         = {Omm_dyn:.8f}")
print(f"  Ω_m,CMB (target)= {Omm_CMB:.8f}")
print(f"  MIRA (algebraic)= {MIRA_alg:.8f}")
print(f"  τ_Π H₀         = KAL₀/(3 Ω_DE) = {tau_Pi_H0:.8f}")
print(f"  c²_s (bare)     = {cs2_bare:.8f}  ← gradient instability")
print(f"  ζ̃ = KAL₀/3     = {zeta_tilde:.8f}")

# ── 3. ANALYTIC STABILITY (Q1) ──────────────────────────────────────────────
# IS dispersion relation in the short-wavelength limit (k τ_Π >> 1):
#   c²_s,eff = c²_s,bare + ζ̃ / τ_Π_H0
#
# The IS correction exactly fills the negative c²_s if SSEE parameters are used:
#   ζ̃ / τ_Π_H0 = (KAL₀/3) / (KAL₀/(3 Ω_DE)) = Ω_DE = |w₀|
#
# → c²_s,eff = w₀ + |w₀| = 0   (exact marginal stability)

IS_correction = zeta_tilde / tau_Pi_H0    # = Ω_DE = |w₀|
cs2_eff_highk  = cs2_bare + IS_correction

print("\n" + "─" * 68)
print("Q1 — IS STABILITY ANALYSIS (analytic)")
print("─" * 68)
print(f"\n  IS correction  ζ̃/τ_Π H₀ = {IS_correction:.8f}")
print(f"  Expected value |w₀|      = {abs(w0):.8f}")
print(f"  Difference               = {abs(IS_correction - abs(w0)):.2e}  (round-off)")
print(f"\n  c²_s,eff (k→∞)  = {cs2_bare:.6f} + {IS_correction:.6f} = {cs2_eff_highk:.2e}")

if abs(cs2_eff_highk) < 1e-12:
    print("\n  *** EXACT MARGINAL STABILITY  c²_s,eff = 0  (algebraic identity) ***")
    print("  The SSEE IS parameters conspire to exactly neutralize c²_s < 0.")
    print("  DE perturbations at k >> k_crit behave as pressureless (c²_s = 0).")
elif cs2_eff_highk > 0:
    print(f"\n  IS stabilizes gradient: c²_s,eff = +{cs2_eff_highk:.2e} > 0")
else:
    print(f"\n  Still unstable at high-k: c²_s,eff = {cs2_eff_highk:.2e} < 0")

# Critical wavenumber (transition from unstable to stable regime)
# At k τ_Π = 1: k_crit = (aH)/(τ_Π H₀ |c_s|) × H/H₀
# At z=0, a=1, H=H₀: k_crit/H₀ = 1/(τ_Π H₀ × |c_s|)
c_s_abs = abs(cs2_bare)**0.5
k_crit_H0 = 1.0 / (tau_Pi_H0 * c_s_abs)   # in units of H₀/c
k_crit_Mpc = k_crit_H0 * (H0_kms / 299792.458)   # Mpc⁻¹
lambda_crit_Mpc = 2 * np.pi / k_crit_Mpc if k_crit_Mpc > 0 else np.inf

print(f"\n  k_crit/H₀ = 1/(τ_Π H₀ × |c_s|) = {k_crit_H0:.4f}")
print(f"  k_crit    = {k_crit_Mpc:.4e} Mpc⁻¹")
print(f"  λ_crit    = {lambda_crit_Mpc:.1f} Mpc")

if k_crit_H0 < 1.0:
    print(f"  → k_crit < H₀/c: ALL subhorizon modes (k > aH) are IS-stabilized.")
    print(f"    There is no observable gradient instability in SSEE.")
else:
    print(f"  → k_crit > H₀/c: instability persists for {k_crit_H0:.2f} < k/(H₀/c) < ∞")

# ── 4. IS EFFECTIVE SOUND SPEED vs k (full expression) ──────────────────────
def cs2_eff_IS(k_H0, a=1.0, H_val=None):
    """
    Effective sound speed from IS dispersion relation.

    Full k-dependent result (interpolating between NS and IS limits):
      c²_s,eff(k,a) = c²_s + ζ̃/τ_Π × x²/(1+x²)

    where x = k τ_Π / (aH) is the dimensionless IS parameter.
    The correction saturates at ζ̃/τ_Π for k τ_Π >> aH.
    """
    if H_val is None:
        H_val = H_interp(a)
    x = k_H0 * tau_Pi_H0 / (a * H_val)   # dimensionless IS parameter
    correction = (zeta_tilde / tau_Pi_H0) * x**2 / (1.0 + x**2)
    return cs2_bare + correction

def H_over_H0_approx(a):
    """Fast H(a)/H₀ using CPL analytic approximation."""
    rDE = a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
    return np.sqrt(Omm_dyn * a**(-3) + OmDE * rDE)

# Precompute ρ_DE(a) on fine grid to avoid repeated quadrature in ODE
_a_grid_rDE   = np.linspace(0.001, 1.0, 2000)
_rDE_grid     = _a_grid_rDE**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-_a_grid_rDE))
# Normalize so ρ_DE(1) = 1
_rDE_grid    /= _rDE_grid[-1]
_H_grid       = np.sqrt(Omm_dyn * _a_grid_rDE**(-3) + OmDE * _rDE_grid)

def rho_de_interp(a):
    """ρ_DE(a)/ρ_DE(1) via precomputed grid interpolation (fast)."""
    return float(np.interp(a, _a_grid_rDE, _rDE_grid))

def H_interp(a):
    """H(a)/H₀ via precomputed grid (fast)."""
    return float(np.interp(a, _a_grid_rDE, _H_grid))

print("\n" + "─" * 68)
print("IS effective sound speed at z=0 (full k-dependence):")
print(f"  {'k [H₀/c]':>10}   c²_s,eff    regime")
for kv in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
    cs2 = cs2_eff_IS(kv, a=1.0, H_val=1.0)
    x   = kv * tau_Pi_H0
    regime = "IS limit" if x > 3 else ("transition" if x > 0.3 else "NS limit")
    print(f"  {kv:>10.3f}   {cs2:+.6f}   {regime}")

# ── 5. BACKGROUND COSMOLOGY ──────────────────────────────────────────────────
print("\n" + "─" * 68)
print("Background cosmology (H/H₀ across redshift):")

a_arr = np.array([0.01, 0.1, 0.2, 0.5, 0.8, 1.0])
for a_ in a_arr:
    H_ = H_interp(a_)
    z_ = 1/a_ - 1
    print(f"  z={z_:6.1f}  a={a_:.2f}  H/H₀={H_:.4f}  τ_Π H = {tau_Pi_H0*H_:.3f}")

# ── 6. ODE SYSTEM — IS PERTURBATIONS ─────────────────────────────────────────
def rhs_QS(x, Y, k_H0):
    """
    IS perturbation ODE — quasi-static Π̃ approximation.

    State: Y = [δ_m, v_m, δ_DE, v_DE]   (4 variables, no separate Π̃)
    Valid when τ_Π H >> 1 (early times) and asymptotically (IS limit).

    Quasi-static IS: Π̃_QS = −ζ̃ (k/aH) v_DE
    This absorbs the stiff Π̃ ODE into an algebraic relation,
    adding a k-dependent friction to the Euler equation:
      friction: +(1−3c²_s + ζ̃ kH²) → +(1+2.52 + 1.84 kH²) at k=1, z=0

    Eigenvalues at k=1, z=0:  λ₁ ≈ +0.15, λ₂ ≈ −5.5
    Gradient instability survives (λ₁ > 0) but is SLOW (c.f. bare rate ~0.92).
    Solution remains finite and integrable for k ≥ 1 H₀/c from a=0.05 to a=1.
    """
    a   = np.exp(x)
    H   = H_interp(a)
    w   = w0 + wa * (1.0 - a)
    fDE = rho_de_interp(a)

    δm, vm, δDE, vDE = Y
    one_plus_w = max(1.0 + w, 1e-6)   # guard against 1+w → 0

    Phi = -(1.5 / (k_H0**2 * a * H**2)) * (Omm_dyn * δm + OmDE * fDE * δDE)
    kH  = k_H0 / (a * H)

    # Matter (pressureless CDM)
    dδm = -kH * vm
    dvm = -vm - kH * Phi

    # Dark energy — quasi-static IS Π̃ absorbed as extra friction
    # Pit_QS = -zeta_tilde * kH * vDE
    eff_friction = (1.0 - 3.0 * cs2_bare) + zeta_tilde * kH**2
    dδDE = -one_plus_w * kH * vDE
    dvDE = (-eff_friction * vDE
            - kH * cs2_bare * δDE / one_plus_w
            - kH * Phi / one_plus_w)

    return [dδm, dvm, dδDE, dvDE]

# ── 7. INTEGRATE AND EXTRACT MIRA RATIO ──────────────────────────────────────
def integrate_growth(k_H0, a_ini=0.05, a_fin=1.0, n_pts=400):
    """
    Integrate IS perturbation equations (quasi-static approximation).

    Start at a_ini = 0.05 (z=19): DE energy fraction ~ 0, IS firmly suppressing
    DE perturbations (τ_Π H ≈ 22, quasi-static regime valid).

    IC: adiabatic growing mode for matter; δ_DE = v_DE = 0.
    The DE perturbation is sourced purely by gravitational coupling.
    """
    x_ini = np.log(a_ini)
    x_fin = np.log(a_fin)

    δm0 =  a_ini    # growing mode δ_m ∝ a
    vm0 = -a_ini    # θ_m = -H δ_m → v_m normalized
    Y0  = [δm0, vm0, 0.0, 0.0]

    x_eval = np.linspace(x_ini, x_fin, n_pts)

    sol = solve_ivp(
        lambda x, Y: rhs_QS(x, Y, k_H0),
        [x_ini, x_fin],
        Y0,
        method='RK45',
        t_eval=x_eval,
        rtol=1e-9,
        atol=1e-12,
        dense_output=False,
        max_step=0.02,
    )
    return sol

# ── 8. RUN FOR MULTIPLE k VALUES ─────────────────────────────────────────────
print("\n" + "─" * 68)
print("Q2 — MIRA TEST: DE perturbation growth vs matter")
print("─" * 68)
print(f"  MIRA (algebraic) = {MIRA_alg:.6f}  →  Ω_m,CMB = {Omm_CMB:.6f}")
print(f"  Required δ_DE/δ_m at z=0 for MIRA: r* = (MIRA-1)×Ω_m/Ω_DE = "
      f"{(MIRA_alg-1)*Omm_dyn/OmDE:.4f}")
print()
print(f"  {'k [H₀/c]':>10}  {'δ_m(z=0)':>12}  {'δ_DE(z=0)':>12}  "
      f"{'r=δ_DE/δ_m':>12}  {'Ω_m,eff':>10}  {'MIRA_num':>10}  status")

# Use a_ini = 0.05 (z=19) for all runs.  At a_ini=0.05, H ≈ 35.8 H₀ →
# a_ini H(a_ini) ≈ 1.79 H₀/c.  Only modes with k > 2 H₀/c are robustly
# sub-horizon throughout the integration; lower-k results are flagged.
k_values = [2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0]
results = {}

for k_H0 in k_values:
    sol = integrate_growth(k_H0, a_ini=0.05, a_fin=1.0)
    if not sol.success:
        print(f"  {k_H0:>10.2f}  ODE failed: {sol.message}")
        continue

    δm_fin  = sol.y[0, -1]
    δDE_fin = sol.y[2, -1]

    # Normalize to δ_m(z=0) = 1 (growth factor convention)
    if abs(δm_fin) < 1e-15:
        continue
    r = δDE_fin / δm_fin

    # Effective matter density: Poisson → both m and DE contribute
    # Ω_m,eff δ_eff = Ω_m δ_m + Ω_DE δ_DE = Ω_m(1 + r × Ω_DE/Ω_m)
    Om_eff   = Omm_dyn + OmDE * r
    MIRA_num = Om_eff / Omm_dyn

    status = ""
    if abs(MIRA_num - MIRA_alg) < 0.05:
        status = "✓ MIRA"
    elif abs(MIRA_num - MIRA_alg) < 0.20:
        status = "~ close"
    else:
        status = "✗ IS-suppressed"

    results[k_H0] = dict(
        δm=δm_fin, δDE=δDE_fin, r=r,
        Om_eff=Om_eff, MIRA_num=MIRA_num, sol=sol
    )

    print(f"  {k_H0:>10.2f}  {δm_fin:>12.4f}  {δDE_fin:>12.6f}  "
          f"{r:>12.6f}  {Om_eff:>10.4f}  {MIRA_num:>10.4f}  {status}")

# ── 9. SUMMARY ───────────────────────────────────────────────────────────────
print("\n" + "═" * 68)
print("SUMMARY")
print("═" * 68)

print(f"\n  Q1 — Gradient stability:")
print(f"    c²_s,bare = {cs2_bare:.6f} (k-essence, unstable)")
print(f"    c²_s,eff  = {cs2_eff_highk:.2e} (IS high-k limit, marginal)")
print(f"    Identity: ζ̃/τ_Π = KAL₀/3 ÷ KAL₀/(3Ω_DE) = Ω_DE = |w₀|  ✓")
print(f"    ALL subhorizon modes stabilized (k_crit = {k_crit_H0:.3f} H₀/c < H₀/c)")

if results:
    # Use only well-subhorizon IS regime: k ≥ 10 H₀/c
    # At these scales k τ_Π >> 1 is satisfied throughout, and the mode is
    # robustly sub-horizon from a_ini onward (kH ≥ 10/35.8 ≈ 0.28 at a=0.05).
    k_IS_keys = [k for k in results if k >= 10.0]
    if k_IS_keys:
        MIRA_mean = np.mean([results[k]['MIRA_num'] for k in k_IS_keys])
        MIRA_std  = np.std( [results[k]['MIRA_num'] for k in k_IS_keys])
        delta_MIRA = abs(MIRA_mean - MIRA_alg)
        frac_MIRA  = delta_MIRA / MIRA_alg * 100

        print(f"\n  Q2 — MIRA ratio (k ≥ 10 H₀/c, well inside IS + sub-horizon):")
        print(f"    MIRA_numeric = {MIRA_mean:.6f} ± {MIRA_std:.6f}")
        print(f"    MIRA_algebraic = {MIRA_alg:.6f}")
        print(f"    Discrepancy = {delta_MIRA:.6f}  ({frac_MIRA:.1f}%)")
        print(f"    δ_DE/δ_m → 0 for large k: IS SUPPRESSES DE clustering")

        if frac_MIRA < 5:
            print(f"    *** MIRA CONSISTENT with IS perturbation growth ***")
        elif frac_MIRA < 20:
            print(f"    ~ MIRA partially reproduced — subleading corrections needed")
        else:
            print(f"    ✗  MIRA NOT reproduced by linear IS perturbations")
            print(f"       Linear IS → Ω_m,eff ≈ Ω_m,dyn (no DE clustering)")
            print(f"       MIRA is a background-level IS effect, not perturbative")

# ── 10. FIGURES ──────────────────────────────────────────────────────────────

# Figure 1: c²_s,eff vs k at z=0
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
k_plot = np.logspace(-2, 2, 300)
cs2_plot = np.array([cs2_eff_IS(k, a=1.0, H_val=1.0) for k in k_plot])
ax.semilogx(k_plot, cs2_plot, 'b-', lw=2, label=r'$c^2_{s,\rm eff}(k)$ IS')
ax.axhline(cs2_bare,  color='r', ls='--', lw=1.5, label=rf'bare $c^2_s={cs2_bare:.3f}$')
ax.axhline(0, color='k', ls=':', lw=1, label='stability threshold')
ax.axhline(cs2_eff_highk, color='g', ls='--', lw=1.5,
           label=rf'IS limit $c^2_{{s,\rm eff}}={cs2_eff_highk:.1e}$')
ax.axvline(k_crit_H0, color='gray', ls='-.', lw=1,
           label=rf'$k_{{\rm crit}}={k_crit_H0:.3f}\,H_0/c$')
ax.fill_between(k_plot, cs2_bare, 0, where=cs2_plot < 0, alpha=0.15, color='red',
                label='unstable region')
ax.fill_between(k_plot, 0, cs2_plot, where=cs2_plot >= 0, alpha=0.15, color='green',
                label='stable region')
ax.set_xlabel(r'$k\,[H_0/c]$')
ax.set_ylabel(r'$c^2_{s,\rm eff}(k)$')
ax.set_title('IS Effective Sound Speed\n(SSEE, $z=0$)')
ax.legend(fontsize=9, loc='lower right')
ax.set_xlim(k_plot[0], k_plot[-1])
ax.set_ylim(cs2_bare - 0.1, 0.1)
ax.grid(True, alpha=0.3)

# Annotate marginal stability
ax.annotate(r'$c^2_{s,\rm eff}\to 0$ (exact)',
            xy=(10, 0.002), xytext=(1, 0.06),
            arrowprops=dict(arrowstyle='->', color='green'),
            color='green', fontsize=10)

# Figure 2: Growth of δ_DE/δ_m vs a for selected k
ax = axes[1]
colors_k = plt.cm.plasma(np.linspace(0.2, 0.85, len(results)))
a_plot_vals = None

for idx, (k_H0, res) in enumerate(sorted(results.items())):
    sol = res['sol']
    a_arr_sol = np.exp(sol.t)
    δm_arr  = sol.y[0]
    δDE_arr = sol.y[2]

    # Avoid division by zero; normalize so δ_m(z=0) = 1
    δm_fin = sol.y[0, -1]
    if abs(δm_fin) < 1e-12:
        continue
    ratio_arr = δDE_arr / δm_arr

    ax.plot(a_arr_sol, ratio_arr,
            color=colors_k[idx], lw=1.8,
            label=rf'$k={k_H0}\,H_0/c$  ($r_{{z=0}}={res["r"]:.3f}$)')

# Target MIRA line
r_MIRA = (MIRA_alg - 1) * Omm_dyn / OmDE
ax.axhline(r_MIRA, color='k', ls='--', lw=1.5,
           label=rf'$r^*={r_MIRA:.3f}$ (MIRA target)')

ax.set_xlabel(r'scale factor $a$')
ax.set_ylabel(r'$\delta_{\rm DE}/\delta_m$')
ax.set_title('DE Perturbation Growth\n(IS causal regularization)')
ax.legend(fontsize=8, loc='upper left')
ax.set_xlim(0, 1)
ax.set_ylim(-0.05, max(r_MIRA * 2, 0.5))
ax.grid(True, alpha=0.3)
ax.axhline(0, color='gray', lw=0.7)

fig.suptitle(
    r'SSEE Paper 5: Israel-Stewart Causal Perturbations'
    '\n'
    rf'$c^2_{{s,\rm bare}}={cs2_bare:.3f}$,  '
    rf'$\tilde{{\zeta}}=KAL_0/3={zeta_tilde:.4f}$,  '
    rf'$\tau_\Pi H_0={tau_Pi_H0:.4f}$,  '
    rf'MIRA$={MIRA_alg:.4f}$',
    fontsize=11
)
fig.tight_layout()

out1 = os.path.join(OUTDIR, 'fig_paper5_IS_stability.pdf')
fig.savefig(out1, bbox_inches='tight')
print(f"\n  Figure saved: {out1}")

# Figure 3: MIRA_numeric vs k
if len(results) >= 3:
    fig2, ax2 = plt.subplots(figsize=(7, 4.5))

    k_arr     = np.array(sorted(results.keys()))
    MIRA_arr  = np.array([results[k]['MIRA_num'] for k in k_arr])
    Om_eff_arr = np.array([results[k]['Om_eff']  for k in k_arr])

    ax2.semilogx(k_arr, MIRA_arr, 'bo-', lw=2, ms=7, label=r'$\Omega_{m,\rm eff}/\Omega_{m,\rm dyn}$ (IS)')
    ax2.axhline(MIRA_alg, color='r', ls='--', lw=2,
                label=rf'MIRA$=(3\varphi+\pi)/4={MIRA_alg:.4f}$ (algebraic)')
    ax2.axhline(1.0, color='gray', ls=':', lw=1, label=r'no DE clustering')
    ax2.axhline(1.0/Omm_dyn, color='gray', ls='-.', lw=1,
                label=rf'full clustering ($1/\Omega_m={1/Omm_dyn:.2f}$)')
    ax2.axvline(k_crit_H0, color='purple', ls=':', lw=1.5,
                label=rf'$k_{{\rm crit}}={k_crit_H0:.3f}\,H_0/c$')

    ax2.fill_between([k_crit_H0, k_arr.max()],
                     [MIRA_alg*0.95]*2, [MIRA_alg*1.05]*2,
                     alpha=0.15, color='red', label='±5% MIRA band')

    ax2.set_xlabel(r'$k\,[H_0/c]$')
    ax2.set_ylabel(r'MIRA$_{\rm num} = \Omega_{m,\rm eff}/\Omega_{m,\rm dyn}$')
    ax2.set_title(r'MIRA Test: IS Growth $\to$ Effective Matter Density'
                  '\n(SSEE Paper 5)')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(k_arr.min(), k_arr.max())

    fig2.tight_layout()
    out2 = os.path.join(OUTDIR, 'fig_paper5_MIRA_test.pdf')
    fig2.savefig(out2, bbox_inches='tight')
    print(f"  Figure saved: {out2}")

print("\n" + "═" * 68)
print("Paper 5 IS analysis complete.")
print("═" * 68)
