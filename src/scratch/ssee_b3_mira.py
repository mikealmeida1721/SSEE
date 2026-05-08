"""
ssee_b3_mira.py — B3: MIRA physical mechanism
Algebraic duality Ω_m,dyn / Ω_m,CMB and IS connection.

Key results verified here:
  MIRA = (3φ+π)/4 = 1.99892...
  Ω_m,dyn  = (π−φ)/[2(π+φ)] = 0.16005
  Ω_m,CMB  = (π−φ)(3φ+π)/[8(π+φ)] = 0.31993
  MIRA     = Ω_m,CMB / Ω_m,dyn  (exact algebraic identity)
  τ_Π H₀  = KAL₀/(3 Ω_DE)      ≈ 2.191  (IS causal scale, Paper 1 App. A)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp

# ── SSEE algebraic constants ──────────────────────────────────────────────────
phi   = (1.0 + 5.0**0.5) / 2.0          # golden ratio  ≈ 1.61803
pi_   = np.pi                             # π
Omega = pi_ + phi                         # ≈ 4.75962
beta  = (pi_ + phi) / 2.0               # ≈ 2.37981
KAL0  = beta + pi_                        # ≈ 5.52140
P_sc  = Omega + phi                       # ≈ 6.37765
Kv    = phi + pi_ + Omega                # ≈ 9.51924
Tr    = 3.0 * (phi + beta)              # ≈ 11.9935
Mv    = phi + pi_ + Kv                   # ≈ 14.2788

w0    = -Tr / Mv                         # ≈ −0.83993
Om_dyn = 1.0 + w0                       # = 1 − |w0| = 0.16007  (≡ Ω_m,dyn)
Om_DE  = Tr / Mv                         # ≈ 0.83993

MIRA_genesis = (phi + beta) / 2.0       # Genesis 5.12 definition (= AURA/2)
MIRA_alg     = (3.0*phi + pi_) / 4.0   # algebraic form

# ── Verify exact algebraic identity ──────────────────────────────────────────
Om_dyn_exact = (pi_ - phi) / (2.0 * (pi_ + phi))  # from Tr/Mv expansion
Om_CMB_exact = (pi_ - phi) * (3.0*phi + pi_) / (8.0 * (pi_ + phi))
MIRA_ratio   = Om_CMB_exact / Om_dyn_exact          # should equal MIRA_alg

# IS relaxation time (Paper 1, Appendix A, eq:tau_IS)
tau_PI = KAL0 / (3.0 * Om_DE)   # τ_Π H₀ (dimensionless)

# Paper 4 formula for comparison
Om_CMB_paper4 = (pi_ - phi) / (pi_ + phi)          # ≈ 0.32001 (approx)

Planck_Om = 0.3153
Planck_sig = 0.0073

print("=" * 60)
print("B3 — MIRA Algebraic Duality Verification")
print("=" * 60)
print(f"\n  φ              = {phi:.10f}")
print(f"  π              = {pi_:.10f}")
print(f"  KAL₀ = β+π     = {KAL0:.10f}")
print(f"  Tr             = {Tr:.10f}")
print(f"  Mv             = {Mv:.10f}")
print()
print("── Sector 1: Dynamic (BAO/DE) ───────────────────────────")
print(f"  Ω_m,dyn = 1−|w₀|        = {Om_dyn:.10f}")
print(f"  Ω_m,dyn = (π−φ)/2(π+φ) = {Om_dyn_exact:.10f}")
print(f"  Difference               = {abs(Om_dyn - Om_dyn_exact):.2e}")
print()
print("── MIRA constant ────────────────────────────────────────")
print(f"  MIRA (Genesis 5.12)  (φ+β)/2     = {MIRA_genesis:.10f}")
print(f"  MIRA (algebraic)     (3φ+π)/4    = {MIRA_alg:.10f}")
print(f"  Difference                         = {abs(MIRA_genesis - MIRA_alg):.2e}")
print()
print("── Sector 2: CMB / observational ───────────────────────")
print(f"  Ω_m,CMB = Ω_m,dyn × MIRA          = {Om_dyn * MIRA_alg:.10f}")
print(f"  Ω_m,CMB = (π−φ)(3φ+π)/8(π+φ)     = {Om_CMB_exact:.10f}")
print(f"  Ω_m,CMB (Paper 4 approx)           = {Om_CMB_paper4:.10f}")
print(f"  Planck 2018                         = {Planck_Om:.4f} ± {Planck_sig:.4f}")
tension = abs(Om_CMB_exact - Planck_Om) / Planck_sig
print(f"  Tension                             = {tension:.2f}σ")
print()
print("── Algebraic identity check ─────────────────────────────")
print(f"  MIRA = Ω_m,CMB / Ω_m,dyn           = {MIRA_ratio:.10f}")
print(f"  MIRA (algebraic)                    = {MIRA_alg:.10f}")
print(f"  Identity holds to                   = {abs(MIRA_ratio - MIRA_alg):.2e}")
print()
print("── Israel-Stewart causal scale (Paper 1, App. A) ────────")
print(f"  τ_Π H₀ = KAL₀/(3Ω_DE)             = {tau_PI:.6f}")
print(f"  MIRA                                = {MIRA_alg:.6f}")
print(f"  τ_Π / MIRA                          = {tau_PI / MIRA_alg:.6f}  (≈ 1.10)")
print(f"  Both are ≈ 2 H₀⁻¹; differ at 9%")
print()
print("── Near-identity: 3φ+π ≈ 8 ─────────────────────────────")
val = 3.0*phi + pi_
print(f"  3φ+π   = {val:.8f}")
print(f"  8      = 8.00000000")
print(f"  Δ/8    = {abs(val - 8)/8 * 100:.4f}%  (<0.06%)")
print(f"  → Ω_m,CMB ≈ (π−φ)/(π+φ)  to 0.054%")
print("=" * 60)


# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 4.5))
gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.38)

# Panel (a): algebraic duality diagram (bar chart)
ax0 = fig.add_subplot(gs[0])
labels  = [r'$\Omega_{m,\mathrm{dyn}}$', r'$\Omega_{m,\mathrm{CMB}}$',
           r'$\Omega_m^{\,\mathrm{Planck}}$']
values  = [Om_dyn, Om_CMB_exact, Planck_Om]
errors  = [0.0,    0.0,          Planck_sig]
colors  = ['#4477BB', '#EE6644', '#AAAAAA']
x = np.arange(len(labels))
bars = ax0.bar(x, values, width=0.55, color=colors, zorder=3)
ax0.errorbar(x, values, yerr=errors, fmt='none', color='k', capsize=5, zorder=4)
ax0.axhline(Planck_Om, color='#AAAAAA', lw=1, ls='--', zorder=2)
ax0.set_xticks(x)
ax0.set_xticklabels(labels, fontsize=12)
ax0.set_ylabel(r'$\Omega_m$', fontsize=12)
ax0.set_ylim(0, 0.42)
ax0.set_title('(a) Two-sector matter density', fontsize=11)
ax0.grid(axis='y', lw=0.4, zorder=1)
# annotate MIRA arrow
ax0.annotate('',
    xy=(1, Om_CMB_exact), xytext=(0, Om_dyn),
    arrowprops=dict(arrowstyle='->', color='#4477BB', lw=1.5))
ax0.text(0.52, (Om_dyn + Om_CMB_exact)/2 + 0.003,
         r'$\times\,\mathrm{MIRA}$', fontsize=10, color='#4477BB')
ax0.text(1.0, Om_CMB_exact + 0.013,
         f'{Om_CMB_exact:.4f}\n({tension:.2f}σ)', ha='center', fontsize=9)

# Panel (b): IS scale comparison
ax1 = fig.add_subplot(gs[1])
names  = [r'$\mathrm{MIRA}=\frac{3\varphi+\pi}{4}$',
          r'$\tau_\Pi H_0 = \frac{\mathrm{KAL}_0}{3\Omega_\mathrm{DE}}$']
vals2  = [MIRA_alg, tau_PI]
cols2  = ['#EE6644', '#44AA77']
x2 = np.arange(len(names))
ax1.bar(x2, vals2, width=0.55, color=cols2, zorder=3)
ax1.axhline(2.0, color='k', lw=1, ls='--', zorder=2, label='$y=2$')
for xi, vi in zip(x2, vals2):
    ax1.text(xi, vi + 0.04, f'{vi:.4f}', ha='center', fontsize=10)
ax1.set_xticks(x2)
ax1.set_xticklabels(names, fontsize=11)
ax1.set_ylabel(r'Value (dimensionless)', fontsize=12)
ax1.set_ylim(0, 2.7)
ax1.set_title(r'(b) MIRA vs IS relaxation scale ($\approx 2\,H_0^{-1}$)', fontsize=11)
ax1.grid(axis='y', lw=0.4, zorder=1)
ax1.legend(fontsize=10)

fig.suptitle('SSEE B3 — MIRA Algebraic Duality\n'
             r'$\Omega_{m,\mathrm{CMB}}/\Omega_{m,\mathrm{dyn}} = '
             r'(3\varphi+\pi)/4 = \mathrm{MIRA}$ (exact)',
             fontsize=12, y=1.01)

plt.savefig('results/figures/fig_b3_mira.pdf', bbox_inches='tight')
plt.savefig('results/figures/fig_b3_mira.png', dpi=150, bbox_inches='tight')
print("\nFigures saved → results/figures/fig_b3_mira.{pdf,png}")
plt.close()
