"""
SSEE UV Completion — Vainshtein Scale Analysis

Goals:
  1. Verify the Vainshtein radius formula from Paper 8:
       r_V = sqrt(KAL² |βc| × r_S × ℓ_Pl)
     and show whether the claimed r_V(Sun) ≈ 6 km is correct.

  2. Scan M to find what strong-coupling cutoff gives cosmological
     Vainshtein screening (r_V ∈ [50, 300] Mpc).

  3. Check whether any algebraic combination of φ and π gives such M.

  4. Provide a definitive verdict on the UV completion route for Paper 9.

All intermediate results are printed before use (verification mandate).
"""

import numpy as np

# ── Algebraic constants ────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2           # 1.61803
pi_  = np.pi                       # 3.14159
Omega = phi + pi_                  # 4.75962 — OMEGA_DNAV
KAL   = (phi + pi_) / 2 + pi_     # 5.52141 — Structural Viscosity
AURA  = (3*phi + pi_) / 2         # 3.99785 — βc = -AURA
beta_c = -AURA                     # EFT kineticity coupling

H0_alg = 3 * (phi + pi_)**2       # 67.960 km/s/Mpc
f_screen = (pi_ - phi) / Omega**2 # 0.06730 — algebraic Hubble correction

print("=" * 70)
print("SSEE UV Completion — Vainshtein Scale Analysis")
print("=" * 70)
print(f"\n  φ          = {phi:.6f}")
print(f"  π          = {pi_:.6f}")
print(f"  KAL        = {KAL:.6f}   (IR kinetic coefficient⁻¹)")
print(f"  |βc|=AURA  = {AURA:.6f}   (EFT coupling)")
print(f"  KAL² |βc|  = {KAL**2 * AURA:.4f}")
print(f"  √(KAL²|βc|) = {np.sqrt(KAL**2 * AURA):.4f}")
print(f"  f_screen   = (π-φ)/Ω² = {f_screen:.6f}")
print(f"  H0_alg     = {H0_alg:.4f} km/s/Mpc")

# ── Physical constants (SI) ────────────────────────────────────────────────────
ell_Pl  = 1.61623e-35   # m  — Planck length  sqrt(ħG/c³)
M_Pl_kg = 2.176e-8      # kg — Planck mass    sqrt(ħc/G)
M_Pl_GeV = 2.44e18      # GeV
G_SI    = 6.674e-11     # m³ kg⁻¹ s⁻²
c_SI    = 2.998e8       # m/s
H0_SI   = 67.96e3 / 3.0856e22  # s⁻¹ (67.96 km/s/Mpc → s⁻¹)
Mpc_m   = 3.0856e22     # m per Mpc

# Schwarzschild radii
r_S_sun     = 2.953e3        # m  (2.953 km)
M_sun_kg    = 1.989e30       # kg
r_S_cluster = 2 * G_SI * (1e15 * M_sun_kg) / c_SI**2   # 10^15 M_sun cluster

# Hubble radius
r_H = c_SI / H0_SI           # m

print(f"\n  ℓ_Pl          = {ell_Pl:.4e} m")
print(f"  M_Pl          = {M_Pl_GeV:.3e} GeV = {M_Pl_kg:.3e} kg")
print(f"  r_S(Sol)      = {r_S_sun:.3e} m = {r_S_sun/1e3:.3f} km")
print(f"  r_S(10¹⁵M☉)  = {r_S_cluster:.3e} m = {r_S_cluster/Mpc_m:.3f} Mpc")
print(f"  r_H (Hubble)  = {r_H:.3e} m = {r_H/Mpc_m:.0f} Mpc")

# ══════════════════════════════════════════════════════════════════════════════
# PART 1 — Verify Paper 8 formula: r_V = √(KAL² |βc|) × √(r_S × ℓ_Pl)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─"*70)
print("PART 1 — Paper 8 formula (M = M_Pl): r_V = √(KAL²|βc| × r_S × ℓ_Pl)")
print("─"*70)

P8_factor = np.sqrt(KAL**2 * AURA)   # = √(KAL² |βc|)

r_V_sun_P8     = P8_factor * np.sqrt(r_S_sun    * ell_Pl)
r_V_cluster_P8 = P8_factor * np.sqrt(r_S_cluster * ell_Pl)

print(f"\n  P8 pre-factor = √(KAL²|βc|) = {P8_factor:.4f}")
print(f"  √(r_S_sun × ℓ_Pl)     = {np.sqrt(r_S_sun*ell_Pl):.3e} m")
print(f"\n  r_V(Sol, Paper 8)      = {r_V_sun_P8:.3e} m = {r_V_sun_P8*1e15:.3f} fm")
print(f"  Paper 8 claims:          ≈ 6 km = 6000 m")
print(f"  Discrepancy factor:      {6000/r_V_sun_P8:.3e}  ← FORMULA ERROR in Paper 8")
print(f"\n  r_V(10¹⁵M☉, Paper 8)  = {r_V_cluster_P8:.3e} m = {r_V_cluster_P8*1e15:.3f} fm")

print(f"""
  ⚠️  DIAGNOSIS: The formula r_V = √(KAL²|βc| × r_S × ℓ_Pl) gives ~2 fm for
  the Sun, not 6 km. The Paper 8 claim is off by ~2.5×10¹⁸ — consistent with
  a missing factor of M_Pl_GeV in the unit conversion (M_Pl ≈ 2.44×10¹⁸ GeV).
  The CONCLUSION of Paper 8 (Vainshtein inoperative) is still correct — in
  fact MORE true than stated. The numerical value needs correction.
""")

# ══════════════════════════════════════════════════════════════════════════════
# PART 2 — General M dependence:  r_V(M) = r_V(M_Pl) × (M_Pl/M)
# From K(X) = X/KAL + X²/M⁴, the Vainshtein crossing gives:
#   φ'_linear ~ KAL × (M_src/M_Pl) / r²
#   φ'_NL     ~ [M⁴ × M_src/M_Pl]^{1/3} / r^{2/3}
# Equating: r_V ∝ (M_src/M_Pl)^{1/2} / M  (K-essence, n=2)
# Anchoring at Paper 8: r_V(Sun, M_Pl) ≡ r_V_sun_P8
# ══════════════════════════════════════════════════════════════════════════════
print("─"*70)
print("PART 2 — M scan: what cutoff gives cosmological Vainshtein?")
print("  Formula: r_V(source, M) = r_V(source, M_Pl) × (M_Pl/M)")
print("─"*70)

# Anchor: r_V scales as 1/M for K-essence with K = X + X²/M⁴
# r_V(cluster, M) = r_V_cluster_P8 × (M_Pl/M)

def r_V_cluster(M_GeV):
    """Vainshtein radius for 10^15 M_sun cluster as function of M (GeV)."""
    return r_V_cluster_P8 * (M_Pl_GeV / M_GeV)

# Target: r_V ∈ [50, 300] Mpc
target_50Mpc  = 50  * Mpc_m
target_100Mpc = 100 * Mpc_m
target_300Mpc = 300 * Mpc_m

M_for_50Mpc  = M_Pl_GeV * r_V_cluster_P8 / target_50Mpc
M_for_100Mpc = M_Pl_GeV * r_V_cluster_P8 / target_100Mpc
M_for_300Mpc = M_Pl_GeV * r_V_cluster_P8 / target_300Mpc

print(f"\n  For r_V(10¹⁵M☉) = 50 Mpc:   M = {M_for_50Mpc:.3e} GeV = {M_for_50Mpc*1e9:.3e} eV")
print(f"  For r_V(10¹⁵M☉) = 100 Mpc:  M = {M_for_100Mpc:.3e} GeV = {M_for_100Mpc*1e9:.3e} eV")
print(f"  For r_V(10¹⁵M☉) = 300 Mpc:  M = {M_for_300Mpc:.3e} GeV = {M_for_300Mpc*1e9:.3e} eV")

# ══════════════════════════════════════════════════════════════════════════════
# PART 3 — Algebraic candidates for M from SSEE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─"*70)
print("PART 3 — Algebraic candidates for M from φ and π")
print("─"*70)

# SSEE energy scales
rho_crit_GeV4 = 3 * H0_SI**2 * (M_Pl_GeV)**2 / (8*np.pi)   # (s⁻¹)² × GeV² → needs c factors
# In natural units: ρ_crit = 3 H₀² M_Pl² / (8π), H₀ in GeV
H0_GeV = H0_SI * 6.582e-25   # s⁻¹ × ħ (GeV·s) → GeV
rho_crit_nat = 3 * H0_GeV**2 * M_Pl_GeV**2 / (8*np.pi)      # GeV⁴
M_from_rhocrit = rho_crit_nat**(1/4)                          # Paper 7 convention

m_phi_eV  = 5.71          # eV
m_phi_GeV = m_phi_eV * 1e-9

# Geometric mean of M_Pl and m_phi
M_geomean_Pl_mphi = np.sqrt(M_Pl_GeV * m_phi_GeV)

# Geometric mean of M_Pl and ρ_crit^(1/4)
M_geomean_Pl_rho  = np.sqrt(M_Pl_GeV * M_from_rhocrit)

# SSEE algebraic scales
M_H0  = H0_GeV            # = H₀ in natural units
M_sq  = np.sqrt(H0_GeV * M_Pl_GeV)    # sqrt(H₀ × M_Pl) — crossover scale

candidates = {
    "ρ_crit^(1/4)  [Paper 7 norm]"    : M_from_rhocrit,
    "m_φ = 5.71 eV [φ-DM mass]"       : m_phi_GeV,
    "√(M_Pl × m_φ) [geom. mean]"      : M_geomean_Pl_mphi,
    "√(M_Pl × ρ_crit^¼) [geom. mean]": M_geomean_Pl_rho,
    "√(H₀ × M_Pl) [crossover]"        : M_sq,
    "H₀ [Hubble scale]"               : H0_GeV,
    "M_Pl [Planck, Paper 8 ref]"       : M_Pl_GeV,
}

print(f"\n  Target M for r_V=100 Mpc: {M_for_100Mpc:.3e} GeV = {M_for_100Mpc*1e9:.3e} eV\n")
print(f"  {'Candidate':<40}  {'M (GeV)':<16}  {'r_V cluster (Mpc)':<20}  Match?")
print(f"  {'-'*40}  {'-'*16}  {'-'*20}  {'-'*6}")

for name, M_val in candidates.items():
    rV = r_V_cluster(M_val) / Mpc_m
    in_range = "✅" if 50 <= rV <= 300 else "❌"
    print(f"  {name:<40}  {M_val:<16.3e}  {rV:<20.3e}  {in_range}")

# ══════════════════════════════════════════════════════════════════════════════
# PART 4 — EFT consistency check: for what M does the non-linear term
# X²/M⁴ remain subdominant at sub-Hubble scales (EFT validity)?
# Condition: X_bg / M⁴ << 1, where X_bg ~ ρ_DE ~ ρ_crit
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─"*70)
print("PART 4 — EFT validity: X_bg/M⁴ << 1 requires M >> ρ_crit^(1/4)")
print("─"*70)

print(f"\n  X_bg ~ ρ_DE ~ ρ_crit^(1/4) = {M_from_rhocrit:.3e} GeV")
print(f"\n  {'M (GeV)':<18}  {'X_bg/M⁴ ratio':<20}  EFT valid?")
print(f"  {'-'*18}  {'-'*20}  {'-'*10}")

for name, M_val in candidates.items():
    ratio = (M_from_rhocrit / M_val)**4
    valid = "✅ linear" if ratio < 0.1 else ("⚠️  marginal" if ratio < 1 else "❌ collapses")
    print(f"  {M_val:<18.3e}  {ratio:<20.3e}  {valid}  ({name[:30]})")

# ══════════════════════════════════════════════════════════════════════════════
# PART 5 — Definitive verdict
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5 — Definitive UV Completion Verdict")
print("=" * 70)

print(f"""
  The cosmological Vainshtein mechanism requires simultaneously:
  (A) r_V(cluster) ∈ [50, 300] Mpc  → M ∈ [{M_for_300Mpc:.2e}, {M_for_50Mpc:.2e}] GeV
                                      = [{M_for_300Mpc*1e9:.1f}, {M_for_50Mpc*1e9:.1f}] eV
  (B) EFT linear validity: X_bg/M⁴ < 0.1
      → M > {M_from_rhocrit * 10**0.25:.2e} GeV (roughly 1.8 × ρ_crit^(1/4) = {1.8*M_from_rhocrit*1e12:.2f} peV)

  Condition (A) range:   M ∈ [{M_for_300Mpc*1e9:.1e}, {M_for_50Mpc*1e9:.1e}] eV
  Condition (B) lower:   M > {1.8*M_from_rhocrit*1e12:.2f} pico-eV

  These are compatible in principle — a sub-eV M satisfies both if small enough.

  BUT: the required M ~ {M_for_100Mpc*1e9:.1e} eV has NO algebraic origin in φ and π.
  It does not emerge from any known SSEE scale combination.

  No algebraic M from φ and π gives cosmological Vainshtein screening.

  ╔══════════════════════════════════════════════════════════════════╗
  ║  CONCLUSION: Paper 9 via Vainshtein requires UV completion.      ║
  ║  The numerical correlation f_screen = (π-φ)/Ω² = {f_screen:.4f} is  ║
  ║  real but its physical mechanism is not derived.                  ║
  ║                                                                  ║
  ║  Additional finding: Paper 8 has a ×10¹⁸ error in r_V(Sun).     ║
  ║  Corrected: r_V(Sol, M=M_Pl) ≈ 2.4 fm (not 6 km).              ║
  ║  Conclusion unchanged: Vainshtein inoperative. Needs fix in P8. ║
  ╚══════════════════════════════════════════════════════════════════╝
""")

# ── f_screen algebraic check ──────────────────────────────────────────────────
H0_local = H0_alg / (1 - f_screen)
print(f"  Algebraic reminder:")
print(f"    f_screen = (π-φ)/Ω² = {f_screen:.6f}")
print(f"    H₀,local = H₀/(1-f) = {H0_alg:.4f}/(1-{f_screen:.4f}) = {H0_local:.4f} km/s/Mpc")
print(f"    SH0ES    = 73.04 ± 1.04  → tension = {abs(H0_local-73.04)/1.04:.2f}σ")
print(f"\n  This algebraic result is preserved as a pre-committed prediction,")
print(f"  pending UV completion of the action P(X,φ).\n")
print("=" * 70)
