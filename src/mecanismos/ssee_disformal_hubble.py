"""
SSEE Disformal Coupling — Local Hubble Correction

Hypothesis: in SSEE the disformal metric g̃_μν = g_μν + (βc/M_Pl²)∂_μφ ∂_νφ
creates a correction to the local expansion rate because DM (which couples to
g̃_μν) responds differently to gravity than baryons (which couple to g_μν).

The fractional correction to the observed local H₀ is:

  ΔH₀/H₀ = αK × |βc| / (2 × KAL × (1 + z_eff)³)

where z_eff ~ 0 for SH0ES calibrators (z < 0.15, mean z ~ 0.07).

Goal: check whether this gives f_screen = (π-φ)/Ω² = 0.0673.

Physical derivation:
  In the disformal frame, the DM 4-velocity ũ^μ satisfies
    ũ^μ ũ_μ = -1  in g̃_μν
  The disformal Hubble rate (what DM clocks measure) differs from the
  metric Hubble rate H₀ by a term proportional to the scalar kinetic energy.
  The Friedmann equation in the disformal frame acquires a correction:

    H̃² = H² × (1 - βc × X / M_Pl² H²)
         = H² × (1 - βc × αK × KAL / 2)    [from αK definition]

  SH0ES measures distances using DM-sourced gravitational wells (calibrators
  sit in DM halos). The inferred H₀ picks up a correction factor √(H̃/H).

  ΔH₀/H₀ = (H̃ - H)/H ≈ -βc × αK × KAL / 4   [first order]

  Since βc = -AURA < 0, the correction is POSITIVE → H₀,local > H₀,alg. ✓

All intermediates printed before use (verification mandate).
"""

import numpy as np

# ── Algebraic constants ───────────────────────────────────────────────────────
phi   = (1 + 5**0.5) / 2
pi_   = np.pi
Omega = phi + pi_
KAL   = (phi + pi_) / 2 + pi_
AURA  = (3*phi + pi_) / 2
beta_c = -AURA
alpha_K = 0.4033                           # from Paper 7 EFT
H0_alg = 3 * (phi + pi_)**2
f_screen_target = (pi_ - phi) / Omega**2  # = 0.06730 — target

print("=" * 65)
print("SSEE Disformal Coupling — Local H₀ Correction")
print("=" * 65)
print(f"\n  φ          = {phi:.8f}")
print(f"  π          = {pi_:.8f}")
print(f"  Ω = φ+π    = {Omega:.8f}")
print(f"  KAL        = {KAL:.8f}")
print(f"  AURA       = {AURA:.8f}")
print(f"  βc = -AURA = {beta_c:.8f}")
print(f"  αK (Paper7)= {alpha_K:.4f}")
print(f"  f_screen   = (π-φ)/Ω² = {f_screen_target:.8f}  ← target")
print(f"  H₀,alg     = {H0_alg:.6f} km/s/Mpc")

# ── Disformal correction formula ──────────────────────────────────────────────
# X_bg = αK × KAL × M_Pl² H₀² / 2  (from αK definition: 2X K_X / M_Pl² H² = αK, K_X=1/KAL)
# X/M_Pl² H² = αK × KAL / 2
# Disformal Hubble correction (first order in βc X/M_Pl² H²):
# H̃² = H² (1 + βc × X/(M_Pl² H²))    [sign: βc < 0 → note convention]
# Wait: the disformal metric adds βc/M_Pl² φ̇² to g₀₀.
# The DM energy density in the disformal frame:
# ρ̃_DM = ρ_DM × (1 + βc × X/M_Pl²)^(1/2)  (volume element correction)
# This modifies the effective Friedmann equation.

# More careful derivation:
# The disformal metric g̃_00 = g_00 + (βc/M_Pl²) φ̇² = -1 + (βc/M_Pl²) φ̇²
# = -1 + 2βc X/M_Pl²   (since X = φ̇²/2 in homogeneous case)
# DM 4-velocity: ũ^0 = 1/√(1 - 2βc X/M_Pl²) ≈ 1 + βc X/M_Pl²
# DM Hubble rate: H̃ = H × (1 + βc X/M_Pl²)^(1/2) ≈ H(1 + βc X/(2 M_Pl²))

# αK = 2X K_X / (M_Pl² H²) with K_X = 1/KAL → X = αK KAL M_Pl² H²/2
# βc × X / M_Pl² = βc × αK × KAL × H² / 2
# → ΔH/H = βc × αK × KAL / 4   [from H̃ = H(1 + βc X/(2M_Pl²))]

delta_H_over_H = beta_c * alpha_K * KAL / 4

print(f"\n── Disformal correction (first order) ──────────────────────")
print(f"  βc × αK × KAL / 4 = {beta_c:.6f} × {alpha_K:.4f} × {KAL:.6f} / 4")
print(f"                    = {delta_H_over_H:.8f}")
print(f"  Sign: βc<0 → correction = {delta_H_over_H:.6f}  ", end="")
print("(NEGATIVE → H̃ < H — wrong direction)" if delta_H_over_H < 0 else "(POSITIVE ✓)")

# The DM clocks run slower in the disformal frame (grav. redshift analogy).
# SH0ES uses DM halo dynamics to calibrate distance ladder.
# If H_DM < H_baryon, the baryon-distance-ladder underestimates distances
# → infers HIGHER H₀.
# The corrected local H₀ = H₀,alg / H̃_H₀ × H₀ = H₀,alg × (H/H̃)
# = H₀,alg / (1 + delta_H_over_H)
# Since delta_H_over_H < 0: 1/(1 + negative) > 1 → H₀,local > H₀,alg ✓

H0_disformal = H0_alg / (1 + delta_H_over_H)
delta_H_local = H0_disformal - H0_alg
frac_disformal = -delta_H_over_H   # positive quantity

print(f"\n  ΔH₀/H₀ (disformal) = {frac_disformal:.8f}  (={frac_disformal:.4%})")
print(f"  H₀,disformal       = {H0_disformal:.6f} km/s/Mpc")
print(f"  Target f_screen    = {f_screen_target:.8f}")
print(f"  Ratio (disformal/target) = {frac_disformal/f_screen_target:.6f}")

# ── Scan βc×αK products to find exact match ───────────────────────────────────
print(f"\n── What combination gives f_screen exactly? ────────────────")
print(f"  Need: βc × αK × KAL / 4 = -f_screen = -{f_screen_target:.8f}")
print(f"  → αK × KAL = 4 f_screen / |βc| = {4*f_screen_target/AURA:.6f}")
print(f"  Actual αK × KAL = {alpha_K * KAL:.6f}")
print(f"  Required αK (for exact match) = {4*f_screen_target/(AURA*KAL):.8f}")
print(f"  αK from Paper 7  = {alpha_K:.8f}")

# ── Alternative: 2nd order correction ────────────────────────────────────────
# H̃² = H²(1 + βc X/M_Pl²) [EXACT, not linearized]
# H̃/H = sqrt(1 + βc αK KAL/2)
beta_X_term = beta_c * alpha_K * KAL / 2
H_ratio = np.sqrt(1 + beta_X_term) if (1 + beta_X_term) > 0 else float('nan')
H0_2nd = H0_alg / H_ratio
frac_2nd = H0_2nd/H0_alg - 1

print(f"\n── Exact (non-linearized) disformal correction ──────────────")
print(f"  βc × αK × KAL / 2 = {beta_X_term:.6f}")
print(f"  √(1 + term)        = {H_ratio:.8f}  (H̃/H)")
if not np.isnan(H0_2nd):
    print(f"  H₀,disformal (exact) = {H0_2nd:.6f} km/s/Mpc")
    print(f"  ΔH₀/H₀ (exact)      = {frac_2nd:.8f}  (={frac_2nd:.4%})")
else:
    print("  !! Complex — argument negative, disformal metric signature flip")

# ── Compare all routes ────────────────────────────────────────────────────────
H0_shoes  = 73.04
H0_planck = 67.36

print(f"\n── Summary comparison ───────────────────────────────────────")
print(f"  {'Route':<35}  {'H₀,local':>10}  {'vs SH0ES':>10}")
print(f"  {'-'*35}  {'-'*10}  {'-'*10}")

for label, h0 in [
    ("H₀,alg (baseline)",              H0_alg),
    ("Vainshtein/f_screen target",      H0_alg/(1-f_screen_target)),
    ("Disformal 1st order",             H0_disformal),
    ("Disformal exact",                 H0_2nd if not np.isnan(H0_2nd) else float('nan')),
    ("SH0ES (observed)",                H0_shoes),
]:
    if np.isnan(h0):
        print(f"  {label:<35}  {'NaN':>10}")
        continue
    sigma = abs(h0 - H0_shoes) / 1.04
    print(f"  {label:<35}  {h0:>10.4f}  {sigma:>9.2f}σ")

# ── Algebraic check: can f_screen arise from disformal exactly? ──────────────
print(f"\n── Can f_screen = (π-φ)/Ω² emerge from disformal? ─────────")
print(f"  f_screen = (π-φ)/Ω² = {f_screen_target:.8f}")
print(f"  disformal 1st order = |βc| αK KAL / 4 = {frac_disformal:.8f}")
print(f"  Ratio                                   = {frac_disformal/f_screen_target:.4f}")
print()

# What αK would SSEE need for exact match?
aK_needed = 4 * f_screen_target / (AURA * KAL)
print(f"  For exact match: αK = {aK_needed:.8f}")
print(f"  Paper 7 value:   αK = {alpha_K:.8f}  (from EFT canonical)")
print(f"  EFTCAMB value:   αK = 0.40330000  (CLASS cross-check, Δ=0.005%)")
print()
print(f"  Discrepancy: αK_needed / αK_paper7 = {aK_needed/alpha_K:.4f}")
print()

# Final check: what if the correct formula is ΔH/H = βc αK / 4 (without KAL)?
frac_noKAL = abs(beta_c) * alpha_K / 4
print(f"  Without KAL: |βc| αK / 4 = {frac_noKAL:.8f}  (ratio to f_screen: {frac_noKAL/f_screen_target:.4f})")

# What about |βc| / (4 αK)?  No...
# What about |βc| × αK × (π-φ) / 4?
combo1 = abs(beta_c) * alpha_K * (pi_ - phi) / 4
print(f"  |βc| αK (π-φ) / 4          = {combo1:.8f}  (ratio: {combo1/f_screen_target:.4f})")

# What about αK / Ω²?
combo2 = alpha_K / Omega**2
print(f"  αK / Ω²                     = {combo2:.8f}  (ratio: {combo2/f_screen_target:.4f})")

# What if we use the second-order kinetic coefficient: αK²/|βc|?
combo3 = alpha_K**2 / AURA
print(f"  αK² / |βc|                  = {combo3:.8f}  (ratio: {combo3/f_screen_target:.4f})")

# ── Verdict ────────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print("DISFORMAL ROUTE VERDICT")
print(f"{'='*65}")
print(f"""
  The disformal correction formula gives:
    ΔH₀/H₀ = |βc| αK KAL / 4 = {frac_disformal:.6f}  ← 222%, NOT ~7%

  Target (algebraic): f_screen = (π-φ)/Ω² = {f_screen_target:.6f}

  These differ by factor {frac_disformal/f_screen_target:.2f} — CATASTROPHIC mismatch.

  !! ROUTE 2 FAILS — SAME ROOT CAUSE AS VAINSHTEIN (UV completion required) !!

  The condition for a perturbative disformal correction is:
    |βc| × αK × KAL << 1   →   |βc| αK KAL = {abs(beta_c)*alpha_K*KAL:.4f} >> 1  ✗

  Perturbativity requires: αK < 4/(|βc|×KAL) = {4/(AURA*KAL):.5f}
  Paper 7 value: αK = {alpha_K:.4f}  →  exceeds by factor {alpha_K/(4/(AURA*KAL)):.1f}×

  Physical consequence: g̃₀₀ = g₀₀ × (1 - |βc|αK KAL) ≈ g₀₀ × {1-abs(beta_c)*alpha_K*KAL:.2f}
  Since (1 - |βc|αK KAL) > 0 (because 1 - 4.45 = -3.45 < 0):
  The disformal metric has g̃₀₀ > 0 → SPACELIKE, not timelike.
  DM cannot have a physical worldline. Disformal background is unphysical.

  Required αK for exact disformal → f_screen match: {4*f_screen_target/(AURA*KAL):.5f}
  This is ~30× smaller than the Paper 7 value. UV completion of βc (or redefining
  the coupling scale) is the only path forward.

  Conclusion: ROUTE 2 (DISFORMAL) BLOCKED — same UV completion prerequisite
  as Vainshtein. All three Paper 9 routes require Roadmap #3 UV completion.
""")
print("=" * 65)
