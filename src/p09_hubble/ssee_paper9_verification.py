"""
SSEE Paper 9 — Algebraic Identity Verification
f_screen = alpha_K / (3 * MIRA) = (pi - phi) / (phi + pi)^2

Verifies, step by step:
  Step 1: alpha_K algebraic expression from K(X)=X/KAL + Friedmann
  Step 2: MIRA = AURA/2 from Paper 8 disformal null geodesic
  Step 3: AURA cancels in alpha_K / (3 * MIRA)  →  (pi - phi) / Omega^2

Also computes both H0 corrections (additive and multiplicative) and
compares against SH0ES.
"""

import numpy as np

# ── Fundamental constants ────────────────────────────────────────────────────
import os as _reloc_os, sys as _reloc_sys  # reloc: anclar src/
_reloc_sys.path.insert(0, _reloc_os.path.dirname(_reloc_os.path.dirname(_reloc_os.path.abspath(__file__))))
from ssee_core import (
    PHI as phi, PI as pi, OMEGA as Omega, AURA, KAL0 as KAL,
    W0 as w0, OMEGA_DE as Omega_DE, H0_ALG as H0_alg, MIRA,
)

# ── Paper 1 background quantities (importados de ssee_core) ──────────────────
# Identidades verificadas: w₀ = -AURA/Ω = -Tr/Mv ;  Ω_DE = AURA/Ω ;
#                          H₀^alg = 3(φ+π)²  [km/s/Mpc]

# ── Step 1: exact alpha_K from Friedmann ────────────────────────────────────
# Friedmann: X/KAL = rho_phi*(1+w0)/2  →  alpha_K = 3*Omega_DE*(1+w0)
one_plus_w0    = (pi - phi) / (2 * Omega)   # exact identity Ω - AURA = (π-φ)/2
alpha_K_exact  = 3 * Omega_DE * (1 + w0)   # from background equations
alpha_K_formula = 3 * AURA * (pi - phi) / (2 * Omega**2)  # substituted

# Cross-check: should match EFTCAMB value from Paper 7
alpha_K_eftcamb = 0.40330         # Paper 7 CLASS/EFTCAMB result

# ── Step 2: MIRA from Paper 8 disformal null geodesic ───────────────────────
# theta_E^SSEE / theta_E^GR = sqrt(beta_c) = MIRA,  MIRA = AURA/2 (importado)

# ── Step 3: cancellation ─────────────────────────────────────────────────────
# alpha_K / (3*MIRA) = [3*AURA*(pi-phi)/(2*Omega^2)] / [3*AURA/2]
#                    = (pi - phi) / Omega^2
f_screen_identity = alpha_K_formula / (3 * MIRA)
f_screen_direct   = (pi - phi) / Omega**2

# Percentage discrepancy between formula and direct
discrepancy_pct = abs(f_screen_identity - f_screen_direct) / f_screen_direct * 100

# ── H0 predictions ───────────────────────────────────────────────────────────
H0_multiplicative = H0_alg / (1 - f_screen_direct)   # H₀/(1−f)
H0_additive       = H0_alg * (1 + f_screen_direct)   # H₀(1+f)
H0_SH0ES          = 73.04
sigma_SH0ES       = 1.04

tension_mult = abs(H0_multiplicative - H0_SH0ES) / sigma_SH0ES
tension_add  = abs(H0_additive       - H0_SH0ES) / sigma_SH0ES

# ── Print verification ───────────────────────────────────────────────────────
sep = "=" * 68
print(sep)
print("SSEE Paper 9 — Algebraic Identity Verification")
print(sep)

print(f"\n── Fundamental constants ─────────────────────────────────────────")
print(f"  φ      = {phi:.10f}")
print(f"  π      = {pi:.10f}")
print(f"  Ω=φ+π  = {Omega:.10f}")
print(f"  AURA   = (3φ+π)/2 = {AURA:.10f}")
print(f"  MIRA   = AURA/2   = {MIRA:.10f}")
print(f"  KAL    = (φ+π)/2+π = {KAL:.10f}")

print(f"\n── Step 1: exact alpha_K ─────────────────────────────────────────")
print(f"  w₀               = −AURA/Ω         = {w0:.10f}")
print(f"  1 + w₀           = (π−φ)/(2Ω)      = {one_plus_w0:.10f}")
print(f"  Ω_DE             = AURA/Ω           = {Omega_DE:.10f}")
print(f"  αK (background)  = 3·Ω_DE·(1+w₀)  = {alpha_K_exact:.8f}")
print(f"  αK (substituted) = 3AURA(π−φ)/2Ω² = {alpha_K_formula:.8f}")
print(f"  αK (EFTCAMB P7)  =                   {alpha_K_eftcamb:.8f}")
print(f"  Δ(formula/EFTCAMB)                 = {abs(alpha_K_formula-alpha_K_eftcamb)/alpha_K_eftcamb*100:.4f}%")

print(f"\n── Step 2: MIRA ──────────────────────────────────────────────────")
print(f"  MIRA = AURA/2 = (3φ+π)/4 = {MIRA:.10f}")
print(f"  3·MIRA         = 3AURA/2  = {3*MIRA:.10f}")

print(f"\n── Step 3: cancellation ──────────────────────────────────────────")
print(f"  αK / (3·MIRA)  = {f_screen_identity:.10f}")
print(f"  (π−φ)/Ω²       = {f_screen_direct:.10f}")
print(f"  Discrepancy    = {discrepancy_pct:.8f}%  ← AURA cancels exactly")

print(f"\n── H₀ predictions ───────────────────────────────────────────────")
print(f"  H₀,alg (Paper 1)    = {H0_alg:.6f} km/s/Mpc")
print(f"  f_screen            = {f_screen_direct:.6f}  ({f_screen_direct:.4%})")
print(f"")
print(f"  H₀,local (mult.)    = H₀/(1−f)  = {H0_multiplicative:.4f} km/s/Mpc  →  {tension_mult:.2f}σ from SH0ES")
print(f"  H₀,local (add.)     = H₀·(1+f)  = {H0_additive:.4f} km/s/Mpc  →  {tension_add:.2f}σ from SH0ES")
print(f"  H₀,SH0ES (Riess22)  =             {H0_SH0ES:.4f} km/s/Mpc  ± {sigma_SH0ES:.2f}")
print(f"")
diff_mult_add = H0_multiplicative - H0_additive
print(f"  Difference (mult − add) = {diff_mult_add:.4f} km/s/Mpc  ≈ f²·H₀  [{diff_mult_add:.4f}]")

print(f"\n── Key identity (exact) ─────────────────────────────────────────")
print(f"  f_screen = αK/(3·MIRA) = (π−φ)/(φ+π)² = {f_screen_direct:.8f}")
print(f"  AURA = {AURA:.6f} cancels in the ratio → result is pure φ,π")

print(f"\n── Cross-check: 1+w₀ identity ───────────────────────────────────")
# Verify Ω - AURA = (pi-phi)/2
lhs = Omega - AURA
rhs = (pi - phi) / 2
print(f"  Ω − AURA           = {lhs:.10f}")
print(f"  (π − φ)/2          = {rhs:.10f}")
print(f"  Difference         = {abs(lhs-rhs):.2e}  (machine precision)")

print(f"\n{'='*68}")
print("VERDICT")
print(f"{'='*68}")
print(f"""
  The identity f_screen = αK/(3·MIRA) = (π−φ)/(φ+π)² holds
  with {discrepancy_pct:.1e}% discrepancy — limited only by floating-point
  precision. AURA cancels analytically; the result depends only on φ and π.

  Both H₀ corrections (mult: {H0_multiplicative:.2f}, add: {H0_additive:.2f} km/s/Mpc)
  are consistent with SH0ES ({H0_SH0ES} ± {sigma_SH0ES}):
    Multiplicative: {tension_mult:.2f}σ
    Additive:       {tension_add:.2f}σ
""")
print("=" * 68)
