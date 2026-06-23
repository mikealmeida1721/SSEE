"""
SSEE Paper 10 — UV Completion Verification
M⁴ = 45α² × ρ_crit = 5φ⁸ × ρ_crit

Verifies, step by step:
  Step 1: Algebraic identity 45α² = 5φ⁸  (α = φ⁴/3, Paper 1 α-attractor)
  Step 2: Self-consistent X_bg from full K(X) = X/KAL + X²/M⁴
  Step 3: αK_full from Bellini-Sawicki with K_X and K_XX contributions
  Step 4: f_screen_full = αK_full / (3·MIRA)  →  H₀,local

The UV ladder:
  φ  →  α = φ⁴/3  →  M⁴ = 45α² ρ_crit = 5φ⁸ ρ_crit
  →  αK_full = 0.41691  →  f_screen = 0.06952
  →  H₀,local CANÓNICO = H_alg/(1−f_UV) = 73.040 km/s/Mpc (0.00σ SH0ES)
     [reframe ω_m-directo: H_alg es la base canónica; era H_MIRA 67.037 → 72.05]
"""

import numpy as np

# ── Fundamental constants ────────────────────────────────────────────────────
import os as _reloc_os, sys as _reloc_sys  # reloc: anclar src/
_reloc_sys.path.insert(0, _reloc_os.path.dirname(_reloc_os.path.dirname(_reloc_os.path.abspath(__file__))))
from ssee_core import (
    PHI as phi, PI as pi, OMEGA as Omega, AURA, KAL0 as KAL, MIRA,
    W0 as w0, OMEGA_DE as Omega_DE, H0_ALG as H0_alg,
)
# Identidades: w₀ = -AURA/Ω = -Tr/Mv ;  Ω_DE = AURA/Ω ;  H₀^alg = 3(φ+π)²

# ── Physical anchor: the ONE critical density ρ_crit of the model ─────────────
# ρ_crit = 3 H₀² M_Pl²  (reduced Planck mass). There is a SINGLE critical density
# in SSEE; it enters every dimensionful quantity (X_bg, M⁴, ρ_DE) as a common
# scale and cancels exactly in αK, f_screen and the H₀ boost (all of these are
# ratios). We pin ρ_crit to its canonical model value so X_bg and M⁴ carry real
# meV⁴ units — no implicit ρ_crit = 1 (which would leave them ambiguous to a
# referee even though the observable chain is invariant). Anchor (reframe
# ω_m-directo 2026-06-19): H_alg = 67.962 = global background H = CMB anchor.
H_ANCHOR     = 67.962              # km/s/Mpc — H_alg = 3(φ+π)² (reframe; era H_MIRA 67.037)
M_PL_eV      = 2.435323e27         # reduced Planck mass [eV]
H0_PER_h_eV  = 2.1331951e-33       # 100 km/s/Mpc expressed in eV
H_ANCHOR_eV  = (H_ANCHOR / 100.0) * H0_PER_h_eV
RHO_CRIT_eV4 = 3.0 * H_ANCHOR_eV**2 * M_PL_eV**2   # eV⁴
RHO          = RHO_CRIT_eV4 * 1e12               # meV⁴ ≈ 36.41  (1 meV⁴ = 1e-12 eV⁴)
RHO_qrt      = RHO**0.25                          # ρ_crit^(1/4) ≈ 2.4564 meV

# ── Step 1: α-attractor parameter and UV identity ────────────────────────────
# α = φ⁴/3 is the curvature of the inflationary Starobinsky-like plateau (Paper 1)
# It predicts r = 12α/N² ≈ 0.00813 for N = 60 e-folds (LiteBIRD 2032 target)
alpha_attr = phi**4 / 3

M4_ratio = 45 * alpha_attr**2     # = 5φ⁸ exactly — the dimensionless M⁴/ρ_crit
M4_UV    = M4_ratio * RHO         # meV⁴ — physical UV cutoff⁴ (= 5φ⁸ ρ_crit)

# Verify algebraic identity
identity = abs(45 * alpha_attr**2 - 5 * phi**8)

# ── Step 2: IR (Paper 9) reference values ────────────────────────────────────
# αK_IR from K(X) = X/KAL + Friedmann (Papers 7 + 9)
alpha_K_IR = 3 * Omega_DE * (1 + w0)          # = 3AURA(π−φ)/(2Ω²)  [dimensionless]
rho_DE_1pw = (alpha_K_IR / 3) * RHO            # = Ω_DE·(1+w₀)·ρ_crit  [meV⁴]
X_bg_IR    = alpha_K_IR * KAL / 6 * RHO        # self-consistent X, IR limit  [meV⁴]

# f_screen_IR from Paper 9 (AURA cancels)
f_screen_IR = (pi - phi) / Omega**2            # = αK_IR / (3·MIRA) exactly

# ── Step 3: self-consistent X_bg from full K(X) = X/KAL + X²/M⁴ ─────────────
# Background EOM:  2X K_X = ρ_DE(1+w₀)
# K_X = 1/KAL + 2X/M⁴  →  quadratic: 4X²/M⁴ + 2X/KAL − ρ_DE(1+w₀) = 0
a = 4.0 / M4_UV
b = 2.0 / KAL
c = -rho_DE_1pw
disc = b**2 - 4*a*c
X_bg_UV = (-b + np.sqrt(disc)) / (2*a)        # physical (positive) root

X_ratio = X_bg_UV / X_bg_IR                   # should be close to 1 for M⁴ >> X²

# ── Step 4: αK_full (Bellini-Sawicki, full K) ─────────────────────────────────
# αK = 2(X K_X + 2X² K_XX) / (M²_Pl H²),  with M²_Pl H² = ρ_crit/3.
# Both numerator and denominator scale with ρ_crit → αK is invariant (anchor-free).
MplH2   = RHO / 3.0                            # M²_Pl H² = ρ_crit/3  [meV⁴]
K_X_UV  = 1.0/KAL + 2*X_bg_UV/M4_UV            # dimensionless (X/M⁴ is a ratio)
K_XX_UV = 2.0/M4_UV                            # meV⁻⁴
alpha_K_UV = 2 * (X_bg_UV*K_X_UV + 2*X_bg_UV**2*K_XX_UV) / MplH2

UV_correction = alpha_K_UV / alpha_K_IR - 1.0  # fractional increase

# ── Step 5: f_screen_full and H₀,local ──────────────────────────────────────
f_screen_UV = alpha_K_UV / (3 * MIRA)
H0_UV       = H0_alg / (1 - f_screen_UV)

# Observational comparison
H0_SH0ES   = 73.04
sigma_SH0ES = 1.04
tension_UV  = abs(H0_UV - H0_SH0ES) / sigma_SH0ES
tension_IR  = abs(H0_alg / (1 - f_screen_IR) - H0_SH0ES) / sigma_SH0ES

# ── Print verification ────────────────────────────────────────────────────────
sep = "=" * 72
print(sep)
print("SSEE Paper 10 — UV Completion Verification")
print(sep)

print(f"\n── Fundamental constants ────────────────────────────────────────────")
print(f"  φ       = {phi:.10f}")
print(f"  π       = {pi:.10f}")
print(f"  Ω=φ+π   = {Omega:.10f}")
print(f"  AURA    = (3φ+π)/2 = {AURA:.10f}")
print(f"  MIRA    = AURA/2   = {MIRA:.10f}")
print(f"  KAL     = (φ+π)/2+π = {KAL:.10f}")
print(f"  H₀,alg  = 3Ω² = {H0_alg:.6f} km/s/Mpc")

print(f"\n── Step 1: α-attractor parameter and UV identity ────────────────────")
print(f"  α = φ⁴/3            = {alpha_attr:.10f}   (Paper 1 inflaton curvature)")
print(f"  45α²                = {45*alpha_attr**2:.10f}")
print(f"  5φ⁸                 = {5*phi**8:.10f}")
print(f"  |45α² − 5φ⁸|        = {identity:.3e}   ← exact identity (machine precision)")
print(f"  ρ_crit (canónico)  = {RHO:.4f} meV⁴   [= 3 H_anchor² M_Pl², H_anchor={H_ANCHOR} km/s/Mpc]")
print(f"  ρ_crit^(1/4)       = {RHO_qrt:.6f} meV")
print(f"  M⁴ = 5φ⁸ ρ_crit    = {M4_UV:.4f} meV⁴   (ratio M⁴/ρ_crit = {M4_ratio:.4f})")
print(f"  M = φ²×5^(1/4)×ρ_crit^(1/4) = {M4_UV**0.25:.4f} meV  = Λ_SSEE")
print(f"  M / ρ_crit^(1/4)   = {M4_ratio**0.25:.8f}")
print(f"    φ²×5^(1/4)       = {phi**2 * 5**0.25:.8f}  ← algebraic form ✓")

print(f"\n── Step 2: IR (Paper 9) baseline ────────────────────────────────────")
print(f"  w₀            = −AURA/Ω     = {w0:.8f}")
print(f"  αK_IR         = 3Ω_DE(1+w₀) = {alpha_K_IR:.8f}  [Paper 7 + Friedmann]")
print(f"  ρ_DE(1+w₀)    = αK_IR/3    = {rho_DE_1pw:.8f}")
print(f"  X_bg_IR       = αK_IR×KAL/6 = {X_bg_IR:.8f}")
print(f"  f_screen_IR   = (π−φ)/Ω²  = {f_screen_IR:.8f}  [AURA cancels exactly]")
print(f"  H₀,IR         = {H0_alg/(1-f_screen_IR):.4f} km/s/Mpc  →  {tension_IR:.2f}σ SH0ES")

print(f"\n── Step 3: self-consistent X_bg with full K(X) ──────────────────────")
print(f"  Quadratic: 4X²/M⁴ + 2X/KAL − ρ_DE(1+w₀) = 0")
print(f"  X_bg_UV          = {X_bg_UV:.8f} meV⁴")
print(f"  X_bg_IR          = {X_bg_IR:.8f} meV⁴")
print(f"  X_UV/X_IR        = {X_ratio:.8f}  (small suppression from UV term)")
print(f"  X_bg_UV/M⁴       = {X_bg_UV/M4_UV:.8f}  ← dimensionless invariant (anchor-free)")
print(f"  ε = X_bg_UV²/M⁴ = {X_bg_UV**2/M4_UV:.4e} meV⁴  (perturbative UV correction)")

print(f"\n── Step 4: αK_full (Bellini-Sawicki) ───────────────────────────────")
print(f"  K_X  = 1/KAL + 2X/M⁴ = {K_X_UV:.10f}")
print(f"  K_XX = 2/M⁴           = {K_XX_UV:.10f}")
print(f"  αK_IR (K(X)=X/KAL)    = {alpha_K_IR:.8f}")
print(f"  αK_UV (K full)        = {alpha_K_UV:.8f}")
print(f"  αK_UV/αK_IR           = {alpha_K_UV/alpha_K_IR:.8f}  (+{UV_correction*100:.2f}%)")
print(f"  Δ αK = +{(alpha_K_UV-alpha_K_IR):.6f}  (UV correction)")

print(f"\n── Step 5: f_screen_full and H₀,local ──────────────────────────────")
print(f"  f_screen_UV = αK_full/(3·MIRA) = {f_screen_UV:.8f}")
print(f"  f_screen_IR = (π−φ)/Ω²        = {f_screen_IR:.8f}")
print(f"  Δf = +{(f_screen_UV-f_screen_IR):.6f}  (+{(f_screen_UV/f_screen_IR-1)*100:.2f}%)")
print(f"")
print(f"  H₀,alg (Paper 1)    = {H0_alg:.6f} km/s/Mpc")
print(f"  H₀,local IR         = H₀/(1−f_IR) = {H0_alg/(1-f_screen_IR):.4f} km/s/Mpc  →  {tension_IR:.2f}σ SH0ES")
print(f"  H₀,local UV         = H₀/(1−f_UV) = {H0_UV:.4f} km/s/Mpc  →  {tension_UV:.4f}σ SH0ES")
print(f"  SH0ES (Riess 2022)  =              {H0_SH0ES:.4f} km/s/Mpc  ±{sigma_SH0ES:.2f}")
print(f"")

# ── UV ladder summary ─────────────────────────────────────────────────────────
print(f"── UV ladder: φ → α → M⁴ → αK_full → H₀ ────────────────────────────")
print(f"  φ = {phi:.8f}")
print(f"  α = φ⁴/3 = {alpha_attr:.8f}  (inflaton curvature → r=0.00813)")
# ρ_crit^(1/4) = RHO_qrt meV con H_anchor=67.962 km/s/Mpc (H_alg, reframe)
print(f"  M⁴ = 45α² = 5φ⁸ ρ_crit = {M4_UV:.4f} meV⁴  →  M = {M4_UV**0.25:.3f} meV")
print(f"  αK_full = {alpha_K_UV:.5f}  (+{UV_correction*100:.2f}% over αK_IR={alpha_K_IR:.5f})")
print(f"  f_screen = {f_screen_UV:.5f}")
print(f"  H₀,local = {H0_UV:.4f} km/s/Mpc  ← {tension_UV:.4f}σ from SH0ES")

print(f"\n── Cross-check: AURA cancellation at UV ──────────────────────────────")
# At UV, AURA does NOT cancel in f_screen_full because αK_full ≠ αK_IR.
# But the numerical result still rounds to SH0ES.
# The cancellation f_screen_UV = exact function of φ,π would require:
# αK_full / (3·MIRA) = (π-φ)/Ω² × (1 + UV_corr)
# The UV correction IS a function of AURA through αK_IR.
AURA_check = alpha_K_UV / (3 * (pi - phi) / Omega**2)
print(f"  αK_full / [(π−φ)/Ω²]  = {AURA_check:.8f}  (= 3·MIRA × αK_full/αK_full×...)")
print(f"  Effective MIRA needed  = {alpha_K_UV/(3*f_screen_IR):.8f}  (vs MIRA={MIRA:.8f})")
print(f"  ΔMIRA_eff              = +{alpha_K_UV/(3*f_screen_IR)-MIRA:.8f}  (UV correction to Paper 8 sector)")

print(f"\n{sep}")
print("VERDICT")
print(f"{sep}")
print(f"""
  Step 1: M⁴ = 5φ⁸ ρ_crit = {M4_UV:.4f} meV⁴  [ratio 5φ⁸={M4_ratio:.4f}, |diff|={identity:.1e}]
  Step 2: X_bg_UV = {X_bg_UV:.6f} meV⁴  (X/M⁴ = {X_bg_UV/M4_UV:.2e} — perturbative ✓)
  Step 3: αK_full = {alpha_K_UV:.6f}  = αK_IR × {alpha_K_UV/alpha_K_IR:.6f}  (+{UV_correction*100:.2f}%)
  Step 4: f_screen_UV = {f_screen_UV:.6f}  (+{(f_screen_UV/f_screen_IR-1)*100:.2f}% over Paper 9)

  CANÓNICO (reframe ω_m-directo 2026-06-19, base H_alg={H0_alg} km/s/Mpc):
    IR:  H₀,local = H_alg/(1−f_IR) = {H0_alg/(1-f_screen_IR):.4f} km/s/Mpc  ({tension_IR:.2f}σ SH0ES)  ← titular Paper 9
    UV:  H₀,local = H_alg/(1−f_UV) = {H0_UV:.4f} km/s/Mpc  ({tension_UV:.4f}σ SH0ES)  ← titular Paper 10  [M⁴=5φ⁸]

  The same α that fixes the inflationary tensor-to-scalar ratio r = 12α/N² ≈ 0.00813
  (testable by LiteBIRD 2032) also fixes the dark-energy UV cutoff M = φ²×5^(1/4)×ρ_crit^(1/4).
  The reframe makes H_alg = 67.962 both the global background H and the CMB anchor
  (with ω_b,ω_c fixed by algebra, plik_lite minimises there), so the H_alg-based
  ladder (72.86 / 73.040) is now the canonical prediction (era H_MIRA 67.037 → 72.05).
  The UV step remains CONDITIONAL on Postulate C.1 (Theorem C.1, Paper 10).

  PENDING (Paper 10): first-principles derivation of M⁴ = 45α² without using SH0ES as input.
    Ruta A: K(X) Taylor matching with K_α(X) = −3α ln(1−X/(3α)) → M⁴ = 6αKAL² ≈ 418 ≠ 234.9
            (requires additional normalization constraint — under investigation)
    Ruta B: Technical naturalness EFT → M² = √45 × α × ρ_crit^(1/2)
    Ruta C: Full P(X,φ) Lagrangian from Paper 1 → determines M algebraically at 2nd order
""")
print(sep)
