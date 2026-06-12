"""
Forward derivation: next term X³/M⁸ in the K(X) UV ladder.

K(X) = X/KAL + X²/M⁴ + X³/M⁸  (coefficient 1 each — NO new free parameter,
same structure as Paper 10's X²/M⁴ truncation).

Physics-first: re-solve the self-consistent background EOM 2X·K_X = ρ_DE(1+w₀)
(now a CUBIC in X), recompute αK_full, f_screen, and H_local from the MIRA base.
We do NOT peek at any target; we report whatever boost comes out.
"""
import numpy as np
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssee_core import (
    PHI as phi, PI as pi, OMEGA as Omega, AURA, KAL0 as KAL, MIRA,
    W0 as w0, OMEGA_DE as Omega_DE, H0_ALG as H0_alg,
)

# ── anchor (same as Paper 10) ─────────────────────────────────────────────────
H_MIRA      = 67.037
H0_SH0ES    = 73.04
sigma_SH0ES = 1.04
M_PL_eV     = 2.435323e27
H0_PER_h_eV = 2.1331951e-33
H_MIRA_eV   = (H_MIRA / 100.0) * H0_PER_h_eV
RHO         = 3.0 * H_MIRA_eV**2 * M_PL_eV**2 * 1e12   # meV⁴
MplH2       = RHO / 3.0

alpha_attr = phi**4 / 3
M4_UV      = 5 * phi**8 * RHO        # = 45α² ρ_crit

alpha_K_IR = 3 * Omega_DE * (1 + w0)
rho_DE_1pw = (alpha_K_IR / 3) * RHO
f_screen_IR = (pi - phi) / Omega**2

def solve_X(coeffs_KX, rhs):
    """Solve 2X·K_X(X) = rhs for the physical positive root.
    coeffs_KX: list of (power_of_X, coeff/M^k) building K_X."""
    # Build polynomial 2X·K_X - rhs = 0 numerically via np.roots.
    # K_X = 1/KAL + 2X/M⁴ + 3X²/M⁸  → 2X·K_X = 2X/KAL + 4X²/M⁴ + 6X³/M⁸
    pass

# ── order-2 (Paper 10 reference): K_X = 1/KAL + 2X/M⁴ ─────────────────────────
# 2X·K_X = 2X/KAL + 4X²/M⁴ = rho_DE_1pw
a2 = 4.0 / M4_UV; b2 = 2.0 / KAL; c2 = -rho_DE_1pw
X2 = (-b2 + np.sqrt(b2**2 - 4*a2*c2)) / (2*a2)
KX2  = 1/KAL + 2*X2/M4_UV
KXX2 = 2/M4_UV
aK2  = 2*(X2*KX2 + 2*X2**2*KXX2) / MplH2

# ── order-3: K = X/KAL + X²/M⁴ + X³/M⁸ ────────────────────────────────────────
# K_X  = 1/KAL + 2X/M⁴ + 3X²/M⁸
# 2X·K_X = 2X/KAL + 4X²/M⁴ + 6X³/M⁸ = rho_DE_1pw
# cubic: 6/M⁸·X³ + 4/M⁴·X² + 2/KAL·X − rho_DE_1pw = 0
coeffs = [6.0/M4_UV**2, 4.0/M4_UV, 2.0/KAL, -rho_DE_1pw]
roots = np.roots(coeffs)
real_pos = [r.real for r in roots if abs(r.imag) < 1e-6 and r.real > 0]
X3 = min(real_pos)          # physical perturbative root
KX3  = 1/KAL + 2*X3/M4_UV + 3*X3**2/M4_UV**2
KXX3 = 2/M4_UV + 6*X3/M4_UV**2
aK3  = 2*(X3*KX3 + 2*X3**2*KXX3) / MplH2

def report(label, aK):
    f = aK / (3*MIRA)
    Hm = H_MIRA / (1 - f)
    Ha = H0_alg / (1 - f)
    tm = abs(Hm - H0_SH0ES)/sigma_SH0ES
    ta = abs(Ha - H0_SH0ES)/sigma_SH0ES
    print(f"  {label:22s} αK={aK:.6f}  f_screen={f:.6f}  "
          f"H_local(MIRA)={Hm:.4f} ({tm:.2f}σ)  H_local(alg)={Ha:.4f} ({ta:.2f}σ)")

print("="*96)
print("K(X) UV ladder — does the X³/M⁸ term add useful boost?  (base anchor = MIRA 67.037)")
print("="*96)
print(f"  M⁴/ρ_crit = 5φ⁸ = {5*phi**8:.4f}")
print(f"  X₂/M⁴ = {X2/M4_UV:.4e}   X₃/M⁴ = {X3/M4_UV:.4e}   (perturbative ratio)")
print(f"  X₃/X₂ = {X3/X2:.8f}\n")
report("IR (X/KAL only)", alpha_K_IR)
report("UV order-2 (X²/M⁴)", aK2)
report("UV order-3 (+X³/M⁸)", aK3)
print()
print(f"  Δ(order2→order3) in αK     = {aK3-aK2:+.6e}  ({(aK3/aK2-1)*100:+.4f}%)")
print(f"  Δ(order2→order3) in f_screen = {aK3/(3*MIRA)-aK2/(3*MIRA):+.6e}")
print(f"  Δ(order2→order3) in H_local(MIRA) = "
      f"{H_MIRA/(1-aK3/(3*MIRA)) - H_MIRA/(1-aK2/(3*MIRA)):+.5f} km/s/Mpc")
print("="*96)
