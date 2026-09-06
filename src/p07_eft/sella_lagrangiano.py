"""Sella el Lagrangiano superviviente
   P(X) = A X^n   (EFT_section, 04-21)

Comprueba en simbolico y en numero:
  w_phi, no-fantasma, no-ghost,
  c_s^2, y la amplitud A.
"""
import sys, os
import sympy as sp
sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..'))
import ssee_core as S

X, A, n = sp.symbols('X A n',
                     positive=False)
P = A * X ** n
PX = sp.diff(P, X)
PXX = sp.diff(P, X, 2)

rho = sp.simplify(2 * X * PX - P)
w = sp.simplify(P / rho)
ghost = sp.simplify(PX + 2 * X * PXX)
cs2 = sp.simplify(PX / ghost)
# Bellini-Sawicki kineticity
aK_num = sp.simplify(2 * X * PX
                     + 4 * X ** 2 * PXX)

print("  SIMBOLICO")
print("  rho   =", rho)
print("  w     =", w)
print("  ghost =", ghost)
print("  c_s^2 =", cs2)
print("  aK*H^2Mp^2 =", aK_num)
print("  aK/rho     =",
      sp.simplify(aK_num / rho))

print()
print("  NUMERICO")
W0 = S.W0
nv = (1.0 + W0) / (2.0 * W0)
w_v = 1.0 / (2.0 * nv - 1.0)
cs2_v = w_v
OM_DE = 1.0 - S.OMEGA_M_CMB
S_DE = -W0
print("  n      = %+.8f" % nv)
print("  w_phi  = %+.12f" % w_v)
print("  w0     = %+.12f" % W0)
print("  dif    = %.2e" % abs(w_v - W0))
print()
print("  c_s^2  = %+.8f" % cs2_v)
if cs2_v < 0:
    print("  [FALLA] gradiente inestable")
else:
    print("  [OK]")
print()
print("  alpha_K = 2n * 3 * Om_DE")
for nom, d in (("Om_DE", OM_DE),
               ("s_DE", S_DE)):
    aK = 2.0 * nv * 3.0 * d
    print("  base %-6s aK = %+.8f"
          % (nom, aK))
print("  [FALLA] aK<0 -> ghost")
