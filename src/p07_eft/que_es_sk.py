"""Que cantidad fisica es
   s_K = 3(-w0)(1+w0) = 0.403302 ?

Candidato: la tasa de dilucion de la
PRESION, no de la energia.
  rho ~ a^{-3(1+w)}
  p = w rho
  dp/dN = -3w(1+w) rho
"""
import sys, os
import sympy as sp
sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..'))
import ssee_core as S

N, w = sp.symbols('N w')
r0 = sp.symbols('rho0', positive=True)
rho = r0 * sp.exp(-3 * (1 + w) * N)
p = w * rho
drho = sp.diff(rho, N) / rho
dp = sp.diff(p, N) / rho

print("  SIMBOLICO")
print("  (drho/dN)/rho =",
      sp.simplify(drho))
print("  (dp/dN)/rho   =",
      sp.simplify(dp))
print()

W0 = S.W0
v_rho = float(-3 * (1 + W0))
v_p = float(3 * W0 * (1 + W0))
print("  NUMERICO en w0=%.6f" % W0)
print("  -(drho/dN)/rho = %+.12f"
      % -v_rho)
print("  -(dp/dN)/rho   = %+.12f"
      % -v_p)
print("  s_K            = %+.12f"
      % S.S_K)
print("  dif con -dp/dN = %.2e"
      % abs(-v_p - S.S_K))
