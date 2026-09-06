"""Verifica la normalizacion de
zeta_0 = KAL_0 * H / (8 pi G).

Friedmann fija 8 pi G:
  rho_crit = 3 H^2 / (8 pi G)
  => H^2/(8 pi G) = rho_crit/3
No hay convencion libre.
"""
import sys, os
sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..'))
import ssee_core as S

PHI = S.PHI
PI = S.PI
KAL = S.KAL0
W0 = S.W0
S_DE = -W0
OM_DE = 1.0 - S.OMEGA_M_CMB

print("  NORMALIZACION DE zeta_0")
print("  KAL_0    = %.6f" % KAL)
print("  w0       = %.6f" % W0)
print("  Om_DE    = %.6f" % OM_DE)
print("  s_DE     = %.6f" % S_DE)
print()

# Pi_0 = -3 zeta_0 H
#      = -3 KAL H^2/(8 pi G)
#      = -3 KAL rho_crit/3
#      = -KAL rho_crit
pi_sobre_crit = -KAL
print("  Pi_0/rho_crit = %+.6f" % pi_sobre_crit)
print()

print("  Pi_0/rho_DE segun que use de rho_DE")
print("  base          Pi_0/rho_DE   objetivo")
print("  " + "-" * 42)
for nom, d in (("Om_DE", OM_DE), ("s_DE", S_DE)):
    r = pi_sobre_crit / d
    print("  %-12s  %+11.6f   %+.6f"
          % (nom, r, W0))
print("  " + "-" * 42)
print()

# que coeficiente c en zeta = c H/(8 pi G)
# daria w_eff = w0 con p_phi = 0 ?
print("  c requerido si zeta = c*H/(8 pi G)")
print("  base        c_req    KAL_0    razon")
print("  " + "-" * 42)
for nom, d in (("Om_DE", OM_DE), ("s_DE", S_DE)):
    c = S_DE * d
    print("  %-10s %8.6f %8.6f %7.3f"
          % (nom, c, KAL, KAL / c))
print("  " + "-" * 42)

print()
print("  DOBLE CONTEO con la seccion P(X)")
n = (1.0 + W0) / (2.0 * W0)
w_phi = 1.0 / (2.0 * n - 1.0)
print("  n        = %.8f" % n)
print("  w_phi    = %.12f  (ya es w0)" % w_phi)
r = pi_sobre_crit / OM_DE
w_tot = w_phi + r
print("  Pi/rho_DE= %+.6f" % r)
print("  w_eff tot= %+.6f" % w_tot)
print("  objetivo = %+.6f" % W0)
print()
if abs(w_tot - W0) > 0.01:
    print("  [FALLA] regimen constante")
else:
    print("  [OK]")

print()
print("  CONTROL: y el regimen (1-a)?")
IG = S.IGNIS
G_FIX = S.P_SC / (3.0 * IG)
print("  Gamma_fix= %+.8f" % G_FIX)
print("  Pi/rho_DE en a->0 = %+.8f"
      % (-3.0 * G_FIX))
print("  w_a resultante    = %+.8f"
      % (-S.P_SC / IG))
