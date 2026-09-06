"""OP-23: existe K(X) con
   w_phi = -T_r/M_v  y  c_s^2 > 0 ?

Forma minima de dos terminos:
   K(X) = c1 X + c2 X^2
"""
import sys, os
import sympy as sp
sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..'))
import ssee_core as S

X, c1, c2 = sp.symbols('X c1 c2')
K = c1 * X + c2 * X ** 2
KX = sp.diff(K, X)
KXX = sp.diff(K, X, 2)

rho = sp.simplify(2 * X * KX - K)
w = sp.simplify(K / rho)
gh = sp.simplify(KX + 2 * X * KXX)
cs2 = sp.simplify(KX / gh)

u = sp.symbols('u')
sub = {c2: u * c1 / X}
wu = sp.simplify(w.subs(sub))
cu = sp.simplify(cs2.subs(sub))
print("  con u = c2 X/c1")
print("  w     =", wu)
print("  c_s^2 =", cu)

wsym = sp.symbols('w')
usol = sp.solve(sp.Eq(wu, wsym), u)[0]
csw = sp.simplify(cu.subs(u, usol))
print()
print("  u(w)     =", sp.simplify(usol))
print("  c_s^2(w) =", csw)

print()
print("  NUMERICO en SSEE")
W0 = S.W0
TR, MV = S.T_R, S.M_V
uv = (1 - W0) / (3 * W0 - 1)
cv = (1 + W0) / (5 - 3 * W0)
print("  w0     = %+.12f" % W0)
print("  u      = %+.12f" % uv)
print("  c_s^2  = %+.12f" % cv)
print()
alg = (MV - TR) / (5 * MV + 3 * TR)
print("  forma algebraica")
print("  (M_v-T_r)/(5M_v+3T_r)")
print("       = %.12f" % alg)
print("  dif  = %.2e" % abs(alg - cv))
print()
print("  1+w0 = %.12f" % (1 + W0))
print("  (M_v-T_r)/M_v = %.12f"
      % ((MV - TR) / MV))
print()
print("  CONDICIONES DE SIGNO")
print("  1+3u = %+.8f" % (1 + 3 * uv))
print("  1+6u = %+.8f" % (1 + 6 * uv))
print("  1+2u = %+.8f" % (1 + 2 * uv))
print("  rho>0 y sin ghost piden c1<0")
print()
aK = (5 - 3 * W0) * 3 * (1 - S.OMEGA_M_CMB)
print("  alpha_K = (5-3w)*3*Om_DE")
print("          = %+.6f" % aK)
print("  P7 dice   +0.403300")

print()
print("  VERIFICACION DIRECTA")
c1v = -1.0 / S.KAL0
Xv = 1.0
c2v = uv * c1v / Xv
Kv = c1v * Xv + c2v * Xv ** 2
KXv = c1v + 2 * c2v * Xv
KXXv = 2 * c2v
rv = 2 * Xv * KXv - Kv
pv = Kv
ghv = KXv + 2 * Xv * KXXv
csv = KXv / ghv
ok = True
pruebas = (
    ("rho > 0", rv, rv > 0),
    ("rho+p > 0 (no fantasma)",
     rv + pv, rv + pv > 0),
    ("K_X+2XK_XX > 0 (no ghost)",
     ghv, ghv > 0),
    ("c_s^2 > 0", csv, csv > 0),
    ("w = w0", pv / rv,
     abs(pv / rv - W0) < 1e-12),
)
for nom, val, cond in pruebas:
    est = "[OK]" if cond else "[FALLA]"
    ok = ok and cond
    print("  %-28s %+12.8f %s"
          % (nom, val, est))
print()
print("  VEREDICTO:",
      "EXISTE" if ok else "NO EXISTE")
