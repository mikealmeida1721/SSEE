"""
A_s — ¿puede SSEE derivar su propia amplitud escalar, o es el abismo de escala?
==============================================================================
Mike (2026-06-01): el CMB de SSEE pide prestado A_s=2.1e-9 a Planck. A_s es
ADIMENSIONAL, así que la intuición dice "SSEE debería poder sacarlo solo".

Test honesto. En α-attractor (el potencial que SSEE YA usa para n_s y r):
        A_s = V₀ / (24·π²·M_Pl⁴·ε)
SSEE ya tiene LIMPIO:
        ε = r/16 = φ⁻¹⁰/16     (de OP-2, sin pedir nada prestado)
Despejando, A_s queda fijo SI Y SOLO SI conocemos la altura del potencial:
        V₀/M_Pl⁴ = 24·π²·ε·A_s
Eso es un número ~1e-10. PREGUNTA DISCRIMINANTE:
  (1) ¿V₀/M_Pl⁴ sale como monomio (φ,π) de grado chico? → A_s derivable.
  (2) ¿o necesita φ^(−N) con N enorme? → es la ESCALA de inflación disfrazada,
      el mismo abismo 10⁻⁶¹ que el Postulado D admite que nadie cruza.
  (3) control transmutación: ¿ln(M_Pl⁴/V₀) = f(φ,π) limpio de grado chico?
Cero maquillaje: si solo salen flukes de alta densidad, A_s NO es de SSEE.
"""
import numpy as np
import itertools
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import PHI, PI, OMEGA, BETA, KAL0, MIRA, AURA, T_R, M_V

A_S   = 2.1e-9
r     = PHI**-10
eps   = r / 16.0
ratio = 24.0 * PI**2 * eps * A_S          # V₀/M_Pl⁴ requerido

print("=" * 74)
print(" A_s — ¿derivable de (φ,π) o escala de inflación disfrazada?")
print("=" * 74)
print(f"  r = φ⁻¹⁰            = {r:.6e}")
print(f"  ε = r/16 (limpio)  = {eps:.6e}")
print(f"  A_s objetivo       = {A_S:.4e}")
print(f"  ⟹ V₀/M_Pl⁴ exigido = {ratio:.6e}   (log10 = {np.log10(ratio):.3f})")

# ── paso 1: ¿qué φ⁻N daría ese ratio? exponente chico = limpio, enorme = numerología
n_phi = -np.log(ratio) / np.log(PHI)
print(f"\n── paso 1: lineal — ¿qué potencia de φ da V₀/M_Pl⁴? ──")
print(f"  necesitaría φ^(−{n_phi:.1f})  "
      f"({'razonable' if n_phi < 12 else 'ENORME ⟹ escala disfrazada, no limpio'})")

# ── paso 2: monomio (φ,π) directo a <0.2% + densidad combinatoria ──────────
print("\n── paso 2: ¿V₀/M_Pl⁴ es monomio (φ,π) grado≤6 a <0.2%? ──")
basis = {'φ':PHI,'π':PI,'Ω':OMEGA,'β':BETA,'KAL':KAL0,'MIRA':MIRA,'AURA':AURA}
names = list(basis); vals = np.array([basis[n] for n in names])
hits = []; dens = 0
for combo in itertools.product(range(-6, 7), repeat=len(names)):
    s = sum(abs(c) for c in combo)
    if s == 0 or s > 6:
        continue
    val = np.prod([vals[i]**combo[i] for i in range(len(names))])
    if val <= 0:
        continue
    err = abs(val - ratio) / ratio
    if err < 0.05:
        dens += 1
    if err < 0.002:
        expr = "·".join(f"{names[i]}^{combo[i]}" for i in range(len(names)) if combo[i])
        hits.append((err, expr, val, s))
hits.sort()
if hits:
    print(f"  candidatos a <0.2%:")
    for err, expr, val, s in hits[:6]:
        print(f"     {expr:<30} = {val:.4e}  ({err*100:+.3f}%, grado {s})")
    print(f"  densidad (#en ±5%) = {dens}")
    if dens > 20:
        print(f"  ⚠ densidad ALTA → match esperado por combinatoria, NO señal")
else:
    print(f"  NINGÚN monomio (φ,π) grado≤6 cae a <0.2% de {ratio:.3e}")
    print(f"  densidad (#en ±5%) = {dens}")

# ── paso 3: transmutación — ¿ln(M_Pl⁴/V₀) limpio? (como Λ_QCD) ────────────
L = -np.log(ratio)        # = ln(M_Pl⁴/V₀)
print(f"\n── paso 3: transmutación — ln(M_Pl⁴/V₀) = {L:.4f} ──")
hitsL = []; densL = 0
for combo in itertools.product(range(-3, 4), repeat=len(names)):
    s = sum(abs(c) for c in combo)
    if s == 0 or s > 4:
        continue
    val = np.prod([vals[i]**combo[i] for i in range(len(names))])
    if val <= 0:
        continue
    err = abs(val - L) / L
    if err < 0.05:
        densL += 1
    if err < 0.002:
        expr = "·".join(f"{names[i]}^{combo[i]}" for i in range(len(names)) if combo[i])
        hitsL.append((err, expr, val, s))
hitsL.sort()
if hitsL:
    for err, expr, val, s in hitsL[:6]:
        print(f"     {expr:<30} = {val:.4f}  ({err*100:+.3f}%, grado {s})")
    print(f"  densidad (#en ±5%) = {densL}")
    if densL > 30:
        print(f"  ⚠ densidad ALTA → match esperado por combinatoria, NO señal")
else:
    print(f"  NINGÚN monomio (φ,π) grado≤4 a <0.2% de L={L:.4f}")
    print(f"  densidad (#en ±5%) = {densL}")

print("\n" + "=" * 74)
print(" LECTURA")
print("=" * 74)
print(f"""
  Si paso 1 da φ^(−N) con N enorme y paso 2/3 solo dan flukes de alta
  densidad → A_s NO es de SSEE: es la escala de inflación V₀/M_Pl⁴
  disfrazada de adimensional. Mismo abismo 10⁻⁶¹ del Postulado D.
  Honesto: A_s queda como INPUT (o predicción condicional-a-M_Pl si algún
  día V₀ se fija con un ancla). NO es "0 parámetros".

  Si por milagro paso 2 o 3 da algo LIMPIO y de baja densidad → A_s pasa a
  ser derivado, y el CMB de SSEE deja de pedirle A_s prestado a Planck.
""")
