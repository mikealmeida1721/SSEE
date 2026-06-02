#!/usr/bin/env python3
"""Look-elsewhere sobre la FAMILIA SSEE REAL (diccionario cerrado).

El conteo previo barria cualquier (a*phi+b*pi)/(c*phi+d*pi) — espacio abierto,
denominador equivocado. La familia SSEE real es el diccionario FIJO de
constantes con nombre de sandbox_unificado/lib/ssee.py (Genesis 5.12).

Aqui el test honesto: de TODAS las razones A/B entre constantes con nombre,
cuantas caen en 0.840 (=w0) y 0.670 (=wa). Ese es el look-elsewhere fiel al
modelo: no inventamos vocabulario, usamos el que ya existia.
"""
import math
from itertools import permutations

PHI = 1.618033988749895
PI  = 3.141592653589793
OMEGA = PHI + PI

BIAL    = (PHI + PI) / 2
KAL     = BIAL + PI
SOLAR   = BIAL + KAL
MAR     = OMEGA + PI
VITA    = PI + KAL
ANMA    = BIAL + VITA
PYROS   = OMEGA + PHI
IGNIS   = PI + PYROS
KRYSTOS = PHI + PI + OMEGA
MIKA    = KRYSTOS + PHI
AURA    = PHI + BIAL
MIRA    = AURA / 2
DUAL    = AURA * 2
TRIAL   = AURA * 3
CUARTAL = AURA * 4
MIKAEL_V = PHI + PI + KRYSTOS
BUFFER  = MIKAEL_V - TRIAL
KRYSTOS_V = 2 * OMEGA

FAMILY = {
    "PHI": PHI, "PI": PI, "OMEGA": OMEGA, "BIAL": BIAL, "KAL": KAL,
    "SOLAR": SOLAR, "MAR": MAR, "VITA": VITA, "ANMA": ANMA, "PYROS": PYROS,
    "IGNIS": IGNIS, "KRYSTOS": KRYSTOS, "MIKA": MIKA, "AURA": AURA,
    "MIRA": MIRA, "DUAL": DUAL, "TRIAL": TRIAL, "CUARTAL": CUARTAL,
    "MIKAEL_V": MIKAEL_V, "BUFFER": BUFFER, "KRYSTOS_V": KRYSTOS_V,
}

W0 = TRIAL / MIKAEL_V      # 0.83995
WA = PYROS / KRYSTOS       # 0.66997
print(f"w0 = TRIAL/MIKAEL_V = {W0:.9f}")
print(f"wa = PYROS/KRYSTOS  = {WA:.9f}")
print(f"constantes en la familia: {len(FAMILY)}\n")

# Razones unicas A/B (A != B por valor). Dedup por valor redondeado.
names = list(FAMILY)
ratios = []          # (valor, "A/B")
seen = set()
for a, b in permutations(names, 2):
    va, vb = FAMILY[a], FAMILY[b]
    if abs(vb) < 1e-9:
        continue
    r = va / vb
    if r <= 0 or r > 5:
        continue
    key = round(r, 6)
    if key in seen:
        continue
    seen.add(key)
    ratios.append((r, f"{a}/{b}"))

total = len(ratios)
print(f"razones A/B distintas en (0,5]: {total}\n")

for label, tgt in (("w0", W0), ("wa", WA)):
    for tol in (0.01, 0.005, 0.001):
        hits = sorted((r, s) for (r, s) in ratios if abs(r - tgt) <= tol)
        print(f"  {label}={tgt:.6f}  tol +-{tol:<6}: {len(hits):3d} razones")
        for r, s in sorted(hits, key=lambda x: abs(x[0]-tgt))[:5]:
            print(f"        {r:.6f}  {s:<22} |d|={abs(r-tgt):.2e}")
    print()
