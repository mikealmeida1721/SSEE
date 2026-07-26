#!/usr/bin/env python3
"""Look-elsewhere / numerologia audit para w0 y wa.

Pregunta de Mike (2026-06-01): ¿cuantas combinaciones de phi y pi dan
0.840 y 0.670? Si hay muchas -> el match a DESI es cherry-picking.
Si son escasas -> robusto.

No opinamos: contamos. Buscamos sobre razones (a*phi+b*pi)/(c*phi+d*pi)
con coeficientes enteros pequenos, y reportamos cuantas caen cerca de
cada objetivo a varias tolerancias.
"""
import math
from fractions import Fraction

PHI = (1.0 + math.sqrt(5.0)) / 2.0
PI  = math.pi

# Objetivos canonicos (magnitudes; los signos son fisica, no numerologia)
OMEGA = PHI + PI
BETA  = (PHI + PI) / 2.0
P_SC  = OMEGA + PHI
K_V   = PHI + PI + OMEGA
T_R   = 3.0 * (PHI + BETA)
M_V   = PHI + PI + K_V
W0_TARGET = T_R / M_V          # 0.8399...
IGNIS = PI + P_SC              # π+PYROS — denominador de wₐ (R21)
WA_TARGET = P_SC / IGNIS       # 0.6699...

print(f"objetivo w0 = {W0_TARGET:.9f}")
print(f"objetivo wa = {WA_TARGET:.9f}\n")

# --- Espacio de busqueda: razones de combos lineales enteros de phi, pi ---
# numerador  = a*phi + b*pi
# denominador= c*phi + d*pi
# coeficientes en COEF (incluye 0). Excluimos combos nulos / degenerados.
def run(coef_max, tols, targets):
    coefs = range(-coef_max, coef_max + 1)
    seen_values = []   # lista de (valor, (a,b,c,d)) unicos por valor redondeado
    seen_keys = set()
    for a in coefs:
        for b in coefs:
            num = a*PHI + b*PI
            if abs(num) < 1e-9:
                continue
            for c in coefs:
                for d in coefs:
                    den = c*PHI + d*PI
                    if abs(den) < 1e-9:
                        continue
                    val = num / den
                    if val <= 0:           # buscamos magnitudes positivas
                        continue
                    if val > 5:            # rango razonable
                        continue
                    key = round(val, 6)
                    if key in seen_keys:
                        continue
                    seen_keys.add(key)
                    seen_values.append((val, (a, b, c, d)))
    total = len(seen_values)
    print(f"--- coef en [-{coef_max},{coef_max}] : {total} valores distintos en (0,5] ---")
    for name, tgt in targets:
        for tol in tols:
            hits = [(v, q) for (v, q) in seen_values if abs(v - tgt) <= tol]
            print(f"  {name}={tgt:.6f}  tol +-{tol:<6}: {len(hits):4d} combos")
            # mostrar los 3 mas cercanos
            hits_sorted = sorted(hits, key=lambda x: abs(x[0]-tgt))[:3]
            for v, (a, b, c, d) in hits_sorted:
                print(f"        {v:.6f}  ({a}*phi+{b}*pi)/({c}*phi+{d}*pi)  |d|={abs(v-tgt):.2e}")
        print()

tols = [0.01, 0.005, 0.001]
targets = [("w0", W0_TARGET), ("wa", WA_TARGET)]

run(3, tols, targets)
run(5, tols, targets)
