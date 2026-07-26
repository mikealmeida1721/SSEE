#!/usr/bin/env python3
"""[SUPERSEDIDO 2026-07-25] conteo conjunto por barrido de enteros.

Este script barre razones (a·φ+b·π)/(c·φ+d·π) con coeficientes ENTEROS —
el espacio abierto. El PRD §2.3 NO usa este conteo y explícitamente lo
descarta: «a stronger statement than a generic (aφ+bπ)/(cφ+dπ) scan, which
by construction inflates the denominator with expressions the model never
uses». El conteo VIGENTE es sobre el diccionario CERRADO de 55 constantes
con nombre → `look_elsewhere_full.py` (55/25/490, 1-de-490 a ±0.0005, R27).
Se conserva como historia de cómo se llegó al denominador correcto.


El conteo suelto (look_elsewhere_w0wa.py) pregunta cuantas razones dan 0.840
y, por separado, cuantas dan 0.670. Eso sobreestima el azar: trata w0 y wa
como dos elecciones independientes.

SSEE NO elige por separado. Construye UNA familia: un esqueleto comun de
constantes derivadas de phi y pi, del cual salen los DOS numeros a la vez.
Esqueleto NULO (álgebra convencional): M_v = 3*Omega, K_v = 2*Omega, con Omega =
phi+pi. (En álgebra SSEE el linaje es K_v = KRYSTOS_V = phi+pi+Omega, padres {φ,π,Ω},
NO 2Ω; aquí el múltiplo entero es la hipótesis nula contra la que se mide.) Los dos
denominadores son el MISMO Omega por un entero.

Aqui contamos familias: fijamos una base compartida D = (c*phi + d*pi) y
preguntamos si CON ESA MISMA BASE existen numeradores enteros que den
~0.840 (sobre k1*D) y ~0.670 (sobre k2*D). Solo cuentan las bases que
aciertan AMBOS objetivos. Ese es el denominador correcto del look-elsewhere.
"""
import math

PHI = (1.0 + math.sqrt(5.0)) / 2.0
PI  = math.pi

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
print(f"objetivo wa = {WA_TARGET:.9f}")
print(f"check: M_v/Omega = {M_V/OMEGA:.6f} (=3?),  K_v/Omega = {K_V/OMEGA:.6f} (=2?)\n")


def joint(coef_max, num_max, kmax, tol):
    """Cuenta familias con base compartida D=(c*phi+d*pi) que aciertan AMBOS.

    coef_max : rango de la BASE compartida (c,d)
    num_max  : rango de los numeradores enteros (a,b)
    kmax     : multiplo entero del esqueleto (1..kmax), p.ej. M_v=3*Omega usa k=3
    """
    base_coefs = range(-coef_max, coef_max + 1)
    num_coefs  = range(-num_max, num_max + 1)
    ks = range(1, kmax + 1)

    joint_bases = []   # (D, base_coefs) que aciertan ambos
    seen_base = set()

    for c in base_coefs:
        for d in base_coefs:
            D = c * PHI + d * PI
            if abs(D) < 1e-9:
                continue
            keyD = round(abs(D), 6)
            if keyD in seen_base:
                continue
            seen_base.add(keyD)

            hit_w0 = False
            hit_wa = False
            for k in ks:
                den = abs(k * D)
                if den < 1e-9:
                    continue
                for a in num_coefs:
                    for b in num_coefs:
                        num = a * PHI + b * PI
                        if num <= 0:
                            continue
                        val = num / den
                        if abs(val - W0_TARGET) <= tol:
                            hit_w0 = True
                        if abs(val - WA_TARGET) <= tol:
                            hit_wa = True
                    if hit_w0 and hit_wa:
                        break
                if hit_w0 and hit_wa:
                    break
            if hit_w0 and hit_wa:
                joint_bases.append((abs(D), (c, d)))

    total_bases = len(seen_base)
    print(f"--- base coef [-{coef_max},{coef_max}], numer [-{num_max},{num_max}], "
          f"k in 1..{kmax}, tol +-{tol} ---")
    print(f"    bases compartidas distintas probadas : {total_bases}")
    print(f"    bases que aciertan AMBOS (familias)   : {len(joint_bases)}")
    # destacar la base canonica (Omega = phi+pi -> c=1,d=1)
    omega_in = any(abs(D - OMEGA) < 1e-6 for D, _ in joint_bases)
    print(f"    incluye la base canonica Omega=phi+pi : {omega_in}")
    for D, (c, d) in sorted(joint_bases)[:8]:
        print(f"        D={D:.4f}  ({c}*phi+{d}*pi)")
    print()


# Tolerancia realista: el match a DESI es sub-porcentual.
for tol in (0.005, 0.001):
    joint(coef_max=3, num_max=3, kmax=3, tol=tol)
    joint(coef_max=5, num_max=5, kmax=4, tol=tol)
