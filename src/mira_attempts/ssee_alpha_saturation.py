"""
ssee_alpha_saturation.py — Paso A: derivación simbólica de las ecuaciones de
movimiento del sistema φ + materia + radiación con acoplamiento conformal.

Objetivo (Pasos A→B→C, este script implementa solo A):

  A. Escribir Friedmann + Klein-Gordon + conservación de materia con la acción
     de Paper 7 + acoplamiento conformal C(φ) = exp(α φ/M_pl), sin asumir el
     valor de α.

  B. Reducir a sistema autónomo en variables adimensionales (x, y, z, Ω_m).
     Encontrar puntos fijos (atractor tardío).

  C. Linealizar alrededor del atractor. Imponer Re(autovalor) < 0 para todos.
     Eso da una desigualdad explícita en α. La SATURACIÓN de esa desigualdad
     fija α — y solo entonces lo comparamos con expresiones algebraicas de φ,π.

Este script ejecuta SOLO el Paso A: deriva las ecuaciones y las imprime.
No asume α saturado; no asume nada de SSEE específico. Pura mecánica.

Convención: signatura (-,+,+,+); X = (1/2) g^μν ∂_μφ ∂_νφ (signo cosmológico
positivo para φ̇²/2 en background FRW).
"""

import sympy as sp

# ──────────────────────────────────────────────────────────────────
# Símbolos
# ──────────────────────────────────────────────────────────────────
# Funciones del tiempo cósmico t
t = sp.Symbol('t', real=True, positive=True)
a = sp.Function('a')(t)              # factor de escala
phi = sp.Function('phi')(t)          # campo escalar
rho_m = sp.Function('rho_m')(t)      # densidad de materia (DM dominante)
rho_r = sp.Function('rho_r')(t)      # densidad de radiación (desacoplada)

# Constantes
M_pl = sp.Symbol('M_pl', positive=True)   # masa de Planck reducida
alpha = sp.Symbol('alpha', real=True)     # coeficiente del acoplamiento (LO QUE BUSCAMOS)
KAL = sp.Symbol('KAL', positive=True)     # escala k-essence lineal (SSEE: KAL_0)
M4 = sp.Symbol('M^4', positive=True)      # escala k-essence cuadrática (SSEE: M)

# Potencial V(φ) — dejamos forma genérica simbólica
V = sp.Function('V')(phi)
Vp = sp.diff(V, phi)                 # V'(φ)

# Hubble
H = sp.diff(a, t) / a

# Derivadas del campo
phi_dot = sp.diff(phi, t)
phi_ddot = sp.diff(phi, t, 2)

# X = φ̇²/2  (homogéneo, signo positivo en FRW)
X = phi_dot**2 / 2

# K-essence SSEE: K(X) = X/KAL + X²/M⁴
K = X/KAL + X**2/M4
K_X = sp.diff(K, X) if False else sp.Rational(1,1)/KAL + 2*X/M4  # ∂K/∂X
K_XX = sp.Rational(2,1)/M4                                        # ∂²K/∂X²

# Densidad y presión del campo (k-essence)
#   ρ_φ = 2X K_X − K + V    p_φ = K − V
rho_phi = 2*X*K_X - K + V
p_phi   = K - V

# Acoplamiento conformal: C(φ) = exp(α φ / M_pl) — Wetterich/Amendola
# El intercambio energético entre φ y materia es:
#   Q = (α / M_pl) ρ_m φ̇   (signo según convención Amendola)
Q = (alpha / M_pl) * rho_m * phi_dot

# ──────────────────────────────────────────────────────────────────
# Ecuaciones de movimiento
# ──────────────────────────────────────────────────────────────────

print("="*72)
print("PASO A — Derivación simbólica de las ecuaciones de movimiento")
print("="*72)

# (1) Friedmann
#     3 M_pl² H² = ρ_φ + ρ_m + ρ_r
Friedmann = sp.Eq(3*M_pl**2 * H**2, rho_phi + rho_m + rho_r)
print("\n(1) Ecuación de Friedmann:")
sp.pprint(Friedmann, use_unicode=True)

# (2) Conservación de materia con acoplamiento
#     ρ̇_m + 3H ρ_m = Q
matter_cons = sp.Eq(sp.diff(rho_m, t) + 3*H*rho_m, Q)
print("\n(2) Conservación de materia (acoplada):")
sp.pprint(matter_cons, use_unicode=True)

# (3) Conservación de radiación (desacoplada)
#     ρ̇_r + 4H ρ_r = 0
rad_cons = sp.Eq(sp.diff(rho_r, t) + 4*H*rho_r, 0)
print("\n(3) Conservación de radiación (desacoplada):")
sp.pprint(rad_cons, use_unicode=True)

# (4) Klein-Gordon acoplado
#     d/dt[2X K_X − K] − 3H · 2X·K_X ... más correcto: de la conservación
#     del 4-momento total ∇_μ T^μν_total = 0 con Q como fuente:
#
#     ρ̇_φ + 3H(ρ_φ + p_φ) = -Q
#
# Que en términos del campo se reduce a:
#
#     (K_X + 2X K_XX) φ̈ + 3H K_X φ̇ + V'(φ) = -(α/M_pl) ρ_m
#
# (con signo opuesto de Q porque φ pierde lo que materia gana)
KG = sp.Eq(
    (K_X + 2*X*K_XX) * phi_ddot + 3*H*K_X*phi_dot + Vp,
    -(alpha/M_pl) * rho_m
)
print("\n(4) Klein-Gordon acoplado (φ pierde lo que materia gana):")
sp.pprint(KG, use_unicode=True)

# ──────────────────────────────────────────────────────────────────
# Variables adimensionales (Copeland-Liddle-Wands 1998)
# ──────────────────────────────────────────────────────────────────
#
# Para el sistema autónomo introducimos:
#   x  ≡ φ̇ / (√6 H M_pl)        (variable cinética)
#   y  ≡ √V / (√3 H M_pl)        (variable potencial)
#   z  ≡ √ρ_r / (√3 H M_pl)      (variable de radiación)
#   N  ≡ ln a                     (número de e-folds como tiempo)
#
# Entonces:
#   Ω_φ_canon = x² + y²        (canónico — se modifica con K(X) no trivial)
#   Ω_r       = z²
#   Ω_m       = 1 − x² − y² − z²
#
# Para K(X) = X/KAL + X²/M⁴, el peso cinético es modificado. La definición
# correcta de Ω_φ involucra (2X K_X − K)/(3 H² M_pl²). Para el Paso B
# trabajaremos primero el caso canónico (KAL→∞-equivalente) para reproducir
# Amendola 1999, y luego añadiremos el término X²/M⁴.

print("\n" + "="*72)
print("PASO A — Variables adimensionales (CLW 1998 + extensión SSEE)")
print("="*72)
print("""
  x ≡ φ̇ / (√6 H M_pl)     [cinética]
  y ≡ √V / (√3 H M_pl)     [potencial]
  z ≡ √ρ_r / (√3 H M_pl)   [radiación]
  N ≡ ln a                  [tiempo]

  Ω_m = 1 − (peso_K(x)) − y² − z²

  λ ≡ −M_pl · V'/V          [pendiente del log-potencial]
""")

print("\n" + "="*72)
print("NOTA importante sobre lo que viene")
print("="*72)
print("""
Próximo paso (B): pasar las ecuaciones (1)-(4) a las variables (x, y, z),
escribir el sistema autónomo dx/dN, dy/dN, dz/dN, y resolver puntos fijos.

Para Amendola 1999 (caso canónico K=X), el resultado conocido es que el
atractor de materia con acoplamiento conformal existe y es estable solo si:

      α² < 3/2         (esta es la cota canónica)

La SATURACIÓN canónica daría α = √(3/2) ≈ 1.2247.

Para SSEE — con K(X) = X/KAL + X²/M⁴ — esa cota se modifica. La modificación
depende del régimen (lineal X/KAL vs cuadrático X²/M⁴ dominante). Eso es lo
que el Paso B tiene que calcular limpio.

Honestidad: si la cota saturada SSEE coincide con alguna expresión algebraica
de φ,π → veta-2 funciona. Si NO coincide → SSEE no tiene mecanismo de
acoplamiento por saturación, y se sabe limpiamente.
""")
