"""
OP-14 — Σm_ν Phenomenological Derivation Audit
================================================

Hipótesis a testar:
  (1) FRAGILIDAD: R = {4·KAL} es coincidencia "casi-entera" sin contenido físico
  (2) SATURACIÓN: existe cota física con valor entero (22) que R mide como exceso
  (3) FORMA ALGEBRAICA: R tiene representación φ,π más natural que 4·KAL-22

Estrategia:
  - Test de robustez: perturbar φ, π → ¿colapsa la predicción?
  - Scan ~150 monomios algebraicos en {φ, π, Ω, β, KAL, P_sc, K_v, T_r, M_v}
  - Test de la hipótesis DoF: ¿22 = conteo del SM (15 fermiones + 4 gauge + 1 Higgs + 2 grav)?
  - Reformulación: ¿se puede expresar Σm_ν directamente sin pasar por R?
"""
import numpy as np
from itertools import product

# Constantes algebraicas SSEE
PHI  = (1 + np.sqrt(5)) / 2
PI   = np.pi
OMEGA = PI + PHI                # 4.7596
BETA  = (PI + PHI) / 2          # 2.3798
KAL   = BETA + PI               # 5.5214 = (φ + 3π)/2
P_SC  = OMEGA + PHI             # 6.3776
K_V   = PHI + PI + OMEGA        # 9.5192
T_R   = 3 * (PHI + BETA)        # 11.9935
M_V   = PHI + PI + K_V          # 14.2789

# Valores observacionales
R_observed     = 4 * KAL - 22           # 0.0856...
Sigma_mnu_obs  = 0.0824                  # eV (Paper 4)
Omega_b_h2     = 0.02237
H0_alg         = 3 * (PHI + PI)**2 / 100 # Mpc^-1 dimensionless h
tau_Pi_H0      = KAL / (3 * 0.840)       # 2.191

print("=" * 72)
print("OP-14 AUDIT — Σm_ν phenomenological derivation")
print("=" * 72)
print(f"\nConstantes algebraicas:")
print(f"  φ            = {PHI:.10f}")
print(f"  π            = {PI:.10f}")
print(f"  KAL          = (φ+3π)/2 = {KAL:.10f}")
print(f"  4·KAL        = 2φ + 6π = {4*KAL:.10f}")
print(f"  R = 4·KAL−22 = {R_observed:.10f}")
print(f"  Σm_ν target  = R·Ωbh²·94.07/τΠH₀ = {R_observed * Omega_b_h2 * 94.07 / tau_Pi_H0:.6f} eV")
print(f"  Σm_ν Paper 4 = {Sigma_mnu_obs} eV")

# ============================================================
# HIPÓTESIS 1 — FRAGILIDAD: ¿el resultado depende del "casi-entero"?
# ============================================================
print("\n" + "=" * 72)
print("HIPÓTESIS 1 — Test de fragilidad")
print("=" * 72)

# Si perturbamos φ o π por una fracción pequeña, ¿qué le pasa a R?
# Una predicción "estructural" debería ser estable; una "coincidencia
# numérica" colapsa porque 4·KAL podría caer del otro lado del entero.

for delta_frac in [1e-6, 1e-4, 1e-3, 1e-2]:
    phi_pert = PHI * (1 + delta_frac)
    pi_pert  = PI  * (1 + delta_frac)
    KAL_pert = (phi_pert + 3*pi_pert) / 2
    R_pert   = 4*KAL_pert - 22
    Sigma_pert = R_pert * Omega_b_h2 * 94.07 / tau_Pi_H0
    drift = (Sigma_pert / Sigma_mnu_obs - 1) * 100
    print(f"  δ(φ,π)/x = {delta_frac:.0e}  → R={R_pert:.5f}  Σm_ν={Sigma_pert:.5f} eV  ({drift:+.1f}%)")

# Conclusión esperada: el resultado deriva linealmente con φ, π — NO es
# infinitamente frágil. La razón: 4·KAL es ~22.086, está suficientemente
# lejos del entero como para que perturbaciones pequeñas no crucen.
# Pero verifiquemos.

# ============================================================
# HIPÓTESIS 2 — SATURACIÓN: ¿22 = conteo físico SM?
# ============================================================
print("\n" + "=" * 72)
print("HIPÓTESIS 2 — ¿22 es un conteo de DoF físico?")
print("=" * 72)

# SM bosonic + fermionic DoF (sin contar antipartícula):
#   Fermiones: 6 quarks × 2(L,R) × 3(color) = 36, + 6 leptones × 2 = 12 → 48? No, demasiado
#   Gauge: 8 gluones + 3 W + 1 B + 1 Higgs = 13 bosones
#   Más conservador: contar familias completas: 3 gen × ...
#
# Otro conteo: SM tiene 28 parámetros libres totales. Demasiado.
#
# Conteo de generadores: SU(3)×SU(2)×U(1) → 8+3+1 = 12 generadores gauge.
# No coincide con 22.
#
# 11 = dim M-theory; 22 = 2 × 11. ¿Coincidencia?
# Otra opción: 22 = rank de E11 (álgebra de Kac-Moody) — sí, rank(E11) = 11,
# pero 2·rank no es estándar.

candidates_22 = {
    "2 × 11 (M-theory dim doblada)": 22,
    "rank(E11) × 2": 22,
    "DoF fermiónicos quarks (6×2×3 sin color factor)": 36,  # no match
    "Generadores SM (8+3+1)": 12,                            # no match
    "Bosones SM (8 gluones+W±+Z+γ+H+graviton+...)": None,
    "Higgs doublet + 3 gauge fields acoplados": None,
    "(N_f + N_g) algún conteo": None,
}

# El conteo "22 = 2·11" es sugerente pero no constructivo. Sin un argumento
# físico independiente de por qué SSEE debería verse en M-theory o E11, esto
# es post-hoc.

print("  22 = 2 × 11 → consistente con dim M-theory (11) doblada, pero post-hoc")
print("  22 ≠ Generadores SM (12), ni 28 (params libres SM)")
print("  Veredicto preliminar: no hay match limpio con conteo SM canónico")

# ============================================================
# HIPÓTESIS 3 — SCAN ALGEBRAICO
# ============================================================
print("\n" + "=" * 72)
print("HIPÓTESIS 3 — Scan: ¿existe forma algebraica limpia para R = 0.0856?")
print("=" * 72)

target = R_observed
TOL_ACCEPT = 1e-3  # 0.1% match para considerar candidato
constants = {
    "φ": PHI, "π": PI, "Ω": OMEGA, "β": BETA,
    "KAL": KAL, "Psc": P_SC, "Kv": K_V, "Tr": T_R, "Mv": M_V,
}

candidates = []

# Categoría A: x^n con n ∈ [-10, 10]
for name, val in constants.items():
    for n in range(-10, 11):
        if n == 0: continue
        v = val ** n
        if abs(v - target) / target < TOL_ACCEPT * 100:  # 10% para shortlist
            err = (v - target) / target * 100
            candidates.append((f"{name}^{n}", v, err))

# Categoría B: 1/(a+b·c) con a, b, c en constants
for (n1, v1), (n2, v2) in product(constants.items(), repeat=2):
    for k in [1, 2, 3]:
        for op_name, op_val in [("+", v1+v2), ("-", v1-v2), ("·", v1*v2)]:
            if op_val == 0: continue
            v_test = 1 / (k * op_val)
            if abs(v_test) < 1.0 and abs(v_test - target) / target < 0.1:
                err = (v_test - target) / target * 100
                candidates.append((f"1/({k}·({n1}{op_name}{n2}))", v_test, err))

# Categoría C: (a−b)/c (típico SSEE)
for (n1, v1), (n2, v2), (n3, v3) in product(constants.items(), repeat=3):
    if v3 == 0: continue
    for k_pow in [1, 2, 3]:
        denom = v3 ** k_pow
        if denom == 0: continue
        v_test = (v1 - v2) / denom
        if 0 < v_test < 1 and abs(v_test - target) / target < 0.05:
            err = (v_test - target) / target * 100
            candidates.append((f"({n1}−{n2})/{n3}^{k_pow}", v_test, err))

# Categoría D: φ^n / π^m
for n in range(-6, 7):
    for m in range(-6, 7):
        if n == 0 and m == 0: continue
        v = PHI**n * PI**m
        if 0 < v < 1 and abs(v - target) / target < 0.05:
            err = (v - target) / target * 100
            candidates.append((f"φ^{n}·π^{m}", v, err))

# Ordenar por error absoluto
candidates.sort(key=lambda x: abs(x[2]))

print(f"\n  Mejores candidatos (≤5% de error):")
print(f"  {'forma':<40} {'valor':<14} {'err %':>8}")
print(f"  {'-'*40} {'-'*14} {'-'*8}")
seen = set()
shown = 0
for name, val, err in candidates:
    if name in seen: continue
    seen.add(name)
    if abs(err) <= 5.0:
        print(f"  {name:<40} {val:<14.7f} {err:>+8.4f}")
        shown += 1
    if shown >= 15: break

if shown == 0:
    print("  → NINGÚN match dentro del 5%. La forma 'R = 4·KAL−22' no tiene")
    print("    representación algebraica obvia más natural.")

# ============================================================
# REFORMULACIÓN — ¿se puede saltar R y derivar Σm_ν directamente?
# ============================================================
print("\n" + "=" * 72)
print("REFORMULACIÓN — Atacar Σm_ν directamente (sin pasar por R)")
print("=" * 72)

# Σm_ν = 0.0824 eV
# H_0 = 67.96 km/s/Mpc → en eV: H₀ ≈ 1.435 × 10⁻³³ eV
# Σm_ν / H_0 ≈ 0.0824 / 1.435e-33 = 5.74e31 (adimensional gigante)
# No es prometedor.

# ¿Σm_ν en unidades naturales del SSEE?
# m_φ = Σm_ν × H_0^alg (Paper 6) → m_φ = 5.60 eV
# Esto es solo un cambio de variable; no resuelve nada.

# Conclusión: si OP-9 ya admite que m_φ es ansatz, y m_φ = Σm_ν·H_0^alg,
# entonces Σm_ν tampoco es independiente. El problema real es el factor
# numérico 5.60 eV / H_0^alg → necesita una escala de masa fundamental.

m_phi_obs = 5.60  # eV (Paper 6 ansatz)
print(f"  Σm_ν observado    = {Sigma_mnu_obs} eV")
print(f"  m_φ = Σm_ν · H_alg = {m_phi_obs} eV (Paper 6)")
print(f"  Ratio m_φ/Σm_ν    = {m_phi_obs/Sigma_mnu_obs:.4f}")
print(f"  H_0 in 1/Mpc (alg) = {H0_alg:.4f}")
print()
print("  → m_φ y Σm_ν tienen el mismo grado de libertad fenomenológico.")
print("  → Derivar uno deriva el otro; no hay independencia.")

# ============================================================
# CONCLUSIÓN
# ============================================================
print("\n" + "=" * 72)
print("VEREDICTO PRELIMINAR")
print("=" * 72)
print("""
  H1 (fragilidad):    R depende linealmente de φ,π — no es infinitamente frágil,
                      pero tampoco es estructural — depende del "lugar" donde
                      4·KAL cae en la recta real respecto a enteros.

  H2 (saturación):    22 = 2·11 sugiere conexión M-theory/E11 pero es POST-HOC.
                      Sin argumento físico independiente, es numerología.

  H3 (scan algebr.):  [Ver tabla arriba — si vacía, ningún match limpio]

  Reformulación:      Σm_ν y m_φ comparten el mismo DoF. OP-14 y OP-9 NO son
                      problemas independientes — son el mismo problema con
                      dos caras.

  RECOMENDACIÓN:
    (a) Si scan H3 NO produce match: aceptar Σm_ν como input experimental SM
        y degradar el claim de Paper 4 a "consistencia con dato observado",
        no "derivación".
    (b) Si scan H3 SÍ produce match: investigar si la nueva forma tiene
        significado físico (saturación, conteo, etc.) o sigue siendo
        numerología en disfraz mejor.
    (c) Atacar OP-10 (V(φ) unificador): si V tiene mínimo a 5.60 eV
        derivado de su forma algebraica, entonces m_φ → corolario,
        OP-9 cae, y Σm_ν queda como predicción del modelo (no input).
""")
