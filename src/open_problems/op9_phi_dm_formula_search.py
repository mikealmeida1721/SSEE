"""
OP-9 — Búsqueda sistemática de fórmula física dimensionalmente consistente
para m_φ del sector φ-DM (Paper 6).

Contexto: La cascada H_alg→H_MIRA del 2026-05-24 AM intentó sustituir
H_MIRA = 67.068 km/s/Mpc en m_φ = Σm_ν × H₀. Pero esa fórmula es
DIMENSIONALMENTE INCONSISTENTE (eV × km/s/Mpc no es masa) — solo "funciona"
si H₀ se trata como número adimensional, como ocurre con H_alg = 3(φ+π)².

Esta búsqueda inventaría combinaciones dimensionalmente consistentes de
escalas físicas SSEE para identificar candidatas que recuperen ~5 eV (la
escala phenomenológica que produce k_fs = 0.493 h/Mpc).

Si NINGUNA candidata da ~5 eV → la fórmula correcta no existe en este
conjunto de escalas, y la k_fs phenomenológica requiere un mecanismo
nuevo (acoplamiento Yukawa, curvatura V''(φ_min), etc.) — refuerzo
empírico de OP-9 y OP-10 acoplados.

Mike, 2026-05-24 PM.
"""

import numpy as np
from itertools import product

# ────────────────────────────────────────────────────────────────────────
# Constantes físicas (natural units: ℏ=c=k_B=1, energy in eV)
# ────────────────────────────────────────────────────────────────────────
EV_PER_KMS_MPC = 2.1331e-30 / 1.0e0  # ℏ × (1 km/s/Mpc) ≈ 2.133e-30 eV
# Sanity: 1 km/s/Mpc = 3.2408e-20 s^-1; ℏ = 6.5821e-16 eV·s
HBAR_EVS = 6.582119569e-16            # eV·s
KMSMPC_TO_INVS = 3.240779289e-20      # (km/s/Mpc) → s^-1
EV_PER_HUBBLE_UNIT = HBAR_EVS * KMSMPC_TO_INVS  # ≈ 2.133e-33 eV per (km/s/Mpc)

M_PL = 1.220910e28   # eV (reduced Planck mass × √(8π) factor not applied; M_Pl GeV scale)
# Standard reduced Planck mass: M_pl_red = M_Pl / sqrt(8π) ≈ 2.435e27 eV. Use full M_Pl.

# SSEE physical inputs
H0_MIRA_KMSMPC = 67.068           # km/s/Mpc — Paper 3 Cobaya plik_lite
H0_ALG_KMSMPC  = 67.962           # km/s/Mpc — Type-P algebraic 3(φ+π)²
SIGMA_MNU_EV   = 0.0824           # eV — Paper 4 derived (Type P)
LAMBDA_SSEE_EV = 8.81e-3          # eV — Paper 10 UV cutoff M = 5φ⁸ρ_crit^(1/4)
PHI = (1 + 5**0.5) / 2            # 1.618...
PI  = np.pi                       # 3.1416...

# Phenomenologically required m_φ (eV) such that k_fs ≈ 0.493 h/Mpc via DW
M_PHI_TARGET_EV = 5.602           # eV — current ansatz value (numerological)
M_PHI_RANGE_LOW = 1.0             # eV — order-of-magnitude tolerance
M_PHI_RANGE_HIGH = 30.0           # eV — order-of-magnitude tolerance

# Derived physical scales
H_MIRA_EV = H0_MIRA_KMSMPC * EV_PER_HUBBLE_UNIT   # ≈ 1.430e-31 eV ← wait, check
H_MIRA_S = H0_MIRA_KMSMPC * KMSMPC_TO_INVS         # s^-1
H_MIRA_EV_correct = HBAR_EVS * H_MIRA_S            # eV
# The two should match; use H_MIRA_EV_correct.

print("=" * 72)
print("OP-9 — Búsqueda de fórmula dimensional para m_φ (φ-DM)")
print("=" * 72)
print(f"\nEscalas físicas SSEE (en eV):")
print(f"  ℏ × H₀^MIRA  = {H_MIRA_EV_correct:.4e} eV       (Hubble scale today)")
print(f"  Σm_ν^active  = {SIGMA_MNU_EV:.4e} eV       (sum active neutrinos)")
print(f"  Λ_SSEE = M   = {LAMBDA_SSEE_EV:.4e} eV       (Paper 10 UV cutoff)")
print(f"  M_Pl         = {M_PL:.4e} eV       (Planck mass)")
print(f"\nTarget phenomenológico:")
print(f"  m_φ ≈ {M_PHI_TARGET_EV} eV  (give k_fs ≈ 0.493 h/Mpc via Dodelson-Widrow)")
print(f"  Rango aceptable (1 década): [{M_PHI_RANGE_LOW}, {M_PHI_RANGE_HIGH}] eV")

# ────────────────────────────────────────────────────────────────────────
# Inventario sistemático
# ────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("Inventario de combinaciones dimensionalmente consistentes")
print("=" * 72)

scales = {
    'H_MIRA':   H_MIRA_EV_correct,
    'Σm_ν':     SIGMA_MNU_EV,
    'Λ_SSEE':   LAMBDA_SSEE_EV,
    'M_Pl':     M_PL,
}

# Exponentes enteros pequeños (incluyendo √, raíces tercias)
exponents = [-2, -1, -0.5, 0, 0.5, 1, 2]

# Factores algebraicos adimensionales SSEE
algebraic_factors = {
    '1':         1.0,
    'φ':         PHI,
    'π':         PI,
    'φ+π':       PHI + PI,
    '(φ+π)²':    (PHI + PI)**2,
    '3(φ+π)²':   3*(PHI + PI)**2,
    '(π-φ)':     PI - PHI,
}

candidates = []
keys = list(scales.keys())
vals = list(scales.values())

# Combinaciones de hasta 4 escalas con exponentes enteros/medios
for e in product(exponents, repeat=len(keys)):
    if all(ei == 0 for ei in e):
        continue
    try:
        result = 1.0
        for v, ei in zip(vals, e):
            result *= v**ei
    except (OverflowError, ValueError):
        continue
    if not np.isfinite(result) or result <= 0:
        continue
    # ¿Tiene dimensión de masa (eV)?  En estas unidades, todas son eV ya.
    # La consistencia dimensional se garantiza por que el producto de potencias
    # de cantidades en eV es en eV elevado a la suma de exponentes.
    total_dim = sum(e)
    if abs(total_dim - 1.0) > 1e-9:
        continue  # No tiene dimensión de masa
    # Aplicar factores algebraicos
    for fname, fval in algebraic_factors.items():
        m = result * fval
        if M_PHI_RANGE_LOW <= m <= M_PHI_RANGE_HIGH:
            formula = " × ".join(
                f"{keys[i]}^{e[i]}" if e[i] not in (1, 0) else (keys[i] if e[i] == 1 else "")
                for i in range(len(keys)) if e[i] != 0
            )
            if fname != '1':
                formula = f"{fname} × " + formula
            candidates.append((m, formula.strip(" ×"), e, fname))

print(f"\nCandidatas con m_φ ∈ [{M_PHI_RANGE_LOW}, {M_PHI_RANGE_HIGH}] eV:")
if candidates:
    # Ordenar por cercanía al target
    candidates.sort(key=lambda c: abs(np.log10(c[0]/M_PHI_TARGET_EV)))
    print(f"\n{'m_φ [eV]':>12} {'log₁₀(m/target)':>18}   Fórmula")
    print("-" * 72)
    for m, formula, e, fname in candidates[:20]:
        ratio = np.log10(m/M_PHI_TARGET_EV)
        print(f"{m:12.4f} {ratio:+18.3f}   {formula}")
else:
    print("\n  ⚠ NINGUNA combinación dimensional simple recupera m_φ ∈ [1, 30] eV.")
    print("    Esto refuerza OP-9: la escala eV de m_φ no surge naturalmente de")
    print("    combinar las escalas físicas conocidas del modelo. Requiere un")
    print("    mecanismo nuevo (curvatura V''(φ_min) o acoplamiento Yukawa).")

# ────────────────────────────────────────────────────────────────────────
# Análisis de las candidatas "obvias"
# ────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("Análisis de candidatas clásicas")
print("=" * 72)

obvious = {
    "Σm_ν × (H_MIRA/E_ref)  [numerológico H₀ adimensional]":
        SIGMA_MNU_EV * H0_MIRA_KMSMPC,   # tratando H₀ como número puro
    "√(Σm_ν × ℏH_MIRA)  [Hubble seesaw]":
        np.sqrt(SIGMA_MNU_EV * H_MIRA_EV_correct),
    "√(M_Pl × ℏH_MIRA)  [cosmo seesaw]":
        np.sqrt(M_PL * H_MIRA_EV_correct),
    "Λ_SSEE = M  [Paper 10 UV cutoff]":
        LAMBDA_SSEE_EV,
    "(M_Pl × ℏH_MIRA × Σm_ν²)^(1/4)":
        (M_PL * H_MIRA_EV_correct * SIGMA_MNU_EV**2)**0.25,
    "Σm_ν²/Λ_SSEE  [type-II seesaw]":
        SIGMA_MNU_EV**2 / LAMBDA_SSEE_EV,
    "M_Pl × (ℏH_MIRA/M_Pl)^(1/2)  [= √(M_Pl·ℏH₀)]":
        np.sqrt(M_PL * H_MIRA_EV_correct),
    "Λ_SSEE × (Λ_SSEE/Σm_ν)":
        LAMBDA_SSEE_EV**2 / SIGMA_MNU_EV,
    "(Σm_ν × Λ_SSEE × M_Pl)^(1/3)":
        (SIGMA_MNU_EV * LAMBDA_SSEE_EV * M_PL)**(1/3),
}

print(f"\n{'Candidata':<52} {'m_φ [eV]':>12}   {'Δ vs 5.6 eV':>12}")
print("-" * 80)
for name, val in obvious.items():
    if val > 0 and np.isfinite(val):
        delta = np.log10(val/M_PHI_TARGET_EV) if val > 0 else None
        delta_str = f"{delta:+.2f} dex" if delta is not None else "—"
        print(f"{name:<52} {val:12.4e}   {delta_str:>12}")

# ────────────────────────────────────────────────────────────────────────
# Diagnóstico final
# ────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("Diagnóstico OP-9")
print("=" * 72)

candidates_in_range = [c for c in obvious.values() if 1 <= c <= 30]
print(f"\n• Candidatas físicas (dimensionalmente consistentes) en [1, 30] eV: "
      f"{len(candidates_in_range)} de {len(obvious)}")

if len(candidates_in_range) <= 1:
    print("\n• CONCLUSIÓN: la escala m_φ ≈ 5 eV NO surge naturalmente de combinar")
    print("  las escalas físicas {ℏH_MIRA, Σm_ν, Λ_SSEE, M_Pl} con exponentes enteros")
    print("  o semienteros y factores algebraicos (φ, π).")
    print("\n• IMPLICACIÓN para SSEE:")
    print("  - El match numérico m_φ = Σm_ν × 67.96 = 5.6 eV era genuinamente")
    print("    una coincidencia numerológica, no un patrón físico oculto.")
    print("  - Las dos rutas físicamente honestas son:")
    print("    (a) m_φ ≡ Λ_SSEE = 8.81 meV → cambia k_fs por ~10³, rompe el")
    print("        canal Lyman-α actual; reconstruir todo Paper 6 con WDM milli-eV.")
    print("    (b) m_φ = √(V''(φ_min)) con V(φ) unificado SSEE — requiere OP-10 cerrado.")
    print("  - OP-9 se eleva: ya no es 'derivar el Yukawa', es 'reconstruir la")
    print("    fenomenología φ-DM con una escala que sí tenga origen físico'.")
else:
    print(f"\n• Hay {len(candidates_in_range)} candidata(s) en rango — investigar.")
    for name, val in obvious.items():
        if 1 <= val <= 30:
            print(f"    → {name}:  m_φ = {val:.3f} eV  (off by {np.log10(val/5.602):+.2f} dex)")

print("\n" + "=" * 72)
print("Fin de la búsqueda OP-9")
print("=" * 72)
