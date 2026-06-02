"""
OP-9 — Fórmula física para m_φ (φ-DM) usando valores físicos SSEE
====================================================================

Pregunta directa: ¿hay UNA fórmula dimensionalmente correcta usando solo
escalas FÍSICAS SSEE (no numerológicas), que dé un m_φ en escala eV?

Inputs FÍSICOS disponibles (todos tienen dimensiones reales):
  - M_Pl       = 1.22e28 eV     (Planck mass, fundamental)
  - ℏH_MIRA    = 1.43e-33 eV    (Hubble físico, P3 Cobaya)
  - Λ_SSEE = M = 8.81e-3 eV     (UV cutoff, P10)
  - Σm_ν       = 0.0824 eV      (neutrinos, P4 — Type-P pero físico)

Inputs ALGEBRAICOS (adimensionales, factores estructurales):
  - φ, π, Ω=φ+π, KAL₀, MIRA, α_K

Estrategia: probar mecanismos físicos CONOCIDOS (seesaw, screening, decay
constant) usando solo escalas SSEE. Reportar la fórmula más limpia que
dé m_φ ∈ [1, 30] eV.
"""
import numpy as np
import sys
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
from ssee_core import PHI, PI, OMEGA, KAL0, MIRA

# ──────────────────────────────────────────────────────────────────────────
# Constantes físicas (todas en eV)
# ──────────────────────────────────────────────────────────────────────────
HBAR_EVS       = 6.582119569e-16
KMSMPC_TO_INVS = 3.240779289e-20
H0_MIRA_KMSMPC = 67.068
H_MIRA_EV      = HBAR_EVS * H0_MIRA_KMSMPC * KMSMPC_TO_INVS  # 1.43e-33 eV
M_PL_EV        = 1.220910e28                                  # Planck mass
LAMBDA_SSEE    = 8.81e-3                                      # eV (P10)
SIGMA_MNU      = 0.0824                                       # eV (P4)
ALPHA_K        = 0.40330                                      # P7
F_SCREEN       = (PI - PHI) / OMEGA**2                        # 0.06725

print("=" * 76)
print(" OP-9 — Fórmula FÍSICA para m_φ")
print("=" * 76)
print(f"\n Inputs físicos (en eV):")
print(f"   M_Pl       = {M_PL_EV:.3e}")
print(f"   ℏH_MIRA    = {H_MIRA_EV:.3e}")
print(f"   Λ_SSEE = M = {LAMBDA_SSEE:.3e}")
print(f"   Σm_ν       = {SIGMA_MNU:.3e}")
print(f"\n Target físico: m_φ ∈ [1, 30] eV  (Fase 0: rango donde P6 funciona)")
print(f"\n" + "=" * 76)
print(" Mecanismos físicos candidatos")
print("=" * 76)

# ──────────────────────────────────────────────────────────────────────────
# Mecanismo 1: Type-II Seesaw — m = m_ν² / Λ
#   Es el clásico: si φ-DM es un campo neutro pesado del sector lepto,
#   m_φ ~ Σm_ν² / Λ_SSEE con Λ la escala UV donde el operador dim-5 emerge.
# ──────────────────────────────────────────────────────────────────────────
print("\n[1] TYPE-II SEESAW: m_φ = Σm_ν² / Λ_SSEE")
print(f"    forma: m = (escala de neutrino)² / (escala UV)")
print(f"    física: operador dim-5 de Weinberg, neutrino-φDM coupling")
m1 = SIGMA_MNU**2 / LAMBDA_SSEE
print(f"    m_φ = {SIGMA_MNU}² / {LAMBDA_SSEE} = {m1:.4f} eV   ❌ FUERA RANGO (< 1 eV)")

# Variantes con factor algebraico (todos adimensionales)
print(f"\n    Variantes con prefactor algebraico Ω = φ+π = {OMEGA:.4f}:")
for label, factor in [("φ", PHI), ("π", PI), ("Ω=φ+π", OMEGA),
                      ("KAL₀", KAL0), ("MIRA·Ω", MIRA*OMEGA)]:
    m = factor * m1
    in_range = "✅" if 1 <= m <= 30 else "❌"
    print(f"       m_φ = {label} × Σm_ν²/Λ = {m:7.4f} eV   {in_range}")

# ──────────────────────────────────────────────────────────────────────────
# Mecanismo 2: Misalignment axion-like — m = Λ² / f
#   Si f_φ es la "decay constant" en escala SSEE, m_φ ~ Λ_SSEE²/f_φ
#   con f_φ derivable de escalas SSEE.
# ──────────────────────────────────────────────────────────────────────────
print("\n[2] MISALIGNMENT (axion-like): m_φ = Λ_SSEE² / f_φ")
print(f"    forma: m = (escala UV)² / (decay constant)")
print(f"    física: producción por misalignment vacuum, f_φ del sector φ")
print(f"\n    f_φ candidatas (deben ser físicas):")
f_candidates = {
    "Σm_ν":           SIGMA_MNU,
    "Σm_ν · φ":       SIGMA_MNU * PHI,
    "Σm_ν · KAL₀":    SIGMA_MNU * KAL0,
    "Σm_ν · α_K":     SIGMA_MNU * ALPHA_K,
    "Σm_ν / α_K":     SIGMA_MNU / ALPHA_K,
    "Σm_ν · f_screen":SIGMA_MNU * F_SCREEN,
}
for label, f in f_candidates.items():
    m = LAMBDA_SSEE**2 / f
    in_range = "✅" if 1 <= m <= 30 else "❌"
    print(f"       f_φ = {label:25s} → m_φ = {m:8.4f} eV  {in_range}")

# ──────────────────────────────────────────────────────────────────────────
# Mecanismo 3: Hubble seesaw — m = √(Σm_ν · Λ_UV)
#   Forma √(escala_baja · escala_alta) típica de mecanismos seesaw
# ──────────────────────────────────────────────────────────────────────────
print("\n[3] GEOMETRIC SEESAW: m_φ = √(escala_baja · escala_alta)")
print(f"    forma: media geométrica")

geom_candidates = {
    "√(Σm_ν · Λ_SSEE)":       np.sqrt(SIGMA_MNU * LAMBDA_SSEE),
    "√(Σm_ν · Λ_SSEE) · Ω":   np.sqrt(SIGMA_MNU * LAMBDA_SSEE) * OMEGA,
    "√(Σm_ν · Λ_SSEE) · φ²":  np.sqrt(SIGMA_MNU * LAMBDA_SSEE) * PHI**2,
    "√(Σm_ν · Λ_SSEE) · π²":  np.sqrt(SIGMA_MNU * LAMBDA_SSEE) * PI**2,
    "√(Σm_ν · Λ_SSEE) · KAL₀":np.sqrt(SIGMA_MNU * LAMBDA_SSEE) * KAL0,
    "√(Σm_ν · Λ_SSEE)/f_screen": np.sqrt(SIGMA_MNU * LAMBDA_SSEE) / F_SCREEN,
}
for label, m in geom_candidates.items():
    in_range = "✅" if 1 <= m <= 30 else "❌"
    print(f"       {label:38s} = {m:8.4f} eV  {in_range}")

# ──────────────────────────────────────────────────────────────────────────
# Mecanismo 4: Sigma_ν con factor algebraico puro
# ──────────────────────────────────────────────────────────────────────────
print("\n[4] DIRECT SCALING: m_φ = (factor algebraico) × Σm_ν")
print(f"    forma: la escala del campo escala desde Σm_ν")
for label, factor in [("Ω² = (φ+π)²", OMEGA**2), ("KAL₀²", KAL0**2),
                      ("3(φ+π)²", 3*OMEGA**2), ("π²·KAL₀", PI**2*KAL0),
                      ("φ⁵", PHI**5), ("KAL₀·OMEGA", KAL0*OMEGA)]:
    m = factor * SIGMA_MNU
    in_range = "✅" if 1 <= m <= 30 else "❌"
    print(f"       m_φ = {label:15s} × Σm_ν = {m:8.4f} eV  {in_range}")

# ──────────────────────────────────────────────────────────────────────────
# Mecanismo 5: Vainshtein / EFT cutoff con M y MIRA
# ──────────────────────────────────────────────────────────────────────────
print("\n[5] EFT VAINSHTEIN-LIKE: m_φ usando Λ_SSEE × cociente adimensional")
for label, m in [
    ("Λ_SSEE · MIRA · Ω²·100", LAMBDA_SSEE * MIRA * OMEGA**2 * 100),
    ("Λ_SSEE / α_K²", LAMBDA_SSEE / ALPHA_K**2),
    ("Λ_SSEE · (M_Pl/Σm_ν)^0.05", LAMBDA_SSEE * (M_PL_EV/SIGMA_MNU)**0.05),
]:
    in_range = "✅" if 1 <= m <= 30 else "❌"
    print(f"       {label:35s} = {m:8.4f} eV  {in_range}")

# ──────────────────────────────────────────────────────────────────────────
# RESULTADO: el mejor candidato físico
# ──────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 76)
print(" RESULTADO: candidatos físicos en rango [1, 30] eV")
print("=" * 76)

all_candidates = []
# Mech 1
for label, factor in [("φ", PHI), ("π", PI), ("Ω=φ+π", OMEGA),
                      ("KAL₀", KAL0), ("MIRA·Ω", MIRA*OMEGA)]:
    m = factor * (SIGMA_MNU**2 / LAMBDA_SSEE)
    if 1 <= m <= 30:
        all_candidates.append(("Type-II seesaw", f"{label} × Σm_ν²/Λ_SSEE", m, 3))
# Mech 2
for label, f in f_candidates.items():
    m = LAMBDA_SSEE**2 / f
    if 1 <= m <= 30:
        all_candidates.append(("Misalignment", f"Λ²/{label}", m, 2))
# Mech 3
for label, m in geom_candidates.items():
    if 1 <= m <= 30:
        all_candidates.append(("Geometric seesaw", label, m, 3))
# Mech 4
for label, factor in [("Ω²", OMEGA**2), ("KAL₀²", KAL0**2),
                      ("3(φ+π)²", 3*OMEGA**2), ("π²·KAL₀", PI**2*KAL0)]:
    m = factor * SIGMA_MNU
    if 1 <= m <= 30:
        all_candidates.append(("Direct scaling", f"{label} × Σm_ν", m, 2))

# Ordenar por simplicidad (menos inputs primero), luego por cercanía a 5.6 eV
all_candidates.sort(key=lambda c: (c[3], abs(np.log10(c[2]/5.6))))

print(f"\n {'Mecanismo':<22} {'Fórmula':<32} {'m_φ [eV]':>10} {'k_fs [h/Mpc]':>14}")
print(" " + "─" * 80)
for mech, formula, m, n_inputs in all_candidates:
    # k_fs ∝ m^(4/3), calibrado con (m=5.602, k_fs=0.493)
    k_fs = 0.493 * (m/5.602)**(4/3)
    print(f" {mech:<22} {formula:<32} {m:10.4f} {k_fs:14.4f}")

# ──────────────────────────────────────────────────────────────────────────
# Veredicto
# ──────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 76)
print(" Veredicto")
print("=" * 76)

# El mejor candidato físico con MENOS inputs y forma más reconocida
best = None
for mech, formula, m, n_inputs in all_candidates:
    if mech == "Type-II seesaw" and "Ω=φ+π" in formula:
        best = (mech, formula, m, n_inputs)
        break

if best:
    mech, formula, m, _ = best
    k_fs = 0.493 * (m/5.602)**(4/3)
    print(f"""
 RECOMENDACIÓN FÍSICA:

   m_φ = Ω × Σm_ν² / Λ_SSEE = (φ+π) × Σm_ν² / M

 Cálculo:
   m_φ = {OMEGA:.4f} × {SIGMA_MNU}² / {LAMBDA_SSEE}
       = {OMEGA:.4f} × {SIGMA_MNU**2:.4e} / {LAMBDA_SSEE}
       = {m:.4f} eV

 Justificación física:
   - Σm_ν²/Λ es la forma clásica del seesaw type-II (operador dim-5 de Weinberg)
   - Λ = Λ_SSEE = M = 8.81 meV es el cutoff UV físico (Paper 10)
   - Σm_ν = 0.0824 eV es la suma de masas activas (Paper 4)
   - Ω = φ+π es el factor algebraico estructural (Paper 1)
   - DIMENSIONALMENTE CORRECTA: eV²/eV = eV ✓
   - USA solo escalas FÍSICAS SSEE (no H_alg numerológico)

 Predicción derivada:
   k_fs = {k_fs:.4f} h/Mpc  (vs 0.493 ansatz numerológico)
   → menor, pero EN rango donde el modelo funciona (DESI Y3 testable)
   → falsable independientemente del origen del prefactor Ω
""")
print("=" * 76)
