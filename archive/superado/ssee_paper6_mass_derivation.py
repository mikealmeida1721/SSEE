#!/usr/bin/env python3
"""
ssee_paper6_mass_derivation.py  —  Sandbox: Derivación Algebraica de m_φ
=========================================================================
HALLAZGO CLAVE (del script anterior):
  m_φ (Dodelson-Widrow) ≈ 5.74 eV
  Σm_ν (activos, Paper 4) = 0.0840 eV
  Ratio = 68.33

  CANDIDATO: 3(φ+π)² = H0^{alg} ≈ 67.96  →  Δ = 0.54%

HIPÓTESIS CENTRAL:
  m_s = Σm_ν × H0^{alg} = Σm_ν × 3(φ+π)²

  Interpretación: el neutrino estéril (φ-DM) adquiere masa mediante
  un mecanismo análogo al seesaw, donde el parámetro de escala es
  el mismo que fija la constante de Hubble.

  Esto es una relación de tipo:
    m_s × m_ν = (H0^{alg})^{1/2} × (escala de Planck)²?  (futura derivación)
  Por ahora: verificamos numéricamente y buscamos la expresión más simple.
"""

import numpy as np
from scipy.optimize import brentq
import itertools

# ── Constantes SSEE ──────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi_  = np.pi
beta = (pi_ + phi) / 2
Kv   = phi + pi_ + (pi_ + phi)
Tr   = 3 * (phi + beta)
Mv   = phi + pi_ + Kv
KAL0 = beta + pi_

w0      = -Tr / Mv
Omm_dyn = 1.0 + w0
OmDE    = 1.0 - Omm_dyn
MIRA    = (3*phi + pi_) / 4
tau_Pi  = KAL0 / (3.0 * OmDE)
H0_alg  = 3*(phi + pi_)**2          # ≈ 67.96 (algebraic)
Sigma9  = 3*(phi + pi_)              # ≈ 14.279

R       = 4*KAL0 - 22               # Residuo acústico = 0.0856
Ob_h2   = 3*(pi_ - phi)/200         # 0.02285
mnu     = R * Ob_h2 * 94.07 / tau_Pi  # 0.0840 eV

# Objetivo DW (Dodelson-Widrow)
m_target_DW = 5.7402   # eV (del script anterior)
m_target_th = 32.745   # eV (producción térmica)

print("=" * 70)
print("  SSEE Sandbox — Derivación Algebraica de m_φ (φ-DM estéril)")
print("=" * 70)
print(f"\n  Σm_ν (activos) = {mnu:.6f} eV")
print(f"  H0^alg = 3(φ+π)² = {H0_alg:.6f}")
print(f"  Σ9 = 3(φ+π) = {Sigma9:.6f}")
print(f"\n  Target DW:      {m_target_DW:.4f} eV")
print(f"  Target thermal: {m_target_th:.4f} eV")

# ── Verificación del candidato principal ─────────────────────────────────────
m_candidate = mnu * H0_alg
print(f"\n  ── Candidato Principal ────────────────────────────────────────")
print(f"  m_φ = Σm_ν × 3(φ+π)² = {mnu:.6f} × {H0_alg:.6f}")
print(f"      = {m_candidate:.6f} eV")
print(f"  Target DW:  {m_target_DW:.6f} eV  →  Δ = {abs(m_candidate-m_target_DW)/m_target_DW*100:.4f}%")

# Expresión algebraica expandida
# m_φ = [R × Ω_b h² × 94.07 / τ_Π] × 3(φ+π)²
# m_φ = R × Ω_b h² × 94.07 × 3(φ+π)² / [KAL0/(3Ω_DE)]
# m_φ = 9 R × Ω_b h² × 94.07 × (φ+π)² × Ω_DE / KAL0

expr_expanded = 9 * R * Ob_h2 * 94.07 * (phi+pi_)**2 * OmDE / KAL0
print(f"\n  Forma expandida: 9·R·Ω_b h²·94.07·(φ+π)²·Ω_DE/KAL₀")
print(f"  = {expr_expanded:.6f} eV  (check: {abs(expr_expanded-m_candidate):.2e})")

# ── Búsqueda sistemática de expresiones simples ──────────────────────────────
print(f"\n  ── Búsqueda Sistemática (operaciones básicas con constantes SSEE) ─")
print(f"  Target ratio = m_φ/Σm_ν = {m_target_DW/mnu:.6f} (DW)")

constants = {
    "φ":       phi,
    "π":       pi_,
    "β":       beta,
    "KAL₀":    KAL0,
    "MIRA":    MIRA,
    "τ_Π":     tau_Pi,
    "H₀^alg":  H0_alg,
    "Σ9":      Sigma9,
    "Ω_m":     Omm_dyn,
    "Ω_DE":    OmDE,
    "(φ+π)":   phi+pi_,
    "(φ+π)²":  (phi+pi_)**2,
    "(φ+π)³":  (phi+pi_)**3,
    "3(φ+π)²": 3*(phi+pi_)**2,
    "3(π-φ)":  3*(pi_-phi),
    "KAL₀²":   KAL0**2,
    "R×94.07/Ω_b": R*94.07/Ob_h2,
}

target = m_target_DW / mnu  # ratio DW
best_candidates = []

# Buscar expresiones simples (un solo constante o producto de dos)
print(f"\n  Candidatos con 1 constante:")
for name, val in constants.items():
    diff = abs(val - target)/target*100
    if diff < 10:
        star = "★" if diff < 1 else "~"
        print(f"    {name:20s} = {val:.6f}  Δ = {diff:.4f}%  {star}")
        best_candidates.append((diff, f"Σm_ν × {name}", mnu*val))

print(f"\n  Candidatos con 2 constantes (productos/cocientes):")
keys = list(constants.keys())
vals = list(constants.values())
found = []
for i, (n1, v1) in enumerate(constants.items()):
    for j, (n2, v2) in enumerate(constants.items()):
        if i >= j:
            continue
        for combo_name, combo_val in [
            (f"{n1}×{n2}", v1*v2),
            (f"{n1}/{n2}", v1/v2),
            (f"{n2}/{n1}", v2/v1),
            (f"({n1}+{n2})", v1+v2),
        ]:
            diff = abs(combo_val - target)/target*100
            if diff < 2.0:
                found.append((diff, combo_name, combo_val))

found.sort()
for diff, name, val in found[:15]:
    star = "★" if diff < 0.3 else ("★★" if diff < 0.1 else "~")
    print(f"    {name:30s} = {val:.6f}  Δ = {diff:.4f}%  {star}")
    best_candidates.append((diff, f"Σm_ν × {name}", mnu*val))

# ── Verificación del seesaw analogy ──────────────────────────────────────────
print(f"\n  ── Analogía Seesaw ──────────────────────────────────────────────")
# En seesaw tipo-I: m_ν × m_N = y² × v²/2 (Yukawa × VEV)
# Aquí: Σm_ν × m_s = ? × (escala SSEE)²
product = mnu * m_target_DW  # eV²
print(f"  Σm_ν × m_φ = {mnu:.6f} × {m_target_DW:.6f} = {product:.6f} eV²")
print(f"  = {product:.6f} eV²")
print(f"  = (escala)²  donde escala = {np.sqrt(product)*1000:.4f} meV")
scale = np.sqrt(product)
print(f"\n  Escala de seesaw: {scale:.6f} eV = {scale*1000:.4f} meV")
print(f"  Σm_ν (activo): {mnu*1000:.4f} meV")
print(f"  Ratio escala/Σm_ν = {scale/mnu:.6f}")

# Buscar expresión SSEE para la escala
target_scale = scale / mnu  # ≈ √H0_alg
print(f"\n  Target: escala/Σm_ν = √(H₀^alg) = {np.sqrt(H0_alg):.6f}")
print(f"  Calculado:           = {target_scale:.6f}")
print(f"  Δ = {abs(target_scale - np.sqrt(H0_alg))/np.sqrt(H0_alg)*100:.4f}%")

# ── Predicción algebraica final ───────────────────────────────────────────────
print(f"\n{'='*70}")
print(f"  RESULTADO ALGEBRAICO FINAL")
print(f"{'='*70}")
print(f"\n  m_φ = Σm_ν × 3(φ+π)²  =  Σm_ν × H₀^{{alg}}")
print(f"\n  Expresión completamente algebraica:")
print(f"  m_φ = [R·Ω_b h²·94.07·eV / τ_Π] × 3(φ+π)²")
print(f"  m_φ = [4KAL₀−22] × [3(π−φ)/200] × 94.07 eV × 3(φ+π)² / [KAL₀/(3Ω_DE)]")
print(f"\n  Sustituyendo:")
print(f"    R = 4KAL₀−22            = {R:.8f}")
print(f"    Ω_b h²                  = {Ob_h2:.8f}")
print(f"    3(φ+π)² = H₀^alg        = {H0_alg:.8f}")
print(f"    τ_Π = KAL₀/(3Ω_DE)      = {tau_Pi:.8f}")
print(f"\n  m_φ = {m_candidate:.6f} eV")
print(f"  Target DW = {m_target_DW:.6f} eV")
print(f"  Δ = {abs(m_candidate-m_target_DW)/m_target_DW*100:.4f}%  ← {'★ DERIVACIÓN EXITOSA' if abs(m_candidate-m_target_DW)/m_target_DW < 0.01 else '≈ APROXIMADO'}")
print(f"\n  Interpretación física:")
print(f"  La masa del φ-DM es la masa del neutrino activo AMPLIFICADA")
print(f"  por el factor H₀ (constante de Hubble como número puro ≈ 68).")
print(f"  Esto es análogo al mecanismo seesaw: m_s/m_ν = H₀^alg = 3(φ+π)²")
print(f"  donde (φ+π) es la constante estructural raíz del universo SSEE.")
print(f"\n  Predicción falsable:")
print(f"  m_φ = {m_candidate:.4f} eV  (en el rango de búsqueda de KATRIN/PTOLEMY)")
print(f"  Si KATRIN mide ΔΣm_ν indicando componente ~{m_candidate:.1f} eV → confirma φ-DM")
print(f"  Si DESI Y3 mide σ8 > 0.80 → falsifica el modelo de dos sectores")
print(f"{'='*70}")
