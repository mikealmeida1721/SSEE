"""
P6 HONEST MATRIX — separando predicciones derivadas vs ajustes a datos
========================================================================
Tras descubrir que el "σ₈_eff = 0.737" de Paper 6 viene de FITEAR un
parámetro libre (alpha_WDM=1.6561 h/Mpc) para que CLASS reproduzca
el observable KiDS (0.737±0.020), esta matriz separa CLARAMENTE:

  - PREDICCIONES derivadas (sin parámetros fiteados a obs)
  - AJUSTES con parámetros libres fiteados a datos

Sin interpretación. Solo los números honestos en cada categoría.
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import OMEGA_DE, OMEGA_M_DYN, MIRA

# Constantes verificadas (de scripts existentes, no inventadas)
SIGMA8_LCDM       = 0.811        # Planck 2018 baseline
G_single_P5       = 0.8661       # P5 growth factor (Ω_m=0.160, ODE)
G_2s_P6           = 0.9788       # P6 two-sector growth factor (Ω_m=0.320, ODE)

# Observables (con errores)
KIDS_sig8         = 0.737;   KIDS_sig8_err = 0.020
KIDS_S8           = 0.766;   KIDS_S8_err   = 0.020
DES_S8            = 0.759;   DES_S8_err    = 0.023

# Corrección baryonic (HMcode-2020 OP-5, multiplicativa)
B_sigma8_OP5      = 0.9956       # = 0.4% supresión

# ─── PREDICCIONES DERIVADAS (sin fittear a datos) ───────────────────────
print("=" * 80)
print(" P6 HONEST MATRIX — sin interpretación, números separados")
print("=" * 80)

print(f"\n{'─'*80}")
print(" BLOQUE 1 — PREDICCIONES DERIVADAS (sin parámetros fiteados a obs)")
print(f"{'─'*80}")

# P5 single-sector — derivado de ODE de crecimiento
sig8_P5 = SIGMA8_LCDM * G_single_P5
S8_P5_single  = sig8_P5 * np.sqrt(OMEGA_M_DYN / 0.3)         # con Ω_m=0.160 puro
S8_P5_MIRA    = sig8_P5 * np.sqrt(0.320 / 0.3)                # con Ω_m,CMB=0.320

print(f"\n[A] Single-sector P5 (ODE growth, Ω_m,dyn=0.160)")
print(f"    σ₈ = 0.811 × G_1s = 0.811 × {G_single_P5:.4f} = {sig8_P5:.4f}")
print(f"    S₈ (con Ω_m=0.160): {S8_P5_single:.4f}")
print(f"    S₈ (con Ω_m,CMB=0.320, como reporta P5): {S8_P5_MIRA:.4f}")

# P6 two-sector — derivado de ODE de crecimiento (sin WDM, sin fittear)
sig8_P6_noWDM = SIGMA8_LCDM * G_2s_P6
S8_P6_noWDM   = sig8_P6_noWDM * np.sqrt(0.320 / 0.3)

print(f"\n[B] Two-sector P6 RAW (ODE growth, Ω_m,growth=0.320, NO WDM)")
print(f"    σ₈ = 0.811 × G_2s = 0.811 × {G_2s_P6:.4f} = {sig8_P6_noWDM:.4f}")
print(f"    S₈ = {S8_P6_noWDM:.4f}")

# ─── AJUSTES CON PARÁMETROS FITEADOS A DATOS ────────────────────────────
print(f"\n{'─'*80}")
print(" BLOQUE 2 — AJUSTES con parámetro libre fiteado a observación")
print(f"{'─'*80}")

print(f"""
[C] Two-sector P6 + WDM transfer
    Método: alpha_WDM CALIBRADO para que σ₈ = 0.737 (KiDS observado)
    Script:  class_ssee/calibrate_wdm_alpha.py
    Parámetro libre fiteado: alpha_WDM = 1.6561 h/Mpc
    Resultado: σ₈ = 0.7370 (PORQUE el target ERA 0.737)
               S₈ = 0.737 × √(0.32/0.3) = 0.7616

    ⚠ Esta no es predicción — es ajuste.
""")

sig8_OP5  = 0.737 * B_sigma8_OP5
S8_OP5    = sig8_OP5 * np.sqrt(0.320 / 0.3)

print(f"[D] Two-sector P6 + WDM (fitted) + HMcode baryonic feedback")
print(f"    σ₈ = 0.737 × {B_sigma8_OP5:.4f} = {sig8_OP5:.4f}")
print(f"    S₈ = {S8_OP5:.4f}")
print(f"    Este es el '0.758' reportado en OP-5 — built ENCIMA del fit.")

# ─── BLOQUE 3 — Comparación lado a lado con observaciones ──────────────
print(f"\n{'─'*80}")
print(" BLOQUE 3 — Comparación con observaciones (con estado epistemológico)")
print(f"{'─'*80}")

scenarios = [
    ("[A] Single-sector P5 (DERIVADO)", sig8_P5, S8_P5_MIRA, "derivado"),
    ("[B] Two-sector P6 RAW (DERIVADO)", sig8_P6_noWDM, S8_P6_noWDM, "derivado"),
    ("[C] P6 + WDM (FITEADO a KiDS)",  0.737, 0.7616, "ajustado"),
    ("[D] P6 + WDM + baryonic",         sig8_OP5, S8_OP5, "ajustado+corr"),
]

print(f"\n{'Escenario':<40} {'σ₈':>7} {'S₈':>7} {'Tipo':<14}")
print(" " + "─"*78)
for name, sig, S8, tipo in scenarios:
    print(f" {name:<40} {sig:>7.4f} {S8:>7.4f} {tipo:<14}")

print(f"\n{'Tensiones vs observaciones':<40}")
print(f"{'Escenario':<40} {'σ₈ vs KiDS':>12} {'S₈ vs KiDS':>12} {'S₈ vs DES':>12}")
print(" " + "─"*88)
for name, sig, S8, tipo in scenarios:
    t1 = (sig - KIDS_sig8) / KIDS_sig8_err
    t2 = (S8  - KIDS_S8)   / KIDS_S8_err
    t3 = (S8  - DES_S8)    / DES_S8_err
    print(f" {name:<40} {t1:>+12.2f} {t2:>+12.2f} {t3:>+12.2f}")

print(f"""
{'─'*80}
 INTERPRETACIÓN HONESTA
{'─'*80}

PREDICCIONES DERIVADAS:
  [A] Single-sector: σ₈=0.702, S₈=0.725  → ~2σ tensión KiDS/DES
  [B] Two-sector RAW: σ₈=0.794, S₈=0.820 → ~2.7σ tensión KiDS/DES

  Estas son las PREDICCIONES VERDADERAS de SSEE. Ni una ni otra
  matchea KiDS/DES limpiamente — ambas en ~2σ.

AJUSTES (NO predicciones):
  [C] P6 + WDM fitted: σ₈=0.737 (= target = KiDS observado)
      → 0σ tensión POR CONSTRUCCIÓN, no por predicción
      → alpha_WDM = 1.6561 es parámetro libre fiteado

  [D] P6 + WDM fitted + baryonic: S₈=0.758, 0.06σ DES
      → "0.06σ" es PUBLICITADO en Paper 6 pero viene del ajuste,
         no de la física

QUÉ DICE ESTO HONESTAMENTE:
  • Paper 6 reporta σ₈=0.737 como si fuera predicción del modelo,
    cuando en realidad es resultado de fitear alpha_WDM a KiDS.
  • El "0.06σ DES" comparable con observación es producto del
    parámetro libre fiteado, no de derivación.
  • Las predicciones REALES son 0.702 (single) o 0.794 (two-sector
    RAW), ambas con tensión moderada ~2σ.

CONCLUSIÓN: Paper 6 tiene un parámetro libre (alpha_WDM) además
            del m_φ numerológico ya identificado. Eso lo hace
            doblemente dependiente de calibración numerológica.
""")
print("=" * 80)
