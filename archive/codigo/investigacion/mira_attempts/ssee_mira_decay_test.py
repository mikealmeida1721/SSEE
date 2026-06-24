#!/usr/bin/env python3
"""
ssee_mira_decay_test.py — ¿Existe τ algebraico de φ,π que ubique el decaimiento
de la componente extra MIRA (0.160) entre z~1100 y z~2?

Modelo: ρ_extra(a) = ρ₀ · a⁻³ · exp(-t(a)/τ)
Criterio de paso:
  f(z=1100) > 0.95  →  componente presente en recombinación (Planck ve 0.320)
  f(z=2)    < 0.10  →  componente ausente en surveys tardíos (DESI ve 0.160)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssee_core import (PHI as phi, PI as pi, OMEGA, BETA, KAL0, P_SC,
                       K_V, M_V, T_R, W0, WA, OMEGA_M_DYN)

import numpy as np
from scipy.integrate import quad

Om  = OMEGA_M_DYN        # 0.1601
Ode = 1.0 - Om           # 0.8399


def f_de(a):
    return a ** (-3 * (1 + W0 + WA)) * np.exp(-3 * WA * (1 - a))


def E2(a):
    return Om * a**-3 + Ode * f_de(a)


def t_hubble(z):
    """Tiempo cósmico en unidades de H₀⁻¹."""
    a_end = 1.0 / (1.0 + z)
    val, _ = quad(lambda a: 1.0 / (a * np.sqrt(E2(a))), 1e-7, a_end, limit=300)
    return val


print("=" * 70)
print("TEST CAMINO A — ¿Decaimiento algebraico para MIRA?")
print("=" * 70)
print("Calculando tiempos cósmicos SSEE...")

t_rec   = t_hubble(1100)
t_z2    = t_hubble(2.0)
t_z05   = t_hubble(0.5)
t_now   = t_hubble(0.0)

GYR = 14.4   # H₀⁻¹ ≈ 14.4 Gyr para H₀=67.96

print(f"  t(z=1100) = {t_rec:.5f} H₀⁻¹  = {t_rec*GYR*1e3:.1f} Myr")
print(f"  t(z=2)    = {t_z2:.5f} H₀⁻¹  = {t_z2*GYR:.2f} Gyr")
print(f"  t(z=0.5)  = {t_z05:.5f} H₀⁻¹  = {t_z05*GYR:.2f} Gyr")
print(f"  t(z=0)    = {t_now:.5f} H₀⁻¹  = {t_now*GYR:.2f} Gyr")

# Candidatos algebraicos (en unidades de H₀⁻¹)
candidates = [
    ("1/π",         1/pi),
    ("1/φ",         1/phi),
    ("1/Ω",         1/OMEGA),
    ("φ",           phi),
    ("π",           pi),
    ("φ²",          phi**2),
    ("β=(φ+π)/2",   BETA),
    ("φ³",          phi**3),
    ("Ω=φ+π",       OMEGA),
    ("π²/φ²",       pi**2/phi**2),
    ("KAL₀",        KAL0),
    ("P_sc",        P_SC),
    ("φ⁴",          phi**4),
    ("π²",          pi**2),
    ("φ⁵",          phi**5),
    ("K_V",         K_V),
    ("φ⁶",          phi**6),
    ("T_R",         T_R),
    ("M_V=Σ₉",      M_V),
    ("π³",          pi**3),
    ("φ⁷",          phi**7),
    ("Ω²",          OMEGA**2),
    ("φ⁸",          phi**8),
]

PASS_REC  = 0.95
PASS_DESI = 0.10

print(f"\n{'Candidato τ':<18} {'τ[H₀⁻¹]':>9}  {'f(z=1100)':>10}  "
      f"{'f(z=2)':>8}  {'f(z=0.5)':>9}  {'?':>6}")
print("-" * 70)

passing = []
for name, tau in candidates:
    f_r  = np.exp(-t_rec  / tau)
    f_d2 = np.exp(-t_z2   / tau)
    f_d5 = np.exp(-t_z05  / tau)
    ok   = (f_r > PASS_REC) and (f_d2 < PASS_DESI)
    flag = "PASA ✓" if ok else ""
    print(f"  {name:<16} {tau:>9.4f}  {f_r:>10.4f}  {f_d2:>8.4f}  {f_d5:>9.4f}  {flag}")
    if ok:
        passing.append((name, tau, f_r, f_d2, f_d5))

print("-" * 70)

# Búsqueda del τ óptimo numérico
tau_grid = np.logspace(-4, 2, 50000)
best_score, tau_opt = -1, tau_grid[0]
for tau in tau_grid:
    f_r = np.exp(-t_rec / tau)
    f_d = np.exp(-t_z2  / tau)
    if f_r > 0.5 and f_d < 0.5:
        s = f_r * (1 - f_d)
        if s > best_score:
            best_score, tau_opt = s, tau

f_opt_rec = np.exp(-t_rec  / tau_opt)
f_opt_d2  = np.exp(-t_z2   / tau_opt)
f_opt_d5  = np.exp(-t_z05  / tau_opt)

print(f"\nτ óptimo numérico (máxima separación f_rec vs f_z2):")
print(f"  τ_opt = {tau_opt:.5f} H₀⁻¹  = {tau_opt*GYR:.3f} Gyr")
print(f"  f(z=1100)={f_opt_rec:.4f}   f(z=2)={f_opt_d2:.4f}   f(z=0.5)={f_opt_d5:.4f}")

print(f"\n  Candidatos algebraicos más cercanos a τ_opt:")
for name, tau in candidates:
    ratio = tau / tau_opt
    if 0.3 < ratio < 3.0:
        print(f"    {name:<16} = {tau:.5f}  (τ_alg/τ_opt = {ratio:.3f})")

print("\n" + "=" * 70)
if passing:
    print(f"VEREDICTO: {len(passing)} candidato(s) algebraico(s) PASAN el corte.")
    for name, tau, fr, fd2, fd5 in passing:
        print(f"  τ = {name} = {tau:.4f} H₀⁻¹")
        print(f"    f(z=1100)={fr:.4f}  f(z=2)={fd2:.4f}  f(z=0.5)={fd5:.4f}")
    print("  → Camino A ABIERTO: existe τ algebraico con la época correcta.")
else:
    print("VEREDICTO: Ningún candidato algebraico pasa el corte exacto.")
    print(f"  τ_opt numérico = {tau_opt:.4f} H₀⁻¹ ({tau_opt*GYR:.2f} Gyr)")
    print("  Verificar si τ_opt corresponde a alguna combinación no listada.")
    print("  → Camino A NECESITA revisión o está cerrado.")
print("=" * 70)
