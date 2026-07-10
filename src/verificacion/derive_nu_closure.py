#!/usr/bin/env python3
"""Demostración desde primeros principios de la constante de cierre de neutrinos C
en la relación  Ω_ν h² = Σm_ν / C.

Objetivo (auditoría 2026-07-10): probar —no afirmar— que 94.07 eV y 93.14 eV son
EL MISMO objeto físico (C), y por qué difieren ~1%. Un referí hostil puede correr
esto y comprobar cada paso.

Cadena:
  1. n_γ hoy desde T_γ (Bose-Einstein).
  2. n_ν por especie (ν+ν̄) = (3/11) n_γ  [factor fermiónico 3/4 × (T_ν/T_γ)³=(4/11)].
  3. ρ_crit,0 (para h=1) en eV/cm³.
  4. Ω_ν h² = n_ν Σm_ν / ρ_crit  →  C ≡ Σm_ν/(Ω_ν h²) = ρ_crit / n_ν.
  Esto da C_instantáneo ≈ 94.1 eV (desacople instantáneo, N_eff=3 exacto).
  El valor preciso 93.14 eV (Mangano+2005; usado por Planck) baja C ~1% por el
  calentamiento no-instantáneo de los ν en la aniquilación e⁺e⁻ (mismo efecto que
  N_eff=3.046). Ese último paso NO se re-deriva aquí (requiere integrar la ec. de
  Boltzmann del desacople); se cita. Lo demostrable es que AMBOS son la misma C.
"""
import math

# ── Constantes físicas (CODATA / PDG) ────────────────────────────────────────
k_B   = 8.617333262e-5      # eV/K
hbar  = 6.582119569e-16     # eV·s
c     = 2.99792458e10       # cm/s
zeta3 = 1.2020569           # ζ(3)
G     = 6.67430e-8          # cm³ g⁻¹ s⁻² (cgs)
T_gamma = 2.7255            # K (COBE/FIRAS)
eV_per_g = 5.609588604e32   # eV / g  (=c²/e en unidades apropiadas)

# ── 1. Densidad numérica de fotones hoy ──────────────────────────────────────
# n_γ = (2 ζ(3)/π²) (k_B T / ħ c)³
T_eV = k_B * T_gamma                        # eV
n_gamma = (2 * zeta3 / math.pi**2) * (T_eV / (hbar * c))**3   # cm⁻³
print(f"T_γ            = {T_gamma} K = {T_eV:.4e} eV")
print(f"n_γ            = {n_gamma:.2f} cm⁻³   (esperado ≈ 411)")

# ── 2. Densidad numérica de neutrinos por especie (ν + ν̄) ────────────────────
# n_ν = (3/4)·(g_ν/g_γ)·(T_ν/T_γ)³·n_γ ; g_ν=2, g_γ=2, (T_ν/T_γ)³=4/11  → 3/11
n_nu = (3.0/11.0) * n_gamma
print(f"n_ν/especie    = {n_nu:.2f} cm⁻³   (esperado ≈ 112; ν+ν̄)")

# ── 3. Densidad crítica hoy para h=1 ─────────────────────────────────────────
H0_h1 = 100 * 1e5 / (3.0857e24)             # h=1: 100 km/s/Mpc → s⁻¹  (1 Mpc=3.0857e24 cm)
rho_crit_g = 3 * H0_h1**2 / (8 * math.pi * G)   # g/cm³
rho_crit_eV = rho_crit_g * eV_per_g             # eV/cm³
print(f"ρ_crit (h=1)   = {rho_crit_g:.4e} g/cm³ = {rho_crit_eV:.4e} eV/cm³")

# ── 4. Constante de cierre C = ρ_crit / n_ν  (Ω_ν h² = Σm_ν/C) ────────────────
C_derived = rho_crit_eV / n_nu
print(f"\nC (derivada, desacople instantáneo) = ρ_crit/n_ν = {C_derived:.2f} eV")
print(f"  → coincide con el valor de libro  94.07 eV  (dif {abs(C_derived-94.07)/94.07*100:.2f}%)")

# ── 5. El valor preciso citado y el ~1% ──────────────────────────────────────
C_instant, C_precise = 94.07, 93.14
print(f"\nDos valores de la MISMA C:")
print(f"  C_instantáneo (N_eff=3.000)  = {C_instant} eV   [derivado arriba]")
print(f"  C_preciso     (N_eff=3.046)  = {C_precise} eV   [Mangano+2005, usado por Planck]")
print(f"  ratio = {C_instant/C_precise:.4f}  (~1%: calentamiento e⁺e⁻ no-instantáneo)")

# ── 6. Impacto en Σm_ν y m_φ (misma fórmula SSEE, sólo cambia C) ─────────────
phi = (1+5**0.5)/2; pi = math.pi
Om = phi+pi; KAL = (phi+pi)/2 + pi; TRIAL = 3*(phi+(phi+pi)/2)
R2 = Om/(KAL*TRIAL); wb = (pi-phi)/(3*(phi+pi)**2); tau = KAL/(3*(3*(phi+(phi+pi)/2)/(3*(phi+pi))))
mult = 594.28
print(f"\nImpacto (Σm_ν = R₂·ω_b·C/τ_Π H₀,  m_φ = Σm_ν·{mult}):")
for C in (C_instant, C_precise):
    S = R2 * wb * C / 2.191
    print(f"  C={C}: Σm_ν={S:.5f} eV | m_φ={S*mult:.2f} eV | ω_ν=Σ/93.14={S/93.14:.6f}")
print("\nσ(Σm_ν) esperado DESI DR5/Euclid ≈ 0.015 eV (22%): el ~1% de C es RUIDO frente")
print("a la barra de error. Pero como es la MISMA C, lo sólido es usar la PRECISA (93.14)")
print("en todo, y así no hay dos valores que un referí pueda picar.")
