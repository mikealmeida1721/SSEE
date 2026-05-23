#!/usr/bin/env python3
"""
ssee_mira_projection_test.py — Test Camino B: ¿proyección 3D→2D?

Si MIRA es un efecto de proyección, la materia real es Ω_m=0.320 y DESI
la "ve" como 0.160 por geometría. Test: calcular los observables BAO de DESI
con Ω_m=0.320 y Ω_m=0.160 y ver cuál se ajusta mejor a los datos reales.

Si 0.320 falla en DESI (especialmente en D_H=c/H(z), que es radial, no proyección)
→ Camino B cerrado: la materia real ES 0.160, no 0.320 proyectado.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssee_core import W0, WA, H0_ALG

import numpy as np
from scipy.integrate import quad

# ── DESI DR2 datos (Abdul-Karim et al. 2025, arXiv:2503.14738) ──────────────
# Exactamente los mismos que ssee_paper2_mcmc.py (fuente canónica)
DESI_Z    = [0.295, 0.510, 0.510, 0.706, 0.706, 0.930, 0.930,
             1.317, 1.317, 1.491, 1.491, 2.330, 2.330]
DESI_TYPE = ["DV_rd","DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd",
             "DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd"]
DESI_OBS  = np.array([7.93, 13.62, 20.08, 16.85, 19.50, 21.71, 17.88,
                      27.79, 13.82, 30.21, 13.23, 39.71,  8.52])
DESI_SIGMA = np.array([0.15, 0.25, 0.60, 0.32, 0.55, 0.28, 0.35,
                       0.69, 0.42, 0.79, 0.55, 0.94, 0.17])

c_kms    = 299792.458
OB_H2    = 0.02237   # Ωb h² fijo (Planck 2018)

def sound_horizon_rd(Om, H0=H0_ALG, ob_h2=OB_H2):
    """r_d auto-consistente con Ω_m (Eisenstein & Hu 1998, igual que Paper 2)."""
    om_h2 = Om * (H0/100)**2
    return 147.27 * (om_h2/0.1432)**(-0.255) * (ob_h2/0.02237)**(-0.134)

def f_de(a, w0=W0, wa=WA):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def H_over_H0(z, Om, w0=W0, wa=WA):
    a = 1/(1+z)
    Ode = 1 - Om
    return np.sqrt(Om * a**-3 + Ode * f_de(a, w0, wa))

def D_M(z, Om, w0=W0, wa=WA):
    """Distancia comóvil en Mpc."""
    DH0 = c_kms / H0_ALG
    val, _ = quad(lambda zp: 1.0 / H_over_H0(zp, Om, w0, wa), 0, z, limit=200)
    return DH0 * val

def D_H(z, Om, w0=W0, wa=WA):
    """Distancia de Hubble c/H(z) en Mpc."""
    return c_kms / (H0_ALG * H_over_H0(z, Om, w0, wa))

def D_V(z, Om, w0=W0, wa=WA):
    """Distancia volumétrica [z * DM² * DH]^(1/3)."""
    dm = D_M(z, Om, w0, wa)
    dh = D_H(z, Om, w0, wa)
    return (z * dm**2 * dh)**(1.0/3.0)

# ── Calcular χ² para cada caso ───────────────────────────────────────────────
print("=" * 68)
print("TEST CAMINO B — ¿Ω_m=0.320 (3D real) vs 0.160 (dinámico)?")
print("TEST usa r_d AUTO-CONSISTENTE con cada Ω_m (igual que Paper 2)")
print("=" * 68)
print(f"  Parámetros fijos: w₀={W0:.4f}, wₐ={WA:.4f}, H₀={H0_ALG:.4f}, Ωb h²={OB_H2}")

cases = [("Ω_m=0.160 (dinámico, Paper 2)", 0.1601),
         ("Ω_m=0.320 (CMB/MIRA, 3D real)", 0.3199)]

results = {}
for label, Om in cases:
    chi2 = 0.0
    n    = 0
    rows = []
    r_d  = sound_horizon_rd(Om)
    for i, (z, typ, obs, sigma) in enumerate(zip(DESI_Z, DESI_TYPE, DESI_OBS, DESI_SIGMA)):
        if typ == "DM_rd":
            th = D_M(z, Om) / r_d
        elif typ == "DH_rd":
            th = D_H(z, Om) / r_d
        elif typ == "DV_rd":
            th = D_V(z, Om) / r_d
        pull = (th - obs) / sigma
        chi2 += pull**2
        n    += 1
        rows.append((z, typ, th, obs, sigma, pull))
    results[label] = (chi2, n, rows)

# ── Imprimir tabla comparativa ───────────────────────────────────────────────
for label, Om in cases:
    chi2, n, rows = results[label]
    r_d = sound_horizon_rd(Om)
    print(f"\n  {label}   [r_d={r_d:.2f} Mpc]")
    print(f"  {'z':>5}  {'tipo':>6}  {'teórico':>10} {'observado':>10} {'pull':>8}")
    print("  " + "-" * 50)
    for z, typ, th, obs, sigma, pull in rows:
        print(f"  {z:.3f}  {typ:>6}  {th:>10.3f} {obs:>10.2f} {pull:>+8.2f}σ")
    print(f"  χ²_total = {chi2:.3f}  (n={n})   χ²_r = {chi2/n:.4f}")

# ── Veredicto ────────────────────────────────────────────────────────────────
chi2_160, n_160, _ = results["Ω_m=0.160 (dinámico, Paper 2)"]
chi2_320, n_320, _ = results["Ω_m=0.320 (CMB/MIRA, 3D real)"]

delta = chi2_320 - chi2_160
print("\n" + "=" * 68)
print(f"  Δχ² (0.320 − 0.160) = {delta:+.3f}")
print(f"  χ²_r(0.160) = {chi2_160/n_160:.4f}   χ²_r(0.320) = {chi2_320/n_320:.4f}")

# Camino B hipótesis: DESI mediría 0.160 pese a que la realidad es 0.320.
# Para validarse, Ω_m=0.320 debería fallar DESI tanto como 0.160.
# Si en cambio 0.320 PASA y 0.160 FALLA → DESI mide directamente 0.320 (no proyección).

if delta > 10:
    # 0.320 es mucho peor — DESI mide 0.160 directamente → proyección innecesaria
    print("  VEREDICTO: CAMINO B INNECESARIO")
    print("  Ω_m=0.320 falla DESI. DESI mide 0.160 directamente (H(z) radial).")
    print("  No hace falta proyección: la materia real ES 0.160.")
elif delta < -10:
    # 0.160 es mucho peor — DESI mide 0.320, NO proyecta 0.320→0.160
    print("  VEREDICTO: CAMINO B FALSIFICADO — DESI MIDE 0.320, NO 0.160")
    print(f"  Con los w₀,wₐ de SSEE, Ω_m=0.320 reproduce DESI (χ²_r={chi2_320/n_320:.2f})")
    print(f"  y Ω_m=0.160 falla catastrófico (χ²_r={chi2_160/n_160:.1f}).")
    print("  H(z) radial de DESI es sensible directamente a Ω_m — no es proyectable.")
    print("  DESI y Planck coinciden: la materia es 0.320, no una proyección de 0.320.")
    print("  → La identidad Ω_m,dyn=1+w₀=0.160 NO es el Ω_m cosmológico del fondo.")
else:
    print("  VEREDICTO: AMBIGUO (Δχ² pequeño — ambos Ω_m comparables)")
print("=" * 68)
