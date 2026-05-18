#!/usr/bin/env python3
"""
ssee_core.py — Fuente única de verdad para SSEE-V3.6
=====================================================
TODO script de SSEE (src/ y class_ssee/) DEBE importar de aquí sus constantes
algebraicas, sus datos observacionales y la función de fondo de Friedmann.

NUNCA redefinir localmente phi, w0, los datos fsigma8, etc. La divergencia entre
tres copias de la ODE de crecimiento (verification / plot_calibrated / mcmc) y el
bug E²(z=0)=1.16 de Paper 6 ocurrieron exactamente por copias locales divergentes.

Uso desde src/:           from ssee_core import PHI, W0, friedmann_E2, load_fsigma8
Uso desde class_ssee/:    import sys; sys.path.insert(0, '../src'); import ssee_core

Verificación rápida:      python3 src/ssee_core.py
  (imprime todas las constantes y corre todos los asserts de cordura).
"""
import math
import os

# ── Constantes fundamentales ─────────────────────────────────────────────────
PHI = (1.0 + math.sqrt(5.0)) / 2.0          # razón áurea
PI  = math.pi

# ── Registros estructurales derivados (Paper 1) ──────────────────────────────
OMEGA = PHI + PI                # Stability Metric        ≈ 4.7596
BETA  = (PHI + PI) / 2.0        # Base Coupling Scalar    ≈ 2.3798
KAL0  = BETA + PI               # Structural Viscosity    ≈ 5.5214
P_SC  = OMEGA + PHI             # Dynamical Evolution     ≈ 6.3776
K_V   = PHI + PI + OMEGA        # Structural Constraint   ≈ 9.5192
T_R   = 3.0 * (PHI + BETA)      # 3D Saturation Horizon   ≈ 11.9935
M_V   = PHI + PI + K_V          # Maximal Dim. Invariant  ≈ 14.2788

# ── Ecuación de estado y densidades del sector dinámico ──────────────────────
W0          = -T_R / M_V        # ≈ -0.840
WA          = -P_SC / K_V       # ≈ -0.670
OMEGA_DE    = T_R / M_V         # ≈ 0.840
OMEGA_M_DYN = 1.0 + W0          # ≈ 0.160  (sector dinámico: BAO, H(z), Paper 1/2)

# ── MIRA / AURA ──────────────────────────────────────────────────────────────
MIRA = (3.0*PHI + PI) / 4.0     # ≈ 1.9989  (hipótesis auxiliar, no derivada)
AURA = (3.0*PHI + PI) / 2.0     # ≈ 3.9978  (= 2*MIRA = PHI+BETA)

# ── Ω_m para el CMB — DOS definiciones, ~0.05% de diferencia ─────────────────
# NO unificar en silencio: la dualidad es un problema abierto (CLAUDE.md Roadmap pt.2).
OMEGA_M_CMB_MIRA      = MIRA * OMEGA_M_DYN          # ≈ 0.31983  (vía MIRA)
OMEGA_M_CMB_GEOMETRIC = (PI - PHI) / (PI + PHI)     # ≈ 0.32010  (geométrico)

# ── Otros observables algebraicos ────────────────────────────────────────────
H0_ALG     = 3.0 * OMEGA**2          # ≈ 67.96 km/s/Mpc  = 3(phi+pi)^2
N_S        = 1.0 - PHI**(-7)         # ≈ 0.96556
R_TENSOR   = PHI**(-10)              # ≈ 0.00813
OMEGA_B_H2 = (PI - PHI) / (3.0 * OMEGA**2)   # ≈ 0.02242  (OP-1)

# ── Valores ΛCDM de referencia (Planck 2018, A&A 641 A6) ─────────────────────
LCDM_OMEGA_M = 0.3153
LCDM_H0      = 67.36
LCDM_SIGMA8  = 0.8111
LCDM_NS      = 0.9649

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', 'data', 'raw')


# ── Fondo de Friedmann — ÚNICA implementación válida ─────────────────────────
def friedmann_E2(a, Omega_m, Omega_de=None, w0=W0, wa=WA, Omega_r=0.0):
    """E²(a) = H²(a)/H₀² para un fondo con energía oscura CPL.

    Si Omega_de se omite, se fija a 1 - Omega_m - Omega_r, de modo que
    Σ Ω = 1 por construcción y E²(a=1) == 1 SIEMPRE. El bug de Paper 6
    (E²(z=0)=1.16) vino de pasar Omega_m=0.320 y Omega_de=0.840 a la vez:
    nunca pasar ambos salvo que su suma + Omega_r sea exactamente 1.
    """
    if Omega_de is None:
        Omega_de = 1.0 - Omega_m - Omega_r
    f_de = a**(-3.0*(1.0 + w0 + wa)) * math.exp(-3.0*wa*(1.0 - a))
    return Omega_r*a**-4 + Omega_m*a**-3 + Omega_de*f_de


# ── Cargadores de datos observacionales (con procedencia en el archivo) ──────
def load_fsigma8():
    """Devuelve la lista canónica fσ₈(z): dicts z, fsigma8, sigma, survey, ref."""
    path = os.path.join(_DATA_DIR, 'fsigma8_rsd.csv')
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = line.split(',')
                continue
            vals = line.split(',')
            row = dict(zip(header, vals))
            rows.append({
                'z':       float(row['z_eff']),
                'fsigma8': float(row['fsigma8']),
                'sigma':   float(row['sigma']),
                'survey':  row['survey'],
                'ref':     row['reference'],
            })
    return rows


# ── Asserts de cordura — corren al importar el módulo ────────────────────────
def _sanity_checks():
    assert abs(PHI**2 - (PHI + 1.0)) < 1e-12, "identidad áurea phi^2=phi+1 rota"
    assert abs(AURA - 2.0*MIRA) < 1e-12,      "AURA debe ser 2*MIRA"
    assert abs(OMEGA_M_DYN + OMEGA_DE - 1.0) < 1e-12, "Omega_m,dyn + Omega_DE != 1"
    assert round(W0, 3) == -0.840,  f"w0 fuera de rango: {W0}"
    assert round(WA, 3) == -0.670,  f"wa fuera de rango: {WA}"
    assert round(OMEGA_M_DYN, 3) == 0.160, f"Omega_m,dyn fuera de rango: {OMEGA_M_DYN}"
    assert 0.0 < N_S < 1.0,         f"n_s no físico: {N_S}"
    assert 66.0 < H0_ALG < 69.0,    f"H0_alg fuera de rango: {H0_ALG}"
    assert abs(friedmann_E2(1.0, OMEGA_M_DYN) - 1.0) < 1e-9, "E^2(a=1) != 1"
    assert abs(friedmann_E2(1.0, OMEGA_M_CMB_MIRA) - 1.0) < 1e-9, "E^2(a=1) != 1"


_sanity_checks()


def check():
    """Imprime todas las constantes y confirma los asserts. Paso de auditoría."""
    print("=" * 64)
    print("ssee_core.py — constantes canónicas SSEE-V3.6")
    print("=" * 64)
    for name in ('PHI', 'PI', 'OMEGA', 'BETA', 'KAL0', 'P_SC', 'K_V', 'T_R',
                 'M_V', 'W0', 'WA', 'OMEGA_DE', 'OMEGA_M_DYN', 'MIRA', 'AURA',
                 'OMEGA_M_CMB_MIRA', 'OMEGA_M_CMB_GEOMETRIC', 'H0_ALG', 'N_S',
                 'R_TENSOR', 'OMEGA_B_H2'):
        print(f"  {name:24s} = {globals()[name]:.10f}")
    fs = load_fsigma8()
    print(f"\n  fsigma8 canónico: {len(fs)} puntos "
          f"(z={fs[0]['z']}..{fs[-1]['z']}, surveys {sorted(set(r['survey'] for r in fs))})")
    _sanity_checks()
    print("\n  Todos los asserts de cordura: OK")


if __name__ == '__main__':
    check()
