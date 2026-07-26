#!/usr/bin/env python3
"""
ssee_core.py — Fuente única de verdad para SSEE
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

# ── Ecuación de estado y densidades ──────────────────────────────────────────
W0          = -T_R / M_V        # ≈ -0.840
IGNIS       = PI + P_SC         # IGNIS = π+PYROS         ≈ 9.5192
WA          = -P_SC / IGNIS     # ≈ -0.670  (denominador IGNIS, NUNCA el
                                #  scaffold K_V: mismo valor 2Ω, otra entidad)
OMEGA_DE    = T_R / M_V         # ≈ 0.840  (= |w0|; es la EoS, NO la fracción de densidad)
#
# ⚠️ REGLA DE Ω_m (2026-07-09, tras el hallazgo del χ²=726, V-L4-DESI):
#   La GEOMETRÍA de fondo — E(z), distancias BAO, H(z), r_d, cualquier E²(z) —
#   usa SIEMPRE la materia TOTAL:  OMEGA_M_TOTAL = ω_m/h² = 0.30889.
#   El 0.160 NO es una densidad de materia: es el SECTOR frío dinámico (1+w0),
#   componente de perturbaciones (Paper 6: 0.160 + φ-DM 0.149 = 0.308) y factor
#   de la fórmula EFT α_K = 3|w0|·0.160. NUNCA va en un E(z).
#   Meter el sector (0.160) en la geometría fue el bug que dio χ²=726 en DESI DR2.
OMEGA_CDM_SECTOR = 1.0 + W0      # ≈ 0.160  sector frío dinámico (Paper 6 + α_K); NO es materia total
OMEGA_M_DYN      = OMEGA_CDM_SECTOR   # [ALIAS DEPRECADO] usar OMEGA_CDM_SECTOR o, si es geometría, OMEGA_M_TOTAL

# ── MIRA / AURA ──────────────────────────────────────────────────────────────
# MIRA persiste como ENTIDAD (= AURA/2, valor 1.9989); conserva su rol en
# f_screen y la identidad AURA = 2·MIRA. Lo que se reasignó (reframe 2026-06-17)
# es su ROL de factor-materia → pasa a π/φ. Ver project_halg_pifi_investigation.
MIRA = (3.0*PHI + PI) / 4.0     # ≈ 1.9989  (entidad; f_screen; = AURA/2)
AURA = (3.0*PHI + PI) / 2.0     # ≈ 3.9978  (= 2*MIRA = PHI+BETA)

# ── Observables algebraicos de fondo ─────────────────────────────────────────
H0_ALG     = 3.0 * OMEGA**2          # ≈ 67.96 km/s/Mpc  = 3(phi+pi)^2
H0_GLOBAL  = H0_ALG                  # CANÓNICO 2026-06-17: H global de fondo = H_alg
H0_MIRA    = 67.037                  # ancla CMB-fit del escenario VIEJO (cascada pendiente re-run)
N_S        = 1.0 - PHI**(-7)         # ≈ 0.96556
R_TENSOR   = PHI**(-10)              # ≈ 0.00813
OMEGA_B_H2 = (PI - PHI) / (3.0 * OMEGA**2)   # ≈ 0.02242  (OP-1)
# Σm_ν^active = R₂·ω_b·C_ν/(τ_Π·H₀) con C_ν=93.14 eV PDG (N_eff=3.046).
# 2026-07-25: era 0.06902, residuo de C_ν=94.07. El cambio C_ν 94.07→93.14
# (2026-07-10) propagó a m_φ (41.02→40.70) pero NO a esta constante ni a sus
# consumidores. Σm_ν es UNA sola cantidad: la misma que da m_φ y la que da ω_ν.
# Coherencia vigilada por V-L2-11a/b/c en ssee_verify.py.
SUM_MNU_EV = 0.06849                          # Σm_ν activos (canónico; = 40.70/594.28)

# ── Densidad de materia CMB — ω_m DIRECTO (reframe 2026-06-18) ────────────────
# OP-8 CERRADO: no hay "factor materia" que derivar. Ω_m,CMB es el observable
# físico estándar ω_m/h², con cada pieza algebraica de SSEE:
#   ω_b = (π-φ)/(3Ω²)              (OP-1)
#   ω_c = KAL₀·ω_b·n_s            (identidad forward, ya en Paper 1; 0.41σ)
#   ω_ν = Σm_ν/93.14 eV           (neutrinos activos)
# Ω_m,CMB = ω_m/h² es DERIVADO, no postulado. El viejo factor π/φ (y antes MIRA)
# queda RETIRADO: Ω_m,dyn=0.160 (DESI) y ω_m (CMB) son dos predicciones
# independientes, ya no ligadas por un factor. Ver project_halg_pifi_investigation.
OMEGA_C_H2 = KAL0 * OMEGA_B_H2 * N_S           # ≈ 0.11952  (forward, Paper 1)
OMEGA_NU_H2 = SUM_MNU_EV / 93.14               # ≈ 0.000735
OMEGA_M_H2 = OMEGA_B_H2 + OMEGA_C_H2 + OMEGA_NU_H2   # ≈ 0.14268  (ω_m físico)
# MATERIA TOTAL — la ÚNICA densidad que entra en cualquier geometría de fondo:
OMEGA_M_TOTAL = OMEGA_M_H2 / (H0_GLOBAL/100.0)**2    # ≈ 0.30889  (CANÓNICO, derivado, ω_m/h²)
OMEGA_M_CMB   = OMEGA_M_TOTAL                         # [ALIAS] mismo número; usar OMEGA_M_TOTAL
# Identidades históricas (RETIRADAS como factor-materia, conservadas para trazar):
OMEGA_M_CMB_PIPHI     = (PI / PHI) * OMEGA_CDM_SECTOR # ≈ 0.31069  (factor π/φ — superado)
OMEGA_M_CMB_MIRA      = MIRA * OMEGA_CDM_SECTOR       # ≈ 0.31983  (identidad MIRA — histórico)
OMEGA_M_CMB_GEOMETRIC = (PI - PHI) / (PI + PHI)      # ≈ 0.32010  (geométrico)

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
    assert abs(OMEGA_CDM_SECTOR + OMEGA_DE - 1.0) < 1e-12, "sector frío + Omega_DE(|w0|) != 1"
    assert round(W0, 3) == -0.840,  f"w0 fuera de rango: {W0}"
    assert round(WA, 3) == -0.670,  f"wa fuera de rango: {WA}"
    assert round(OMEGA_CDM_SECTOR, 3) == 0.160, f"sector frío fuera de rango: {OMEGA_CDM_SECTOR}"
    assert OMEGA_M_DYN == OMEGA_CDM_SECTOR, "alias OMEGA_M_DYN roto"
    assert 0.0 < N_S < 1.0,         f"n_s no físico: {N_S}"
    assert 66.0 < H0_ALG < 69.0,    f"H0_alg fuera de rango: {H0_ALG}"
    assert abs(friedmann_E2(1.0, OMEGA_CDM_SECTOR) - 1.0) < 1e-9, "E^2(a=1) != 1 (sector)"
    assert abs(friedmann_E2(1.0, OMEGA_M_TOTAL) - 1.0) < 1e-9, "E^2(a=1) != 1 (total)"
    assert round(OMEGA_M_TOTAL, 4) == 0.3089, f"Omega_m,total (ω_m/h²) fuera de rango: {OMEGA_M_TOTAL}"
    assert OMEGA_M_CMB == OMEGA_M_TOTAL, "alias OMEGA_M_CMB roto"
    assert abs(OMEGA_M_TOTAL - OMEGA_M_H2 / (H0_GLOBAL/100.0)**2) < 1e-12, "Omega_m,total debe ser ω_m/h²"
    assert round(OMEGA_M_H2, 4) == 0.1427, f"ω_m algebraico fuera de rango: {OMEGA_M_H2}"
    # Ω_m,CMB reproducible en 3 s: los TRES términos son necesarios.
    # Omitir ω_ν devuelve 0.30729 (no 0.30889) — el ν NO es opcional (Sealed §two-omega).
    _h2 = (H0_GLOBAL / 100.0) ** 2
    _om_sin_nu = (OMEGA_B_H2 + OMEGA_C_H2) / _h2
    assert round(_om_sin_nu, 4) == 0.3073, \
        f"Ω_m,CMB sin ω_ν debe dar ~0.30729 (ν no opcional): {_om_sin_nu}"
    assert round((OMEGA_B_H2 + OMEGA_C_H2 + OMEGA_NU_H2) / _h2, 4) == 0.3089, \
        "Ω_m,CMB con los tres términos debe dar 0.30889"
    assert abs(OMEGA_C_H2 - KAL0 * OMEGA_B_H2 * N_S) < 1e-12, "ω_c debe ser KAL₀·ω_b·n_s"


_sanity_checks()


def check():
    """Imprime todas las constantes y confirma los asserts. Paso de auditoría."""
    print("=" * 64)
    print("ssee_core.py — constantes canónicas SSEE")
    print("=" * 64)
    for name in ('PHI', 'PI', 'OMEGA', 'BETA', 'KAL0', 'P_SC', 'K_V', 'T_R',
                 'M_V', 'W0', 'WA', 'OMEGA_DE', 'OMEGA_CDM_SECTOR', 'OMEGA_M_TOTAL',
                 'OMEGA_M_DYN', 'MIRA', 'AURA',
                 'OMEGA_B_H2', 'OMEGA_C_H2', 'OMEGA_NU_H2', 'OMEGA_M_H2',
                 'OMEGA_M_CMB', 'OMEGA_M_CMB_PIPHI', 'OMEGA_M_CMB_MIRA',
                 'OMEGA_M_CMB_GEOMETRIC', 'H0_ALG', 'H0_GLOBAL', 'N_S',
                 'R_TENSOR', 'SUM_MNU_EV'):
        print(f"  {name:24s} = {globals()[name]:.10f}")
    fs = load_fsigma8()
    print(f"\n  fsigma8 canónico: {len(fs)} puntos "
          f"(z={fs[0]['z']}..{fs[-1]['z']}, surveys {sorted(set(r['survey'] for r in fs))})")
    _sanity_checks()
    print("\n  Todos los asserts de cordura: OK")


if __name__ == '__main__':
    check()
