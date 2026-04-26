"""
SSEE-V3.6 — Verificación numérica del horizonte de sonido (task 2A)

Compara r_d y theta_* de SSEE contra los valores MEDIDOS por Planck 2018
(Planck Collaboration 2020, A&A 641, A6, Tabla 2).

No se usa ΛCDM como árbitro. El árbitro es la observación directa.

Casos evaluados:
  A) SSEE + MIRA (Ωb h² estándar, Planck 2018)
  B) SSEE + MIRA + IS (Ωb h² derivado algebraicamente en Paper 4)
  C) SSEE dinámico puro (Ωm = 0.160, sin MIRA) — muestra por qué falla
"""

import numpy as np

# ---------------------------------------------------------------------------
# Constantes algebraicas SSEE (todas de φ y π, sin ajuste)
# ---------------------------------------------------------------------------
phi   = (1 + 5**0.5) / 2          # 1.61803…
pi    = np.pi
Omega = pi + phi                   # 4.75963…  Stability Metric
beta  = (pi + phi) / 2             # 2.37981…  Base Coupling Scalar
KAL0  = beta + pi                  # 5.52141…  Structural Viscosity
P_sc  = Omega + phi                # 6.37766…  Dynamical Evolution Scalar
Kv    = phi + pi + Omega           # 9.51926…  Structural Constraint
Tr    = 3 * (phi + beta)           # 11.99353… 3D Saturation Horizon
Mv    = phi + pi + Kv              # 14.27889… Maximal Dimensional Invariant

w0    = -Tr / Mv                   # -0.84027 (sector dinámico)
wa    = -P_sc / Kv                 # -0.67030
OmDE  = Tr / Mv                    # 0.84027
Omm_dyn = 1.0 - OmDE              # 0.15973 (sector dinámico BAO/cúmulos)

BIAL  = (phi + pi) / 2
AURA  = phi + BIAL
MIRA  = AURA / 2                   # 1.99892 (Genesis 5.12, pre-data 2026-01-28)
Omm_cmb = Omm_dyn * MIRA          # 0.31939 (sector CMB/alto-z)

# H0 algebraico (Paper 4)
H0_ssee = 3 * (phi + pi)**2        # ≈ 67.96 km/s/Mpc
H0_run  = 66.66                    # valor MCMC Paper 2 (best-fit posterior)

# Ωb h² — dos casos:
#   Estándar: Planck 2018 prior (Paper 2)
#   IS: derivado algebraicamente en Paper 4 §sec:baryons
#     Ωb = 100(π−φ)/[6(φ+π)⁴] = 0.0495 (densidad fraccional)
#     Ωb h² = 3(π−φ)/200 = 0.02285 (con H0 algebraico = 3(φ+π)² ≈ 67.96)
#     Tensión con Planck: |0.02285−0.02237|/0.00015 = 3.2σ — documentado en Paper 4 Tabla
#   NOTA: corregido en sesión 8 — todos los docs ahora usan 0.02285
Ombh2_std = 0.02237
Ombh2_IS  = 3 * (pi - phi) / 200               # 0.02285 (Paper 4 §sec:baryons)

# n_s algebraico (Paper 4): 1 − φ⁻⁷
ns_ssee = 1.0 - phi**(-7)          # 0.96556

# Observables Planck 2018 MEDIDOS (Planck 2020, Tabla 2, TT+TE+EE+lowE+lensing)
# Estos son los árbitros — no ΛCDM
RD_OBS   = 147.09   # Mpc  (valor central)
RD_SIGMA =   0.26   # Mpc  (1σ)
THETA_OBS   = 0.59668  # grados (100θ_* = 1.04110 → θ_* = 0.010411 rad = 0.59668°)
THETA_SIGMA = 0.00046  # grados


def run_camb_case(H0_val, ombh2, Omm_val, w0_val, wa_val, ns_val,
                  label, mnu=0.06):
    """
    Calcula r_d, θ* y z_drag via CAMB.

    r_d es independiente de la reionización (z_drag~1060 >> z_reion~8),
    por eso desactivamos Reionization para evitar fallos de convergencia
    con parámetros IS (Ωb h² ligeramente distinto al prior de Planck).
    """
    try:
        import camb
    except ImportError:
        raise ImportError("camb no instalado: pip install camb")

    h      = H0_val / 100.0
    omch2  = Omm_val * h**2 - ombh2
    if omch2 < 0:
        print(f"  [{label}] omch2 = {omch2:.5f} < 0 — no computable")
        return None

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0_val, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=0)
    pars.Reion.Reionization = False   # r_d insensible a reionización
    pars.set_dark_energy(w=w0_val, wa=wa_val, dark_energy_model="ppf")
    pars.InitPower.set_params(ns=ns_val, As=np.exp(3.044) * 1e-10)
    pars.set_for_lmax(2500, lens_potential_accuracy=0)
    pars.Want_CMB = True
    pars.WantTensors = False

    results  = camb.get_results(pars)
    d        = results.get_derived_params()

    r_d      = d["rdrag"]
    theta_s  = d.get("thetastar", None)
    z_drag   = d.get("zdrag", None)
    # CAMB devuelve thetastar como 100×θ* (puro, ≈1.041) — no en radianes
    theta_deg = np.degrees(theta_s / 100) if theta_s is not None else None

    return {
        "r_d":       r_d,
        "theta_deg": theta_deg,
        "z_drag":    z_drag,
        "omch2":     omch2,
    }


def tension(val, obs, sigma):
    return abs(val - obs) / sigma


def print_comparison(label, res):
    if res is None:
        print(f"  [{label}] No se pudo calcular (omch2 < 0)")
        return

    r_d   = res["r_d"]
    t_rd  = tension(r_d, RD_OBS, RD_SIGMA)
    print(f"\n  [{label}]")
    print(f"    r_d         = {r_d:.3f} Mpc   (obs: {RD_OBS:.2f} ± {RD_SIGMA:.2f} Mpc,  tensión: {t_rd:.2f}σ)")

    if res["theta_deg"] is not None:
        t_th = tension(res["theta_deg"], THETA_OBS, THETA_SIGMA)
        print(f"    θ*          = {res['theta_deg']:.5f}°  (obs: {THETA_OBS:.5f} ± {THETA_SIGMA:.5f}°, tensión: {t_th:.2f}σ)")

    if res["z_drag"] is not None:
        print(f"    z_drag      = {res['z_drag']:.1f}")
    if res["omch2"] is not None:
        print(f"    Ωc h²       = {res['omch2']:.5f}")


def main():
    print("=" * 65)
    print("SSEE-V3.6 — Verificación numérica del horizonte de sonido r_d")
    print("Referencia observacional: Planck 2018, Tabla 2")
    print(f"  r_d^obs   = {RD_OBS} ± {RD_SIGMA} Mpc")
    print(f"  θ*^obs    = {THETA_OBS} ± {THETA_SIGMA} °")
    print("=" * 65)

    # -----------------------------------------------------------------------
    # Caso A: SSEE + MIRA, Ωb h² estándar (Planck 2018 prior)
    # -----------------------------------------------------------------------
    A = run_camb_case(
        H0_val=H0_run, ombh2=Ombh2_std, Omm_val=Omm_cmb,
        w0_val=w0, wa_val=wa, ns_val=ns_ssee,
        label="A: SSEE+MIRA, Ωb h²=0.02237 (Planck prior)")

    # -----------------------------------------------------------------------
    # Caso B: SSEE + MIRA, Ωb h² IS (derivado algebraicamente en Paper 4)
    # -----------------------------------------------------------------------
    B = run_camb_case(
        H0_val=H0_run, ombh2=Ombh2_IS, Omm_val=Omm_cmb,
        w0_val=w0, wa_val=wa, ns_val=ns_ssee,
        label=f"B: SSEE+MIRA+IS, Ωb h²={Ombh2_IS:.5f} (Paper 4 algebraic)")

    # -----------------------------------------------------------------------
    # Caso C: SSEE dinámico puro (Ωm = 0.160, sin MIRA) — falsificación
    # -----------------------------------------------------------------------
    C = run_camb_case(
        H0_val=H0_run, ombh2=Ombh2_std, Omm_val=Omm_dyn,
        w0_val=w0, wa_val=wa, ns_val=ns_ssee,
        label="C: SSEE dinámico (Ωm=0.160, sin MIRA) — caso naive")

    print("\n--- Resultados vs observación directa Planck 2018 ---")
    print_comparison("A: SSEE+MIRA, Ωb h²=0.02237", A)
    print_comparison(f"B: SSEE+MIRA+IS, Ωb h²={Ombh2_IS:.5f}", B)
    print_comparison("C: SSEE dinámico (Ωm=0.160)", C)

    print("\n--- Resumen ---")
    print(f"  Ωb h² IS algebraico  = (π−φ)/[6(φ+π)⁴] = {Ombh2_IS:.5f}")
    print(f"  n_s algebraico       = 1 − φ⁻⁷ = {ns_ssee:.5f}")
    print(f"  Ωm,dyn               = {Omm_dyn:.5f}")
    print(f"  MIRA                 = {MIRA:.6f}")
    print(f"  Ωm,CMB               = {Omm_cmb:.5f}")
    print(f"  w₀                   = {w0:.5f}")
    print(f"  wₐ                   = {wa:.5f}")

    if A is not None and B is not None:
        delta_rd = A["r_d"] - B["r_d"]
        print(f"\n  Δr_d (A−B, IS correction) = {delta_rd:+.3f} Mpc")
        print(f"  → IS Ωb h² shift moves r_d {'toward' if B['r_d'] < A['r_d'] and A['r_d'] > RD_OBS else 'away from'} observation")
    print("=" * 65)


if __name__ == "__main__":
    main()
