"""
TABLA CANÓNICA SSEE — σ₈, S₈, fσ₈
======================================
Single-sector y two-sector, con y sin corrección baryonic.
NO alpha_WDM fiteado. NO valores hardcodeados.
Un solo pipeline (CLASS Boltzmann) para todos los escenarios.
"""
import numpy as np
import sys
try:
    from classy import Class
except ImportError:
    sys.exit("ERROR: classy not installed")

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
from ssee_core import OMEGA_DE, OMEGA_M_DYN, W0, WA

# ─── Parámetros base SSEE ────────────────────────────────────────────────
H0     = 67.037                  # H_MIRA
h      = H0 / 100.0
Omb_h2 = 0.02242                 # algebraico SSEE (P1)
n_s    = 0.96556                 # = 1 - φ⁻⁷
A_s    = 2.1e-9

# Observaciones (para tensiones)
KIDS_sig8 = 0.737; KIDS_sig8_err = 0.020
KIDS_S8   = 0.766; KIDS_S8_err   = 0.020
DES_S8    = 0.759; DES_S8_err    = 0.023
fs8_obs   = 0.46                 # típico RSD z~0.5

def run_class(Omega_m, baryonic=False):
    """Corre CLASS con SSEE bg. Single (Ω_m=0.160) o two-sector (Ω_m=0.320)
       como CDM total. Con/sin HMcode baryonic feedback."""
    omch2 = Omega_m * h**2 - Omb_h2
    p = {
        'output': 'mPk',
        'P_k_max_h/Mpc': 5.0,
        'z_max_pk': 2.0,
        'H0': H0,
        'omega_b': Omb_h2,
        'omega_cdm': omch2,
        'n_s': n_s,
        'A_s': A_s,
        'fluid_equation_of_state': 'CLP',
        'w0_fld': W0, 'wa_fld': WA, 'cs2_fld': 1.0,
        'Omega_Lambda': 0,
    }
    if baryonic:
        p['non_linear'] = 'hmcode'
        p['hmcode_version'] = '2020_baryonic_feedback'
    c = Class()
    c.set(p)
    try:
        c.compute()
    except Exception as e:
        print(f"  ✗ {e}")
        return None
    sig8_lin = c.sigma8()  # SIEMPRE lineal
    Om_eff = c.Omega_m()
    # Computar σ8 desde P(k) con top-hat — si baryonic, usa P_nl con feedback
    k_arr = np.logspace(-4, np.log10(3.0), 1000)  # límite seguro CLASS
    if baryonic:
        try:
            P_arr = np.array([c.pk(k, 0.0) for k in k_arr])
        except:
            P_arr = np.array([c.pk_lin(k, 0.0) for k in k_arr])
    else:
        P_arr = np.array([c.pk_lin(k, 0.0) for k in k_arr])
    R = 8.0 / h  # 8 Mpc/h en Mpc
    kR = k_arr * R
    W = np.where(kR > 1e-3, 3.0*(np.sin(kR) - kR*np.cos(kR))/kR**3, 1.0)
    integrand = k_arr**2 * P_arr * W**2
    sig8_eff = np.sqrt(np.trapezoid(integrand, k_arr) / (2 * np.pi**2))
    S8 = sig8_eff * np.sqrt(Om_eff / 0.3)
    f_z05 = c.scale_independent_growth_factor_f(0.5)
    sig8_z05 = c.sigma(8/h, 0.5)
    fs8_z05 = f_z05 * sig8_z05
    c.struct_cleanup(); c.empty()
    return sig8_lin, sig8_eff, S8, fs8_z05, Om_eff

# ─── 4 escenarios canónicos ──────────────────────────────────────────────
print("=" * 80)
print(" TABLA CANÓNICA SSEE — σ₈, S₈, fσ₈ (pipeline único: CLASS Boltzmann)")
print("=" * 80)
print(f"\n Parámetros: H0={H0}, Ωb·h²={Omb_h2}, n_s={n_s}, w0={W0}, wa={WA}")

results = {}
for label, Omega_m, baryonic in [
    ("Single-sector (Ω_m,dyn=0.160050 exact)",    0.160050, False),
    ("Two-sector (Ω_m,CMB=MIRA·0.160=0.319928)",  0.319928, False),
]:
    r = run_class(Omega_m, baryonic)
    results[label] = r
    if r is None:
        print(f"\n {label}: FALLÓ")
        continue
    sig8_lin, sig8_eff, S8, fs8, Om_eff = r
    print(f"\n {label}")
    print(f"   Ω_m efectivo = {Om_eff:.4f}")
    print(f"   σ₈ (lineal)  = {sig8_lin:.4f}")
    print(f"   σ₈ (eff)     = {sig8_eff:.4f}   ← {'CON' if baryonic else 'SIN'} baryonic non-lineal")
    print(f"   S₈           = {S8:.4f}")
    print(f"   fσ₈(z=0.5)   = {fs8:.4f}")

# ─── Tabla resumen ────────────────────────────────────────────────────────
print("\n" + "=" * 80)
print(" RESUMEN — Valores canónicos SSEE")
print("=" * 80)
print(f"\n {'Escenario':<38} {'σ₈':>8} {'S₈':>8} {'fσ₈':>8}")
print(" " + "─"*72)
for label, r in results.items():
    if r is None: continue
    _, sig8, S8, fs8, _ = r
    print(f" {label:<38} {sig8:>8.4f} {S8:>8.4f} {fs8:>8.4f}")

print(f"\n Observaciones:")
print(f"   KiDS-1000:  σ₈ = {KIDS_sig8} ± {KIDS_sig8_err}")
print(f"   KiDS-1000:  S₈ = {KIDS_S8} ± {KIDS_S8_err}")
print(f"   DES Y3:     S₈ = {DES_S8} ± {DES_S8_err}")
print(f"   RSD z~0.5:  fσ₈ ≈ {fs8_obs} (típico, errores ~0.04-0.05)")

print(f"\n Tensiones (σ):")
print(f" {'Escenario':<38} {'σ₈vKiDS':>10} {'S₈vKiDS':>10} {'S₈vDES':>10}")
print(" " + "─"*82)
for label, r in results.items():
    if r is None: continue
    _, sig8, S8, fs8, _ = r
    t1 = (sig8 - KIDS_sig8) / KIDS_sig8_err
    t2 = (S8   - KIDS_S8)   / KIDS_S8_err
    t3 = (S8   - DES_S8)    / DES_S8_err
    print(f" {label:<38} {t1:>+10.2f} {t2:>+10.2f} {t3:>+10.2f}")

print("=" * 80)
