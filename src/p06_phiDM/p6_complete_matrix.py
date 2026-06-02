"""
P6 complete matrix — σ₈/S₈ con y sin correcciones, single y two-sector
========================================================================
Llena las celdas que faltaban en la tabla comparativa de Mike.
Usa CLASS uniformemente para todos los escenarios.
"""
import numpy as np
import sys
try:
    from classy import Class
except ImportError:
    print("ERROR: classy not installed"); sys.exit(1)

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
from ssee_core import PHI, PI, OMEGA_DE, OMEGA_M_DYN, MIRA

# ── Parámetros base SSEE ─────────────────────────────────────────────────────
H0_ssee   = 67.068                # H_MIRA (más reciente)
Omega_b   = 0.0492
h         = H0_ssee / 100.0
Omb_h2    = Omega_b * h**2
n_s       = 0.96556
w0, wa    = -OMEGA_DE, -0.670

# Observaciones
KIDS_sig8 = 0.737; KIDS_sig8_err = 0.020
KIDS_S8   = 0.766; KIDS_S8_err   = 0.020
DES_S8    = 0.759; DES_S8_err    = 0.023

# Para fσ8 comparación — del paper 5/6 verificado ya
fs8_obs_mean = 0.46   # típico z~0.5

# ── Función auxiliar ─────────────────────────────────────────────────────────
def run_class(scenario_name, Omega_m_total, use_ncdm=False, use_baryonic=False, m_phi_eV=5.602):
    Omega_cdm = Omega_m_total - Omega_b
    params = {
        'output': 'mPk',
        'P_k_max_h/Mpc': 5.0,
        'z_max_pk': 2.0,
        'H0': H0_ssee,
        'omega_b': Omb_h2,
        'omega_cdm': Omega_cdm * h**2,
        'n_s': n_s,
        'A_s': 2.1e-9,
        'fluid_equation_of_state': 'CLP',
        'w0_fld': w0, 'wa_fld': wa, 'cs2_fld': 1.0,
        'Omega_Lambda': 0,
    }
    if use_ncdm:
        # Two-sector: φ-DM como ncdm warm
        # Pero aquí queremos partition CDM/φDM
        Om_phi = Omega_m_total - OMEGA_M_DYN
        Om_cdm_only = OMEGA_M_DYN - Omega_b
        params['omega_cdm'] = Om_cdm_only * h**2
        params['N_ncdm'] = 1
        params['m_ncdm'] = m_phi_eV
        params['T_ncdm'] = 0.716
        params['omega_ncdm'] = Om_phi * h**2
    if use_baryonic:
        params['non_linear'] = 'hmcode'
        params['hmcode_version'] = '2020_baryonic_feedback'

    c = Class()
    c.set(params)
    try:
        c.compute()
    except Exception as e:
        print(f"   ✗ {scenario_name}: {e}")
        c.struct_cleanup(); c.empty()
        return None, None, None
    sigma8 = c.sigma8()
    Omega_m_eff = c.Omega_m()
    S8 = sigma8 * np.sqrt(Omega_m_eff / 0.3)
    fs8_z05 = c.scale_independent_growth_factor_f(0.5) * c.sigma(8/h, 0.5)
    c.struct_cleanup(); c.empty()
    return sigma8, S8, fs8_z05

# ─── 4 ESCENARIOS DE LA TABLA ──────────────────────────────────────────────
print("=" * 78)
print(" Matriz completa SSEE — σ₈, S₈, fσ₈(z=0.5)  con y sin correcciones")
print("=" * 78)

results = {}

print("\n[A] Single-sector lineal (Ω_m=0.160, NO correcciones)")
sig, S8, fs8 = run_class("A", 0.160, use_ncdm=False, use_baryonic=False)
results['A'] = (sig, S8, fs8)
print(f"    σ₈ = {sig:.4f}   S₈ = {S8:.4f}   fσ₈(0.5) = {fs8:.4f}")

print("\n[B] Single-sector + HMcode baryonic (Ω_m=0.160, CON baryonic)")
sig, S8, fs8 = run_class("B", 0.160, use_ncdm=False, use_baryonic=True)
results['B'] = (sig, S8, fs8)
if sig: print(f"    σ₈ = {sig:.4f}   S₈ = {S8:.4f}   fσ₈(0.5) = {fs8:.4f}")

print("\n[C] Two-sector lineal todo-CDM (Ω_m=0.320, NO correcciones, NO WDM)")
sig, S8, fs8 = run_class("C", 0.320, use_ncdm=False, use_baryonic=False)
results['C'] = (sig, S8, fs8)
print(f"    σ₈ = {sig:.4f}   S₈ = {S8:.4f}   fσ₈(0.5) = {fs8:.4f}")

print("\n[D] Two-sector + HMcode baryonic todo-CDM (Ω_m=0.320, CON baryonic, NO WDM)")
sig, S8, fs8 = run_class("D", 0.320, use_ncdm=False, use_baryonic=True)
results['D'] = (sig, S8, fs8)
if sig: print(f"    σ₈ = {sig:.4f}   S₈ = {S8:.4f}   fσ₈(0.5) = {fs8:.4f}")

print("\n[E] Two-sector + WDM φ-DM (Ω_m=0.320, NO baryonic, CON WDM ncdm)")
sig, S8, fs8 = run_class("E", 0.320, use_ncdm=True, use_baryonic=False)
results['E'] = (sig, S8, fs8)
if sig: print(f"    σ₈ = {sig:.4f}   S₈ = {S8:.4f}   fσ₈(0.5) = {fs8:.4f}")

print("\n[F] Two-sector + WDM + baryonic (FULL P6 + OP-5)")
sig, S8, fs8 = run_class("F", 0.320, use_ncdm=True, use_baryonic=True)
results['F'] = (sig, S8, fs8)
if sig: print(f"    σ₈ = {sig:.4f}   S₈ = {S8:.4f}   fσ₈(0.5) = {fs8:.4f}")

# ─── REPORTE FINAL ─────────────────────────────────────────────────────────
print("\n" + "=" * 78)
print(" TABLA FINAL (tensiones vs observaciones)")
print("=" * 78)
print(f"\n{'Escenario':<55} {'σ₈':>7} {'S₈':>7} {'fσ₈':>7}")
print(" " + "─"*76)
labels = {
    'A': "[A] Single-sector lineal (no corr.)",
    'B': "[B] Single-sector + baryonic feedback",
    'C': "[C] Two-sector lineal todo-CDM (no WDM, no bar)",
    'D': "[D] Two-sector + baryonic (no WDM)",
    'E': "[E] Two-sector + WDM φ-DM (no baryonic)",
    'F': "[F] Two-sector + WDM + baryonic (full P6+OP-5)",
}
for k, lbl in labels.items():
    sig, S8, fs8 = results[k]
    if sig is None:
        print(f" {lbl:<55} {'FAIL':>7} {'FAIL':>7} {'FAIL':>7}")
        continue
    print(f" {lbl:<55} {sig:>7.4f} {S8:>7.4f} {fs8:>7.4f}")

print(f"\n Observaciones:")
print(f"   KiDS:  σ₈={KIDS_sig8}±{KIDS_sig8_err},  S₈={KIDS_S8}±{KIDS_S8_err}")
print(f"   DES:   S₈={DES_S8}±{DES_S8_err}")
print(f"   fσ₈(z~0.5): ~0.46 (típico RSD)")

print(f"\n Tensiones σ vs observaciones:")
print(f"{'Escenario':<55} {'σ₈vsKiDS':>10} {'S₈vsKiDS':>10} {'S₈vsDES':>10}")
print(" " + "─"*86)
for k, lbl in labels.items():
    sig, S8, fs8 = results[k]
    if sig is None:
        continue
    t_sig = (sig - KIDS_sig8) / KIDS_sig8_err
    t_S8K = (S8 - KIDS_S8) / KIDS_S8_err
    t_S8D = (S8 - DES_S8) / DES_S8_err
    print(f" {lbl:<55} {t_sig:>+10.2f} {t_S8K:>+10.2f} {t_S8D:>+10.2f}")

print("\n" + "=" * 78)
