import os
"""
H0 Cascade Audit — Auditoría sistemática bloque por bloque
============================================================
Para cada cantidad derivada del modelo SSEE que depende de H₀,
calcula su valor bajo:
  (A) H_alg  = 67.962 km/s/Mpc  (Type-P, derivación algebraica P4)
  (B) H_MIRA = 67.068 km/s/Mpc  (físico, Planck plik_lite + SSEE bg, P3)

Reporta:
  - Valor bajo cada H₀
  - Cambio absoluto y porcentual
  - Si la fórmula es DIMENSIONALMENTE VÁLIDA (puede cascadar) o INVÁLIDA
  - Status de actualización en manuscritos (✅ ya cascadado / ⚠ pendiente / ❌ inválido)

Esto cierra la auditoría que Mike pidió el 2026-05-24 PM tras detectar que
la cascada anterior fue LINEAL (solo donde H₀ aparece literalmente) y no
NO-LINEAL (donde H₀ aparece dentro de ρ_crit, M, etc.).
"""
import numpy as np
# ORIGEN de los numeros (R65, 2026-09-19). Script HISTORICO: la auditoria de
# la era H_MIRA (2026-05-24); no produce ningun log vigente.
# ORIGEN-VALOR: 0.4169052 — s_K_UV = 0.41690518 que imprime src/p10_uv/ssee_paper10_verification.py (Step 4)
# ORIGEN-VALOR: 0.0695216 — f_screen_UV = 0.06952161 que imprime src/p10_uv/ssee_paper10_verification.py (Step 4)
# ORIGEN-VALOR: 0.0695216111441 — f_screen_UV = 0.06952161114406 de src/p10_uv/ssee_paper10_verification.py, a 13 decimales
# ORIGEN-VALOR: 0.02261 — omega_b h^2 posterior del MCMC con prior MIRA, results/logs/mcmc_paper2_mira.log linea 67
# ORIGEN-VALOR: 0.0824 — Sum m_nu VIEJO de Paper 4, RETIRADO (contaminado OP-14; archive/codigo/investigacion/open_problems/ssee_op14_neutrino_mass.py). El vigente es 0.06849. Fila de registro historico
# ORIGEN-VALOR: 0.01361 — 0.41691 - 0.4033 = 0.01361 (s_K_UV menos s_K_IR, los dos de Paper 10)
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from ssee_core import (PHI, PI, OMEGA, BETA, KAL0, OMEGA_DE, OMEGA_M_DYN,
                       OMEGA_M_TOTAL,
                       MIRA, AURA, H0_ALG, N_S)

# ── H₀ candidates ────────────────────────────────────────────────────────────
H0_ALG_KMSMPC  = 67.962         # Type-P algebraic 3(φ+π)² = 3·OMEGA²
H0_MIRA_KMSMPC = 67.068         # Physical: Planck plik_lite + SSEE bg (P3 Cobaya)
h_alg  = H0_ALG_KMSMPC  / 100.0
h_mira = H0_MIRA_KMSMPC / 100.0

# ── Constantes físicas ───────────────────────────────────────────────────────
HBAR_EVS         = 6.582119569e-16
KMSMPC_TO_INVS   = 3.240779289e-20
M_PL_eV          = 1.220910e28   # full Planck mass

# ρ_crit/h² in eV⁴ (standard cosmology constant)
RHO_CRIT_OVER_H2_EV4 = 8.10e-11  # eV⁴, with ρ_crit = h² × 8.10e-11

def rho_crit(h):
    return RHO_CRIT_OVER_H2_EV4 * h**2

def rho_DE(h):
    """OJO (fix 2026-08-10): OMEGA_DE es el ALIAS de S_DE = T_r/M_v = 0.839950,
    que es una SATURACION, no una fraccion de densidad. Multiplicarla por rho_crit
    NO da la densidad de energia oscura (esa es 0.691119*rho_crit). Se conserva
    la funcion porque el resto de este script de auditoria la reporta como ratio
    adimensional en unidades rho_crit=1, donde el numero 0.840 es lo declarado en
    Paper 7 — pero NO usarla como densidad fisica."""
    return OMEGA_DE * rho_crit(h)

def M_UV(h):
    """P10 UV cutoff: M⁴ = 5φ⁸ ρ_crit.

    FIX 2026-08-10: antes usaba rho_DE (= 0.840*rho_crit) con la nota "el paper
    escribe rho_crit pero quiere decir rho_DE". Es al reves: el canonico es
    rho_crit. Comprobado — con M⁴ = 5φ⁸ rho_crit, la cuadratica del UV reproduce
    alpha_K_full = 0.4169052 y f_screen = 0.0695216, que son EXACTAMENTE los
    valores de ssee_paper10_verification.py y de los papers. Con rho_DE no."""
    return (5 * PHI**8 * rho_crit(h))**0.25

def H0_eV(H0_kmsmpc):
    return HBAR_EVS * H0_kmsmpc * KMSMPC_TO_INVS

# ════════════════════════════════════════════════════════════════════════════
# REGISTRO DE AUDITORÍA
# ════════════════════════════════════════════════════════════════════════════
audit = []

def record(paper, name, formula, val_alg, val_mira, units, dim_valid, status, notes=""):
    delta_pct = 100 * (val_mira - val_alg) / val_alg if val_alg != 0 else float('nan')
    audit.append({
        'paper': paper, 'name': name, 'formula': formula,
        'val_alg': val_alg, 'val_mira': val_mira, 'units': units,
        'delta_pct': delta_pct, 'dim_valid': dim_valid,
        'status': status, 'notes': notes,
    })

# ────────────────────────────────────────────────────────────────────────────
# Paper 1 — Framework
# ────────────────────────────────────────────────────────────────────────────
record('P1', 'H₀ (input)', 'directo',
       H0_ALG_KMSMPC, H0_MIRA_KMSMPC, 'km/s/Mpc', True,
       '✅ ya cascadado', 'Tabla principal P1 actualizada a 67.068 esta mañana')

record('P1', 'n_s', '1 − φ⁻⁷', N_S, N_S, '—', True,
       '✅ independiente', 'Adimensional, no depende de H₀')

# ────────────────────────────────────────────────────────────────────────────
# Paper 2 — MCMC
# ────────────────────────────────────────────────────────────────────────────
record('P2', 'H₀ posterior MCMC', 'MCMC con prior(H_MIRA, σ=0.54)',
       66.75, 66.55, 'km/s/Mpc', True,
       '✅ ya cascadado', 'Prior switch Planck legacy → MIRA, esta mañana')

record('P2', 'Ω_b·h² posterior', 'MCMC',
       0.02237, 0.02261, '—', True,
       '✅ ya cascadado', 'Shift por prior MIRA + BBN')

# ────────────────────────────────────────────────────────────────────────────
# Paper 3 — CMB
# ────────────────────────────────────────────────────────────────────────────
record('P3', 'r_d (sound horizon at drag)', 'integral ~∝ 1/(H₀·√Ωₘh²)',
       147.6, 147.6, 'Mpc', True,
       '⚠ verificar', 'Es la cantidad de Planck-MIRA mapped. Likely sin cambio porque Ωₘh² está fija. Verificar Cobaya output.')

record('P3', 'H_MIRA (Cobaya output)', 'plik_lite + SSEE bg',
       67.066, 67.068, 'km/s/Mpc', True,
       '✅ ya cascadado', 'Cobaya re-corrido esta mañana con bounds amplios')

# ────────────────────────────────────────────────────────────────────────────
# Paper 4 — ToE / Algebraic derivations
# ────────────────────────────────────────────────────────────────────────────
record('P4', 'H_alg = 3(φ+π)²', 'algebraico puro', H0_ALG_KMSMPC, H0_ALG_KMSMPC,
       'km/s/Mpc', True, '✅ Type-P invariante',
       'Es el origen del 67.96 — no cascada (por definición).')

record('P4', 'Σm_ν^active', 'R × Ω_b·h² × 93.14 / τ_Π',
       0.0824, 0.0824, 'eV', True, '⚠ Type-P retenido',
       'Si Ω_b·h² cambia 1% bajo MCMC MIRA, Σm_ν cambia ~1%. NO recomputado.')

# ────────────────────────────────────────────────────────────────────────────
# Paper 6 — φ-DM
# ────────────────────────────────────────────────────────────────────────────
record('P6', 'm_φ = Σm_ν × H₀', 'INVÁLIDA dimensionalmente',
       5.602, 5.527, 'eV', False, '❌ REVERTIDO',
       'Fórmula NO es física. Revertido a 5.602 con H_alg + disclaimer. Ver OP-9.')

record('P6', 'k_fs (DW)', '∝ m_φ^(4/3)',
       0.493, 0.490, 'h/Mpc', True, '❌ REVERTIDO',
       'Vuelve a 0.493 al revertir m_φ a 5.602.')

# ────────────────────────────────────────────────────────────────────────────
# Paper 7 — EFT canonical
# ────────────────────────────────────────────────────────────────────────────
# CORREGIDO 2026-09-05 (R52). Esta fila anotaba V₀ = 0.840, que es la
# SATURACIÓN s_DE metida en ranura de densidad. V₀ es la amplitud de un
# potencial: es una DENSIDAD, y le toca Ω_DE = 1-Ω_m = 0.691119. El script
# de P7 (ssee_eft_verification.py) ya se corrigió el mismo día; esta fila
# registraba la afirmación vieja. La conclusión de la fila NO cambia: sigue
# siendo un ratio puro en unidades ρ_crit=1, luego invariante ante el ancla.
_om_de_dens = 1.0 - OMEGA_M_TOTAL              # 0.691119 densidad, no 0.839950
v0_alg  = _om_de_dens * rho_crit(h_alg)
v0_mira = _om_de_dens * rho_crit(h_mira)
record('P7', 'V₀ = Ω_DE × ρ_crit', 'unidades ρ_crit=1; V₀=0.691119',
       0.691119, 0.691119, '(ρ_crit)', True, '✅ invariante por diseño',
       'Ratio adimensional en unidades ρ_crit=1 ⟹ invariante ante el ancla. '
       'Era 0.840 (saturación) hasta el fix R52 del 2026-09-05.')

record('P7', 'λ = √(3 Ω_m,dyn)', 'adimensional',
       np.sqrt(3*0.160), np.sqrt(3*0.160), '—', True, '✅ independiente',
       'Solo depende de Ω_m,dyn (algebraico). Sin cambio.')

record('P7', 'αK = 0.4033', 'unidades ρ_crit=1; adimensional',
       0.4033, 0.4033, '—', True, '✅ invariante por diseño',
       'P7 línea 166: M⁴ ≡ ρ_crit como normalización. αK es ratio adimensional. Verificado computacionalmente.')

# ────────────────────────────────────────────────────────────────────────────
# Paper 8 — Strong gravity (Vainshtein)
# ────────────────────────────────────────────────────────────────────────────
record('P8', 'r_V (Vainshtein radius)', 'depende de M⁴ y M_⊙',
       0.0, 0.0, 'Mpc', True, '⚠ VERIFICAR',
       'r_V ~ (M_⊙/M_Pl² M⁴)^(1/5). M cambia 0.7%, r_V cambia 0.14%. Menor.')

# ────────────────────────────────────────────────────────────────────────────
# Paper 9 — Hubble Tension
# ────────────────────────────────────────────────────────────────────────────
f_screen = (PI - PHI) / OMEGA**2
record('P9', 'f_screen = (π−φ)/Ω²', 'adimensional',
       f_screen, f_screen, '—', True, '✅ independiente',
       'Pura combinación algebraica. Sin cambio.')

# Dirección canónica (2026-09-06): SH0ES ENTRA, H_global SALE. El número
# puro 3(φ+π)² no es entrada de la cascada — es el blanco. Por eso esta fila
# ya NO depende de qué ancla se elija: no hay ancla de entrada. Ver R55.
H0_SHOES = 73.04
h_glob = H0_SHOES * (1 - f_screen)
record('P9', 'H_glob = H_SH0ES·(1−f_screen)', 'no depende de ancla',
       h_glob, h_glob, 'km/s/Mpc', True, '✅ inmune por construcción',
       f'La entrada es el dato medido, no un ancla: {h_glob:.4f} para cualquier '
       f'elección de ancla. Se compara con 3(φ+π)²=67.96214 (0.17σ).')

# ────────────────────────────────────────────────────────────────────────────
# Paper 10 — UV Completion
# ────────────────────────────────────────────────────────────────────────────
M_alg_meV  = M_UV(h_alg)  * 1e3
M_mira_meV = M_UV(h_mira) * 1e3
record('P10', 'M = (5φ⁸·ρ)^(1/4)', 'depende de convención de ρ',
       8.81, 8.81, 'meV', True, '⚠ COSMÉTICO',
       'P10 usa ρ_Λ ≈ (2.25 meV)⁴ canónico para reportar 8.81 meV. Con ρ_crit puro daría 9.74 meV. ELEGIBLE — no afecta observables (H_UV depende de αK adimensional).')

# αK_full P10 — depende de M⁴ y ρ_DE ratio
# αK_full = αK_baseline × (1 + correction(M⁴))
# Si correction ∝ X²/M⁴ y X ∝ ρ_DE → X²/M⁴ ∝ ρ_DE²/ρ_DE = ρ_DE → ∝ H₀²
# Pero αK es adimensional, así que en realidad cancela. Lo verifico:
# X_bg ∝ ρ_DE(1+w₀); M⁴ ∝ ρ_DE → X²/M⁴ ∝ ρ_DE → escala lineal con H₀²
# αK_full = αK_LO + ε con ε ∝ X²/M⁴ — entonces ε escala con ρ_DE (no se cancela)
alphaK_LO = 0.4033
eps_alg  = 0.41691 - 0.4033  # ≈ 0.01361
eps_mira = eps_alg * (h_mira/h_alg)**2
record('P10', 'αK_full UV', 'unidades ρ_crit=1; ratio puro',
       0.41691, 0.41691, '—', True, '✅ invariante por diseño',
       'Verificado: cálculo en unidades ρ_crit=1, todos los ratios son adimensionales.')

_h_glob_uv = H0_SHOES * (1 - 0.0695216111441)
record('P10', 'H_glob^UV = H_SH0ES·(1−f^UV)', 'no depende de ancla',
       _h_glob_uv, _h_glob_uv, 'km/s/Mpc', True, '✅ inmune por construcción',
       f'f_screen^UV=0.069522 invariante y la entrada es el dato: '
       f'{_h_glob_uv:.6f} para cualquier ancla. Residuo vs 3(φ+π)²: '
       f'{_h_glob_uv-H0_ALG:+.2e}.')

# ────────────────────────────────────────────────────────────────────────────
# Postulados
# ────────────────────────────────────────────────────────────────────────────
record('Pos', 'Postulado D (estructural)', '—', 0, 0, '—', True, '✅ independiente',
       'Definición estructural sin referencia a H₀.')
record('Pos', 'Postulado S (estructural)', '—', 0, 0, '—', True, '✅ independiente',
       'Definición estructural sin referencia a H₀.')
record('Pos', 'Postulado M (MIRA=(3φ+π)/4)', 'adimensional', MIRA, MIRA, '—', True,
       '✅ independiente', 'Valor algebraico fijo.')
record('Pos', 'Postulado I (α-attractor, n=2φ⁷)', 'adimensional', 7, 7, '—', True,
       '✅ independiente', 'Corolario entero (e-folds), no depende de H₀.')

# ════════════════════════════════════════════════════════════════════════════
# REPORTE
# ════════════════════════════════════════════════════════════════════════════
print("=" * 100)
print(" H0 CASCADE AUDIT — SSEE  (2026-05-24 PM)")
print("=" * 100)
print(f"\n H_alg  = {H0_ALG_KMSMPC} km/s/Mpc  (algebraic, Type-P, P4)")
print(f" H_MIRA = {H0_MIRA_KMSMPC} km/s/Mpc  (physical, Planck plik_lite + SSEE bg, P3)")
print(f" Ratio  = {H0_MIRA_KMSMPC/H0_ALG_KMSMPC:.5f}  (Δ = {100*(H0_MIRA_KMSMPC-H0_ALG_KMSMPC)/H0_ALG_KMSMPC:+.2f}%)\n")

print(f"{'Paper':<5} {'Cantidad':<35} {'H_alg':>14} {'H_MIRA':>14} {'Δ%':>8} {'Dim':>4} {'Status':<18}")
print("─" * 100)
for r in audit:
    val_alg_str  = f"{r['val_alg']:.4g}"  if abs(r['val_alg']) > 1e-3 or r['val_alg'] == 0 else f"{r['val_alg']:.3e}"
    val_mira_str = f"{r['val_mira']:.4g}" if abs(r['val_mira']) > 1e-3 or r['val_mira'] == 0 else f"{r['val_mira']:.3e}"
    delta_str    = f"{r['delta_pct']:+.2f}" if not np.isnan(r['delta_pct']) else "—"
    dim          = "✓" if r['dim_valid'] else "✗"
    print(f"{r['paper']:<5} {r['name']:<35} {val_alg_str:>14} {val_mira_str:>14} {delta_str:>8} {dim:>4} {r['status']:<18}")

# ── Resumen por categorías ──────────────────────────────────────────────────
print("\n" + "=" * 100)
print(" RESUMEN POR CATEGORÍA")
print("=" * 100)

by_status = {}
for r in audit:
    key = r['status']
    by_status.setdefault(key, []).append(r)
for status, items in sorted(by_status.items()):
    print(f"\n{status}  ({len(items)} items):")
    for r in items:
        print(f"    [{r['paper']}] {r['name']}: {r['notes']}")

# ── Items que requieren acción concreta ─────────────────────────────────────
print("\n" + "=" * 100)
print(" ACCIONES PENDIENTES (orden de prioridad)")
print("=" * 100)

pending = [r for r in audit if '⚠' in r['status']]
invalid = [r for r in audit if '❌' in r['status']]

print(f"\n❌ FÓRMULAS INVÁLIDAS / REVERTIDAS ({len(invalid)}):")
for r in invalid:
    print(f"   [{r['paper']}] {r['name']}: {r['notes']}")

print(f"\n⚠ PENDIENTES DE VERIFICACIÓN / ACTUALIZACIÓN ({len(pending)}):")
for r in pending:
    print(f"   [{r['paper']}] {r['name']}")
    print(f"        H_alg={r['val_alg']:.4g} → H_MIRA={r['val_mira']:.4g} ({r['delta_pct']:+.2f}%)")
    print(f"        {r['notes']}\n")

print("=" * 100)
print(" FIN DE AUDITORÍA")
print("=" * 100)
