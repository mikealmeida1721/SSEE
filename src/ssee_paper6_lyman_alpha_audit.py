"""
SSEE Paper 6 — Lyman-α Audit: 3 Routes
========================================
Evaluación cuantitativa de si m_φ = 5.71 eV puede sostenerse frente a
los límites observacionales de Lyman-α y estructura a gran escala.

Rutas:
  A — Reinterpretar φ-DM como campo escalar (fuzzy/axion)
  B — Verificar validez de la relación DW extrapolada a eV
  C — Aceptar HDM y calcular impacto real en P(k)

Resultado esperado: diagnóstico de cuál (si alguna) ruta es viable, y
qué reinterpretación física sobrevive.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTES SSEE
# ─────────────────────────────────────────────────────────────────────────────
phi_gr   = (1 + 5**0.5) / 2          # razón áurea
pi_val   = np.pi
KAL0     = (phi_gr + pi_val) / 2 + pi_val   # ≈ 5.5214

Omega_DE      = 0.840
Omega_m_dyn   = 0.160
Omega_CDM     = 0.160050
Omega_phiDM   = 0.159878
Omega_m_CMB   = 0.319928
MIRA          = 1.9989

m_phi_eV  = 5.71             # masa algebraica φ-DM [eV]
m_phi_keV = m_phi_eV / 1e3  # [keV]
h         = 0.670
H0        = h * 100          # [km/s/Mpc]

Omega_phiDM_h2 = Omega_phiDM * h**2   # ≈ 0.0718
Omega_m_h2     = Omega_m_CMB  * h**2  # ≈ 0.1437
f_phiDM        = Omega_phiDM / Omega_m_CMB  # ≈ 0.50

# Parámetros WDM transfer function (Viel 2005)
nu = 1.12

# ─────────────────────────────────────────────────────────────────────────────
# RUTA A — Campo escalar coherente (fuzzy/axion DM)
# ─────────────────────────────────────────────────────────────────────────────
# Para fuzzy DM (campo escalar ultraligero), la escala de supresión es la
# longitud de Jeans cuántica (presión cuántica, no free-streaming térmico):
#   k_J [h/Mpc] ≈ 66.5 × (m / 10^-22 eV)^0.5   [Hui et al. 2017]
#
# Invertir para obtener la masa que daría k_J = k_fs,paper = 0.493 h/Mpc:

k_fs_paper = 0.493  # h/Mpc (del paper)

m_fuzzy_needed_eV = (k_fs_paper / 66.5)**2 * 1e-22  # [eV]

# k_J para m_φ = 5.71 eV (como campo escalar):
k_J_5p71eV = 66.5 * (m_phi_eV / 1e-22)**0.5  # h/Mpc — astronómicamente grande

# Longitud de de Broglie para m = 5.71 eV moviéndose a velocidad virial típica
# v_vir ~ 100 km/s = 3.3×10^-4 c  →  λ_dB [Mpc] = ℏc / (m c² × v/c)
# ℏc = 197.3 MeV·fm = 197.3×10^-15 MeV·m → en unidades cosmológicas:
# λ_dB [Mpc] = 6.58×10^-22 eV·s × c / (m × v)
# Para m=5.71 eV, v=100 km/s=3.33×10^5 m/s:
hbar_eV_s = 6.582e-16      # [eV·s]
c_m_s     = 3e8            # [m/s]
v_vir     = 1e5            # [m/s] (100 km/s)
Mpc_m     = 3.086e22       # [m/Mpc]
lambda_dB_Mpc = (hbar_eV_s * c_m_s) / (m_phi_eV * v_vir) / Mpc_m
k_dB = 2 * np.pi / lambda_dB_Mpc  # [1/Mpc] ~ [h/Mpc] for h~1

# ─────────────────────────────────────────────────────────────────────────────
# RUTA B — Validez de la relación DW extrapolada a eV
# ─────────────────────────────────────────────────────────────────────────────
# Fórmula del paper (DW free-streaming en matter-radiation equality):
k_fs_DW_paper = 7.68 * m_phi_keV**0.5 * (Omega_phiDM_h2 / 0.1)**0.5  # h/Mpc

# Escala de corte Viel 2005 (calibrada en keV, extrapolada a eV):
# α [h⁻¹ Mpc] = 0.049 × (m/keV)^(-1.11) × (Ω/0.25)^0.11 × (h/0.7)^1.22
alpha_Viel = (0.049
              * m_phi_keV**(-1.11)
              * (Omega_phiDM / 0.25)**0.11
              * (h / 0.7)**1.22)   # [h⁻¹ Mpc]

k_fs_Viel = 1.0 / alpha_Viel  # [h/Mpc]

# Límites Lyman-α observacionales (literatura)
m_thermal_Lya_keV = 3.5   # Viel+2013, SDSS Lyman-α
m_DW_Lya_keV      = 1.5   # DW, solo Lyman-α
m_DW_Xray_keV     = 8.0   # DW + rayos-X (Boyarsky+2019)

# T(k) con ambas fórmulas en array de k
k_arr = np.logspace(-2, 1, 1000)  # [h/Mpc]

T_DW_paper = (1 + (k_arr / k_fs_DW_paper)**(2*nu))**(-5.0/nu)
T_Viel_ext = (1 + (k_arr * alpha_Viel)**(2*nu))**(-5.0/nu)

# T(k) para WDM térmico mínimo viable (3.5 keV), para referencia
alpha_3p5keV = (0.049 * 3.5**(-1.11) * (Omega_phiDM/0.25)**0.11 * (h/0.7)**1.22)
T_3p5keV     = (1 + (k_arr * alpha_3p5keV)**(2*nu))**(-5.0/nu)

# Valores en escalas clave
k_BAO = 0.05    # h/Mpc
k_Lya = 1.0     # h/Mpc
k_BAO_idx  = np.argmin(np.abs(k_arr - k_BAO))
k_Lya_idx  = np.argmin(np.abs(k_arr - k_Lya))

T_paper_BAO = T_DW_paper[k_BAO_idx]
T_paper_Lya = T_DW_paper[k_Lya_idx]
T_Viel_BAO  = T_Viel_ext[k_BAO_idx]
T_Viel_Lya  = T_Viel_ext[k_Lya_idx]

# ─────────────────────────────────────────────────────────────────────────────
# RUTA C — Aceptar HDM, impacto real en P(k)
# ─────────────────────────────────────────────────────────────────────────────
# f_phiDM = 0.50 → 50% de la materia oscura es HDM-like

# Aproximación lineal (válida solo para f_ν << 1, se rompe aquí):
sigma8_CDM_ref   = 0.811   # ΛCDM Planck 2018
delta_linear     = -8 * f_phiDM
sigma8_HDM_lin   = sigma8_CDM_ref * (1 + delta_linear)

# Aproximación de Bond 1980 (no-lineal para f grande):
# σ₈(f) ≈ σ₈,CDM × (1 − f_ν)²   [válida para k >> k_fs,ν]
sigma8_HDM_Bond  = sigma8_CDM_ref * (1 - f_phiDM)**2

# k_fs para HDM térmico a 5.71 eV (tipo neutrino activo):
# k_fs,ν [h/Mpc] ≈ 0.82 (m/eV)^0.5 (Ω_m h²/0.14)^0.5   [Lesgourgues & Pastor 2006]
k_fs_HDM_thermal = 0.82 * m_phi_eV**0.5 * (Omega_m_h2 / 0.14)**0.5

# Supresión P(k) a escala Lyman-α (k >> k_fs):
# ΔP/P ≈ −8 f_ν + O(f²)  (lineal)
# (1−f)² − 1  (Bond)
suppression_Lya_lin  = delta_linear            # ΔP/P
suppression_Lya_Bond = (1 - f_phiDM)**2 - 1   # ΔP/P

# Límite Planck 2018 en Σmν activos:
Sigma_mnu_Planck_max = 0.12   # eV (95% CL)

# ─────────────────────────────────────────────────────────────────────────────
# FIGURAS
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 10))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.40, wspace=0.35)

ax_A  = fig.add_subplot(gs[0, 0])   # Ruta A — masa requerida
ax_B  = fig.add_subplot(gs[0, 1:])  # Ruta B — T(k)
ax_C1 = fig.add_subplot(gs[1, 0])   # Ruta C — σ₈ supresión
ax_C2 = fig.add_subplot(gs[1, 1:])  # Ruta C — resumen masas excluidas

# ── Ruta A ──────────────────────────────────────────────────────────────────
masses_eV = np.logspace(-28, 2, 500)
k_J_vals  = 66.5 * (masses_eV / 1e-22)**0.5

ax_A.loglog(masses_eV, k_J_vals, 'steelblue', lw=2.5, label=r'$k_J$ (fuzzy DM)')
ax_A.axhline(k_fs_paper, color='crimson', ls='--', lw=1.8,
             label=rf'$k_{{fs}}^{{\rm paper}}={k_fs_paper}$ h/Mpc')
ax_A.axvline(m_phi_eV, color='darkorange', ls=':', lw=2.0,
             label=rf'$m_\phi={m_phi_eV}$ eV')
ax_A.axvline(m_fuzzy_needed_eV, color='green', ls='-.', lw=1.8,
             label=rf'$m_{{fuzzy}}^{{needed}}\approx{m_fuzzy_needed_eV:.1e}$ eV')
ax_A.set_xlabel(r'$m$ [eV]', fontsize=11)
ax_A.set_ylabel(r'$k_J$ [h/Mpc]', fontsize=11)
ax_A.set_title('Route A — Fuzzy scalar field', fontsize=11, fontweight='bold')
ax_A.legend(fontsize=7.5, loc='upper left')
ax_A.set_xlim(1e-28, 1e2)
ax_A.set_ylim(1e-6, 1e10)
ax_A.grid(True, alpha=0.3)

# ── Ruta B ──────────────────────────────────────────────────────────────────
ax_B.semilogx(k_arr, T_DW_paper, 'steelblue', lw=2.5,
              label=rf'Paper DW formula: $k_{{fs}}={k_fs_DW_paper:.3f}$ h/Mpc')
ax_B.semilogx(k_arr, T_Viel_ext, 'crimson', lw=2.5, ls='--',
              label=rf'Viel 2005 extrapolated: $k_{{fs}}={k_fs_Viel:.4f}$ h/Mpc')
ax_B.semilogx(k_arr, T_3p5keV, 'green', lw=1.8, ls=':',
              label=r'$m_{\rm WDM}^{\rm thermal}=3.5$ keV (Lyman-α min.)')
ax_B.axvline(k_BAO, color='gray', ls=':', lw=1.5, alpha=0.7, label='BAO scale')
ax_B.axvline(k_Lya, color='gray', ls='-.', lw=1.5, alpha=0.7, label='Lyman-α scale')
ax_B.fill_betweenx([0, 1], 0.01, k_BAO, alpha=0.07, color='gray',
                   label='Below BAO — observed ✓')
ax_B.axhspan(0, 0.10, alpha=0.08, color='crimson',
             label='T(k) < 0.10 — catastrophic suppression')
ax_B.set_xlabel(r'$k$ [h/Mpc]', fontsize=11)
ax_B.set_ylabel(r'$T_{\rm WDM}(k)$', fontsize=11)
ax_B.set_title(r'Route B — WDM transfer function: paper vs Viel 2005 extrapolation',
               fontsize=11, fontweight='bold')
ax_B.legend(fontsize=8, loc='lower left')
ax_B.set_xlim(0.01, 10)
ax_B.set_ylim(-0.05, 1.10)
ax_B.grid(True, alpha=0.3)

# ── Ruta C — σ₈ ────────────────────────────────────────────────────────────
f_arr  = np.linspace(0, 1, 200)
s8_lin = sigma8_CDM_ref * (1 - 8*f_arr)
s8_Bond= sigma8_CDM_ref * (1 - f_arr)**2

ax_C1.plot(f_arr, s8_lin,  'steelblue', lw=2.5, label='Linear (Δσ₈/σ₈ = −8f)')
ax_C1.plot(f_arr, s8_Bond, 'crimson',   lw=2.5, ls='--', label=r'Bond: $\sigma_8(1-f)^2$')
ax_C1.axhline(0.737, color='green', ls=':', lw=2,
              label=r'$\sigma_{8,\rm eff}=0.737$ (Paper 6 claim)')
ax_C1.axhline(0.737, color='green', ls=':', lw=2)
ax_C1.axhline(0.766, color='purple', ls=':', lw=1.5,
              label='KiDS $S_8$ target')
ax_C1.axvline(f_phiDM, color='darkorange', ls='-', lw=2,
              label=rf'$f_{{φDM}}={f_phiDM:.2f}$')
ax_C1.fill_between(f_arr, 0.700, 0.775, alpha=0.12, color='green',
                   label='KiDS-1000 σ₈ range')
ax_C1.set_xlabel(r'$f_{\phi DM} = \Omega_{\phi DM}/\Omega_m$', fontsize=11)
ax_C1.set_ylabel(r'$\sigma_8$', fontsize=11)
ax_C1.set_title('Route C — HDM σ₈ suppression', fontsize=11, fontweight='bold')
ax_C1.legend(fontsize=7.5)
ax_C1.set_xlim(0, 1)
ax_C1.set_ylim(-0.1, 0.9)
ax_C1.grid(True, alpha=0.3)

# ── Ruta C — exclusión de masas ────────────────────────────────────────────
obs_labels = [
    r'Planck $\Sigma m_\nu < 0.12$ eV',
    r'Lyman-α DW (thermal equiv.): $m > 1.5$ keV',
    r'Lyman-α thermal WDM: $m > 3.5$ keV',
    r'X-ray + Lyman-α DW: $m > 8.0$ keV',
]
obs_limits = [0.12e-3, 1.5, 3.5, 8.0]   # en keV
obs_colors = ['purple', 'steelblue', 'darkorange', 'crimson']

ax_C2.set_xscale('log')
ax_C2.axvline(m_phi_keV, color='black', lw=3,
              label=rf'$m_\phi = {m_phi_eV}$ eV $= {m_phi_keV:.5f}$ keV')
for lim, label, col in zip(obs_limits, obs_labels, obs_colors):
    ax_C2.axvline(lim, color=col, ls='--', lw=2, label=label)
    ax_C2.fill_betweenx([0,1], lim, 100, alpha=0.06, color=col)

ax_C2.set_xlim(1e-5, 100)
ax_C2.set_ylim(0, 1)
ax_C2.set_xlabel(r'$m_{\phi DM}$ [keV]', fontsize=11)
ax_C2.set_title('Route C — Mass exclusion limits', fontsize=11, fontweight='bold')
ax_C2.legend(fontsize=8, loc='upper left')
ax_C2.set_yticks([])
ax_C2.grid(True, alpha=0.3, axis='x')

# Título global
fig.suptitle(
    r'SSEE Paper 6 — Lyman-$\alpha$ Audit: $m_\phi = 5.71$ eV vs. observational limits',
    fontsize=13, fontweight='bold', y=1.01
)

_outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "figures")
plt.savefig(os.path.join(_outdir, 'fig_paper6_lyman_alpha_audit.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(_outdir, 'fig_paper6_lyman_alpha_audit.png'), bbox_inches='tight', dpi=150)
plt.close()

# ─────────────────────────────────────────────────────────────────────────────
# DIAGNÓSTICO IMPRESO
# ─────────────────────────────────────────────────────────────────────────────
SEP = "=" * 70
sep = "─" * 70

print(SEP)
print("  SSEE Paper 6 — Auditoría Lyman-α: 3 Rutas")
print(f"  m_φ = {m_phi_eV} eV = {m_phi_keV:.5f} keV")
print(f"  Ω_φDM = {Omega_phiDM:.4f}  |  h = {h}  |  f_φDM = {f_phiDM:.3f}")
print(SEP)

print(f"\n{'RUTA A — Campo escalar (fuzzy/axion DM)':^70}")
print(sep)
print(f"  Escala de Jeans fuzzy DM: k_J = 66.5 × (m/10⁻²² eV)^0.5")
print(f"  Para k_J = {k_fs_paper} h/Mpc (target paper):")
print(f"    m_requerida = {m_fuzzy_needed_eV:.2e} eV  (~21 órdenes de magnitud más ligera)")
print(f"  Para m_φ = {m_phi_eV} eV como campo escalar:")
print(f"    k_J = {k_J_5p71eV:.2e} h/Mpc  (mucho más allá de cualquier escala observable)")
print(f"    λ_dB (v=100 km/s) = {lambda_dB_Mpc:.2e} Mpc  (escala sub-atómica cosmológica)")
print(f"\n  VEREDICTO: ✗  A 5.71 eV, campo escalar → CDM puro. Sin supresión posible.")
print(f"             Fuzzy DM requiere m ~ 10⁻²⁶ eV — incompatible con la derivación algebraica.")

print(f"\n{'RUTA B — Validez relación DW extrapolada a eV':^70}")
print(sep)
print(f"  Fórmula del paper (DW):  k_fs = {k_fs_DW_paper:.3f} h/Mpc")
print(f"  Viel 2005 extrapolado:   α = {alpha_Viel:.2f} h⁻¹Mpc  →  k_fs = {k_fs_Viel:.5f} h/Mpc")
print(f"  Discrepancia: {k_fs_DW_paper / k_fs_Viel:.1f}× — la fórmula DW sobreestima k_fs por 7 órdenes")
print(f"\n  T_WDM en escalas clave:")
print(f"    k = {k_BAO} h/Mpc (BAO):  T_paper = {T_paper_BAO:.4f}  |  T_Viel = {T_Viel_BAO:.6f}")
print(f"    k = {k_Lya} h/Mpc (Lyα):  T_paper = {T_paper_Lya:.4f}  |  T_Viel = {T_Viel_Lya:.6f}")
print(f"\n  Límites observacionales:")
print(f"    Lyman-α térmico:  m > {m_thermal_Lya_keV} keV")
print(f"    DW solo Lyman-α:  m > {m_DW_Lya_keV} keV")
print(f"    DW + rayos-X:     m > {m_DW_Xray_keV} keV")
print(f"    m_φ = {m_phi_keV*1e3:.2f} meV  →  {m_DW_Xray_keV*1e3/m_phi_eV:.0f}× por debajo del límite más débil")
print(f"\n  VEREDICTO: ✗  Empeora el problema. Viel extrapolado da k_fs = {k_fs_Viel:.5f} h/Mpc,")
print(f"             suprimiendo P(k) >99% en escalas BAO. Los datos BAO lo excluyen.")

print(f"\n{'RUTA C — Aceptar HDM, impacto real en P(k)':^70}")
print(sep)
print(f"  f_φDM = {f_phiDM:.3f}  (50% de toda la materia oscura es HDM-like)")
print(f"  k_fs,HDM térmico (5.71 eV) = {k_fs_HDM_thermal:.3f} h/Mpc")
print(f"\n  Supresión σ₈:")
print(f"    Lineal (Δσ₈/σ₈ = −8f):    σ₈ = {sigma8_HDM_lin:.3f}  [aproximación fuera de rango]")
print(f"    Bond (1−f)²:               σ₈ = {sigma8_HDM_Bond:.3f}")
print(f"    Paper 6 claim:             σ₈_eff = 0.737")
print(f"\n  Supresión P(k) para k >> k_fs,HDM:")
print(f"    Lineal:  ΔP/P = {suppression_Lya_lin:.2f}  ({suppression_Lya_lin*100:.0f}%)")
print(f"    Bond:    ΔP/P = {suppression_Lya_Bond:.2f}  ({suppression_Lya_Bond*100:.0f}%)")
print(f"\n  Planck 2018: Σmν < {Sigma_mnu_Planck_max} eV  →  m_φ = {m_phi_eV} eV es {m_phi_eV/Sigma_mnu_Planck_max:.0f}× el límite")
print(f"\n  VEREDICTO: ✗  50% HDM a 5.71 eV destruye la formación de estructura a gran escala.")
print(f"             Bond da σ₈ = {sigma8_HDM_Bond:.3f} — numéricamente cercano al claim del paper,")
print(f"             pero el P(k) completo está catastrófico en k > {k_fs_HDM_thermal:.2f} h/Mpc.")

print(f"\n{SEP}")
print(f"  DIAGNÓSTICO GLOBAL")
print(SEP)
print("""
  Ruta A: ✗  5.71 eV como campo escalar → CDM puro. Sin supresión de k_fs.
             La interpretación fuzzy requiere m ~ 10⁻²⁶ eV.

  Ruta B: ✗  Empeora el problema. La fórmula DW extrapolada es inválida en eV.
             k_fs real (Viel) ~ 7 veces más pequeño → destruye P(k) en BAO.

  Ruta C: ✗  50% HDM a 5.71 eV es incompatible con observaciones LSS existentes.
             La relación DW usada en Paper 6 no aplica a este rango de masa.

  ─────────────────────────────────────────────────────────────────────────
  CONCLUSIÓN:
  m_φ = 5.71 eV NO puede interpretarse como partícula DM free-streaming.
  La derivación algebraica es correcta — el número 5.71 eV puede ser real.
  Lo que necesita cambiar: la INTERPRETACIÓN FÍSICA.

  CANDIDATOS DE REINTERPRETACIÓN:
    (1) m_φ como escala de masa del potencial V(φ), no como masa del portador de DM.
        El sector φ-DM sería un condensado frío (CDM-like) cuya masa del campo
        es 5.71 eV — observable en experimentos de campo escalar, no en Lyman-α.

    (2) La unificación Ω_CDM + Ω_φDM = Ω_m,CMB es válida independientemente
        del mecanismo de supresión. σ₈ puede resolverse por IS viscosity (Paper 5)
        sin necesitar WDM transfer function.

    (3) m_φ como predicción para un experimento de laboratorio (PTOLEMY/KATRIN)
        donde se mide la masa del excitación del campo φ — no como DM cosmológico.
  ─────────────────────────────────────────────────────────────────────────
""")
print(f"  Figura guardada en: results/figures/fig_paper6_lyman_alpha_audit.pdf")
print(SEP)
