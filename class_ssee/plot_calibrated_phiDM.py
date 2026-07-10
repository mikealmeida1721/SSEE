"""
SSEE-V3.6 — Figuras φ-DM FORWARD (Paper 6 / Unified Journal).

Reescrito 2026-06-25 (reframe ω_m-directo + partícula SOLAR²·KRYSTOS):
ELIMINA el fit α_WDM (que estaba calibrado a KiDS, output/alpha_calibrated.npy)
y CONSUME el forward canónico — cero parámetros libres:

  • P(k):  salida CLASS real con la partícula m_φ=40.70 eV
           (output/can_part__pk.dat = dos sectores; output/can_cold__pk.dat = techo frío)
           generada por ../src/ssee_paper6_canonical_particle.py
  • fσ₈:   ODE de crecimiento con fondo PLANO Ω_m=0.30889 (Ω_DE=1−Ω_m=0.69111),
           idéntica a ../src/p06_phiDM/ssee_paper6_verification.py

Amplitud σ₈ del two-sector = 0.747 (CON free-streaming), usada para fσ₈ Y S₈
(CORREGIDO 2026-06-29):
  σ₈ = 0.747 → tanto S₈/lensing como fσ₈/RSD. La supresión de free-streaming
       muerde DENTRO de la ventana σ₈(R=8) que mide RSD (k_half≈0.35 h/Mpc),
       no solo el lensing. El uso previo de σ₈=0.814 (sin supresión) para fσ₈
       suponía "k≪k_fs, φ-DM agrupa como frío" — VERIFICADO FALSO.
  σ₈_single = 0.811·G_2s = 0.814 → baseline single-sector (SIN free-streaming, 0.70σ).
  k_fs = 0.754 h/Mpc (analítico DW, falseable DESI/Euclid)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ── Constantes (φ,π) — idénticas a verification.py ────────────────────────────
phi   = (1 + 5**0.5) / 2
pi_   = np.pi
beta  = (pi_ + phi) / 2
Kv    = phi + pi_ + (pi_ + phi)
Tr    = 3 * (phi + beta)
Mv    = phi + pi_ + Kv
KAL0  = beta + pi_
w0    = -Tr / Mv                       # -0.840
wa    = -(pi_ + phi + phi) / Kv        # -0.670
Om_DNAV = pi_ + phi

# Fondo ω_m-directo (sin factor materia)
H0_GLOBAL = 67.962; h_ = H0_GLOBAL/100.0
n_s     = 1.0 - phi**-7
omega_b = (pi_ - phi) / (3 * Om_DNAV**2)
omega_c = KAL0 * omega_b * n_s
R2      = Om_DNAV / (KAL0 * Tr)
mnu     = R2 * 0.960318
omega_nu= mnu / 93.14
omega_m = omega_b + omega_c + omega_nu
Om_total = omega_m / h_**2             # 0.30889 (Ω_m,CMB, clustering)
Om_CDM   = 1.0 + w0                    # 0.160 (sector dinámico)
m_phi    = 40.70
k_fs     = 0.754
sig8_Planck = 0.811

# Cada fondo PLANO por separado: Ω_DE = 1 − Ω_m (corrige el artefacto fσ₈ 0.82σ)
OmDE_total = 1.0 - Om_total            # 0.69111
OmDE_cdm   = 1.0 - Om_CDM              # 0.84000

# ── Ecuación de estado φ (CPL) normalizada f_DE(a=1)=1 ────────────────────────
def f_DE_raw(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
_norm = f_DE_raw(1.0)
def f_DE(a): return f_DE_raw(a) / _norm

def E2_total(a):
    return max(Om_total*a**(-3) + OmDE_total*f_DE(a), 1e-30)

def growth_rhs_twosector(lna, Y):
    a = np.exp(lna); e2 = E2_total(a)
    h = 0.5*(-3*Om_total*a**(-3) + OmDE_total*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a)/e2
    Omm = Om_total*a**(-3)/e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm*D]

def growth_lcdm(lna, Y, Om_lcdm=0.3153):
    a = np.exp(lna); e2 = max(Om_lcdm*a**(-3) + (1-Om_lcdm), 1e-30)
    h = 0.5*(-3*Om_lcdm*a**(-3))/e2
    Omm = Om_lcdm*a**(-3)/e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm*D]

lna_pts = np.linspace(np.log(0.01), 0.0, 600)
IC = [0.01, 0.01]
sol_2s   = solve_ivp(growth_rhs_twosector, [lna_pts[0],0.0], IC, t_eval=lna_pts, rtol=1e-10, atol=1e-13)
sol_lcdm = solve_ivp(growth_lcdm,          [lna_pts[0],0.0], IC, t_eval=lna_pts, rtol=1e-10, atol=1e-13)

G_2s = sol_2s.y[0][-1] / sol_lcdm.y[0][-1]
sig8_growth = sig8_Planck * G_2s       # amplitud RSD (fσ₈)
sig8_eff    = 0.747                     # free-streaming CLASS (S₈/lensing)
S8_eff      = sig8_eff * np.sqrt(Om_total/0.3)

D_2s  = sol_2s.y[0]/sol_2s.y[0][-1];   Dp_2s  = sol_2s.y[1]/sol_2s.y[0][-1]
D_lc  = sol_lcdm.y[0]/sol_lcdm.y[0][-1]; Dp_lc = sol_lcdm.y[1]/sol_lcdm.y[0][-1]

# ── P(k) forward desde CLASS (cero fit) ───────────────────────────────────────
def load(f):
    d = np.loadtxt(f, comments="#"); return d[:,0], d[:,1]
k_p, P_p = load("output/can_part__pk.dat")   # dos sectores (partícula m_φ=40.70)
k_c, P_c = load("output/can_cold__pk.dat")   # frío (techo σ₈=0.8335)
k_l, P_l = load("output/lcdm_planck2018__pk.dat")

def W_tophat(k, R=8.0):
    x = k*R; return 3.0*(np.sin(x) - x*np.cos(x))/x**3
def sigma8_real(k, P):
    W = W_tophat(k); return np.sqrt(np.trapezoid(k**2*P*W**2/(2*np.pi**2), k))
s8_part = sigma8_real(k_p, P_p)
s8_cold = sigma8_real(k_c, P_c)

# k_half: donde el cociente partícula/frío cae a 0.5 (borrado medido por Boltzmann)
kg = np.logspace(np.log10(0.05), np.log10(2.0), 400)
ratio = np.clip(np.interp(kg, k_p, P_p)/np.interp(kg, k_c, P_c), 1e-6, 1.2)
k_half = float(np.interp(0.5, ratio[::-1], kg[::-1])) if ratio.min() < 0.5 else np.nan

print("="*64)
print("  SSEE φ-DM FORWARD — figuras Paper 6 / Unified")
print("="*64)
print(f"  m_φ        = {m_phi} eV")
print(f"  Ω_m,CMB    = {Om_total:.5f}   Ω_φDM = {Om_total-Om_CDM:.5f}")
print(f"  G_2s       = {G_2s:.4f}")
print(f"  σ₈_single  = {sig8_growth:.4f}  (single-sector baseline, SIN free-streaming)")
print(f"  σ₈_eff     = {sig8_eff:.4f}  (two-sector, CON free-streaming → S₈ Y fσ₈)")
print(f"  S₈_eff     = {S8_eff:.4f}")
print(f"  σ₈(part)   = {s8_part:.4f}   σ₈(frío,techo) = {s8_cold:.4f}  [CLASS top-hat]")
print(f"  k_fs       = {k_fs} h/Mpc   k_half(CLASS) = {k_half:.4f} h/Mpc")

# ── Datos fσ₈ (set canónico Paper 5) ──────────────────────────────────────────
Z_RSD   = np.array([0.067, 0.150, 0.380, 0.510, 0.610, 1.480])
surveys = ['6dFGRS','SDSS MGS','BOSS DR12','BOSS DR12','BOSS DR12','eBOSS DR16']
fo      = np.array([0.423, 0.490, 0.497, 0.458, 0.436, 0.462])
fe      = np.array([0.055, 0.145, 0.045, 0.038, 0.034, 0.045])

def get_fsig8(z, D, Dp, s8):
    lna = np.log(1/(1+z))
    f_t = float(np.interp(lna, lna_pts, Dp/D))
    D_t = float(np.interp(lna, lna_pts, D))
    return f_t * s8 * D_t

# two-sector usa σ₈=0.747 (CON free-streaming); single-sector usa 0.814 (baseline)
fs_2s     = np.array([get_fsig8(z, D_2s, Dp_2s, sig8_eff)    for z in Z_RSD])
fs_single = np.array([get_fsig8(z, D_2s, Dp_2s, sig8_growth) for z in Z_RSD])
fs_lcdm   = np.array([get_fsig8(z, D_lc, Dp_lc, sig8_Planck) for z in Z_RSD])
t_2s      = (fo - fs_2s)/fe
t_single  = (fo - fs_single)/fe
print(f"\n  fσ₈ tensión media two-sector (σ₈=0.747) = {np.mean(np.abs(t_2s)):.4f}σ  (<1σ)")
print(f"  fσ₈ tensión media single-sector (σ₈={sig8_growth:.4f}) = {np.mean(np.abs(t_single)):.4f}σ  (baseline)")

# ═══════════════════════════════════════════════════════════════════════════════
# Figura 1 — P(k) forward
# ═══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 1, figsize=(8, 10), gridspec_kw={'height_ratios':[2,1]})
ax1 = axes[0]
ax1.loglog(k_l, P_l, 'k-',  lw=2,   label=r'$\Lambda$CDM (Planck 2018)', zorder=3)
ax1.loglog(k_c, P_c, 'b--', lw=1.6, label=rf'SSEE frío (techo, $\sigma_8={s8_cold:.3f}$)', alpha=0.85)
ax1.loglog(k_p, P_p, 'r-',  lw=2.2, label=rf'SSEE $\varphi$-DM forward ($m_\varphi=41.0$ eV, $\sigma_8={s8_part:.3f}$)', zorder=4)
ax1.axvline(k_fs, color='orange', ls='--', lw=1.5, label=rf'$k_{{fs}}={k_fs}$ h/Mpc (DW, falseable)')
if np.isfinite(k_half):
    ax1.axvline(k_half, color='green', ls=':', lw=1.5, label=rf'$k_{{1/2}}={k_half:.3f}$ h/Mpc (CLASS)')
ax1.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
ax1.set_ylabel(r'$P(k)$ [$(h^{-1}$ Mpc$)^3$]', fontsize=13)
ax1.set_title(r'SSEE-V3.6: $P(k)$ forward $\varphi$-DM — salida CLASS, cero fit', fontsize=13)
ax1.legend(fontsize=9.5, loc='lower left')
ax1.set_xlim(1e-3, 2); ax1.grid(True, alpha=0.3)

# Subpanel: cociente partícula/frío (transferencia REAL del Boltzmann)
ax2 = axes[1]
ax2.semilogx(kg, ratio, 'r-', lw=2, label=r'$P_{\varphi\mathrm{DM}}/P_\mathrm{frío}$ (CLASS, forward)')
ax2.axvline(k_fs, color='orange', ls='--', lw=1.5)
if np.isfinite(k_half):
    ax2.axvline(k_half, color='green', ls=':', lw=1.5)
ax2.axhline(0.5, color='gray', ls=':', lw=1)
ax2.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
ax2.set_ylabel(r'$P_{\varphi\mathrm{DM}}/P_\mathrm{frío}$', fontsize=12)
ax2.set_title(r'Supresión free-streaming medida por CLASS (no impuesta)', fontsize=12)
ax2.legend(fontsize=9.5); ax2.set_xlim(1e-2, 3); ax2.set_ylim(-0.05, 1.15); ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("output/ssee_phiDM_Pk_calibrated.pdf", dpi=150, bbox_inches='tight')
plt.savefig("output/ssee_phiDM_Pk_calibrated.png", dpi=150, bbox_inches='tight')
print("\n  Figura guardada: output/ssee_phiDM_Pk_calibrated.{pdf,png}")

# ═══════════════════════════════════════════════════════════════════════════════
# Figura 2 — fσ₈ forward
# ═══════════════════════════════════════════════════════════════════════════════
fig2, ax = plt.subplots(figsize=(8, 5))
z_th = np.linspace(0.05, 1.6, 100)
fs_lcdm_th = np.array([get_fsig8(z, D_lc, Dp_lc, sig8_Planck) for z in z_th])
fs_2s_th   = np.array([get_fsig8(z, D_2s, Dp_2s, sig8_eff)    for z in z_th])
fs_sgl_th  = np.array([get_fsig8(z, D_2s, Dp_2s, sig8_growth) for z in z_th])
ax.plot(z_th, fs_lcdm_th, 'k-',  lw=2,   label=rf'$\Lambda$CDM ($\sigma_8={sig8_Planck:.3f}$, $\bar\sigma={np.mean(np.abs((fo-fs_lcdm)/fe)):.2f}$)')
ax.plot(z_th, fs_sgl_th,  'b--', lw=1.6, label=rf'SSEE single-sector ($\sigma_8={sig8_growth:.3f}$, $\bar\sigma={np.mean(np.abs(t_single)):.2f}$)')
ax.plot(z_th, fs_2s_th,   'r-',  lw=2.2, label=rf'SSEE two-sector forward ($\sigma_8={sig8_eff:.3f}$, $\bar\sigma={np.mean(np.abs(t_2s)):.2f}$)')
ax.errorbar(Z_RSD, fo, yerr=fe, fmt='ko', capsize=4, ms=6, label='Datos RSD (6 encuestas)')
for i, z in enumerate(Z_RSD):
    ax.annotate(f'{t_2s[i]:+.1f}$\\sigma$', (z, fo[i]), textcoords="offset points",
                xytext=(5, 8), fontsize=8, color='red')
ax.set_xlabel(r'Redshift $z$', fontsize=13)
ax.set_ylabel(r'$f\sigma_8(z)$', fontsize=13)
ax.set_title(rf'SSEE-V3.6: $f\sigma_8$ forward two-sector — tensión media {np.mean(np.abs(t_2s)):.2f}$\sigma$ ($<1\sigma$; single-sector {np.mean(np.abs(t_single)):.2f}$\sigma$, $\Lambda$CDM 0.73$\sigma$)', fontsize=11)
ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1.6); ax.set_ylim(0.25, 0.65)
plt.tight_layout()
plt.savefig("output/ssee_fsig8_tension_phiDM.pdf", dpi=150, bbox_inches='tight')
plt.savefig("output/ssee_fsig8_tension_phiDM.png", dpi=150, bbox_inches='tight')
print("  Figura guardada: output/ssee_fsig8_tension_phiDM.{pdf,png}")
print("="*64)
