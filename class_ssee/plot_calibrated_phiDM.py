"""
SSEE-V3.6 — Fase 2c/3: Figura P(k) calibrado + análisis fσ₈ con CLASS
Alpha calibrado: 1.6561 h/Mpc → sigma_8_eff = 0.737 (top-hat R=8 Mpc/h)
Nota Fase 3: alpha anterior (0.1597) usaba sigma8_proxy sin filtro top-hat.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ── Constantes físicas SSEE ──────────────────────────────────────────────────
phi = (1 + 5**0.5) / 2
pi  = np.pi
Omega_m_dyn  = 0.160
Omega_phiDM  = 0.159878
f_phiDM      = Omega_phiDM / (Omega_m_dyn + Omega_phiDM)  # ≈ 0.5
k_fs         = 0.493   # h/Mpc — Paper 6
alpha_cal    = np.load("output/alpha_calibrated.npy").item()  # recalibrado Fase 3
mu           = 1.12

# ── Cargar P(k) de CLASS ─────────────────────────────────────────────────────
ssee_pk   = np.loadtxt("output/ssee_v36__pk.dat",           comments="#")
lcdm_pk   = np.loadtxt("output/lcdm_planck2018__pk.dat",    comments="#")
two_pk    = np.loadtxt("output/ssee_v36_twosector__pk.dat", comments="#")

k_s, P_s  = ssee_pk[:,0], ssee_pk[:,1]
k_l, P_l  = lcdm_pk[:,0], lcdm_pk[:,1]
k_2, P_2  = two_pk[:,0],  two_pk[:,1]

# pk_cb = CDM+baryon (excluye NCDM) — esto es el sigma_8 medido por lensing
two_pk_cb = np.loadtxt("output/ssee_v36_twosector__pk_cb.dat", comments="#")
k_2cb, P_2cb = two_pk_cb[:,0], two_pk_cb[:,1]

# Interpolar ΛCDM a grilla SSEE
P_l_i = np.interp(k_s, k_l, P_l)

# ── Transferencia WDM calibrada ──────────────────────────────────────────────
def T_wdm(k, alpha):
    return (1 + (alpha * k)**(2*mu))**(-5/mu)

T2_cal = T_wdm(k_s, alpha_cal)
P_2sector_cal = P_s * ((1 - f_phiDM) + f_phiDM * T2_cal)**2

# ── Sigma_8 estándar con filtro top-hat R=8 Mpc/h (Peebles 1980) ─────────────
def W_tophat(k, R=8.0):
    x = k * R
    return 3.0 * (np.sin(x) - x * np.cos(x)) / x**3

def sigma8_real(k, P):
    W = W_tophat(k)
    return np.sqrt(np.trapezoid(k**2 * P * W**2 / (2 * np.pi**2), k))

sig_lcdm  = sigma8_real(k_s, P_l_i)
sig_ssee  = sigma8_real(k_s, P_s)
sig_cal   = sigma8_real(k_s, P_2sector_cal)
sig_two   = sigma8_real(k_s, np.interp(k_s, k_2,   P_2))
sig_twocb = sigma8_real(k_s, np.interp(k_s, k_2cb, P_2cb))

sigma8_planck = 0.811
sigma8_ssee5  = 0.702   # Paper 5 baseline (single-sector ODE)
sigma8_eff    = sig_cal  # top-hat directo desde CLASS P(k) calibrado

print("=" * 60)
print("SSEE-V3.6 — Fase 2c: Resultados calibración φ-DM CLASS")
print("=" * 60)
print(f"alpha_WDM calibrado : {alpha_cal:.6f} h/Mpc")
print(f"f_phiDM             : {f_phiDM:.4f}")
print(f"T_WDM(k_fs)         : {T_wdm(k_fs, alpha_cal):.4f}  ({(1-T_wdm(k_fs,alpha_cal)**2)*100:.1f}% supresión)")
print()
print(f"sigma_8 (top-hat R=8 Mpc/h):")
print(f"  ΛCDM                    : {sig_lcdm:.4f}  ratio=1.0000  (ref)")
print(f"  SSEE+MIRA               : {sig_ssee:.4f}  ratio={sig_ssee/sig_lcdm:.4f}")
print(f"  SSEE+MIRA+φDM cal       : {sig_cal:.4f}  ratio={sig_cal/sig_lcdm:.4f}  (target=0.9088)")
print(f"  CLASS two-sector (total): {sig_two:.4f}  ratio={sig_two/sig_lcdm:.4f}")
print(f"  CLASS two-sector (cb)   : {sig_twocb:.4f}  ratio={sig_twocb/sig_lcdm:.4f}")
print()
print(f"sigma_8_eff (CLASS top-hat directo): {sigma8_eff:.4f}")
print(f"sigma_8 phi-DM (Paper 6 analítico) : 0.737")

# ── Análisis fσ₈ — ODE idéntica a ssee_paper6_verification.py ───────────────
from scipy.integrate import solve_ivp

# Parámetros SSEE
w0, wa    = -0.840, -0.670
Om_CDM    = 0.160050    # sector CDM activo (Paper 6)
Om_phiDM  = 0.159878    # sector φ-DM (solo k < k_fs)
Om_total  = Om_CDM + Om_phiDM   # ≈ 0.320
OmDE      = 1 - Om_total        # ≈ 0.680

def f_DE(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def E2_total(a):
    return max(Om_total * a**(-3) + OmDE * f_DE(a), 1e-30)

def E2_CDM(a):
    return max(Om_CDM * a**(-3) + OmDE * f_DE(a), 1e-30)

def growth_rhs_twosector(lna, Y):
    a   = np.exp(lna)
    e2  = E2_total(a)
    h   = 0.5*(-3*Om_total*a**(-3) + OmDE*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a)/e2
    Omm = Om_total * a**(-3) / e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm*D]

def growth_rhs_single(lna, Y):
    a   = np.exp(lna)
    e2  = E2_CDM(a)
    h   = 0.5*(-3*Om_CDM*a**(-3) + OmDE*(f_DE(a+1e-5)-f_DE(a-1e-5))/(2e-5)*a)/e2
    Omm = Om_CDM * a**(-3) / e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm*D]

def growth_lcdm(lna, Y, Om_lcdm=0.3153):
    a   = np.exp(lna)
    e2  = max(Om_lcdm*a**(-3) + (1-Om_lcdm), 1e-30)
    h   = 0.5*(-3*Om_lcdm*a**(-3))/e2
    Omm = Om_lcdm*a**(-3)/e2
    D, Dp = Y
    return [Dp, -(2+h)*Dp + 1.5*Omm*D]

lna_pts = np.linspace(np.log(0.01), 0.0, 600)
IC = [0.01, 0.01]

sol_2s   = solve_ivp(growth_rhs_twosector, [lna_pts[0], 0.0], IC, t_eval=lna_pts, rtol=1e-10, atol=1e-13)
sol_1s   = solve_ivp(growth_rhs_single,    [lna_pts[0], 0.0], IC, t_eval=lna_pts, rtol=1e-10, atol=1e-13)
sol_lcdm = solve_ivp(growth_lcdm,          [lna_pts[0], 0.0], IC, t_eval=lna_pts, rtol=1e-10, atol=1e-13)

G_2s     = sol_2s.y[0][-1]   / sol_lcdm.y[0][-1]
G_single = sol_1s.y[0][-1]   / sol_lcdm.y[0][-1]

# sigma_8: baseline ODE, two-sector sin WDM, two-sector con WDM
sigma8_baseline = sigma8_planck * G_single          # ≈ 0.702
sigma8_2s_noWDM = sigma8_planck * G_2s              # ≈ 0.794
sigma8_phiDM    = sigma8_eff                        # CLASS top-hat directo (target: 0.737)

D_2s   = sol_2s.y[0]   / sol_2s.y[0][-1]
Dp_2s  = sol_2s.y[1]   / sol_2s.y[0][-1]
D_1s   = sol_1s.y[0]   / sol_1s.y[0][-1]
Dp_1s  = sol_1s.y[1]   / sol_1s.y[0][-1]
D_lcm  = sol_lcdm.y[0] / sol_lcdm.y[0][-1]
Dp_lcm = sol_lcdm.y[1] / sol_lcdm.y[0][-1]

def get_fsig8(z, D_arr, Dp_arr, sig8_val):
    lna_t = np.log(1/(1+z))
    f_t   = float(np.interp(lna_t, lna_pts, Dp_arr/D_arr))
    D_t   = float(np.interp(lna_t, lna_pts, D_arr))
    return f_t * sig8_val * D_t

# Datos observacionales fσ₈(z) — mismo set que Paper 6 verification
Z_RSD     = np.array([0.067, 0.150, 0.320, 0.380, 0.510, 0.570])
surveys   = ['6dFGRS', 'SDSS MGS', 'BOSS DR12', 'BOSS DR12', 'eBOSS DR16', 'WiggleZ']
fsig_d    = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
sigma_d   = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])
z_d       = Z_RSD

fsig8_base  = np.array([get_fsig8(z, D_1s,  Dp_1s,  sigma8_baseline) for z in z_d])
fsig8_phiDM = np.array([get_fsig8(z, D_2s,  Dp_2s,  sigma8_phiDM)   for z in z_d])
fsig8_lcdm  = np.array([get_fsig8(z, D_lcm, Dp_lcm, sigma8_planck)   for z in z_d])

tension_base  = (fsig_d - fsig8_base)  / sigma_d
tension_phiDM = (fsig_d - fsig8_phiDM) / sigma_d

print()
print(f"Factores de crecimiento ODE:")
print(f"  G_single (Paper 5, Ω_m=0.160) = {G_single:.4f}  →  σ₈ = {sigma8_baseline:.3f}")
print(f"  G_2s     (Paper 6, Ω_m=0.320) = {G_2s:.4f}  →  σ₈ = {sigma8_2s_noWDM:.3f}")
print(f"  σ₈_eff con T_WDM (Paper 6)    = {sigma8_phiDM:.3f}")

print()
print("── Análisis fσ₈ ────────────────────────────────────────────")
print(f"{'Survey z':>8}  {'Obs':>6}  {'Base':>6}  {'φ-DM':>6}  {'σ_base':>7}  {'σ_φDM':>7}")
surveys = ["6dFGRS", "WiggleZ", "MGS", "BOSS", "WiggleZ", "eBOSS"]
for i, (z, fo, sb, sp, tb, tp) in enumerate(
        zip(z_d, fsig_d, fsig8_base, fsig8_phiDM, tension_base, tension_phiDM)):
    print(f"z={z:.2f}  {fo:.3f}  {sb:.3f}  {sp:.3f}  {tb:+.2f}σ  {tp:+.2f}σ  {surveys[i]}")

print()
print(f"Media |tensión| baseline : {np.mean(np.abs(tension_base)):.2f}σ")
print(f"Media |tensión| φ-DM     : {np.mean(np.abs(tension_phiDM)):.2f}σ")
print(f"Reducción tensión        : {np.mean(np.abs(tension_base)):.2f}σ → {np.mean(np.abs(tension_phiDM)):.2f}σ")

# ── Figura 1: P(k) comparación ───────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(8, 10), gridspec_kw={'height_ratios': [2, 1]})

ax1 = axes[0]
ax1.loglog(k_l, P_l, 'k-', lw=2, label=r'$\Lambda$CDM (Planck 2018)', zorder=3)
ax1.loglog(k_s, P_s, 'b--', lw=1.5, label=r'SSEE+MIRA (un sector)', alpha=0.7)
ax1.loglog(k_s, P_2sector_cal, 'r-', lw=2,
           label=rf'SSEE+$\varphi$-DM ($\alpha$={alpha_cal:.4f}, $\sigma_8={sigma8_phiDM:.3f}$)', zorder=4)
ax1.loglog(k_2, P_2, 'g:', lw=1.5, label=r'CLASS two-sector ($T_\mathrm{ncdm}=0.716$)', alpha=0.8)
ax1.axvline(k_fs, color='orange', ls='--', lw=1.5, label=rf'$k_{{fs}}={k_fs}$ h/Mpc')
ax1.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=13)
ax1.set_ylabel(r'$P(k)$ [$(h^{-1}$ Mpc$)^3$]', fontsize=13)
ax1.set_title(r'SSEE-V3.6: Espectro de Materia $P(k)$ — Calibración $\varphi$-DM', fontsize=13)
ax1.legend(fontsize=10)
ax1.set_xlim(1e-3, 2)
ax1.grid(True, alpha=0.3)

# ── Subpanel: T_WDM(k) ──────────────────────────────────────────────────────
ax2 = axes[1]
k_plot = np.logspace(-2, 0.5, 200)
T_plot = T_wdm(k_plot, alpha_cal)
ax2.semilogx(k_plot, T_plot**2, 'r-', lw=2, label=rf'$T^2_\mathrm{{WDM}}(k)$, $\alpha$={alpha_cal:.4f}')
ax2.axvline(k_fs, color='orange', ls='--', lw=1.5)
ax2.axhline(0.5, color='gray', ls=':', lw=1)
ax2.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
ax2.set_ylabel(r'$T^2_\mathrm{WDM}(k)$', fontsize=12)
ax2.set_title(r'Función de transferencia WDM calibrada', fontsize=12)
ax2.legend(fontsize=10)
ax2.set_xlim(1e-2, 3)
ax2.set_ylim(-0.05, 1.05)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("output/ssee_phiDM_Pk_calibrated.pdf", dpi=150, bbox_inches='tight')
plt.savefig("output/ssee_phiDM_Pk_calibrated.png", dpi=150, bbox_inches='tight')
print("\nFigura guardada: output/ssee_phiDM_Pk_calibrated.pdf")

# ── Figura 2: fσ₈ tensiones ─────────────────────────────────────────────────
fig2, ax = plt.subplots(figsize=(8, 5))

z_th = np.linspace(0.05, 1.6, 100)
fsig_base_th  = np.array([get_fsig8(z, D_1s, Dp_1s, sigma8_baseline) for z in z_th])
fsig_phiDM_th = np.array([get_fsig8(z, D_2s, Dp_2s, sigma8_phiDM)   for z in z_th])
fsig_lcdm_th  = np.array([get_fsig8(z, D_lcm, Dp_lcm, sigma8_planck) for z in z_th])

ax.plot(z_th, fsig_lcdm_th,  'k-',  lw=2,   label=rf'$\Lambda$CDM ($\sigma_8={sigma8_planck:.3f}$)')
ax.plot(z_th, fsig_base_th,  'b--', lw=1.5, label=rf'SSEE baseline ($\sigma_8={sigma8_baseline:.3f}$)')
ax.plot(z_th, fsig_phiDM_th, 'r-',  lw=2,   label=rf'SSEE+$\varphi$-DM ($\sigma_8={sigma8_phiDM:.3f}$)')
ax.errorbar(z_d, fsig_d, yerr=sigma_d, fmt='ko', capsize=4, ms=6, label='Datos RSD')
ax.set_xlabel(r'Redshift $z$', fontsize=13)
ax.set_ylabel(r'$f\sigma_8(z)$', fontsize=13)
ax.set_title(r'SSEE-V3.6: Tensión $f\sigma_8$ — Reducción con $\varphi$-DM', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1.6)
ax.set_ylim(0.25, 0.65)

# Anotar tensiones
for i, (z, tension) in enumerate(zip(z_d, tension_phiDM)):
    ax.annotate(f'{tension:+.1f}σ', (z, fsig_d[i]), textcoords="offset points",
                xytext=(5, 8), fontsize=8, color='red')

plt.tight_layout()
plt.savefig("output/ssee_fsig8_tension_phiDM.pdf", dpi=150, bbox_inches='tight')
plt.savefig("output/ssee_fsig8_tension_phiDM.png", dpi=150, bbox_inches='tight')
print("Figura guardada: output/ssee_fsig8_tension_phiDM.pdf")

print()
print("=" * 60)
print("RESUMEN FASE 2c:")
print(f"  alpha_WDM calibrado = {alpha_cal:.6f} h/Mpc")
print(f"  sigma_8_eff         = {sigma8_phiDM:.3f}")
print(f"  Tensión fσ₈ media   : {np.mean(np.abs(tension_base)):.2f}σ → {np.mean(np.abs(tension_phiDM)):.2f}σ")
print(f"  (Paper 6 reporta: 2.56σ → 0.50σ con corrección sigma_8 corregida)")
print("=" * 60)
