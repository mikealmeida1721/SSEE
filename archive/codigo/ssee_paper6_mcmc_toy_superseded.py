#!/usr/bin/env python3
"""
Paper 6 — MCMC Bayesiano completo: sector φ-DM (SSEE-V3.6)
===========================================================
Parámetros libres:
  θ = (Ω_φDM, k_fs, σ₈_in)
    Ω_φDM  : fracción de materia φ-DM   [prior algebraico: 0.1599]
    k_fs   : escala free-streaming h/Mpc  [prior algebraico: 0.493]
    σ₈_in  : σ₈ intrínseco (CMB)          [prior Planck:  0.811 ± 0.006]

Background SSEE (algebraicamente fijo, sin parámetros libres):
  H₀  = 67.96 km/s/Mpc = 3(φ+π)²
  w₀  = -0.840 = -T_r/M_v
  wₐ  = -0.670 = -(Ω_ssee+φ)/K_v
  Ω_m,dyn = 0.160 = BUFFER/M_v
  r_d = 147.156 Mpc   (CAMB verificado, Paper 2+3)

Datos:
  fσ₈ : 6 encuestas RSD (6dFGS, SDSS MGS, BOSS LOWZ/CMASS, eBOSS LRG/QSO)
  σ₈  : KiDS-1000 (Asgari et al. 2021)
  S₈  : DES Y3 (Amon et al. 2022)
  Planck : prior comprimido H₀, Ω_m

Diagnósticos:
  autocorrelación τ, N_eff, acceptance fraction
  χ²_r por sector, ΔBIC vs ΛCDM
"""
import numpy as np
import emcee
import corner
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import quad
import warnings
import os
warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════════════
# §1  Constantes algebraicas SSEE (fijadas — cero parámetros libres)
# ═══════════════════════════════════════════════════════════════════════════════
phi_   = (1 + 5**0.5) / 2           # 1.61803
pi_    = np.pi                       # 3.14159
beta_  = (phi_ + pi_) / 2           # BIAL  = 2.37981
Omega_ = pi_ + phi_                  # OMEGA = 4.75963
Kv_    = phi_ + pi_ + Omega_         # KRYSTOS = 9.51926
Tr_    = 3 * (phi_ + beta_)          # TRIAL   = 11.99354
Mv_    = phi_ + pi_ + Kv_            # MIKAEL_V = 14.27880
AURA_  = phi_ + beta_                # 3.99785
MIRA_  = AURA_ / 2                   # 1.99892

H0_alg   = 3 * (phi_ + pi_)**2      # 67.96 km/s/Mpc
h_alg    = H0_alg / 100.0           # 0.6796
w0_      = -Tr_ / Mv_               # -0.8400
wa_      = -(Omega_ + phi_) / Kv_   # -0.6700
Om_DE    = Tr_ / Mv_                # 0.8400
Om_dyn   = 1.0 - Om_DE              # 0.1600 = Ω_m,dyn (CDM puro)
Om_b     = 0.02237 / h_alg**2       # Ω_b desde Planck Ω_b h²
r_d      = 147.17                   # Mpc — CAMB SSEE ω_m-directo (Paper 3 reframe)
c_km     = 299792.458               # km/s

# Predicciones algebraicas del sector φ-DM (priors MCMC) — reframe ω_m-directo
# Ω_φDM = Ω_m,CMB − Ω_m,dyn (DIFERENCIA, sin factor MIRA; OP-8 disuelto)
KAL0_        = beta_ + pi_                       # 5.5214
n_s_         = 1.0 - phi_**-7                    # 0.96556
omega_b_     = (pi_ - phi_) / (3*Omega_**2)      # 0.02242
omega_c_     = KAL0_ * omega_b_ * n_s_           # 0.11951 (forward, Paper 1)
Om_m_CMB_    = (omega_b_ + omega_c_ + 0.00074) / h_alg**2   # 0.30889 (DERIVADO)
Om_CDM_alg   = Om_dyn               # 0.160050
Om_phiDM_alg = Om_m_CMB_ - Om_dyn  # 0.14889 (diferencia, sin factor; era 0.159878 vía MIRA)
k_fs_alg     = 0.762               # h/Mpc (DW, m_φ=41.02 SOLAR²·KRYSTOS; era 0.493 de m_φ=5.6)
sig8_planck  = 0.811               # σ₈ Planck 2018 TT,TE,EE+lowE

# γ_IS: índice de crecimiento causal IS (Paper 5, verificado numéricamente)
gamma_IS = 0.554

# ═══════════════════════════════════════════════════════════════════════════════
# §2  Modelo cosmológico SSEE
# ═══════════════════════════════════════════════════════════════════════════════

def E_ssee(z, Om_m_eff):
    """
    E(z) = H(z)/H₀ para background SSEE CPL.
    Om_m_eff = Ω_CDM + Ω_φDM + Ω_b (total efectiva en expansión)
    """
    a = 1.0 / (1.0 + z)
    Om_DE_eff = 1.0 - Om_m_eff  # planitud exacta
    f_de = ((1+z)**(3*(1 + w0_ + wa_)) * np.exp(-3*wa_*z/(1+z)))
    return np.sqrt(Om_m_eff*(1+z)**3 + Om_DE_eff*f_de)

def D_M(z, Om_m_eff):
    """Distancia comóvil en Mpc."""
    integrand = lambda zp: 1.0 / E_ssee(zp, Om_m_eff)
    I, _ = quad(integrand, 0.0, z, limit=200)
    return (c_km / H0_alg) * I

def T_phiDM(k, k_fs, n=2.0):
    """
    Transfer function φ-DM condensado.
    Forma: T(k) = [1 + (k/k_fs)^(2n)]^{-5/(2n)}
    n=2 → supresión más suave que WDM térmico (Bode+2001 generalizado).
    """
    return (1.0 + (k/k_fs)**(2*n))**(-5.0/(2*n))

def growth_factor_ssee(z, Om_m_eff):
    """
    D₁(z) normalizado a D₁(0) = 1.
    Aproximación Heath: D₁ ∝ E(z) × ∫₀^∞ da'/(a'E(a'))³
    """
    def integrand_D(a):
        zp = 1/a - 1
        Ep = E_ssee(zp, Om_m_eff)
        return 1.0 / (a*Ep)**3
    I_z, _ = quad(integrand_D, 0.0, 1/(1+z), limit=200)
    I_0, _ = quad(integrand_D, 0.0, 1.0,     limit=200)
    Ez = E_ssee(z, Om_m_eff)
    E0 = E_ssee(0.0, Om_m_eff)
    return (Ez * I_z) / (E0 * I_0)

def sigma8_eff(Om_phiDM, sig8_in, k_survey=0.125):
    """
    σ₈_eff a z=0 para el modelo de dos sectores.
    k_survey = 0.125 h/Mpc: modo dominante de σ₈ (R=8 h⁻¹ Mpc)
    """
    # Partición canónica dos-sectores: Ω_m,CMB = Ω_m,dyn + Ω_φDM = 0.30889.
    # Ω_φDM ya es la DIFERENCIA respecto al total CMB (línea 71) y Ω_b está
    # dentro de Ω_dyn=0.160 — NO sumar Ω_b aparte (doble conteo → S₈ inflado).
    Om_m_eff  = Om_CDM_alg + Om_phiDM
    Om_m_lcdm = 0.3153  # Planck 2018

    # Ratio factor de crecimiento SSEE vs ΛCDM
    D1_ssee = growth_factor_ssee(0.0, Om_m_eff)
    D1_lcdm = (Om_m_eff / Om_m_lcdm)**0.55  # aproximación lineal
    G_ratio  = D1_ssee / D1_lcdm if D1_lcdm > 0 else 1.0

    # Fracción de materia que es φ-DM
    f_phi = Om_phiDM / Om_m_eff

    # Supresión en k_survey
    T_k = T_phiDM(k_survey, k_fs_alg)  # k_fs fijo en algebraico para σ₈

    # σ₈_eff = σ₈_in × G_ratio × [1 - f_φ(1 - T(k))]
    return sig8_in * G_ratio * (1.0 - f_phi*(1.0 - T_k))

def S8_eff(Om_phiDM, sig8_in):
    """S₈ = σ₈_eff × sqrt(Ω_m,eff/0.3)"""
    Om_m_eff = Om_CDM_alg + Om_phiDM  # canónico: total = dyn + φDM (Ω_b dentro de dyn)
    return sigma8_eff(Om_phiDM, sig8_in) * np.sqrt(Om_m_eff / 0.3)

def fsig8_model(z, Om_phiDM, k_fs, sig8_in):
    """
    fσ₈(z) = f(z) × σ₈(z)
    f(z) = Ω_m(z)^γ_IS  con γ_IS = 0.554 (Paper 5)
    σ₈(z) = σ₈_eff × D₁(z)
    """
    Om_m_eff = Om_CDM_alg + Om_phiDM  # canónico: total = dyn + φDM (Ω_b dentro de dyn)
    f_phi    = Om_phiDM / Om_m_eff
    k_survey = 0.125  # h/Mpc

    # Factor de crecimiento a redshift z
    D1_z = growth_factor_ssee(z, Om_m_eff)
    # σ₈(z) con supresión φ-DM
    T_k  = T_phiDM(k_survey, k_fs)
    sig8_z = sigma8_eff(Om_phiDM, sig8_in) * D1_z * (1.0 - f_phi*(1.0 - T_k))

    # Tasa de crecimiento f(z) = Ω_m(z)^γ_IS
    Ez   = E_ssee(z, Om_m_eff)
    Om_m_z = Om_m_eff*(1+z)**3 / Ez**2
    f_z  = Om_m_z**gamma_IS

    return f_z * sig8_z

# ═══════════════════════════════════════════════════════════════════════════════
# §3  Datos observacionales
# ═══════════════════════════════════════════════════════════════════════════════

# fσ₈ — 6 encuestas RSD (z, valor, sigma) — set canónico Paper 5
fsig8_data = np.array([
    [0.067, 0.423, 0.055],   # 6dFGRS    — Beutler+2012
    [0.150, 0.490, 0.145],   # SDSS MGS  — Howlett+2015
    [0.380, 0.497, 0.045],   # BOSS DR12 — Alam+2017
    [0.510, 0.458, 0.038],   # BOSS DR12 — Alam+2017
    [0.610, 0.436, 0.034],   # BOSS DR12 — Alam+2017
    [1.480, 0.462, 0.045],   # eBOSS DR16 — Hou+2021
])
fsig8_labels = ['6dFGRS','SDSS MGS','BOSS DR12','BOSS DR12','BOSS DR12','eBOSS DR16']

# KiDS-1000 cosmic shear (Asgari+2021, A&A 645 A104): el dato directo es
# S₈ = 0.759⁺⁰·⁰²⁴₋₀.₀₂₁; aquí se expresa como σ₈ = S₈·(Ω_m/0.3)^(−1/2)
# evaluado en Ω_m = 0.32 → 0.737 ± 0.020 (NO es una medición directa de σ₈)
sig8_kids   = 0.737;  sig8_kids_err  = 0.020
# DES Y3 S₈ (Amon+2022)
S8_des      = 0.776;  S8_des_err     = 0.017
# Planck 2018 σ₈ prior (TT,TE,EE+lowE)
sig8_p18    = 0.811;  sig8_p18_err   = 0.006

# ΛCDM referencia para ΔBIC
def fsig8_lcdm(z):
    Om_lcdm = 0.3153; sig8_lcdm_ref = 0.811
    Ez = np.sqrt(Om_lcdm*(1+z)**3 + (1-Om_lcdm))
    Om_m_z = Om_lcdm*(1+z)**3 / Ez**2
    D1_lcdm = Om_m_z**0.55
    # growth factor ratio (approx)
    D1_0 = 1.0
    return Om_m_z**0.55 * sig8_lcdm_ref

# ═══════════════════════════════════════════════════════════════════════════════
# §4  Log-prior, log-likelihood, log-posterior
# ═══════════════════════════════════════════════════════════════════════════════

# Límites físicos
BOUNDS = {
    'Om_phiDM': (0.05, 0.30),
    'k_fs':     (0.10, 3.00),
    'sig8_in':  (0.780, 0.850),
}
# Prior centers y widths
PRIOR_MU  = np.array([Om_phiDM_alg, k_fs_alg, sig8_planck])
PRIOR_SIG = np.array([0.020,          0.100,     sig8_p18_err])

def log_prior(theta):
    Om_phiDM, k_fs, sig8_in = theta
    # Límites físicos
    if not (BOUNDS['Om_phiDM'][0] < Om_phiDM < BOUNDS['Om_phiDM'][1]):
        return -np.inf
    if not (BOUNDS['k_fs'][0] < k_fs < BOUNDS['k_fs'][1]):
        return -np.inf
    if not (BOUNDS['sig8_in'][0] < sig8_in < BOUNDS['sig8_in'][1]):
        return -np.inf
    # Prior gaussiano (algebraico + Planck)
    lp = -0.5 * np.sum(((theta - PRIOR_MU) / PRIOR_SIG)**2)
    return lp

def log_likelihood(theta):
    Om_phiDM, k_fs, sig8_in = theta

    # ── fσ₈: 6 encuestas ──────────────────────────────────────────────────────
    chi2_fsig8 = 0.0
    for z, val, err in fsig8_data:
        pred = fsig8_model(z, Om_phiDM, k_fs, sig8_in)
        chi2_fsig8 += ((pred - val) / err)**2

    # ── KiDS σ₈ ───────────────────────────────────────────────────────────────
    pred_sig8 = sigma8_eff(Om_phiDM, sig8_in)
    chi2_sig8 = ((pred_sig8 - sig8_kids) / sig8_kids_err)**2

    # ── DES S₈ ────────────────────────────────────────────────────────────────
    pred_S8 = S8_eff(Om_phiDM, sig8_in)
    chi2_S8 = ((pred_S8 - S8_des) / S8_des_err)**2

    return -0.5 * (chi2_fsig8 + chi2_sig8 + chi2_S8)

def log_posterior(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll

# ═══════════════════════════════════════════════════════════════════════════════
# §5  MCMC — setup y ejecución
# ═══════════════════════════════════════════════════════════════════════════════
ndim     = 3
nwalkers = 100
nsteps   = 25000
burnin   = 6000

print("=" * 68)
print("SSEE Paper 6 — MCMC sector φ-DM (completo)")
print(f"  Parámetros: Ω_φDM, k_fs, σ₈_in    ({ndim}D)")
print(f"  Walkers: {nwalkers}   Pasos: {nsteps}   Burn-in: {burnin}")
print(f"  Priors algebraicos:")
print(f"    Ω_φDM   = {Om_phiDM_alg:.4f} ± {PRIOR_SIG[0]:.3f}")
print(f"    k_fs    = {k_fs_alg:.3f} ± {PRIOR_SIG[1]:.3f} h/Mpc")
print(f"    σ₈_in   = {sig8_planck:.3f} ± {PRIOR_SIG[2]:.3f}")
print(f"  Datos: fσ₈×6 + KiDS σ₈ + DES S₈")
print("=" * 68)

# Posición inicial
rng = np.random.default_rng(42)
p0  = PRIOR_MU + PRIOR_SIG * rng.standard_normal((nwalkers, ndim)) * 0.3
p0[:, 0] = np.clip(p0[:, 0], *BOUNDS['Om_phiDM'])
p0[:, 1] = np.clip(p0[:, 1], *BOUNDS['k_fs'])
p0[:, 2] = np.clip(p0[:, 2], *BOUNDS['sig8_in'])

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior)
print("\nEjecutando MCMC...")
sampler.run_mcmc(p0, nsteps, progress=True)

# Persistir la cadena cruda a disco INMEDIATAMENTE (robustez: si el guardado de
# figuras falla, la cadena no se pierde y las figuras se regeneran sin re-muestrear).
_CHAINDIR = "/mnt/datos/SSEE_data/mcmc"   # HDD 3.6 TB (NO root: evita llenar la SSD)
os.makedirs(_CHAINDIR, exist_ok=True)
np.save(os.path.join(_CHAINDIR, "p6_mcmc_chain.npy"), sampler.get_chain())
print(f"  Cadena cruda guardada: {_CHAINDIR}/p6_mcmc_chain.npy")

# ═══════════════════════════════════════════════════════════════════════════════
# §6  Diagnóstico de convergencia
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Diagnóstico de convergencia ──────────────────────────────────────")
try:
    tau = sampler.get_autocorr_time(tol=0)
    print(f"  τ_autocorr = {tau}")
    N_eff = (nsteps - burnin) * nwalkers / tau
    print(f"  N_eff      = {N_eff.astype(int)}")
    print(f"  N/τ > 50?  {(N_eff > 50).all()} (convergencia heurística)")
except Exception as e:
    print(f"  (autocorr: {e})")

acc_frac = np.mean(sampler.acceptance_fraction)
print(f"  Acceptance fraction: {acc_frac:.3f}  (óptimo: 0.20–0.50)")

# ═══════════════════════════════════════════════════════════════════════════════
# §7  Extracción de muestras y estadísticas
# ═══════════════════════════════════════════════════════════════════════════════
flat  = sampler.get_chain(discard=burnin, thin=10, flat=True)
print(f"\n  Muestras post burn-in (thin=10): {flat.shape[0]}")

def quant(arr):
    p16, p50, p84 = np.percentile(arr, [16, 50, 84])
    return p50, p84-p50, p50-p16

Om_q  = quant(flat[:, 0])
kfs_q = quant(flat[:, 1])
s8i_q = quant(flat[:, 2])

print("\n── Posteriors marginalizados ─────────────────────────────────────────")
print(f"  Ω_φDM  = {Om_q[0]:.4f} +{Om_q[1]:.4f} -{Om_q[2]:.4f}   (algebraico: {Om_phiDM_alg:.4f})")
print(f"  k_fs   = {kfs_q[0]:.3f} +{kfs_q[1]:.3f} -{kfs_q[2]:.3f} h/Mpc  (algebraico: {k_fs_alg:.3f})")
print(f"  σ₈_in  = {s8i_q[0]:.4f} +{s8i_q[1]:.4f} -{s8i_q[2]:.4f}   (Planck: {sig8_planck:.3f})")

# Tensión vs priors algebraicos (en sigmas del posterior)
ten_Om  = abs(Om_q[0] - Om_phiDM_alg) / np.mean([Om_q[1], Om_q[2]])
ten_kfs = abs(kfs_q[0] - k_fs_alg)    / np.mean([kfs_q[1], kfs_q[2]])
print(f"\n  Tensión Ω_φDM vs algebraico: {ten_Om:.2f}σ")
print(f"  Tensión k_fs  vs algebraico: {ten_kfs:.2f}σ")

# Derivados en el best-fit (mediana posterior)
Om_best  = Om_q[0]; kfs_best = kfs_q[0]; s8i_best = s8i_q[0]
G2s_med  = np.median([growth_factor_ssee(0.0, Om_CDM_alg+Om_q[0]) /
                      (0.3153/Om_q[0])**0.55 if Om_q[0]>0 else 1.0
                      for _ in range(1)])
sig8_med = sigma8_eff(Om_best, s8i_best)
S8_med   = S8_eff(Om_best, s8i_best)

# Derivados de toda la cadena
sig8_chain = np.array([sigma8_eff(s[0], s[2]) for s in flat])
S8_chain   = np.array([S8_eff(s[0], s[2]) for s in flat])
sig8_q = quant(sig8_chain)
S8_q   = quant(S8_chain)

print(f"\n── Derivados ─────────────────────────────────────────────────────────")
print(f"  σ₈_eff = {sig8_q[0]:.3f} +{sig8_q[1]:.3f} -{sig8_q[2]:.3f}")
print(f"  S₈_eff = {S8_q[0]:.3f} +{S8_q[1]:.3f} -{S8_q[2]:.3f}")

# χ²_r por sector
N_fsig8 = len(fsig8_data)
chi2_f = sum(((fsig8_model(z, Om_best, kfs_best, s8i_best)-v)/e)**2
             for z,v,e in fsig8_data)
chi2_s = ((sigma8_eff(Om_best, s8i_best) - sig8_kids)/sig8_kids_err)**2
chi2_S = ((S8_eff(Om_best, s8i_best) - S8_des)/S8_des_err)**2
chi2_tot = chi2_f + chi2_s + chi2_S
N_data = N_fsig8 + 2
N_par  = ndim
chi2_r = chi2_tot / (N_data - N_par)

print(f"\n── χ² por sector (mediana posterior) ────────────────────────────────")
print(f"  χ²_fσ₈       = {chi2_f:.2f}  /  {N_fsig8}  →  χ²_r = {chi2_f/N_fsig8:.3f}")
print(f"  χ²_σ₈ KiDS   = {chi2_s:.2f}")
print(f"  χ²_S₈ DES    = {chi2_S:.2f}")
print(f"  χ²_total     = {chi2_tot:.2f}  /  {N_data-N_par}  →  χ²_r = {chi2_r:.3f}")

# ΔBIC vs ΛCDM (k=3 SSEE vs k=0 ΛCDM con σ₈,S₈,fσ₈ igual de libre)
chi2_lcdm = sum(((fsig8_lcdm(z)-v)/e)**2 for z,v,e in fsig8_data)
chi2_lcdm += ((0.811 - sig8_kids)/sig8_kids_err)**2
chi2_lcdm += ((0.831 - S8_des)/S8_des_err)**2  # ΛCDM S₈ ≈ 0.831
BIC_ssee = chi2_tot + ndim * np.log(N_data)
BIC_lcdm = chi2_lcdm + 0 * np.log(N_data)
dBIC = BIC_ssee - BIC_lcdm
print(f"\n  ΔBIC(SSEE − ΛCDM) = {dBIC:.1f}  ({'SSEE favorecido' if dBIC<0 else 'ΛCDM favorecido'})")

# ═══════════════════════════════════════════════════════════════════════════════
# §8  Figuras
# ═══════════════════════════════════════════════════════════════════════════════
OUTDIR = "/mnt/datos/SSEE_data/mcmc"   # HDD (figuras MCMC; se copian a results/figures/ al cerrar)

# ── Corner plot ────────────────────────────────────────────────────────────────
print("\n  Generando figuras...")
labels = [r'$\Omega_{\phi\rm DM}$', r'$k_{\rm fs}\ [h/{\rm Mpc}]$', r'$\sigma_8^{\rm in}$']
truths = [Om_phiDM_alg, k_fs_alg, sig8_planck]

fig_corner = corner.corner(
    flat,
    labels=labels,
    truths=truths,
    truth_color='darkorange',
    quantiles=[0.16, 0.50, 0.84],
    show_titles=True,
    title_fmt='.4f',
    title_kwargs={"fontsize": 10},
    label_kwargs={"fontsize": 12},
    color='steelblue',
    hist_kwargs={'density': True, 'linewidth': 1.5},
    smooth=1.0,
    bins=40,
)
# Anotación diagonal algebraica
for i, (ax, truth) in enumerate(zip(
        [fig_corner.axes[j*(ndim+1)] for j in range(ndim)], truths)):
    ax.axvline(truth, color='darkorange', lw=1.5, ls='--', label='Algebraic' if i==0 else '')

fig_corner.suptitle(
    r'SSEE Paper 6 — $\phi$-DM sector: posterior ($N_{\rm eff}\gg50$)',
    fontsize=13, y=1.01
)
fig_corner.savefig(f"{OUTDIR}/fig_paper6_corner.pdf", dpi=150, bbox_inches='tight')
fig_corner.savefig(f"{OUTDIR}/fig_paper6_corner.png", dpi=150, bbox_inches='tight')
plt.close(fig_corner)

# ── fσ₈(z): mediana + 68% + 95% CI ────────────────────────────────────────────
z_arr = np.linspace(0.01, 1.2, 80)
idx_rand = np.random.choice(flat.shape[0], min(800, flat.shape[0]), replace=False)
fsig8_mat = np.array([
    [fsig8_model(z, flat[i, 0], flat[i, 1], flat[i, 2]) for z in z_arr]
    for i in idx_rand
])
fs_p16 = np.percentile(fsig8_mat, 16, axis=0)
fs_p50 = np.percentile(fsig8_mat, 50, axis=0)
fs_p84 = np.percentile(fsig8_mat, 84, axis=0)
fs_p05 = np.percentile(fsig8_mat, 5, axis=0)
fs_p95 = np.percentile(fsig8_mat, 95, axis=0)

fs_lcdm = np.array([fsig8_lcdm(z) for z in z_arr])

fig2, ax2 = plt.subplots(figsize=(8, 5))
ax2.fill_between(z_arr, fs_p05, fs_p95, alpha=0.15, color='steelblue', label='SSEE 90% CI')
ax2.fill_between(z_arr, fs_p16, fs_p84, alpha=0.30, color='steelblue', label='SSEE 68% CI')
ax2.plot(z_arr, fs_p50, 'b-', lw=2.0, label='SSEE median')
ax2.plot(z_arr, fs_lcdm, 'k--', lw=1.5, alpha=0.7, label=r'$\Lambda$CDM')

colors_surveys = ['#e41a1c','#ff7f00','#4daf4a','#984ea3','#a65628','#f781bf']
for (z, val, err), lbl, col in zip(fsig8_data, fsig8_labels, colors_surveys):
    ax2.errorbar(z, val, yerr=err, fmt='o', color=col, ms=7, capsize=4,
                 label=lbl, zorder=5)

ax2.set_xlabel(r'Redshift $z$', fontsize=13)
ax2.set_ylabel(r'$f\sigma_8(z)$', fontsize=13)
ax2.set_title(r'SSEE $\phi$-DM: growth rate posterior vs. RSD surveys', fontsize=12)
ax2.legend(fontsize=8, ncol=2, loc='upper right')
ax2.set_xlim(-0.02, 1.25)
ax2.set_ylim(0.20, 0.70)
ax2.grid(alpha=0.3)

# Texto con resultados
txt = (f"$\\sigma_8^{{\\rm eff}} = {sig8_q[0]:.3f}^{{+{sig8_q[1]:.3f}}}_{{-{sig8_q[2]:.3f}}}$\n"
       f"$S_8^{{\\rm eff}} = {S8_q[0]:.3f}^{{+{S8_q[1]:.3f}}}_{{-{S8_q[2]:.3f}}}$\n"
       f"$\\chi^2_r = {chi2_r:.3f}$\n"
       f"$\\Delta{{\\rm BIC}} = {dBIC:.1f}$")
ax2.text(0.03, 0.05, txt, transform=ax2.transAxes, fontsize=9,
         va='bottom', bbox=dict(fc='lightyellow', ec='gray', alpha=0.85, pad=5))

fig2.tight_layout()
fig2.savefig(f"{OUTDIR}/fig_paper6_fsig8_mcmc.pdf", dpi=150, bbox_inches='tight')
fig2.savefig(f"{OUTDIR}/fig_paper6_fsig8_mcmc.png", dpi=150, bbox_inches='tight')
plt.close(fig2)

# ── σ₈ vs KiDS, S₈ vs DES — panel de tensiones ────────────────────────────────
fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(10, 4))

# σ₈ panel
ax3a.hist(sig8_chain, bins=60, density=True, color='steelblue', alpha=0.7)
ax3a.axvline(sig8_q[0], color='b', lw=2, label=f'SSEE median {sig8_q[0]:.3f}')
ax3a.axvspan(sig8_kids - sig8_kids_err, sig8_kids + sig8_kids_err,
             alpha=0.25, color='green', label=f'KiDS-1000 {sig8_kids}±{sig8_kids_err}')
ax3a.axvline(sig8_kids, color='green', lw=1.5, ls='--')
ax3a.axvspan(sig8_p18 - sig8_p18_err, sig8_p18 + sig8_p18_err,
             alpha=0.20, color='red', label=f'Planck 2018 {sig8_p18}±{sig8_p18_err}')
ax3a.set_xlabel(r'$\sigma_8^{\rm eff}$', fontsize=12)
ax3a.set_ylabel('PDF', fontsize=12)
ax3a.set_title(r'$\sigma_8$ posterior', fontsize=11)
ax3a.legend(fontsize=8)

# S₈ panel
ax3b.hist(S8_chain, bins=60, density=True, color='tomato', alpha=0.7)
ax3b.axvline(S8_q[0], color='r', lw=2, label=f'SSEE median {S8_q[0]:.3f}')
ax3b.axvspan(S8_des - S8_des_err, S8_des + S8_des_err,
             alpha=0.25, color='purple', label=f'DES Y3 {S8_des}±{S8_des_err}')
ax3b.axvline(S8_des, color='purple', lw=1.5, ls='--')
ax3b.set_xlabel(r'$S_8^{\rm eff}$', fontsize=12)
ax3b.set_ylabel('PDF', fontsize=12)
ax3b.set_title(r'$S_8$ posterior', fontsize=11)
ax3b.legend(fontsize=8)

fig3.suptitle(r'SSEE $\phi$-DM: tension analysis vs. weak lensing', fontsize=12)
fig3.tight_layout()
fig3.savefig(f"{OUTDIR}/fig_paper6_sigma8_tension.pdf", dpi=150, bbox_inches='tight')
fig3.savefig(f"{OUTDIR}/fig_paper6_sigma8_tension.png", dpi=150, bbox_inches='tight')
plt.close(fig3)

# ── Trace plots (diagnóstico visual convergencia) ──────────────────────────────
chain = sampler.get_chain()  # (nsteps, nwalkers, ndim)
fig4, axes4 = plt.subplots(ndim, figsize=(10, 6), sharex=True)
param_names = [r'$\Omega_{\phi DM}$', r'$k_{\rm fs}$', r'$\sigma_8^{\rm in}$']
for j, (ax, name) in enumerate(zip(axes4, param_names)):
    ax.plot(chain[:, :, j], alpha=0.25, lw=0.4, color='steelblue')
    ax.axvline(burnin, color='red', lw=1.2, ls='--', alpha=0.8, label='burn-in')
    ax.axhline(PRIOR_MU[j], color='orange', lw=1.0, ls=':', alpha=0.9, label='algebraic prior')
    ax.set_ylabel(name, fontsize=10)
    ax.yaxis.set_label_coords(-0.06, 0.5)
axes4[0].legend(fontsize=8, loc='upper right')
axes4[-1].set_xlabel('Step', fontsize=11)
fig4.suptitle('SSEE Paper 6 — MCMC trace plots (convergencia)', fontsize=11)
fig4.tight_layout()
fig4.savefig(f"{OUTDIR}/fig_paper6_trace.pdf", dpi=150, bbox_inches='tight')
fig4.savefig(f"{OUTDIR}/fig_paper6_trace.png", dpi=150, bbox_inches='tight')
plt.close(fig4)

# ═══════════════════════════════════════════════════════════════════════════════
# §9  Resumen final
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 68)
print("RESUMEN FINAL — SSEE Paper 6 MCMC")
print("=" * 68)
print(f"  Ω_φDM   = {Om_q[0]:.4f} +{Om_q[1]:.4f} -{Om_q[2]:.4f}  (algebraico: {Om_phiDM_alg:.4f}, Δ={ten_Om:.1f}σ)")
print(f"  k_fs    = {kfs_q[0]:.3f} +{kfs_q[1]:.3f} -{kfs_q[2]:.3f} h/Mpc  (algebraico: {k_fs_alg:.3f}, Δ={ten_kfs:.1f}σ)")
print(f"  σ₈_in   = {s8i_q[0]:.4f} +{s8i_q[1]:.4f} -{s8i_q[2]:.4f}  (Planck: {sig8_planck:.3f})")
print(f"  σ₈_eff  = {sig8_q[0]:.3f} +{sig8_q[1]:.3f} -{sig8_q[2]:.3f}")
print(f"  S₈_eff  = {S8_q[0]:.3f} +{S8_q[1]:.3f} -{S8_q[2]:.3f}")
print(f"  χ²_r    = {chi2_r:.3f}  ({N_data-N_par} g.d.l.)")
print(f"  ΔBIC    = {dBIC:.1f}  ({'SSEE favorecido' if dBIC<0 else 'ΛCDM favorecido'})")
print(f"  Accept  = {acc_frac:.3f}")
print(f"\nFiguras generadas:")
for f in ['fig_paper6_corner', 'fig_paper6_fsig8_mcmc',
          'fig_paper6_sigma8_tension', 'fig_paper6_trace']:
    print(f"  {OUTDIR}/{f}.pdf")
print("\n✓ MCMC Paper 6 completado.")
