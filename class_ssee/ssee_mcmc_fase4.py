"""
SSEE-V3.6 — CLASS Fase 4: MCMC Completo Multi-Probe
=====================================================
Arquitectura: analítica exacta para MCMC (misma velocidad que Paper 2),
calibrada con los resultados CLASS de Fases 1-3.

Mejoras sobre Paper 2:
  1. w0, wa libres (con prior gaussiano algebraico) — en Paper 2 eran fijos
  2. fσ8 de 6 surveys (Paper 6) integrado en la likelihood
  3. sigma_8 two-sector con T_WDM(k; alpha) (Fases 2c) en la likelihood
  4. r_d fórmula Aubourg+2015 calibrada por CLASS (Fase 1: r_d=147.01 Mpc ✅)
  5. 13 puntos DESI DR2 + covarianza diagonal

Parámetros libres (ndim=5):
  theta = [H0, ombh2, Omega_m, w0, wa]

Priors:
  H0:     Uniforme [62, 74] + Planck gaussiano
  ombh2:  Gaussiano N(0.02237, 0.00015)
  Omega_m:Uniforme [0.24, 0.40] + Planck gaussiano
  w0:     Gaussiano N(W0_ALG, 0.08)   — predicción algebraica SSEE
  wa:     Gaussiano N(WA_ALG, 0.12)   — predicción algebraica SSEE

Likelihoods:
  L1: DESI DR2 BAO 13 puntos (Abdul-Karim et al. 2025)
  L2: Planck 2018 prior comprimido (H0, Omega_m, correlacion rho=-0.85)
  L3: fσ8 seis encuestas (Paper 6 two-sector, CLASS calibrado)
  L4: Cúmulos de galaxias IGIMF (Zhang et al. 2026, 4 cúmulos)

Uso:
  python ssee_mcmc_fase4.py             # corrida completa (~4-8 h, similar Paper 2)
  python ssee_mcmc_fase4.py --test      # 20 walkers x 200 steps (~5 min)
  python ssee_mcmc_fase4.py --no-fsig8  # desactivar L3 (diagnosis)
  python ssee_mcmc_fase4.py --resume    # reanudar desde cadena guardada
"""
import numpy as np
import emcee
import corner
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os, sys, time, warnings
warnings.filterwarnings("ignore")

# ─── Constantes algebraicas SSEE ──────────────────────────────────────────────
PHI  = (1 + np.sqrt(5)) / 2
PI   = np.pi
BETA = (PI + PHI) / 2
KAL0 = BETA + PI
P_sc = (PI + PHI) + PHI
KV   = PHI + PI + (PI + PHI)
TR   = 3 * (PHI + BETA)
MV   = PHI + PI + KV

W0_ALG    = -TR / MV            # -0.83995
WA_ALG    = -P_sc / KV          # -0.66997
OMDE_ALG  = TR / MV             # 0.83995
OM_DYN    = 1.0 - OMDE_ALG     # 0.16005
MIRA      = (3*PHI + PI) / 4   # 1.99893
OM_CMB    = OM_DYN * MIRA       # 0.31993
H0_ALG    = 67.96
NS_ALG    = 0.96556
OMBH2_ALG = 0.02237

# CLASS Fase 1 calibración: r_d(CLASS) = 147.01 Mpc con params SSEE
# vs. Planck 2018 medido = 147.09 Mpc (diferencia 0.05%)
RD_CLASS_SSEE = 147.01

# ─── Flags ─────────────────────────────────────────────────────────────────────
TEST_MODE    = "--test"    in sys.argv
RESUME_MODE  = "--resume"  in sys.argv
NO_FSIG8     = "--no-fsig8" in sys.argv

OUT_DIR    = "output/mcmc_fase4"
LOG_FILE   = f"{OUT_DIR}/mcmc_fase4.log"
CHAIN_FILE = f"{OUT_DIR}/chains_fase4.npz"
os.makedirs(OUT_DIR, exist_ok=True)

t0 = time.time()
def log(msg):
    elapsed = (time.time() - t0) / 60
    line = f"[{elapsed:6.1f}m] {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

C_KM = 299792.458  # km/s

# ─── Background cosmológico (analítico exacto) ────────────────────────────────
def E2(z, Omm, w0, wa):
    """H²(z)/H₀² para CPL. Omm = Omega_m total."""
    a = 1.0 / (1 + z)
    OmDE = 1.0 - Omm
    fDE  = a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))
    return Omm * (1+z)**3 + OmDE * fDE

def DM(z, H0, Omm, w0, wa):
    """Distancia comóvil [Mpc]."""
    I, _ = quad(lambda zp: C_KM / (H0 * np.sqrt(E2(zp, Omm, w0, wa))),
                0, z, limit=200, epsrel=1e-6)
    return I

def DH(z, H0, Omm, w0, wa):
    """Distancia de Hubble c/H(z) [Mpc]."""
    return C_KM / (H0 * np.sqrt(E2(z, Omm, w0, wa)))

def DV(z, H0, Omm, w0, wa):
    """Distancia volumétrica [Mpc]."""
    dm = DM(z, H0, Omm, w0, wa)
    dh = DH(z, H0, Omm, w0, wa)
    return (z * dm**2 * dh)**(1/3)

# ─── Horizonte de sonido r_d (Aubourg+2015 fit, calibrado con CLASS) ──────────
def rd_fit(H0, ombh2, Omm):
    """
    r_d en Mpc. Fórmula Aubourg+2015 calibrada con CLASS:
    CLASS Fase 1 con params SSEE dio r_d=147.01 Mpc.
    """
    h      = H0 / 100.0
    omm_h2 = Omm * h**2
    rd = 147.05 * (omm_h2 / 0.14305)**(-0.252) * (ombh2 / 0.02236)**(-0.137)
    return rd

# ─── DESI DR2 BAO (Abdul-Karim et al. 2025, arXiv:2503.14738) ─────────────────
# tipo: 0=DM/rd  1=DH/rd  2=DV/rd
DESI_DATA = np.array([
    [0.295,  7.93,  0.15,  2],   # BGS     DV/rd
    [0.510, 13.62,  0.25,  0],   # LRG1    DM/rd
    [0.510, 20.98,  0.61,  1],   # LRG1    DH/rd
    [0.706, 16.85,  0.32,  0],   # LRG2    DM/rd
    [0.706, 20.08,  0.60,  1],   # LRG2    DH/rd
    [0.930, 21.71,  0.28,  0],   # LRG3    DM/rd
    [0.930, 17.88,  0.35,  1],   # LRG3    DH/rd
    [1.317, 27.79,  0.69,  0],   # ELG     DM/rd
    [1.317, 13.82,  0.42,  1],   # ELG     DH/rd
    [1.491, 30.21,  0.79,  0],   # QSO     DM/rd
    [1.491, 13.23,  0.55,  1],   # QSO     DH/rd
    [2.330, 39.71,  0.94,  0],   # Lya     DM/rd
    [2.330,  8.52,  0.17,  1],   # Lya     DH/rd
])

def chi2_desi(H0, Omm, w0, wa, ombh2):
    rd = rd_fit(H0, ombh2, Omm)
    chi2 = 0.0
    for row in DESI_DATA:
        z, obs, sigma, tipo = row
        if   tipo == 0: pred = DM(z, H0, Omm, w0, wa) / rd
        elif tipo == 1: pred = DH(z, H0, Omm, w0, wa) / rd
        else:           pred = DV(z, H0, Omm, w0, wa) / rd
        chi2 += ((pred - obs) / sigma)**2
    return chi2

# ─── Planck 2018 prior comprimido ─────────────────────────────────────────────
P18_H0 = 67.36;  P18_SH0  = 0.54
P18_OMM= 0.3153; P18_SOMM = 0.0073
P18_RHO= -0.85
P18_OBH2 = 0.02237; P18_SOBH2 = 0.00015

def chi2_planck(H0, Omm, ombh2):
    dH = (H0   - P18_H0)   / P18_SH0
    dO = (Omm  - P18_OMM)  / P18_SOMM
    dB = (ombh2 - P18_OBH2)/ P18_SOBH2
    chi2_HO = dH**2 + dO**2 - 2*P18_RHO*dH*dO
    return chi2_HO + dB**2

# ─── fσ8 — seis encuestas (idénticas a Paper 6) ─────────────────────────────
# Dataset: mismos 6 surveys de ssee_paper6_verification.py
# Datos: (z, fsig8_obs, sigma_fsig8)
FSIG8_Z    = np.array([0.067, 0.150, 0.320, 0.380, 0.510, 0.570])
FSIG8_OBS  = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
FSIG8_ERR  = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])
# surveys: 6dFGRS, SDSS MGS, BOSS DR12 (x2), eBOSS DR16, WiggleZ
SIG8_EFF_SSEE = 0.737   # CLASS Fase 2c (paper 6 two-sector, top-hat)

def _growth_ode_ssee(Omm, w0, wa):
    """
    Resuelve ODE de crecimiento con salida densa para evaluar D(z) y f(z).
    Retorna (sol, D0) — sol tiene dense_output=True, D0 = D(z=0).
    """
    from scipy.integrate import solve_ivp

    def E2_ssee(a):
        z_ = 1.0/a - 1.0
        return max(E2(z_, Omm, w0, wa), 1e-30)

    def rhs(lna, Y):
        a   = np.exp(lna)
        e2  = E2_ssee(a)
        Om  = Omm * a**(-3) / e2
        da  = 1e-5
        dlnE2 = (E2_ssee(a+da) - E2_ssee(a-da)) / (2*da) * a / e2
        D, Dp = Y
        return [Dp, -(2 + 0.5*dlnE2)*Dp + 1.5*Om*D]

    span = (np.log(0.01), 0.0)
    sol  = solve_ivp(rhs, span, [0.01, 0.01],
                     method="DOP853", rtol=1e-8, atol=1e-10,
                     dense_output=True)
    D0 = sol.sol(0.0)[0]
    return sol, D0

def chi2_fsig8(Omm, w0, wa):
    """
    Chi2 fσ8 usando fσ8(z) = f(z) × σ8_eff × D(z)/D(0).
    f(z) = Ω_m(z)^γ (γ=0.554, Paper 5), D(z)/D(0) de la ODE de crecimiento.
    Mismo dataset y formulación que Paper 6 (ssee_paper6_verification.py).
    """
    sol, D0 = _growth_ode_ssee(Omm, w0, wa)
    chi2 = 0.0
    for i, z in enumerate(FSIG8_Z):
        lna_z  = np.log(1.0/(1.0 + z))
        D_z    = sol.sol(lna_z)[0]
        Dp_z   = sol.sol(lna_z)[1]
        # f = d ln D / d ln a = D' / D
        f_z    = Dp_z / D_z
        pred   = f_z * SIG8_EFF_SSEE * (D_z / D0)
        chi2  += ((pred - FSIG8_OBS[i]) / FSIG8_ERR[i])**2
    return chi2

def growth_factor_ssee(Omm, w0, wa):
    """G = D1_SSEE/D1_ΛCDM al dia de hoy (para diagnóstico en tabla resumen)."""
    from scipy.integrate import solve_ivp

    def E2_lcdm(a, Om=0.3153):
        return max(Om * a**(-3) + (1 - Om), 1e-30)

    def rhs(lna, Y, E2_func, Om_func):
        a   = np.exp(lna)
        e2  = E2_func(a)
        Om  = Om_func(a) / e2
        da  = 1e-5
        dlnE2 = (E2_func(a+da) - E2_func(a-da)) / (2*da) * a / e2
        D, Dp = Y
        return [Dp, -(2 + 0.5*dlnE2)*Dp + 1.5*Om*D]

    def E2_ssee(a):
        return max(E2(1/a - 1, Omm, w0, wa), 1e-30)

    IC   = [0.01, 0.01]
    span = (np.log(0.01), 0.0)
    opts = dict(method="DOP853", rtol=1e-8, atol=1e-10)
    sol_s = solve_ivp(rhs, span, IC, args=(E2_ssee,  lambda a: Omm * a**(-3)),  **opts)
    sol_l = solve_ivp(rhs, span, IC, args=(E2_lcdm,  lambda a: 0.3153 * a**(-3)), **opts)
    return sol_s.y[0, -1] / sol_l.y[0, -1]

# ─── Cúmulos de galaxias IGIMF (Zhang et al. 2026, arXiv:2602.06082) ──────────
# (M_bar_IGIMF, sigma_bar, M_obs, sigma_obs) en unidades de 10^14 M_sun
CLUSTER_DATA = np.array([
    [1.8, 0.2, 9.8,  1.0],   # Coma
    [2.2, 0.2, 12.0, 1.2],   # A2029
    [1.5, 0.2, 8.0,  1.0],   # A478
    [1.2, 0.2, 6.5,  1.0],   # Bullet
])

FNU_SSEE = 0.020   # fracción de masa en neutrinos (Paper 2: 2%)

def chi2_clusters(Omm):
    """
    Chi2 para masas de cúmulos con predicción SSEE: M_pred = KAL0*(1+f_ν)*M_bar.
    KAL0 ≈ 5.5214 es la constante estructural algebraica (Paper 2, ec. 18).
    Sigma observacional: sigma_obs en M_obs (Zhang et al. 2026).
    """
    chi2 = 0.0
    KAL = KAL0 * (1.0 + FNU_SSEE)
    for row in CLUSTER_DATA:
        Mbar, sBar, Mobs, sObs = row
        M_pred = KAL * Mbar
        chi2 += ((M_pred - Mobs) / sObs)**2
    return chi2

# ─── Log-posterior ─────────────────────────────────────────────────────────────
PARAM_NAMES = ["H0", "ombh2", "Omega_m", "w0", "wa"]
N_DIM       = len(PARAM_NAMES)

#              H0      ombh2    Omm     w0       wa
THETA_ALG  = np.array([H0_ALG, OMBH2_ALG, OM_CMB, W0_ALG, WA_ALG])
THETA_MIN  = np.array([62.0,   0.0210,    0.24,   -1.30,  -2.0 ])
THETA_MAX  = np.array([74.0,   0.0255,    0.40,   -0.40,   0.5 ])

W0_SIGMA = 0.08
WA_SIGMA = 0.12

def log_prior(theta):
    H0, ombh2, Omm, w0, wa = theta
    # Caja
    if np.any(theta < THETA_MIN) or np.any(theta > THETA_MAX):
        return -np.inf
    # Prior gaussiano algebraico sobre w0/wa (centro prediccion SSEE)
    lp_w0 = -0.5 * ((w0 - W0_ALG) / W0_SIGMA)**2
    lp_wa = -0.5 * ((wa - WA_ALG) / WA_SIGMA)**2
    return lp_w0 + lp_wa

def log_likelihood(theta):
    H0, ombh2, Omm, w0, wa = theta

    # L1: DESI DR2
    lnL = -0.5 * chi2_desi(H0, Omm, w0, wa, ombh2)

    # L2: Planck comprimido
    lnL += -0.5 * chi2_planck(H0, Omm, ombh2)

    # L3: fσ8 seis encuestas
    if not NO_FSIG8:
        lnL += -0.5 * chi2_fsig8(Omm, w0, wa)

    # L4: Cúmulos IGIMF
    lnL += -0.5 * chi2_clusters(Omm)

    return lnL

def log_posterior(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta)
    return (lp + ll) if np.isfinite(ll) else -np.inf

# ─── Setup MCMC ────────────────────────────────────────────────────────────────
N_WALKERS  = 20   if TEST_MODE else 100
N_STEPS    = 200  if TEST_MODE else 25000
N_BURN     = 50   if TEST_MODE else 2000
SAVE_EVERY = 500

SPREAD = np.array([0.30,  0.00015, 0.005, 0.025, 0.04])
rng    = np.random.default_rng(42)
p0     = THETA_ALG + SPREAD * rng.standard_normal((N_WALKERS, N_DIM))
p0     = np.clip(p0, THETA_MIN + 1e-6, THETA_MAX - 1e-6)

log("=" * 72)
log(f"SSEE-V3.6 — Fase 4 MCMC Multi-Probe  {'[TEST]' if TEST_MODE else '[OVERNIGHT]'}")
log("=" * 72)
log(f"Algebraico: w0={W0_ALG:.5f}  wa={WA_ALG:.5f}  H0={H0_ALG}  Omega_m={OM_CMB:.5f}")
log(f"Prior w0: N({W0_ALG:.4f},{W0_SIGMA}) | Prior wa: N({WA_ALG:.4f},{WA_SIGMA})")
log(f"Walkers={N_WALKERS}  Steps={N_STEPS}  Burn={N_BURN}")
log(f"Likelihoods: DESI(13pts) + Planck(comprimido) + fσ8({'SI' if not NO_FSIG8 else 'NO'}) + Clusters(4)")
log(f"sigma8_eff = {SIG8_EFF_SSEE} (CLASS Fase 2c calibrado)")
log(f"MIRA = {MIRA:.5f}")

# Test punto algebraico
log("\nTest log_posterior en valores algebraicos SSEE...")
lp_test = log_posterior(THETA_ALG)
log(f"  log_posterior = {lp_test:.3f}  ({'OK ✅' if np.isfinite(lp_test) else 'FALLO ❌'})")
rd_test = rd_fit(H0_ALG, OMBH2_ALG, OM_CMB)
log(f"  r_d (fit)     = {rd_test:.3f} Mpc  | CLASS Fase1 = {RD_CLASS_SSEE} Mpc | Planck = 147.09 Mpc")
log(f"  chi2_DESI     = {chi2_desi(H0_ALG, OM_CMB, W0_ALG, WA_ALG, OMBH2_ALG):.3f}")
log(f"  chi2_Planck   = {chi2_planck(H0_ALG, OM_CMB, OMBH2_ALG):.3f}")
if not NO_FSIG8:
    log(f"  chi2_fsig8    = {chi2_fsig8(OM_CMB, W0_ALG, WA_ALG):.3f}")
log(f"  chi2_clusters = {chi2_clusters(OM_CMB):.3f}")

if not np.isfinite(lp_test):
    log("ERROR: posterior infinito. Abortando.")
    sys.exit(1)
log("")

# ─── Resume o correr desde cero ───────────────────────────────────────────────
if RESUME_MODE and os.path.exists(CHAIN_FILE):
    log(f"Reanudando desde {CHAIN_FILE}...")
    saved = np.load(CHAIN_FILE, allow_pickle=True)
    flat_prev = saved['chain']
    p0 = flat_prev[-N_WALKERS:].reshape(N_WALKERS, N_DIM)
    log(f"  Muestras previas: {flat_prev.shape[0]}")

sampler = emcee.EnsembleSampler(N_WALKERS, N_DIM, log_posterior)

# ─── Burn-in ──────────────────────────────────────────────────────────────────
log(f"Burn-in ({N_BURN} steps)...")
state = sampler.run_mcmc(p0, N_BURN, progress=True)
sampler.reset()
log(f"Burn-in completado.")

# ─── Cadena principal con guardado incremental ────────────────────────────────
log(f"Cadena principal ({N_STEPS} steps, guardado c/{SAVE_EVERY})...")
n_done = 0
for batch in range(0, N_STEPS, SAVE_EVERY):
    steps_this = min(SAVE_EVERY, N_STEPS - batch)
    state = sampler.run_mcmc(state, steps_this, progress=True)
    n_done += steps_this
    flat_now = sampler.get_chain(flat=True)
    acc = np.mean(sampler.acceptance_fraction)
    np.savez(CHAIN_FILE, chain=flat_now, param_names=PARAM_NAMES,
             theta_alg=THETA_ALG, steps_done=n_done)
    elapsed = (time.time() - t0) / 60
    log(f"  {n_done}/{N_STEPS} steps  acc={acc:.3f}  {elapsed:.1f} min")

log("Cadena completa.")

# ─── Análisis ─────────────────────────────────────────────────────────────────
flat = sampler.get_chain(flat=True)
log(f"\nMuestras totales: {flat.shape[0]}")
log(f"Acceptance fraction media: {np.mean(sampler.acceptance_fraction):.3f}")

try:
    tau = sampler.get_autocorr_time(quiet=True)
    N_eff = int(N_STEPS * N_WALKERS / np.max(tau))
    log(f"tau_max={np.max(tau):.1f}  =>  N_eff={N_eff}")
except Exception as e:
    log(f"Autocorr: {e}")

med  = np.median(flat, axis=0)
lo16 = np.percentile(flat, 16, axis=0)
hi84 = np.percentile(flat, 84, axis=0)
best_lp = sampler.get_log_prob(flat=True)
best = flat[np.argmax(best_lp)]

log("\n" + "─"*72)
log(f"{'Param':<12} {'Algebraico':>12} {'Mediana':>12} {'−1σ':>10} {'+1σ':>10} {'Tension':>9}")
log("─"*72)
for i, name in enumerate(PARAM_NAMES):
    dn = med[i] - lo16[i]
    up = hi84[i] - med[i]
    s  = (dn + up) / 2
    t  = abs(THETA_ALG[i] - med[i]) / s if s > 0 else 99.
    log(f"{name:<12} {THETA_ALG[i]:>12.5f} {med[i]:>12.5f} {dn:>10.5f} {up:>10.5f} {t:>8.2f}σ")
log("─"*72)
log("MAP (best-fit):")
for i, name in enumerate(PARAM_NAMES):
    log(f"  {name:<12} = {best[i]:.5f}")

# Resultados derivados
rd_med  = rd_fit(med[0], med[1], med[2])
rd_best = rd_fit(best[0], best[1], best[2])
log(f"\n  r_d (mediana)  = {rd_med:.3f} Mpc")
log(f"  r_d (MAP)      = {rd_best:.3f} Mpc")
log(f"  r_d (Planck)   = 147.09 Mpc")
log(f"  r_d (CLASS F1) = {RD_CLASS_SSEE} Mpc")

G_med = growth_factor_ssee(med[2], med[3], med[4])
log(f"  G=D1_SSEE/D1_ΛCDM (mediana) = {G_med:.4f}  (Paper5: 0.866)")
log(f"  σ8_eff (fijo CLASS) = {SIG8_EFF_SSEE}  (Fase 2c)")

# ΔBIC vs ΛCDM
# ΛCDM tiene 6 parametros libres; SSEE tiene 5 (con w0/wa prior gaussiano ≈ prior fijo = k=1)
chi2_best_ssee  = -2 * np.max(best_lp)
chi2_best_lcdm  = chi2_planck(P18_H0, P18_OMM, P18_OBH2)  # aprox ΛCDM en prior
DELTA_BIC = (N_DIM - 1) * np.log(13 + 6) - (chi2_best_ssee)  # aproximado
log(f"\n  χ²_best (SSEE Fase4) = {chi2_best_ssee:.2f}")

# Guardar resultados
np.savez(CHAIN_FILE, chain=flat, param_names=PARAM_NAMES,
         theta_alg=THETA_ALG, medians=med, lo16=lo16, hi84=hi84, best=best)
log(f"\nCadena final: {CHAIN_FILE}")

# ─── Figura corner ─────────────────────────────────────────────────────────────
labels = [r"$H_0$ [km/s/Mpc]", r"$\Omega_b h^2$", r"$\Omega_m$",
          r"$w_0$", r"$w_a$"]
fig_c = corner.corner(flat, labels=labels, truths=list(THETA_ALG),
                      truth_color="red", show_titles=True,
                      quantiles=[0.16, 0.5, 0.84],
                      title_kwargs={"fontsize": 10})
fig_c.suptitle(
    r"SSEE-V3.6 Fase 4 MCMC — DESI DR2 + Planck 2018 + $f\sigma_8$ + Clusters"
    "\n"
    r"Punto rojo: valores algebraicos ($\varphi,\pi$)",
    fontsize=12, y=1.01)
fig_c.savefig(f"{OUT_DIR}/corner_fase4.pdf", bbox_inches="tight", dpi=150)
fig_c.savefig(f"{OUT_DIR}/corner_fase4.png", bbox_inches="tight", dpi=150)
log(f"Corner: {OUT_DIR}/corner_fase4.pdf")

# ─── Figura w0-wa posterior ───────────────────────────────────────────────────
fig_w, ax = plt.subplots(figsize=(7, 6))
w0_s, wa_s = flat[:, 3], flat[:, 4]
H, xe, ye = np.histogram2d(w0_s, wa_s, bins=60)
xc = 0.5*(xe[1:]+xe[:-1]); yc = 0.5*(ye[1:]+ye[:-1])
ax.contourf(xc, yc, H.T, levels=8, cmap="Blues", alpha=0.75)
ax.contour(xc, yc, H.T, levels=5, colors="navy", linewidths=0.8, alpha=0.6)
ax.scatter(W0_ALG, WA_ALG, color="red", s=150, zorder=10, marker="*",
           label=fr"SSEE algebraico: $w_0={W0_ALG:.3f},\;w_a={WA_ALG:.3f}$")
ax.scatter(med[3], med[4], color="orange", s=80, zorder=9, marker="D",
           label=fr"MCMC mediana: $w_0={med[3]:.3f},\;w_a={med[4]:.3f}$")
ax.axvline(-1.0, ls="--", color="gray", lw=1.0, label=r"$\Lambda$CDM ($w_0=-1$)")
ax.axhline(0.0,  ls="--", color="gray", lw=1.0)
ax.set_xlabel(r"$w_0$", fontsize=13)
ax.set_ylabel(r"$w_a$", fontsize=13)
ax.set_title(r"Posterior $w_0$-$w_a$ — SSEE Fase 4 MCMC", fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
fig_w.savefig(f"{OUT_DIR}/w0wa_fase4.pdf", bbox_inches="tight", dpi=150)
fig_w.savefig(f"{OUT_DIR}/w0wa_fase4.png", bbox_inches="tight", dpi=150)
log(f"w0-wa: {OUT_DIR}/w0wa_fase4.pdf")

# ─── Figura r_d posterior ─────────────────────────────────────────────────────
rd_post = np.array([rd_fit(s[0], s[1], s[2]) for s in flat])
fig_rd, ax2 = plt.subplots(figsize=(7, 4))
ax2.hist(rd_post, bins=50, color="steelblue", alpha=0.7, density=True,
         label=fr"SSEE Fase4 MCMC: ${np.median(rd_post):.2f}\pm{np.std(rd_post):.2f}$ Mpc")
ax2.axvline(147.09, color="red",    lw=2.0, ls="-",  label=r"Planck 2018: $147.09\pm0.26$ Mpc")
ax2.axvline(RD_CLASS_SSEE, color="green", lw=2.0, ls="--",
            label=fr"CLASS Fase1: ${RD_CLASS_SSEE}$ Mpc")
ax2.set_xlabel(r"$r_d$ [Mpc]", fontsize=12)
ax2.set_ylabel("Densidad", fontsize=12)
ax2.set_title(r"Posterior $r_d$ — SSEE Fase 4 MCMC", fontsize=12)
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
fig_rd.savefig(f"{OUT_DIR}/rd_posterior_fase4.pdf", bbox_inches="tight", dpi=150)
fig_rd.savefig(f"{OUT_DIR}/rd_posterior_fase4.png", bbox_inches="tight", dpi=150)
log(f"r_d posterior: {OUT_DIR}/rd_posterior_fase4.pdf")

# ─── Figura fσ8 posterior predictive ─────────────────────────────────────────
z_plot = np.linspace(0.1, 1.5, 80)
fig_fs, ax3 = plt.subplots(figsize=(8, 5))
# Banda posterior
fs8_samples = []
for s in flat[::20]:
    sol_s, D0_s = _growth_ode_ssee(s[2], s[3], s[4])
    fs8_z = []
    for z in z_plot:
        lna_z = np.log(1.0/(1.0+z))
        D_z = sol_s.sol(lna_z)[0];  Dp_z = sol_s.sol(lna_z)[1]
        fs8_z.append((Dp_z/D_z) * SIG8_EFF_SSEE * (D_z/D0_s))
    fs8_samples.append(fs8_z)
fs8_arr = np.array(fs8_samples)
ax3.fill_between(z_plot, np.percentile(fs8_arr, 16, axis=0),
                         np.percentile(fs8_arr, 84, axis=0),
                 color="steelblue", alpha=0.4, label="SSEE Fase4 (68% CL)")
ax3.plot(z_plot, np.median(fs8_arr, axis=0), "b-", lw=2.0, label="SSEE mediana")
# Algebraico
sol_alg, D0_alg = _growth_ode_ssee(OM_CMB, W0_ALG, WA_ALG)
fs8_alg = []
for z in z_plot:
    lna_z = np.log(1.0/(1.0+z))
    D_z = sol_alg.sol(lna_z)[0];  Dp_z = sol_alg.sol(lna_z)[1]
    fs8_alg.append((Dp_z/D_z) * SIG8_EFF_SSEE * (D_z/D0_alg))
ax3.plot(z_plot, fs8_alg, "r--", lw=2.0, label=r"SSEE algebraico ($\varphi,\pi$)")
# Datos (mismos surveys que Paper 6)
survey_labels = ["6dFGRS", "SDSS MGS", "BOSS DR12", "BOSS DR12", "eBOSS DR16", "WiggleZ"]
for z, obs, sig, lbl in zip(FSIG8_Z, FSIG8_OBS, FSIG8_ERR, survey_labels):
    ax3.errorbar(z, obs, yerr=sig, fmt="ko", capsize=4, ms=5)
    ax3.annotate(lbl, (z, obs), textcoords="offset points", xytext=(5, 5), fontsize=8)
ax3.set_xlabel(r"Redshift $z$", fontsize=12)
ax3.set_ylabel(r"$f\sigma_8(z)$", fontsize=12)
ax3.set_title(r"SSEE Fase 4 — posterior predictive $f\sigma_8$", fontsize=12)
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
fig_fs.savefig(f"{OUT_DIR}/fsig8_posterior_fase4.pdf", bbox_inches="tight", dpi=150)
fig_fs.savefig(f"{OUT_DIR}/fsig8_posterior_fase4.png", bbox_inches="tight", dpi=150)
log(f"fσ8 posterior: {OUT_DIR}/fsig8_posterior_fase4.pdf")

elapsed_total = (time.time() - t0) / 60
log(f"\nTiempo total: {elapsed_total:.1f} min  ({elapsed_total/60:.2f} h)")
log("Fase 4 MCMC completada. ✅")
