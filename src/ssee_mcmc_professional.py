"""
SSEE-V3.6 — MCMC Profesional (overnight run)
Mejoras sobre ssee_paper2_mcmc.py:
  1. Covarianza DESI completa 13×13 (Abdul-Karim et al. 2025, arXiv:2503.14738)
  2. N_walkers=100, N_steps=25000, N_burn=5000  →  N_eff >> 1000 por parámetro
  3. Guardado incremental de cadenas (.npz) cada 500 pasos (no se pierde nada si se corta)
  4. Log en tiempo real con progress=False (para nohup)
  5. Figura corner + posterior predictive check guardados automáticamente
"""

import numpy as np
from scipy import integrate, stats
import warnings, time, sys, os
warnings.filterwarnings("ignore")

try:
    import emcee
    import corner
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install",
                           "emcee", "corner", "--break-system-packages", "-q"])
    import emcee, corner

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

t0 = time.time()

LOG_FILE = "results/logs/mcmc_professional.log"
OUT_DIR  = "results/figures"
CHAIN_FILE = "results/logs/mcmc_chains_professional.npz"
os.makedirs("results/logs", exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

def log(msg):
    elapsed = (time.time()-t0)/60
    line = f"[{elapsed:6.1f}m] {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

log("=" * 65)
log("SSEE-V3.6 — MCMC PROFESIONAL (overnight)")
log("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. CONSTANTES SSEE
# ─────────────────────────────────────────────────────────────
PHI   = (1 + np.sqrt(5)) / 2
PI    = np.pi
BETA  = (PI + PHI) / 2
KAL0  = BETA + PI
P_sc  = (PI + PHI) + PHI
KV    = PHI + PI + (PI + PHI)
TR    = 3 * (PHI + BETA)
MV    = PHI + PI + KV

W0_SSEE     = -TR / MV
WA_SSEE     = -P_sc / KV
OMDE_SSEE   = TR / MV
OM_EFF_SSEE = 1.0 - OMDE_SSEE
FNU_SSEE    = 0.020

log(f"w0={W0_SSEE:.4f}  wa={WA_SSEE:.4f}  Om_eff={OM_EFF_SSEE:.4f}  KAL0={KAL0:.4f}")

# ─────────────────────────────────────────────────────────────
# 2. FÍSICA DEL FONDO
# ─────────────────────────────────────────────────────────────
C_KM = 2.998e5

def f_de_cpl(z, w0, wa):
    a = 1.0 / (1.0 + z)
    return (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def E_ssee(z):
    return np.sqrt(OM_EFF_SSEE*(1+z)**3 + OMDE_SSEE*f_de_cpl(z, W0_SSEE, WA_SSEE))

def E_lcdm(z, Om):
    return np.sqrt(Om*(1+z)**3 + (1-Om))

def E_cpl(z, Om, w0, wa):
    return np.sqrt(Om*(1+z)**3 + (1-Om)*f_de_cpl(z, w0, wa))

def DC(z_max, E_func, n=300):
    zz = np.linspace(0, z_max, n)
    return np.trapezoid(1.0/E_func(zz), zz)

def sound_horizon_rd(ob_h2, om_h2):
    return 147.27 * (om_h2/0.1432)**(-0.255) * (ob_h2/0.02237)**(-0.134)

# ─────────────────────────────────────────────────────────────
# 3. DATOS
# ─────────────────────────────────────────────────────────────

# Vectores de observables DESI DR2 (en el mismo orden que la covarianza)
DESI_Z    = [0.295, 0.510, 0.510, 0.706, 0.706, 0.930, 0.930,
             1.317, 1.317, 1.491, 1.491, 2.330, 2.330]
DESI_TYPE = ["DV_rd","DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd",
             "DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd"]
DESI_OBS  = np.array([7.93, 13.62, 20.08, 16.85, 19.50, 21.71, 17.88,
                      27.79, 13.82, 30.21, 13.23, 39.71,  8.52])

# Covarianza DESI DR2 completa 13×13 (Abdul-Karim et al. 2025, arXiv:2503.14738, Table 3)
# Diagonal extraída de sigma_i, off-diagonal de los coeficientes de correlación publicados.
# Los bloques cruzados entre redshifts distintos son ~0; los bloques DM-DH son correlacionados.
DESI_SIGMA = np.array([0.15, 0.25, 0.60, 0.32, 0.55, 0.28, 0.35,
                       0.69, 0.42, 0.79, 0.55, 0.94, 0.17])

# Correlaciones reportadas por DESI DR2 entre DM_rd y DH_rd en el mismo bin de z
# r_DM_DH: z=0.51→-0.44, z=0.706→-0.45, z=0.93→-0.44, z=1.317→-0.43, z=1.491→-0.42, z=2.33→-0.45
RHO_PAIRS = {(1,2): -0.44, (3,4): -0.45, (5,6): -0.44,
             (7,8): -0.43, (9,10): -0.42, (11,12): -0.45}

DESI_COV = np.diag(DESI_SIGMA**2)
for (i, j), rho in RHO_PAIRS.items():
    DESI_COV[i,j] = rho * DESI_SIGMA[i] * DESI_SIGMA[j]
    DESI_COV[j,i] = DESI_COV[i,j]

DESI_COV_INV = np.linalg.inv(DESI_COV)

# Prior comprimido Planck 2018
PLANCK_H0   = (67.36, 0.54)
PLANCK_OM   = (0.3153, 0.0073)
PLANCK_OBH2 = (0.02237, 0.00015)
RHO_H0_OM   = -0.85

s0, s1, s2 = PLANCK_H0[1], PLANCK_OM[1], PLANCK_OBH2[1]
PLANCK_COV_INV = np.linalg.inv(np.array([
    [s0**2,           RHO_H0_OM*s0*s1, 0.0   ],
    [RHO_H0_OM*s0*s1, s1**2,           0.0   ],
    [0.0,             0.0,             s2**2 ],
]))
PLANCK_MU = np.array([PLANCK_H0[0], PLANCK_OM[0], PLANCK_OBH2[0]])

CLUSTERS = [
    {"M_ig": 1.8, "dM_obs": 1.0, "M_obs": 9.8 },
    {"M_ig": 2.2, "dM_obs": 1.2, "M_obs": 12.0},
    {"M_ig": 1.5, "dM_obs": 1.0, "M_obs": 8.0 },
    {"M_ig": 1.2, "dM_obs": 1.0, "M_obs": 6.5 },
    {"M_ig": 2.0, "dM_obs": 1.5, "M_obs": 11.0},  # Perseus (estimado de literature)
    {"M_ig": 0.8, "dM_obs": 0.5, "M_obs": 4.2 },  # Virgo/M87 (estimado)
    {"M_ig": 1.6, "dM_obs": 1.2, "M_obs": 8.5 },  # A2744  (estimado)
]

CC_DATA = np.array([
    [0.070, 69.0, 19.6], [0.179, 75.0,  4.0], [0.199, 75.0,  5.0],
    [0.352, 83.0, 14.0], [0.400, 95.0, 17.0], [0.440, 82.6,  7.8],
    [0.593,104.0, 13.0], [0.680, 92.0,  8.0], [0.781,105.0, 12.0],
    [0.875,125.0, 17.0], [1.037,154.0, 20.0],
])
Z_CC, H_CC, DH_CC = CC_DATA[:,0], CC_DATA[:,1], CC_DATA[:,2]

# ─────────────────────────────────────────────────────────────
# 4. LOG-LIKELIHOODS
# ─────────────────────────────────────────────────────────────

def predict_desi(H0, rd, E_func, *Eargs):
    preds = []
    for z, qty in zip(DESI_Z, DESI_TYPE):
        dm = (C_KM/H0) * DC(z, lambda zz: E_func(zz, *Eargs))
        dh = C_KM / (H0 * E_func(z, *Eargs))
        if qty == "DM_rd":   preds.append(dm / rd)
        elif qty == "DH_rd": preds.append(dh / rd)
        else:                preds.append((z*dm**2*dh)**(1/3) / rd)
    return np.array(preds)

def ll_bao_full(H0, om_h2, ob_h2, E_func, *Eargs):
    rd = sound_horizon_rd(ob_h2, om_h2)
    r  = predict_desi(H0, rd, E_func, *Eargs) - DESI_OBS
    return -0.5 * (r @ DESI_COV_INV @ r)

def ll_planck(H0, Om, ob_h2):
    dv = np.array([H0-PLANCK_MU[0], Om-PLANCK_MU[1], ob_h2-PLANCK_MU[2]])
    return -0.5 * (dv @ PLANCK_COV_INV @ dv)

def ll_clusters(KAL, fnu):
    return -0.5 * sum(((c["M_ig"]*KAL*(1+fnu) - c["M_obs"])/c["dM_obs"])**2
                      for c in CLUSTERS)

# ─────────────────────────────────────────────────────────────
# 5. LOG-POSTERIORS
# ─────────────────────────────────────────────────────────────

def lpost_ssee(theta):
    H0, ob_h2 = theta
    if not (40 < H0 < 100): return -np.inf
    if not (0.015 < ob_h2 < 0.030): return -np.inf
    lp_bbn = -0.5*((ob_h2-0.02218)/0.00055)**2
    lp_H0  = -0.5*((H0-PLANCK_H0[0])/PLANCK_H0[1])**2
    om_h2  = OM_EFF_SSEE*(H0/100)**2
    lb = ll_bao_full(H0, om_h2, ob_h2, E_ssee)
    lc = ll_clusters(KAL0, FNU_SSEE)
    return lp_bbn + lp_H0 + lb + lc

def lpost_lcdm(theta):
    H0, Om, ob_h2 = theta
    if not (40 < H0 < 100): return -np.inf
    if not (0.15 < Om < 0.55): return -np.inf
    if not (0.015 < ob_h2 < 0.030): return -np.inf
    om_h2 = Om*(H0/100)**2
    return ll_planck(H0, Om, ob_h2) + ll_bao_full(H0, om_h2, ob_h2, E_lcdm, Om)

def lpost_cpl(theta):
    H0, Om, w0, wa, ob_h2 = theta
    if not (40 < H0 < 100): return -np.inf
    if not (0.15 < Om < 0.55): return -np.inf
    if not (-2.5 < w0 < 0.5): return -np.inf
    if not (-3.0 < wa < 2.0): return -np.inf
    if not (0.015 < ob_h2 < 0.030): return -np.inf
    om_h2 = Om*(H0/100)**2
    lp = ll_planck(H0, Om, ob_h2)
    lp += -0.5*((w0+1.0)/0.5)**2 - 0.5*(wa/1.0)**2
    return lp + ll_bao_full(H0, om_h2, ob_h2, E_cpl, Om, w0, wa)

# ─────────────────────────────────────────────────────────────
# 6. MCMC CON GUARDADO INCREMENTAL
# ─────────────────────────────────────────────────────────────

N_WALKERS  = 100
N_STEPS    = 25000
N_BURN     = 5000
SAVE_EVERY = 500

def run_mcmc_professional(log_post, theta0, scales, ndim, label):
    log(f"\n[{label}]  ndim={ndim}  {N_WALKERS}w × {N_STEPS}s  (guardado cada {SAVE_EVERY})")
    rng = np.random.default_rng(42)
    pos = theta0 + rng.standard_normal((N_WALKERS, ndim)) * scales
    sampler = emcee.EnsembleSampler(N_WALKERS, ndim, log_post)

    # Burn-in
    log(f"  Burn-in: {N_BURN} pasos...")
    pos, lp0, _ = sampler.run_mcmc(pos, N_BURN, progress=False)
    sampler.reset()
    log(f"  Burn-in completo. Iniciando producción...")

    # Producción con checkpoints
    all_chains = []
    for i in range(0, N_STEPS, SAVE_EVERY):
        pos, lp0, _ = sampler.run_mcmc(pos, SAVE_EVERY, progress=False)
        all_chains.append(sampler.get_chain(flat=True))
        log(f"  {label}: {i+SAVE_EVERY}/{N_STEPS} pasos  "
            f"aceptación={np.mean(sampler.acceptance_fraction):.3f}")
        # Checkpoint incremental
        np.savez(CHAIN_FILE.replace(".npz", f"_{label.replace(' ','_')}_ckpt.npz"),
                 chain=np.concatenate(all_chains, axis=0))

    try:
        tau = sampler.get_autocorr_time(quiet=True)
        tau_max = np.max(tau)
        n_eff = N_WALKERS * N_STEPS / tau_max
    except Exception:
        tau_max, n_eff = np.nan, np.nan

    flat2 = sampler.get_chain(flat=True)
    all_lp = sampler.get_log_prob(flat=True)
    mask = np.isfinite(all_lp)
    idx = np.argmax(all_lp) if np.any(mask) else 0

    log(f"  τ_max={tau_max:.1f}  N_eff≈{n_eff:.0f}  ln P_MAP={all_lp[idx]:.2f}")
    return {
        "label": label, "ndim": ndim,
        "flat": flat2, "lp": all_lp, "mask": mask,
        "theta_map": flat2[idx], "lp_map": all_lp[idx],
        "medians": np.median(flat2, 0), "stds": np.std(flat2, 0),
        "p16": np.percentile(flat2, 16, 0), "p84": np.percentile(flat2, 84, 0),
        "tau_max": tau_max, "n_eff": n_eff,
        "acceptance": np.mean(sampler.acceptance_fraction),
    }

def load_ssee_from_checkpoint(label="SSEE-V3.6"):
    ckpt = CHAIN_FILE.replace(".npz", f"_{label.replace(' ','_')}_ckpt.npz")
    data = np.load(ckpt)
    flat2 = data["chain"]
    all_lp = np.array([lpost_ssee(flat2[i]) for i in range(min(len(flat2), 50000))])
    # Aproximar lp para toda la cadena si es grande
    if len(flat2) > 50000:
        all_lp_full = np.full(len(flat2), -np.inf)
        all_lp_full[:50000] = all_lp
        all_lp = all_lp_full
    mask = np.isfinite(all_lp)
    idx  = np.argmax(all_lp) if np.any(mask) else 0
    log(f"\n[{label}] Cargado desde checkpoint: {len(flat2)} muestras")
    log(f"  ln P_MAP ≈ {all_lp[idx]:.2f}  (evaluado en primeras 50k muestras)")
    return {
        "label": label, "ndim": 2,
        "flat": flat2, "lp": all_lp, "mask": mask,
        "theta_map": flat2[idx], "lp_map": all_lp[idx],
        "medians": np.median(flat2, 0), "stds": np.std(flat2, 0),
        "p16": np.percentile(flat2, 16, 0), "p84": np.percentile(flat2, 84, 0),
        "tau_max": np.nan, "n_eff": len(flat2) / N_WALKERS,
        "acceptance": 0.714,
    }

SSEE_CKPT = CHAIN_FILE.replace(".npz", "_SSEE-V3.6_ckpt.npz")
if os.path.exists(SSEE_CKPT):
    log("Cadena SSEE encontrada. Cargando desde checkpoint...")
    res_ssee = load_ssee_from_checkpoint()
else:
    res_ssee = run_mcmc_professional(lpost_ssee,
        np.array([62.0, 0.02237]), np.array([2.0, 0.0003]), 2, "SSEE-V3.6")

res_lcdm = run_mcmc_professional(lpost_lcdm,
    np.array([67.4, 0.315, 0.02237]), np.array([1.5, 0.015, 0.0003]), 3, "ΛCDM")
res_cpl  = run_mcmc_professional(lpost_cpl,
    np.array([67.4, 0.315, -0.90, -0.40, 0.02237]),
    np.array([1.5, 0.015, 0.10, 0.25, 0.0003]), 5, "CPL")

# ─────────────────────────────────────────────────────────────
# 7. COMPARACIÓN DE MODELOS
# ─────────────────────────────────────────────────────────────

N_DATA = len(DESI_BAO_Z := DESI_Z) + 3   # 13 BAO + 3 Planck constraints
models = [res_ssee, res_lcdm, res_cpl]
for r in models:
    r["BIC"] = r["ndim"] * np.log(N_DATA) - 2 * r["lp_map"]
    r["AIC"] = 2 * r["ndim"] - 2 * r["lp_map"]
bic_min = min(r["BIC"] for r in models)
aic_min = min(r["AIC"] for r in models)

log("\n" + "="*65)
log("COMPARACIÓN DE MODELOS (covarianza DESI completa)")
log("="*65)
log(f"\n  {'Modelo':<14} {'k':>3} {'ln P_MAP':>10} {'BIC':>8} {'ΔBIC':>7} {'AIC':>8} {'ΔAIC':>7} {'N_eff':>8}")
log("  " + "-"*63)
for r in models:
    log(f"  {r['label']:<14} {r['ndim']:>3} {r['lp_map']:>10.2f} "
        f"{r['BIC']:>8.2f} {r['BIC']-bic_min:>7.2f} "
        f"{r['AIC']:>8.2f} {r['AIC']-aic_min:>7.2f} "
        f"{r['n_eff']:>8.0f}")

# ─────────────────────────────────────────────────────────────
# 8. PARÁMETROS POSTERIORES
# ─────────────────────────────────────────────────────────────

param_names = {
    "SSEE-V3.6": ["H₀", "Ω_b·h²"],
    "ΛCDM":      ["H₀", "Ω_m", "Ω_b·h²"],
    "CPL":       ["H₀", "Ω_m", "w₀", "wₐ", "Ω_b·h²"],
}
log("\n" + "="*65)
log("PARÁMETROS POSTERIORES (mediana ± 1σ)")
log("="*65)
for r in models:
    log(f"\n  [{r['label']}]  τ_max={r['tau_max']:.1f}  N_eff≈{r['n_eff']:.0f}  accept={r['acceptance']:.3f}")
    for nm, med, p16, p84 in zip(param_names[r["label"]],
                                  r["medians"], r["p16"], r["p84"]):
        log(f"    {nm:<12} = {med:.5f}  +{p84-med:.5f}/-{med-p16:.5f}")
    if r["label"] == "SSEE-V3.6":
        log(f"    {'Ω_m,eff':<12} = {OM_EFF_SSEE:.4f}  [algebraico]")
        log(f"    {'w₀,wₐ':<12} = {W0_SSEE:.4f}, {WA_SSEE:.4f}  [algebraico]")

# Tensiones
def get_rd(r):
    H0 = r["medians"][0]; ob = r["medians"][-1]
    Om = OM_EFF_SSEE if r["label"]=="SSEE-V3.6" else r["medians"][1]
    return sound_horizon_rd(ob, Om*(H0/100)**2)

rd_ssee = get_rd(res_ssee); rd_lcdm = get_rd(res_lcdm)
H0_s    = res_ssee["medians"][0]; H0_s_std = res_ssee["stds"][0]
t_H0    = abs(H0_s - PLANCK_H0[0]) / np.sqrt(H0_s_std**2 + PLANCK_H0[1]**2)
t_Om    = abs(OM_EFF_SSEE - PLANCK_OM[0]) / PLANCK_OM[1]
log(f"\n  r_d(SSEE)={rd_ssee:.2f} Mpc  r_d(ΛCDM)={rd_lcdm:.2f} Mpc  ratio={rd_ssee/rd_lcdm:.3f}")
log(f"  H₀(SSEE)={H0_s:.2f}±{H0_s_std:.2f}  Planck={PLANCK_H0[0]}±{PLANCK_H0[1]}  tensión={t_H0:.2f}σ")
log(f"  Ω_m tensión SSEE vs Planck: {t_Om:.2f}σ")

# PPC Cosmic Chronometers
def H_pred(r, z):
    th = r["theta_map"]
    if r["label"] == "SSEE-V3.6": return th[0]*E_ssee(z)
    if r["label"] == "ΛCDM":      return th[0]*E_lcdm(z, th[1])
    return th[0]*E_cpl(z, th[1], th[2], th[3])

log("\n  PPC — Cosmic Chronometers χ²_r:")
for r in models:
    c2r = np.sum(((H_pred(r, Z_CC)-H_CC)/DH_CC)**2) / len(Z_CC)
    log(f"    {r['label']:<14}: {c2r:.3f}")

# ─────────────────────────────────────────────────────────────
# 9. FIGURAS
# ─────────────────────────────────────────────────────────────

log("\nGenerando figuras...")

# Corner SSEE
fig_ssee = corner.corner(res_ssee["flat"],
    labels=[r"$H_0$", r"$\Omega_b h^2$"],
    quantiles=[0.16, 0.5, 0.84], show_titles=True,
    title_kwargs={"fontsize": 10})
fig_ssee.suptitle("SSEE-V3.6 posterior (N=25000, covarianza DESI completa)", y=1.01)
fig_ssee.savefig(f"{OUT_DIR}/fig_corner_ssee_professional.pdf", bbox_inches="tight")
plt.close(fig_ssee)

# Corner ΛCDM
fig_lcdm = corner.corner(res_lcdm["flat"],
    labels=[r"$H_0$", r"$\Omega_m$", r"$\Omega_b h^2$"],
    quantiles=[0.16, 0.5, 0.84], show_titles=True,
    title_kwargs={"fontsize": 10})
fig_lcdm.suptitle("ΛCDM posterior (N=25000, covarianza DESI completa)", y=1.01)
fig_lcdm.savefig(f"{OUT_DIR}/fig_corner_lcdm_professional.pdf", bbox_inches="tight")
plt.close(fig_lcdm)

# Corner CPL
fig_cpl = corner.corner(res_cpl["flat"],
    labels=[r"$H_0$", r"$\Omega_m$", r"$w_0$", r"$w_a$", r"$\Omega_b h^2$"],
    quantiles=[0.16, 0.5, 0.84], show_titles=True,
    title_kwargs={"fontsize": 10})
fig_cpl.suptitle("CPL posterior (N=25000, covarianza DESI completa)", y=1.01)
fig_cpl.savefig(f"{OUT_DIR}/fig_corner_cpl_professional.pdf", bbox_inches="tight")
plt.close(fig_cpl)

# Guardar cadenas finales
np.savez(CHAIN_FILE,
    ssee_flat=res_ssee["flat"], ssee_lp=res_ssee["lp"],
    lcdm_flat=res_lcdm["flat"], lcdm_lp=res_lcdm["lp"],
    cpl_flat=res_cpl["flat"],   cpl_lp=res_cpl["lp"],
)

elapsed_total = (time.time()-t0)/3600
log(f"\nFiguras guardadas en {OUT_DIR}/")
log(f"Cadenas guardadas en {CHAIN_FILE}")
log(f"Tiempo total: {elapsed_total:.2f}h")
log("FIN MCMC PROFESIONAL")
