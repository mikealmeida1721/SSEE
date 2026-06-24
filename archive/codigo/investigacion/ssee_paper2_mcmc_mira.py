"""
SSEE Paper 2 — MCMC PRODUCCIÓN bajo prior MIRA
================================================
Variante del MCMC oficial de Paper 2 con un único cambio crítico:

  Prior H₀:  Planck (67.36 ± 0.54)  →  MIRA (67.04 ± 0.54)

Donde 67.04 es el H₀ que sale de Paper 3 al aplicar MIRA a Planck CMB
(Ω_m,CMB = MIRA × Ω_m,dyn). Es decir, "el H₀ que SSEE deriva de Planck
vía MIRA" en vez del "H₀ algebraico puro" o "H₀ Planck estándar".

Tamaño producción: 100 walkers × 25000 steps (= original Paper 2).
Solo corre SSEE (no ΛCDM ni CPL — ya están bien establecidos).

Hipótesis a falsar:
  Si posterior H₀ ≈ 67.04 ± ~0.4 → MIRA es un prior robusto, refuerza a
  MIRA como pivote físico del modelo.
  Si posterior H₀ se aleja → MIRA-from-Planck no encaja con DESI tardío.
"""
import numpy as np
import time, os, sys, warnings
warnings.filterwarnings("ignore")
import emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner
import os as _reloc_os, sys as _reloc_sys  # reloc: anclar src/
_reloc_sys.path.insert(0, _reloc_os.path.dirname(_reloc_os.path.dirname(_reloc_os.path.abspath(__file__))))
from ssee_core import (
    PHI, PI, BETA, KAL0, P_SC, K_V, T_R, M_V,
    W0, WA, OMEGA_DE, OMEGA_M_DYN,
    OMEGA_M_CMB_MIRA, H0_ALG,
)

t0 = time.time()
LOG = "results/logs/mcmc_paper2_mira.log"
CKPT = "results/logs/mcmc_paper2_mira_ckpt.npz"
OUT = "results/figures"
os.makedirs("results/logs", exist_ok=True)
os.makedirs(OUT, exist_ok=True)

def log(msg):
    line = f"[{(time.time()-t0)/60:6.1f}m] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f: f.write(line + "\n")

log("=" * 70)
log("SSEE-V3.6 — MCMC PRODUCCIÓN bajo prior MIRA (67.04)")
log("=" * 70)
log(f"  Ω_m,dyn (exacto) = {OMEGA_M_DYN:.10f}")
log(f"  Ω_m,CMB (MIRA)   = {OMEGA_M_CMB_MIRA:.10f}")
log(f"  w0 = {W0:.10f},  wa = {WA:.10f}")
log(f"  H0_alg = {H0_ALG:.6f}")

C_KM = 2.998e5
FNU_SSEE = 0.020

# ─── Datos DESI DR2 ───
DESI_Z    = [0.295, 0.510, 0.510, 0.706, 0.706, 0.930, 0.930,
             1.317, 1.317, 1.491, 1.491, 2.330, 2.330]
DESI_TYPE = ["DV_rd","DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd",
             "DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd"]
DESI_OBS  = np.array([7.93, 13.62, 20.08, 16.85, 19.50, 21.71, 17.88,
                      27.79, 13.82, 30.21, 13.23, 39.71,  8.52])
DESI_SIGMA = np.array([0.15, 0.25, 0.60, 0.32, 0.55, 0.28, 0.35,
                       0.69, 0.42, 0.79, 0.55, 0.94, 0.17])
RHO_PAIRS = {(1,2): -0.44, (3,4): -0.45, (5,6): -0.44,
             (7,8): -0.43, (9,10): -0.42, (11,12): -0.45}
DESI_COV = np.diag(DESI_SIGMA**2)
for (i, j), rho in RHO_PAIRS.items():
    DESI_COV[i,j] = rho * DESI_SIGMA[i] * DESI_SIGMA[j]
    DESI_COV[j,i] = DESI_COV[i,j]
DESI_COV_INV = np.linalg.inv(DESI_COV)

CLUSTERS = [
    {"M_ig": 1.8, "dM_obs": 1.0, "M_obs": 9.8 },
    {"M_ig": 2.2, "dM_obs": 1.2, "M_obs": 12.0},
    {"M_ig": 1.5, "dM_obs": 1.0, "M_obs": 8.0 },
    {"M_ig": 1.2, "dM_obs": 1.0, "M_obs": 6.5 },
]

# ─── PRIOR MIRA EXACTO ───
# H₀ derivado de re-correr P3 (Cobaya plik_lite + CAMB) con Ω_m,CMB exacto
# = MIRA × Ω_m,dyn = 0.319928188 (no el 0.320 redondeado).
# Resultado óptimo del scan: H_MIRA = 67.037 km/s/Mpc (SSEE@Σm_ν=0.069
# self-consistente; era 67.068 con mnu=0.06 baseline — 2026-06-09).
# σ = 0.54 (error Planck H₀ propagado; el del minimize_scalar es subestimado
# porque mantiene ns/As/τ/w0/wa fijos en vez de marginalizar).
MIRA_H0 = (67.037, 0.54)
BBN_OBH2 = (0.02218, 0.00055)

def f_de_cpl(z, w0, wa):
    a = 1.0/(1.0+z)
    return (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def E_ssee(z):
    return np.sqrt(OMEGA_M_DYN*(1+z)**3 + (1.0-OMEGA_M_DYN)*f_de_cpl(z, W0, WA))

def DC(z_max, n=300):
    zz = np.linspace(0, z_max, n)
    return np.trapezoid(1.0/E_ssee(zz), zz)

def sound_horizon_rd(ob_h2, om_h2):
    return 147.27 * (om_h2/0.1432)**(-0.255) * (ob_h2/0.02237)**(-0.134)

def predict_desi(H0, rd):
    preds = []
    for z, qty in zip(DESI_Z, DESI_TYPE):
        dm = (C_KM/H0)*DC(z)
        dh = C_KM/(H0*E_ssee(z))
        if qty == "DM_rd":   preds.append(dm/rd)
        elif qty == "DH_rd": preds.append(dh/rd)
        else:                preds.append((z*dm**2*dh)**(1/3)/rd)
    return np.array(preds)

def ll_bao(H0, om_h2, ob_h2):
    rd = sound_horizon_rd(ob_h2, om_h2)
    r  = predict_desi(H0, rd) - DESI_OBS
    return -0.5 * (r @ DESI_COV_INV @ r)

LL_CLUSTERS = -0.5 * sum(((c["M_ig"]*KAL0*(1+FNU_SSEE) - c["M_obs"])/c["dM_obs"])**2
                         for c in CLUSTERS)

def lpost(theta):
    H0, ob_h2 = theta
    if not (40 < H0 < 100): return -np.inf
    if not (0.015 < ob_h2 < 0.030): return -np.inf
    lp_H0  = -0.5*((H0-MIRA_H0[0])/MIRA_H0[1])**2
    lp_bbn = -0.5*((ob_h2-BBN_OBH2[0])/BBN_OBH2[1])**2
    om_h2  = OMEGA_M_DYN*(H0/100)**2
    return lp_H0 + lp_bbn + ll_bao(H0, om_h2, ob_h2) + LL_CLUSTERS

# ─── MCMC ───
N_W, N_S, N_B, SAVE = 100, 25000, 5000, 500
rng = np.random.default_rng(42)
pos = np.array([62.0, 0.02237]) + rng.standard_normal((N_W, 2)) * np.array([2.0, 0.0003])
sampler = emcee.EnsembleSampler(N_W, 2, lpost)

log(f"\nBurn-in {N_B} steps...")
pos, _, _ = sampler.run_mcmc(pos, N_B, progress=False)
sampler.reset()
log("Burn-in OK. Producción...")

all_chains = []
for i in range(0, N_S, SAVE):
    pos, _, _ = sampler.run_mcmc(pos, SAVE, progress=False)
    all_chains.append(sampler.get_chain(flat=True))
    log(f"  {i+SAVE}/{N_S} steps  accept={np.mean(sampler.acceptance_fraction):.3f}")
    np.savez(CKPT, chain=np.concatenate(all_chains, axis=0))

flat = sampler.get_chain(flat=True)
lp   = sampler.get_log_prob(flat=True)
mask = np.isfinite(lp)
idx  = np.argmax(lp) if np.any(mask) else 0

try:
    tau = sampler.get_autocorr_time(quiet=True)
    n_eff = N_W * N_S / np.max(tau)
except Exception:
    tau, n_eff = np.array([np.nan]), np.nan

H0_med = np.median(flat[:,0]); H0_std = np.std(flat[:,0])
H0_p16, H0_p84 = np.percentile(flat[:,0], [16, 84])
H0_map = flat[idx, 0]
ob_med = np.median(flat[:,1])

log("\n" + "=" * 70)
log("RESULTADO MCMC PRODUCCIÓN — Prior MIRA")
log("=" * 70)
log(f"  H₀     = {H0_med:.4f}  +{H0_p84-H0_med:.4f} / -{H0_med-H0_p16:.4f}  km/s/Mpc")
log(f"  H₀ MAP = {H0_map:.4f}")
log(f"  Ω_b·h² = {ob_med:.5f}")
log(f"  ln P_MAP = {lp[idx]:.3f}")
log(f"  τ_max = {np.max(tau):.1f}    N_eff ≈ {n_eff:.0f}")
log(f"  acceptance = {np.mean(sampler.acceptance_fraction):.3f}")

# BIC vs documentado en CLAUDE.md (SSEE k=2)
N_data = 13 + 2  # 13 BAO + 2 priors (H0 + BBN)
BIC = 2 * np.log(N_data) - 2 * lp[idx]
log(f"\n  BIC (k=2, N={N_data}): {BIC:.3f}")

# Comparación
log(f"\nComparación con resultados previos:")
log(f"  CLAUDE.md doc (prior Planck, 100w×25k): H₀ = 67.756 ± 0.442")
log(f"  Exploratorio (prior Planck, 50w×10k):   H₀ = 66.751 ± 0.445")
log(f"  Exploratorio (prior MIRA,   50w×10k):   H₀ = 66.559 ± 0.440")
log(f"  ESTE (prior MIRA, 100w×25k):            H₀ = {H0_med:.3f} ± {H0_std:.3f}")

# Tensión con DESI puro (resultado exploratorio H0=65.530)
desi_pure_H0 = 65.530
delta = abs(H0_med - desi_pure_H0)
log(f"\n  Distancia a DESI-puro (65.530): {delta:.3f} km/s/Mpc")

# Figura
fig = corner.corner(flat,
    labels=[r"$H_0$ [km/s/Mpc]", r"$\Omega_b h^2$"],
    quantiles=[0.16, 0.5, 0.84], show_titles=True,
    title_kwargs={"fontsize": 10})
fig.suptitle(f"SSEE-V3.6 posterior — Prior MIRA (67.04)\nH₀={H0_med:.3f}±{H0_std:.3f}, MAP={H0_map:.3f}", y=1.02)
fig.savefig(f"{OUT}/fig_corner_ssee_mira_prior.pdf", bbox_inches="tight")
plt.close(fig)
log(f"\nFigura: {OUT}/fig_corner_ssee_mira_prior.pdf")
log(f"Cadena: {CKPT}")
log(f"Tiempo total: {(time.time()-t0)/60:.1f} min")
