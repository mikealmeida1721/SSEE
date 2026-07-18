"""
SSEE Paper 2 — MCMC PRODUCCIÓN bajo prior H_alg (reframe ω_m-directo)
=====================================================================
Variante de producción del MCMC de Paper 2 para el reframe ω_m-DIRECTO
(2026-06-18, OP-8 cerrado; prior MIRA 67.037 RETIRADO):

  Prior H₀:  H_MIRA (67.037, retirado)  →  H_alg = 3(φ+π)² = 67.962

El ancla 67.962 NO es ad hoc: con ω_b=(π−φ)/(3Ω²) y ω_c=KAL₀·ω_b·n_s
FIJOS por álgebra, la verosimilitud CMB plik_lite se MINIMIZA en
H₀=67.962 (scan results/logs/p3_h0anchor_reframe.log; χ²=1005.50 mín).
Es decir, el H global de fondo y el ancla CMB coinciden: 67.962.

Tamaño producción: 100 walkers × 25000 steps (= original Paper 2).
Solo corre SSEE (no ΛCDM ni CPL — ya están bien establecidos).

Hipótesis a falsar:
  Posterior H₀ = compromiso entre el ancla CMB (67.962) y DESI tardío.
  El sector DESI usa Ω_m,dyn=0.160 (INALTERADO por el reframe).
"""
import numpy as np
import time, os, sys, warnings
_SSEE_DATA = os.environ.get("SSEE_DATA_DIR") or ("/mnt/datos/SSEE_data" if os.path.isdir("/mnt/datos") else "results/data")  # portable: HDD si existe, si no results/ local
warnings.filterwarnings("ignore")
import emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner
import os as _reloc_os, sys as _reloc_sys  # reloc: anclar src/
_reloc_sys.path.insert(0, _reloc_os.path.dirname(_reloc_os.path.dirname(_reloc_os.path.abspath(__file__))))
from ssee_core import (
    PHI, PI, BETA, KAL0, P_SC, K_V, T_R, M_V,
    W0, WA, OMEGA_DE, OMEGA_M_TOTAL, OMEGA_CDM_SECTOR,
    OMEGA_M_CMB_MIRA, H0_ALG,
)

t0 = time.time()
LOG = "results/logs/mcmc_paper2_reframe.log"
# Regla de disco: la cadena (~1 GB) va al HDD, NO al SSD root (91% lleno)
CKPT = _SSEE_DATA + "/mcmc/paper2_reframe/mcmc_paper2_reframe_ckpt.npz"
OUT = "results/figures"
os.makedirs("results/logs", exist_ok=True)
os.makedirs(OUT, exist_ok=True)

def log(msg):
    line = f"[{(time.time()-t0)/60:6.1f}m] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f: f.write(line + "\n")

log("=" * 70)
log("SSEE — MCMC PRODUCCIÓN bajo prior H_alg (67.962, reframe ω_m-directo)")
log("=" * 70)
log(f"  Ω_m,total (geometría) = {OMEGA_M_TOTAL:.8f}  |  Ω_cdm,sector = {OMEGA_CDM_SECTOR:.8f}")
log(f"  Ω_m,CMB (ω_m/h², reframe) = 0.30889  (sin factor; OP-8 cerrado)")
log(f"  w0 = {W0:.10f},  wa = {WA:.10f}")
log(f"  H0_alg = {H0_ALG:.6f}")

C_KM = 2.998e5
FNU_SSEE = 0.020

# ─── Datos DESI DR2 (2503.14738 Tabla 4) — FUENTE ÚNICA data/raw/desi_dr2_bao.csv ───
# NO hardcodear (evita drift DR1/DR2). src/ ya está en path (línea ~29).
from desi_dr2_data import load_desi_dr2, desi_covariance  # noqa: E402
_DESI        = load_desi_dr2()
DESI_Z       = list(_DESI["z"])
_QSHORT      = {"DM_over_rd": "DM_rd", "DH_over_rd": "DH_rd", "DV_over_rd": "DV_rd"}
DESI_TYPE    = [_QSHORT[q] for q in _DESI["quantity"]]
DESI_OBS     = _DESI["value"]
DESI_COV_INV = np.linalg.inv(desi_covariance(_DESI))   # bloque-diagonal (r_MH oficiales DR2)

CLUSTERS = [
    {"M_ig": 1.8, "dM_obs": 1.0, "M_obs": 9.8 },
    {"M_ig": 2.2, "dM_obs": 1.2, "M_obs": 12.0},
    {"M_ig": 1.5, "dM_obs": 1.0, "M_obs": 8.0 },
    {"M_ig": 1.2, "dM_obs": 1.0, "M_obs": 6.5 },
]

# ─── PRIOR H_alg EXACTO (reframe ω_m-directo) ───
# Ancla CMB del reframe: con ω_b y ω_c FIJOS por álgebra SSEE, la
# verosimilitud plik_lite TTTEEE se minimiza en H₀ = 3(φ+π)² = 67.962
# (scan: results/logs/p3_h0anchor_reframe.log; χ²=1005.50 mín, ΔBIC=−23.93).
# El H global de fondo y el ancla CMB coinciden — no son dos números.
# σ = 0.54 (error Planck H₀ propagado, conservador).
# (prior MIRA 67.037 RETIRADO: usaba Ω_m,CMB=MIRA×Ω_m,dyn, factor disuelto OP-8)
MIRA_H0 = (67.962, 0.54)
BBN_OBH2 = (0.02218, 0.00055)

def f_de_cpl(z, w0, wa):
    a = 1.0/(1.0+z)
    return (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def E_ssee(z):
    return np.sqrt(OMEGA_M_TOTAL*(1+z)**3 + (1.0-OMEGA_M_TOTAL)*f_de_cpl(z, W0, WA))  # geometría: materia total 0.30889

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
    om_h2  = OMEGA_M_TOTAL*(H0/100)**2   # r_d con materia total (geometría)
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
log("RESULTADO MCMC PRODUCCIÓN — Prior H_alg (67.962, reframe)")
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
log(f"  Retirado (prior MIRA 67.037, 100w×25k): H₀ = 66.53 ± 0.44")
log(f"  ESTE (prior H_alg 67.962, 100w×25k):    H₀ = {H0_med:.3f} ± {H0_std:.3f}")

# Tensión con DESI puro (resultado exploratorio H0=65.530)
desi_pure_H0 = 65.530
delta = abs(H0_med - desi_pure_H0)
log(f"\n  Distancia a DESI-puro (65.530): {delta:.3f} km/s/Mpc")

# Figura
fig = corner.corner(flat,
    labels=[r"$H_0$ [km/s/Mpc]", r"$\Omega_b h^2$"],
    quantiles=[0.16, 0.5, 0.84], show_titles=True,
    title_kwargs={"fontsize": 10})
fig.suptitle(f"SSEE posterior — Prior H_alg (67.962, reframe)\nH₀={H0_med:.3f}±{H0_std:.3f}, MAP={H0_map:.3f}", y=1.02)
fig.savefig(f"{OUT}/fig_corner_ssee_halg_prior.pdf", bbox_inches="tight")
plt.close(fig)
log(f"\nFigura: {OUT}/fig_corner_ssee_halg_prior.pdf")
log(f"Cadena: {CKPT}")
log(f"Tiempo total: {(time.time()-t0)/60:.1f} min")
