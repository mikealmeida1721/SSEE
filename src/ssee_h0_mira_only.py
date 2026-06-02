"""
SSEE — MCMC corto: solo prior MIRA (H₀ = 67.08 ± 0.54)
========================================================
Aprovecha resultados previos para Planck/SSEE-alg/plano (ssee_h0_prior_experiment.py).
Esta corrida solo añade el prior MIRA-from-Planck que faltaba.

Datos previos (de ssee_h0_prior_experiment.py, 50w × 10k, seed 42):
  Prior Planck (67.36)  →  H₀ = 66.751 ± 0.445
  Prior SSEE (67.96)    →  H₀ = 67.159 ± 0.445
  Prior plano U(50,90)  →  H₀ = 65.530 ± 0.761
"""
import numpy as np
import time, os, sys, warnings
warnings.filterwarnings("ignore")
import emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from ssee_core import (PHI, PI, KAL0, W0, WA, OMEGA_M_DYN, H0_ALG)

# ─── Setup idéntico al script previo ───
C_KM = 2.998e5
FNU_SSEE = 0.020

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

def f_de_cpl(z, w0, wa):
    a = 1.0 / (1.0 + z)
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
        dm = (C_KM/H0) * DC(z)
        dh = C_KM / (H0 * E_ssee(z))
        if qty == "DM_rd":   preds.append(dm / rd)
        elif qty == "DH_rd": preds.append(dh / rd)
        else:                preds.append((z*dm**2*dh)**(1/3) / rd)
    return np.array(preds)

def ll_bao(H0, om_h2, ob_h2):
    rd = sound_horizon_rd(ob_h2, om_h2)
    r  = predict_desi(H0, rd) - DESI_OBS
    return -0.5 * (r @ DESI_COV_INV @ r)

LLC_CLUSTERS = -0.5 * sum(((c["M_ig"]*KAL0*(1+FNU_SSEE) - c["M_obs"])/c["dM_obs"])**2
                          for c in CLUSTERS)

PRIOR_MIRA = (67.08, 0.54)

def lpost_mira(theta):
    H0, ob_h2 = theta
    if not (40 < H0 < 100): return -np.inf
    if not (0.015 < ob_h2 < 0.030): return -np.inf
    lp_H0  = -0.5*((H0-PRIOR_MIRA[0])/PRIOR_MIRA[1])**2
    lp_bbn = -0.5*((ob_h2-0.02218)/0.00055)**2
    om_h2  = OMEGA_M_DYN*(H0/100)**2
    return lp_H0 + lp_bbn + ll_bao(H0, om_h2, ob_h2) + LLC_CLUSTERS

N_WALKERS, N_STEPS, N_BURN = 50, 10000, 1000

print("=" * 78)
print(f"MCMC corto — Prior MIRA: N({PRIOR_MIRA[0]}, {PRIOR_MIRA[1]})")
print(f"  Comparar contra (ya corridos):")
print(f"    Planck (67.36) → 66.751 ± 0.445")
print(f"    SSEE   (67.96) → 67.159 ± 0.445")
print(f"    Plano  (50,90) → 65.530 ± 0.761")
print("=" * 78)

t0 = time.time()
rng = np.random.default_rng(42)
pos = np.array([65.0, 0.02237]) + rng.standard_normal((N_WALKERS, 2)) * np.array([3.0, 0.0005])
sampler = emcee.EnsembleSampler(N_WALKERS, 2, lpost_mira)
pos, _, _ = sampler.run_mcmc(pos, N_BURN, progress=False)
sampler.reset()
sampler.run_mcmc(pos, N_STEPS, progress=False)
elapsed = time.time() - t0

flat = sampler.get_chain(flat=True)
lp   = sampler.get_log_prob(flat=True)
idx_map = np.argmax(lp)

H0_med = np.median(flat[:,0]); H0_std = np.std(flat[:,0])
H0_p16 = np.percentile(flat[:,0], 16); H0_p84 = np.percentile(flat[:,0], 84)
H0_map = flat[idx_map, 0]
ob_map = flat[idx_map, 1]
ll_b_map = ll_bao(H0_map, OMEGA_M_DYN*(H0_map/100)**2, ob_map)

print(f"\nResultado Prior MIRA:")
print(f"  H₀ = {H0_med:.3f} +{H0_p84-H0_med:.3f}/-{H0_med-H0_p16:.3f}  km/s/Mpc")
print(f"  MAP = {H0_map:.3f}    ln L_BAO = {ll_b_map:.2f}    accept = {np.mean(sampler.acceptance_fraction):.3f}")
print(f"  tiempo = {elapsed:.1f}s")

# ─── Tabla comparativa completa ───
print("\n" + "=" * 78)
print("TABLA COMPLETA — 4 priors")
print("=" * 78)
print(f"  {'Prior':<22} {'H₀ mediana':>12} {'±':>8} {'distancia DESI(65.53)':>22}")
print("  " + "-"*70)
priors_resumen = [
    ("Planck (67.36)",  66.751, 0.445),
    ("MIRA-P3 (67.08)", H0_med, H0_std),
    ("SSEE-alg (67.96)", 67.159, 0.445),
    ("Plano U(50,90)",   65.530, 0.761),
]
for name, med, std in priors_resumen:
    dist = abs(med - 65.530)
    print(f"  {name:<22} {med:>12.3f} {std:>8.3f} {dist:>22.3f} km/s/Mpc")

# Guardar cadena MIRA
np.savez("results/logs/h0_mira_only.npz", mira=flat, mira_lp=lp)

# Figura 4 priors
print("\nGenerando figura comparativa 4 priors...")
prev = np.load("results/logs/h0_three_priors.npz")
fig, ax = plt.subplots(figsize=(11, 5))
chains = {
    "Prior Planck (67.36)": (prev["planck"][:,0], "#1f77b4"),
    "Prior MIRA-P3 (67.08)": (flat[:,0],          "#ff7f0e"),
    "Prior SSEE-alg (67.96)": (prev["ssee"][:,0],  "#d62728"),
    "Prior plano (50-90)":    (prev["flat"][:,0],  "#2ca02c"),
}
for lab, (data, col) in chains.items():
    ax.hist(data, bins=80, density=True, alpha=0.5, color=col, label=lab)
    ax.axvline(np.median(data), color=col, ls="--", lw=1.5)
for x, lab in [(67.36,"Planck"), (67.08,"MIRA"), (67.962,"SSEE")]:
    ax.axvline(x, color="gray", ls=":", alpha=0.5)
    ax.text(x, ax.get_ylim()[1]*0.95 if hasattr(ax,'get_ylim') else 1.0, lab,
            rotation=90, fontsize=8, alpha=0.6)
ax.set_xlabel(r"$H_0$ [km/s/Mpc]")
ax.set_ylabel("posterior density")
ax.set_title("Posterior H₀ bajo 4 priors — DESI DR2 + BBN + cúmulos (SSEE fondo)")
ax.legend(fontsize=9, loc="upper left")
ax.set_xlim(63, 70)
plt.tight_layout()
os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/fig_h0_four_priors.pdf", bbox_inches="tight")
print("Figura: results/figures/fig_h0_four_priors.pdf")
