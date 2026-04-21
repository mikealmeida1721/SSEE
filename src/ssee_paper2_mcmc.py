"""
SSEE-V3.6 — Paper 2: Bayesian MCMC Validation (v2)
Modelos: SSEE vs ΛCDM vs CPL
Likelihoods: DESI DR2 BAO + Planck 2018 prior comprimido + Masas de cúmulos
Nota diagnóstica: SSEE fija Ω_m,eff = 0.160, lo que modifica r_d y altera
la escala H0 preferida respecto a ΛCDM — tensión cuantificada explícitamente.
"""

import numpy as np
from scipy import integrate, stats
import warnings
warnings.filterwarnings("ignore")

try:
    import emcee
    import corner
except ImportError:
    import subprocess, sys
    subprocess.check_call([sys.executable, "-m", "pip", "install",
                           "emcee", "corner", "--break-system-packages", "-q"])
    import emcee, corner

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("=" * 65)
print("SSEE-V3.6 — PAPER 2: BAYESIAN MCMC VALIDATION (v2)")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. CONSTANTES SSEE (algebraicas — cero parámetros libres)
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

print(f"\n  Constantes SSEE fijas:")
print(f"    w0={W0_SSEE:.4f}, wa={WA_SSEE:.4f}")
print(f"    Ω_DE={OMDE_SSEE:.4f}, Ω_m,eff={OM_EFF_SSEE:.4f}, KAL0={KAL0:.4f}")

# ─────────────────────────────────────────────────────────────
# 2. FÍSICA DEL FONDO
# ─────────────────────────────────────────────────────────────

C_KM = 2.998e5

def f_de_cpl(z, w0, wa):
    a = 1.0 / (1.0 + z)
    return (1 + z)**(3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

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

DESI_BAO = [
    (0.295, "DV_rd",  7.93, 0.15),
    (0.510, "DM_rd", 13.62, 0.25),
    (0.510, "DH_rd", 20.08, 0.60),
    (0.706, "DM_rd", 16.85, 0.32),
    (0.706, "DH_rd", 19.50, 0.55),
    (0.930, "DM_rd", 21.71, 0.28),
    (0.930, "DH_rd", 17.88, 0.35),
    (1.317, "DM_rd", 27.79, 0.69),
    (1.317, "DH_rd", 13.82, 0.42),
    (1.491, "DM_rd", 30.21, 0.79),
    (1.491, "DH_rd", 13.23, 0.55),
    (2.330, "DM_rd", 39.71, 0.94),
    (2.330, "DH_rd",  8.52, 0.17),
]

# Prior comprimido Planck 2018 (arXiv:1807.06209)
PLANCK_H0   = (67.36, 0.54)
PLANCK_OM   = (0.3153, 0.0073)
PLANCK_OBH2 = (0.02237, 0.00015)
RHO_H0_OM   = -0.85

s0, s1, s2 = PLANCK_H0[1], PLANCK_OM[1], PLANCK_OBH2[1]
PLANCK_COV = np.array([
    [s0**2,        RHO_H0_OM*s0*s1, 0.0  ],
    [RHO_H0_OM*s0*s1, s1**2,        0.0  ],
    [0.0,          0.0,              s2**2],
])
PLANCK_COV_INV = np.linalg.inv(PLANCK_COV)
PLANCK_MU = np.array([PLANCK_H0[0], PLANCK_OM[0], PLANCK_OBH2[0]])

CLUSTERS = [
    {"M_ig": 1.8, "dM_obs": 1.0, "M_obs": 9.8 },
    {"M_ig": 2.2, "dM_obs": 1.2, "M_obs": 12.0},
    {"M_ig": 1.5, "dM_obs": 1.0, "M_obs": 8.0 },
    {"M_ig": 1.2, "dM_obs": 1.0, "M_obs": 6.5 },
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

def ll_bao(H0, om_h2, ob_h2, E_func, *Eargs):
    rd = sound_horizon_rd(ob_h2, om_h2)
    chi2 = 0.0
    for (z, qty, obs, sig) in DESI_BAO:
        dm = (C_KM/H0) * DC(z, lambda zz: E_func(zz, *Eargs))
        dh = C_KM / (H0 * E_func(z, *Eargs))
        if qty == "DM_rd":   pred = dm / rd
        elif qty == "DH_rd": pred = dh / rd
        else:                pred = (z*dm**2*dh)**(1/3) / rd
        chi2 += ((pred - obs)/sig)**2
    return -0.5 * chi2

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
    lb = ll_bao(H0, om_h2, ob_h2, E_ssee)
    lc = ll_clusters(KAL0, FNU_SSEE)
    return lp_bbn + lp_H0 + lb + lc

def lpost_lcdm(theta):
    H0, Om, ob_h2 = theta
    if not (40 < H0 < 100): return -np.inf
    if not (0.15 < Om < 0.55): return -np.inf
    if not (0.015 < ob_h2 < 0.030): return -np.inf
    om_h2 = Om*(H0/100)**2
    return ll_planck(H0, Om, ob_h2) + ll_bao(H0, om_h2, ob_h2, E_lcdm, Om)

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
    return lp + ll_bao(H0, om_h2, ob_h2, E_cpl, Om, w0, wa)

# ─────────────────────────────────────────────────────────────
# 6. MCMC
# ─────────────────────────────────────────────────────────────

N_WALKERS = 32
N_STEPS   = 25000
N_BURN    = 1500

def run_mcmc(log_post, theta0, scales, ndim, label):
    print(f"\n  [{label}]  ndim={ndim}, {N_WALKERS}w × {N_STEPS}s")
    rng = np.random.default_rng(42)
    pos = theta0 + rng.standard_normal((N_WALKERS, ndim)) * scales
    sampler = emcee.EnsembleSampler(N_WALKERS, ndim, log_post)
    sampler.run_mcmc(pos, N_STEPS, progress=True)
    flat = sampler.get_chain(discard=N_BURN, flat=True)
    lp   = sampler.get_log_prob(discard=N_BURN, flat=True)
    mask = np.isfinite(lp)
    flat, lp = flat[mask], lp[mask]
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        tau_max = np.max(tau)
        n_eff = N_WALKERS*(N_STEPS-N_BURN)/tau_max
    except Exception:
        tau_max, n_eff = np.nan, np.nan
    idx = np.argmax(lp)
    print(f"    τ_max={tau_max:.1f}  N_eff≈{n_eff:.0f}  ln P_MAP={lp[idx]:.2f}")
    return {
        "label": label, "ndim": ndim, "flat": flat, "lp": lp,
        "theta_map": flat[idx], "lp_map": lp[idx],
        "medians": np.median(flat,0), "stds": np.std(flat,0),
        "p16": np.percentile(flat,16,0), "p84": np.percentile(flat,84,0),
    }

print("\nIniciando MCMC...")
res_ssee = run_mcmc(lpost_ssee,
                    np.array([62.0, 0.02237]),
                    np.array([2.0,  0.0003]), 2, "SSEE-V3.6")
res_lcdm = run_mcmc(lpost_lcdm,
                    np.array([67.4, 0.315, 0.02237]),
                    np.array([1.5,  0.015, 0.0003]), 3, "ΛCDM")
res_cpl  = run_mcmc(lpost_cpl,
                    np.array([67.4, 0.315, -0.90, -0.40, 0.02237]),
                    np.array([1.5,  0.015,  0.10,  0.25, 0.0003]), 5, "CPL")

# ─────────────────────────────────────────────────────────────
# 7. COMPARACIÓN DE MODELOS
# ─────────────────────────────────────────────────────────────

N_DATA = len(DESI_BAO) + 3
models = [res_ssee, res_lcdm, res_cpl]
for r in models:
    r["BIC"] = r["ndim"]*np.log(N_DATA) - 2*r["lp_map"]
    r["AIC"] = 2*r["ndim"] - 2*r["lp_map"]
bic_min = min(r["BIC"] for r in models)
aic_min = min(r["AIC"] for r in models)

print("\n" + "="*65)
print("COMPARACIÓN DE MODELOS")
print("="*65)
print(f"\n  {'Modelo':<14} {'k':>3} {'ln P_MAP':>10} {'BIC':>8} {'ΔBIC':>7} {'AIC':>8} {'ΔAIC':>7}")
print("  "+"-"*57)
for r in models:
    print(f"  {r['label']:<14} {r['ndim']:>3} {r['lp_map']:>10.2f} "
          f"{r['BIC']:>8.2f} {r['BIC']-bic_min:>7.2f} "
          f"{r['AIC']:>8.2f} {r['AIC']-aic_min:>7.2f}")

# ─────────────────────────────────────────────────────────────
# 8. PARÁMETROS Y DIAGNÓSTICO DE TENSIÓN
# ─────────────────────────────────────────────────────────────

param_names = {
    "SSEE-V3.6": ["H₀", "Ω_b·h²"],
    "ΛCDM":      ["H₀", "Ω_m", "Ω_b·h²"],
    "CPL":       ["H₀", "Ω_m", "w₀", "wₐ", "Ω_b·h²"],
}
print("\n"+"="*65)
print("PARÁMETROS POSTERIORES")
print("="*65)
for r in models:
    print(f"\n  [{r['label']}]")
    for nm, med, p16, p84 in zip(param_names[r["label"]],
                                  r["medians"], r["p16"], r["p84"]):
        print(f"    {nm:<12} = {med:.5f}  +{p84-med:.5f}/-{med-p16:.5f}")
    if r["label"] == "SSEE-V3.6":
        print(f"    {'Ω_m,eff':<12} = {OM_EFF_SSEE:.4f}  [algebraico]")
        print(f"    {'w₀,wₐ':<12} = {W0_SSEE:.4f}, {WA_SSEE:.4f}  [algebraico]")

# r_d predicho
def get_rd(r):
    H0 = r["medians"][0]; ob = r["medians"][-1]
    Om = OM_EFF_SSEE if r["label"] == "SSEE-V3.6" else r["medians"][1]
    return sound_horizon_rd(ob, Om*(H0/100)**2)

rd_ssee = get_rd(res_ssee)
rd_lcdm = get_rd(res_lcdm)
H0_ssee = res_ssee["medians"][0]; H0_ssee_s = res_ssee["stds"][0]
H0_lcdm = res_lcdm["medians"][0]
t_H0 = abs(H0_ssee-PLANCK_H0[0])/np.sqrt(H0_ssee_s**2+PLANCK_H0[1]**2)
t_Om = abs(OM_EFF_SSEE-PLANCK_OM[0])/PLANCK_OM[1]

print(f"""
  DIAGNÓSTICO DE TENSIÓN SSEE vs PLANCK 2018:
    r_d(SSEE) = {rd_ssee:.2f} Mpc  vs  r_d(ΛCDM) = {rd_lcdm:.2f} Mpc  (×{rd_ssee/rd_lcdm:.3f})
    H₀(SSEE)  = {H0_ssee:.2f} ± {H0_ssee_s:.2f}  vs  Planck = {PLANCK_H0[0]} ± {PLANCK_H0[1]}
    Tensión H₀: {t_H0:.2f}σ  |  Tensión Ω_m: {t_Om:.2f}σ

  La tensión es consecuencia directa de Ω_m,eff=0.160 (SSEE) frente a
  Ω_m=0.315 (Planck), que alarga r_d en ×{rd_ssee/rd_lcdm:.3f} y requiere H₀ más
  bajo para preservar los ratios D/r_d. Esto es la manifestación BAO de
  la tensión Ω_DE (Sección 3). La reconciliación requiere el análisis
  CMB completo con la ecuación de Friedmann modificada (Paper 3).
""")

# ─────────────────────────────────────────────────────────────
# 9. POSTERIOR PREDICTIVE CHECK
# ─────────────────────────────────────────────────────────────

def H_pred(r, z):
    th = r["theta_map"]
    if r["label"] == "SSEE-V3.6": return th[0]*E_ssee(z)
    if r["label"] == "ΛCDM":      return th[0]*E_lcdm(z, th[1])
    return th[0]*E_cpl(z, th[1], th[2], th[3])

chi2_cc = {r["label"]: np.sum(((H_pred(r,Z_CC)-H_CC)/DH_CC)**2)/len(Z_CC)
           for r in models}
print("  PPC — Cosmic Chronometers (χ²/N):")
for nm, v in chi2_cc.items():
    print(f"    {nm:<14}: {v:.3f}")

# ─────────────────────────────────────────────────────────────
# 10. FIGURAS
# ─────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":"serif","font.size":11,"axes.labelsize":12,
    "figure.dpi":150,"savefig.dpi":300,"savefig.bbox":"tight",
    "xtick.direction":"in","ytick.direction":"in",
})
colors = {"SSEE-V3.6":"#E6002B","ΛCDM":"#2166AC","CPL":"#4DAF4A"}
ls_map  = {"SSEE-V3.6":"-","ΛCDM":"--","CPL":"-."}

# Fig 5: Corner SSEE
fig5 = corner.corner(res_ssee["flat"], labels=[r"$H_0$",r"$\Omega_b h^2$"],
    quantiles=[0.16,0.50,0.84], show_titles=True,
    title_kwargs={"fontsize":11}, color="#E6002B", plot_datapoints=False,
    truths=[PLANCK_H0[0], PLANCK_OBH2[0]], truth_color="#2166AC")
fig5.suptitle("SSEE-V3.6 posteriors (azul = Planck 2018)", y=1.02, fontsize=11)
fig5.savefig("fig5_corner_ssee.pdf"); fig5.savefig("fig5_corner_ssee.png")
plt.close(fig5); print("\nFig 5 guardada")

# Fig 6: Corner ΛCDM
fig6 = corner.corner(res_lcdm["flat"],
    labels=[r"$H_0$",r"$\Omega_m$",r"$\Omega_b h^2$"],
    quantiles=[0.16,0.50,0.84], show_titles=True,
    title_kwargs={"fontsize":11}, color="#2166AC", plot_datapoints=False)
fig6.suptitle("ΛCDM posteriors (DESI DR2 + Planck 2018)", y=1.02, fontsize=11)
fig6.savefig("fig6_corner_lcdm.pdf"); fig6.savefig("fig6_corner_lcdm.png")
plt.close(fig6); print("Fig 6 guardada")

# Fig 7: H(z)
z_plot = np.linspace(0, 1.5, 200)
fig7, ax7 = plt.subplots(figsize=(8,5.5))
ax7.errorbar(Z_CC, H_CC, yerr=DH_CC, fmt="o", color="black",
             ms=5, capsize=3, zorder=6, label="Cosmic Chronometers")
for r in models:
    H_v = np.array([H_pred(r, z) for z in z_plot])
    ax7.plot(z_plot, H_v, color=colors[r["label"]], lw=2.3,
             ls=ls_map[r["label"]], label=f"{r['label']} (MAP, χ²/N={chi2_cc[r['label']]:.2f})")
idx_s = np.random.default_rng(0).integers(0, len(res_ssee["flat"]), 300)
H_band = np.array([res_ssee["flat"][i][0]*E_ssee(z_plot) for i in idx_s])
ax7.fill_between(z_plot, np.percentile(H_band,16,0), np.percentile(H_band,84,0),
                 color="#E6002B", alpha=0.15, label="SSEE 68% posterior")
ax7.set_xlabel("Redshift $z$"); ax7.set_ylabel(r"$H(z)$ [km s$^{-1}$ Mpc$^{-1}$]")
ax7.set_title("Posterior Predictive Check — $H(z)$")
ax7.legend(fontsize=9); ax7.set_xlim(0,1.5)
fig7.tight_layout()
fig7.savefig("fig7_Hz_comparison.pdf"); fig7.savefig("fig7_Hz_comparison.png")
plt.close(fig7); print("Fig 7 guardada")

# Fig 8: Tensiones — H0 posteriors + r_d
fig8, (ax8a, ax8b) = plt.subplots(1, 2, figsize=(11,4.5))
for r in models:
    ax8a.hist(r["flat"][:,0], bins=50, density=True,
              color=colors[r["label"]], alpha=0.55, label=r["label"])
ax8a.axvline(PLANCK_H0[0], color="black", ls="--", lw=1.5,
             label=f"Planck: {PLANCK_H0[0]}±{PLANCK_H0[1]}")
ax8a.fill_betweenx([0,10], PLANCK_H0[0]-PLANCK_H0[1], PLANCK_H0[0]+PLANCK_H0[1],
                   color="black", alpha=0.12)
ax8a.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]")
ax8a.set_ylabel("Densidad posterior")
ax8a.set_title(r"Posterior $H_0$ — comparación de modelos")
ax8a.legend(fontsize=9); ax8a.set_ylim(0,None)

rd_vals = [get_rd(r) for r in models]
y_pos   = np.arange(len(models))
bars = ax8b.barh(y_pos, rd_vals,
                 color=[colors[r["label"]] for r in models],
                 alpha=0.75, edgecolor="black", linewidth=0.6)
ax8b.axvline(147.27, color="black", ls="--", lw=1.5,
             label=r"$r_d^{\rm fid}=147.27$ Mpc")
ax8b.set_yticks(y_pos)
ax8b.set_yticklabels([r["label"] for r in models])
ax8b.set_xlabel("$r_d$ [Mpc]")
ax8b.set_title("Horizonte de sonido $r_d$ (MAP)")
ax8b.legend(fontsize=9)
for bar, v in zip(bars, rd_vals):
    ax8b.text(v+0.5, bar.get_y()+bar.get_height()/2,
              f"{v:.1f}", va="center", fontsize=10)
fig8.tight_layout()
fig8.savefig("fig8_tension_summary.pdf"); fig8.savefig("fig8_tension_summary.png")
plt.close(fig8); print("Fig 8 guardada")

# ─────────────────────────────────────────────────────────────
# 11. RESUMEN FINAL
# ─────────────────────────────────────────────────────────────

print("\n"+"="*65)
print("RESUMEN EJECUTIVO PAPER 2")
print("="*65)
print(f"""
  SSEE (2 params libres):
    H₀ = {H0_ssee:.2f} ± {H0_ssee_s:.2f} km/s/Mpc
    r_d = {rd_ssee:.1f} Mpc  (ΛCDM: {rd_lcdm:.1f} Mpc, factor ×{rd_ssee/rd_lcdm:.3f})
    Tensión H₀ vs Planck: {t_H0:.2f}σ  |  Tensión Ω_m: {t_Om:.2f}σ

  ΛCDM (3 params):
    H₀ = {H0_lcdm:.2f} ± {res_lcdm['stds'][0]:.2f} km/s/Mpc
    Ω_m = {res_lcdm['medians'][1]:.4f} ± {res_lcdm['stds'][1]:.4f}

  ΔBIC (respecto al mejor):
    SSEE={res_ssee['BIC']-bic_min:+.2f}  ΛCDM={res_lcdm['BIC']-bic_min:+.2f}  CPL={res_cpl['BIC']-bic_min:+.2f}

  PPC H(z) (χ²/N):
    SSEE={chi2_cc['SSEE-V3.6']:.3f}  ΛCDM={chi2_cc['ΛCDM']:.3f}  CPL={chi2_cc['CPL']:.3f}
""")
print("Figuras 5–8 guardadas. MCMC completado.")
print("="*65)
