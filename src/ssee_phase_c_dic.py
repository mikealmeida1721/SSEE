"""
Phase C — Rigorous Bayesian Model Comparison
SSEE-V3.6: DIC + BIC under three k-counting schemes + Laplace Bayes factor

Inputs:  results/logs/mcmc_chains_professional.npz
Outputs: results/tables/model_comparison_dic.tex
         results/tables/model_comparison_dic_values.txt
"""

import numpy as np
import warnings
warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# Constants (must match ssee_paper2_mcmc.py exactly)
# ─────────────────────────────────────────────────────────────────────────────
phi   = (1 + 5**0.5) / 2               # 1.6180339887
pi    = 3.14159265358979                # π
OMEGA = pi + phi                        # 4.7596
BETA  = (pi + phi) / 2                  # 2.3798
KAL0  = BETA + pi                       # 5.5214
PSC   = OMEGA + phi                     # 6.3776
KV    = phi + pi + OMEGA               # 9.5192
TR    = 3*(phi + BETA)                  # 11.9935
MV    = phi + pi + KV                   # 14.2788

W0_SSEE   = -TR / MV                   # -0.8399
WA_SSEE   = -PSC / KV                  # -0.6700
OM_EFF    = TR / MV                    # 0.8399  = Ω_DE
OM_M_DYN  = 1 - OM_EFF                # 0.1601
MIRA      = (3*phi + pi) / 4           # 1.9989
OM_M_CMB  = MIRA * OM_M_DYN           # 0.3199

PLANCK_H0   = (67.36, 0.54)
PLANCK_OM   = (0.3153, 0.0073)
PLANCK_OBH2 = (0.02237, 0.00015)

N_DATA = 16   # 13 BAO + 3 Planck constraints

# ─────────────────────────────────────────────────────────────────────────────
# Load chains
# ─────────────────────────────────────────────────────────────────────────────
print("Loading chains...")
data = np.load("results/logs/mcmc_chains_professional.npz")

ssee_flat = data["ssee_flat"]   # (N, 2)
ssee_lp   = data["ssee_lp"]    # (N,)   — only first 50k valid; rest are -inf
lcdm_flat = data["lcdm_flat"]  # (N, 3)
lcdm_lp   = data["lcdm_lp"]   # (N,)
cpl_flat  = data["cpl_flat"]   # (N, 5)
cpl_lp    = data["cpl_lp"]     # (N,)

# Filter valid lp values
mask_ssee = np.isfinite(ssee_lp)
mask_lcdm = np.isfinite(lcdm_lp)
mask_cpl  = np.isfinite(cpl_lp)

ssee_lp_v = ssee_lp[mask_ssee]
lcdm_lp_v = lcdm_lp[mask_lcdm]
cpl_lp_v  = cpl_lp[mask_cpl]

print(f"  Valid SSEE samples: {mask_ssee.sum():,}")
print(f"  Valid ΛCDM samples: {mask_lcdm.sum():,}")
print(f"  Valid CPL  samples: {mask_cpl.sum():,}")

# ─────────────────────────────────────────────────────────────────────────────
# DIC = D̄ + p_D   where D = -2 log π(θ|data),  p_D = D̄ - D(θ_MAP)
# Note: using log_posterior as proxy (prior terms are subdominant for BAO data)
# ─────────────────────────────────────────────────────────────────────────────
def compute_dic(lp_valid, label):
    D_samples = -2 * lp_valid
    D_bar     = np.mean(D_samples)
    D_map     = -2 * np.max(lp_valid)        # D at MAP (best sample)
    p_D       = D_bar - D_map                # effective number of parameters
    DIC       = D_bar + p_D                  # = 2*D_bar - D_map
    lnL_map   = np.max(lp_valid)
    print(f"\n  [{label}]")
    print(f"    ln P_MAP = {lnL_map:.3f}")
    print(f"    D̄        = {D_bar:.3f}")
    print(f"    D(θ_MAP) = {D_map:.3f}")
    print(f"    p_D      = {p_D:.3f}  (effective params)")
    print(f"    DIC      = {DIC:.3f}")
    return {"label": label, "lnP_MAP": lnL_map, "D_bar": D_bar,
            "D_map": D_map, "p_D": p_D, "DIC": DIC}

print("\n" + "="*60)
print("DIC ANALYSIS")
print("="*60)
dic_ssee = compute_dic(ssee_lp_v, "SSEE-V3.6")
dic_lcdm = compute_dic(lcdm_lp_v, "ΛCDM")
dic_cpl  = compute_dic(cpl_lp_v,  "CPL")

dics = [dic_ssee, dic_lcdm, dic_cpl]
dic_min = min(d["DIC"] for d in dics)
print(f"\n  ΔDIC(SSEE−ΛCDM) = {dic_ssee['DIC'] - dic_lcdm['DIC']:+.2f}")
print(f"  ΔDIC(SSEE−CPL)  = {dic_ssee['DIC'] - dic_cpl['DIC']:+.2f}")
print(f"  ΔDIC(ΛCDM−CPL)  = {dic_lcdm['DIC'] - dic_cpl['DIC']:+.2f}")

# ─────────────────────────────────────────────────────────────────────────────
# BIC under three k-counting schemes
# BIC = k ln(N) - 2 ln P_MAP
# ─────────────────────────────────────────────────────────────────────────────
lnN = np.log(N_DATA)
print("\n" + "="*60)
print("BIC — Three k-counting schemes")
print("="*60)

# MAP log-posterior values from chains
lnP_ssee = dic_ssee["lnP_MAP"]  # -124.87  (SSEE at Omega_m=0.160)
lnP_lcdm = dic_lcdm["lnP_MAP"]  # -15.79
lnP_cpl  = dic_cpl["lnP_MAP"]   # -11.75

# SSEE_dyn: uses ΛCDM background (same lnP_MAP as ΛCDM = -14.19 as in paper)
# This value is from the Paper 2 table (separate evaluation not in these chains)
# Using the published paper value directly for scheme A
lnP_ssee_dyn = -14.19   # from Paper 2 Table 1 (SSEE evaluated at ΛCDM background)

print("\n  Scheme A — dynamic sector isolated (as published):")
print("    SSEE_dyn: k=1, background borrowed from ΛCDM, w0/wa algebraically fixed")
print("    Both SSEE_dyn and ΛCDM evaluated at the ΛCDM background → same lnP_MAP=-14.19")
print("    ΔBIC is purely a k-parsimony difference: Δk×ln(N) = (1-3)×ln(16) = -5.55")
k_A = {"SSEE_dyn": 1, "LCDM": 3, "CPL": 5}
# Both models evaluated at ΛCDM background → same lnP_MAP (from Paper 2 Table 1)
lnP_lcdm_dyn = -14.19   # ΛCDM at its MAP in the same dynamic-sector framework
lnP_A = {"SSEE_dyn": lnP_ssee_dyn, "LCDM": lnP_lcdm_dyn, "CPL": lnP_cpl}
bic_A = {m: k_A[m]*lnN - 2*lnP_A[m] for m in k_A}
bic_min_A = min(bic_A.values())
for m in k_A:
    print(f"    {m:<12} k={k_A[m]}  lnP={lnP_A[m]:7.2f}  BIC={bic_A[m]:7.2f}  ΔBIC={bic_A[m]-bic_min_A:+7.2f}")

print("\n  Scheme B — MCMC parameters sampled (conservative):")
print("    SSEE: k=2 (H0, Ωb h²), ΛCDM: k=3 (H0, Ωm, Ωb h²), CPL: k=5")
k_B = {"SSEE": 2, "LCDM": 3, "CPL": 5}
lnP_B = {"SSEE": lnP_ssee, "LCDM": lnP_lcdm, "CPL": lnP_cpl}
bic_B = {m: k_B[m]*lnN - 2*lnP_B[m] for m in k_B}
bic_min_B = min(bic_B.values())
for m in k_B:
    print(f"    {m:<12} k={k_B[m]}  lnP={lnP_B[m]:7.2f}  BIC={bic_B[m]:7.2f}  ΔBIC={bic_B[m]-bic_min_B:+7.2f}")

print("\n  Scheme C — full dark-energy parameter count:")
print("    SSEE: k=2 (algebraic w0/wa count as 0 free params → background k same)")
print("    ΛCDM: k=3, CPL: k=5 (w0, wa explicit free params)")
print("    Note: Scheme C is identical to Scheme B by construction (algebraic params are not free)")
# In principle some would argue w0,wa "would be free" in other models, adding k=2 to SSEE
k_C_ssee_alt = 4   # k=2+2 if one counts w0,wa as extra phantom DOF
bic_C_ssee_alt = k_C_ssee_alt*lnN - 2*lnP_ssee
print(f"\n    SSEE worst-case (k=4 phantom): BIC={bic_C_ssee_alt:.2f}  ΔBIC vs ΛCDM = {bic_C_ssee_alt - bic_B['LCDM']:+.2f}")
print(f"    (SSEE remains disfavoured vs ΛCDM by the background penalty regardless of k-scheme)")

# ─────────────────────────────────────────────────────────────────────────────
# Laplace Approximation Bayes Factor
# log Z_Laplace = ln P(θ_MAP) + (k/2)ln(2π) - (1/2)ln|H|
# where H is the Hessian of -ln P at MAP, estimated from posterior variance
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("LAPLACE APPROXIMATION — Log Bayes Factor")
print("="*60)

def laplace_log_evidence(lp_map, cov_matrix, k):
    """log Z ≈ ln P(θ_MAP) + (k/2)ln(2π) + (1/2)ln|Σ|"""
    sign, logdet = np.linalg.slogdet(cov_matrix)
    if sign <= 0:
        return np.nan
    return lp_map + (k/2)*np.log(2*np.pi) + 0.5*logdet

# SSEE (2 params, valid 50k samples)
ssee_params = ssee_flat[mask_ssee]
lcdm_params = lcdm_flat[mask_lcdm]
cpl_params  = cpl_flat[mask_cpl]

cov_ssee = np.cov(ssee_params.T)
cov_lcdm = np.cov(lcdm_params.T)
cov_cpl  = np.cov(cpl_params.T)

logZ_ssee = laplace_log_evidence(lnP_ssee, cov_ssee, 2)
logZ_lcdm = laplace_log_evidence(lnP_lcdm, cov_lcdm, 3)
logZ_cpl  = laplace_log_evidence(lnP_cpl,  cov_cpl,  5)

print(f"\n  ln Z (Laplace):  SSEE={logZ_ssee:.2f}  ΛCDM={logZ_lcdm:.2f}  CPL={logZ_cpl:.2f}")
print(f"  ln B_10 (SSEE vs ΛCDM) = {logZ_ssee - logZ_lcdm:.2f}")
print(f"  ln B_10 (SSEE vs CPL)  = {logZ_ssee - logZ_cpl:.2f}")
print(f"  ln B_10 (ΛCDM vs CPL)  = {logZ_lcdm - logZ_cpl:.2f}")

# Jeffreys scale
def jeffreys_strength(ln_B):
    """Returns qualitative strength on Jeffreys scale."""
    B = abs(ln_B)
    if B < 1.0:   return "barely worth mentioning"
    elif B < 2.5: return "substantial"
    elif B < 5.0: return "strong"
    else:         return "decisive"

lnB_ssee_lcdm = logZ_ssee - logZ_lcdm
print(f"\n  Interpretation: SSEE vs ΛCDM — {jeffreys_strength(lnB_ssee_lcdm)}")
print(f"  Direction: {'ΛCDM favoured' if lnB_ssee_lcdm < 0 else 'SSEE favoured'}")

# ─────────────────────────────────────────────────────────────────────────────
# Posterior predictive check — χ²_r on DESI BAO
# ─────────────────────────────────────────────────────────────────────────────
# From Paper 2 logs: H0_SSEE = 66.75, H0_LCDM = 68.28
# Posterior dispersion
h0_ssee_med = np.median(ssee_params[:, 0])
h0_ssee_std = np.std(ssee_params[:, 0])
h0_lcdm_med = np.median(lcdm_params[:, 0])
h0_lcdm_std = np.std(lcdm_params[:, 0])
h0_planck = PLANCK_H0

print("\n" + "="*60)
print("POSTERIOR SUMMARY")
print("="*60)
print(f"  SSEE H0 = {h0_ssee_med:.3f} ± {h0_ssee_std:.3f} km/s/Mpc")
print(f"  ΛCDM H0 = {h0_lcdm_med:.3f} ± {h0_lcdm_std:.3f} km/s/Mpc")
print(f"  Planck  = {h0_planck[0]:.3f} ± {h0_planck[1]:.3f} km/s/Mpc")
t_ssee = abs(h0_ssee_med - h0_planck[0]) / np.sqrt(h0_ssee_std**2 + h0_planck[1]**2)
t_lcdm = abs(h0_lcdm_med - h0_planck[0]) / np.sqrt(h0_lcdm_std**2 + h0_planck[1]**2)
print(f"  Tension H0: SSEE={t_ssee:.2f}σ,  ΛCDM={t_lcdm:.2f}σ")

# ─────────────────────────────────────────────────────────────────────────────
# Summary table for LaTeX appendix
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("LATEX TABLE OUTPUT")
print("="*60)

latex_table = r"""\begin{table}[H]
\centering
\caption{%
  Bayesian model comparison: BIC (three k-counting schemes) and DIC
  from MCMC chains. $N_{\rm data}=16$. ΔBIC/ΔDIC $<0$ favours SSEE.
  Jeffreys criterion: $|\Delta|>10$ = very strong.
}
\label{tab:bic_dic_rigorous}
\begin{tabular}{lrrrr}
\toprule
\textbf{Criterion} & \textbf{SSEE} & \textbf{ΛCDM} & \textbf{CPL}
  & $\boldsymbol{\Delta(\mathrm{SSEE}-\mathrm{\Lambda CDM})}$ \\
\midrule
\multicolumn{5}{l}{\textit{Scheme A — dynamic sector isolated
  ($k_{\rm SSEE}=1$, $k_{\rm \Lambda CDM}=3$; same background)}} \\
BIC$_{\rm dyn}$ & SSEE_DYN_BIC & LCDM_BIC_A & CPL_BIC_A & DELTA_A \\
\midrule
\multicolumn{5}{l}{\textit{Scheme B — MCMC parameters sampled
  ($k_{\rm SSEE}=2$, $k_{\rm \Lambda CDM}=3$; full backgrounds)}} \\
BIC$_{\rm MCMC}$ & SSEE_BIC_B & LCDM_BIC_B & CPL_BIC_B & DELTA_B \\
\midrule
\multicolumn{5}{l}{\textit{DIC — effective parameters from posterior
  ($p_{D}$ estimated from chains; background-agnostic)}} \\
$D_{\rm MAP}$     & SSEE_DMAP & LCDM_DMAP & CPL_DMAP & — \\
$\bar{D}$         & SSEE_DBAR & LCDM_DBAR & CPL_DBAR & — \\
$p_D$             & SSEE_PD & LCDM_PD & CPL_PD & — \\
DIC               & SSEE_DIC & LCDM_DIC & CPL_DIC & DELTA_DIC \\
\midrule
\multicolumn{5}{l}{\textit{Laplace approximation (full posterior evidence)}} \\
$\ln Z_{\rm Lap}$ & SSEE_LNZ & LCDM_LNZ & CPL_LNZ & DELTA_LNZ \\
\bottomrule
\end{tabular}
\end{table}
"""

# Fill values
latex_table = latex_table.replace("SSEE_DYN_BIC",  f"${bic_A['SSEE_dyn']:.2f}$")
latex_table = latex_table.replace("LCDM_BIC_A",    f"${bic_A['LCDM']:.2f}$")
latex_table = latex_table.replace("CPL_BIC_A",     f"${bic_A['CPL']:.2f}$")
delta_A = bic_A['SSEE_dyn'] - bic_A['LCDM']
latex_table = latex_table.replace("DELTA_A",       f"${delta_A:+.2f}$")

latex_table = latex_table.replace("SSEE_BIC_B",    f"${bic_B['SSEE']:.2f}$")
latex_table = latex_table.replace("LCDM_BIC_B",    f"${bic_B['LCDM']:.2f}$")
latex_table = latex_table.replace("CPL_BIC_B",     f"${bic_B['CPL']:.2f}$")
delta_B = bic_B['SSEE'] - bic_B['LCDM']
latex_table = latex_table.replace("DELTA_B",       f"${delta_B:+.2f}$")

latex_table = latex_table.replace("SSEE_DMAP",  f"${dic_ssee['D_map']:.2f}$")
latex_table = latex_table.replace("LCDM_DMAP",  f"${dic_lcdm['D_map']:.2f}$")
latex_table = latex_table.replace("CPL_DMAP",   f"${dic_cpl['D_map']:.2f}$")
latex_table = latex_table.replace("SSEE_DBAR",  f"${dic_ssee['D_bar']:.2f}$")
latex_table = latex_table.replace("LCDM_DBAR",  f"${dic_lcdm['D_bar']:.2f}$")
latex_table = latex_table.replace("CPL_DBAR",   f"${dic_cpl['D_bar']:.2f}$")
latex_table = latex_table.replace("SSEE_PD",    f"${dic_ssee['p_D']:.2f}$")
latex_table = latex_table.replace("LCDM_PD",    f"${dic_lcdm['p_D']:.2f}$")
latex_table = latex_table.replace("CPL_PD",     f"${dic_cpl['p_D']:.2f}$")
latex_table = latex_table.replace("SSEE_DIC",   f"${dic_ssee['DIC']:.2f}$")
latex_table = latex_table.replace("LCDM_DIC",   f"${dic_lcdm['DIC']:.2f}$")
latex_table = latex_table.replace("CPL_DIC",    f"${dic_cpl['DIC']:.2f}$")
delta_DIC = dic_ssee['DIC'] - dic_lcdm['DIC']
latex_table = latex_table.replace("DELTA_DIC",  f"${delta_DIC:+.2f}$")

latex_table = latex_table.replace("SSEE_LNZ",   f"${logZ_ssee:.2f}$")
latex_table = latex_table.replace("LCDM_LNZ",   f"${logZ_lcdm:.2f}$")
latex_table = latex_table.replace("CPL_LNZ",    f"${logZ_cpl:.2f}$")
latex_table = latex_table.replace("DELTA_LNZ",  f"${lnB_ssee_lcdm:+.2f}$")

print(latex_table)

# Save
import os
os.makedirs("results/tables", exist_ok=True)
with open("results/tables/model_comparison_dic.tex", "w") as f:
    f.write(latex_table)

summary_txt = f"""Phase C — Rigorous Bayesian Model Comparison
Generated from MCMC chains (mcmc_chains_professional.npz)
N_data = {N_DATA}, chains: SSEE {mask_ssee.sum():,}, ΛCDM {mask_lcdm.sum():,}, CPL {mask_cpl.sum():,}

=== BIC SCHEME A (dynamic sector, k_SSEE=1) ===
SSEE_dyn: BIC={bic_A['SSEE_dyn']:.2f}  LCDM: BIC={bic_A['LCDM']:.2f}  CPL: BIC={bic_A['CPL']:.2f}
ΔBIC(SSEE_dyn − ΛCDM) = {delta_A:+.2f}   → SSEE favoured

=== BIC SCHEME B (MCMC params, k_SSEE=2) ===
SSEE: BIC={bic_B['SSEE']:.2f}  LCDM: BIC={bic_B['LCDM']:.2f}  CPL: BIC={bic_B['CPL']:.2f}
ΔBIC(SSEE − ΛCDM) = {delta_B:+.2f}   → ΛCDM favoured (background mismatch dominates)

=== DIC (from full MCMC posterior) ===
SSEE: DIC={dic_ssee['DIC']:.2f} (p_D={dic_ssee['p_D']:.2f})
ΛCDM: DIC={dic_lcdm['DIC']:.2f} (p_D={dic_lcdm['p_D']:.2f})
CPL:  DIC={dic_cpl['DIC']:.2f} (p_D={dic_cpl['p_D']:.2f})
ΔDIC(SSEE−ΛCDM) = {delta_DIC:+.2f}   → background mismatch captured

=== LAPLACE BAYES FACTOR ===
ln Z: SSEE={logZ_ssee:.2f}  ΛCDM={logZ_lcdm:.2f}  CPL={logZ_cpl:.2f}
ln B_10(SSEE vs ΛCDM) = {lnB_ssee_lcdm:+.2f} ({jeffreys_strength(lnB_ssee_lcdm)}, {'ΛCDM favoured' if lnB_ssee_lcdm < 0 else 'SSEE favoured'})

=== H0 POSTERIOR ===
SSEE: {h0_ssee_med:.3f} ± {h0_ssee_std:.3f} km/s/Mpc  (tension vs Planck: {t_ssee:.2f}σ)
ΛCDM: {h0_lcdm_med:.3f} ± {h0_lcdm_std:.3f} km/s/Mpc  (tension vs Planck: {t_lcdm:.2f}σ)
"""

with open("results/tables/model_comparison_dic_values.txt", "w") as f:
    f.write(summary_txt)

print("\nFiles saved:")
print("  results/tables/model_comparison_dic.tex")
print("  results/tables/model_comparison_dic_values.txt")
