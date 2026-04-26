# SSEE-V3.6 — External Audit Guide

**Author:** Mike Edison Almeida Vallejo  
**ORCID:** 0009-0008-2195-7836  
**Repo:** https://github.com/mikealmeida1721/SSEE  
**Date:** 2026-04-25

---

## What this framework claims

SSEE-V3.6 (Structural Self-Energy Expansion) is a **zero-free-parameter** dark energy model.
All cosmological predictions are derived algebraically from two constants: the golden ratio φ
and π. No fitting to data is performed to obtain the central predictions.

The four falsifiable predictions — all pre-data (timestamped Zenodo: 10.5281/zenodo.19679049, 2026-01-28):

| Observable | SSEE prediction | Observed | Status |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.827, −0.75) | 0.05σ |
| CMB peak ℓ₁ | 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 | Planck 2018: 0.3153 | 0.63σ |
| n_s | 1 − φ⁻⁷ = 0.96556 | Planck 2018: 0.9649 | 0.2σ |

**Future falsifier (LiteBIRD ~2032):** tensor-to-scalar ratio r = φ⁻¹⁰ = 0.0083.

---

## Repo structure

```
SSEE/
├── manuscript/                     — LaTeX source files (authoritative)
│   ├── SSEE_Paper1_Framework.tex   — Paper 1: zero-parameter framework + EFT
│   ├── SSEE_EFT_section.tex        — EFT appendix (input'd by Paper 1)
│   ├── SSEE_Paper2_draft.tex       — Paper 2: MCMC validation
│   ├── SSEE_Paper3_draft.tex       — Paper 3: CMB confrontation
│   ├── SSEE_appendix_Friedmann.tex — Appendix
│   ├── ssee_paper3.bib
│   ├── abstracts_arXiv.txt
│   └── cover_letter_paper*.txt
├── sandbox_unificado/              — Paper 4 LaTeX + SSEE web compendium (submodule)
│   └── SSEE_Paper4_ToE.tex
├── src/
│   ├── ssee_paper2_mcmc.py         — MCMC: SSEE vs ΛCDM vs CPL (emcee)
│   ├── ssee_paper2_analysis.py     — Analytical w₀-wₐ plane
│   ├── ssee_paper2_figures.py      — Figure generation
│   └── ssee_paper3_cmb.py          — CAMB CMB spectrum (TT+TE+EE+lensing)
├── data/raw/                       — Observational data (CSV + Planck PR4 spectra)
├── results/
│   ├── figures/                    — All paper figures (PDF+PNG)
│   ├── logs/                       — MCMC chains (.npz) and run logs
│   └── *.txt                       — Diagnostic outputs
├── docs/                           — Compiled PDFs for review (4 papers only)
│   ├── SSEE_Paper1_Framework.pdf
│   ├── SSEE_Paper2_draft.pdf
│   ├── SSEE_Paper3_draft.pdf
│   └── SSEE_Paper4_ToE.pdf
├── archive/                        — Superseded versions (see archive/README.md)
├── submission_packages/            — arXiv-ready .tar.gz for each paper
└── AUDIT.md                        — This file
```

---

## How to reproduce the results

### Prerequisites
```bash
pip install camb emcee scipy numpy matplotlib
```

### Paper 2 — MCMC validation
```bash
python3 src/ssee_paper2_mcmc.py
```
Produces posterior chains for SSEE, ΛCDM, CPL against DESI DR2 + Planck 2018 priors + clusters.
Runtime: ~30–60 min (N_eff ≈ 637,500 for SSEE, 100 walkers).

Expected output:
- H₀ = 66.75⁺⁰·⁴⁴₋₀.₄⁴ km/s/Mpc
- χ²_r clusters = 0.122 (4 clusters, IGIMF-corrected)
- χ²_2D (w₀-wₐ vs DESI) = 0.080 → 0.05σ

### Paper 3 — CMB power spectrum
```bash
python3 src/ssee_paper3_cmb.py
```
Uses local Planck PR4 data in `data/raw/`. Runs CAMB, computes χ² and BIC.

Expected output:
```
TT: SSEE χ²_r=1.062  |  ΛCDM χ²_r=1.043  (N=1971)
TE: SSEE χ²_r=1.053  |  ΛCDM χ²_r=1.040  (N=1967)
EE: SSEE χ²_r=1.040  |  ΛCDM χ²_r=1.039  (N=1967)
PP: SSEE χ²_r=0.730  |  ΛCDM χ²_r=0.757  (N=9)  [lensing — SSEE favoured]
ΔBIC(TT, k=1 vs k=6) = −40.3  [Cobaya/plik_lite]
ΔBIC(TT, k=0 vs k=6) = −6.9   [diagonal covariance]
ΔBIC(TT+TE+EE+PP, upper bound) = +13.7
Peak positions: ℓ = 221, 538, 815
```

---

## Known limitations — full disclosure

### 1. H(z) tension (Paper 2, §Results)
SSEE χ²_r(H(z)) = 1.861 vs ΛCDM = 0.458.
Genuine tension (4× worse), not a calibration artefact. Paper 2 states this explicitly.

### 2. ΔBIC = +218 (Paper 2, §BIC) — background penalty
Applies BIC within the standard Friedmann framework, penalising SSEE's modified background
against ΛCDM's 6 fitted parameters. Paper 2 has an explicit paragraph distinguishing this
from the dynamic-sector test (ΔBIC = −5.55, SSEE favoured). These are two different physical
questions and should not be conflated.

### 3. Two-sector Ωm (Papers 2 & 3)
Ωm,dyn = 0.160 ≠ Ωm,CMB = 0.3199. Bridged by algebraic MIRA = 1.9989.
Full derivation from first principles deferred (Level 3 work).

### 4. MIRA — not post-hoc
MIRA defined in Genesis 5.12 (Zenodo: 10.5281/zenodo.19679049, 2026-01-28),
83 days before the CMB analysis (2026-04-19).

### 5. Diagonal CMB likelihood (Paper 3)
χ²_r uses diagonal covariance. Off-diagonal terms required for PRD/PRL submission.
The Cobaya/plik_lite test (ΔBIC = −40.3) uses the proper likelihood.

### 6. Eckart viscosity in EFT section
Eckart formulation violates causality (Hiscock-Lindblom 1985). Israel-Stewart (IS)
formulation is implemented for Ωc h² derivation; full perturbation predictions pending.

### 7. Ωc h² — IS-corrected to −0.6σ
Static Eckart: 3.7σ. IS derivation: KAL₀ × Ωb h² × n_s = 0.11926 → −0.6σ from Planck.
Physical justification: n_s = 1 − φ⁻⁷ acts as IS inflationary modulator (Paper 4).

### 8. ΔBIC(TT+TE+EE+PP) = +13.7 — upper bound
Four spectra combined ignoring correlations. Full joint test requires Planck CMB likelihood.
TT-only ΔBIC ≈ −6.9 (diagonal) / −40.3 (plik_lite) is the primary result.

---

## What would falsify SSEE

1. CMB peak ℓ₁ outside 221 ± 3 by more than 2σ in a new measurement
2. DESI DR2+ requiring Ωm > 0.20 in the dynamic BAO sector
3. Ratio Ωm,CMB/Ωm,dyn significantly different from MIRA = 1.9989 in independent CMB experiments
4. χ²_r(CMB) > 2 in future datasets (CMB-S4, Simons Observatory)
5. Tensor-to-scalar ratio r ≠ φ⁻¹⁰ = 0.0083 measured by LiteBIRD (~2032)

---

## Predictive Register — anti-selection-bias defence

Chronological record of what was derived before what data:

| Quantity | Algebraic formula | Value | Test dataset | Status |
|---|---|---|---|---|
| w₀ | −Tr/Mv | −0.8399 | DESI DR2 (2025) | Pre-data |
| wₐ | −Psc/Kv | −0.6699 | DESI DR2 (2025) | Pre-data |
| MIRA | (φ+β)/2 | 1.9989 | CMB/BAO ratio | Pre-data |
| Ωm | (π−φ)/(π+φ) | 0.3201 | Planck 2018 | Retrodiction |
| n_s | 1−φ⁻⁷ | 0.96556 | Planck 2018 | Retrodiction |
| H₀ | 3(φ+π)² | 67.962 | Planck 2018 | Retrodiction |
| r | φ⁻¹⁰ | 0.0083 | LiteBIRD (~2032) | Future prediction |

Genesis 5.12 commit: https://github.com/mikealmeida1721/SSEE_UNIFICADO (2026-01-28)

---

## Reference to companion repo

Genesis 5.12 (prior definitions of all algebraic constants including MIRA):
https://github.com/mikealmeida1721/SSEE_UNIFICADO

Zenodo DOI: https://doi.org/10.5281/zenodo.19679049
