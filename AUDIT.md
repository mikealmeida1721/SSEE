# SSEE-V3.6 — External Audit Guide

**Author:** Mike Edison Almeida Vallejo  
**ORCID:** 0009-0008-2195-7836  
**Repo:** https://github.com/mikealmeida1721/SSEE  
**Date:** 2026-04-20

---

## What this framework claims

SSEE-V3.6 (Structural Self-Energy Expansion) is a **zero-free-parameter** dark energy model.
All cosmological predictions are derived algebraically from two constants: the golden ratio φ
and π. No fitting to data is performed to obtain the central predictions.

The three falsifiable predictions are:

| Observable | SSEE prediction | Observed | Status |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.827, −0.75) | 0.05σ separation |
| CMB peak ℓ₁ | 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 | Planck 2018: 0.3153 | 0.63σ |

---

## Repo structure

```
SSEE/
├── manuscript/          — LaTeX source files
│   ├── SSEE_Paper2_draft.tex   — Paper 2: MCMC validation (submission-ready)
│   ├── SSEE_Paper3_draft.tex   — Paper 3: CMB confrontation
│   ├── ssee_paper2.bib
│   └── ssee_paper3.bib
├── src/
│   ├── ssee_paper2_mcmc.py     — MCMC analysis (emcee, N_s=6000)
│   ├── ssee_paper2_analysis.py — Analytical w₀-wₐ plane
│   └── ssee_paper3_cmb.py      — CAMB CMB spectrum (TT+TE+EE)
├── data/raw/                   — Planck PR4 TT/TE/EE spectra (downloaded)
├── results/figures/            — All paper figures (PDF)
├── docs/                       — Compiled PDFs for review
│   ├── SSEE_Paper1_Framework_v3.6.pdf
│   ├── SSEE_Paper2_MCMC_Validation_draft.pdf  (outdated — use .tex)
│   └── SSEE_Paper3_CMB_Confrontation_v1.pdf
└── archive/                    — Obsolete drafts (see archive/README.md)
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
Runtime: ~10 min (N_s=6000, N_w=32).

Expected output:
- H₀ = 66.66⁺⁰·⁴⁷₋₀.₄₆ km/s/Mpc
- χ²_r clusters = 0.131 (7 clusters, IGIMF-corrected)
- χ²_2D (w₀-wₐ vs DESI) = 0.080 → 0.05σ

### Paper 3 — CMB power spectrum
```bash
python3 src/ssee_paper3_cmb.py
```
Downloads Planck PR4 TT/TE/EE data, runs CAMB, computes χ² and BIC.

Expected output:
```
TT: SSEE χ²_r=1.062  |  ΛCDM χ²_r=1.043  (N=1971)
TE: SSEE χ²_r=1.053  |  ΛCDM χ²_r=1.040  (N=1967)
EE: SSEE χ²_r=1.040  |  ΛCDM χ²_r=1.039  (N=1967)
PP: SSEE χ²_r=0.730  |  ΛCDM χ²_r=0.757  (N=9)  [lensing, SSEE favoured]
ΔBIC(TT, k=0 vs k=6) = -6.9
ΔBIC(TT+TE+EE+PP combined, upper bound) = +13.7
Peak positions: ℓ = 221, 538, 815
Lensing data: Planck 2018 MV bandpowers (arXiv:1807.06209, Table 1)
              sourced from CobayaSampler/planck_supp_data_and_covmats
```

---

## Known limitations — full disclosure

These are documented honestly in the papers. Reviewers should verify they are addressed:

### 1. H(z) tension (Paper 2, §Results)
SSEE χ²_r(H(z)) = 1.861 vs ΛCDM = 0.458.
This is a **genuine tension** (4× worse), not a calibration artefact.
Paper 2 states this explicitly. It reflects the fixed H₀=66.66 and modified background.

### 2. ΔBIC = +206 (Paper 2, §BIC)
This applies BIC within the standard Friedmann framework, where SSEE's modified
background is penalised against ΛCDM's 6 fitted parameters.
Paper 2 includes an explicit "Scope of this penalty" paragraph distinguishing this
from model-independent falsification. The CMB sector ΔBIC = −6.9 (Paper 3) provides
the complementary test under the SSEE background.

### 3. Two-sector Ωm (Papers 2 & 3)
Ωm,dyn = 0.160 (BAO/clusters) ≠ Ωm,CMB = 0.160 × MIRA = 0.3199.
The physical justification is discussed in Paper 3 §2.3: MIRA encodes the ratio of
observational depth to dynamical depth across the full expansion history.
A full derivation from first principles is deferred (Level 3 work).

### 4. MIRA — not post-hoc
MIRA was defined in SSEE Genesis 5.12 (commit f8728c3, 2026-01-28):
https://github.com/mikealmeida1721/SSEE_UNIFICADO
The CMB analysis was performed 2026-04-19 — **83 days after** the commit.
Paper 3 §3.2 documents this timeline explicitly with the commit hash.

### 5. No Lagrangian action (Paper 1)
SSEE is currently a phenomenological framework. A formal Lagrangian connecting
it to known EFT of dark energy is identified as Level 3 work. The papers frame
SSEE as a "zero-parameter phenomenological framework", not a theory from first principles.

### 6. ΔBIC(TT+TE+EE+PP) = +13.7 — upper bound
When TT, TE, EE, PP are combined ignoring correlations, ΔBIC ≈ +13.7 (ΛCDM favoured).
This is an upper bound: the four spectra are not independent and a proper
joint test requires the full Planck CMB likelihood. The TT-only result ΔBIC ≈ −7
is the clean primary result. Both are reported in Paper 3 §4.3.
Note: in the lensing sector alone, SSEE has χ²_r=0.730 vs ΛCDM 0.757 (N=9 bins).

---

## What would falsify SSEE

From Paper 3 §5:

1. CMB peak ℓ₁ outside 221 ± 3 by more than 2σ in a new measurement
2. DESI DR2+ requiring Ωm > 0.20 in the dynamic BAO sector
3. Ratio Ωm,CMB/Ωm,dyn significantly different from MIRA=1.9989 in independent CMB experiments
4. χ²_r(CMB) > 2 in future datasets (CMB-S4, Simons Observatory)

---

## Reference to companion repo

Genesis 5.12 (prior definitions of all algebraic constants including MIRA):
https://github.com/mikealmeida1721/SSEE_UNIFICADO
