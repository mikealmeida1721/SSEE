# SSEE-V3.6 — External Audit Guide

**Author:** Mike Edison Almeida Vallejo  
**ORCID:** 0009-0008-2195-7836  
**Repo:** https://github.com/mikealmeida1721/SSEE  
**Date:** 2026-05-17 (10-paper suite — bibliography audit + Phase 1/2 hardening; Paper 1 Postulates D & S; OPEN_PROBLEMS OP-2/3/4/6 resolved, OP-1/5 partial)

---

## What this framework claims

SSEE-V3.6 (Structural Self-Energy Expansion) is a **minimal-parameter** dark energy model
(~3 effective parameters vs. 6 for $\Lambda$CDM). The shape of the background sector is derived
algebraically from two constants: the golden ratio φ and π. No fitting to data is performed to
obtain the central predictions of that core; the absolute vacuum scale enters as a single measured
input (Postulate D), and the dark-matter-origin extensions are tracked as open problems.

Falsifiable predictions — all pre-data (timestamped Zenodo: 10.5281/zenodo.19679049, 2026-01-28):

| Observable | SSEE prediction | Observed | Status |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.838±0.057, −0.631±0.226) | 0.05σ / 0.17σ |
| CMB peak ℓ₁ | 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 (algebraic via MIRA) | Planck 2018: 0.3153 | 0.63σ |
| n_s | 1 − φ⁻⁷ = 0.96556 | Planck 2018: 0.9649 | 0.2σ |
| m_φ (φ-DM mass) | 5.60 eV algebraic | gravitationally-produced scalar — probed via k_fs | Future prediction |
| k_fs | 0.493 h/Mpc algebraic | DESI Y3/Euclid P(k): 2026–2028 | Future prediction |
| r (tensor-to-scalar) | φ⁻¹⁰ = 0.00813 | LiteBIRD (~2032) | Future prediction |

---

## Repo structure

```
SSEE/
├── manuscript/                     — LaTeX sources (authoritative)
│   ├── SSEE_Paper1_Framework.tex … SSEE_Paper10_UVCompletion.tex  — the 10 papers
│   ├── SSEE_EFT_section.tex         — EFT section (\input by Paper 1)
│   ├── SSEE_appendix_Friedmann.tex  — Friedmann-derivation appendix
│   ├── SSEE_Unified_Journal.tex     — consolidated journal paper (Papers 1–7)
│   ├── SSEE_Endorser_Summary.tex    — 2-page arXiv endorser brief
│   └── *.bib                        — ssee_paper3/5/6, SSEE_Paper4, ssee_unified
├── src/                             — 42 Python scripts (analysis & verification)
│   ├── ssee_paper2…paper10_*.py     — per-paper analysis, MCMC, figures
│   ├── ssee_op1…op6_*.py            — OPEN_PROBLEMS resolution scripts (OP-1..OP-6)
│   ├── ssee_paperB_*.py             — Paper B groundwork (baryogenesis DW, N_*)
│   ├── ssee_audit_consistency.py    — cross-paper parameter consistency audit
│   └── scratch/                     — exploratory scripts (not load-bearing)
├── class_ssee/                      — CLASS Boltzmann fork — SSEE .ini configs + plots
│   ├── ssee_v36.ini                 — SSEE MIRA sector (Ω_m=0.3199)
│   ├── ssee_v36_nomira.ini          — SSEE dynamic sector only (Ω_m=0.160)
│   ├── ssee_v36_twosector.ini       — φ-DM two-sector (ncdm m=5.60 eV)
│   └── ssee_v36_IS.ini              — IS viscosity (cs2_fld=0.001)
├── data/                            — observational data (DESI DR2, Planck PR4, clusters)
├── results/                         — generated figures, tables, logs
├── docs/                            — compiled PDFs: 10 papers + endorser + unified journal
├── notebooks/ · submission_packages/ · archive/
├── README.md · CHANGELOG.md · OPEN_PROBLEMS.md · CITATION.cff
└── AUDIT.md                         — this file
```
*(`eftcamb_ssee/` — EFTCAMB fork — is not versioned; clone separately.)*

---

## How to reproduce the results

### Prerequisites
```bash
pip install camb emcee scipy numpy matplotlib corner cobaya classy getdist astropy
```

### Paper 2 — MCMC validation
```bash
python3 src/ssee_paper2_mcmc.py
```
Runtime: ~30–60 min (N_eff ≈ 637,500 for SSEE, 100 walkers × 25,000 steps).

Expected output:
```
H₀ = 66.75 ± 0.44 km/s/Mpc
χ²_r clusters = 0.122  (4 clusters, IGIMF-corrected, MCMC)
χ²_r clusters = 0.126  (7 clusters, analytic sample)
χ²_2D (w₀-wₐ vs DESI) = 0.080 → 0.05σ
ΔBIC (dynamic sector, k_SSEE=1 vs k_ΛCDM=3) = −5.55  [SSEE favoured]
ΔBIC (full background, k=0 vs k_ΛCDM=6) = +206        [framework penalty]
```

### Paper 3 — CMB power spectrum
```bash
python3 src/ssee_paper3_cmb.py
```

Expected output:
```
TT: SSEE χ²_r=1.062  |  ΛCDM χ²_r=1.043  (N=1971)
TE: SSEE χ²_r=1.053  |  ΛCDM χ²_r=1.040  (N=1967)
EE: SSEE χ²_r=1.040  |  ΛCDM χ²_r=1.039  (N=1967)
PP: SSEE χ²_r=0.730  |  ΛCDM χ²_r=0.757  (N=9)    [lensing — SSEE favoured]
ΔBIC (diagonal, TT+TE+EE+PP, k=2 vs k=6) = +31.1   [N amplifies small Δχ²_r]
ΔBIC (plik_lite Cobaya, TTTEEE, k=2 vs k=6)  = −31.3 [SSEE decisively favoured]
ΔBIC (conservative k=4 vs k=6) = −13.8              [SSEE strongly favoured]
Peak positions: ℓ = 221, 538, 815
```

CLASS cross-check (αK Bellini-Sawicki):
```bash
python3 src/ssee_paper3_hiclass_check.py
```
Expected: αK(0) = 0.4033, Δ = 0.005% vs algebraic prediction.

### Paper 4 — Press-Schechter δc comparison
```bash
python3 src/ssee_press_schechter.py
```
Expected output:
```
δc(SSEE) = 1.6284  (= δc,EdS × n_s)
δc(ΛCDM) = 1.6865
M [M☉]      σ_M(z=10)    n_SSEE/n_ΛCDM
1.00e+11      1.0103       1.061
3.00e+11      0.7267       1.159
1.00e+12      0.5064       1.406
```

### Paper 5 — Israel-Stewart causal perturbations
```bash
python3 src/ssee_paper5_IS_perturbations.py
```
Expected output:
```
c²_s,eff = 0 (exact algebraic)
MIRA_num (k≥10) = 0.989 ± 0.017   [background effect, not perturbative]
γ_IS = 0.554 ± 0.001               [≈ γ_ΛCDM = 0.55]
G = D₁_SSEE/D₁_ΛCDM = 0.866      [13.4% growth suppression]
σ₈_SSEE = 0.702 ± 0.005
S₈_SSEE = 0.725 ± 0.005            [1.96σ KiDS / 2.84σ DES]
fσ₈(z=0.5) = 0.341                 [mean 2.67σ, 6 surveys — baseline for Paper 6]
```

### Paper 6 — φ-DM two-sector
```bash
python3 src/ssee_paper6_verification.py
python3 src/ssee_paper6_mcmc.py    # 64 walkers × 15000 steps
```
Expected output (verification):
```
m_φ = 5.60 eV (algebraic, zero free parameters)
k_fs = 0.493 h/Mpc
Ω_CDM = 0.160050
Ω_φDM = 0.159878
Ω_total = 0.319928 ≈ Ω_m,CMB
σ₈_eff = 0.737  (S₈ = 0.761, 2.29σ DES)
Mean fσ₈ tension: 0.50σ (was 2.56σ, single-sector baseline corrected)
```
Expected output (MCMC):
```
Ω_φDM = 0.161 ± 0.011
ΔBIC = −12.1  [SSEE favoured]
χ²_r = 0.497
```

### Paper 7 — Canonical EFT βc plateau test
```bash
python3 src/ssee_eft_verification.py
```
Expected output:
```
βc = −3.990 ± 0.001  (8 initial conditions)
AURA = (3φ+π)/2 = 3.997847
Gap = 0.199%  [attributed to Ωb+Ωr excluded from background]
αT = αM = αB = 0 (exact)
αK(0) = 0.4033 (algebraic)
```

### CLASS Boltzmann — MIRA test and two-sector P(k)

```bash
# Build CLASS (requires gfortran + python3-dev):
cd class_ssee && make
# Run MIRA vs no-MIRA comparison:
./class ssee_v36.ini && ./class ssee_v36_nomira.ini
python3 plot_ssee_cmb.py
```

Expected output (peak positions):
```
SSEE+MIRA:   ℓ₁=220  ℓ₂=535  ℓ₃=810   RMS vs ΛCDM = 1.4%
SSEE no-MIRA: ℓ₁=240  ℓ₂=597  ℓ₃=922   RMS vs ΛCDM = 31.5%
```
*Without MIRA, all three CMB peaks shift ~10% and RMS rises 22×.*

Two-sector P(k) calibration:
```bash
./class ssee_v36_twosector.ini
python3 calibrate_wdm_alpha.py
```
Expected: α_WDM = 1.6561 h/Mpc, σ₈_eff = 0.737 (top-hat integral).

IS viscosity test:
```bash
./class ssee_v36_IS.ini && python3 plot_IS_viscosity.py
```
Expected: cs² effect on σ₈ = 0.03% (negligible); G = D₁_SSEE/D₁_ΛCDM = 0.866.

### OP-5 — HMcode-2020 baryonic feedback (requires CLASS + classy)

`src/ssee_op5_hmcode.py` needs CLASS compiled with the Python wrapper and
HMcode-2020 enabled. A bare `pip install classy` does **not** enable
HMcode-2020 — the local CLASS build is required:
```bash
cd class_ssee && make
cd python && python setup.py install
# Verify HMcode-2020 is available:
python -c "from classy import Class; c=Class(); c.set({'non_linear':'hmcode','hmcode_version':'2020_baryonic_feedback'})"
```

### MCMC Fase 4 — multi-probe background

5 parameters (H₀, w₀, wₐ, Ωm, r_d); 100 walkers × 11,000 steps; N_eff = 4,402.

Expected output:
```
H₀ = 67.34 ± 0.77  (algebraic: 67.96, tension: 0.81σ)
w₀ = −0.812 ± 0.046  (algebraic: −0.840, tension: 0.61σ)
wₐ = −0.659 ± 0.113  (algebraic: −0.670, tension: 0.10σ)
Ωm = 0.163 ± 0.012   (algebraic: 0.160, tension: 0.25σ)
r_d = 147.44 ± 0.89  (algebraic: 147.156, tension: 0.32σ)
Mean tension (5 params): 0.36σ
```

---

## Known limitations — full disclosure

### 1. H(z) tension (Paper 2, §Results)
SSEE χ²_r(H(z)) = 1.861 vs ΛCDM = 0.458.
Genuine tension (4× worse), not a calibration artefact. Paper 2 states this explicitly.
Physical explanation: SSEE's modified expansion history shifts H(z) predictions systematically.

### 2. ΔBIC = +206 (Paper 2, §BIC) — background penalty
Applies BIC within the standard Friedmann framework, penalising SSEE's modified background
against ΛCDM's 6 fitted parameters. Paper 2 explicitly distinguishes this from the
dynamic-sector test (ΔBIC = −5.55, SSEE favoured). These address different physical questions.

### 3. CMB pipeline — canonical result and H0-anchoring sensitivity (Paper 3)

**Canonical CMB result (Paper 3 §5):** H0 scanned freely over [66.5, 67.5]; minimum at
H0=67.066 km/s/Mpc. Δχ²=+3.793, ΔBIC=−31.3 (k=2 vs k=6, N=6469). Decisively favours SSEE.

**Illustrative H0-anchored comparison** (`results/tables/planck_fulllike_results.txt`):
When H0 is fixed to the DESI BAO value (66.66, calibrated at z~0.3–2.3), Δχ²=+144.8,
ΔBIC=+92.1. This is NOT a contradiction — it decomposes as: +141 from the 0.74σ
H0 offset between DESI and CMB preferred values, +3.5 from the intrinsic w0/wa effect.
The file is explicitly labelled as illustrative; the primary result is ΔBIC=−31.3.

**k count fixed** (BC2): SSEE uses k=2 in all CMB comparisons (H0 + Ωb h² prior borrowed
from Planck). k=0 counting appears only in the illustrative file with matching caveat.

**Parameter table** (authoritative):

| Parameter | Status in SSEE | k contribution |
|---|---|---|
| w₀, wₐ | algebraic (Zenodo pre-data) | 0 |
| n_s | algebraic (1−φ⁻⁷) | 0 |
| Ωm,CMB | algebraic (MIRA hypothesis) | 0 |
| H₀ | scanned to CMB minimum | 1 |
| Ωb h² | fixed to Planck prior | 1 (conservative) |
| τ, As | fixed to Planck priors | 0 (conservative option: +2) |
| **Total SSEE** | | **k=2** (k=4 hyper-conservative) |
| **ΛCDM baseline** | | **k=6** |

Diagonal covariance gives ΔBIC = +31.1 (N amplifies small Δχ²_r ≈ 0.01).
Full plik/CamSpec likelihood with official Planck covariance matrix remains a blocker for PRD/PRL.

### 4. Two-sector Ωm (Papers 2, 3, 5, 6)
Ωm,dyn = 0.160 ≠ Ωm,CMB = 0.3199. Bridged by algebraic MIRA = (3φ+π)/4 = 1.9989.
Full derivation of MIRA from first-principles field equations is deferred (Level 3 work).

### 5. MIRA — status and provenance (BC3)

**Formal status:** MIRA = (3φ+π)/4 = 1.998924 is a **phenomenological auxiliary hypothesis**
(not yet a derived quantity). It bridges Ωm,dyn=0.160 and Ωm,CMB=0.3199 algebraically,
but its derivation from first-principles field equations (Blocker B3) remains open.
It is NOT simultaneously a free parameter and a derived quantity — it is a fixed hypothesis
with documented pre-data provenance:

MIRA defined in Genesis 5.12 (Zenodo: 10.5281/zenodo.19679049, 2026-01-28),
83 days before the CMB analysis (2026-04-19). Anti-post-hoc guarantee.

This status (hypothesis pending formal derivation) is stated explicitly in Papers 1, 3, 5, 6.
All MIRA-dependent claims are sector-specific and conditional on this hypothesis.

### 6. Israel-Stewart perturbations (Papers 5, 6)
Full causal IS (Hiscock-Lindblom 1985) perturbation theory is implemented for background sector.
B-mode LiteBIRD forecasts requiring IS tensor perturbations remain future work.

### 7. Ωc h² (Paper 4)
Static Eckart: 3.7σ tension. IS derivation: KAL₀ × Ωb h² × n_s = 0.11926 → −0.6σ.
Note: Ωb h² algebraic = (π−φ)/(3Ω²) = 0.02242 (0.32σ from Planck — see OPEN_PROBLEMS OP-1; the earlier 3(π−φ)/200 form is superseded).

### 8. φ-DM mass scale and field content (Paper 6)
m_φ = 5.60 eV is not keV-scale WDM; φ-DM is modelled as a real scalar field (explicit
Lagrangian in Paper 6 §4) populated by gravitational particle production. The
Dodelson-Widrow (sterile-neutrino) route is excluded — its mixing angle has no
zero-parameter SSEE form. φ-DM has no Standard-Model portal, so it is not accessible
to neutrino-mass experiments; free-streaming suppression enters observation via
k_fs = 0.493 h/Mpc, not direct m_φ. Lyman-α bounds apply to thermal relics; the
scalar fraction f_φ ≈ 0.50 gives an effective ~0.1–0.5 keV equivalent — within
observational bounds. Quantified in Paper 6 §Lyman-α. The ab-initio relic abundance
from the gravitational-production integral is deferred to Paper B.

---

## What would falsify SSEE

1. CMB peak ℓ₁ outside 221 ± 3 by more than 2σ in a new measurement
2. DESI DR2+ requiring Ωm > 0.20 in the dynamic BAO sector
3. S₈ measured above 0.85 by Euclid weak-lensing (two-sector φ-DM predicts S₈ ≈ 0.76)
4. k_fs cutoff absent or at significantly different scale in Euclid/DESI Y3 P(k) (2026–28)
5. Tensor-to-scalar ratio r ≠ φ⁻¹⁰ = 0.00813 measured by LiteBIRD (~2032)
6. |w₀ + 0.840| > 3σ confirmed by DESI DR5 or Euclid

---

## Predictive Register — anti-selection-bias defence

Chronological record of what was derived before what data:

| Quantity | Algebraic formula | Value | Test dataset | Status |
|---|---|---|---|---|
| w₀ | −Tr/Mv | −0.8399 | DESI DR2 (2025) | Pre-data |
| wₐ | −Psc/Kv | −0.6699 | DESI DR2 (2025) | Pre-data |
| MIRA | (3φ+π)/4 | 1.9989 | CMB/BAO ratio | Pre-data (Zenodo 2026-01-28) |
| Ωm | (π−φ)/(π+φ) | 0.3201 | Planck 2018 | Retrodiction |
| n_s | 1 − φ⁻⁷ | 0.96556 | Planck 2018 | Retrodiction |
| H₀ | 3(φ+π)² | 67.962 | Planck 2018 | Retrodiction |
| r_d (MIRA) | CAMB, Ωm,CMB=0.3199 | 147.156 Mpc | Planck 2018: 147.09 ± 0.26 Mpc | 0.25σ (retrodiction) |
| m_φ | Σm_ν × H₀^alg | 5.60 eV | Euclid/DESI Y3 P(k) via k_fs 2026–28 | Future prediction |
| k_fs | free-streaming (m_φ) | 0.493 h/Mpc | DESI Y3/Euclid 2026–28 | Future prediction |
| r | φ⁻¹⁰ | 0.00813 | LiteBIRD (~2032) | Future prediction |

Genesis 5.12 commit: https://github.com/mikealmeida1721/SSEE_UNIFICADO (2026-01-28)
Zenodo DOI: https://doi.org/10.5281/zenodo.19679049

---

## Open blockers (path to PRD/PRL)

| Blocker | Description | Status |
|---|---|---|
| **B1** | plik_lite TTTEEE done (ΔBIC=−31.3); full CamSpec off-diagonal covariance pending | Weeks |
| **B2** | Full causal IS tensor perturbations (Hiscock-Lindblom 1985) — B-mode forecasts | Months |
| **B3** | MIRA: geometric derivation given in Paper 8 (disformal geodesic); full field-equation closure partial | Months |
| **B4** | arXiv endorsement — deferred by choice (journal-first strategy, see README) | External |
