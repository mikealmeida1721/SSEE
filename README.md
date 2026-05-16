# SSEE — Structural Self-Energy Expansion

**A zero-free-parameter dark energy model derived exclusively from φ (golden ratio) and π, tested against DESI DR2 BAO, Planck 2018 CMB (TT+TE+EE+lensing), galaxy cluster masses, and large-scale structure growth.**

---

## Falsifiable predictions — current status

| Observable | SSEE (algebraic) | Observed | Separation |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.838±0.057, −0.631±0.226) | 0.05σ / 0.17σ |
| CMB first peak | ℓ₁ = 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 (algebraic via MIRA) | Planck 2018: 0.3153 | 0.63σ |
| n_s | 1 − φ⁻⁷ = 0.96556 | Planck 2018: 0.9649 | 0.2σ |
| αT (GW speed) | 0 exact | GW170817: \|αT\| < 10⁻¹⁵ | exact match |
| αK (kineticity, z=0) | 0.4033 algebraic | Euclid forecast: < 0.1 | testable 2026–2028 |
| m_φ (φ-DM mass) | 5.71 eV algebraic | KATRIN/PTOLEMY: 2027–2030 | falsifiable |
| k_fs (free-streaming) | 0.493 h/Mpc algebraic | DESI Y3/Euclid P(k): 2026–2028 | falsifiable |

**Future falsifiers:** r = φ⁻¹⁰ = 0.00813 (LiteBIRD ~2032); k_fs cutoff in matter power spectrum (Euclid ~2028).

---

## Papers

| # | Title | Pages | Status | PDF |
|---|---|---|---|---|
| 1 | Zero-Parameter Framework (φ, π → w₀, wₐ, EFT) + Predictive Register | 22 | arXiv-ready | [docs/](docs/SSEE_Paper1_Framework.pdf) |
| 2 | Bayesian MCMC Validation — DESI DR2 + Planck + clusters | 21 | arXiv-ready | [docs/](docs/SSEE_Paper2_MCMC.pdf) |
| 3 | CMB Confrontation — Planck PR4 TT+TE+EE+lensing | 18 | arXiv-ready | [docs/](docs/SSEE_Paper3_CMB.pdf) |
| 4 | Algebraic Derivation of the CMB Background from φ and π | 15 | Preprint | [docs/](docs/SSEE_Paper4_ToE.pdf) |
| 5 | Israel-Stewart Causal Perturbation Theory — Gradient Stability, MIRA Origin, S₈ Tension Characterization | 24 | Preprint | [docs/](docs/SSEE_Paper5_IS.pdf) |
| 6 | φ-Dark Matter in SSEE-V3.6: Algebraic Mass Derivation and Resolution of the fσ₈ Tension | 15 | Preprint | [docs/](docs/SSEE_Paper6_phiDM.pdf) |
| 7 | Canonical EFT of SSEE-V3.6: Action, β_c = −AURA, and Bellini-Sawicki α-Functions | 14 | Preprint | [docs/](docs/SSEE_Paper7_EFT.pdf) |
| 8 | Strong Gravity Regime — Disformal Geodesics, MIRA Lensing Emergence, Vainshtein Screening | 13 | Preprint | [docs/](docs/SSEE_Paper8_StrongGravity.pdf) |
| 9 | Hubble Tension via Algebraic Screening Fraction: f_screen = αK/(3·MIRA) | 10 | Draft | [docs/](docs/SSEE_Paper9_HubbleTension.pdf) |
| 10 | UV Completion of SSEE: K(X) = X/KAL + X²/M⁴, M = φ²·5^(1/4)·ρ_crit^(1/4) = 8.81 meV | ~15 | Preprint | [docs/](docs/SSEE_Paper10_UVCompletion.pdf) |
| — | **Unified Journal Paper** (consolidation of Papers 1–7 + CLASS + MCMC Fase 4) | 14 | Journal submission candidate | [manuscript/](manuscript/SSEE_Unified_Journal.tex) |

---

## Repository structure

```
SSEE/
├── src/
│   ├── ssee_paper2_mcmc.py               # Bayesian MCMC: SSEE vs ΛCDM vs CPL (Paper 2)
│   ├── ssee_paper2_analysis.py           # Analytical w₀-wₐ plane (Paper 2)
│   ├── ssee_paper3_cmb.py                # CAMB CMB spectrum TT+TE+EE+lensing (Paper 3)
│   ├── ssee_verify_rd.py                 # CAMB r_d verification vs Planck 2018
│   ├── ssee_press_schechter.py           # Press-Schechter δc=1.6284 vs 1.6865 (Paper 4)
│   ├── ssee_inflation_connection.py      # α-attractor: α=φ⁴/3 exact, N=2φ⁷≈58.07 (Paper 1)
│   ├── ssee_paper5_IS_perturbations.py   # IS causal perturbations: γ_IS, σ₈, fσ₈ (Paper 5)
│   ├── ssee_paper6_verification.py       # φ-DM two-sector: fσ₈ tensions, σ₈_eff (Paper 6)
│   ├── ssee_paper6_sterile_neutrino.py   # k_fs derivation, Dodelson-Widrow, m_φ scan (Paper 6)
│   ├── ssee_paper6_mcmc.py               # Full MCMC 3-param φ-DM: 64 walkers × 15000 steps (Paper 6)
│   ├── ssee_eft_verification.py          # β_c=-AURA plateau test, 8 ICs (Paper 7)
│   └── ssee_audit_consistency.py         # Cross-paper parameter consistency audit (all papers)
├── data/raw/
│   ├── planck_pr4_lensing.txt     # Planck 2018 MV lensing bandpowers (14 bins)
│   └── ...                        # Planck PR4 TT/TE/EE spectra
├── results/figures/               # All generated figures (PDF/PNG)
├── manuscript/
│   ├── SSEE_Paper1_Framework.tex  # Paper 1 source
│   ├── SSEE_EFT_section.tex       # EFT appendix (included by Paper 1)
│   ├── SSEE_Paper2_MCMC.tex       # Paper 2 source
│   ├── SSEE_Paper3_CMB.tex        # Paper 3 source + ssee_paper3.bib
│   ├── SSEE_Paper5_IS.tex         # Paper 5 source + ssee_paper5.bib
│   ├── SSEE_Paper4_ToE.tex        # Paper 4 source + SSEE_Paper4.bib
│   ├── SSEE_Paper6_phiDM.tex      # Paper 6 source + ssee_paper6.bib
│   ├── SSEE_Paper7_EFT.tex        # Paper 7 source
│   ├── SSEE_Unified_Journal.tex   # Unified journal paper (Papers 1–7 + CLASS + MCMC Fase 4)
│   ├── ssee_unified.bib           # Bibliography for unified journal
│   └── SSEE_Endorser_Summary.tex  # 2-page arXiv endorser brief
├── class_ssee/                    # CLASS Boltzmann code (fork of class_public)
│   ├── ssee_v36.ini               # SSEE parameters — MIRA sector (Ω_m=0.3199)
│   ├── ssee_v36_nomira.ini        # SSEE parameters — dynamic sector only (Ω_m=0.160)
│   ├── ssee_v36_twosector.ini     # φ-DM two-sector (ncdm m=5.71 eV)
│   ├── ssee_v36_IS.ini            # IS viscosity (cs2_fld=0.001)
│   ├── plot_ssee_cmb.py           # CMB TT comparison CLASS vs Planck
│   ├── plot_ssee_twosector_pk.py  # P(k) two-sector + WDM transfer
│   ├── plot_IS_viscosity.py       # IS viscosity effect on σ₈
│   └── calibrate_wdm_alpha.py     # WDM α calibration (sigma8 top-hat correct)
├── submission_packages/           # arXiv-ready .tar.gz for each paper
├── docs/                          # Compiled PDFs (7 papers + endorser brief)
└── AUDIT.md                       # Full reproducibility guide + known limitations
```

---

## How to reproduce

```bash
pip install camb emcee scipy numpy matplotlib
```

**Paper 2 — MCMC validation** (100-walker, N_eff ≈ 637,500 for SSEE):
```bash
python src/ssee_paper2_mcmc.py
```

**Paper 3 — CMB power spectrum** (Cobaya MCMC, TT+TE+EE+lowl):
You can set the `COBAYA_PACKAGES_PATH` environment variable to point to your local Planck 2018 data:
```bash
export COBAYA_PACKAGES_PATH=/path/to/your/cobaya_packages
python src/ssee_paper3_cobaya_unified.py
```

See [AUDIT.md](AUDIT.md) for expected outputs and known limitations.

---

## Key results

### Paper 2 (DESI DR2 + Planck 2018 + clusters)

| Metric | Value |
|---|---|
| χ²_2D (w₀-wₐ vs DESI DR2) | 0.080 → 0.05σ |
| χ²_r clusters (7 clusters, analytic) | 0.126 |
| χ²_r clusters (4 clusters, MCMC) | 0.122 |
| H₀ SSEE (MCMC) | 66.75 ⁺⁰·⁴⁴₋₀.₄₄ km/s/Mpc |
| ΔBIC (dynamic sector, k=1 vs ΛCDM k=3) | −5.55 (SSEE favoured) |
| ΔBIC (full background, k=0 vs ΛCDM k=6) | +206 (ΛCDM favoured — framework penalty) |

### Paper 3 (Planck PR4 CMB)

| Spectrum | SSEE χ²_r | ΛCDM χ²_r | N |
|---|---|---|---|
| TT | 1.062 | 1.043 | 1971 |
| TE | 1.053 | 1.040 | 1967 |
| EE | 1.040 | 1.039 | 1967 |
| PP (lensing) | **0.730** | 0.757 | 9 |
| ΔBIC (Combined, k=2 vs k=6) | **-31.3** (SSEE favoured) | — | — |

*Methodology Note on parameter counting ($k$): SSEE is evaluated with $k=2$ ($H_0$ and $\Omega_b h^2$ free) against the full official `plik_lite` TTTEEE + lowl covariance matrix via Cobaya. This is compared to the standard $\Lambda$CDM model which uses $k=6$. Even with a hyper-conservative $k=4$ penalty (treating fixed priors $\tau, A_s$ as free), $\Delta\mathrm{BIC} = -13.8$ continues to favour SSEE.*

**Growth structure (Paper 3 §5.4–5.5):**

| Metric | SSEE | ΛCDM |
|---|---|---|
| Growth index γ_IS (Paper 1 App.A) | 0.657 ± 0.002 | 0.55 |
| S8 (IS numerical) | 0.818 | 0.830 |
| fσ8 χ²/N (13 RSD surveys) | **0.524** | 0.596 |

### Paper 5 (Israel-Stewart Causal Perturbations)

| Result | Value | Status |
|---|---|---|
| c²_s,eff | 0 (exact algebraic) | Q1: all modes stable |
| k_crit / (H₀/c) | 0.498 < 1 | Sub-Hubble stability window |
| MIRA (numerical, k≥10) | 0.989 ± 0.017 | Background IS origin confirmed |
| γ_IS | 0.554 ± 0.001 | ≈ γ_ΛCDM = 0.55 |
| G = D₁_SSEE/D₁_ΛCDM | 0.866 | 13.4% growth suppression |
| σ₈_SSEE | 0.702 ± 0.005 | — |
| **S₈_SSEE** | **0.725 ± 0.005** | **1.96σ KiDS / 2.84σ DES (IS sector explored, tension remains)** |
| fσ₈(z=0.5) | 0.341 | 2.67σ vs BOSS (structural tension, Ω_m,dyn=0.160) |

**Diagnostic:** SSEE predicts *amplitude* (S₈) lower than Planck but tension with DES/KiDS remains.
Split points to kinetic braiding (αB in Bellini-Sawicki) as the natural Paper 6 extension.

### Paper 6 (φ-Dark Matter, Two-Sector Model)

| Result | Value | Status |
|---|---|---|
| Ω_CDM | 0.160050 | Active at all k |
| Ω_φDM = (MIRA−1)×Ω_m,dyn | 0.159878 | Active for k < k_fs only |
| Ω_total (two-sector) | 0.319928 ≈ Ω_m,CMB | Zero-parameter unification |
| m_φ = Σm_ν × H₀^alg | 5.71 eV | Algebraic — no fitting |
| k_fs (Dodelson-Widrow) | 0.493 h/Mpc | From m_φ algebraic |
| T_WDM(k=0.125 h/Mpc) | 0.8175 | Soft cutoff (scalar condensate, not thermal WDM) |
| σ₈_eff | **0.737** (CLASS+WDM, α=1.6561 h/Mpc) | 0.00σ KiDS-1000 |
| S₈_eff | **0.761** | 2.29σ DES (KiDS-DES internal tension) |
| **Mean fσ₈ tension (6 surveys)** | **0.50σ** | **Resolved from 2.56σ (single-sector baseline, corrected); fσ₈ only** |
| ΛCDM fσ₈ mean tension | 0.51σ | Reference |

> **Note:** The fσ₈ (growth rate) tension is resolved at 0.50σ. The S₈ (weak lensing amplitude) tension remains at ~2.3σ DES / ~3σ KiDS; it is not claimed to be solved.

**Lyman-α compatibility:** φ-DM is a non-thermal condensate. f_φ≈0.50 → effective bound ~0.1–0.5 keV; k_fs=0.493 h/Mpc is the correct observable (not m_φ directly); ΔP/P ≈ −f²(k/k_fs)² — quantified falsifiable prediction.

**Falsifiable predictions:** m_φ = 5.71 eV testable by KATRIN/PTOLEMY (2027–2030); k_fs = 0.493 h/Mpc via DESI Y3/Euclid P(k) power spectrum (2026–2028). Zero free parameters.

### Paper 7 (Canonical EFT)

| Result | Value | Status |
|---|---|---|
| Action | S = ∫d⁴x√(−g)[−X + V₀e^{β_c φ}] | Canonical, minimal coupling |
| β_c | −AURA = −3.9978 | Algebraic exact |
| β_c (plateau test, 8 ICs) | −3.98991 ± 0.00001 | Δ = 0.199% (systematic: Ω_b+Ω_r excluded) |
| αT | 0 exact | GW170817 \|αT\| < 10⁻¹⁵ satisfied ✓ |
| αM | 0 exact | Euclid forecast < 0.05 satisfied ✓ |
| αB | 0 exact | Euclid forecast < 0.05 satisfied ✓ |
| αK(z=0) | 3·Ω_DE·Ω_m,dyn = 0.4033 | Algebraic (Euclid will constrain < 0.1) |
| G₂_s (running) | 0.979 | Two-sector unification ✓ |

### CLASS Boltzmann Validation (Fases 1–3)

| Test | SSEE+MIRA | SSEE sin MIRA | ΛCDM | Significance |
|---|---|---|---|---|
| CMB peak 1 (ℓ) | **220** | 240 | 221 | MIRA necessary |
| CMB peak 2 (ℓ) | **535** | 597 | 537 | MIRA necessary |
| CMB peak 3 (ℓ) | **810** | 922 | 814 | MIRA necessary |
| RMS vs ΛCDM | **1.4%** | 31.5% | — | 22× degradation without MIRA |
| σ₈ (MIRA+WDM, CLASS) | **0.737** | — | ~0.82 | 0.00σ KiDS-1000 |
| IS cs² effect on σ₈ | 0.03% | — | — | Negligible ✓ |
| G = D₁_SSEE/D₁_ΛCDM (ODE) | **0.866** | — | 1.0 | 13.4% suppression confirmed |

*CLASS confirms MIRA is physically necessary: without it, all three CMB peaks shift ~10% and RMS error rises 22×.*

### MCMC Fase 4 (Multi-probe background)

| Parameter | Algebraic | Posterior (median ±1σ) | Tension |
|---|---|---|---|
| H₀ (km/s/Mpc) | 67.96 | 67.34 ± 0.77 | **0.81σ** ✅ |
| w₀ | −0.840 | −0.812 ± 0.046 | **0.61σ** ✅ |
| wₐ | −0.670 | −0.659 ± 0.113 | **0.10σ** ✅ |
| Ωm,dyn | 0.160 | 0.163 ± 0.012 | **0.25σ** ✅ |
| r_d (Mpc) | 147.156 | 147.44 ± 0.89 | **0.32σ** ✅ |
| **Mean tension (5 params)** | — | — | **0.36σ** ✅ |

*100 walkers × 11,000 steps; N_eff = 4,402. All algebraic predictions within 1σ of joint posterior.*

### Paper 4 (Algebraic ToE)

| Observable | SSEE algebraic | Planck 2018 | Tension |
|---|---|---|---|
| n_s | 1 − φ⁻⁷ = 0.96556 | 0.9649 ± 0.0042 | 0.2σ |
| H₀ | 3(φ+π)² = 67.96 km/s/Mpc | 67.36 ± 0.54 | 1.1σ |
| Ωm | (π−φ)/(π+φ) = 0.3201 | 0.3153 ± 0.0073 | 0.66σ |
| Ωb h² | 3(π−φ)/200 = 0.02285 | 0.02237 ± 0.00015 | 3.2σ |
| Ωc h² (IS) | KAL₀ × Ωb h² × n_s = 0.11926 | 0.1200 ± 0.0012 | **−0.6σ** |
| Y_p (BBN) | AlterBBN(Ωb h²=0.02285) = 0.2476 | 0.2449 ± 0.0040 | 0.7σ |
| δc | δc,EdS × n_s = 1.6284 | 1.6865 (EdS) | — |

**Press-Schechter halo-count enhancement** at z=10 (`src/ssee_press_schechter.py`):

| Halo mass | σ_M(z=10) | n_SSEE/n_ΛCDM |
|---|---|---|
| 10^11 M☉ | 1.010 | ×1.06 |
| 3×10^11 M☉ | 0.727 | ×1.16 |
| 10^12 M☉ | 0.506 | **×1.41** |

*Partially alleviates the JWST z>10 galaxy excess (ΛCDM deficit: ~10–100×); residual tension points to IS perturbation physics.*

**Inflationary embedding (Paper 1 App.A §A.5):**

| Quantity | SSEE value | Relation |
|---|---|---|
| α-attractor parameter | φ⁴/3 = 2.285 | exact (0.00e+00) |
| e-folds N | 2φ⁷ ≈ 58.07 | from n_s = 1 − φ⁻⁷ |
| Kähler curvature R | −φ⁻⁴ ≈ −0.146 | R = −2/(3α) |
| IS growth index γ_IS | 0.657 ± 0.002 | Paper 1 App.A |
| IS relaxation time τ_Π H₀ | KAL₀/(3Ω_DE) ≈ 2.191 | algebraic |
| MIRA algebraic identity | (3φ+π)/4 = 1.9989 | exact (0.00e+00) |

---

## Known limitations and open problems

Disclosed honestly in the papers. Editorial limitations in [AUDIT.md](AUDIT.md).
**Physics gaps that future work must address:** [OPEN_PROBLEMS.md](OPEN_PROBLEMS.md) (6 items, OP-1..OP-6).

### Editorial/pipeline limitations
1. **Diagonal CMB likelihood** (Paper 3): χ²_r uses diagonal covariance. Off-diagonal terms required for PRD/PRL.
2. **Full causal IS** (Paper 5): IS growth index γ_IS=0.657 is derived analytically. Full Hiscock-Lindblom 1985 treatment for B-mode predictions remains a blocker for LiteBIRD forecasts.
3. **ΔBIC = +206** (Paper 2 full model): applies ΛCDM Friedmann background to SSEE parameters — acknowledged as a framework-internal constraint.
4. **H(z) tension** (Paper 2): SSEE χ²_r = 1.861 vs ΛCDM 0.458 on cosmic chronometers.
5. **CMB ΔBIC status:** Diagonal approximation (TT+TE+EE+PP, k=2 vs k=6) yields $\Delta\mathrm{BIC}=+31.1$ favouring $\Lambda$CDM. Full Cobaya plik\_lite (TTTEEE) gives $\Delta\mathrm{BIC}=-31.3$ decisively favouring SSEE.

### Open physics problems (see [OPEN_PROBLEMS.md](OPEN_PROBLEMS.md))
| ID | Problem | Severity |
|----|---------|----------|
| OP-1 | Factor 200 in Ω_b h² = 3(π−φ)/200 — no first-principles derivation | High |
| OP-2 | n_s = 1−φ⁻⁷ — exponent 7 motivated by hierarchy count, not V(φ) | High |
| OP-3 | UV-IR separability conjecture (φ/π factorization) — unproven at loop level | Medium |
| OP-4 | Vainshtein radius for Sun r_V > r_Hubble — EFT internal inconsistency | High |
| OP-5 | S₈ weak-lensing tension 2.29σ DES (fσ₈ resolved in Paper 6 at 0.50σ) | Medium |
| OP-6 | Screening form (multiplicative vs additive) not derived from EFT Lagrangian | Medium |

---

## Provenance

- **Genesis 5.12** (Zenodo DOI: [10.5281/zenodo.19679049](https://doi.org/10.5281/zenodo.19679049)): documents that MIRA = (3φ+π)/4 was defined on 2026-01-28, prior to Planck 2018 Ωm comparison.

---

## Roadmap

- [x] Paper 1: Algebraic framework (φ, π → w₀, wₐ, EFT) — compiles clean
- [x] Paper 2: MCMC validation (100-walker, DESI 13×13 block-diagonal covariance)
- [x] Paper 3: CMB TT+TE+EE+lensing (Planck PR4)
- [x] Paper 4: Algebraic derivation of CMB observables from φ and π only
- [x] Genesis 5.12 — Zenodo DOI: 10.5281/zenodo.19679049
- [x] IS growth index γ_IS=0.657 ± 0.002, τ_Π H₀=2.191, S8=0.837 (Paper 1 App.A)
- [x] MIRA algebraic identity: (3φ+π)/4 = 1.9989 exact — Paper 3 §2.4
- [x] α-attractor embedding: α=φ⁴/3 exact, N=2φ⁷≈58.07, R_Kähler=−φ⁻⁴ — Paper 1 App.A §A.5
- [x] Task 2B: Press-Schechter δc=1.6284 vs 1.6865 — quantitative table in Paper 4 §sec:deltac
- [x] Peer-review audit round 1 — all C1–C3, M1–M5, m1–m4 fixes applied across 4 papers
- [x] Zenodo manuscripts suite (v3) — DOI: 10.5281/zenodo.19932301 (all 4 PDFs archived)
- [x] arXiv endorsement emails sent — Dr. Shafieloo (KASI), Dr. Lee (SKKU), Dr. Di Valentino (Sheffield)
- [x] External audit round 2 (PRD/JCAP level): C1–C4, M1–M5, Mod1–Mod6 — all resolved
  - C1–C2: AURA/MIRA glosario corregido (AURA=β+φ, MIRA=AURA/2)
  - C3: H₀=3(φ+π)² reclasificado Type P† con footnote dimensional
  - M1: App D "Algebraic Prior Space" añadido a Paper 1 (ahora 22 pp.)
  - M5: Phantom crossing a*=0.761 documentado explícitamente (Paper 1 §2)
  - Mod2: Bullet Cluster hedge — masa proyectada, no fit κ(θ)
  - Mod3: c²_s<0 caveat en EFT — deferred full IS treatment a Paper 4
- [x] Paper 5: IS causal perturbation theory — Q1+Q2+Q3+fσ₈, 24 pp
- [x] Paper 6: φ-DM two-sector model — m_φ=5.71 eV algebraic, fσ₈ 2.56σ→0.50σ (corrected baseline), 15 pp (incl. Lyman-α defense)
- [x] Paper 7: Canonical EFT — β_c=−AURA exact, αT=αM=αB=0, αK=0.4033 algebraic, 14 pp
- [x] Cross-paper parameter consistency audit — no critical inconsistencies (`src/ssee_audit_consistency.py`)
- [x] Paper 6 MCMC — 64 walkers × 15000 steps: Ω_φDM=0.161±0.011, ΔBIC=−12.1, χ²_r=0.497 ✓
- [x] hi_class/CLASS cross-check — αK(0)=0.4033, Δ=0.005% vs algebraic prediction ✓
- [x] CLASS Fase 1A — MIRA test: SSEE sin MIRA → RMS 31.5% (MIRA physically necessary confirmed)
- [x] CLASS Fase 1B — P(k) two-sector: WDM α=1.6561 h/Mpc, σ₈_eff=0.737, S₈=0.761 ✓
- [x] CLASS Fase 2c — sigma8 top-hat fix: correct σ₈ integral with W²(kR) filter ✓
- [x] CLASS Fase 3 — IS viscosity: cs² effect 0.03%, G=0.866 confirmed ✓
- [x] MCMC Fase 4 — 5-param multi-probe: max tension 0.81σ, mean 0.36σ ✓
- [x] Unified journal paper — 884 lines, full audit, Irsic+2017 citation corrected ✓
- [x] Paper 8: Strong gravity regime — disformal geodesics, MIRA lensing emergence, Vainshtein (13 pp, commit d7f3b38) ✓
- [x] Paper 9: Hubble tension via f_screen=αK/(3·MIRA)=(π−φ)/Ω²=0.06725, H₀,local=72.86 km/s/Mpc (0.17σ SH0ES) — draft 10 pp ✓
- [x] Editorial audit 10-paper suite (2026-05-15) — overclaims, nomenclature, BIC corrections across Papers 1–9 ✓
- [x] Paper 4 nomenclature: mythological names (SOLAR, IGNIS, PYROS…) → algebraic symbols (φ, π, Ω, β, KAL, P_sc, K_v); bibliography 10→40 entries ✓
- [x] OPEN_PROBLEMS.md — 6 physics gaps catalogued (OP-1..OP-6) with severity and resolution paths ✓
- [x] Paper 10: UV completion K(X)=X/KAL+X²/M⁴; M⁴=5φ⁸ρ_crit exact; M=8.81 meV=Λ_SSEE; αK_full=0.41691; Conditional Theorem C.1; H₀^UV=73.040 km/s/Mpc ✓
- [ ] OP-4: Recompute solar r_V using Λ_SSEE=M=8.81 meV from Paper 10 (k-mouflage vs Galileon) — in progress
- [ ] OP-6: Derive screening form H₀^local from SSEE Lagrangian perturbatively (after OP-4)
- [ ] OP-3: Prove Postulate C.1 unconditionally — Jacobian ∂φ/∂χ|_transition = φ²√(5/2) (after OP-6)
- [ ] OP-1+OP-2: Derive Ω_b h² factor 200 and n_s exponent 7 from UV structure (after OP-3)
- [ ] Zenodo v2 — archive Papers 1–10 + OPEN_PROBLEMS.md + CLASS scripts with DOI chaining
- [ ] Journal submission to JCAP (pending Zenodo v2 + endorser response)

---

## Author

Mike Edison Almeida Vallejo — mike.almeida1721@gmail.com
ORCID: [0009-0008-2195-7836](https://orcid.org/0009-0008-2195-7836)

## License

Apache 2.0 — see [LICENSE](LICENSE).
