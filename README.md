# SSEE — Structural Self-Energy Expansion

**A minimal-parameter dark energy framework derived from φ (golden ratio) and π. The background sector $(w_0, w_a, \Omega_\mathrm{DE}, \Omega_{m,\mathrm{dyn}})$ carries zero fitted dimensionless parameters; the framework rests on 4 postulates (D, S fundamentals + M, I auxiliary register-level) plus open problems tracked in [`OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md). Tested against DESI DR2 BAO, Planck 2018 CMB (TT+TE+EE+lensing), galaxy cluster masses, and large-scale structure growth.**

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
| m_φ (φ-DM mass) | 36.95 eV algebraic | gravitationally-produced scalar — observable via k_fs, not neutrino experiments | falsifiable |
| k_fs (free-streaming) | 0.659 h/Mpc algebraic | DESI Y3/Euclid P(k): 2026–2028 | falsifiable |

**Future falsifiers:** r = φ⁻¹⁰ = 0.00813 (LiteBIRD ~2032); k_fs cutoff in matter power spectrum (Euclid ~2028).

---

## Papers

| # | Title | Pages | Status | PDF |
|---|---|---|---|---|
| 1 | Minimal-Parameter Framework (φ, π → w₀, wₐ, EFT) + Predictive Register + Two-Ω_m Criterion | 28 | arXiv-ready | [docs/](docs/SSEE_Paper1_Framework.pdf) |
| 2 | Bayesian MCMC Validation — DESI DR2 + Planck + clusters | 26 | arXiv-ready | [docs/](docs/SSEE_Paper2_MCMC.pdf) |
| 3 | CMB Confrontation — Planck PR4 TT+TE+EE+lensing | 24 | arXiv-ready | [docs/](docs/SSEE_Paper3_CMB.pdf) |
| 4 | Algebraic Derivation of the CMB Background from φ and π | 16 | Preprint | [docs/](docs/SSEE_Paper4_ToE.pdf) |
| 5 | Israel-Stewart Causal Perturbation Theory — Gradient Stability, MIRA Origin, S₈ Tension Characterization | 25 | Preprint | [docs/](docs/SSEE_Paper5_IS.pdf) |
| 6 | φ-Dark Matter in SSEE-V3.6: Algebraic Mass Derivation and Resolution of the fσ₈ Tension | 24 | Preprint | [docs/](docs/SSEE_Paper6_phiDM.pdf) |
| 7 | Canonical EFT of SSEE-V3.6: Action, β_c = −AURA, and Bellini-Sawicki α-Functions | 16 | Preprint | [docs/](docs/SSEE_Paper7_EFT.pdf) |
| 8 | Strong Gravity Regime — Two-limit analysis (alt MOND-like vs canonical EFT B-S) | 20 | Preprint | [docs/](docs/SSEE_Paper8_StrongGravity.pdf) |
| 9 | Hubble Tension via Algebraic Screening Fraction: f_screen = αK/(3·MIRA) | 18 | Preprint | [docs/](docs/SSEE_Paper9_HubbleTension.pdf) |
| 10 | UV Completion of SSEE: K(X) = X/KAL + X²/M⁴, M = φ²·5^(1/4)·ρ_crit^(1/4) = 9.62 meV | 14 | Preprint | [docs/](docs/SSEE_Paper10_UVCompletion.pdf) |
| — | **Unified Journal Paper** (consolidation of Papers 1–7 + CLASS + MCMC Fase 4) | 19 | Journal submission candidate | [docs/](docs/SSEE_Unified_Journal.pdf) |
| ★ | **Sealed Journal** — consolidated late-universe dark-energy paper (φ → w₀, wₐ; closed-dictionary look-elsewhere; two-stage H₀; honest accounting of ~3 vs 6 parameters) | 10 | **Sealed — external-audit candidate** | [docs/](docs/SSEE_Sealed_Journal.pdf) |

---

## Repository structure

```
SSEE/
├── manuscript/                     # LaTeX sources, bibliographies, cover letters
│   ├── SSEE_Paper1_Framework.tex … SSEE_Paper10_UVCompletion.tex  # the 10 papers
│   ├── SSEE_EFT_section.tex        # EFT section (\input by Paper 1)
│   ├── SSEE_Unified_Journal.tex    # consolidated journal paper (Papers 1–10)
│   ├── SSEE_Endorser_Summary.tex   # 2-page arXiv endorser brief
│   ├── *.bib                       # ssee_paper3/4/5/6 + ssee_unified bibliographies
│   └── cover_letter_*.txt, abstracts_arXiv.txt
├── src/                            # 42 Python scripts — analysis & verification
│   ├── ssee_paper2_*.py … ssee_paper10_*.py  # per-paper analysis, MCMC, figures
│   ├── ssee_op1_*.py … ssee_op6_*.py         # OPEN_PROBLEMS resolution (OP-1..OP-6)
│   ├── ssee_paperB_DW.py, ssee_paperB_Nstar.py        # Paper B groundwork (baryogenesis, N_*)
│   ├── ssee_phase_c_dic.py, ssee_phase_d_savage_cv.py # Bayesian model-selection phases
│   ├── ssee_verify_rd.py, ssee_press_schechter.py, …  # standalone physics checks
│   ├── ssee_audit_consistency.py   # cross-paper parameter consistency audit
│   └── scratch/                    # exploratory scripts (not load-bearing)
├── class_ssee/                     # CLASS Boltzmann fork — SSEE .ini configs + plot scripts
├── data/                           # observational data (DESI DR2, Planck PR4, clusters)
├── results/                        # generated figures, tables, logs
├── notebooks/                      # Jupyter exploration
├── docs/                           # compiled PDFs — 10 papers + endorser + unified journal
├── submission_packages/            # arXiv-ready .tar.gz bundles per paper
├── archive/                        # superseded drafts (historical)
├── notes/                          # internal work-notes & attack plans (not load-bearing)
├── build_arxiv_packages.py         # regenerates submission_packages/
├── requirements.txt · environment.yml   # reproducible Python environment
├── OPEN_PROBLEMS.md                # physics gaps OP-1..OP-6 with status
├── AUDIT.md                        # reproducibility guide + known limitations
├── CHANGELOG.md · CITATION.cff · LICENSE
└── (eftcamb_ssee/ — EFTCAMB fork, not versioned: clone separately)
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
| k_crit / (H₀/c) | 0.456 < 1 | Sub-Hubble stability window |
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
| Σm_ν = R₂ × 0.960318 eV | 0.0690 eV | R₂ = Ω/(KAL·TRIAL) = 0.071875 |
| m_φ = Σm_ν × (Ω⁴_DNAV + AURA·KAL) | 36.95 eV | Forward-prediction — no fitting (multiplier 535.28 is a pure number) |
| α (Viel fit to particle/cold P(k) ratio) | 1.243 Mpc/h | CLASS output — not imposed |
| k_fs (free-streaming) | 0.659 h/Mpc | From m_φ, CLASS-derived |
| σ₈_eff (two-sector particle) | 0.742 | — |
| **S₈ (two-sector particle)** | **0.766** | **0.01σ KiDS-1000 (0.766±0.020) — resolves S₈ lensing tension** |
| Single-sector linear (cold source) | σ₈=0.820, S₈=0.847 | 3.9σ — open before two-sector |
| Mean fσ₈ tension (6 surveys) | 0.74σ | from 3.9σ baseline; ties ΛCDM |

> **Note:** With the forward-predicted particle (m_φ = 36.95 eV) and its own relic temperature, the two-sector model **resolves** the S₈ lensing tension: S₈ = 0.766 sits 0.01σ from KiDS-1000, down from the 3.9σ single-sector baseline. The fσ₈ (growth-rate) sector is unchanged and ties ΛCDM (0.74σ across six RSD surveys).

**Lyman-α compatibility:** φ-DM is a non-thermal condensate. The correct observable is k_fs=0.659 h/Mpc in the matter power spectrum (not m_φ directly); ΔP/P ≈ −f²(k/k_fs)² — quantified falsifiable prediction.

**φ-DM is a gravitationally-produced scalar** (explicit Lagrangian in Paper 6 §4): no Standard-Model portal, so it is *not* accessible to neutrino-mass experiments. The Dodelson-Widrow (sterile-neutrino) route is excluded — its mixing angle has no zero-parameter SSEE form.

**Falsifiable prediction:** k_fs = 0.659 h/Mpc — the free-streaming imprint of m_φ = 36.95 eV on the matter power spectrum — testable via the Lyman-α forest and DESI Y3/Euclid P(k) (2026–2028). Zero free parameters.

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
| Ωb h² | (π−φ)/(3Ω²) = 0.02242 | 0.02237 ± 0.00015 | 0.32σ |
| Ωc h² (IS) | KAL₀ × Ωb h² × n_s = 0.11926 | 0.1200 ± 0.0012 | **−0.6σ** |
| Y_p (BBN) | AlterBBN(Ωb h²=0.02242) = 0.2476 | 0.2449 ± 0.0040 | 0.7σ |
| δc | δc,EdS × n_s = 1.6284 | 1.6865 (EdS) | — |

**Press-Schechter halo-count enhancement** at z=10 (`src/ssee_press_schechter.py`):

| Halo mass | σ_M(z=10) | n_SSEE/n_ΛCDM |
|---|---|---|
| 10^11 M☉ | 1.010 | ×1.01 |
| 3×10^11 M☉ | 0.727 | ×1.04 |
| 10^12 M☉ | 0.506 | ×1.09 |

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
| ID | Problem | Status (2026-05-16) |
|----|---------|---------------------|
| OP-1 | Baryon density Ω_b h² | **Partial** — formula (π−φ)/(3Ω²) = 0.32σ Planck; ab-initio baryogenesis → Paper B |
| OP-2 | n_s = 1−φ⁻⁷ exponent | **Resolved** (conditional) — α-attractor universality + N_*=2φ⁷; new prediction r=φ⁻¹⁰ |
| OP-3 | UV-IR separability | **Resolved** — EFT hierarchy (H₀/M)²≈10⁻⁶²; KALeff=φ²√(5/2) unique |
| OP-4 | Solar Vainshtein radius | **Resolved** — k-mouflage (not Galileon) + αB=αM=αT=0 EFT suppression |
| OP-5 | S₈ weak-lensing tension | **Partial** — HMcode-2020 → S₈=0.758 (0.06σ DES); full N-body → future |
| OP-6 | Screening form (mult. vs add.) | **Resolved** — separate-universe k-essence + identity 1+w₀=Ω_m |

**Foundational postulates (Paper 1 §2.4, Postulates D & S):** the dimensional scale of
H₀ is an explicit anchor input (zero *dimensionless* fitted parameters, like ΛCDM); the
saturation correspondence (Ω_DE, w₀) = (s, −s) interpolates matter ↔ de Sitter. These are
pre-registered axioms, not derived theorems — see OPEN_PROBLEMS.md and Paper 1.

---

## Provenance

- **Genesis 5.12** (Zenodo DOI: [10.5281/zenodo.19679049](https://doi.org/10.5281/zenodo.19679049)): documents that MIRA = (3φ+π)/4 was defined on 2026-01-28, prior to Planck 2018 Ωm comparison.

---

## Roadmap

**Status (2026-05-17):** all 10 papers complete and compile clean (0 LaTeX errors,
0 orphan bibitems, 0 undefined citations). Foundational postulates D & S formalised
in Paper 1; OPEN_PROBLEMS OP-2/3/4/6 resolved, OP-1/5 partial. Full development
history in [CHANGELOG.md](CHANGELOG.md).

**Done**
- [x] Papers 1–10 — algebraic framework, Bayesian MCMC, CMB confrontation, algebraic
      CMB derivation, IS causal perturbations, φ-DM two-sector, canonical EFT,
      strong-gravity regime, Hubble-tension screening, UV completion
- [x] CLASS Boltzmann validation — MIRA necessity (RMS 1.4% vs 31.5%), σ₈, IS viscosity
- [x] Multi-probe MCMC (DESI DR2 + Planck + clusters) — mean tension 0.36σ
- [x] Bibliography brought to JCAP/PRD standard — all papers 36–42 refs, 0 orphans
- [x] OPEN_PROBLEMS OP-2/3/4/6 resolved; OP-1/5 partial
- [x] Paper 1 Postulates D (dimensional anchor) & S (saturation correspondence)
- [x] Hostile-referee overclaim sweep across all 10 papers
- [x] Zenodo v1 — Papers 1–7 archived (DOI 10.5281/zenodo.20093447)

**Pending**
- [ ] External audits → green light
- [ ] Zenodo v2 — archive Papers 1–10 + OPEN_PROBLEMS.md + CLASS scripts
- [ ] Journal submission to JCAP (after external green light)
- [ ] Paper B — ab-initio baryogenesis (OP-1 closure) + φ-DM relic abundance
- [ ] OP-5 closure — full N-body S₈ (BAHAMAS / IllustrisTNG-SSEE)

---

## Author

Mike Edison Almeida Vallejo — mike.almeida1721@gmail.com
ORCID: [0009-0008-2195-7836](https://orcid.org/0009-0008-2195-7836)

## License

Apache 2.0 — see [LICENSE](LICENSE).
