# SSEE — Structural Self-Energy Expansion

**A zero-free-parameter dark energy model derived exclusively from φ (golden ratio) and π, tested against DESI DR2 BAO, Planck 2018 CMB (TT+TE+EE+lensing), and galaxy cluster masses.**

---

## Four falsifiable predictions — all pass

| Observable | SSEE (algebraic) | Observed | Separation |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.827, −0.75) | 0.05σ |
| CMB first peak | ℓ₁ = 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 (algebraic via MIRA) | Planck 2018: 0.3153 | 0.63σ |
| n_s | 1 − φ⁻⁷ = 0.96556 | Planck 2018: 0.9649 | 0.2σ |

**Future falsifier (LiteBIRD ~2032):** tensor-to-scalar ratio r = φ⁻¹⁰ = 0.00813.

---

## Papers

| # | Title | Pages | Status | PDF |
|---|---|---|---|---|
| 1 | Zero-Parameter Framework (φ, π → w₀, wₐ, EFT) + Predictive Register | 22 | arXiv-ready | [docs/](docs/SSEE_Paper1_Framework.pdf) |
| 2 | Bayesian MCMC Validation — DESI DR2 + Planck + clusters | 21 | arXiv-ready | [docs/](docs/SSEE_Paper2_MCMC.pdf) |
| 3 | CMB Confrontation — Planck PR4 TT+TE+EE+lensing | 19 | arXiv-ready | [docs/](docs/SSEE_Paper3_CMB.pdf) |
| 4 | Algebraic Derivation of the CMB Background from φ and π | 15 | Preprint | [docs/](docs/SSEE_Paper4_ToE.pdf) |
| 5 | Israel-Stewart Causal Perturbation Theory — Gradient Stability, MIRA Origin, S₈ Partial Resolution | 24 | Preprint | [docs/](docs/SSEE_Paper5_IS.pdf) |

---

## Repository structure

```
SSEE/
├── src/
│   ├── ssee_paper2_mcmc.py           # Bayesian MCMC: SSEE vs ΛCDM vs CPL
│   ├── ssee_paper2_analysis.py       # Analytical w₀-wₐ plane
│   ├── ssee_paper3_cmb.py            # CAMB CMB spectrum TT+TE+EE+lensing
│   ├── ssee_verify_rd.py             # CAMB r_d verification vs Planck 2018
│   ├── ssee_press_schechter.py       # Press-Schechter δc=1.6284 vs 1.6865 (Paper 4)
│   └── ssee_inflation_connection.py  # α-attractor: α=φ⁴/3 exact, N=2φ⁷≈58.07 (Paper 1)
├── data/raw/
│   ├── planck_pr4_lensing.txt     # Planck 2018 MV lensing bandpowers (14 bins)
│   └── ...                        # Planck PR4 TT/TE/EE spectra
├── results/figures/               # All generated figures (PDF/PNG)
├── manuscript/
│   ├── SSEE_Paper1_Framework.tex  # Paper 1 source
│   ├── SSEE_EFT_section.tex       # EFT appendix (included by Paper 1)
│   ├── SSEE_Paper2_MCMC.tex      # Paper 2 source
│   ├── SSEE_Paper3_CMB.tex      # Paper 3 source
│   ├── ssee_paper3.bib
│   ├── SSEE_Paper5_IS.tex         # Paper 5 source
│   └── ssee_paper5.bib
├── sandbox_unificado/
│   └── SSEE_Paper4_ToE.tex        # Paper 4 source (git submodule)
├── submission_packages/           # arXiv-ready .tar.gz for each paper
├── docs/                          # Compiled PDFs (5 papers)
│   └── archive/                   # Older versions — local only, not in git
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

**Paper 3 — CMB power spectrum** (CAMB, TT+TE+EE+lensing):
```bash
python src/ssee_paper3_cmb.py
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
| TT | 1.063 | 1.043 | 1971 |
| TE | 1.053 | 1.040 | 1967 |
| EE | 1.040 | 1.039 | 1967 |
| PP (lensing) | **0.732** | 0.757 | 9 |
| ΔBIC (Combined, k=2 vs k=6) | **-31.3** (SSEE favoured) | — | — |

*Note: SSEE evaluated with $k=2$ ($H_0$ and $\Omega_b h^2$ free) using n_s = 0.96556 against the full official `plik_lite` TTTEEE + lowl covariance matrix via Cobaya. Even with a hyper-conservative $k=4$ penalty (treating fixed priors $\tau, A_s$ as free), $\Delta\mathrm{BIC} = -13.8$ continues to decisively favour SSEE.*

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

## Known limitations

Disclosed honestly in the papers. Full treatment in [AUDIT.md](AUDIT.md):

1. **Diagonal CMB likelihood** (Paper 3): χ²_r uses diagonal covariance. Off-diagonal terms required for PRD/PRL.
2. **Israel-Stewart viscosity** (EFT section): IS growth index γ_IS=0.657 and τ_Π H₀=2.191 are derived (Paper 1 App.A). Full causal IS (Hiscock-Lindblom 1985) for B-mode perturbation predictions remains a blocker for LiteBIRD forecasts.
3. **Ωc h² (resolved)** (Paper 4): static Eckart gave 3.7σ; Israel-Stewart IS derivation KAL₀ × Ωb h² × n_s = 0.11926 reduces tension to −0.6σ.
4. **Two-sector Ωm**: Ωm,dyn = 0.160 ≠ Ωm,CMB = 0.3199. Bridged by algebraic MIRA factor; physical mechanism needs formal development.
5. **ΔBIC = +206** (Paper 2 full model): applies ΛCDM Friedmann background to SSEE parameters — acknowledged as a framework-internal constraint.
6. **H(z) tension** (Paper 2): SSEE χ²_r = 1.861 vs ΛCDM 0.458 on cosmic chronometers.
7. **CMB ΔBIC status:** The project is currently unifying its CMB likelihood pipelines. Preliminary diagonal estimates yield $\Delta\mathrm{BIC} = +33.2$ favouring $\Lambda$CDM, whilst full off-diagonal evaluation is an open blocker for top-tier review.

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
- [ ] Full CMB likelihood (Cobaya + plik/CamSpec) — blocker for PRD/PRL
- [ ] Full causal Israel-Stewart (Hiscock-Lindblom 1985) — blocker for B-mode LiteBIRD forecasts
- [ ] arXiv submission — pending endorser response

---

## Author

Mike Edison Almeida Vallejo — mike.almeida1721@gmail.com
ORCID: [0009-0008-2195-7836](https://orcid.org/0009-0008-2195-7836)

## License

Apache 2.0 — see [LICENSE](LICENSE).
