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

**Future falsifier (LiteBIRD ~2032):** tensor-to-scalar ratio r = φ⁻¹⁰ = 0.0083.

---

## Papers

| # | Title | Status | PDF |
|---|---|---|---|
| 1 | Zero-Parameter Framework (φ, π → w₀, wₐ, EFT) | arXiv-ready | [docs/](docs/SSEE_Paper1_Framework.pdf) |
| 2 | Bayesian MCMC Validation — DESI DR2 + Planck + clusters | arXiv-ready | [docs/](docs/SSEE_Paper2_draft.pdf) |
| 3 | CMB Confrontation — Planck PR4 TT+TE+EE+lensing | arXiv-ready | [docs/](docs/SSEE_Paper3_draft.pdf) |
| 4 | Theory of Everything — algebraic derivation from φ and π only | Preprint | [docs/](docs/SSEE_Paper4_ToE.pdf) |

---

## Repository structure

```
SSEE/
├── src/
│   ├── ssee_paper2_mcmc.py        # Bayesian MCMC: SSEE vs ΛCDM vs CPL
│   ├── ssee_paper2_analysis.py    # Analytical w₀-wₐ plane
│   └── ssee_paper3_cmb.py         # CAMB CMB spectrum TT+TE+EE+lensing
├── data/raw/
│   ├── planck_pr4_lensing.txt     # Planck 2018 MV lensing bandpowers (14 bins)
│   └── ...                        # Planck PR4 TT/TE/EE spectra
├── results/figures/               # All generated figures (PDF/PNG)
├── manuscript/
│   ├── SSEE_Paper1_Framework.tex  # Paper 1 source
│   ├── SSEE_EFT_section.tex       # EFT appendix (included by Paper 1)
│   ├── SSEE_Paper2_draft.tex      # Paper 2 source
│   ├── SSEE_Paper3_draft.tex      # Paper 3 source
│   └── ssee_paper3.bib
├── sandbox_unificado/
│   └── SSEE_Paper4_ToE.tex        # Paper 4 source (git submodule)
├── submission_packages/           # arXiv-ready .tar.gz for each paper
├── docs/                          # Compiled PDFs (4 papers only)
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
| χ²_r clusters (4 clusters, IGIMF-corrected) | 0.122 |
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
| ΔBIC (TT, k=0 vs k=6) | **−6.9** (SSEE favoured) | — | — |

*Note: χ²_r computed with diagonal covariance. Full off-diagonal likelihood (Cobaya/plik) is a known open item — see AUDIT.md.*

### Paper 4 (Algebraic ToE)

| Observable | SSEE algebraic | Planck 2018 | Tension |
|---|---|---|---|
| n_s | 1 − φ⁻⁷ = 0.96556 | 0.9649 ± 0.0042 | 0.2σ |
| H₀ | 3(φ+π)² = 67.96 km/s/Mpc | 67.36 ± 0.54 | 1.1σ |
| Ωm | (π−φ)/(π+φ) = 0.3201 | 0.3153 ± 0.0073 | 0.66σ |
| Ωb h² | (π−φ)/[6(φ+π)⁴] = 0.02260 | 0.02237 ± 0.00015 | 1.5σ |
| Ωc h² | (φ+π)/[3(3φ+π)] = 0.1245 | 0.1200 ± 0.0012 | **3.7σ** |

*Ωc h² 3.7σ tension is the principal open problem of Paper 4.*

---

## Known limitations

Disclosed honestly in the papers. Full treatment in [AUDIT.md](AUDIT.md):

1. **Diagonal CMB likelihood** (Paper 3): χ²_r uses diagonal covariance. Off-diagonal terms required for PRD/PRL.
2. **Eckart viscosity** (EFT section): Eckart formulation violates causality (Hiscock-Lindblom 1985). Israel-Stewart formulation needed for perturbation predictions.
3. **Ωc h² 3.7σ** (Paper 4): cold dark matter density is the model's main quantitative open problem.
4. **Two-sector Ωm**: Ωm,dyn = 0.160 ≠ Ωm,CMB = 0.3199. Bridged by algebraic MIRA factor; physical mechanism needs formal development.
5. **ΔBIC = +206** (Paper 2 full model): applies ΛCDM Friedmann background to SSEE parameters — acknowledged as a framework-internal constraint.
6. **H(z) tension** (Paper 2): SSEE χ²_r = 1.861 vs ΛCDM 0.458 on cosmic chronometers.

---

## Provenance

- **Genesis 5.12** (Zenodo DOI: [10.5281/zenodo.19679049](https://doi.org/10.5281/zenodo.19679049)): documents that MIRA = (3φ+π)/4 was defined on 2026-01-28, prior to Planck 2018 Ωm comparison.

---

## Roadmap

- [x] Paper 1: Algebraic framework (φ, π → w₀, wₐ, EFT) — compiles clean
- [x] Paper 2: MCMC validation (100-walker, full DESI 13×13 covariance)
- [x] Paper 3: CMB TT+TE+EE+lensing (Planck PR4)
- [x] Paper 4: Algebraic derivation of CMB observables from φ and π only
- [x] Genesis 5.12 — Zenodo DOI: 10.5281/zenodo.19679049
- [ ] Full CMB likelihood (Cobaya + plik/CamSpec) — blocker for PRD/PRL
- [ ] Israel-Stewart viscosity formulation — blocker for perturbation predictions
- [ ] Spanish translations of all 4 papers
- [ ] arXiv submission (Papers 1–3 ready, Paper 4 pending Ωc h² resolution)

---

## Author

Mike Edison Almeida Vallejo — mike.almeida1721@gmail.com
ORCID: [0009-0008-2195-7836](https://orcid.org/0009-0008-2195-7836)

## License

Apache 2.0 — see [LICENSE](LICENSE).
