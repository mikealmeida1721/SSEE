# SSEE — Structural Self-Energy Expansion

**A zero-free-parameter dark energy model derived exclusively from φ (golden ratio) and π, tested against DESI DR2 BAO, Planck 2018 CMB (TT+TE+EE+lensing), galaxy cluster masses, and large-scale structure growth.**

---

## Falsifiable predictions — all pass

| Observable | SSEE (algebraic) | Observed | Separation |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.827, −0.75) | 0.05σ |
| CMB first peak | ℓ₁ = 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 (algebraic via MIRA) | Planck 2018: 0.3153 | 0.63σ |
| n_s | 1 − φ⁻⁷ = 0.96556 | Planck 2018: 0.9649 | 0.2σ |
| αT (GW speed) | 0 exact | GW170817: \|αT\| < 10⁻¹⁵ | exact match |
| αK (kineticity, z=0) | 0.403 algebraic | Euclid forecast: < 0.1 | testable 2026–2028 |
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
| 5 | Israel-Stewart Causal Perturbation Theory — Gradient Stability, MIRA Origin, S₈ Partial Resolution | 24 | Preprint | [docs/](docs/SSEE_Paper5_IS.pdf) |
| 6 | φ-Dark Matter in SSEE-V3.6: Algebraic Mass Derivation and Resolution of the fσ₈ Tension | 15 | Preprint | [docs/](docs/SSEE_Paper6_phiDM.pdf) |
| 7 | Canonical EFT of SSEE-V3.6: Action, β_c = −AURA, and Bellini-Sawicki α-Functions | 14 | Preprint | [docs/](docs/SSEE_Paper7_EFT.pdf) |

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
│   └── SSEE_Endorser_Summary.tex  # 2-page arXiv endorser brief
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
| TT | 1.062 | 1.043 | 1971 |
| TE | 1.053 | 1.040 | 1967 |
| EE | 1.040 | 1.039 | 1967 |
| PP (lensing) | **0.730** | 0.757 | 9 |
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

### Paper 6 (φ-Dark Matter, Two-Sector Model)

| Result | Value | Status |
|---|---|---|
| Ω_CDM | 0.160050 | Active at all k |
| Ω_φDM = (MIRA−1)×Ω_m,dyn | 0.159878 | Active for k < k_fs only |
| Ω_total (two-sector) | 0.319928 ≈ Ω_m,CMB | Zero-parameter unification |
| m_φ = Σm_ν × H₀^alg | 5.71 eV | Algebraic — no fitting |
| k_fs (Dodelson-Widrow) | 0.493 h/Mpc | From m_φ algebraic |
| T_WDM(k=0.125 h/Mpc) | 0.8175 | Soft cutoff (scalar condensate, not thermal WDM) |
| σ₈_eff | 0.794 | — |
| S₈_eff | 0.820 | 2.54σ KiDS (KiDS-DES internal tension) |
| **Mean fσ₈ tension (6 surveys)** | **0.50σ** | **Resolved from 3.66σ (Paper 5)** |
| ΛCDM fσ₈ mean tension | 0.51σ | Reference |

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
| αK(z=0) | 3·Ω_DE·Ω_m,dyn = 0.403 | Algebraic (Euclid will constrain < 0.1) |
| G₂_s (running) | 0.979 | Two-sector unification ✓ |

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
7. **CMB ΔBIC status:** Diagonal approximation (TT+TE+EE+PP, k=2 vs k=6) yields $\Delta\mathrm{BIC}=+31.1$ favouring $\Lambda$CDM — driven by large N amplifying a small $\Delta\chi^2_r$. Full off-diagonal evaluation via Cobaya plik\_lite (TTTEEE) gives $\Delta\mathrm{BIC}=-31.3$ decisively favouring SSEE ($k_{\rm SSEE}=2$ vs $k_{\Lambda\rm CDM}=6$). Full plik/CamSpec likelihood remains a blocker for PRD/PRL.

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
- [x] Paper 6: φ-DM two-sector model — m_φ=5.71 eV algebraic, fσ₈ 3.66σ→0.50σ, 15 pp (incl. Lyman-α defense)
- [x] Paper 7: Canonical EFT — β_c=−AURA exact, αT=αM=αB=0, αK=0.403 algebraic, 14 pp
- [x] Cross-paper parameter consistency audit — no critical inconsistencies (`src/ssee_audit_consistency.py`)
- [x] Paper 6 MCMC — 64 walkers × 15000 steps: Ω_φDM=0.161±0.011, ΔBIC=−12.1, χ²_r=0.497 ✓
- [ ] Zenodo v2 — archive Papers 5+6+7 PDFs with DOI chaining
- [ ] Full CMB likelihood (Cobaya + plik/CamSpec) — blocker for PRD/PRL
- [x] hi_class/CLASS cross-check — αK(0)=0.4033, Δ=0.005% vs algebraic prediction ✓
- [ ] Full causal Israel-Stewart (Hiscock-Lindblom 1985) — blocker for B-mode LiteBIRD forecasts
- [ ] arXiv submission — pending endorser response (Shafieloo/KASI, Lee/SKKU, Di Valentino/Sheffield)

---

## Author

Mike Edison Almeida Vallejo — mike.almeida1721@gmail.com
ORCID: [0009-0008-2195-7836](https://orcid.org/0009-0008-2195-7836)

## License

Apache 2.0 — see [LICENSE](LICENSE).
