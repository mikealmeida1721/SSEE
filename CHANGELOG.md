# Changelog

## V3.6 — Session 7: Repo audit & documentation (2026-04-25)

- Repo structure: `docs/` reduced to 4 canonical PDFs (one per paper)
- Consolidated `docs/archive/` → `archive/` at root (eliminated nested archive)
- Fixed `manuscript/SSEE_Paper2_draft.pdf` removed from git tracking
- Cleaned root: removed `texput.log`; moved `mcmc_run.log` to `results/logs/`
- Updated AUDIT.md: Paper 4 section, IS Ωc h², Predictive Register table, ΔBIC two-question distinction
- Updated CHANGELOG, CITATION.cff, archive/README.md, data/README.md

## V3.6 — Session 6: Peer-review fixes 1B, 1F, 1G, 2C (2026-04-25)

- **Paper 2**: added "Two ΔBIC values — two physical questions" paragraph (ΔBIC=+218 background vs ΔBIC=−5.55 dynamic sector)
- **Paper 3**: added asymmetry paragraph (k=1 SSEE predictions vs k=6 ΛCDM fits; Zenodo timestamp)
- **Paper 1**: added Predictive Register table (w₀, wₐ, MIRA pre-DESI; Ωm, n_s, H₀ retrodictions; r future)
- **Paper 4**: Y_p BBN = 0.2473 (AlterBBN, Ωb h²=0.02260, 0.7σ); Ωc h² IS = 0.11926 (−0.6σ); δc = 1.6284
- **EFT section**: fixed 3 broken refs (eq:Pi→eq:zeta, DESI2025DR2→AbdulKarim2025, sec:constants→sec:axioms)
- **Paper 4 bib**: added Pisanti2008 (AlterBBN, doi:10.1016/j.cpc.2007.11.013)

## V3.6 — Paper 4: Algebraic ToE (2026-04-22)

- Added Paper 4: algebraic derivation of all CMB background observables from φ and π only
- Nine Sovereignties: 9 independent algebraic paths to 3(φ+π)
- Ωb h² derived: (π−φ)/[6(φ+π)⁴] = 0.02260 (1.5σ)
- Ωc h² IS: KAL₀ × Ωb h² × n_s = 0.11926 (−0.6σ)
- δc SSEE: δc,EdS × n_s = 1.6284 (linked to JWST early galaxy excess)
- Y_p BBN: 0.2473 (AlterBBN, 0.7σ)
- Title revised to "Algebraic Derivation of CMB Background from φ and π"

## V3.6 — Paper 3: CMB Confrontation (2026-04-20)

- **MIRA hypothesis confirmed**: Ωm,CMB = Ωm,dyn × MIRA = 0.160 × 1.9989 = 0.3199 (within 0.63σ of Planck 0.3153)
- MIRA = (φ+β)/2, pre-existing in Genesis 5.12 as "Observational Frequency"
- CMB TT peaks with SSEE+MIRA: ℓ=221, 538, 815 (Planck: 220, ~540, ~810)
- χ²_r = 1.062 vs ΛCDM 1.043 over 1971 multipoles; ΔBIC(TT, k=1 vs k=6) = −40.3 (plik_lite)
- TT+TE+EE+PP lensing all computed; SSEE favoured in lensing (χ²_r=0.730 vs 0.757)
- Added `src/ssee_paper3_cmb.py` — full CAMB implementation with IS growth suppression (γ=φ⁻¹)
- Added `manuscript/SSEE_Paper3_draft.tex` — complete paper (Sections 1–7, 17 pages)

## V3.6 — Paper 2 (2026-04-20)

- Bayesian MCMC validation against DESI DR2 BAO + Planck 2018 + cluster masses (IGIMF-corrected)
- N_eff ≈ 637,500 (SSEE), 63,789 (ΛCDM), 43,793 (CPL); 100 walkers
- Full 13×13 DESI DR2 covariance matrix implemented
- Dynamic sector (w₀, wₐ): 0.05σ from DESI best-fit; ΔBIC = −5.55 (SSEE favoured)
- Full-model ΔBIC = +218 (background mismatch, not falsification; documented explicitly)
- H₀ = 66.75⁺⁰·⁴⁴₋₀.₄⁴ km/s/Mpc; χ²_r clusters = 0.122 (4 clusters)
- Added `manuscript/SSEE_Paper2_draft.tex` — 18 pages, submission-ready

## V3.6 — Paper 1 (2026-04-15)

- First derivation of SSEE algebraic constants from φ and π
- EFT section: Israel-Stewart viscosity framework, Ωc h² IS derivation
- Lagrangian action (phenomenological), σ8/S8 prediction via IS growth index γ=φ⁻¹=0.618
- Qualitative comparison with ΛCDM background

## V3.5 — Paper 2 draft

- Initial analytic comparison on w₀-wₐ plane
- Ω_DE chi-squared analysis vs DESI

## V3.0 — Paper 1 initial

- First derivation of SSEE algebraic constants from φ and π
- Qualitative comparison with ΛCDM background
