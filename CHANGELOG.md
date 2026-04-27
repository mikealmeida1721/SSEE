# Changelog

## V3.6 — Session 13: ΔBIC consistency fix +218 → +206 (2026-04-26)

- **ΔBIC audit**: Detected discrepancy between script output (+206) and paper text (+218)
  - Root cause: paper table used lp_MAP=-124.87 from less-converged prior run; current run (N_eff=637,500) gives lp_MAP=-118.67 → ΔBIC=+206.19
  - The ΛCDM MAP (-14.19) was unchanged; full change came from improved SSEE MCMC convergence
  - Three ΔBIC values now clearly distinguished in script and docs:
    - [1] MCMC (k_SSEE=2 vs k_ΛCDM=3): +206 ← updated everywhere
    - [2] Framework (k_SSEE=0 vs k_ΛCDM=6): +192 ← noted in script (not in paper prose, avoids confusion)
    - [3] Dynamic sector (k=1 w₀,wₐ): −5.55 ← unchanged, SSEE favoured
  - Updated: Paper 2 table + all prose, README, AUDIT.md, CLAUDE.md, cover letters ×2
  - Paper 2 recompiled: 20 p → docs/SSEE_Paper2_draft.pdf
- **Script improvement**: `ssee_paper2_mcmc.py` now prints all three ΔBIC values with labels
- **README fix**: ΔBIC +206 → +218 stale value corrected; limitation item #7 added (ΔBIC TT+TE+EE+PP=+13.7)

## V3.6 — Session 12: Task 2B — Press-Schechter quantitative comparison (2026-04-26)

- **Task 2B completed**: Paper 4 §sec:deltac replaces "deferred to Paper 5" with computed values
  - New script `src/ssee_press_schechter.py`: numerically integrates D(z)/D(0), computes σ_M(M,z)
    via Eisenstein-Hu fit, evaluates PS ratio n_SSEE/n_ΛCDM = (δc_SSEE/δc_ΛCDM) × exp[…]
  - D(z=10)/D(0) = 0.1153 (Planck 2018 cosmology, flat ΛCDM)
  - Correction: σ_M≈0.5 at z=10 corresponds to M~10^12 M☉ (not 10^11 M☉ as stated before)
  - Results: M=10^11 M☉ → ×1.06; M=10^12 M☉ → ×1.41; M=3×10^12 M☉ → ×2.0
  - Honest conclusion: SSEE partially alleviates JWST tension but ~10–100× ΛCDM deficit remains
  - New figures: `results/figures/fig_press_schechter.{pdf,png}` (ratio vs σ_M and vs z)
  - New reference in Paper 4.bib: Boylan-Kolchin2023 (Nature Astronomy 7, 728)
  - Paper 4: 13 pages (was 12; +1 for Table 8 + eq:ps_ratio)
- **Internal audit session 11** (commits 7fbe5bc + 3e2918d): fixes C3–C7 + bonus macros (see session 11 entry)

## V3.6 — Session 11: Internal audit C3–C7 + bonus Paper 3 macros (2026-04-26)

- **C3**: Paper 1 r=0.0083 → 0.00813 (abstract + Predictive Register table); AUDIT.md ×3; README.md ×1
- **C4**: Paper 1 cluster table caption R_500 → R_200 (coherent with Paper 2 M1/M2 mass corrections)
- **C5**: Paper 3 χ²_r=0.131 → 0.122 (abstract + conclusions + abstracts_arXiv.txt)
- **C6**: Paper 2 prose "χ²_r=12.95" → 14.91 (KAL₀ only; coherent with sensitivity table)
- **C7**: Paper 2 fig2 caption Δχ²≈+85 → ≈+97 (coherent with table)
- **Bonus**: Paper 3 \SSEE and \LCDM macros undefined → added to preamble (silent pre-existing error)
- Papers recompiled: Paper 1 (11 p), Paper 2 (20 p), Paper 3 (17 p) → docs/
- All 4 papers audited; Paper 4 clean (r=0.00813 correct, M3 applied earlier)

## V3.6 — Session 9: Task 1C — Nine Sovereignties explicit table + Independence Appendix (2026-04-25)

- **Task 1C completed** (Paper 4 Theorem 2): replaced `= \ldots` ellipsis with explicit table of all 9 paths
  - Table~6 (`tab:sovereignties`): all 9 sovereignty names, named-constant formulas, expanded φ,π forms, values
  - Appendix A (`app:independence`): formal statement of independence — hierarchy-level vs formula-level;
    documents SOLAR=MAR and IGNIS=KRYSTOS algebraic equalities; explains why redundant pairs are retained
- **Fixed**: ICEBERG formula in `ssee-data.json` corrected from "MIRA + PYROS" → "MAR + PYROS"
- **Fixed**: Paper 4 params table — Ωb h² tension 2.1σ → 3.2σ; Y_p 0.2473 → 0.2476, tension 0.7σ → 1.7σ
- **Fixed**: Conclusion section — Y_p 0.2473 → 0.2476, updated tensions
- Paper 4 PDF recompiled → 12 pages (was 11) → `docs/SSEE_Paper4_ToE.pdf`

## V3.6 — Session 8: r_d numerical verification + Paper 3 §2.3 fix (2026-04-25)

- **Task 2A completed**: CAMB numerical verification of SSEE sound horizon vs Planck 2018 measured r_d
  - Case A (SSEE+MIRA, Ωb h²=0.02237): r_d = 147.156 Mpc → **0.25σ** from Planck (147.09 ± 0.26 Mpc) ✅
  - Case B (SSEE+MIRA+IS, Ωb h²=0.02285): r_d = 146.749 Mpc → 1.31σ ✅
  - Case C (SSEE dynamic, Ωm=0.160, no MIRA): r_d = 171.650 Mpc → 94.46σ ❌ (expected; MIRA essential)
  - θ* residual tension: 4.46σ in Case A (requires full D_A geometry match — linked to CMB likelihood blocker)
- **Script added**: `src/ssee_verify_rd.py` — standalone CAMB verification, arbiter = Planck observation
- **Paper 3 §2.3 fixed**: removed erroneous "r_d,eff = r_d × |w₀|" claim; MIRA resolves r_d directly
- **Paper 3 §4.1 fixed**: replaced "0.03% of ΛCDM reference" with "0.25σ of Planck 2018 measured r_d"
- **CHANGELOG/AUDIT fixed**: Ωb h² IS 0.02260 → 0.02285; wrong formula corrected to 3(π−φ)/200

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
- **Paper 4**: Y_p BBN = 0.2473 (AlterBBN, Ωb h²=0.02285, 0.7σ); Ωc h² IS = 0.11926 (−0.6σ); δc = 1.6284
- **EFT section**: fixed 3 broken refs (eq:Pi→eq:zeta, DESI2025DR2→AbdulKarim2025, sec:constants→sec:axioms)
- **Paper 4 bib**: added Pisanti2008 (AlterBBN, doi:10.1016/j.cpc.2007.11.013)

## V3.6 — Paper 4: Algebraic ToE (2026-04-22)

- Added Paper 4: algebraic derivation of all CMB background observables from φ and π only
- Nine Sovereignties: 9 independent algebraic paths to 3(φ+π)
- Ωb h² derived: 3(π−φ)/200 = 0.02285 (3.2σ; Paper 4 §sec:baryons)
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
