# SSEE — Structural Self-Energy Expansion

**A zero-free-parameter dark energy model derived from φ (golden ratio) and π, tested against DESI DR2 BAO, Planck 2018 CMB (TT+TE+EE+lensing), and galaxy cluster masses.**

---

## Three falsifiable predictions — all pass

| Observable | SSEE (algebraic) | Observed | Separation |
|---|---|---|---|
| (w₀, wₐ) | (−0.840, −0.670) | DESI DR2: (−0.827, −0.75) | 0.05σ |
| CMB peak ℓ₁ | 221 | Planck PR4: ~220 | Δℓ = 1 |
| Ωm,CMB | 0.3199 (via MIRA) | Planck 2018: 0.3153 | 0.63σ |

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
│   ├── SSEE_Paper2_draft.tex      # Paper 2: MCMC validation
│   ├── SSEE_Paper3_draft.tex      # Paper 3: CMB confrontation
│   ├── ssee_paper2.bib
│   └── ssee_paper3.bib
├── docs/                          # Compiled PDFs
├── AUDIT.md                       # Full reproducibility guide + known limitations
└── archive/                       # Obsolete drafts
```

---

## How to reproduce

```bash
pip install camb emcee scipy numpy matplotlib
```

**Paper 2 — MCMC validation** (~10 min, N_s=6000):
```bash
python src/ssee_paper2_mcmc.py
```

**Paper 3 — CMB power spectrum**:
```bash
python src/ssee_paper3_cmb.py
```

See [AUDIT.md](AUDIT.md) for expected outputs and known limitations.

---

## Key results

### Paper 2 (DESI DR2 + Planck 2018 + clusters)

| Metric | Value |
|---|---|
| χ²_r clusters (7 clusters, IGIMF-corrected) | 0.131 |
| χ²_2D (w₀-wₐ vs DESI DR2) | 0.080 → 0.05σ |
| H₀ SSEE | 66.66⁺⁰·⁴⁷₋₀.₄₆ km/s/Mpc |
| ΔBIC (dynamic sector isolated) | −13.5 (SSEE favoured) |

### Paper 3 (Planck PR4 CMB)

| Spectrum | SSEE χ²_r | ΛCDM χ²_r | N |
|---|---|---|---|
| TT | 1.062 | 1.043 | 1971 |
| TE | 1.053 | 1.040 | 1967 |
| EE | 1.040 | 1.039 | 1967 |
| PP (lensing) | **0.730** | 0.757 | 9 |
| ΔBIC (TT, k=0 vs k=6) | **−6.9** (SSEE favoured) | — | — |

---

## Known limitations

Documented honestly in the papers. Full disclosure in [AUDIT.md](AUDIT.md):

1. **H(z) tension**: SSEE χ²_r=1.861 vs ΛCDM 0.458 — genuine tension, stated explicitly.
2. **ΔBIC = +206** (Paper 2): applies ΛCDM's Friedmann framework to an SSEE background — acknowledged as a framework-internal penalty.
3. **Two-sector Ωm**: Ωm,dyn=0.160 ≠ Ωm,CMB=0.3199. Physical justification via MIRA in Paper 3 §2.3.
4. **No Lagrangian action**: SSEE is a phenomenological framework. Formal EFT connection is Level 3 work.

---

## Roadmap

- [x] Paper 1: Algebraic framework (φ, π → w₀, wₐ)
- [x] Paper 2: MCMC validation against DESI DR2 + Planck 2018 + clusters
- [x] Paper 3: CMB confrontation TT+TE+EE+lensing (Planck PR4)
- [ ] Zenodo deposit of Genesis 5.12 (DOI in preparation)
- [ ] arXiv submission

---

## Author

Mike Edison Almeida Vallejo — mike.almeida1721@gmail.com  
ORCID: 0009-0008-2195-7836

## License

Apache 2.0 — see [LICENSE](LICENSE).
