# SSEE — Scalar Structural Evolution Equations

**A parameter-free dark energy model derived from algebraic constants (φ, π), tested against DESI DR2 BAO, Planck 2018, and galaxy cluster masses.**

---

## What this repo shows

SSEE V3.6 has two faces:

- **Dynamic sector** (w₀, wₐ, Ωₘ_eff): competitive with CPL and ΛCDM on DESI BAO + cluster data — validated in Paper 2.
- **Background sector** (CMB H₀ scale): tension with Planck 2018 is quantified explicitly — open problem for Paper 3.

The model fixes *all* cosmological parameters algebraically from φ and π. There are no free parameters fitted from data.

---

## What remains open

- CMB background tension via H₀ scale (SSEE fixes Ωₘ_eff = 0.160, shifting r_d relative to ΛCDM).
- Full CMB power spectrum comparison (Paper 3).
- Folding Symmetry interpretation in a quantum gravity context.

---

## Repository structure

```
SSEE/
├── src/              # Analysis scripts (pipeline reproducible)
├── notebooks/        # Exploratory notebooks
├── data/
│   ├── raw/          # Observational data (DESI, Planck, clusters)
│   └── processed/    # Derived inputs
├── results/
│   ├── figures/      # All generated figures
│   ├── tables/       # Statistical results
│   └── logs/         # MCMC logs
├── manuscript/       # Paper 2 LaTeX draft
└── docs/             # Reference documents
```

---

## Installation

```bash
git clone https://github.com/mikealmeida1721/SSEE.git
cd SSEE
pip install -r requirements.txt
```

---

## How to reproduce Paper 2

Run the three scripts in order:

```bash
python src/ssee_paper2_analysis.py   # w0-wa plane, sigma deviations, sensitivity
python src/ssee_paper2_mcmc.py       # Bayesian MCMC: SSEE vs ΛCDM vs CPL
python src/ssee_paper2_figures.py    # Generate all figures → results/figures/
```

---

## Key results (Paper 2)

| Observable | SSEE (algebraic) | Observed | Tension |
|---|---|---|---|
| w₀ | −0.840 | −0.838 ± 0.096 | < 0.03σ |
| wₐ | −0.670 | −0.62 ± 0.39 | < 0.13σ |
| Ω_DE | 0.840 | 0.685 ± 0.019 | ~8σ (background) |

---

## What's hypothesis vs. result

- **Algebraic (hypothesis):** All SSEE constants (w₀, wₐ, Ωₘ_eff, etc.) are derived from φ and π — no fitting.
- **Result:** Bayesian MCMC comparison shows SSEE dynamic sector is statistically competitive with CPL on DESI + cluster data.
- **Open:** Background CMB tension is a quantified discrepancy, not a resolved result.

---

## Roadmap

- [x] Paper 2: Dynamic sector validation (DESI BAO + clusters)
- [ ] Paper 3: CMB background tension resolution

---

## License

Apache 2.0 — see [LICENSE](LICENSE).

## Contact

Mike Almeida — mike.almeida1721@gmail.com
