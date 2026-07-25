# submission_packages/ — ⚠️ STALE, regenerate before submitting

**These `.tar.gz` arXiv bundles are dated 2026-06-12 and are SUPERSEDED.**

They predate three sets of major corrections and MUST NOT be submitted as-is:

- **ω_m-direct reframe** (2026-06-18/19) — OP-8 dissolved; the bundles still carry the
  old matter-factor narrative.
- **Canonical φ-DM particle lineage** — the Paper 6 bundle still contains
  **m_φ = 36.95 eV** (the retired AURA·KAL era, two reframes before the canonical
  **40.70 eV** SOLAR²·KRYSTOS with the ν-closure constant C=93.14).
- **V-L4-DESI geometry fix** (2026-07-09) — H₀ posterior 66.41 → 67.95; luego **67.79 ± 0.35** (R25, ω_m fijo),
  r_d, θ*, χ²_BAO all corrected.
- **Full 10-paper + consolidation audit** (2026-07-10) to 10/10.

For canonical values see `../CANONICAL_VALUES.yaml`, `../README.md`,
`../VERIFICATION_LEDGER.md`.

## How to regenerate (do this at submission time, after the papers are final)

```bash
python3 archive/codigo/build_arxiv_packages.py    # generator (moved to archive/ in the zenodo cleanup)
```

The generator bundles each paper's current `.tex` + `.bbl` + `figures/` from
`manuscript/` and `results/figures/`. Rebuild only once the manuscript suite is
frozen for submission, then verify each bundle against `CANONICAL_VALUES.yaml`
(no retired numbers) before uploading.
