# Texto para Zenodo V7 — copiar/pegar

## Título (sin cambios)
SSEE-V3.6: Structural Self-Energy Expansion — A Minimal-Parameter Dark Energy Framework (Papers 1–10 + Consolidated Documents)

## Descripción (campo "Description")

Version 7 of the SSEE-V3.6 preprint suite (10 papers + consolidated journal
documents). This is a **major revision**: the framework is no longer presented
as a "zero free parameter" model. It is now framed honestly as a
**minimal-parameter framework** (~3 effective parameters vs 6 for ΛCDM), with
all background-sector quantities derived algebraically from φ (golden ratio)
and π, and phenomenological dark-matter-sector extensions explicitly tracked
as open problems (OPEN_PROBLEMS.md, included).

**Pre-registered falsifiable prediction (key item of this version):**
the canonical φ-DM particle, m_φ = 36.95 eV
(= Σm_ν^active × (Ω⁴ + AURA·KAL) = 0.0690 eV × 535.28, forward prediction,
zero fitting), implies a free-streaming imprint at
**k_fs = 0.659 h/Mpc** in the matter power spectrum, testable by
**DESI Y3 and Euclid (2026–2028)**. This deposit timestamps the prediction
before those data releases.

Main results (canonical values, see VERIFICATION_LEDGER.md):
- Background: w₀ = −0.840, wₐ = −0.670 derived algebraically; 0.05σ vs DESI DR2 BAO.
- CMB: Planck PR4 TT+TE+EE+lensing, χ²_r ≈ 1.04 (parity with ΛCDM); ΔBIC = −32.2 in favor of SSEE (k=2 vs k=6, plik_lite/Cobaya).
- MCMC: H₀ = 66.53 ± 0.44 km/s/Mpc (DESI DR2 + Planck prior + clusters).
- S₈ tension: two-sector φ-DM gives S₈_eff = 0.766 (0.01σ from KiDS) — resolved within the analyzed regime.
- Hubble tension: screening fraction f_screen = (π−φ)/Ω² = 0.06725 → H₀,local = 71.87 km/s/Mpc (1.12σ from SH0ES; UV completion: 72.05 km/s/Mpc, 0.96σ, conditional).

## Changelog V7 (campo "Notes" o al final de la descripción)

**V7 (2026-06-12) — major revision:**
1. Honest reframing: "0 free parameters" → minimal-parameter framework (~3 effective); explicit postulate accounting (D, S fundamental + M, I auxiliary) in Paper 1.
2. Canonical φ-DM particle: m_φ = 36.95 eV (forward prediction, zero fitting; free scalar Lagrangian); supersedes the retired 5.60 eV ansatz. k_fs = 0.659 h/Mpc pre-registered for DESI Y3/Euclid.
3. Self-consistent Hubble cascade: Σm_ν = 0.0690 eV → H_MIRA = 67.037 km/s/Mpc → local 71.87 (IR) / 72.05 (UV) km/s/Mpc.
4. Paper 8 revised: lensing prediction now unconditional (mass fixed by Paper 6); two-limit analysis retained.
5. Paper 9: new Appendix P + quantitative comparison vs EDE/SIDR/late-DE.
6. Full hostile-referee audit: figure-level sweep (pdftotext) — all retired numbers purged from text AND figures; verification guardian 100/100 green.
7. OPEN_PROBLEMS.md and VERIFICATION_LEDGER.md included: all open physical gaps (OP-1…OP-14) and canonical-value provenance documented.

## Archivos a subir
- SSEE_V7_papers.zip (13 PDFs + OPEN_PROBLEMS.md + VERIFICATION_LEDGER.md)
- SSEE_V7_code.zip (10 arXiv source tarballs + build script + requirements)
  (o subir los PDFs sueltos si prefieres que se vean individualmente)

## Metadatos
- Resource type: Preprint
- Creator: Almeida Vallejo, Mike Edison
- License: CC-BY-4.0 (igual que versiones previas)
- Keywords: dark energy, cosmology, DESI, Hubble tension, S8 tension, golden ratio, minimal-parameter, free-streaming dark matter
- Related identifiers: "New version of" → DOI previo (10.5281/zenodo.20093447)
