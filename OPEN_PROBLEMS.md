# SSEE-V3.6 — Open Problems and Known Limitations

This document records physics gaps and open theoretical questions in the SSEE framework
that are **acknowledged but not resolved** in the current preprint suite (Papers 1–9).
These are not editorial issues; they are genuine scientific limitations that future work
must address before the framework can claim completeness.

Last updated: 2026-05-15

---

## OP-1 — Factor 200 in the Baryon Density (Paper 4)

**Location:** Paper 4, §3.2 (baryon buffer constant)

**Statement:**
$$\Omega_b h^2 = \frac{3(\pi - \varphi)}{200} \approx 0.02237$$

The factor 200 in the denominator has no derivation from first principles. The value matches
Planck 2018 to four significant figures, but it is inserted by dimensional matching, not
derived from the algebraic hierarchy φ, π, Ω, β, KAL, etc.

**Why it matters:** Without a first-principles derivation, this is a numerical coincidence,
not a prediction. A skeptical referee will flag this immediately.

**What would resolve it:** A derivation showing that 200 emerges from the SSEE constant
hierarchy (e.g., as a ratio of dimensionless invariants) or from a UV completion
(string landscape, large-field inflation). Until then, this should be presented explicitly
as a phenomenological fit in the paper text, not as an algebraic derivation.

---

## OP-2 — Spectral Index Exponent ns = 1 − φ⁻⁷ (Paper 4)

**Location:** Paper 4, §3.3 (inflationary sector)

**Statement:**
$$n_s = 1 - \varphi^{-7} \approx 0.9649$$

The exponent 7 is motivated in the paper by the count of independent SSEE constants in the
hierarchy (levels 1–3), not by a derivation from an inflationary Lagrangian V(φ).

**Why it matters:** The prediction matches Planck 2018 to 0.1%, but exponent counting is
not a physical argument. The standard derivation of ns requires a slow-roll calculation
from a specific potential V(φ_inflaton). The SSEE framework does not provide this.

**What would resolve it:** An α-attractor inflation model (e.g., α = φ⁴/3 as in Paper 1
Appendix A.5) that yields ε, η, and thus ns analytically. This is the most tractable path;
the Starobinsky limit already gives ns ≈ 1 − 2/N with N~60 e-folds, which is close but
not identical to 1−φ⁻⁷.

---

## OP-3 — UV-IR Separability Conjecture (Paper 10 / future work)

**Location:** Conjectured in multi-paper discussions; not yet written as a standalone section

**Statement:** Conjecture C.1 — the φ-sector (IR, large-scale structure) and the π-sector
(UV, quantum corrections) decouple algebraically, so that SSEE constants factorize as
$f(\varphi, \pi) = g(\varphi) \cdot h(\pi)$ in the relevant observables.

**Why it matters:** If true, this would justify treating SSEE as a consistent EFT with a
well-defined decoupling limit. If false, φ and π mix at loop level and the hierarchy is
not radiatively stable.

**Current status:** Not proven. The algebraic coincidences (e.g., Ω = φ + π mixes both)
suggest the factorization does NOT hold in general. Routes A/B/C in the M⁴ derivation
give consistent results only at tree level.

**What would resolve it:** A one-loop Coleman-Weinberg calculation in the SSEE scalar
potential showing that φ-π cross-terms are suppressed by a large mass ratio.

---

## OP-4 — Vainshtein Radius Exceeds Observable Universe (Paper 8)

**Location:** Paper 8, §4.2 (solar system Vainshtein screening)

**Statement:** The Vainshtein radius for the Sun computed from the SSEE EFT parameters gives
$$r_V \sim 1.84 \times 10^{44}\ \text{m} \gg r_{\rm Hubble} \sim 1.3 \times 10^{26}\ \text{m}$$

**Why it matters:** If r_V is larger than the Hubble radius, the Vainshtein mechanism
screens modifications everywhere, including at cosmological scales where SSEE predicts
deviations from ΛCDM. This would make the EFT of Paper 7 (αK = 0.4033) internally
inconsistent with the screening model of Paper 8.

**Current status:** Flagged as a tension in Paper 8 footnote. The resolution depends on
whether the effective EFT cutoff Λ_SSEE is identified with M_Pl or with a lower scale
set by the SSEE hierarchy.

**What would resolve it:** A careful identification of the EFT cutoff and the effective
Planck mass M_eff in the SSEE Lagrangian, followed by recomputation of r_V with the
correct M_eff. Alternatively, a transition from Vainshtein to k-mouflage screening
(which has a different radius formula) may be more appropriate given αK ≠ 0.

---

## OP-5 — fσ₈ Residual Tension with BOSS (Papers 5–6)

**Location:** Paper 5, Table 3; Paper 6, Table 2

**Statement:** In the single-sector (CDM only) SSEE background with Ω_m,dyn = 0.160, the
predicted fσ₈(z = 0.5) = 0.341, which is 2.67σ below the BOSS measurement.

**Current mitigation:** Paper 6 introduces the φ-DM two-sector model, which reduces the
mean tension across six RSD surveys from 2.56σ to 0.50σ (σ₈_eff = 0.737, S₈ = 0.761).
The BOSS point at z = 0.5 specifically is reduced but not eliminated.

**Why the tension persists:** Ω_m,dyn = 0.160 is a fixed algebraic prediction. Any model
with this effective matter density will produce suppressed growth at intermediate redshift
relative to ΛCDM. The φ-DM suppression partially compensates but acts primarily on small
scales (k > k_fs = 0.493 h/Mpc), not at the BAO/RSD scale.

**What would resolve it:** N-body simulations with baryonic feedback (AGN, supernova) that
modify the effective σ₈ on the scales probed by BOSS. If baryonic suppression brings the
ΛCDM prediction down by ~10%, the residual tension with SSEE narrows. This requires
supercomputer time and is beyond the current analytical framework.

---

## OP-6 — Screening Form Ambiguity in Hubble Tension Resolution (Paper 9)

**Location:** Paper 9, §3 (f_screen derivation)

**Statement:** The f_screen = αK/(3·MIRA) = (π−φ)/Ω² = 0.06725 is derived from a
multiplicative screening ansatz:
$$H_0^{\rm local} = H_0^{\rm SSEE} \times (1 + f_{\rm screen})$$

An additive form H₀^local = H₀^SSEE + ΔH₀ with ΔH₀ derived from the same constant
gives a different central value.

**Why it matters:** The choice between multiplicative and additive screening is physically
motivated by the Vainshtein/k-mouflage geometry, but Paper 9 does not derive which form
is correct from the SSEE Lagrangian of Paper 7. The agreement H₀,local = 72.86 km/s/Mpc
(0.17σ from SH0ES) is dependent on the multiplicative form.

**What would resolve it:** A full perturbative calculation of the fifth-force contribution
to the local distance ladder in the SSEE EFT, analogous to the Vainshtein calculation of
Paper 8. This would fix the functional form of the screening correction unambiguously.

---

## Summary Table

| ID | Paper | Problem | Severity | Path to Resolution |
|----|-------|---------|----------|--------------------|
| OP-1 | P4 | Factor 200 in Ω_b h² | High | UV completion or hierarchy derivation |
| OP-2 | P4 | n_s exponent 7 not derived from V(φ) | High | α-attractor inflation calculation |
| OP-3 | P10 | UV-IR separability unproven | Medium | One-loop CW effective potential |
| OP-4 | P8 | r_V > r_Hubble for Vainshtein | High | Identify Λ_SSEE, recompute r_V |
| OP-5 | P5-6 | fσ₈ residual 2.67σ (single-sector) | Medium | N-body + baryonic feedback |
| OP-6 | P9 | Screening form ambiguity | Medium | Full EFT perturbative calculation |

**Severity legend:** High = referee would likely request resolution before acceptance;
Medium = requires acknowledgment and discussion; Low = cosmetic or presentational.

---

## Note on Methodology

These problems are documented here rather than concealed because scientific integrity
requires pre-registration of known limitations. Referees and collaborators should be
directed to this document when evaluating the strength of the SSEE predictions.

None of these problems falsify SSEE at the current observational precision — they define
the boundary of what has been rigorously established versus what remains as working
hypotheses. The resolution of OP-1 through OP-6 constitutes the research agenda for
SSEE-V4.0.
