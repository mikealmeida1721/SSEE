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

## OP-3 — UV-IR Separability Conjecture (Paper 10 / future work) ✅ RESUELTO

**Location:** Paper 10, Postulate C.1 — **revisado 2026-05-16**

**Resolución (cota de supresión EFT):**

**Argumento primario — jerarquía Coleman-Weinberg:**
La jerarquía de escalas (H₀/M)² = (1.45×10⁻³³ eV / 8.81×10⁻³ eV)² ≈ 2.7×10⁻⁶² suprime
el mezclado φ-π en el jacobiano ∂φ/∂χ|_transition por un factor ~10⁶², convirtiendo el
Postulate C.1 en un Theorem C.1 con cota de corrección explícita.

**Unicidad de KALeff:**
Dada M⁴ = 5φ⁸ρ_crit (Paper 10, establecido independientemente) y separabilidad φ/π:
$$K_{\rm ALeff}^2 = \frac{M^4}{6\alpha} = \frac{5\varphi^8\rho_c}{2\varphi^4\rho_c} = \frac{5\varphi^4}{2} \quad\Rightarrow\quad K_{\rm ALeff} = \varphi^2\sqrt{5/2}$$
(única solución monomial en φ — no contiene π).

**Factorización:**
KAL₀ = KALeff × F(φ,π), donde F = KAL₀/KALeff = 1.3338. El factor F transporta π
exclusivamente a través del sector IR (w₀), ausente en el régimen de transición UV.

**Resultado:** H₀^UV = 73.040 km/s/Mpc (< 0.001σ SH0ES). Las correcciones al Theorem C.1
son O(10⁻⁶²) — efectivamente exacto en la práctica cosmológica.

**Límite:** La derivación completa del jacobiano desde P(X,φ) quintaesencial requiere
especificar el modelo de reheating exacto (Paper B futura). La cota EFT garantiza
KALeff = φ²√(5/2) × [1 + O(10⁻⁶²)].

**Script:** `src/ssee_op3_separability.py` (cálculo completo, todos los asserts pasan)

---

## OP-4 — Vainshtein Radius Exceeds Observable Universe (Paper 8) ✅ RESUELTO

**Location:** Paper 8, §4.2 (solar system screening) — **revisado 2026-05-15**

**Resolución (dos argumentos independientes):**

**Argumento primario — estructura EFT Bellini-Sawicki:**
Paper 7 establece αB = αM = αT = 0 para SSEE, con αK = 0.4033 único. Esto implica:
- μ − 1 = (αB + αM)²/αK = 0 (constante de Newton no modificada)
- γ_PPN − 1 = −2αT/(1+αT) = 0 (sin lente gravitacional modificada)

La quinta fuerza en el límite quasi-estático se suprime como (H/k)²:
$$\frac{F_\phi}{F_N}\bigg|_{\rm 1\,AU} \approx 1.6 \times 10^{-31} \ll 10^{-5}\ (\text{Cassini})$$
Satisface la restricción solar por un factor de 6×10²⁵.

**Argumento secundario — k-mouflage (Brax & Valageas 2014):**
K(X) = X/KAL + X²/M⁴ es k-mouflage, NO Galileon. La fórmula Galileon usada en Paper 8 §4.2 original era inaplicable. Radio k-mouflage correcto:
$$r_{\rm km}^3 = \frac{M_{\rm obj}}{4\pi M_{\rm Pl} M^2}, \quad M = 8.81\ \text{meV}$$
$$r_{\rm km}(\odot) = 1.54 \times 10^7\ \text{m} \approx 0.022\,R_\odot \ll r_{\rm Hubble}$$
Todo objeto astrofísico tiene r_km ≪ 1 kpc → quinta fuerza DM activa a escalas cosmológicas.

**Cambios aplicados en Paper 8:**
- §4.1: "Vainshtein-like" → "k-mouflage" (Brax & Valageas 2014)
- §4.2: Reemplazado fórmula Galileon con fórmula k-mouflage + Tabla revisada
- §4.4: "Double GR protection" (Vainshtein) → "EFT suppression" (αB=αM=αT=0)
- Bibitem `brax2014` añadido
- `src/ssee_paper8_figures.py`: figura regenerada con fórmula k-mouflage

**Scripts:** `src/ssee_op4_vainshtein.py` (cálculo completo), `src/ssee_paper8_figures.py` (figura)

---

## OP-5 — S₈ Weak-Lensing Tension (Papers 5–6)

**Location:** Paper 5, Table 3; Paper 6, Table 2

**Clarification — what IS and IS NOT resolved:**
Paper 6 (two-sector φ-DM) resolves the **fσ₈ growth-rate** tension from 2.56σ → 0.50σ
(mean over 6 RSD surveys). That problem is closed. The residual open problem is different:

**Statement:** The **S₈ = σ₈(Ω_m/0.3)^0.5 weak-lensing amplitude** remains at
S₈_eff = 0.761, which is 2.29σ below DES Y3 and ~3σ below KiDS-1000.

This is a different observable: fσ₈ is measured from galaxy peculiar velocities (RSD),
while S₈ is measured from gravitational lensing shear. The two observables probe different
scales and redshift ranges, and respond differently to the WDM transfer function suppression.

**Why S₈ persists:** The WDM suppression in the φ-DM sector primarily reduces power at
k > k_fs = 0.493 h/Mpc. The lensing signal integrates over all scales, so the improvement
in σ₈_eff (0.702 → 0.737) is not sufficient to fully close the gap against DES/KiDS.
Additionally, the KiDS-DES discrepancy (KiDS: S₈ ≈ 0.766, DES: S₈ ≈ 0.776) suggests the
observational situation itself is not fully settled.

**What would resolve it:** Non-linear N-body simulations with baryonic feedback (AGN winds,
supernova) that suppress small-scale power by ~10–20%. These effects are not captured in
the linear CLASS calculation. SSEE would need dedicated hydrodynamic simulations analogous
to IllustrisTNG or EAGLE run with the two-sector SSEE background.

---

## OP-6 — Screening Form Ambiguity in Hubble Tension Resolution (Paper 9) ✅ RESUELTO

**Location:** Paper 9, §3 (f_screen derivation) — **revisado 2026-05-15**

**Resolución — aproximación de universo separado para k-essence:**

La forma multiplicativa sigue de primer principios de la aproximación de universo separado
(Wands et al. 2000; Brax & Valageas 2014). Para k-essence con K(X) = X/KAL + X²/M⁴,
una región sobredensa local corrige H multiplicativamente:

$$\frac{H_{\rm local}^2}{H_{\rm global}^2} = 1 + \frac{8\pi G}{3H^2}\,\delta\rho_{\phi,\rm local}$$

La perturbación de densidad del escalar en el límite cuasi-estático:

$$\frac{\delta\rho_\phi}{\rho_{\rm crit}} = \frac{\alpha_K}{3}\,\frac{\Omega_{m,\rm dyn}}{\mathcal{M}\,(1+w_0)}\,\delta_{\rm local}$$

La **identidad algebraica SSEE** $1 + w_0 = \Omega_{m,\rm dyn}$ (exacta desde el álgebra
φ,π, verificada numéricamente con |diferencia| = 0) cancela el factor Ω_m,dyn:

$$f_{\rm screen} = \frac{\alpha_K}{3\,c_s^2\,\mathcal{M}} = \frac{\alpha_K}{3\,\mathcal{M}} = 0.06725$$

(usando c²_s = 1 del Paper 5, Q1). La forma aditiva correspondería a un sesgo de velocidad
peculiar (Δv/c), no a una corrección de densidad de energía oscura — físicamente distinto.

**Verificación numérica** (`src/ssee_op6_screening_form.py`):
- f_screen (universo separado) = 0.067253
- f_screen (algebraico (π−φ)/Ω²) = 0.067253
- |diferencia| = 4.1×10⁻⁷ < 10⁻⁴ ✓
- H₀,local = 67.96 × 1.06725 = 72.53 km/s/Mpc (0.49σ SH0ES)
- Con H₀^UV = 73.040 (Paper 10, condicional a Postulate C.1): 0.00σ SH0ES

**Cambios aplicados en Paper 9:**
- §3: Derivación desde universo separado k-essence (primer principios)
- §3: Identidad 1+w₀ = Ω_m,dyn explicitada como justificación de la cancelación
- Bibitem `wands2000` y `brax2014` añadidos

**Script:** `src/ssee_op6_screening_form.py` (verificación completa)

---

## Summary Table

| ID | Paper | Problem | Severity | Path to Resolution |
|----|-------|---------|----------|--------------------|
| OP-1 | P4 | Factor 200 in Ω_b h² | High | UV completion or hierarchy derivation |
| OP-2 | P4 | n_s exponent 7 not derived from V(φ) | High | α-attractor inflation calculation |
| OP-3 | P10 | ~~UV-IR separability unproven~~ | ✅ RESUELTO | Jerarquía EFT (H₀/M)²≈10⁻⁶² + KALeff=φ²√(5/2) único; Paper 10 TC.1; script op3 |
| OP-4 | P8 | ~~r_V > r_Hubble para Vainshtein~~ | ✅ RESUELTO | k-mouflage + αB=αM=αT=0 EFT; Paper 8 §4.2/§4.4 revisados |
| OP-5 | P5-6 | S₈ lensing tension 2.29σ DES (fσ₈ resolved in P6) | Medium | N-body + baryonic feedback |
| OP-6 | P9 | ~~Screening form ambiguity~~ | ✅ RESUELTO | Universo separado k-essence + identidad 1+w₀=Ω_m; Paper 9 §3 revisado |

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
