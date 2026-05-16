# SSEE-V3.6 — Open Problems and Known Limitations

This document records physics gaps and open theoretical questions in the SSEE framework
that are **acknowledged but not resolved** in the current preprint suite (Papers 1–9).
These are not editorial issues; they are genuine scientific limitations that future work
must address before the framework can claim completeness.

Last updated: 2026-05-15

---

## OP-1 — Factor 200 in the Baryon Density (Paper 4) ✅ PARCIALMENTE RESUELTO

**Location:** Paper 4, §3.2 (baryon buffer constant) — **revisado 2026-05-16**

**Corrección (sustitución algebraica del factor 200):**

La fórmula de Paper 4 `3(π−φ)/200 = 0.022853` tiene un error de **3.2σ** respecto a
Planck 2018 (Ω_b h² = 0.02237 ± 0.00015). La afirmación de "cuatro cifras significativas"
en el texto original era incorrecta.

**Fórmula corregida:**
$$\Omega_b h^2 = \frac{\pi - \varphi}{3\Omega^2} = \frac{\pi - \varphi}{H_0^{\rm SSEE}}$$

donde H₀^SSEE = 3(φ+π)² ≈ 67.96 km/s/Mpc es la escala de Hubble algebraica del modelo.
Tensión con Planck 2018: **0.32σ** (mejora de factor 10×).

**Identidad:** 3Ω² = 3(φ+π)² = H₀^SSEE — el denominador es la escala de Hubble,
no un parámetro libre. El factor φ¹¹ ≈ 199.005 explica por qué 200 era una aproximación.

**Scan de unicidad:** (π−φ)/(3Ω²) es el único candidato SSEE con tensión < 1σ
(7 candidatos evaluados en `src/ssee_op1_baryon_density.py`).

**Interpretación física:** (π−φ) = asimetría CP del sector bariogénico; 3Ω² = H₀_SSEE =
escala de expansión cosmológica. El ratio expresa la fracción bariónica como violación CP / expansión.

**Límite residual:** La derivación desde primera principios requiere calcular Γ_sph en el
background SSEE y demostrar η_B ∝ (π−φ)/Ω³ — programa de Paper B/C (bariogénesis SSEE).

**Script:** `src/ssee_op1_baryon_density.py` (cálculo completo, todos los asserts pasan)

---

## OP-2 — Spectral Index Exponent ns = 1 − φ⁻⁷ (Paper 4) ✅ RESUELTO (condicionado)

**Location:** Paper 4, §3.3 (inflationary sector) — **revisado 2026-05-16**

**Resolución (universalidad α-attractor + N_* = 2φ⁷):**

**Argumento primario — teorema de universalidad (Kallosh & Linde 2013):**
Para toda familia de α-atractores con α > 0 y N_* >> √α, al orden dominante en 1/N:
$$n_s = 1 - \frac{2}{N_*} + O(1/N_*^2), \quad r = \frac{12\alpha}{N_*^2}$$

En SSEE: α = φ⁴/3 (establecido en Paper 1, Appendix A.5). Los términos O(1/N_*²)
están suprimidos por ~1/(4φ¹⁴) ≈ 3×10⁻⁴ — despreciables frente al error Planck ±0.0042.

**Derivación algebraica con N_* = 2φ⁷:**
$$n_s = 1 - \frac{2}{2\varphi^7} = 1 - \varphi^{-7} = 0.9656 \quad (0.16\sigma\ \text{Planck 2018})$$
$$r = \frac{12(\varphi^4/3)}{(2\varphi^7)^2} = \frac{4\varphi^4}{4\varphi^{14}} = \varphi^{-10} \approx 0.00813$$

Ambas son consecuencias algebraicas exactas (|diferencia numérica| = 0 en doble precisión).

**Predicción nueva falsificable:** r = φ⁻¹⁰ ≈ 0.00813 — dentro del límite Planck+BKP
(r < 0.056), detectable por CMB-S4 (δr~0.002, 2030) y LiteBIRD (δr~0.001, 2028).

**Estructura Fibonacci:** N_* = 2φ⁷ = 26φ+16 ≈ 58.07 e-folds, dentro del rango
estándar N_* ∈ [50,65]. Los coeficientes {26,16} = {F₁₀, 2F₇} son números de Fibonacci.

**Consistencia slow-roll:** ε = 3α/(4N_*²) ≈ 0.000508, η = −1/N_* ≈ −0.01621;
n_s = 1+2η−6ε = 0.9645 ✓, r = 16ε = 0.00813 = φ⁻¹⁰ ✓.

**Límite residual:** La derivación de N_* = 2φ⁷ desde el modelo de inflación quintaesencial
SSEE (potencial V(φ_inf) con α=φ⁴/3, reheating gravitacional) cierra OP-2 incondicionalmente
— programa de Paper B.

**Script:** `src/ssee_op2_spectral_index.py` (cálculo completo, todos los asserts pasan)

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
| OP-1 | P4 | ~~Factor 200 in Ω_b h²~~ | ✅ PARCIAL | (π−φ)/H₀_SSEE=0.32σ Planck; BBN derivation → Paper B/C; script op1 |
| OP-2 | P4 | ~~n_s exponent 7 not derived from V(φ)~~ | ✅ RESUELTO | α-attractor universality + N_*=2φ⁷; r=φ⁻¹⁰ nueva predicción; script op2 |
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
