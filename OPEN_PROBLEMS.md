# SSEE-V3.6 — Open Problems and Known Limitations

This document records physics gaps and open theoretical questions in the SSEE framework
that are **acknowledged but not resolved** in the current preprint suite (Papers 1–9).
These are not editorial issues; they are genuine scientific limitations that future work
must address before the framework can claim completeness.

Last updated: 2026-05-16

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

**Argumento de bariogénesis Sakharov (refuerzo formal — 2026-05-16):**

El script `src/ssee_op1_baryogenesis.py` establece la estructura Sakharov que sustenta
la fórmula (π−φ)/H₀_SSEE:

**Condición 1 — Violación de número bariónico:** Esfalerón electroweak con tasa
Γ_sph ~ α_W⁵ T⁴ (Arnold & McLerran 1987); activo para T > T_sph ≈ 131.7 GeV.

**Condición 2 — Violación CP:** δ_CP = (π−φ)/Ω = 0.3201. Este es el parámetro SSEE
que mide la asimetría entre el sector trascendental π (gauge boson loops) y el sector
algebraico φ (campo escalar). A temperatura T_EW: H_EW = (π/√90)×√g*×T_EW²/M_Pl ≈ 2.8×10⁻¹⁵ GeV.

**Condición 3 — No-equilibrio térmico:** Inflación quintaesencial → reheating gravitacional
con T_rh ≪ T_EW. El factor de dilución requerido f_dil ≈ 1.095×10⁻¹⁸ implica T_rh ~ 10⁻⁴ GeV,
consistente con el rango T_rh ~ 10⁻² – 10⁴ GeV para producción gravitacional de partículas.

**Estructura algebraica de η_B:**
$$\eta_B \sim \frac{\delta_{\rm CP} \times \Gamma_{\rm sph}}{T_{\rm EW}^4 H_{\rm EW}} \times f_{\rm dil}^{-1} \propto \frac{\pi - \varphi}{\Omega} \cdot \frac{\alpha_W^5}{H_{\rm EW}/T_{\rm EW}^4} \cdot f_{\rm dil}^{-1}$$

Cuando se normaliza con el denominador cosmológico H₀_SSEE = 3Ω², la dependencia estructural
η_B ∝ (π−φ)/Ω³ ∝ (π−φ)/(3Ω²) × (1/Ω) reproduce la fórmula empírica Ω_b h² = (π−φ)/H₀_SSEE
a nivel de conteo de potencias en Ω.

**ESTATUS HONESTO (auditoría 2026-05-16):** lo anterior NO es una derivación de η_B.
El script `ssee_op1_baryogenesis.py` calcula un η_B^naive ~ 5×10⁸ (no físico: η_B ≤ 1)
y **retro-calcula** el factor de dilución f_dil ~ 10⁻¹⁸ exigiendo que el producto
iguale el η_B observado. Es un *consistency check* (verifica que el T_rh implicado,
~10⁻⁴ GeV, cae en el rango de reheating gravitacional), NO una predicción. El scan
del paso [7] del script no discrimina entre candidatos δ_CP — todos "funcionan" con
algún f_dil. La fórmula Ω_b h²=(π−φ)/(3Ω²) es un ansatz algebraico (Type-P, Postulado D
de Paper 1); el mecanismo Sakharov motiva su FORMA, no deriva su valor. Paper 4 §3.2
revisado en consecuencia (commit de la sesión).

**Límite residual de OP-1 — la derivación genuina (programa Paper B/C):**
1. Calcular T_rh exacto desde V(φ_inf) con α = φ⁴/3 (quintessential inflation)
2. Integrar g*(T) desde T_rh hasta T_EW para obtener el factor de dilución exacto
3. Evaluar Γ_sph(T_EW)/H(T_EW) en el background SSEE (no ΛCDM)
4. Resolver la ecuación de Boltzmann para η_B sin f_dil retro-calculado y demostrar
   que el producto reproduce η_B = 6.12×10⁻¹⁰ (BBN observacional) ab initio

**Scripts:** `src/ssee_op1_baryon_density.py` + `src/ssee_op1_baryogenesis.py`

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
SSEE (potencial V(φ_inf) con α=φ⁴/3) cierra OP-2 incondicionalmente — programa de Paper B.

**Resultado numérico Paper B (ssee_paperB_Nstar.py):**

Script `src/ssee_paperB_Nstar.py` verifica la Conjetura B.1 numéricamente:
- V_end^(1/4) = 2.27×10¹⁶ GeV  (φ_end = 1.382 Mpl, ε(φ_end) = 1.000 ✓)
- T_rh que produce N_* = 2φ⁷ exacto: **9.345×10¹⁵ GeV**
- Este T_rh corresponde a ρ_rh ≈ V_end → T_rh ≈ (30/π²g*)^(1/4) × V_end^(1/4) ≈ 0.41 × V_end^(1/4)
- Elegante: si ρ_end/ρ_rh = 3 exacto (físicamente ρ_end = 3V_end con ε=1 y ρ_rh = V_end):
  N_* = 58.25 − ln(3)/6 = 58.067 ≈ 2φ⁷ = 58.069 → Δ = 0.002 e-folds = O(1/N_*)

**Nota:** T_rh(inflación) ≈ 9.4×10¹⁵ GeV es la temperatura de reheating para N_*;
distinta de T_bary ~ 10⁻⁴ GeV del argumento Sakharov (OP-1), que corresponde al
epoch de bariogénesis. Estas son dos temperaturas físicamente distintas.

**Script:** `src/ssee_op2_spectral_index.py` (n_s, r) + `src/ssee_paperB_Nstar.py` (T_rh completo)

**Resultado numérico Paper B (ssee_paperB_DW.py) — RESULTADO NEGATIVO:**

Script `src/ssee_paperB_DW.py` evalúa el segundo problema de Paper B: el mecanismo
de producción de φ-DM (m_φ=5.60 eV) que reproduce Ω_φDM h²=0.0739.

- Mecanismo Dodelson-Widrow (mezcla activo-estéril): el ángulo requerido es
  sin²(2θ)_DW = 4.7846×10⁻⁵ (fórmula Boyarsky-Ruchayskiy-Shaposhnikov 2009).
- Scan de ~80 combinaciones algebraicas de φ,π,Ω,MIRA,AURA: el mejor candidato es
  φ⁻²¹ = 4.086×10⁻⁵, con **Δ = 14.6%** — supera el umbral del 10% para "limpio".
- **Conclusión:** el ángulo de mezcla DW NO admite forma algebraica SSEE. Si se
  adoptara DW, sin²(2θ) sería un parámetro libre — violando el principio de cero
  parámetros. Es el primer resultado del modelo que no entra en el rango esperado.
- Nota de rigor: la fórmula DW (constante C_DW) tiene ~15–20% de incertidumbre
  teórica propia, por lo que perseguir un match sub-1% contra ese blanco carece
  de sentido — `Ω×10⁻⁵` (mantisa 0.5%) es numerología: el factor 10⁻⁵ no es
  algebraico (φ⁻ⁿ no genera 10⁻⁵ con n entero; requeriría n=23.92).
- **Ruta correcta (abierta):** producción gravitacional (Parker; Chung-Kolb-Riotto)
  — sin ángulo de mezcla, depende solo de m_φ y V(φ), ya fijos en SSEE. La
  estimación de orden de magnitud del script sale corta por un factor grande
  (solo prefactor); requiere la integral de producción completa con α=φ⁴/3.

**Script:** `src/ssee_paperB_DW.py` (DW scan + producción gravitacional preliminar)

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

## OP-5 — S₈ Weak-Lensing Tension (Papers 5–6) ✅ PARCIALMENTE RESUELTO (Nivel 1)

**Location:** Paper 5, Table 3; Paper 6, Table 2 — **revisado 2026-05-16**

**Clarification — what IS and IS NOT resolved:**
Paper 6 (two-sector φ-DM) resolves the **fσ₈ growth-rate** tension from 2.56σ → 0.50σ
(mean over 6 RSD surveys). That problem is closed. The residual open problem is different:

**Statement:** The **S₈ = σ₈(Ω_m/0.3)^0.5 weak-lensing amplitude** tiene como baseline
S₈_eff = 0.761 (Paper 6), que era 2.29σ por debajo de DES Y3 y ~3σ de KiDS-1000.

**Nivel 1 — HMcode-2020 baryonic feedback (CLASS, laptop) — COMPLETADO 2026-05-16:**

Script `src/ssee_op5_hmcode.py` implementa retroalimentación bariónica AGN via
HMcode-2020_baryonic_feedback en CLASS (Mead et al. 2020, log10T_heat=7.8):

Resultados CLASS HMcode-2020 con parámetros SSEE (H₀=66.75, Ω_m=0.3199, w₀=−0.840, wₐ=−0.670):

| k [h/Mpc] | B(k) = P_bar/P_hm |
|---|---|
| 0.1 | 0.9974 |
| 0.3 | 0.9878 |
| 0.5 | 0.9765 |
| 1.0 | 0.9560 |
| 2.0 | 0.9203 |

B_eff (lensing k=0.03–2 h/Mpc, peso k) = **0.9447** (supresión 5.53% en P(k))
B_sigma8 (top-hat integral, k<2 h/Mpc) = **0.9956** (supresión 0.44% en σ₈_eff)

Tensiones S₈ (aplicando supresión al baseline Paper 6 S₈=0.761):

| Escenario | S₈ | DES Y3 | KiDS-1000 |
|---|---|---|---|
| Paper 6 baseline (φ-DM + WDM) | 0.761 | 0.09σ | −0.25σ |
| + HMcode-2020 baryonic (Mead+20) | 0.758 | −0.06σ | −0.42σ |
| DES Y3 (observado) | 0.759 | 0.00σ | — |
| KiDS-1000 (observado) | 0.766 | — | 0.00σ |

**El baseline Paper 6 ya está dentro de 1σ DES** (0.09σ). HMcode añade Δσ = 0.03σ de mejora.

**Por qué el baseline está tan bien:** El two-sector φ-DM (m_φ=5.60 eV, k_fs=0.493 h/Mpc)
ya suprime P(k) en k > k_fs, combinado con MIRA (Ω_m=0.320). El HMcode añade supresión bariónica
suave adicional, principalmente a k > 0.5 h/Mpc.

**Nivel 2 — N-body completo (proyección):**

HMcode-2020 captura ~60–70% de la supresión bariónica real (McCarthy+2017, Chisari+2019).
Rango adicional N-body: ΔB_sigma8 ~ 0.03–0.07, llevando S₈^N-body ≈ 0.705–0.735.

Tensión DES proyectada:
- Optimista: (0.705 − 0.759)/0.023 = −2.36σ
- Conservador: (0.735 − 0.759)/0.023 = −1.05σ

**Falsificación:** Si N-body produce S₈ < 0.785 → OP-5 resuelto (<1.2σ DES).

**Recursos Nivel 2:** BAHAMAS-SSEE: ~5,000–10,000 CPU-horas (~USD 500–1,000);
IllustrisTNG-SSEE: ~10,000–20,000 CPU-horas (~USD 1,000–2,000).

**Script:** `src/ssee_op5_hmcode.py` (HMcode-2020 completo en CLASS, todos los pasos documentados)

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
- H₀,local = 67.96 / (1 − 0.06725) = 72.86 km/s/Mpc (0.17σ SH0ES)
- Con H₀^UV = 73.040 (Paper 10, condicional a Postulate C.1): 0.00σ SH0ES

**Cambios aplicados en Paper 9:**
- §3: Derivación desde universo separado k-essence (primer principios)
- §3: Identidad 1+w₀ = Ω_m,dyn explicitada como justificación de la cancelación
- Bibitem `wands2000` y `brax2014` añadidos

**Script:** `src/ssee_op6_screening_form.py` (verificación completa)

---

## OP-7 — QFT Derivation of Genesis Role Assignments (Transversal) — PARCIALMENTE RESUELTO

**Location:** Transversal — Papers 4, 7, 8 principalmente; también Papers 1, 9, 10

**Añadido:** 2026-05-19  
**Actualizado:** 2026-05-21

**Resolución parcial (2026-05-21):**

El argumento de unicidad estructural para $\beta_c = -\mathrm{AURA}$ está ahora
formalizado en dos lugares:

1. **Paper 7 §5.2** (Eq.~`eq:phi_pi_duality`, añadida en esta sesión): el párrafo
   "Discrete duality $\varphi\leftrightarrow\pi$" muestra que KAL₀ y AURA son imágenes
   duales bajo la $\mathbb{Z}_2$ discreta, y que por lo tanto $\beta_c = -\mathrm{AURA}$
   es la única elección consistente con esa simetría — no un parámetro libre.
2. **Paper 1 §5.3** (Eqs.~`eq:id_kal`–`eq:id_aura`): las identidades algebraicas exactas
   están enunciadas formalmente con nota explícita "The open question of whether this
   discrete symmetry is an exact symmetry of the full $P(X,\varphi)$ action is catalogued
   as OP-7."

El argumento de unicidad a nivel EFT está cerrado. Lo que permanece abierto es el nivel
QFT más profundo (puntos 1–3 abajo).

**Statement:**

SSEE Genesis 5.12 (sistema algebraico adimensional, preexistente a todos los papers)
asigna roles funcionales a las constantes algebraicas:

- **AURA = (3φ+π)/2** → "Portal / Contenedor dimensional" — límite umbral 1
- **KAL₀ = (φ+3π)/2** → "Retención / Anclaje Local" — sector cinético escalar

Los papers cosmológicos dan peso dimensional a estos roles vía observables:

- Paper 7: βc = −AURA para el acoplamiento disformal del fotón
- Paper 8: geodésicas disformales del fotón usan AURA
- Papers 1–6: el sector cinético escalar usa KAL₀ como constante de viscosidad estructural

**La dualidad φ↔π (establecida en esta sesión, 2026-05-19):**

La intercambiabilidad exacta es verificable numéricamente con error = 0 en doble precisión:

```
KAL₀|_{φ↔π} = AURA      (error = 0.00e+00)
AURA|_{φ↔π} = KAL₀      (error = 0.00e+00)
```

Los dos sectores son exactamente duales: el escalar φ usa la combinación π-dominada (KAL₀);
el fotón usa la combinación φ-dominada (AURA). La cadena cierra:

```
Ω_m,CMB = MIRA × Ω_m,dyn = (AURA/2) × (π−φ)/(2Ω) = AURA(π−φ)/(4Ω) = 0.31993
```

(0.63σ de Planck 2018: 0.3153±0.0073)

**El gap:**

La dualidad es estructuralmente consistente y numéricamente verificada, pero **no existe un
teorema que derive por qué la teoría cuántica de campos debe reproducir estas asignaciones
de roles** desde primeros principios.

Específicamente, falta demostrar:

1. **Por qué el sector disformal (fotón) debe acoplarse con AURA y no con KAL₀:** El
   argumento actual es la dualidad φ↔π — como el sector cinético usa la combinación
   π-dominada (KAL₀), el sector fotónico debe usar la φ-dominada (AURA). Esta es una
   restricción estructural plausible, pero no un teorema de QFT.

2. **Por qué βc = −AURA es la única solución:** Paper 7 establece αT=αM=αB=0 y deriva βc
   desde el Lagrangian P(X)=X/KAL₀+X²/M⁴. La derivación D=−P_{XX}/P_X²|_{X→0}=−2KAL₀²/M⁴
   requiere que βc = −AURA para satisfacer la dualidad — pero la dualidad misma no está
   derivada desde un principio de simetría de la acción.

3. **La conexión escala Planck → escala cosmológica:** El sistema Genesis 5.12 es
   adimensional. La dimensional weight se asigna vía cosmología (Papers 1–10). Falta el
   teorema que muestre que la física de QFT en la escala de Planck reproduce inevitablemente
   los roles de Genesis 5.12. Esta es la "Derivación Cuántica Dimensional" de la arquitectura.

**Consecuencia si permanece abierto:**

La predicción βc=−AURA es correcta (verificada por CAMB, CLASS, DESI a <0.2%) pero el
argumento de unicidad descansa en la dualidad estructural φ↔π, no en una simetría explícita
de la acción. Un referee podría aceptarlo como condición de consistencia pero no como
derivación. Es el gap más profundo del modelo.

**Programa de cierre (largo plazo):**

Identificar la simetría en el espacio de parámetros de KAL₀ / AURA bajo φ↔π como una
simetría discreta de la acción P(X,φ), y mostrar que esa simetría discreta — cuando se
impone a nivel de QFT — fuerza βc=−AURA como valor único del acoplamiento disformal.
Esto conectaría el sistema Genesis 5.12 con la física de campos desde primeros principios.

**Severidad:** Alta — es el gap conceptual más profundo. No falsifica SSEE a la precisión
observacional actual, pero es el requisito para que el modelo tenga derivación completa
desde primeros principios (necesario para aspirar a nivel de premio).

---

## OP-8 — MIRA Dynamical Mechanism (Papers 3, 5, 7, 8, 9) — ABIERTO

**Location:** Transversal — referenced in Papers 3, 5, 6, 7, 8, 9 as "the MIRA bridge."

**Problem:** The MIRA factor $(3\varphi+\pi)/4 \approx 1.9989$ is **algebraically exact**
and **observationally required** (CLASS validation 2026-05-09 showed CMB peaks shift
10% wrong without MIRA, RMS 31.5% vs 1.4%). However, **no dynamical mechanism derives
its value from first principles.** Seven candidate mechanisms have been tested and
ruled out (VERIFICATION_LEDGER §V-L3-mira, V-L3-saturacion, derivative coupling,
forward φ-MDE, etc.):

1. $c_s^2$ clustering — fails by magnitude (~100× too small)
2. $\mu>1$ Poisson modification — fails structurally (P7 forced $\alpha_B=\alpha_M=0$)
3. Disformal P8 — fails axiomatically (requires $\psi_{DM}$, P1 originally prohibited)
4. Conformal coupling $\beta_c = -$AURA — fails magnitude+sign+timing
5. Conformal coupling $\beta_c = \pm\alpha_{\rm sat}$ (veta-2) — fails connectivity
6. Derivative coupling $L_{\rm int} = (X/M^4) L_{DM}$ — fails UV scale
7. Forward integration from φ-MDE — fails matter drainage

**Current status:** MIRA enters the model as an **empirical input with exact algebraic
value**, not as a derived consequence. The framework's working interpretation
(two-sector dark matter, Paper 6) provides a phenomenological structure that
reproduces MIRA's numerical value via $\Omega_{m,{\rm CMB}}/\Omega_{m,{\rm dyn}} \approx 2$,
but this is a parametrization, not a dynamical derivation.

**Path to resolution:** A first-principles QFT/UV derivation showing why MIRA must
take this specific value. Possible candidates not yet tested: k-mouflage with
matter-coupled velocity-dependent screening, two-scale UV completion with discrete
self-similarity, holographic principle constraint.

**Severity:** High — this is the central open problem of SSEE-V3.6. Until resolved,
the model has ~3 effective free parameters (vs 6 in $\Lambda$CDM); when resolved,
it returns to ≤2 (only $H_0$ and $\Omega_b h^2$ as observation-tunable).

---

## OP-9 — Phenomenological Mass m_φ = 5.60 eV (Paper 6) — PARCIAL

**Location:** Paper 6, Eq. (mass_algebraic), §3.2.

**Problem:** The mass relation $m_\phi = \Sigma m_\nu^{\rm active} \times H_0^{\rm alg}
= 0.0824\,{\rm eV} \times 67.96 = 5.60\,{\rm eV}$ is **dimensionally inconsistent**
when treated as eV (the product has units [eV·s⁻¹·Mpc⁻¹]). A proper dimensional analysis
with $\hbar=c=1$ would yield $\sqrt{m_\nu \cdot H_0^{\rm alg}} \approx 10^{-17}$ eV
(Fuzzy DM regime), not 5.60 eV.

**Paper 6 honestly admits this (L371-377)** as a "numerological ansatz whose physical
origin (a Yukawa coupling through the IS sector) is to be derived in future work."

**Current status:** m_φ = 5.60 eV is adopted **phenomenologically** as the value that
gives the required free-streaming scale $k_{\rm fs} = 0.493\,h/$Mpc. The numerical
coincidence with Σm_ν × H_0 is presented as evidence of an underlying algebraic
structure to be derived.

**Path to resolution:** Derive m_φ from a proper field-theory mechanism — most likely
a Yukawa coupling between χ and the IS-relaxation operator, or from the curvature
of a unified potential $V(\phi)$ at its minimum (see OP-10).

**Severity:** Medium — the falsifiable prediction is $k_{\rm fs} = 0.493\,h/$Mpc
in DESI Y3 / Euclid, which is testable **independently of the mass formula's origin**.

---

## OP-10 — Unification of φ and χ into a Single Field (Papers 6, 7) — ABIERTO

**Location:** Paper 6 introduces χ (the DM scalar) as **distinct** from Paper 7's φ
(the DE k-essence scalar). L394-395: "We model φ-dark matter as a single real scalar
field $\chi$ (distinct from the SSEE background k-essence scalar $\phi$ of Paper 7)."

**Problem:** This historically arose post-hoc as a fix for the fσ8 tension in Paper 5.
The current potential $V(\phi) = V_0 \exp(-\alpha\phi)$ of Paper 7 has **no minimum**,
so the same field cannot also produce matter-mode oscillations. Hence χ was
introduced as a separate species with its own mass term.

**Consequence:** The "single field" philosophical claim of SSEE (one φ, derived from
$\varphi,\pi$) is technically violated — there are two cosmological scalar fields.
Both have parameters constrained by $(\varphi,\pi)$ but they are not the same entity.

**Path to resolution:** Find a potential $V(\phi)$ with:
- A slow-roll region (produces DE behavior at low energy)
- A stable minimum (produces matter-mode oscillation, hence DM)
- Both regions parametrized by $\varphi,\pi$ only

If such a $V(\phi)$ exists, χ becomes the matter-phase of φ itself, m_φ emerges as
the curvature at the minimum (resolving OP-9), and the "zero free parameters"
philosophy is restored.

**Severity:** Medium-High — does not affect observational predictions of SSEE-V3.6
but is the natural next step toward genuine zero-parameter status.

---

## OP-11 — Free Non-Minimal Coupling ξ (Paper 6) — ABIERTO

**Location:** Paper 6, Eq.~\eqref{eq:phiDM_lagrangian}, L406-407.

**Problem:** The action of χ contains a non-minimal coupling $\xi R \chi^2/2$ where ξ
is **explicitly a free continuous parameter** (P6 L416-417: "no free continuous
parameter other than ξ"). No algebraic constraint from $(\varphi,\pi)$ is provided
to fix ξ.

**Current status:** ξ is treated as a free parameter, partially constrained by the
production mechanism (gravitational particle production, OP-12).

**Path to resolution:** Either derive ξ from the underlying field theory (e.g., a
conformal symmetry requirement giving $\xi=1/6$, or a specific algebraic value), or
absorb χ into φ (OP-10) so that ξ becomes a structural piece of $V(\phi)$ rather
than a free input.

**Severity:** Medium — contributes one free parameter to the current model count.

---

## OP-12 — Relic Abundance Ω_φ-DM h² Not Computed Ab Initio (Paper 6) — ABIERTO

**Location:** Paper 6, §4.2 (gravitational particle production), L460-468.

**Problem:** P6 admits explicitly:
> "We have **not** computed the full gravitational-production integral that would
> yield $\Omega_{\phi DM} h^2$ from $m_\phi$ and $V(\phi)$ ab initio; an order-of-magnitude
> estimate is consistent with the required abundance, but the complete calculation
> (including the ξ-dependence) is deferred to a dedicated study (Paper B)."

**Current status:** Ω_φ-DM h² ≈ 0.074 is **imposed** to match the (MIRA−1)·Ω_dyn = 0.160
algebraic identity. The gravitational-production calculation that should derive this
abundance from inflationary parameters and m_φ is deferred.

**Path to resolution:** Carry out the full Parker-Kolb-Riotto calculation with the
SSEE α-attractor inflation potential and $m_\phi = 5.60$ eV. Should give
Ω_φ-DM h² ≈ 0.074 without additional tuning.

**Severity:** Medium — the relic abundance match is currently algebraic; the dynamical
calculation would close the gap.

---

## OP-13 — Inconsistencia interna Paper 8: §3-4 (factor √AURA) vs §4.5 (B-S) — ✅ RESUELTO 2026-05-23 (Opción A)

**Location:** Paper 8, §3-4 ("Disformal null geodesic" + "MIRA emerges in lensing") vs
§4.5 ("EFT suppression of fifth-force corrections") vs §6.2 item 3 (falsifiability).

**Diagnóstico (revisión 2026-05-23):**

P8 contiene dos derivaciones que asumen físicas incompatibles:

| Sección | Asunción física | Predicción |
|---|---|---|
| §3-4 | "DM efectiva sourced enteramente por baryonic seed" (escenario MOND-like, sin DM real) | $\theta_E^{\rm SSEE}/\theta_E^{\rm GR-bary} = \sqrt{\beta_c} \approx 2$ |
| §4.5 | DM real existe (Ω_CDM=0.16 + Ω_φDM=0.16), EFT con α_B=α_M=α_T=0 | $F_\phi/F_N \sim (H/k)^2 \sim 10^{-15}$ a escala kpc |
| §6.2 item 3 | (saca de §3-4) | $M_{\rm dyn}/M_{\rm lens} \approx 4$ |

**El propio paper admite** (P8 L363-370) que §3-4 es "working phenomenological estimate"
y difiere "full solution to future N-body simulations". Sin embargo §6.2 item 3 lo
presenta como predicción rígida falsable.

**Datos confrontados (revisión literatura 2026-05-23):**
- SLACS γ = 2.078 ± 0.027 (Auger 2010; Sonnenfeld 2015)
- M_E/M_dyn ≈ 1.207 (SIS, Cao et al. 2018, arXiv:1803.00819) — i.e., M_lens excede M_dyn ~20%
- Power-law profiles: M_E/M_dyn within 1σ de 1
- NO desviaciones factor 2 ni factor 4 reportadas en ningún sample

**Confrontación:**
- §3-4 predice $M_{\rm dyn}/M_{\rm lens}\approx 4$ → datos dan ≈ 0.83 (factor 5 wrong direction)
- §4.5 predice $M_{\rm dyn}/M_{\rm lens}\approx 1$ → datos compatibles (offset ~20% atribuible
  a contaminación línea-de-vista y errores de modelo)

**Veredicto:**
1. La física correcta para el modelo real (two-sector con DM existente) es §4.5
2. §3-4 describe un escenario MOND-like que NO corresponde al modelo SSEE canónico
3. La identidad "$\sqrt{\beta_c}\approx$ MIRA al 0.03%" es una near-coincidence
   numérica entre dos cantidades de regímenes incompatibles, no derivación física

**Implicaciones cross-paper:**
- P8 título/abstract: "MIRA emerges in lensing" pierde fuerza si §3-4 no es la física canónica
- P1 §1.4 / P8 introducción: la "unificación CMB↔lensing via MIRA" descansa en
  near-coincidence numérica, no en identidad estructural
- P6 y P7: no afectados directamente (no usan factor √β_c en lensing)

**Path to resolution:**
1. **Opción A (defensiva)**: Reformular P8 §3-4 como "alternative limit scenario", marcar
   §4.5 como predicción canónica. §6.2 item 3 ajustado a $M_{\rm dyn}/M_{\rm lens}\approx 1$
   (consistente con datos). Reescribir narrativa MIRA↔lensing como "near-coincidence
   numérica" en lugar de "emergencia estructural".
2. **Opción B (constructiva)**: Derivar consistentemente el lensing en el escenario
   two-sector real (con DM existente + EFT B-S), determinar si hay alguna firma
   observable distinta de ΛCDM, reformular P8 alrededor de eso.
3. **Opción C (radical)**: Aceptar que P8 §3-4 fue overclaim y reescribir el paper
   sin la afirmación factor 2.

**Severity:** **High** — afecta la narrativa central de P8 ("MIRA emerges in lensing")
y la cadena de argumentos de unificación. NO afecta la solidez del background
SSEE (Papers 1-5, 7, 9, 10), solo el régimen de gravedad fuerte (Paper 8).

**Cómo este OP escapó de la auditoría previa:** la auditoría hostil chequeó consistencia
**dentro** de cada sección y entre paper↔scripts↔ssee_core, pero no consistencia
**entre secciones del mismo paper** ni paper↔datos publicados. Patch metodológico:
añadir Capa 8 "consistencia inter-sección" al guardián.

---

### Resolución aplicada (Opción A — 2026-05-23)

**Insight clave (correctión de terminología señalada por Mike):** el factor
derivado en P8 §3-4 NO es $\MIRA$ — es $\sqrt{\AURA}$. Son dos operaciones
algebraicas distintas sobre $\AURA$ (división por 2 vs raíz cuadrada). Coinciden
numéricamente al 0.03% porque $\AURA\approx 4$ hace que ambas rondan 2. Drafts
previos de P8 conflated los dos bajo la etiqueta "$\MIRA$ emerges in lensing";
esa identificación se retira.

**Cambios aplicados (commit de la sesión):**

1. **P8 abstract reescrito**: framing "dos límites" explícito (alternative MOND-like
   vs canonical two-sector); $\sqrt{\AURA}$ en lugar de "√β_c ≈ MIRA"; canonical
   limit prediction declarada como $\thetaE^{\rm SSEE}\approx\thetaE^{\rm GR-with-DM}$.

2. **P8 §1 (Introduction) reescrito**: dos limits identificados explícitamente;
   $\MIRA$ aparece estructuralmente solo en CMB y LSS, no en lensing; aclaración
   $\sqrt{\AURA}\neq\MIRA$ algebraicamente.

3. **P8 §3-4 (Disformal geodesic)** ahora titulado "Alternative MOND-like limit";
   tcolorbox al inicio aclarando scope; sección retenida como pedagógica/histórica,
   no canónica.

4. **P8 §5 (formerly "MIRA emerges in lensing")** re-titulado "Lensing factor
   $\sqrt{\AURA}$ in the Alternative Limit (and Its Numerical Near-Coincidence with
   $\MIRA$)"; tcolorbox naming clarification explícito al inicio.

5. **P8 §6.2 items 2 y 3 (Falsifiable predictions)** reformulados a predicción
   canónica: $\thetaE^{\rm SSEE}\approx\thetaE^{\rm GR-with-DM}$ y
   $M_{\rm dyn}/M_{\rm lens}\approx 1$, consistente con SLACS/BELLS observados
   ($\gamma=2.078\pm 0.027$; $M_E/M_{\rm dyn}\approx 1.21$ SIS). Footnote retira
   explícitamente la predicción "$\sim 30\%$ deficit" previa.

6. **P8 §8.1 (Discussion)** reescrita: $\MIRA$ aparece exactly en 2 sectores
   (no 3); la aparición en lensing era name-conflation de $\sqrt{\AURA}$ con
   $\MIRA$ vía near-coincidence numérica, retirada en esta versión.

7. **Unified Journal §2 + tabla predicciones**: framing dos-límites; entry de
   lensing marcado como "Conditional prediction (not canonical; OP-13)".

8. **Endorser Summary §IV (Paper 8 bullet)**: framing dos-límites; SLACS/BELLS
   citado como evidencia del canonical limit.

**PDFs regenerados:** P8 (18 pp, +1), Unified (21 pp), Endorser (2 pp).
Guardián: VERDE 102/102.

**Estado del modelo tras resolución:**
- $\MIRA$ aparece exactly en 2 sectores estructurales: CMB horizon mapping
  (Paper 3) y two-sector mass ratio (Paper 6). Ambos: $\MIRA=\AURA/2$ exacto.
- $\sqrt{\AURA}$ aparece en lensing alternative-limit (P8 §3-4) — preservado como
  derivación válida pero marcado explícitamente como no canónico.
- Canonical SSEE lensing prediction: $\thetaE^{\rm SSEE}\approx\thetaE^{\rm GR-with-DM}$
  (de EFT B-S structure P8 §4.5), consistente con datos SLACS/BELLS.
- La "unificación CMB↔DM↔lensing via $\MIRA$" se degrada a "unificación CMB↔DM
  via $\MIRA$ exacta; lensing es near-coincidence numérica con $\sqrt{\AURA}$".

**Lo que NO cambia:**
- Postulado M ($\MIRA=(3\phiG+\pi)/4$) sigue siendo central del framework.
- Papers 1-7, 9, 10 no afectados.
- EFT canonical (P7), Hubble tension (P9), UV completion (P10): intactos.

**Lo que cambia (honesto):**
- P8 ya no es "MIRA-in-lensing paper" sino "two-limit analysis paper": el factor
  $\sqrt{\AURA}$ del alternative MOND-like limit es preservado como derivación
  pedagógica; la predicción canónica observable es GR-with-DM.
- El claim "MIRA aparece en 3 sectores" se retira de la narrativa cross-paper.

---

## OP-14 — Σm_ν Phenomenological Derivation (Paper 4) — ABIERTO (atacado 2026-05-23, blocked-by-OP-10)

**Location:** Paper 4, §"Neutrino Mass Sum", L675-700.

**Formula:**
$$\Sigma m_\nu^{\rm active} = \mathcal{R}\,\Omega_b h^2 \cdot \frac{94.07~\mathrm{eV}}{\tau_\Pi H_0} = 0.0824~\mathrm{eV}$$

con $\mathcal{R} = 4\cdot\mathrm{KAL} - 22 = 0.0856$.

**Problem:** The paper itself classifies this as **Type P (phenomenological)** at L691-697,
explicitly admitting it "has no independent φ,π derivation." Two weak elements:

1. **Integer subtraction $4\cdot\mathrm{KAL} - 22$**: el offset 22 no se deriva — se ajusta
   numéricamente para reproducir la cota cosmológica $\Sigma m_\nu \lesssim 0.12$ eV.
2. **94.07 eV** es input del Modelo Estándar (relación entre densidad relíquica y masa
   de neutrino), no de SSEE.

**Why this matters (cascada):** Σm_ν alimenta directamente a $m_\phi = \Sigma m_\nu \cdot H_{\rm alg}$
en P6, así que OP-9 depende de OP-14. **OP-14 es el eslabón más débil de la cadena
de derivaciones.** Si se deriva $\mathcal{R}$ desde primeros principios, OP-9 cae casi
automáticamente.

**Path to resolution:**
- Derivar el offset 22 desde teoría de Genesis Roles (¿conteo de grados de libertad?
  ¿índice topológico?).
- Reformular como cota saturada de una desigualdad física estándar (análogo al método
  veta-2 de saturación que dio α_sat = √(3/(φ+3π)) — ver memoria
  [[project-veta2-saturation]]).
- Alternativa: aceptar Σm_ν como input experimental SM y solo derivar la **razón**
  $m_\phi/\Sigma m_\nu = H_{\rm alg}$.

**Severity:** Medium-High — bloquea el camino a cero parámetros porque rompe la
narrativa "todo emerge de φ,π" en P6 y P4 simultáneamente.

### Ataque ejecutado (2026-05-23) — script `src/ssee_op14_neutrino_mass.py`

Tres hipótesis testadas:

**H1 — Fragilidad estructural: CONFIRMADA**

| Perturbación δ(φ,π)/x | Drift en Σm_ν |
|---|---|
| 10⁻⁶ | −0.2% |
| 10⁻⁴ | +2.4% |
| 10⁻³ | **+25.5%** |
| 10⁻² | +257% |

La forma $\mathcal{R} = 4\cdot\text{KAL} - 22$ es una diferencia entre números casi
iguales (22.086 vs 22.000). Pequeñas perturbaciones en (φ,π) se amplifican
relativamente. Una predicción genuinamente estructural debería ser estable bajo
perturbaciones de orden 10⁻³ — esta no lo es.

**H2 — ¿22 = conteo físico de DoF?: SIN MATCH**

- SM tiene 12 generadores gauge (≠22), 28 parámetros libres (≠22), 13 bosones (≠22).
- "22 = 2 × 11 (dim M-theory)" o "22 = 2·rank(E11)" son post-hoc sin argumento físico
  independiente.

**H3 — Scan algebraico (~150 monomios): SIN IDENTIDAD EXACTA**

Mejor candidato: $\mathcal{R} \approx 1/[3(\text{KAL}-\varphi)] = 2/[3(3\pi-\varphi)]$
con error de **−0.27%**. Más limpio que "$4\text{KAL}-22$" pero sigue siendo
aproximación numérica, no identidad. Si fuera estructural el error sería 0%.

### Conclusión: OP-14 colapsa con OP-9

Reformulación crítica: $\Sigma m_\nu$ y $m_\varphi$ comparten el **mismo grado de
libertad fenomenológico** vía la identidad de Paper 6:
$$m_\varphi = \Sigma m_\nu \cdot H_0^{\rm alg}$$

OP-14 y OP-9 NO son problemas independientes — son **dos caras del mismo
problema**. Atacar OP-14 directamente no produce derivación porque no hay
ataque local: la única salida es derivar $m_\varphi$ desde la curvatura
de un potencial $V(\varphi)$ fundamental.

### Nueva cadena de dependencias (post-ataque)

```
OP-10  V(φ) unificador DE + DM
   │
   ├──► OP-9   (m_φ = curvatura V en mínimo)
   │       │
   │       └──► OP-14 (Σm_ν = m_φ / H_alg)    ← invertido: ya no flecha OP-14→OP-9
   │
   ├──► OP-11  (ξ acoplamiento → funcional de V)
   │
   └──► OP-12  (Ω_φDM h² desde dinámica de V)

OP-8   MIRA mecanismo (independiente, eslabón duro)
```

Si OP-10 cae, **4 OPs caen juntas** (9, 11, 12, 14). El modelo pasa de 4 a 3
postulados (D, S, I; M sigue abierto vía OP-8).

### Acción inmediata sobre Paper 4

El §"Neutrino Mass Sum" ya admite Type P en L704-711. **No requiere edición
urgente** — el lenguaje actual es honesto. La acción es **catalogar OP-14 como
blocked-by-OP-10** en el registro de OPs y dirigir el ataque siguiente a OP-10
(Fase 1: catálogo de candidatos $V(\varphi)$).

---

## Summary Table

| ID | Paper | Problem | Severity | Path to Resolution |
|----|-------|---------|----------|--------------------|
| OP-1 | P4 | ~~Factor 200 in Ω_b h²~~ | ✅ PARCIAL | (π−φ)/H₀_SSEE=0.32σ Planck; BBN derivation → Paper B/C; script op1 |
| OP-2 | P4 | ~~n_s exponent 7 not derived from V(φ)~~ | ✅ RESUELTO | α-attractor universality + N_*=2φ⁷; r=φ⁻¹⁰ nueva predicción; script op2 |
| OP-3 | P10 | ~~UV-IR separability unproven~~ | ✅ RESUELTO | Jerarquía EFT (H₀/M)²≈10⁻⁶² + KALeff=φ²√(5/2) único; Paper 10 TC.1; script op3 |
| OP-4 | P8 | ~~r_V > r_Hubble para Vainshtein~~ | ✅ RESUELTO | k-mouflage + αB=αM=αT=0 EFT; Paper 8 §4.2/§4.4 revisados |
| OP-5 | P5-6 | ~~S₈ 2.29σ DES (fσ₈ resuelto P6)~~ | ✅ PARCIAL | HMcode-2020 CLASS: S₈=0.758 (0.06σ DES); N-body full → Paper B/ext |
| OP-6 | P9 | ~~Screening form ambiguity~~ | ✅ RESUELTO | Universo separado k-essence + identidad 1+w₀=Ω_m; Paper 9 §3 revisado |
| OP-7 | P4/7/8 | QFT derivation of Genesis role assignments | ✅ PARCIAL | EFT uniqueness formalizado P7 §5.2 + P1 §5.3; QFT desde primeros principios → largo plazo |
| OP-8 | Transv. | **MIRA dynamical mechanism** | High | 7 mechs ruled out (Ledger §V-L3); future: k-mouflage matter-coupled, holographic |
| OP-9 | P6 | m_φ=5.60 eV is numerological (dim. inconsistent) | Medium | Yukawa/IS-coupling derivation; or absorb into OP-10 |
| OP-10 | P6/P7 | Unify χ into φ via richer V(φ) | Medium-High | V(φ) with slope (DE) + minimum (DM matter-mode); restores zero-param status |
| OP-11 | P6 | ξ (non-minimal coupling) is free parameter | Medium | Algebraic constraint or absorb via OP-10 |
| OP-12 | P6 | Ω_φ-DM h² not computed ab initio | Medium | Parker-Kolb-Riotto with α-attractor + m_φ |
| OP-13 | P8 | ~~Contradicción interna §3-4 vs §4.5~~ | ✅ RESUELTO | Opción A aplicada: framing dos-límites, retirado claim "MIRA en lensing", $\sqrt{\AURA}$ ≠ $\MIRA$ aclarado, canonical prediction = GR-with-DM (2026-05-23) |
| OP-14 | P4 | Σm_ν phenomenological (Type P admitted L691-697); offset 22 ad hoc | Medium-High | **Atacado 2026-05-23** (script op14): forma frágil (1‰→25% drift), sin DoF match, scan sin identidad exacta. **Blocked-by-OP-10**: $\Sigma m_\nu$ y $m_\varphi$ comparten DoF, derivar V(φ) los resuelve juntos |

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
SSEE-V4.0; OP-8 through OP-13 constitute the agenda for SSEE-V5.0 (full unification
of the dark sector).

## Parameter-Count Status (2026-05-22)

SSEE-V3.6 currently has **~3 effective free parameters** vs **6 for $\Lambda$CDM**:

| Parameter | Status | Comment |
|---|---|---|
| $H_0$ | Sampled in MCMC | Algebraic prediction 67.96 vs posterior 66.75 |
| $\Omega_b h^2$ | Sampled in MCMC | Algebraic 0.02242 vs Planck 0.02237 |
| MIRA factor | Algebraic value $(3\varphi+\pi)/4$ but no derived mechanism | OP-8 |
| $m_\phi$ | Phenomenological ansatz | OP-9 (paper admits) |
| $\xi$ | **Free parameter** (P6 L416) | OP-11 |
| $\Sigma m_\nu$ | **Phenomenological** (P4 L691-697 admite Type P) | offset 22 ad hoc — OP-14 |
| $\alpha$ (attractor) | $\varphi^4/3$ derived | from $\varphi$ |
| $V(\phi)$ form | Adopted exponential | locks DE-only behavior — OP-10 |

**Path to zero parameters:** Resolve OP-8 (MIRA derivation) + OP-10 (unify χ into φ
via richer $V(\phi)$) + OP-14 (Σm_ν offset 22). This would also resolve OP-9 (m_φ =
curvature of V at min) and OP-11 (ξ becomes part of V structure). End state: $H_0$
and $\Omega_b h^2$ as the only observation-tunable inputs.

**Dependency chain (cascada):** OP-14 → OP-9 → (parte de) OP-10. Atacar OP-14
primero es el camino más corto: si se deriva $\mathcal{R}=4\cdot\mathrm{KAL}-22$,
$m_\phi$ pasa de ansatz a predicción automáticamente.
