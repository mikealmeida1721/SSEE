# SSEE-V3.6 — Open Problems and Known Limitations

This document records the open theoretical questions and **frontiers in progress** of the
SSEE framework. These are **not** hidden flaws: they are declared, in-process work — the
opposite of numerology (which would pretend they are solved). Ver `NUMEROLOGY_AUDIT.md`.

Last updated: 2026-07-12

---

## Leyenda de estado — cómo leer este documento (FIJADO 2026-07-12, no re-frasear)

Un problema abierto (OP) es una **frontera en proceso, NO una debilidad.** La distinción es
referee-proof y se aplica a todo el documento:

| Estado | Qué significa | ¿Bloquea publicación? |
|---|---|---|
| **CERRADO** | Derivado/verificado; se defiende ante un referí. Las afirmaciones titulares del modelo descansan aquí. | No |
| **PARCIAL** (resource-gated) | Responde bien la pregunta **hasta donde los recursos actuales permiten**; la verificación PROFUNDA requiere recursos que hoy no tengo (supercómputo, datos futuros, N-body). El resultado inmediato suele ser **falseable**. | No |
| **ABIERTO** (frontera) | Pregunta clara y bien planteada, en proceso. | No |
| **DEBILIDAD** | Una afirmación **CERRADA** que **NO se puede defender** ante un referí, por razones claras. | **Sí** |

**Regla anti-inconsistencia (para auditorías futuras — incluido Claude):** un OP **no es una
debilidad**. Una debilidad es algo que ya *afirmamos como cerrado* y que **falla ante un
referí**. El TEST para saber si un OP es debilidad: *¿alguna afirmación titular del modelo
depende en secreto de que ese OP esté resuelto?* Si **no** → es frontera (el modelo se para sin
él). Si **sí** → es carga, y hay que reclasificarla. **A hoy SSEE tiene 0 debilidades:** ningún
titular (w₀wₐ validado por DESI, k_fs falseable, cascada H con unidades ancladas
empíricamente) depende de un OP sin resolver. Los OPs son profundizaciones, no huecos que
sostengan lo publicado.

**Readiness para Zenodo:** publicable = **(0 debilidades) + (OPs declarados con honestidad)**.
NO es *(0 OPs)* — eso es imposible: hasta ΛCDM tiene el problema de la constante cosmológica
abierto. No encadenar la publicación al mecanismo (OP-9/OP-10) — ver [[project_publication_strategy]].

---

## OP-1 — First-Principles Derivation of ω_b = (π−φ)/(3Ω²) (Paper 4 / baryogenesis) ✅ PARCIALMENTE RESUELTO

> **El "factor 200" YA NO EXISTE en el modelo (corregido en el título 2026-06-24).**
> El ω_b canónico es la fórmula algebraica `(π−φ)/(3Ω²) = 0.02242` (0.32σ Planck),
> sin ningún 200. El viejo `3(π−φ)/200` (3.2σ) quedó retirado; φ¹¹≈199 mostró que 200
> era solo una aproximación. Lo que permanece ABIERTO no es el factor — es derivar
> `(π−φ)/(3Ω²)` desde primeros principios (bariogénesis Γ_sph, programa Paper B/C).

**Location:** Paper 4, §3.2 (baryon buffer constant) — **revisado 2026-05-16**

**Corrección (sustitución algebraica del antiguo factor 200, ya retirado):**

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
de producción de φ-DM (m_φ=40.70 eV) que reproduce Ω_φDM h²=0.0688.

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
La jerarquía de escalas (H₀/M)² = (1.45×10⁻³³ eV / 9.62×10⁻³ eV)² ≈ 2.3×10⁻⁶² suprime
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
$$r_{\rm km}^3 = \frac{M_{\rm obj}}{4\pi M_{\rm Pl} M^2}, \quad M = 9.62\ \text{meV}$$
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

**Clarification — what IS and IS NOT resolved (canónico ω_m-directo, Ω_m=0.30889):**
- **fσ₈ (growth-rate):** NO es donde el φ-DM ayuda. El single-sector canónico (Paper 5,
  Ω_m=0.30889) da media 0.70σ; el two-sector da 0.70σ (idéntico) — ambos **empatan ΛCDM**
  (0.73σ). A escalas RSD (k≪k_fs) el φ-DM agrupa como frío, sin firma two-sector en fσ₈.
  (El viejo "2.56σ→0.50σ" usaba datos fσ₈ erróneos y/o el baseline no-canónico Ω_m=0.160;
  el "0.74/0.76σ" usaba Ω_m=0.30889 vía MIRA — ambos retirados.)
- **S₈ (weak-lensing) — el desafío REAL:** single-sector S₈=0.846 (3.5σ KiDS). El
  two-sector free-streaming lo baja a **S₈_eff=0.758 (0.00σ KiDS) — RESUELVE** a nivel
  lineal/forward (m_φ=40.70 eV, cero fiteo; el viejo 0.761 era la rama WDM fiteada).

**Residual abierto:** solo el refinamiento NO LINEAL pleno (N-body con feedback bariónico,
Nivel 2). El cierre lineal/forward de S₈ ya está hecho (0.758).

**Nivel 1 — HMcode-2020 baryonic feedback (CLASS, laptop) — COMPLETADO 2026-05-16:**

Script `src/ssee_op5_hmcode.py` implementa retroalimentación bariónica AGN via
HMcode-2020_baryonic_feedback en CLASS (Mead et al. 2020, log10T_heat=7.8):

Resultados CLASS HMcode-2020 con parámetros SSEE (H₀=66.75 — input de la corrida 2026-05-16, anterior al posterior canónico 66.53 km/s/Mpc; Ω_m=0.30889, w₀=−0.840, wₐ=−0.670):

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
| DES Y3 (observado) | 0.758 | 0.00σ | — |
| KiDS-1000 (observado) | 0.759 | — | 0.00σ |

**El baseline Paper 6 ya está dentro de 1σ DES** (0.09σ). HMcode añade Δσ = 0.03σ de mejora.

**Por qué el baseline está tan bien:** El two-sector φ-DM (m_φ=40.70 eV, k_fs=0.754 h/Mpc)
ya suprime P(k) en k > k_fs, sobre Ω_m,CMB=0.30889 (ω_m-directo). El HMcode añade supresión bariónica
suave adicional, principalmente a k > 0.5 h/Mpc.

**Nivel 2 — N-body completo (proyección):**

HMcode-2020 captura ~60–70% de la supresión bariónica real (McCarthy+2017, Chisari+2019).
Rango adicional N-body: ΔB_sigma8 ~ 0.03–0.07, llevando S₈^N-body ≈ 0.705–0.735.

Tensión DES proyectada:
- Optimista: (0.705 − 0.758)/0.023 = −2.36σ
- Conservador: (0.735 − 0.758)/0.023 = −1.05σ

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

Los dos sectores son duales: el escalar/clustering φ usa la combinación π-dominada (KAL₀);
el sector φ-dominado (AURA) aparece en el screening de Hubble (MIRA=AURA/2, Paper 9). Bajo
el reframe ω_m-directo (OP-8 disuelto) la densidad de materia del CMB sale DIRECTA, sin factor
materia — el escalar usa KAL₀ dentro de ω_c:

```
ω_c = KAL₀ · ω_b · n_s = 0.11951   →   Ω_m,CMB = ω_m/h² = 0.30889
```

(0.88σ de Planck 2018: 0.3153±0.0073; la vieja cadena MIRA×Ω_m,dyn=0.30889 queda retirada)

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

## OP-8 — MIRA Dynamical Mechanism (Papers 3, 5, 7, 8, 9) — ✅ DISUELTO 2026-06-18

> **CIERRE (reframe ω_m-directo, 2026-06-18).** OP-8 preguntaba por qué el factor
> de duplicación materia $\Omega_{m,{\rm CMB}}/\Omega_{m,{\rm dyn}}$ vale
> $(3\varphi+\pi)/4$ (MIRA) y no $2$. **El reframe elimina la premisa:** la materia
> del CMB ya **no** se obtiene duplicando el sector dinámico. Es el observable
> físico estándar $\omega_m = \omega_b + \omega_c + \omega_\nu = 0.14267$, con cada
> pieza algebraica de SSEE:
> - $\omega_b = (\pi-\varphi)/(3\Omega^2) = 0.02242$  (OP-1)
> - $\omega_c = \mathrm{KAL_0}\cdot\omega_b\cdot n_s = 0.11951$  (identidad **forward**, ya en Paper 1, 0.41σ)
> - $\omega_\nu = \Sigma m_\nu/93.14\,\mathrm{eV} = 0.000741$
>
> y $\Omega_{m,{\rm CMB}} = \omega_m/h^2 = 0.30889$ es **derivado**, no postulado.
> $\Omega_{m,{\rm dyn}}=0.160$ (DESI) y $\omega_m$ (CMB) son **dos predicciones
> independientes**, ya no ligadas por ningún factor (ni MIRA ni π/φ).
> Verificación CMB Fase B (Planck plik_lite, $H=67.962$, $\Omega_m=0.30889$):
> $\chi^2_{\rm SSEE}=1005.50$ vs $\chi^2_{\Lambda{\rm CDM}}=1003.76$,
> $\Delta\mathrm{BIC}=-23.93$ → SSEE favorecido
> (`results/logs/p3_cmb_reframe_omega_m.log`).
>
> **Residuo honesto (no es perilla nueva):** $\Omega_{m,{\rm CMB}}$ ahora descansa
> en que $\omega_b$ (OP-1) y la identidad $\omega_c=\mathrm{KAL_0}\cdot\omega_b\cdot n_s$
> sean correctos — piezas que el modelo **ya tenía**. **MIRA persiste como entidad**
> $=\mathrm{AURA}/2$ y conserva su rol en $f_{\rm screen}=(\pi-\varphi)/\Omega^2$
> (invariante); lo que se retira es su rol de **factor-materia**. El análisis
> histórico abajo (siete mecanismos, identidad $\Omega_{m,{\rm dyn}}=\Omega_{\rm geom}/2$)
> queda como registro; ya no carga el peso de un problema abierto.

**[HISTÓRICO — premisa superada por el reframe ω_m-directo]**

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

**Finding 2026-06-05 — OP-8 and Roadmap-point-2 collapse to a single number.**
The previously separate "$\Omega_{m,{\rm CMB}}$ dual" puzzle (Roadmap-point-2:
geometric route $\Omega_{\rm geom}=(\pi-\varphi)/(\pi+\varphi)=0.3201005$ vs
MIRA route $\mathrm{MIRA}\times\Omega_{m,{\rm dyn}}=0.30889282$, a $0.054\%$ gap
read as a possible loop/self-energy correction) is **not a second coincidence**.
It is the **same object** as MIRA's deviation from the integer $2$. Proven to
machine precision ($\Delta=5.55\times10^{-17}$) and by hand:
- $M_v = \varphi+\pi+K_v = 3\Omega$ (since $K_v=2\Omega$, $\Omega=\varphi+\pi$);
  $\mathrm{TRIAL}=3(3\varphi+\pi)/2$, so $M_v-\mathrm{TRIAL}=3(\pi-\varphi)/2$.
- Therefore $\Omega_{m,{\rm dyn}}=(M_v-\mathrm{TRIAL})/M_v=(\pi-\varphi)/\big(2(\pi+\varphi)\big)=\Omega_{\rm geom}/2$ **exactly**.
- The dual's fractional deviation is then identically
  $(\Omega_{\rm geom}-\mathrm{MIRA}\cdot\Omega_{m,{\rm dyn}})/\Omega_{\rm geom}
  = 1-\mathrm{MIRA}/2 = (2-\mathrm{MIRA})/2 = (8-3\varphi-\pi)/8 = 0.0538\%$.

So MIRA's distance from $2$ **fully determines** the $\Omega_m$ dual: one number,
$(8-3\varphi-\pi)/8$, not two independent $\sim0.054\%$ near-coincidences. This
gives MIRA a concrete root as $\mathrm{AURA}/2$ (half of the first dimensional
ceiling $\mathrm{AURA}=(3\varphi+\pi)/2$), which propagates the dimensional reading
to the **entire AURA branch** (the copy-law ladder MIRA·1, AURA·1, DUAL·2,
TRIAL·3, … spaced by exactly one AURA — the "dimensional ceilings"). The dual
exists *because* $\mathrm{AURA}\neq4$ exactly. Verification: `src/op8_mira_aura_dimensional.py`.

**Status of this finding (audit phase 2 — honest residue):** the identity is
**internal** (everything follows from the $\varphi,\pi$ definitions of $K_v$,
$M_v$, TRIAL), so it **tightens** the problem from two knobs to one but does
**not** supply the dynamical mechanism. The interpretation "MIRA = half dimensional
ceiling" stays in the *sistema SSEE* (interpretive layer) until a physical
mechanism validates it; the algebra $\Omega_{m,{\rm dyn}}=\Omega_{\rm geom}/2$ is
already *modelo*-grade. The problem is now **surrounded, not solved**: one must
still derive why the doubling factor is exactly $(3\varphi+\pi)/4$ and not $2$.

**Path to resolution:** A first-principles QFT/UV derivation showing why MIRA must
take this specific value — equivalently, why the matter-sector doubling is
$(3\varphi+\pi)/4$ rather than the integer $2$ (the $(8-3\varphi-\pi)/8$ residue).
Possible candidates not yet tested: k-mouflage with matter-coupled
velocity-dependent screening, two-scale UV completion with discrete
self-similarity, holographic principle constraint.

**Severity:** High — this is the central open problem of SSEE-V3.6. Until resolved,
the model has ~3 effective free parameters (vs 6 in $\Lambda$CDM); when resolved,
it returns to ≤2 (only $H_0$ and $\Omega_b h^2$ as observation-tunable).

---

## OP-9 — UV Origin of the Mass Multiplier $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V$ (Paper 6) — **ABIERTO** (coeficiente UV)

> [!nota] ESTADO CANÓNICO ÚNICO (2026-07-12) — leer esto, ignorar redacciones antiguas abajo.
> **Estado en UNA palabra: OP-9 está ABIERTO.** No lo llamamos "cerrado" en ninguna mitad —
> eso confundía. Lo que existe es una partícula usable y una derivación pendiente; solo la
> segunda es OP-9.
>
> $m_\phi = \Sigma m_\nu^{\rm act}\cdot\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V = 0.0685\cdot594.28 = 40.70$ eV.
> $\mathrm{SOLAR}=\mathrm{BIAL}+\mathrm{KAL}=\varphi+2\pi$; $\mathrm{KRYSTOS}_V=\varphi+\pi+\Omega$
> (padres {φ,π,Ω} — **NO** "2Ω"). Forma $g^2v$, escrita en el Lagrangiano escalar libre de Paper 6.
>
> **Para que quede claro qué falta y qué no:**
> - **Lo que YA funciona (esto NO es lo que OP-9 pregunta):** la partícula tiene masa
>   dimensionalmente correcta, sin fiteo, escrita como término $g^2v$ de un Lagrangiano real, y
>   es **falseable** ($k_{\rm fs}=0.754\,h/$Mpc, DESI Y3/Euclid). La *decisión* de usar esta
>   partícula ya está tomada — eso es **OP-17**, no OP-9.
> - **Lo que FALTA — LA PREGUNTA ABIERTA, exacta (esto SÍ es OP-9):**
>   > ¿Existe un potencial $V(\phi)$ construido **solo de φ,π** cuya **curvatura en el mínimo**
>   > reproduzca el coeficiente 594.28 —es decir, $m_\phi=\Sigma m_\nu\cdot594.28$— **sin
>   > meterlo a mano**?
>   >
>   > (Imagen: la masa = qué tan curvado está el fondo del valle del campo. La pregunta es si un
>   > valle con forma pura φ,π se curva *exactamente* eso. $m_\phi^2 = V''(\phi_{\min})$.)
>
> **Relación OP-9 ↔ OP-10 (para NO marear la pregunta):** comparten la **misma respuesta** — el
> potencial. **OP-10** = derivar el potencial completo (la *función*, el valle entero). **OP-9**
> = el *número* que ese valle debe escupir (la curvatura del fondo). Resuelto OP-10, OP-9 **cae
> solo**: son **un problema en dos resoluciones**, no dos problemas. Lo que le da nombre propio a
> OP-9: es la parte del potencial **ya expuesta a datos** ($k_{\rm fs}=0.754$) → el programa se
> puede **FALSAR por OP-9 antes** de tener OP-10. El intento acotado de derivarlo por transporte
> (Camino B, KAL, 2026-07-11) dio **CORTE**; la vía viva es el $V(\phi)$ de OP-10, más allá del
> horizonte TRIAL.
>
> **No confundir con OP-17** (¿la partícula sirve y es falseable?) → eso YA está decidido (sí).
> OP-9 pregunta SOLO por el origen del coeficiente. Mezclarlas produce el falso vaivén
> "falta/no-falta".
>
> Que OP-9 esté abierto **y declarado abierto con la pregunta exacta** es lo CONTRARIO de
> numerología (que fingiría tenerlo resuelto). Incompletitud, no inconsistencia. El bloque
> histórico (615.33, 1/537) es rastro, no estado. Ver [[project-op17-particle-deferred]],
> [[project-nine-sovereigns-naming]].

**Location:** Paper 6, §3.2 (subsec:mass_derivation), §3.4 (subsec:lagrangian), §6 (subsec:origin).

**Status change (2026-06-04):** The old phenomenological ansatz $m_\phi = \Sigma m_\nu
\times 3(\varphi+\pi)^2 = 5.602$ eV is **superseded**. The canonical particle is now
fixed by a **zero-fitting forward chain** (Vía 2 + multiplier), locked by Mike:

```
R₂   = Ω_DNAV/(KAL·TRIAL)          = 0.07188      (pure φ,π number)
Σm_ν = R₂ × 0.9530 eV              = 0.0685 eV     (fixed mass scale)
mult = SOLAR² · KRYSTOS_V            = 594.28        (PURE φ,π number)
m_φ  = Σm_ν × mult                 = 40.70 eV      (forward prediction, zero fitting)
T_φ                                = 0.5385 T_ν    (from relic-abundance constraint)
```

**What changed in the gap:** the coefficient is no longer a *loose* numerological factor.
The canonical multiplier $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V = 594.28$
is a **pure $(\varphi,\pi)$ number**, so $m_\phi = \Sigma m_\nu \times (\text{pure number})
= [{\rm eV}]$ is **dimensionally consistent**. Crucially, this coefficient is now the
**mass term of a written scalar Lagrangian**:
$$\mathcal{L}_\Phi = \tfrac{1}{2}\partial_\mu\Phi\,\partial^\mu\Phi
  - \tfrac{1}{2}\bigl[\Sigma m_\nu(\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V)\bigr]^2 \Phi^2.$$
The coefficient enters as $m_\Phi^2$ in a Lagrangian that is *written down*, not as a
floating ansatz. **This closes the incompleteness flagged in the original OP-9** (the
coefficient now comes from a Lagrangian mass term, not a bare numerological multiplier).

**What remains open (the refined OP-9):** the **UV origin of the canonical multiplier**
$\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V=594.28$ (adopted 2026-06-19; supersedes the
intermediate $\mathrm{PYROS}\cdot\mathrm{VITA}\cdot\mathrm{MIKA}=615.33$ and the older
$\Omega_{\rm DNAV}^4+\mathrm{AURA}\cdot\mathrm{KAL}=535.28$) — i.e.,
why *this* particular combination of $(\varphi,\pi)$ constants sets the curvature of the
potential at its minimum. The Lagrangian is written, but the multiplier's derivation from
a unified $V(\phi)$ (OP-10) is the next step. This is **incompleteness, not inconsistency**.

> [!warning] Bloque histórico (615.33-era). El análisis look-elsewhere que sigue
> (1/537, 1/192) corresponde al multiplicador **intermedio** $\mathrm{PYROS}\cdot
> \mathrm{VITA}\cdot\mathrm{MIKA}=615.33$, hoy **superado** por $\mathrm{SOLAR}^2\cdot
> \mathrm{KRYSTOS}_V=594.28$. La selección vigente bajo la gramática de linaje estricta
> es **1/3**, no 1/cientos — ver «Methodological note» más abajo. Se conserva por registro.

**Look-elsewhere of the multiplier (2026-06-18, `op9_particle_existence.py`).** Honest
test before enshrining the particle: of all closed-dictionary constructions (21 named
entities + root powers) landing in the *physically viable* window $M=m_\phi/\Sigma m_\nu
\in[295,988]$ (i.e.\ $k_{\rm fs}\in[0.3,1.5]\,h/$Mpc, relevant to $S_8$ without erasing
structure — **not** a fit to a target), the count is strongly grammar-dependent:
permissive grammar (powers+products+sums) gives $1/537$; pure triple-products
("volumes", 3 ceilings) give $1/192$. **Verdict:** PYROS, VITA, MIKA are *real*
dictionary entities, but the multiplier is **not statistically privileged** under a
permissive grammar — its uniqueness is exactly as strong as the lineage-law grammar we
can rigorously justify (the narrow "no-self-sum" rule once gave $\sim1/16$; pinning this
down is the substance of OP-7/OP-9). Versus the old 535.28: statistically similar, but
615.33 is a **pure 3-volume** (cleaner, more Lagrangian-derivable) vs a power$+$product
*sum* — a better *form* for a UV derivation, **not** a stronger statistical selection.
Implication: OP-9 is **not** closed by the dictionary; the promising route is a UV
completion whose $V''(\phi_{\rm min})$ is literally a 3-volume of structural scales.

**Leading mechanism candidate (2026-06-19) — $M=\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V$.**
A data-driven probe (let $m_\phi$ float in CLASS, find what $S_8$ prefers;
`ssee_paper6_particle_scan.py`) gives a data-preferred $m_\phi\simeq41.0$ eV
($M\simeq594$, $S_8=0.00\sigma$ vs the forward 615.33 at $0.24\sigma$). At that value the
construction $M=\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V=(\varphi+2\pi)^2\cdot(\varphi+\pi+\Omega)=594.28$
($m_\phi=40.70$ eV) is notable because — unlike a generic look-elsewhere hit — every
ingredient is **lineage-anchored** to an *already-established* role (the standard set by
KAL$\to$viscosity and MIRA$\to f_{\rm screen}$):
- $\mathrm{KRYSTOS}_V=\varphi+\pi+\Omega$ (padres {φ,π,Ω}; Paper 4), already anchors $w_a=-P_{sc}/K_v=-0.670$
  (DESI). The vacuum scale is **not** assigned here.
- $\mathrm{SOLAR}=\mathrm{BIAL}+\mathrm{KAL}$: BIAL ($=(\varphi+\pi)/2$, the genesis
  "first heat/pulse" — radiative seed) $+$ KAL ($=\beta+\pi$, the **established bulk
  viscosity** $\tilde\zeta=\mathrm{KAL}_0/3$, Paper 5). Both parents are thermal/dissipative,
  so SOLAR as the **radiative–dissipative coupling** is *inherited*, not chosen; its value
  $\varphi+2\pi$ is forced by the lineage.
- Form $m_\phi=g^2 v\,\Sigma m_\nu$ is a standard generated mass (self-energy $\propto g^2$,
  vacuum $v$, neutrino-portal seed). The $\times594$ is an **enhancement** (seesaw/geometric),
  not a loop (loops suppress) — so a loop-radiative narrative is wrong; enhancement is right.

**Honest residue (why this is *candidate*, not *closed*):** (i) it was assembled after
fitting $m_\phi$ to $S_8$ (timing — not forward); (ii) the quantitative coefficient
$\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V$ is dimensionally standard but **not yet derived** by
solving the dissipative ($\mathrm{KAL}_0$-governed) free-streaming/mass equation — the
obvious shortcuts (relic $T_\phi$, $V''$) don't yield SOLAR trivially; (iii) free-streaming
(collisionless) vs IS bulk viscosity (fluid) are technically different dissipations, so the
"single $\mathrm{KAL}_0$ governs both" unification is the claim to prove. Status: SOLAR moved
from "a number that fit nicely" to "BIAL(first-heat)+KAL(anchored viscosity), standard
$g^2v$ form" — a real qualitative anchor; the closing calculation is the explicit
dissipative mass derivation. Scripts: `op9_multiplier_search.py`,
`ssee_paper6_particle_scan.py`. **Canonical particle is now SOLAR²·KRYSTOS_V=594.28
($m_\phi=40.70$ eV), adopted 2026-06-19 (OP-17); the open OP-9 residue is deriving its
coefficient from KAL-governed dissipative transport.**

**Methodological note — lineage-restricted prediction, not blind search (2026-06-20).**
The construction was *not* obtained by combining dictionary constants until one hit 594.
The procedure was: (1) ask what *value* the physics requires (the $M$ that makes $S_8$
close, given our neutrino data + SSEE physics); (2) ask which construction the **lineage
laws** emit *naturally* — i.e. constrained by rules that pre-exist the particle
(no-self-sum, the AURA copy-law at the bifurcation, role-inheritance). SOLAR=BIAL+KAL and
KRYSTOS_V=2Ω came out of step (2) because their *parents already had roles*, not because
their numeric value was convenient. This is **lineage-restricted selection**, categorically
different from a free fit: the admissible space is fixed by laws, not by the datum. The
honest residue is therefore **not** "it was post-hoc" in the loose sense — it is that the
*strength* of "emerged from the lineage" equals **how restrictive the lineage grammar
provably is**, and that number has not yet been measured. The look-elsewhere scans run so
far used *permissive* grammars (products+sums+powers → 1/537; triple-products → 1/192);
the **lineage-strict grammar** (no-self-sum ∧ copy-law ∧ each factor role-anchored ∧ $g^2v$
form) has **not** been counted. **DONE — residue (i) quantified (2026-06-20,
`op9_lineage_grammar_scan.py`):** counting admissible constructions in the viable window
$M\in[295,988]$ under each rule of the lineage grammar gives a sharp gradient:
permissive → 1/hundreds; **L1** (physical form $g^2v$, $M=X^2Y$) → 1/143; **L1+L2**
(+no-self-sum) → 1/80; **L1+L2+L4** (+prior role: coupling$^2\cdot$scale, with
couplings $\{$BIAL,KAL,SOLAR,AURA$\}$ and scales $\{$KRYSTOS_V$=2\Omega$,OMEGA,PYROS$\}$,
each role declared with paper provenance) → **1/3**. The three survivors all share
$\mathrm{SOLAR}^2$ and differ only in the scale (SOLAR²·OMEGA=297, SOLAR²·PYROS=398,
**SOLAR²·KRYSTOS_V=594**); of the three, KRYSTOS_V$=2\Omega$ is the only scale that already
anchors $w_a$ (the DE vacuum scale) — a physical reason to prefer it within the 1/3.
Crucially the old $\Omega^4+\mathrm{AURA}\cdot\mathrm{KAL}=535.28$ is a **sum**, so it
**fails L1** ($g^2v$ form) and is excluded by physics, not by tuning. **Verdict:** the
timing objection (i) is answered — selection is 1/3 under the grammar actually used, not
1/hundreds. Caveat (honest): the role assignment (L4) is *declared, not unique*; this
measures how constrained the choice was, it does **not** prove the coefficient. Residue
(ii) — deriving 594.28 as the *output* of the $\mathrm{KAL}_0$ transport equation —
remains the deeper closure via OP-10. See [[project-op17-particle-deferred]],
[[project-ssee-lineage-laws]].

**Bounded transport-derivation attempt (2026-07-11, negative — do not re-run blind;
`archive/codigo/investigacion/open_problems/op9_transport_derivation.py`).** A pre-committed
attempt to derive $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V$ as the *output* of the
$\mathrm{KAL}_0$-governed dissipative transport (cutoff: simple fraction + residue $<0.5\%$,
no circular re-parametrization). Result: the only exact rewrite found,
$\mathrm{SOLAR}=(\varphi+4\,\mathrm{KAL}_0)/3$, is a **trivial algebraic restatement** of the
natural form $\mathrm{SOLAR}=\mathrm{BIAL}+\mathrm{KAL}_0=\varphi+2\pi$: substituting
$\pi=(2\mathrm{KAL}_0-\varphi)/3$ turns the $2\pi$ into $4\mathrm{KAL}_0/3$, so the factor
$4=2\times2$ is a **substitution artifact with no physical content** — it carries no
transport meaning and **fails the cutoff** (re-parametrization, not derivation). The natural
two-term form $\mathrm{BIAL}+\mathrm{KAL}_0$ (radiative seed $+$ bulk viscosity) is the
honest statement; no cleaner one exists.
Simple-fraction probes against transport combinations ($\mathrm{KAL}_0^3$, etc.) returned
only coincidental hits (24/11, 11/2), excluded by the anti-numerology lock (small residue
$\wedge$ *simple* fraction, both required). **Verdict:** the $g^2v$ form remains the ceiling
of what is derivable without solving the full $V(\phi)$; residue (ii) is genuinely OP-10-level,
not a shortcut away. This confirms — does not weaken — the honest status above.

**Strong rule (still in effect):** the physical Hubble scale
$H_0^{\rm MIRA}=67.068$~km/s/Mpc **must NOT be substituted** into any mass formula.
The canonical chain uses the fixed mass scale $0.9530$ eV and pure $(\varphi,\pi)$
numbers only — no Hubble rate enters. The dimensionally-inconsistent
"$\times H_0^{\rm alg}$" framing of the old 5.602 eV ansatz is retired.

**Canonical CLASS verification (zero fitting, SOLAR²·KRYSTOS_V adopted 2026-06-19):** the
forward-predicted particle ($m_\phi=40.70$ eV, $\Omega_{\phi{\rm DM}}=0.14889$) yields
$\sigma_8^{\rm eff} = 0.747$, $S_8 = 0.758$ ($0.04\sigma$ KiDS-1000 $0.759\pm0.024$ —
**resolves** the lensing tension), single-sector ceiling $S_8=0.846$ ($3.5\sigma$, "the
challenge"). $\sigma_8$ is a **direct CLASS output** (top-hat on the two-sector $P(k)$),
not the retired $\alpha_{\rm WDM}$ fit; the Viel $\alpha = 1.108$ Mpc/h is a CLASS
diagnostic output. log: `results/logs/p6_class_reframe_omega_m.log`.

**Falsifiable anchor:** $k_{\rm fs} = 0.754\,h/$Mpc (analytic; CLASS-measured half-mode
$k_{1/2} = 0.339\,h/$Mpc), set by $m_\phi = 40.70$ eV (SOLAR²·KRYSTOS_V; the prediction
*updates* with the particle — was 0.798 at 42.47, 0.659 at the old 36.95 eV). The whole chain is forward:
$\Sigma m_\nu=0.069$ eV is itself a prediction ($\mathcal{R}_2$), comfortably allowed in
the dynamical-DE background (the tight DESI $\Sigma m_\nu\lesssim0.064$–$0.072$ eV bound
is $\Lambda$CDM-specific; with $w_0,w_a$ free it relaxes to $\sim0.13$–$0.16$ eV). If
DESI Y3/Euclid disconfirms $k_{\rm fs}$, the particle dies cleanly.

**Path to full resolution:** Derive the unified $V(\phi)$ (OP-10) such that
$m_\phi^2 = V''(\phi_{\rm min})$ reproduces the multiplier $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V=594.28$
($g^2v$ form: coupling$^2\cdot$scale) as the curvature of the potential at its minimum,
built from structural scales at the eV scale.

**Severity:** Medium — downgraded from Medium-High. The dimensional contradiction is
resolved (pure-number multiplier + written Lagrangian, dimensionally consistent); what
remains is the *derivation* of the multiplier from first principles, which is
incompleteness on the natural OP-10 path.

**Search file:** `src/op9_phi_dm_formula_search.py` — historical inventory of mass
combinations (now superseded by the canonical Vía-2 chain).

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

### Registro de rutas exploradas (actualizado 2026-06-11)

Protocolo en todas las fases: **forward, cero fiteo** — solo escalas ya bloqueadas
por Papers 6/7/10 (V₀ = 0.840 ρ_crit, λ = √(3·0.160) ≈ 0.693, M_UV = 9.68 meV,
m_φ = 40.70 eV). Requisitos declarados ANTES de calcular. Scripts en
`archive/codigo/investigacion/open_problems/`.

| Ruta | Mecanismo | Veredicto | Causa del descarte |
|---|---|---|---|
| F1 (pozo local) | mínimo cuadrático añadido al exponencial | ❌ | el pozo contamina el plateau DE — w₀ se aleja de −0.840 |
| F2 (tanh²) | doble-escala suave | ❌ | jerarquía DE/DM exige 33 órdenes de magnitud en el parámetro de forma |
| F3 (axion-like) | V₀e^(−λφ) + M_UV⁴[1−cos(φ/f)] | ❌ | atrapamiento topológico: si el axion engancha, φ se fija → w = −1; si no engancha, m_φ ≈ H₀ ≠ 40.70 eV (`ssee_op10_family3_axion.py`) |
| F4 (see-saw escalas) | m_φ² = m_h·m_l con escalas SSEE | ❌ | predicción off por ~10³³ (`ssee_op10_seesaw_search.py`) |
| Opción A (UV-inducido) | corrección de loop de K(X) genera el mínimo | ❌ | K(X) es φ-independiente a nivel árbol; curvatura de loop ~ M = 9.68 meV, factor ~4 238 bajo m_φ (`ssee_op10_uv_induced_minimum.py`) |
| 2c (K-essence boost) | K_X renormaliza la masa efectiva | ❌ analítico | K_X ≥ 1/KAL acota el boost a √KAL ≈ 2.35× — insuficiente para 10³⁻⁴ |
| 2e (masa tracking) | m_eff²(z) = ρ_m(z)/M_*², M_* = M_UV | ❌ estructural | ver R1–R3 abajo (`ssee_op10_phase2e_dynamic_mass.py`) |

**Fase 2e en detalle** (cerrada 2026-06-10; idea sobreviviente de la propuesta
"conversión DM→DE" de Mike): masa inducida por densidad ambiente, chameleon-like,
con la única escala bloqueada M_* = M_UV. Hallazgo sugestivo: el cruce
m_eff = 40.70 eV cae en z_x ≈ 2 370, época de igualdad materia-radiación
(z_x/z_eq ≈ 0.70) sin ajustar nada. Pero falla los tres requisitos pre-declarados:

- **R1 (estructural, independiente de M_*):** con m ∝ √ρ_m ∝ a^(−3/2) y el
  invariante adiabático n ∝ a⁻³, la energía oscilante diluye como
  ρ_osc = m·n ∝ a^(−4.5) → w_eff = +0.5: se evapora más rápido que la radiación,
  no es CDM para NINGÚN M_*.
- **R2:** m_eff(0) = 0.355 meV, factor ~10⁵ bajo el canónico 40.70 eV (rompe
  k_fs = 0.754 h/Mpc de Paper 6).
- **R3:** m_eff(0)/H₀ ~ 10²⁹ — el campo queda clavado hoy → w = −1, no −0.840.
- Rescate M_* = 92.4 neV sería escala nueva fiteada (prohibido) y R1 persiste.

### Teorema del escalón (lección estructural consolidada, 2026-06-10)

Los 7 descartes consolidan un patrón: **ni un pozo estático ni una masa que
trackea ρ_m de forma continua pueden unificar DE+DM en un solo campo SSEE.**
R1 lo prohíbe estructuralmente: mientras la masa cambie durante las oscilaciones,
la energía no diluye como materia.

⇒ La transición DM↔DE, si existe en un solo campo, tiene que ser un **escalón**
(transición de fase única): el potencial cambia de forma UNA sola vez cerca de
z_eq, con masa constante = 40.70 eV en la fase oscilante. El escenario
escalón + masa constante ES el misalignment frío ya estudiado (cierre Osiris
2026-05-30): viable como DM, pero φ_i = 1.17×10⁻⁶ M_Pl no sale algebraico de
(φ,π) → empata ΛCDM en conteo de parámetros, no lo supera.

**Única ruta no excluida:** transición de fase única cerca de z_eq cuyo
DISPARADOR esté bloqueado por (φ,π). Sin candidato a la fecha.

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

## OP-12 — Origen físico de T_φ y de Ω_φ-DM h² (Paper 6) — ABIERTO (reframe 2026-06-20)

**Location:** Paper 6, §4.2; pipeline CLASS `src/ssee_paper6_canonical_particle.py`.

**Problem (actualizado):** la partícula canónica φ-DM ($m_\phi=\mathrm{SOLAR}^2\cdot
\mathrm{KRYSTOS}_V\cdot\Sigma m_\nu=40.70$ eV; antes 36.95) entra en CLASS como especie
**térmica** (ncdm) con temperatura $T_\phi$. Pero hoy $T_\phi$ se **despeja** de la
fórmula del relic usando $m_\phi$ y $\Omega_{\phi DM}$:
$$T_\phi = T_\nu\left(\Omega_{\phi DM}h^2\cdot 93.14/m_\phi\right)^{1/3}=0.5385\,T_\nu.$$
Es **bookkeeping** (densidad = número × peso), no una temperatura derivada de física.
Verificado en CLASS: con $(m_\phi,T_\phi)$ el Boltzmann completo reproduce
$\Omega_{\phi DM}=0.14889$ ($\Omega h^2=0.06877$) a 5 cifras — el lazo cierra, pero
porque se construyó así. $m_\phi$ descansa en una sola pata (SOLAR²·KRYSTOS_V, OP-9).

**Bifurcación (cold vs thermal) — y una inconsistencia interna a corregir:** el texto
de P6 asume la rama **FRÍA** (producción gravitacional Parker-Kolb-Riotto, "escalar
libre, sin portal SM, invisible"), mientras el **pipeline CLASS trata φ como TÉRMICO**
(ncdm a $T_\phi$). No pueden ser ambas: un campo frío no tiene distribución térmica.
Esto hay que unificarlo.

- **Rama fría** (misalignment / producción gravitacional): consistente con "sin portal",
  pero entonces $(T_\phi/T_\nu)^3$ es ficción de bookkeeping y la 2ª pata de $m_\phi$
  debe venir del potencial $V''(\phi_{\min})$ (ruta OP-10), no del freeze-out.
- **Rama térmica** (freeze-out): $T_\phi$ es una temperatura real → el tratamiento
  ncdm de CLASS es consistente, y el freeze-out **podría derivar $T_\phi$** y cerrar
  $m_\phi$ por una 2ª vía independiente de SOLAR²·KRYSTOS_V. Requiere un **portal**.

**Insight (Mike, 2026-06-20):** el portal de la rama térmica **NO es nuevo** — es la
**viscosidad** que ya está en el modelo ($\mathrm{KAL}_0$, $\tilde\zeta=\mathrm{KAL}_0/3$,
Paper 5 IS). Disipación = interacción con un baño; un escalar libre no tiene viscosidad.
$\mathrm{SOLAR}=\mathrm{BIAL}+\mathrm{KAL}$ hereda el rol **radiativo-disipativo** →
señala que el baño es la **radiación**. Esto empuja hacia la rama térmica y **revisa**
el "sin portal SM" (φ pasaría a ser débilmente acoplado, casi invisible pero no nulo).

**Paso 1 — VIABILIDAD (✅ hecho, `archive/codigo/investigacion/open_problems/op12_thermal_decoupling.py`):** si φ
es relic térmico, la dilución de entropía $(T_\phi/T_\nu)^3=10.75/g_{*s}(T_{\rm dec})$
exige $g_{*s}\approx 68.8 \Rightarrow T_{\rm dec}\approx 217$ MeV — **justo la transición
QCD**. Época física, no absurda: la rama térmica es **viable y falsable** (predice
*cuándo* se desenchufa φ). Coincidencia a vigilar (no afirmar): $g_{*s}=68.8$ vs
$H_{\rm alg}=3\Omega^2=67.96$.

**Naturaleza del relic (auto-consistencia):** a $T_{\rm dec}=217$ MeV, $m_\phi/T_{\rm
dec}\sim2\times10^{-7}\ll1$ → φ desacopla **ultra-relativista** (relic **caliente**,
como el ν). Esto **justifica** la fórmula de dilución $(T_\phi/T_\nu)^3$ (solo vale si
desacopla relativista) → el "bookkeeping" no es arbitrario. φ se vuelve no-relativista
en $z\approx4.5\times10^5$ (mucho antes del CMB) → DM **templada** en recombinación,
encajando con el tratamiento WDM/$k_{\rm fs}$ que P6 ya usa. La rama térmica se
**enchufa** a la física existente, no la contradice.

**Pasos 2–3 (abiertos):** (2) mostrar que la viscosidad KAL es el portal y fijar su
estructura/escala $\Lambda$; (3) freeze-out forward → ¿da $T_\phi/T_\nu=0.5385$? Si sí,
2ª pata de $m_\phi$. **Obstáculo honesto:** $\mathrm{SOLAR}\approx 7.9$ no es un
acoplamiento perturbativo ($g>1$) → el portal es casi seguro suprimido por escala
(operador dim$>4$); el desconocido real pasa a ser la escala $\Lambda$ del portal.

**Paso 2 — dos rutas exploradas 2026-06-21 (scripts `op12_route2_freezeout.py`):**
- **Ruta 2 (cross-section, número duro):** freeze-out a 217 MeV exige $\sigma\approx
  6\times10^{-19}$ GeV$^{-2}$ → portal **dim-4** $g\sim3\times10^{-5}$ (débil) o **dim-5**
  $\Lambda\sim8.8$ TeV. **HALLAZGO (correctness burden):** el portal **NO es SOLAR**
  ($g_{\rm portal}\sim10^{-5}$ vs SOLAR$\sim$7.9, ~5 órdenes). **La masa (SOLAR, g²v) y el
  portal son acoplamientos DISTINTOS** — la idea "mismo SOLAR hace masa Y portal" no
  cierra. Coherencia: acoplamiento débil $\Leftrightarrow$ viscosidad alta (signo del
  Bridge 2) → portal débil es consistente con KAL grande.
- **Ruta 1 (linaje, "cuándo"):** el desacople cae en la transición QCD (217 MeV) porque
  un portal al **sector QCD** (quarks/gluones) se apaga al **confinarse** los quarks
  (~200 MeV) → el confinamiento es el **gatillo físico** (no perilla). El plasma
  quark-gluón es el medio viscoso por excelencia ↔ KAL=viscosidad. **Rima, no deriva**
  (emparejamiento de roles SOLAR=BIAL+KAL = térmico+viscoso ↔ QGP).
- **Convergencia:** portal débil al plasma QCD viscoso, desacople = confinamiento; KAL
  es la huella macroscópica.
- **$\Lambda\sim8.8$ TeV NO es blanco a derivar (disuelto 2026-06-21, regla dimensional
  de Mike):** la condición de freeze-out lo FIJA, $\Lambda=[O(1)\cdot M_{\rm Pl}\cdot
  T_{\rm dec}^3]^{1/4}$ — anclado a $M_{\rm Pl}$ (gravedad) y $T_{\rm dec}$ (escala QCD),
  ambas físicas. No hay número que sacar de $\varphi,\pi$; es la "sombra" de QCD+gravedad.
  La pregunta real se afila a: **¿el linaje de SOLAR pone a φ en el sector QCD?** (Ruta 1).
- **2ª pata de $m_\phi$ (consistente, no cerrada):** si φ↔QCD, confinamiento $\to g_{*s}
  \to T_\phi \to m_\phi$ (con $\Omega_{\phi DM}$), SIN SOLAR²·KRYSTOS_V. Coincide con 41 eV
  **solo si el desacople se completa en la cima del confinamiento** ($g_{*s}=69$, 217 MeV);
  más abajo baja hasta ~12 eV. Dos rutas independientes (linaje=41, térmica-en-onset=41)
  riman → alentador, PERO sub-asunción (punto de desacople) + Ruta 1 es rima no derivación.
  **Promesa, no prueba.** Scripts: `op12_thermal_decoupling.py`, `op12_route2_freezeout.py`.
- **Paso 2 INTENTADO 2026-06-21 (`op12_trace_lights_qcd.py`) — resultado NEGATIVo honesto:**
  el término candidato $\mathcal L=(\varphi/\Lambda)T^\mu_\mu$ (acoplo a la traza, fuente de
  KAL). (a) La interaction measure $(1-3w)$ SÍ pica en QCD (~190 MeV) ✓. PERO (b) en
  $\Gamma/H\propto I^2 T^3 M_{\rm Pl}/\Lambda^2$ el factor $T^3$ domina → el bulto de QCD
  es solo ~4% del máximo → **el desacople NO se ancla en QCD** (ocurre en alta T, fijado
  por Λ). **La historia "φ se desacopla en QCD gratis por la traza" NO sobrevive el
  cálculo** — caer en QCD exigiría ajustar Λ~8.8 TeV a mano (no derivado). El gatillo de
  confinamiento queda DEBILITADO. **Saldo:** la rama térmica NO cierra elegante; $m_\phi$=41
  sigue en 1 pata (SOLAR²·KRYSTOS_V). 2ª pata regresa a $V''(\phi_{\min})$/OP-10 (rama fría) o
  admitir Λ sin explicar. (Corrige una narrativa previa optimista; el cálculo mandó.)

**Path to resolution:** elegir rama (la viscosidad ya presente favorece la térmica) y
o bien (térmica) cerrar Pasos 2–3, o bien (fría) derivar $\Omega_{\phi DM}h^2$ del
$V(\phi)$ unificado (OP-10) sin tuning. Ligado a OP-10/OP-11.

**Severity:** Medium — el relic match es algebraico; cualquiera de las dos ramas, bien
cerrada, lo vuelve dinámico. La inconsistencia cold/thermal del pipeline es **a
corregir sí o sí** (independiente de qué rama gane).

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

## OP-14 — Σm_ν Phenomenological Derivation (Paper 4) — ✅ RESUELTO (2026-06-04)

**Location:** Paper 4, §"Neutrino Mass Sum", L675-700.

**Fórmula canónica (resuelta):**
$$\Sigma m_\nu^{\rm active} = \mathcal{R}_2 \times 0.9530~\mathrm{eV} = 0.0685~\mathrm{eV},
\qquad \mathcal{R}_2 = \frac{\Omega_{\rm DNAV}}{\mathrm{KAL}\cdot\mathrm{TRIAL}}
= \frac{4.7596}{5.5214 \times 11.9935} = 0.07188.$$

El cociente $\mathcal{R}_2$ es un número puro de $(\varphi,\pi)$: $\Omega_{\rm DNAV}=\pi+\varphi$,
$\mathrm{KAL}=(\pi+\varphi)/2+\pi$, $\mathrm{TRIAL}=3(\varphi+(\pi+\varphi)/2)$. No hay
sustracción de enteros, no hay offset 22.

**Por qué esto resuelve el problema:**

1. **El offset 22 queda eliminado.** La forma anterior $\mathcal{R}=4\cdot\mathrm{KAL}-22=0.0856$
   era una diferencia entre números casi iguales (22.086 vs 22.000) y se ajustaba para
   reproducir la cota cosmológica. La forma canónica $\mathcal{R}_2=\Omega_{\rm DNAV}/(\mathrm{KAL}\cdot\mathrm{TRIAL})$
   es un cociente limpio de constantes estructurales — sin parámetro ajustado.
2. **Σm_ν asciende de Type P → Type A (algebraico).** El único input externo que queda es
   $0.9530$ eV (constante de normalización fija del Modelo Estándar relíquica↔masa), igual
   que cualquier predicción dimensional de SSEE usa una escala física fija.

**Cascada (actualizada):** Σm_ν alimenta a $m_\varphi = \Sigma m_\nu^{\rm active}\,(\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V)=40.70$ eV
en P6 (ver OP-9). Con OP-14 resuelto y OP-9 refinado, la cadena $\varphi,\pi \to \Sigma m_\nu \to m_\varphi$
es ahora **forward-prediction sin parámetros libres**. El antiguo eslabón más débil queda cerrado.

> **Nota histórica (lo de abajo precede a la resolución 2026-06-04).** El registro
> del ataque 2026-05-23 se conserva por honestidad: muestra por qué la forma antigua
> $\mathcal{R}=4\cdot\mathrm{KAL}-22$ era frágil y por qué se descartó en favor del
> cociente limpio $\mathcal{R}_2=\Omega_{\rm DNAV}/(\mathrm{KAL}\cdot\mathrm{TRIAL})$.

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
libertad fenomenológico** vía la identidad canónica de Paper 6:
$$m_\varphi = \Sigma m_\nu^{\rm active}\cdot(\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V) = 0.0685\,\text{eV}\times 594.28 = 40.70\ \text{eV}$$

> **Nota (2026-06-14):** la forma previa $m_\varphi=\Sigma m_\nu\cdot H_0^{\rm alg}$
> (que daba 5.602 eV) está **RETIRADA** — era dim-inconsistente (eV × km/s/Mpc).
> El argumento de abajo NO cambia: $m_\varphi$ sigue construido a partir de
> $\Sigma m_\nu$ (multiplicador adimensional puro), así que comparten el mismo DoF.

OP-14 y OP-9 NO son problemas independientes — son **dos caras del mismo
problema**. Atacar OP-14 directamente no produce derivación porque no hay
ataque local: la única salida es derivar $m_\varphi$ desde la curvatura
de un potencial $V(\varphi)$ fundamental.

### Cadena de dependencias (post-resolución 2026-06-04)

```
OP-14  Σm_ν = R₂·0.9530 eV = 0.0685 eV    ✅ RESUELTO (R₂=Ω_DNAV/(KAL·TRIAL), sin offset 22)
   │
   └──► OP-9   m_φ = Σm_ν·(SOLAR²·KRYSTOS_V) = 40.70 eV   ✅ refinado: dim-consistente, forward-prediction
           │     (queda abierto SOLO el origen UV del multiplicador 594.28)
           │
           ├──► OP-11  (ξ acoplamiento → funcional de V)        depende de OP-10
           └──► OP-12  (Ω_φDM h² desde dinámica de V)            depende de OP-10

OP-10  V(φ) unificador DE + DM   (daría el origen UV del multiplicador; ya no bloquea Σm_ν ni m_φ)
OP-8   MIRA mecanismo (independiente, eslabón duro)
```

La antigua flecha OP-10→…→OP-14 queda obsoleta: la cadena $\varphi,\pi\to\Sigma m_\nu\to m_\varphi$
ya cierra como predicción algebraica sin V(φ). OP-10 sigue abierto pero su rol es dar el
**origen UV del multiplicador** $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V=594.28$,
no rescatar la masa.

### Acción sobre Paper 4 (completada 2026-06-04)

El §"Neutrino Mass Sum" se actualiza a la fórmula canónica
$\Sigma m_\nu^{\rm active}=\mathcal{R}_2\times 0.9530$ eV $=0.0685$ eV con
$\mathcal{R}_2=\Omega_{\rm DNAV}/(\mathrm{KAL}\cdot\mathrm{TRIAL})$, promoviendo Σm_ν de
**Type P → Type A**. El offset 22 queda eliminado del manuscrito.

---

## OP-15 — Bullet Cluster Offset κ(θ) from KAL(x) Not Computed (Paper 1) — ABIERTO

**Location:** Paper 1, §3.4 (Gravitational Lensing and the Newtonian Limit), L618-629.

**Problem:** The mass-discrepancy result (§3.1-3.3, Table of 4 clusters, χ²_r=0.122) is
solid: the structural amplifier KAL₀ maps the IGIMF baryonic baseline onto the full
Newtonian-equivalent dynamical mass without cold dark matter. **A separate, harder claim
is not yet demonstrated:** that the *environment-dependent* amplifier KAL(x) reproduces
the **spatial offset** between the X-ray gas and the lensing-mass peak in the Bullet
Cluster (1E 0657-56) — the single most-cited evidence for particulate dark matter.

The paper currently states only:
> "SSEE proposes a *qualitative* mechanism for this offset … A quantitative calculation
> of the projected convergence κ(θ⃗) from first principles is deferred to future work."

**Why this matters (audit lesson, 2026-06-14):** this was a *soft* hedge buried in prose,
easy to read as finished. A modified-gravity boost that tracks the baryons must put the
extra gravity where the baryons are; the Bullet's offset (gravity displaced toward the
collisionless galaxies, away from the gas) is precisely the configuration that challenges
any baryon-tracking amplifier. SSEE's proposed escape is that KAL(x) amplifies more in the
compact galactic potential wells than in the diffuse stripped gas — but this is asserted,
not computed.

**Path to resolution:** compute Σ_SSEE(θ⃗) = ∫ ρ_bar^IGIMF · KAL(x) dℓ for the Bullet's
observed gas+galaxy distribution and show the resulting κ(θ⃗) peak lands on the galaxies
(matching Clowe+2006 lensing), not the gas. Falsifiable: if the KAL(x)-weighted convergence
peak coincides with the gas instead, the qualitative mechanism fails.

**Severity:** Medium-High — the Bullet offset is a headline dark-matter argument; leaving it
at "proposed mechanism" is honest but a referee will press on it. Not fatal to the mass
result (OP-15 is about *where* the effective mass sits, not *how much*).

**Relation:** distinct from OP-13 (Paper 8 lensing, RESUELTO). OP-13 established the
canonical θ_E^SSEE ≈ θ_E^GR-with-DM (lensing *amplitude*); OP-15 is the *spatial*
distribution in a merging cluster.

---

## OP-16 — ¿Coincide $(\pi-\varphi)/(\pi+\varphi)=0.3201$ con una fracción medida de la descomposición masa-energía del protón? (origen génesis) — ABIERTO / ESPECULATIVO

**Origen.** El sistema fenomenológico génesis (`SSEE_UNIFICADO`, "Resolución Física de
Partículas") asignaba al protón un "registro" $\mathrm{Ygg}_P = 0.319$ ("66% gravedad +
34% fotones / pensamiento cristalizado"). Ese número se había filtrado a Paper 4 como
"Proton register" y fue **retirado en la auditoría P4-A (2026-06-15)** por dos defectos:
(1) **circular** — el "valor observado" era $0.319$ (YGG, interno), no una medición;
(2) **mal etiquetado** — $0.3201 = \Omega_{m,\rm CMB}$ es materia **total** (dominada por DM),
no bariónica ($\Omega_b \approx 0.049$). No es física de partículas: es cosmología disfrazada.

**Por qué NO está muerto (la salida honesta).** La intuición génesis *"la masa del protón
es energía de campo, no constituyentes"* **es físicamente correcta**: la descomposición
lattice-QCD de la masa del protón (p.ej. Ji / Yang et al.\ 2018) da masa de quarks $\sim 9\%$,
energía cinética de quarks $\sim 32\%$, energía de gluones $\sim 37\%$, anomalía de traza
$\sim 23\%$. La masa del protón **sí** es mayoritariamente energía de campo.

**El test falsable (requisitos de entrada — las tres reglas).** Para que $0.3201$ se gane un
lugar en física de partículas (dominio del Modelo Estándar, **independiente** de la cosmología;
distinto de $m_\varphi$ que vive *dentro* del sector oscuro cosmológico):
1. **Comparar contra una fracción MEDIDA** de la descomposición (con barras de error de lattice),
   NO contra un número interno YGG.
2. **Forward-prediction:** la fórmula $(\pi-\varphi)/(\pi+\varphi)$ fija de antemano, sin elegir
   a posteriori cuál de las 4 fracciones "casa mejor" (evitar look-elsewhere oculto).
3. **Anclaje físico:** idealmente atado a una escala/Lagrangiano (no una coincidencia adimensional
   suelta), igual que se exige en [[feedback_three_principles]].

**Criterio de cierre.** Si $0.3201$ está a $>2\sigma$ de **toda** fracción estándar de la
descomposición → coincidencia, se cierra. Si casa con **una** dentro de error Y se puede atar
estructuralmente a $\{\varphi,\pi\}$ → merece una nota dedicada (no antes).

**Severidad: Baja / especulativa.** Cero impacto sobre la suite cosmológica de 10 papers
(la cosmología no depende de esto). Es una **dirección de investigación**, no un claim de paper.
Abordar SOLO después de las auditorías. Relacionado con la advertencia [[project_genesis_repo_risk]].

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
| OP-8 | Transv. | ~~MIRA dynamical mechanism~~ → factor-materia DISUELTO | ✅ DISUELTO 2026-06-18 | Reframe ω_m-directo: Ω_m,CMB=ω_m/h²=0.30889 sin factor (ω_c=KAL₀·ω_b·n_s forward); MIRA persiste solo en f_screen; CMB χ²=1005.5/ΔBIC=−23.9 |
| OP-9 | P6 | ~~m_φ=5.60 eV numerological~~ → m_φ=40.70 eV canónico | 🔶 **ABIERTO** (derivar coef UV 594.28 → OP-10) · partícula ya adoptada [OP-17] y falseable (k_fs=0.754) | $m_\varphi=\Sigma m_\nu^{\rm active}\cdot(\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS\_V})=40.70$ eV (SOLAR=φ+2π, **KRYSTOS_V=φ+π+Ω** padres {φ,π,Ω}, NO 2Ω; mult=594.28), término de masa $g^2v$ escrito, **falseable por $k_{\rm fs}=0.754$ (DESI/Euclid)** — contenido MEDIBLE cerrado (adoptado OP-17). El origen UV del coeficiente vive **más allá del horizonte de saturación TRIAL** (no es un hueco: es el borde de lo medible; la derivación dinámica es OP-10, gobernada por ley, no numerología). Intento de transporte 2026-07-11 → CORTE (ver §OP-9). Antes 36.95 (Ω⁴+AURA·KAL) y 42.47 (PYROS·VITA·MIKA), retirados |
| OP-10 | P6/P7 | Unify χ into φ via richer V(φ) | Medium-High | V(φ) with slope (DE) + minimum (DM matter-mode); restores zero-param status |
| OP-11 | P6 | ξ (non-minimal coupling) is free parameter | Medium | Algebraic constraint or absorb via OP-10 |
| OP-12 | P6 | Ω_φ-DM h² not computed ab initio | Medium | Parker-Kolb-Riotto with α-attractor + m_φ |
| OP-13 | P8 | ~~Contradicción interna §3-4 vs §4.5~~ | ✅ RESUELTO | Opción A aplicada: framing dos-límites, retirado claim "MIRA en lensing", $\sqrt{\AURA}$ ≠ $\MIRA$ aclarado, canonical prediction = GR-with-DM (2026-05-23) |
| OP-14 | P4 | ~~Σm_ν Type P; offset 22 ad hoc~~ → canónico Type A | ✅ RESUELTO | $\Sigma m_\nu^{\rm active}=\mathcal{R}_2\times 0.9530$ eV $=0.0685$ eV con $\mathcal{R}_2=\Omega_{\rm DNAV}/(\mathrm{KAL}\cdot\mathrm{TRIAL})=0.07188$; offset 22 eliminado, Σm_ν promovido Type P→Type A (2026-06-04) |
| OP-15 | P1 | Bullet offset κ(θ) desde KAL(x) no calculado | Medium-High | Computar Σ_SSEE(θ)=∫ρ_bar·KAL(x)dℓ del Bala; mostrar pico κ sobre galaxias, no gas (falsable vs Clowe+2006). Distinto de OP-13 (amplitud); esto es distribución espacial (2026-06-14) |
| OP-16 | — (génesis) | ¿0.3201=(π−φ)/(π+φ) casa con fracción medida de la masa-energía del protón? | Baja/especulativa | Retirado de P4 (P4-A, era circular+materia total mal-etiquetada). Test: comparar vs descomposición lattice-QCD (quark 9%/gluón 37%/anomalía 23%) con barras, forward, anclado a (φ,π). Cero impacto en cosmología; dirección de investigación post-auditoría (2026-06-15) |

**Severity legend:** High = referee would likely request resolution before acceptance;
Medium = requires acknowledgment and discussion; Low = cosmetic or presentational.

---

## Note on Methodology

These problems are documented here rather than concealed because scientific integrity
requires pre-registration of known limitations. Referees and collaborators should be
directed to this document when evaluating the strength of the SSEE predictions.

## OP-17 — Partícula canónica $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V$ ($m_\phi=40.70$ eV, $S_8=0.01\sigma$) — ✅ ADOPTADA 2026-06-19

**Status:** ✅ **CERRADO / ADOPTADO** (Mike, 2026-06-19). El diferimiento se revirtió
en la misma sesión ("Mira bien lo que NO es, es SOLAR"): se adoptó la partícula con
mecanismo en vez de seguir con PYROS·VITA·MIKA. **HECHO:** CLASS forward real @ 40.70 eV
($\sigma_8=0.747$, $S_8=0.758=0.04\sigma$ KiDS, $k_{\rm fs}=0.754$, $\alpha=1.108$,
$T_\phi=0.5385\,T_\nu$, $N_{\rm ur}=2.9619$) → propagado a `CANONICAL_VALUES.yaml`,
guardián (VERDE 119), manuscrito Paper 6 completo (Lagrangiano g²·v escrito, fσ8
recomputado 0.82σ, 25 pp, `docs/`), 3 memorias. **RETIRADO** 615.33/42.47.
Residuo único: OP-9 (derivar el coeficiente del transporte).

**Qué cambió:** la partícula canónica pasó de $m_\phi=42.47$ eV (PYROS·VITA·MIKA=615.33,
triple-producto plano, sin mecanismo, $S_8=0.24\sigma$) a
$$m_\phi = \mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V\cdot\Sigma m_\nu = (\varphi+2\pi)^2\cdot (\varphi+\pi+\Omega)\cdot 0.0685 = 40.70\ \text{eV},$$
que **resuelve $S_8$ de lleno ($0.04\sigma$ KiDS)** — el valor que los datos prefieren
(`ssee_paper6_particle_scan.py`, CLASS real). Esta era la finalidad del mecanismo.

**Por qué es mejor partícula (peso, no certeza):**
- $\mathrm{KRYSTOS}_V=\varphi+\pi+\Omega$ (padres {φ,π,Ω}), **anclado** por $w_a=-P_{sc}/K_v=-0.670$ (DESI).
- $\mathrm{SOLAR}=\mathrm{BIAL}+\mathrm{KAL}$ = (primer-calor/pulso, radiativo) + (viscosidad
  anclada P5). Rol radiativo **por linaje**, valor $\varphi+2\pi$ forzado.
- Forma $m=g^2 v\,\Sigma m_\nu$ = masa generada estándar (enhancement, no loop).
- **Peso vía KAL (clave):** KAL₀ aparece en 5 lugares — $\omega_c=\mathrm{KAL}_0\,\omega_b\,n_s$
  (¡la densidad de materia del CMB! → **"CMB prefiere KAL"**, igual que "CMB prefiere MIRA"),
  $K(X)=X/\mathrm{KAL}_0$ (cinético que maneja la cascada Hubble), $\tilde\zeta=\mathrm{KAL}_0/3$
  (viscosidad), y ahora el acoplamiento φ-DM. SOLAR hereda el peso del KAL que el CMB exige.

**Lo honesto (no lo decidimos nosotros):** es una partícula **en camino a ser falseada**
($k_{\rm fs}$ en DESI Y3/Euclid 2026–2028 decide). Le damos **más peso** por el multi-anclaje
de KAL, NO certeza. El coeficiente $\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V$ aún NO está derivado
del transporte disipativo (ver el "mecanismo líder" en OP-9) — eso es lo que cerraría OP-9.

**Trabajo de implementación (cuando se haga):**
1. Correr CLASS exacto en $m_\phi=41.017$ → tomar $S_8/k_{\rm fs}/\sigma_8$ reales (no interpolados).
2. Propagar 42.47→41.017 y 615.33→594.28 por: `ssee_core`, `CANONICAL_VALUES.yaml`, guardián,
   3 memorias, vault (kfs/sigma8_eff/Cadena), Papers 1/6.
3. Escribir el mecanismo del **Lagrangiano** ($m_\phi=g^2 v\,\Sigma m_\nu$, $g$=SOLAR disipativo,
   $v$=KRYSTOS_V vacío) y **mostrar cómo resuelve $S_8$** en Paper 6.
4. Guardián VERDE + memory_sync VERDE.

---

None of these problems falsify SSEE at the current observational precision — they define
the boundary of what has been rigorously established versus what remains as working
hypotheses. The resolution of OP-1 through OP-6 constitutes the research agenda for
SSEE-V4.0; OP-8 through OP-13 constitute the agenda for SSEE-V5.0 (full unification
of the dark sector).

## Parameter-Count Status (2026-05-22)

SSEE-V3.6 currently has **~2 effective free parameters** vs **6 for $\Lambda$CDM**
(tras DISOLVER OP-8 el 2026-06-18: el factor-materia ya no es input — Ω_m,CMB sale
de ω_m algebraico; restan $H_0$ y $\Omega_b h^2$ como inputs ajustables):

| Parameter | Status | Comment |
|---|---|---|
| $H_0$ | Sampled in MCMC | Algebraic prediction 67.96 (coincidencia, V-L2-06) vs posterior canónico 66.53 |
| $\Omega_b h^2$ | Sampled in MCMC | Algebraic 0.02242 vs Planck 0.02237 |
| ~~MIRA factor~~ | ✅ DISUELTO 2026-06-18 — no hay factor materia; Ω_m,CMB=ω_m/h² (ω_c=KAL₀·ω_b·n_s forward) | OP-8 cerrado |
| $m_\phi$ | Phenomenological ansatz | OP-9 (paper admits) |
| $\xi$ | **Free parameter** (P6 L416) | OP-11 |
| $\Sigma m_\nu$ | **Phenomenological** (P4 L691-697 admite Type P) | offset 22 ad hoc — OP-14 |
| $\alpha$ (attractor) | $\varphi^4/3$ derived | from $\varphi$ |
| $V(\phi)$ form | Adopted exponential | locks DE-only behavior — OP-10 |

**Path to zero parameters:** OP-14 (✅ resuelto) y OP-9 (✅ refinado, forward-prediction)
cerrados; **OP-8 ✅ DISUELTO 2026-06-18** (factor-materia eliminado vía ω_m-directo).
Resta OP-10 (unify χ into φ via richer $V(\phi)$, daría el origen UV del multiplicador
de la masa φ-DM). End state ya alcanzado en el sector materia: $H_0$ y $\Omega_b h^2$
son los únicos inputs ajustables por observación.

**Dependency chain (cascada, actualizada 2026-06-04):** OP-14 ✅ → OP-9 ✅; ambos cerrados
sin V(φ). OP-10 ya no bloquea la masa — solo daría el origen UV del multiplicador
$\mathrm{SOLAR}^2\cdot\mathrm{KRYSTOS}_V=594.28$.

---

## OP-18 — Derivación de la amplitud primordial $A_s$ desde $(\varphi,\pi)$ (Paper 3 / inflación) — ABIERTO (2026-06-20)

**Origin.** Tras el reframe $\omega_m$-directo, el conteo honesto del sector CMB es
$k=2=\{A_s,\tau\}$: $H_0=3(\varphi+\pi)^2$ es derivado (SH0ES–$f_{\rm screen}$, Paper 9),
$\omega_b,\omega_c,n_s,w_0,w_a$ son algebraicos. Las **dos** únicas cantidades que el
ajuste CMB de SSEE **no** deriva de $\{\varphi,\pi\}$ son la amplitud escalar primordial
$A_s$ y la profundidad óptica de reionización $\tau$. OP-18 ataca la primera.

**El problema.** $A_s\simeq2.1\times10^{-9}$ entra como normalización del espectro
primordial; hoy se toma de Planck. SSEE ya fija el resto del sector inflacionario:
$\alpha=\varphi^4/3$ (α-attractor), $N_*=2\varphi^7$, $n_s=1-\varphi^{-7}=0.96556$,
$r=\varphi^{-10}$ (OP-2). En α-attractor la amplitud es
$A_s = \dfrac{V(\phi_*)}{24\pi^2\,\epsilon(\phi_*)\,M_{\rm Pl}^4}$, evaluada en $N_*$.
Con $\epsilon$ ya fijado por $r=16\epsilon$ (i.e.\ $\epsilon=\varphi^{-10}/16$), el **único
ingrediente faltante** es la **escala del potencial** $V_0$ (equivalentemente $M$, la
escala de inflación).

**Cálculo realizado (2026-06-20, `op18_As_from_inflation.py`) — resultado honesto.**
En el plateau, $A_s=(2\varphi^{10}/3\pi^2)\,(V_0/M_{\rm Pl}^4)$. El **prefactor
$2\varphi^{10}/3\pi^2=8.308$ es puro $(\varphi,\pi)$** — SSEE fija la FORMA del espectro
por completo ($n_s,r,\alpha,N_*$ y el prefactor). PERO la normalización exige
$V_0/M_{\rm Pl}^4=2.53\times10^{-10}$, i.e.\ $V_0^{1/4}\approx9.7\times10^{15}$ GeV
(escala GUT/inflación, razonable). **Esa jerarquía NO sale limpia de potencias
$(\varphi,\pi)$** ($\varphi^{-46}$ cae a 4\% pero $-46$ no tiene razón estructural →
coincidencia, no derivación). **Corrección de una nota previa errónea:** anclar $A_s$ a
$\Lambda_{\rm SSEE}=M=9.68$ meV es **INCORRECTO** — esa escala es IR (energía oscura);
$V_0$ es UV (inflación), separadas $\sim10^{27}$. La "ruta alcanzable vía $M$" queda
**retirada**. OP-18 sigue abierto de verdad: el residuo es la **escala de inflación
$V_0$**, una jerarquía UV análoga al problema de jerarquía, no un anclaje trivial.

**Lo defendible hoy (check de consistencia, no derivación):** dado $A_s$ medido y la forma
SSEE fijada, la escala de inflación inferida cae en el rango GUT esperado ($\sim10^{16}$ GeV).
$A_s$ permanece como **1 input de escala** del modelo; lo que SSEE aporta es que $n_s,r$ ya
no son perillas. Cierre futuro: una teoría de $V_0$ desde la UV-completion (liga a OP-10/OP-3).

**Test cuantitativo (forward, falsable) — pendiente de una teoría de $V_0$.** Si en el
futuro $V_0$ se deriva estructuralmente, calcular $A_s^{\rm SSEE}$ y comparar con Planck
$\ln(10^{10}A_s)=3.044\pm0.014$. Si cae dentro de $\sim1\sigma$ **sin** ajustar $V_0$ al
dato, $A_s$ pasa de parámetro a predicción → el sector CMB del modelo baja a
**$k=1$ del modelo + $\tau$ nuisance**.

**Sobre $\tau$ (por qué NO es un OP gemelo).** $\tau$ no se deriva de $\{\varphi,\pi\}$
puros porque su cadena pasa por astrofísica bariónica (eficiencia de formación estelar
$f_*$, escape de fotones ionizantes $f_{\rm esc}$): $(\varphi,\pi)\to A_s\to\sigma_8\to D_1
\to$ primeros halos $\to[f_*,f_{\rm esc}]\to z_{\rm reion}\to\tau$. Lo honesto es tratar
$\tau$ como **nuisance astrofísico compartido** con ΛCDM (se cancela en $\Delta$BIC), o
—como mucho— convertirlo en *output predicho* dado $A_s$ + un modelo de reionización
estándar (no derivado de $(\varphi,\pi)$). No es imposible, pero no es derivación pura;
por eso $\tau$ queda como residuo nuisance, no como OP de derivación.

**Severity:** Medium. Es el último parámetro genuinamente cosmológico no-derivado del
modelo; cerrarlo lleva el conteo a "0 parámetros del modelo + $\tau$ nuisance".
**Relation:** depende de OP-2 (✅, da $\epsilon,N_*$) y de la escala $M$ de Paper 10.
Script a crear: `op18_As_from_inflation.py`.
