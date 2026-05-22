# Registro de Verificación — Modelo SSEE-V3.6

> Reconstrucción verificada del modelo. No se reescribe desde cero: cada elemento
> se **verifica** — qué es, qué función cumple, dónde se usa, dónde debe estar, y
> con qué conecta — antes de declararse resuelto y, finalmente, sellarse.

## Por qué existe este registro

La revisión árbitro hostil de los 11 documentos (2026-05-21) encontró que varias
"resoluciones" marcadas `✅ RESUELTO` en `CLAUDE.md` nunca fueron verificadas de
verdad. Caso confirmado por `git`: la fórmula k-mouflage de Paper 8
(`r_km³ = M_obj/(4π M_pl M²)`) es **dimensionalmente inconsistente** y se
introdujo en el commit `295ed6e` ("OP-4 resolved"). Es decir: una corrección
creó un error nuevo, se commiteó como "fix", se escribió como resuelto, y pasó
**siete auditorías posteriores** sin que nadie lo detectara — porque esas
auditorías revisaban strings, citas y retórica, nunca la física.

Este registro existe para que eso no vuelva a pasar: nada se da por resuelto sin
pasar las seis comprobaciones, y nada se sella sin re-verificarse.

## Ciclo de vida de un elemento

1. **pendiente** — identificado, sin verificar.
2. **verificado** — pasó las 6 comprobaciones.
3. **resuelto** — confirmado por Mike; se le asigna el nuevo estado resuelto.
4. **re-verificado** — el estado resuelto se comprueba otra vez, de forma
   independiente. El elemento queda **aprobado**.

## Las 6 comprobaciones

Para cada elemento se verifica y se registra:

1. **Numérica** — el valor o la identidad se calcula y coincide.
2. **Dimensional** — las unidades cierran (esto es lo que faltó en Paper 8).
3. **Derivación** — el elemento se sigue de sus premisas; no es post-hoc ni un fit.
4. **Rol** — qué función cumple en el modelo está claro y es necesario.
5. **Ubicación** — está en el lugar correcto y se usa de forma consistente en
   todos los papers donde aparece.
6. **Conexiones** — de qué depende y qué depende de él; el grafo cierra.

## Sello (por paper)

Un paper se **sella** solo cuando *todos* sus elementos están `re-verificado` **y**
pasa un chequeo completo final. El sello es invisible — no aparece en el PDF: se
registra aquí el **commit git + `sha256` del `.tex`** al momento del sellado. Si
el archivo cambia después, el `sha256` deja de coincidir → el sello se rompe → el
paper vuelve a verificación.

## Capas de verificación (orden ascendente — bottom-up)

| Capa | Contenido | Estado |
|------|-----------|--------|
| **L1** | Axiomas y constantes algebraicas | re-verificado ✓ |
| **L2** | Parámetros cosmológicos derivados | re-verificado ✓ |
| **L3** | Mecanismos y derivaciones (OP-1..OP-7, EFT, IS, dos-sectores, k-mouflage, f_screen, m_φ, K(X) UV) | en progreso (OP-1..OP-7 ✓; faltan EFT, IS, dos-sectores, m_φ, K(X)) |
| **L4** | Confrontaciones con datos (MCMC, CMB, ΔBIC, fσ₈, S₈) | pendiente |
| **L5** | Papers (sellado) | pendiente |

Regla: un elemento de una capa no pasa de `verificado` si sus insumos de capas
inferiores no están al menos `verificado`.

---

# Capa 1 — Axiomas y constantes algebraicas

Verificación numérica: `python3` (2026-05-21). Todas las constantes son
adimensionales → comprobación dimensional trivialmente ✓. Derivación: V-L1-01/02
son axiomas; el resto son definiciones algebraicas de φ y π.

| ID | Constante | Definición | Valor verificado | Rol en el modelo | Estado |
|----|-----------|-----------|------------------|------------------|--------|
| V-L1-01 | φ | (1+√5)/2 | 1.6180339887 | Axioma generador | re-verificado |
| V-L1-02 | π | — | 3.1415926536 | Axioma generador | re-verificado |
| V-L1-03 | Ω | φ+π | 4.7596266423 | Métrica de estabilidad; base de H₀^alg=3Ω² | re-verificado |
| V-L1-04 | β | (φ+π)/2 | 2.3798133212 | Escalar de acoplamiento base | re-verificado |
| V-L1-05 | AURA | (3φ+π)/2 | 3.9978473099 | Acoplamiento EFT βc; genera MIRA | re-verificado |
| V-L1-06 | MIRA | AURA/2 | 1.9989236550 | Razón Ω_m,cosm / Ω_m,dyn | re-verificado |
| V-L1-07 | KAL₀ | (φ+3π)/2 | 5.5214059748 | Viscosidad estructural; fija τ_Π | re-verificado |
| V-L1-08 | T_r | 3(φ+β) | 11.9935419298 | Horizonte de saturación 3D; numerador de w₀ | re-verificado |
| V-L1-09 | K_v | 2(φ+π) | 9.5192532847 | Invariante de restricción estructural | re-verificado |
| V-L1-10 | M_v | φ+π+K_v | 14.2788799270 | Invariante dimensional máximo; denominador de w₀ | re-verificado |

### Conexiones L1 (grafo de dependencias)

- φ, π → todo.
- Ω = φ+π → H₀^alg (L2), K_v.
- β = Ω/2 → T_r, KAL₀ (vía β+π... ver nota).
- AURA → MIRA → Ω_m,cosm (L2), βc (L2/L3).
- T_r, M_v → w₀ (L2).
- KAL₀ → τ_Π, viscosidad IS (L3).

### Notas de verificación L1 (hallazgos)

- **V-L1-07 KAL₀** — definición canónica `(φ+3π)/2`. En varios papers aparece
  también como `β+π`; es equivalente (`β+π = (φ+π)/2+π = (φ+3π)/2`) ✓, pero
  conviene unificar a una sola forma. Comprobación de ubicación: **pendiente de
  barrido cross-paper.**
- **V-L1-09 K_v** — Paper 1 la define con dos formas (`φ+π+Ω` y `2(φ+π)`),
  numéricamente idénticas (Ω=φ+π). Unificar a `2(φ+π)`. Mismo caso M_v
  (`φ+π+K_v` vs `3(φ+π)`). Comprobación de ubicación: **observación abierta.**
- Nomenclatura: AURA, MIRA, KAL₀ se conservan como identificadores físicos
  establecidos. Los nombres mitológicos (BIAL, KRYSTOS, TRIAL, MIKAEL_V, PYROS,
  Ω_DNAV) deben eliminarse — limpieza en curso (Clase A).

**Estado Capa 1:** 10/10 `re-verificado` (confirmado por Mike 2026-05-21;
re-verificado de forma independiente por el guardián `ssee_verify.py`).
Pendiente: barrido de la comprobación 5 (ubicación) en los 11 `.tex`.

---

# Capa 2 — Parámetros cosmológicos derivados

Verificación numérica: `python3` (2026-05-21), todas las identidades recomputadas
desde las constantes de Capa 1. La columna **Dim.** indica si la comprobación
dimensional pasa.

| ID | Parámetro | Definición | Valor verificado | Dim. | Estado |
|----|-----------|-----------|------------------|------|--------|
| V-L2-01 | w₀ | −T_r/M_v | −0.8399497713 | ✓ | verificado |
| V-L2-02 | wₐ | −P_sc/K_v  (P_sc=Ω+φ=6.3776617) | −0.6699748857 | ✓ | verificado |
| V-L2-03 | Ω_DE | T_r/M_v | 0.8399497713 | ✓ | verificado |
| V-L2-04 | Ω_m,dyn | 1+w₀ | 0.1600502287 | ✓ | verificado |
| V-L2-05 | Ω_m,cosm | MIRA·Ω_m,dyn | 0.3199281880 | ✓ | verificado |
| V-L2-06 | H₀^alg | 3Ω² | 67.9621373234 | ✗ | **ABIERTO** |
| V-L2-07 | n_s | 1−φ⁻⁷ | 0.9655581463 | ✓ | verificado |
| V-L2-08 | αK | 3·Ω_DE·Ω_m,dyn | 0.4033024589 | ✓ | verificado |
| V-L2-09 | βc | −AURA | −3.9978473099 | ✓ | verificado |
| V-L2-10 | m_φ | Σm_ν^act·H₀^alg (0.0824·67.96) | 5.6001 eV | ✗ | **ABIERTO** |
| V-L2-11 | k_fs | free-streaming de m_φ (vía DW) | 0.493 h/Mpc | — | pendiente L3 |
| V-L2-12 | r | 12α/N²  (α=φ⁴/3, N=2φ⁷) | 0.00813062 | ✓ | verificado |
| V-L2-13 | f_screen | αK/(3·MIRA) = (π−φ)/Ω² | 0.0672532703 | ✓ | verificado |

### Cross-checks de identidad (todas pasan)

- **w₀**: dos rutas coinciden — `−T_r/M_v` y `1/(2n−1)` con `n=(T_r−M_v)/(2T_r)`.
- **f_screen**: dos fórmulas coinciden exactamente — `αK/(3·MIRA)` y `(π−φ)/Ω²`.
- **Ω_m,dyn + Ω_DE = 1** ✓; **Ω_m,cosm = MIRA·Ω_m,dyn** ✓.

### Problemas ABIERTOS detectados en Capa 2

- **V-L2-06 H₀^alg = 3Ω²** — numéricamente da 67.962, pero `3Ω²` es
  **adimensional** mientras que H₀ tiene unidades km/s/Mpc. No es una
  derivación dimensional; es una coincidencia numérica (Postulado D en P1).
  El modelo lo admite a medias. **No puede pasar a `resuelto` como derivación.**
- **V-L2-10 m_φ = Σm_ν·H₀^alg** — numéricamente da 5.60 eV, pero el producto
  es `[eV]·[km/s/Mpc]`, que **no es [eV]**. Solo da "5.60 eV" si H₀ se trata
  como el número puro 67.96. Es numerología dimensional (hallazgo crítico del
  árbitro de P6). Hereda además el problema de H₀^alg. **ABIERTO.**

### Dependencias hacia Capa 3 (derivación — comprobación 3)

Estos parámetros pasan numérica y dimensionalmente, pero su **derivación** se
apoya en mecanismos de Capa 3 aún no verificados — no pueden pasar de
`verificado` a `resuelto` hasta que Capa 3 los sostenga:

- **n_s** (V-L2-07): el exponente 7 depende de OP-2 (N_*=2φ⁷ + α-attractor).
- **βc** (V-L2-09): depende de OP-7 (unicidad EFT; el árbitro lo llamó fit al 0.2 %).
- **αK** (V-L2-08): depende de la acción EFT de P7.
- **r** (V-L2-12): depende de α=φ⁴/3 y N=2φ⁷ (OP-2).
- **f_screen** (V-L2-13): la *identidad algebraica* está verificada, pero el
  *mecanismo* de screening depende de P9 (el árbitro lo llamó circular).
- **k_fs** (V-L2-11): depende de m_φ (ABIERTO) y de la relación DW (L3).

**Estado Capa 2:** 10/13 `verificado` (numérica + dimensional + identidades);
2 `ABIERTO` (H₀^alg, m_φ — dimensionalmente inconsistentes); 1 `pendiente L3`
(k_fs). Ninguno pasa a `resuelto` todavía: la comprobación 3 (derivación) de
varios alcanza la Capa 3.

# Capa 3 — Mecanismos y derivaciones

Aquí vive la "nueva física". Cada OP-1..OP-7 está marcada `✅ RESUELTO` en
CLAUDE.md — se re-verifican una por una, paso a paso. Estado por elemento:
`verificado` (la derivación cierra), `ABIERTO` (coincidencia/conjetura vestida
de derivación), o `PARCIAL` (mezcla — partes verificadas, partes abiertas).

## V-L3-OP2 — n_s = 1 − φ⁻⁷ (índice espectral) — **PARCIAL**

*Claim CLAUDE.md:* "✅ RESUELTO — α-attractor universality + N_*=2φ⁷".

Cadena de derivación, paso a paso:

1. **✓** Universalidad α-attractor: `n_s = 1 − 2/N_*` (Kallosh & Linde 2013).
   Física estándar correcta.
2. **⚠ pendiente** `α = φ⁴/3` — insumo de Paper 1, aún sin verificar
   (elemento V-L3-alpha, pendiente en esta misma capa).
3. **✓** Álgebra exacta: con `N_* = 2φ⁷` → `n_s = 1−2/(2φ⁷) = 1−φ⁻⁷`, y
   `r = 12(φ⁴/3)/(2φ⁷)² = φ⁻¹⁰`. Ambas identidades exactas — el guardián
   las recomputa.
4. **✗** `N_* = 2φ⁷` — **NO derivado.** Es la Conjecture B.1. El script
   `ssee_paperB_Nstar.py` es honesto: *invierte* la fórmula para hallar el
   T_rh que da 2φ⁷; la cuasi-coincidencia ρ_end/ρ_rh≈3 necesita un ajuste
   δ~O(1/N) en la constante 58.25 de la fórmula estándar para cerrar; y el
   puente físico (eficiencia de reheating gravitacional con α=φ⁴/3 → ρ_rh=V_end)
   está **ausente**. El argumento alterno "contar 7 constantes SSEE" es
   racionalización post-hoc — el conteo depende de cómo se agrupen.

**Veredicto:** dado N_*=2φ⁷, todo cierra exacto. Pero N_*=2φ⁷ es una conjetura
no probada. OP-2 NO está "RESUELTO": es **condicional a la Conjecture B.1**.

## V-L3-OP7 — βc = −AURA (acoplamiento EFT) — **PARCIAL**

*Claim CLAUDE.md:* "PARCIALMENTE RESUELTO — unicidad EFT vía dualidad Z₂".

1. **✓** Extracción numérica: integrar el fondo EFT (P7 §6, shooting con
   Ω_φ(a=1)=Ω_DE) da `βc ≈ −3.990`, sin parámetros libres.
2. **✗** Identificación `βc = −AURA`: −AURA = −3.99785; el valor extraído
   −3.990 está a **0.2 %** (|Δ|/βc = 0.196 %). La ecuación `\boxed{βc=−AURA}`
   de P7 lo presenta como exacto — no lo es. El origen del 0.2 % se atribuye
   a la aproximación de shooting, pero el árbitro halló que el plateau test
   muestra que NO viene de las condiciones iniciales. "Bariones+radiación lo
   llevarían a −AURA dentro de 0.01 %" es una predicción no ejecutada.
3. **✓** Dualidad Z₂: `KAL₀=(φ+3π)/2 ↔ AURA=(3φ+π)/2` bajo φ↔π es una
   identidad algebraica exacta — el guardián la recomputa.
4. **✗** Pero la dualidad **no genera** βc. Es una relación entre dos
   constantes hechas con la misma plantilla, no una simetría de la acción EFT
   (V₀e^{αφ}, K(X) no son φ↔π-invariantes). "βc queda determinado por KAL a
   través de la dualidad" es un overclaim lógico: una coincidencia notacional
   no es una derivación.

**Veredicto:** el resultado numérico βc≈−3.990 es sólido; `βc=−AURA` es una
coincidencia numérica al 0.2 %, no una derivación; la dualidad Z₂ es álgebra
real pero descriptiva, no generativa. **ABIERTO.**

## V-L3-alpha — α = φ⁴/3 (parámetro α-attractor) — **verificado (como consecuencia de axiomas)**

Al verificar de dónde sale α=φ⁴/3 se descubre la lógica real de P1
(§Inflationary Embedding) — y es más honesta de lo que CLAUDE.md sugiere:

1. **P1 declara explícitamente** que `n_s=1−φ⁻⁷` y `r=φ⁻¹⁰` son **AXIOMAS**
   (postulados algebraicos de {φ,π}): *"Neither is derived from the
   α-attractor framework"* (P1 EFT, L498). Son predicciones falsables.
2. **✓** Dados esos dos axiomas: `N = 2/φ⁻⁷ = 2φ⁷`, y
   `α = r·N²/12 = φ⁻¹⁰·(2φ⁷)²/12 = φ⁴/3` — consecuencia algebraica exacta
   (guardián la recomputa). Identidad Fibonacci: φ⁴=3φ+2 ⇒ α=φ+2/3.
3. El resultado genuino: los dos axiomas son **mutuamente consistentes** —
   seleccionan una única geometría α-attractor. Eso es honesto y correcto.

**Hallazgo nuevo — error aritmético en P1 (corregido):** la curvatura de
Kähler. P1 escribía `R = −2/(3α) = −2/φ⁴ = −φ⁻⁴ ≈ −0.146`; el último paso
perdía un factor 2. Correcto: `−2/φ⁴ = −2φ⁻⁴ ≈ −0.292`. Corregido en
`SSEE_EFT_section.tex` (eq:kahler_curvature + tabla). El guardián lo comprueba.

**Contradicción de framing detectada:** P1 dice que n_s es un AXIOMA;
`ssee_op2_spectral_index.py` y CLAUDE.md dicen "exponent 7 → derivado" y
marcan OP-2 "✅ RESUELTO". No pueden ser ambas ciertas. La honesta es la de
P1: n_s y r son postulados. El error está en CLAUDE.md y el script, no en P1.

**Veredicto:** α=φ⁴/3 `verificado` como consecuencia exacta de los axiomas
n_s, r — no es derivación independiente, y P1 nunca afirmó que lo fuera.

## V-L3-OP4 — radio de screening k-mouflage (P8) — **ABIERTO (fórmula rota)**

*Claim CLAUDE.md:* "OP-4 ✅ RESUELTO 2026-05-15 — k-mouflage + αB=αM=0".

`git blame` confirma que la fórmula `r_km³ = M_obj/(4π·M_pl·M²)` (P8
eq:rkm) nació en el commit `295ed6e` ("OP-4 resolved — replace Galileon
Vainshtein with k-mouflage") y pasó 7 auditorías posteriores sin detectarse.

**✗ Comprobación dimensional (falla):** el lado derecho tiene dimensión
`[GeV]/([GeV]·[GeV²]) = GeV⁻²`. Entonces `r_km = (RHS)^{1/3}` tiene dimensión
`GeV^{−2/3}` — **no es una longitud** (una longitud es GeV⁻¹). La fórmula es
dimensionalmente inconsistente. Toda la §4–5 de P8 (Tabla 2, Fig. 2, ejemplo
A1689) descansa sobre ella.

**Forma dimensionalmente correcta:** del propio Lagrangiano de P8
(K(X)=X/KAL+X²/M⁴, cruce X=M⁴/KAL, cierre de gradiente
∇φ=2·KAL·βc·M_pl·∇Φ_N) el cruce da `r_km⁴ ∝ M_obj²/(M_pl²·M⁴)` — raíz
cuarta con escala M_obj², no la raíz cúbica de M_obj¹ del paper.

**Veredicto:** OP-4 NO está resuelto. La "resolución" introdujo un error de
derivación. Requiere re-derivar el radio k-mouflage desde cero y propagar a
P8 §4–5 (tabla, figura, ejemplo A1689). **ABIERTO** — remediación grande,
pendiente de sesión dedicada.

## V-L3-OP1 — Ω_b h² = (π−φ)/(3Ω²) (densidad bariónica) — **ABIERTO (coincidencia escaneada)**

*Claim CLAUDE.md:* "OP-1 PARCIAL — (π−φ)/H₀_SSEE = 0.32σ".

1. **✓ numérico:** (π−φ)/(3Ω²) = 1.52356/67.962 = 0.022418; Planck
   0.02237±0.00015 → tensión **0.32σ**. La identidad pura (número/número)
   es dimensionalmente correcta.
2. **✗ no derivado:** la fórmula sale de un **scan** — `ssee_op1_baryon_density.py`
   Paso 5 prueba 7 expresiones {(π−φ)/(3Ω²), 3(π−φ)/200, (π−φ)/φ¹¹, …} y
   elige la de menor tensión. Eso es ajuste a un objetivo conocido.
3. **✗ interpretación rota:** el script reescribe (π−φ)/(3Ω²) como
   "(π−φ)/H₀_SSEE" — pero con H₀ en km/s/Mpc eso tiene unidades, no es
   adimensional. La "derivación desde la tasa de esfalerón" se difiere a
   Paper B/C.

**Veredicto:** una coincidencia numérica decente (0.32σ) hallada por scan,
no una derivación. **ABIERTO.**

## V-L3-OP3 — separabilidad UV-IR / KALeff = φ²√(5/2) — **ABIERTO (cota asertada, dimensión confusa)**

*Claim CLAUDE.md:* "OP-3 RESUELTO — jerarquía EFT (H₀/M)²≈10⁻⁶²".

1. **✓** La jerarquía es real: (H₀/M)² ≈ 2.7×10⁻⁶² (M=8.81 meV ≫ H₀).
2. **✗** Que esa jerarquía *pruebe* la separabilidad (Postulate C.1 →
   Theorem C.1) es una **aserción**: `ssee_op3_separability.py` (L176–179)
   admite que el jacobiano ∂φ/∂χ que mide el mezclado real se difiere a
   Paper B. "RESUELTO (cota EFT)" asevera la cota y aplaza el cálculo.
3. **✗ dimensión confusa:** KALeff² = M⁴/(6α). Con M⁴=5φ⁸ρ_crit y 6α=2φ⁴
   da KALeff² = (5/2)φ⁴·ρ_crit — el script escribe "5φ⁴/2", **dropeando
   ρ_crit**. KALeff resulta dimensional (∝GeV²) pero se factoriza
   `KAL₀=KALeff·F` contra KAL₀=(φ+3π)/2, que es adimensional.
4. El script contiene una auto-corrección sin resolver
   (√(6α)=φ² → "corrección: =φ²√2") — la derivación no está limpia.

**Veredicto:** la jerarquía (H₀/M)²~10⁻⁶² es un hecho; la "prueba" de
separabilidad no lo es (jacobiano diferido), y KALeff arrastra un ρ_crit
dropeado. OP-3 NO está "RESUELTO". **ABIERTO.**

## V-L3-OP5 — tensión S₈ weak-lensing / HMcode bariónico — **ABIERTO (anclado en rama secundaria)**

*Claim CLAUDE.md:* "OP-5 PARCIAL — HMcode-2020 CLASS: S₈=0.758 (0.06σ DES)".

1. **✓ definición:** S₈ = σ₈(Ω_m/0.3)^½ cierra para ambas ramas — rama
   secundaria σ₈=0.737 → S₈=0.761; rama titular σ₈=0.794 → S₈=0.820.
2. **✗ rama elegida:** `ssee_op5_hmcode.py` aplica la supresión bariónica
   sobre `S8_Paper6 = 0.761` — la **rama WDM CLASS, marcada "secundario, NO
   titular"** en CLAUDE.md. El resultado titular del Paper 6 es S₈_eff=0.820.
3. **✗ magnitud del efecto:** HMcode da B_σ₈ ≈ 0.996 (supresión ~0.4% en
   σ₈, porque R=8 Mpc/h cae en k~0.13 h/Mpc, donde el feedback AGN apenas
   actúa). Aplicado a 0.761 → 0.758 (0.06σ); aplicado al titular 0.820 →
   ~0.817 (~2.6σ). La "resolución" depende de anclar en la rama baja.
4. **✗ N-body diferido:** el cierre real (<1σ) requiere simulaciones N-body
   bariónicas SSEE (Nivel 2, ~5k–20k CPU-h) — diferido.

**Veredicto:** HMcode reduce la tensión solo escogiendo el S₈ más bajo de
dos valores internamente inconsistentes. El titular sigue a ~2.6σ. **ABIERTO.**

## V-L3-OP6 — forma de screening f_screen / universo separado — **PARCIAL (forma derivada, valor con insumo)**

*Claim CLAUDE.md:* "OP-6 ✅ RESUELTO — universo separado k-essence + identidad 1+w₀=Ω_m".

1. **✓ valor:** f_screen = α_K/(3·MIRA) = (π−φ)/Ω² = 0.067253 — álgebra
   exacta, ya verificada en V-L2-13 y en la identidad cruzada de Capa 2.
   H₀,local = H₀^alg/(1−f_screen) = 72.86 km/s/Mpc cierra numéricamente.
2. **✓ forma:** que la corrección sea **multiplicativa** sí sigue de la
   aproximación de universo separado para k-essence (Wands 2000; Brax &
   Valageas 2014) — ese paso es una derivación legítima.
3. **✗ paso δρ_φ asertado:** `ssee_op6_screening_form.py` (L68) escribe
   δρ_φ/ρ_crit = (α_K/3)(Ω_m,dyn/MIRA)δ_local/(1+w₀) sin derivarla; el
   factor 1/MIRA se justifica con un argumento de plausibilidad, no un
   cálculo. El c²_s aparece y desaparece entre L65 y L73.
4. **✗ insumo δ_local=2:** el valor f_screen=α_K/(3·MIRA) exige fijar
   δ_local=2 (sobredensidad del Grupo Local) para cancelar el δ_local/2.
   Es un insumo astrofísico razonable, **no derivado de φ,π**.

**Veredicto:** la forma multiplicativa está derivada; el valor f_screen es
álgebra exacta *condicionada* a δ_local=2 y a una expresión δρ_φ asertada.
No es "RESUELTO" pleno. **PARCIAL.**

## V-L3-mphi — masa del campo φ-DM, m_φ = 5.60 eV — **ABIERTO (numerología + dimensión)**

*Claim CLAUDE.md:* "m_φ = Σm_ν^active × H₀^alg = 5.60 eV — cero parámetros libres".

1. **✓ numérico:** la cadena cierra — R=4·KAL₀−22=0.085624, Σm_ν^active =
   R·Ω_b h²·94.07/τ_Π = 0.082407 eV, m_φ = Σm_ν^active·H₀^alg = 5.6005 eV.
2. **✗ R = 4·KAL₀−22:** se resta el **entero pelado 22** a 4·KAL₀=22.0856
   para extraer el remanente 0.0856. No hay derivación de por qué 4, ni de
   por qué 22. Es ingeniería inversa de un número pequeño.
3. **✗ dimensional:** m_φ = [eV]·H₀^alg. H₀^alg=3Ω²=67.96 se usa como
   número adimensional (si fuese km/s/Mpc, m_φ no sería [eV] — V-L2-10).
   Multiplicar una energía por 67.96 sin mecanismo físico no da una masa.
4. **✗ insumo heredado:** Ω_b h² es la coincidencia escaneada de OP-1
   (**ABIERTO**). m_φ se construye encima de un input no derivado.

**Veredicto:** la cadena aritmética cierra a 5.60 eV, pero descansa en una
resta de entero pelado, un producto dimensionalmente injustificado y el
input abierto OP-1. **ABIERTO.**

## V-L3-2sec — modelo dos sectores φ-DM — **PARCIAL (identidad sí, split físico no)**

*Claim CLAUDE.md:* "Ω_total (dos sectores) = 0.319928 ≈ Ω_m,CMB — unificación algebraica".

1. **✓ identidad:** Ω_CDM + Ω_φDM = Ω_m,dyn + (MIRA−1)·Ω_m,dyn =
   MIRA·Ω_m,dyn = 0.319928. Diferencia con V-L2-05 = 0 exacto. Es una
   **re-partición algebraica** de Ω_m,cosm en dos mitades casi iguales.
2. **✗ split físico:** que un sector (φ-DM) free-streame para k>k_fs y el
   otro (CDM) no, depende de m_φ (**ABIERTO**, ver arriba) y de k_fs
   (**pendiente** en Capa 2). El mecanismo Dodelson-Widrow que fija k_fs
   parte de m_φ.

**Veredicto:** la suma Ω_total es un re-enunciado exacto de V-L2-05; el
modelo físico de dos sectores hereda la apertura de m_φ y k_fs. **PARCIAL.**

## Pendientes de Capa 3

La acción EFT (Paper 7); las perturbaciones Israel-Stewart; la completación
UV K(X) (Paper 10); el bug dos-Ω_m (uso de Ω_m,dyn vs Ω_m,cosm en E(z),
Poisson y CMB).

**Estado Capa 3:** re-verificados OP-1..OP-7, α=φ⁴/3, m_φ y dos sectores
(10 elementos). Resultado: α=φ⁴/3 **verificado**; OP-2, OP-6, OP-7, dos
sectores **PARCIAL** (álgebra/forma cierra, insumo físico no); OP-1, OP-3,
OP-4, OP-5, m_φ **ABIERTO**. 1 bug aritmético corregido de paso (curvatura
de Kähler, P1). Patrón confirmado en 10/10: ninguna "✅ RESUELTO" de
CLAUDE.md lo estaba sin reservas.

# Capa 4 — Confrontaciones con datos — *pendiente*

MCMC DESI+Planck, espectros CMB, ΔBIC, fσ₈, S₈, tensión de Hubble.

# Capa 5 — Sellado de papers — *pendiente*

| Paper | Todos los elementos re-verificados | Chequeo completo | Sello (commit + sha256) |
|-------|-----------------------------------|------------------|-------------------------|
| P1–P10, Unified | — | — | — |

---

*Registro iniciado 2026-05-21 tras la revisión árbitro hostil de los 11 documentos.*
