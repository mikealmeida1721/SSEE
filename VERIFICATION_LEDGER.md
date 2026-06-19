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
| **L3** | Mecanismos y derivaciones (OP-1..OP-7, EFT, IS, dos-sectores, k-mouflage, f_screen, m_φ, K(X) UV, T_μν, disformal, retención MIRA) | re-verificado ✓ (17 elementos; 1 verif., 7 PARCIAL, 9 ABIERTO — dos-Ω_m es el central) |
| **L4** | Confrontaciones con datos (MCMC, CMB, ΔBIC, fσ₈, S₈) | re-verificado ✓ (pipelines re-corridos; CMB+MCMC reproducen, r_d/H₀ derivaron) |
| **L5** | Papers (sellado) | pendiente |

Regla: un elemento de una capa no pasa de `verificado` si sus insumos de capas
inferiores no están al menos `verificado`.

---

# Valores Canónicos del Modelo SSEE

**Esta es la fuente canónica.** "Canónico" no significa *permanente* —
significa **el valor que refleja el estado real más actual del modelo**.
Cuando un pipeline se re-corre y el número cambia, se actualiza **aquí
primero** y luego se propaga a todo lo que lo use.

**Regla estructural.** Si un resultado cambia, cambia en *todas* partes
que lo usan. Cualquier script, paper o cálculo que use uno de estos
números debe reflejar el valor de esta tabla. Un valor distinto sin
justificación es un error que debe detectarse de inmediato.

## A. Constantes algebraicas — invariantes (derivadas de φ, π)

Fuente única: `src/ssee_core.py`. Todo script importa de ahí. El guardián
(sección «Fuente canónica») re-computa cada una y verifica `ssee_core`
contra esa recomputación — si el módulo se edita mal, el guardián → ROJO.

| Símbolo | Valor | Identidad | Rol |
|---|---|---|---|
| φ | 1.6180339887 | (1+√5)/2 | Axioma generador |
| π | 3.1415926536 | — | Axioma generador |
| Ω | 4.7596266423 | φ+π | Métrica de Estabilidad |
| β | 2.3798133212 | (φ+π)/2 | Escalar de Acoplamiento Base |
| KAL₀ | 5.5214059748 | β+π | Viscosidad Estructural |
| P_sc | 6.3776606311 | Ω+φ | Escalar de Evolución Dinámica |
| K_v | 9.5192532847 | φ+π+Ω | Restricción Estructural |
| T_r | 11.9935419298 | 3(φ+β) | Horizonte de Saturación 3D |
| M_v | 14.2788799270 | φ+π+K_v | Invariante Dimensional Máximo |
| AURA | 3.9978473099 | (3φ+π)/2 | = 2·MIRA = φ+β |
| MIRA | 1.9989236550 | AURA/2 | Frecuencia de Observación |
| w₀ | −0.8399497713 | −T_r/M_v | Ecuación de estado hoy |
| wₐ | −0.6699748857 | −P_sc/K_v | Evolución de la EoS |
| Ω_DE | 0.8399497713 | T_r/M_v | Densidad de energía oscura |
| Ω_m,dyn | 0.1600502287 | 1+w₀ | Sector dinámico (fija w₀; **no** en E(z)) |
| Ω_m,cosm | 0.3199281880 | MIRA·Ω_m,dyn | Fondo gravitacional (E(z)/Poisson/CMB) |
| H₀^alg | 67.9621373234 | 3(φ+π)² | H₀ algebraico — **invariante** |
| n_s | 0.9655581463 | 1−φ⁻⁷ | Índice espectral |
| α_K | 0.4033024589 | 3·Ω_DE·Ω_m,dyn | Kineticity EFT |
| Ω_b h² (alg) | 0.0224177568 | (π−φ)/(3Ω²) | Densidad bariónica OP-1 (**ABIERTO**) |

## B. Valores de pipeline — dependientes de estado (script + datos + fecha)

Estos **no** se derivan de φ,π — los calcula un pipeline (CAMB/CLASS/emcee).
Cambian si el script o los datos cambian. Cada uno lleva su **procedencia**.

| Cantidad | Valor canónico | Fuente | Re-anclado |
|---|---|---|---|
| H₀ anchor/prior (H_MIRA, CMB-óptimo) | 67.037 km/s/Mpc | `ssee_paper3_cobaya_unified.py` (SSEE@Σm_ν=0.069 self-consistente; era 67.068 con mnu=0.06 baseline) | 2026-06-09 |
| H₀ MCMC posterior (prior MIRA 67.037, corregido por DESI) | 66.533 ± 0.442 km/s/Mpc | `ssee_paper2_mcmc_mira.py` (100w×25k, seed 42, prior 67.037) → `results/logs/mcmc_paper2_mira.log` | 2026-06-09 |
| ΔBIC MCMC (ΛCDM−SSEE) | +7.91 (SSEE favorecido) | `ssee_paper2_mcmc.py` | 2026-05-22 |
| Ω_b h² (posterior MCMC) | 0.02183 ± 0.00048 | `ssee_paper2_mcmc.py` | 2026-05-22 |
| r_d,SSEE (crudo) | 175.16 Mpc | `ssee_paper2_mcmc.py` | 2026-05-22 |
| r_d,eff (CAMB, en anchor 67.037, mnu=0.069) | 146.73 Mpc — **1.38σ** | `ssee_verify_rd.py` (re-run 2026-06-09) | 2026-06-09 |
| r_d,eff (CAMB, en posterior 66.533, mnu=0.069) | 147.30 Mpc — **0.80σ** | `ssee_verify_rd.py` (re-run 2026-06-09) | 2026-06-09 |
| χ²_r CMB TT (SSEE) | 1.044 | `ssee_paper3_cmb.py` (Σm_ν=0.0690) → `results/logs/paper3_cmb_canonical.log` | 2026-06-14 |
| ΔBIC CMB diagonal (SSEE−ΛCDM) | −28.0 (SSEE favorecido) | `ssee_paper3_cmb.py` (Σm_ν=0.0690) → `results/logs/paper3_cmb_canonical.log` | 2026-06-14 |
| ΔBIC CMB plik_lite Cobaya (k=2) | −32.2 (ΛCDM@0.06·SSEE@0.069; H₀_opt=67.037) | `ssee_paper3_cobaya_unified.py` | 2026-06-09 |
| θ* (CAMB, en anchor 67.037, mnu=0.069) | 0.59638° — **0.66σ** | `ssee_verify_rd.py` (re-run 2026-06-09) | 2026-06-09 |
| θ* (CAMB, en posterior 66.533, mnu=0.069) | 0.59423° — **5.33σ** (sensibilidad de θ* a H₀; ver V-L4-θ*) | `ssee_verify_rd.py` (re-run 2026-06-09) | 2026-06-09 |
| σ₈ / S₈ SSEE (Paper 5 single-sector, "el desafío") | 0.820 / 0.847 — **3.9σ KiDS** | `ssee_paper5_IS_perturbations.py` (G=1.011, Ω_m,CMB) | 2026-06-13 (canónico; cadena 0.702/0.725 con G=0.866 retirada) |
| σ_LS (Paper 6, amplitud RSD k<k_fs, fσ₈) | 0.811·G_2s = 0.794 | `ssee_paper6_verification.py` | 2026-06-04 (canónico) |
| σ_eff / S₈ (Paper 6 titular lensing, R=8 cruza k_fs) | 0.742 / 0.766 — **0.01σ KiDS** | `ssee_paper6_verification.py` | 2026-06-04 (canónico) |
| m_φ (Paper 6 forward-pred, dim. consistente) | 36.9463 eV | Σm_ν^act·(Ω⁴+AURA·KAL₀), `ssee_core` | 2026-06-04 (canónico) |
| k_fs (Paper 6, output CLASS) | 0.659 h/Mpc | CLASS (`calibrate_wdm_alpha.py`, α=1.243 Mpc/h) | 2026-06-04 (canónico) |

**Historial de deriva de H₀ MCMC** (para entender por qué cambió): el MCMC
es determinista (semilla fija 42) — *mismo script → mismo número*. La
deriva corresponde a **ediciones del script / cambios de prior**, no a azar:
- 66.75 ± 0.44 — prior Planck-ΛCDM legacy, 50w×10k.
- 67.756 ± 0.442 — prior Planck-ΛCDM, 100w×25k (commit `4892b53`, corrección
  dos-Ω_m). **OBSOLETO desde el switch de prior 2026-05-24.**
- **66.533 ± 0.442 — prior MIRA self-consistent, 100w×25k (CANÓNICO ACTUAL).**
  `ssee_paper2_mcmc_mira.py`, ref `results/logs/mcmc_paper2_mira.log` L78.
  *Re-run 2026-06-09:* corrido con el prior canónico nuevo 67.037 (mnu=0.069
  self-consistente). El posterior pasó de 66.553 (prior 67.068) a **66.533**
  (corrimiento −0.020, = 0.06σ del prior, como se predijo). N_eff≈77709,
  acceptance 0.715. Este es el valor canónico actual.

**Anatomía del H₀ (un solo número, dos etapas):** el modelo tiene UN H₀.
(1) **Anchor/prior** = H_MIRA = 67.037 km/s/Mpc: el H₀ que minimiza la
tensión CMB (Cobaya `minimize_scalar` sobre Planck plik_lite con fondo SSEE,
con Σm_ν=0.069 SSEE self-consistente; era 67.068 con mnu=0.06 baseline).
(2) **Posterior** = 66.533 ± 0.442: el mismo H₀ tras dejar que el MCMC ajuste
DESI DR2 BAO encima del prior. DESI lo baja ~0.5 km/s/Mpc (≈1.2σ). Esto es
el clásico **split BAO–CMB de los modelos w₀wₐ**, no un error.

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
| V-L2-10 | m_φ | Σm_ν^act·(Ω⁴+AURA·KAL₀) — forward-pred | 36.9463 eV | ✓ | verificado (dim.) |
| V-L2-11 | k_fs | free-streaming de m_φ (output CLASS) | 0.659 h/Mpc | — | pendiente L3 |
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
  **Pregunta abierta CERRADA (2026-06-10):** se probó la hipótesis "H_alg es la
  tasa física del fondo (la sábana desnuda) y MIRA es el observable curvado por
  la materia". REFUTADA por dos vías independientes: (1) metiendo H_alg=67.962
  como tasa de expansión real en el CMB da θ*=7.80σ y r_d=5.35σ fuera de Planck
  → H_alg NO puede ser la tasa de fondo. (2) El signo está invertido: el límite
  de Sitter (sábana solo-DE) da H_dS=H_MIRA·√Ω_DE=61.44 < H_MIRA, es decir
  quitar materia BAJA H (no lo sube a 67.96); Friedmann ata H₀ al contenido
  total, no hay "H de sábana vacía" separado a z=0. **Veredicto: la dualidad
  H_alg≈SH0ES es COINCIDENCIA** de un número adimensional en la ventana de SH0ES.
  La veta real abierta NO es el frenado por materia sino el origen dimensional
  de H₀ (escala Mpc ↔ Planck, roadmap #1).
- **V-L2-10 m_φ = Σm_ν^active · (Ω⁴+AURA·KAL₀)** — forward-prediction canónica:
  `[eV]·(número puro)=[eV]`, dimensionalmente **consistente**. Con
  Σm_ν^active = (Ω/(KAL₀·T_r))·0.960318 eV = 0.069023 eV y multiplicador
  Ω⁴+AURA·KAL₀ = 535.2795 → m_φ = 36.9463 eV (cero fiteo). Reemplaza la vieja
  cadena numerológica `Σm_ν·H₀^alg = 5.60 eV` (RETIRADA — sí era `[eV]·[km/s/Mpc]`).
  Lo que queda **ABIERTO** es el Lagrangiano φ-DM que justifique el multiplicador
  (OP-9), no la dimensión.

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
- **k_fs** (V-L2-11): depende de m_φ (canónico 36.95 eV) y del cómputo CLASS (L3).

**Estado Capa 2:** 11/13 `verificado` (numérica + dimensional + identidades);
1 `ABIERTO` (H₀^alg — adimensional vs km/s/Mpc); 1 `pendiente L3` (k_fs).
m_φ pasa a `verificado (dim.)` tras la cadena forward-prediction canónica; su
derivación del multiplicador alcanza la Capa 3 (OP-9). Ninguno pasa a
`resuelto` todavía: la comprobación 3 (derivación) de varios alcanza la Capa 3.

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

1. **✓** La jerarquía es real: (H₀/M)² ≈ 2.3×10⁻⁶² (M=9.62 meV ≫ H₀).
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

*Claim CLAUDE.md (canónico 2026-06):* titular two-sector S₈_eff=0.766 (0.01σ KiDS).

1. **✓ definición:** S₈ = σ₈(Ω_m/0.3)^½ con Ω_m,CMB=0.31993 (√(Ω_m/0.3)=1.0327).
2. **✓ single-sector (el desafío):** σ₈=0.820 → S₈=0.847 — **3.9σ KiDS**.
   Es el baseline que el modelo debe resolver.
3. **✓ two-sector φ-DM (TITULAR, forward):** el free-streaming en k_fs=0.659
   h/Mpc (de m_φ=36.95 eV, cero fiteo) baja σ₈_eff a 0.742 → **S₈_eff=0.766
   = 0.01σ KiDS-1000**. RESUELVE la tensión S₈, sin parámetros libres.
4. **○ refinamiento no-lineal (Nivel 2, diferido):** el cierre no-lineal pleno
   con feedback bariónico (N-body SSEE, ~5k–20k CPU-h) queda pendiente; HMcode-2020
   da una corrección ~0.4% (B_σ₈≈0.996). No altera el resultado lineal forward.

**Veredicto:** la tensión S₈ la **resuelve el two-sector lineal forward**
(0.766, 0.01σ). Las ramas viejas σ₈=0.737/0.794 → S₈=0.761/0.820 (HMcode,
internamente inconsistentes) y 0.702/0.725 (G=0.866, Ω_m,dyn) están **retiradas**.
Sólo el refinamiento no-lineal Nivel 2 queda ABIERTO.

## V-L3-OP6 — forma de screening f_screen / universo separado — **PARCIAL (forma derivada, valor con insumo)**

*Claim CLAUDE.md:* "OP-6 ✅ RESUELTO — universo separado k-essence + identidad 1+w₀=Ω_m".

1. **✓ valor:** f_screen = α_K/(3·MIRA) = (π−φ)/Ω² = 0.067253 — álgebra
   exacta, ya verificada en V-L2-13 y en la identidad cruzada de Capa 2.
   **Canónico (espeja Paper 9, Eq. H0local + Tabla L397):**
   H₀,local = H₀^MIRA/(1−f_screen) = 67.037/0.93275 = **71.87 km/s/Mpc**
   (1.12σ SH0ES). La ruta Type-P H₀^alg/(1−f_screen) = 72.86 km/s/Mpc es
   **comparación** (Paper 9 Tabla L398), NO el valor canónico.
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

## V-L3-mphi — masa del campo φ-DM, m_φ = 36.95 eV — **PARCIAL (cadena dim. consistente, Lagrangiano abierto OP-9)**

*Claim CLAUDE.md (canónico 2026-06-04):* "m_φ = Σm_ν^active × (Ω⁴+AURA·KAL₀)
= 36.95 eV — forward-prediction, cero fiteo".

1. **✓ numérico:** la cadena cierra — Σm_ν^active = (Ω/(KAL₀·T_r))·0.960318 eV
   = 0.071875·0.960318 = 0.069023 eV; multiplicador Ω⁴+AURA·KAL₀ = 535.2795;
   m_φ = 0.069023·535.2795 = 36.9463 eV.
2. **✓ dimensional:** `[eV]·(número puro)=[eV]`. El multiplicador es
   combinación de constantes adimensionales (Ω, AURA, KAL₀). Reemplaza la vieja
   cadena `Σm_ν·H₀^alg` (RETIRADA — esa sí era `[eV]·[km/s/Mpc]`).
3. **✓ cero entero pelado:** la razón Σm_ν usa R₂=Ω/(KAL₀·T_r), cociente puro
   de constantes — ya no la resta `4·KAL₀−22` de la versión numerológica.
4. **✗ Lagrangiano no cerrado:** falta la acción P(X,φ) que produzca el
   multiplicador Ω⁴+AURA·KAL₀ desde primeros principios. Es forward-prediction
   estructural, no derivación de mecanismo (**OP-9 ABIERTO**).

**Veredicto:** la cadena es dimensionalmente consistente y sin fiteo (cero
parámetros libres), pero su justificación desde un Lagrangiano sigue abierta
(OP-9). **PARCIAL** (era ABIERTO bajo la cadena 5.60 eV retirada).

## V-L3-2sec — modelo dos sectores φ-DM — **PARCIAL (identidad sí, split físico no)**

*Claim CLAUDE.md:* "Ω_total (dos sectores) = 0.319928 ≈ Ω_m,CMB — unificación algebraica".

1. **✓ identidad:** Ω_CDM + Ω_φDM = Ω_m,dyn + (MIRA−1)·Ω_m,dyn =
   MIRA·Ω_m,dyn = 0.319928. Diferencia con V-L2-05 = 0 exacto. Es una
   **re-partición algebraica** de Ω_m,cosm en dos mitades casi iguales.
2. **⚠ split físico:** que un sector (φ-DM) free-streame para k>k_fs y el
   otro (CDM) no, descansa en m_φ=36.95 eV (canónico, cadena dim. consistente
   — V-L3-mphi PARCIAL) y en k_fs=0.659 h/Mpc (output CLASS). El cierre del
   Lagrangiano que justifique el split sigue en OP-9.

**Veredicto:** la suma Ω_total es un re-enunciado exacto de V-L2-05; el
modelo físico de dos sectores hereda la apertura del Lagrangiano de m_φ
(OP-9). **PARCIAL.**

## V-L3-EFT — acción EFT canónica (Paper 7) — **PARCIAL (parámetros sí, M⁴ inconsistente)**

*Claim CLAUDE.md:* "Paper 7: EFT canónico — λ/V₀/M/g² bloqueados".

1. **✓ λ, α_pot, V₀:** son consecuencias algebraicas exactas de constantes
   ya verificadas — λ²=3·Ω_m,dyn (λ=0.6929), α_pot=λ/√KAL₀ (=0.2949),
   V₀=Ω_DE·ρ_crit (=0.8400). No son parámetros libres: re-enuncian
   Ω_m,dyn (V-L2-04), KAL₀ (V-L1-07) y Ω_DE (V-L2-03).
2. **✗ M⁴ inconsistente entre papers:** `ssee_eft_verification.py` (L70)
   fija **M⁴ = ρ_crit = 1**; `ssee_paper10_verification.py` (L30) fija
   **M⁴ = 5φ⁸·ρ_crit = 234.9**. Factor ~235 de diferencia en el mismo
   término X²/M⁴ del mismo Lagrangiano K(X). Con M⁴=1 el término UV no es
   perturbativo; con M⁴=234.9 sí. Es una incoherencia cruzada P7↔P10.

**Veredicto:** los parámetros λ, α_pot, V₀ están algebraicamente fijados;
M⁴ tiene dos valores incompatibles según el paper. **PARCIAL.**

## V-L3-KX — completación UV K(X) (Paper 10) — **ABIERTO (M⁴ calibrado a SH0ES)**

*Claim CLAUDE.md:* "Paper 10: M⁴=5φ⁸ρ_crit exacto; H₀^UV canónico=72.05 (MIRA=67.037, 0.96σ; condicional C.1). Type-P 73.040 es coincidencia numérica".

1. **✓ identidad 45α² = 5φ⁸:** exacta a precisión de máquina (α=φ⁴/3 →
   45α²=45φ⁸/9=5φ⁸=234.89). El *valor numérico* de M⁴/ρ_crit cierra.
2. **✗ normalización física no derivada:** el propio
   `ssee_paper10_verification.py` (VERDICT, L167–175) admite — *"This UV
   result is CONDITIONAL on Postulate C.1, which is calibrated to SH0ES —
   not an independent prediction"* y *"PENDING: first-principles derivation
   of M⁴ = 45α² without using SH0ES as input"*. La Ruta A da M⁴≈418 ≠ 234.9.
   M⁴=5φ⁸ es el número que hace falta para llegar a H₀=73.04, expresado en
   φ — mismo patrón que OP-1 (ajuste a objetivo conocido).
3. αK_full=0.41691 y H₀^UV (canónico 72.05 vía MIRA=67.037; Type-P 73.040) son aguas
   abajo de este M⁴ **ABIERTO**.

**Veredicto:** la forma algebraica 5φ⁸ es exacta, pero su identificación
como el cutoff UV físico está calibrada a SH0ES, no derivada. El script de
P10 ya lo declara honestamente. **ABIERTO.**

## V-L3-IS — perturbaciones Israel-Stewart (Paper 5) — **PARCIAL (c²_s,eff=0 sí, mecanismo asertado)**

*Claim CLAUDE.md:* "Paper 5: c²_s,eff = 0 (exacto algebraico) — todos los modos estables".

1. **✓ c²_s,eff = 0:** la corrección IS es ζ̃/τ_Π = (KAL₀/3)/(KAL₀/(3Ω_DE))
   = Ω_DE, y c²_s,eff = w₀ + Ω_DE = −0.8399 + 0.8399 = 0. Cierra exacto.
2. **✗ el mecanismo IS no hace trabajo:** el factor KAL₀/3 **se cancela**
   — ζ̃=KAL₀/3 y τ_Π=KAL₀/(3Ω_DE) comparten KAL₀/3, así que ζ̃/τ_Π=Ω_DE
   para *cualquier* ζ̃. El resultado se reduce a la identidad ya verificada
   w₀ = −Ω_DE (V-L2-01/03). La hipótesis ζ̃=KAL₀/3 es decorativa; la
   derivación de τ_Π (estado estacionario IS) está **asertada, no mostrada**
   en el script. Si τ_Π se deriva de verdad independientemente, el resultado
   es no-trivial; tal como está, es la tautología w₀+|w₀|=0.
3. **✓ Q2 reportado con honestidad:** el test MIRA perturbativo da
   MIRA_num=0.989 (k≥10) vs MIRA_alg=1.999 — **no coinciden**. Paper 5
   concluye correctamente que MIRA es un efecto de fondo, no perturbativo.
   Es un resultado negativo bien reportado, no un problema.

**Veredicto:** c²_s,eff=0 es cierto pero se reduce a w₀=−Ω_DE; el aparato
IS (ζ̃, τ_Π) está construido para reproducir esa identidad y su parte
no-trivial (derivación de τ_Π) no se muestra. **PARCIAL.**

## V-L3-cs2 — extracción del T_μν efectivo del sector k-essence — **ABIERTO (resultado clave)**

*Origen:* el autor pidió extraer el tensor de energía efectivo de la
acción para ver si el sector geométrico puede ser la materia que el CMB
exige (la «0.320»). Extracción hecha 2026-05-22.

**La extracción (k-essence estándar, Garriga–Mukhanov 1999).** Para la
acción de SSEE P(X,φ) = X/KAL₀ + X²/M⁴ − V(φ):
- ρ_φ = 2X·P_X − P = X/KAL₀ + 3X²/M⁴ + V
- p_φ = P = X/KAL₀ + X²/M⁴ − V
- w_φ = p_φ/ρ_φ
- **c_s² = P_X/(P_X + 2X·P_XX) = (A + 2BX)/(A + 6BX)**, con A=1/KAL₀>0,
  B=1/M⁴>0.

**Resultado clave (exacto, analítico, parámetro-independiente).** c_s² es
monótona decreciente en X: vale 1 en X=0, tiende a 1/3 en X→∞. Por tanto
**c_s² ∈ [1/3, 1] para todo X≥0, cualquier M⁴>0, cualquier KAL₀>0.**
Depende solo de la *forma* K = A·X + B·X² (suma de coeficientes
positivos), no de los valores. La inconsistencia M⁴ P7↔P10 no lo afecta.

**Qué significa — respuesta a «¿la geometría tiene peso?».**
- **Peso de fondo: SÍ.** El sector tiene una densidad ρ_φ genuina; su
  w_φ de fondo puede transitar (kinético→potencial).
- **Peso de agrupamiento: NO.** Un fluido se agrupa como materia fría
  solo si c_s² ≈ 0. La k-essence de SSEE **nunca baja de 1/3** — se
  difumina, no forma pozos de potencial. No puede agruparse como CDM.

**Veredicto.** El CMB necesita materia que se *agrupe* (que forme los
pozos donde oscila el fluido fotón-barión). El sector k-essence de SSEE,
tal como está la acción, **no puede hacerlo** — su c_s² es estructuralmente
demasiado alto. **El mecanismo MIRA no está en la acción vigente.** Es un
resultado negativo limpio: descarta «la k-essence es la materia geométrica»
y dice qué hay que buscar — un sector cuyo c_s² pueda anularse a alto z
(la viscosidad IS pretendía eso, pero V-L3-IS la halló decorativa).
**ABIERTO.**

## V-L3-disf — mecanismo disformal de Paper 8 (Ruta B) — **ABIERTO (inconsistencia interna P1↔P8)**

*Origen:* tras el resultado negativo de V-L3-cs2 (la k-essence no se
agrupa), se exploró la «Ruta B»: que el 0.320 no sea sustancia sino
gravedad amplificada. Paper 8 ya invoca «MIRA emergence» vía geodésicas
disformales — se auditó si ese mecanismo deriva MIRA sin materia oscura.
Auditoría hecha 2026-05-22.

**Lo que dice Paper 8.** La acción (P8 eq.1) es
S = ∫√−g[M_pl²R/2 + P(X,φ)] + **S_DM[g̃_μν; ψ_DM]** + S_b[g_μν; ψ_b].
Contiene un **campo de materia oscura ψ_DM** acoplado a la métrica
disformal g̃_μν = g_μν + (2/M⁴)∂_μφ∂_νφ. Todo el mecanismo se sostiene
en ρ_DM: la ecuación del escalar (1/KAL)∇²φ = (β_c/M_pl)·ρ_DM (P8 L221),
el «régimen dominado por materia oscura» (L238), «la fuerza sobre una
partícula de materia oscura» (L271).

**La contradicción.** Paper 1 L51: el modelo se construye «without
introducing free parameters or exotic dark-matter particles». Paper 1
L278: SSEE «is refuted by direct detection of a collisionless
dark-matter particle». **Paper 8 presupone exactamente lo que Paper 1
declara falsador.** No es una extensión — es una violación del postulado
fundacional.

**Además, MIRA no «emerge».** El propio Paper 8 (L69) admite que
√β_c/MIRA = 1.00030 es «a near-coincidence, **not an algebraic
identity**». La «emergencia» de MIRA en el lensing es una coincidencia
numérica al 0.03 %, no una derivación.

**Veredicto.** La Ruta B tal como está escrita en Paper 8 **no deriva
MIRA sin materia oscura** — la introduce. Con V-L3-cs2 (sustancia que se
agrupa: descartada) y α_B=α_M=0 de Paper 7 (Poisson modificada μ>1:
imposible, μ=1 exacto), **los tres mecanismos posibles para el 0.320 en
un marco sin materia oscura fallan.** El 0.320/MIRA no está derivado en
ninguno de los 10 papers.

**Lo rescatable.** La ecuación (1/KAL)∇²φ = (β_c/M_pl)·ρ, con ρ → materia
bariónica + dinámica real (0.160) en vez de ρ_DM, sería un modified
gravity genuino sin materia oscura: el campo geométrico φ desarrolla
perfiles alrededor de materia ordinaria y amplifica su gravedad. Pero
(a) Paper 8 da G_eff/G ≈ 177 en cúmulos, no ≈2=MIRA en el fondo
cosmológico — régimen distinto, sin verificar; (b) amplificar el
agrupamiento no mueve las posiciones de los picos del CMB, fijadas por
r_d/D_A del fondo. Rescate posible en principio, no hecho. **ABIERTO.**

**Bloquea:** sellado de Paper 8 (inconsistencia con Paper 1); y deja a
P3/P6/P9 sin un mecanismo válido para el 0.320.

## V-L3-mira — test del mecanismo de retención conformal para MIRA — **ABIERTO (mecanismo candidato FALLA)**

*Origen:* tras V-L3-cs2 y V-L3-disf se construyó y probó la «Ruta B» en su
forma más concreta — MIRA como retención de energía vía un acoplamiento
conformal de quintaesencia entre φ y la materia real (el L_int de Paper 7
re-apuntado de materia oscura a materia real, como exige Paper 1).
Cálculo: `src/mira_attempts/ssee_mira_mechanism.py`, 2026-05-22.

**El mecanismo (física estándar, quintaesencia acoplada Amendola 2000).**
Un acoplamiento conformal da, en el fondo:
- ρ_m = ρ_m,0 a⁻³·exp[−β_c(φ−φ_0)]  ⟹  R(a) = exp[−β_c·Δφ(a)]
- KG con fuente: (P_X+2X P_XX)φ̈ + 3H P_X φ̇ + V_φ = β_c ρ_m
- Tasa de retención Γ = |β_c φ̇|; desacople en Γ = H.
Para R_temprano = MIRA se requiere AURA·Δφ_total = ln(MIRA) = 0.693, i.e.
una excursión del campo Δφ ≈ 0.173 M_pl.

**El test (β_c = −AURA fijo, sin ajustar nada).** Integración del fondo
acoplado Friedmann+KG hacia atrás desde hoy hasta z=1100:

| | Necesario | Obtenido |
|---|---|---|
| R(hoy) | 1 | 1.00000 ✓ (maquinaria OK) |
| R(z=1100) | MIRA ≈ 1.999 | ≈ 0 ✗ |
| AURA·Δφ_total | +0.693 | **−12.3** ✗ |
| signo | R>1 (carga materia) | R<1 (drena materia) ✗ |
| timing | retención activa temprano | activa tarde (Γ/H≈3 hoy) ✗ |

**Veredicto — negativo limpio.** β_c = AURA ≈ 4 es un acoplamiento
cosmológico enorme (los límites realistas de quintaesencia acoplada dan
β ≲ 0.1; AURA es ~40× ese techo). Con esa fuerza el término β_c ρ_m de la
KG **golpea el campo** y lo hace rodar Δφ ≈ −3.1 (×18 lo necesario) con
signo invertido: en vez de cargar el sector materia lo **drena** —
R(z=2) ≈ 0.016, la materia casi desaparece. Además la KG empuja la
k-essence a un régimen sin solución real de Friedmann en el 44 % de la
trayectoria. Tres fallas independientes (magnitud, signo, timing) ⟹ no es
«casi» — es el mecanismo equivocado.

**Estado del problema MIRA.** Cuatro mecanismos probados para el «0.320»,
cuatro negativos: sustancia que se agrupa (V-L3-cs2), Poisson modificada
μ>1 (α_B=α_M=0, Paper 7), disformal (V-L3-disf, exige ψ_DM), y retención
conformal (esta entrada). **MIRA no tiene derivación en el marco vigente
de SSEE por ninguno de los cuatro mecanismos naturales.** El «0.320» de
P3/P6/P8/P9 está insertado, no derivado. El núcleo w₀wₐ (DESI 0.5σ) no se
ve afectado. Ver [[feedback-impact-analysis]] y V-L3-2Om. **ABIERTO.**

## V-L3-saturacion — α por saturación φ-MDE (veta-2, Amendola-rescaled) — **PARCIAL: álgebra OK, dinámica FALLA**

*Origen:* tras V-L3-mira (β_c=−AURA falla con tres modos independientes), se
abrió la búsqueda de un **principio físico** que fije α en el L_int conformal
sin recurrir a fits ni a numerología. Hipótesis: α toma el valor *saturado*
de una desigualdad física estándar (veta-2 del programa de reconstrucción).
Cálculo: `src/mira_attempts/ssee_alpha_saturation_stepD2.py`, 2026-05-22.

**La desigualdad física (Amendola 2000, canonical).** En quintaesencia
acoplada con K=X y conformal coupling, el φ-MDE (fixed-point y=0 durante era
de materia) está en x = −α√6/3 y tiene Ω_φ = 2α²/3. La condición Ω_φ ≤ 1
(la era de materia debe existir, *the universe must pass through matter
domination*) impone:

> α² ≤ 3/2 → α_sat,canonical = √(3/2) ≈ 1.2247

**Traducción a SSEE.** Con K(X) = X/KAL + X²/M⁴, en el límite σ→0 (régimen
X/KAL dominante, baja energía / σ = H²M_pl²/M⁴ ≪ 1), la renormalización
natural del campo χ = φ/√KAL traduce la cota canónica a:

> **α_sat = √(3/(2·KAL₀)) = √(3/(φ+3π)) ≈ 0.5212**

Equivalente: **α² · 2·KAL₀ = 3**, o (α·√KAL₀)² = 3/2. Identidad estructural,
no fit.

**Verificación numérica.**

| | Esperado (analítico) | Obtenido (numérico) | Desvío |
|---|---|---|---|
| α_sat canónico (σ=0, KAL=1) | √(3/2) = 1.22474487 | 1.22474487 | 1.4·10⁻⁷ % |
| x_φMDE canónico (α=1.2) | −α√6/3 = −0.97980 | −0.97980 | exacto |
| α_sat SSEE (σ→0, KAL=KAL₀) | √(3/(φ+3π)) = 0.521220 | 0.520674 | 0.10 % |

**Las 6 comprobaciones:**

1. **Numérica:** OK, 0.10 % de error al límite σ→0 (la diferencia es
   suppresión por σ finito; en σ=10⁻⁹ ya converge).
2. **Dimensional:** α adimensional (sale en exp(α·φ/M_pl)); 2·KAL₀ adimensional;
   √(3/(2·KAL₀)) adimensional. OK.
3. **Derivación:** Amendola 2000 (literatura física estándar, no ad-hoc),
   aplicado a K(X) de Paper 7 con normalización canónica del campo.
   *No es post-hoc, no es fit.*
4. **Rol:** fija el coeficiente del L_int conformal en el sector φ-materia.
   *Candidato* para mecanismo MIRA — pendiente test dinámico.
5. **Ubicación:** *no aplicado todavía* en ningún paper. Esta entrada lo registra
   como candidato derivado, no como resultado de un paper.
6. **Conexiones:** depende de KAL₀ (V-L1) y de la estructura k-essence de Paper 7.
   Sería insumo de V-L3-mira si reemplaza a β_c=−AURA.

**Test dinámico (2026-05-22).** `src/mira_attempts/ssee_mira_saturated_diagnostic.py`:
matriz 2×2 con β_c = ±α_sat y u₀ = ±|u_w₀|, integrando hacia atrás desde el
estado de hoy (fijado por w_φ(hoy)=w₀=−0.840) hasta z=1100, régimen lineal
puro K=X/KAL₀ (donde la derivación aplica):

| β_c | u₀ signo | resultado |
|---|---|---|
| +α_sat | + | integración rompe ~z<1 (runaway u) |
| +α_sat | − | rompe |
| −α_sat | + | rompe |
| −α_sat | − | rompe |

Las cuatro fallan. **No existe trayectoria suave** desde el atractor tardío
(estado de hoy, u₀≈+1.5) hasta el φ-MDE (u≈−2α_sat≈−1.04) con esta familia
de modelos y esta normalización. Cualquier signo se vuelve runaway antes de
cruzar la transición de régimen.

**Interpretación física.** La cota Amendola garantiza que el φ-MDE *exista*
como FP cuando α² < 3/(2KAL₀), pero **no garantiza** que la trayectoria
backward desde el atractor tardío de DE alcance ese FP suavemente. Son dos
preguntas distintas: (a) ¿hay FP físico? — sí. (b) ¿el atractor tardío
desciende a él? — no, al menos no con la α y la inicialización canónicas
del modelo. La veta-2 produjo un número honesto, pero la dinámica acoplada
en SSEE no lo aprovecha.

**Lo que sí queda como producto positivo:**

- **Primer número derivado por saturación** en SSEE: α_sat=√(3/(φ+3π))
  cumple las 6 comprobaciones de la entrada y no es numerología.
- **Confirmación estructural:** Amendola 2000 generaliza limpiamente a
  SSEE bajo la normalización inducida por KAL₀ — el método funciona, lo
  que falla es la aplicación al mecanismo MIRA específico.
- **Evidencia adicional** para V-L3-2Om: ni el coupling AURA-grande
  (V-L3-mira) ni la cota Amendola-saturada conectan el background tardío
  con z=1100. Quinto mecanismo descartado para MIRA.

**Estado.** Álgebra **verificada**. Dinámica acoplada **falla** numéricamente
en las cuatro variantes de signo — el bound es real pero no derivado en una
forma constructiva del mecanismo MIRA. **PARCIAL.**

### Extensión #6+#7 (2026-05-22) — derivativo + forward desde φ-MDE

Tras V-L3-saturacion, se cazaron dos mecanismos adicionales:

**Mecanismo #6 — acoplamiento derivativo** `L_int=(X/M⁴)·L_DM`.
Script `src/mira_attempts/ssee_mira_derivative_test.py`. Coupling β_eff=u/√M⁴ depende
de la *velocidad* del campo (no del valor), encendido en matter era y
apagado hoy (lo que MIRA "necesita"). Barrido en M⁴_code:
- M⁴ ∈ {0.01, 1, 10}: integración OK, **R(z=1100)∈[1.32, 1.41]** — ~70%
  del log de MIRA, no llega.
- M⁴ ∈ {100, 462 (físico M≈9.62 meV canónico; el 8.81 meV era normalización ρ_crit previa), 10⁴}: integración rompe en z≈1–2.

Diagnóstico: el mecanismo funciona como mecanismo, pero requiere M
muchísimo más chico que Λ_SSEE para producir MIRA — UV-incompatible.

**Mecanismo #7 — integración forward desde φ-MDE hacia hoy.** Script
`src/mira_attempts/ssee_mira_phimde_forward_test.py`. Estado inicial Ω_m≈0.99 en
matter era, integrar hasta N=0. Barrido α∈{0.295, 0.521, 0.8, 1.0, 1.5}:

| α | w_φ(0) | Ω_m(0) | resultado |
|---|---|---|---|
| 0.295 (Copeland exp puro) | −0.88 ✓ | **0.00** ✗ | drenaje total |
| 0.521 (α_sat veta-2)     | −0.70   | **0.00** ✗ | drenaje |
| 0.8–1.5                  | crece   | **0.00** ✗ | drenaje |

**Hallazgo central:** las trayectorias forward llegan a N=0, pero el
coupling conformal drena *toda* la materia hacia el scalar en ~7 e-folds.
Ningún α reproduce simultáneamente w₀=−0.840 y Ω_m=0.32.

**Cierre matrix backward+forward.** Combinando #5 (backward) y #7 (forward):
no existe coupling conformal canónico que conecte φ-MDE con el estado de
hoy preservando materia. **6 mecanismos descartados para MIRA**, dos
cualitativamente distintos (topología #5/#7, escala #6).

Las dos lecturas siguen abiertas:
- **A**: continuar caza (disformal velocity-dependent, screening, ...)
- **B**: pivot a MIRA-como-etiqueta de transición, no engranaje.

## V-L3-2Om — mecanismo MIRA / transición Ω_m(z) — **ABIERTO (problema central del modelo)**

*Corrección de criterio (2026-05-22):* una versión anterior de esta entrada
marcaba esto «verificado (regla)». **Era una sobreafirmación.** Lo único
verificado es el *álgebra*; la *regla física* de uso no está derivada.

1. **✓ álgebra:** Ω_m,dyn = 1+w₀ = 0.16005 **sí** se deriva de φ,π.
   Ω_m,cosm = MIRA·Ω_m,dyn = 0.31993 — el número es algebraico, pero
   **MIRA es hipótesis auxiliar no derivada** (lo dice `ssee_core.py` L41).
2. **✗ la regla NO está derivada.** «Ω_m,dyn fija w₀ y no entra en E(z);
   Ω_m,cosm va en E(z)/Poisson/CMB» es una **aserción**, no un teorema.
   No existe:
   - una derivación de *por qué* MIRA (≈2) amplifica la materia gravitante;
   - una función Ω_m,eff(z) que conecte el régimen tardío (BAO: 0.160
     funciona, 0.5σ DESI, cero parámetros libres) con el temprano (CMB:
     0.320 encaja). El modelo usa 0.320 y 0.160 como dos regímenes sin
     puente — **no se sabe qué pasa en la transición**.
3. **✗ riesgo de fondo:** poner 0.320 como término ∝(1+z)³ en E(z) es,
   fenomenológicamente, **materia oscura** — y el postulado fundacional de
   SSEE es *sin materia oscura, solo geometría y viscosidad*. Usar 0.320 en
   el fondo sin un mecanismo derivado puede estar contradiciendo el axioma
   que el modelo dice no necesitar.

**Por qué importa.** El resultado más fuerte de SSEE (DESI 0.5σ, cero
parámetros) se construyó con 0.160 en el fondo. El commit `4892b53`
(«corrección dos-Ω_m») cambió E(z) a 0.320 y, con el prior Planck-ΛCDM legacy,
movió H₀ a 67.756 — que producía r_d/θ* a 4–5σ. **El switch de prior MIRA
(2026-05-24) reancló H₀ a 66.533 posterior / 67.037 anchor**, y con esos H₀
r_d queda a ≤1.5σ y θ* a 0.66σ *en el anchor* (ver V-L4-rd/θ* re-corridas
2026-06-09). Lo que queda abierto NO es ya un 5σ duro en r_d, sino: (a) la
ausencia de Ω_m,eff(z) que puentee 0.160↔0.320, y (b) el split BAO–CMB de
1.2σ en H₀, que se manifiesta como tensión en θ* solo si se usa el posterior
en el observable CMB. Sigue siendo un cambio de configuración con
consecuencias, no una «corrección» sin costo.

**Veredicto.** Qué es MIRA físicamente, por qué la materia gravitante es 2×
la dinámica sin partícula de materia oscura, qué papel cumple AURA (=2·MIRA),
y cuál es Ω_m(z) entre CMB y BAO — **ése es el problema central abierto del
modelo**. Todo lo demás (r_d, θ*, H₀) son síntomas de esta brecha. **ABIERTO.**

## Estado final de Capa 3

Re-verificados **15 elementos**: OP-1..OP-7, α=φ⁴/3, m_φ, dos sectores,
EFT, K(X), IS, c_s² (T_μν), dos-Ω_m.

| Veredicto | Elementos |
|---|---|
| **verificado** | α=φ⁴/3 |
| **PARCIAL** (álgebra/forma cierra, insumo físico no) | OP-2, OP-6, OP-7, dos sectores, EFT, IS, m_φ (cadena dim. consistente; Lagrangiano OP-9) |
| **ABIERTO** | OP-1, OP-3, OP-4, OP-5, K(X), c_s² (T_μν), **dos-Ω_m (central)** |

2 bugs corregidos/detectados de paso: curvatura de Kähler (P1, **corregido**)
y M⁴ inconsistente P7↔P10 (**detectado**, pendiente de re-derivación).
Patrón confirmado en 14/14: ninguna "✅ RESUELTO" de CLAUDE.md lo estaba
sin reservas. **Ninguno es regresión** — todos son brechas preexistentes,
ahora rastreadas por el guardián (67 comprobaciones, 15 ABIERTO).

# Capa 4 — Confrontaciones con datos — *en progreso*

Aquí se verifica la **aritmética** que conecta cantidades reportadas:
tensiones (model−obs)/σ, S₈=σ₈√(Ω_m/0.3), χ²_r=χ²/N, ΔBIC. Los χ² de CMB
y los posteriores MCMC en sí salen de CAMB/CLASS/emcee — el guardián
verifica que los números encajen entre ellos, no re-corre los pipelines.

## V-L4-S8 — cadena S₈ weak-lensing (canónico 2026-06) — **verificado (aritmética)**

Usa Ω_m,cosm=0.31993 → √(Ω_m,cosm/0.3)=1.03268 (S₈ es amplitud gravitacional).

1. **✓ single-sector ("el desafío"):** σ₈ = σ₈_LCDM·G = 0.811·1.011 = 0.820;
   G=1.011 es el ODE de crecimiento con fuente Poisson Ω_m,CMB=0.320.
   S₈ = 0.820·1.03268 = 0.847 → **3.9σ KiDS-1000** (DES-Y3 ≈ 3.9σ).
2. **✓ two-sector φ-DM (TITULAR, forward):** σ₈_eff = 0.742 (free-streaming
   CLASS, k_fs=0.659 de m_φ=36.95 eV, cero fiteo). S₈_eff = 0.742·1.03268 =
   **0.766 → 0.01σ KiDS-1000**. RESUELVE la tensión.

**Veredicto:** la cadena S₈ es aritméticamente correcta y usa la Ω_m correcta.
El titular es el two-sector (0.766, 0.01σ). La cadena vieja G=0.866 →
σ₈=0.7023 → S₈=0.7253 (fuente Ω_m,dyn) está **retirada**. **Verificado.**

## V-L4-DES — referencia DES-Y3 inconsistente entre scripts — **ABIERTO**

`ssee_paper5_IS_perturbations.py` usa **S₈_DES = 0.776±0.017** (DES-Y3
3×2pt, Abbott et al. 2022). `ssee_op5_hmcode.py` usa **S₈_DES = 0.759±0.023**
(DES-Y3 cosmic shear, Amon et al. 2022). Son dos análisis DES-Y3 reales y
distintos, pero la suite debería fijar **una** referencia para que las
tensiones reportadas sean comparables entre papers. Un árbitro lo marcaría.

## V-L4-CMB — espectro CMB Planck PR4 (Paper 3) — **verificado (re-run CAMB 2026-05-22)**

Re-corrido `ssee_paper3_cmb.py` con CAMB 1.6.5. Reproduce **exactamente**
lo reportado en CLAUDE.md:

| Espectro | SSEE χ²_r | ΛCDM χ²_r | N |
|---|---|---|---|
| TT | 1.047 | 1.043 | 1971 |
| TE | 1.041 | 1.040 | 1967 |
| EE | 1.041 | 1.039 | 1967 |
| PP | 0.837 | 0.757 | 9 |
| **ΔBIC(SSEE−ΛCDM)** | **−20.8** | — | N=5914 |

Picos TT en ℓ = 220, 536, 812 — también reproducidos. La aritmética
χ²_r→ΔBIC solo acota a [−22.9, −11.1] por el redondeo de χ²_r a 3
decimales; el −20.8 del pipeline cae dentro. **El test de datos central
de Paper 3 es reproducible. Verificado.**

## V-L4-rd — horizonte de sonido r_d — **r_d sano en ambos anclajes**

Re-corrido con CAMB (2026-06-01) tras el switch de prior MIRA. El "4.47σ"
previo era **artefacto del H₀ obsoleto 67.756** (prior Planck-ΛCDM legacy,
superado). Con los H₀ canónicos actuales:

| H₀ usado | r_d resultante | tensión Planck (147.09±0.26) |
|---|---|---|
| 66.533 (posterior MCMC, prior MIRA, mnu=0.069) | 147.30 Mpc | **0.80σ ✓** |
| 67.037 (anchor H_MIRA, CMB-óptimo, mnu=0.069) | 146.73 Mpc | **1.38σ ✓** |
| 67.756 (obsoleto, prior Planck legacy) | 145.93 Mpc | 4.47σ — *superado* |

**Mecanismo:** la parametrización SSEE fija Ω_m,cosm=0.31993 (fracción) y
deriva ω_c = Ω_m,cosm·h² − ω_b. Al subir H₀, sube h², sube ω_c, la igualdad
materia-radiación se adelanta y r_d **cae**. r_d es genuinamente sensible a
H₀. Con el H₀ real del modelo (66.533 posterior / 67.037 anchor), r_d está a
≤1.5σ de Planck. **El "✅" de r_d se restituye.**

## V-L4-θ* — escala acústica angular θ* — **ABIERTO (sensibilidad extrema a H₀)**

θ* es el observable CMB *más preciso* (σ≈0.08%) y por eso es hipersensible a
H₀ vía D_A. Re-run CAMB 2026-06-01, vs Planck 2018 **0.59668±0.00046°**:

| H₀ usado | θ* resultante | tensión |
|---|---|---|
| 67.037 (anchor H_MIRA, CMB-óptimo, mnu=0.069) | 0.59638° | **0.66σ ✓** |
| 66.533 (posterior MCMC, corregido por DESI, mnu=0.069) | 0.59423° | **5.33σ** |
| 67.756 (obsoleto) | 0.59927° | 5.62σ — *superado* |

**Crux honesto:** en el **anchor CMB-óptimo (67.037)**, θ* (0.66σ) y r_d
(1.38σ) están **ambos sanos**. La tensión de 5.33σ aparece SOLO si se inyecta
el H₀ posterior corregido-por-DESI (66.533) en el observable CMB — es decir,
es la manifestación en θ* del **split BAO–CMB de 1.2σ en H₀** de los modelos
w₀wₐ. Un control ΛCDM (w=−1) al mismo H₀ da θ* casi idéntico: **el tirón en
θ* lo causa H₀, NO la energía oscura w₀wₐ de SSEE.**

**Lectura para el documento de journal:** anclar el CMB en H₀=67.037 (donde
r_d y θ* ≤1.5σ), reportar el posterior BAO 66.533±0.442, y declarar el split
BAO–CMB (1.2σ) como feature conocido de w₀wₐ — no como falla. Ningún H₀ único
satisface BAO-posterior + r_d + θ* todos a <2σ a la vez, pero 67.037 deja los
tres a ≤1.5σ.

## V-L4-MCMC — MCMC DESI+Planck (Paper 2) — **re-run 2026-05-22; H₀ derivó**

Re-corrido `ssee_paper2_mcmc.py` (100 walkers × 25000 pasos × 3 modelos,
1.52 h). Posteriores:

| Modelo | k | H₀ | ln P_MAP | BIC | ΔBIC |
|---|---|---|---|---|---|
| **SSEE** | 2 | 67.756 ± 0.442 | −13.22 | 31.98 | **0.00** |
| ΛCDM | 3 | 68.283 ± 0.380 | −15.79 | 39.89 | +7.91 |
| CPL | 5 | 67.301 ± 0.523 | −11.74 | 37.35 | +5.37 |

1. **✓ aritmética BIC:** BIC = k·ln(16) − 2·lnP_MAP se recomputa exacto —
   SSEE 31.98, ΛCDM 39.89, **ΔBIC=+7.91 a favor de SSEE**. La tensión
   H₀ vs Planck (67.36±0.54) es 0.57σ — también recomputada y cierra.
2. **✗ H₀ derivó:** el MCMC da **H₀ = 67.76 ± 0.44**; CLAUDE.md registra
   **66.75 ± 0.44**. Deriva ~1 km/s/Mpc — más de 2σ del propio ancho del
   posterior. El valor registrado está **obsoleto**.
3. **✗ Ω_b h² tira más bajo:** el posterior da Ω_b h²=0.02183±0.00048;
   la fórmula algebraica OP-1 da 0.02242. El dato prefiere ~1.2σ **menos**
   barión que (π−φ)/(3Ω²). Es una señal en contra de la coincidencia OP-1.
4. r_d(SSEE)=175.16 Mpc en el MCMC — es el horizonte *crudo* (pre-mapeo
   MIRA), consistente con la narrativa r_d,SSEE≈175 → r_d,eff≈147.
5. El MCMC reporta "Ω_m tensión SSEE vs Planck: 21.3σ" — es la huella
   estructural dos-Ω_m (Ω_m,dyn=0.160 vs Planck 0.315), ya documentada en
   V-L3-2Om, no un bug nuevo.

**Veredicto:** la aritmética estadística (BIC, ΔBIC, tensión H₀) cierra y
**SSEE sigue favorecido (ΔBIC=+7.91)**. Pero el headline H₀ registrado
(66.75) está obsoleto — el valor vivo es 67.76. Y el dato bariónico
empuja en contra de OP-1.

## Estado final de Capa 4

Re-corridos los tres pipelines (CAMB r_d, CAMB CMB, emcee MCMC) el
2026-05-22.

| Veredicto | Elemento |
|---|---|
| **verificado / reproducido** | S₈ (P5), χ²+ΔBIC del CMB (P3), BIC+ΔBIC del MCMC (P2) |
| **ABIERTO — tensión grave (enmascarada)** | r_d 4.47σ, θ* 5.62σ — al usar el H₀ canónico |
| **ABIERTO — deriva de valor (resuelta)** | H₀ MCMC 66.75→67.76 — re-anclado a 67.756 |
| **ABIERTO — tensión física** | Ω_b h² −1.2σ vs OP-1 |
| **ABIERTO — inconsistencia de referencia** | DES-Y3 (0.776 vs 0.759) |

**Lo que reprodujo es sólido**: el ajuste al CMB (χ²_r) y la preferencia
estadística por SSEE (ΔBIC negativo en P2 y P3) se sostienen al re-correr.

**Lo que la campaña destapó** — y es lo más importante de toda Capa 4: el
valor stale H₀=66.75 estaba **enmascarando** una tensión de 4.47σ en r_d
y 5.62σ en θ*. Al re-anclar H₀ al canónico 67.756 y propagarlo, las
tensiones aparecieron. El "r_d ✅ 0.25σ" era un artefacto. **Esta es la
brecha más grave del modelo y bloquea el sellado de Papers 2, 3 y 9**
hasta que se encare (revisar la parametrización Ω_m,cosm, o aceptar la
tensión y reportarla). Pendiente menor: fσ₈ (Papers 5–6).

# Capa 5 — Sellado de papers — *pendiente*

| Paper | Todos los elementos re-verificados | Chequeo completo | Sello (commit + sha256) |
|-------|-----------------------------------|------------------|-------------------------|
| P1–P10, Unified | — | — | — |

---

*Registro iniciado 2026-05-21 tras la revisión árbitro hostil de los 11 documentos.*
