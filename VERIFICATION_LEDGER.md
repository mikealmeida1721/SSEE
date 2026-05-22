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
| **L2** | Parámetros cosmológicos derivados | pendiente |
| **L3** | Mecanismos y derivaciones (OP-1..OP-7, EFT, IS, dos-sectores, k-mouflage, f_screen, m_φ, K(X) UV) | pendiente |
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

# Capa 3 — Mecanismos y derivaciones — *pendiente*

Incluye la re-verificación de OP-1..OP-7 (todas marcadas `✅ RESUELTO` en
CLAUDE.md — **a tratar como sospechosas hasta verificarlas**), la acción EFT, las
perturbaciones Israel-Stewart, el modelo de dos sectores φ-DM, el screening
k-mouflage / Vainshtein (fórmula rota confirmada), f_screen, la masa m_φ, y la
completación UV K(X).

# Capa 4 — Confrontaciones con datos — *pendiente*

MCMC DESI+Planck, espectros CMB, ΔBIC, fσ₈, S₈, tensión de Hubble.

# Capa 5 — Sellado de papers — *pendiente*

| Paper | Todos los elementos re-verificados | Chequeo completo | Sello (commit + sha256) |
|-------|-----------------------------------|------------------|-------------------------|
| P1–P10, Unified | — | — | — |

---

*Registro iniciado 2026-05-21 tras la revisión árbitro hostil de los 11 documentos.*
