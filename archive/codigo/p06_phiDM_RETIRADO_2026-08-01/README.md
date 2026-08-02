# Cajón φ-DM — RETIRADO el 2026-08-01

Aquí vive el código del **sector φ-DM de dos componentes y su partícula**
(`m_φ = 40.70 eV`, `k_fs = 0.754 h/Mpc`, `α = 1.117 Mpc/h`, `Ω_φDM = 0.14889`,
`σ₈ = 0.747`, `S₈ = 0.758`). **Nada de esto es vigente.** Se conserva como
registro de una hipótesis que se construyó, se probó y se retiró — no como
código a reutilizar.

## Por qué cayó — dos golpes independientes, cada uno suficiente

**1. La resta que definía la partícula no era física.**

    Ω_φDM = Ω_m,CMB − Ω_m,dyn = 0.308881 − 0.160050

El primer término es una densidad medida. El segundo es `1+w₀`: menos la
ecuación de estado de la energía oscura, una identidad algebraica del sector
de fondo. Ambos son adimensionales, así que la resta está **bien formada
aritméticamente** — pero restarle un parámetro de ecuación de estado a una
densidad no produce una densidad. La partícula construida sobre ese residuo
**no tenía de qué estar hecha**.

**2. La tensión S₈ que la motivaba no existía.**

El «3.5σ» se medía (a) contra el estadístico **comprimido** S₈, cuya tubería de
reducción asume un fondo ΛCDM, y (b) con `A_s` fijado al valor de Planck —
o sea, importando la tensión Planck–KiDS a un modelo que por sí solo no la
tiene. `A_s` es uno de los dos parámetros libres del modelo a nivel CMB.

Medido contra los **225 puntos crudos** de ξ± de KiDS-1000, con **un solo
sector** y `A_s` libre (MCMC convergido, R−1 = 0.019, N_eff = 4.2×10⁴):

    σ₈ = 0.7446 ± 0.0189      S₈ = 0.7555 ± 0.0192      →  0.11σ

**No hay tensión S₈**, y por lo tanto no había nada que el segundo sector
tuviera que hacer.

## Qué se perdió y qué se ganó

**Perdido:** una predicción pre-registrada y genuinamente falsable — el escalón
de free-streaming en P(k) que DESI Y3 / Euclid habrían confirmado o excluido.
No la refutó el dato: se retira porque la construcción que la generaba no era
física.

**Ganado:** un postulado menos, y cuatro problemas abiertos (OP-9, OP-10,
OP-11, OP-12) cerrados **por disolución** — ninguno se resolvió; todos dejaron
de ser preguntas cuando se retiró el objeto del que trataban.

## Dos lecciones durables

1. **La falsabilidad es el requisito mínimo para que una hipótesis sea
   científica, no evidencia de que la entidad exista.** Publicar una cota sobre
   esta partícula le habría dado existencia por la puerta de atrás.
2. **Una cadena puede ser dimensionalmente impecable y aun así no significar
   nada.** Verificar unidades no sustituye verificar que cada término sea *la
   clase de cosa* que dice ser.

## Dónde está lo vigente

- Código de crecimiento en uso: `src/p06_growth/`
- Paper reescrito: `manuscript/SSEE_Paper6_Growth.tex`
- Paper viejo archivado: `archive/manuscript_superseded/SSEE_Paper6_phiDM_TWOSECTOR_RETIRED.tex`
- Resultado canónico: `results/logs/growth_2026-07/R3_ssee_kids_S8.json`
- Valores: `CANONICAL_VALUES.yaml` (`S8_kids_mcmc`, `sigma8_kids_mcmc`)
