# El techo de σ₈ que publican los papers está 2.3% alto — le faltaban los neutrinos

**Corrida:** `class_ssee/techo_ssee_canonico.ini` + su control
`techo_lcdm_referencia.ini`, con CLASS v3.3.4. Terminadas 2026-09-08 11:38.
**Coste medido:** menos de 1 minuto cada una.
**Estado:** ✅ completa, **control PASA** con criterio escrito antes de correr.

> Esta corrida **supersede** el `σ₈ = 0.8335 / S₈ = 0.846` que citan Papers 2,
> 3 y 5, el Sealed Journal, el README y `CANONICAL_VALUES.yaml`. Nada se ha
> propagado: este informe es para que Mike decida.

## 1. Qué se preguntó, y por qué se volvió a preguntar

El techo contesta **cuánta estructura da el fondo de SSEE si la amplitud
primordial `A_s` se hereda de Planck en vez de ajustarse**. No es una
predicción: `A_s` es uno de los dos libres del sector CMB, así que fijarlo
importa una inferencia hecha dentro de ΛCDM.

Se volvió a correr porque el número publicado **no se podía certificar**. Salía
de `class_ssee/output/can_cold__pk.dat`, un fichero que:

- **no tenía archivo de configuración** — no constaba con qué se corrió;
- **estaba fuera del repositorio** — `class_ssee/output/` está en la lista de
  ignorados y nunca se commiteó;
- era **del mismo minuto** (26 jul, 20:49) que `can_part__pk.dat`, la corrida
  de dos sectores con la partícula, retirada el 2026-08-01. El techo nació
  como su término de comparación.

Y la sospecha concreta, que Mike planteó: **si esa corrida no llevaba
neutrinos masivos, no era el techo del modelo canónico.** El fondo canónico sí
los lleva, Σm_ν = 0.06849 eV.

## 2. Qué está fijo y qué libre

| | SSEE | control ΛCDM |
|---|---|---|
| ω_b | 0.0224177568 (álgebra) | 0.02237 (Planck) |
| ω_c | 0.1195144084 (álgebra) | 0.1200 (Planck) |
| h | 0.6796213732 (álgebra) | 0.6736 (Planck) |
| n_s | 0.9655581463 (álgebra) | 0.9649 (Planck) |
| w₀ / wₐ | −0.8399 / −0.6700 (álgebra) | −1 / 0 |
| **Σm_ν** | **0.06849 eV** | **0.06 eV** |
| A_s | **2.1e−9, clavado a Planck** | 2.100e−9 |
| libres | **ninguno** | ninguno |

Los ingredientes de SSEE se leyeron del núcleo, no se re-tecleraron. Misma
precisión, mismo z, mismo k máximo y el mismo integrador en los dos.

## 3. El control, con su criterio escrito antes

El criterio quedó escrito en el propio archivo de configuración **antes** de
correr, para que no se pudiera acomodar después: el mismo integrador sobre la
línea base de Planck 2018 tiene que dar σ₈ = 0.8111 ± 0.006.

| | σ₈ |
|---|---|
| criterio (Planck 2018) | 0.8111 ± 0.006 |
| medido | **0.810851** |
| desvío | 0.04σ |

**El método queda calibrado.** El intento anterior daba 0.8221, un 1.4% alto, y
la causa era exactamente ésta: las corridas ΛCDM guardadas tampoco llevaban
neutrinos masivos.

## 4. El número

| | σ₈ | S₈ |
|---|---|---|
| publicado en los papers | 0.8335 | 0.846 |
| **recomputado con el fondo canónico** | **0.814854** | **0.826827** |
| desvío del publicado | +2.3% | +2.3% |

Razón SSEE sobre ΛCDM con la misma receta: **1.0049**. SSEE da un 0.5% más de
estructura, no el 1.4% que salía antes de incluir los neutrinos.

## 5. Qué cambia en las tensiones

| | S₈ publicado 0.846 | S₈ recomputado 0.8268 |
|---|---|---|
| KiDS-1000 (0.759 ± 0.024) | 3.5σ | **2.74σ** |
| DES-Y3 (0.776 ± 0.017) | 3.9σ | **2.82σ** |
| Planck 2018 (0.832 ± 0.013) | 1.1σ | **0.36σ** |

**El «3.5σ» del techo era, en parte, neutrinos que faltaban.** Ya sabíamos que
era artefacto de fijar `A_s`; ahora se ve que además el número estaba mal
calculado dentro de esa condición.

Detalle que conviene mirar: el recomputado **0.8268 coincide con el
`S₈ = 0.8256` que Paper 5 obtiene por su vía independiente**, la de crecimiento
lineal IS, que ya reportaba 2.7σ. Las dos estimaciones convergen al incluir los
neutrinos. La separación entre «techo» y «estimación IS» era, en buena medida,
el neutrino ausente.

## 6. Qué NO establece esta corrida

- No dice nada sobre el resultado vivo. Con `A_s` libre contra el ξ± crudo de
  KiDS-1000, el sector único da **S₈ = 0.7555 ± 0.0192, 0.11σ** (Paper 6, R3).
  Ese sigue siendo el número que se publica.
- No he podido comparar contra la corrida vieja al detalle, porque la vieja no
  tiene configuración. Sé que las dos difieren un 2.3% y que la nueva pasa un
  control que la vieja no pasaba; no puedo demostrar cuál fue exactamente el
  ingrediente que faltaba en la vieja, sólo que con neutrinos sale esto.
- La barra `± 0.006` que arrastran los papers no la he recomputado. Viene de
  antes y no la he tocado.

## 7. Lo que habría que propagar, si Mike aprueba

Son **6 documentos y 2 archivos de datos**, todos con el mismo cambio de par de
números y de las tres tensiones asociadas.

| dónde | qué dice hoy |
|---|---|
| `manuscript/SSEE_Paper5_IS.tex` | 15 sitios: la ecuación con caja, la tabla de tensiones, la fila de comparación, el resumen |
| `manuscript/SSEE_Paper3_CMB.tex` | 3 sitios, tabla σ₈ |
| `manuscript/SSEE_Paper2_MCMC.tex` | 2 sitios |
| `manuscript/SSEE_Sealed_Journal.tex` | 3 sitios |
| `manuscript/SSEE_Unified_Journal.tex` | 3 sitios |
| `README.md` | 4 sitios |
| `CANONICAL_VALUES.yaml` | `sigma8_single_ceiling` |
| `VERIFICATION_LEDGER.md` | fila del techo, ya reescrita como procedencia incompleta |

**Recomendación.** Propagarlo. El número nuevo tiene configuración versionada,
control que pasa con criterio previo, y coincide con la vía independiente de
Paper 5. El viejo no tiene ni configuración. Y el cambio va en la dirección
prudente: baja una tensión que ya sabíamos artificial, así que no hay riesgo de
estar mejorando el modelo por conveniencia. El coste es una tarde de
propagación y recompilar cinco documentos.

**Riesgo de no tocarlo:** los papers publican un número que no reproduce ni su
propia herramienta, con tres tensiones citadas un 30% más altas de lo que salen.
Es exactamente el hallazgo que un referee hostil encontraría al primer intento
de reproducción.
