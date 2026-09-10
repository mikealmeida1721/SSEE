# PREDICCIÓN REGISTRADA — cola #24, antes de mirar la rejilla

**Escrito:** 2026-09-09, con los controles corriendo y la rejilla sin lanzar.
**Corrida:** `src/p06_growth/particula_que_prefiere_kids.py`
**Autor de la pregunta:** Mike. Yo la anoto y la mido.
**No toca ningún paper.**

## Los símbolos, antes de usarlos

| símbolo | qué es | ¿libre? |
|---|---|---|
| `m_x` | masa de la partícula, en eV | **sí, se barre** |
| `ω_x` | cuánta densidad lleva (adimensional, como `ω_c`) | **sí, se barre** |
| `ξ` | su temperatura dividida por la de los neutrinos | no: sale de `m_x` y `ω_x` |
| `ΔN_eff` | energía relativista que aporta, en «neutrinos». Vale `ξ⁴` | no: sale de `ξ` |
| `ω_c` | materia oscura fría de siempre. Se le **resta** `ω_x` | fijo por álgebra |
| `logA` | `ln(10¹⁰ A_s)` | **CLAVADO** en 3.0448 |

Materia total constante: `ω_c → ω_c − ω_x`, así que `ω_m` no se mueve.
Lo pidió Mike explícitamente: *«de ω_m sea el 100% y de este le quitas el
2.08%»*.

## Lo que ya sabemos, y es la vara

| | χ² |
|---|---|
| fondo SSEE clavado, `A_s` impuesto, sin partícula | **283.69** |
| lo mejor que consiguió la familia inventada de la #19 | 265.32 |
| soltar `A_s` del todo (cadena R3) | 265.44 |

O sea: hay **~18 de χ²** sobre la mesa, y la #19 los conseguía con una supresión
indistinguible de plana.

## Corrección a un número que ya te di

Le dije a Mike que la partícula lleva el **2.08%** de la materia. Ese número
salía de mi regla `ΔP/P ≈ −8·f`, aplicada al `A_sup = 0.1660 ± 0.0417` de la
#19.

**CAMB dice que esa regla se queda corta.** Con `ω_x = 0.002960` (el 2.08%) la
supresión real que sale no es del 16.6% sino del **23% al 27%**, según la masa.
Hay un efecto que mi regla ignoraba por completo: mientras la partícula va
rápida **añade radiación**, y eso retrasa el momento en que la materia manda,
lo que suprime el crecimiento **una segunda vez**.

Por eso `ω_x` va **libre** en esta corrida y no clavado en 0.002960. Si la
supresión que hace falta es del 16.6%, la densidad que la produce es más bien
**~1.4%**, no 2.08%. Que lo diga el dato, no mi regla.

## La predicción, en números

| # | predicción | qué la falsaría |
|---|---|---|
| 1 | el χ² **baja** desde 283.69 hacia ~265 con alguna combinación | que la rejilla entera se quede por encima de 280 |
| 2 | la mejor `ω_x` sale **por debajo** de 0.00296, entre 0.0015 y 0.0025 | que prefiera 0.004 o más |
| 3 | la mejor `m_x` sale **ligera**, por debajo de 2 eV | que prefiera 10 eV o más |
| 4 | esa `m_x` ligera arrastra un `ΔN_eff` **por encima de 0.3**, o sea que **Planck la prohíbe** | que salga con `ΔN_eff < 0.3`, y entonces la partícula sobrevive a las dos sondas |
| 5 | `halo_A` y `A_IA` vuelven a ~2.60 y ~0.55 en los puntos buenos, como en la #19 | que se queden deformados incluso en el mejor punto |

**La predicción 4 es la que importa y es la que espero que salga MAL para la
partícula.** Mis cuentas a mano de esta mañana daban una ventana de 0.68 a
0.95 eV, pero la verificación con CAMB mostró que mi fórmula de la malla
también estaba mal: a 1 eV el corte real cae en 44.8 Mpc, **dentro** de la
ventana de KiDS, no fuera. Con eso los dos bordes se cruzan y la ventana se
cierra.

**Si la rejilla lo confirma, el resultado es que la partícula NO existe**, y
eso es un resultado limpio que hay que decir igual de fuerte que el contrario.

## Los tres controles, y por qué C2 es el que manda

| | qué comprueba | criterio |
|---|---|---|
| **C0** | toqué `run_camb`, física validada. Sin partícula debe dar lo mismo bit a bit | idéntico |
| **C1** | la materia total no se mueve al meter la partícula | < 1e−6 |
| **C2** | una partícula **muy pesada y muy fría** es materia oscura fría normal, así que debe devolver el χ² de sin partícula | \|Δχ²\| < 1.0 |

**Estado: los tres PASAN.** C2 da 0.0003 contra un criterio de 1.0.

**C2 falló primero, por 1.72, y el fallo era mío.** Al habilitar el canal en
`run_camb` pasé también `num_massive_neutrinos=3`. Este pipeline llevaba **un
solo** autoestado masivo cargando todo Σm_ν (degeneración 1.0147); con el 3 lo
convertí en **tres** neutrinos ligeros repartiéndose la misma masa. Tres
ligeros escapan más que uno pesado, y eso movía el P(k) un **0.6% plano en
todas las escalas** — la forma delató que no era una supresión física sino un
desajuste de amplitud. Quitado el argumento, C2 pasa.

**Lo que esto dice del control:** C2 detectó un cambio del 0.6% en un
ingrediente que no tenía nada que ver con lo que la corrida mide. Sin él, la
rejilla entera habría salido con ese sesgo dentro y yo no lo habría visto.

## Lo que esta corrida NO hace

- **No reabre la partícula φ-DM** retirada el 2026-08-01. Aquella salía de
  restar una densidad menos una ecuación de estado. Ésta sale de un déficit
  medido en el χ².
- **No decide que la partícula exista.** Mide qué prefiere KiDS. Lo que Planck
  permite se cruza después.
- **Ninguna cifra entra en ningún paper.**
