# ΔBIC de Paper 3 con τ ajustado: −22.59, no mejora — y el porqué es interesante

**Corrida:** `src/p03_cmb/dbic_tau_ajustado.py`, terminada 2026-09-08 ~12:31.
**Coste medido:** 2612 s = **43 min**.
**Estado:** ✅ completa, **control PASA**.

## 1. Qué se preguntó

Paper 3 publica **ΔBIC = −24.02** sobre el dato de Planck, comparando SSEE
contra ΛCDM. Ese número se obtuvo con `τ` (profundidad óptica, el parámetro que
mide cuánta niebla dejó la reionización) **puesto al valor que prefiere ΛCDM**,
no al que prefiere SSEE.

Prestarle a un modelo el parámetro de otro no es neutral: puede favorecerlo o
perjudicarlo, y no se sabe cuál sin medirlo. Aquí se ajusta `τ` **en los dos**.

## 2. Qué está fijo y qué libre

| | SSEE | ΛCDM |
|---|---|---|
| ω_b (densidad de bariones) | álgebra, 0.0224178 | **libre** |
| ω_c (densidad de materia oscura fría) | álgebra, 0.1195144 | **libre** |
| H₀ (constante de Hubble) | álgebra, 67.9621 | **libre** |
| n_s (índice espectral) | álgebra, 0.9655581 | **libre** |
| w₀, wₐ (ecuación de estado) | álgebra, −0.8399 / −0.6700 | fijos en −1 / 0 |
| **A_s (amplitud primordial)** | **libre** | **libre** |
| **τ (profundidad óptica)** | **libre** | **libre** |
| **cuenta de libres, k** | **2** | **6** |

Dato: Planck `plik_lite` TTTEEE + lowT + lowE, **N = 271** bandas.
`BIC = χ²_min + k·ln(N)`, con `ln(271) = 5.6021`.

## 3. Los números

| | χ²_min | k | BIC |
|---|---|---|---|
| SSEE | **1003.586** | 2 | 1014.790 |
| ΛCDM | **1003.769** | 6 | 1037.381 |
| | | | **ΔBIC = −22.59** |

Publicado con `τ` prestado: **−24.02**. Negativo favorece a SSEE en los dos
casos.

## 4. Por qué NO mejoró, que es lo que se preguntaba

Se esperaba que mejorara porque SSEE gana al recuperar su propio `τ`:
**1005.41 → 1003.586, gana 1.82**. Pero ΛCDM también se movió, y la resta se
come casi toda la ganancia. El resultado neto es **1.43 peor** que el publicado.

Lo interesante no es el ΔBIC sino la fila de arriba: **SSEE ajusta el dato
igual de bien que ΛCDM (1003.586 contra 1003.769, diferencia 0.18 sobre 271
puntos) usando cuatro perillas menos.** Todo el ΔBIC sale de esa diferencia de
perillas, no de que un modelo ajuste mejor. Eso es más honesto de contar que un
número grande.

## 5. El control

Criterio escrito en el script **antes** de correr: el mínimo de ΛCDM tiene que
caer a menos de 2σ de la línea base de Planck 2018 en sus cuatro parámetros de
fondo. Si el optimizador aterrizara lejos, estaría midiendo su propia
convergencia y no el modelo.

| | encontrado | Planck 2018 | desvío |
|---|---|---|---|
| ω_b | 0.022349 | 0.02237 | 0.14σ |
| ω_c | 0.120314 | 0.1200 | 0.26σ |
| H₀ | 67.1208 | 67.36 | 0.44σ |
| n_s | 0.964678 | 0.9649 | 0.05σ |

**PASA con holgura.** El optimizador reencuentra Planck por su cuenta, así que
la comparación es limpia.

## 6. Qué NO establece

- No es un ΔBIC bayesiano completo: usa el mínimo de χ², no la evidencia
  integrada. Es la misma aproximación que usa el número publicado, así que la
  comparación entre los dos es válida; el valor absoluto arrastra esa
  aproximación en los dos casos.
- El `N = 271` es una cuenta declarada (215 bandas de la parte de alta escala
  más 28 y 28 de las dos de baja). Si esa cuenta cambiara, el ΔBIC cambia por
  su `ln(N)`. Está anotada en el log para que se pueda discutir.

## 7. Qué habría que propagar, si Mike aprueba

| dónde | qué dice hoy | qué diría |
|---|---|---|
| `manuscript/SSEE_Paper3_CMB.tex` | ΔBIC = −24.02 | −22.59, con `τ` ajustado en los dos |
| `CLAUDE.md`, tabla de Paper 3 | −24.02 | ídem |
| `VERIFICATION_LEDGER.md` | fila con log `p3_cmb_reframe_nu_fix.log` | añadir esta corrida al lado |

**Recomendación: cambiarlo, y además cambiar cómo se cuenta.** El número nuevo
es más defendible (nadie le presta nada a nadie) y va en dirección
desfavorable, lo cual quita cualquier sospecha de haberlo buscado. Y el
titular debería dejar de ser el ΔBIC y pasar a ser: *el mismo ajuste con cuatro
parámetros menos*. Un referee acepta eso mucho mejor que un número de
selección de modelos que depende de cómo cuentes N.

**Riesgo de no tocarlo:** el número publicado se obtuvo prestándole a SSEE el
`τ` de ΛCDM, y eso es justo la clase de detalle que un referee pregunta.
