# ¿El segundo campo es π? La prueba dice NO — pero mata una pareja, no la idea

**Corrida:** `src/p07_eft/dos_campos_phi_pi.py` · **Log:** `results/logs/eft_dos_campos_phi_pi.json`
**Control PASA.** Coste: segundos. **Nada de esto toca ningún paper todavía.**

## 1. Qué se preguntó, exactamente

Tu corazonada: si un campo solo no puede cruzar la barrera fantasma, quizá hacen
falta dos, uno con la ley de φ y otro con la ley de π.

Para que la prueba pudiera **fallar**, no se ajustó nada al resultado. Se
**despejó**:

| pieza | quién la fija |
|---|---|
| el total | obligado a la forma CPL con w₀ y wₐ **algebraicos** |
| componente 1 | el campo de Paper 7, con su `u` fijado por el álgebra |
| componente 2 | **lo que sobra** de la resta. Su ecuación de estado sale despejada |
| f₁, cuánto aporta hoy el campo 1 | **el único número libre** |

Y `f₁` se fija con la condición más simple que puede cumplir un segundo campo:
que su w sea **constante**. Eso es una exigencia sobre todo el rango, no sobre un
punto, así que puede no tener solución.

## 2. El control, primero (R24)

Se fabricó un total sintético que **es**, por construcción, el campo 1 con
f₁ = 0.6 más un componente de w = −1.3. El despeje tiene que devolver esos dos.

```
despejado  f1 = 0.600000   w2 = -1.300000   dispersion = 4.5e-08   -> PASA
```

Nota de proceso: en la primera versión este control **falló**, y no por la
física. El mínimo es una ranura estrechísima entre una meseta y un acantilado, y
el buscador acotado se iba al borde. Se cambió por barrido denso más refinado
local. Es exactamente el caso que justifica R24: sin control primero, habría
publicado el resultado de una máquina rota.

## 3. El resultado

| | |
|---|---|
| f₁, fracción del campo de Paper 7 hoy | **0.0001** (el borde inferior) |
| w₂ medio del sobrante | −1.074658 |
| dispersión de w₂ (0 = constante) | **0.1360** |
| w₂ en a = 0.30 / 0.65 / 1.00 | −1.3106 / −1.0739 / −0.8399 |

**El sobrante no es un campo de w constante.** Se mueve 0.47 de punta a punta y
cruza −1 él solo, que era justo el problema de partida.

## 4. Lo importante no es el NO, es la FORMA del no

El `f₁` se fue al borde. Verifiqué si eso era un fallo del buscador o la
respuesta, y es la respuesta:

| f₁ | dispersión |
|---|---|
| 0.0001 | 0.136 |
| 0.01 | 0.165 |
| 0.05 | 0.448 |
| ≥ 0.08 | **la densidad sobrante se vuelve negativa** |

La dispersión **crece** desde el borde. Es decir: el campo de Paper 7 no sólo no
ayuda, **estorba**, y pasado el 8% exigiría que el otro componente tuviera
densidad negativa. El mejor valor de f₁ es «que no esté».

## 5. Qué queda en pie y qué no

**Muere** la pareja concreta *(campo de Paper 7) + (campo de w constante)*. El
lagrangiano k-esencia de Paper 7 no puede ser una de las dos mitades.

**No muere** tu idea de dos opuestos. Lo que muere es un candidato para el primer
opuesto. La prueba nunca examinó una pareja de dos campos genéricos, uno normal y
uno fantasma, ambos de w constante — ése es el siguiente paso natural y es
barato.

**No se tocó** el acuerdo con DESI. El par (w₀, wₐ) sigue siendo algebraico y
sigue ajustando a 0.24σ. Esta prueba es sobre si el modelo concuerda **consigo
mismo**, no con el dato.

## 6. Deuda que deja (va a la cola hoy)

- La pareja «dos campos genéricos de w constante» no está probada.
- Sigue sin resolverse que en la suite convivan **tres** valores de wₐ:
  −0.669975 (álgebra, el que se compara con DESI), 0 (lo que se le entrega a
  hi_class en Paper 7) y +0.4135 (lo que de verdad hace el campo de Paper 7).
- **No se buscó nada en el diccionario**, a propósito: primero el número, después
  la búsqueda. Aquí no hay número que buscar, porque no hubo solución.
