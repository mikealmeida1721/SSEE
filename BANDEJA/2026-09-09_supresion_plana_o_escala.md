# Cola #19 — la supresión que le falta a KiDS es indistinguible de PLANA

**Corrida:** `src/p06_growth/supresion_plana_o_con_escala.py` · fondo SSEE ·
~1 h · **Salida en el `.out`, NO en `.json`** (ver §5) ·
**No toca ningún paper.**

**Qué pregunta.** Con el fondo clavado y solo la amplitud suelta, KiDS pide
~0.18 menos de `logA` que su propio fondo cósmico. Ya está descartado que lo
fabriquen los ingredientes (`fuga3 fondo`: se mueven <0.6σ), que lo fabrique el
fondo de SSEE (`lcdmfijo`: el de Planck da el mismo agujero) y que sea la regla
de medir (el volumen es `ln det F`, y KiDS no lo lleva). Queda **lo que se
observa**, y eso deja huella de escala.

## 1. Los dos controles, antes de leer nada

| | qué comprueba | criterio | medido | |
|---|---|---|---|---|
| **C1** | apagada la supresión, el χ² no depende de `k_c` | < 1e−6 | 283.6872 = 283.6872 | **PASA** |
| **C2** | plana ⟹ equivale a bajar A_s hasta el que midió la cadena | < 2σ | **1.47σ** | **PASA** |

C2 en detalle: la mejor supresión plana es `A_sup` = 0.225, que da un A_s
efectivo de **2.7899** contra el **2.8633 ± 0.0500** de la cadena.

> **Defecto declarado de C2** (dicho antes de ver el resultado): compara un A_s
> de **perfil** contra el **marginal** de la cadena — el mismo desemparejamiento
> de reglas de la cola #22. El criterio de 2σ es holgado y vale como
> comprobación de que la maquinaria suprime lo que dice, **no** como acuerdo de
> precisión.

## 2. La rejilla

Familia de dos parámetros que contiene las tres hipótesis:

```
P(k) → P(k) · [ 1 − A_sup · x²/(1+x²) ],   x = k/k_c
```

`k_c`→0 es plana (amplitud) · `k_c` medio es fuga por escala · `k_c` grande son
bariones. **La amplitud va IMPUESTA en la del CMB (3.0448)**: toda la reducción
la hace la supresión, no hay ningún A_s que ayude.

| `k_c` \ `A_sup` | 0.000 | 0.075 | 0.150 | **0.225** | 0.300 | 0.375 | 0.450 | 0.525 | 0.600 |
|---|---|---|---|---|---|---|---|---|---|
| 1e−6 *(plana)* | 283.7 | 272.9 | 267.4 | **266.6** | 269.6 | 276.8 | 289.8 | 308.3 | 331.3 |
| 0.020 | 283.7 | 272.7 | 266.8 | 265.5 | 267.7 | 274.0 | 286.2 | 303.9 | 326.7 |
| **0.042** | 283.7 | 273.0 | 267.1 | **265.3** ← | 266.8 | 272.4 | 283.6 | 300.3 | 322.6 |
| 0.089 | 283.7 | 273.7 | 268.2 | 265.9 | 266.3 | 270.7 | 280.5 | 295.7 | 316.4 |
| 0.189 | 283.7 | 275.2 | 270.5 | 268.1 | 267.4 | 270.6 | 278.6 | 291.6 | 309.5 |
| 0.400 | 283.7 | 277.5 | 274.3 | 272.3 | 271.2 | 273.2 | 279.1 | 289.0 | 302.9 |
| 0.846 | 283.7 | 280.0 | 278.1 | 277.4 | 276.7 | 277.0 | 279.9 | 285.4 | 293.6 |
| 1.789 | 283.7 | 281.7 | 280.6 | 280.3 | 280.7 | 280.7 | 280.9 | 281.8 | 283.9 |
| 3.783 | 283.7 | 282.7 | 281.8 | 281.3 | 280.9 | 280.7 | 280.8 | 281.1 | 281.6 |
| 8.000 | 283.7 | 283.1 | 282.5 | 282.0 | 281.5 | 281.1 | 280.7 | 280.3 | 280.0 |

## 3. Las tres lecturas, con las restas hechas

| | χ² | |
|---|---|---|
| amplitud del CMB, **sin** supresión | 283.7 | columna `A_sup`=0 |
| **mejor de la rejilla** (`k_c`=0.042, `A_sup`=0.225) | **265.3** | |
| mejor **plana** (`k_c`→0) | 266.6 | primera fila |
| A_s **libre**, sin supresión | 265.44 | cadena R3 |
| mejor con corte en escalas pequeñas (`k_c`≥1.8) | 280.0 | tres últimas filas |

**(a) La supresión hace el trabajo entero.** 283.7 → 265.3: gana **18.4**, y
queda a **0.14** de lo que consigue soltar A_s del todo. Una supresión de la
materia y una bajada de la amplitud son, para KiDS, la misma cosa.

**(b) La forma NO se distingue.** Plana 266.6 contra con-escala 265.3:
**Δχ² = 1.3 por un parámetro extra**. No es significativo por ningún criterio.
KiDS no tiene poder para decidirlo. El corte preferido, 0.042 h/Mpc, cae por
debajo de las escalas donde KiDS mide — o sea, indistinguible de plana.

**(c) Los bariones quedan EXCLUIDOS.** Un corte confinado a las escalas
pequeñas se queda en 280.0, **14.7 peor** que el óptimo. Es lo único contundente
de la rejilla, y descarta una de las tres ramas.

## 4. Cómo se lee esto (criterio corregido por Mike, registrado antes)

El criterio no es *«¿hace falta lo mismo en los dos fondos?»* — eso era mío y era
tosco. Es **comparar `k_c`, no `A_sup`**: la amplitud depende del presupuesto de
materia de cada modelo y **tiene** que diferir; la escala es de la física y
**tiene** que coincidir.

**Esa mitad está pendiente:** la corrida de ΛCDM abortó y hay que relanzarla.

## 5. Dos fallos míos en esta corrida

**El `.json` nunca se escribió.** `TypeError: Object of type bool is not JSON
serializable` — un booleano de numpy en la rama de aborto. **Mismo fallo que ya
arreglé hoy en otro script**: la rama de fallo no se prueba nunca y revienta el
día que por fin se usa. Los números están íntegros en el `.out`; son los de
arriba. Arreglado en `8f49b48`.

**El control C1 de ΛCDM falló por mi culpa, no por física.** Dio 290.0161 contra
290.0160 con un criterio de 1e−6. Mi perfil arranca tibio desde la solución
anterior y con tope de 40 iteraciones no converge a esa cifra: **el criterio
quedaba por debajo del ruido de mi propio minimizador**. Corregido para que C1
compruebe la **función** con las molestias fijas, que es exacto — y hecho **con
la rejilla de ΛCDM aún sin ver**.

## 6. Lo que este informe NO autoriza

- **No dice que haya materia que no se agrupa.** Dice que la supresión que KiDS
  necesita es indistinguible de una bajada de amplitud, y que **no** puede estar
  confinada a escalas pequeñas.
- **No cierra la #19**, porque falta ΛCDM, que es la mitad que distingue física
  de parche.
- **Ninguna cifra entra en ningún paper.**
