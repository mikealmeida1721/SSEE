# Dos campos tampoco bastan **si sólo conviven**. El álgebra pide que se INTERCAMBIEN

**Corrida:** `src/p07_eft/dos_campos_genericos.py` (cola #12) ·
**Log:** `results/logs/eft_dos_campos_genericos.json` · **Control PASA** ·
coste: segundos. **No toca ningún paper.**

## 1. Qué se probó

La corrida anterior mató una pareja concreta (el campo de Paper 7 más un
compañero). Esta prueba la versión limpia de tu corazonada, sin ese campo de por
medio: **un campo normal (w > −1) y uno fantasma (w < −1), los dos de w
constante**. Es la forma estándar de un modelo *quintom*.

Tres números libres — w_A, w_B y cuánto aporta hoy cada uno — contra 400 puntos
de la curva algebraica. Podía fallar.

**Control primero (R24):** se fabricó un total que **es** dos fluidos con
w_A = −0.70, w_B = −1.40 y fracción 0.35. El ajuste los devolvió con desvío
**0.0e+00**. La máquina mide.

## 2. El resultado, y lo raro que tiene

| | |
|---|---|
| w_A ajustado | −1.074441 |
| w_B ajustado | −1.074441 |
| desvío máximo \|Δw\| | **0.234** (criterio previo: < 0.01 igual, < 0.05 aproxima) |

**Los dos salieron iguales.** El ajuste no encontró una pareja: colapsó a **un
solo fluido** con el w promedio. Eso no es «ajusta mal». Eso es la máquina
diciendo que el segundo campo no le sirve para nada.

## 3. Por qué. Y esto es lo que vale de la corrida

No es cuestión de afinar valores. Es el **sentido del movimiento**.

Cuando dos fluidos de w fijo comparten el universo, el que se diluye más despacio
acaba mandando, y ése es siempre el **más negativo**. Así que la mezcla se corre
hacia él: **su w total siempre BAJA con el tiempo.** Comprobado con cuatro
parejas muy distintas, las cuatro bajan.

Y el álgebra de SSEE pide lo contrario:

| | en a = 0.30 | hoy |
|---|---|---|
| lo que pide SSEE | −1.3089 | **−0.8399** |

**SSEE SUBE.** Empieza fantasma y va saliendo de ahí. Es lo que ya intuiste
cuando preguntaste si «choca y regresa»: eso es exactamente lo que hace la curva.

Ningún par de valores arregla un problema de dirección.

## 4. La puerta que esto abre — y es una sola

Si la mezcla no puede subir mientras cada pieza guarda lo suyo, entonces las
piezas **no pueden guardar lo suyo**. Las opciones quedan reducidas a dos, y son
la misma vista de dos lados:

1. uno de los dos campos **no** tiene w constante, o
2. los dos **se pasan energía** el uno al otro.

Dicho en tu lenguaje: tus dos opuestos no bastan con **coexistir**. Tienen que
**tocarse**. El término que hoy no está en el lagrangiano de Paper 7 es el de
intercambio entre los dos sectores.

Eso ya no es una corazonada suelta: es la única forma que sobrevive después de
haber cerrado las otras con prueba y control.

## 5. Lo que NO dice

- No toca el acuerdo con DESI. El par (w₀, wₐ) sigue siendo algebraico y sigue
  ajustando el dato. Esto mide si el modelo concuerda **consigo mismo**.
- No dice que exista tal acoplamiento. Dice dónde hay que buscarlo, y que si no
  aparece, la energía oscura de SSEE no es «un campo escalar» en el sentido
  simple en que Paper 7 lo dice hoy.
- **No se buscó nada en el diccionario.** No hay número que buscar: no hubo par.

## 6. Deuda que deja

Va a la cola: probar un par **acoplado**, con un término de intercambio
Q entre los dos sectores, y ver si el sentido se invierte. Ésa es la primera
prueba que puede salir que **sí**.
