# KiDS **sí** mide A_s. BOSS **no**. Y KiDS lo mide gracias a que tu fondo es rígido

**Corrida:** `src/p06_growth/mide_As_o_producto.py` ·
**Log:** `results/logs/growth_2026-07/As_medido_o_producto.json` ·
**Control PASA** · sobre cadenas ya corridas, no se corrió nada nuevo.
**No toca ningún paper.**

## 1. La pregunta, que la puso Mike

A_s es primordial, no cambia con la época. El fondo cósmico lo mide directo en
la altura de los picos. Entonces, ¿por qué las sondas tardías ven menos
amplitud en su propia época? Dos respuestas, y había que separarlas:

- **(a)** el crecimiento tardío realmente es menor, y hay física que contar;
- **(b)** la sonda no puede separar A_s de su acompañante, así que su A_s baja
  «no porque realmente baje, sino porque no puede verlo».

## 2. Cómo se separa, que es medida y no opinión

De la misma cadena salen dos anchuras para cada parámetro. La **marginal** es
lo que la sonda sabe de A_s sin ayuda de nadie. La **condicional** es lo que
sabría si todo lo demás estuviera clavado. Su cociente **D** dice cuánto de su
anchura viene de no poder separarlo.

D cerca de 1 significa medido por sí mismo. D grande significa que lo que se
mide es una combinación, y que el valor central depende de por dónde se
deslizó el ajuste.

**Control, primero:** una gaussiana con una degeneración puesta a mano de 0.99
debe dar D = 7.09 en los enredados y 1.00 en el libre. Salió 7.08, 7.08 y 1.00.

## 3. El resultado

| cadena | D de logA | enredado con | veredicto |
|---|---|---|---|
| **KiDS con fondo SSEE** | **1.16** | `halo_A`, r = −0.39 | **lo mide por sí mismo** |
| BOSS con fondo SSEE | 4.01 | `b1_z1NGC`, r = −0.89 | ve una combinación |
| KiDS con fondo ΛCDM libre | 4.73 | `ω_c`, r = −0.95 | ve una combinación |
| BOSS con fondo ΛCDM libre | 3.94 | `b1_z1NGC`, r = −0.88 | ve una combinación |

## 4. Las dos lecturas, y la segunda no me la esperaba

**BOSS no mide A_s.** Lo que mide es `−0.95·logA − 0.24·b1`, con la amplitud y
el sesgo de galaxia pegados a −0.89. Su 2.7636 ± 0.0981 no es una medición de
la amplitud: es la proyección de una dirección que el dato no distingue. La
sospecha de Mike era correcta **para BOSS**.

**KiDS sí la mide, y la mide porque el fondo de SSEE está clavado por álgebra.**
La misma sonda, con fondo ΛCDM libre, pasa a D = 4.73 enredada con la densidad
de materia oscura a −0.95. O sea: la rigidez del modelo es lo que convierte a
KiDS en un medidor de amplitud. No es un accidente del ajuste, es una
consecuencia de tener el fondo fijo.

## 5. Qué cambia esto

**El promedio 2.8418 es peor de lo que ya habíamos dicho.** No solo estrecha la
barra artificialmente: mezcla **una medición real** (KiDS) con **un número
dominado por degeneración** (BOSS), y les da peso por su barra, que es
justamente lo que la degeneración distorsiona.

**La deriva que hay que probar es KiDS contra el fondo cósmico: 3.59σ.** No los
4.50σ del promedio.

**Y a la pregunta de por qué ven menos:** BOSS no ve menos, ve borroso. KiDS sí
ve menos, y eso sigue pidiendo explicación.

## 6. Deuda que deja

- Rehacer el barrido de ingredientes clavando **el A_s de KiDS**, no el
  promedio. El de BOSS pierde sentido como valor a clavar, aunque sirve como
  comprobación de consistencia.
- El `halo_A` de KiDS (r = −0.39 con logA) es el único acompañante que queda.
  Conviene ver cuánto de la deriva sobrevive si se le pone el prior de
  retroalimentación bariónica más estrecho.
