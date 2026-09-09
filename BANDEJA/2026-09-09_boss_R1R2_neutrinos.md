# BOSS R1/R2 con la masa de neutrino corregida — y el A_s de BOSS se movió 1.85σ

**Corrida:** `src/p06_growth/boss_lpt_R1R2.py` (cola #10) ·
**Log:** `results/logs/growth_2026-07/R1R2_boss_lpt_kmax0.200.json` ·
7 min. **No toca ningún paper.**

## 1. El resultado que se buscaba: fσ₈ de BOSS con LPT

| | χ² | dof | logA |
|---|---|---|---|
| SSEE | 197.438 | 185 | 2.94479 ± 0.12385 |
| ΛCDM | 197.036 | 185 | 2.95468 |

**Δχ² (SSEE − ΛCDM) = +0.402, con los MISMOS libres.** Sobre 222 puntos eso es
cero estadístico: **empatan**. Es la fila que Paper 6 tenía pendiente para
fσ₈, y sale limpia — con modelado LPT de un lazo, no con Kaiser, que ya se
había descartado por depender del corte en k.

## 2. Lo que no se buscaba, y es más importante

El `logA` de BOSS pasó de **2.7636 ± 0.0981** a **2.9448 ± 0.1238**.

**Se movió +0.1811, que son 1.85σ de su propia barra vieja.**

Y con eso su tensión con el fondo cósmico se desploma:

| | logA | deriva contra el CMB (3.0448) |
|---|---|---|
| antes | 2.7636 ± 0.0981 | **2.87σ** |
| ahora | 2.9448 ± 0.1238 | **0.81σ** |

## 3. Por qué esto CONFIRMA el diagnóstico de anoche

Ayer se midió, sobre las cadenas, que **BOSS no mide A_s**: su factor de
degeneración es 4.01 y está enredado con el sesgo de galaxia a −0.89. Lo que
mide es la combinación `−0.95·logA − 0.24·b1`, no la amplitud.

**Un número dominado por degeneración desliza su centro cuando cambia el
método.** Es exactamente lo que acaba de pasar: 1.85σ de desplazamiento. No es
que antes estuviera «mal» y ahora «bien»; es que **ese número nunca estuvo
determinado por el dato**, y cada método lo deja en un punto distinto de la
misma raya.

## 4. Aviso de honestidad: hay DOS cambios a la vez

No se puede atribuir el desplazamiento a una sola causa:

- la **masa de neutrino** dejó de estar prestada (SSEE usa la suya, 0.06849;
  ΛCDM la de Planck, 0.06);
- el **método** cambió: la vieja era una cadena MCMC, ésta es un perfil por
  minimización a k_max = 0.200.

Separar las dos causas pediría una tercera corrida. **No lo he hecho**, y por
eso no digo cuál pesó más.

## 5. Qué NO cambia

**El veredicto de anoche sobre el A_s no se mueve.** Aquel se construyó sobre
**KiDS sola** (3.46σ), precisamente porque ya se había establecido que BOSS no
mide la amplitud. BOSS estaba fuera de esa cuenta antes de que este número se
moviera. La decisión de sacarlo queda vindicada por una vía que no se buscó.

## 6. Lo que SÍ queda tocado

El `logA` promedio **2.8418** que clava los barridos `fuga2` y `fuga3` es la
combinación de KiDS con **el BOSS viejo**. Esa mitad acaba de moverse 1.85σ.
Ese promedio ya estaba desaconsejado por mezclar una medición con un número
degenerado; ahora además está desactualizado. **La corrida que vale es la #17,
con el A_s de KiDS solo.**
