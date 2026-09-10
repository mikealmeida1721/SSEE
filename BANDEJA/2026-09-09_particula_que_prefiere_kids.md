# Cola #24 — qué partícula prefiere KiDS con el fondo de SSEE clavado

**Corrida:** `src/p06_growth/particula_que_prefiere_kids.py` · 0.77 h · 4 núcleos
**Salida:** `results/logs/growth_2026-07/particula_que_prefiere_kids.json`
**Predicción registrada ANTES:** `BANDEJA/2026-09-09_prediccion_particula_kids.md`
**No toca ningún paper.**

---

## 0. Los símbolos

| símbolo | qué es | en esta corrida |
|---|---|---|
| `m_x` | masa de la partícula, en eV | **libre**, 9 valores |
| `ω_x` | densidad que lleva, adimensional como `ω_c` | **libre**, 5 valores |
| `ξ` | su temperatura ÷ la de los neutrinos | **no libre**, sale de `m_x` y `ω_x` |
| `ΔN_eff` | energía relativista que aportó, en «neutrinos» = `ξ⁴` | **no libre**, sale de `ξ` |
| `halo_A` | apantallamiento por bariones. Bajo = suprime más | **libre**, rejilla de 5 sobre el prior entero U(2.0, 3.13) |
| `A_IA` | alineamiento intrínseco de las galaxias | **libre**, minimizado |
| `delta_c` | punto cero de la cizalla | **libre**, minimizado |

**Clavado en las 45 casillas, idéntico:** `ω_b` = 0.0224178 · `h` = 0.679621 ·
`n_s` = 0.965558 · `logA` = **3.044834** · `Σm_ν` = 0.06849 · `w₀`, `wₐ` de SSEE.
Y `ω_c` = 0.1195144 **menos** `ω_x`, para que `ω_m` no se mueva.

**Sonda:** KiDS-1000, 225 puntos ξ±. **Modelo:** SSEE, fondo rígido.

---

## 1. Los controles, antes de leer nada

| | qué comprueba | criterio | medido | |
|---|---|---|---|---|
| **C0** | toqué `run_camb`; sin partícula debe dar lo mismo | idéntico | 291.8870612400 = 291.8870612400 | **PASA** |
| **C1** | la materia total no se mueve | < 1e−6 | 2.78e−17 | **PASA** |
| **C2** | una partícula pesada y fría es materia oscura fría normal | \|Δχ²\| < 1.0 | **0.0003** | **PASA** |

C2 falló primero por 1.72 y **el fallo era mío**: al habilitar el canal pasé
`num_massive_neutrinos=3` y convertí el único neutrino masivo del pipeline en
tres ligeros. La forma lo delató, un 0.6% **plano en todas las escalas**, que
no es una supresión física. Sin ese control la rejilla entera sale sesgada.

---

## 2. LA REJILLA — χ², todas de SSEE, todas contra KiDS

`m_x` \ `ω_x` →

| `m_x` | 0.0010 | 0.0020 | 0.0030 | 0.0040 | 0.0050 |
|---|---|---|---|---|---|
| 0.3 eV | 269.66 | **265.82** | 270.82 | 286.48 | 308.91 |
| 0.5 eV | 271.41 | 266.45 | 266.33 | 272.89 | 285.96 |
| **0.8 eV** | 272.92 | 267.70 | **265.73** ← | 267.49 | 274.50 |
| 1.3 eV | 274.47 | 269.42 | 266.46 | **265.75** | 268.24 |
| 2.2 eV | 276.16 | 271.49 | 268.37 | 266.39 | **265.80** ⚠ |
| 4.0 eV | 278.04 | 274.11 | 271.31 | 268.96 | **267.22** ⚠ |
| 10.0 eV | 280.58 | 278.47 | 276.48 | 274.80 | **273.37** ⚠ |
| 50.0 eV | 282.21 | 282.24 | 282.25 | 282.24 | 282.24 |
| 3000 eV | 282.17 | 282.17 | 282.17 | 282.17 | 282.17 |

⚠ = el mínimo de esa fila cae en el **borde de mi rejilla**, no en un mínimo
real. Esas filas son cotas superiores.

Las dos últimas filas son el **límite frío**: la partícula se vuelve
indistinguible de materia oscura fría normal y el χ² vuelve exactamente a la
referencia. Eso es C2 visto en la rejilla entera.

---

## 3. LOS OCHO MEJORES, con sus parámetros en la misma fila

| `m_x` | `ω_x` | % de `ω_m` | `ξ` | `ΔN_eff` | `halo_A` | `A_IA` | **χ²** |
|---|---|---|---|---|---|---|---|
| **0.8 eV** | 0.0030 | 2.10% | 0.7066 | **0.249** | 2.900 | 0.5748 | **265.73** |
| 1.3 eV | 0.0040 | 2.80% | 0.6615 | 0.191 | 3.130 ⚠ | 0.5703 | 265.75 |
| 2.2 eV | 0.0050 | 3.50% | 0.5979 | 0.128 | 3.130 ⚠ | 0.5634 | 265.80 |
| 0.3 eV | 0.0020 | 1.40% | 0.8559 | **0.537** | 2.900 | 0.5777 | 265.82 |
| 0.5 eV | 0.0030 | 2.10% | 0.8264 | **0.466** | 3.130 ⚠ | 0.5431 | 266.33 |
| 2.2 eV | 0.0040 | 2.80% | 0.5551 | 0.095 | 2.900 | 0.6020 | 266.39 |
| 0.5 eV | 0.0020 | 1.40% | 0.7219 | 0.272 | 2.600 | 0.5948 | 266.45 |
| 1.3 eV | 0.0030 | 2.10% | 0.6010 | 0.131 | 2.600 | 0.5834 | 266.46 |

**Las varas:**

| | χ² | `halo_A` |
|---|---|---|
| sin partícula, `A_s` clavado | **282.17** | **2.000** ⚠ borde de abajo |
| soltar `A_s` del todo (cadena R3) | 265.44 | 2.568 ± 0.305 |
| mejor partícula | **265.73** | 2.900 |

---

## 4. Los cinco puntos de la predicción, con su veredicto

| # | lo que predije | qué salió | |
|---|---|---|---|
| 1 | el χ² baja hacia ~265 | 282.17 → **265.73**, gana 16.45 | **ACIERTO** |
| 2 | `ω_x` sale por debajo de 0.00296, entre 0.0015 y 0.0025 | **0.00300**, el 2.10% | **FALLO** |
| 3 | `m_x` sale por debajo de 2 eV | **0.80 eV** | acierto flojo (ver §5) |
| 4 | esa masa arrastra `ΔN_eff` > 0.3 y **Planck la prohíbe** | **`ΔN_eff` = 0.249, Planck la PERMITE** | **FALLO** |
| 5 | `halo_A` y `A_IA` vuelven a ~2.60 / ~0.55 | `A_IA` = 0.575 sí; `halo_A` = 2.900, a 1.1σ de la cadena | acierto parcial |

**El fallo del 4 es el resultado.** Aposté a que la partícula se moría entre
las dos sondas y **no se muere**. Yo predije en contra y el dato fue al revés.

**Y el fallo del 2 también importa:** dije que el 2.08% estaba mal calculado
porque mi regla `ΔP/P ≈ −8f` se quedaba corta. La rejilla, con física de CAMB
y sin mi regla, cae en **2.10%**. La regla llegaba al sitio correcto por un
camino que yo había declarado roto. Lo digo porque me equivoqué al desconfiar
de ella, no solo al confiar.

---

## 5. LO QUE ESTO NO MIDE — y es la mitad del resultado

**KiDS NO determina la masa.** Hay **8 puntos dentro de 1 de χ²** del mejor,
con masas de **0.3 a 2.2 eV** y `ΔN_eff` de **0.095 a 0.537**. El valle es
plano. El «0.80 eV» es la casilla más baja de una rejilla, **no una medición**,
y el paso de la rejilla en masa es un factor 1.6 entre casillas.

Lo que KiDS mide es una **combinación** de masa y densidad, no cada una.

**Y aquí es donde el CMB sí decide.** A lo largo del valle:

| tramo del valle | `ΔN_eff` | Planck (< 0.30) |
|---|---|---|
| 0.3 – 0.5 eV | 0.47 – 0.54 | **prohibido** |
| 0.8 eV | 0.249 | permitido |
| 1.3 eV | 0.191 | permitido |
| 2.2 eV | 0.128 | permitido |

**La pinza funciona pero no cierra.** `N_eff` corta el extremo ligero del valle
y deja vivo el resto. La lectura conjunta es **`m_x` ≳ 0.7 eV**, sin cota
superior desde estos datos.

---

## 6. Defectos declarados de esta corrida

1. **La rejilla se pega al borde en `ω_x` = 0.005** para `m_x` ≥ 2.2 eV. El
   valle continúa fuera de lo que barrí. **Hay que extenderlo.**
2. **Tres de los ocho mejores se pegan a `halo_A` = 3.130**, el borde de
   ARRIBA del prior. Piden menos apantallamiento del permitido. El mejor punto
   (2.900) no se pega, pero esos tres son cotas.
3. **La referencia sin partícula se pega a `halo_A` = 2.000**, el borde de
   ABAJO: sin la partícula KiDS pide más apantallamiento del que hay.
4. **La parte no lineal es una extrapolación.** HMcode-2015 está ajustado sobre
   simulaciones con materia oscura fría y neutrinos. Usarlo a otra temperatura
   afecta sobre todo a la profundidad. Declarado antes de leer.
5. **Solo se probó la familia térmica.** Una partícula producida de otra manera
   (por desintegración tardía, por ejemplo) podría ir rápida sin pagar en
   `N_eff`. Esta corrida no dice nada de esa rama.

---

## 7. Lo que este informe NO autoriza

- **No dice que la partícula exista.** Dice que KiDS la prefiere a no tenerla
  por 16.45 de χ², y que en el tramo `m_x` ≳ 0.7 eV Planck no la prohíbe.
- **No mide su masa.** El valle es plano en 8 casillas.
- **No reabre la partícula φ-DM** retirada el 2026-08-01. Aquella salía de una
  resta sin contenido físico; ésta de un déficit medido en el χ².
- **Ninguna cifra entra en ningún paper.**

## 8. Lo que hay que correr después

1. **Extender `ω_x`** hasta 0.010 y mapear el valle entero, que ahora mismo se
   sale por el borde.
2. **Cruzar con el CMB de verdad**, no con la cota de `N_eff` a mano: correr
   Planck con la misma partícula dentro y ver el χ² conjunto.
3. **`halo_A` al minimizador**, no en rejilla, para que los railes de arriba y
   de abajo dejen de ser artefactos míos.
