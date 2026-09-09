# ΔBIC de Paper 3 con τ ajustado: **−26.21**, mejora — y el titular cambia

**Corrida:** `src/p03_cmb/dbic_tau_ajustado.py`, re-corrida 2026-09-08 tras
corregir el conteo de datos. **Coste:** 44 min. **Control PASA.**

> ⚠️ **Este informe ANULA la versión anterior**, que decía −22.59 y «no mejora».
> Aquélla estaba mal por dos cosas, las dos mías, y las dos las destapó Mike al
> pedir que se verificara antes de decidir sobre el título.

## 1. Los dos errores de la versión anterior

**(a) Dije que ΛCDM también mejoraba y se comía la ganancia. Falso.** ΛCDM va de
1003.760 a 1003.769: nueve milésimas, ruido del optimizador. **El único que
mejora es SSEE**, y mejora 1.82 al recuperar su propio `τ`.

**(b) Conté 271 puntos donde hay 669.** Escribí «215 bandas de plik_lite
TTTEEE», pero 215 son las de TT **solo**. Contadas en el propio `.clik`
(`cl_cmb_plik_v22.dat`): TTTEEE tiene **613**. Con lowT (28) y lowE (28),
**N = 669**. El BIC lleva `4·ln(N)`: 22.41 declarado contra **26.02** real.
El script ahora **mide** N del likelihood y aborta si no lo encuentra.

## 2. Qué está fijo y qué libre

| | SSEE | ΛCDM |
|---|---|---|
| ω_b · ω_c · H₀ · n_s | álgebra (leídos del núcleo) | **libres** |
| w₀ · wₐ | álgebra (−0.8399 / −0.6700) | fijos en −1 / 0 |
| A_s · τ | **libres** | **libres** |
| **k** | **2** | **6** |

## 3. Los números

| | publicado (τ prestado) | ahora (τ ajustado) |
|---|---|---|
| χ² SSEE, k=2 | 1005.409 | **1003.586** |
| χ² ΛCDM, k=6 | 1003.760 | **1003.769** |
| diferencia | +1.649 | **−0.183** |
| ΔBIC con N=669 | −24.37 | **-26.21** |

**Mejora 1.84.** Y la fila de la diferencia **cambia de signo**: con el `τ`
prestado SSEE perdía por 1.65; con el suyo **gana por 0.18**.

Aviso de honestidad: 0.18 sobre 669 puntos es **cero estadístico**. Lo correcto
es decir que **ajustan igual**, no que SSEE ajuste mejor. Lo que no es cero es
la cuenta de perillas.

## 4. Hallazgo colateral: el −24.02 publicado tampoco está bien contado

Su χ² de 1005.409 incluye **las tres piezas** del likelihood — medido en el
punto nuevo: plik_lite 584.170 + lowT 23.206 + lowE 396.210 = 1003.586 — pero
su N contaba **sólo la grande**. Con el conteo completo, el número publicado
debería haber dicho **−24.37**, no −24.02.

## 5. El control

Criterio escrito en el script **antes** de correr: el mínimo de ΛCDM debe caer a
menos de 2σ de Planck 2018 en sus cuatro parámetros de fondo.

| | encontrado | Planck | desvío |
|---|---|---|---|
| ω_b | 0.022349 | 0.02237 | 0.14σ |
| ω_c | 0.120314 | 0.1200 | 0.26σ |
| H₀ | 67.1208 | 67.36 | 0.44σ |
| n_s | 0.964678 | 0.9649 | 0.05σ |

**PASA con holgura.** El optimizador reencuentra Planck por su cuenta.

## 6. Qué propagar y qué decidir

**Propagar (no es decisión, es lo que toca):** −24.02 → **−26.21** en Paper 3,
`CLAUDE.md` y el Registro, con la nota de que el viejo estaba mal contado.

**Decidir (esto sí es tuyo):** el titular. Hoy Paper 3 vende el ΔBIC. La
alternativa es vender la fila de arriba: **el mismo ajuste con cuatro
parámetros menos**. Recomiendo el cambio, porque el ΔBIC depende de cómo
cuentes N —hoy mismo se ha visto— y la comparación de χ² a igual dato no
depende de nada.

**Riesgo de no tocarlo:** el número publicado se obtuvo prestándole a SSEE el
`τ` de ΛCDM y contando 613 puntos donde el χ² usaba 669. Las dos cosas las
encuentra un referee que intente reproducirlo.
