# Cola de corridas

**Abierto 2026-09-07.** Nace de una regla de Mike: *«dejarlo anotado no
significa que se vaya a hacer».* Una corrida que no está en esta cola con su
coste no existe — se olvida. Y como la máquina no puede con todas a la vez,
hay que **decidir** cuál corre, no lanzarlas según van saliendo.

## ORDEN DE ATAQUE (reordenado 2026-09-09 — lo pidió Mike)

> El número de cada corrida es un **identificador**, no un puesto en la fila:
> es el orden en que se abrió, no el orden en que se hace. Se reordena cuando
> cambia lo que sabemos, y hoy cambió.

**Corriendo ahora:** barrido de combinaciones (viejo, 12/15) · cadena KiDS
ΛCDM-fijo · control `fuga3 fondo`.

**Siguiente, por orden real de valor:**

1. **#19 — ¿la supresión es plana o depende de la escala?** Es la que decide si
   el déficit del 8.7% es amplitud o materia que se fuga. Todo lo demás de la
   línea del A_s depende de su respuesta.
2. **#18 — apretar `halo_A`.** Única escapatoria observacional que queda en KiDS.
   Barata y puede cerrar la pregunta antes que la #19.
3. **#17 — barrido con el A_s de KiDS. SUBE, lanzada 2026-09-09.** Los dos
   resultados de hoy apuntan a ella: el barrido viejo clavaba el promedio
   2.8418, cuya mitad de BOSS acaba de moverse 1.85σ. El único A_s que se
   sostiene como medición es el de KiDS.
4. **#14 — par de campos acoplado.** Línea distinta (el wₐ), sigue viva.
5. Lo demás (#4, #5, #7, #9, #10, #11, #13, #16) sin cambio de prioridad.

---

## Reglas de la casa

1. **Nada se lanza sin `preflight.py` verde en la misma línea de comando.**
2. **Nada se lanza sin `OMP_NUM_THREADS=1`.** Medido el 2026-09-07: al
   relanzar `cobaya_kids` sin esa variable, cada uno de los 4 procesos MPI
   tomó ~2.3 núcleos en vez de 1 y la carga subió a **52 sobre 12 núcleos**.
   Las librerías de álgebra abren hilos por su cuenta y se pelean entre sí.
3. **`setsid nohup ... < /dev/null &`**, nunca `nohup` a secas. Se perdieron
   dos corridas el 2026-09-07 porque colgaban del shell y el shell murió.
4. **Presupuesto: 12 núcleos.** No pasar de ~10 ocupados. Una corrida MPI de
   4 procesos gasta 4; el resto queda para lo interactivo.
5. Al terminar una corrida: **escribir su resultado en el Registro antes de
   lanzar la siguiente.** Un resultado sin propagar es un resultado perdido.

## Corriendo ahora

| corrida | núcleos | lanzada | coste medido/estimado | qué desbloquea |
|---|---|---|---|---|
| `cobaya_kids lcdmfijo` | 4 (mal, ver regla 2) | 07-09 19:43 | **15 h y sigue**; `R−1` 0.226 (06:51) → **0.165** (10:57), para en 0.03 | la casilla que falta de la tabla 2×2: ¿ΛCDM con fondo fijo también muestra la tensión en `A_s`? |
| `punto_de_fuga2` | 1 | 08-09 17:32 | est. ~1 h | con la amplitud clavada en el valor tardío, qué ingrediente del fondo pide moverse y cuánto recupera cada uno |

## Terminadas 2026-09-08 (madrugada, sin Mike delante)

| corrida | coste medido | resultado |
|---|---|---|
| `analiza_lcdm_R4` | 10 h 02 min | reproduce lo publicado (`S₈=0.757105±0.019354`, `χ²=262.7462/212`); cierra R35 |
| `ajuste_conjunto_wc_ns` | ~1 h | SSEE y ΛCDM ajustan el CMB igual (χ² 1002.77 vs 1002.75); el `ω_c` que pide el CMB queda a 0.3% del que da la identidad |
| `perfil_wc_boss` | 18 min | **la brecha de amplitud NO es `ω_c` disfrazado**; traspasado a Paper 6 |
| `perfil_wc_boss_lcdm` | 16 min | el desplazamiento es del DATO; el `ω_c` algebraico cae más cerca del dato que el de Planck |

## En cola — ordenadas por (valor / coste)

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| ~~1~~ | ~~**ΔBIC de Paper 3 con `τ` ajustado**~~ ✅ **HECHA 08-09, 43 min** — **ΔBIC = −22.59** (publicado −24.02): NO mejora, empeora 1.43. SSEE sí gana al recuperar su `τ` (1005.41 → **1003.586**) pero ΛCDM también se mueve y la resta se lo come. Lo que queda es mejor titular: **SSEE ajusta igual que ΛCDM (1003.586 vs 1003.769 sobre 271 puntos) con cuatro perillas menos**. Control PASA: el mínimo de ΛCDM reencuentra Planck a ≤0.44σ en los cuatro. Informe: `BANDEJA/2026-09-08_dbic_tau_ajustado.md` | 43 min | — |
| ~~2~~ | ~~**perfil de `w_c`, control ΛCDM**~~ ✅ **HECHA** — `results/logs/cmb_perfil_wc.json`. SSEE `ω_c = 0.119334 ± 0.000246`, ΛCDM `0.119748 ± 0.000252`; **razón de anchuras 0.977** ⟹ la barra la pone el DATO, no el álgebra. De paso mostró que el `0.119534` publicado en Paper 8 estaba sesgado por la rejilla (paso 0.0020 sobre un tramo donde χ² sube 900: mandan las alas sobre el vértice) | ~40 min | la mitad que faltó al morir con el reinicio. Sin él no se puede afirmar que el `±0.000248` de SSEE lo pone el dato y no el álgebra | — |
| ~~3~~ | ~~**pendiente `d ln w_c / d ln n_s`**~~ ✅ **HECHA** — `results/logs/cmb_ns_forzado.json`, pendiente medida **−0.042** (la identidad predice +1). Y su gemela `cmb_wb_forzado.json`: pendiente **+0.430**, con el mínimo de χ² justo en el `ω_b` algebraico. Forzar un ingrediente algebraico saca al modelo de sí mismo, así que esto mide la verosimilitud, no la fórmula | ~1 h | complementa la #2 de la cola de arriba: impone `n_s` y mide si `w_c` responde con pendiente +1 | — |
| 4 | **BOSS ΛCDM con fondo LIBRE** | **cientos de horas** | la única versión publicable de R1/R2: hoy ΛCDM corrió con el fondo fijo, lo que le impide mostrar sus propias tensiones. Requiere reconstruir las tablas LPT por muestra (~18 s/llamada) | rediseño previo: emulador o templates precalculados |
| 5 | **`b1_*` de Paper 3** | horas–día | las 4 figuras rancias (41 días) que R36 marca | — |
| ~~6~~ | **`fuga2` — punto de fuga con `τ` libre** 🔄 **LANZADA 08-09 17:32** | ~1 h | con `A_s` clavado donde lo ponen KiDS+BOSS, qué ingrediente quiere moverse. La v1 tenía `τ` congelado y dio una base absurda (Δχ²=13 618) | — |

## Descartadas, con su razón

- **Barrido Kaiser de `fσ₈`** — fue sondeo, no resultado: el `Δχ²` cambia de
  signo según dónde se corte en `k` (+0.8 a −11.3 entre `k=0.06` y `0.12`).
  Kaiser no publica. Superada por LPT.
- **N-body para `S₈`** — 5 000–20 000 CPU-horas. Fuera de esta máquina.

## Añadido 2026-09-07 (tarde)

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 7 | **R7 — prueba de transporte de ingredientes** (diseño de M. Almeida, Paper 6 §sec:r7, tabla `tab:r7`) | ajuste conjunto ~1 día + las dos dedicadas | ajusta las 3 sondas JUNTAS con la amplitud libre, y evalúa ese mismo juego de ingredientes contra cada sonda por separado sin re-ajustar. Para SSEE `Δχ²≡0` por construcción (no hay ingrediente que transportar); para ΛCDM `Δχ²>0` se **mide**. Es la formalización falsable de todo lo visto hoy | la celda BOSS necesita la #4 (fondo libre); la celda KiDS sale de `lcdmfijo` menos R4 |

**Nota sobre la celda KiDS de R7:** cuando termine `cobaya_kids lcdmfijo`, su
`χ²_min` menos el de R4 (`262.746`, fondo libre) da directamente el
`Δχ²(KiDS)` de ΛCDM. Es la única de las tres celdas que sale gratis.

## Regla 6 (2026-09-07, tras perder `lcdmfijo` dos veces)

**Toda corrida MPI se lanza desde un script en disco, no desde una línea
inline.** Medido: `setsid nohup mpirun ... &` escrito en la línea de comando
murió en silencio dos veces (a los ~35 min y a los ~2 min), sin dejar error en
el log — el log simplemente corta. Lanzada desde un `.sh` con `exec mpirun`,
sobrevive. El script fija además las tres variables de hilos (regla 2), que es
donde se olvidan.

## Resuelto 2026-09-07

- **`ajuste_conjunto_wc_ns`** ✅ — `n_s` y `ω_c` libres a la vez. El CMB elige
  `n_s=0.96703`, `ω_c=0.119332`; la identidad con ese mismo `n_s` pide
  `0.119698` → **−0.31%**. Sin `n_s` el error sería `+3.73%`: el índice
  espectral lleva la predicción de 3.7% a 0.3%, factor 12, con `n_s` elegido
  por el dato y no impuesto. Control ΛCDM: `+1.35%`, cuatro veces peor.
  χ² gemelos (1002.772 vs 1002.750) con los mismos 4 libres.
  Log: `results/logs/cmb_ajuste_conjunto_wc_ns.json`.

## Reprioridad 2026-09-07 — la #4 sube a lo más alto

**Razón (observación de Mike):** todas las pruebas de la identidad
`ω_c = KAL₀·ω_b·n_s` hechas hoy viven en **una sola época**. `ω_c` sólo se
mide con precisión en el CMB (`±1%`); en KiDS sale con `±21%` y además
deslizándose por la degeneración con `A_s`, así que no es una medición. Con un
solo punto no se puede distinguir *«la identidad es una ley»* de *«la
identidad vale en `z=1100`»*.

La salida existe: `ω_m` deja huella tardía en la **forma** del espectro, a
través de la escala de igualdad materia-radiación `k_eq ∝ ω_m`, que un
análisis de forma completa de BOSS con el fondo libre sí mide, a `z≈0.5`.

⟹ **la corrida #4 deja de ser «la versión publicable de R1/R2» y pasa a ser la
única prueba capaz de decidir si la identidad es ley o coincidencia de época.**
Sigue bloqueada por el coste (reconstruir las tablas LPT por muestra, ~18 s por
llamada); el trabajo previo es un emulador o una rejilla de templates. Ese
rediseño pasa a ser la tarea de mayor valor de la cola.

**Anotado aparte, sin perseguir:** en ese mismo ajuste de ΛCDM a KiDS con el
fondo libre, `H₀ = 73.15 ± 4.96`, que cae a **0.02σ de SH0ES** (73.04±1.04) y
a 1.16σ de Planck. Con esa barra cabe cualquiera de los dos, así que es una
coincidencia sugerente, no una medición — pero conviene volver a mirarla si
alguna vez se estrecha.

## Añadido 2026-09-08 (tanda autónoma)

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| ~~8~~ | ~~**techo σ₈ con el fondo canónico y sus neutrinos**~~ ✅ **HECHA** — sale **0.814854 / S₈=0.826827**, un 2.3% por debajo del `0.8335 / 0.846` publicado; control ΛCDM pasa a 0.04σ del criterio previo. Informe: `BANDEJA/2026-09-08_techo_sigma8_neutrinos.md` | 1 min | si el techo publicado es el del modelo canónico o el de una variante sin neutrinos masivos | — |
| 9 | **auditoría de las configuraciones rescatadas** | ~1 h | los `.ini` de `config/class/` rescatados el 08-09 **no llevan neutrinos masivos**: medido, `ncdm` aparece en 0 de 5 (`ssee_v36`, `_canonical`, `_IS`, `_nomira`, `lcdm_planck2018_ref`). El único con `ncdm` es `_twosector`, y ahí es la **partícula retirada**, no un neutrino. Sus salidas alimentan números publicados: los picos `ℓ = 220, 535, 811` del Unified Journal y el **31.5% sin MIRA** (Unified L110/504/534). Falta ver si el `ℓ=220, 536, 813` de Paper 3 sale de CLASS o de CAMB | la #8, que ya dio el método y el control |

## Regla 7 (2026-09-08) — la configuración se versiona o el número no existe

`class_ssee/` es un fork de CLASS y trae **el `.gitignore` de CLASS**, que
ignora `*.ini` y `output/`. Consecuencia medida: **ninguna configuración de
CLASS de este proyecto estaba en el repositorio**, así que ningún número salido
de CLASS se podía certificar. No era un caso aislado, era el estado normal.
Desde hoy viven en `config/class/`, versionadas, y se corren desde ahí.

## Añadido 2026-09-08 (tarde) — lo vio Mike preguntando por σ₈

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 10 ✅ | ~~re-correr R1/R2 de BOSS con la masa de neutrino correcta~~ **CERRADA 2026-09-09** (`BANDEJA/2026-09-09_boss_R1R2_neutrinos.md`): Δχ²=+0.402, empatan; y el logA de BOSS se movió **1.85σ**, confirmando que no mide A_s** | ~1 h (LPT, 222 pts) | `boss_lpt_R1R2.py` tenía `MNU = 0.06` suelto **usado para los dos modelos**. 0.06 eV es el fiducial de Planck; la de SSEE es 0.06849 eV. Mismo patrón del `τ` prestado: un modelo evaluado con el ingrediente del otro. **Medido antes de arreglar** (mismo fondo, CAMB, z=0.51): σ₈(0) −0.276%, fσ₈(0.51) −0.247% = **0.058σ** de la barra. Pequeño pero con signo, no ruido. El código ya está arreglado y lee del núcleo; falta re-correr | — |

**Comprobado al mismo tiempo, y sale limpio:** `cobaya_kids.py`, que produce el
S₈ canónico de Paper 6 (R3 y su control R4), **sí** pasa `mnu = S.SUM_MNU_EV`
para SSEE y `0.06` para ΛCDM. El titular `S₈ = 0.7555 ± 0.0192` **no** está
afectado. El evaluador del CMB `cmb_eval.py` también lo lee del núcleo. El
agujero de los neutrinos era **sólo** de las corridas de CLASS y de esta línea
de BOSS.

## Añadido 2026-09-08 (noche) — deuda de la propagación del techo

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 11 | **recomputar la barra ±0.006 del techo de σ₈** | ~30 min | al propagar el techo (0.8335→0.814854) se recomputó el **valor** pero **no su incertidumbre**: el `±0.006` que imprimen Papers 2, 3, 5, los dos Journals, el PRD y el README viene de antes y nadie sabe de dónde sale. Con esa barra se calculan las tres tensiones publicadas (2.7σ / 2.8σ / 0.4σ), así que **no es cosmética**. Hay que decidir de qué es esa barra (¿precisión del integrador? ¿banda de A_s? ¿de Ω_m?) y medirla, o retirarla | la #8, que ya dejó el `.ini` y el método |

**Por qué está aquí y no como nota (regla de Mike, 2026-09-08):** *«dejar una
nota por ahí sólo hace que te des cuenta cuando vuelves a pasar por ahí; si no
lo hacemos en meses, queda anotado sin que nadie lo vea».* Toda deuda que deje
una propagación entra en esta cola el mismo día, o no existe.

## Añadido 2026-09-08 (noche) — deuda de la prueba de dos campos

Informe: `BANDEJA/2026-09-08_dos_campos_phi_pi.md`.

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 12 ✅ | ~~dos campos genéricos de w constante~~ CERRADA 2026-09-08 | minutos | la prueba de hoy mató la pareja *(campo de Paper 7) + (w constante)*: el campo de Paper 7 no ayuda, estorba, y pasado el 8% pediría densidad negativa. Queda sin probar la pareja **genérica**: un campo normal (w>−1) más uno fantasma (w<−1), los dos de w constante. Tres números libres contra 400 puntos, así que **puede fallar** | ninguna |
| 14 | **par de campos ACOPLADO (término de intercambio Q)** | ~1 h | la #12 cerró la vía sin acoplar: dos fluidos de w constante siempre mueven su w **hacia abajo**, y SSEE lo pide **hacia arriba** (−1.309 → −0.840). Es fallo de sentido, no de precisión, así que ningún par de valores lo arregla. Sobrevive una sola puerta: que los dos campos **se pasen energía**. Hay que probar un Q y ver si el sentido se invierte. Es la primera prueba que puede salir que SÍ | la #12, ya cerrada |
| 13 | **los TRES wₐ que conviven en la suite** | lectura | −0.669975 (álgebra, el que se compara con DESI), 0 (lo que se le entrega a hi_class en Paper 7) y +0.4135 (lo que de verdad hace el campo de Paper 7 al evolucionar). No es un error de cálculo: son tres cosas distintas mal etiquetadas con el mismo símbolo. Hay que decidir qué nombre lleva cada uno y dónde se dice | la #12 |

## Añadido 2026-09-08 (noche) — lo vio Mike preguntando por qué el fondo fijo tarda

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 15 ✅ | ~~separar rápido/lento en `cobaya_kids.py`~~ **HECHO 2026-09-08** (commit c201945) | 1 h de código | **7 de los 9 libres NO tocan CAMB** (A_IA, dz1..dz5, delta_c): sólo `logA` y `halo_A` entran en `_camb_cached`. Pero la verosimilitud es UNA función externa con los 9 parámetros en un solo bloque, así que Cobaya recomputa CAMB también cuando mueve los 7 baratos. **Medido: 8.965 s el paso que toca CAMB, 0.491 s el que no — 18x.** La caché ya existe y acierta el 100% en los pasos baratos; lo único que falta es declararle a Cobaya el `blocking` con velocidades para que los oversamplee | ninguna |

**Cómo se cerró, y por qué se relanzó después de todo.** Mi primera
recomendación fue dejarla correr, y el argumento era **coste hundido**: «ya
esperamos 25 h». Mike lo tumbó con el criterio correcto — *¿esperar da datos que
no se verían de otro modo?* No: las dos cadenas persiguen la MISMA distribución,
así que esperar no compra información, sólo horas. Diagnóstico por parámetro que
lo decidió: los caros ya iban convergidos (logA 0.0044, halo_A 0.0087) y los que
frenaban eran los baratos (dz1/dz2 0.0235). Se relanzó con `blocking` 1/18 y
`oversample_power` 0.7, **sembrada con la covmat aprendida** por la corrida
vieja, que se archivó entera en `chains_p6/kids/sin_reparto_2026-09-08/`.
Medido tras relanzar: **5531 pasos/h contra 379 = 14.6x** (ventana de 3 min).

**`cobaya_boss.py` NO lo necesita, comprobado:** sus tablas LPT se calculan una
vez por conjunto y viven en `st['lpt'][name]['L']`; no hay CAMB en el bucle y
todos sus libres (logA y los sesgos) recorren el mismo camino. No hay caro y
barato que separar.

**Lección para el guardián:** una verosimilitud externa de un solo bloque
esconde su propia estructura rápido/lento. Cada vez que se escriba una, hay que
preguntar qué parámetros entran de verdad en la parte cara.

## Añadido 2026-09-08 (noche) — lo señaló Mike leyendo el barrido a medias

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 16 | **recuperación DENTRO de rango (±3σ de Planck)** | ~4 h | el barrido de combinaciones deja a cada ingrediente libre entre cotas **numéricas** (`LIM`), no físicas, así que el optimizador se va a −21.3σ en ω_c y +20.5σ en H₀ para recuperar el 58.4%. **Una recuperación comprada a 20σ no dice que el dato prefiera nada**: dice que hay una degeneración por la que A_s se puede canjear. Hay que repetirlo acotando cada ingrediente a ±3σ de Planck y volver a medir. Sólo ese número contesta «cuánto del castigo se absorbe **sin salirse de lo que el propio dato permite**» | el barrido en curso, para comparar |

**Por qué está aquí (lo dijo Mike, 2026-09-08):** *«lo único que señala es que
A_s está absorbiendo la discrepancia de dos o más ingredientes; y si esos
tuvieran que ser corregidos fuera de su rango, no estaría señalando algo más que
la preferencia del dato»*. Yo leí la columna de porcentajes y pasé por alto la
de σ, que estaba en la misma tabla. **Regla que deja: un porcentaje de
recuperación no se lee nunca sin el desplazamiento que lo compró.**

## Añadido 2026-09-08 (noche) — tras medir si BOSS y KiDS miden A_s

Informe: `BANDEJA/2026-09-08_As_medido_o_producto.md`.

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 17 | **barrido de ingredientes clavando el A_s de KiDS (2.8627)** | ~8 h | medido: KiDS **sí** mide A_s (D=1.16) y BOSS **no** (D=4.01, enredado con el sesgo a −0.89). Así que el valor a clavar es el de KiDS, no el promedio ni el de BOSS. Script listo: `src/p03_cmb/fuga3_por_As.py kids`, con las dos cotas en la misma pasada | el control `fondo`, en marcha |
| 18 | **¿cuánto de la deriva de KiDS sobrevive con un prior estrecho en `halo_A`?** | ~4 h | `halo_A` (retroalimentación bariónica) es el único acompañante que le queda a logA en KiDS, r = −0.39. Es débil, pero es el que hay. Si con un prior estrecho la deriva de 3.59σ no se mueve, entonces es física | la #17 |

**Lo que esto retira:** clavar el A_s de BOSS deja de tener sentido como
medición (sigue valiendo como comprobación de consistencia), y el promedio
2.8418 queda desaconsejado por mezclar una medición con un número dominado por
degeneración. La deriva canónica a probar pasa de **4.50σ (promedio)** a
**3.59σ (KiDS sola)**.

## Añadido 2026-09-09 — la hipótesis de Mike: materia que está pero no se agrupa

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 19 | **¿la supresión es plana o depende de la escala?** | ~1 día | **la prueba que decide, y no hay que postular nada.** Bajar A_s suprime TODAS las escalas por igual. Que parte de la materia no se agrupe suprime **solo por debajo de su escala de fuga**, y deja las grandes intactas. KiDS mide ξ± en un rango de escalas, así que **el dato puede distinguirlas**. Se ajusta KiDS dos veces: (a) A_s libre, supresión plana — es lo que ya está; (b) A_s **clavado en el del fondo cósmico** (3.0448) más una supresión con escala de corte libre. Gana quien dé menor χ². Si gana (b), la hipótesis de Mike es correcta y además queda medida la escala | nada; los ingredientes ya están |

**Los números que la motivan (medidos 2026-09-09):** el modelo YA cuenta materia
que no se agrupa — los neutrinos, 0.52% del presupuesto, que ya suprimen σ₈ un
2.06%. Para cubrir el 8.7% que falta haría falta que **otro 2.18%** no se
agrupara, 4.2 veces más que los neutrinos actuales. Si fueran neutrinos serían
Σm_ν = 0.358 eV, muy por encima del límite habitual de 0.12 eV: **no cabe**.

**AVISO, y va en mayúsculas.** Esta hipótesis tiene **la misma forma** que la
partícula φ-DM retirada el 2026-08-01. Lo que ha cambiado es que la tensión que
la motivaba **ahora sí está medida** (3.46σ, con las dos escapatorias cerradas)
en vez de venir de un estadístico comprimido. Lo que NO ha cambiado es la
lección que dejó su retirada: la densidad tiene que salir de algo real, no de
una resta. La #19 no postula ninguna partícula: pregunta si la **forma** de la
supresión es la de una fuga de escala o la de una amplitud plana. Eso lo decide
el dato, y si sale que no, se cierra otra vez y ya está.

## Añadido 2026-09-09 — ROJO en Paper 8, hallado al reclasificar la deuda

Informe: `BANDEJA/2026-09-09_p8_radio_kmouflage.md`.

| # | corrida | coste | qué contesta | depende de |
|---|---|---|---|---|
| 20 🔴 | **derivar r_km desde K(X) y rehacer tabla+figura de P8** | días | la tabla `tab:vainshtein` publica tres radios con **tres defectos**: el exponente 1/3 sobre un corchete GeV⁻² da GeV^−0.667 (no es longitud); los números no salen de la propia fórmula (factor ~10³ con cualquier convención de M_pl); y al reparar el exponente a 1/2 **dos de las tres filas cambian de veredicto** (el Sol pasa de «dentro del cuerpo estelar» a 823 AU, la Vía Láctea de «≪1 kpc» a 4 kpc). Es el contenido real de OP-4 | ninguna |
| 21 ✅ | ~~¿la predicción de lente de P8 depende del radio k-mouflage?~~ **NO, verificado 2026-09-09** | lectura | Comprobado en las tres menciones: L621 dice que los bariones están protegidos por acoplamiento selectivo **«regardless of r_km»**, y L697 que el radio es **«supplementary»** frente a la supresión EFT, que es **«the dominant mechanism»**. El titular se sostiene sobre `α_B=α_M=0`. La #20 queda acotada a §4.2 | — |

**Medida inmediata recomendada (decisión de Mike):** retirar tabla y figura del
paper y dejar el radio como pendiente declarado, mientras se deriva el bueno. Un
paper no puede sostener tres números que no salen de su propia ecuación.
