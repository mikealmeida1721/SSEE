# Cola de corridas

**Abierto 2026-09-07.** Nace de una regla de Mike: *«dejarlo anotado no
significa que se vaya a hacer».* Una corrida que no está en esta cola con su
coste no existe — se olvida. Y como la máquina no puede con todas a la vez,
hay que **decidir** cuál corre, no lanzarlas según van saliendo.

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
| `dbic_tau_ajustado` | 1 | 08-09 11:47 | est. ~1–2 h (Nelder-Mead 6D, cada χ² es un CAMB) | si el `ΔBIC=−24.02` de Paper 3 cambia al ajustar `τ` en los dos modelos en vez de prestarle a SSEE el fiducial de ΛCDM |

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
| ~~1~~ | ~~**ΔBIC de Paper 3 con `τ` ajustado**~~ 🔄 **LANZADA 08-09 11:47** | ~1 h | si el `ΔBIC=−24.02` mejora al no prestarle el `τ` fiducial de ΛCDM. Medido aparte: SSEE gana 1.84, ΛCDM 0.47 ⟹ podría mejorar ~1.4, pero hay que rehacerlo con la precisión de lente y el `Σm_ν` del paper, no con los míos | — |
| ~~2~~ | ~~**perfil de `w_c`, control ΛCDM**~~ ✅ **HECHA** — `results/logs/cmb_perfil_wc.json`. SSEE `ω_c = 0.119334 ± 0.000246`, ΛCDM `0.119748 ± 0.000252`; **razón de anchuras 0.977** ⟹ la barra la pone el DATO, no el álgebra. De paso mostró que el `0.119534` publicado en Paper 8 estaba sesgado por la rejilla (paso 0.0020 sobre un tramo donde χ² sube 900: mandan las alas sobre el vértice) | ~40 min | la mitad que faltó al morir con el reinicio. Sin él no se puede afirmar que el `±0.000248` de SSEE lo pone el dato y no el álgebra | — |
| ~~3~~ | ~~**pendiente `d ln w_c / d ln n_s`**~~ ✅ **HECHA** — `results/logs/cmb_ns_forzado.json`, pendiente medida **−0.042** (la identidad predice +1). Y su gemela `cmb_wb_forzado.json`: pendiente **+0.430**, con el mínimo de χ² justo en el `ω_b` algebraico. Forzar un ingrediente algebraico saca al modelo de sí mismo, así que esto mide la verosimilitud, no la fórmula | ~1 h | complementa la #2 de la cola de arriba: impone `n_s` y mide si `w_c` responde con pendiente +1 | — |
| 4 | **BOSS ΛCDM con fondo LIBRE** | **cientos de horas** | la única versión publicable de R1/R2: hoy ΛCDM corrió con el fondo fijo, lo que le impide mostrar sus propias tensiones. Requiere reconstruir las tablas LPT por muestra (~18 s/llamada) | rediseño previo: emulador o templates precalculados |
| 5 | **`b1_*` de Paper 3** | horas–día | las 4 figuras rancias (41 días) que R36 marca | — |
| 6 | **`fuga2` — punto de fuga con `τ` libre** | ~1 h | con `A_s` clavado donde lo ponen KiDS+BOSS, qué ingrediente quiere moverse. La v1 tenía `τ` congelado y dio una base absurda (Δχ²=13 618) | — |

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
| 9 | **auditoría de las configuraciones rescatadas** | ~1 h | los `.ini` de `config/class/` rescatados el 08-09 (`ssee_v36*`) **no llevan neutrinos masivos**. Cualquier número citado que salga de ellos hereda el mismo sesgo del 2.3% que se acaba de encontrar. Hay que ver cuáles alimentan un número publicado | la #8, que ya dio el método y el control |

## Regla 7 (2026-09-08) — la configuración se versiona o el número no existe

`class_ssee/` es un fork de CLASS y trae **el `.gitignore` de CLASS**, que
ignora `*.ini` y `output/`. Consecuencia medida: **ninguna configuración de
CLASS de este proyecto estaba en el repositorio**, así que ningún número salido
de CLASS se podía certificar. No era un caso aislado, era el estado normal.
Desde hoy viven en `config/class/`, versionadas, y se corren desde ahí.
