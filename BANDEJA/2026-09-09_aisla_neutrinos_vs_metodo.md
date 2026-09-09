# NO fueron los neutrinos: el 98% del desplazamiento del A_s de BOSS es el MÉTODO

**Corrida:** `src/p06_growth/boss_aisla_neutrinos.py` ·
**Log:** `results/logs/growth_2026-07/boss_aisla_neutrinos.json` · 55 min ·
**Pedida por Mike:** *«cuando algo cambia así, incluso si son buenas noticias,
siempre se reverifica antes de pasarlo al paper»*. **No toca ningún paper.**

## 1. El control, que va primero (R24)

| | logA | χ² |
|---|---|---|
| A · canónico, m_ν = 0.06849 | 2.9447867806 | 197.438 |
| C · repetición idéntica de A | 2.9447867806 | 197.559 → *(mismo logA)* |

**|A − C| = 0.00e+00**, criterio < 1e−5. Cero **exacto**: el perfil es
determinista, no hay ruido de método que restar. Cualquier diferencia medida
después es señal entera.

## 2. La descomposición

Total a explicar: el logA de BOSS pasó de **2.7636** (MCMC viejo) a **2.9448**
(perfil nuevo) = **+0.1811**.

| causa | cuánto | fracción | en σ de la barra vieja (0.0981) |
|---|---|---|---|
| masa de neutrino (A − B, método fijo) | **+0.0030** | **2%** | **0.03σ** |
| método (lo que queda) | **+0.1782** | **98%** | **1.82σ** |

Con m_ν prestada (0.06): logA = 2.94181 ± 0.12166, χ² = 197.559.
Con la suya (0.06849): logA = 2.94479 ± 0.12385, χ² = 197.438.

**La masa de neutrino mueve el A_s de BOSS tres milésimas.** Es
indistinguible de nada.

## 3. Corrección de mi informe de ayer — y de la lectura que hicimos los dos

El informe `2026-09-09_boss_R1R2_neutrinos.md` §4 decía: *«no se puede atribuir
el desplazamiento a una sola causa… no digo cuál pesó más»*. Eso era correcto
entonces. Ahora está medido y **la causa es el método, casi por entero**.

La lectura de que «BOSS, al tener en cuenta el peso de los neutrinos, deja de
estar en conflicto con el CMB» **no se sostiene**. Los neutrinos no hicieron
nada. El conflicto no se resolvió: **se leyó con otra regla**.

## 4. Y qué es exactamente ese «método» — verificado, no supuesto

Los dos usan **el mismo dato, el mismo k_max = 0.20, el mismo modelado LPT y el
mismo fondo**. Comprobado: `cobaya_boss.py` L71–73 importa `boss_lpt_R1R2` y le
presta su propio KMAX. **Lo único que cambia es el estimador.**

| | qué reporta | logA | mejor χ² alcanzado |
|---|---|---|---|
| viejo | **media marginal** de una cadena de 104 272 muestras | 2.7636 | **207.638** |
| nuevo | **mínimo del perfil** | 2.9448 | **197.438** |

Y ahí está lo que no esperaba: **la cadena vieja nunca encontró el fondo.** Su
mejor punto está **10.2 unidades de χ² por encima** del mínimo real, con la
misma verosimilitud. No es que las dos lecturas señalen dos puntos legítimos de
un mismo valle: es que una de ellas **no llegó al valle**.

La explicación está en la estructura del problema: a logA fijo, la
verosimilitud de BOSS **factoriza** en seis problemas independientes de 6
molestias cada uno. El perfil explota esa factorización; la cadena no, y se
pasea por 36 dimensiones sin tocar fondo en 104 mil muestras.

## 5. Qué se salva, y esto es lo importante

**El empate en fσ₈ SÍ se salva, y sale reforzado por dos vías.**

| | Δχ² (SSEE − ΛCDM) |
|---|---|
| método viejo (MCMC) | +0.062 |
| método nuevo (perfil) | +0.402 |

Ambos son cero estadístico sobre 222 puntos. **La conclusión no depende del
método**, porque es una *diferencia* entre dos modelos medidos con la misma
regla, y los errores de regla se cancelan.

Y ahora además sabemos que **no depende del préstamo de la masa**: cambiarla
mueve el χ² 0.12 unidades. El empate no lo compraba la m_ν prestada.

## 6. Lo que NO se puede pasar a Paper 6

**El «0.81σ contra el CMB» no es un resultado, es una comparación mal
emparejada.** Enfrenta el **mínimo del perfil** de BOSS contra la **media
marginal** de Planck. Son dos estadísticos distintos; para un parámetro
degenerado —y el de BOSS lo es, D = 4.01— se separan justamente por lo que aquí
se midió: 1.8σ.

El «2.87σ» viejo tampoco vale, porque su cadena no había llegado al mínimo.

**Ninguno de los dos números de tensión sirve.** Para comparar hay que medir
ambos lados con la misma regla: o perfil contra perfil, o cadena contra cadena.

## 7. Veredicto

| afirmación | estado |
|---|---|
| fσ₈ de BOSS: SSEE y ΛCDM empatan con LPT | ✅ **medido y robusto al método** — pasa a Paper 6 |
| el empate no depende de la masa prestada | ✅ **medido** (Δχ² = 0.12) |
| el A_s de BOSS se movió por los neutrinos | ❌ **falso** — 2%, 0.03σ |
| el A_s de BOSS ya concuerda con el del fondo | ⛔ **indeterminado** — comparación desemparejada |
| el A_s de BOSS mide algo | ❌ ya se sabía (D = 4.01); esto lo confirma otra vez |

## 8. Deuda que deja (va a `COLA_CORRIDAS.md`)

**Nueva #22 — perfil contra perfil.** Medir el A_s del fondo cósmico por
perfil, no por marginal, para que la tensión BOSS–CMB tenga las dos mitades en
la misma unidad. Sin eso no hay número de tensión que citar.

**Toca a #17 y a la cadena de KiDS.** Ambas clavan un A_s de BOSS; la cifra de
2.7636 que usan viene de la cadena que no llegó al fondo. Hay que decidir si se
reemplaza por 2.9448 — pero eso es exactamente lo que la #22 tiene que zanjar
primero.
