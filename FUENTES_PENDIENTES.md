# Afirmaciones con fuente pendiente de localizar

**Abierto 2026-07-27.** Nace de un error de método: al no encontrar la fuente de
una cifra del Paper 1, la **borré**. Mike lo corrigió:

> «No que elimines si no encuentras su origen, sino que lo marques. Puede que
> esté en otra parte donde no buscaste, o puede que no se sostenga — ¿pero cómo
> saberlo si lo eliminaste? Si mañana se borra la fuente de H, no vamos a
> eliminar H: volvemos a revisar de dónde salió y creamos la fuente.»

**Regla.** Una afirmación sin fuente localizada **no se borra**: se anota aquí
con qué dice, dónde vive, qué se buscó y qué hipótesis quedan abiertas. Sale de
esta lista de una de dos formas — se encuentra/reconstruye la prueba y se cita,
o se demuestra que es falsa y **entonces** se corrige. Nunca por no encontrarla.

Nada entra al modelo «porque sí». Si algo se sostiene sólo en intuición, no se
publica; el día que se pueda probar, entra **con su prueba**.

---

## ✅ FP-1 · «DR1→DR2 shrank the uncertainties by ~40%» — CERRADA 2026-07-27

- **Dónde:** `manuscript/SSEE_Paper1_Framework.tex:125` (§Predictive Register)
- **Estado:** ✅ **CERRADA** — fuente primaria localizada y citada
- **Buscado en:** `results/logs/`, `CANONICAL_VALUES.yaml`, `data/raw/` (no hay
  datos DR1), `archive/`, bibliografía de la suite, y literatura pública.

**Reconstrucción numérica.** Con el DR2 que ya cita Paper 2 (σ_w0 = 0.055,
σ_wa = 0.208) y los valores DR1 de literatura (σ_w0 ≈ 0.063, σ_wa ≈ 0.29):

| métrica | reducción |
|---|---|
| σ(w₀) sola | 13% |
| σ(wₐ) sola | 28% |
| **área del contorno 2D** | **37% ≈ ~40%** |

La cifra encaja con la **reducción del área del contorno (w₀,wₐ)**, que además
es la métrica correcta para la frase que la acompaña («el punto siguió dentro
del 68%»): lo que importa es el área, no cada eje por separado.

**Pendiente:** fijar los valores DR1 desde la fuente primaria (DESI 2024,
arXiv:2404.03002) y, o bien citar el área explícitamente, o bien reescribir la
frase como «el área del contorno 68% se redujo ~40%», que es lo que la cifra
mide. **No borrar hasta entonces.**

---

## ✅ FP-2 · «DR3 ~1.5× más preciso → ~0.5σ» — CERRADA 2026-07-27

- **Dónde:** `manuscript/SSEE_Paper1_Framework.tex:128`
- **Estado:** ✅ **CERRADA** — se retiró el factor inventado, no la conclusión

Las dos cifras no se deducen una de la otra. Con el punto inmóvil la tensión
escala con el inverso del error:

    1.5 × 0.24σ = 0.36σ        (no 0.5σ)
    2.0 × 0.24σ = 0.48σ ≈ 0.5σ  ✓

Y la literatura DESI describe el salto DR1→DR2 como una mejora de **factor ~2**
en las constricciones. Si DR2→DR3 repite ese factor, **0.5σ es la cifra
correcta y el «1.5×» es el intruso** — al revés de lo que supuse al principio,
cuando iba a «corregir» el 0.5σ.

**Pendiente:** confirmar el factor esperado para DR3 y dejar UNA sola cifra con
su aritmética a la vista. **No borrar hasta entonces.**

---

## Cómo se cierra una entrada

1. Se localiza o reconstruye la prueba → se cita en el `.tex` y se registra en
   `VERIFICATION_LEDGER.md` con su procedencia.
2. Se demuestra que la afirmación es falsa → se corrige el texto **citando la
   demostración**, no el silencio.

En ambos casos la entrada se archiva con la fecha y el desenlace.


---

## Desenlace de FP-1 y FP-2 (2026-07-27)

**FP-1.** Fuente primaria localizada: DESI 2024 VI (arXiv:2404.03002),
DESI+CMB+PantheonPlus da w₀ = −0.827 ± 0.063, wₐ = −0.75 (+0.29/−0.25). Contra
DR2 (w₀ = −0.838 ± 0.055, wₐ = −0.617 ± 0.208):

| métrica | reducción |
|---|---|
| σ(w₀) | 13% |
| σ(wₐ) | 23% |
| **área del contorno 68%** | **33–37%** |

La cifra era real; le faltaba la cita y el nombre de la métrica. El texto ahora
da **ambos valores con su referencia** y dice «~35% en el ÁREA del contorno»,
que es lo que la cifra mide y lo que la frase necesita. Si la hubiera borrado,
habríamos perdido un dato correcto.

**FP-2.** El «~1.5× tighter» no tenía fuente y era el número equivocado; el
«~0.5σ» era el defendible (la literatura describe DR1→DR2 como factor ~2, y
2 × 0.24σ = 0.48σ). Reescrito sin inventar factor: se enuncia la ley
—«con el punto inmóvil la tensión escala con el inverso del error»— y se da el
rango 0.4–0.5σ para un salto comparable al de DR1→DR2. La fecha 2027 se
conserva; la proyección numérica sin fuente, no.

**Lección de método:** en los dos casos mi primer impulso —borrar el 40%,
«corregir» el 0.5σ— habría empeorado el paper. Investigar antes de tocar no es
lentitud: es la diferencia entre corregir y romper.

---

## ✅ FP-3 · ¿Cuál de los dos α_K es el falsador de Euclid? — CERRADA 2026-07-27

- **Dónde:** Paper 1 Tabla 1 (fila α_K^field), Paper 7 §504/§600, PRD §4.1.
- **Estado:** ✅ **CERRADA — no había contradicción. La alarma fue mía.**

**Desenlace.** Leí `σ(α_{K,0}) < 0.1` como una **cota superior** cuando es la
**precisión** del forecast de Euclid. Con esa lectura correcta:

| valor | frente a σ=0.1 | lectura |
|---|---|---|
| α_K^field = 0.072567 | por debajo de la barra de error | indistinguible de cero, **no testeable aún** ✓ |
| α_K^eff = 0.403302 | ~4σ sobre cero | sería **detección**, no exclusión |

Nada quedaba «ya excluido». Paper 1 y Paper 7 **coinciden**: ambos ponen
α_K^field frente a Euclid, y la frase de Paper 1 («el valor usado para el
cross-check de hi\_class/EFTCAMB») es literalmente cierta — 0.403302 es lo que
se le *pasa* a esos códigos, no lo que ellos devuelven.

Queda en pie lo bueno de la investigación: la **unificación de nomenclatura**
(la misma cantidad se llamaba de tres formas y la otra de dos → ahora
α_K^field / α_K^eff en toda la suite).

→ Lo que la investigación SÍ destapó, en el Paper 7 y sólo ahí, se anota como
**FP-4** y se trata **cuando la lectura llegue al Paper 7**, no antes.

---

## FP-4 · Dos sobreafirmaciones en Paper 7 — ATENDER AL LLEGAR AL PAPER 7

- **Dónde:** `manuscript/SSEE_Paper7_EFT.tex` líneas 553 y 584 (y sus ecos en la
  Tabla `tab:alpha`, notas al pie $^a$/$^b$).
- **Estado:** 🟡 **anotado, NO perseguir** — detectado de paso mientras se leía
  Paper 1 pág. 4. Verificado que **no está en Paper 1** (su redacción es exacta).
  Verificar si aparece en Papers 2–6 al pasar por ellos.
- **No toca ningún número del modelo.** Son afirmaciones *sobre la evidencia*.

**(a) La falsación invierte el test — línea 584.**
Dice: «una detección de $\alpha_K > 0.1$ por Euclid falsificaría el EFT de SSEE».
Pero 0.1 es la **precisión** σ, no un techo. Si el que Euclid reconstruye fuera
α_K^eff = 0.403, entonces α_K > 0.1 sería justo **lo que SSEE predice** — el
paper estaría ofreciendo como falsador lo que sería su confirmación.

**(b) «hi\_class confirma independientemente al 0.005%» — línea 553.**
El script `src/p03_cmb/ssee_paper3_hiclass_check.py` calcula:

```
línea  57:  aK_alg = 3 * Om_DE   * (1 + w0_)
línea 212:  aK_z   = 3 * Om_DE_z * (1 + weff_z)     # "desde CLASS"
```

En z=0 son **la misma cuenta**: `weff_z[0] = w0_` y `Om_DE_z[0] = Om_DE/E²(0)`
con E²(0)=1 por construcción. El 0.005% es el residuo numérico de que CLASS no
devuelve E²(0) exactamente 1. **No es confirmación independiente.** Además no es
hi\_class: usa CLASS con fluido `w0_fld/wa_fld`, sin el módulo `eft_smg` — y el
propio Paper 7 (línea 850) lista implementar hi\_class como **trabajo futuro**.

Igual con EFTCAMB: `'RPHkineticity0': alphaK_ssee` — recibe α_K como **input**,
no lo mide.

**Corrección propuesta (nada se borra):**
0. *(añadido 07-27, ver FP-5)* revisar junto con la tensión P7↔P9 sobre la
   definición de $\alpha_K$ — puede que (a) se resuelva sola al fijarla.
1. reescribir la falsación en términos de σ: a cuántas σ mediría Euclid cada
   valor, sin invertir el signo del test;
2. degradar «confirma independientemente» a lo que el script sí prueba —
   **consistencia interna de la evolución del fondo**, que es cierto;
3. conservar EFTCAMB por lo que demuestra (**ghost-free, gradient-stable**), que
   es un resultado real, y no por lo que no demuestra (el valor).

**El problema.** La suite es consistente en que existen dos kineticidades, cada
una válida en su marco (Paper 7 lo explica bien). Lo que NO queda resuelto es
cuál de las dos confronta el dato:

| | valor | qué es |
|---|---|---|
| `α_K^field` | 0.072567 | campo desnudo, antes del acoplamiento conformal |
| `α_K^eff` | 0.403302 | fondo efectivo con el w₀ **observado**; confirmado por hi_class (0.005%) y EFTCAMB |

- Paper 1 Tabla 1 y Paper 7 §600 ponen **α_K^field** frente al umbral de Euclid
  (σ(α_K,0) < 0.1) → cómodo, 0.073 < 0.1.
- Pero Paper 7 §508 dice, textual, que «the observationally relevant effective
  value is α_K^eff = 0.403302».

Si lo observacionalmente relevante es el efectivo y el umbral de Euclid es 0.1,
**0.403 > 0.1** y la predicción se leería como ya excluida. Si Euclid mide el
campo desnudo, no hay problema.

**No se resuelve con grep ni redondeos:** depende de qué cantidad reconstruye el
weak lensing de Euclid. Hay que fijarlo desde la definición de Bellini–Sawicki y
declararlo en UN solo sitio, del que los demás documentos citen.

**Riesgo si no se cierra:** es una fila marcada «Prediction» en el Predictive
Register. Un referee que tire de ese hilo encuentra dos respuestas en el mismo
paper.

**Hecho mientras tanto (nomenclatura, sin tocar la física):** la misma cantidad
se llamaba de tres formas (`α_K^field` ×7, `α_K^bare` ×2, `α_K` a secas ×7) y la
otra de dos (`α_K^eff` ×5, `α_K^fluid` ×4). Unificado a los nombres de la fuente
—**α_K^field** y **α_K^eff**— y retirada de Unified la frase «using w₀
overestimates α_K by ≈5.5», que contradecía el marco de «dos frames
autoconsistentes» de Paper 7 al presentar uno de los dos como error.

---

## FP-5 · La misma definición de α_K da dos valores en P7 y P9 — AL LLEGAR AL PAPER 7/9

- **Dónde:** `SSEE_Paper7_EFT.tex:525` vs `SSEE_Paper9_HubbleTension.tex:190,254`
- **Estado:** 🟡 **anotado, NO perseguir** — detectado desde Paper 1 pág. 6.
  **No bloquea Paper 1:** ahí la aritmética decide sin ambigüedad (ver abajo).

Los dos papers parten de la definición Bellini–Sawicki y llegan a números
distintos por un factor ≈5.6:

```
P7 §525:  α_K = 2X /(H² M_Pl²)        con X = ρ_φ(1+w_φ)/2  →  3Ω_DE(1+w_φ) = 0.072567
P9 §254:  α_K = 2X /(KAL H² M_Pl²)    →  3·AURA(π−φ)/(2(φ+π)²) = 0.403302
```

La diferencia aparente es el **K_X = 1/KAL** que P9 incluye y P7 (con K(X)=X,
K_X=1) no. Puede ser perfectamente consistente —son acciones distintas— pero
**los dos papers escriben la misma ecuación con el mismo símbolo**. Hay que
fijar cuál acción corresponde a cuál valor y decirlo en UN sitio.

**Relación con FP-4:** si esto se resuelve, probablemente resuelva también la
afirmación (b) de FP-4 — porque aclararía qué es lo que hi\_class debería
reproducir y qué significa el 0.005%.

**Ya hecho, y es independiente de cómo se resuelva:** la fórmula
`f_screen = α_K/(3·MIRA) = (π−φ)/Ω²` se citaba con el símbolo **pelado** en 12
sitios de 6 documentos. La aritmética no admite duda —sólo α_K^eff = 0.403302
da 0.067253; α_K^field da 0.012101, factor 5.6 fuera— así que el símbolo se
marcó `α_K^{\rm eff}` en P1(4), P4(1), P7(1), P8(2), Sealed(2), PRD(2).
**Paper 9 y 10 NO se tocaron**: son los que *derivan* la fórmula, con el símbolo
definido localmente, y se revisan al llegar a ellos.
