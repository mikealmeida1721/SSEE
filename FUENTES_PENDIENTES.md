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

## FP-1 · «DR1→DR2 shrank the uncertainties by ~40%»

- **Dónde:** `manuscript/SSEE_Paper1_Framework.tex:125` (§Predictive Register)
- **Estado:** 🟡 hipótesis fuerte, falta la cita exacta
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

## FP-2 · «DR3 ~1.5× más preciso → ~0.5σ»

- **Dónde:** `manuscript/SSEE_Paper1_Framework.tex:128`
- **Estado:** 🟡 inconsistencia interna — probablemente sobra el «1.5×», no el «0.5σ»

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
