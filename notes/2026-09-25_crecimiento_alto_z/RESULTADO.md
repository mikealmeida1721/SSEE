# Crecimiento a alto z: Press-Schechter SSEE vs ΛCDM — 2026-09-25

Pedido de Mike en chat (~06:40 ET): cotejar los objetos masivos tempranos
(estrella de agujero negro MoM-BH*-1, cuásares de Euclid) contra el modelo.
Pregunta: ¿el fondo algebraico de SSEE acomoda halos masivos a z~7–10 mejor que ΛCDM?

## Método

Script existente `src/p02_mcmc/ssee_press_schechter.py` (Task 2B), corrido en la VM de Nova.
- Factor de crecimiento lineal D(z) por integración ODE con f = Ω(a)^γ:
  γ_SSEE = 0.5504 (Paper 5, medido), γ_ΛCDM = 0.55.
- Fondo SSEE: (w₀,wₐ) = (−0.8399,−0.6699), Ω_m = 0.308881, H₀ = 67.962.
- Umbral de colapso: δc_SSEE = δc_EdS × n_s = 1.6284 vs δc_ΛCDM = 1.6865 (−3.44%).
  ⚠️ El escalado δc × n_s es el ansatz del script — verificar su justificación
  en el manuscrito antes de citarlo como predicción del modelo.

## Corrección aplicada al script (2026-09-25)

σ₈ ancla: 0.7446 → **0.8153**. El 0.7446 era el MCMC R3 contra KiDS-1000
(A_s libre), superado el 2026-09-19 por KiDS-Legacy. Vigente:
`sigma8_ssee_unif = 0.8153` (CANONICAL_VALUES.yaml) — predicción del modelo
unificado con A_s fijo al CMB. Para esta pregunta el ancla correcta es la
predicción propia, no un dato recalibrado.

## Resultados

D(z)/D(0) — SSEE ≈ ΛCDM en todo el rango alto (diferencia ≤0.2%, SSEE
marginalmente MENOR): z=10 → 0.11513 vs 0.11535; z=15 → 0.07915 vs 0.07932.
Esperable: a z~10 ambos universos están dominados por materia; la DE no pinta.

Enhancement n_SSEE/n_ΛCDM a z=10 (viene casi todo del δc menor, no del crecimiento):
| M [M☉] | ratio |
|---|---|
| 1e11 | 1.05× |
| 3e11 | 1.13× |
| 1e12 | 1.34× |
| 3e12 | 1.80× |
| 10^10.8 (régimen JWST, Boylan-Kolchin 2023) | 1.03× |

## Veredicto honesto

Empujón en la dirección correcta, NO una resolución. El déficit ΛCDM de
Boylan-Kolchin 2023 era ~10–100× en densidad numérica; SSEE da 1.03× en ese
régimen. Además esa tensión se ha suavizado desde 2023 (interlopers,
contaminación AGN — irónicamente los propios LRDs). El modelo nunca queda
peor que ΛCDM aquí, pero este diagnóstico no convierte el "too big too soon"
en una victoria. Chequeo de consistencia con barras anchas, no kill-or-confirm.

Figura: `results/figures/fig_press_schechter.pdf` (+ .png).
