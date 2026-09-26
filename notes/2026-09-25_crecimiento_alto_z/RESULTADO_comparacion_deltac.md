# Comparación δc: Ruta 1 (derivado) vs postulado n_s — RESULTADO

Fecha: 2026-09-25. Script: `comparar_deltac_ruta1_vs_postulado.py` (este directorio).
Figura: `results/figures/fig_deltac_comparacion.pdf/.png`.

## Diseño

Press–Schechter con el FONDO SSEE en ambos casos (w₀=−0.83995, wₐ=−0.66997,
Ωm=0.308881, σ₈=0.8153 vs 0.811, D(z) con γ=0.5504/0.55). Lo ÚNICO que cambia
es δc:

- **Caso A (Ruta 1)**: δc_S(z), δc_L(z) del colapso esférico derivado.
  A z=10: 1.68647 vs 1.68646 (indistinguibles).
- **Caso B (postulado)**: δc_S=1.62839=1.68647×n_s, δc_L=1.68647.

## Tabla: n_SSEE/n_ΛCDM

| z  | M [M☉]  | Caso A (Ruta 1) | Caso B (n_s) | B/A = lo que compra n_s |
|----|---------|-----------------|--------------|------------------------|
| 10 | 3×10¹⁰  | 0.998           | 1.009        | 1.01                   |
| 10 | 10¹¹    | 0.990           | 1.052        | 1.06                   |
| 10 | 3×10¹¹  | 0.976           | 1.133        | 1.16                   |
| 10 | 10¹²    | 0.946           | 1.335        | 1.41                   |
| 10 | 3×10¹²  | 0.894           | 1.798        | 2.01                   |
| 12 | 10¹¹    | 0.984           | 1.085        | 1.10                   |
| 12 | 10¹²    | 0.922           | 1.513        | 1.64                   |
| 15 | 10¹¹    | 0.973           | 1.149        | 1.18                   |
| 15 | 10¹²    | 0.881           | 1.898        | 2.16                   |
| 15 | 3×10¹²  | 0.778           | 3.548        | 4.56                   |

## Lectura honesta

1. **La apuesta queda cuantificada.** El factor n_s vale: 1% a 3×10¹⁰ M☉,
   6% a 10¹¹, 41% a 10¹², ×2 a 3×10¹² (z=10). Esa curva masa-dependiente es
   la HUELLA que la Ruta 2 tendría que producir si algún día se deriva el
   mecanismo — y el listón observacional para distinguirla (abundancias a
   z~10 con precisión ~10–20% a M>10¹² M☉, era Euclid-futura, no actual).
2. **En el régimen JWST (10^10.8 M☉) la conjetura da ~1.03×**: ni siquiera el
   postulado resuelve el "too big too soon" (déficit 10–100×). Los factores
   grandes solo aparecen en objetos extremadamente raros (3×10¹²).
3. **Caso A ≈ 1.00×, sin diferencia significativa.** El leve déficit a masas
   extremas (0.89 a 3×10^12, z=10) viene de los anclas σ₈/M₈/D(z), no de δc,
   y está dentro de las incertidumbres de esos anclas. No se reclama como
   predicción.
4. **Lo que la diferencia NO dice.** La forma de la curva es la sensibilidad
   exponencial genérica de los picos raros a δc — no apunta a ningún
   mecanismo. Cuantifica la apuesta, no enseña a ganarla.

## Status

Conjetura motivada (n_s vía ω_c, OP-19) con apuesta observable cuantificada,
mecanismo pendiente → OP-27 (borrador en `borrador_op27_deltac.md`).
