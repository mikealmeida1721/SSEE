# δc y los halos tempranos — cheatsheet para futuro Mike

**La pregunta**: ¿SSEE predice más halos gigantes a z~10 que ΛCDM
(los "monstruos tempranos" de JWST/Euclid)?

## La intuición (la dinámica)

- Los halos nacen de **picos raros** de densidad. Todo es una carrera entre
  dos cosas: qué tan ALTA está la barra (δc ≈ 1.686, el umbral de colapso)
  vs. qué tan altos son los picos típicos (σ(M,z), que crece con el tiempo).
- Los halos masivos son picos rarísimos (eventos de 3–4σ). En la cola de la
  gaussiana todo es exponencial: e^(−ν²/2) con ν = δc/σ. Bajar la barra 3.4%
  a ν=4 multiplica los halos ×3; al mismo 3.4% a ν=1 casi no hace nada.
  **Por eso el efecto crece con la masa**: no es magia, es estadística de colas.
- A alto z, σ es pequeño (las fluctuaciones aún no han crecido) → hasta masas
  modestas son "picos raros" → el efecto también crece con z.

## Los dos casos

- **Derivado — Ruta 1 (vigente)**: la barra está a la misma altura en SSEE y
  ΛCDM (1.68647 vs 1.68646 a z=10). A esa época ambos universos están
  dominados por materia y la gravedad es la misma. Resultado: **≈1.00×**,
  sin diferencia.
- **Conjetura n_s (OP-27)**: la barra baja 3.4% (×n_s). Resultado: 1.06× en
  10¹¹ M☉ → 2.01× en 3×10¹² M☉ (z=10). Esa curva es su huella.

## Qué esperar ver en el cielo

- Si los conteos futuros salen como ΛCDM → la derivación queda, la conjetura
  sigue en OP. Nada cambia.
- Si sale un **exceso que crece con la masa** siguiendo la curva → evidencia
  para la conjetura; la Ruta 2 (derivar el mecanismo) se vuelve urgente.
- **HOY: no hay test posible.** Las incertidumbres son de factor ~varios y la
  señal de la conjetura a masas JWST es 3% → invisible. Sin prisa.

## Status (2026-09-25)

Conjetura motivada (n_s hace trabajo estructural vía ω_c = KAL_0·ω_b·n_s),
mecanismo pendiente. Borrador OP-27 + cambios Paper 4/5 en
`borrador_op27_deltac.md`, pendiente de tu revisión. Manuscrito no tocado.
