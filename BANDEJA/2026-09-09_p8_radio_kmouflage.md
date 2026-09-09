# ROJO — la tabla del radio k-mouflage de Paper 8 tiene TRES defectos

**Dónde:** `manuscript/SSEE_Paper8_StrongGravity.tex`, ecs. (rkm_sun / rkm_mw /
rkm_cluster), Tabla `tab:vainshtein` y Figura `fig:vainshtein`.
**Publicado en:** `docs/SSEE_Paper8_StrongGravity.pdf` (preprint, y en Zenodo).
**Encontrado:** 2026-09-09, al reclasificar la deuda por dependencia (Mike).
**Nada corregido todavía. Esto es un informe, no un arreglo.**

## Lo que dice el paper hoy

```
r_km = [ M_obj / (4π M_pl M²) ]^(1/3)      con M = 9.68 meV
```

| Objeto | M_obj | r_km publicado | «Status» publicado |
|---|---|---|---|
| Sol | 1 M_☉ | 1.45×10⁷ m ≈ 0.021 R_☉ | Inside stellar body |
| Vía Láctea | 10¹² M_☉ | 1.45×10¹¹ m ≈ 1 AU | ≪ 1 kpc |
| Cúmulo | 10¹⁵ M_☉ | 1.45×10¹² m ≈ 10 AU | ≪ 1 Mpc |

Y el pie de tabla concluye: *«All values satisfy r_km ≪ 1 kpc: the DM fifth
force is active at galactic and cosmological scales.»*

## Defecto 1 — las unidades no cierran

En unidades naturales el corchete vale
`GeV / (GeV · GeV²) = GeV⁻²`. Elevado a **1/3** da **GeV^−0.667**.
Una longitud es GeV⁻¹. **El exponente 1/3 no puede dar un radio.**
El único exponente que sí da longitud es **1/2**.

## Defecto 2 — los números no salen de su propia fórmula

Evaluando la fórmula publicada tal cual, con el M declarado:

| convención de M_pl | lo que da la fórmula (Sol) | lo que publica el paper |
|---|---|---|
| reducida, 2.435×10¹⁸ GeV | 1.44×10⁴ m | 1.45×10⁷ m |
| completa, 1.221×10¹⁹ GeV | 8.42×10³ m | 1.45×10⁷ m |

Falla por un factor ~10³ con cualquiera de las dos convenciones. El
**escalado** sí es el de 1/3 (masa ×10¹² ⟹ radio ×10⁴, comprobado), así que la
tabla es internamente coherente consigo misma pero **no con la ecuación que
la encabeza**. Los números vienen de otro sitio que el paper no dice.

## Defecto 3 — al reparar el exponente, la conclusión se rompe

Con **1/2**, que es el único que da longitud, y la misma M_pl reducida:

| Objeto | r_km reparado | ¿sobrevive su «Status»? |
|---|---|---|
| Sol | 1.23×10¹⁴ m = **823 AU** | **NO** — R_☉ = 7×10⁸ m; queda 176 000 radios solares fuera, no «inside stellar body» |
| Vía Láctea | 3.99×10³ pc ≈ **4 kpc** | **NO** — no es ≪ 1 kpc, es cuatro veces mayor |
| Cúmulo | 1.26×10⁵ pc = 126 kpc | **sí** — sigue ≪ 1 Mpc |

**Dos de las tres filas cambian de veredicto.** Y con ellas la frase del pie de
tabla, que es la que sostiene el argumento de que la quinta fuerza opera a
escala galáctica.

## Por qué esto es ROJO y no naranja

Por el criterio de Mike: rojo es *aquello de lo que dependen resultados ya
publicados*. Aquí no es una brecha por derivar — es una **fórmula publicada con
las unidades rotas**, unos **números que no salen de ella**, y una **conclusión
que cambia** al repararla. Un árbitro lo ve en la primera lectura, y la
verificación son tres líneas de aritmética.

## Lo que NO afirmo

Que el exponente correcto sea 1/2. Es la reparación mínima que cierra unidades
conservando el corchete, pero **el radio k-mouflage depende de la forma de
K(X)** y derivarlo bien es el contenido de OP-4. Lo que sí está establecido sin
ninguna hipótesis es que **la fórmula publicada no puede ser un radio**.

## Qué NO toca

La predicción de lente de Paper 8 se apoya en `ω_c = KAL₀·ω_b·n_s` y en
`α_B = α_M = 0` de Paper 7, no en este radio (CLAUDE.md, entrada de Paper 8).
Hay que **comprobarlo explícitamente** antes de asegurar que el titular del
paper no se mueve. No lo he comprobado todavía.

## Decisión que es de Mike

1. **Derivar r_km bien desde K(X)** (cierra OP-4 de verdad) y rehacer tabla y
   figura, o
2. **retirar la tabla y la figura** del paper y dejar el radio como pendiente
   declarado, que es más rápido y más honesto si la derivación va a tardar.

Recomiendo la 2 como medida inmediata y la 1 como trabajo de fondo: el paper
no puede quedarse con tres números que no salen de su propia ecuación mientras
se deriva la buena.
