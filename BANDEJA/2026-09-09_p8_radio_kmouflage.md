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

## De dónde salió — la arqueología (lo pidió Mike)

**No dejó de ser coherente al cambiar otra cosa. Nació roto, en un solo
commit**: `295ed6e`, del 2026-05-15, titulado *«fix(paper8): OP-4 resolved —
replace Galileon Vainshtein with k-mouflage»*.

Lo que había antes era **peor y evidente**: la fórmula de Galileon daba
`r_V ≈ 1.8×10⁴⁴ m`, mayor que el radio de Hubble. Un número absurdo a simple
vista. Eso era OP-4.

El arreglo lo sustituyó por la fórmula k-mouflage con **exponente 1/3**, y el
propio mensaje del commit la llama *«Correct radius»* y marca OP-4 como
RESUELTO. `CLAUDE.md` todavía dice «OP-4 ✅ RESUELTO 2026-05-15».

**Y ahí está la lección, que es incómoda:** el arreglo de un problema de
unidades introdujo otro problema de unidades. Y lo empeoró en un sentido
concreto — cambió un número **absurdo** (10⁴⁴ m, que cualquiera caza) por uno
**plausible** (10⁷ m, que nadie mira dos veces). Un número absurdo se detecta
solo; uno plausible sobrevive meses.

## El exponente correcto SÍ se puede derivar, y es 1/2

No hace falta buscar en la literatura. El radio k-mouflage es donde el término
cinético no lineal alcanza al lineal, o sea donde `X ~ M⁴`, es decir
`φ' ~ M²`. En régimen lineal, alrededor de una fuente,
`φ' = M_obj/(4π M_pl r²)`. Igualando:

```
M_obj / (4π M_pl r²) = M²   ⟹   r² = M_obj/(4π M_pl M²)
```

Es decir **exponente 1/2**, que es justo el que cierra unidades. El paper puso
un cubo donde va un cuadrado.

## Qué NO toca — VERIFICADO 2026-09-09

**El titular de Paper 8 NO depende de este radio.** Comprobado leyendo los tres
sitios que lo usan:

- L621: *«Baryonic matter is independently protected by selective coupling
  **regardless of $r_{\rm km}$**»* — el cumplimiento con GR no cuelga de él.
- L697: el radio *«provides **supplementary** non-linear screening… but the EFT
  suppression … is the **dominant** mechanism for all solar-system tests»*.

O sea que la predicción de lente se sostiene sobre la pata EFT
(`α_B = α_M = 0`), y esa pata está intacta. **Lo que dice `CLAUDE.md` sobre
Paper 8 queda confirmado, ya no es una afirmación sin comprobar.**

**Lo que sí se cae** son tres cosas de §4.2, todas locales:
1. los tres números de la tabla y la figura;
2. la frase *«r_km ≪ 1 kpc para todos los objetos»* — con 1/2, la Vía Láctea da
   4 kpc;
3. la frase de que el radio da apantallamiento *«dentro de los cuerpos
   estelares»* — con 1/2 son 823 AU, muy fuera.

Nótese la dirección: el error corregido da **más** apantallamiento, no menos.
Eso hace **más fácil** pasar las pruebas del sistema solar, no más difícil. Lo
que se debilita es el argumento secundario de que la quinta fuerza opera a
escala galáctica.

## Reclasificación tras verificar

Sigue siendo **ROJO, pero acotado a §4.2**: hay tres números publicados que
están mal y dos veredictos de tabla que se invierten. No es rojo a nivel de
modelo: el titular del paper no se mueve.

## Decisión que es de Mike

1. **Derivar r_km bien desde K(X)** (cierra OP-4 de verdad) y rehacer tabla y
   figura, o
2. **retirar la tabla y la figura** del paper y dejar el radio como pendiente
   declarado, que es más rápido y más honesto si la derivación va a tardar.

Recomiendo la 2 como medida inmediata y la 1 como trabajo de fondo: el paper
no puede quedarse con tres números que no salen de su propia ecuación mientras
se deriva la buena.
