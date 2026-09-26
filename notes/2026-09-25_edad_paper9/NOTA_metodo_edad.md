# Edad del universo — Paper 9: qué cambia con Ω_m=0.308881

Fecha: 2026-09-25. Corregida 2026-09-25 noche (reclamo de Mike — ver "Corrección"
abajo). Estado: **investigación — NO propagado al manuscrito.**
Pregunta de Mike: con s_m sale un universo más viejo, con Ω_m uno más joven — ¿por cuánto?

## Framework H del modelo (fijado por Mike 2026-09-25)

SSEE tiene **dos** H, conectadas por el f_screen completo: **H_local** (medida
SH0ES, 73.04) y **H_global** (inferida, 67.962). La edad del universo es una
cantidad **global** → se calcula con **H_global**. Punto. No hay menú de H.

## Método (reproducible, córrelo tú mismo)

```
t_0 = (1/H_0) ∫_0^∞ dz / [(1+z)·E(z)]
E²(z) = Ω_m(1+z)³ + Ω_DE·f_DE(z)      (plano: Ω_DE = 1 − Ω_m)
f_DE(z) = (1+z)^{3(1+w_0+w_a)} · exp(−3w_a·z/(1+z)),  (w_0,w_a) = (−0.840,−0.670)
1/H_0 [Gyr] = 977.8 / H_0[km/s/Mpc]
```

Código: `edad_paper9.py` (numpy + scipy.integrate.quad). Radiación omitida: ~0.01% en t_0.

## Resultados — modelo vs modelo, cada uno con sus parámetros propios

| Configuración | Ω_m en E(z) | H_0 | t_0 |
|---|---|---|---|
| Paper 9 original (s_m como densidad) | 0.160050 | 73.04 (local) | 15.52 (paper) / 15.28 (repro) |
| **SSEE corregida — predicción del modelo** | 0.308881 | **67.962 (H_global)** | **13.73** |
| ΛCDM (parámetros Planck) | 0.315 | 67.4 | 13.80 |

Diagnóstico intermedio (**no** es predicción del modelo): si se corrige solo la
densidad y se mantiene el H del paper → 12.78 Gyr. Sirve únicamente para aislar
el bug de s_m en E(z) (−2.74 Gyr del error de densidad).

Diagnóstico de forma (**no** es edad de ningún modelo): a parámetros fijos, la
dinámica (w_0,w_a)≠(−1,0) aporta −0.02 Gyr frente a ΛCDM. Irrelevante.

## Lectura honesta

1. **SSEE 13.73 vs ΛCDM 13.80: indistinguibles (0.07 Gyr).** El "+1.72 Gyr más
   viejo" del paper muere por completo.
2. El paper mezclaba **dos** inconsistencias: densidad equivocada (s_m como Ω_m)
   **y** H_local en una cantidad global.
3. Con la edad corregida no hay tensión con cúmulos globulares (≳12.5 Gyr) en
   ningún sentido: el selling point "more accommodating" se retira, no se invierte.
4. Brecha repro 15.28 vs 15.52 (1.5%): detalle de integración/redondeo de (w_0,w_a).
   No mueve la conclusión.

## Corrección 2026-09-25 noche (reclamo de Mike, fundado)

La primera versión de esta nota presentaba "SSEE 12.78 (H_0=73.04)" como
candidata a la edad del modelo y evaluaba ΛCDM a H_0=73.04 en la sensibilidad.
Ambos son cruces de categoría: el modelo tiene H_local y H_global (f_screen
completo) y la edad global usa H_global; ΛCDM se evalúa con sus parámetros
Planck, no con H ajenos. Corregido arriba: la predicción del modelo es **13.73**.

## Qué falta

- El "adelante" de Mike para editar `manuscript/SSEE_Paper9_HubbleTension.tex`
  l.924–946 (retirar/reescribir el párrafo de 15.52 Gyr y el argumento de
  cúmulos globulares).
- Radiación Ω_r (~0.01%), curvatura (plano, como el paper) y neutrinos:
  despreciables para t_0, ya revisado.
