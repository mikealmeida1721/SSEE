# Edad del universo — Paper 9: qué cambia con Ω_m=0.308881

Fecha: 2026-09-25. Estado: **investigación — NO propagado al manuscrito.**
Pregunta de Mike: con s_m sale un universo más viejo, con Ω_m uno más joven — ¿por cuánto? ¿y cómo se calcula, para ver si falta algo?

## Método (reproducible, córrelo tú mismo)

```
t_0 = (1/H_0) ∫_0^∞ dz / [(1+z)·E(z)]
E²(z) = Ω_m(1+z)³ + Ω_DE·f_DE(z)      (plano: Ω_DE = 1 − Ω_m)
f_DE(z) = (1+z)^{3(1+w_0+w_a)} · exp(−3w_a·z/(1+z)),  (w_0,w_a) = (−0.840,−0.670)
1/H_0 [Gyr] = 977.8 / H_0[km/s/Mpc]
```

Código: `edad_paper9.py` (numpy + scipy.integrate.quad). Radiación omitida: aporta ~0.01% a t_0.

## Resultados

Comparación honesta: **cada modelo con sus parámetros propios** (la fila híbrida
era solo análisis de sensibilidad — ver punto 2).

| Configuración | Ω_m en E(z) | H_0 | t_0 |
|---|---|---|---|
| Paper 9 actual (s_m como densidad) | 0.160050 | 73.04 | 15.52 (paper) / 15.28 (repro) |
| **SSEE corregida** (parámetros propios) | 0.308881 | 73.04 | **12.78** |
| SSEE corregida, H_global del framework | 0.308881 | 67.962 | **13.73** |
| ΛCDM (parámetros Planck: Ω_m=0.315) | 0.315 | 67.4 | 13.80 |

Sensibilidad (no es comparación entre modelos): a H_0=73.04 y Ω_m=0.308881 fijos,
ΛCDM da 12.80 Gyr vs SSEE 12.78 Gyr — la dinámica (w_0,w_a)≠(−1,0) aporta **−0.02 Gyr**,
irrelevante. Todo el "universo más viejo" era el artefacto de s_m en E(z).

## Lectura honesta

1. **Modelo vs modelo, cada uno con lo suyo:** ΛCDM 13.80 vs SSEE 12.78 (H_0=73.04)
   → SSEE **1.02 Gyr más joven**; con H_global=67.962 → 13.73, indistinguible de ΛCDM.
   El "+1.72 Gyr más viejo" del paper muere en ambos casos.
2. La comparación original del paper mezclaba **dos** inconsistencias: densidad
   equivocada (s_m como Ω_m) **y** H_0 distintos entre modelos (73.04 vs 67.4).
3. Con H_0=73.04, la edad corregida (12.78) deja el argumento de cúmulos globulares
   (≳12.5 Gyr) en un margen de 0.3 Gyr: **el selling point se invierte**.
4. Con H_global=67.962 (el consistente con el framework SSEE) da 13.73 Gyr, sin
   tensión — pero es puro escalado 1/H_0, no física nueva.
5. Brecha repro 15.28 vs 15.52 (1.5%): detalle de integración/redondeo de (w_0,w_a).
   No mueve la conclusión.

## Qué NO está incluido (revisado)

- Radiación Ω_r: ~0.01% en t_0, despreciable. Curvatura: plano, como el paper. Neutrinos: despreciables para t_0.
- Lo que **sí** falta por decidir (no es cálculo, es modelado): **¿qué H_0 debe entrar en la edad del universo dentro del framework SSEE?** El paper usó el local SH0ES (73.04) sin justificarlo; el global del framework es 67.962. La edad es una cantidad global. Decide Mike.

## Cambio neto si se corrige

- t_0: 15.52 → 12.78 Gyr (**−2.74 Gyr**).
- vs ΛCDM: de +1.72 (más viejo) a −1.02 (más joven) a H_0=73.04; a H_0 común la diferencia es −0.02 Gyr (nula).
- El párrafo "more accommodating of the oldest globular clusters" debe reescribirse o retirarse.
