# H_alg vs H_MIRA — ¿Coincidencia o estructura?

**Fecha:** 2026-05-24 PM
**Pregunta:** El modelo tiene dos valores para H₀:
- $H_{\rm alg} = 3(\varphi+\pi)^2 = 67.962$ km/s/Mpc (Type-P algebraico)
- $H_{\rm MIRA} = 67.068$ km/s/Mpc (Planck plik_lite + SSEE bg)

Gap: $\Delta = 0.894$ km/s/Mpc, ratio = 1.01334 (exceso 1.334%).

¿Es coincidencia numerológica o predicción estructural escondida?

---

## TL;DR visual

```
              VALOR ALGEBRAICO            VALOR FÍSICO
              H_alg = 67.962              H_MIRA = 67.068
                  │                            │
                  │                            │
                  ▼                            ▼
            ┌─────────────┐            ┌──────────────┐
            │  numerología │           │  Planck data  │
            │  3(φ+π)²    │ ←Δ=0.894→ │  + SSEE bg    │
            └─────────────┘            └──────────────┘
                  │                            │
                  ▼                            ▼
          "casualidad afortunada"      "valor canónico para
           (Type P coincidence)         predicciones físicas"

                        ⇕  ¿1.334%?

         BÚSQUEDA: ¿alguna identidad SSEE
         da exactamente 1.334%?
                       │
                       ▼
                ❌ NO se encontró

         Mejor candidato: α_K²/(3·MIRA²) = 1.355%
                          (0.02% off, sugiere pista, no prueba)
```

---

## Lo que el barrido encontró

### Test 1 — ¿es 3(φ+π)² la única combinación que llega cerca de H_phys?

| Combinación | Valor | σ vs Planck |
|---|---|---|
| **3(φ+π)² = OMEGA²·3** | **67.962** | **1.66σ** ⭐ |
| AURA·KAL₀·3 | 66.221 | 1.57σ |
| Σ_Sov/π·15 | 68.177 | 2.05σ |
| π²·7 | 69.087 | 3.74σ |
| 4(π+φ)² | 90.616 | 43.6σ |
| (φ+π)·5φ | 38.506 | 52.9σ |

**Hallazgo:** las combinaciones simples de (φ, π) que dan valores en el rango cosmológico [60-80] son **escasas**. La mayoría dan valores muy fuera de rango. Que existan 2-3 a <2σ es algo, pero no extraordinario.

Probabilidad de coincidencia aleatoria (H en [60,80] con tolerancia 1.5%): ~10%. Significativa pero no decisiva.

### Test 2 — ¿el ratio 1.01334 coincide con algo conocido?

| Candidata | Predicción | Δ del exceso |
|---|---|---|
| **α_K²/(3·MIRA²) ≈ α_K²/12** | **1.01355** | **1.6%** ⭐ |
| 1 + (α_K · α_K^UV)/12 | 1.01401 | 5.0% |
| 1 + (π−φ)/Ω³ | 1.01413 | 6.0% |
| 1 + 1/(KAL₀·15) | 1.01207 | 9.5% |
| 1 + 1/π⁴ | 1.01027 | 23% |
| 1 + ε (UV corr.) | 1.00057 | 96% |
| 1 + f_screen | 1.06725 | 405% |

**Hallazgo:** ningún candidato exacto. El más cercano es $\alpha_K^2/(3\,{\rm MIRA}^2) \approx 0.0136$, off por 1.6%. Esto **sugiere** que el exceso podría tener relación con $\alpha_K$ (la kineticidad EFT) pero no es prueba: 1.6% es ruido típico de coincidencias numerológicas.

### Test 3 — ¿la diferencia 0.894 km/s/Mpc match alguna escala física?

Nada significativo. Las escalas naturales del modelo (H·f_screen ≈ 4.6, H·ε ≈ 0.04, H·(MIRA−2) ≈ 0.07) están todas muy lejos de 0.894 km/s/Mpc.

---

## Tres interpretaciones posibles

### Interpretación (A) — Coincidencia numerológica afortunada

H_alg es **solo numerología**: 3(φ+π)² da casualmente un número en el rango cosmológico correcto, a 1.3% del valor físico. Esto es estadísticamente plausible (~10% por azar) pero notable: solo 3 combinaciones simples de (φ,π) caen en rango [60,80].

**Consecuencia:** H_alg no es una predicción del modelo, sino una **pista heurística** que sugiere que la combinación φ,π,Ω, etc., genera escalas del orden correcto. El valor predictivo está en H_MIRA.

### Interpretación (B) — Relación estructural pendiente de identificar

El exceso 1.334% está cerca de α_K²/(3·MIRA²) = 1.355%. Si esta relación fuera exacta, indicaría que **H_alg encapsula a H_MIRA más una corrección dinámica controlada por la kineticidad α_K**.

Pero el desajuste de 1.6% es demasiado grande para reclamar identidad — probablemente sería un argumento esponjoso ante referee. La hipótesis queda **especulativa**.

### Interpretación (C) — Predicción todavía no entendida

El gap es un observable real del modelo que apunta a física que aún no hemos identificado: tal vez un efecto de loop quántico, screening sub-leading, o corrección de back-reaction. Esta interpretación es **prematura** sin teorema que la respalde.

---

## Mi recomendación honesta

**Adoptar la interpretación (A) en las publicaciones** y reportar H_alg como Type-P coincidence:

> *"H_alg = 3(φ+π)² = 67.96 km/s/Mpc coincide numéricamente con el valor físico H_MIRA = 67.07 km/s/Mpc a 1.3%. Esta proximidad es notable dado que pocas combinaciones algebraicas simples de (φ,π) caen en el rango cosmológico; sin embargo, no constituye una predicción estructural derivada y se reporta como una coincidencia numérica (Type-P). El valor predictivo del modelo está en H_MIRA, obtenido autoconsistentemente del fit Planck + background SSEE."*

Esto:
- ✅ Es honesto: no sobre-vende H_alg como predicción
- ✅ Reconoce el patrón sin reclamar mecanismo
- ✅ Posiciona H_MIRA correctamente como cantidad física
- ✅ Sobrevive auditoría hostil

**No adoptar** la interpretación (B) hasta que el gap colapse a <0.1% con alguna identidad limpia, lo cual no sucede con las combinaciones probadas.

---

## Implicación para el Roadmap

Esto resuelve un punto del Roadmap Hacia la Excelencia (CLAUDE.md):

> *"H₀ ≈ 67.96 km/s/Mpc es un éxito numérico, pero falta el teorema que conecte la escala macroscópica (Megaparsec) con unidades fundamentales (escala de Planck)."*

**Conclusión post-investigación:** ese teorema **no existe** en el barrido actual. H_alg es coincidencia numerológica, no teorema dimensional. El reto real sigue siendo derivar **H_MIRA** desde principios fundamentales — H_alg es solo una pista numérica.

Esto baja una de las "vulnerabilidades pendientes" de "esperanza de derivación" a "honestidad declarada".

---

## Archivos

- `src/h_alg_vs_mira_investigation.py` — script ejecutable de auditoría
- Esta nota: `H_ALG_VS_MIRA_INVESTIGATION.md`
