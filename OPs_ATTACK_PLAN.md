# Plan de Ataque a Open Problems — SSEE-V3.6

**Fecha:** 2026-05-23
**Estado del modelo:** 4 postulados (D, S, M, I); cross-paper consistente; sin contradicciones.
**OPs abiertas:** OP-8, OP-9, OP-10, OP-11, OP-12, OP-13, OP-14.

---

## Mapa de dependencias (cadenas)

```
OP-14 (Σm_ν Type P, offset 22 ad hoc)
    │
    └──► OP-9 (m_φ = Σm_ν · H_alg ansatz)
              │
              └──► OP-12 (relic abundance Ω_φDM)

OP-10 (unify χ into φ via richer V(φ))
    │
    ├──► OP-9  (m_φ = curvature V at minimum) ◄── resuelto por OP-10
    ├──► OP-11 (ξ becomes part of V structure) ◄── resuelto por OP-10
    └──► OP-12 (relic abundance from V dynamics) ◄── resuelto por OP-10

OP-8 (MIRA dynamical mechanism — 7 mech ruled out)
    │
    └──► Papers 5,6,7,8,9 dejan de ser "condicionales a MIRA"

OP-13 (θ_E factor 2 in SLACS/BELLS — observational)
    │
    └──► Self-contained: falsifica o confirma Paper 8
```

## Análisis de impacto vs esfuerzo

| OP | Tipo | Esfuerzo estimado | Impacto si se resuelve | Notas |
|----|------|---|---|---|
| **OP-13** | Observacional | 1-2 sesiones (literatura) | Bajo si confirma, alto si falsifica P8 | Self-contained |
| **OP-14** | Matemático | 2-4 sesiones | Resuelve OP-9 automáticamente | Eslabón más débil |
| **OP-10** | Teórico | 5-10 sesiones | Resuelve OP-9 + OP-11 + OP-12 | Más ambicioso |
| **OP-8**  | Físico | Indef. (7 mech ya fallaron) | Sella todos los condicionales | Más difícil |
| OP-9      | (derivado) | Sale como corolario de OP-14 u OP-10 | — | No atacar directamente |
| OP-11     | (derivado) | Sale como corolario de OP-10 | — | No atacar directamente |
| OP-12     | (derivado) | Sale como corolario de OP-10 | — | No atacar directamente |

## Orden recomendado de ataque

### Fase 1 — Verificación observacional (corto)
**OP-13: θ_E SLACS/BELLS**
- Compilar literatura: SLACS, BELLS, H0LiCOW analyses
- Buscar análisis que reporten desviaciones sistemáticas en θ_E vs GR
- Paper 8 predice factor 2 en régimen disformal-unscreened; si los datos publicados son consistentes con GR puro a unos pocos % → falsificación parcial de P8 (mecanismo de screening k-mouflage debe activarse antes)
- **Salida:** confirmación, calibración del screening, o falsificación del factor 2

### Fase 2 — Derivación matemática (medio)
**OP-14: Σm_ν Type P → derivar R = 4·KAL − 22**

Hipótesis a probar para el offset 22:
1. **Conteo de grados de libertad** — ¿22 = 2 × 11 corresponde a algún conteo de DoF SM (15 fermiones + 4 gauge + 1 Higgs + 2 grav = 22)?
2. **Índice topológico** — ¿22 sale de algún teorema de Atiyah-Singer o invariante de Witten en una variedad SSEE?
3. **Saturación de cota** — análogo a veta-2 (α_sat = √(3/(φ+3π))). Buscar cota física que sature a R = 4·KAL − 22.
4. **Otra forma algebraica** — quizás R no es 4·KAL − 22 sino una expresión más natural; testar 80+ monomios en φ, π, KAL como en `mira_attempts/`.

**Si OP-14 se deriva:**
- $m_\phi = \Sigma m_\nu \cdot H_{\rm alg}$ pasa de ansatz a corolario (OP-9 cae)
- $\Omega_{\phi\rm DM} h^2$ se puede recomputar ab initio (OP-12 cae si se acopla con OP-10)

### Fase 3 — Construcción teórica (largo)
**OP-10: V(φ) unificador DE + DM**

Buscar potencial $V(\phi)$ tal que:
- Plateau DE en $\phi$ tardío (reproduce $w_0 = -0.840$, $w_a = -0.670$)
- Mínimo con curvatura $m_\phi^2$ en $\phi$ post-rolling (reproduce $\phi$-DM con $m_\phi$ algebraico)
- $\alpha$-attractor compatible (preserva $n_s = 1-\varphi^{-7}$, $r = \varphi^{-10}$)

Candidatos a explorar:
- $V(\phi) = V_0 e^{-\lambda\phi} + V_1 (\phi - \phi_0)^2$ — exponencial + harmónico
- $V(\phi) = V_0 \tanh^2(\phi/f)$ con $f = \varphi^k$ — α-attractor + binding mínimo
- Combinaciones $\varphi^n + \pi^m$ con coeficientes algebraicos

**Si OP-10 se construye:** OP-9, OP-11, OP-12 caen juntas. El modelo pasa de 4 postulados a 3 (D, S, I; M sigue siendo postulado hasta resolver OP-8).

### Fase 4 — Mecanismo MIRA (más difícil, dejar para último)
**OP-8: MIRA dynamical mechanism**

7 mecanismos ya descartados (ver `src/mira_attempts/README.md`). Ideas residuales:
- Mecanismo holográfico (relación área/volumen en superficie de Sachs-Wolfe)
- k-mouflage matter-coupled velocity-dependent
- Transición topológica entre dos vacuums

**Status pragmático:** la resolución de OP-8 sería extraordinaria pero el modelo es publicable sin ella (MIRA queda como Postulate M, valor algebraico exacto). Atacar después de OP-10/14 cuando haya más estructura teórica.

## Decisión recomendada

**Empezar por OP-13** (Fase 1, una sesión) porque:
- Es observacional, no requiere construcción teórica
- Resultado binario: confirma o falsifica → cierra puerta
- No bloquea las otras OPs

**Después OP-14** porque:
- Eslabón más débil de la cadena
- Resuelve OP-9 automáticamente
- Es matemático puro (no requiere observación nueva)
- Métodos: scan algebraico tipo `mira_attempts/`, o saturación tipo veta-2

**OP-10 y OP-8 quedan para el largo plazo** — son los retos teóricos profundos. No bloquean publicación: las OPs están abiertas y catalogadas.
