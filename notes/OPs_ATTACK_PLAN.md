# Plan de Ataque a Open Problems — SSEE-V3.6

**Fecha:** 2026-05-23 (revisado tras ataque a OP-14)
**Estado del modelo:** 4 postulados (D, S, M, I); cross-paper consistente; sin contradicciones.
**OPs abiertas:** OP-8, OP-9, OP-10, OP-11, OP-12, OP-14.
**OPs cerradas recientemente:** OP-13 ✅ (Opción A, P8 → DRAFT pending OP-9).

---

## Mapa de dependencias (REVISADO — 2026-05-23 tras ataque a OP-14)

```
ANTES (hipótesis original):

  OP-14 (Σm_ν Type P, offset 22 ad hoc)
      │
      └──► OP-9 (m_φ = Σm_ν · H_alg ansatz)
                │
                └──► OP-12 (relic abundance Ω_φDM)


AHORA (post-ataque OP-14, ver OPEN_PROBLEMS.md §OP-14):

  OP-10  V(φ) unificador DE + DM
     │
     ├──► OP-9   (m_φ = curvatura V en mínimo)
     │      │
     │      └──► OP-14 (Σm_ν = m_φ / H_alg)    ← invertido
     │
     ├──► OP-11  (ξ acoplamiento → funcional de V)
     │
     └──► OP-12  (Ω_φDM h² desde dinámica de V)

  OP-8   MIRA mecanismo (independiente)

  OP-13  ✅ RESUELTO (Opción A: P8 → DRAFT pending OP-9)
```

**Por qué cambió la flecha OP-14 → OP-9:** El ataque directo a OP-14 (script
`src/ssee_op14_neutrino_mass.py`, 2026-05-23) reveló que:
- La forma $\mathcal{R}=4\cdot\text{KAL}-22$ es **estructuralmente frágil**
  (perturbación 10⁻³ en φ,π → 25% drift en Σm_ν).
- Scan algebraico de ~150 monomios no halla identidad exacta (mejor candidato
  $1/[3(\text{KAL}-\varphi)]$ con error −0.27%, sigue siendo aproximación).
- $\Sigma m_\nu$ y $m_\varphi$ comparten **el mismo grado de libertad**:
  derivar uno deriva el otro vía $m_\varphi = \Sigma m_\nu \cdot H_0^{\text{alg}}$.

→ La única salida es derivar $m_\varphi$ desde la curvatura de un $V(\varphi)$
fundamental, que es lo que OP-10 ya pretendía atacar como Fase 3 de medio plazo.
**OP-10 absorbe OP-9, OP-11, OP-12 Y OP-14 simultáneamente.**

## Análisis de impacto vs esfuerzo (REVISADO 2026-05-23)

| OP | Tipo | Esfuerzo estimado | Impacto si se resuelve | Notas |
|----|------|---|---|---|
| ~~**OP-13**~~ | ~~Observacional~~ | ~~1-2 sesiones~~ | ✅ RESUELTO (Opción A, 2026-05-23) | P8 → DRAFT pending OP-9 |
| ~~**OP-14**~~ | ~~Matemático~~ | ~~Atacado 2026-05-23~~ | No derivable directamente | Blocked-by-OP-10 |
| **OP-10** | Teórico | 5-10 sesiones | Resuelve OP-9 + OP-11 + OP-12 + **OP-14** | **Máxima prioridad** ahora |
| **OP-8**  | Físico | Indef. (7 mech ya fallaron) | Sella todos los condicionales | Más difícil; dejar para último |
| OP-9      | (derivado) | Sale como corolario de OP-10 | — | No atacar directamente |
| OP-11     | (derivado) | Sale como corolario de OP-10 | — | No atacar directamente |
| OP-12     | (derivado) | Sale como corolario de OP-10 | — | No atacar directamente |
| OP-14     | (derivado) | Sale como corolario de OP-10 | — | No atacar directamente (ataque 2026-05-23 falló) |

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

## Decisión recomendada (REVISADA 2026-05-23)

**OP-13 ✅ resuelto** (Opción A: P8 → DRAFT pending OP-9, 2026-05-23, commit 89710b3).

**OP-14 ✅ atacado y archivado** (script `src/ssee_op14_neutrino_mass.py`,
2026-05-23): no derivable directamente, blocked-by-OP-10. Ver OPEN_PROBLEMS.md
§OP-14 para evidencia (fragilidad H1, no DoF match H2, scan sin identidad H3).

**Próximo objetivo: OP-10 — Fase 1 (catálogo $V(\varphi)$)**

Plan de ataque OP-10 en 6 fases:

| Fase | Contenido | Esfuerzo | Salida |
|---|---|---|---|
| 1 | Catálogo de 3-5 familias $V(\varphi)$ | 1 sesión | Lista de candidatos viables |
| 2 | Test consistencia background ($w_0, w_a$) | 1-2 sesiones | Subconjunto que pasa Fase 2 |
| 3 | Mínimo y curvatura $m_\varphi$ | 1 sesión | Resuelve OP-9 o descarta |
| 4 | α-attractor compatibility ($n_s$) | 1-2 sesiones | Conexión con OP-2 |
| 5 | $\xi$ acoplamiento como funcional de $V$ | 1 sesión | Resuelve OP-11 |
| 6 | Paper 10b o reescritura P6 §two-sector | 1 sesión | Cierre OP-9, 11, 12, 14 |

Candidatos a explorar en Fase 1:
- $V_1(\varphi) = V_0 e^{-\lambda\varphi/M_{pl}} + \frac{1}{2}m^2(\varphi-\varphi_0)^2$
- $V_2(\varphi) = V_0\tanh^2(\varphi/f) + (m^2/2)\text{sech}^2(\varphi/f)$
- $V_3(\varphi) = M^4[1 - \cos(\varphi/f)] + V_0 e^{-\lambda\varphi}$
- $V_4(\varphi) = V_0(1 - e^{-\sqrt{2/3\alpha}\,\varphi/M_{pl}})^2 + V_{\min}$ (α-attractor puro Starobinsky)

**Constraints duros** (de papers anteriores, no negociables):
- $\alpha = \varphi^4/3$ (Paper 1)
- $\lambda = 0.693$ (Paper 7)
- $V_0$ bloqueado (Paper 7)
- $M = 8.81$ meV (Paper 10)

→ El espacio de búsqueda es pequeño: la forma del mínimo y la curvatura son
casi los únicos grados de libertad libres. Si existe solución, debe ser única.

**OP-8 sigue para el final** — la resolución sería extraordinaria pero el modelo
es publicable sin ella (MIRA como Postulate M con valor algebraico exacto).
