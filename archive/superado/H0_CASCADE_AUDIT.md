# H₀ Cascade Audit — SSEE-V3.6 (REVISADO)

> ⚠️ **SUPERADO (2026-06-04) — snapshot histórico, no leer como estado actual.**
> La parte de m_φ (`5.602 eV con disclaimer dimensional`) está **obsoleta**: el
> canónico es **m_φ = 36.95 eV** = Σm_ν·(Ω⁴_DNAV+AURA·KAL), forward-prediction
> dim-consistente. La cascada H₀ (P9 IR 71.90 / P10 UV 72.077) sigue vigente.
> Fuente autoritativa: `VERIFICATION_LEDGER.md` §Valores Canónicos. Rastro de auditoría.

**Fecha:** 2026-05-24 PM
**Motivo:** Tras revertir P6 (m_φ) por inconsistencia dimensional, auditoría bloque
por bloque para asegurar que NINGUNA otra cascada superficial ocultara errores.

**Conclusión:** Tras lectura profunda de las derivaciones, el modelo es notablemente
robusto. P5, P7, P10 están **diseñados en unidades naturales (ρ_crit = 1)** donde
los observables son ratios adimensionales — invariantes bajo cualquier H₀.

---

## TL;DR — el panorama real

```
              H_alg = 67.962                  H_MIRA = 67.068
              (Type-P)                        (físico Planck+SSEE bg)
                  │                                  │
                  └─────────── Δ = −1.32% ───────────┘

  IMPACTO REAL DE LA CASCADA H_alg → H_MIRA:

  ✅ INVARIANTE POR DISEÑO       ✅ YA CASCADADO          ❌ NO VÁLIDO
  ─────────────────────────      ──────────────────       ─────────────
  • Postulados D, S, M, I        • H₀ posterior (P2)      • m_φ (P6):
  • n_s, λ, f_screen             • H_local (P9)             fórmula
  • α_K, ε, M⁴ (unidades         • H_UV (P10)               dimensional-
    naturales P7/P10)            • Ω_b·h² posterior         mente invá-
  • τ_Π·H₀, γ_IS (P5)              (P2 MCMC)                lida → re-
  • r_V en metros (P8)                                      vertido a
                                                            5.602 eV
                                                            con H_alg
                                                            + disclaimer

  ⚠ COSMÉTICO / PEQUEÑO   ⚠ POR VERIFICAR
  ──────────────────────  ──────────────────
  • M en meV (P10):       • Σm_ν (P4): podría
    8.81 con convención     cambiar ~1% si re-
    ρ_Λ ≈ (2.25 meV)⁴;     computas con Ω_b·h²
    el OBSERVABLE H_UV     del MCMC MIRA (cos-
    no depende de M(meV)   mético — P6 numerol.
                           ya revertido)
                          • r_d (P3): Cobaya output
                            verificar (likely sin
                            cambio)
```

---

## Diagrama del flujo correcto

```
                          H_alg = 67.962
                                  │
       ┌──────────────────────────┼──────────────────────────┐
       │                          │                          │
       ▼                          ▼                          ▼
  ADIMENSIONALES            LINEALES en H₀            COSMÉTICOS
  (rho_crit=1)              (cascada trivial)         (unidades meV/eV)

  α_K=0.4033              H_local=H₀/(1−f_screen)    M(meV) según
  f_screen=0.06725                                    convención de ρ
  α_K^full=0.41691        H_UV=H₀/(1−f_screen^UV)    elegida (ρ_crit
  ε=5.7e-4                                            puro vs ρ_Λ
  V₀=0.840                                            canónico)
  M⁴=234.89 (norm.)
  γ_IS=0.554
  z̃, τ_Π·H₀
  KAL₀=5.5214
       │
       ▼
  NO CAMBIAN al cascadar H₀
  (son ratios adimensionales)
```

---

## La razón física por la que casi nada cambia

P5, P7 y P10 fueron construidos en **unidades naturales donde ρ_crit ≡ 1**.
Cita textual de P7 (línea 166):

> *"$M^{4} \equiv \rho_{\rm crit}$ as a unit normalisation ($M^4 = 1$ in units of $\rho_{\rm crit}$)"*

En ese sistema:
- $V_0 = \Omega_{DE} \cdot \rho_{\rm crit} = 0.840 \cdot 1 = 0.840$ ← adimensional
- $M^4 = 5\phi^8 \cdot \rho_{\rm crit} = 234.89 \cdot 1 = 234.89$ ← adimensional
- $\epsilon = X^2/M^4 = 5.7 \times 10^{-4}$ ← adimensional
- $\alpha_K^{\rm full} = \alpha_K^{\rm IR} + 24X^2/M^4 = 0.41691$ ← adimensional
- $f_{\rm screen}^{\rm UV} = \alpha_K^{\rm full}/(3\,{\rm MIRA}) = 0.06952$ ← adimensional

Todo es adimensional → **invariante bajo qué H₀ uses**.

El ÚNICO momento donde H₀ entra es en la cascada lineal final:
```
H_local^UV = H_0 / (1 − f_screen^UV)
       └────────┘   └────────────┘
       depende      adimensional, invariante
       de qué H₀    bajo cascada H₀
       elijas
```

Esta es la razón por la que la cascada AM (73.040 → 72.077) **es correcta** —
solo se necesita cambiar el H₀ de entrada.

---

## El "issue" de M = 8.81 meV es solo cosmético

El manuscript dice "$M = \phi^2 \times 5^{1/4} \times \rho_{\rm crit}^{1/4} \approx 8.81$ meV" con $\rho_{\rm crit}^{1/4} \approx 2.25$ meV.

Pero $(2.25\text{ meV})^4 = 2.56 \times 10^{-11}$ eV⁴, que NO es $\rho_{\rm crit}$
puro para $h = 0.6796$ ($3.74 \times 10^{-11}$ eV⁴). Es más cercano a $\rho_\Lambda$
canónico de ΛCDM.

**Convenciones posibles para reportar M numéricamente:**
| Convención | ρ usado | M [meV] |
|---|---|---|
| ρ_crit puro (h=0.6796, H_alg) | 3.74×10⁻¹¹ | 9.74 |
| ρ_crit puro (h=0.67068, H_MIRA) | 3.65×10⁻¹¹ | 9.65 |
| ρ_DE-SSEE = 0.840×ρ_crit (H_alg) | 3.14×10⁻¹¹ | 9.27 |
| ρ_DE-SSEE = 0.840×ρ_crit (H_MIRA) | 3.06×10⁻¹¹ | 9.21 |
| ρ_Λ canónico ΛCDM ≈ (2.25 meV)⁴ | 2.56×10⁻¹¹ | **8.81** ← lo que P10 usa |

**TODOS estos valores predicen el mismo H_UV** (porque H_UV solo depende de
α_K^full, que es adimensional). El número en meV es solo cómo presentamos M
en unidades convencionales.

**Acción editorial recomendada:** clarificar en P10 §3.2 que el "2.25 meV"
es la escala de Λ canónica, no ρ_crit puro. Pequeño cambio de redacción,
no afecta ningún resultado.

---

## Tabla de acciones REVISADA

| Paper | Acción real | Impacto |
|---|---|---|
| **P6** | YA REVERTIDO (m_φ → 5.602 con H_alg + disclaimer dimensional) | ✅ hecho |
| **P9** | YA CASCADADO (H_local = 71.90, 1.10σ SH0ES) | ✅ hecho |
| **P10** | YA CASCADADO (H_UV = 72.077, 0.93σ SH0ES) — cascada lineal correcta | ✅ hecho |
| **P10** | Editorial: clarificar convención de ρ en M(meV) (8.81 ↔ 9.74 según ρ_crit puro/ρ_Λ) | ⚠ pendiente menor |
| **P4** | Recomputar Σm_ν si quieres usar Ω_b·h² del MCMC MIRA (cambio ~+1%) | ⚠ cosmético |
| **P3** | Verificar r_d output de Cobaya con MCMC MIRA-recalibrado | ⚠ likely sin cambio |

**Casos sin acción (mi alarma inicial era falsa, mea culpa):**
- P10 αK_full: ✅ adimensional, invariante
- P10 H_UV: ✅ cascada lineal AM ya correcta
- P7 V₀: ✅ unidades naturales, invariante
- P5 γ_IS, z̃, τ_Π·H₀: ✅ adimensionales
- P8 r_V: ✅ en metros vía G físico, invariante
- Postulados D, S, M, I: ✅ adimensionales

---

## Implicación para OP-9 (búsqueda de m_φ físico)

**Puedes usar Λ_SSEE = M en OP-9 sin miedo a cascadas escondidas.** El valor
de M en meV es cosmético — lo que importa es la ESCALA RELATIVA respecto a
otras escalas físicas (M_Pl, ℏH₀, Σm_ν).

Las rutas físicas para m_φ que propuse antes siguen válidas:
- Misalignment (m_φ = Λ²/f_φ con Λ = M)
- Producción gravitacional inflacionaria
- Curvatura V''(φ_min) de potencial unificado

---

## Lección epistemológica

Este audit confirma algo importante sobre el modelo SSEE:
**la mayoría de los observables están en unidades naturales/adimensionales.**
Eso significa que las cascadas H_alg ↔ H_MIRA son sólo cambios cosméticos en
la presentación, NO en la física. El único lugar donde la cascada importa
físicamente es:
- Donde H₀ aparece como anchor dimensional (H_local, H_UV)
- Donde una fórmula está mal construida dimensionalmente (m_φ P6, ya revertido)

Eso da confianza para futuras decisiones de cascada: si una fórmula está en
unidades naturales (ρ_crit=1), es invariante. Si tiene H₀ explícito como
dimensional, hay que cascadar linealmente.

---

## Anexo: Cómo correr la auditoría

```bash
python src/h0_cascade_audit.py
```

Devuelve la tabla completa con clasificaciones por paper, cantidad y dimensión.
