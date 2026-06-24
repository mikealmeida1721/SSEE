# SSEE Constants Audit — Agua vs Aceite

> ⚠️ **SUPERADO (2026-06-04) — snapshot histórico, no leer como estado actual.**
> Anterior a la partícula canónica φ-DM. `m_φ = 5.602 eV` y el framing
> "dimensionalmente inválido" están **obsoletos**. Canónico: **m_φ = 36.95 eV** =
> Σm_ν·(Ω⁴_DNAV+AURA·KAL), forward-prediction dim-consistente. Fuente autoritativa:
> `VERIFICATION_LEDGER.md` §Valores Canónicos + Paper 6. Se conserva como rastro de auditoría.

**Fecha:** 2026-05-25
**Motivo:** Clasificar HONESTAMENTE cada "constante" de SSEE por origen real,
no por la unidad escrita al lado. Evitar usar valores con etiqueta física
cuando en realidad son numerología disfrazada.

**Leyenda:**
- 🟢 **PHYSICAL**: derivado de observación o constante fundamental
- 🔵 **ALGEBRAIC**: combinación pura de (φ, π), adimensional
- 🔴 **NUMEROLOGICAL**: algebraico con unidades que no le corresponden, o con factor ad hoc
- 🟠 **HYBRID**: mezcla — fórmula tiene parte algebraica + parte físico-externa
- ⚪ **UNVERIFIED**: claim de Paper que no he verificado contra su derivación

---

## Tabla maestra

### Fundamentales
| Símbolo | Valor | Origen | Clasificación |
|---|---|---|---|
| M_Pl | 1.221×10²⁸ eV | Constante fundamental | 🟢 PHYSICAL |
| ℏ | 6.582×10⁻¹⁶ eV·s | Constante fundamental | 🟢 PHYSICAL |
| c | 1 (natural) | Constante fundamental | 🟢 PHYSICAL |

### Algebraicas puras (adimensionales)
| Símbolo | Definición | Valor | Clasificación |
|---|---|---|---|
| φ | (1+√5)/2 | 1.6180 | 🔵 ALGEBRAIC |
| π | — | 3.1416 | 🔵 ALGEBRAIC |
| Ω = φ+π | suma | 4.7596 | 🔵 ALGEBRAIC |
| β = (φ+π)/2 | media | 2.3798 | 🔵 ALGEBRAIC |
| KAL₀ = β+π | combinación | 5.5214 | 🔵 ALGEBRAIC |
| P_sc, K_v, T_r, M_v | combinaciones de φ,π | varios | 🔵 ALGEBRAIC |
| MIRA = (3φ+π)/4 | postulado auxiliar OP-8 | 1.9989 | 🔵 ALGEBRAIC (mecanismo: OP-8 abierto) |
| AURA = 2·MIRA | derivado | 3.9978 | 🔵 ALGEBRAIC |
| w₀ = -T_r/M_v | algebraico | -0.840 | 🔵 ALGEBRAIC |
| wₐ = -P_sc/K_v | algebraico | -0.670 | 🔵 ALGEBRAIC |
| Ω_DE = T_r/M_v | algebraico | 0.840 | 🔵 ALGEBRAIC |
| Ω_m,dyn = 1+w₀ | algebraico | 0.160 | 🔵 ALGEBRAIC |
| n_s = 1−φ⁻⁷ | algebraico (Paper 1, OP-2 partial) | 0.9656 | 🔵 ALGEBRAIC |
| r_tensor = φ⁻¹⁰ | algebraico | 0.00813 | 🔵 ALGEBRAIC |
| f_screen = (π−φ)/Ω² | algebraico (Paper 9) | 0.06725 | 🔵 ALGEBRAIC |
| α_K = 0.4033 | calculado en unidades ρ_crit=1 (Paper 7) | 0.4033 | 🔵 ALGEBRAIC |
| α_K^full = 0.41691 | UV correction adimensional (P10) | 0.41691 | 🔵 ALGEBRAIC |
| ε = X²/M⁴ | adimensional en unidades ρ_crit=1 | 5.7×10⁻⁴ | 🔵 ALGEBRAIC |
| Ω_m,CMB = MIRA·Ω_m,dyn | algebraico (predice Planck) | 0.31983 | 🔵 ALGEBRAIC |
| Ω_b·h² = (π−φ)/(3·Ω²) | algebraico (Paper 1, OP-1 partial) | 0.02242 | 🔵 ALGEBRAIC (vs 0.02237 obs) |

### Físicas (calibradas con datos)
| Símbolo | Valor | Origen | Clasificación |
|---|---|---|---|
| H_MIRA | 67.068 km/s/Mpc | Cobaya plik_lite + SSEE bg | 🟢 PHYSICAL |
| ℏH_MIRA | 1.43×10⁻³³ eV | derivado de H_MIRA | 🟢 PHYSICAL |
| ρ_crit(H_MIRA) | 3.65×10⁻¹¹ eV⁴ | = 3H²M_Pl²/(8π) con H_MIRA | 🟢 PHYSICAL |
| ρ_DE = Ω_DE·ρ_crit | 0.840 × ρ_crit | algebraico × físico | 🟢 PHYSICAL (efectivo) |

### 🔴 NUMEROLÓGICAS (algebraico vestido de físico)
| Símbolo | Claim Paper | Realidad | Por qué es numerológico |
|---|---|---|---|
| **H_alg** = 3(φ+π)² | "67.962 km/s/Mpc" (Paper 1, 4) | 🔴 algebraico × unidad arbitraria | El "km/s/Mpc" es etiqueta puesta a mano. Por eso H_alg/H_MIRA tiene 1.34% de gap sin explicación. **Ya documentado como Type-P coincidence hoy** |

### 🟠 HÍBRIDAS (verificadas hoy)
| Símbolo | Valor reportado | Lo que descubrí | Status |
|---|---|---|---|
| **Λ_SSEE = M** | 8.81 meV (Paper 10) | 🟠 La fórmula M⁴=5φ⁸ρ es dim-válida. PERO el valor numérico 8.81 meV usa **ρ_Λ canónico ΛCDM ≈ (2.25 meV)⁴**, no ρ_crit puro SSEE. Con ρ_crit puro daría 9.74 meV (h=0.6796) o 9.65 meV (h=0.6707). | Footnote agregada hoy. Convención calibrada a Λ-ΛCDM externamente |
| **Σm_ν^active** | 0.0824 eV (Paper 4, OP-14) | 🟠 Fórmula: R·Ω_b·h²·94.07/(τ_Π·H₀). Pero R = 4·KAL−22 con "**−22 ad hoc**" admitido en P4 L707-711 como "calibrated ansatz". El paper lo clasifica como Type-P explícitamente. | OP-14 abierto. **NO usar como input físico limpio** |
| **m_φ = 5.602 eV** | (Paper 6) | 🔴 Ya revertido hoy: fórmula `Σm_ν × H_alg` es dimensionalmente inválida. | Revertido + disclaimer dimensional |

---

## Lo que esto significa para construir fórmulas

### Inputs limpios disponibles (🟢 PHYSICAL puro)
```
M_Pl ≈ 1.22×10²⁸ eV
ℏH_MIRA ≈ 1.43×10⁻³³ eV
ρ_crit(H_MIRA) ≈ 3.65×10⁻¹¹ eV⁴
ρ_DE(H_MIRA) ≈ 3.06×10⁻¹¹ eV⁴   (= Ω_DE × ρ_crit; Ω_DE algebraico)
```

### Inputs algebraicos limpios (🔵 — adimensionales)
```
φ, π, Ω, β, KAL₀, MIRA, w₀, wₐ, Ω_DE, Ω_m,dyn,
n_s, f_screen, α_K, α_K^full, ε, Ω_m,CMB, Ω_b·h²
```

### Inputs CONTAMINADOS (NO usar sin clarificar)
```
🔴 H_alg = 67.962 km/s/Mpc      → solo etiqueta. Usar H_MIRA físico.
🟠 Λ_SSEE = 8.81 meV            → calibrado externamente. Usar ρ_crit o ρ_DE directo.
🟠 Σm_ν = 0.0824 eV             → Type-P con −22 ad hoc. NO usar hasta cerrar OP-14.
```

---

## Implicación para OP-9 (m_φ del φ-DM)

**Lo que estaba haciendo MAL:**  
Construir m_φ usando Λ_SSEE (🟠 calibrado externamente) y Σm_ν (🟠 con offset ad hoc). Aunque la fórmula `Ω × Σm_ν²/Λ_SSEE` cierre dimensionalmente, los inputs arrastran numerología.

**Lo que se PUEDE construir con inputs limpios:**  
Combinaciones de {M_Pl, ℏH_MIRA, ρ_crit, ρ_DE} × {φ, π, Ω, KAL₀, ...} con factores algebraicos puros.

**Lo que ya sabemos del barrido op9_phi_dm_formula_search.py (mañana):**  
NINGUNA combinación dimensional simple de {M_Pl, ℏH_MIRA, Σm_ν, Λ_SSEE} dio m_φ en [1, 30] eV con mecanismo físico claro. La mejor candidata (φ × Σm_ν² × Λ_SSEE / (ℏH·M_Pl) = 5.54 eV) usaba **3 inputs híbridos**.

Si restringimos a SOLO inputs 🟢 + 🔵 limpios, hay que repetir el barrido sin Σm_ν ni Λ_SSEE como inputs. **Pendiente verificar si existe alguna combinación.**

---

## Acciones derivadas de esta auditoría

1. **NO usar Λ_SSEE = 8.81 meV** como anchor físico en fórmulas para OP-9 sin antes:
   - Decidir si usamos ρ_crit puro SSEE → M ≈ 9.65 meV bajo H_MIRA
   - O reconocer que P10 usa Λ-canónico ΛCDM como calibrador externo

2. **NO usar Σm_ν = 0.0824 eV** como input físico para construir m_φ en OP-9:
   - OP-14 está abierto. El offset −22 es ad hoc.
   - Si Σm_ν entra en la fórmula de m_φ, m_φ hereda la contaminación.

3. **Rehacer op9_phi_dm_formula_search.py** usando SOLO {M_Pl, ℏH_MIRA, ρ_crit_MIRA, ρ_DE_MIRA} + factores algebraicos. Reportar si existe alguna fórmula limpia que dé m_φ ∈ [1, 30] eV.

4. **Si NO existe formula limpia** que dé m_φ ∈ eV-scale: aceptar honestamente que SSEE actual NO tiene mecanismo físico para producir m_φ ∈ eV. OP-9 queda abierto hasta resolver OP-14 y/o cerrar la convención de Λ_SSEE.

5. **Si SÍ existe**: esa es la predicción física que buscábamos.

---

## Disculpa metodológica

Esta auditoría debió ser el PRIMER paso de OP-9, no el último. Estaba construyendo fórmulas con ingredientes contaminados y reportándolas como "físicas" porque la fórmula cerraba dimensionalmente. Mike tuvo que insistir tres veces para que finalmente verificara los inputs. El error fue mío.

**Regla para sesiones futuras:** antes de usar cualquier "constante" SSEE en una derivación, verificar su origen contra esta tabla. Si está marcada 🟠 o 🔴, NO usarla sin clarificar primero.
