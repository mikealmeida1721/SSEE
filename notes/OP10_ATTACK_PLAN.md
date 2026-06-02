# OP-10 Attack Plan — V(φ) unificado DE+DM

**Objetivo:** Construir V(φ) cuyos parámetros estén DERIVADOS de la estructura SSEE,
no inventados. Si tiene mínimo con curvatura V''(φ_min), entonces:

```
m_φ² = V''(φ_min)   ← predicción derivada limpia (cumple Regla 3)
```

Si esta predicción cae en rango eV: P6 se rescata, alpha_WDM se deriva, σ₈/S₈ tensión se reduce honestamente.
Si no cae en eV: aceptar lo que dé SSEE, ajustar fenomenología P6.

---

## Lo que YA exploramos (sesión 2026-05-25)

```
src/op10_principled_search.py — 3 familias probadas:
─────────────────────────────────────────────────────
C1: ΔV = c·V₀·(αφ)²·exp(-αφ)         → m_φ ≈ 1-2 μeV   (loop corr.)
C2: V = V₀·tanh²(βφ/√6) α-attractor   → m_φ ≈ 1.7 μeV
C3: V = -μ²φ² + λφ⁴/4 Higgs-like     → m_φ ≈ 10⁻³³ eV  (Hubble seesaw)

TODAS dan μeV o menos. NINGUNA llega a eV con principios estrictos.
```

---

## Mecanismos que FALTAN explorar

### Mecanismo 2 — Yukawa neutrino-φDM
```
L_int = y · φ · ν̄ν

m_φ² emerge de loop con neutrino:
  m_φ² ~ y² · Σm_ν² · log(Λ/m_ν)

Si y derivado de (φ,π): podría dar eV
Inputs limpios: y, Σm_ν (Type-P), Λ (UV cutoff)

Riesgo: Σm_ν está contaminada (OP-14 con -22 ad hoc)
Salida: limpiar OP-14 primero, o usar mecanismo que no use Σm_ν
```

### Mecanismo 3 — Axion-like misalignment
```
L = (Λ_ax⁴/f²) · cos(φ/f)
V(φ) = Λ_ax⁴ · (1 - cos(φ/f))

m_φ = Λ_ax² / f

Si Λ_ax y f derivados de SSEE:
  Λ_ax ~ Λ_SSEE = M ≈ 8.81 meV (con disclaimer convención)
  f ~ ? (la incógnita crítica)
  
Si f ~ M_Pl: m_φ ~ (8.81 meV)²/M_Pl ~ 6×10⁻³³ eV (muy chico)
Si f ~ √(Λ·M_Pl): m_φ ~ Λ_ax ~ meV
Si f ~ Σm_ν: depende de OP-14
Si f ~ algebraico·ρ_crit^(1/4): cálculo pendiente

→ Necesita f DERIVADO desde primeros principios SSEE
```

### Mecanismo 4 — Gravitational production
```
m_φ desde producción durante inflación:
  m_φ² ~ ξ · H_inflation²
  
Donde:
  ξ = acoplamiento no-mínimo a curvatura (OP-11 libre)
  H_inflation = escala Hubble durante inflación (derivable)

Si ξ ~ O(1): m_φ ~ H_inflation
SSEE α-attractor con N_*=58 e-folds → H_inf ~ 10¹³ GeV
  m_φ ~ 10¹³ GeV = 10²² eV (muy alto)

Necesitaríamos ξ << 1, pero ξ debe ser DERIVADO no fitted
```

### Mecanismo 5 — Symmetry breaking (Higgs-like 2do orden)
```
V(φ) = ((φ²-v²)²/4) · λ

Con v ~ algebraico·M_Pl o v ~ √(M_Pl·ℏH₀):
  v = √(M_Pl·ℏH₀) ≈ 4 meV (cosmo seesaw)
  λ = adimensional (φ,π)
  
m_φ² = 2λv² → m_φ ~ √λ · v
       ~ √λ · 4 meV

Para llegar a eV: λ ~ 60000 (no físico)
Para m_φ ~ meV: λ ~ 1 (natural)
→ Predice meV no eV
```

### Mecanismo 6 — Two-field hierarchy
```
Si φ-DE y χ-DM son DOS campos SEPARADOS:
  L = K_φ(X_φ) - V_φ(φ) + K_χ(X_χ) - V_χ(χ) + L_int(φ, χ)

V_φ(φ) ya está en P7
V_χ(χ) puede tener forma diferente, escala diferente

Si V_χ usa escala UV diferente: m_χ podría ser eV-scale natural
PERO: introduce dos campos, conteo de parámetros sube
→ Viola minimalismo SSEE — pero podría ser la realidad

Test: ¿puede SSEE postular V_χ con todos los parámetros derivados
desde (φ,π) y M_Pl/ρ_crit?
```

---

## Estrategia sistemática para mañana

### Fase 1 — Inventario de escalas eV-naturales

Verificar qué combinaciones de escalas físicas SSEE limpias
{M_Pl, ℏH_MIRA, ρ_crit^(1/4), ρ_DE^(1/4)} con factores algebraicos
{φ, π, KAL₀, MIRA, T_r, M_v, ...} dan exactamente eV scale.

Si ninguna combinación natural da eV: SSEE V3.6 NO predice
eV-scale DM, OP-10 cierra negativamente y P6 debe reformularse
a otro régimen (μeV fuzzy o meV milli-eV).

### Fase 2 — Test cada mecanismo 2-6

Para cada uno:
1. Aplicar 3 principios (checklist)
2. Identificar parámetros libres
3. Buscar SI alguno puede fijarse algebraicamente
4. Computar m_φ resultante
5. Reportar honestamente

### Fase 3 — Decisión

```
Si Fase 2 encuentra un mecanismo con m_φ ∈ eV derivado limpiamente:
  → OP-9 cerrada, P6 reescrita con derivación física
  → alpha_WDM emerge desde m_φ, no fiteado
  → σ₈/S₈ tensión reducible honestamente

Si Fase 2 no encuentra ninguno:
  → SSEE V3.6 confirma: m_φ está en μeV/meV scale
  → P6 debe reformularse: φ-DM es fuzzy o milli-eV
  → k_fs es muy chico, no afecta σ₈/S₈ via free-streaming
  → Tensión σ₈/S₈ 5σ aceptada como prediction honesta
  → O P6 elimina del modelo principal
```

---

## El criterio claro

```
Para que OP-10 "rescate" P6 sin trucos:

  1. V(φ) con minimo derivable (no inventado)
  2. m_φ² = V''(φ_min) ∈ [1, 30] eV
  3. Todos los parámetros desde {(φ,π), M_Pl, ρ_crit, ℏH_MIRA}
  4. NO usar Σm_ν (contaminado OP-14)
  5. NO usar Λ_SSEE=8.81 meV (convención canónica-ΛCDM, ya marcado)

Si esto existe: rescatamos P6 limpio.
Si no existe: aceptar prediction real ~5σ tensión, declarar honestamente.
```

---

## Conexión con otras OPs

```
OP-10 cierra ↘
              → resuelve OP-9 (m_φ desde V'')
              → resuelve OP-11 (ξ desde V(φ) forma)
              → aporta predicción k_fs física
              → reduce dependencia de alpha_WDM
              
OP-10 abierta → P6 declarado fenomenológico
              → alpha_WDM declarado parámetro libre
              → σ₈/S₈ con tensión real ~5σ
```

---

## Próxima sesión — primer paso

1. Construir `src/op10_systematic_search.py`:
   - Inventario completo de escalas y combinaciones
   - Test cada mecanismo 2-6
   - Reporte estructurado: m_φ, parámetros usados, status principios

2. Decidir basado en resultados

3. Aplicar conclusiones a P6 cleanup (ver P6_CLEANUP_NOTES.md)

Ver también: P6_CLEANUP_NOTES.md, SSEE_CONSTANTS_AUDIT.md, OPEN_PROBLEMS.md
