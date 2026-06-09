# Paper 6 — Notes para limpieza honesta

> ⚠️ **SUPERADO (2026-06-04) — snapshot histórico, no leer como estado actual.**
> Anterior a la partícula canónica φ-DM. `m_φ = 5.602 eV` (numerológico) y
> `alpha_WDM fiteado` están **obsoletos**. Canónico: **m_φ = 36.95 eV** =
> Σm_ν·(Ω⁴_DNAV+AURA·KAL), forward-prediction cero fiteo; **alpha = 1.243 Mpc/h**
> es OUTPUT de CLASS, no fiteo. Fuente autoritativa: `VERIFICATION_LEDGER.md` + Paper 6.
> Se conserva como rastro de auditoría.

**Fecha:** 2026-05-25 (sesión nocturna)
**Motivo:** Tras descubrir múltiples capas de fitting numerológico en P6,
documentar los números CANÓNICOS reales del modelo y los ediciones necesarias.

---

## 1. Valores canónicos verificados (CLASS Boltzmann completo)

```
Pipeline: src/p6_canonical_table.py
Asunciones: H₀=67.068, Ωb·h²=0.02242, n_s=0.96556, w0=-0.840, wa=-0.670,
            A_s=2.1e-9 (Planck), CLASS Boltzmann completo, SIN fittear nada

ESCENARIO              Ω_m         σ₈       S₈       fσ₈(z=0.5)
─────────────────────────────────────────────────────────────────
Single-sector          0.160050    0.4423   0.3230   0.2151
Two-sector (MIRA)      0.319928    0.8344   0.8617   0.4870

Tensiones vs observaciones:
                       σ₈ vs KiDS  S₈ vs KiDS   S₈ vs DES
Single-sector          -14.7σ      -22.2σ       -19.0σ
Two-sector             +4.9σ       +4.8σ        +4.5σ
```

**Conclusión:** two-sector es claramente preferible. 5σ vs 20σ.

---

## 2. fσ₈ es donde two-sector brilla limpiamente

```
fσ₈(z=0.5) CANÓNICO:
  Single-sector:  0.215   (catastróficamente bajo)
  Two-sector:     0.487   (excelente, obs ~0.46)
  
SIN fittear nada. Sin alpha_WDM. Sin trucos.
```

Este es el observable donde el modelo dos-sectores TIENE PREDICCIÓN REAL.
Reportar como fortaleza principal de P6.

---

## 3. Descubrimientos críticos del día

### 3.1 σ₈ = 0.737 reportado en P6 es FIT, no predicción

**Evidencia:** `class_ssee/calibrate_wdm_alpha.py` línea 21:
```python
sigma8_target = 0.737   # Paper 6 two-sector
```

El script literalmente busca alpha_WDM tal que CLASS reproduzca 0.737
(= valor observado por KiDS). Resultado: alpha_WDM = 1.6561 h/Mpc fiteado.

### 3.2 m_φ = 5.602 eV es numerológico

Fórmula `m_φ = Σm_ν × H_alg` es dimensionalmente inválida (eV × km/s/Mpc).
Ya revertido en la sesión de la mañana con disclaimer dimensional explícito.

### 3.3 Paper 6 tiene DOS capas de numerología
1. m_φ = 5.602 eV (numerológico, ya documentado en boxed eq)
2. alpha_WDM = 1.6561 h/Mpc (fiteado a observable KiDS, NO declarado)

### 3.4 Verificación scripts (P5/P6) usan aproximación

`ssee_paper6_verification.py` reporta σ₈ = 0.794, pero asume:
- σ₈ ≈ σ₈_LCDM × G_growth (aproximación que asume T_LCDM)
- Esto da número diferente del CLASS Boltzmann real (0.835)
- Para Ω_m=0.32 la aproximación es razonable (~5% diferencia)
- Para Ω_m=0.16 la aproximación falla feo (0.702 vs 0.442 CLASS)

---

## 4. Ediciones requeridas en `manuscript/SSEE_Paper6_phiDM.tex`

### Cambios CRÍTICOS

| Línea aprox | Estado actual | Cambio requerido |
|---|---|---|
| Abstract (~89) | "σ₈=0.794" | Reportar σ₈=0.835 canónico CLASS (raw, no fit) |
| Sección 4.1 (~110) | "σ₈_eff = 0.794 (2.6σ KiDS)" | Actualizar a 0.835 (4.9σ KiDS) raw |
| §6 (origen) | falta declarar alpha_WDM | Agregar declaración: alpha_WDM=1.6561 es parámetro libre fiteado a KiDS, pendiente derivación física vía OP-9/OP-10 |
| §7 (falsificación) | "k_fs = 0.493" | Mantener pero clarificar: depende de m_φ numerológico Y de alpha_WDM fiteado, ambos pendientes |

### Sección nueva sugerida — "Status epistemológico de los parámetros"

```latex
\section{Parámetros libres y derivados — accounting honesto}

El sector φ-DM postula la existencia de una segunda especie de materia
oscura con $\Omega_{\phi DM} = (\text{MIRA}-1)\Omega_{m,dyn} = 0.160$,
predicción algebraica del modelo.

Sin embargo, la caracterización física del sector requiere dos cantidades
adicionales:

\begin{itemize}
  \item \textbf{$m_\phi$}: actualmente fijado por ansatz numerológico
    $m_\phi = \Sigma m_\nu \times H_0^{\rm alg} = 5.602$ eV. Tracked como
    OP-9; derivación física pendiente vía $V(\phi)$ unificado (OP-10).
  
  \item \textbf{$\alpha_{\rm WDM}$}: parámetro de free-streaming suppression.
    En la presente versión, calibrado a $\alpha_{\rm WDM} = 1.6561\,h/$Mpc
    mediante fit a $\sigma_8 = 0.737$ observado por KiDS-1000. \textbf{Este
    es un parámetro libre fiteado, no una predicción.} Derivable desde
    $m_\phi$ y mecanismo de production una vez que OP-9 se resuelva.
\end{itemize}

\textbf{Predicciones genuinas del modelo (sin fittear):}
\begin{align}
  \sigma_8^{\rm RAW} &= 0.835 \quad (+4.9\sigma\text{ vs KiDS}) \\
  S_8^{\rm RAW} &= 0.862 \quad (+4.8\sigma\text{ vs KiDS}, +4.5\sigma\text{ vs DES}) \\
  f\sigma_8(z=0.5)^{\rm RAW} &= 0.487 \quad (\text{excelente match obs}\approx 0.46)
\end{align}

\textbf{Predicción con WDM-suppression (alpha\_WDM como parámetro libre):}
\begin{align}
  \sigma_8^{\rm WDM} &= 0.737 \quad (\text{= KiDS por construcción del fit})
\end{align}

La tensión real del modelo en $\sigma_8/S_8$ es de orden $5\sigma$, comparable
con la tensión $\sigma_8/S_8$ que aflige al $\Lambda$CDM (3$\sigma$ vs KiDS).
La predicción \emph{robusta} de $f\sigma_8$ es independiente de
$\alpha_{\rm WDM}$ y constituye el observable principal donde el sector dos-
sectores muestra su valor sin trucos numerológicos.
```

---

## 5. Reviews pendientes — papers 5, 7, 8, 9

### Paper 5 (IS perturbations)
- ⚠ Usa aproximación σ₈ ≈ σ₈_LCDM × G (no Boltzmann real)
- ⚠ Asume T(k) = T_ΛCDM en transfer function
- Revisar con CLASS canónico para ver si σ₈ = 0.702 es robusto

### Paper 7 (EFT)
- ✓ V₀, αK en unidades ρ_crit=1, invariantes
- ⚠ Verificar que λ = √(3 Ω_m,dyn) es genuinamente derivado, no fiteado
- ⚠ Verificar derivación de αK = 0.4033 desde Bellini-Sawicki

### Paper 8 (Strong gravity)
- ⚠ OP-9 ya marcada como abierta (interpretación m_φ)
- ⚠ Verificar cálculos r_V Vainshtein vs observación SLACS/BELLS (OP-13)
- ⚠ Buscar otros fits ocultos similar a alpha_WDM

### Paper 9 (Hubble tension)
- ✓ f_screen = (π-φ)/Ω² es algebraico puro
- ✓ Cascada H_local lineal verificada
- ⚠ Comparación con EDE/SIDR/late-DE: verificar que la tabla comparativa no esconda asunciones

---

## 6. Resumen del estado P6 tras hoy

```
SÓLIDO (sin trucos):
  ✓ Ω_φDM = 0.160 (predicción algebraica MIRA)
  ✓ fσ₈(z=0.5) = 0.487 vs obs ~0.46 (excelente)
  ✓ Estructura two-sector (mejor que single-sector por mucho)
  ✓ Lagrangiano formal del campo χ (P6 §5)

CON NUMEROLOGÍA (declarar):
  ⚠ m_φ = 5.602 eV (ansatz, OP-9 abierto) — ya declarado
  ⚠ alpha_WDM = 1.6561 h/Mpc (fit, OP-9 abierto) — falta declarar

TENSIÓN REAL (sin fittear):
  ⚠ σ₈/S₈ +4.9σ vs KiDS/DES — comparable con S₈ tension ΛCDM
```

---

## 7. Scripts producidos hoy (referencias)

```
src/p6_canonical_table.py     ← canonical CLASS Boltzmann (este doc)
src/p6_honest_matrix.py        ← separación derivado vs fiteado
src/p6_complete_matrix.py      ← intento CLASS completo (ncdm pendiente)
src/h0_cascade_audit.py        ← cascade H_alg → H_MIRA
src/h_alg_vs_mira_investigation.py  ← coincidencia o estructura
src/op9_physical_formula.py    ← búsqueda m_φ físico
src/op10_principled_search.py  ← búsqueda V(φ) candidatos
src/op9_phi_dm_formula_search.py  ← barrido dimensional combinaciones
```

---

## 8. Próxima sesión — POA

1. **Aplicar las ediciones P6** propuestas en §4
2. **Recompilar PDF P6** y verificar consistencia
3. **Empezar review Paper 5** (similar audit a P6)
4. **Continuar OP-10** (Mecanismos 2-5 que no exploré: Yukawa, misalignment, axion, gauge)

Cuando OP-10 cierre: alpha_WDM emerge naturalmente desde m_φ derivado,
y la tensión σ₈/S₈ se reduce honestamente.
