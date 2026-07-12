# NUMEROLOGY_AUDIT — ¿Es SSEE numerología? (respuesta pre-registrada)

> **Propósito.** Este documento existe para que la pregunta *"¿tu modelo no será
> numerología?"* NO se re-litigue en cada auditoría. Operacionaliza la acusación en
> criterios falsables, los aplica al estado real del repo, y da un veredicto con su
> evidencia. Si un revisor (o un Claude futuro) plantea la duda, **la respuesta es este
> archivo**, no una nueva derivación desde cero. Actualizar SOLO si cambia la evidencia
> (nuevos datos, nueva derivación, o un fallo encontrado).
>
> Autor del reto: Mike Almeida — *"pruébalo; si puedes demostrar que es numerología, lo
> dejamos ahí. Pero necesito pruebas."* Este es el intento honesto de probarlo.

---

## 1. Definición operacional (no vibes)

"Numerología" no es "usa números bonitos como φ y π". Es un **modo de fallo epistémico**
con propiedades verificables. Siguiendo a Popper, su rasgo **definitorio** es la
**infalsabilidad + acomodo post-hoc**. Operacionalizado en 4 tests:

| # | Criterio de numerología | Cómo se prueba |
|---|---|---|
| **N1** | **Sobreajuste / look-elsewhere grande.** El diccionario es tan flexible que acertar los blancos observados no sorprende. | Contar cuántas expresiones simples de (φ,π) caen en el blanco. Si son muchas → azar. |
| **N2** | **Post-hoc.** Todas las "predicciones" son retro-dicciones de valores YA conocidos. | ¿Existe alguna predicción sobre un dato que AÚN no se ha medido? Si no → post-hoc. |
| **N3** | **Sin anclaje dimensional.** Un número puro se iguala a una cantidad con unidades sin escala física. | Revisar cada cantidad dimensional: ¿pasa por una escala física, o es número=unidad por fiat? |
| **N4** | **Infalsable.** No hay ninguna medición futura que pueda matar el modelo. | ¿Existe un criterio de falsación pre-comprometido? Si no → no es ciencia. |

**Para PROBAR que SSEE es numerología basta con demostrar UNO de N1–N4.** Vamos.

---

## 2. El test aplicado (referí hostil)

### N1 — Look-elsewhere: ¿el diccionario es demasiado flexible?

**Contra la acusación:**
- **Conteo individual:** cada blanco (w₀, wₐ, …) es raro entre las expresiones simples de
  (φ,π): **1/317** razones distintas (`src/estadistica/look_elsewhere_full.py`).
- **Conteo CONJUNTO (el que importa):** SSEE NO ajusta w₀ y wₐ por separado. El **mismo
  esqueleto** (Ω=φ+π y sus múltiplos enteros) debe dar **los DOS a la vez**: −0.840 y
  −0.670. El conteo de bases que aciertan AMBOS con la misma base es mucho menor
  (`look_elsewhere_joint.py`). Un ajustador con 2 perillas independientes NO tiene esta
  rigidez; SSEE la tiene por construcción.
- **Diccionario ACOTADO:** ley de no-auto-suma + techo de cierre M_v=3Ω, alcanzado por
  **exactamente 20 rutas** (conjunto finito, enumerable — Paper 4 App. independencia). No
  es un espacio abierto donde "todo cabe".

**A favor de la acusación (dónde el referí anota — y el modelo lo CONCEDE):**
- El coeficiente de masa de la partícula φ-DM (594.28, OP-9) tiene look-elsewhere
  **dependiente de la gramática**: 1/537 permisiva, 1/192 "volúmenes", ~1/16 estricta. Para
  ESA cantidad la selección **no es estadísticamente abrumadora**. Por eso **OP-9 está
  etiquetado ABIERTO**, no cerrado (ver §OP-9 en OPEN_PROBLEMS.md, mitad UV congelada).

**Veredicto N1:** No se puede probar numerología en el **núcleo** (sector DE: w₀,wₐ,H₀
con look-elsewhere conjunto restrictivo). El coeficiente de OP-9 es débilmente seleccionado
— **y el modelo lo declara**, no lo esconde. → N1 **no demostrado** para el núcleo.

### N2 — ¿Todo es post-hoc?

- **w₀wₐ:** SSEE **no reclama prioridad temporal** sobre DESI/Planck públicos (Paper 3
  §B1 explícito). El match a w₀≈−0.84, wₐ≈−0.6 es **consistencia**, no pre-dicción probada.
  Honesto: aquí NO hay crédito de predicción.
- **Predicciones genuinamente FORWARD (dato aún inexistente):**
  - **k_fs = 0.754 h/Mpc** (imprint de m_φ en P(k)) → **DESI Y3 / Euclid, 2026–2028**.
  - **r = φ⁻¹⁰** (tensor-a-escalar) → CMB-S4 / LiteBIRD.
  - Estas **no pueden ser post-hoc**: el dato no existe todavía. Un sistema numerológico
    **no apuesta** sobre cantidades no medidas.

**Veredicto N2:** No se puede probar N2. El modelo tiene predicciones forward sobre datos
no medidos. → **Evidencia anti-numerología más fuerte.**

### N3 — ¿Número puro = cantidad dimensional por fiat?

- **Auto-corrección documentada:** la vieja m_φ = Σm_ν × H_alg (eV × km/s/Mpc) **era**
  dimensionalmente inválida — y el modelo la **detectó y RETIRÓ** (reframe). La actual
  m_φ = Σm_ν[eV] × (número puro φ,π) → eV. Válida.
- **Masas y densidades** pasan por una escala física (Σm_ν, ρ_crit), no número=unidad.
- **Punto débil genuino (referí anota):** H₀ = 3(φ+π)² "= 67.96 km/s/Mpc". Defensa del
  modelo: las unidades entran por la inversión SH0ES–f_screen (dato dimensional), no por
  fiat. Pero **el teorema del puente dimensional Planck→Mpc sigue FALTANDO** (roadmap #1).
  Etiquetado Type-P, abierto.

**Veredicto N3:** No se puede probar N3 en general. Un punto débil real (puente H₀),
**etiquetado**, no oculto. → N3 **no demostrado**.

### N4 — ¿Es infalsable? (el criterio definitorio)

- SSEE tiene **tabla de falsación pre-comprometida** (Paper 1): k_fs, r, S₈, evolución
  w₀wₐ, todas con criterio de muerte cuantitativo.
- Una medición de 2027 (k_fs de DESI Y3) puede **matar el modelo**.

**Veredicto N4:** **Imposible probar N4.** El modelo es falsable por diseño. Esto es lo
opuesto exacto a la numerología.

---

## 3. Veredicto: la prueba que se pidió

> **NO se puede probar que SSEE es numerología.**

El intento de prueba **falla en el criterio definitorio (N4)**: la numerología es
infalsable y post-hoc; SSEE **apuesta predicciones falsables (k_fs=0.754, r=φ⁻¹⁰) sobre
datos que aún no existen**. Un sistema que una medición de 2027 puede matar **no es
numerología** — es una **hipótesis física** (posiblemente falsa, pero hipótesis).

**El rasgo decisivo:** la numerología **esconde su libertad** y reclama completitud. SSEE
**enumera su libertad** — OPEN_PROBLEMS.md cataloga 23 problemas abiertos; los scripts de
look-elsewhere **reportan los conteos desfavorables** (1/537, dependencia de gramática) en
vez de ocultarlos. **Catalogar tus propias debilidades es la firma anti-numerología.**

---

## 4. La otra mitad honesta: qué NO está probado

"No es numerología" **≠** "es correcto". Ubicación en el espectro:

```
  NUMEROLOGÍA ─────────[ SSEE: falsable, rígido, forward ]──────────── FÍSICA CONFIRMADA
  (infalsable,          ↑ está aquí                                    (predicción
   post-hoc)            pendiente de DESI Y3 / Euclid                   confirmada)
```

SSEE está **del lado de la ciencia** de la línea de numerología (falsable + rígido +
predicciones forward), **corto de confirmado**. Lo que lo movería a la derecha: medir
k_fs. Lo que lo mataría: k_fs ≠ 0.754, o descartar la evolución w₀wₐ.

**Soft-spots abiertos, catalogados (no ocultos):**
1. **OP-9** — coeficiente 594.28 débilmente seleccionado; mitad UV congelada (→OP-10).
2. **Puente dimensional** Planck→Mpc para H₀ (roadmap #1, Type-P).
3. **A_s** no derivado (OP-18).
4. Otros en OPEN_PROBLEMS.md (23 rastreados).

---

## 5. Criterios de muerte (pre-registrados)

El modelo se declara FALSO si:
- **DESI Y3 / Euclid** miden k_fs fuera de [0.70, 0.81] h/Mpc → mata la partícula φ-DM.
- La evolución **w₀wₐ** se descarta (phantom crossing excluido a >3σ) → mata el sector DE.
- **S₈** two-sector se aleja >3σ de KiDS/DES con feedback bariónico controlado.
- **r** (tensor-a-escalar) se mide y excluye φ⁻¹⁰.

Si alguno ocurre, este documento se cierra con veredicto FALSO — que es lo que hace un
modelo científico y **no** un sistema numerológico.

---

*Última verificación: guardián VERDE 139; álgebra convencional y álgebra SSEE consistentes
(linaje KRYSTOS_V=φ+π+Ω, M_v=φ+π+K_v restaurado en código, prosa y PDFs).*
