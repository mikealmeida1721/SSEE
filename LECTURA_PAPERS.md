# Lectura de la suite — Papers 1→10, en orden

**Abierta 2026-07-27.** Antes se leyó el PRD (consolidado). Mike corrigió la
prioridad: **el PRD y el Sealed son resúmenes DERIVADOS; los Papers son el
árbol**. Si un Paper y un consolidado dicen cosas distintas, quien compare
concluye que el modelo se contradice a sí mismo — y eso pesa más que cualquier
acierto. Se riega el árbol, no la fruta.

## El método (no cambia)

Se lee **en orden de lectura**, se para en la **primera afirmación**, y no se
avanza hasta que las tres capas coincidan:

1. **la afirmación** — qué dice el texto
2. **la evidencia** — qué la respalda (script, log, figura, cita)
3. **la trazabilidad** — de dónde sale cada número y cada entidad

Reglas que Mike fijó para esta lectura:

- **Ningún número sin origen localizable.** Si no se puede señalar de dónde
  sale, no se afirma.
- **Las figuras cuentan como afirmación.** Una figura rancia miente igual que
  un número mal.
- **Cada MCMC debe apuntar a UNA pregunta declarada.** Es la mayor fuente de
  ambigüedad: unos dejan parámetros libres, otros los fijan, y confundirlos
  invalida la conclusión. Sólo en Paper 2 hay ~5.
- **Una pregunta abierta a la vez.** Lo lateral se anota aquí y NO se persigue
  hasta cerrar la actual.

## Estado

| Paper | Páginas | Estado |
|---|---|---|
| **1 — Framework** | 33 | 🔵 EN CURSO |
| 2 — MCMC | 26 | ⬜ |
| 3 — CMB | 24 | ⬜ |
| 4 — ToE | 16 | ⬜ |
| 5 — IS | 25 | ⬜ |
| 6 — φ-DM | 24 | ⬜ |
| 7 — EFT | 16 | ⬜ |
| 8 — Strong gravity | 20 | ⬜ |
| 9 — Hubble | 18 | ⬜ |
| 10 — UV | 14 | ⬜ |

## Bitácora — Paper 1

### Afirmación 1 · §1 pág. 1 — la lista de las siete constantes fundacionales
**Dice:** Ω≈4.7596, β≈2.3798, KAL₀≈5.5214, P≈6.3776, K_v≈9.5192, T_r≈11.9935,
M_v≈14.2788.
**Verificado contra `ssee_core`:** cuatro correctas. **TRES TRUNCADAS**, no
redondeadas — en la primera página del paper que define el diccionario:

| | exacto | decía | correcto |
|---|---|---|---|
| P (PYROS) | 6.377661 | 6.3776 | **6.3777** |
| K_v | 9.519253 | 9.5192 | **9.5193** |
| M_v | 14.278880 | 14.2788 | **14.2789** |

Corregidas en **14 sitios de 5 documentos** (P1, P4, P6, P9, Unified) — regar el
árbol, no sólo donde se vio. Ancladas en R30 (12 cantidades × 18 documentos).
→ **CERRADA**

### Afirmación 2 · §3 ec. (w₀,wₐ) — la caja destacada
**Dice:** `w₀ = −(3φ+π)/(2(φ+π)) = −0.8399`, `wₐ = −(2φ+π)/(2(φ+π)) = −0.6699`.
**Dos hallazgos:**
1. **`−0.6699` está MAL.** Exacto −0.669974886 → a 4 decimales **−0.6700**
   (el 5.º dígito es 7). Error 7.5e−5, por encima de la media unidad. Corregido
   en P1, P4, P5 y P6. (`w₀ = −0.8399` sí es correcto: su 5.º dígito es 4.)
2. **La entidad estaba colapsada a su reducción.** La forma `(2φ+π)/(2(φ+π))`
   es algebraicamente correcta —equivale exacto a P/I_g— pero aparecía SOLA.
   Los dos denominadores dan el mismo número, 9.5193, y son entidades
   DISTINTAS: M_v = φ+π+K_v vs I_g = π+P. Ningún chequeo numérico las separa;
   sólo la construcción. Restaurado el linaje delante de la reducción, con nota
   explícita de por qué el orden importa.
→ **CERRADA**

## Anotado al pasar (NO perseguir hasta cerrar la pregunta en curso)

- Propagar a Paper 3 los ΔBIC verificados hoy: plik_lite −24.02 y plik full
  −25.77 (hoy sólo viven en PRD/Sealed y en el Registro).
- Los Papers muestran menos decimales que los consolidados (0.4033 vs 0.403302).
  No es error —son redondeos correctos— pero R30 ya vigila las 9 anclas en los
  18 documentos.
