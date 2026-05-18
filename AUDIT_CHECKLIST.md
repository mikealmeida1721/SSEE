# SSEE — Checklist de Auditoría Escalonada

**Propósito:** procedimiento escrito y reproducible que se corre **después de
consolidar un paper o hacer un cambio**. No depende de la atención de quien audita
ese día — se siguen los pasos, en orden, y se marca cada casilla.

Motivación: las auditorías previas solo comprobaban `script ↔ PDF` (consistencia
interna de un bucle cerrado). Un bucle cerrado **bendice un error si el dato de
entrada ya era malo**. Así sobrevivieron el bug de datos fσ₈ y el bug E²(z=0)=1.16
de Paper 6. Esta checklist añade los invariantes que el modelo no puede falsear
consigo mismo.

---

## Regla de oro — fuente única de verdad

- Toda constante algebraica, todo dato observacional y la función de fondo de
  Friedmann viven en **`src/ssee_core.py`** y en **`data/raw/*.csv`**.
- Ningún script redefine `phi`, `w0`, los datos `fsigma8`, etc. localmente.
  Si un script lo hace, **es un hallazgo de auditoría** aunque el número sea correcto:
  una copia local correcta hoy es una copia divergente mañana.
- Cada archivo de datos en `data/raw/` lleva cabecera `#` con la línea `Source:`
  y la referencia primaria de cada dato (valor **y** error).

---

## Cómo elegir el tier

| El cambio fue… | Tier |
|---|---|
| Redacción, tipografía, bibliografía, comentarios | **Tier 0** |
| Un resultado numérico cambió en un paper (tabla, abstract, figura) | **Tier 1** |
| Cambió física, una ecuación, un script de cálculo, o datos de entrada | **Tier 2** |

Cada tier **incluye** los anteriores. Ante la duda, sube de tier.

---

## Tier 0 — Editorial / bibliográfico

- [ ] El `.tex` compila con **0 errores**.
- [ ] **0** `Citation undefined` y **0** `Reference undefined` en el `.log`
      (papers con `\bibliography{}` requieren `bibtex`; ver [feedback-audit-method]).
- [ ] Ningún bibitem huérfano. Para P1, el grep de citas **debe incluir**
      `SSEE_EFT_section.tex` (entra por `\input`).
- [ ] Sin tocar números ni ecuaciones → no se sube de tier.

## Tier 1 — Cambio numérico (Tier 0 +)

- [ ] **Procedencia externa.** Cada dato observacional usado traza a su referencia
      citada y coincide con ella en **valor y error**. No basta que script y PDF
      concuerden — el dato mismo debe ser el de la fuente.
- [ ] **Cross-documento.** La cantidad cambiada tiene el **mismo valor** en TODOS
      los sitios donde aparece: el `.tex` del paper, `CLAUDE.md`, `AUDIT.md`,
      `README.md`, `CHANGELOG.md`, `SSEE_Unified_Journal.tex`,
      `SSEE_Endorser_Summary.tex`, y las memorias.
- [ ] **Script ↔ PDF.** El número del PDF es la salida actual del script que lo
      produce (re-correr el script, no confiar en la salida cacheada).
- [ ] El número y su σ/tensión son **autoconsistentes** (p.ej. una tensión Nσ
      coincide con `|valor − obs| / σ`).

## Tier 2 — Cambio de física / script (Tier 1 +)

- [ ] **`python3 src/ssee_core.py`** corre y todos los asserts de cordura pasan.
- [ ] El script importa sus constantes y datos de `ssee_core` — **no** los redefine.
- [ ] **Cordura física.** El script tiene `assert` explícitos para sus invariantes:
      `E²(a=1) ≈ 1`, `0 ≤ Ω_m,source ≤ 1`, `0 < G < 1.5`, `σ₈` en rango razonable,
      densidades suman 1. Que reviente ruidosamente, no que imprima `G=6.15`.
- [ ] **Cross-implementación.** Si otro script calcula la misma cantidad física,
      ambos concuerdan dentro de tolerancia. Lo ideal: una sola función importada.
      (Caso de origen: tres ODEs de crecimiento divergentes en Paper 6.)
- [ ] El script se **re-corrió de cero** y su salida se copió al paper.
- [ ] Si el cambio invalida resultados de otros papers, esos papers entran a
      auditoría también (el cambio se propaga).

---

## Cierre de toda auditoría

- [ ] Hallazgos verificados **uno por uno contra los archivos reales** antes de
      corregir — nunca a ciegas (los auditores externos suelen traer PDFs viejos;
      ver `project_ssee_state`).
- [ ] Correcciones aplicadas en **todos** los lugares afectados, no solo donde se
      detectó.
- [ ] Si la auditoría es por emoción (un resultado mejoró), aplicar el protocolo
      anti-sesgo: actuar como referee hostil, comprobar que no empeoró otro
      observable (ver memoria `feedback_emotional_bias`).
- [ ] Commit con mensaje que nombre el tier y los archivos tocados.

---

## Pendientes conocidos para este sistema

- `src/ssee_audit_consistency.py` es el auditor viejo basado en regex con valores
  canónicos hardcodeados (algunos ya obsoletos, p.ej. busca `m_φ=5.71`). Debe
  reescribirse para importar de `ssee_core` y comparar, o retirarse.
- `AUDIT.md` arrastra números fσ₈ desactualizados (Paper 6); su pasada Tier 1
  depende de cerrar antes el bug E²/dos-sectores de Paper 6.
- Los scripts existentes aún no importan de `ssee_core`; migrarlos es trabajo de
  la fase de reinicio, uno por uno, cada uno con su auditoría Tier 2.
