# Bitácora de Archivado — SSEE

Registro de **qué está archivado, cuándo y por qué**. Nada en `archive/` es
vigente: es el historial del proceso. Lo vivo está en los cajones vigentes
(ver Mapa de Vigencia). **Regla:** cuando algo deja de ser canónico, se mueve
aquí **con una entrada en esta bitácora** (fecha + razón + qué lo reemplaza).
Nada obsoleto se queda en un cajón vivo sin marcar; nada se archiva sin entrada.

**Organización (reorg 2026-06-24):** `archive/` se ordena **por tipo**, no por fecha.
Cada cajón contiene solo su tipo; las fechas y razones viven en esta bitácora.

```
archive/
  pdfs/      — PDFs de documentos viejos (drafts, versiones previas, enmiendas)
  figuras/   — figuras no usadas en los papers finales (.pdf + .png)
  chains/    — cadenas MCMC legacy (datos pesados, pre-reframe)
  codigo/    — scripts y fuentes (.py/.tex/.bib) superados
  superado/  — docs de trabajo superados, plegados al Registro canónico
  zenodo_dictionary.zip — empaquetado histórico del diccionario
```

---

## Mapa de Vigencia (qué cajón debe estar siempre actual)

| Cajón | Estado | Auditado por |
|---|---|---|
| `manuscript/` (12 fuentes .tex) | **VIGENTE** | guardián + compilación + lectura hostil |
| `docs/` (PDFs compilados) | **VIGENTE** | recompilan de `manuscript/` |
| `src/` (módulos `pXX_*`, `verificacion/`, `estadistica/`) | **VIGENTE** | guardián (procedencia) |
| `results/logs/` (logs de procedencia) | **VIGENTE** | guardián R9 |
| `VERIFICATION_LEDGER.md`, `CANONICAL_VALUES.yaml` | **VIGENTE (fuente de verdad)** | guardián |
| `CLAUDE.md`, `SSEE-Vault/` | **VIGENTE (memorias)** | guardián Capa Memorias |
| `zenodo_dictionary/` | **VIGENTE** (registro publicado) | scripts corren limpio |
| `AUDIT.md`, `RIGOR_CHECKLIST.md` | **VIGENTE** (sistema de auditoría) | — |
| `archive/` | **HISTÓRICO** — no se edita, no se cita como vigente | esta bitácora |
| `sandbox_unificado/`, `notes/`, `class_ssee/output/`, `eftcamb_ssee/` | **TRABAJO** — no canónico, no auditado | — |

> `AUDIT_V7_PREFLIGHT.md` y `HALG_PIFI_CHANGEMAP.md` (antes en la raíz) se archivaron
> en `superado/` el 2026-06-24 (snapshots ya ejecutados; llevan banner interno).

---

## Catálogo de Archivado

### `pdfs/` — versiones anteriores de papers
**Archivados:** 2026-04 a 2026-05 · **Razón:** drafts / versiones intermedias superadas.
**Reemplazados por:** `docs/SSEE_PaperN_*.pdf` (recompilan de `manuscript/`).
- `SSEE_Paper1_Framework_v3.6.pdf`, `SSEE_Paper2_MCMC_Validation*.pdf` (draft/v1/v2),
  `SSEE_Paper2_Summary_ES_v1.pdf`, `SSEE_Paper2_draft.pdf`,
  `SSEE_Paper3_CMB_Confrontation*.pdf` (v1/v2), `SSEE_Paper3_Ch1_Perturbations_v1.pdf`,
  `SSEE_Paper3_Conclusion_PlanckPR4_v1.pdf`, `SSEE_Paper3_draft.pdf`,
  `SSEE_Paper3_CMB_backup.pdf` (backup pre-edición).
- `SSEE_Paper6_Holographic_Amendment.pdf` (2026-05-19) · enmienda holográfica explorada,
  no adoptada en el Paper 6 canónico (antes en `experimental/`).

### `figuras/` — figuras no usadas en los papers finales
**Archivadas:** 2026-06-04 · **Razón:** generadas en exploración; no entran en los
PDFs finales. Conservadas como registro visual del proceso (corner plots, scans,
diagnósticos B1/B3/IS, hi_class, Lyman-α audit, MCMC fase 4, rutas P(k)). 48 archivos (.pdf + .png).

### `chains/` — cadenas B1 CMB legacy pre-reframe
**Archivadas:** 2026-05-08 · **Razón:** cadenas Cobaya MCMC SSEE (`ssee_cmb.*`) de la tabla
B1 de Paper 3 con ingredientes **pre-reframe** (Ω_m,CMB=0.31993 vía MIRA, ω_b=0.02237,
H₀≈67.04). Superadas por el reframe ω_m-directo (ω_b=0.02242, ω_c=0.11951 forward fijo,
Σm_ν=0.0690, ancla H_alg=67.962); re-corrida con `src/p03_cmb/ssee_paper3_b1_mcmc.py`.
Las cadenas ΛCDM (`lcdm_cmb.*`) NO se archivan: sus ingredientes no cambiaron.

### `codigo/` — scripts y fuentes superados
**Archivados:** 2026-04 a 2026-06 · **Razón:** reemplazados por los módulos `pXX_*` vigentes
de `src/`, o experimentos exploratorios de un solo uso.
- **Reemplazados por módulos vigentes:** `ssee_inflation_connection.py` → `src/pB_inflation/`;
  `ssee_uv_completion.py` → `src/p10_uv/`; `ssee_paper3_cobaya.py` →
  `src/p03_cmb/ssee_paper3_cobaya_unified.py`; `ssee_paper2_mcmc_legacy.py` (2026-04-30) →
  `src/p02_mcmc/`.
- **Fuente obsoleta:** `SSEE_EFT_Fundamental_obsolete.tex` (2026-05-07) → Paper 7;
  `SSEE_appendix_Friedmann.tex` (2026-04-21) → integrado en Paper 1/2;
  `ssee_paper3_docs_duplicate.bib`, `ssee_paper6_docs_duplicate.bib` (duplicados de `manuscript/`).
- **Pre-canónico m_φ:** `p6_honest_matrix_PRECANONICAL_2026-06-02.py` (2026-05-25) → Paper 6
  canónico (2026-06-04).
- **Scratch exploratorio:** sensibilidad de cúmulos, scans α_M/θ*/σ₈, IS growth, candidatos;
  ningún resultado canónico depende de ellos.

### `superado/` — docs de trabajo plegados al Registro canónico
**Archivados:** 2026-06-08 (lote inicial) y 2026-06-24 (snapshots de raíz) · **Razón:** su
contenido se consolidó en `VERIFICATION_LEDGER.md`, `CANONICAL_VALUES.yaml` y los papers canónicos.
- `H0_CASCADE_AUDIT.md` → cascada H₀ en P9/P10 + Registro §B.
- `SSEE_CONSTANTS_AUDIT.md` → `CANONICAL_VALUES.yaml` + guardián.
- `P6_CLEANUP_NOTES.md` → Paper 6 canónico.
- `OPEN_PROBLEMS_LAGRANGIAN_MAP.md` → `OPEN_PROBLEMS.md`.
- `SEALED_STATUS.md` → estado del Sealed (vigente en CLAUDE.md).
- `ssee_paper6_mass_derivation.py` (2026-06-04) → derivación m_φ previa; superada por la
  forward-prediction canónica.
- `AUDIT_V7_PREFLIGHT.md` (archivado 2026-06-24) → snapshot de auditoría pre-vuelo del
  2026-06-14; lleva banner interno "pre-reframe ω_m-directo".
- `HALG_PIFI_CHANGEMAP.md` (archivado 2026-06-24) → changemap del H_alg/π-φ ya ejecutado.
