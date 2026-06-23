# Bitácora de Archivado — SSEE

Registro de **qué está archivado, cuándo y por qué**. Nada en `archive/` es
vigente: es el historial del proceso. Lo vivo está en los cajones vigentes
(ver Mapa de Vigencia). **Regla:** cuando algo deja de ser canónico, se mueve
aquí **con una entrada en esta bitácora** (fecha + razón + qué lo reemplaza).
Nada obsoleto se queda en un cajón vivo sin marcar; nada se archiva sin entrada.

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
| `AUDIT.md`, `RIGOR_CHECKLIST.md`, `AUDIT_V7_PREFLIGHT.md` | **VIGENTE** (sistema de auditoría) | — |
| `archive/` | **HISTÓRICO** — no se edita, no se cita como vigente | esta bitácora |
| `sandbox_unificado/`, `notes/`, `class_ssee/output/`, `eftcamb_ssee/` | **TRABAJO** — no canónico, no auditado | — |

---

## Catálogo de Archivado

### PDFs de papers — versiones anteriores
**Archivados:** 2026-04 a 2026-05 · **Razón:** versiones intermedias/drafts superados.
**Reemplazados por:** `docs/SSEE_PaperN_*.pdf` (recompilan de `manuscript/`).
- `SSEE_Paper1_Framework_v3.6.pdf`, `SSEE_Paper2_MCMC_Validation*.pdf` (draft/v1/v2),
  `SSEE_Paper2_Summary_ES_v1.pdf`, `SSEE_Paper2_draft.pdf`,
  `SSEE_Paper3_CMB_Confrontation*.pdf` (v1/v2), `SSEE_Paper3_Ch1_Perturbations_v1.pdf`,
  `SSEE_Paper3_Conclusion_PlanckPR4_v1.pdf`, `SSEE_Paper3_draft.pdf`,
  `SSEE_Paper3_CMB_backup.pdf` (backup pre-edición).

### `superseded_2026-06-04/` — docs de trabajo plegados al Registro canónico
**Archivados:** 2026-06-08 · **Razón:** su contenido se consolidó en
`VERIFICATION_LEDGER.md` y los papers canónicos durante la reorganización canónica.
- `H0_CASCADE_AUDIT.md` → cascada H₀ en P9/P10 + Registro §B.
- `SSEE_CONSTANTS_AUDIT.md` → `CANONICAL_VALUES.yaml` + guardián.
- `P6_CLEANUP_NOTES.md` → Paper 6 canónico (m_φ=36.95, 2026-06-04).
- `OPEN_PROBLEMS_LAGRANGIAN_MAP.md` → `OPEN_PROBLEMS.md`.
- `SEALED_STATUS.md` → estado del Sealed (vigente en CLAUDE.md).

### `superseded_2026-06-04/ssee_paper6_mass_derivation.py` — m_φ pre-canónico
**Archivado:** 2026-06-04 · **Razón:** derivación previa de m_φ; superada por la
forward-prediction canónica (36.95 eV, cero fiteo). Ver `p6_honest_matrix_PRECANONICAL_2026-06-02.py`.

### `unused_figures_2026-06-04/` — figuras no usadas en los papers finales
**Archivadas:** 2026-06-04 · **Razón:** generadas en exploración; no entran en los
PDFs finales. Conservadas como registro visual del proceso (corner plots, scans,
diagnósticos B1/B3/IS, hi_class, Lyman-α audit, MCMC fase 4, rutas P(k)).

### `codigo_viejo/` — scripts superados por módulos actuales de `src/`
**Archivados:** 2026-04 a 2026-05 · **Razón:** reemplazados por los módulos `pXX_*` vigentes.
- `ssee_inflation_connection.py` → `src/pB_inflation/`.
- `ssee_uv_completion.py` → `src/p10_uv/`.
- `ssee_paper3_cobaya.py` → `src/p03_cmb/ssee_paper3_cobaya_unified.py`.

### Scripts sueltos obsoletos
- `ssee_paper2_mcmc_legacy.py` (2026-04-30) → `src/p02_mcmc/` (versión MIRA actual).
- `SSEE_EFT_Fundamental_obsolete.tex` (2026-05-07) → Paper 7 (`SSEE_Paper7_EFT.tex`).
- `p6_honest_matrix_PRECANONICAL_2026-06-02.py` (2026-05-25) → Paper 6 canónico (2026-06-04).
- `SSEE_appendix_Friedmann.tex` (2026-04-21) → integrado en Paper 1/2.

### `scratch/` — experimentos exploratorios de un solo uso
**Razón:** pruebas puntuales (sensibilidad de cúmulos, scans α_M, θ*, σ₈, IS growth,
candidatos). Conservadas como registro del proceso; ningún resultado canónico depende de ellas.

### `experimental/` — direcciones exploradas y no adoptadas
- `SSEE_Paper6_Holographic_Amendment.pdf` (2026-05-19) · enmienda holográfica explorada,
  no adoptada en el Paper 6 canónico.
