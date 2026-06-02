# Código viejo — cuarentena (NO borrar hasta comprobar)

Scripts movidos aquí porque quedaron **superados por una versión más nueva y
canónica**. No se borran: pueden contener algo valioso bajo un nombre engañoso.
Para revertir cualquiera: `git mv src/codigo_viejo/<archivo> src/`.

| Archivo | Fecha | Reemplazado por | Motivo |
|---|---|---|---|
| `ssee_paper3_cobaya.py` | 2026-05-18 | `ssee_paper3_cobaya_unified.py` (2026-05-24) | El viejo evaluaba la verosimilitud Planck en un punto fijo; el nuevo escanea H₀ y es el que CLAUDE.md cita para H_MIRA=67.068. |
| `ssee_uv_completion.py` | 2026-05-15 | `ssee_paper10_verification.py` (2026-05-18) | "First calculation" de K(X); la verificación paso a paso del Paper 10 lo subsume (M⁴=5φ⁸ρ_crit). |
| `ssee_inflation_connection.py` | 2026-04-27 | `ssee_op2_spectral_index.py` + `ssee_paperB_Nstar.py` | Pregunta exploratoria temprana (¿n_s deriva de inflación?); resuelta por α-attractor N*=2φ⁷. |

## Revisados a fondo (2026-06-02) — VEREDICTO: NO mover, quedan en src/

Cada candidato se leyó aplicando dos criterios: (a) ¿algún paper lo necesita
(figura/resultado citado)? (b) ¿contiene algo importante sin verificar? Quien
cumple cualquiera de los dos NO va a cuarentena. Resultado: ninguno de los 15
califica — a diferencia de los 3 de arriba (reemplazo canónico claro + fecha
vieja), todos estos son de mayo 24–25 y atados a hilos vivos.

**Provenance de figura publicada (criterio a):**
- `regenerate_fig_corner_p2.py` → genera `fig_corner_ssee_professional.pdf`, que Paper 2 incluye (L551). Es el script que produjo el corner committeado con la cadena prior-MIRA. KEEP.

**Caza de mecanismo OP-8/OP-9 — diagnósticos vivos sin cerrar (criterio b):**
- `ssee_disformal_hubble.py`, `ssee_mira_disformal_test.py`, `ssee_void_hubble.py`, `ssee_vainshtein_scale.py`, `ssee_paper6_kinetic_braiding.py` — candidatos a mecanismo de MIRA/H₀/screening/αB. NO superados por nada. Si se confirman descartados, van a `src/mira_attempts/`, no aquí.

**Review activo de Paper 6 (criterio b — atados a P6_CLEANUP_NOTES.md):**
- `p6_canonical_table.py`, `p6_complete_matrix.py`, `p6_honest_matrix.py` — matrices honestas fit-vs-predicción (todas 2026-05-25, untracked, scratch vivo). `honest_matrix` documenta el fit de alpha_WDM. KEEP hasta cerrar el review de P6.
- `ssee_paper6_lyman_alpha_audit.py` — defensa de m_φ=5.60 eV vs Lyman-α (OP-9 vivo).
- `ssee_nonlinear_s8.py` — estimación Halofit de S₈ (territorio OP-5 Nivel 2).

**Calibración H₀ reciente (2026-05-24, provenance de números canónicos):**
- `ssee_h0_prior_experiment.py`, `ssee_h0_mira_only.py` — produjeron `fig_h0_three_priors.pdf`/`fig_h0_four_priors.pdf` (untracked, recién generados) y el hallazgo de comparación de priors (qué prior prefiere DESI). Sustentan la narrativa H₀=66.55.
- `ssee_paper2_mcmc_lcdm_baseline.py` — baseline ΛCDM para ΔBIC limpio vs SSEE-MIRA.

**Tooling de verificación:**
- `ssee_precision_audit.py` — auditoría de propagación de error de redondeo a H₀/r_d/Ω_m. Pregunta distinta a la del guardián `ssee_verify.py` (matching de constantes). KEEP.

## Exonerados — alimentan figuras de manuscritos publicados (provenance, KEEP en src/)
Verificado por cruce `\includegraphics` ↔ `savefig`:
- `ssee_eftcamb_validation.py` → `ssee_eftcamb_CMB_TT.pdf` + `ssee_eftcamb_Pk.pdf` (pese a la nota "EFTCAMB never installed" en CLAUDE.md).
- `ssee_paper4_toe.py` → `fig_toe_cmb_TT.pdf` + `fig_toe_derivations.pdf` (Paper 4).
- `ssee_paper8_figures.py` → `fig_paper8_vainshtein.pdf` + `fig_paper8_lensing_ratio.pdf` (Paper 8).
- `ssee_paper9_figures.py` → `fig_paper9_fscreen_z.pdf` + `fig_paper9_h0_tension.pdf` (Paper 9).
- `ssee_paper10_figures.py` → `fig_paper10_alphaK_vs_alpha.pdf` + `fig_paper10_KX_profile.pdf` (Paper 10).
- `ssee_eft_verification.py` → `fig_eft_verification.pdf` + `fig_eft_beta_convergence.pdf` (Paper 7).
