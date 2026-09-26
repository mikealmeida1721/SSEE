# Paquete de trabajo para Claude — revisión pre-Zenodo SSEE
**Fecha:** 2026-09-25 · **Objetivo:** dejar el repo y los 10 papers listos para la
nueva versión de Zenodo (DOI 10.5281/zenodo.20093447). Mike ya decidió: la revisión
va ANTES de Zenodo, y NO se espera a DESI DR3.
**Contexto de cómputo:** la PC de Mike está encendida y con internet — puedes correr
cadenas y recálculos pesados ahí.

---

## 1. Terreno verificado (Nova lo comprobó el 2026-09-25 ~20:35 ET; confírmalo con los 3 comandos)

| Elemento | Estado verificado |
|---|---|
| Repo | `github.com/mikealmeida1721/SSEE` (rama `main` en remoto = `de122d4`) |
| Rama de trabajo | **`fix/etiquetas-h0-ir`** — existe SOLO en local, **no está en el remoto**. HEAD = `25c0b89` "fix: s_m/s_DE relabel + Δ≡SSEE−ΛCDM convention" |
| Ramas remotas que NO tocas | `reframe-omega-m-direct`, `trabajo/beta-c-lagrangiano-2026-09` |
| Cambios SIN commitear (¡no perder!) | `src/p02_mcmc/ssee_press_schechter.py` (corrección σ8 **0.7446 → 0.8153** del 2026-09-25, con comentario de justificación) + `results/figures/fig_press_schechter.{pdf,png}` regeneradas + `auditoria_papers_zenodo_2026-09-20.md` sin trackear |
| Auditoría vigente | `notes/2026-09-25_auditoria_papers_contradicciones/REPORTE.md` — **19 hallazgos (6 contradicciones A1–A6, 12 residuos A7–A17+A19, 1 sospecha A18)**; ninguno cubierto por la rama |
| Borrador δc (NO aplicado) | `notes/2026-09-25_crecimiento_alto_z/borrador_op27_deltac.md` |
| Nota edad Paper 9 (NO aplicada) | `notes/2026-09-25_edad_paper9/NOTA_metodo_edad.md` + `edad_paper9.py` |
| Fuente canónica de números | `CANONICAL_VALUES.yaml` (manda sobre logs y sobre el reporte de auditoría si discrepan) |

**Confirma antes de actuar:**
```
git status -sb && git log --oneline -3 && git branch -a
```
Si tu clon no muestra `fix/etiquetas-h0-ir` con HEAD `25c0b89`, o los cambios sin
commitear no están, DETENTE y avisa a Mike — no reconstruyas nada a ciegas.

---

## 2. REGLAS DURAS — léelas dos veces

1. **Lo no commiteado se commitea primero.** El cambio σ8 (0.7446→0.8153) y las
   figuras son avance real con justificación escrita. `git add` + commit en
   `fix/etiquetas-h0-ir` ANTES de cualquier otra edición. No lo reescribas, no lo
   "mejores", no lo pierdas.
2. **δc (A1–A3) y edad de Paper 9 (A4): PREPARAR sí, APLICAR no.** El postulado
   δc=1.6284 está falsificado bajo los supuestos corridos, pero su retiro del
   manuscrito y la reescritura de Paper 9 requieren el "adelante" explícito de
   Mike. Deja los parches listos en notas, no en `manuscript/`.
3. **Nada a `main` sin aprobación de Mike.** Trabaja en `fix/etiquetas-h0-ir`,
   haz push de LA RAMA al remoto como respaldo cuando esté verde. El merge a
   `main` + tag pre-Zenodo lo autoriza Mike.
4. **El manuscrito manda sobre los logs.** Los logs (`results/logs/`) contienen
   cifras intermedias que el texto ya corrigió. Antes de citar un número de un
   log, verifícalo en `CANONICAL_VALUES.yaml` o en el .tex.
5. **Convención Δ unificada:** Δ≡SSEE−ΛCDM en toda la suite; **negativo favorece
   SSEE**; se declara UNA vez por documento. Canónicos: `deltaBIC_lcdm_minus_ssee:
   6.43` (>0 favorece SSEE) ≡ ΔBIC(SSEE−ΛCDM)=−6.43; `deltaDIC_ssee_minus_lcdm:
   −5.66`. Si el reporte de auditoría y el YAML discrepan en redacción, **manda
   el YAML** y escala la discrepancia.
6. **No toques** `polux-automatizaciones`, `polux.online`, ni las ramas remotas
   ajenas. Nova lleva el checklist de Zenodo y la auditoría cruzada; tú no
   dupliques su trabajo: tu carril es revisión consolidada, MCMC/figuras,
   `model_comp`, sincronización y propagaciones **autorizadas**.

---

## 3. Trabajos, en orden

### T0. Asegurar el terreno (primero, 15 min)
1. Commit de lo no commiteado (regla 1).
2. `git push -u origin fix/etiquetas-h0-ir` (respaldo de la rama; NO es merge a main).
3. `git fetch --all` y confirma que tu local ve lo mismo que el remoto.
- **Done:** `git status` limpio, la rama existe en `origin`, nada sin respaldo.

### T1. Verificar el diff de Nova + tabla `model_comp` única (núcleo del encargo)
1. Revisa el diff `main..fix/etiquetas-h0-ir`: ~45 ediciones en 7 .tex
   (re-etiquetado s_m/s_DE + convención Δ). Marca cada hunk como verificado o
   repórtalo.
2. **Reconstruye UNA tabla `model_comp` para Paper 2** con columnas: fuente (log),
   N_data, prior, parámetros, convención. Datos:
   - Corrida canónica: `results/logs/mcmc_paper2_reframe.log` (100w×25k, seed 42,
     N_eff=78170, accept=0.715) → BIC_SSEE=16.49, BIC_ΛCDM=22.92, **ΔBIC=−6.43**
     (convención Δ≡SSEE−ΛCDM).
   - Los números **17.24 / 5.68 / 4.91** son de OTRA corrida
     (`mcmc_paper2_3models_om308.log`, Ω_m congelado) — van a una nota al pie o
     fuera de la tabla, nunca mezclados.
   - DIC: 14.95 vs 20.61 → ΔDIC=−5.66 (`dic_from_chains.log`).
3. Unifica los comentarios de convención en `CANONICAL_VALUES.yaml` (regla 5).
- **Done:** tabla reconciliada aritméticamente (22.92−16.49=6.43 ✓), sin mezcla
  de corridas, convención declarada una vez.

### T2. A5+A6 — signo ΔBIC en la suite
Unifica a la convención de la regla 5 en: `SSEE_Paper1_Framework.tex` (l.812 y
menciones), `SSEE_Paper2_MCMC.tex` (l.967 + tabla), `SSEE_Endorser_Summary.tex`
(l.62), `SSEE_Sealed_Journal.tex` (l.644), `SSEE_Unified_Journal.tex` (l.330).
- **Done:** el mismo estadístico, el mismo signo, la misma convención declarada.

### T3. A1–A3 (δc) — preparar, NO aplicar
Con `notes/2026-09-25_crecimiento_alto_z/borrador_op27_deltac.md` como base,
prepara el parche exacto para `SSEE_Paper4_ToE.tex` (§deltac l.478–498, bullet
l.772) y `SSEE_Paper5_IS.tex` (§JWST l.1402–1440, tabla de ratios 1.051–1.887 →
≈1.00×). Guárdalo como nota, **no edites `manuscript/`**.
- **Done:** parche listo para el "adelante" de Mike, cero ediciones al manuscrito.

### T4. A4 (edad Paper 9) — EN ESPERA
No tocar `SSEE_Paper9_HubbleTension.tex` l.924–946 hasta que Mike decida el
encuadre y qué H₀ representa la edad global. La investigación está en
`notes/2026-09-25_edad_paper9/`.

### T5. Residuos mecánicos A7–A17, A19
Fixes directos según el REPORTE.md: etiquetas IR en 68.13 (A7, A8, A9, A15),
re-etiquetado s_DE/s_m en Paper 3 (A10, A11, A12), citas al criterio retirado
(A13, A14, A17), label `eq:sDE` (A16), falsa precisión 0.160050σ→0.16σ (A19).
- **Done:** cada hallazgo marcado con su fix y verificado en el PDF compilado.

### T6. A18 — verificar referencia cruzada
`SSEE_Paper1_Framework.tex` l.817 cita "Paper 2 §5.3" hardcodeado. Verifica que
el label existe y apunta al Ly-α audit; si está colgada, corrige el puntero.

### T7. Cierre pre-Zenodo (solo con aprobación de Mike)
1. Compila los 10 papers + Sealed + Unified + Endorser; cero warnings nuevos.
2. Presenta a Mike: diff final `main..fix/etiquetas-h0-ir`, lista A1–A19 con
   estado, y la decisión pendiente (T3/T4).
3. Solo con su "adelante": merge a `main`, tag `vX.Y-zenodo`, push.
4. El checklist de publicación en Zenodo lo cierra Nova.

---

## 4. Qué le entregas a Mike al terminar cada trabajo
Una línea por trabajo: qué hiciste, qué verificaste, qué quedó pendiente de él.
Sin sorpresas, sin avances perdidos, sin ediciones fuera de tu carril.
