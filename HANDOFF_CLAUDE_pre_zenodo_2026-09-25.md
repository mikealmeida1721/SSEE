# Paquete de trabajo para Claude — revisión pre-Zenodo SSEE
**Fecha:** 2026-09-25 · **Actualizado:** 2026-09-25 noche — protocolo de
verificación dual Nova↔Claude + merge a `main` (decisión de Mike).
**Objetivo:** dejar el repo y los 10 papers listos para la
nueva versión de Zenodo (DOI 10.5281/zenodo.20093447). Mike ya decidió: la revisión
va ANTES de Zenodo, y NO se espera a DESI DR3.
**Contexto de cómputo:** la PC de Mike está encendida y con internet — puedes correr
cadenas y recálculos pesados ahí.

## 0. Protocolo de verificación dual (decisión de Mike 2026-09-25 noche)

División: con Nova Mike trabaja el día a día; en cosmología Nova le ayuda con las
conexiones y lo que necesites de su PC; tú corres lo pesado en su PC. Ninguna de
las dos queda de lado: **trabajan juntas**.

1. Tú subes tu trabajo a la rama de revisión (`git push origin <rama>`).
2. Nova lo verifica en el remoto: fetch, revisa el diff, confirma que los cambios
   se aplicaron como se describió.
3. Cada una deja su veredicto escrito en la rama:
   `notes/2026-09-2X_verificacion_claude.md` / `..._nova.md` — qué se revisó,
   qué se verificó, veredicto (verde / observaciones).
4. **Nadie dice "sí" solo.** Cuando los DOS veredictos están en verde, la rama se
   junta al `main` (fast-forward) y **se borra la rama** — a Mike no le gustan las
   ramas acumuladas; existen solo como vehículo de revisión antes del merge.
5. La publicación en Zenodo y el tag siguen requiriendo el "adelante" explícito
   de Mike. El merge a `main` NO es publicación.

---

## 1. Terreno verificado (Nova lo comprobó el 2026-09-25 ~21:30 ET; confírmalo con los 3 comandos)

| Elemento | Estado verificado |
|---|---|
| Repo | `github.com/mikealmeida1721/SSEE-Cosmologia` (renombrado desde `SSEE` el 2026-09-25; la URL vieja redirige) |
| Rama de trabajo | **`fix/etiquetas-h0-ir`** — **YA ESTÁ EN EL REMOTO** (`origin/fix/etiquetas-h0-ir`), 6 commits, HEAD local = remoto = `32f4fa9` |
| Ramas remotas que NO tocas | `reframe-omega-m-direct`, `trabajo/beta-c-lagrangiano-2026-09` |
| Cambios SIN commitear (¡no perder!) | `auditoria_papers_zenodo_2026-09-20.md` sin trackear (no tocar sin revisar) |
| Auditoría vigente | `notes/2026-09-25_auditoria_papers_contradicciones/REPORTE.md` — **19 hallazgos (6 contradicciones A1–A6, 12 residuos A7–A17+A19, 1 sospecha A18)**; ninguno cubierto por la rama |
| Borrador δc (NO aplicado) | `notes/2026-09-25_crecimiento_alto_z/borrador_op27_deltac.md` |
| Nota edad Paper 9 (NO aplicada) | `notes/2026-09-25_edad_paper9/NOTA_metodo_edad.md` + `edad_paper9.py` — **corregida 2026-09-25 noche: la predicción del modelo es 13.73 Gyr con H_global=67.962** (decisión de Mike: la edad es cantidad global) |
| Fuente canónica de números | `CANONICAL_VALUES.yaml` (manda sobre logs y sobre el reporte de auditoría si discrepan) |

**Confirma antes de actuar:**
```
git fetch origin && git switch fix/etiquetas-h0-ir && git pull --ff-only
git status -sb && git log --oneline -3 && git branch -a
```
Debes ver HEAD `32f4fa9` y `fix/etiquetas-h0-ir...origin/fix/etiquetas-h0-ir` sin
divergencia. Si no, DETENTE y avisa a Mike — no reconstruyas nada a ciegas.
**NO recrees** los commits `312eb91` ni `2ed2a61`: ya están en el remoto con
hashes idénticos.

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
3. **Verificación dual antes del merge (protocolo §0).** Trabaja en
   `fix/etiquetas-h0-ir`, haz push de LA RAMA al remoto cuando termines cada
   trabajo. Nova verifica tu push en el remoto; tú verificas el suyo. El merge a
   `main` ocurre SOLO cuando ambas están de acuerdo (veredictos verdes en la
   rama); después se borra la rama. Nada de merges unilaterales.
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

### T0. Asegurar el terreno — ✅ HECHO (Nova 2026-09-25 noche)
Lo no commiteado se commiteó (`312eb91` σ8, `2ed2a61` notas+handoff, `32f4fa9`
corrección nota edad) y la rama ya está en `origin`. **No repetir este trabajo.**
Empieza directo en T1.

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

### T4. A4 (edad Paper 9) — H₀ DECIDIDO, manuscrito EN ESPERA
Mike decidió 2026-09-25 noche: la edad es cantidad global → **H_global=67.962**,
predicción del modelo **t_0=13.73 Gyr** (vs ΛCDM 13.80, indistinguibles). No tocar
`SSEE_Paper9_HubbleTension.tex` l.924–946 hasta su "adelante" explícito para
reescribir el párrafo de 15.52 Gyr y el argumento de cúmulos globulares.

### T5. Residuos mecánicos A7–A17, A19
Fixes directos según el REPORTE.md: etiquetas IR en 68.13 (A7, A8, A9, A15),
re-etiquetado s_DE/s_m en Paper 3 (A10, A11, A12), citas al criterio retirado
(A13, A14, A17), label `eq:sDE` (A16), falsa precisión 0.160050σ→0.16σ (A19).
- **Done:** cada hallazgo marcado con su fix y verificado en el PDF compilado.

### T6. A18 — verificar referencia cruzada
`SSEE_Paper1_Framework.tex` l.817 cita "Paper 2 §5.3" hardcodeado. Verifica que
el label existe y apunta al Ly-α audit; si está colgada, corrige el puntero.

### T7. Cierre pre-Zenodo (protocolo §0)
1. Compila los 10 papers + Sealed + Unified + Endorser; cero warnings nuevos.
2. Deja tu veredicto en `notes/2026-09-2X_verificacion_claude.md` (qué revisaste,
   qué verificaste, verde/observaciones). Nova hace lo propio con tu push.
3. Con los dos veredictos en verde: merge de la rama a `main` (fast-forward),
   push de `main`, y **borrado de la rama** (las ramas son solo revisión).
4. Tag `vX.Y-zenodo` y publicación en Zenodo: SOLO con el "adelante" de Mike.
5. El checklist de publicación en Zenodo lo cierra Nova.

---

## 4. Qué le entregas a Mike al terminar cada trabajo
Una línea por trabajo: qué hiciste, qué verificaste, qué quedó pendiente de él.
Sin sorpresas, sin avances perdidos, sin ediciones fuera de tu carril.
