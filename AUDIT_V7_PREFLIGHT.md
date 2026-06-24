# Auditoría Pre-Vuelo V7 — 12 documentos (10 papers + Unified + Sealed)

**Fecha:** 2026-06-14 · **Estándar:** revista / premio (referee hostil)
**Regla:** marcar lo NUEVO (no re-litigar lo ya señalado, p.ej. cronología).
**GATE:** nada corregido todavía · **Zenodo NO se toca hasta luz verde de Mike.**

> ⚠️ **SNAPSHOT histórico pre-reframe ω_m-directo (2026-06-19).** Este documento
> es una foto de la auditoría al 2026-06-14 y lista como "canónicos" varios valores
> que **fueron superados** por el reframe ω_m-directo (OP-8 disuelto): H_local IR/UV
> 71.87/72.05 → **72.86/73.040**; ΔBIC CMB −32.2 → **−33.0** (full plik); m_φ 36.95 →
> **41.02 eV** (SOLAR²·KRYSTOS); k_fs 0.659 → **0.762**; S₈/σ₈ two-sector 0.766/0.742
> → **0.759/0.748**. No se reescribe (es registro de lo que se auditó entonces). Para
> los valores vigentes manda `CANONICAL_VALUES.yaml` + `VERIFICATION_LEDGER.md`.

> Protocolo: (1) marcar aquí → (2) Mike revisa → (3) corregir lo aprobado →
> (4) solo entonces actualizar Zenodo V7 (nombre, descripción, conteo, resultados).

---

## Capa objetiva — completada

| Check | Resultado |
|---|---|
| Guardián `ssee_verify.py` | ✅ VERDE 102/102 (constantes φ,π + 3 memorias) |
| Compilación de los 12 | ✅ 0 errores, 0 citas/refs sin resolver (el "undefined" de P1/P2 era warning de fuente) |
| Deriva de valores obsoletos en papers | ✅ limpio (m_φ 5.60, k_fs 0.493, H₀ 66.55/67.068/67.756, σ8 0.702/0.737, cascada 71.90/72.077 — ninguno presente como valor vigente) |
| Conteo de papers stale ("seven-paper") | ✅ ninguno |
| Headlines canónicos (71.87 / 72.05 / 0.766 / −32.2 / m_φ 36.95 / k_fs 0.659) | ✅ consistentes entre docs |

---

## HALLAZGOS

> **F1 y F2 — RESUELTOS 2026-06-14** (con luz verde de Mike). Detalle abajo.

### F1 — Error en el Registro canónico (NO en los papers) · Severidad MEDIA · ✅ RESUELTO
**Qué:** `VERIFICATION_LEDGER.md` líneas 119–120 dan, como canónico re-anclado
2026-06-09 (Σm_ν=0.069):
- ΔBIC CMB diagonal = **−24.7**
- χ²_r CMB TT (SSEE) = **1.045**

**Verdad objetiva (re-corrida CAMB hoy, misma config Σm_ν=0.069, H₀=67.04):**
- ΔBIC diagonal = **−28.0**
- χ²_r TT = **1.044** (TE 1.041, EE 1.040, PP 0.866, combinado 1.041)
- Log: `/tmp/paper3_rerun.log` · script `src/p03_cmb/ssee_paper3_cmb.py`

**Veredicto:** los **papers (P1, P3, Sealed) están CORRECTOS** con −28.0 / 1.044.
El **Registro está mal**: −24.7 y 1.045 no aparecen en ningún log/script y **no se
reproducen** bajo su propia procedencia declarada. CLAUDE.md (−28.0/1.044) correcto.
CANONICAL_VALUES.yaml no contiene estos valores (no arrastra el error).

**Acción propuesta (tras luz verde):** corregir `VERIFICATION_LEDGER.md` L119–120
→ −28.0 / 1.044, con nota de que la re-corrida 2026-06-09 NO cambió estos valores
respecto a Σm_ν=0.06 (la entrada −24.7/1.045 fue mala anotación). **No tocar papers.**

### F2 — Grieta en el sistema de verificación · Severidad MEDIA (proceso)
**Qué:** el guardián solo recomputa las constantes algebraicas (sección A) y las 3
memorias; **no verifica los valores de pipeline de la sección B** del Registro. Por
eso F1 sobrevivió sin detección. Es el mismo tipo de punto ciego que dejó pasar la
cronología.
**Acción propuesta:** añadir al guardián un check que, para cada valor de sección B,
exija (a) un log de procedencia existente y (b) coincidencia con ese log; o al menos
marque "sin log → no verificable". Convierte el Registro en algo auditável, no
declarativo.

---

## Pendiente — capa subjetiva (referee hostil por paper)

La capa objetiva (números, compilación, consistencia) está hecha. Falta la lectura
hostil documento-por-documento (física, derivaciones, coherencia argumental,
overclaim residual, presentación nivel-revista). Plan al retomar:

- [ ] P1 Framework · P2 MCMC · P3 CMB · P4 ToE · P5 IS
- [ ] P6 φDM · P7 EFT · P8 StrongGravity · P9 Hubble · P10 UV
- [ ] Unified Journal · Sealed Journal (las dos caras consolidadas → Zenodo)
- [ ] Pasada conjunta: coherencia inter-paper (mismo número, misma definición, mismas unidades en los 12)

---

## Bitácora — lectura hostil por paper

### Paper 1 (Framework) — en curso
**Hallazgos corregidos (2026-06-14):**
- **TÍTULO** "A Zero-Parameter Framework" → **"A Minimal-Parameter Framework"** (overclaim
  vs reframe aprobado y vs su propio abstract). +cita del título en P2 actualizada. +comentario L1.
- Abstract "extensiones en Papers~3--7" → **"3--10"** (stale; el cuerpo ya decía 3--10).
- Cronología residual L1452 "deposited 2026-01-28, prior to all Planck comparisons" →
  reframe honesto + DOI diccionario limpio (la encontró la nueva Capa R1 automática).
**Verificados como NO-error:** ΔBIC −5.55 (consistente P1/P2, test distinto del +7.91 MCMC);
claim dark:baryon (loose en abstract pero el cuerpo L1175 es preciso: KAL₀·Ω_b h²·n_s=0.41σ).
**Cuerpo §1.1–1.5 (corregidos 2026-06-14):**
- **P1-A (HIGH):** register (l.132) marcaba H₀=3(φ+π)²=67.96 como "Retrodiction" mientras §1.4
  lo llama "Type-P numerical coincidence" → contradicción interna + overclaim. Status → "Type-P coincidence".
- **P1-B (HIGH):** §1.3 enumera "tres postulados físicos" y §1.4 "cuatro (D,S,M,I)" sin reconciliar
  → confusión de conteo. Añadido puente conservador (solo afirma P1→S, que el paper ya hace).
- **P1-B' (presentación, 2026-06-14):** la **colisión de vocabulario** que causaba la pregunta "¿3 o 4?":
  §1.3 llamaba a los tres principios "physical postulates" y §1.4 al presupuesto "four postulates" — misma
  palabra, dos ejes distintos. Reservada "postulate" SOLO para el presupuesto D/S/M/I; §1.3 → "three physical
  **principles**" (+refs cruzadas "Postulate~1"→"Principle~1" en L202, L252). Cero cambio de física ni de
  conteo honesto (sigue 2 fundamentales D,S + 2 auxiliares M,I); compila 30 pp 0 err; guardián VERDE 111.
- **P1-C (menor):** quitada nota editorial "(Paper 5 corrected)" del register (l.136).
- **Notado, NO tocado (defensible):** título de §1.4 "Scope of the Zero-Parameter Claim" + el conteo
  honesto (1 dim H₀ + 0 adim vs ΛCDM 1+5) — el contenido es honesto y fuerte; el "zero-parameter"
  ahí está scoped y definido. Opcional harmonizar el título.
**Cuerpo §1.5–2 (corregidos):**
- **P1-D (MEDIUM, R7):** nombre mitológico "YGG\_PROTON" (l.438, 450) → "geometric form"; codename
  "MIKAEL\_V" quitado de la tabla de notación (l.1396). P1 limpio de YGG/MIKAEL.
- **Flag cross-paper (c_s²):** P1/P7 dicen k-essence c_s²=1; Paper 5 dice c²_s,eff=0 (IS). Verificar
  coherencia en la pasada conjunta (¿dos regímenes o tensión?).

### Paper 2 (MCMC) — en curso (lectura hostil 2026-06-14)
- **P2-A (HIGH) — ✅ CORREGIDO 2026-06-14:** L666 afirmaba que H₀^alg=67.962 y H₀^MIRA=67.037
  "agree at **0.10σ** … striking … supports interpreting MIRA as a physical bridge". **FALSO:** la
  diferencia es 0.925 km/s/Mpc = **1.71σ** (Planck-error 0.54), no 0.10σ (verificado). El "0.10σ" es stale
  (correspondería a H_MIRA~67.91). Doble problema: (1) número falso; (2) overclaim — usa el "acuerdo
  asombroso" como argumento pro-MIRA, cuando H₀^alg es **coincidencia Type-P** (igual que el register de P1),
  no predicción. Las otras σ de la sección verifican (0.74σ H_MCMC↔H_MIRA ✓; ~1.85σ H_MCMC↔H_alg ✓). **Fix
  propuesto:** reportar el gap real (~1.3%, 1.7σ), marcar H₀^alg como Type-P, y quitar la frase "striking…
  physical bridge" (el soporte de MIRA vive en el ajuste CMB + two-sector, no en esta falsa coincidencia).
  NO aplicado aún — cambia un argumento.
- **P2-B (menor, notación):** L338 escribe $\Omega_b h^2=(\pi-\varphi)/H_0^{\rm SSEE}=0.02242$, mezclando con
  un H₀ dimensional; P1 §5 usa la forma limpia $(\pi-\varphi)/[3(\varphi+\pi)^2]$. Armonizar a la forma algebraica
  (mismo número, evita el aire de "dividir adimensional por km/s/Mpc"). Alimenta la confusión Type-P/estructural.
- **P2-C (menor):** L678 "~1.2 km/s/Mpc difference (2.7%)" — los dos no cuadran (1.2 km/s/Mpc ≈1.8%; 2.7%≈1.8 km/s/Mpc).
  Revisar a qué par se refiere y corregir el % o el valor.
- **P2-D (stale, MEDIA) — ✅ CORREGIDO 2026-06-14:** L935 decía "Paper 3 reports χ²_r=**1.062** (TT)"; Paper 3
  reporta **1.044** en todas partes (L67/L597) + log canónico F1 + el propio P2 L962. → 1.044. (Mismo tipo que F1.)
- **P2-E (consistencia con P2-A) — ✅ CORREGIDO 2026-06-14:** L984 "SSEE **predicts** H₀^alg… consistent with Planck
  1.1σ" contradecía el reframe Type-P de L666. → "algebraic value… Type-P coincidence, not a derived prediction…
  lies within 1.1σ"; +L991 "prediction"→"coincidence". P2 ahora coherente en las 3 menciones de H_alg.
- **P2-F (inconsistencia interna, MEDIA-ALTA) — ⏸ ESPERA OK MIKE (verificado contra la cadena):** la mediana del
  posterior CPL se reporta DISTINTA: Tabla §6 (w₀=−0.800, wₐ=−0.583) vs Apéndice B (w₀=−0.731, wₐ=−0.919), Δwₐ~0.8σ.
  **Resuelto con la cadena pública** `results/logs/mcmc_chains_professional_CPL_ckpt.npz` (63.75M, cols [H0,Om,w0,wa,Obh2]):
  mediana real **w₀=−0.731⁺⁰·¹⁰⁷₋₀·¹⁰³, wₐ=−0.923⁺⁰·⁴⁴³₋₀·⁴⁷⁴** → **el Apéndice B es correcto, la Tabla §6 está STALE**
  en w₀,wₐ (H0/Om/Obh2 OK dentro de rounding: 67.30/0.316/0.02237). Además L585-587 atribuye el 0.05σ "a ese
  posterior CPL" — FALSO: el 0.05σ (χ²_2D=0.080) es vs el **best-fit de DESI DR2**, no vs el posterior CPL; vs la
  mediana CPL real SSEE está a ~1.2σ (1.0σ w₀, 0.55σ wₐ), aún consistente. **Fix propuesto:** (1) Tabla 2 columna CPL
  w₀→−0.731, wₐ→−0.923 (+H0→67.30, Om→0.316); (2) reescribir L585-587 separando "0.05σ vs DESI best-fit" de
  "~1.2σ vs posterior CPL". El titular 0.05σ (vs DESI) SE MANTIENE válido.
- **P2-F': sub-corrección pendiente + HALLAZGO degeneración (lab Mike 2026-06-14):** mi propio reframe de P2-F
  puso "~1.2σ jointly" — **MAL** (cuadratura por-eje ignorando la correlación ρ=−0.95). El propio (Mahalanobis 2D)
  es **1.74σ** al prior Planck. El laboratorio (`notebooks/lab_w0wa_degeneracy.ipynb`) reveló además algo real:
  la distancia SSEE↔CPL **depende del prior de H₀** porque SSEE vive SOBRE la degeneración w₀-w_a-H₀:
  reponderando, CPL@H_MIRA=2.15σ, @Planck=1.74σ, @H_alg=0.99σ — sube H₀ → CPL se desliza hacia SSEE.
  El H₀ que clava SSEE ≈ **68.8** (NO SH0ES 73; a 73 sobrepasa a fantasma extrema). PERO reponderar muere
  (ESS≈0) más allá de ~68 → **rerun real lanzado** (`src/p02_mcmc/rerun_cpl_h0anchors.py`, 5 anclas incl. 73).
  **RESUELTO 2026-06-15:** (a) convención reconciliada — bajo la estándar (σ-equivalente 2dof, DESI/Planck)
  SSEE vs CPL BAO = **1.2σ** (mi "1.74σ" era Mahalanobis crudo); el "1.2σ" original del paper ERA correcto, solo
  se limpió la redacción ("comfortably within 1σ"→"1.2σ joint 2D significance"). (b) **rerun REAL de anclas
  (`rerun_cpl_h0anchors.py`) corrigió el espejismo de la reponderación:** NINGÚN H₀ aterriza el ajuste libre en SSEE
  — w₀ sigue a H₀ pero wₐ queda clavado en ~−0.92 a todo H₀ (SSEE quiere −0.67); SH0ES(72.76)→(−0.94,−0.945), MÁS
  lejos. El "H0*=68.8 clava SSEE" era artefacto ESS≈0. SSEE NO es SH0ES (confirmado por muestreo directo).
  (c) **Subsección robustez AÑADIDA a §9** (honesta, sin claim H₀-clava-SSEE): Fase 4 full excluye ΛCDM 2.4σ,
  SSEE 0.68σ, ρ=−0.92 + **figura `fig_w0wa_degeneracy.pdf`**. Mitiga R4 (likelihood completa vs prior comprimido).
  El futuro de SSEE se juega en wₐ (DESI DR3/Euclid): predicción falsable limpia. Compila 28pp 0err, guardián VERDE 111.
- **Verificado OK (no-error):** Ω_b h²=0.02218 es el prior BBN Cooke2018 (no stale, distinto del algebraico 0.02242,
  declarado a 0.4σ); 21.3σ→0.63σ two-sector ✓; r_d 175.6×0.840=147.6 ✓; H_alg=M_v·Ω=3(φ+π)²=67.962 ✓;
  χ²_2D=0.080 ✓; CMB-only w₀ 6.3σ tratado honestamente vía r_d.

### Paper 3 (CMB) — lectura hostil 2026-06-14/15
- **Verificado sólido:** χ²_r TT=1.044 (=F1 canónico), ΔBIC diagonal −28.0, Cobaya −32.2 (Δχ²=2.875, ln6469=8.775),
  k=4 conserv. −14.7, n_s=1−φ⁻⁷=0.96556, MIRA=(φ+β)/2=(3φ+π)/4, peaks 220/536/813, r_d≈147, H₀local 71.87+Type-P
  72.86 "for comparison only" (coherente con reframe P2-A). §k=2 maneja overclaim honestamente. Falsabilidad pre-comprometida.
- **P3-B (símbolo, MENOR) — ✅ CORREGIDO:** L839 γ_eff=0.55+0.05(1+**w₀**)=0.542 — el valor requiere w(z=1)=−1.175
  (CPL en pivote), no w₀=−0.84 (que daría 0.558). Símbolo corregido a w(z=1) + Linder formula. Número OK.
- **P3-D (stale S₈, MEDIA-ALTA) — ✅ CORREGIDO:** L933 "two-sector reduces S₈ from 3.9σ to **2.6σ**" contradecía el
  canónico. Triple-verificado (Paper 6 L112 + CANONICAL_VALUES + CLAUDE.md = **0.766, 0.01σ, RESUELVE**). El barrido
  automático no lo cazó ("2.6" genérico). Corregido a "resolves, S₈_eff=0.766, 0.01σ KiDS". (Stale era MÁS conservador
  → fix mejora S₈ pero es alinear-al-canónico, no inventar.)
- **P3-A (presentación, MENOR — notado):** el ΔBIC "canónico" −32.2 usa k=2 (point-estimate, τ/A_s fijos); el B1 MCMC
  riguroso usa k=3 (−23.6). Ambos reportados y decisivos (<−10). Llamar "canónico" al k=2 (más favorable) es defensible
  (τ/A_s NO se re-fitean en ese escenario) pero un referee podría preguntar. Opcional liderar con k=3.
- **P3-C (presentación, MEDIA — notado):** la tabla fσ8 usa γ_IS=0.657/σ8=0.792 (pre-Paper 5, etiquetado) y el canónico
  (0.74σ) vive solo en el caveat box. Honesto pero confuso (¿cuál es la predicción?). Considerar actualizar la tabla al
  canónico γ=0.554/σ8=0.820 o marcarla claramente como "upper-bound cross-check".
- **P3-E (contradicción cross-paper, ALTA) — ✅ RESUELTO 2026-06-15:** §hiclass decía αK^SSEE=0.073 (bare-field) y que
  0.4032 "**does NOT represent the SSEE prediction**", contradiciendo a Paper 7 (autoridad EFT) que mantiene AMBOS
  válidos y declara 0.4033 como el observacional. **Forensia git:** el texto nació CORRECTO el 05-07 (commit 1fa4ee4,
  junto con el script hi_class: 0.4033 = SSEE) → el commit 05-08 `fc30b6e` ("Integrate B1 full MCMC + update EFT")
  lo sobre-corrigió al bare-field y degradó 0.4033 → la pasada de reconciliación 05-15 `d64da64` SOLO tocó Paper 7,
  dejando a P3 fosilizado en el estado intermedio. Papers 8/10 ya citan 0.4033 universalmente. **Fix aplicado:**
  L961-967 reescrito — conserva ambos valores (0.073 bare / 0.4033 effective), quita "does not represent", declara
  0.4033 como el observacionalmente relevante (idéntico framing a P7). Compila 25pp/0err, guardián VERDE 111.

### Hallazgos para auditorías de OTROS papers (marcados, no tocados)
- **P4-A (SUSTANTIVO) — ✅ RESUELTO 2026-06-15:** "Proton/Neutron register" (L278-279, 385) era numerología
  del sistema génesis SSEE_UNIFICADO (Ygg_P=0.319 "66% gravedad + 34% fotones / pensamiento cristalizado")
  filtrada al paper serio. DOS defectos: (1) circular — "observado 0.319 (YGG)" es valor INTERNO, no medición;
  (2) error físico — 0.3201=Ω_m,CMB es materia TOTAL (DM-dominada), NO bariónica (Ω_b=0.049 ya está en tabla).
  Inflaba el overclaim "7 ratios / 3 dominios / astronómicamente pequeña". **Confirmado por Mike: viene de su
  trabajo génesis fenomenológico (8 dominios), la cosmología es la única rama que se graduó a física falsable.**
  **Fix:** quitadas filas protón/neutrón de tabla cross-domain; §matter reescrito (0.3201 = geométrico vs 0.31993
  = MIRA two-sector canónico, las dos rutas a materia total, coinciden 0.054%); intro L66 limpiada; overclaim
  probabilístico → look-elsewhere honesto (Paper 1). 0 restos génesis. Compila 17pp/0err, guardián VERDE 111.
  ⚠️ Pendiente desacople pre-submission: SSEE_UNIFICADO público + H₀=73.5 (67.4+5.52 aditivo) CONTRADICE
  cascada 71.87 de los papers — desvincular de Zenodo V7/arXiv (ver memoria project_genesis_repo_risk).
  → 0.319 preservado como OP-16 (test falsable lattice-QCD, severidad baja, post-auditoría).
- **P4-B (overclaim) — ✅ RESUELTO 2026-06-15:** título ya de-escalado a "Two-Axiom Cosmology" pero cuerpo
  usaba "SSEE-ToE" (Theory of Everything) 5× en captions/L682 — fosilizado. Cambiado a "\SSEE". Mismo patrón
  título-corregido/cuerpo-fosilizado que P3-E. Compila 17pp/0err.
- **P4-C (desacople génesis) — ✅ RESUELTO 2026-06-15:** subtítulo L30 decía "SSEE-V3.6 / Genesis 5.12" —
  enlace directo al sistema multidominio de riesgo. Quitado "/ Genesis 5.12" (el origen se preserva en el repo
  SSEE_UNIFICADO + git + memoria de linaje, NO en el paper de arXiv). Mike: "ten en cuenta que ese es el origen"
  → respetado: separar ≠ borrar.
- **P4 nota (no tocado):** filename SSEE_Paper4_ToE.tex aún tiene "ToE" pero NO es visible al referee (el PDF
  compilado tiene título limpio); renombrar cascadearía refs/docs → diferido. "Nine Sovereignties"/sovereignty
  = renombrado-diccionario (todo #10), no error.
- **Paper 7:** codename "MIKAEL\_V" en tabla de notación (l.630) → quitar en su turno.
- **Notación cross-paper:** `Ω_DNAV` (subíndice del scaffold Ω) en 7 papers — jerga establecida pero
  redundante (Ω basta). Limpieza coordinada futura, no urgente.

### Pendiente de automatizar (cuando esté limpio)
- **R7 mitológicos** (YGG, LUCIFER, OSIRIS, MIKAEL, MIKE…): automatizar como check duro DESPUÉS de
  limpiar P4/P7 (si se activa antes, bloquea por esas instancias reales).

**Cuerpo §3–4 (lectura hostil 2026-06-14) — NUEVOS hallazgos, NO tocados (gate):**
- **P1-E (MEDIUM, error objetivo de física) — ✅ CORREGIDO 2026-06-14:** §3.2 daba $\rho_\nu\approx23.2$ eV/cm³.
  Era **3× alto** (usaba $n_{\rm total}=336$ cm⁻³ × Σm_ν en vez de $n_{\rm especie}=112$ cm⁻³ × Σm_ν). Corregido a
  **7.7 eV/cm³** con la fórmula explícita $\rho_\nu=n_\nu^{\rm species}\sum m_\nu$. **NO toca f_ν ni la tabla:**
  f_ν≈0.018–0.022 es un cociente de espacio de fases (fracción bajo v_esc), independiente de la normalización ρ_ν;
  χ²_r=0.122 intacto. Compila 30 pp 0 err; guardián VERDE 111.
- **P1-F (contradicción real) — ✅ CORREGIDO 2026-06-14 con el dato fuente:** se consultó el abstract de
  Zhang2026 (arXiv:2602.06082, WebFetch): el 88% es vs la **masa dinámica MOND** (hidrostática bajo Milgrom),
  NO la newtoniana. Verificado en los 4 cúmulos: M_MOND≈M_bar/0.88, boost newtoniano/MOND≈4.8×, KAL₀ implícito
  5.33–5.45 ≈ 5.5214. Conclusión: 88% (marco MOND) y ×5.52 (mapeo SSEE bariones→masa newtoniana total) son **marcos
  distintos y NO se componen** (la frase vieja implicaba 88%×5.52=486%, absurdo). §3.1 reescrito para separar ambas
  masas explícitamente y citar Zhang2026 en el punto correcto. Compila 30 pp 0 err; guardián VERDE 111.
- **P1-G (MEDIUM, consistencia interna) — ✅ CORREGIDO 2026-06-14:** añadida nota tras Eq.(mdyn) — KAL₀ es la asíntota
  deep-MOND de KAL(x), adoptada como amplificador efectivo de masa encerrada en R₂₀₀; proyección KAL(x)-pesada diferida
  a trabajo futuro. Caveat honesto, sin cambiar ningún número (caveat de modelado, no de valor).
- **Verificado honesto/bien acotado (NO-error):** acción bg Eq.(action) sólo background, perturbaciones diferidas a
  EFT App.A; NEC violación reconocida y diferida a P4; phantom-crossing z*≈0.31 marcado como artefacto CPL lineal;
  H₀=3(φ+π)² L660 "numerical coincidence" (consistente con el fix del register). §3.4 lensing "qualitative mechanism…
  deferred to future work" — hedge correcto, sin overclaim.

**Flag cross-paper (pasada conjunta):** L482 k-essence $c_s^2=1$ (P1/P7) vs $c^2_{s,\rm eff}=0$ (P5 IS) — verificar
que son dos objetos distintos (escalar DE propagante vs sound speed efectivo del fluido viscoso), no contradicción.

**§3.4 lensing — hallazgo de rastreo (2026-06-14):**
- **OP-15 (Medium-High, NUEVO OP creado):** el offset del Cúmulo Bala vía κ(θ) desde KAL(x) estaba como
  hedge suave "deferred to future work" en §3.4, **sin OP numerado** → podía leerse como terminado. Formalizado:
  entrada completa en OPEN_PROBLEMS.md + fila en tabla índice + hedge de §3.4 reescrito visible apuntando a OP-15
  (distingue: masa-magnitud §3 SÓLIDA vs distribución-espacial Bala ABIERTA vs amplitud P8). Causa: hedge en prosa
  sin marca formal. Candidato a trabajar tras la auditoría (Mike, 2026-06-14).
- **Stale formula en OPEN_PROBLEMS.md L880 (corregido):** la "conclusión" de OP-14 mostraba
  $m_\varphi=\Sigma m_\nu\cdot H_0^{\rm alg}$ (retirada, 5.602 eV) como "identidad de Paper 6" presente, sin marca.
  → actualizada a la canónica (×535.28=36.95 eV) + nota de retiro; el argumento (mismo DoF) se conserva.
  **Gap de cobertura detectado:** OPEN_PROBLEMS.md no está entre las 3 memorias que escanea memory_sync, y la forma
  era simbólica (sin el número "5.602" en esa línea) → ni el scan de patrones retirados la habría cazado. Anotado
  para el perfeccionamiento del sistema (cobertura OPEN_PROBLEMS + patrones simbólicos, no solo numéricos).

**§A.5 EFT + §5 Domain/Falsification (lectura hostil 2026-06-14):**
- **Verificado algebraicamente (sólido):** n=(T_r−M_v)/2T_r=−0.09527 con 1/(2n−1)=w₀=−0.840; α=φ⁴/3=(3φ+2)/3=2.285;
  η/ε=3−8φ³=−5−16φ=−30.9; τ_Π H₀=KAL₀/(3Ω_DE)=2.191; N=2φ⁷=58.07; T_r=3(φ+β)=3(3φ+π)/2 y M_v=φ+π+K_v=3(φ+π)
  coinciden con §4. Ω_b h²=(π−φ)/[3(φ+π)²]=0.02242 (0.32σ Planck), clasificado **retrodicción** (dimensionless → forma
  estructural post-hoc legítima, distinto del Type-P de H₀ que sí tiene problema de unidades). Coherente con OP-1.
- **Cross-paper c_s² (FLAG para pasada conjunta, importante):** ahora son TRES contextos de c_s²=1
  — (a) k-essence escalar gradient-stable (§2.5, P7/P10), (b) señal bulk-viscosa en el límite de causalidad IS
  (Eq IS_causality §A.5, fija τ_Π) — más (c) Paper 5 c²_{s,eff}=0 (clustering). Plausiblemente son objetos
  distintos (velocidad del modo escalar vs señal viscosa vs sound speed efectivo de clustering), NO contradicción,
  pero un referee lo va a preguntar. Verificar definiciones y enunciarlas claras en la pasada conjunta P5↔P7↔P1.
  NO resuelto a ciegas (requiere leer P5/P7 a fondo).
- **S₈ multiplicidad (presentación, menor):** §A.5 sola cita S₈ = 0.818 / 0.837 / 0.847 / 0.766 / 0.59 (y σ₈ 0.742/0.820)
  — todos etiquetados por método/sector y consistentes con la divulgación two-Ω_m, pero visualmente "cherry-picky".
  Opcional: una tabla/frase consolidadora. NO es error (cada uno es un cálculo distinto, el titular 0.766 manda).
- **R7 menor:** "Nine Sovereignties closure" (L1099, tabla de falsación) — jerga inventada en el cuerpo; idealmente
  "the nine-constant closure" o referir al apéndice de linaje. No urgente.

**Paper 1 — estado:** §1–5 + §A.5 EFT auditados. Hallazgos objetivos corregidos (título, cronología, register,
ρ_ν, §3.1 MOND, KAL₀ caveat, OP-15). Pendiente menor: barrido de valores stale en apéndices B (Master Derivation
Chain) / glosario / prior space — tablas, scan rápido. Luego P2.

### Mejoras al sistema de auditoría (esta sesión)
El guardián pasó de **102 → 109** comprobaciones. Capas nuevas, una corrida las cubre todas:
- **R9 Procedencia** (valores de pipeline vs log committeado; atrapó/selló F1).
- **R1 Cronología** automatizada (atrapó violación residual en P1).
- **R2 Multidominio** automatizada.
- **R10 Overclaim** "zero-parameter" + **R11 Conteo de serie** stale (atraparon 5 instancias
  en P1/P2/P3/P8/P10; refinadas para 0 falsos positivos).

## Para la actualización Zenodo V7 (cuando haya luz verde)

Cambia respecto a la versión previa (a confirmar con Mike):
- Nombre de presentación del depósito
- Descripción (10 papers + 2 journals; framing minimal-parameter; postdicción honesta)
- Conteo de documentos y resultados headline (w₀wₐ, S₈=0.766, H₀ cascada 71.87/72.05, m_φ=36.95)
- DOI del diccionario limpio (concept 10.5281/zenodo.20684908) como procedencia del registro de constantes

---

## Paper 5 (IS) — auditado 2026-06-16

**Bandera conjunta c²_s (P1↔P5↔P7↔P10↔Unified) — RESUELTA (redacción, 0 física):**
- Diagnóstico: el mismo campo tenía 3 c²_s citados sin marcar régimen. No es error de física —
  P10 ya lo distinguía bien (c²_s,ad≈0.967 modo campo vs c²_s≈0 modo IS, mismo KAL₀).
- P5 L182: "canonical scalar field k-essence" → "canonical-kinetic (constant-w) k-essence **proxy**".
- P5 §instability: +párrafo — campo K(X)=X/KAL₀ estable de por sí (P7/P10); IS = descripción de
  fluido; el mismo KAL₀ fija ζ̃=KAL₀/3 → cero parámetros nuevos. IS NO es mecanismo aparte.
- P1 L482: "c²_s=1; Papers 7 y 10" (se citaba en contra) → "c²_s>0; =1 IR, ~0.60–0.95 con UV".
- Unified L679/696: "c²_s=1 exact" → "=1 IR; [0.60,0.95] completo; c²_s>0".
- Compila: P5 26pp · P1 30pp · Unified 23pp · 0 err. Guardián VERDE 112.

**Sección de crecimiento (γ_IS, σ₈, S₈, fσ8) — LIMPIA, valores canónicos verificados:**
- γ_IS=0.5503±0.0003 ✓ · G=D₁/D₁=1.0109 ✓ · σ₈=0.811×1.0109=0.820 ✓
- S₈=0.820×√(0.3199/0.3)=0.8466 ✓ · tensiones 1.0σ Planck / 3.85σ KiDS / 3.90σ DES (aritmética confirmada)
- fσ8 media=0.74σ vs ΛCDM 0.73σ ✓ (refleja corrección 2026-06-13; el 2.67σ/0.50σ viejos ausentes)

**Framing S₈ (anti-sesgo, modo referee) — EJEMPLAR:**
- Single-sector S₈=0.847 presentado HONESTAMENTE como "the honest baseline tension / upper bound /
  single-sector challenge", 3.9σ KiDS, "marginally worse than ΛCDM". NO overclaim.
- Resolución two-sector atribuida correctamente a P6: S₈=0.766 (0.01σ KiDS), σ₈=0.742,
  baryonic 0.758 (0.06σ DES); m_φ=36.95 eV "forward prediction, no parameter fitted to lensing".
- Criterio de falsabilidad explícito (S₈>0.86 favorece SSEE; S₈<0.80 desfavorece ambos).

**Veredicto P5:** limpio. Único hallazgo = la bandera c²_s (redacción), ya corregida. Sin números tocados.
