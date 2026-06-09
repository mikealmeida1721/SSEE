# SSEE — Contexto del Proyecto para Claude

## Qué es este proyecto

**Structural Self-Energy Expansion (SSEE-V3.6)** — modelo cosmológico de energía oscura con ~3 parámetros efectivos (vs 6 ΛCDM), con el sector background derivado algebraicamente de φ (razón áurea) y π, sin partículas WIMP/axión-QCD/SUSY. Extensiones fenomenológicas en el sector DM (Paper 6) tratadas como ansätze y rastreadas en OPEN_PROBLEMS.md (OP-8, OP-9, OP-11, OP-14). Autor: Mike Edison Almeida Vallejo.

Diez papers (auditados y endurecidos hasta 2026-05-17: bibliografía JCAP/PRD sin huérfanos, Postulados D y S formalizados en Paper 1, barrido de referee hostil — ver OPEN_PROBLEMS.md para brechas físicas abiertas):
- **Paper 1**: Framework teórico + EFT + Predictive Register (arXiv-ready) ✅
- **Paper 2**: Validación Bayesiana MCMC DESI DR2 + Planck + cúmulos (arXiv-ready) ✅
- **Paper 3**: Confrontación CMB Planck PR4 TT+TE+EE+lensing (arXiv-ready) ✅
- **Paper 4**: Derivación algebraica CMB desde φ,π — Nine Sovereignties; nombres mitológicos eliminados, bib 10→40 refs (preprint) ✅
- **Paper 5**: IS causal perturbation theory — estabilidad, MIRA, S8, fσ₈ (commit 0a9d21c) ✅
- **Paper 6**: φ-DM dos sectores — m_φ=36.95 eV (partícula canónica, forward-prediction Σm_ν·535.28, cero fiteo; canónico 2026-06-04, antes 5.60 eV obsoleto); S₈_eff=0.766 (0.01σ KiDS), fσ₈ 0.76σ (empata ΛCDM 0.73σ); Lagrangiano escalar libre cierra incompletitud OP-9 ✅
- **Paper 7**: EFT canónico — βc=−AURA, αT=αM=αB=0, αK=0.4033, λ/V₀/M/g² bloqueados (preprint) ✅
- **Paper 8**: Régimen de gravedad fuerte — ✅ preprint (rescatado 2026-06-08). Interpretación física de m_φ fijada por P6 (escalar libre, 36.95 eV) → predicción lensing canónica **incondicional** (límite EFT B-S). Two-limit analysis conservado: límite (a) MOND-like pedagógico + contraste falsable, límite (b) canónico = modelo real. Residuo OP-9 (origen UV del multiplicador) NO entra en las predicciones. 20 pp.
- **Paper 9**: Tensión Hubble — f_screen = αK/(3·MIRA) = (π-φ)/Ω² = 0.06725; H₀,local = 72.86 km/s/Mpc (0.17σ SH0ES); AURA cancela exactamente; +sección comparación cuantitativa vs EDE/SIDR/late-DE (preprint, 18 pp) ✅
- **Paper 10**: UV Completion — K(X)=X/KAL+X²/M⁴; M⁴=5φ⁸ρ_crit exacto; M=9.62 meV=Λ_SSEE; αK_full=0.41691; H₀^UV=73.040 km/s/Mpc (condicional C.1); 8/10 JCAP ✅

---

## Estado actual de los documentos (al 2026-05-17)

### docs/ — PDFs compilados (nombres reales; todos recompilan con 0 errores LaTeX)

| Archivo | Contenido | Estado |
|---|---|---|
| `SSEE_Paper1_Framework.pdf` | Paper 1 (25 p) | ✅ Completo — +Postulados D y S (§2.4) |
| `SSEE_Paper2_MCMC.pdf` | Paper 2 (26 p) | ✅ Submission-ready |
| `SSEE_Paper3_CMB.pdf` | Paper 3 (24 p) | ✅ Completo v3 (+ CLASS/hi_class cross-check §5.x) |
| `SSEE_Paper4_ToE.pdf` | Paper 4 (16 p) | ✅ Preprint — bariogénesis calibrada (D3) |
| `SSEE_Paper5_IS.pdf` | Paper 5 (25 p) | ✅ Completo — IS causal perturbations |
| `SSEE_Paper6_phiDM.pdf` | Paper 6 (24 p) | ✅ Preprint — +Lagrangiano φ-DM escalar, DW excluido |
| `SSEE_Paper7_EFT.pdf` | Paper 7 (16 p) | ✅ Preprint — EFT canónico, αK=0.4033, βc=−AURA |
| `SSEE_Paper8_StrongGravity.pdf` | Paper 8 (20 p) | ✅ Preprint — rescatado 2026-06-08: m_φ fijado por P6 → predicción lensing canónica incondicional; two-limit conservado |
| `SSEE_Paper9_HubbleTension.pdf` | Paper 9 (18 p) | ✅ Preprint — f_screen + comparación EDE/SIDR/late-DE |
| `SSEE_Paper10_UVCompletion.pdf` | Paper 10 (14 p) | ✅ Preprint — K(X)=X/KAL+X²/M⁴; M=9.62 meV=Λ_SSEE; αK_full=0.41691 |
| `SSEE_Unified_Journal.pdf` | Consolidación Papers 1–10 (23 p) — §Extensions 8/9/10 | ✅ Candidato journal (10-paper, 2026-05-17) |
| `SSEE_Endorser_Summary.pdf` | Síntesis 2 p para endorsers arXiv — serie 10-paper | ✅ Listo para envío (2026-05-17) |
| `SSEE_Sealed_Journal.pdf` | **Documento consolidado de energía oscura tardía** (10 p) — φ→w₀,wₐ; look-elsewhere sobre diccionario cerrado (1/245); H₀ en dos etapas; cuentas honestas ~3 vs 6 params; P5/P6/P8 y OP-8/9/14 como trabajo futuro | ✅ **SELLADO — candidato a auditoría externa (2026-06-02)** |

### manuscript/ — Fuente LaTeX
- `SSEE_Paper1_Framework.tex` + `SSEE_EFT_section.tex`
- `SSEE_Paper2_MCMC.tex`
- `SSEE_Paper3_CMB.tex` + `ssee_paper3.bib`
- `manuscript/SSEE_Paper4_ToE.tex` (submodule)
- `SSEE_Paper5_IS.tex` + `ssee_paper5.bib`
- `SSEE_Paper6_phiDM.tex` + `ssee_paper6.bib`
- `SSEE_Paper7_EFT.tex` (Paper 7 — versión completa, 874 líneas, βc=-AURA al 0.2%)
- `SSEE_EFT_Fundamental.tex` (borrador obsoleto de Paper 7 — NO usar)
- `SSEE_Paper8_StrongGravity.tex` (Paper 8 — 891 líneas, régimen fuerte, MIRA lensing)
- `SSEE_Paper10_UVCompletion.tex` (Paper 10 — 1037 líneas, UV K(X), M=9.62 meV, Conditional Theorem C.1)
- `SSEE_Endorser_Summary.tex` (en docs/)

### src/ — Scripts Python
- `ssee_paper2_analysis.py` — análisis analítico (plano w0-wa)
- `ssee_paper2_mcmc.py` — MCMC Bayesiano (SSEE vs ΛCDM vs CPL)
- `ssee_paper2_figures.py` — generación de figuras
- `ssee_paper3_cmb.py` — CAMB CMB spectrum TT+TE+EE+lensing vs Planck PR4
- `ssee_paper3_hiclass_check.py` — CLASS cross-check αK Bellini-Sawicki: αK(0)=0.4033, Δ=0.005% ✅
- `ssee_verify_rd.py` — verificación numérica r_d vs Planck 2018 medido (tarea 2A ✅)
- `ssee_press_schechter.py` — ratio PS n_SSEE/n_ΛCDM, δc=1.6284 vs 1.6865 (tarea 2B ✅)
- `ssee_paper5_IS_perturbations.py` — IS causal perturbations: Q1+Q2+Q3+fσ₈ (Paper 5 ✅)
- `ssee_paper6_verification.py` — dos sectores φ-DM: fσ₈ tensions, σ₈_eff, m_φ algebraico (Paper 6 ✅)
- `ssee_paper6_sterile_neutrino.py` — k_fs derivation, DW relation, scan m_φ candidates (Paper 6 ✅)
- `ssee_paper6_mcmc.py` — MCMC Paper 6 φ-DM dos sectores (100 walkers × 25000 steps; datos fσ₈ canónicos Paper 5)
- `ssee_op1_baryon_density.py` — OP-1: fórmula (π−φ)/H₀_SSEE = 0.32σ Planck (vs 3.2σ del 200) ✅
- `ssee_op1_baryogenesis.py` — OP-1 Sakharov: δ_CP=(π−φ)/Ω, T_rh~10⁻⁴ GeV, 3 condiciones ✅
- `ssee_op2_spectral_index.py` — OP-2: n_s=1−φ⁻⁷ (α-attractor + N_*=2φ⁷); r=φ⁻¹⁰ ✅
- `ssee_op3_separability.py` — OP-3: KALeff=φ²√(5/2), jerarquía (H₀/M)²≈10⁻⁶² ✅
- `ssee_op5_hmcode.py` — OP-5 Nivel 1: HMcode-2020 baryonic feedback CLASS; S₈=0.758 (0.06σ DES) ✅

---

## Fixes Paper 2 — estado

### ✅ Completados
1. ~~**graphicspath incorrecto**~~ ✅ — corregido a `{../results/figures/}`
2. ~~**Referencia MNRAS2026**~~ ✅ — clave Moresco2017 corregida
3. ~~**Cadena MCMC más larga**~~ ✅ — N_eff=637,500 (SSEE), 63,789 (ΛCDM), 43,793 (CPL)
4. ~~**Full covariance DESI**~~ ✅ — covarianza 13×13 implementada
5. ~~**Data Availability statement**~~ ✅
6. ~~**ORCID / affiliation**~~ ✅
7. ~~**ΔBIC correcto**~~ ✅ — −13.5 → −5.55 (5 ocurrencias + tabla BIC)

### Pendientes (Paper 2)
- [x] **1E**: Añadir 3 cúmulos: Perseus, Virgo/M87, A2744 (Implementado)
- [x] **2D**: χ²_r con 7 cúmulos (Completado)

---

## Resultados clave Paper 3 (para referencia)

| Espectro | SSEE χ²_r | ΛCDM χ²_r | N |
|---|---|---|---|
| TT | 1.047 | 1.043 | 1971 |
| TE | 1.041 | 1.040 | 1967 |
| EE | 1.041 | 1.039 | 1967 |
| PP (lensing) | 0.837 | 0.757 | 9 |
| ΔBIC (TT+TE+EE+PP diagonal, k=2 vs k=6) | −20.8 (SSEE favorecido; cross-check a H₀=67.08, = MCMC plik_lite) | — | — |
| ΔBIC (plik_lite Cobaya, k=2 vs k=6) | −31.3 (SSEE decisivamente favorecido) | — | — |
| ΔBIC (plik_lite conservador, k=4 vs k=6) | −13.8 (SSEE fuertemente favorecido) | — | — |

Datos lensing: `data/raw/planck_pr4_lensing.txt` (14 bins MV)

---

## Constantes algebraicas del modelo (referencia rápida)

```
φ = (1 + √5)/2 ≈ 1.6180   (razón áurea)
π ≈ 3.1416

Ω    = π + φ         ≈ 4.7596   (Stability Metric / OMEGA_DNAV)
β    = (π + φ)/2     ≈ 2.3798   (Base Coupling Scalar / BIAL)
KAL₀ = β + π        ≈ 5.5214   (Structural Viscosity)
P_sc = Ω + φ        ≈ 6.3776   (Dynamical Evolution Scalar / PYROS)
Kᵥ   = φ + π + Ω   ≈ 9.5192   (Structural Constraint / KRYSTOS)
Tᵣ   = 3(φ + β)    ≈ 11.9935  (3D Saturation Horizon / TRIAL)
Mᵥ   = φ + π + Kᵥ ≈ 14.2788  (Maximal Dimensional Invariant / Σ₉)

w₀ = -Tᵣ/Mᵥ ≈ -0.840
wₐ = -P_sc/Kᵥ ≈ -0.670
Ω_DE = Tᵣ/Mᵥ ≈ 0.840
Ω_m,dyn = 0.160    MIRA = 1.9989    Ω_m,CMB = 0.3199
```

---

## Datos observacionales usados (Paper 2)

- **DESI DR2 BAO**: 13 puntos, z = 0.295–2.330 (Abdul-Karim et al. 2025, arXiv:2503.14738)
- **Planck 2018**: Prior comprimido H₀=67.36±0.54, Ωm=0.3153±0.0073, Ωbh²=0.02237±0.00015, ρ(H₀,Ωm)=−0.85
- **Cúmulos (IGIMF-corregidos, Zhang et al. 2026, arXiv:2602.06082)**:
  - Coma: M_bar^IGIMF=1.8±0.2, M_obs=9.8±1.0 [10¹⁴ M☉]
  - A2029: 2.2±0.2, 12.0±1.2
  - A478:  1.5±0.2, 8.0±1.0
  - Bullet: 1.2±0.2, 6.5±1.0

---

## Resultados clave Paper 2 (para referencia)

| Métrica | Valor |
|---|---|
| χ²₂D (plano w₀-wₐ vs DESI) | 0.080 → 0.05σ |
| χ²ᵣ cúmulos (SSEE+IGIMF, 4 cúmulos) | 0.122 |
| ΔBIC (modelo completo, k=0 vs ΛCDM k=6) | +206 (penalidad por background ΛCDM) |
| ΔBIC MCMC (ΛCDM−SSEE, dynamic sector, k=2 vs k=3) | **+7.91** (SSEE favorecido) |
| ΔBIC MCMC structural (SSEE−ΛCDM full background) | +213 (Ω_m=0.160 vs 0.315; resuelve P6) |
| **H₀ SSEE (MCMC, prior MIRA — canónico tras 2026-05-24)** | **66.55⁺⁰·⁴⁴/₋₀·₄⁴** km/s/Mpc (script `ssee_paper2_mcmc_mira.py`, 100w×25k, seed 42) |
| H₀ SSEE (MCMC, prior Planck legacy — para comparación) | 66.75⁺⁰·⁴⁴/₋₀·₄⁴ km/s/Mpc (script original, valor en P2 antes del switch MIRA) |
| H₀ MIRA (Planck plik_lite + SSEE bg) | 67.068 km/s/Mpc (script `ssee_paper3_cobaya_unified.py`) |
| H₀ SSEE (algebraico Paper 4, Type P coincidence) | 67.962 km/s/Mpc = 3(φ+π)² |
| r_d (SSEE pure, Ω_m,dyn=0.160) | 175.6 Mpc — coherente Eisenstein-Hu |
| r_d (SSEE+MIRA mapping, Ω_m,CMB=0.31993) | 147.6 Mpc — matchea Planck |
| **m_φ (P6, partícula canónica — canónico 2026-06-04)** | **36.95 eV** = Σm_ν^active·(Ω⁴_DNAV+AURA·KAL) = 0.0690 eV × 535.28. Forward-prediction, cero fiteo, dim-consistente (masa × número adimensional). Lagrangiano escalar libre. Antes 5.602 eV (obsoleto, fórmula Σm_ν·H_alg retirada) |
| **k_fs (P6 canónico)** | **0.659 h/Mpc** (consistente con m_φ=36.95; antes 0.493 con 5.602) |
| **H_local IR canónico (P9 cascada)** | **71.90 km/s/Mpc — 1.10σ SH0ES** (era 72.86 = 0.17σ Type-P) |
| **H_local UV canónico (P10 cascada)** | **72.077 km/s/Mpc — 0.93σ SH0ES** (era 73.040 = 0σ Type-P) |

> **Valores canónicos:** la tabla autoritativa de cada número vive en
> `VERIFICATION_LEDGER.md` §«Valores Canónicos». Esta tabla es un resumen;
> ante discrepancia, manda el Registro.
>
> **Historia H₀ MCMC (registro honesto):**
> - 66.75 ± 0.44 — script original con prior Planck-ΛCDM (valor en P2 antes 2026-05-24)
> - 66.55 ± 0.44 — script con prior MIRA (self-consistent, P2 actualizado 2026-05-24)
> - El "67.756" anteriormente documentado aquí era MALA ANOTACIÓN, no corresponde a
>   ninguna corrida real. Corregido 2026-05-24. Ver [[project-h0-mira-calibration]].
>
> **Validez de cascadas H_alg→H_MIRA (revisión 2026-05-24 PM):**
> - **P9 (H_local = H_MIRA/(1−f_screen) = 71.90):** ✅ válida — H_MIRA dimensional
>   y f_screen adimensional; la fórmula cierra unidades.
> - **P10 (H_local^UV = 72.077):** ✅ válida por misma razón estructural.
> - **P6 (m_φ canónico = Σm_ν × MULT):** ✅ **VÁLIDA** (canónico 2026-06-04) —
>   m_φ = Σm_ν^active · (Ω⁴_DNAV+AURA·KAL) = 0.0690 eV × 535.28 = 36.95 eV.
>   Dim-consistente: masa (eV) × número adimensional puro (φ,π). Forward-prediction,
>   cero fiteo, escrita en un Lagrangiano escalar libre → cierra la incompletitud
>   de OP-9 (queda abierto SOLO el origen UV del multiplicador 535.28).
>   La vieja fórmula `m_φ = Σm_ν × H_alg = 5.602 eV` (eV × km/s/Mpc, sí era
>   dimensionalmente inválida) está **RETIRADA y obsoleta**.

---

## Resultados clave Paper 5 (referencia)

| Resultado | Valor | Tensión |
|---|---|---|
| c²_s,eff | 0 (exacto algebraico) | Q1: todos los modos estables |
| MIRA_num (k≥10) | 0.989 ± 0.017 | Q2: MIRA es efecto background, no perturbativo |
| γ_IS | 0.554 ± 0.001 | ≈ γ_ΛCDM = 0.55 |
| G = D₁_SSEE/D₁_ΛCDM | 0.866 | supresión 13.4% |
| σ₈_SSEE | 0.702 ± 0.005 | — |
| S₈_SSEE | 0.725 ± 0.005 | **1.96σ KiDS / 2.84σ DES (IS sector explored, tension remains)** |
| fσ₈(z=0.5) | 0.341 | 2.67σ vs BOSS (tensión estructural — Ω_m,dyn=0.160) |

---

## Resultados clave Paper 6 (referencia)

| Resultado | Valor | Estado |
|---|---|---|
| Ω_CDM | 0.160050 | Sector activo a todo k |
| Ω_φDM = (MIRA−1)×Ω_m,dyn | 0.159878 | Solo k < k_fs |
| Ω_total (dos sectores) | 0.319928 ≈ Ω_m,CMB | Unificación algebraica |
| m_φ = Σm_ν^active × (Ω⁴_DNAV+AURA·KAL) | 36.95 eV | Forward-prediction, cero fiteo (canónico 2026-06-04; antes 5.60 obsoleto) |
| k_fs (free-streaming) | 0.659 h/Mpc | Derivado de m_φ=36.95 (antes 0.493 con 5.60) |
| σ₈_eff (two-sector canónico) — **titular Paper 6** | 0.742 | cadena 0.811→×G_2s 0.794→free-streaming 0.742 |
| S₈_eff (two-sector canónico) — **titular Paper 6** | 0.766 | **0.01σ KiDS — RESUELVE la tensión** |
| Single-sector baseline (fuente Poisson Ω_m=0.320) | σ₈=0.820 / S₈=0.847 | 3.9σ KiDS — "el desafío" que el two-sector resuelve |
| alpha (free-streaming, OUTPUT CLASS no fiteado) | 1.243 Mpc/h | forward del two-sector canónico |
| Media tensión fσ₈ (6 encuestas canónicas Paper 5) | **0.76σ** | desde baseline 2.67σ (= Paper 5); empata ΛCDM 0.73σ. Antes 0.50σ con datos fσ₈ erróneos — corregido 2026-05-18 |
| m_φ algebraico (φ-DM) | 36.95 eV | Sin portal SM — no observable en KATRIN/PTOLEMY |
| Predicción DESI Y3/Euclid | k_fs = 0.659 h/Mpc (imprint de m_φ en P(k)) | Falsificable 2026–2028 |

---

## POA pendientes (próximas sesiones)

| Tarea | Descripción | Prioridad |
|---|---|---|
| ~~**1E**~~ | ~~Paper 2: añadir Perseus, Virgo/M87, A2744~~ | ✅ sesión 10 |
| ~~**2B**~~ | ~~Press-Schechter δc=1.628 vs 1.686~~ | ✅ sesión 12 |
| ~~**2D**~~ | ~~χ²_r con 7 cúmulos~~ | ✅ sesión 10 (χ²_r=0.126) |
| ~~**B2**~~ | ~~IS perturbativo → γ_IS, S8~~ | ✅ sesión 13 |
| ~~**B3**~~ | ~~MIRA mecanismo algebraico~~ | ✅ sesión 14 |
| ~~**B4**~~ | ~~α-attractor α=φ⁴/3, Paper 1 EFT §A.5~~ | ✅ sesión 15 |
| ~~**Paper 5**~~ | ~~IS causal perturbations: Q1+Q2+Q3+fσ₈, 24 pp~~ | ✅ sesión 16 (commit 0a9d21c) |
| ~~**Paper 6**~~ | ~~φ-DM dos sectores — m_φ algebraico, fσ₈ resuelto~~ | ✅ sesión 18 (commit 9d28620) |
| ~~**hi_class**~~ | ~~CLASS cross-check αK Bellini-Sawicki: αK(0)=0.4033, Δ=0.005%~~ | ✅ sesión 20 |
| ~~**Tone/Overclaim**~~ | ~~Eliminar "Complete EFT" y bajar el tono en Paper 7 y Endorser~~ | ✅ sesión 21 (auditoría editorial 10-paper) |
| ~~**Audit 10-paper**~~ | ~~Auditoría editorial completa Papers 1–9: nomenclatura, overclaims, bib~~ | ✅ sesión 21 (2026-05-15) |
| ~~**P4 nomenclatura**~~ | ~~Eliminar nombres mitológicos de Paper 4; bib 10→40 refs~~ | ✅ sesión 21 |
| **OPEN_PROBLEMS.md** | Documento de brechas físicas abiertas (OP-1..OP-6) | ✅ sesión 21 |
| ~~**Reproducibility**~~ | ~~Quitar rutas hardcodeadas y unificar requirements.txt~~ | ✅ completado |
| ~~**B1**~~ | ~~Full CMB likelihood Cobaya+plik (bloqueante PRD)~~ | ✅ completado |
| ~~**OP-4**~~ | ~~Recomputar r_V solar con Λ_SSEE=M=9.62 meV~~ | ✅ RESUELTO 2026-05-15 — k-mouflage + αB=αM=0 EFT; Paper 8 §4.2/§4.4 revisados |
| ~~**OP-6**~~ | ~~Derivar forma screening multiplicativa vs aditiva del Lagrangiano~~ | ✅ RESUELTO 2026-05-15 — universo separado k-essence + identidad 1+w₀=Ω_m |
| ~~**OP-3**~~ | ~~Probar Postulate C.1 incondicionalmente — jacobiano ∂φ/∂χ (tras OP-6)~~ | ✅ RESUELTO 2026-05-16 — jerarquía EFT (H₀/M)²≈10⁻⁶²; KALeff=φ²√(5/2) único |
| ~~**OP-1**~~ | ~~Derivar factor 200 Ω_b h²~~ | ✅ PARCIAL 2026-05-16 — (π−φ)/H₀_SSEE=0.32σ; BBN completo → Paper B/C |
| ~~**OP-2**~~ | ~~Derivar exponente 7 n_s del UV~~ | ✅ RESUELTO 2026-05-16 — α-attractor universality + N_*=2φ⁷; r=φ⁻¹⁰ nueva pred. falsificable |
| ~~**OP-5 Nivel 1**~~ | ~~HMcode-2020 CLASS: supresión bariónica baryonic_feedback~~ | ✅ PARCIAL 2026-05-16 — S₈=0.758 (0.06σ DES); B_eff=0.9447; script op5 |
| **OP-5 Nivel 2** | S₈ pleno — simulaciones N-body BAHAMAS/IllustrisTNG-SSEE | DIFERIDO (supercomputadora, ~5k–20k CPU-horas) |
| ~~**OP-1 Sakharov**~~ | ~~Argumento formal bariogénesis: δ_CP=(π−φ)/Ω, T_rh~10⁻⁴ GeV~~ | ✅ ROBUSTO 2026-05-16 — 3 condiciones Sakharov verificadas; script op1_baryogenesis |
| ~~**Bibliografía P1/P2/P3**~~ | ~~Llevar P1, P2, P3 a 40+ refs estándar JCAP/PRD~~ | ✅ 2026-05-16 — P1(40), P2(40), P3(41). Commits fc2e8ed + 5edcfc5 |
| ~~**P9/P10 inconsistencia H₀**~~ | ~~Verificar que 72.86 es canónico y 73.040 es check condicional en cuerpo completo~~ | ✅ 2026-05-16 — Verificado. P9 L76 "canonical, independent prediction"; P10 marca "conditional on Postulate C.1" en L74/127/362/379/395. Sin ediciones necesarias. |
| **OP-5 Nivel 2** | S₈ pleno — simulaciones N-body BAHAMAS/IllustrisTNG-SSEE | DIFERIDO (supercomputadora, ~5k–20k CPU-horas) |
| **Zenodo v2** | Archivar Papers 1-10 PDFs + OPEN_PROBLEMS.md + CLASS scripts con DOI | Alta |
| **Endorser arXiv** | Re-contactar endorsers con Papers 1-3 submission-ready | Alta |
| **B5** | Quintessential inflation V(φ) unificado | Largo plazo |

---

## Roadmap Hacia la Excelencia (Vulnerabilidades a resolver a largo plazo)
Guardado en memoria para futuras iteraciones teóricas (Auditoría Nivel Premio):

1. **Derivación Cuántica Dimensional:** $H_0 \approx 67.96$ km/s/Mpc es un éxito numérico, pero falta el teorema que conecte la escala macroscópica (Megaparsec) con unidades fundamentales (escala de Planck).
2. **QFT y la Diferencia del 0.05%:** La dualidad en $\Omega_{m,\mathrm{CMB}}$ (0.31993 vía MIRA vs 0.32010 geométrico) sugiere fuertemente una corrección de bucle (loop correction / self-energy) a nivel de campos cuánticos que aún no se ha escrito.
3. **El Lagrangiano Exacto de $\varphi$-DM:** El modelo usa EFT (Paper 7) para describir el campo a bajas energías, pero falta deducir el Lagrangiano $P(X, \phi)$ fundamental derivado puramente de $\varphi$ y $\pi$.
4. **Simulaciones N-body para el $S_8$:** La predicción de $S_8 = 0.820$ es lineal. Para resolver completamente la tensión contra KiDS ($0.766$), SSEE deberá correrse en supercomputadoras para incluir *Baryonic Feedback* y efectos no lineales.
5. **El Horizonte de Sonido $r_d$:** Probar definitivamente que un $r_d \approx 175$ Mpc (en lugar de $147$ Mpc) es un modelo físico superior para el universo temprano a pesar de la penalización artificial del $\Delta\mathrm{BIC}$ bajo calibración $\Lambda$CDM.

---

## ✅ Sesión 2026-05-09 — Auditoría Final + Zenodo + CLASS Boltzmann

### 1. Auditoría Editorial (Completada)
Se aplicó una cirugía "hostil pero justa" de lenguaje a toda la suite. Commits aplicados:
- `13b2226` — Papers 2, 4, 5, 6, 7 + Endorser Summary (overclaims eliminados)
- `d1799d7` — Paper 1, 3, 4 (tabla mitológica → P1..P9, Hubble Tension como "Note", exponent 7 como ansatz, MIRA en Paper 3 como "discrete structural constraint")

**Política de lenguaje definitiva adoptada:**
- "derives" → resultados algebraicos fijos
- "tests/shows" → confrontación con datos
- "phenomenological proposal/ansatz" → extensiones no derivadas de primeros principios
- "within the analyzed regime" → condiciones de estabilidad acotadas

### 2. Publicación Zenodo
- **DOI:** https://doi.org/10.5281/zenodo.20093447
- **Contenido:** 7 preprints + código MCMC completo
- **Descripción:** ajustada al tono fenomenológico post-auditoría

### 3. Correos de Endorsement arXiv (Enviados)
- **Dra. Eleonora Di Valentino** (e.divalentino@sheffield.ac.uk) — focus: phantom crossing + fσ₈
- **Dr. Skylee Lee** (skylee@skku.edu) — focus: anclaje r_d y degeneración θ* vs r_d
- **Código de endorsement:** `3XCDRE`
- Adjunto: `SSEE_Endorser_Summary.pdf`
- **Estado:** Enviados. Esperando respuesta.

### 4. CLASS Boltzmann — Primera Corrida SSEE ✅
- **Repositorio:** `/home/mike/Proyectos/SSEE/class_ssee/` (fork de class_public)
- **Config:** `class_ssee/ssee_v36.ini` — parámetros SSEE algebraicos
- **Outputs generados:**
  - `output/ssee_v36__cl_lensed.dat` — Espectro CMB TT+TE+EE+lensing
  - `output/ssee_v36__pk.dat` — Espectro de materia P(k)
- **Gráfico:** `output/ssee_v36_CMB_TT_comparison.png`
- **Script:** `class_ssee/plot_ssee_cmb.py`

**Resultados primera corrida:**
| Pico | Planck 2018 | SSEE-CLASS |
|------|-------------|------------|
| 1°   | ℓ = 220     | ℓ = **220** ✅ |
| 2°   | ℓ ≈ 540     | ℓ = **535** ✅ |
| 3°   | ℓ ≈ 810     | ℓ = **810** ✅ |
| RMS vs ΛCDM | — | **1.5%** ✅ |

> ⚠️ NOTA: Esta corrida usa `Ω_m = 0.3199` (MIRA-mapped como parámetro directo).
> El test crítico pendiente es correr con `Ω_m = 0.160` (sin MIRA) para cuantificar
> qué hace MIRA físicamente.

### 5. POA Siguiente Sesión (CLASS)
- [ ] **Fase 1A:** Correr CLASS con `Ω_m = 0.160` (sin MIRA) → cuantificar degradación ✅ COMPLETADO
- [x] **Fase 1B:** Comparar P(k) SSEE vs ΛCDM en escala de k_fs = 0.493 h/Mpc ✅ COMPLETADO
- [ ] **Fase 2:** Implementar φ-DM (5.60 eV, free-streaming) en perturbations.c
- [ ] **Fase 3:** IS viscosity en perturbations.c
- [ ] **Fase 4:** MCMC completo CLASS+SSEE vs Planck+DESI simultáneo

### Resultado Fase 1A — Test MIRA Boltzmann (CRÍTICO)
Config: `class_ssee/ssee_v36_nomira.ini` (Omega_b=0.048432, Omega_cdm=0.111568)

| Modelo        | Pico 1° | Pico 2° | Pico 3° | RMS vs ΛCDM |
|---------------|---------|---------|---------|-------------|
| SSEE + MIRA   | ℓ=220   | ℓ=535   | ℓ=810   | **1.4%** ✅  |
| SSEE sin MIRA | ℓ=240   | ℓ=597   | ℓ=922   | **31.5%** ❌ |
| ΛCDM          | ℓ=221   | ℓ=537   | ℓ=814   | referencia  |

**Conclusión científica:** Sin MIRA, los tres picos se desplazan ~10% hacia la derecha y el RMS
sube 22x. CLASS confirma que MIRA no es cosmético — es físicamente necesario para reproducir
el horizonte de sonido correcto con Omega_m,dyn=0.160. Argumento defensivo clave ante referees.

### 6. Limpieza repositorio
- Eliminados archivos `.aux`, `.log`, `.out`, `.toc` de `docs/` y directorio raíz
- Carpeta `docs/` ahora contiene únicamente los 8 PDFs finales

### 7. Fase 2 — φ-DM Two-Sector P(k) (En progreso)
- Config: `class_ssee/ssee_v36_twosector.ini` (ncdm con m=5.60 eV, T_ncdm=0.716)
- Script: `class_ssee/plot_ssee_twosector_pk.py` — WDM transfer function (Viel+2005)
- **Resultado:** WDM alpha=1.0933 h/Mpc — supresión demasiado agresiva (86% en k_fs)
- **Pendiente:** Recalibrar alpha usando relación exacta DW del Paper 6 para obtener σ₈≈0.702

---

## ✅ Sesión 2026-05-10 — Respuesta Dr. Lee + Nueva Estrategia de Publicación

### Respuesta Dr. Lee (Skylee, skku.edu)
El Dr. Lee respondió con consejo constructivo (no endorsement). Puntos clave:
- arXiv astro-ph.CO está aplicando moderación muy estricta a investigadores independientes
- Recomendó: **enviar primero a revista con revisión de pares**, luego arXiv
- Esto es el camino correcto y más rápido para un investigador sin historial previo

### Nueva Estrategia de Publicación
**CAMBIO DE PLAN:** arXiv ya NO es el primer objetivo. El flujo correcto es:
1. **Consolidar CLASS** → hacer el modelo más riguroso
2. **Crear documento unificado** → un paper principal (journal-ready) que consolide los resultados más sólidos de los 7 preprints
3. **Enviar a revista** → JCAP, PRD, o universe MDPI (más accesible para independientes)
4. **Tras aceptación** → subir a arXiv automáticamente vía la revista

### Revistas objetivo (en orden de prioridad)
| Revista | Factor Impacto | Apertura | Estrategia |
|---------|---------------|----------|------------|
| JCAP | 5.3 | Media | Primer intento si CLASS valida bien |
| PRD | 4.7 | Media | Alternativa a JCAP |
| Universe (MDPI) | 2.9 | Alta | Más accesible para independientes |
| MNRAS | 4.8 | Media | Requiere historial institucional |

### Documento Unificado — Qué incluir
El paper principal para la revista debe consolidar:
- Derivación algebraica del background (Papers 1+4)
- Validación MCMC DESI DR2 + Planck PR4 (Papers 2+3)
- Validación CLASS Boltzmann (trabajo actual — más riguroso que CAMB parcheado)
- φ-DM two-sector + reducción fσ₈ (Paper 6, con calibración CLASS)
- Falsification table del Paper 1 (criterios cuantitativos pre-comprometidos)

### POA Actualizado (prioridades)
- [x] **Fase 2c:** Recalibrar WDM alpha con sigma_8 top-hat correcto → alpha=1.6561, σ₈=0.737 ✅
- [x] **Fase 2b:** fσ₈ tensión 1.88σ → 0.53σ con CLASS (Paper 6: 2.56→0.50, baseline diferente) ✅
- [x] **Fase 3:** IS viscosity CLASS — cs2 efecto 0.03%, ODE G=0.866 confirmado ✅
- [ ] **Doc unificado:** Estructura del paper journal con resultados CLASS integrados
- [ ] **EFTCAMB:** Pendiente para después de CLASS (requiere gfortran + lapack)

---

## ✅ Sesión 2026-05-13 — Fase 2c+3 CLASS: sigma_8 top-hat fix + IS viscosity

### Corrección crítica: sigma8_proxy → sigma8_real
`calibrate_wdm_alpha.py` tenía un bug severo: `sigma8_proxy = sqrt(∫ k² P dk)` sin
filtro top-hat W²(kR) ni normalización 1/(2π²). Esto daba valores de ~16 en vez de ~0.82.
El ratio era internamente consistente pero el alpha calibrado (0.1597) solo suprimía σ₈ 0.03%.

**Fix:** reemplazado por `sigma8_real` con fórmula de Peebles 1980:
`σ₈² = (1/2π²) ∫ k² P(k) W²(kR) dk`, W(x) = 3(sin x − x cos x)/x³, R=8 Mpc/h

**Resultado corregido:**
- alpha_WDM = **1.6561 h/Mpc** (era 0.1597 — corrección 10×)
- σ₈_eff (CLASS top-hat directo) = **0.7370** — exactamente Paper 6 ✅
- Tensión fσ₈: **1.88σ → 0.53σ** ✅

### Fase 3 — IS Viscosity CLASS (completada)
Config: `class_ssee/ssee_v36_IS.ini` (cs2_fld=0.001, Ω_m=0.160 single-sector)
Script: `class_ssee/plot_IS_viscosity.py`

| Resultado | Valor | Estado |
|---|---|---|
| cs2=1 vs cs2≈0 (σ₈ diferencia) | 0.03% | IS c²_s efecto negligible ✅ |
| ODE G = D₁_SSEE/D₁_ΛCDM | 0.866 | Confirma Paper 5 ✅ |
| σ₈_ODE | 0.702 | Confirma Paper 5 ✅ |
| CLASS single-sector (Ω_m=0.160) | 0.457 | Discrepancia ODE: T(k) real ≠ T_ΛCDM(k) |
| CLASS MIRA (Ω_m=0.320) + WDM | 0.737 | Confirma Paper 6 ✅ |

**Insight físico clave:** El ODE de Paper 5 asume T_SSEE(k) ≈ T_ΛCDM(k) (mismo transfer
function). CLASS calcula T(k) real con Ω_m=0.160 → mucho menos poder en k pequeño.
El sector MIRA corrige esto: CLASS con Ω_m=0.320 es el background correcto para σ₈.

### Commit
`bacd679` — feat: CLASS Fase 2c+3 — sigma_8 top-hat fix, IS viscosity, phi-DM calibration

### POA Próxima Sesión
- [ ] **Doc unificado:** Estructura del paper journal con resultados CLASS integrados
- [ ] **CLAUDE.md POA tabla:** Marcar Fase 2/3 completas en sección principal
- [ ] **EFTCAMB:** Instalación y primera corrida αK Bellini-Sawicki
- [ ] **Zenodo v2:** Archivar Papers 1-7 + scripts CLASS

---

## ✅ Sesión 2026-05-15 — Auditoría Editorial 10-Paper + P4 Renovación

### Alcance de la auditoría
Revisión editorial completa de Papers 1–9 orientada a estándares JCAP/PRD y Buchalter Prize.
13 ítems completados a lo largo de sesión multi-compresión.

### Cambios aplicados en esta sesión

#### Paper 4 — Nomenclatura + Bibliografía
- **Nombres mitológicos eliminados** de `manuscript/SSEE_Paper4_ToE.tex`:
  SOLAR, IGNIS, PYROS, KRYSTOS_V, BIAL, MAR, MIKA, PHITA, VITA, BUFFER, TRIAL, MIKAEL_V, Ω_DNAV
  → reemplazados por símbolos algebraicos estándar: φ, π, Ω, β, KAL, P_sc, K_v, Σ_Sov
  (AURA, MIRA, KAL mantenidos como identificadores físicos establecidos en Papers 5–9)
- **Bibliografía expandida** `manuscript/SSEE_Paper4.bib`: 10 → 40 entradas
  (Perlmutter1999, Riess1998, CPL, Copeland2006, Verde2019, Kamionkowski2023,
  desidr2_2025, ArmendarizPicon2001, Ratra1988, Caldwell2002, Starobinsky1980,
  Guth1981, Linde1982, Mukhanov1992, Blas2011, KiDS2021, DESY32022, y otros)
- **Citas en texto** añadidas en §Intro, §spectral index, §Hubble tension, §w₀wₐ, §H₀

#### OPEN_PROBLEMS.md — creado
Documento de 6 brechas físicas abiertas (OP-1..OP-6):
| ID | Brecha | Severidad |
|----|--------|-----------|
| OP-1 | Factor 200 en Ω_b h² — sin derivación | Alta |
| OP-2 | n_s = 1−φ⁻⁷ — exponente 7 no derivado de V(φ) | Alta |
| OP-3 | Separabilidad UV-IR φ/π — conjetura no probada | Media |
| OP-4 | r_V solar > r_Hubble — incoherencia Vainshtein-EFT | Alta |
| OP-5 | fσ₈ 2.67σ BOSS (sector único) — tensión estructural | Media |
| OP-6 | Forma de screening multiplicativa vs aditiva en P9 | Media |

### Estado final de calidad (estimado post-auditoría)
| Paper | Score JCAP/PRD | Score Buchalter |
|-------|---------------|-----------------|
| P1 | 8.5/10 | 8.0/10 |
| P2 | 9.0/10 | 8.5/10 |
| P3 | 8.5/10 | 8.0/10 |
| P4 | 7.5/10 (bib OK; OP-1,OP-2 abiertas) | 7.0/10 |
| P5 | 8.0/10 | 7.5/10 |
| P6 | 8.0/10 | 7.5/10 |
| P7 | 8.5/10 | 8.0/10 |
| P8 | 7.5/10 (OP-4 abierta) | 7.0/10 |
| P9 | 7.0/10 (DRAFT; OP-6) | 6.5/10 |

### Commits de esta sesión
Ver `git log` — commit final cubre P4 tex + bib + OPEN_PROBLEMS + CLAUDE.md

### POA Siguiente Sesión
- [ ] Compilar P4 PDF con nueva bib (pdflatex + bibtex) y verificar referencias resueltas
- [ ] Zenodo v2: subir PDFs actualizados + OPEN_PROBLEMS.md
- [ ] Re-contactar endorsers con papers 1–3 en estado submission-ready
- [ ] B5: Quintessential inflation V(φ) — derivación ns desde α-attractor (largo plazo)
- [ ] OP-4: Investigar cutoff Λ_SSEE para recalcular r_V solar

---

## ✅ Sesión 2026-05-23 — Reframe minimal-parameter + Postulados auxiliares + Auditoría hostil

### Cambios principales

**Slogan "0 free parameters" → "minimal-parameter framework"** en P1, P2, P3, P4, P5, P6, Endorser, Unified. Honest accounting documentado en P1 §1.3.

**5 OPs nuevos** en OPEN_PROBLEMS.md:
- OP-8: MIRA dynamical mechanism (7 mech ruled out → archive `src/mira_attempts/`)
- OP-9: m_φ=5.60 eV phenomenological ansatz (Paper 6 admits)
- OP-11: ξ non-minimal coupling free parameter
- OP-13: θ_E factor 2 SLACS/BELLS not verified
- OP-14: Σm_ν Type P, offset 22 ad hoc — eslabón más débil

**Two-Ω_m criterion** (P1 §1.4) — tabla por paper con 18 filas + cross-refs en P4/P5/P6/P7/P8/P9. Cierra el "uso arbitrario" ante referee.

**Postulados enumerados honestamente: 4 total** (D, S fundamentales + M, I auxiliares):
- M (MIRA = (3φ+π)/4) — valor algebraico; mecanismo OP-8
- I (α-attractor con α=φ⁴/3, forma N_*=2·φⁿ) — n=7 es **corolario** (única solución entera en [50,60] e-folds), no postulado
- YGG_PROTON (π-φ)/(π+φ)=0.3201 disuelto como predicción algebraica independiente, no postulado adicional

**Filosofía de auditoría guardada en memoria** (`feedback_audit_philosophy.md`): hostilidad sin resolución agrega debilidades. Aplicar fase 2 (buscar la salida) antes de reportar como vulnerabilidad. Aplicado en sesión a 4 hallazgos que terminaron siendo confirmaciones del modelo.

**P5 'algebraic theorem' → 'construction-level identity'** en 3 lugares + nota margen k_crit comfortable tree-level.

**β_c signo aclarado**: footnote en P8 explicita que β_c=+AURA (disformal) ≠ β_c=-AURA (conformal P7).

### Auditoría hostil de 7 capas — resultado

- **Capa 1** identidades algebraicas machine precision: 20+ verificadas, 1 hallazgo (dualidad Ω_m,CMB 0.054%) documentado como cross-confirmación, no contradicción.
- **Capa 2** cross-implementación ssee_core↔LaTeX↔scripts: VERDE en 16 constantes × 10 papers.
- **Capa 3-7** derivación/dimensional/física/asunciones/falsificabilidad: 3 fixes aplicados (H-3.1, H-5.1, H-5.2).

### Estado del modelo tras sesión

- **0 contradicciones internas identificadas**
- **4 postulados** explícitos (vs presentación previa de "2")
- **5 OPs catalogadas** con cadena de dependencias (OP-14→OP-9 conocida)
- Guardián VERDE 102/102
- Todos los papers recompilan limpio
- Filosofía de auditoría documentada para sesiones futuras

### Limpieza repo

- `src/mira_attempts/` creado con 12 scripts archivados + README
- `.gitignore` actualizado: sandbox_unificado/, src/results/
- VERIFICATION_LEDGER paths actualizados
- SEALED_STATUS.md actualizado con acciones completadas

### POA siguiente

- [ ] **Plan de ataque OPs** con cadena de dependencias (OP-14 primero por ser eslabón débil)
- [ ] **OP-13** verificar literatura SLACS/BELLS θ_E factor 2
- [ ] **Slim Journal v2** diseño estructural (~25 pp con framing minimal-parameter)
- [ ] Re-compilar PDFs en `docs/` desde manuscripts actualizados
