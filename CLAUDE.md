# SSEE — Contexto del Proyecto para Claude

## Qué es este proyecto

**Structural Self-Energy Expansion (SSEE-V3.6)** — modelo cosmológico de energía oscura de cero parámetros libres, derivado algebraicamente de φ (razón áurea) y π. Autor: Mike Edison Almeida Vallejo.

Cuatro papers:
- **Paper 1**: Framework teórico + EFT + Predictive Register (arXiv-ready)
- **Paper 2**: Validación Bayesiana MCMC DESI DR2 + Planck + cúmulos (arXiv-ready)
- **Paper 3**: Confrontación CMB Planck PR4 TT+TE+EE+lensing (arXiv-ready)
- **Paper 4**: Derivación algebraica CMB desde φ,π — Nine Sovereignties (preprint)

---

## Estado actual de los documentos (al 2026-04-26)

### docs/ — PDFs compilados (nombres reales)

| Archivo | Contenido | Estado |
|---|---|---|
| `SSEE_Paper1_Framework.pdf` | Paper 1 (10 p) | ✅ Completo |
| `SSEE_Paper2_draft.pdf` | Paper 2 (20 p) | ✅ Submission-ready |
| `SSEE_Paper3_CMB_Confrontation_v2.pdf` | Paper 3 (17 p) | ✅ Completo v2 |
| `SSEE_Paper4_ToE.pdf` | Paper 4 (13 p) | ✅ Preprint actualizado |

### manuscript/ — Fuente LaTeX
- `SSEE_Paper1_Framework.tex` + `SSEE_EFT_section.tex`
- `SSEE_Paper2_draft.tex`
- `SSEE_Paper3_draft.tex` + `ssee_paper3.bib`
- `sandbox_unificado/SSEE_Paper4_ToE.tex` (submodule)

### src/ — Scripts Python
- `ssee_paper2_analysis.py` — análisis analítico (plano w0-wa)
- `ssee_paper2_mcmc.py` — MCMC Bayesiano (SSEE vs ΛCDM vs CPL)
- `ssee_paper2_figures.py` — generación de figuras
- `ssee_paper3_cmb.py` — CAMB CMB spectrum TT+TE+EE+lensing vs Planck PR4
- `ssee_verify_rd.py` — verificación numérica r_d vs Planck 2018 medido (tarea 2A ✅)
- `ssee_press_schechter.py` — ratio PS n_SSEE/n_ΛCDM, δc=1.6284 vs 1.6865 (tarea 2B ✅)

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
- [ ] **1E**: Añadir 3 cúmulos: Perseus (Simionescu 2011), A2142 (Tchernin 2016), A1689
- [ ] **2D**: χ²_r con 7 cúmulos (tras completar 1E)

---

## Resultados clave Paper 3 (para referencia)

| Espectro | SSEE χ²_r | ΛCDM χ²_r | N |
|---|---|---|---|
| TT | 1.062 | 1.043 | 1971 |
| TE | 1.053 | 1.040 | 1967 |
| EE | 1.040 | 1.039 | 1967 |
| PP (lensing) | 0.730 | 0.757 | 9 |
| ΔBIC (TT, k=0 vs k=6) | −6.9 (SSEE favorecido) | — | — |
| ΔBIC (TT+TE+EE+PP, upper bound) | +13.7 | — | — |
| ΔBIC (plik_lite, k=1 vs k=6) | −40.3 | — | — |

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
| ΔBIC (sector dinámico aislado, k=1) | −5.55 (SSEE favorecido) |
| H₀ SSEE (MCMC best-fit) | 66.75⁺⁰·⁴⁴₋₀.₄₄ km/s/Mpc |
| H₀ SSEE (algebraico Paper 4) | 67.96 km/s/Mpc = 3(φ+π)² |
| r_d (SSEE+MIRA, CAMB) | 147.156 Mpc → **0.25σ** de Planck 2018 ✅ |

---

## POA pendientes (próximas sesiones)

| Tarea | Descripción | Prioridad |
|---|---|---|
| **1E** | Paper 2: añadir Perseus, A2142, A1689 | Alta |
| ~~**2B**~~ | ~~Press-Schechter δc=1.628 vs 1.686~~ | ✅ sesión 12 |
| **2D** | χ²_r con 7 cúmulos (tras 1E) | Depende de 1E |
| **B1** | Full CMB likelihood Cobaya+plik (bloqueante PRD) | Semanas |
| **B2** | IS perturbativo completo → cs², S8 | Semanas |
| **B3** | MIRA mecanismo físico formal | Largo plazo |
