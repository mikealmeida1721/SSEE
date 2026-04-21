# SSEE — Contexto del Proyecto para Claude

## Qué es este proyecto

**Structural Self-Energy Expansion (SSEE-V3.6)** — modelo cosmológico de energía oscura de cero parámetros libres, derivado algebraicamente de φ (razón áurea) y π. Autor: Mike Edison Almeida Vallejo.

Trilogía de papers:
- **Paper 1**: Framework teórico (publicado como preprint)
- **Paper 2**: Validación Bayesiana MCMC (draft en progreso → submission-ready)
- **Paper 3**: Confrontación CMB Planck PR4 (✅ completo — TT+TE+EE+lensing)

---

## Estado actual de los documentos (al 2026-04-21)

### docs/ — Documentos fuente

| Archivo | Contenido | Estado |
|---|---|---|
| `SSEE_Paper1_Framework_v3.6.pdf` | Paper 1 completo (5p) | ✅ Completo |
| `SSEE_Paper2_MCMC_Validation_draft.pdf` | Paper 2 draft (15p) | 🟡 Draft — ver fixes pendientes abajo |
| `SSEE_Paper3_CMB_Confrontation_v1.pdf` | Paper 3 completo (TT+TE+EE+lensing) | ✅ Completo — pendiente recompilar v2 |

### manuscript/ — Fuente LaTeX
- `SSEE_Paper2_draft.tex` — fuente completa del Paper 2
- `SSEE_Paper3_draft.tex` — fuente completa del Paper 3 (TT+TE+EE+PP)

### src/ — Scripts Python
- `ssee_paper2_analysis.py` — análisis analítico (plano w0-wa, sigma, tabla sensibilidad)
- `ssee_paper2_mcmc.py` — MCMC Bayesiano (SSEE vs ΛCDM vs CPL)
- `ssee_paper2_figures.py` — generación de figuras → results/figures/
- `ssee_paper3_cmb.py` — CAMB CMB spectrum TT+TE+EE+lensing vs Planck PR4

---

## Fixes pendientes en Paper 2 (.tex)

### Críticos (rompen compilación o son errores factuales):
1. **graphicspath incorrecto**: está `{../figures/}` debe ser `{../results/figures/}` (reorganizamos el repo)
2. **Referencia MNRAS2026**: clave engañosa — el paper de Moresco & Marulli es de 2017, no 2026

### Mejoras para submission-ready:
3. **Tabla H(z) CC**: la Fig 7 compara con Moresco 2022 pero no hay tabla con los datos reales de Cosmic Chronometers
4. **Más cúmulos**: solo tiene 4 (Coma, A2029, A478, Bullet de Zhang 2026). Candidatos: Perseus, Virgo/M87, A2744
5. ~~**Cadena MCMC más larga**~~ ✅ COMPLETADO — N_eff=637,500 (SSEE), 63,789 (ΛCDM), 43,793 (CPL)
6. ~~**Full covariance DESI**~~ ✅ COMPLETADO — covarianza 13×13 implementada
7. ~~**Data Availability statement**~~ ✅ YA EXISTÍA
8. ~~**ORCID / affiliation**~~ ✅ YA EXISTÍA

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

Datos lensing: `data/raw/planck_pr4_lensing.txt` (14 bins MV, fuente: Cobaya planck_supp_data_and_covmats)

## Pendientes Paper 3

- [ ] Recompilar PDF → `docs/SSEE_Paper3_CMB_Confrontation_v2.pdf`
- [x] Subir Genesis 5.12 a Zenodo → DOI incluido (10.5281/zenodo.19679049)

---

## Constantes algebraicas del modelo (referencia rápida)

```
φ = (1 + √5)/2 ≈ 1.6180   (razón áurea)
π ≈ 3.1416

Ω    = π + φ         ≈ 4.7596   (Stability Metric)
β    = (π + φ)/2     ≈ 2.3798   (Base Coupling Scalar)
KAL₀ = β + π        ≈ 5.5214   (Structural Viscosity)
P_sc = Ω + φ        ≈ 6.3776   (Dynamical Evolution Scalar)
Kᵥ   = φ + π + Ω   ≈ 9.5192   (Structural Constraint)
Tᵣ   = 3(φ + β)    ≈ 11.9935  (3D Saturation Horizon)
Mᵥ   = φ + π + Kᵥ ≈ 14.2788  (Maximal Dimensional Invariant)

w₀ = -Tᵣ/Mᵥ ≈ -0.840
wₐ = -P_sc/Kᵥ ≈ -0.670
Ω_DE = Tᵣ/Mᵥ ≈ 0.840
Ω_m,eff = 0.160
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
| χ²ᵣ cúmulos (SSEE+IGIMF) | 0.122 |
| ΔBIC (modelo completo) | +206 (penalidad por background) |
| ΔBIC (sector dinámico aislado) | −13.5 (favorecido) |
| H₀ SSEE | 66.66⁺⁰·⁴⁷₋₀.₄₆ km/s/Mpc |
| Tensión Ωm | 21.3σ (background mismatch, no falsificación) |
| r_d SSEE | 175.6 Mpc vs 147.6 Mpc (ΛCDM) |
| r_d,eff (Folding Symmetry) | 147.2 Mpc (dentro de 1σ de Planck) |

---

## Próximos pasos (2026-04-22)

1. Zenodo deposit Genesis 5.12 → DOI → actualizar bib
2. Abstracts arXiv (250 palabras × 3 papers)
3. Recompilar Paper 3 PDF
4. Fixes pendientes Paper 2 (full covariance DESI, Data Availability, ORCID)
