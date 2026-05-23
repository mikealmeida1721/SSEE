# SEALED_STATUS — Auditoría paper por paper

**Creado:** 2026-05-22
**Última actualización:** 2026-05-23
**Propósito:** registrar para cada paper SSEE qué partes están **selladas**
(resuelven cosa sin depender de derivar MIRA), qué partes son **condicionales
a MIRA** (necesitan MIRA como input observacional sin mecanismo derivado), y
qué partes requieren **edición** antes de ser publicables.

## ✅ ACCIONES EDITORIALES CRÍTICAS — COMPLETADAS 2026-05-23

### 1. Two-Ω_m Criterion — escrito y referenciado
- **P1 §1.4** "The Two-Ω_m Criterion" creada con tabla por paper (los 10 papers).
- Cross-references añadidas en P4, P5, P6, P7, P8, P9 a la primera mención de Ω_m.
- Inversión de signo del paréntesis para el referee: "uso arbitrario" → "regla pre-registrada física".

### 2. Reframe del slogan "0 parámetros libres"
- P1, P2, P4, Endorser, Unified: "zero-free-parameter" → "minimal-parameter framework"
- Honest accounting: ~3 efectivos vs 6 ΛCDM, documentado en P1 §1.3.
- 5 OPs nuevos: OP-8 (MIRA mech), OP-9 (m_φ ansatz), OP-11 (ξ libre), OP-13 (θ_E SLACS), OP-14 (Σm_ν Type P).

### 3. Postulados auxiliares enumerados (Auditoría hostil → resolución)
- P1 §1.3 ahora declara honestamente **4 postulados**: D, S (fundamentales) + M, I (auxiliares).
- YGG_PROTON $(\pi-\varphi)/(\pi+\varphi)$ NO es postulado adicional: es predicción algebraica independiente que coincide con MIRA·Ω_m,dyn a 0.054% (ambas dentro de Planck 1σ).
- Postulate I: forma `2·φⁿ` es postulado; **n=7 es corolario** (única solución entera en ventana 50-60 e-folds).

### 4. P5 'theorem' → 'construction-level identity'
- 3 lugares matizados: el c²_s,eff=0 es identidad de construcción dentro de SSEE, no teorema sobre IS arbitrarios.
- Margen k_crit documentado como cómodo tree-level.

### 5. β_c signo aclarado P8 vs P7/P9/P10
- Footnote en P8 explicita que ahí β_c = +AURA (disformal), distinto del β_c = -AURA (conformal) de P7.

### 6. Scripts MIRA experimentales archivados
- 12 scripts movidos a `src/mira_attempts/` con README explicativo de los 7 mecanismos descartados.
- Paths en VERIFICATION_LEDGER actualizados.

---

## Estado global del modelo (recordatorio)

- **MIRA = (3φ+π)/4 ≈ 1.9989** es algebraicamente exacto.
- **MIRA es observacionalmente REQUERIDO**: sin él, CLASS muestra RMS CMB 31.5%
  vs 1.4% con MIRA (sesión 2026-05-09). El modelo no es viable sin MIRA.
- **MIRA NO está derivado dinámicamente**: 7 mecanismos descartados (ver
  V-L3-mira, V-L3-saturacion, V-L3-2Om en `VERIFICATION_LEDGER.md`).
- **Postura adoptada**: lectura B explícita. MIRA entra como input observacional
  con valor algebraico exacto; el mecanismo dinámico es **open problem**.

## Frame ontológico sobre materia oscura

Los Papers 5–9 contienen una densidad oscura no-bariónica Ω_dark ≈ 0.110.
Frame propuesto consistente:

> "SSEE postula un solo campo escalar φ. En régimen de baja energía (z<2)
> φ se comporta como dark energy. En régimen sub-galáctico o k > k_fs,
> oscila coherentemente y se comporta como materia oscura. No hay partícula
> exótica nueva (WIMP, axión QCD); hay un solo campo con dos fases dinámicas.
> Esto es scalar field dark matter / ultra-light DM, no exotic particle DM."

Pertinente para resolver aparente contradicción P1 ↔ P6/P7.

## Plantilla por paper

```
### Paper N — <título>
- **Numerical**: PASS / FAIL [refs ledger]
- **Ontological (DM frame)**: COMPATIBLE / NEEDS REWORDING [palabras problemáticas]
- **MIRA dependency**: NONE / EMPIRICAL INPUT / DERIVED CLAIM
- **Status**: SEALED / NEEDS EDIT / EXPLORATORY
- **Notas**: [hallazgos clave]
```

## Estado actual

### Paper 1 — Framework  *(auditado 2026-05-22)*
- **Numerical**: PASS (guardián 102/102 OK: W0, WA, KAL0, MIRA, AURA, n_s, Ω_m,dyn, OMEGA_B_H2, H0_ALG todos OK)
- **Ontological (DM frame)**: ANTI-DM-PARTICLE STRICT
  - L51 abstract: "without introducing free parameters or exotic dark-matter particles"
  - L364: clusters resuelven "without cold dark matter" vía IGIMF + KAL + ν
  - L413-421: Bullet cluster vía KAL(x), "without collisionless particles"
  - L423-429: MOND-like μ_SSEE(x) para escala galáctica
- **MIRA dependency**: EMPIRICAL INPUT (no derivado mecánicamente)
  - L122 tabla: "Pre-data" ratio Ω_m,CMB/Ω_m,dyn
  - L836-845 falsificación: incluye MIRA two-sector crossover
  - L955 ya admite "in early internal documents MIRA was loosely approximated"
- **Status**: **SEALED (numérico + álgebra + MIRA-as-input)** + **PENDING (frame DM)**
  - Lo que se sella ahora:
    1. Identidades algebraicas (w₀, wₐ, KAL₀, MIRA, AURA, n_s, Ω_b h², H₀_alg)
    2. MIRA tratado como input observacional con falsificación explícita
    3. Cluster χ²_r=0.122 vía IGIMF+KAL+ν (internamente coherente)
  - Lo que queda pendiente hasta auditar P6/P7/P8:
    1. ¿Se mantiene "no exotic DM particles" estricto? (Si P6 con m_φ=5.60 eV
       resulta sólido y falsificable, P1 podría flexibilizar L51)
    2. ¿La interpretación de Ω_dark=0.110 final?
    3. Decisión Opción A / B sobre frame ontológico
- **Hallazgo crítico**: P1 ≠ frame "scalar field two-phase DM" que propuse
  inicialmente. P1 usa **modificación de gravedad efectiva** vía KAL(x), no
  campo escalar DM con masa. Esto **contradice directamente Paper 6** si
  P6 mantiene φ-DM con m=5.60 eV como sector CDM-like.
- **Acción pendiente sobre MIRA**: en próxima edición, agregar §"MIRA Status"
  haciendo explícito que es input observacional, no derivado, con los 7
  mecanismos descartados listados (V-L3 del Registro).

### Paper 2 — MCMC DESI+Planck+Cúmulos  *(auditado 2026-05-22)*
- **Numerical**: PASS para invariantes (ΔBIC=+7.91, BIC=31.99, tensión 0.57σ)
- **Resolución del drift (auditado 2026-05-22)**:
  - Backup `mcmc_chains_professional_SSEE-V3.6_ckpt_Om160_backup.npz`
    (abril 21, anterior al commit 4892b53) contiene:
    - H₀ = 66.7505 ± 0.4431  (63.25M muestras post-burn-in)
    - Ω_b h² = 0.02269 ± 0.00048
  - **Coincide exactamente** con el texto del Paper 2 (L68/L481/L507)
  - **Veredicto**: el chain canónico de Paper 2 es Om160_backup
    (Ω_m=0.160 uniforme). El commit 4892b53 (dos-Ω_m, H₀=67.76) fue
    un experimento que **no se canonizó en el paper**. No es drift
    paper→script; es script que se adelantó al paper.
  - **HECHO (2026-05-22)**: revertido `ssee_paper2_mcmc.py` a OM_EFF_SSEE
    uniforme en E(z) y r_d. Comentario explícito agregado: "MIRA NO entra
    al fondo gravitacional (input observacional no derivado)".
  - **Razón histórica del drift**: experimento dos-Ω_m introducido por
    Claude en sesión anterior (commit 4892b53). El experimento agregaba
    ruido sin justificación (metía MIRA al fondo cuando MIRA no está
    derivado). Lectura B explícita resuelve: MIRA solo aparece donde
    entra empíricamente (P3 CMB peaks, P9 H₀ screening), NO en E(z).
- **Ontological**: CONTRADICCIÓN INTERNA
  - L118: "no exotic particles are introduced"
  - L606-610: "two dark matter sectors, Ω_CDM=0.16005 + Ω_φ-DM=0.15988"
    con m_φ=5.60 eV, k_fs=0.493 h/Mpc
  - Misma contradicción que P1↔P6, aquí dentro del mismo paper
- **MIRA dependency**: EMPIRICAL INPUT (acknowledged open)
  - L75, L83: "MIRA two-sector mechanism established in Paper 3"
  - L614: **honesto** — "a first-principles QFT derivation of the
    MIRA factor remains open"
  - L979: "the MIRA bridge is load-bearing"
- **Status**: **SEALED-NUMERICAL** + **PENDING (frame DM)**
  - Sellado ahora:
    1. H₀ = 66.75 ± 0.44 (Om160_backup) — coincide con el paper
    2. ΔBIC, BIC, χ²_r cúmulos, comparación 3 modelos
    3. Tratamiento de MIRA como "load-bearing pero open" (L614, L83)
  - Pendiente (frame DM):
    - Decisión global P1 sobre exoticidad de m_φ=5.60 eV
    - Unificar L118 ("no exotic particles") con L606-610 ("Ω_CDM + Ω_φ-DM")
  - Acción pendiente menor:
    - Revertir `ssee_paper2_mcmc.py` al setup viejo o documentar
      el _legacy como canónico
- **Lo que SÍ se sella ya**:
  1. Pipeline MCMC, posteriors, χ²_r cúmulos, ΔBIC dyn-sector −5.55
  2. Tratamiento de MIRA como "load-bearing pero open" (correcto)
  3. Comparación SSEE vs ΛCDM vs CPL como tres modelos

### Paper 3 — CMB Planck PR4  *(auditado 2026-05-22)*
- **Numerical**: PASS
  - (π−φ)(3φ+π)/[8(φ+π)] = 0.319928 ✓ (L55)
  - Tensión vs Planck = 0.63σ ✓ (L273)
  - ΔBIC = −20.8 ✓ (L666, guardián V-L4-05)
  - H₀ = 67.08±0.03 (plik_lite Cobaya MCMC propio de P3,
    diferente del 66.75 de P2 — son MCMCs de datos distintos, OK)
- **Ontological**: CONFLICTO HEREDADO de P6
  - L299: "without the need for cold dark matter" (en §methodology)
  - L59: "two-sector **dark matter** structure with $m_\phi=5.60$\,eV"
    (en intro, dependiendo de P6)
  - Mismo patrón que P2 L118 ↔ L606
- **MIRA framing**: **EL MEJOR DE LA SUITE**
  - L59: "discrete structural constraint, **rather than a continuously
    fitted parameter**"
  - L59: "**cannot be sourced** by classical Israel--Stewart fluid mechanics"
  - L74: explícito "MIRA-corrected claim" — condicional
  - L188-192: sección dedicada "The MIRA Two-Sector Solution"
  - Este es el framing que debería propagarse a P1, P2 cuando se editen
- **Status**: **SEALED-NUMERICAL** + **SEALED-MIRA-FRAMING** + **PENDING (frame DM)**
  - Si Paper 6 sobrevive auditoría con frame SFDM, P3 ya está alineado
  - Si Paper 6 cae, P3 mantiene resultado de fitting pero L59 debe reformularse

### Paper 4 — ToE algebraico  *(auditado 2026-05-22)*
- **Numerical**: PASS (5/5 identidades verificadas)
  - (π−φ)/(π+φ) = 0.320100 ✓ (L383: Ω_m geométrico)
  - KAL·Ω_b h²·n_s = 0.11926 ✓ (L395-397: Ω_c h²)
  - |w₀|·β = 1.998924 = MIRA ✓ (L842)
  - (π−φ)/Ω² = 0.06725 ✓ (L633: f_screen)
  - Dualidad Ω_m,CMB: MIRA·Ω_dyn=0.319928 vs (π−φ)/(π+φ)=0.320100,
    discrepancia 0.017% (más pequeño que el 0.05% recordado del Roadmap)
- **Ontological**: AGNOSTIC (más cuidadoso que P1, P2)
  - L48, L549: "no free parameter introduced" — no afirma sobre partículas
  - L386 menciona "cold dark matter" pero como **etiqueta observacional**
    de Ω_c h² que el paper **deriva algebraicamente**
  - L688 referencia φ-DM mass relation (Paper 6) — heredando frame de P6
  - **No contiene** "no exotic particles" (a diferencia de P1, P2)
- **MIRA dependency**: STRUCTURAL IDENTITY (no derived mechanism)
  - L211, L439: MIRA en tablas de identidades
  - L842: MIRA = |w₀|·β (identidad algebraica, no derivación mecanística)
- **Status**: **SEALED**
  - Más limpio ontológicamente que P1/P2. No requiere reformulación
    cualquiera sea la decisión global sobre frame DM.
  - Si en sesión futura se decide endurecer P1 a "no exotic particles",
    P4 mantiene su frame agnóstico sin contradicción.
  - Si se suaviza P1 a "SFDM permitido", P4 sigue válido.

### Paper 5 — IS perturbations  *(auditado 2026-05-22)*
- **Numerical**: PASS-PARCIAL + DRIFT-METODOLÓGICO
  - γ_IS = 0.5503 ✓ ≈ γ_ΛCDM
  - MIRA_num = 0.989 ± 0.017 ✓ (NO 1.9989 — confirma que IS perturbations
    no sourcean MIRA)
  - fσ₈ tensión media = 0.74σ ≈ ΛCDM (0.73σ) ✓
  - **DRIFT**: el texto del paper (L102-107) usa 𝒢=1.0109, σ₈=0.821,
    S₈=0.848 (two-sector Poisson, Ω_m=0.32). El guardián usa G=0.866,
    σ₈=0.702, S₈=0.725 (single-sector Poisson, Ω_m=0.16). Ambos
    verifican con la fórmula, son **convenciones distintas**.
  - **En frame suave adoptado**: el paper text es el correcto (matter-mode
    de φ contribuye en Poisson en k<k_fs). El guardián refleja convención
    single-sector vieja.
  - Tensión S₈: paper text da 3.85σ vs KiDS / 3.90σ vs DES; Paper 6
    luego corrige a 2.6σ con WDM transfer function.
- **Ontological**: AGNOSTIC (no contiene "no exotic particles")
  - Usa two-sector φ-DM citando P6 — compatible con frame suave
  - L1069: paragrafo grande sobre two-sector, m_φ=5.60 eV, condicional
    a falsificación en DESI Y3/Euclid
- **MIRA framing**: **EL PAPER QUE DESCARTA EL MECANISMO IS** — referencia
  para todos los otros papers
  - L75-98: Q2 — IS no puede sourcear MIRA, demostrado numéricamente
  - L1206: "MIRA factor does not emerge from classical IS fluid mechanics"
  - L1227: "first-principles QFT derivation remaining open"
  - Honestidad exemplar — debería ser modelo para el §"MIRA Status"
- **Status**: **SEALED-CONCEPTUAL** + **PENDING (resolver drift σ₈/S₈)**
  - Lo sellado: γ_IS, MIRA_num=0.989, fσ₈, framing MIRA-no-IS
  - Pendiente menor: actualizar guardián para reflejar convención
    two-sector del paper (σ₈=0.821, S₈=0.848), o vice-versa
  - **Acción**: agregar entrada V-L4 al Registro: "P5 σ₈/S₈
    convención two-sector vs single-sector"

### Paper 6 — φ-DM dos sectores  *(auditado 2026-05-22)*
- **Numerical**: PASS para identidades algebraicas pero ANSATZ admitido para m_φ
  - Ω_φ-DM = (MIRA−1)·Ω_dyn = 0.160 ✓ (algebraico)
  - σ₈_eff = 0.811 × G_2s = 0.811 × 0.979 = 0.794 ✓ (L109)
  - S₈_eff = 0.820, tensión 2.6σ KiDS ✓ (L109)
  - k_fs = 0.493 h/Mpc — derivación Dodelson-Williams pendiente verificar
  - **m_φ = 5.60 eV — PROBLEMA DIMENSIONAL**:
    - Fórmula L353: m_φ = Σm_ν × H_alg = 0.0824 [eV] × 67.96 [km/s/Mpc]
    - Producto numérico = 5.60, pero unidades [eV·s⁻¹·Mpc⁻¹], NO eV
    - Análisis dimensional propio (ℏ=c=1): √(m_ν · H_alg) ≈ 10⁻¹⁷ eV
      (régimen Fuzzy DM), NO 5.60 eV
    - **Paper 6 admite L371-377**: "numerological ansatz whose physical
      origin... is to be derived in future work"
- **Ontological**: DOS CAMPOS ESCALARES (no "un campo dos fases")
  - L394-395: χ "distinct from the SSEE background k-essence scalar φ
    of Paper 7" — **dos species separadas**
  - L398-404: acción explícita para χ con masa m_φ y coupling ξ·R·χ²
  - **Implicación para frame suave**: la analogía que adoptamos debe ser
    "dos campos escalares cosmológicos (ninguno WIMP/axión)", no
    "un campo dos fases". El resultado defensivo es el mismo —
    sigue siendo no-WIMP, no-axión-QCD, no-SUSY.
- **MIRA dependency**: STRUCTURAL — Paper 6 ES la "physical basis"
  citada por P3, P5, P7 para MIRA two-sector
- **Status**: **NEEDS HONEST REFRAMING**
  - Lo que se sostiene:
    1. Ω_φ-DM = 0.160 algebraico
    2. σ₈_eff = 0.794, S₈_eff = 0.820 (G_2s + WDM transfer)
    3. k_fs = 0.493 h/Mpc como predicción falsificable (DESI Y3/Euclid)
    4. m_φ = 5.60 eV como CHOICE FENOMENOLÓGICA (no derivación)
  - Lo que hay que ser honesto sobre:
    1. m_φ = 5.60 eV es ansatz, no derivación (paper ya lo admite, pero
       el lenguaje "algebraic mass" debería suavizarse a "phenomenological
       mass coincident with Σm_ν × H_0 in mixed units")
    2. La acción L398-404 introduce un segundo campo χ — explícito
  - Editorial work:
    - Reescribir L350-355 para que sea más honesto sobre el carácter
      numerológico del producto m_ν × H_alg
    - Reescribir L394 para alinear con frame suave: "we model the
      dark matter sector as a second cosmological scalar field χ,
      distinct from but of the same fundamental type as φ"
    - Énfasis en k_fs como predicción primaria falsable (ya está bien
      framed en abstract L113-115)

### Paper 7 — EFT canónico  *(auditado 2026-05-22)*
- **Numerical**: PASS
  - βc = −AURA = −(3φ+π)/2 = −3.997847 ✓ (L63)
  - Derivación numérica βc ≈ −3.990 al 0.2% (L348)
  - αK_eff = 3·Ω_DE·(1+w₀) = 0.40320 ✓ (L544, paper dice 0.4033)
- **Ontological**: **MÁXIMO COMPROMISO CON FRAME DM**
  - L148: `L_DM` componente separada del Lagrangiano total
  - L273-285: sección "Covariant φ-DM Coupling" con C(φ) = exp(βc·φ/Mpl)
  - L283: $g^{\rm DM}_{\mu\nu} = C(\phi) g_{\mu\nu}$ métrica conformal a DM
  - L299: $\nabla_\mu T^{\mu\nu}_{\rm DM} = -Q^\nu$ DM con su propio T^μν
  - L305-312: Q ∝ (βc/Mpl)·ρ_DM·φ̇ acopla ρ_DM (no ρ_total)
  - L725: "Ω_m = Ω_m,dyn = 0.160 (**active CDM only**)"
  - L731: ρ_CDM,0 en derivación algebraica de βc
  - **No hay lectura agnóstica posible** — es coupled quintessence
    estándar con dos species (φ + CDM fluid)
- **MIRA dependency**: ESTRUCTURAL via two-sector
  - G_2s = D₁^SSEE(Ω_m=0.320) / D₁^ΛCDM (L473)
  - αK_eff usa Ω_m,dyn=0.160 pero la justificación cosmológica
    requiere two-sector
  - L859: "Ω_m=0.320, Paper 6 — testable with Euclid"
- **Status**: **TIGHTLY COUPLED TO P6**
  - **Hallazgo crítico**: si P6 cae, P7 no es find-replace — el EFT
    está construido sobre la asunción de un fluido DM acoplable
  - La derivación βc=−AURA al 0.2% **requiere ρ_CDM como fluido
    separado** para cerrar el sistema dinámico
  - Si en futura sesión se endurece P1 a "no exotic DM",
    P7 hay que **rederivar desde cero**
  - Posible salvataje: reinterpretar ρ_DM como ρ_b (solo bariones),
    pero entonces βc=−AURA al 0.2% se rompe (Ω_b=0.05, no 0.16)
- **Lo que sí se sella ya**:
  - Numérico (identidades algebraicas)
  - Framing βc = −AURA como cierre del EFT
  - αK_eff = 0.4033 cross-checked con EFTCAMB

### Paper 8 — Strong gravity  *(auditado 2026-05-22)*
- **Numerical**: PASS
  - θ_E^SSEE/θ_E^GR = √βc = √AURA = 1.999462 ✓ (L66)
  - MIRA = (3φ+π)/4 = 1.998924 ✓
  - Diferencia √AURA vs MIRA: 0.0269% (paper L119 dice 0.03%, OK)
  - **Honestidad**: L69 explícito — "near-coincidence, not an algebraic
    identity" — admite el carácter no-derivado de √βc ≈ MIRA
- **Ontological**: ALTO COMPROMISO CON DM (similar a P7)
  - L163: $S_{DM}[\tilde{g}_{\mu\nu}; \psi_{DM}]$ — DM con acción propia
  - L175: "Baryons couple to physical metric; **dark matter couples
    to disformal metric**" — disformal coupling explícito DM↔φ
  - L221: $\nabla^2\phi = (\bc/M_{pl})\rho_{DM}$ — DM sourcea φ
  - L72-76: "disformal coupling selects dark matter as the sole source
    of φ" — solar system tiene ρ_DM≈0 → screening natural
  - **No contiene** "no exotic particles" — agnóstico
  - Pero requiere ψ_DM como fluido físico real, igual que P7
- **MIRA framing**: PARCIAL HONESTO
  - L69: admite coincidencia, no identidad
  - L124: usa MIRA = Ω_m,CMB/Ω_m,dyn (definición operacional)
  - Predicción "factor 2 en θ_E lensing" es observable y falsable
- **Status**: **TIGHTLY COUPLED CON P6/P7** (frame DM requerido)
  - Sello sujeto a:
    1. Verificación observacional de θ_E lensing factor 2 (probably partial
       — k-mouflage screening puede reducir efecto en algunos regímenes)
    2. Resolución del frame DM global (necesita P6 reframe)
  - Lo que se sella ya:
    - Estructura matemática (acción disformal, geodésica nula)
    - √AURA ≈ MIRA al 0.03% como coincidencia honesta
    - K-mouflage screening solar
  - Acción pendiente: verificar literatura SLACS/BELLS sobre θ_E
    factor 2 (¿está observacionalmente refutado o no?)

### Paper 9 — H₀ tension  *(auditado 2026-05-22)*
- **Numerical**: PASS
  - f_scr = αK/(3·MIRA) = 0.06725 ✓ (L130)
  - (π−φ)/Ω² = 0.06725 ✓ (L120, identidad algebraica equivalente)
  - H₀_local = H₀_alg/(1−f_scr) = 67.9621/(1−0.06725) = 72.862 ✓ (L73: 72.86)
  - Tensión vs SH0ES = 0.171σ ✓ (paper claims 0.17σ)
- **Ontological**: COMPROMISO CON DM (similar a P7, P8)
  - L128: "null geodesic of **dark matter**"
  - L157: "conformally coupled to **dark matter**"
  - L524: cita disformal de Paper 8 → ψ_DM con su propia métrica
  - L546: cita φ-DM de Paper 6 + IS de Paper 5
  - **No contiene** "no exotic particles" — agnóstico operacional
- **MIRA framing**: USA P8's disformal MIRA derivation
  - L213-216: MIRA = βc/2 = (3φ+π)/4
  - L218: cita "Paper 8 derives MIRA = Ω_DE^CMB/Ω_DE^dyn"
  - L67: "AURA = 2·MIRA cancel exactly in the ratio" — cancelación
    estructural produce f_scr
- **Status**: **SEALED-NUMERICAL** + **CONDICIONAL en P8 disformal**
  - Lo que se sella:
    1. Identidad algebraica f_scr = αK/(3·MIRA) = (π−φ)/Ω² = 0.06725
    2. H₀_local = 72.86 con tensión SH0ES 0.17σ (resultado clave)
    3. AURA cancellation en el ratio
  - Lo que es condicional:
    1. Validez del frame disformal de Paper 8 (necesario para
       interpretar f_scr físicamente como screening)
    2. Frame DM global (heredado de P6/P7/P8)

### Paper 10 — UV completion  *(auditado 2026-05-22)*
- **Numerical**: PASS
  - 5·φ⁸ = 234.8936 ✓ (paper L204-207: 234.89)
  - 45·(φ⁴/3)² = 234.8936 ✓ identidad algebraica
  - M = φ²·5^(1/4)·(ρ_Λ)^(1/4) = 8.8085 meV ≈ 8.81 ✓
  - **No aparece "8.87"** en ninguna línea — falsa memoria del 8.81
- **Ontological**: NEUTRAL (no menciona DM exotic/WIMP/axion en cuerpo;
  L1022 solo cita Amendola 2000 en bibliografía)
- **MIRA dependency**: EMPIRICAL INPUT
  - f_sc = α_K/(3·MIRA) usa MIRA estructuralmente (L103-105, L347, L423, L455)
  - L815 explícito: "M only modifies α_K; background MIRA unchanged" —
    MIRA es input no derivado de la UV completion
- **Status**: **SEALED-NUMERICAL** + **NEEDS-MICRO-EDIT** (notación)
- **Micro-edit pendiente**: L192 dice "ρ_c^(1/4) ≈ 2.25 meV is the
  cosmological constant scale today". Debe ser **ρ_Λ** o **ρ_DE**, no
  ρ_c (densidad crítica total tiene 2.4 meV, no 2.25). Resultado M=8.81
  no cambia; es claridad notacional.
