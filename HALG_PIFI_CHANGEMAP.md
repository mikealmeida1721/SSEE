# Mapa de Cambios — Reframe MIRA→π/φ (factor materia) + H_alg global

**Estado:** EN APLICACIÓN 1-por-1 (Mike adoptó 2026-06-17).
- ✅ **Cambio #1 APLICADO — capa código + 3 memorias** (Ω_m,CMB π/φ=0.31076 + H global=H_alg=67.962):
  · código: `ssee_core.py` (PI_OVER_PHI, OMEGA_M_CMB, H0_GLOBAL; MIRA/AURA entidades intactas)
  · guardián: `ssee_verify.py` (V-L2-05 canónico + V-L2-05b identidad MIRA histórica)
  · fuente única: `CANONICAL_VALUES.yaml` (Omega_m_CMB, pi_over_phi_factor, H_alg titular)
  · CLAUDE.md: tabla constantes + nota dualidad
  · Obsidian nodos-definición: `Constantes/MIRA.md` (rol reasignado a π/φ), `Constantes/H0.md` (H global titular)
  Guardián 114 VERDE · memory_sync VERDE. Estado mixto aceptado: partición two-sector
  (MIRA·dyn) sigue como identidad verdadera; nodos vault de RESULTADO (Papers/Resultados)
  recomputan en su cajón. NO retirado aún el 0.31993 (beacon se activa al propagar cajones).

- ✅ **FASE A completa — los 3 inputs FIJOS en las 3 memorias** (álgebra pura, verde):
  · H global = H_alg = 67.962
  · Ω_m,CMB = (π/φ)·0.160 = 0.31076
  · m_φ = Σm_ν·(PYROS·VITA·MIKA) = 0.06902·615.33 = **42.47 eV** (mult 535.28→615.33)
  Propagado a: ssee_core, guardián (V-L2-10/V-L3-mphi → 42.47), CANONICAL_VALUES.yaml,
  CLAUDE.md, vault (MIRA.md/H0.md/kfs.md), memoria p6_canonical_particle + MEMORY.md.
  Guardián 115 VERDE · memory_sync VERDE.

  **DEPENDIENTES marcados PENDIENTE Fase B** (NO se escriben sin correr código; pueden dar ROJO):
  k_fs, alpha, σ8/S8 single+two-sector, fσ8, cascada Hubble IR/UV, r_d, posterior MCMC, χ²_CMB.
  Flag explícito en guardián: track_open "REFRAME-FaseB".

  FASE B (cajón scripts) — por dueño del script:
  - ✅ **P1 ingredientes** (H titular, Ω_m→π/φ, m_φ→42.47, Postulado M, Two-Ω_m, dual retirada) — compila 30pp
  - ✅ **Cascada Hubble (P9/P10)** VERIFICADO corriendo scripts: IR=72.86 (0.17σ), UV=73.040 (0.00σ)
       con H global=67.962. Flip retired↔canonical en yaml (71.87/72.05 ahora viejos). Propagado:
       CLAUDE.md, vault H0.md + 7 nodos espejo. Guardián 115 VERDE, memory_sync VERDE.
  - ✅ **P9/P10 punta a punta**: manuscritos .tex (compilan 19/14pp, cascada invertida a H_alg global,
       Apéndice P reframeado como confirmación del ancla) + figuras regeneradas (72.86/73.040) + scripts figura.
  - ⏳ PENDIENTE (regla: por script, propagar COMPLETO antes del siguiente): P3 CMB (χ²/ΔBIC/r_d con
       Ω_m=0.31076 — posible rojo leve) · P6 CLASS (k_fs/alpha/σ8/S8 con m_φ=42.47) · P5 fσ8 · P2 MCMC ·
       luego P2/P4/P7/P8 + Unified/Endorser/Sealed → limpieza/archive → re-auditoría raíz.
**Origen:** investigación 2026-06-17, memoria `project_halg_pifi_investigation`.

## RAÍZ (un solo cambio conceptual)
MIRA **deja de ser el factor-materia**; ese rol pasa a **π/φ=1.9416**.
MIRA **persiste** como AURA/2=1.9989 (sigue en f_screen, sigue cumpliendo AURA=2·MIRA).

```
Ω_m,CMB = factor × Ω_m,dyn
   VIEJO: MIRA(1.9989) × 0.160 = 0.31993
   NUEVO: (π/φ)(1.9416) × 0.160 = 0.31069
H-global (CMB-fit) sigue al factor:  H_MIRA 67.037 → H_alg 67.962
```

## CADENA — qué se mueve (Tier 1) — VERIFICADO vs PENDIENTE
| Valor | Viejo | Nuevo | Estado |
|---|---|---|---|
| Ω_m,CMB | 0.31993 | 0.3107 | álgebra |
| H global | 67.037 | 67.962 | VERIF (CMB-fit) |
| CMB χ² | ~1019 | 1056 | VERIF cobaya |
| σ8 single | 0.820 | 0.821 | VERIF CAMB |
| S8 single | 0.847 | 0.836 | VERIF CAMB |
| σ8 two-sector | 0.742 | 0.731* | *con m_φ=36.95 FIJO — NO es el target |
| S8 two-sector | 0.766 | 0.744* | *PARTÍCULA SUELTA: 36.95 da 1.11σ; escenario nuevo PIDE m_φ≈42.5 eV para S8=0.766 (0σ). NO resuelto — abierto. ¿SSEE da ~42.5? (forward=36.95, ratio 1.15) |
| H_local IR | 71.87 | 72.86 (0.17σ) | cascada |
| r_d | ~147 | ~invariante | VERIF CAMB |
| H_local UV (P10) | 72.05 | ~73.0 | PENDIENTE escala |
| G_growth | 1.0126 | 1.0078 | VERIF (≈igual, ambos ~1%) |
| fσ8 single (delta) | ref | +0.04σ | VERIF (≈igual; absoluto pipeline-roto, two-sector pendiente) |
| posterior MCMC | 66.53 | ~67.4 estimado | PENDIENTE corrida real |

## VEREDICTO (2026-06-17, verificado): NUEVO 0.311 vs VIEJO 0.320
- 🟢 MEJOR: Hubble (0.17σ vs 1.12σ), Parsimonia (2 vs 3 incógnitas)
- 🟡 leve peor: CMB χ² (1056 vs 1019), S8 two-sec (1.11σ vs 0.42σ, AMBOS resuelven)
- ⚪ IGUAL: r_d, σ8 single, G_growth, fσ8, m_φ, f_screen
- Net: intercambio favorable. NO rompe nada. Falta solo MCMC real.

## NO SE TOCA (invariantes — Tier 2/3)
- **m_φ = 36.95 eV** (usa AURA, no Ω_m) — INTACTO
- **f_screen = 0.06725** (usa MIRA-valor, retenido) — INTACTO
- **MIRA = 1.9989** (sigue existiendo como AURA/2) — INTACTO
- AURA, w0=−0.840, wa, Ω_m,dyn=0.160, H_alg, αK=0.4033, KAL — INTACTOS

## ARCHIVOS afectados (orden de aplicación, SOLO si se adopta)
1. `src/ssee_core.py` — OMEGA_M_CMB: nuevo factor (CUIDADO: NO tocar MIRA ni AURA)
2. `CANONICAL_VALUES.yaml` — Ω_m, H, σ8/S8, H_local (fuente única)
3. Guardián `src/verificacion/ssee_verify.py` — recalibrar checks
4. `CLAUDE.md` — tabla valores canónicos
5. Memorias — sync via `memory_sync.py`
6. Papers: P2(MCMC), P3(CMB), P5(IS/fσ8), P6(two-sector), P9(Hubble), P10(UV)
7. Unified, Sealed (PUERTA UN SOLO SENTIDO — decisión aparte)

## COMPUERTAS antes de aplicar a papers
1. fσ8 sanity (no rompa) — single-sector, sin partícula
2. MCMC posterior
3. Decisión PI: ¿adoptar reframe sin mecanismo de π/φ? (OP-8 sigue abierto)
