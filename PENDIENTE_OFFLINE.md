# Corridas pendientes — se pueden lanzar SIN internet (2026-07-02)

Todo ya está wireado al vector DESI DR2 corregido (V-L4-DESI). Cada corrida es
un solo comando desde la raíz del repo. Los resultados quedan en disco para
analizarlos con Claude cuando vuelva el internet.

## 1. MCMC FULL con DR2 (la grande, ~13–17 h, dejar la PC prendida una noche)
```bash
bash src/mcmc_full/lanzar_mcmc_full.sh
```
- Lanza 4 cadenas Cobaya independientes (12 cores). Sobreviven al cierre de la
  terminal (`setsid`). Output → `/mnt/datos/SSEE_data/mcmc/mcmc_full/` (HDD).
- Monitoreo: `tail -f /mnt/datos/SSEE_data/mcmc/mcmc_full/logs/cadena_1.log`
- Al terminar, con Claude: combinar con getdist, actualizar tabla Unified §8
  (marcada con `% TODO-DR2` en el tex — hoy usa las chains DR1 viejas).

## 2. MCMC 3-modelos Paper 2 (SSEE vs ΛCDM vs CPL, ~1.5 h)
```bash
python3 src/p02_mcmc/ssee_paper2_mcmc.py > results/logs/mcmc_paper2_3models_dr2.log 2>&1 &
```
- Regenera el ΔBIC (+7.91 era DR1) y alimenta fig7_Hz/fig8_tension.

## 3. ΛCDM baseline (~30 min) — necesario para la comparación H₀ de Paper 2
```bash
python3 src/p02_mcmc/ssee_paper2_mcmc_lcdm_baseline.py > results/logs/lcdm_baseline_dr2.log 2>&1 &
```

## 4. Secundarios (minutos–1h cada uno, en cualquier orden)
```bash
python3 src/p02_mcmc/rerun_cpl_h0anchors.py        > results/logs/cpl_h0anchors_dr2.log 2>&1
python3 src/p09_hubble/ssee_h0_prior_experiment.py  > results/logs/h0_prior_exp_dr2.log 2>&1
python3 src/estadistica/ssee_phase_d_savage_cv.py   > results/logs/savage_cv_dr2.log 2>&1
```

## Verificación después de cualquier corrida
```bash
python3 src/verificacion/ssee_verify.py   # debe decir VERDE
```

## Estado al cierre de hoy (commit b6dab36)
- Corregido y committeado: fuente única DR2 + 8 scripts + guardián R14 +
  titulares (README/AUDIT/Paper2/Unified/CANONICAL/Ledger) + fig1 nueva.
- Ya re-corridos HOY con DR2: análisis w₀wₐ (0.16σ Pantheon+ / 1.2–1.5σ resto),
  MCMC P2 reframe (H₀=66.41±0.39), Fase 4 flat ciego (media 0.72σ, máx 1.23σ).
- PENDIENTE con Claude (necesita internet): propagar a los papers restantes
  (P1/P3/P4/P5/P6 menciones DESI), CLAUDE.md tablas, cover letters, lectura
  prosa 10 papers, tabla Unified §8 con el mcmc_full DR2 nuevo, CITATION.cff.
