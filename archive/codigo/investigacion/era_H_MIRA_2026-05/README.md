# Era H_MIRA (mayo 2026) — scripts archivados el 2026-09-19

Archivados por decisión de Mike («si ya no los usa, muévelos») tras el
rastreo de orígenes R65. **Ningún log vigente, paper ni número canónico sale
de ellos**: se comprobó con grep sobre `src/`, `manuscript/`, los `.md` y
`PROPAGACION.yaml` antes de moverlos.

## `h0_cascade_audit.py` (antes en `src/p09_hubble/`)

Auditoría bloque por bloque del 2026-05-24: cada cantidad que depende de H₀,
evaluada bajo H_alg = 67.962 y H_MIRA = 67.068. **Superada dos veces**: el
factor MIRA salió de la cascada con el reframe ω_m-directo (2026-06-18), y la
dirección de la cascada se invirtió el 2026-09-06 (SH0ES entra, H_global sale).
Contiene además `rho_DE = OMEGA_DE * rho_crit`, la saturación metida como
densidad (el caso que motivó una regla del guardián). Su fila de Σm_ν usa
0.0824 eV, valor retirado (contaminado OP-14; el vigente es 0.06849).
Los orígenes de sus números quedaron anotados (`# ORIGEN-VALOR`) antes de
moverlo.

## `ssee_verify_rd.py` (antes en `src/verificacion/`)

Verificación CAMB de r_d y θ* contra Planck (tarea 2A, 2026-05). Evalúa
configuraciones **retiradas**: SSEE+MIRA (Ωc h² = 0.1214) y el sector
Ω_m = 0.160 (el bug del χ²=726). El r_d/θ* canónico sale de
`src/p03_cmb/run_p3_rd_reframe.py` → `results/logs/p3_rd_reframe_omega_m.log`
(r_d 0.32σ, θ* 1.00σ).

**Bug hallado al archivar (R65, 2026-09-19) y corregido antes de mover:**
`THETA_OBS = 0.59668°` equivalía a 100θ* = 1.04140, no al 1.04110 de Planck
2018 que decía su comentario, y `THETA_SIGMA = 0.00046°` inflaba el σ ×2.6
(el de Planck es 0.00031 en 100θ*, 0.000178°). Ahora los grados se calculan
desde 100θ*. Efecto en lo que imprime, sólo en configuraciones retiradas:
θ* caso A 0.97σ → 1.53σ; caso B 1.18σ → 2.09σ. Nadie citaba esos números.
