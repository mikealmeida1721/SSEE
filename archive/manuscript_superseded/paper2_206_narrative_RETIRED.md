# Paper 2 — narrativa "+206 / r_d=175.6 / two-Ω_m" RETIRADA (2026-07-09)

**Motivo:** V-L4-DESI — el bug de geometría. La densidad de materia de fondo
(E(z), r_d, distancias BAO) debe ser la materia **total** gravitante
Ω_m = ω_m/h² = 0.30889, no el sector frío 1+w₀ = 0.160 (que es la ecuación de
estado, no una densidad). Meter 0.160 en E(z) producía SEIS tensiones falsas
que Paper 2 documentaba honestamente como limitaciones reales:

| Observable | Con bug (0.160 en E(z)) | Corregido (0.30889) |
|---|---|---|
| χ²_BAO | 726 | ~11 |
| Ω_m tensión | 21.3σ | 0.88σ |
| r_d | 175.6 Mpc (1.19×) | 148.2 (1.002×) |
| θ* | 6.66σ | 0.91σ |
| H₀ | 66.41 (2.9σ bajo anchor) | 67.95 (0.04σ del anchor) |
| χ²_r(H(z)) | 1.86 (3σ) | 0.482 (≈ΛCDM 0.459) |

## Narrativa retirada (git preserva el detalle; commit anterior a a942c9f+)

1. **"Two ΔBIC domains" (+206 y −5.55):** el +206 era el "full-model penalty
   under standard Friedmann" al forzar SSEE con Ω_m=0.160 contra ΛCDM. Con la
   total, ese escenario no existe. ΔBIC canónico = **+5.68 (ΛCDM), +5.59 (CPL)
   favoreciendo SSEE** (3-modelos, mismo prior, tabla limpia).

2. **Tabla vieja (L894-913):** SSEE-V3.6 k=2 BIC 242.89 ΔBIC +206.19; SSEE_dyn
   k=1 BIC 31.15 ΔBIC −5.55; ΛCDM ref; CPL +0.51. → Reemplazada por la tabla
   del 3-modelos (SSEE ref, ΛCDM +5.68, CPL +5.59).

3. **Apéndice k-scheme (+206/+216 "k-scheme-independent finding"):** defendía
   el +206 bajo tres conteos de k. Retirado — no hay +206 que defender.

4. **Sección H(z) "most direct limitation, needs new gravitational physics
   (non-minimal coupling / Kinetic Braiding)":** era el 0.160 en H(z). Con la
   total, χ²_r(H(z))=0.482, indistinguible de ΛCDM. Sin necesidad de extensión.

5. **Figura "dual sound horizon" (175.6 vs 147.17):** no hay dual; r_d físico
   147.17 (CAMB) / 148.2 (BAO run) ≈ Planck.

**Fuente de los valores nuevos:** results/logs/mcmc_paper2_3models_om308.log,
mcmc_paper2_reframe.log, ssee_phase_c_dic.py, ssee_phase_d_savage_cv.py,
run_p3_rd_reframe.py. Guardián R14 canario-geometría vigila la no-recurrencia.
