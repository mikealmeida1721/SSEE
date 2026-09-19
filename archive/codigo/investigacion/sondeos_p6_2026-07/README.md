# Sondeos de Paper 6 (julio 2026) — archivados

## `lcdm_kids.py` (archivado 2026-09-19)

Sondeo de la época en que renació Paper 6 (commit `fc3b594`): ajuste de KiDS-1000 con
un fondo ΛCDM-Planck **clavado**, minimizando sólo la amplitud y las molestias.

**Por qué se archiva y no se arregla:**
- **No puede correr.** Importa `kids_shear` desde el scratchpad de una sesión de
  julio (`/tmp/claude-1000/.../238cf748-.../scratchpad`), que ya no existe.
- **Nadie lo usa.** Ningún script lo importa y ningún log vigente sale de él. La
  vara ΛCDM canónica es R4 (`cobaya_kids.py lcdm`, S₈ = 0.7571 ± 0.0194), y la de
  fondo fijo, `cobaya_kids.py lcdmfijo`.
- **Sus cuatro números de fondo no tienen fuente identificada.** El comentario
  dice «media ponderada de las 4 cadenas» (ω_b = 0.02238396, ω_c = 0.11995244,
  h = 0.6737330790, n_s = 0.96530732) sin decir qué cadenas. Lo detectó R65 al
  rastrear el origen de cada número; no se le inventa uno.

Se conserva como registro de lo que se intentó.
