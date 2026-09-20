# Cajón: lo que KiDS-Legacy dejó atrás (2026-09-20)

Aquí vive código de la época de **KiDS-1000** cuyo motivo desapareció cuando
KiDS-Legacy (Wright et al. 2025, A&A 703 A158, arXiv:2503.19441) recalibró las
distribuciones de redshift y la tensión Planck–cizalla se disolvió en el dato.

No está mal: está **superado**. Se conserva porque es la procedencia de
resultados que sí se publicaron.

## `s8_barra_kids.py`

**Qué hacía.** Perfil de χ² sobre A_s para dar σ(S₈) del lado KiDS bajo SSEE,
homólogo del ±0.024 que KiDS publicaba bajo ΛCDM. Servía para expresar la
tensión «en la receta de la literatura» sin suponer barras.

**Por qué se archiva, tres razones independientes:**

1. **Ya no hay tensión que expresar.** Con KiDS-Legacy el hueco mide 0.49σ.
2. **No corría.** Importaba desde `/tmp/claude-1000/.../238cf748-.../scratchpad`,
   una ruta de sesión efímera que ya no existe.
3. **Sus números eran de julio**, anteriores al reframe: usaba S₈ del CMB =
   0.8264 cuando con el fondo de hoy son 0.8273.

**Lo que destapó al archivarlo, que es la lección que se queda:** tenía
`S8_cmb = 0.82639` tecleado a mano. Ese valor **no corresponde a ningún log**:
`precio_cmb_de_la_particula.json` da 0.8264027 y
`growth_2026-07/base_sin_particula.json` da 0.8263335. Es un error de tecleo de
0.82640 — una diferencia de 1e-5 que nadie habría visto nunca, porque **un
número escrito a mano no deja registro de dónde salió**. Se corrigió al valor
con origen antes de archivar, para que la historia quede bien.

Regla que deja: ver `memory/feedback_no_hand_typed_numbers.md`.
