# Veredicto de Nova — push de Claude (`37175ad`, `866c1be`)

Fecha: 2026-09-25 noche. Rama: `fix/etiquetas-h0-ir`.
Protocolo §0 del handoff (decisión de Mike): verifico el push de Claude en el remoto.

## Qué revisé

- Diff completo `32f4fa9..866c1be` en `CANONICAL_VALUES.yaml` y los 8 .tex.
- Aritmética independiente (no lectura) de los números que Claude declara
  "comprobados por aritmética".

## Verificación

- **A5/A6 signo ΔBIC/AIC:** el YAML tenía las dos convenciones mezcladas, como
  dice el mensaje. Con Δχ²=+2.65, Δk=4, N=225: AIC = 2.65−8 = **−5.35** ✓;
  BIC = 2.65−4·ln(225) = 2.65−21.664 = **−19.01** ✓. Claves corregidas a la
  convención Δ≡SSEE−ΛCDM declarada una vez. Propagación a P2/P6/Sealed/Unified
  con signo consistente (muestreo: Sealed/Unified −19.0 "favouring SSEE" ✓).
- **A6 `model_comp`:** tabla reconstruida de una sola corrida; reconcilia:
  22.92−16.49=6.43 ✓, DIC 20.61−14.95=5.66 ✓, ΔDIC vs CPL −4.02 ✓. El caption
  documenta la corrida vieja y su valor. ✓
- **H_global único (`37175ad`):** `H0_glob_km_s_Mpc=67.96214` canónico;
  68.13 degradado a `H0_glob_regimen_IR_historico` (no canónico);
  `num_alg_3Omega2` documentado como blanco sin unidades (nunca entrada).
  f_screen = 0.41691/(3×1.99892) = 0.069522 ✓;
  73.04×(1−0.069522) = 67.9621 ✓. Paper 9 usa `H_0^{glob,IR}` en los 9 sitios
  de régimen. ✓
- **A7–A9, A15:** 68.13 → `H_0^{glob,IR}` (muestreo P9 ✓).
- **A10:** P3 reescrito ("lo que s_DE/s_m NO son") + labels `eq:sDE`/`eq:sm`. ✓
- **A11/A12:** confirmados cerrados por el commit `25c0b89` de Nova (Paper 3:
  abstract l.44 y l.957–959 ya dicen `s_m`). ✓
- **A18:** confirmada la sospecha; P1 ahora cita la sección por nombre
  ("Posterior predictive check: H(z)") en vez del §5.3 hardcodeado. ✓
- **A19:** 0.160050σ → 0.157σ; (0.965558−0.9649)/0.0042 = 0.1567 ✓.
- **No tocado:** A1–A3 (δc) y A4 (edad P9) — intactos en `manuscript/`, como
  exige el handoff. ✓
- **Control cruzado de edad:** Claude reprodujo la nota de Nova por su cuenta:
  13.733 Gyr (vs 13.73) y ΛCDM 13.796 (vs 13.797 publicado) — el control cierra
  a 0.01%. ✓
- 9 PDFs recompilados (binarios actualizados en el commit).

## Observaciones (no bloquean)

- Los commits son de las 21:19/21:27, anteriores al protocolo §0 (21:51): por eso
  no traen archivo de veredicto `notes/..._verificacion_claude.md`; los mensajes
  de commit hacen esa función con creces.
- P1 suma una nota de retiro de β_c=−AURA (fuera de la lista A1–A19): es
  documentación de un retiro, no cambio de física; queda registrado aquí.

## Veredicto: VERDE

13 de 19 cerrados por Claude + 2 (A11/A12) ya cerrados por Nova = 15/19.
Faltan A1–A4: requieren el "adelante" de Mike (parche δc y reescritura edad P9).
Nada pendiente de Claude en este ciclo.

---

## Adenda — `3e929f4` "Verificación de los 4 que faltan" (verificado por Nova)

- **Cuenta corregida:** 4 faltantes (A1–A4), no 6. ✓ (A11/A12 ya eran míos.)
- **A1–A3 δc — reproducción independiente:** sus números son idénticos a los de
  mi Ruta 1: z_c=0 → EdS 1.68646, ΛCDM 1.67599, SSEE 1.67634; z_c=10 →
  1.68647/1.68646/1.68647. El control EdS pasa contra el analítico. ✓
- **Corrección a mi script:** la justificación "DE suave: consistente con
  c²_s,eff=0" estaba invertida (c_s²=0 anula el horizonte sonoro, lo que
  permitiría agrupamiento). Corrigió el comentario a la fricción viscosa IS de
  Paper 5; los números no cambian. Aceptado.
- **Refinamiento JWST:** con δc derivado el efecto no es ≈1.00× sino <1
  (0.894 a z=10, 3×10¹² M☉; 0.778 a z=15): SSEE predice MENOS halos masivos
  tempranos que ΛCDM (causa: Ω_m menor + fondo CPL → menos crecimiento a z
  alto). Físicamente plausible; refina mi nota.
- **Radio de δc más ancho — CONFIRMADO por grep:** el postulado vive en 7
  sitios, no 3: Paper 4 (l.487, l.772), Paper 5 (l.1411, l.1441), README.md
  (l.314, la portada), AUDIT.md, CHANGELOG.md y el código vivo
  `src/p02_mcmc/ssee_press_schechter.py` (mi fix σ8 de `312eb91` dejó el δc
  falsificado dentro — catch válido de Claude). Todo requiere el "adelante".
- **A4 tres defectos apilados — CONFIRMADO en el texto** (P9 l.943–971):
  (1) Ω_m=Ω_{m,dyn}=0.160050 en E(z); (2) H_0^SH0ES=73.04 local para cantidad
  global; (3) cita como autoridad viva el Two-Ω_m Criterion que Paper 1 declara
  retirado. Más el cuarto: el criterio precomprometido "t_0>14 Gyr would favour
  SSEE" construido sobre el bug — con 13.733 el modelo se reprobaría a sí mismo.
  Radio CONTENIDO en Paper 9 (no llegó a P1/Sealed/Unified/Endorser). ✓
- **No editó `manuscript/`** en A1–A4: verificación, no propagación. ✓

Veredicto adenda: VERDE. Estado final: 15/19 cerrados; A1–A4 verificados y
acotados, esperando el "adelante" de Mike.
