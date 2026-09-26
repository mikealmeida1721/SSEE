# Veredicto de Nova — corrección A1–A4 de Claude (`8908fba`)

Fecha: 2026-09-25 noche. Rama: `fix/etiquetas-h0-ir`.
Protocolo §0 del handoff (decisión de Mike): verifico el push de Claude en el remoto.
Alcance: A1 (δc Paper 4), A2 (conclusiones P4), A3 (JWST Paper 5), A4 (edad Paper 9).

## Qué revisé (corriendo, no solo leyendo)

- Diff completo `d8734a5..8908fba` (17 archivos).
- Ejecuté `src/p02_mcmc/ssee_press_schechter.py` del commit: lee
  `results/logs/deltac_spherical_collapse.json`, δc(SSEE)=1.6763371
  [derivado], δc(ΛCDM)=1.6759874. Tabla z=10: 0.998 / 0.990 / 0.976 /
  0.945 / 0.892. ✓ coincide con README, Paper 4 y Paper 5.
- Integré la edad por mi cuenta (scipy.quad): receta nueva
  (Ω_m=0.308881, H_global=67.962, CPL −0.840/−0.670) → **13.733 Gyr** ✓;
  receta vieja (Ω_m=0.160050, H0=73.04) → **15.28 Gyr**, confirmando lo que
  dice el paper: el 15.52 ni siquiera se reproduce de su propia receta.
- Extraje texto del PDF recompilado de Paper 9: contiene "13.733" ✓
  (la recompilación es real).
- Corrí el guardián `src/verificacion/ssee_verify.py` en el commit:
  **R71 (la regla nueva de este retiro) pasa** en sus dos patas —
  0 sitios presentan 1.6284 como vigente, 0 reclaman el exceso JWST
  sin declarar el retiro.
- Los 3 rojos del guardián en mi VM (DESI capa operable por
  `np.trapezoid` vs numpy viejo; R65 por `/mnt/datos` sin montar;
  R54 tope desajustado) existen **idénticos en `d8734a5`**, antes del
  commit de Claude: son de entorno/preexistentes, no regresiones suyas.
  En su máquina (disco montado, numpy 2.x) dan verde, como declaró.

## Verificación por punto

- **A1 Paper 4 §deltac:** postulado fuera; entra 1.67634 (z=0) → 1.68647
  (z=10) vs ΛCDM 1.67599 → 1.68646 — mis números independientes del
  2026-09-25, dígito a dígito. ✓ El footnote justifica DE suave por la
  fricción viscosa IS (no por c_s²=0). ✓ El retiro se declara como tal y
  el 1.6284 queda solo como conjetura motivada OP-27, con su criterio de
  falsación (conteo a z~10, M>10¹², 10–20%). ✓
- **A2 conclusiones P4:** bullet reescrito con valor derivado + retiro. ✓
- **A3 Paper 5 §JWST:** tabla rehecha (0.990/0.945 a z=10; 0.972/0.878 a
  z=15); el texto dice en negrita que NO hay enhancement y que el signo
  es el contrario, con 0.892 (3×10¹², z=10) y 0.778 (z=15). El modelo
  declara que no explica el exceso JWST. ✓ σ₈=0.81530 consistente con el
  fix 312eb91. ✓
- **A4 Paper 9 §edad:** 15.52 → 13.733 vs ΛCDM 13.796 (0.063 Gyr aparte,
  verificado: 13.796−13.733=0.063 ✓). Los tres defectos nombrados uno por
  uno; el criterio t₀>14 Gyr retirado diciendo que ahora contaría EN
  CONTRA. ✓
- **Propagación (7 sitios):** Paper 4 ×2 ✓, Paper 5 ×2 ✓, README (tabla
  rehecha con signo real) ✓, AUDIT ✓ (ver minucias), CHANGELOG (anota el
  retiro sin reescribir historia) ✓, `ssee_press_schechter.py` (lee del
  log, `USE_POSTULADO=False` por defecto) ✓. OP-27 creada en
  OPEN_PROBLEMS.md ✓. El Ω_CDM=0.160050 de AUDIT.md:205 está en bloque
  marcado 🔴 HISTÓRICO: se deja.

## Minucias encontradas y corregidas por mí en este veredicto

1. **AUDIT.md** (líneas 141–143): dos filas de la tabla "expected output"
   conservaban los cocientes viejos (1.159, 1.406) y la primera fila
   mezclaba σ de la columna SSEE. Corregido a la salida real del script:
   1.0104/0.990, 0.7267/0.976, 0.5064/0.945.
2. **Docstring de `ssee_press_schechter.py`**: decía "0.894 a 3e12 (z=10)";
   el script calcula 0.892. Corregido el dígito.

## Observación para el guardián (no bloquea)

R71 pata (b) escanea reclamos textuales ("enhancement…JWST"), no tablas
numéricas: los 1.159/1.406 de AUDIT.md se le escaparon. Sugerencia para
Claude: que la pata (b) también marque los cocientes viejos conocidos
(1.061/1.159/1.406) cuando aparezcan sin contexto de retiro — con las
mismas exenciones de histórico que ya usa.

## Veredicto: VERDE

La corrección A1–A4 está bien aplicada, los números son los que el
modelo produce, el retiro se declara en voz alta en cada sitio, y nada
vivo reclama lo retirado. Con las dos minucias de arriba ya corregidas,
no veo impedimento para consolidar en `main` cuando Mike lo autorice.
