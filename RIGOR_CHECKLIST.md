# SSEE — Lista de Rigurosidad (Estándar de Auditoría)

> Cada regla de esta lista nació de algo que **se pasó por alto** en auditorías
> previas. Es el estándar de revisión obligatorio antes de enviar cualquier
> manuscrito. Cuando se detecte un nuevo tipo de fallo, se **marca aquí** para que
> el sistema de auditoría no vuelva a dejarlo pasar.
>
> Leyenda: 🔴 letal (hunde el paper) · 🟠 grave · 🟡 de presentación
> Estado: ⛔ se pasó por alto en el pasado · ✅ regla activa

---

## R1 — Prioridad temporal: verificar fechas SIEMPRE 🔴 ⛔
**Regla:** nunca afirmar "prior to / predating / before / committed before [dataset]"
sin verificar que la fecha del registro (DOI, commit, Zenodo) es **objetivamente
anterior** a la fecha de publicación del dataset.
**Por qué (caso real):** los papers afirmaban que MIRA/w₀ fueron "fixed prior to
DESI DR2", citando Zenodo 19679049 (**2026-01-28**). DESI DR2 es **2025-03-19**.
El registro es ~10 meses POSTERIOR. Una línea decía literalmente "2026 January 28,
prior to the DESI DR2 release on 2025 March 19" — cronológicamente imposible.
**Verificación:** `grep -rniE "predat|prior to|before (desi|the release)|committed before|timestamp"` en todos los `.tex`; cruzar cada fecha contra la del dato.
**Marco correcto:** el match con datos ya públicos es **postdicción libre de
parámetros**; la fortaleza es la **rigidez estructural**, no la cronología. La
prioridad temporal genuina solo se reclama sobre datos NO publicados (DESI DR3, Euclid).

## R2 — Enlaces a material multidominio / numerológico 🔴 ⛔
**Regla:** ningún manuscrito, README, ni Zenodo citado puede enlazar a material que
aplique el sistema a dominios no-cosmológicos (medicina, química, paz, ADN…) ni a
versiones numerológicas crudas. Activa el detector "crackpot" del referee.
**Por qué (caso real):** el repo público `SSEE_UNIFICADO` (Génesis 5.12) aplica φ,π
a 8 dominios y tiene H₀ crudo ("73.483 = H_global + KAL") que **contradice** la
versión sofisticada de los papers (f_screen algebraico).
**Verificación:** seguir cada enlace/DOI de los papers hasta su destino final;
confirmar que solo expone estructura algebraica + física cosmológica.
**Resuelto (2026-06-14):** creado registro limpio dedicado
`SSEE-Constant-Dictionary` (Zenodo concept DOI `10.5281/zenodo.20684908`), solo
diccionario + predicciones, sin material multidominio. Todas las citas del
diccionario (P1, P2, P3, P4, Unified, Sealed + bibs + AUDIT.md) reapuntadas del
viejo `19679049` (SSEE_UNIFICADO multidominio) al DOI limpio. Disponibilidad de
**código** sigue apuntando al archivo legítimo `20093447` (serie de papers + CLASS).

## R3 — Look-elsewhere sobre el diccionario COMPLETO 🟠 ⛔
**Regla:** el conteo look-elsewhere se hace sobre TODAS las constantes del diccionario
fuente, nunca un subconjunto. Reportar la **curva de sensibilidad a la tolerancia** y
la **robustez a extensiones futuras**.
**Por qué (caso real):** el conteo usaba 21 de 31 constantes → acusable de "subset a
conveniencia". Corregido a 31 (317 razones, 1/317 a ±0.001); verificado robusto a las
copias QUINTAL–DECAL. Script: `src/estadistica/look_elsewhere_full.py`.

## R4 — Likelihoods reales, no comprimidas (modelos no-ΛCDM) 🟠 ⛔
**Regla:** no presentar distance priors comprimidos de Planck como equivalentes a la
likelihood completa de los Cℓ, especialmente en un modelo no-ΛCDM (los distance priors
se calibran asumiendo ΛCDM → posible sesgo).
**Por qué (caso real):** el MCMC Fase 4 usaba prior comprimido; un referee objeta el
sesgo. Mitigación: MCMC full con `plik_lite_native` + lensing + DESI + fσ₈ (en curso).

## R5 — Benchmark = ancla canónica, no conjetura obsoleta 🟠 ⛔
**Regla:** comparar posteriores contra el valor canónico vigente, no contra valores
heredados de iteraciones previas.
**Por qué (caso real):** Fase 4 comparaba H₀ posterior contra la conjetura Type-P
67.96 (tensión 0.81σ) en vez del ancla canónica H_MIRA=67.037 (0.39σ), porque se corrió
antes de canonizar H_MIRA. Drift no detectado.

## R6 — Conteo y valores exactos contra la fuente 🟡 ⛔
**Regla:** todo número citado (nº de constantes, razones, valores) debe regenerarse
desde el script/fuente, no copiarse de iteraciones previas.
**Por qué (caso real):** "21 constantes" / "29" cuando el diccionario fuente tiene 31;
valores con ruido en el último decimal copiados de un compendio viejo.

## R7 — Jerga inventada / nombres mitológicos 🟡
**Regla:** la cara del paper usa notación algebraica; los nombres de linaje (IGNIS,
SOLAR, OSIRIS…) van solo en un apéndice de linaje, nunca en ecuaciones del cuerpo.
**Por qué:** activan el detector crackpot. Pendiente de limpieza transversal en el Sealed.

## R8 — Drift de estado (notas/memorias/papers) 🟡 ⛔
**Regla:** antes de enviar, verificar que notas, memorias y papers reflejan el estado
REAL del trabajo (no un estado congelado).
**Por qué (caso real):** el vault marcaba "§8 pendiente / MCMC corriendo" cuando ya
estaba escrito, commiteado y compilado desde mayo.

## R9 — Valor de pipeline sin log que lo reproduzca 🔴 ⛔
**Regla:** ningún valor de pipeline (sección B del Registro: ΔBIC, χ²_r, H₀ posterior,
r_d, θ*, σ₈…) puede ser canónico si no existe un **log committeado** en `results/logs/`
que lo reproduzca. El valor del Registro debe **aparecer en ese log**, no estar
escrito a mano. La fuente de verdad es el log, no la tabla.
**Por qué (caso real, F1 2026-06-14):** el Registro daba ΔBIC=−24.7 y χ²_r TT=1.045
(supuesto re-anclaje Σm_ν=0.069) — valores que **no aparecían en ningún log** y que una
re-corrida CAMB desmintió (reales: −28.0 / 1.044, = los papers). El "source of truth"
estaba mal y casi se "corrigen" los papers a un valor falso. Causa: el guardián
verificaba constantes algebraicas (sección A) y memorias, pero **no** los valores de
pipeline (sección B).
**Verificación (automatizada):** la **Capa Procedencia** del guardián
(`src/verificacion/ssee_verify.py`) extrae el valor de cada fila de la sección B y exige
que aparezca en su log committeado; si no coincide → **ROJO** (probado con prueba
negativa); sin log → **ABIERTO** (grieta visible, no bloquea). Para añadir un valor de
pipeline: guardar su log en `results/logs/` y registrarlo en `PIPELINE_PROVENANCE`.

## R10 — Overclaim "zero-parameter" global 🟠 ⛔ (automatizada)
**Regla:** ningún paper afirma ser un "zero-parameter framework/model" de forma global.
El claim honesto es **scoped**: "el sector background tiene cero parámetros ajustados".
**Por qué (caso real, 2026-06-14):** el TÍTULO de Paper 1 seguía siendo "A Zero-Parameter
Framework" (+cita en P2), contradiciendo el reframe aprobado minimal-parameter y su propio
abstract. La lectura hostil lo encontró; la automatización lo blindó.
**Verificación:** Capa Manuscritos del guardián, patrón asertivo
`(achieves|is|provides…) a zero-parameter (framework|model|theory)`. NO marca claims
scoped, recantaciones ("described as… however"), ni citas de títulos.

## R11 — Conteo de la serie congelado 🟠 ⛔ (automatizada)
**Regla:** las referencias al alcance de la serie deben reflejar los **10 papers + 2
journals** actuales, no un estado anterior.
**Por qué (caso real, 2026-06-14):** el abstract de P1 decía "extensiones en Papers~3--7"
mientras el cuerpo ya decía "3--10". Stale interno.
**Verificación:** Capa Manuscritos, patrón `papers 3--[789]` y `(seven|eight|nine)-paper
series`. NO marca referencias correctas de un paper a sus previos (P8 «1--7», P10 «1--9»).

## R12 — Vigente vs archivado, con bitácora 🟠 ⛔ (automatizada)
**Regla:** cada cajón está en una de dos clases — **VIGENTE** (siempre actual, auditado:
`manuscript/`, `docs/`, `src/`, Registro, memorias…) o **ARCHIVO** (histórico, no se edita
ni se cita como vigente). Lo obsoleto se mueve a `archive/` **con entrada en la Bitácora de
Archivado** (`archive/README.md`: qué, cuándo, por qué, qué lo reemplaza). Nada obsoleto se
queda en un cajón vivo sin marcar; nada se archiva sin entrada.
**Por qué (caso real, 2026-06-14):** el `archive/` mezclaba PDFs viejos, código superado y
figuras sin un "cuándo/por qué" completo; el README de archivo era parcial e inexacto. Sin
un mapa de vigencia, lo viejo se confunde con lo vivo.
**Verificación:** **Capa Archivo** del guardián — toda subcarpeta de `archive/` debe estar
documentada en la bitácora; cajón sin entrada → **ROJO** (probado con prueba negativa).
Mapa de Vigencia completo en `archive/README.md`.

---

## Protocolo de uso
1. Antes de cada envío, recorrer R1–R9 con sus comandos de verificación.
2. Cualquier hallazgo nuevo de "esto se pasó por alto" → **añadir una regla R-n aquí**.
3. La auditoría referee-hostil de los 10 papers usa esta lista como base.
4. Correr el guardián (`src/verificacion/ssee_verify.py`) — debe dar VERDE — antes de
   cualquier sello, commit de resultados o actualización de Zenodo.

## R13 — No vaciar el hueco al quitar un nombre 🟠
**Regla:** en los papers, un nombre del sistema (MIKAEL_V, KRYSTOS, DNAV…) codifica una
**función/ley** del sistema SSEE. Al quitarlo NO se deja un hueco: o (a) se conserva en la
columna de linaje/rol-de-sistema claramente etiquetada como andamiaje heurístico (como hace
la tabla de notación de P1), o (b) si sale de la prosa, se reemplaza por su etiqueta neutra
**ya establecida** (p.ej. "geometric form"), nunca por nada. La cara es álgebra; la función
vive en la capa de linaje + el diccionario citado.
**Por qué (caso real, 2026-06-14):** quité "MIKAEL_V" de la tabla de notación dejando a M_v
como única fila sin nombre de sistema → inconsistente y con el hueco vacío. Revertido. Las
leyes (copia / no-auto-suma) se pierden si se borra el nombre sin reemplazo que las preserve.
