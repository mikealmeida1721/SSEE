# Auditoría de contradicciones entre los 10 papers SSEE

**Fecha:** 2026-09-25 · **Método:** solo lectura, sin editar `manuscript/`
**Alcance:** `SSEE_Paper1..10` (.tex) + `SSEE_Sealed_Journal.tex`, `SSEE_Unified_Journal.tex`,
`SSEE_Endorser_Summary.tex`, `SSEE_Cover_Letter.tex`
**Referencias:** `CANONICAL_VALUES.yaml`, `RETRACCIONES.yaml`,
`auditoria_papers_zenodo_2026-09-20.md` (15 hallazgos previos — no se duplican, se verifica su estado)
**Estado del árbol:** el checkout está SOBRE la rama `fix/etiquetas-h0-ir`
(commits `3a9115c` + `78fa756`). "Cubierto por la rama" = ya corregido en el texto actual.

## Resultado ejecutivo

**19 hallazgos nuevos o persistentes: 6 CONTRADICCIÓN · 12 RESIDUO · 1 SOSPECHA.**
**0 de los 19 están cubiertos por la rama** — la rama corrigió 4 de los 15 del 9/20
(#1 tabla k_fs Paper 1; #8 tabla Paper 9 → `H_0^{glob,IR}`; #9 etiquetas Paper 9;
#14 Endorser → `H_0^{glob,IR}`), pero ninguno de los de abajo.

Los 4 bloques que Mike debe ver antes de la versión Zenodo / envío a revista:
1. **δc falsificado sigue vivo** en Paper 4 §deltac y Paper 5 §JWST (A1–A3).
2. **Paper 9 calcula la edad del universo con el bug retirado** (0.160 en E(z)) y cita
   como justificación el criterio two-Ω_m que Paper 1 ya retiró (A4).
3. **El signo del ΔBIC sigue opuesto** entre Paper 1 (−6.43) y Paper 2 (+6.43) (A5),
   y la tabla `model_comp` de Paper 2 no reconcilia aritméticamente (A6).
4. **Residuos de etiquetas** 68.13 sin IR y notación Ω_{m,dyn}/Ω_DE=0.839950 como
   densidades en Paper 3 y Paper 4 (A7–A17).

---

## CONTRADICCIÓN (6)

### A1. Paper 4 §deltac presenta el postulado falsificado como predicción viva
- **Archivo:** `manuscript/SSEE_Paper4_ToE.tex`, l.478–498 (ecuación l.486–488, `\label{eq:deltac}`)
- **Cita:** `\boxed{\delta_c^{\rm SSEE} = \delta_{c,\rm EdS}\times n_s = 1.6865\times 0.96556 = 1.6284}` …
  "this gives a halo-count enhancement of order ×1.2–1.5 relative to ΛCDM,
  providing a qualitative explanation for the anomalously massive galaxies
  detected by JWST at z≳10."
- **Por qué es contradicción:** el 2026-09-25 se falsificó la ruta dinámica —
  el colapso top-hat con DE suave en GR da δc_SSEE = 1.676 (z=0) → 1.68647 (z=10),
  ≈ ΛCDM; el factor ×n_s no se deriva. El enhancement JWST se evapora (ratio honesto ≈1.00×).
- **Rama:** NO cubierto (el diff de la rama en Paper 4 no toca §deltac).
- **Fix:** retirar §deltac / mover a OP-27 según `notes/2026-09-25_crecimiento_alto_z/borrador_op27_deltac.md`.

### A2. Paper 4 conclusiones listan δc=1.6284 como resultado
- **Archivo:** `manuscript/SSEE_Paper4_ToE.tex`, l.772
- **Cita:** `\item $\delta_c = \delta_{c,\rm EdS}\times n_s = 1.6284$ (qualitative JWST link)`
- **Rama:** NO cubierto. **Fix:** mismo que A1 (quitar el bullet o marcarlo retirado).

### A3. Paper 5 §JWST atribuye el enhancement al δc falsificado
- **Archivo:** `manuscript/SSEE_Paper5_IS.tex`, l.1402–1440 (`\label{subsec:JWST}`; l.1408, l.1437–1438)
- **Cita:** "The enhancement survives and grows with redshift, **driven by the modified
  collapse threshold** $\delta_c^{\mathrm{SSEE}}=1.6284$ against the EdS $1.686$."
  + tabla de ratios 1.051–1.887.
- **Por qué es contradicción:** el motor del enhancement (δc=1.6284) está falsificado;
  con el δc dinámico el ratio honesto es ≈1.00× en todo el rango. (La subsección sí
  corrigió honestamente el otro bug, el 0.160050 como densidad — ese fix se conserva.)
- **Rama:** NO cubierto (la rama no toca Paper 5). **Fix:** rehacer la tabla con δc dinámico o retirar la subsección a OP-27.

### A4. Paper 9 calcula la edad del universo metiendo 0.160050 en E(z)
- **Archivo:** `manuscript/SSEE_Paper9_HubbleTension.tex`, l.924–946
- **Cita:** "the dynamical-sector matter density $\Omega_m = \Omega_{m,\rm dyn} = 0.160050$
  **(Paper~1, Sec.~1.4, Two-$\Omega_m$ Criterion: late-time $E(z)$ in CPL form uses
  the dynamical value)**" → `E^2(z) = Ω_m(1+z)^3 + …` → `t_0^{SSEE} = 15.52 Gyr`.
  "The SSEE universe is 1.72 Gyr older than ΛCDM… A measurement of t_0 > 14 Gyr
  would favour SSEE over ΛCDM at this level."
- **Por qué es contradicción:** es exactamente el bug retirado — Paper 1 §two_omega_m
  (rama actual): "Ω_{m,dyn} was never a density… enters the Friedmann equations
  through the pressure term, **never as a matter density**"; Paper 2 §"Friedmann
  background (ω_m-direct)": "There is one matter density, Ω_m=0.308881". Además cita
  como justificación viva el criterio que Paper 1 declara **withdrawn**.
- **Rama:** NO cubierto (la rama en Paper 9 solo tocó 4 líneas de etiquetas H0).
- **Fix:** recalcular t_0 con Ω_m=0.308881 o retirar el párrafo completo.

### A5. Signo del ΔBIC (DESI) sigue opuesto entre documentos [persiste 9/20 #1–#3, #11, #13, #15]
- **Archivos:**
  - `SSEE_Paper1_Framework.tex` l.812 (+ l.81, 810, 877, 894, 996, 1008, 1323):
    `$\Delta\mathrm{BIC}=-6.43$ favouring SSEE`
  - `SSEE_Paper2_MCMC.tex` l.967 (caption Tab. model_comp):
    "SSEE is the reference; **positive Δ favours SSEE**" → `+6.43`
  - `SSEE_Endorser_Summary.tex` l.62: `$\Delta\text{BIC}=-6.43$ (SSEE favoured…)`
  - `SSEE_Sealed_Journal.tex` l.644 y `SSEE_Unified_Journal.tex` l.330: `+6.43`
- **Por qué es contradicción:** mismo estadístico, signos opuestos según el documento.
  El canónico (`CANONICAL_VALUES.yaml`) fija la convención: ΔBIC(ΛCDM−SSEE)=+6.43, >0 favorece SSEE.
- **Rama:** NO cubierto. **Fix:** unificar a +6.43 (convención canónica) en toda la suite y declararla una sola vez.

### A6. Tab. model_comp de Paper 2 no reconcilia aritméticamente [persiste 9/20 #4]
- **Archivo:** `manuscript/SSEE_Paper2_MCMC.tex`, l.968–978 (tabla), l.983, l.1301, l.1324
- **Cita:** columna BIC: SSEE 17.24, ΛCDM 22.92 → 22.92−17.24=**5.68** ≠ ΔBIC **+6.43**
  citado. Columna DIC: SSEE **15.69**, pero el apéndice (l.1324) lista DIC 14.95/20.61/18.97
  con +5.66, y el texto (l.983, l.1301) cita ΔDIC=−5.66 (= 20.61−14.95, no 20.61−15.69=4.92).
- **Rama:** NO cubierto. **Fix:** reconciliar las columnas BIC/DIC con los valores
  canónicos (+6.43 / −5.66) o añadir nota que explique el estimador.

---

## RESIDUO (12)

### A7. Paper 4 l.63–67: 68.13 con etiqueta desnuda [persiste 9/20 #6]
- `$H_0^{\rm glob} = H_0^{\rm SH0ES}(1-f_{\rm scr}) = 68.13$` en las predicciones
  enumeradas, sin subíndice IR. Los hunks de la rama en Paper 4 no cubren estas líneas.
- **Fix:** `H_0^{\rm glob,IR}`.

### A8. Paper 4 l.249: fila de tabla con etiqueta desnuda [9/20 #7 parcial]
- La rama corrigió la nota al pie ("Canonical IR-cascade value"), pero la fila sigue
  como `$H_0^{\rm glob}$ [km/s/Mpc] … 68.13`. **Fix:** IR en la etiqueta de la fila.

### A9. Paper 4 l.773 (conclusiones): 68.13 desnudo
- `\item $H_0^{\rm glob} = H_0^{\rm SH0ES}(1-f_{\rm scr}) = 68.13$ …`. **Fix:** subíndice IR.

### A10. Paper 3 l.150–159: define Ω_DE=0.839950 y Ω_{m,dyn}=0.160050 como densidades
- `\Omega_{\mathrm{DE}} = T_r/M_v = 0.839950` (`\label{eq:OmDE}`),
  `\Omega_{m,\mathrm{dyn}} = 1 - \Omega_{\mathrm{DE}} = 0.160050` (`\label{eq:Omm}`),
  presentadas como "density parameters [that] follow algebraically", sin marca de
  "naive" en la definición (el paper lo reenmarca como naive después, l.989).
  Contradice la doctrina de la rama (s_DE/s_m: cantidades del sector EoS, no densidades).
- **Rama:** NO cubierto (solo tocó 1 bloque de H0 en Paper 3). **Fix:** reetiquetar
  como s_DE/s_m con nota "naive reading" en la definición.

### A11. Paper 3 abstract l.44: símbolo retirado para el sector background
- "the background sector $(w_0, w_a, \Omega_{m,\rm dyn})$ carries zero fitted
  dimensionless parameters". **Fix:** $(w_0, w_a)$ (+ nota s_m si hace falta).

### A12. Paper 3 l.957–959: "the dynamical-sector density"
- "the raw SSEE dynamical background ($\Omega_{m,\rm dyn}=0.160050$)… the
  dynamical-sector **density** is not the CMB matter density" — revive la categoría
  retirada (Paper 1: *never a density*). **Fix:** "the EoS-sector number s_m=0.160050".

### A13. Paper 5 l.131: cita el criterio retirado como regla viva
- "see Paper~1, Sec.~1.4 (Two-$\Omega_m$ Criterion) for the rule that selects this
  value over $\Omega_{m,\rm dyn}=0.160050$ on clustering scales". El valor elegido
  (0.308881) es correcto; la cita apunta a la sección que hoy declara el criterio
  **withdrawn**. **Fix:** citar la sección ω_m-direct o §two_omega_m como retirada.

### A14. Paper 5 Q2: lenguaje "two matter sectors" (auto-corregido pero persistente)
- l.83: "the CMB-to-dynamical matter enhancement $\Omega_{m,CMB}/\Omega_{m,dyn}\approx2$,
  or are the **two matter sectors** independent"; l.111: "$\Omega_{m,eff}/\Omega_{m,dyn}
  = 0.989\pm0.017$"; l.849–853: "the numerical MIRA is $\mathrm{MIRA}_{num} =
  \Omega_{m,eff}/\Omega_{m,dyn}$". La sección se auto-corrige (l.115–120: "That
  construction is retracted"), pero la notación de dos sectores persiste.
- **Fix:** reenmarcar Q2 sin "two sectors" (p.ej. "the ratio once written …").

### A15. Sealed/Unified: 68.13 con etiqueta desnuda [persisten 9/20 #10, #12]
- `SSEE_Sealed_Journal.tex` l.573; `SSEE_Unified_Journal.tex` l.142, 916, 952, 1130.
  Contexto inmediato lo califica (IR/floor), pero la etiqueta invita al mal citado.
  La rama no toca estos archivos. **Fix:** subíndice IR en cada aparición.

### A16. Paper 9 l.196: label rancio
- `\sDE &= -w_0 = … = 0.839950\,. \label{eq:OmDE}` — la ecuación ya usa s_DE pero el
  label sigue diciendo OmDE. **Fix:** `\label{eq:sDE}`.

### A17. Paper 9 l.243–249: lenguaje two-Ω_m residual
- "the CMB-to-dynamical matter ratio (which is $\Omega_{m,CMB}/\Omega_{m,dyn}\approx1.93$,
  a numerical near-coincidence)". La afirmación sustantiva (MIRA ≠ ese ratio) es
  correcta, pero mantiene vivo el ratio entre una densidad y un no-densidad.
  **Fix:** reformular sin el ratio ("…not identified with any matter ratio; there is
  one matter density").

### A19. Paper 4 l.476: falsa precisión
- "Planck 2018: $n_s = 0.9649\pm0.0042$; deviation $= 0.160050\sigma$."
  El valor real es ≈0.157σ → 0.16σ; escribir 0.160050σ afirma 6 decimales que no
  existen (y evoca el número EoS). **Fix:** `$0.16\sigma$`.

---

## SOSPECHA (1)

### A18. Paper 1 l.817: referencia cruzada posiblemente colgada
- "the earlier 1.861 was an artefact of the cold sector 0.160050 in E(z)
  **(Paper~2~§5.3)**". El §5.3 está hardcodeado; la estructura de Paper 2
  (apéndices B.x, secciones renombradas) puede haberlo desplazado.
- **Fix (Claude, Trabajo 1):** verificar que el label existe y apunta al Ly-α audit.

---

## Verificación de falsos positivos (barridos que salieron limpios)

- **k_fs como falsificador vivo:** solo aparece en contextos de retirada
  (P1 l.485/521 "What is lost"; P3 l.273; P8 §withdrawn l.715–730, 904–907). Limpio.
- **β_c:** P8 usa el acoplamiento **disformal** βc≡+AURA con disambiguación explícita
  del conformal retirado de P7 (l.231–238) — uso deliberado, no residuo. P7 §withdrawn intacto.
- **φ-DM / 40.70 eV / 4/7.4/15 eV:** solo en contexto de retracción (P6 l.92–98, 170–216).
  Limpio.
- **H_MIRA=67.037:** solo con marca de superseded (P9 l.86, 547; P10 l.121). Limpio.
- **Paper 2 §"Friedmann background (ω_m-direct)"** (l.204–226): correcto —
  "There is one matter density, Ω_m=0.308881… s_DE/s_m do not enter the background densities".
- **Paper 2 "Scenario A"** (l.1386–1401): presenta 0.160050 en E(z) **como ilustración
  explícita del category error** ("demonstrating the category error directly… the same
  misplacement that produced the spurious background tensions in superseded drafts").
  Pedagógico y marcado — no es contradicción.
- **Paper 5 §JWST:** la corrección del bug 0.160050-como-densidad (2026-09-07) está
  honestamente narrada in situ; solo el motor δc quedó obsoleto (A3).
- **Valores cruzados** w0=−0.840/wa=−0.670, n_s=0.965558, r_d=147.17, H0=67.962,
  S8 KiDS-Legacy: consistentes entre papers (salvo redondeos de titular legítimos).

---

## Mapa rama vs hallazgos

| Hallazgo 9/20 | Estado hoy |
|---|---|
| #1 (P1 tabla k_fs) | ✅ Corregido por la rama |
| #2/#3/#11/#13/#15 (signo ΔBIC) | ❌ Persiste → A5 |
| #4 (P2 aritmética tabla) | ❌ Persiste → A6 |
| #5 (P2 convención mixta) | ❌ Persiste → ver A5 |
| #6 (P4 l.63-67) | ❌ Persiste → A7 |
| #7 (P4 nota l.255) | ⚠️ Parcial (nota OK, fila desnuda) → A8 |
| #8 (P9 tabla l.417) | ✅ Corregido por la rama |
| #9 (P9 etiquetas) | ✅ Corregido por la rama |
| #10 (Sealed l.573) | ❌ Persiste → A15 |
| #12 (Unified) | ❌ Persiste → A15 |
| #14 (Endorser 68.13) | ✅ Corregido por la rama |

**Nuevos de esta auditoría (post-9/20):** A1, A2, A3 (δc, 2026-09-25), A4 (edad P9),
A10, A11, A12, A13, A14, A16, A17, A18, A19.
