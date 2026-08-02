# Acomodo tras la retracción del sector φ-DM

> **Qué es esto.** Lista de trabajo de lo que quedó **descolocado** (no roto) al
> reescribirse Paper 6 el 2026-08-01. Cada ítem lleva **cómo se comprueba**, no
> solo qué hacer — para que nada entre por afirmación. Se tacha lo hecho.
>
> **Regla:** no se marca ✅ sin que la comprobación de su fila esté corrida y
> pegada. «Compila» no es comprobación (los warnings no son errores).
>
> Estado del guardián al abrir la lista: **VERDE 194/194** · memory_sync **0 drifts**.

---

## Bloque 1 — Descolocado y verificado (hallado por barrido, no de memoria)

| # | Qué | Dónde | Cómo se comprueba | Estado |
|---|---|---|---|---|
| ~~1.1~~ | ~~`.bib` con el título viejo de Paper 6~~ | `manuscript/ssee_unified.bib`, `submission_PRD/ssee_unified.bib` | `grep -c "Two-Sector"` → **0 y 0** | ✅ |
| ~~1.2~~ | ~~Endorser cita OP-9/11 como abiertos~~ | `SSEE_Endorser_Summary.tex:33` | `grep -c "OP-9/11"` → **0**; compila 0 err / 0 refs / 2 pp | ✅ |
| ~~1.3~~ | ~~CLAUDE.md con páginas viejas~~ | `CLAUDE.md`, tabla `docs/` | tabla **regenerada desde `pdfinfo`**, no a mano. También corregidas las 3 filas cuyo *estado* había quedado obsoleto (P1 §1.4, P6, P8) | ✅ |
| ~~1.4~~ | ~~`cover_galaxy.jpg` fuera de sitio~~ | — | **FALSO POSITIVO de mi barrido**: vive junto a su `.tex`, la ruta relativa resuelve | ✅ |
| ~~1.5~~ | ~~`ssee_v36_CMB_TT_comparison` fuera de sitio~~ | — | **FALSO POSITIVO**: el `.tex` usa ruta explícita `../class_ssee/output/` | ✅ |

**Barrido de figuras:** `src/verificacion/barrido_figuras.sh` → sale **vacío**.

> ⚠️ **Este barrido mintió DOS veces antes de servir**, y las dos quedan como
> lección:
> 1. Usaba `ls a.pdf a.png`, que **falla si falta cualquiera de los dos** aunque
>    el `.pdf` esté → **10 falsos positivos**. Corregido a `[ -f ]` separados.
> 2. Sólo miraba `results/figures/`, ignorando las rutas relativas al propio
>    `.tex` y las explícitas → **2 falsos positivos más** (1.4, 1.5).
>
> Moraleja anotada: **una herramienta de verificación hay que verificarla
> antes de creerle.** Las tres versiones «funcionaban»; sólo la tercera decía
> la verdad.

---

## Bloque 2 — Contramedidas que quedaron anotadas pero NO implementadas

| # | Qué | Por qué importa | Cómo se comprueba | Estado |
|---|---|---|---|---|
| ~~2.1~~ | ~~R45 al caso simétrico «OP adoptado que fue revertido»~~ | ídem | **implementado + auto-test con el encabezado REAL de OP-17 pre-arreglo** → lo marca; OP-1 y OP-5 exentos. Afinado tras 2 rondas de falsos positivos | ✅ |
| ~~2.2~~ | ~~R20 no cubría los scripts de figuras~~ | ídem | causa hallada: R20 buscaba el nombre EXACTO `kids_s8`; la figura usaba `S8_KIDS`. Ahora cada ancla lleva **alias** + tolera `= (v, err)`. **Auto-test con la línea real del bug** → la marca | ✅ |

---

## Bloque 3 — Físicas abiertas (NO son acomodo; van con R1/R2)

| # | Qué | Estado |
|---|---|---|
| 3.1 | **OP-21** — derivar `ζ̃/(τ_Π H₀) = Ω_DE` del Lagrangiano, o descartarlo. Diana **única**: la hipótesis `ζ̃=KAL₀/3` y el resultado `c²_s=0` son algebraicamente **la misma afirmación** (verificado: diferencia 0.00e+00), así que derivar una de la otra es circular | 🔴 abierto |
| 3.2 | **Veta A_s** — las tres determinaciones (CMB / cizalla / clustering). Requiere R1/R2 | ⏳ bloqueado por corridas |
| 3.3 | **OP-18** — derivar A_s. Bloqueado *por* 3.2: no se puede validar una derivación contra un blanco que tiene dos valores | ⏳ bloqueado |

---

## Lo que NO hay que hacer

- **No** tocar Paper 8 más: su predicción quedó **restaurada a incondicional**
  sobre base mejor (ω_c de OP-8 + α_B=α_M=0 de Paper 7). Verificado que ninguna
  de las dos patas venía de la partícula.
- **No** publicar ni sellar hasta que Mike termine la lectura página por página
  (regla R39, `project_publication_gate`).
- **No** buscar mecanismo para la veta A_s antes de R1/R2.

---

## Bitácora

| Fecha | Qué |
|---|---|
| 2026-08-02 | Lista abierta. Inventario levantado por barrido del repo, no de memoria. Guardián VERDE 194/194. |
| 2026-08-02 | **Bloque 1 cerrado** (5/5). 3 defectos reales + 2 falsos positivos de mi propio barrido. |
| 2026-08-02 | **Bloque 2 cerrado** (2/2). Guardián **VERDE 197/197** (piso 194→197). Hallazgo extra: OP-5 tenía el encabezado desfasado respecto a su fila del resumen — lo destapó el detector nuevo al dar un falso positivo sobre él. |
