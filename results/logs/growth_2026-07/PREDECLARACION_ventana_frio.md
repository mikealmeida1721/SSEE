# PRE-DECLARACIÓN — búsqueda de partícula φ fría y pesada

**Fecha: 2026-07-30.** Escrito ANTES de enumerar candidatos. La validez de esta
búsqueda depende de que este archivo preceda al escaneo.

## Por qué se re-abre la búsqueda

El multiplicador canónico SOLAR²·KRYSTOS_V = 594,28 → m_φ = 40,70 eV queda
**excluido por el dato**, no por preferencia:

1. SSEE de un sector con A_s libre ajusta la cizalla cruda de KiDS-1000:
   S₈ = 0,75693 contra 0,759 ± 0,024 → **−0,08σ**. No hay tensión S₈ que resolver.
2. Luego un sector con free-streaming solo puede **degradar** el ajuste. Medido:
   m_φ = 40,70 eV suprime σ₈ un **10,42%**, contra una incertidumbre de KiDS
   del 3,16%. La supresión NO la absorbe A_s porque cambia la forma de P(k),
   no la amplitud.

## Ventana declarada (fijada ANTES de enumerar)

Fijada Ω_φDM = 0,148890 (no es libre: = ω_m/h² − 0,160), la masa determina la
temperatura y el borde. Leyes medidas con CLASS en m ∈ [40,7 ; 300] eV:

    k₅₀(h/Mpc) = 0,00603 · m^1,087          (m en eV)
    Δσ₈(%)     = 52.088  · m^−2,283

**Cota inferior** — la supresión debe ser invisible para KiDS (< 3,16%):

    m_φ > 70,3 eV        →  M = m_φ/Σm_ν > 1026     (Σm_ν = 0,06849 eV)

**Cota superior** — el borde debe seguir siendo falsable por DESI Y3 / Euclid,
que alcanzan k ≲ 10 h/Mpc:

    m_φ < 915 eV         →  M < 13.360

**VENTANA DECLARADA:  M ∈ [1026 , 13.360]**  (m_φ ∈ [70,3 ; 915] eV)

Nota: la cota inferior es un *mínimo del mínimo* — solo exige amplitud invisible.
El criterio de forma (T² > 0,99 para k < 3 h/Mpc) exige más masa aún. Si la
búsqueda usara el criterio de forma, la ventana sería MÁS restrictiva, nunca menos.

## Gramática declarada — SIN CAMBIOS respecto a la que produjo el 594,28

Se usa exactamente la de `op9_lineage_grammar_scan.py` (2026-06-20), con el mismo
diccionario cerrado de 21 entidades y las mismas asignaciones de rol:

- **L1** forma física g²·v: M = X²·Y
- **L2** no-autosuma (ruta π→KAL)
- **L4** rol previo: X ∈ COUPLINGS = {BIAL, KAL, SOLAR, AURA};
        Y ∈ SCALES = {KRYSTOS, KRYSTOS_V, OMEGA, PYROS}

**Compromiso explícito:** no se añaden entidades, no se añaden roles, no se cambia
la forma (X²·Y). Si ningún candidato cae en la ventana, el resultado es
"no existe candidato bajo la gramática vigente" — y relajar la gramática después
de ver esto sería ajuste post-hoc, no descubrimiento.

## Criterio de éxito, declarado

- **0 candidatos en ventana** → la gramática de linaje NO puede producir una
  partícula fría admisible. El segundo sector queda sin coeficiente derivable, y
  eso se declara como resultado negativo.
- **1 candidato** → selección forzada (1/1); el más fuerte posible.
- **N > 1 candidatos** → look-elsewhere 1/N, a reportar íntegro sin elegir a mano.
