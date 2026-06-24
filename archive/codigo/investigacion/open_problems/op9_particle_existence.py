#!/usr/bin/env python3
"""
OP-9 / ¿EXISTE la partícula φ-DM dentro de SSEE? — look-elsewhere del multiplicador
====================================================================================
Pregunta de Mike (2026-06-18): antes de enshrinar m_φ=42.47 (mult PYROS·VITA·MIKA
=615.33) en el manuscrito, ¿es una construcción PRIVILEGIADA del diccionario
cerrado, o cae por look-elsewhere entre cientos? ¿Es igual o mejor que el viejo
535.28 (Ω⁴+AURA·KAL)? Si es más privilegiada, OP-9 (origen UV) se vuelve más fácil.

Método (honesto, NO circular):
 - El multiplicador es un número puro M tal que m_φ = Σm_ν·M.
 - La VENTANA VIABLE no es "cerca de un target S8" (eso sería fiteo) sino
   "produce un φ-DM físicamente relevante": k_fs in [0.3,1.5] h/Mpc (afecta las
   escalas de S8 sin borrar estructura ni ser inerte) → M in [295, 988].
 - Se enumera el diccionario CERRADO (21 entidades nombradas + potencias de las
   3 raíces) bajo la gramática del sistema (potencias, productos, sumas-de-término)
   y se cuenta cuántas construcciones caen en la ventana. Eso es el look-elsewhere.
"""
import numpy as np
from itertools import combinations, product as iproduct

PHI = 1.618033988749895
PI  = 3.141592653589793
OMEGA = PHI + PI
SMNU = 0.06902

# ── diccionario CERRADO (Genesis 5.12, idéntico a look_elsewhere_family.py) ──
BIAL=(PHI+PI)/2; KAL=BIAL+PI; SOLAR=BIAL+KAL; MAR=OMEGA+PI; VITA=PI+KAL
ANMA=BIAL+VITA; PYROS=OMEGA+PHI; IGNIS=PI+PYROS; KRYSTOS=PHI+PI+OMEGA
MIKA=KRYSTOS+PHI; AURA=PHI+BIAL; MIRA=AURA/2; DUAL=AURA*2; TRIAL=AURA*3
CUARTAL=AURA*4; MIKAEL_V=PHI+PI+KRYSTOS; BUFFER=MIKAEL_V-TRIAL; KRYSTOS_V=2*OMEGA
FAMILY = {
    "PHI":PHI,"PI":PI,"OMEGA":OMEGA,"BIAL":BIAL,"KAL":KAL,"SOLAR":SOLAR,
    "MAR":MAR,"VITA":VITA,"ANMA":ANMA,"PYROS":PYROS,"IGNIS":IGNIS,
    "KRYSTOS":KRYSTOS,"MIKA":MIKA,"AURA":AURA,"MIRA":MIRA,"DUAL":DUAL,
    "TRIAL":TRIAL,"CUARTAL":CUARTAL,"MIKAEL_V":MIKAEL_V,"BUFFER":BUFFER,
    "KRYSTOS_V":KRYSTOS_V,
}
names = list(FAMILY)

# ── ventana viable del multiplicador (k_fs 0.3..1.5 h/Mpc) ──
M_LO, M_HI = 295.0, 988.0
def in_window(x): return M_LO <= x <= M_HI

# ── enumeración bajo la gramática SSEE ──
# (1) productos de 2 y 3 entidades distintas (la "ley de volumen": producto de
#     ceilings ortogonales). (2) potencias de las 3 raíces. (3) suma de un
#     término-potencia + un término-producto (la forma del viejo Ω⁴+AURA·KAL).
constructions = {}   # etiqueta -> valor

# potencias de raíces φ,π,Ω  (n=2..6)
roots = {"PHI":PHI,"PI":PI,"OMEGA":OMEGA}
for rn, rv in roots.items():
    for n in range(2, 7):
        constructions[f"{rn}^{n}"] = rv**n

# productos de 2 y 3 entidades distintas
for combo in combinations(names, 2):
    constructions["·".join(combo)] = FAMILY[combo[0]]*FAMILY[combo[1]]
for combo in combinations(names, 3):
    constructions["·".join(combo)] = (FAMILY[combo[0]]*FAMILY[combo[1]]*FAMILY[combo[2]])

# suma de potencia-de-raíz + producto-de-2 (forma del viejo multiplicador)
for rn, rv in roots.items():
    for n in range(2, 7):
        for combo in combinations(names, 2):
            constructions[f"{rn}^{n}+{combo[0]}·{combo[1]}"] = rv**n + FAMILY[combo[0]]*FAMILY[combo[1]]

# clase gramatical de cada construcción
def gclass(lab):
    if "+" in lab: return "suma"          # potencia + producto-2 (forma vieja)
    if lab.count("·") == 2: return "prod3" # VOLUMEN (3 ceilings)
    if lab.count("·") == 1: return "prod2" # área
    if "^" in lab: return "pot"            # potencia de raíz
    return "simple"

# dedup por valor redondeado, guardando la clase
seen = {}
for lab, val in constructions.items():
    if val <= 0: continue
    key = round(val, 4)
    if key not in seen:
        seen[key] = (val, lab, gclass(lab))
uniq = list(seen.values())

MULT_NEW = PYROS*VITA*MIKA
MULT_OLD = OMEGA**4 + AURA*KAL

print("="*78)
print("  ¿EXISTE la partícula φ-DM en SSEE? — look-elsewhere del multiplicador")
print("="*78)
print(f"  Diccionario cerrado: {len(FAMILY)} entidades nombradas (+ raíces φ,π,Ω)")
print(f"  Ventana viable: M = m_φ/Σm_ν ∈ [{M_LO:.0f}, {M_HI:.0f}]  (k_fs 0.3–1.5 h/Mpc)")
print(f"  NUEVO  PYROS·VITA·MIKA = {MULT_NEW:.3f} (m_φ={SMNU*MULT_NEW:.2f} eV) — clase: VOLUMEN (prod3 puro)")
print(f"  VIEJO  Ω⁴ + AURA·KAL   = {MULT_OLD:.3f} (m_φ={SMNU*MULT_OLD:.2f} eV) — clase: SUMA (pot4 + prod2)")
print("-"*78)

# conteo por TIER de estrictez gramatical (= cuán privilegiado, según la regla)
tiers = {
    "A permisivo  (pot + prod2 + prod3 + suma)": lambda c: True,
    "B productos  (prod2 + prod3, sin sumas/pot)": lambda c: c in ("prod2","prod3"),
    "C volúmenes  (prod3 puro: 3 ceilings)":       lambda c: c == "prod3",
    "D sumas      (forma del viejo: pot4+prod2)":   lambda c: c == "suma",
}
print(f"  {'Gramática (regla de construcción)':46s} {'#total':>7} {'#ventana':>9}  1/N")
for tname, tfilt in tiers.items():
    pool = [u for u in uniq if tfilt(u[2])]
    win  = [u for u in pool if in_window(u[0])]
    n = len(win)
    print(f"  {tname:46s} {len(pool):7d} {n:9d}  1/{n}")
print("-"*78)

# vecinos inmediatos de NUEVO en la clase VOLUMEN (prod3)
vols = sorted((v,l) for (v,l,c) in uniq if c=="prod3" and in_window(v))
idx = min(range(len(vols)), key=lambda i: abs(vols[i][0]-MULT_NEW))
print("  Vecindad de PYROS·VITA·MIKA entre los VOLÚMENES (prod3) de la ventana:")
for i in range(max(0,idx-2), min(len(vols),idx+3)):
    v,l = vols[i]; mark = "  <== NUEVO" if abs(v-MULT_NEW)<1e-3 else ""
    print(f"    M={v:8.3f}  m_φ={SMNU*v:6.2f} eV   {l:28s}{mark}")
print("="*78)
print("  VEREDICTO:")
print(f"  · PYROS, VITA, MIKA son entidades REALES del diccionario (sí 'existe').")
print(f"  · Bajo gramática PERMISIVA es 1-de-cientos -> NO privilegiada por look-elsewhere.")
print(f"  · Su privilegio depende ENTERAMENTE de fijar la gramática legítima (lineage laws).")
print(f"  · vs VIEJO: estadísticamente similar; gramaticalmente más limpia (volumen puro")
print(f"    vs suma pot4+prod2) -> forma más derivable para OP-9, pero NO más única.")
print("="*78)
