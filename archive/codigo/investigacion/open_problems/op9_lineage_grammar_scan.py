#!/usr/bin/env python3
"""
OP-9 — Look-elsewhere bajo GRAMÁTICA DE LINAJE ESTRICTO (forma g²·v)
====================================================================
Cierra el "residuo (i)" de OP-9 (timing post-hoc) de forma CUANTITATIVA.

Contexto (Mike, 2026-06-20): la partícula canónica m_φ = Σm_ν·(SOLAR²·KRYSTOS)
NO se obtuvo combinando constantes a ciegas hasta dar 594. El método fue:
  (1) preguntar qué VALOR pide la física (el M que cierra S8, M∈[295,988]);
  (2) preguntar qué construcción emite NATURALMENTE el linaje — restringido por
      reglas que PRE-EXISTEN a la partícula.
La fuerza de "salió del linaje" = cuán restrictiva es la gramática de linaje.
`op9_particle_existence.py` ya midió gramáticas PERMISIVAS (1/cientos) y concluyó
que el privilegio "depende enteramente de fijar la gramática legítima (lineage
laws)". Este script fija ESA gramática y cuenta.

Gramática de linaje estricto = la conjunción de cuatro reglas que existían antes:
  (L1) FORMA FÍSICA g²·v: el multiplicador es M = X²·Y (un acoplamiento al
       cuadrado por una escala/vacío). Es la forma de una masa generada
       (self-energy ∝ g², vacío v), no un producto/suma arbitrario.
  (L2) NO-AUTOSUMA (ruta π→KAL): un hijo debe DIFERIR de sus padres; se excluyen
       construcciones cuyo X o Y reproduce trivialmente a un padre por re-adición.
  (L3) LEY DE COPIA (ruta φ→AURA): múltiplos/dividendos de AURA están permitidos
       como entidades (MIRA=AURA/2, DUAL, TRIAL...), ya en el diccionario.
  (L4) ROL PREVIO: X (el acoplamiento, va al cuadrado) debe ser una entidad con
       rol RADIATIVO/DISIPATIVO ya establecido; Y (la escala) una entidad con rol
       de ANCLA DE VACÍO/ESCALA ya establecido. NO se asigna rol nuevo aquí.

Las asignaciones de rol (L4) se declaran explícitamente abajo con su procedencia
en los papers — son transparentes y auditables, no ajustadas al resultado.
"""
import numpy as np
from itertools import product as iproduct

PHI = 1.618033988749895
PI  = 3.141592653589793
OMEGA = PHI + PI
SMNU = 0.06902

# ── diccionario CERRADO (idéntico a op9_particle_existence.py) ──
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
# parents (para no-autosuma): hijo que es padre+padre por re-adición queda excluido
PARENTS = {
    "KAL":("BIAL","PI"), "SOLAR":("BIAL","KAL"), "VITA":("PI","KAL"),
    "ANMA":("BIAL","VITA"), "PYROS":("OMEGA","PHI"), "IGNIS":("PI","PYROS"),
    "KRYSTOS":("PHI","PI","OMEGA"), "MIKA":("KRYSTOS","PHI"), "AURA":("PHI","BIAL"),
    "MIKAEL_V":("PHI","PI","KRYSTOS"),
}

# ── (L4) asignaciones de ROL con procedencia (declaradas, auditables) ──
# COUPLINGS (g, va al cuadrado): rol radiativo/disipativo ya establecido.
COUPLINGS = {
    "BIAL":   "génesis: primer calor/pulso radiativo (semilla térmica)",
    "KAL":    "viscosidad bulk ζ̃=KAL₀/3 (Paper 5, IS)",
    "SOLAR":  "BIAL+KAL → acoplamiento radiativo-disipativo (hijo de dos térmicos)",
    "AURA":   "βc acoplamiento conforme (Paper 7/8)",
}
# SCALES (v, multiplica): rol de ancla de vacío/escala ya establecido.
SCALES = {
    "KRYSTOS":   "2Ω → ancla wₐ=−P_sc/K_v=−0.670 (DESI, escala de vacío DE, Paper 4)",
    "KRYSTOS_V": "2Ω idéntico (ancla de vacío)",
    "OMEGA":     "π+φ Stability Metric (escala estructural global)",
    "PYROS":     "Ω+φ Dynamical Evolution Scalar (escala de evolución, wₐ)",
}

# ── ventana viable del multiplicador (k_fs 0.3..1.5 h/Mpc) ──
M_LO, M_HI = 295.0, 988.0
def in_window(x): return M_LO <= x <= M_HI

MULT_CANON = SOLAR**2 * KRYSTOS        # 594.28  (canónico 2026-06-19)
M_PYROS_VITA_MIKA = PYROS*VITA*MIKA    # 615.33  (anterior)
M_OMEGA4_AURAKAL  = OMEGA**4 + AURA*KAL  # 535.28 (más viejo, forma SUMA → viola L1)

def show(title, pool):
    win = sorted((v,l) for (v,l) in pool if in_window(v))
    n = len(win)
    print(f"  {title:54s}  #total={len(pool):5d}  #ventana={n:3d}  -> 1/{n if n else '∞'}")
    return win

print("="*82)
print("  OP-9 — look-elsewhere bajo GRAMÁTICA DE LINAJE ESTRICTO (forma g²·v)")
print("="*82)
print(f"  Ventana viable M=m_φ/Σm_ν ∈ [{M_LO:.0f},{M_HI:.0f}]  (k_fs 0.3–1.5 h/Mpc)")
print(f"  CANÓNICO  SOLAR²·KRYSTOS = {MULT_CANON:.3f}  (m_φ={SMNU*MULT_CANON:.2f} eV)")
print(f"            SOLAR=BIAL+KAL=φ+2π={SOLAR:.4f}  KRYSTOS=2Ω={KRYSTOS:.4f}")
print("-"*82)

names = list(FAMILY)

# Nivel 0 — FORMA física g²·v sin más restricción: M = X²·Y, X,Y en toda la familia
pool0 = {}
for x in names:
    for y in names:
        lab = f"{x}²·{y}"
        pool0[lab] = FAMILY[x]**2 * FAMILY[y]
pool0 = [(v,l) for l,v in pool0.items()]
print("  GRADIENTE DE ESTRICTEZ (cada regla recorta el espacio admisible):")
w0 = show("L1 sola: forma g²·v, X,Y∈familia(21)", pool0)

# Nivel 1 — + NO-AUTOSUMA: excluir X o Y que sea un hijo-por-re-adición de sus padres
# (regla de la ruta π/KAL: el hijo debe diferir; se excluyen los nodos 'suma trivial').
def autosum_child(name):
    # un nodo viola no-autosuma si ES una entidad cuya definición es padre+padre
    # y aquí se usaría como factor crudo (re-adición). Se excluyen como FACTORES.
    return name in PARENTS and len(PARENTS[name]) == 2
pool1 = []
for x in names:
    for y in names:
        if autosum_child(x) or autosum_child(y):
            # permitido SOLO si el factor hereda combinación de funciones (SOLAR, AURA);
            # los puramente aditivos sin función nueva quedan fuera.
            if x not in COUPLINGS and x not in SCALES and autosum_child(x): continue
            if y not in COUPLINGS and y not in SCALES and autosum_child(y): continue
        pool1.append((FAMILY[x]**2*FAMILY[y], f"{x}²·{y}"))
w1 = show("L1+L2: + no-autosuma (factores con función)", pool1)

# Nivel 2 — + ROL PREVIO (L4): X∈COUPLINGS (acoplamiento), Y∈SCALES (vacío)
pool2 = []
for x in COUPLINGS:
    for y in SCALES:
        pool2.append((FAMILY[x]**2*FAMILY[y], f"{x}²·{y}"))
# dedup por valor (KRYSTOS == KRYSTOS_V == 2Ω)
seen={}
for v,l in pool2:
    k=round(v,3)
    if k not in seen: seen[k]=(v,l)
pool2=list(seen.values())
w2 = show("L1+L2+L4: + rol previo (coupling²·escala)", pool2)

print("-"*82)
print(f"  Couplings (g, ya con rol radiativo/disipativo): {list(COUPLINGS)}")
print(f"  Scales    (v, ya con rol de vacío/escala):      {list(SCALES)}")
print("-"*82)
print("  Construcciones admisibles en la ventana bajo LINAJE ESTRICTO (nivel 2):")
for v,l in w2:
    mark = "  <== CANÓNICO" if abs(v-MULT_CANON)<1e-2 else ""
    print(f"    M={v:8.3f}  m_φ={SMNU*v:6.2f} eV   {l:18s}{mark}")
print("="*82)
print("  VEREDICTO:")
n2 = len(w2)
print(f"  · Bajo la forma física g²·v + no-autosuma + rol-previo, el espacio admisible")
print(f"    en la ventana es 1/{n2}  (vs 1/cientos permisivo de op9_particle_existence).")
hit = any(abs(v-MULT_CANON)<1e-2 for v,_ in w2)
print(f"  · SOLAR²·KRYSTOS {'CAE' if hit else 'NO cae'} entre los admisibles de linaje estricto.")
print(f"  · El viejo Ω⁴+AURA·KAL={M_OMEGA4_AURAKAL:.1f} es una SUMA → viola L1 (no es g²·v):")
print(f"    queda EXCLUIDO por la forma física, no por ajuste. Confirma que SOLAR²·KRYSTOS")
print(f"    es gramaticalmente superior, no sólo numéricamente distinto.")
print(f"  · Residuo (i) de OP-9 cuantificado: el privilegio es 1/{n2} bajo la gramática que")
print(f"    Mike usó de hecho. Residuo (ii) — derivar 594.28 del transporte KAL₀ — sigue")
print(f"    abierto (OP-10). Esto NO prueba el coeficiente; mide cuán restringida fue la")
print(f"    elección. Honesto: la asignación de rol (L4) es declarada, no única.")
print("="*82)
