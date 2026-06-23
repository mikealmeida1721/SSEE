#!/usr/bin/env python3
"""
PASO 2 — ¿hay un multiplicador SIMPLE en el diccionario SSEE cerca del valor
que prefieren los datos (M_best = 593.8, S8=0.00σ)?
============================================================================
Misma fórmula m_φ = Σm_ν·MULT; solo cambia el multiplicador. Si un MULT más
SIMPLE que PYROS·VITA·MIKA=615.33 cae cerca de 593.8, es una partícula MEJOR:
data-preferred (0.00σ vs 0.24σ) Y más derivable (OP-9). Ranking por simplicidad.
"""
import numpy as np
from itertools import combinations

PHI=1.618033988749895; PI=3.141592653589793; OMEGA=PHI+PI
SMNU=0.06902
M_BEST=593.81          # data-preferred (S8=KiDS, paso 1)
MULT_CUR=  None        # se calcula abajo

# diccionario cerrado (Genesis 5.12)
BIAL=(PHI+PI)/2; KAL=BIAL+PI; SOLAR=BIAL+KAL; MAR=OMEGA+PI; VITA=PI+KAL
ANMA=BIAL+VITA; PYROS=OMEGA+PHI; IGNIS=PI+PYROS; KRYSTOS=PHI+PI+OMEGA
MIKA=KRYSTOS+PHI; AURA=PHI+BIAL; MIRA=AURA/2; DUAL=AURA*2; TRIAL=AURA*3
CUARTAL=AURA*4; MIKAEL_V=PHI+PI+KRYSTOS; BUFFER=MIKAEL_V-TRIAL; KRYSTOS_V=2*OMEGA
FAMILY={"PHI":PHI,"PI":PI,"OMEGA":OMEGA,"BIAL":BIAL,"KAL":KAL,"SOLAR":SOLAR,
    "MAR":MAR,"VITA":VITA,"ANMA":ANMA,"PYROS":PYROS,"IGNIS":IGNIS,
    "KRYSTOS":KRYSTOS,"MIKA":MIKA,"AURA":AURA,"MIRA":MIRA,"DUAL":DUAL,
    "TRIAL":TRIAL,"CUARTAL":CUARTAL,"MIKAEL_V":MIKAEL_V,"BUFFER":BUFFER,
    "KRYSTOS_V":KRYSTOS_V}
roots={"PHI":PHI,"PI":PI,"OMEGA":OMEGA}
names=list(FAMILY)
MULT_CUR=PYROS*VITA*MIKA   # 615.33

# ── enumeración con SCORE de simplicidad (menor = más simple) ──
# El modelo prefiere: enteros pequeños × potencia (H0=3Ω²), razones de 2 (w0),
# productos cortos. Penaliza: sumas, potencias altas, muchos factores.
cands=[]   # (M, label, score)
def add(M,label,score):
    if M>0: cands.append((M,label,score))

# k · (potencia de raíz):  score = 1(power) + (k>1) + (n-2)*0.5
for rn,rv in roots.items():
    for n in range(2,8):
        for k in range(1,7):
            add(k*rv**n, f"{k}·{rn}^{n}" if k>1 else f"{rn}^{n}", 1 + (k>1) + max(0,n-2)*0.5)
# k · (producto de 2 entidades): score = 2 + (k>1)
for a,b in combinations(names,2):
    for k in range(1,7):
        add(k*FAMILY[a]*FAMILY[b], f"{k}·{a}·{b}" if k>1 else f"{a}·{b}", 2 + (k>1))
# producto de 3 entidades (VOLUMEN): score = 3
for a,b,c in combinations(names,3):
    add(FAMILY[a]*FAMILY[b]*FAMILY[c], f"{a}·{b}·{c}", 3)
# (potencia de raíz) + (producto de 2): forma del viejo Ω⁴+AURA·KAL — score = 4 (suma)
for rn,rv in roots.items():
    for n in range(2,8):
        for a,b in combinations(names,2):
            add(rv**n+FAMILY[a]*FAMILY[b], f"{rn}^{n}+{a}·{b}", 4)

def rep(M): return f"S8 σ-rank via m_φ={SMNU*M:.2f} eV"

print("="*82)
print("  PASO 2 — multiplicador SIMPLE cerca del data-preferred M_best=593.8 (S8=0.00σ)")
print("="*82)
print(f"  ACTUAL  PYROS·VITA·MIKA = {MULT_CUR:.3f} (m_φ={SMNU*MULT_CUR:.2f} eV, 0.24σ), score=3 (volumen)")
print(f"  Buscamos: M cerca de {M_BEST:.1f} con score MENOR (más simple) o igual y más cerca.")
print("-"*82)

# tolerancia: ±1.5% de M_best
tol = 0.015*M_BEST
near=[(M,l,s) for (M,l,s) in cands if abs(M-M_BEST)<=tol]
# dedup por valor
seen={}
for M,l,s in near:
    key=round(M,3)
    if key not in seen or s<seen[key][2]:
        seen[key]=(M,l,s)
near=sorted(seen.values(), key=lambda x:(x[2], abs(x[0]-M_BEST)))   # por simplicidad, luego cercanía

print(f"  Candidatos en M∈[{M_BEST-tol:.1f},{M_BEST+tol:.1f}] (±1.5%), ordenados por SIMPLICIDAD:")
print(f"  {'score':>5} {'M':>9} {'m_φ(eV)':>8} {'|ΔM|':>6}  construcción")
for M,l,s in near[:16]:
    d=abs(M-M_BEST)
    print(f"  {s:5.1f} {M:9.3f} {SMNU*M:8.2f} {d:6.2f}  {l}")
print("-"*82)
# ¿el más simple es más simple que el actual (score 3 volumen)?
if near:
    best=near[0]
    print(f"  MÁS SIMPLE cerca del data-preferred: {best[1]}  (M={best[0]:.3f}, score={best[2]})")
    if best[2] < 3:
        print(f"  → MÁS SIMPLE que PYROS·VITA·MIKA (score 3). Candidato MEJOR: data-preferred + más derivable.")
    elif best[2]==3:
        print(f"  → mismo nivel (volumen) pero ¿más cerca de 593.8 que 615.33 (Δ={abs(MULT_CUR-M_BEST):.1f})?")
    else:
        print(f"  → no hay nada más simple que un volumen cerca; 615.33 sigue compitiendo.")
print("="*82)
