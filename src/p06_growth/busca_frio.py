#!/usr/bin/env python3
"""
Busqueda de particula phi FRIA y PESADA bajo la gramatica de linaje SIN CAMBIOS.
Ventana y reglas PRE-DECLARADAS en PREDECLARACION_ventana_frio.md (leelo primero).

Diccionario, parentescos y roles: copiados literalmente de
archive/codigo/investigacion/open_problems/op9_lineage_grammar_scan.py (2026-06-20)
para que la gramatica sea LA MISMA que produjo el 594.28. Nada anadido.
"""
import numpy as np, sys

PHI = 1.618033988749895
PI = 3.141592653589793
OMEGA = PHI + PI
SMNU = 0.06849                                     # canonico C_nu=93.14

BIAL = (PHI+PI)/2; KAL = BIAL+PI; SOLAR = BIAL+KAL; MAR = OMEGA+PI
VITA = PI+KAL; ANMA = BIAL+VITA; PYROS = OMEGA+PHI; IGNIS = PI+PYROS
KRYSTOS = PHI+PI+OMEGA; MIKA = KRYSTOS+PHI; AURA = PHI+BIAL; MIRA = AURA/2
DUAL = AURA*2; TRIAL = AURA*3; CUARTAL = AURA*4; MIKAEL_V = PHI+PI+KRYSTOS
BUFFER = MIKAEL_V-TRIAL; KRYSTOS_V = 2*OMEGA
FAMILY = {"PHI": PHI, "PI": PI, "OMEGA": OMEGA, "BIAL": BIAL, "KAL": KAL,
          "SOLAR": SOLAR, "MAR": MAR, "VITA": VITA, "ANMA": ANMA,
          "PYROS": PYROS, "IGNIS": IGNIS, "KRYSTOS": KRYSTOS, "MIKA": MIKA,
          "AURA": AURA, "MIRA": MIRA, "DUAL": DUAL, "TRIAL": TRIAL,
          "CUARTAL": CUARTAL, "MIKAEL_V": MIKAEL_V, "BUFFER": BUFFER,
          "KRYSTOS_V": KRYSTOS_V}
PARENTS = {"KAL": ("BIAL", "PI"), "SOLAR": ("BIAL", "KAL"),
           "VITA": ("PI", "KAL"), "ANMA": ("BIAL", "VITA"),
           "PYROS": ("OMEGA", "PHI"), "IGNIS": ("PI", "PYROS"),
           "KRYSTOS": ("PHI", "PI", "OMEGA"), "MIKA": ("KRYSTOS", "PHI"),
           "AURA": ("PHI", "BIAL"), "MIKAEL_V": ("PHI", "PI", "KRYSTOS")}
COUPLINGS = ["BIAL", "KAL", "SOLAR", "AURA"]
SCALES = ["KRYSTOS", "KRYSTOS_V", "OMEGA", "PYROS"]

M_LO, M_HI = 1026.0, 13360.0                       # PRE-DECLARADA
K50 = lambda m: 0.00603 * m**1.087
DS8 = lambda m: 52087.9358 * m**-2.283


def report(titulo, pool):
    win = sorted((v, l) for v, l in pool if M_LO <= v <= M_HI)
    print(f'  {titulo:52s} #total={len(pool):5d}  #en ventana={len(win):4d}')
    return win


names = list(FAMILY)
print('=' * 84)
print('  BUSQUEDA DE PARTICULA FRIA — gramatica de linaje SIN CAMBIOS')
print('=' * 84)
print(f'  ventana PRE-DECLARADA:  M in [{M_LO:.0f}, {M_HI:.0f}]  '
      f'-> m_phi in [{SMNU*M_LO:.1f}, {SMNU*M_HI:.1f}] eV')
print(f'  canonico vigente:  SOLAR^2*KRYSTOS_V = {SOLAR**2*KRYSTOS_V:.2f}  '
      f'(m={SMNU*SOLAR**2*KRYSTOS_V:.2f} eV)  -> EXCLUIDO por el dato\n')

p0 = [(FAMILY[x]**2 * FAMILY[y], f'{x}²·{y}') for x in names for y in names]
w0 = report('L1 sola: forma g²·v, X,Y en familia(21)', p0)


def autosum(n):
    return n in PARENTS and len(PARENTS[n]) == 2


p1 = []
for x in names:
    for y in names:
        if autosum(x) and x not in COUPLINGS and x not in SCALES:
            continue
        if autosum(y) and y not in COUPLINGS and y not in SCALES:
            continue
        p1.append((FAMILY[x]**2 * FAMILY[y], f'{x}²·{y}'))
w1 = report('L1+L2: + no-autosuma', p1)

p2 = [(FAMILY[x]**2 * FAMILY[y], f'{x}²·{y}') for x in COUPLINGS
      for y in SCALES]
w2 = report('L1+L2+L4: + rol previo (acoplamiento²·escala)', p2)

print(f'\n  --- TODOS los candidatos L1+L2+L4, en orden, con su fisica ---')
print(f'  {"construccion":26} {"M":>10} {"m_phi(eV)":>10} {"k50":>8} '
      f'{"ds8":>8}  estado')
for v, l in sorted(p2):
    m = SMNU * v
    est = ('EN VENTANA' if M_LO <= v <= M_HI else
           'suprime demasiado' if v < M_LO else 'borde no observable')
    print(f'  {l:26} {v:10.2f} {m:10.2f} {K50(m):8.3f} {DS8(m):7.2f}%  {est}')

print(f'\n  MAXIMO alcanzable bajo L1+L2+L4 = {max(v for v,_ in p2):.2f}'
      f'   (necesario: >= {M_LO:.0f})')
print(f'\n  RESULTADO:')
for nm, w in (('L1 sola', w0), ('L1+L2', w1), ('L1+L2+L4', w2)):
    if w:
        print(f'    {nm:10} : {len(w)} candidato(s) -> look-elsewhere 1/{len(w)}')
        for v, l in w[:12]:
            print(f'        {l:26} M={v:10.2f}  m={SMNU*v:8.2f} eV  '
                  f'k50={K50(SMNU*v):.3f} h/Mpc')
    else:
        print(f'    {nm:10} : 0 candidatos en ventana')
