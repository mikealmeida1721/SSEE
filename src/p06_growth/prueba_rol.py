#!/usr/bin/env python3
"""
PRUEBA POR ROL de los candidatos en ventana — no por valor.

Objecion de Mike (correcta): listar candidatos por su NUMERO es trampa. La
gramatica original era valida porque exigia que cada factor tuviera ya un ROL
(g^2*v: X = acoplamiento radiativo/disipativo, Y = escala/ancla de vacio). Si al
quitar L4 se aceptan factores sin rol, no se esta "relajando una regla": se esta
quitando la fisica que hacia que la construccion significara algo.

Roles y leyes tomados de zenodo_dictionary/SSEE_Constant_Dictionary.md y
CLAUDE.md — NO asignados aqui.

Tres tests, en orden:
  T1  LEY DE COPIA (Law I): los objetos-copia n*AURA {MIRA, DUAL, TRIAL, CUARTAL}
      son PAREDES dimensionales terminales — "they are not re-combined into one
      another (else DUAL+DUSTAL would explode)". No pueden ser factores.
  T2  ROL DE X (el que va al cuadrado): debe ser acoplamiento radiativo/disipativo.
  T3  ROL DE Y (el que multiplica): debe ser escala / ancla de vacio.
"""
import numpy as np

PHI = 1.618033988749895; PI = 3.141592653589793; OMEGA = PHI+PI
SMNU = 0.06849
BIAL = OMEGA/2; KAL = BIAL+PI; SOLAR = BIAL+KAL; MAR = OMEGA+PI
VITA = PI+KAL; ANMA = BIAL+VITA; PYROS = OMEGA+PHI; IGNIS = PI+PYROS
KRYSTOS_V = PHI+PI+OMEGA; MIKA = KRYSTOS_V+PHI; AURA = PHI+BIAL
MIRA = AURA/2; DUAL = AURA*2; TRIAL = AURA*3; CUARTAL = AURA*4
ATLAS = PHI+PI+KRYSTOS_V                      # = 3*OMEGA = M_v, invariante 4D
PHITA = VITA+PHI; BUFFER = ATLAS-TRIAL
FAMILY = {"PHI": PHI, "PI": PI, "OMEGA": OMEGA, "BIAL": BIAL, "KAL": KAL,
          "SOLAR": SOLAR, "MAR": MAR, "VITA": VITA, "ANMA": ANMA,
          "PYROS": PYROS, "IGNIS": IGNIS, "KRYSTOS_V": KRYSTOS_V,
          "MIKA": MIKA, "AURA": AURA, "MIRA": MIRA, "DUAL": DUAL,
          "TRIAL": TRIAL, "CUARTAL": CUARTAL, "ATLAS": ATLAS,
          "PHITA": PHITA, "BUFFER": BUFFER}

# Law I — copias de AURA: paredes dimensionales TERMINALES (no recombinan)
COPIAS = {"MIRA": "½·AURA (pared ½; rol fisico f_screen)",
          "DUAL": "2·AURA (pared 2D)", "TRIAL": "3·AURA (pared 3D, horizonte de lo medible)",
          "CUARTAL": "4·AURA (pared 4D)"}

# Rol ACOPLAMIENTO (radiativo/disipativo) — con procedencia
ACOPL = {"BIAL": "Base Coupling Scalar; primer pulso radiativo",
         "KAL": "Structural Viscosity; retencion, ζ̃=KAL₀/3 (Paper 5 IS)",
         "SOLAR": "BIAL+KAL: hijo de dos termicos (radiativo-disipativo)",
         "AURA": "βc acoplamiento conforme (Papers 7/8)"}

# Rol ESCALA / ancla de vacio — con procedencia
ESCALA = {"OMEGA": "Stability Metric (escala estructural global)",
          "PYROS": "Dynamical Evolution Scalar (escala de evolucion, wₐ)",
          "KRYSTOS_V": "structuring; ancla de wₐ (escala de vacio DE)",
          "ATLAS": "3Ω = M_v, invariante de integracion 4D; ancla de w₀"}

# Sin rol de acoplamiento NI de escala declarado en ningun paper
SIN_ROL = {"IGNIS": "rol declarado = 'disruptive' (ni acoplamiento ni vacio)",
           "VITA": "invariante derivado (π+KAL); sin rol de acoplamiento",
           "ANMA": "invariante derivado (BIAL+VITA); sin rol",
           "MIKA": "invariante derivado (KRYSTOS_V+φ); sin rol",
           "MAR": "invariante derivado (Ω+π); sin rol de acoplamiento",
           "PHITA": "invariante derivado (VITA+φ); sin rol",
           "BUFFER": "ATLAS−TRIAL (resta); sin rol",
           "PHI": "primordial, no es acoplamiento ni escala de vacio",
           "PI": "primordial, no es acoplamiento ni escala de vacio"}

M_LO, M_HI = 1026.0, 13360.0
K50 = lambda m: 0.00603*m**1.087


def evalua(x, y):
    """devuelve (pasa, motivo del primer fallo)"""
    for f, papel in ((x, 'X'), (y, 'Y')):
        if f in COPIAS:
            return False, f'T1 ley de copia: {f} es {COPIAS[f]} — terminal'
    if x not in ACOPL:
        return False, (f'T2 rol de X: {x} no es acoplamiento — '
                       f'{SIN_ROL.get(x, ESCALA.get(x, "sin rol"))}')
    if y not in ESCALA:
        return False, (f'T3 rol de Y: {y} no es escala/vacio — '
                       f'{SIN_ROL.get(y, ACOPL.get(y, "sin rol"))}')
    return True, 'PASA los tres'


names = list(FAMILY)
en_ventana = []
for x in names:
    for y in names:
        v = FAMILY[x]**2 * FAMILY[y]
        if M_LO <= v <= M_HI:
            en_ventana.append((v, x, y))
en_ventana.sort()

print('=' * 88)
print('  PRUEBA POR ROL de los candidatos en ventana  M in [1026, 13360]')
print('=' * 88)
print(f'  candidatos que caen en ventana solo por su VALOR: {len(en_ventana)}\n')
pasan = []
for v, x, y in en_ventana:
    ok, motivo = evalua(x, y)
    if ok:
        pasan.append((v, x, y))
    print(f'  {x+"²·"+y:26} M={v:9.2f}  m={SMNU*v:7.2f} eV  '
          f'{"✓ PASA" if ok else "✗"}  {motivo if not ok else ""}')

print(f'\n{"="*88}')
print(f'  candidatos en ventana por valor : {len(en_ventana)}')
print(f'  candidatos que PASAN por rol    : {len(pasan)}')
if pasan:
    for v, x, y in pasan:
        print(f'      {x}²·{y}  M={v:.2f}  m={SMNU*v:.2f} eV  '
              f'k50={K50(SMNU*v):.3f} h/Mpc')
else:
    print('  => NINGUNA construccion en ventana tiene el linaje correcto.')

# control positivo: el canonico debe PASAR los tres tests (aunque este fuera de ventana)
ok, mot = evalua('SOLAR', 'KRYSTOS_V')
print(f'\n  CONTROL POSITIVO — el canonico SOLAR²·KRYSTOS_V (M=594.28, fuera de')
print(f'  ventana por el dato): rol {"correcto" if ok else "INCORRECTO"} — {mot}')
print(f'  => el test no esta rechazando por capricho: acepta el linaje valido y')
print(f'     rechaza los que solo aciertan el numero.')
