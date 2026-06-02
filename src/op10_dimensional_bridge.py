"""
OP-10 raíz — ¿EXISTE el puente puro→dimensional?  (test de transmutación)
==========================================================================
Hipótesis de Mike (2026-05-30): cada parámetro libre del modelo (offset −22,
alpha_WDM, m_φ=5.6, H_alg-como-puro) aparece EXACTAMENTE en el cruce
adimensional→dimensional. El puente actual no es generativo: ancla a θ* y no
produce escalas nuevas. Pregunta: ¿hay UN puente que genere TODAS las escalas
desde M_Pl con coeficientes (φ,π) derivados?

Pista física: la forma natural de sacar una escala chica de M_Pl NO es
multiplicar por un número puro chico (numerología, jerarquías 10⁻⁶² no salen de
φ⁻ⁿ con n chico). Es TRANSMUTACIÓN DIMENSIONAL:
        escala_IR = M_Pl · exp(−f(φ,π))   ⟹   f = ln(M_Pl/escala_IR)
Un f O(1)..O(100) en el EXPONENTE genera la jerarquía (así sale Λ_QCD).

TEST DISCRIMINANTE (honesto):
  (1) ¿ln(M_Pl/escala) = f(φ,π) limpio con expresión CHICA?  (no fluke)
  (2) control de densidad: ¿cuántos monomios caen en la tolerancia? Si miles,
      un match no significa nada (lección del escaneo de y≈121 de hoy).
  (3) la prueba REAL: ¿la MISMA forma de f sirve para H₀, M y m_φ a la vez?
      Si cada escala pide un f no relacionado → es fiteo, no puente.

Si (1)+(3) salen NO → el puente exponencial no existe; SSEE es (por ahora)
teoría de RAZONES, no de escalas. Eso es un resultado, no una rendición.
"""
import numpy as np
import itertools
import sys
sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
from ssee_core import PHI, PI, OMEGA, BETA, KAL0, MIRA, AURA, T_R, M_V

P_SC = OMEGA + PHI          # 6.3776
K_V  = PHI + PI + OMEGA     # 9.5192

# ── escalas físicas (eV) ──────────────────────────────────────────────────
M_PL     = 1.220910e28      # Planck no reducida
M_PL_R   = 2.435e27         # Planck reducida
HBAR_H   = 1.4306e-33       # ℏH_MIRA  (H₀ como energía)
M_CUT    = 8.81e-3          # M = cutoff cinético P10 (DE), LIMPIO vía 5φ⁸ρ_crit
M_PHI_56 = 5.60             # ansatz P6 (numerológico)
M_SEESAW = np.sqrt(M_PL*HBAR_H)   # seesaw cosmológico √(M_Pl·ℏH)
SUM_MNU  = 0.0824           # contaminado OP-14

basis = {'φ':PHI,'π':PI,'Ω':OMEGA,'β':BETA,'KAL':KAL0,'MIRA':MIRA,
         'AURA':AURA,'T_R':T_R,'M_V':M_V,'P_sc':P_SC,'K_v':K_V}
names = list(basis); vals = np.array([basis[n] for n in names])

print("=" * 78)
print(" OP-10 RAÍZ — ¿existe el puente puro→dimensional? (transmutación)")
print("=" * 78)
print(f"  M_Pl = {M_PL:.4e} eV   M_Pl_red = {M_PL_R:.4e} eV")
print(f"  seesaw √(M_Pl·ℏH) = {M_SEESAW*1e3:.4f} meV")

targets = [
    ("H₀ (ℏH_MIRA)", HBAR_H),
    ("M = 8.81 meV (cutoff, LIMPIO)", M_CUT),
    ("seesaw 4.18 meV", M_SEESAW),
    ("m_φ = 5.6 eV (ansatz)", M_PHI_56),
    ("Σm_ν = 0.0824 eV (contam)", SUM_MNU),
]

# ════════════════════════════════════════════════════════════════════════
#  Para cada escala: L = ln(M_Pl/escala). Buscar f(φ,π) ≈ L de DOS formas:
#   (a) monomio producto (φ^a·π^b·...) — y reportar DENSIDAD combinatoria
#   (b) k·(producto chico)  con k entero pequeño (factores de loop 3, 8π², …)
# ════════════════════════════════════════════════════════════════════════
def scan_log(L, tol=0.002, max_deg=4, exp_range=range(-2,3)):
    """devuelve (mejores, densidad) — monomios producto ≈ L a <tol."""
    hits = []
    total_near = 0
    for combo in itertools.product(exp_range, repeat=len(names)):
        s = sum(abs(c) for c in combo)
        if s == 0 or s > max_deg:
            continue
        val = np.prod([vals[i]**combo[i] for i in range(len(names))])
        if val <= 0:
            continue
        err = abs(val - L)/L
        if err < 0.05:
            total_near += 1            # densidad: cuántos en ±5%
        if err < tol:
            expr = "·".join(f"{names[i]}^{combo[i]}" for i in range(len(names)) if combo[i])
            hits.append((err, expr, val, s))
    hits.sort()
    return hits, total_near

print("\n" + "─"*78)
print(" TEST 1+2 — ¿ln(M_Pl/escala) es f(φ,π) limpio?  (+ densidad combinatoria)")
print("─"*78)
results = {}
for name, S in targets:
    L = np.log(M_PL/S)
    Lr = np.log(M_PL_R/S)
    hits, dens = scan_log(L)
    results[name] = (L, hits)
    print(f"\n  {name}:")
    print(f"     ln(M_Pl/S)     = {L:.4f}   (M_Pl_red: {Lr:.4f})   log10 = {np.log10(M_PL/S):.3f}")
    if hits:
        # complejidad mínima entre los hits a <0.2%
        best = hits[0]
        print(f"     mejor monomio  = {best[1]}  = {best[2]:.4f}  ({best[0]*100:+.3f}%, grado {best[3]})")
        print(f"     #monomios a <0.2% = {len(hits)}   |   densidad (#en ±5%) = {dens}")
        if dens > 30:
            print(f"     ⚠ densidad ALTA ({dens}) → un match a este target es esperado por combinatoria, NO señal")
    else:
        print(f"     NINGÚN monomio (φ,π) grado≤4 a <0.2%   |   densidad (±5%) = {dens}")

# ════════════════════════════════════════════════════════════════════════
#  TEST 3 — la prueba REAL: ¿una forma COMÚN sirve para varias escalas?
#  Si el puente existe, los L deberían compartir estructura (ratios racionales
#  limpios entre sí, o un patrón k·base). Si son números sin relación → fiteo.
# ════════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" TEST 3 — ¿los exponentes L comparten estructura entre escalas?")
print("─"*78)
Ls = {name: np.log(M_PL/S) for name, S in targets}
print("  Ratios L_i / L_j (¿racionales limpios o (φ,π) limpios?):")
keys = list(Ls)
for i in range(len(keys)):
    for j in range(i+1, len(keys)):
        r = Ls[keys[i]]/Ls[keys[j]]
        # ¿r cerca de racional simple p/q (q≤6) o de constante base?
        best_rat = min(((abs(r-p/q), f"{p}/{q}") for q in range(1,7) for p in range(1,4*q+1)),
                       key=lambda t: t[0])
        flag = " ★ racional limpio" if best_rat[0]/r < 0.005 else ""
        print(f"     {keys[i][:22]:<22} / {keys[j][:22]:<22} = {r:.4f}  (≈{best_rat[1]}{flag})")

# ════════════════════════════════════════════════════════════════════════
#  CONTROL — el puente LINEAL (numerología): escala = f(φ,π)·M_Pl, f diminuto.
#  Mostrar que NO hay f(φ,π) chico que dé ~10⁻⁶¹ (por eso lineal falla).
# ════════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" CONTROL — puente LINEAL escala=f·M_Pl: ¿algún f(φ,π) chico da el ratio?")
print("─"*78)
for name, S in [("H₀", HBAR_H), ("M=8.81meV", M_CUT)]:
    ratio = S/M_PL
    print(f"  {name}: escala/M_Pl = {ratio:.3e}  → log10 = {np.log10(ratio):.2f}")
    # ¿qué exponente n en φ⁻ⁿ?  n = -log(ratio)/log(φ)
    n_phi = -np.log(ratio)/np.log(PHI)
    print(f"     necesitaría φ^(−{n_phi:.1f}) — exponente ENORME ⟹ lineal = numerología")

print("\n" + "=" * 78)
print(" LECTURA")
print("=" * 78)
print("""
  Si TEST 1 da expresiones CHICAS (grado ≤3) con densidad BAJA, y TEST 3
  muestra que las escalas comparten estructura → el puente exponencial EXISTE
  y es el camino a derivar H₀, M y m_φ de una sola ancla (M_Pl). Bloqueador #1
  del premio resuelto.

  Si TEST 1 solo da flukes de alta densidad/alto grado y TEST 3 da ratios sin
  estructura → el puente exponencial NO existe en (φ,π). SSEE es teoría de
  RAZONES; las escalas absolutas requieren input externo (como m_q en el SM).
  Resultado honesto, no rendición: dice DÓNDE están confinados los params libres.
""")
