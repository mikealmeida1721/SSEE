"""
OP-14 — Ataque por el ángulo NUEVO (portal Yukawa), no por el offset 22
========================================================================
El ataque de 2026-05-23 mató la ruta "derivar 4·KAL−22":
  H1 frágil (25% drift @ 1e-3), H2 sin DoF match, H3 mejor candidato
  R≈1/[3(KAL−φ)] a −0.27% (aprox, no identidad).

Hoy el script op10_mechanisms_2to6 reveló una relación FÍSICA m_φ↔m_ν
que NO pasa por la fórmula rota m_φ=Σm_ν·H_alg:

    portal Yukawa:  m_φ ≈ (y/√(8π²))·m_ν·√(log(M_Pl/m_ν))

Pregunta nítida de este ataque:
  (A) ¿el acoplamiento y sale LIMPIO de (φ,π)?  Si sí → m_φ físico, sin offset.
  (B) ¿qué escala de m_φ implica un y limpio O(1)?  (esperado: sub-eV/meV)
  (C) re-prueba R=1/[3(KAL−φ)] y escaneo de identidad exacta para Σm_ν.

Criterio referee-hostil: limpio = error < 0.2% Y exponente/factor pequeño.
"""
import numpy as np
import itertools
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import PHI, PI, OMEGA, BETA, KAL0, MIRA, AURA, T_R, M_V

# ── escalas físicas (eV) ────────────────────────────────────────────────
M_PL       = 1.220910e28
M_NU_HEAVY = 0.05            # neutrino más pesado (osc. atmosférica)
SUM_MNU    = 0.0824          # CONTAMINADO OP-14
R_OP14     = 4*KAL0 - 22     # 0.0856, la forma frágil actual

print("=" * 78)
print(" OP-14 — ATAQUE NUEVO: portal Yukawa  (no el offset 22)")
print("=" * 78)
print(f"  φ={PHI:.5f}  π={PI:.5f}  Ω={OMEGA:.5f}  β={BETA:.5f}")
print(f"  KAL₀={KAL0:.5f}  MIRA={MIRA:.5f}  AURA={AURA:.5f}")
print(f"  T_R={T_R:.5f}  M_V={M_V:.5f}")
print(f"  R actual (4·KAL−22) = {R_OP14:.6f}  ← forma frágil")

# ════════════════════════════════════════════════════════════════════════
# (A) ¿qué y necesita cada m_φ target?  y ¿hay constante SSEE limpia cerca?
# ════════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" (A) Portal Yukawa:  m_φ = (y/√(8π²))·m_ν·√(log(M_Pl/m_ν))")
print("─"*78)
pref = lambda mnu: np.sqrt(np.log(M_PL/mnu)) / np.sqrt(8*PI**2)  # m_φ = y·pref·m_ν
for mnu_name, mnu in [("m_ν=0.05 (heavy)", M_NU_HEAVY), ("Σm_ν=0.0824", SUM_MNU)]:
    k = pref(mnu)*mnu     # m_φ = y·k
    print(f"  {mnu_name}:  m_φ = y·{k*1e3:.4f} meV")
    for target_name, target in [("meV (=k, y=1)", k), ("0.1 eV", 0.1),
                                 ("1 eV", 1.0), ("5.6 eV (P6)", 5.6)]:
        y_need = target / k
        print(f"     m_φ={target_name:<16} → y = {y_need:.3f}")

# inventario de y limpios O(1) y qué m_φ dan (con m_ν=0.05)
print("\n  Si y es una constante SSEE limpia O(1), m_φ resultante (m_ν=0.05):")
k0 = pref(M_NU_HEAVY)*M_NU_HEAVY
y_clean = [("1", 1.0), ("1/AURA", 1/AURA), ("AURA", AURA), ("φ", PHI),
           ("π", PI), ("Ω", OMEGA), ("β", BETA), ("KAL₀", KAL0),
           ("MIRA", MIRA), ("φ²", PHI**2), ("Ω·φ", OMEGA*PHI)]
for yn, y in y_clean:
    m = y*k0
    flag = " ★ eV!" if 1 <= m <= 30 else (" ·sub-eV" if 0.1 <= m < 1 else "")
    print(f"     y={yn:<8} = {y:6.3f} → m_φ = {m*1e3:8.3f} meV{flag}")

# ════════════════════════════════════════════════════════════════════════
# (B) ¿hay un y limpio que cierre el lazo a un m_φ que SEA constante×escala?
#     (el golazo: m_φ derivado dos veces de forma consistente)
# ════════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" (B) ¿qué y daría m_φ = 5.6 eV, y es ese y una identidad (φ,π)?")
print("─"*78)
y_56 = 5.6 / k0
print(f"  y necesario para m_φ=5.6 eV (m_ν=0.05) = {y_56:.4f}")
# escaneo de identidad para ese y
basis = {'φ':PHI,'π':PI,'Ω':OMEGA,'β':BETA,'KAL':KAL0,'MIRA':MIRA,
         'AURA':AURA,'T_R':T_R,'M_V':M_V}
best = []
names = list(basis); vals = np.array([basis[n] for n in names])
for combo in itertools.product(range(-2,3), repeat=len(names)):
    if sum(abs(c) for c in combo) > 4 or all(c==0 for c in combo):
        continue
    val = np.prod([vals[i]**combo[i] for i in range(len(names))])
    if val <= 0:
        continue
    err = abs(val - y_56)/y_56
    if err < 0.01:
        expr = "·".join(f"{names[i]}^{combo[i]}" for i in range(len(names)) if combo[i])
        best.append((err, expr, val))
best.sort()
if best:
    print(f"  candidatos a y={y_56:.3f} con error <1%:")
    for err, expr, val in best[:8]:
        print(f"     {expr:<28} = {val:.4f}  ({err*100:+.3f}%)")
else:
    print(f"  NINGÚN monomio (φ,π) hasta grado 4 cae a <1% de y={y_56:.3f}.")
    print(f"  ⟹ el y de 5.6 eV NO es una constante SSEE limpia.")

# ════════════════════════════════════════════════════════════════════════
# (C) re-prueba R y escaneo de identidad exacta para Σm_ν / R
# ════════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" (C) R = 4·KAL−22 vs candidatos limpios (objetivo: identidad, no aprox)")
print("─"*78)
R_true = R_OP14
cand_R = [("4·KAL−22", 4*KAL0-22), ("1/[3(KAL−φ)]", 1/(3*(KAL0-PHI))),
          ("2/[3(3π−φ)]", 2/(3*(3*PI-PHI))), ("1/(Ω·β)", 1/(OMEGA*BETA)),
          ("φ⁻⁴·... ", PHI**-4*0.5), ("(π−φ)/Ω³", (PI-PHI)/OMEGA**3)]
for cn, cv in cand_R:
    err = abs(cv-R_true)/R_true
    flag = " ★ identidad?" if err < 0.001 else ""
    print(f"     {cn:<16} = {cv:.6f}  (Δ={err*100:+.3f}%){flag}")
# escaneo limpio para R
print(f"\n  Escaneo monomios (φ,π) → ¿alguno = R={R_true:.6f} a <0.2%?")
bestR = []
for combo in itertools.product(range(-3,4), repeat=len(names)):
    if sum(abs(c) for c in combo) > 5 or all(c==0 for c in combo):
        continue
    val = np.prod([vals[i]**combo[i] for i in range(len(names))])
    if val <= 0:
        continue
    err = abs(val - R_true)/R_true
    if err < 0.002:
        expr = "·".join(f"{names[i]}^{combo[i]}" for i in range(len(names)) if combo[i])
        bestR.append((err, expr, val))
bestR.sort()
if bestR:
    for err, expr, val in bestR[:6]:
        print(f"     {expr:<30} = {val:.6f}  ({err*100:+.3f}%)")
else:
    print(f"     NINGUNO a <0.2%. R no es identidad (φ,π) limpia.")

print("\n" + "=" * 78)
print(" VEREDICTO")
print("=" * 78)
print(f"""
  (A) Yukawa con y limpio O(1) → m_φ en decenas de meV. eV exige y~{y_56:.0f}.
  (B) ese y~{y_56:.0f} {'ES' if best else 'NO es'} una constante (φ,π) limpia.
  (C) R=4KAL−22 {'tiene' if bestR else 'no tiene'} identidad (φ,π) a <0.2%.

  Lectura honesta: si (B) y (C) salen NO, OP-14 sigue sin ataque limpio,
  y la física (Yukawa) vota por m_φ sub-eV/meV — NO por 5.6 eV.
  5.6 eV queda en su caja de Schrödinger hasta DESI/Euclid (k_fs).
""")
