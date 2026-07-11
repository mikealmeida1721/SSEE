"""
OP-11 — ¿qué ξ (acoplamiento no-mínimo ξRχ²) arregla la sobreproducción?
=========================================================================
Continúa op9_gravitational_production.py, que predijo Ω_χh²=4140 (×34 500
de más) SIN parámetros libres (m_φ=46.29 meV Yukawa + r=φ⁻¹⁰).

El término ξRχ² YA está en el Lagrangiano de Paper 6 (L428-436). Durante
inflación R=12H_inf² es enorme, así que la masa efectiva del campo χ se
infla:
        m_eff² = m_φ² + 12 ξ H_inf²
Si ξ~O(1)..O(100), m_eff >> H_inf y χ deja de ser ultraligero: NO crece
linealmente (de Sitter), sino que llega a EQUILIBRIO de Starobinsky-Yokoyama:
        ⟨χ²⟩_eq = 3 H_inf⁴ / (8π² m_eff²)   (en vez de (H/2π)²·N)

PREGUNTA HONESTA (regla del proyecto):
  Fijar ξ pidiendo Ω=0.12 es un FIT (un observable fija un parámetro).
  SOLO es golazo si el ξ que sale es LIMPIO en (φ,π). Si no, es OP-11:
  un parámetro libre confinado, no derivado. Lo reportamos sin maquillar.
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

# ── inputs FIJOS (idénticos a op9_gravitational_production) ────────────────
M_PL  = 2.435e18
A_S   = 2.1e-9
N_EF  = 60.0
r_tensor = PHI**-10
H_inf = M_PL * np.sqrt((PI**2/2) * r_tensor * A_S)

OVERPROD = 34500.0   # factor de sobreproducción medido por op9 (Ω/0.12)

print("=" * 74)
print(" OP-11 — ξ que arregla la sobreproducción gravitacional de χ")
print("=" * 74)
print(f"  H_inf = {H_inf:.4e} GeV   (r=φ⁻¹⁰ fijo)")
print(f"  sobreproducción a corregir: ×{OVERPROD:.0f}")

# ── paso 1: factor de supresión de ⟨χ²⟩ que da el acoplamiento ξ ──────────
#   de Sitter (ξ=0):   ⟨χ²⟩ = (H/2π)²·N      = H²·N/(4π²)
#   equilibrio S-Y:    ⟨χ²⟩ = 3H⁴/(8π²m_eff²) con m_eff²≈12ξH²  (m_φ<<H)
#                            = 3H²/(8π²·12ξ)  = H²/(32π²ξ)
#   supresión = eq/deSitter = [1/(32π²ξ)]/[N/(4π²)] = 1/(8ξN)
#   Ω ∝ ⟨χ²⟩, así que pedimos  1/(8ξN) = 1/OVERPROD
print("\n── paso 1: ξ pedido por Ω=0.12 ──")
xi_needed = OVERPROD / (8.0 * N_EF)
print(f"  supresión necesaria = 1/{OVERPROD:.0f}")
print(f"  1/(8·ξ·N) = 1/{OVERPROD:.0f}  →  ξ = {OVERPROD:.0f}/(8·{N_EF:.0f}) = {xi_needed:.3f}")

# verificación de consistencia del régimen: ¿m_eff >> H_inf?
m_eff_over_H = np.sqrt(12*xi_needed)
print(f"  chequeo régimen: m_eff/H_inf = √(12ξ) = {m_eff_over_H:.2f}  "
      f"({'✓ pesado, equilibrio válido' if m_eff_over_H>3 else '⚠ no tan pesado'})")

# ── paso 2: ¿ξ≈{xi_needed:.1f} es LIMPIO en (φ,π)? ────────────────────────
print("\n── paso 2: ¿ξ es un monomio limpio (φ,π)? (referee-hostil <0.2%) ──")
basis = {'φ':PHI,'π':PI,'Ω':OMEGA,'β':BETA,'KAL':KAL0,'MIRA':MIRA,
         'AURA':AURA}
names = list(basis); vals = np.array([basis[n] for n in names])
hits = []; dens = 0
for combo in itertools.product(range(-3,4), repeat=len(names)):
    s = sum(abs(c) for c in combo)
    if s == 0 or s > 5:
        continue
    val = np.prod([vals[i]**combo[i] for i in range(len(names))])
    if val <= 0:
        continue
    err = abs(val - xi_needed)/xi_needed
    if err < 0.05:
        dens += 1
    if err < 0.002:
        expr = "·".join(f"{names[i]}^{combo[i]}" for i in range(len(names)) if combo[i])
        hits.append((err, expr, val, s))
hits.sort()
if hits:
    print(f"  candidatos a ξ={xi_needed:.3f} con error <0.2%:")
    for err, expr, val, s in hits[:8]:
        print(f"     {expr:<26} = {val:.4f}  ({err*100:+.3f}%, grado {s})")
    print(f"  densidad (#monomios en ±5%) = {dens}")
    if dens > 20:
        print(f"  ⚠ densidad ALTA → un match es esperado por combinatoria, NO señal")
else:
    print(f"  NINGÚN monomio (φ,π) grado≤5 cae a <0.2% de ξ={xi_needed:.3f}")
    print(f"  densidad (#monomios en ±5%) = {dens}")

# referencia: potencias de φ cercanas
print(f"\n  referencia potencias φ:  φ⁸={PHI**8:.3f}  φ⁹={PHI**9:.3f}  "
      f"3φ⁷={3*PHI**7:.3f}  Ω·KAL={OMEGA*KAL0:.3f}")

print("\n" + "=" * 74)
print(" LECTURA")
print("=" * 74)
print(f"""
  ξ ≈ {xi_needed:.1f} apaga la sobreproducción (de ×34 500 a 1).
  - Si arriba salió un monomio (φ,π) LIMPIO y de baja densidad → OP-9+OP-11
    cierran juntos sin parámetro libre: χ se puebla y da Ω=0.12 con (m_φ,ξ)
    ambos derivados. Golazo real.
  - Si NO salió limpio → ξ es un PARÁMETRO LIBRE fiteado a Ω=0.12. Honesto:
    SSEE predice la sobreproducción cualitativa (escalar ligero + inflación
    alta) pero el valor exacto de Ω lo fija un número a mano. Eso es OP-11,
    no resuelto. Cero maquillaje.
""")
