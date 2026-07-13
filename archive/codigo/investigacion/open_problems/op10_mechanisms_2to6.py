"""
OP-10 — Mecanismos 2, 3, 6 EJECUTADOS (los que el plan dejó en servilleta)
===========================================================================
op10_systematic_search.py NUNCA implementó estos tres de verdad:
  - su "M2" era inflación α-attractor (10^49 eV), NO el Yukawa neutrino-φDM
  - el axion-like (Λ²/f) con f derivado nunca se computó
  - el two-field nunca se cuantificó

Aquí los corremos con física real y acoplamientos DERIVADOS, no fitted.
No buscamos confirmar nada: buscamos el m_φ que salga, para después meterlo
al test de frío. Si sale meV, ese meV es un RESULTADO, no una rendición.

Inputs limpios: M_Pl, ℏH_MIRA, ρ_crit, ρ_DE, (φ,π) y derivados.
Σm_ν se usa SOLO marcado como contaminado (OP-14) para ver si la ruta eV
DEPENDE de él — eso es información, no trampa.
"""
import numpy as np
import sys
import os as _os
_r = _os.path.dirname(_os.path.abspath(__file__))
while _r != _os.path.dirname(_r) and not _os.path.isdir(_os.path.join(_r, 'src')):
    _r = _os.path.dirname(_r)
sys.path.insert(0, _os.path.join(_r, 'src'))  # portable: sube hasta la raíz del repo
from ssee_core import PHI, PI, OMEGA, BETA, KAL0, MIRA, AURA, T_R, M_V

# ── Escalas físicas limpias (eV) ───────────────────────────────────────
H_MIRA = 1.4306e-33          # ℏH_MIRA
M_PL   = 1.220910e28         # Planck (no reducida)
M_PL_R = 2.435e27            # Planck reducida
h_mira = 0.67037
RHO_CRIT = 8.10e-11 * h_mira**2
RHO_DE   = 0.840 * RHO_CRIT
rho_qr   = RHO_CRIT**0.25
M_CUTOFF = 8.81e-3           # M = cutoff cinético P10 (sector DE, NO m_φ)
SUM_MNU  = 0.0824            # CONTAMINADO OP-14 — solo para diagnóstico
M_NU_HEAVY = 0.05            # masa del neutrino más pesado (osc. atmosférica)

print("=" * 78)
print(" OP-10 MECANISMOS 2, 3, 6 — ejecutados con física real")
print("=" * 78)
print(f"  ρ_crit^(1/4) = {rho_qr*1e3:.3f} meV   |  M_Pl = {M_PL:.3e} eV")
print(f"  ℏH_MIRA = {H_MIRA:.3e} eV  |  Σm_ν = {SUM_MNU} eV (OP-14, contaminado)")
print(f"  TARGET informativo: eV-scale → si aparece, ver de qué depende")

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 2 — Yukawa neutrino-φDM:  L = y φ ν̄ν
#   masa escalar radiativa (Coleman-Weinberg, pieza log tras quitar cuadrática)
#   m_φ² ≈ (y²/8π²) · m_ν² · log(M_Pl/m_ν)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M2] Yukawa neutrino-φDM:  m_φ² = (y²/8π²)·m_ν²·log(Λ/m_ν)")
print("─"*78)
logfac = np.log(M_PL / M_NU_HEAVY)
print(f"  log(M_Pl/m_ν) = {logfac:.2f}")
# con y = O(1) natural
for yname, y in [("1 (natural)", 1.0), ("1/AURA", 1/AURA), ("φ", PHI),
                 ("4π (loop max)", 4*PI)]:
    m = np.sqrt(y**2 / (8*PI**2) * M_NU_HEAVY**2 * logfac)
    flag = " ★ eV!" if 1 <= m <= 30 else ""
    print(f"   y={yname:<14} → m_φ = {m*1e3:.4f} meV{flag}")
# ¿qué y se necesita para eV?
y_eV = np.sqrt(1.0**2 * 8*PI**2 / (M_NU_HEAVY**2 * logfac))   # m_φ=1 eV
y_56 = y_eV * np.sqrt(5.6)
print(f"  → para m_φ=1 eV se necesita y = {y_eV:.2f}")
print(f"  → para m_φ=5.6 eV se necesita y = {y_56:.2f}  (¿constante SSEE limpia? NO)")

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 3 — Axion-like:  V = Λ_ax⁴(1-cos(φ/f)),  m_φ = Λ_ax²/f
#   En régimen perturbativo f > Λ_ax  ⟹  m_φ < Λ_ax SIEMPRE.
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M3] Axion-like:  m_φ = Λ_ax² / f")
print("─"*78)
print(f"  Cota dura: si f > Λ_ax (perturbativo) ⟹ m_φ < Λ_ax.")
for Lname, Lax in [("M=8.81meV", M_CUTOFF), ("ρ_crit^¼", rho_qr),
                   ("ρ_DE^¼", RHO_DE**0.25)]:
    for fname, f in [("M_Pl", M_PL), ("√(Λ·M_Pl)", np.sqrt(Lax*M_PL)),
                     ("M_Pl_red", M_PL_R)]:
        m = Lax**2 / f
        print(f"   Λ_ax={Lname:<10} f={fname:<12} → m_φ = {m:.3e} eV")
    # f que daría eV (régimen anhárónico f<Λ, dudoso)
    f_eV = Lax**2 / 1.0
    print(f"   → para m_φ=1 eV con Λ_ax={Lname}: f={f_eV:.3e} eV (f<<Λ, no físico)")
print(f"  Veredicto: con Λ_ax ≤ meV, axion-like NUNCA llega a eV en régimen sano.")

# ═══════════════════════════════════════════════════════════════════════
# MECANISMO 6 — Two-field:  V_χ(χ) con escala UV propia
#   ¿hay alguna escala UV en SSEE que sea eV-natural sin invocar Σm_ν?
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "─"*78)
print(" [M6] Two-field: ¿existe escala UV eV-natural en SSEE sin Σm_ν?")
print("─"*78)
# Inventario de TODAS las escalas/combinaciones que dan algo cerca de eV
print("  Geometric means y combinaciones (buscando eV sin factor grande):")
combos = [
    ("√(M_Pl·ℏH)",        np.sqrt(M_PL*H_MIRA)),
    ("√(Σm_ν·M_Pl)",      np.sqrt(SUM_MNU*M_PL)),     # contaminado
    ("√(Σm_ν·ℏH)",        np.sqrt(SUM_MNU*H_MIRA)),   # contaminado
    ("(M_Pl·ρ_crit)^¼",   (M_PL*RHO_CRIT)**0.25),
    ("(M_Pl²ρ_crit)^⅙",   (M_PL**2*RHO_CRIT)**(1/6)),
    ("(ρ_crit·M_Pl²)^⅙",  (RHO_CRIT*M_PL**2)**(1/6)),
    ("Σm_ν (directo)",     SUM_MNU),                  # contaminado
    ("ρ_crit^¼",           rho_qr),
]
for name, m in combos:
    near_eV = " ← cerca de eV" if 0.3 <= m <= 30 else ""
    contam = " [CONTAMINADO Σm_ν]" if "Σm_ν" in name else ""
    print(f"   {name:<20} = {m:.4e} eV{near_eV}{contam}")

print("\n" + "=" * 78)
print(" DIAGNÓSTICO — ¿de qué depende la escala eV?")
print("=" * 78)
print("""
  M2 Yukawa : y~O(1) → meV. Para eV hace falta y~25-60 (no es constante SSEE).
  M3 axion  : con Λ_ax≤meV, m_φ<meV SIEMPRE en régimen sano. No llega a eV.
  M6 2-field: la ÚNICA escala SSEE que toca eV es Σm_ν (0.08 eV) — y está
              contaminada (OP-14, offset -22 ad hoc).

  PATRÓN FÍSICO (no rendición — diagnóstico):
    Todo mecanismo limpio aterriza en meV/μeV.
    La escala eV SOLO aparece anclada al sector neutrino (Σm_ν).
    ⟹ La ruta eV NO está muerta: está BLOQUEADA-POR-OP-14.
       Si Σm_ν se deriva limpio, el portal Yukawa/neutrino da eV legítimo.
       Esa es la siguiente vela concreta, no un punto final.
""")
