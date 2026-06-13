"""
⚠️  PRE-CANÓNICO / INVESTIGACIÓN CERRADA — OP-10 se cerró 2026-06-11 (7 rutas
    descartadas: F1-F4, OpA, 2c, 2e). Los valores m_φ=5.60 eV usados aquí son
    la ruta vieja (RETIRADA). Conservado como registro de exploración fallida.

OP-10 Fase 2 — Familia 2: análisis de consistencia jerárquica
================================================================

Potencial:
    V(φ) = V₀ tanh²(φ/f) + (m_φ² f²/2) sech²(φ/f)

Análisis cualitativo (NO requiere integración FRW completa):

CONSTRAINT 1 — Curvatura en mínimo:
    V''(0) = (2V₀ - m_φ²f²)/f² = m_φ²
    →  f · m_φ = √V₀  ≈ 8.25 μeV²

CONSTRAINT 2 — DE observada (w₀ = -0.840):
    Slow-roll requiere f ~ M_pl para que φ ruede lento.
    Sin slow-roll, el campo está dominado por cinética → w → +1 (NO -0.84).

Si la jerarquía de escalas (f vs M_pl, m_φ vs H₀) es incompatible, Familia 2
queda DESCARTADA antes de necesitar integración numérica.
"""
import numpy as np

# Constantes
PHI          = (1 + np.sqrt(5)) / 2
# ρ_crit canónico SSEE: 8.10e-11 eV⁴ es el coeficiente pelado (h=1); la densidad
# física se escala por h_MIRA² con h_MIRA=H_MIRA/100=0.67037 (canónico 2026-06-09).
H_MIRA       = 67.037                    # km/s/Mpc
h_MIRA       = H_MIRA / 100.0            # 0.67037
RHO_CRIT_eV4 = 8.10e-11 * h_MIRA**2     # ρ_crit en eV⁴ (SSEE, ×h²)
V0           = 0.840 * RHO_CRIT_eV4
M_PL_eV      = 2.435e27
H0_eV        = 1.45e-33
M_PHI_TARGET = 5.60   # eV
M_UV         = (5 * PHI**8 * RHO_CRIT_eV4)**0.25   # ≈ 9.62 meV (Paper 10)

print("=" * 76)
print("OP-10 FASE 2 — Familia 2: consistency check jerárquico")
print("=" * 76)

# ============================================================
# Punto 1 — derivar f desde constraint 1
# ============================================================
f_required = np.sqrt(V0) / M_PHI_TARGET
print(f"\n[1] CONSTRAINT 1 (curvatura mínimo = m_φ²):")
print(f"    f · m_φ = √V₀  →  f = √V₀/m_φ = {f_required:.3e} eV  ≈ {f_required*1e6:.2f} μeV")

# ============================================================
# Punto 2 — slow-roll necesita f ~ M_pl
# ============================================================
print(f"\n[2] SLOW-ROLL CONSISTENCY:")
print(f"    Para α-attractor con plateau DE, f debe ser comparable a M_pl.")
print(f"    M_pl              = {M_PL_eV:.3e} eV")
print(f"    f requerido       = {f_required:.3e} eV")
print(f"    Ratio f/M_pl      = {f_required/M_PL_eV:.3e}")
print(f"    → f es {M_PL_eV/f_required:.1e}× MÁS PEQUEÑO que M_pl.")

# ============================================================
# Punto 3 — velocidad de φ en slow-roll con este f
# ============================================================
# En slow-roll: φ̇ ≈ -V'/(3H). Para φ ≈ φ_today >> f:
#   V'(φ_today) ≈ (V₀/f) · tanh(x_today)·sech²(x_today) ≈ (V₀/f)·η
#   con η ≈ 2 e^(-2x_today) para x_today >> 1
#
# Supongamos x_today = 3 (φ_today = 3f):
x_today = 3.0
tanh_x  = np.tanh(x_today)
sech2_x = (1/np.cosh(x_today))**2
eta     = tanh_x * sech2_x
V_prime = (V0 / f_required) * eta
phi_dot = -V_prime / (3 * H0_eV)
phi_dot_sq = phi_dot**2

print(f"\n[3] CINÉTICA DEL CAMPO (slow-roll con f = {f_required:.2e} eV):")
print(f"    Asumiendo x_today = φ/f = 3:  η = tanh(3)·sech²(3) = {eta:.4f}")
print(f"    V'(φ_today)       ≈ V₀·η/f = {V_prime:.3e} eV³")
print(f"    φ̇ slow-roll        ≈ -V'/(3H₀) = {phi_dot:.3e} eV²")
print(f"    φ̇²                 = {phi_dot_sq:.3e} eV⁴")
print(f"    V₀                 = {V0:.3e} eV⁴")
print(f"    Ratio φ̇²/V₀        = {phi_dot_sq/V0:.3e}")

# ============================================================
# Punto 4 — w₀ implicada
# ============================================================
# w = (½φ̇² - V) / (½φ̇² + V)
# Para slow-roll válido: φ̇² << V → w ≈ -1
# Si φ̇² >> V → w → +1 (cinética domina)
KE = 0.5 * phi_dot_sq
PE_approx = V0  # asumiendo plateau, V≈V₀
w_implied = (KE - PE_approx) / (KE + PE_approx)
print(f"\n[4] w_DE IMPLICADA por esta dinámica:")
print(f"    KE = ½φ̇²           = {KE:.3e} eV⁴")
print(f"    PE ≈ V₀            = {PE_approx:.3e} eV⁴")
print(f"    w = (KE-PE)/(KE+PE) = {w_implied:.4f}")
print(f"    Target w₀          = -0.840")

# ============================================================
# Veredicto
# ============================================================
print("\n" + "=" * 76)
print("VEREDICTO FAMILIA 2")
print("=" * 76)

if phi_dot_sq > 10 * V0:
    print("""
  ✗ DESCARTADA: violación de jerarquía.

  El f tan pequeño (≈ 1.47 μeV) requerido por constraint 1 (curvatura mínimo
  = m_φ² = 31 eV²) hace que V'/f sea enorme. En slow-roll, esto produce un
  campo con φ̇² >> V₀, lo cual viola el régimen DE (w → +1 en vez de -0.84).

  Conclusión estructural: Familia 2 NO puede simultáneamente:
    (a) reproducir m_φ = 5.60 eV (requiere f ≈ 1 μeV)
    (b) reproducir w₀ = -0.840 (requiere f ~ M_pl ≈ 10²⁷ eV)

  La discrepancia es de 33 órdenes de magnitud en f.

  CAUSA RAÍZ: en modelos α-attractor canónicamente normalizados, la masa
  del campo en el mínimo escala como m_φ ~ √V₀/f. Para m_φ ~ eV con V₀ ~
  ρ_crit, f debe ser microscópico. Pero un f microscópico hace que el
  potencial sea "muy abrupto" (gradiente alto), violando slow-roll.

  IMPLICACIÓN PARA OP-10:
    - Familia 2 (α-attractor tanh² + binding) DESCARTADA.
    - Familia 3 (axion-like + exp) y otros enfoques necesitan exploración.
    - El problema fundamental: una sola componente φ canónica NO puede
      hospedar dos escalas tan separadas (Λ_DE^(1/4) ~ meV vs m_φ ~ eV).
""")
elif abs(w_implied - (-0.840)) < 0.05:
    print(f"""
  ✓ FAMILIA 2 VIABLE: w₀ ≈ {w_implied:.3f} compatible con target -0.840.
    Proceder a Fase 3 (verificar V'' completo y m_φ exacto).
""")
else:
    print(f"""
  ⚠ AMBIGUO: w_implied = {w_implied:.4f} (no match limpio).
    Investigar con integración FRW completa antes de descartar.
""")

# ============================================================
# Punto 5 — ¿hay rescate? análisis de dos campos
# ============================================================
print("=" * 76)
print("¿RESCATE POSIBLE? — separación de escalas DE/DM")
print("=" * 76)
print("""
  La raíz del problema es asumir que UN solo campo φ porta TANTO la DE
  como la DM. Las escalas de energía son muy distintas:

    V₀^(1/4) (DE)    ≈ {:.2e} eV (≈ meV, escala de Λ)
    m_φ      (DM)    ≈ 5.60 eV (≈ eV, escala neutrino)

  Razón m_φ / V₀^(1/4) ≈ {:.0f}  → 3-4 órdenes de magnitud.

  POSIBLES RESCATES (para futuras fases):

  A) DOS CAMPOS desacoplados: φ_DE rueda lento en V_DE = V₀ exp(-λφ_DE),
     mientras χ_DM oscila en V_DM = ½m²χ² con m = 5.60 eV. Esto es lo que
     Paper 6 actualmente asume (two-sector). Pero entonces seguimos con
     OP-9 (origen de m) sin resolver.

  B) ACOPLAMIENTO NO-MINIMAL ξφ²R: el acoplamiento al Ricci puede
     generar masa efectiva m_φ² = ξR ≈ ξ H₀². Para m_φ = 5.60 eV con
     H₀ ≈ 10⁻³³ eV, ξ ≈ (5.60/H₀)² ≈ 10⁶⁷. ξ astronómico — descartado.

  C) MECANISMO RADIATIVO: m_φ generada por loops electroweak (Coleman-
     Weinberg). Requiere portal SM → contradice premisa SSEE de φ-DM
     sin portal. Conflicto con Paper 6 §"Particle-physics constraints".

  D) AXION-LIKE (Familia 3): m_φ²f² ~ Λ_QCD-like scale. Pero Familia 3
     da φ-DM oscilante, no estático coherente como P6 asume.

  E) ABANDONAR derivación de m_φ desde V(φ): aceptar OP-9 como ansatz
     permanente. SSEE permanece minimal-parameter (3-4 efectivos) pero
     NO cero-parámetros. Esto es lo más honesto pragmáticamente.

  RECOMENDACIÓN: explorar Familia 3 (Fase 2b) antes de aceptar (E).
""".format(V0**0.25, M_PHI_TARGET/V0**0.25))
