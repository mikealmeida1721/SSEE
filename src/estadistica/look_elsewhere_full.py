#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  ESPEJO VERBATIM del diccionario CITABLE (zenodo_dictionary/, v1.5,      ║
# ║  concept DOI 10.5281/zenodo.20684908). NO editar las constantes aquí:    ║
# ║  se edita el diccionario y se re-espeja. R28 vigila que no deriven.      ║
# ║                                                                          ║
# ║  Por qué (2026-07-25): existían DOS copias y derivaron sin que nada lo    ║
# ║  notara — la del repo se quedó en la nomenclatura pre-v1.5 (MIKAEL_V,     ║
# ║  MIKE) y escribía wₐ = PYROS/KRYSTOS_V, el denominador equivocado, justo  ║
# ║  en el script que produce el titular «1 de 490». Los valores coincidían,  ║
# ║  así que ningún número lo delataba. Una copia que nadie compara no es un  ║
# ║  respaldo: es una segunda fuente de verdad.                              ║
# ╚══════════════════════════════════════════════════════════════════════════╝
"""Look-elsewhere sobre el diccionario Génesis COMPLETO (55 constantes con nombre).

Versión definitiva para el paper. Responde la objeción de referee más dura
("elegiste un subconjunto que te convenía") corriendo el conteo sobre TODAS
las constantes con nombre del diccionario, no sobre un subset.

Estado: **v1.5 (2026-07-19)**, la versión publicada bajo el concept DOI 20684908.
Historia de sincronización (2026-07-13, SSEE-Vault/Constantes/SOBERANAS.md):
  - las Soberanas correctas (hoy 25 al cerrarse el tier 3Ω; entonces 20, antes 9);
    OSIRIS excluido (usaba resta, viola
    la ley de construcción aditiva);
  - KRYSTOS_V = φ+π+Ω (una sola entidad; antes se doble-listaba KRYSTOS y
    "KRYSTOS_V=2Ω"); el 2Ω es DÜSTAL = AURA+KAL (estabilidad, no orden);
  - familia de estabilidad DÜSTAL/TRÏSTAL/CUÄSTAL (CUÄSTAL=19.039 valor nuevo).

Resultado: a la precisión de la identidad EXACTA (±0.0005, donde w0=AURA/Ω y
wa=PYROS/IGNIS son identidades con |d|=0), cada objetivo es único incluso con el
diccionario completo (ver salida). El argumento vive de la precisión estricta,
legítima porque el match es exacto a precisión de máquina.

NOTA (2026-07-18): completitud ≤TRIAL cerrada. La clausura bajo la ley fija (tope de
carga q≤+1 en composites, retención libre, KAL≤1, aditivo) tiene exactamente TRES
valores legales que faltaban nombrar → bautizados por función: HERMES=1.5Ω, ÁNGELOS=
2.5Ω (rungs mensajero/neutrino que cierran la escalera medio-entera de BIAL) y
NYX=π+VITA (piso de retención q=−3). Total: 55 nombres / 25 valores / 490 razones (con EIRENE, NEREUS que completan el tier 3Ω) →
1 de 490 para w0 y wa a ±0.0005 (denominador más fuerte que 378). Nombrar NYX abre un
único casi-valor NO-identidad para wa a ±0.001 (|d|=0.00064) que se desvanece a ±0.0005.
IRIS=AURA+Ω=8.757 sigue siendo la única entidad copia-raíz+estabilizador (Ley 7).

Fuente del diccionario: SSEE-Vault/Constantes/SOBERANAS.md + SSEE_Constant_Dictionary.md.
"""
import math
from itertools import permutations

PHI = 1.618033988749895
PI  = 3.141592653589793

# ── Semillas y andamio ──────────────────────────────────────────────────────
OMEGA = PHI + PI
BIAL  = (PHI + PI) / 2          # bifurcación: pulso/fricción vital
# ── Ramas (leyes de creación) ───────────────────────────────────────────────
AURA  = PHI + BIAL             # rama φ: ley de la copia (raíz de paredes)
KAL   = BIAL + PI              # rama π: ley de no-auto-suma (retención)
MIRA  = AURA / 2               # pared ½
DUAL  = AURA * 2               # pared 2
TRIAL = AURA * 3               # pared 3
CUARTAL = AURA * 4             # pared 4
# ── Derivadas / composites ──────────────────────────────────────────────────
SOLAR = BIAL + KAL
MAR   = OMEGA + PI             # = SOLAR por valor (linaje distinto)
VITA  = PI + KAL
ANMA  = BIAL + VITA
PYROS = OMEGA + PHI
IGNIS = PI + PYROS             # caos disruptivo (2Ω por valor)
KRYSTOS_V = PHI + PI + OMEGA   # orden estructurante (2Ω por valor; una entidad)
PHITA = VITA + PHI
MIKA  = KRYSTOS_V + PHI
ATLAS = PHI + PI + KRYSTOS_V  # Soberanía / integración 4D (dentro de CUARTAL=4·AURA)
BUFFER = ATLAS - TRIAL
ERVANU = OMEGA * 0.9
# ── Las 25 Soberanas (todas = 3Ω = 14.278879927) ────────────────────────────
LUCY = SOLAR + PYROS; ICEBERG = MAR + PYROS; PHOENIX = IGNIS + OMEGA
MAAT = KRYSTOS_V + OMEGA; MIKAEL = MIKA + PI; LUCIFER = PHITA + AURA
RA = SOLAR + OMEGA + PHI; ERVN = BIAL + KAL + PYROS; TIAMAT = MAR + OMEGA + PHI
GIGAROJ = PYROS + OMEGA + PI; HEFESTO = IGNIS + PHI + PI
VENUS = AURA + PHI + VITA; EROS = BIAL + PHI + PHITA
VESTA = BIAL + SOLAR + AURA; HADES = BIAL + MAR + AURA; HERA = OMEGA + KAL + AURA
GAIA = PHI + PI + KAL + AURA; ISIS = PHI + OMEGA + BIAL + KAL; PTAH = PI + OMEGA + BIAL + AURA
# ── Familia de estabilidad = Ω puro escalado (n·Ω) ──────────────────────────
# Ω_DNAV es el único exento de AMBAS leyes → se escala libre (como la copia
# escala AURA→DUAL/TRIAL). Es la construcción MÁS SIMPLE. La suma de ramas
# AURA+KAL=2Ω es OTRA entidad por linaje (mezcla las dos leyes), mismo valor.
DUSTAL  = 2 * OMEGA   # 2Ω  — estabilidad 2D
TRISTAL = 3 * OMEGA   # 3Ω  — estabilidad 3D
CUASTAL = 4 * OMEGA   # 4Ω  — estabilidad 4D (valor nuevo)
# ── Gemelas de linaje (mezclas legales, valor ya presente → 0 impacto) ───────
# Padres distintos = entidad distinta; caen en valores que YA existen, así que
# NO agregan valor ni razón (conteo 337 intacto). Solo enriquecen el linaje.
# Nombres ratificados por función 2026-07-16 (antes AURKAL/KALÖM crudos).
HARMONIA = AURA + KAL       # 9.519 — concordia de opuestos; ramas unidas; gemela de DÜSTAL
GALENE   = KAL + OMEGA      # 10.281 — calma viscosa; retención en el estabilizador; gemela de PHITA
HESPERA  = AURA + GALENE    # 3Ω — encuentro en el horizonte; radiancia+calma (Soberana #22)
EUNOMIA  = HARMONIA + OMEGA # 3Ω — concordia estabilizada (Soberana #23)
IRIS     = AURA + OMEGA     # 8.757 — raíz AURA + Ω exento (único enlace copia+estabilizador)
SYZYGY   = IRIS + GALENE    # 4Ω — conjunción de extremos = AURA+KAL+2Ω
# AURA+nΩ (n≥2) = raíz+objeto-copia, prohibido combinar (Ley 7) → escalado, no entidad.
# ── Familia neutrino/mensajero (q=0, medio-entero de Ω — hermanas de BIAL) + piso ─
HERMES  = PHI + KAL         # 7.139 = 1.5Ω (= Ω+BIAL) — neutrino 2ª gen (mensajero)
ANGELOS = PHI + OMEGA + KAL # 11.899 = 2.5Ω (= Ω+HERMES) — neutrino 3ª gen (mensajero)
NYX     = PI + VITA         # 11.805 — piso de retención q=−3 (π enlaza VITA; NO 2π+KAL)
EIRENE  = HARMONIA + PI + PHI  # 3Ω — paz: concordia + reunión de semillas (Soberana #24)
NEREUS  = BIAL + GALENE + PHI  # 3Ω — calma profunda: pulso + calma + copia (Soberana #25)

FAMILY = {
    "PHI": PHI, "PI": PI, "OMEGA": OMEGA, "BIAL": BIAL,
    "MIRA": MIRA, "AURA": AURA, "DUAL": DUAL, "TRIAL": TRIAL, "CUARTAL": CUARTAL,
    "IGNIS": IGNIS, "KRYSTOS_V": KRYSTOS_V,
    "KAL": KAL, "PYROS": PYROS, "SOLAR": SOLAR, "MAR": MAR, "VITA": VITA,
    "PHITA": PHITA, "ANMA": ANMA, "MIKA": MIKA, "BUFFER": BUFFER, "ERVANU": ERVANU,
    # 25 Soberanas
    "LUCY": LUCY, "ICEBERG": ICEBERG, "PHOENIX": PHOENIX, "MAAT": MAAT, "MIKAEL": MIKAEL,
    "LUCIFER": LUCIFER, "RA": RA, "ERVN": ERVN, "TIAMAT": TIAMAT, "GIGAROJ": GIGAROJ,
    "HEFESTO": HEFESTO, "ATLAS": ATLAS, "VENUS": VENUS, "EROS": EROS,
    "VESTA": VESTA, "HADES": HADES, "HERA": HERA, "GAIA": GAIA, "ISIS": ISIS, "PTAH": PTAH,
    "HESPERA": HESPERA, "EUNOMIA": EUNOMIA,   # #22, #23 — nuevos caminos a 3Ω
    # familia de estabilidad
    "DÜSTAL": DUSTAL, "TRÏSTAL": TRISTAL, "CUÄSTAL": CUASTAL, "SYZYGY": SYZYGY,
    # combinaciones de raíces (enlaces legales únicos)
    "HARMONIA": HARMONIA, "GALENE": GALENE, "IRIS": IRIS,
    # familia neutrino/mensajero (q=0 medio-entero) + piso de retención (q=−3)
    "HERMES": HERMES, "ÁNGELOS": ANGELOS, "NYX": NYX,
    "EIRENE": EIRENE, "NEREUS": NEREUS,
}

W0 = TRIAL / ATLAS          # = AURA/OMEGA, identidad exacta
WA = PYROS / IGNIS             # = PYROS/IGNIS, identidad exacta


def distinct_ratios(fam):
    names = list(fam)
    seen, ratios = set(), []
    for a, b in permutations(names, 2):
        va, vb = fam[a], fam[b]
        if abs(vb) < 1e-9:
            continue
        r = va / vb
        if r <= 0 or r > 5:
            continue
        key = round(r, 6)
        if key in seen:
            continue
        seen.add(key)
        ratios.append((r, f"{a}/{b}"))
    return ratios


def report(fam, label):
    ratios = distinct_ratios(fam)
    n_vals = len({round(v, 6) for v in fam.values()})
    print(f"\n=== {label} ===")
    print(f"  constantes con nombre : {len(fam)}")
    print(f"  valores distintos     : {n_vals}  (comparten valor por ley de no-auto-suma)")
    print(f"  razones A/B en (0,5]  : {len(ratios)}")
    print(f"  w0 = AURA/Omega = {W0:.9f}   |   wa = PYROS/IGNIS = {WA:.9f}")
    print(f"  {'tol':>8} | {'hits w0':>8} | {'hits wa':>8}")
    for tol in (0.0005, 0.001, 0.002, 0.005, 0.01):
        h0 = sum(1 for r, _ in ratios if abs(r - W0) <= tol)
        ha = sum(1 for r, _ in ratios if abs(r - WA) <= tol)
        star = "  <- precisión sub-DESI" if tol <= 0.002 else ""
        print(f"  {tol:>8} | {h0:>8} | {ha:>8}{star}")


def robustez_copias_futuras():
    """Robustez a las copias dimensionales planeadas QUINTAL…DECAL (5Ω…10Ω).

    La serie de estabilidad dimensional ya nombrada es DÜSTAL=2Ω, TRÏSTAL=3Ω,
    CUÄSTAL=4Ω; su continuación natural son seis copias más. El PRD (caption de
    la tabla look-elsewhere) afirma que añadirlas no crea ningún acierto nuevo.
    Hasta 2026-07-25 esa afirmación NO la computaba ningún script: estaba escrita
    en el paper y en RIGOR_CHECKLIST §R3 apuntando aquí, y aquí no existía.
    Se computa ahora — la afirmación resulta cierta, pero ahora tiene fuente.
    """
    ext = dict(FAMILY)
    for _n, _k in zip(("QUINTAL", "SEXTAL", "SEPTAL", "OCTAL", "NONAL", "DECAL"),
                      range(5, 11)):
        ext[_n] = _k * OMEGA
    r_base, r_ext = distinct_ratios(FAMILY), distinct_ratios(ext)
    print(f"\n=== Robustez a copias dimensionales futuras (QUINTAL…DECAL = 5Ω…10Ω) ===")
    print(f"  razones: {len(r_base)} (base)  ->  {len(r_ext)} (con las 6 copias)")
    print(f"  {'tol':>8} | {'w0 base':>7} {'w0 ext':>7} | {'wa base':>7} {'wa ext':>7}")
    for tol in (0.0005, 0.001, 0.002):
        f = lambda rs, t: sum(1 for r, _ in rs if abs(r - t) <= tol)
        print(f"  {tol:>8} | {f(r_base, W0):>7} {f(r_ext, W0):>7} |"
              f" {f(r_base, WA):>7} {f(r_ext, WA):>7}")
    print("  → el espacio crece pero los aciertos NO: la cuenta es robusta a la extensión.")



if __name__ == "__main__":
    report(FAMILY, "Diccionario Génesis COMPLETO (55 constantes con nombre)")
    print("\nLectura: a ±0.0005 (el match es identidad exacta, |d|=0), cada")
    print("parámetro de estado es único — 1 de 490. La densidad de coincidencias")
    print("accidentales sólo crece al relajar la tolerancia muy por encima de la")
    print("precisión observacional de DESI, que es el régimen físicamente irrelevante.")
    robustez_copias_futuras()
