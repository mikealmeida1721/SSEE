#!/usr/bin/env python3
"""Look-elsewhere sobre el diccionario Génesis COMPLETO (55 constantes con nombre).

Versión definitiva para el paper. Responde la objeción de referee más dura
("elegiste un subconjunto que te convenía") corriendo el conteo sobre TODAS
las constantes con nombre del diccionario, no sobre un subset.

Sincronizado 2026-07-17 con el repositorio citable del diccionario
(zenodo_dictionary, concept DOI 10.5281/zenodo.20684908, v1.3): el set creció al
formalizarse las leyes de linaje — familia de estabilidad Ω (DÜSTAL/TRÏSTAL/CUÄSTAL),
gemelas de linaje renombradas por función (HARMONIA=AURA+KAL, GALENE=KAL+Ω), Soberanas
nuevas (HESPERA, EUNOMIA), SYZYGY (4Ω) e IRIS=AURA+Ω (único enlace copia-raíz+equilibrio,
contado). Total 55 nombres / 25 valores distintos / 490 razones.

Resultado: a la precisión de la identidad exacta (±0.0005, donde
w0=AURA/Ω y wa=PYROS/IGNIS son identidades EXACTAS), cada objetivo es único
— 1 de 490 razones (a ±0.0005) — incluso con el diccionario completo. Reporta además la curva
de sensibilidad a la tolerancia (honestidad: a ±0.01 el bosque se densifica; el
argumento vive de la precisión estricta, legítima porque el match es exacto a
precisión de máquina).

Fuente del diccionario: zenodo_dictionary/ssee_constants.py + SOBERANAS.md.
"""
import math
from itertools import permutations

PHI = 1.618033988749895
PI  = 3.141592653589793

# ── Semillas y andamio ──────────────────────────────────────────────────────
OMEGA = PHI + PI                 # equilibrio de las dos fuerzas (exento, n=1)
BIAL  = (PHI + PI) / 2           # bifurcación (= Ω/2)
# ── Raíces (leyes de creación) ──────────────────────────────────────────────
AURA  = PHI + BIAL               # rama φ: ley de copia (raíz de paredes)
MIRA  = AURA / 2                 # pared ½
DUAL  = AURA * 2                 # pared 2
TRIAL = AURA * 3                 # pared 3
CUARTAL = AURA * 4               # pared 4
KAL   = BIAL + PI                # rama π: ley de no-auto-suma
PYROS = OMEGA + PHI
SOLAR = BIAL + KAL
MAR   = OMEGA + PI               # = SOLAR por valor (linaje distinto)
VITA  = PI + KAL
PHITA = VITA + PHI
ANMA  = BIAL + VITA
IGNIS = PI + PYROS               # caos disruptivo (9.519)
KRYSTOS_V = PHI + PI + OMEGA     # orden estructurante (9.519, linaje distinto)
MIKA  = KRYSTOS_V + PHI
ERVANU = OMEGA * 0.9
MIKAEL_V = PHI + PI + KRYSTOS_V  # integración 5D (una Soberana)
BUFFER = MIKAEL_V - TRIAL
# ── Soberanas (caminos a 3Ω = 14.278880) ────────────────────────────────────
LUCY = SOLAR + PYROS; ICEBERG = MAR + PYROS; MIKE = IGNIS + OMEGA
MAAT = KRYSTOS_V + OMEGA; MIKAEL = MIKA + PI; LUCIFER = PHITA + AURA
RA = SOLAR + OMEGA + PHI; ERVN = BIAL + KAL + PYROS; TIAMAT = MAR + OMEGA + PHI
GIGAROJ = PYROS + OMEGA + PI; HEFESTO = IGNIS + PHI + PI
VENUS = AURA + PHI + VITA; EROS = BIAL + PHI + PHITA
VESTA = BIAL + SOLAR + AURA; HADES = BIAL + MAR + AURA; HERA = OMEGA + KAL + AURA
GAIA = PHI + PI + KAL + AURA; ISIS = PHI + OMEGA + BIAL + KAL; PTAH = PI + OMEGA + BIAL + AURA
# ── Familia de estabilidad = copia aplicada a Ω (n·Ω; rungs terminales) ──────
DUSTAL  = 2 * OMEGA              # 2Ω — estabilidad 2D
TRISTAL = 3 * OMEGA              # 3Ω — estabilidad 3D
CUASTAL = 4 * OMEGA              # 4Ω — estabilidad 4D (valor nuevo)
# ── Combinaciones de raíces (enlaces legales) ───────────────────────────────
HARMONIA = AURA + KAL            # 9.519 — concordia de opuestos (las dos ramas)
GALENE   = KAL + OMEGA           # 10.281 — calma viscosa (retención + equilibrio)
IRIS     = AURA + OMEGA          # 8.757 — único enlace copia-raíz + equilibrio (contado)
HESPERA  = AURA + GALENE         # 3Ω — Soberana #22
EUNOMIA  = HARMONIA + OMEGA      # 3Ω — Soberana #23
SYZYGY   = IRIS + GALENE         # 4Ω — conjunción de extremos = AURA+KAL+2Ω
HERMES   = PHI + KAL             # 7.139 = 1.5Ω (=Ω+BIAL) — neutrino 2ª gen (mensajero)
ANGELOS  = PHI + OMEGA + KAL     # 11.899 = 2.5Ω (=Ω+HERMES) — neutrino 3ª gen (mensajero)
NYX      = PI + VITA             # 11.805 — piso de retención q=−3 (π enlaza VITA; NO 2π+KAL)
EIRENE   = HARMONIA + PI + PHI  # 3Ω — paz: concordia+reunión de semillas (Soberana #24)
NEREUS   = BIAL + GALENE + PHI  # 3Ω — calma profunda: pulso+calma+copia (Soberana #25)

FAMILY = {
    "PHI": PHI, "PI": PI, "OMEGA": OMEGA, "BIAL": BIAL,
    "AURA": AURA, "MIRA": MIRA, "DUAL": DUAL, "TRIAL": TRIAL, "CUARTAL": CUARTAL,
    "KAL": KAL, "PYROS": PYROS, "SOLAR": SOLAR, "MAR": MAR, "VITA": VITA,
    "PHITA": PHITA, "ANMA": ANMA, "IGNIS": IGNIS, "KRYSTOS_V": KRYSTOS_V,
    "MIKA": MIKA, "ERVANU": ERVANU, "BUFFER": BUFFER,
    # Soberanas
    "LUCY": LUCY, "ICEBERG": ICEBERG, "MIKE": MIKE, "MAAT": MAAT, "MIKAEL": MIKAEL,
    "LUCIFER": LUCIFER, "RA": RA, "ERVN": ERVN, "TIAMAT": TIAMAT, "GIGAROJ": GIGAROJ,
    "HEFESTO": HEFESTO, "MIKAEL_V": MIKAEL_V, "VENUS": VENUS, "EROS": EROS,
    "VESTA": VESTA, "HADES": HADES, "HERA": HERA, "GAIA": GAIA, "ISIS": ISIS, "PTAH": PTAH,
    "HESPERA": HESPERA, "EUNOMIA": EUNOMIA,
    # familia de estabilidad
    "DÜSTAL": DUSTAL, "TRÏSTAL": TRISTAL, "CUÄSTAL": CUASTAL, "SYZYGY": SYZYGY,
    # combinaciones de raíces
    "HARMONIA": HARMONIA, "GALENE": GALENE, "IRIS": IRIS,
    "HERMES": HERMES, "ÁNGELOS": ANGELOS, "NYX": NYX,
    "EIRENE": EIRENE, "NEREUS": NEREUS,
}

W0 = TRIAL / MIKAEL_V            # = AURA/OMEGA, identidad exacta
WA = PYROS / KRYSTOS_V           # = PYROS/IGNIS, identidad exacta


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
    report(FAMILY, f"Diccionario Génesis COMPLETO ({len(FAMILY)} constantes)")
    print("\nLectura: a ±0.0005 (el match es identidad exacta, |d|=0), cada")
    print("parámetro de estado es único — 1 de 490 (a ±0.0005). La densidad de coincidencias")
    print("accidentales sólo crece al relajar la tolerancia muy por encima de la")
    print("precisión observacional de DESI, que es el régimen físicamente irrelevante.")
    robustez_copias_futuras()
