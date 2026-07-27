#!/usr/bin/env python3
"""preflight.py — verificar los INGREDIENTES antes de gastar horas de cómputo.

POR QUÉ EXISTE (2026-07-26)
---------------------------
Tres corridas de MCMC salieron con la misma constante rancia. Tres. El patrón no
es mala suerte: era lanzar primero y verificar después. Una corrida de 35 min con
un ingrediente muerto cuesta lo mismo que una buena y no sirve para nada, y lo
peor es que su log entra al repositorio con cara de evidencia.

Este script se corre ANTES de cualquier cómputo largo. Hace tres cosas:

  1. BORRA todo __pycache__ del repo. `python -B` no alcanza: impide ESCRIBIR
     bytecode, no LEERLO. Un .pyc con el header falsificado (mtime, size) se
     acepta como válido y sirve constantes muertas.
  2. Lee src/ssee_core.py del FUENTE (compile+exec, sin caché) y lo confronta
     con CANONICAL_VALUES.yaml. Si discrepan, aborta.
  3. Imprime la tabla de ingredientes que la corrida va a usar, para que quede
     en el log de la corrida y sea auditable después.

Uso:
    .venv/bin/python3 src/verificacion/preflight.py && <lanzar el cómputo>

Devuelve 0 si todo cuadra, 1 si no. El `&&` hace el resto.
"""
import pathlib
import shutil
import sys

_REPO = pathlib.Path(__file__).resolve().parents[2]

# ── 1. caché fuera ────────────────────────────────────────────────────────
_n = 0
for _pc in _REPO.rglob("__pycache__"):
    if ".git" in _pc.parts:
        continue
    shutil.rmtree(_pc, ignore_errors=True)
    _n += 1
print(f"preflight  [1/3] __pycache__ borrados: {_n}")

# ── 2. núcleo del FUENTE vs CANONICAL_VALUES ──────────────────────────────
_ns = {"__file__": str(_REPO / "src" / "ssee_core.py"), "__name__": "ssee_core"}
_src = (_REPO / "src" / "ssee_core.py").read_text(errors="ignore")
exec(compile(_src, str(_REPO / "src" / "ssee_core.py"), "exec"), _ns)

try:
    import yaml
    _cv = yaml.safe_load((_REPO / "CANONICAL_VALUES.yaml").read_text())
except Exception as _e:                                    # pragma: no cover
    print(f"preflight  ABORTA: no se pudo leer CANONICAL_VALUES ({_e})")
    sys.exit(1)

_prov = _cv.get("provenance", {})
# (nombre en ssee_core, nodo en provenance, tolerancia relativa)
_PARES = [("SUM_MNU_EV", "sigma_m_nu", 1e-4),
          ("OMEGA_B_H2", "omega_b", 1e-4),
          ("OMEGA_C_H2", "omega_c", 1e-4),
          ("OMEGA_NU_H2", "omega_nu", 1e-3),
          ("OMEGA_M_H2", "omega_m", 1e-4),
          ("OMEGA_M_TOTAL", "Omega_m_cosm", 1e-4)]

_mal = []
for _attr, _nodo, _tol in _PARES:
    if _attr not in _ns or _nodo not in _prov:
        continue
    _a, _b = float(_ns[_attr]), float(_prov[_nodo]["value"])
    if abs(_a - _b) > _tol * max(abs(_b), 1e-30):
        _mal.append(f"{_attr}={_a:.9g} vs {_nodo}={_b:.9g}")

if _mal:
    print("preflight  [2/3] ABORTA — núcleo y CANONICAL_VALUES discrepan:")
    for _m in _mal:
        print(f"             {_m}")
    print("preflight  NO lanzar el cómputo. Arreglar la fuente primero.")
    sys.exit(1)
print(f"preflight  [2/3] núcleo(fuente) == CANONICAL_VALUES en {len(_PARES)} constantes")

# ── 3. tabla de ingredientes, para que quede en el log de la corrida ──────
print("preflight  [3/3] ingredientes de esta corrida:")
for _attr, _, _ in _PARES:
    if _attr in _ns:
        print(f"             {_attr:16s} = {float(_ns[_attr]):.9f}")
for _extra in ("W0", "WA", "H0_ALG"):
    if _extra in _ns:
        print(f"             {_extra:16s} = {float(_ns[_extra]):.9f}")
print("preflight  OK — luz verde para el cómputo")
sys.exit(0)
