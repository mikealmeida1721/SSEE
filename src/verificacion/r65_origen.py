#!/usr/bin/env python3
"""Lista los numeros SIN ORIGEN que ve R65, sin correr el guardian entero.

Ejecuta exactamente el bloque de R65 de `ssee_verify.py` (no una copia: si
la regla cambia, esto cambia con ella) y muestra, por script, cada numero y la
linea donde vive.

    python3 src/verificacion/r65_origen.py            # resumen por script
    python3 src/verificacion/r65_origen.py kids       # detalle de los que casan
"""
import importlib.util
import pathlib
import re
import sys

SRC = pathlib.Path(__file__).resolve().parent.parent
_v = (SRC / "verificacion" / "ssee_verify.py").read_text()
_blk = _v[_v.index("# R65 CRECE (2026-09-19)"):
          _v.index("# CONTROL (R53): cada via de origen")]
_sp = importlib.util.spec_from_file_location("_core65", SRC / "ssee_core.py")
_core = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(_core)
g = dict(re=re, pathlib=pathlib, ROOT=SRC, _core63=_core, _ERR_CORE=None,
         _R65_NUM=re.compile(r"(?<![\w.])(\d+\.\d{4,})(?![\w])"),
         track_open=lambda *a, **k: None, print=lambda *a, **k: None,
         check=lambda *a, **k: None, _DEUDA_REAL={}, _DEUDA_MAX={})
exec(_blk, g)
sin = g["_sin65"]
filtro = sys.argv[1] if len(sys.argv) > 1 else ""
print(f"TOTAL {sum(len(v) for v in sin.values())} numeros sin origen "
      f"en {len(sin)} scripts")
for k, v in sorted(sin.items(), key=lambda kv: -len(kv[1])):
    if filtro and filtro not in k:
        continue
    print(f"## {k} ({len(v)})")
    if filtro:
        for x, c in v:
            print(f"    {x:>14} | {c}")
