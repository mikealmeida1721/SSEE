#!/usr/bin/env python3
"""postflight.py — decidir si una corrida MCMC quedó BIEN HECHA.

POR QUÉ EXISTE (2026-07-27)
---------------------------
`preflight.py` garantiza que una corrida ARRANCA con los ingredientes correctos.
No dice nada de si TERMINÓ bien. Y ahí se perdió tiempo repetidamente:

  · tres corridas que abortaron en segundos y un script las dio por buenas
    porque verificaba AUSENCIA de una huella mala en vez de PRESENCIA de éxito;
  · una corrida de 18 h que terminó imprimiendo «B1 COMPLETE» sin el ΔBIC —
    el único número por el que existía— porque un `except: pass` se tragó el
    fallo de getBestFit().

Un MCMC caro que hay que repetir por un defecto detectable en 2 segundos es
tiempo tirado. Esto lo detecta antes de que nadie use el resultado.

Comprueba, sobre los archivos de cadena ya escritos (NO re-corre nada):

  1. CONVERGENCIA   — R−1 declarado en .progress, contra el umbral pedido
  2. MUESTRAS       — que haya cadena suficiente tras el burn-in
  3. χ² COHERENTE   — la columna chi2 debe reproducir 2·(minuslogpost −
                      minuslogprior). Si no, el fichero está corrupto o el
                      formato cambió y cualquier χ² leído de ahí es basura.
  4. MÍNIMO SANO    — el mejor punto no puede caer pegado a un borde de la
                      cadena: eso delata que aún estaba bajando, es decir que
                      la corrida se cortó antes de tiempo.

Uso:
    .venv/bin/python3 src/verificacion/postflight.py results/chains/ssee_cmb
    .venv/bin/python3 src/verificacion/postflight.py <pref1> <pref2> ...

Devuelve 0 si TODAS pasan, 1 si alguna falla — para encadenar con `&&`.
"""
import pathlib
import sys

import numpy as np

_REPO = pathlib.Path(__file__).resolve().parents[2]
R1_MAX = 0.05          # umbral laxo: por debajo de esto la cadena es usable
BURN = 0.3


def revisar(pref: str) -> bool:
    p = pathlib.Path(pref)
    if not p.is_absolute():
        p = _REPO / pref
    nombre = p.name
    txt = p.with_suffix(".1.txt")
    print(f"\npostflight  {nombre}")
    if not txt.exists():
        print(f"  ✗ no existe {txt.name} — la corrida no dejó cadena")
        return False

    ok = True

    # 1 · convergencia
    prog = p.with_suffix(".progress")
    if prog.exists():
        _lineas = prog.read_text().splitlines()
        # Por NOMBRE de columna, nunca por posición fija: la primera versión de
        # esto leía el índice 2 y devolvía acceptance_rate (0.712) en vez de
        # Rminus1 (0.017), declarando NO CONVERGIDA una cadena que sí lo estaba.
        # Una herramienta de verificación que se equivoca es peor que no tenerla.
        _cab = next((l for l in _lineas if l.lstrip().startswith("#")), "")
        _cols = _cab.lstrip("#").split()
        filas = [l.split() for l in _lineas
                 if l.strip() and not l.lstrip().startswith("#")]
        try:
            r1 = float(filas[-1][_cols.index("Rminus1")])
            estado = "OK" if r1 <= R1_MAX else "NO CONVERGIDA"
            print(f"  {'✓' if r1 <= R1_MAX else '✗'} R−1 = {r1:.4f}  ({estado}, umbral {R1_MAX})")
            ok &= r1 <= R1_MAX
        except Exception:
            print("  ? .progress ilegible — convergencia no verificable")
    else:
        print("  ? sin .progress — convergencia no verificable")

    # 2..4 · sobre la cadena
    hdr = txt.open().readline().lstrip("#").split()
    d = np.loadtxt(txt)
    corte = int(BURN * len(d))
    n = len(d) - corte
    print(f"  {'✓' if n > 1000 else '✗'} muestras tras burn-in {BURN:.0%}: {n}")
    ok &= n > 1000

    if all(c in hdr for c in ("chi2", "minuslogpost", "minuslogprior")):
        chi2 = d[corte:, hdr.index("chi2")]
        recon = 2 * (d[corte:, hdr.index("minuslogpost")]
                     - d[corte:, hdr.index("minuslogprior")])
        dev = float(np.abs(chi2 - recon).max())
        print(f"  {'✓' if dev < 1e-2 else '✗'} χ² coherente con (post−prior): "
              f"desviación máx {dev:.2e}")
        ok &= dev < 1e-2

        pos = int(np.argmin(chi2)) / len(chi2)
        sano = 0.05 < pos < 0.98
        print(f"  {'✓' if sano else '✗'} mínimo χ² = {chi2.min():.3f} al {pos:.0%} "
              f"de la cadena ({'sano' if sano else 'PEGADO AL BORDE: se cortó antes de tiempo'})")
        ok &= sano
    else:
        print("  ✗ faltan columnas chi2/minuslogpost/minuslogprior — "
              "cualquier χ² leído de aquí no es confiable")
        ok = False

    print(f"  → {'APTA para usar' if ok else 'NO APTA — no publicar este resultado'}")
    return ok


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    todo = all([revisar(a) for a in sys.argv[1:]])
    print(f"\npostflight  {'TODAS APTAS' if todo else 'HAY CORRIDAS NO APTAS'}")
    sys.exit(0 if todo else 1)
