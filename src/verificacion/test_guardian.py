#!/usr/bin/env python3
"""
test_guardian.py — el vigilante vigilado ("who watches the watcher").

El guardián (`ssee_verify.py`) protege contra el FALSO VERDE (decir VERDE
cuando debería ser ROJO). Pero un check puede dejar de detectar en silencio
tras un refactor → el guardián daría falso verde sin que nadie lo note.

Este test vuelve PERMANENTES las pruebas negativas: para cada check crítico,
perturba un cajón a un estado que DEBE disparar ROJO, corre el guardián, y
exige que (a) salga con código ≠ 0 y (b) aparezca el check correcto. Luego
restaura el archivo (siempre, vía finally). Si un check deja de cazar lo suyo,
ESTE test falla → la grieta se delata en vez de pasar callada.

    .venv/bin/python3 src/verificacion/test_guardian.py
    → TODOS VERDES: el guardián detecta lo que debe.
    → ROJO: un check se debilitó; arreglarlo antes de confiar en un VERDE.
"""
import contextlib
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
GUARD = "src/verificacion/ssee_verify.py"
PY = sys.executable
_passed, _failed = 0, 0


def run_guardian():
    r = subprocess.run([PY, GUARD], cwd=REPO, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


@contextlib.contextmanager
def perturbed(relpath, transform):
    """Perturba un archivo y SIEMPRE lo restaura (bytes originales)."""
    p = REPO / relpath
    original = p.read_bytes()
    try:
        p.write_text(transform(original.decode("utf-8", errors="ignore")),
                     encoding="utf-8")
        yield
    finally:
        p.write_bytes(original)


def expect_red(name, relpath, transform, must_contain):
    global _passed, _failed
    with perturbed(relpath, transform):
        rc, out = run_guardian()
    low = out
    ok = rc != 0 and all(k in low for k in must_contain)
    print(f"  [{'PASS' if ok else 'FALLA'}] {name}")
    if not ok:
        print(f"        esperaba ROJO con {must_contain}; rc={rc}")
    _passed += ok
    _failed += (not ok)


# ── Test positivo: estado limpio → VERDE (exit 0) ───────────────────────
print("Vigilante vigilado — el guardián debe detectar cada falla simulada\n")
rc0, out0 = run_guardian()
_ok0 = rc0 == 0 and "VERDE" in out0
print(f"  [{'PASS' if _ok0 else 'FALLA'}] estado limpio → VERDE (exit 0)")
_passed += _ok0
_failed += (not _ok0)

# ── Tests negativos: cada perturbación DEBE dar ROJO ────────────────────
# F1 — valor de pipeline que NO coincide con su log de procedencia (R9)
expect_red("R9 procedencia: ΔBIC −24.7 no está en el log → ROJO",
           "VERIFICATION_LEDGER.md",
           lambda t: t.replace("−34.9 (SSEE favorecido)", "−24.7 (SSEE favorecido)"),
           ["procedencia", "ROJO"])

# R1 — claim de prioridad temporal
expect_red("R1 cronología: 'predating DESI DR2' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n% test predating DESI DR2\n",
           ["R1", "ROJO"])

# R2 — enlace al compendio multidominio retirado
expect_red("R2 multidominio: DOI 19679049 → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n% test zenodo.19679049\n",
           ["R2", "ROJO"])

# R10 — overclaim 'zero-parameter' global
expect_red("R10 overclaim: 'is a zero-parameter framework' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\nSSEE is a zero-parameter framework.\n",
           ["R10", "ROJO"])

# R11 — conteo de serie congelado
expect_red("R11 conteo serie: 'Papers~3--7' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\nextensions in Papers~3--7.\n",
           ["R11", "ROJO"])


# ── Archivo: cajón sin entrada en la Bitácora → ROJO (caso especial dir) ─
def test_archive_orphan():
    global _passed, _failed
    d = REPO / "archive" / "_test_orphan_dir"
    try:
        d.mkdir(exist_ok=True)
        rc, out = run_guardian()
    finally:
        if d.exists():
            d.rmdir()
    ok = rc != 0 and "archivo" in out
    print(f"  [{'PASS' if ok else 'FALLA'}] R12 archivo: cajón sin bitácora → ROJO")
    _passed += ok
    _failed += (not ok)


test_archive_orphan()

# ── Veredicto ───────────────────────────────────────────────────────────
print(f"\n{'='*55}")
total = _passed + _failed
if _failed:
    print(f"ROJO — {_failed}/{total} pruebas del vigilante FALLARON.")
    print("Un check del guardián se debilitó: arreglarlo antes de confiar en un VERDE.")
    sys.exit(1)
print(f"VERDE — el guardián detecta las {total} fallas simuladas. Vigilante íntegro.")
sys.exit(0)
