"""Prueba de MUTACIÓN del guardián, dirigida por el registro de reglas.

Inyecta un defecto CONOCIDO en un manuscrito, corre el guardián, y exige que
falle **la regla que dice cubrirlo** — no cualquiera. Después restaura el archivo
y verifica que quedó idéntico.

POR QUÉ EXISTE (2026-07-29). Mike: «que el guardián me dé verde no significa que
esté bien». Un VERDE demuestra que nada saltó; NO demuestra que algo habría
saltado. Eso sólo se prueba rompiendo el documento a propósito.

POR QUÉ SE EXIGE LA REGLA CONCRETA. La primera versión sólo miraba si el guardián
enrojecía. Eso deja pasar el peor caso: la regla nueva no funciona, pero otra
salta por casualidad y el resultado parece correcto. Atribuir el fallo demuestra
que la regla añadida es la que trabaja.

Así se destapó el agujero original: «$n_s = 1-\\varphi^{-7}$ & $0.965123$» —un
valor sencillamente MAL, a 4.4e-4 del exacto— pasaba VERDE, porque la ventana de
identificación de R38 lo leía como «esta celda no habla de n_s» y los patrones de
R30 esperan «fórmula = valor» con el valor en la misma celda.

Uso:  python3 src/verificacion/mutacion_guardian.py
"""
import re
import sys
import pathlib
import subprocess

_AQUI = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_AQUI))
import registro_reglas as _reg   # noqa: E402

REPO = _AQUI.parent.parent
TEX = REPO / _reg.TEX_MUTACION
GUARDIAN = _AQUI / "ssee_verify.py"


def corre():
    """(verde?, [nombres de los checks que fallaron])."""
    r = subprocess.run([sys.executable, str(GUARDIAN)],
                       capture_output=True, text=True, cwd=REPO)
    fallos = re.findall(r"^\s*x\s+(\S+)", r.stdout, re.M)
    return ("VERDE" in r.stdout.split("=" * 40)[-1]), fallos


orig = TEX.read_text()
verde, _ = corre()
print(f"línea base: {'VERDE' if verde else 'ROJO'}")
if not verde:
    print("  el guardián ya está en rojo: arréglalo antes de mutar.")
    sys.exit(1)

problemas = []
for regla, info in _reg.REGLAS.items():
    for caso in info["mutacion"]:
        # Un caso puede declarar el archivo que muta. Sin declararlo, es el .tex.
        # Las capas de FÍSICA no viven en un manuscrito: se prueban tocando
        # `ssee_core.py` o el propio guardián, que es donde vive lo que afirman.
        arch, (nombre, viejo, nuevo) = ((info.get("archivo", _reg.TEX_MUTACION), caso)
                                        if len(caso) == 3 else (caso[0], caso[1:]))
        destino = REPO / arch
        base = destino.read_text()
        etiqueta = f"{regla} · {nombre}"
        if viejo not in base:
            problemas.append(f"{etiqueta}: el ancla ya no existe en {arch}")
            print(f"  [ ANCLA?  ] {etiqueta}")
            continue
        destino.write_text(base.replace(viejo, nuevo, 1))
        try:
            _, fallos = corre()
        finally:
            destino.write_text(base)      # restaurar SIEMPRE
        if not fallos:
            problemas.append(f"{etiqueta}: nadie lo detecta (VERDE por vacío)")
            print(f"  [  PASA   ] {etiqueta}")
        elif not any(f.startswith(tuple(info.get("prefijos", [regla]))) for f in fallos):
            problemas.append(f"{etiqueta}: lo detecta {fallos[0]}, no {regla}")
            print(f"  [ OTRA    ] {etiqueta} → lo caza {fallos[0]}")
        else:
            print(f"  [ DETECTA ] {etiqueta}")

# Control negativo: sin defecto, nadie debe disparar. Una regla que se queja del
# documento correcto es tan inútil como una que no se queja del roto.
TEX.write_text(orig)
verde, _ = corre()
print(f"  [{' CONTROL ' if verde else ' FALLA   '}] documento intacto → "
      f"{'VERDE, nadie dispara' if verde else 'ALGUIEN DISPARA SIN DEFECTO'}")
if not verde:
    problemas.append("control negativo: hay ruido sobre el documento correcto")

assert TEX.read_text() == orig, "¡el archivo NO volvió a su estado original!"
print("\narchivo restaurado idéntico ✓")

if problemas:
    print(f"\nMUTACIÓN-ROJO — {len(problemas)} caso(s):")
    for p in problemas:
        print("   x ", p)
    sys.exit(1)
print(f"\nMUTACIÓN-VERDE — {sum(len(v['mutacion']) for v in _reg.REGLAS.values())} "
      f"defectos inyectados, cada uno detectado por SU regla.")
