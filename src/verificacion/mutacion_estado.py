"""Mutaciones de ESTADO: las que no se escriben en un fichero.

POR QUE EXISTE (2026-09-08). Las 23 reglas restantes se probaron con un
reemplazo de texto: se escribe un defecto en un archivo, se corre el guardian y
se exige que enrojezca ESA regla. R31 no se deja: lo que vigila es un `.pyc`
rancio, o sea un estado del sistema de ficheros con el fuente INTACTO byte a
byte. Un caso (archivo, viejo, nuevo) no puede expresarlo.

En vez de declararla «no mutable» —que es como se acumulan los puntos ciegos—
se le da a la suite una segunda forma de caso: un preparador y un restaurador.

EL BUG QUE REPRODUCE es real y esta documentado: un `.pyc` rancio dio 0.06902
donde el fuente decia 0.06849, y el guardian amanecio en ROJO por una causa que
no estaba escrita en ningun sitio. Reproducirlo tiene una sutileza: Python
valida el `.pyc` por (mtime, tamano) del `.py`, asi que hay que compilar el
fuente mutado CON EL MTIME DEL ORIGINAL y con el mismo numero de bytes. Si se
compila y despues se reajusta el mtime, el `.pyc` queda invalido y Python
recompila: el primer intento fallo justo por eso y R31 salio verde sin que
hubiera nada que ver.

Uso:  python3 mutacion_estado.py r31_pyc_rancio prepara|restaura
"""
import os
import pathlib
import py_compile
import re
import shutil
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
CORE = REPO / "src" / "ssee_core.py"
CACHE = REPO / "src" / "__pycache__"

# (regla que debe enrojecer, que hace el caso)
CASOS = {
    "r31_pyc_rancio": ("R31", "un .pyc rancio: el fuente dice 0.06849 y quien "
                              "lo importa recibe 0.06850"),
}


def prepara():
    orig = CORE.read_bytes()
    st = os.stat(CORE)
    _m = re.search(rb"SUM_MNU_EV\s*=\s*0\.06849", orig)
    if not _m:
        raise SystemExit("el ancla SUM_MNU_EV = 0.06849 ya no esta en el nucleo")
    mut = orig.replace(_m.group(0), _m.group(0).replace(b"0.06849", b"0.06850"))
    if len(mut) != len(orig):
        raise SystemExit("la mutacion cambia el tamano: Python invalidaria el .pyc")
    shutil.rmtree(CACHE, ignore_errors=True)
    CORE.write_bytes(mut)
    os.utime(CORE, (st.st_atime, st.st_mtime))   # el .pyc registra ESTE mtime
    py_compile.compile(str(CORE), doraise=True)
    CORE.write_bytes(orig)                       # fuente identico byte a byte...
    os.utime(CORE, (st.st_atime, st.st_mtime))   # ...y con su mtime original


def restaura():
    shutil.rmtree(CACHE, ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in CASOS:
        raise SystemExit(f"uso: {sys.argv[0]} {'|'.join(CASOS)} prepara|restaura")
    {"prepara": prepara, "restaura": restaura}[sys.argv[2]]()
