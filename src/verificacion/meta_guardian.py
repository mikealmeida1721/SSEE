"""META-GUARDIÁN — verifica al verificador.

`ssee_verify.py` comprueba el modelo. Nadie comprobaba a `ssee_verify.py`, y por
eso su VERDE significaba «nada saltó», no «algo habría saltado». Esto cierra ese
hueco: cruza el CÓDIGO del guardián contra el REGISTRO de reglas y contra la
prueba de MUTACIÓN, y falla si aparece cualquiera de las cuatro patologías.

  M1  toda capa del código tiene entrada en el registro (o está en la deuda)
  M2  toda entrada del registro corresponde a una capa que existe
  M3  toda regla registrada declara al menos un caso de mutación
  M4  toda exención nombra quién cubre el caso, o se declara PUNTO CIEGO
  M5  no hay dos reglas con la misma (intencion, ambito)
  M6  la deuda de capas sin cobertura no crece

M4 es la que responde a «reglas contradictorias que provocan que no vea el
error»: una exención sin `cubierta_por` es un agujero declarado, no un olvido.
Se permiten —a veces la exención es correcta y nada más la cubre— pero se
IMPRIMEN, para que sean una decisión visible y no un descubrimiento futuro.

Uso:  python3 src/verificacion/meta_guardian.py [--mutacion]
"""
import re
import sys
import pathlib
import collections

_AQUI = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_AQUI))
import registro_reglas as _reg   # noqa: E402

REPO = _AQUI.parent.parent
GUARDIAN = _AQUI / "ssee_verify.py"
DEUDA_MAX = 12          # medida 2026-07-29; sólo puede BAJAR

_fallos, _avisos = [], []


def chk(nombre, ok, detalle=""):
    print(f"  [{' OK ' if ok else 'FALLA'}] {nombre}" + (f"  — {detalle}" if detalle else ""))
    if not ok:
        _fallos.append(nombre)


print("META-GUARDIÁN — ¿el guardián puede demostrar lo que afirma?\n")

src = GUARDIAN.read_text(errors="ignore")
capas_codigo = [m.group(1) for m in re.finditer(r'print\("\\n(Capa [^"]+)"\)', src)]
registradas = {v["capa"]: k for k, v in _reg.REGLAS.items()}

# ── M1/M2 · el registro y el código dicen lo mismo ────────────────────────────
# Se comparan por PREFIJO: el texto de la capa lleva la descripción larga y basta
# con que empiece igual para identificarla sin obligar a repetirla carácter a
# carácter (que sería otra fuente de deriva).
def _casa(capa_codigo, capa_reg):
    return capa_codigo.replace("Capa ", "").startswith(capa_reg.split(" — ")[0])


sin_entrada = []
for c in capas_codigo:
    if any(_casa(c, r) for r in registradas):
        continue
    if c in _reg.SIN_COBERTURA:
        continue
    sin_entrada.append(c)
chk("M1 toda capa del código está registrada o declarada en deuda",
    not sin_entrada, "; ".join(sin_entrada) if sin_entrada
    else f"{len(capas_codigo)} capas: {len(registradas)} con cobertura, "
         f"{len(_reg.SIN_COBERTURA)} en deuda")

# Una regla se reconoce por CUALQUIERA de sus dos huellas en el código: una capa
# impresa, o comprobaciones con su prefijo. No todas imprimen cabecera —R30 vive
# dentro de otra capa— y exigir sólo lo primero daba un fallo espurio: la regla
# existía y trabajaba. Lo que M2 persigue es la entrada FANTASMA: registrada aquí
# y sin nada en el código que la respalde.
# La segunda huella son los PREFIJOS de sus comprobaciones. No basta con «R\d+»:
# las capas de física rotulan «V-L2», «canon», «DESI», y exigir el formato R-número
# las declaraba fantasma cuando existían y trabajaban.
nombres_check = re.findall(r'check\(\s*(?:f)?"([^"{]+)', src)
huerfanas = []
for cap, r in registradas.items():
    if any(_casa(c, cap) for c in capas_codigo):
        continue
    pref = tuple(_reg.REGLAS[r].get("prefijos", [r]))
    if not any(n.startswith(pref) for n in nombres_check):
        huerfanas.append(r)
chk("M2 toda regla registrada existe en el código",
    not huerfanas, "; ".join(huerfanas) if huerfanas
    else f"sin entradas fantasma ({len(nombres_check)} comprobaciones inspeccionadas)")

# ── M3 · nada entra sin poder probarse ───────────────────────────────────────
sin_mut = [k for k, v in _reg.REGLAS.items() if not v["mutacion"]]
chk("M3 toda regla registrada declara un caso de mutación",
    not sin_mut, "; ".join(sin_mut) if sin_mut
    else f"{sum(len(v['mutacion']) for v in _reg.REGLAS.values())} casos "
         f"para {len(_reg.REGLAS)} reglas")

# ── M4 · las exenciones no pueden abrir agujeros callados ────────────────────
ciegos, mal_cubiertas = [], []
for k, v in _reg.REGLAS.items():
    for motivo, cubre in v["exenciones"]:
        if cubre is None:
            ciegos.append(f"{k}: {motivo}")
        elif cubre not in _reg.REGLAS:
            mal_cubiertas.append(f"{k}: «{motivo}» dice cubrirse con {cubre}, que no está registrada")
chk("M4 ninguna exención apunta a una regla inexistente",
    not mal_cubiertas, "; ".join(mal_cubiertas) if mal_cubiertas
    else f"{len(ciegos)} punto(s) ciego(s) DECLARADO(s), ninguna referencia rota")
for c in ciegos:
    print(f"         · punto ciego declarado — {c}")

# ── M5 · dos reglas no pueden hacer lo mismo ─────────────────────────────────
firmas = collections.defaultdict(list)
for k, v in _reg.REGLAS.items():
    firmas[(v["intencion"], v["ambito"])].append(k)
choques = [f"{'/'.join(ks)} comparten ({i}, {a})"
           for (i, a), ks in firmas.items() if len(ks) > 1]
chk("M5 no hay dos reglas con la misma intención y ámbito",
    not choques, "; ".join(choques) if choques
    else f"{len(firmas)} pares (intención, ámbito) distintos")

# Aviso, no fallo: misma intención en ámbitos distintos es legítimo, pero es la
# señal temprana de la acumulación por parches. Se imprime para que se vea venir.
por_intencion = collections.Counter(v["intencion"] for v in _reg.REGLAS.values())
for inten, n in por_intencion.items():
    if n > 2:
        print(f"         · aviso: {n} reglas comparten la intención «{inten}» "
              f"en ámbitos distintos — candidatas a unificarse")

# ── M6 · la deuda sólo baja ──────────────────────────────────────────────────
chk("M6 la deuda de capas sin cobertura no crece",
    len(_reg.SIN_COBERTURA) <= DEUDA_MAX,
    f"{len(_reg.SIN_COBERTURA)} capas sin demostrar (tope {DEUDA_MAX})")

if "--mutacion" in sys.argv:
    print("\n— prueba de mutación —")
    import subprocess
    r = subprocess.run([sys.executable, str(_AQUI / "mutacion_guardian.py")],
                       capture_output=True, text=True, cwd=REPO)
    print(r.stdout.rstrip())
    if "PASA" in r.stdout.split("CONTROL")[0] or r.returncode:
        _fallos.append("mutación: algún defecto inyectado NO se detecta")

print("\n" + "=" * 60)
if _fallos:
    print(f"META-ROJO — {len(_fallos)} comprobación(es) fallaron:")
    for f in _fallos:
        print("   x ", f)
    sys.exit(1)
print(f"META-VERDE — el guardián demuestra {len(_reg.REGLAS)} de "
      f"{len(_reg.REGLAS) + len(_reg.SIN_COBERTURA)} capas. "
      f"Las {len(_reg.SIN_COBERTURA)} restantes están VERDES SIN DEMOSTRAR.")
