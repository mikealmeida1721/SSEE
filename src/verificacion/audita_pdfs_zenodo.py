#!/usr/bin/env python3
"""AUDITORIA DE LOS PDF QUE VAN A ZENODO (2026-09-19)

Por que existe
--------------
`memory_sync.py` ya sabe distinguir un valor retirado VIVO de uno NARRADO:
para eso estan `context_markers:` e `historical_sections:` en
CANONICAL_VALUES.yaml. Pero su cajon de papers es `manuscript/*.tex`.
A Zenodo NO suben los .tex: suben los PDF de `docs/`. Y entre el .tex y el
PDF hay dos maneras de que se cuele algo:

  1. el PDF es mas viejo que su .tex (se limpio la fuente y no se recompilo);
  2. el PDF no tiene .tex en el cajon barrido — es el caso de
     `submission_PRD/SSEE_PRD.tex`, que nadie estaba mirando.

Este script cierra las dos: lee el TEXTO IMPRESO de cada PDF de docs/ —lo que
un lector de Zenodo va a ver— y le aplica el MISMO criterio vivo-vs-narrado.

LA VENTANA VA EN CARACTERES, NO EN LINEAS — y esto se aprendio midiendo.
`memory_sync` mira la linea del hit y sus vecinas +-1, y eso funciona en un
.tex porque ahi la marca de contexto («superseded», «earlier») cae casi
siempre en la linea de al lado. En el texto impreso de un PDF las lineas son
trozos de parrafo cortados por la justificacion: la frase
«The earlier register used the CMB-mapped value H0=67.037» sale partida en
tres, y la marca «earlier» queda a +-2 o +-3 lineas del numero. Con ventana
+-1 los 7 hits de Paper 9 salian VIVOS estando perfectamente narrados en la
fuente. Por eso aqui se pega el texto de cada pagina en un solo flujo y la
ventana son +-360 caracteres alrededor del hit, que es el tamano de una
frase larga con su clausula de contexto.

R53 (el control del otro lado): dos controles, porque la ventana puede fallar
por los dos lados.
  · que SI aparezcan narrados — si salen 0 de las dos clases, el extractor no
    leyo texto y el «limpio» no significa nada;
  · que una frase de prueba con un valor retirado y SIN marca se clasifique
    VIVA. Sin este, una ventana demasiado ancha se tragaria cualquier marca
    del parrafo entero y todo saldria «narrado» — verde por ceguera.

Salida: results/logs/auditoria_pdfs_zenodo.json
"""
import json, pathlib, re, sys
import yaml, fitz

ROOT = pathlib.Path(__file__).resolve().parents[2]
CANON = ROOT / "CANONICAL_VALUES.yaml"
SALIDA = ROOT / "results" / "logs" / "auditoria_pdfs_zenodo.json"


RADIO = 360          # caracteres a cada lado del hit


def clasifica(flujo_bajo, pos, largo, marcas, hist):
    """¿El hit en `pos` esta NARRADO? Lo esta si en su vecindad de +-RADIO
    caracteres hay una marca de contexto o un encabezado de seccion historica."""
    win = flujo_bajo[max(0, pos - RADIO): pos + largo + RADIO]
    return any(m in win for m in marcas) or any(h in win for h in hist)


def hits(flujo, marcas, hist, retirados):
    bajo = flujo.lower()
    vivos, narrados = [], []
    for it in retirados:
        pat = it["pattern"].lower()
        req = [r.lower() for r in it.get("requires", [])]
        exc = [e.lower() for e in it.get("excludes", [])]
        for m in re.finditer(re.escape(pat), bajo):
            i = m.start()
            cerca = bajo[max(0, i - RADIO): i + len(pat) + RADIO]
            if req and not any(t in cerca for t in req):
                continue
            if exc and any(t in cerca for t in exc):
                continue
            reg = dict(patron=it["pattern"], pos=i,
                       texto=" ".join(flujo[max(0, i - 110): i + 110].split()))
            (narrados if clasifica(bajo, i, len(pat), marcas, hist)
             else vivos).append(reg)
    return vivos, narrados


def main():
    cfg = yaml.safe_load(CANON.read_text(encoding="utf-8"))
    retirados = cfg["retired"]
    marcas = [m.lower() for m in cfg["context_markers"]]
    hist = [h.lower() for h in cfg.get("historical_sections", [])]

    # CONTROL R53-b: una frase inventada con un valor retirado y sin marca
    # alguna TIENE que salir viva. Si sale narrada, la ventana esta ciega.
    sonda = ("The Hubble constant of the model is H0 = 67.037 km/s/Mpc and "
             "the matter density follows from it directly in every sector.")
    cv, cn = hits(sonda, marcas, hist, retirados)
    if not cv:
        print("  [CONTROL R53-b FALLA] la sonda sin marca salio narrada")
        return 2
    print("  control R53-b: la sonda sin marca sale VIVA  (%d)\n" % len(cv))

    informe, vivos_tot, narrados_tot = {}, 0, 0
    for pdf in sorted((ROOT / "docs").glob("*.pdf")):
        doc = fitz.open(pdf)
        flujo = "\n".join(pag.get_text() for pag in doc)
        doc.close()
        vivos, narrados = hits(flujo, marcas, hist, retirados)
        informe[pdf.name] = dict(caracteres=len(flujo), vivos=vivos, narrados=narrados)
        vivos_tot += len(vivos); narrados_tot += len(narrados)
        estado = "LIMPIO" if not vivos else "REVISAR %d" % len(vivos)
        print("  %-38s %7d car  vivos %3d  narrados %3d   %s"
              % (pdf.name, len(flujo), len(vivos), len(narrados), estado), flush=True)

    print("\n  TOTAL  vivos %d  ·  narrados %d" % (vivos_tot, narrados_tot))
    if narrados_tot == 0 and vivos_tot == 0:
        print("  [CONTROL R53 FALLA] ni vivos ni narrados: el extractor no leyo texto")
    SALIDA.write_text(json.dumps(dict(
        patrones=len(retirados), vivos=vivos_tot, narrados=narrados_tot,
        por_documento=informe), indent=1, ensure_ascii=False))
    print("  escrito -> %s" % SALIDA.relative_to(ROOT))
    return 1 if vivos_tot else 0


if __name__ == "__main__":
    sys.exit(main())
