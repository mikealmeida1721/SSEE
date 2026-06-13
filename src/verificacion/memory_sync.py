#!/usr/bin/env python3
"""
memory_sync.py — el sincronizador de las 3 memorias de SSEE.

Las 3 memorias deben concordar SIEMPRE (regla de Mike):
  • Guardián  — VERIFICATION_LEDGER.md (resultados/valores)
  • Obsidian  — /home/mike/SSEE-Vault (conexiones/cadenas)
  • CLAUDE.md — contexto/estado

Este script lee CANONICAL_VALUES.yaml (fuente única de verdad) y revisa que
NINGUNA memoria presente un valor RETIRADO como vigente. Un valor retirado
sólo se permite si en su misma línea hay una marca de contexto
(retirado, viejo, Type-P, coincidencia, sin re-correr…).

    .venv/bin/python3 src/verificacion/memory_sync.py            # las 3 memorias
    .venv/bin/python3 src/verificacion/memory_sync.py --vault    # sólo el vault

  VERDE → memorias sincronizadas.   DRIFT → una memoria quedó desfasada.

Diseñado para ser barato: correr tras cada cambio de valor. También lo invoca
ssee_verify.py (el Guardián) como su capa de coherencia de memorias.
"""
import pathlib
import re
import sys

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
VAULT = pathlib.Path("/home/mike/SSEE-Vault")
CANON = ROOT / "CANONICAL_VALUES.yaml"


def _load():
    with open(CANON, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _targets(vault_only=False):
    """Las memorias a escanear, como (etiqueta, lista de Paths)."""
    out = []
    if not vault_only:
        out.append(("Guardián (Ledger)", [ROOT / "VERIFICATION_LEDGER.md"]))
        out.append(("CLAUDE.md", [ROOT / "CLAUDE.md"]))
    if VAULT.exists():
        out.append(("Obsidian (vault)",
                    sorted(VAULT.rglob("*.md"))))
    elif not vault_only:
        print(f"  [aviso] vault no encontrado en {VAULT} — se omite Obsidian")
    return out


def _marked(lines_low, idx, markers):
    """¿La línea idx (o sus vecinas ±1) contiene una marca de contexto?
    La ventana ±1 cubre el caso de notas que parten valor y marca en
    líneas contiguas (p. ej. «…da 72.86» / «(0.17σ) — coincidencia Type-P»)."""
    for j in (idx - 1, idx, idx + 1):
        if 0 <= j < len(lines_low) and any(mk in lines_low[j] for mk in markers):
            return True
    return False


def scan(vault_only=False):
    """Devuelve (drifts, scanned). drifts = lista de (memoria, archivo, lineno, patrón, texto)."""
    cfg = _load()
    retired = cfg["retired"]
    markers = [m.lower() for m in cfg["context_markers"]]
    hist_heads = [h.lower() for h in cfg.get("historical_sections", [])]
    drifts = []
    scanned = 0

    for label, paths in _targets(vault_only):
        for path in paths:
            if not path.exists():
                continue
            scanned += 1
            lines = path.read_text(encoding="utf-8").splitlines()
            low = [ln.lower() for ln in lines]
            in_hist = False
            hist_level = 0
            for i, raw in enumerate(lines):
                # Seguimiento de sección con NIVEL de encabezado: una sección
                # histórica (##) sigue siéndolo a través de sus sub-encabezados
                # (###) y sólo se cierra con un encabezado hermano o superior.
                stripped = raw.lstrip()
                if stripped.startswith("#"):
                    level = len(stripped) - len(stripped.lstrip("#"))
                    if any(h in low[i] for h in hist_heads):
                        in_hist, hist_level = True, level
                    elif in_hist and level <= hist_level:
                        in_hist, hist_level = False, 0
                if in_hist or _marked(low, i, markers):
                    continue
                for item in retired:
                    pat = item["pattern"]
                    if pat.lower() in low[i]:
                        rel = path.relative_to(VAULT if label.startswith("Obsidian") else ROOT)
                        drifts.append((label, str(rel), i + 1, pat, raw.strip()[:90]))
    return drifts, scanned


def run(vault_only=False, verbose=True):
    """Imprime el informe y devuelve la lista de drifts (vacía = sincronizado)."""
    drifts, scanned = scan(vault_only)
    if verbose:
        print(f"memory_sync — {scanned} archivos de memoria escaneados")
        if not drifts:
            print("  VERDE — las memorias concuerdan con CANONICAL_VALUES.yaml.")
        else:
            print(f"  DRIFT — {len(drifts)} valor(es) retirado(s) sin marcar:")
            for label, rel, ln, pat, txt in drifts:
                print(f"   [{label}] {rel}:{ln}  «{pat}»  →  {txt}")
    return drifts


if __name__ == "__main__":
    vo = "--vault" in sys.argv
    sys.exit(1 if run(vault_only=vo) else 0)
