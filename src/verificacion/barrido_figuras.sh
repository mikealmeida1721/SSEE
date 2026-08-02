#!/bin/bash
# Barrido de figuras: una referencia esta OK si el archivo existe en
#   (a) results/figures/   (graphicspath habitual)
#   (b) el directorio del propio .tex   (ruta relativa simple)
#   (c) la ruta explicita tal cual, resuelta desde el dir del .tex
cd /home/mike/Proyectos/SSEE
for f in manuscript/*.tex submission_PRD/*.tex; do
  d=$(dirname "$f")
  for g in $(grep -oP 'includegraphics(\[[^]]*\])?\{\K[^}]+' "$f" 2>/dev/null); do
    ok=0
    for cand in "results/figures/$(basename "$g")" "$d/$g" "$d/$(basename "$g")"; do
      for ext in "" .pdf .png .jpg; do
        [ -f "${cand}${ext}" ] && ok=1 && break 2
      done
    done
    [ $ok -eq 0 ] && echo "  ${f##*/}: $g"
  done
done | sort -u
