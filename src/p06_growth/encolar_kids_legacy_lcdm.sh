#!/usr/bin/env bash
# Espera turno para la TERCERA configuracion de KiDS-Legacy (LCDM, fondo libre)
# y la lanza solo cuando hay sitio. Nace del 2026-09-19, 08:16: con dos
# corridas MPI de 4 procesos (~8 GB) el kernel se quedo sin memoria y mato a
# VS Code (oom_score_adj 300 lo pone primero). Una tercera corrida de 4
# procesos (~4 GB mas) habria repetido lo mismo, o algo peor.
#
# Condiciones para lanzar, las DOS a la vez:
#   1. hueco: menos de 2 corridas `cobaya_kids_legacy.py` vivas (se cuentan
#      los mpirun, no los hijos). Presupuesto real = 2 corridas de 4 procesos.
#   2. memoria: MemAvailable >= 5.5 GB (4 procesos a ~1 GB + margen).
# Si solo se cumple una, sigue esperando. Revisa cada 5 minutos.
#
# Lanza con el mismo lanzador de siempre (preflight + hilos a 1 + exec mpirun,
# reglas 1, 2 y 6 de COLA_CORRIDAS.md). Todo lo que decide queda en el log.
set -u
ROOT=/home/mike/Proyectos/SSEE
LOG=$ROOT/results/logs/cola_kids_legacy.log
MIN_KB=$((5500 * 1024))

echo "$(date '+%F %T')  en cola: kids_legacy/lcdm (fondo libre, 13 libres)" >> "$LOG"
while true; do
  vivas=$(pgrep -f "^mpirun .*cobaya_kids_legacy.py" | wc -l)
  libre=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  if [ "$vivas" -lt 2 ] && [ "$libre" -ge "$MIN_KB" ]; then
    echo "$(date '+%F %T')  HAY SITIO: $vivas corridas vivas, $((libre/1024)) MB libres -> lanzo lcdm" >> "$LOG"
    setsid nohup "$ROOT/src/p06_growth/lanzar_kids_legacy.sh" lcdm \
        > "$ROOT/results/logs/kids_legacy_lcdm.log" 2>&1 < /dev/null &
    echo "$(date '+%F %T')  lanzada, PID del lanzador $!" >> "$LOG"
    exit 0
  fi
  echo "$(date '+%F %T')  espera: $vivas corridas vivas, $((libre/1024)) MB libres" >> "$LOG"
  sleep 300
done
