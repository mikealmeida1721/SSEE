#!/usr/bin/env bash
# Lanza 4 cadenas Cobaya independientes (sin MPI) para el MCMC FULL de SSEE.
# Cada cadena escribe en su propia subcarpeta c1..c4 (sin conflictos de output).
# 3 threads de CAMB por cadena (4×3 = 12 cores). Persisten al cierre con setsid.
# Al terminar, combinar las 4 con getdist para el R-1 Gelman-Rubin.
set -u

ROOT=/home/mike/Proyectos/SSEE
YAML=$ROOT/src/mcmc_full/ssee_full.yaml
BASE=$ROOT/class_ssee/output/mcmc_full
LOGDIR=$BASE/logs
mkdir -p "$LOGDIR"

export COBAYA_PACKAGES_PATH=/home/mike/cobaya_packages
export OMP_NUM_THREADS=3

for i in 1 2 3 4; do
  OUT=$BASE/c${i}/ssee_full
  mkdir -p "$BASE/c${i}"
  echo "$(date '+%H:%M:%S')  lanzando cadena $i  →  $OUT"
  setsid nohup $ROOT/.venv/bin/cobaya-run "$YAML" --output "$OUT" --force \
     > "$LOGDIR/cadena_${i}.log" 2>&1 &
  echo "  PID=$!  log=$LOGDIR/cadena_${i}.log"
  sleep 20
done

echo "$(date '+%H:%M:%S')  4 cadenas lanzadas."
echo "Monitoreo:  tail -f $LOGDIR/cadena_1.log"
