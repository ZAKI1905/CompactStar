#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}"
cd "$RF_SOURCE"
export PYTHONDONTWRITEBYTECODE=1
PY=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/final-allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task=24 --cpu-bind=cores "$PY" tools/linux/preliminary_pair/inventory.py "$RF_OUTPUT"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores "$PY" tools/linux/preliminary_pair/assess.py "$RF_OUTPUT"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores "$PY" tools/linux/preliminary_pair/supplement.py "$RF_OUTPUT"
