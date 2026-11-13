#!/bin/bash
# Resources selected at submission to match independent stationary-root work.
set -euo pipefail
: "${COMP_SOURCE:?}" "${COMP_OUTPUT:?}"
BASE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
cd "$COMP_SOURCE"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$BASE/bnv-rotochemical-competition-20261112-r1/python-deps"
export MPLCONFIGDIR="$COMP_OUTPUT/matplotlib-cache"
mkdir -p "$MPLCONFIGDIR"
scontrol show job "$SLURM_JOB_ID" > "$COMP_OUTPUT/allocation-$SLURM_JOB_ID.txt"
exec srun --cpu-bind=cores "$BASE/native-p2-plotting-env/bin/python" -B \
  tools/research/bnv_competition/scan.py "$COMP_OUTPUT" --workers "$SLURM_CPUS_PER_TASK" "${@}"
