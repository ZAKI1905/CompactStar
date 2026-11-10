#!/bin/bash
# Resource counts/time are selected explicitly at submission.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_MODE:?}"
cd "$RF_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export RF_STAGE_WORKERS=24
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores "$TC/python/bin/python" -B tools/linux/preliminary_pair/run.py "$RF_MODE" "$RF_OUTPUT" "${SLURM_ARRAY_TASK_ID:-0}"
