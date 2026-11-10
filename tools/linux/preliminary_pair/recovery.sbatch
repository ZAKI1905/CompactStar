#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_MODE:?}"
cd "$RF_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores "$TC/python/bin/python" -B tools/linux/preliminary_pair/recovery.py "$RF_MODE" "$RF_OUTPUT"
