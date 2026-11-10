#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_MODE:?}"
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
export PATH="$TC/bin:$PATH" RF_STAGE_WORKERS=24
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RF_SOURCE"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores "$TC/python/bin/python" -B tools/linux/stage_b/supplement.py "$RF_OUTPUT" "$RF_MODE"
