#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}"
cd "$RF_SOURCE"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores "$TC/python/bin/python" -B tools/linux/stage_b/assess.py "$RF_OUTPUT" "$SLURM_CPUS_PER_TASK"
"$TC/python/bin/python" -B tools/linux/stage_b/chemical_response.py "$RF_OUTPUT"
