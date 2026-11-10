#!/bin/bash
set -euo pipefail
: "${RF_OUTPUT:?}" "${SLURM_ARRAY_TASK_ID:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores /bin/bash "$RF_OUTPUT/direct/task-$SLURM_ARRAY_TASK_ID.sh"
