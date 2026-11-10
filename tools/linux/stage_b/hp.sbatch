#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}"
cd "$RF_SOURCE"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
roles=(source control)
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python -B tools/linux/stage_b/reference.py "$RF_OUTPUT" "${roles[$SLURM_ARRAY_TASK_ID]}"
