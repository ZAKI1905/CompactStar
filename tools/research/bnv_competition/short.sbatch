#!/bin/bash
set -euo pipefail
: "${COMP_SOURCE:?}" "${COMP_OUTPUT:?}" "${COMP_COMMIT:?}"
cd "$COMP_SOURCE"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
scontrol show job "$SLURM_JOB_ID" > "$COMP_OUTPUT/short-allocation-$SLURM_JOB_ID.txt"
exec srun python3 -B tools/research/bnv_competition/run_short.py "$COMP_OUTPUT"
