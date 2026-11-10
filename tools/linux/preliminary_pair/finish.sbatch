#!/bin/bash
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_EXTERNAL_EVIDENCE:?}"
cd "$RF_SOURCE"
export PYTHONDONTWRITEBYTECODE=1
BASE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
PY="$BASE/toolchain/python/bin/python"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/finish-allocation-$SLURM_JOB_ID.txt"
srun --ntasks=1 --cpus-per-task=24 --cpu-bind=cores "$PY" tools/linux/preliminary_pair/summarize.py "$RF_OUTPUT"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores "$BASE/native-p2-plotting-env/bin/python" tools/linux/preliminary_pair/compare.py "$RF_OUTPUT" "$RF_SOURCE/docs/validation/results-first-bnv-pilot/preliminary-pair/evidence"
srun --ntasks=1 --cpus-per-task=24 --cpu-bind=cores "$PY" tools/linux/preliminary_pair/seal.py "$RF_OUTPUT" "$RF_EXTERNAL_EVIDENCE"
