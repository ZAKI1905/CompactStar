#!/bin/bash
#SBATCH --job-name=stage-a-source-24
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=02:00:00
# 24 distinct physical cores, one execution thread per worker. Source A only.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_SNAPSHOTS:?}" "${TC:?}"
cd "$RF_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$RF_OUTPUT/tmp" RF_STAGE_ROLE=source RF_STAGE_WORKERS=24
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation.txt"
srun --ntasks=1 --cpus-per-task=24 --cpu-bind=cores bash -c '
 set -euo pipefail
 "$TC/python/bin/python" -I -B -c "import os; print(sorted(os.sched_getaffinity(0)))" > "$RF_OUTPUT/coordinator-affinity.txt"
 exec "$TC/python/bin/python" -I -B tools/linux/measure.py \
  "$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py stage "$RF_SNAPSHOTS" 1 "$RF_OUTPUT/run"
' > "$RF_OUTPUT/stage.log" 2>&1
