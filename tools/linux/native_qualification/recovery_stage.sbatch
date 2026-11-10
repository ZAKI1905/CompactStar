#!/bin/bash
#SBATCH --job-name=stage-a-recovery
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_SNAPSHOTS:?}" "${TC:?}"
cd "$RF_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$RF_OUTPUT/tmp"
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation.txt"
# This invocation can execute source only. Control needs a separate successful
# actual-output assessment; no shell dependency submits or launches it.
export RF_STAGE_ROLE=source
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
 set -euo pipefail
 "$TC/python/bin/python" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$RF_OUTPUT/affinity.txt"
 exec "$TC/python/bin/python" -I -B tools/linux/measure.py \
  "$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py stage "$RF_SNAPSHOTS" 1 "$RF_OUTPUT/run"
' > "$RF_OUTPUT/stage.log" 2>&1
