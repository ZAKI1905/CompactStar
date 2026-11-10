#!/bin/bash
#SBATCH --job-name=stage-a-recovery-static
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${TC:?}"
cd "$RF_SOURCE"
export TMPDIR="$RF_OUTPUT/tmp" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/static-allocation.txt"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
 set -euo pipefail
 "$TC/python/bin/python" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$RF_OUTPUT/static-affinity.txt"
 exec "$TC/python/bin/python" -B tools/linux/native_qualification/recovery_snapshots.py verify-publish "$RF_OUTPUT"
' > "$RF_OUTPUT/static-stage.log" 2>&1
