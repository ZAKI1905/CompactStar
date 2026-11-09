#!/bin/bash
#SBATCH --job-name=native-checkpoint-diagnosis
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_REPLAY_RETAINED:?}" "${TC:?}"
cd "$RF_SOURCE"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/replay-allocation.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
 set -euo pipefail
 "$TC/python/bin/python" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$RF_OUTPUT/cpu-affinity.txt"
 exec "$TC/python/bin/python" -I -B tools/linux/measure.py "$TC/python/bin/python" -B \
  tools/linux/native_qualification/run_diagnostic.py "$RF_OUTPUT/build" "$RF_REPLAY_RETAINED" "$RF_OUTPUT/replay"
' > "$RF_OUTPUT/replay.log" 2>&1
