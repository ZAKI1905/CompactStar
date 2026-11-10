#!/bin/bash
#SBATCH --job-name=echem-local-feasibility
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_ACTION:?}" "${TC:?}"
cd "$RF_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$RF_OUTPUT/tmp"
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/$RF_ACTION-allocation.txt"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
 set -euo pipefail
 "$TC/python/bin/python" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$RF_OUTPUT/$RF_ACTION-affinity.txt"
 if [[ "$RF_ACTION" == build ]]; then
  exec "$TC/python/bin/python" -I -B tools/linux/measure.py "$TC/python/bin/python" -B tools/linux/native_qualification/build_refinement.py "$TC/../native-p2-regeneration-20261108-r1" "$RF_OUTPUT/build"
 elif [[ "$RF_ACTION" == study ]]; then
  exec "$TC/python/bin/python" -I -B tools/linux/measure.py "$TC/python/bin/python" -B tools/linux/native_qualification/run_refinement.py "$RF_OUTPUT/build" "$TC/../native-p2-qualification-1yr-r2" "$RF_OUTPUT/study"
 else exit 2; fi
' > "$RF_OUTPUT/$RF_ACTION.log" 2>&1
