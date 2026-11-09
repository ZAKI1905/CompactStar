#!/bin/bash
#SBATCH --job-name=native-p2-bounded-retry
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# One authorized Stage A retry; Stage B requires separate assessment of A.
# The retained collector and all its numerical settings remain unchanged.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_SNAPSHOTS:?}" "${RF_YEARS:?}" "${TC:?}"
cd "$RF_SOURCE"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
  set -euo pipefail
  "$TC/python/bin/python" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(\"EXECUTION_CPU_AFFINITY\", a, flush=True); assert len(a)==1" > "$RF_OUTPUT/cpu-affinity.txt"
  exec bash tools/linux/native_pair/stage.sbatch
'
