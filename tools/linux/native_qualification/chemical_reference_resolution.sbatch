#!/bin/bash
#SBATCH --job-name=bnv-d1-resolution
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
: "${D1_SOURCE:?}" "${D1_OUTPUT:?}" "${D1_PYTHON:?}"
cd "$D1_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PATH="$(dirname "$D1_PYTHON")/../../bin:$PATH"
export TMPDIR="$D1_OUTPUT/tmp"
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$D1_OUTPUT/allocation.txt"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
    set -euo pipefail
    "$D1_PYTHON" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$D1_OUTPUT/affinity.txt"
    "$D1_PYTHON" -B tools/linux/native_qualification/test_chemical_reference.py > "$D1_OUTPUT/tests-original.log" 2>&1
    "$D1_PYTHON" -B tools/linux/native_qualification/test_chemical_reference_resolution.py > "$D1_OUTPUT/tests-aligned.log" 2>&1
    exec "$D1_PYTHON" -B tools/linux/native_qualification/resolve_chemical_reference.py --output "$D1_OUTPUT/diagnostic"
' > "$D1_OUTPUT/execution.log" 2>&1
