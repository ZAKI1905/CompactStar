#!/bin/bash
#SBATCH --job-name=bnv-thermo-d1
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
# Compute nodes use the retained cluster toolchain, including its Git launcher.
export PATH="$(dirname "$D1_PYTHON")/../../bin:$PATH"
export TMPDIR="$D1_OUTPUT/tmp"
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$D1_OUTPUT/allocation.txt"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
    set -euo pipefail
    "$D1_PYTHON" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$D1_OUTPUT/affinity.txt"
    "$D1_PYTHON" -B tools/linux/native_qualification/test_chemical_reference.py > "$D1_OUTPUT/tests.log" 2>&1
    args=(--output "$D1_OUTPUT/diagnostic")
    if [[ -n "${D1_ADJUDICATION:-}" ]]; then args+=(--adjudication "$D1_ADJUDICATION"); fi
    if [[ -n "${D1_ADJUDICATION_SHA256:-}" ]]; then args+=(--adjudication-sha256 "$D1_ADJUDICATION_SHA256"); fi
    exec "$D1_PYTHON" -B tools/linux/native_qualification/verify_chemical_reference.py "${args[@]}"
' > "$D1_OUTPUT/execution.log" 2>&1
