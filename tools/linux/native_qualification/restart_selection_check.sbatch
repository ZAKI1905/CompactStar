#!/bin/bash
#SBATCH --job-name=restart-selector-only
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_CHECK_OUTPUT:?}" "${TC:?}" "${CXX:?}"
cd "$RF_SOURCE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="$RF_CHECK_OUTPUT/tmp" PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_CHECK_OUTPUT/allocation.txt"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads \
  "$TC/python/bin/python" -B tools/linux/native_qualification/check_restart_coverage_proposal.py \
  "$RF_CHECK_OUTPUT/check" > "$RF_CHECK_OUTPUT/check.log" 2>&1
