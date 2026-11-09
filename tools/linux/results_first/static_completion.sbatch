#!/bin/bash
#SBATCH --job-name=rf-bnv-static
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=1G
#SBATCH --time=00:05:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_COMMIT:?}" "${TC:?}"
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
exec "$TC/python/bin/python" -I -B "$RF_SOURCE/tools/linux/results_first/static_completion_run.py"
