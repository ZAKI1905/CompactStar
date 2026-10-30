#!/bin/bash
#SBATCH --job-name=ba13-final
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=1024M
#SBATCH --time=00:03:00
set -euo pipefail
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
export PATH="$ROOT/toolchain/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
: "${BA13_FINAL_OUTPUT_ROOT:?explicit empty output root required}"
scontrol show job "$SLURM_JOB_ID" > "$BA13_FINAL_OUTPUT_ROOT/allocation-$SLURM_JOB_ID.txt"
exec "$ROOT/toolchain/python/bin/python" -B "$ROOT/worktrees/ba13-final-decision/tools/linux/ba13_final_decision/run.py" "$1"
