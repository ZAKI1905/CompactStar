#!/bin/bash
#SBATCH --job-name=ba13-evidence
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=1024M
#SBATCH --time=00:05:00
#SBATCH --output=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-evidence-20261026/slurm-%j.out
#SBATCH --error=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-evidence-20261026/slurm-%j.err
set -euo pipefail
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
export PATH="$ROOT/toolchain/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
scontrol show job "$SLURM_JOB_ID" > "$ROOT/ba13-evidence-20261026/allocation-$SLURM_JOB_ID.txt"
exec "$ROOT/toolchain/python/bin/python" -B "$ROOT/worktrees/ba13-evidence/tools/linux/ba13_evidence/$1.py" "${@:2}"
