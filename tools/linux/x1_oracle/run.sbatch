#!/bin/bash
#SBATCH --job-name=x1-oracle-anchor
#SBATCH --partition=phyall
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024M
#SBATCH --time=00:04:00
#SBATCH --output=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/x1-oracle-anchor/slurm-%j.out
#SBATCH --error=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/x1-oracle-anchor/slurm-%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
scontrol show job "$SLURM_JOB_ID" > "$ROOT/x1-oracle-anchor/slurm-allocation.txt"
exec "$ROOT/toolchain/python/bin/python" -I -B "$ROOT/worktrees/x1-oracle-anchor/tools/linux/x1_oracle/run.py"
