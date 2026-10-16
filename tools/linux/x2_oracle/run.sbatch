#!/bin/bash
#SBATCH --job-name=x2-fixed-step
#SBATCH --partition=phyall
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=1024M
#SBATCH --time=00:10:00
#SBATCH --output=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/x2-fixed-step/slurm-%j.out
#SBATCH --error=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/x2-fixed-step/slurm-%j.err
set -euo pipefail
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
export PATH="$ROOT/toolchain/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
scontrol show job "$SLURM_JOB_ID" > "$ROOT/x2-fixed-step/allocation-$SLURM_JOB_ID.txt"
exec "$ROOT/toolchain/python/bin/python" -I -B "$ROOT/worktrees/x2-fixed-step/tools/linux/x2_oracle/run.py" "$1"
