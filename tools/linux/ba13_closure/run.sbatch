#!/bin/bash
#SBATCH --job-name=ba13-closure
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=1024M
#SBATCH --time=00:10:00
#SBATCH --output=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-closure-20261019/slurm-%j.out
#SBATCH --error=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-closure-20261019/slurm-%j.err
set -euo pipefail
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
export PATH="$ROOT/toolchain/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
BA13_OUTPUT=${BA13_CLOSURE_OUTPUT_ROOT:-"$ROOT/ba13-closure-20261019"}
mkdir -p "$BA13_OUTPUT"
scontrol show job "$SLURM_JOB_ID" > "$BA13_OUTPUT/allocation-$SLURM_JOB_ID.txt"
exec "$ROOT/toolchain/python/bin/python" -B "$ROOT/worktrees/ba13-closure/tools/linux/ba13_closure/$1.py" "${@:2}"
