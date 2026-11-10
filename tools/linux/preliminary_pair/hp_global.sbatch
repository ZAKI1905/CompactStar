#!/bin/bash
set -euo pipefail
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/preliminary-pair-20261110-r1
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
role=source
if (( SLURM_ARRAY_TASK_ID >= 5 )); then role=control; fi
index=$((SLURM_ARRAY_TASK_ID % 5))
steps=(1 2 4 8 8)
precision=110
if (( index == 4 )); then precision=80; fi
srun --cpu-bind=cores "$TC/python/bin/python" -B tools/linux/preliminary_pair/high_precision.py global "$RF_OUTPUT" --role "$role" --n "${steps[$index]}" --precision "$precision"
