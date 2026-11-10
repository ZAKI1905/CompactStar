#!/bin/bash
set -euo pipefail
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export PYTHONDONTWRITEBYTECODE=1
: "${RF_ROLE:?}"
name="reference-$SLURM_ARRAY_TASK_ID"
if (( SLURM_ARRAY_TASK_ID == 5 )); then name=main; fi
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python -B tools/linux/preliminary_pair/roundoff.py /mnt/sdd/zaki/CompactStar/linux-transition-20261010/preliminary-pair-20261110-r1 "$RF_ROLE" "$name"
