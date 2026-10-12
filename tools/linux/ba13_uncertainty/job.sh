#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
export PYTHONDONTWRITEBYTECODE=1
exec "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/run.py" "${1:?case}" "${2:-$SLURM_ARRAY_TASK_ID}"
