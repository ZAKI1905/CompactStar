#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
cases=(r320000 e32768 center025)
targets=(0 10 20)
a=${SLURM_ARRAY_TASK_ID:?}
exec "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/run.py" "${cases[$((a / 3))]}" "${targets[$((a % 3))]}"
