#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
cases=(nominal r20000 r40000 r160000)
targets=(5 15)
a=${SLURM_ARRAY_TASK_ID:?}
exec "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/run.py" "${cases[$((a / 2))]}" "${targets[$((a % 2))]}"
