#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
cases=(r20000 r40000 r160000 e2048 e4096 e16384 center2 center05 b5e-12 b5e-13 floor15 floor16)
targets=(0 10 20)
a=${SLURM_ARRAY_TASK_ID:?}
exec "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/run.py" "${cases[$((a / 3))]}" "${targets[$((a % 3))]}"
