#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
export PYTHONDONTWRITEBYTECODE=1
case "${SLURM_ARRAY_TASK_ID:?}" in
 0) mode=source;;
 1) mode=control;;
 2) mode=events;;
 *) exit 2;;
esac
"$TC/python/bin/python" "$SRC/tools/linux/run_bounded.py" --handoff /mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010 --executable "$TC/build/CompactStar-Debug/tests/phase6_breakpoint_qualification" --output "$TC/bounded/$mode-${SLURM_ARRAY_JOB_ID}" --mode "$mode" --timeout 1800
