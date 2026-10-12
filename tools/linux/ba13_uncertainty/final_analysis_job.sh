#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
OUT=${BA13_STUDY_OUT:-/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty}
for name in profile_diagnostic p1_validation; do
 "$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/$name.py" > "$OUT/${name//_/-}.log" 2>&1
done
