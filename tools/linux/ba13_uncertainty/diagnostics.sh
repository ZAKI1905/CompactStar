#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
OUT=${BA13_STUDY_OUT:-/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty}
for name in thermal_diagnostic fit_diagnostic; do
 "$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/$name.py" > "$OUT/$name.log" 2>&1
done
