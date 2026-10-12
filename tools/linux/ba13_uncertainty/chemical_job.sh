#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty
exec "$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/chemical_diagnostic.py" > "$OUT/chemical-diagnostic.log" 2>&1
