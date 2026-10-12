#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty
exec "$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/integral_diagnostic.py" > "$OUT/integral-diagnostic.log" 2>&1
