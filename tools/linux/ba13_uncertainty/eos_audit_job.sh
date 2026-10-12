#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
OUT=${BA13_STUDY_OUT:-/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty}
cmake -S "$SRC" -B "$OUT/build" > "$OUT/eos-audit-configure.log" 2>&1
"$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" cmake --build "$OUT/build" --target ba13_uncertainty_eos_audit -j2 > "$OUT/eos-audit-build.log" 2>&1
"$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$OUT/build/tests/ba13_uncertainty_eos_audit" > "$OUT/eos-ordering.tsv" 2> "$OUT/eos-ordering.log"

"$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/profile_diagnostic.py" > "$OUT/profile-diagnostic.log" 2>&1
