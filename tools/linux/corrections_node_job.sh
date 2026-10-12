#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
OUT="$TC/corrections/node-$(hostname)-${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
for mode in Debug Release; do
 "$TC/python/bin/python" "$SRC/tools/linux/arithmetic_environment.py" --build "$TC/build/CompactStar-$mode" --output "$OUT/runtime-$mode.json"
 "$TC/python/bin/python" "$SRC/tests/linux/arithmetic_oracle.py" --executable "$TC/build/CompactStar-$mode/tests/linux_arithmetic_probe" --output "$OUT/arithmetic-$mode.json" > "$OUT/arithmetic-$mode.log"
done
"$TC/python/bin/python" "$SRC/tests/linux/libm_oracle.py" --libm "$TC/sysroot/usr/lib64/libm.so.6" --output "$OUT/libm-oracle.json"
