#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources
OUT="$TC/corrections/dependencies-${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
for mode in Debug Release; do
 "$TC/python/bin/python" "$SRC/CONFIND/tests/linux/qualify.py" --build "$TC/build/CONFIND-$mode" --output "$OUT/CONFIND-$mode" > "$OUT/confind-$mode.log" 2>&1
 set +e
 ctest --test-dir "$TC/build/Zaki-$mode" --output-on-failure > "$OUT/zaki-$mode.log" 2>&1
 echo "$?" > "$OUT/zaki-$mode.exit"
 ctest --test-dir "$TC/build/CONFIND-$mode" -R '^(ownership_inputs_exceptions|threading_discrimination|historical_export)$' --output-on-failure > "$OUT/confind-native-$mode.log" 2>&1
 echo "$?" > "$OUT/confind-native-$mode.exit"
 set -e
done
"$TC/python/bin/python" "$SRC/CompactStar/tests/linux/libm_oracle.py" --libm "$TC/sysroot/usr/lib64/libm.so.6" --output "$OUT/libm-oracle.json"
for name in zaki_axis_oracle zaki_axis_attrib confind_pow; do
 "$TC/python/bin/python" "$SRC/CompactStar/docs/validation/linux-corrections/independent-review-20261011/experiments/scripts/$name.py" > "$OUT/$name.txt" 2>&1
done
