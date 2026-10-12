#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources
OUT="$TC/corrections/integrals-${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
GSL="$TC/install/gsl-2.7.1"
for mode in Debug Release; do
 "$CXX" -std=c++17 -O0 -ffp-contract=off -I"$TC/install/$mode/Zaki/include" -I"$GSL/include" "$SRC/Zaki/tests/linux/integral_probe.cpp" "$TC/install/$mode/Zaki/lib/libZaki.a" "$GSL/lib/libgsl.a" "$GSL/lib/libgslcblas.a" "$TC/sysroot/usr/lib64/libz.so" -lm -o "$OUT/integral-probe-$mode"
 "$TC/python/bin/python" "$SRC/Zaki/tests/linux/integral_oracle.py" --executable "$OUT/integral-probe-$mode" --output "$OUT/integral-$mode.json" > "$OUT/integral-$mode.log" 2>&1
done
