#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources
OUT="$TC/corrections/cells-${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
GSL="$TC/install/gsl-2.7.1"
for mode in Debug Release; do
 "$CXX" -std=c++17 -O0 -ffp-contract=off -I"$TC/install/$mode/CONFIND/include" -I"$TC/install/$mode/Zaki/include" -I"$GSL/include" "$SRC/CONFIND/tests/linux/cell_probe.cpp" "$TC/install/$mode/CONFIND/lib/libCONFIND.a" "$TC/install/$mode/Zaki/lib/libZaki.a" "$GSL/lib/libgsl.a" "$GSL/lib/libgslcblas.a" "$TC/sysroot/usr/lib64/libz.so" -lm -pthread -o "$OUT/cell-probe-$mode"
 "$TC/python/bin/python" "$SRC/CONFIND/tests/linux/cell_oracle.py" --executable "$OUT/cell-probe-$mode" --output "$OUT/cells-$mode.json" > "$OUT/cells-$mode.log" 2>&1
done
