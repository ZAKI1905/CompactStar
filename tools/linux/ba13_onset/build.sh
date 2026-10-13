#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-onset-pilot
OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-onset-pilot
OLD=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty/build
mkdir -p "$OUT/build"
"$TC/python/bin/python" -B "$SRC/tools/linux/ba13_uncertainty/extract_helpers.py" "$SRC/tests/bnv/frozen_sensitivity.cpp" "$OUT/build/ba13_original_helpers.hpp"
for name in onset_pilot onset_table_audit onset_structure_oracle; do
 "$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" "$CXX" -ffp-contract=off -O3 -DNDEBUG -std=c++17 -I"$OUT/build" -I"$OLD/generated/include" -I"$SRC" -isystem "$TC/install/gsl-2.7.1/include" -isystem "$TC/install/Release/Zaki/include" -isystem "$TC/install/Release/CONFIND/include" "$SRC/tests/bnv/$name.cpp" -o "$OUT/build/$name" "$OLD/libCompactStar.a" "$TC/install/Release/CONFIND/lib/libCONFIND.a" "$TC/install/Release/Zaki/lib/libZaki.a" "$TC/sysroot/usr/lib64/libz.so" "$TC/install/gsl-2.7.1/lib/libgsl.a" "$TC/install/gsl-2.7.1/lib/libgslcblas.a"
done
