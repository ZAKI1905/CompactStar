#!/bin/bash
# Run under Slurm after sourcing eku_environment.sh. Reuses authenticated Linux
# correction Release build/lock; records all commands through the supplied log.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
OUT=${1:?fresh build receipt directory}
mkdir "$OUT"
cd "$SRC"
cmake -S "$SRC" -B "$TC/build/CompactStar-Release" > "$OUT/configure.log" 2>&1
"$TC/python/bin/python" -B tools/linux/measure.py cmake --build "$TC/build/CompactStar-Release" --target phase6a1_frozen_sensitivity phase6a1_ba13_profile_binding phase6a1_ba13_native_certificate_binding phase6a1_frozen_monitor -j4 > "$OUT/build.log" 2>&1
cp "$TC/build/CompactStar-Release/tests/CMakeFiles/phase6a1_frozen_sensitivity.dir/flags.make" "$OUT/flags.make"
cp "$TC/build/CompactStar-Release/tests/CMakeFiles/phase6a1_frozen_sensitivity.dir/link.txt" "$OUT/link.txt"
cp "$TC/build/CompactStar-Release/dependency-provenance.json" "$OUT/dependency-provenance.json"
"$TC/python/bin/python" -B tools/linux/arithmetic_environment.py --build "$TC/build/CompactStar-Release" --output "$OUT/runtime.json"
