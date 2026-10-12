#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty
OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-uncertainty
GSL="$TC/install/gsl-2.7.1"
cmake -S "$SRC" -B "$OUT/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_FLAGS=-ffp-contract=off -DCMAKE_CXX_FLAGS=-ffp-contract=off -DPython3_EXECUTABLE="$TC/python/bin/python" -DGIT_EXECUTABLE="$TC/bin/git" -DGSL_INCLUDE_DIR="$GSL/include" -DGSL_LIBRARY="$GSL/lib/libgsl.a" -DGSL_CBLAS_LIBRARY="$GSL/lib/libgslcblas.a" -DZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" -DZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so" -DCOMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json" -DCOMPACTSTAR_ZAKI_PREFIX="$TC/install/Release/Zaki" -DCOMPACTSTAR_CONFIND_PREFIX="$TC/install/Release/CONFIND" -DCOMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$TC/install/Debug/Zaki" -DCOMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$TC/install/Debug/CONFIND" -DCOMPACTSTAR_EOS_DATA_ROOT=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/compose -DBUILD_TESTING=ON > "$OUT/configure.log" 2>&1
"$TC/python/bin/python" -B "$SRC/tools/linux/measure.py" cmake --build "$OUT/build" --target ba13_uncertainty_pilot -j4 > "$OUT/build.log" 2>&1
