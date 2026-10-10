#!/bin/bash
set -u
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
export PYTHONDONTWRITEBYTECODE=1
export COMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$TC/install/Debug/Zaki"
export COMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$TC/install/Debug/CONFIND"
export COMPACTSTAR_QUALIFICATION_GSL_PREFIX="$TC/install/gsl-2.7.1"
export COMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json"
export ZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" ZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so"
for mode in Debug Release; do
 ctest --test-dir "$TC/build/CompactStar-$mode" -E '^phase5d1_controlled_evolution_regression$' --output-on-failure --output-junit "$TC/logs/CompactStar-$mode-suite.xml" --timeout 900 -j1 > "$TC/logs/CompactStar-$mode-suite.log" 2>&1
 rc=$?
 echo "$rc" > "$TC/logs/CompactStar-$mode-suite.exit"
done
