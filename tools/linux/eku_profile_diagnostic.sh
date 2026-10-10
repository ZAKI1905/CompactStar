#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
HANDOFF=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010
GSL="$TC/install/gsl-2.7.1"
mkdir -p "$TC/profile-diagnostic"
"$CXX" -std=c++17 -O0 -g -ffp-contract=off -I"$SRC" -I"$TC/build/CompactStar-Debug/generated/include" -I"$TC/install/Debug/Zaki/include" -I"$TC/install/Debug/CONFIND/include" -I"$GSL/include" "$SRC/tools/linux/profile_identity.cpp" "$TC/build/CompactStar-Debug/libCompactStar.a" "$TC/install/Debug/CONFIND/lib/libCONFIND.a" "$TC/install/Debug/Zaki/lib/libZaki.a" "$GSL/lib/libgsl.a" "$GSL/lib/libgslcblas.a" "$TC/sysroot/usr/lib64/libz.so" -pthread -lm -o "$TC/profile-diagnostic/profile-identity"
cd "$SRC"
"$TC/profile-diagnostic/profile-identity" "$HANDOFF/scientific/phase6/profile_tree" "$TC/profile-diagnostic/work" "$TC/profile-diagnostic/profile.tsv"
"$TC/build/CompactStar-Debug/tests/phase6a1_frozen_monitor" "$HANDOFF/scientific/phase6/frozen.tsv" > "$TC/logs/frozen-monitor-authenticated.log" 2>&1
