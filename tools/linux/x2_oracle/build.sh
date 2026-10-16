#!/bin/bash
set -euo pipefail
SRC=$(cd "$(dirname "$0")/../../.." && pwd)
ROOT=$(dirname "$(dirname "$SRC")")
TC="$ROOT/toolchain"
OUT="$ROOT/x2-fixed-step/build"
mkdir -p "$OUT"
for mode in extended double; do
 flags=()
 if [[ "$mode" == double ]]; then flags+=(-DX2_DOUBLE); fi
 /usr/bin/time -v -o "$OUT/$mode.build.txt" "$TC/bin/g++" -O3 -std=c++17 -ffp-contract=off -Wall -Wextra "${flags[@]}" -isystem "$TC/install/Release/Zaki/include" "$SRC/tests/bnv/x2_fixed_step.cpp" -o "$OUT/$mode" "$TC/install/Release/Zaki/lib/libZaki.a"
done
"$TC/bin/g++" --version > "$OUT/compiler.txt"
