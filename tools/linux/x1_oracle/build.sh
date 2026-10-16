#!/bin/bash
set -euo pipefail
SRC=$(cd "$(dirname "$0")/../../.." && pwd)
ROOT=$(dirname "$(dirname "$SRC")")
TC="$ROOT/toolchain"
OUT="$ROOT/x1-oracle-anchor/build"
mkdir -p "$OUT"
for name in onset_structure_oracle x1_regularized_oracle; do
 /usr/bin/time -v -o "$OUT/$name.build-receipt.txt" "$TC/bin/g++" \
  -ffp-contract=off -O3 -DNDEBUG -std=c++17 \
  -isystem "$TC/install/gsl-2.7.1/include" \
  -isystem "$TC/install/Release/Zaki/include" "$SRC/tests/bnv/$name.cpp" \
  -o "$OUT/$name" "$TC/install/Release/Zaki/lib/libZaki.a" \
  "$TC/install/gsl-2.7.1/lib/libgsl.a" "$TC/install/gsl-2.7.1/lib/libgslcblas.a"
done
"$TC/bin/g++" --version > "$OUT/compiler.txt"
