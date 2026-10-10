#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources
EOS=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/compose/DS-CMF-1-with-crust
GSL="$TC/install/gsl-2.7.1"
mkdir -p "$TC/historical"
for kind in fixed stress; do
 name=oracle.tsv
 if [ "$kind" = stress ]; then name=oracle-stress.tsv; fi
 gzip -dc "$SRC/Zaki/tests/fp_preservation/fixtures/$name.gz" > "$TC/historical/reference-$kind.tsv"
done
for mode in Debug Release; do
 zaki="$TC/install/$mode/Zaki"
 "$CXX" -std=c++17 -O0 -g -ffp-contract=off -I"$zaki/include" -I"$GSL/include" "$SRC/Zaki/tests/fp_preservation/characterize.cpp" "$zaki/lib/libZaki.a" "$GSL/lib/libgsl.a" "$GSL/lib/libgslcblas.a" "$TC/sysroot/usr/lib64/libz.so" -lm -o "$TC/historical/zaki-$mode"
 for kind in fixed stress; do
  "$TC/historical/zaki-$mode" "$TC/historical/Zaki-$mode-$kind.tsv" "$kind" "$EOS" 'DS(CMF)-1_with_crust.eos' > "$TC/logs/Zaki-$mode-$kind-producer.log" 2>&1
  "$TC/python/bin/python" "$SRC/CompactStar/tools/linux/compare_hex_tables.py" "$TC/historical/reference-$kind.tsv" "$TC/historical/Zaki-$mode-$kind.tsv" --output "$TC/logs/Zaki-$mode-$kind-differences.json"
 done
 for dep in Zaki CONFIND; do
  ar t "$TC/install/$mode/$dep/lib/lib$dep.a" > "$TC/logs/$dep-$mode-archive-members.txt"
  nm --defined-only "$TC/install/$mode/$dep/lib/lib$dep.a" | gzip > "$TC/logs/$dep-$mode-symbols.txt.gz"
 done
done
