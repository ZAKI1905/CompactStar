#!/bin/bash
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
HANDOFF=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010
PY="$TC/python/bin/python"
export PYTHONDONTWRITEBYTECODE=1
"$PY" -c 'import numpy, scipy, mpmath; print(numpy.__version__, scipy.__version__, mpmath.__version__)'
"$PY" "$SRC/tools/linux/record_dependencies.py" --toolchain "$TC" --sources "$(dirname "$SRC")" --output "$TC/dependency-lock-linux.json"
GSL="$TC/install/gsl-2.7.1"
common=(-DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_FLAGS=-ffp-contract=off -DCMAKE_CXX_FLAGS=-ffp-contract=off -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DPython3_EXECUTABLE="$PY" -DGIT_EXECUTABLE="$TC/bin/git" -DGSL_ROOT_DIR="$GSL" -DGSL_INCLUDE_DIR="$GSL/include" -DGSL_LIBRARY="$GSL/lib/libgsl.a" -DGSL_CBLAS_LIBRARY="$GSL/lib/libgslcblas.a" -DZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" -DZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so" -DCOMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json" -DCOMPACTSTAR_EOS_DATA_ROOT="$HANDOFF/scientific/compose" -DBUILD_TESTING=ON -DCOMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$TC/install/Debug/Zaki" -DCOMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$TC/install/Debug/CONFIND")
for mode in Debug Release; do
 build="$TC/build/CompactStar-$mode"
 cmake -S "$SRC" -B "$build" -DCMAKE_BUILD_TYPE="$mode" -DCOMPACTSTAR_ZAKI_PREFIX="$TC/install/$mode/Zaki" -DCOMPACTSTAR_CONFIND_PREFIX="$TC/install/$mode/CONFIND" "${common[@]}" > "$TC/logs/CompactStar-$mode-configure.log" 2>&1
 cmake --build "$build" --parallel 4 -- -k > "$TC/logs/CompactStar-$mode-build.log" 2>&1
 ctest --test-dir "$build" --show-only=json-v1 > "$TC/logs/CompactStar-$mode-inventory.json"
 set +e
 ctest --test-dir "$build" -R '^(chemical_production_contract|chemical_exact_oracles|chemical_curved_gc9|particle_number_analytic|passive_checkpoint_output_contract|phase6a1_direct_energy|phase6a1_thermal_ledger|heat_capacity_v1|cache_contract|cache_thermal_contract|phase5d_component_tolerances|phase6_breakpoint_contract|phase6_currentness_bytes)$' --output-on-failure --output-junit "$TC/logs/CompactStar-$mode-readiness.xml" -j1 > "$TC/logs/CompactStar-$mode-readiness.log" 2>&1
 rc=$?
 set -e
 echo "$rc" > "$TC/logs/CompactStar-$mode-readiness.exit"
done
