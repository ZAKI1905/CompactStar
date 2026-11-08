#!/bin/bash
#SBATCH --job-name=rf-bnv-native-binding
#SBATCH --partition=phy01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${TC:?}"
export TMPDIR="$RF_OUTPUT/tmp"
mkdir -p "$TMPDIR"
test "$(git -C "$RF_SOURCE" rev-parse HEAD)" = "$RF_COMMIT"
test -z "$(git -C "$RF_SOURCE" status --porcelain --untracked-files=no)"
build="$RF_OUTPUT/build"
native="$TC/../ba13-recovery/run1"
cmake -S "$RF_SOURCE" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DPython3_EXECUTABLE="$TC/python/bin/python" \
  -DCOMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json" \
  -DCOMPACTSTAR_ZAKI_PREFIX="$TC/install/Release/Zaki" \
  -DCOMPACTSTAR_CONFIND_PREFIX="$TC/install/Release/CONFIND" \
  -DCOMPACTSTAR_EOS_DATA_ROOT=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/compose \
  > "$RF_OUTPUT/configure.log" 2>&1
/usr/bin/time -v -o "$RF_OUTPUT/build.time" cmake --build "$build" \
  --target phase6a1_native_input_binding -j4 > "$RF_OUTPUT/build.log" 2>&1
/usr/bin/time -v -o "$RF_OUTPUT/native-binding.time" \
  "$build/tests/phase6a1_native_input_binding" \
  "$native/raw/chemical-profiles/target-0" \
  "$native/raw/chemical-profiles/target-0/thermal" \
  "$RF_SOURCE/docs/validation/results-first-bnv-pilot/preparation-registry.json" \
  "$native/raw/coefficients.tsv" "$RF_OUTPUT/native-fixture" \
  > "$RF_OUTPUT/native-binding.log" 2>&1
