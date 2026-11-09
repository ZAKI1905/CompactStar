#!/bin/bash
#SBATCH --job-name=native-p2-build
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${TC:?}"
export TMPDIR="$RF_OUTPUT/tmp"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/allocation-$SLURM_JOB_ID.txt"
test "$(git -C "$RF_SOURCE" rev-parse HEAD)" = "$RF_COMMIT"
test -z "$(git -C "$RF_SOURCE" status --porcelain --untracked-files=no)"
cd "$RF_SOURCE" # retained fixture authenticates provider sources relative to the repository
build="${RF_BUILD:-$RF_OUTPUT/build}"

cmake -S "$RF_SOURCE" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS=-ffp-contract=off \
  -DGSL_ROOT_DIR="$TC/install/gsl-2.7.1" \
  -DZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" \
  -DZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so" \
  -DPython3_EXECUTABLE="$TC/python/bin/python" \
  -DCOMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json" \
  -DCOMPACTSTAR_ZAKI_PREFIX="$TC/install/Release/Zaki" \
  -DCOMPACTSTAR_CONFIND_PREFIX="$TC/install/Release/CONFIND" \
  -DCOMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$TC/install/Debug/Zaki" \
  -DCOMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$TC/install/Debug/CONFIND" \
  -DCOMPACTSTAR_EOS_DATA_ROOT=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/compose \
  > "$RF_OUTPUT/configure.log" 2>&1
python3 tools/linux/native_pair/authenticate.py > "$RF_OUTPUT/authentication.json"
"$TC/python/bin/python" -I -B "$RF_SOURCE/tools/linux/measure.py" cmake --build "$build" \
  --target native_p2_regenerate native_p2_verify native_p2_stage phase6a1_direct_energy phase6a1_thermal_ledger passive_checkpoint_output_contract phase6_breakpoint_contract phase5d_component_tolerances phase6_currentness_bytes -j4 > "$RF_OUTPUT/build.log" 2>&1
