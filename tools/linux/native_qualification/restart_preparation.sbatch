#!/bin/bash
#SBATCH --job-name=restart-preparation
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# Build and static verification only; no stellar regeneration or thermal run.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_COMMIT:?}" "${TC:?}"
cd "$RF_SOURCE"
export TMPDIR="$RF_OUTPUT/tmp" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/build-allocation.txt"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=threads bash -c '
 set -euo pipefail
 "$TC/python/bin/python" -I -B -c "import os; a=sorted(os.sched_getaffinity(0)); print(a); assert len(a)==1" > "$RF_OUTPUT/build-affinity.txt"
 test "$(git rev-parse HEAD)" = "$RF_COMMIT"
 test -z "$(git status --porcelain --untracked-files=no)"
 cmake -S "$RF_SOURCE" -B "$RF_OUTPUT/build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS=-ffp-contract=off \
  -DGSL_ROOT_DIR="$TC/install/gsl-2.7.1" \
  -DZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" -DZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so" \
  -DPython3_EXECUTABLE="$TC/python/bin/python" \
  -DCOMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json" \
  -DCOMPACTSTAR_ZAKI_PREFIX="$TC/install/Release/Zaki" \
  -DCOMPACTSTAR_CONFIND_PREFIX="$TC/install/Release/CONFIND" \
  -DCOMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$TC/install/Debug/Zaki" \
  -DCOMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$TC/install/Debug/CONFIND" \
  -DCOMPACTSTAR_EOS_DATA_ROOT=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/compose > "$RF_OUTPUT/configure.log" 2>&1
 "$TC/python/bin/python" -B tools/linux/native_pair/authenticate.py > "$RF_OUTPUT/authentication.json"
 cmake --build "$RF_OUTPUT/build" --target native_p2_stage native_o2_configuration native_restart_schedule \
  phase6a1_direct_energy phase6a1_thermal_ledger passive_checkpoint_output_contract \
  phase6_breakpoint_contract phase5d_component_tolerances phase6_currentness_bytes -j1 > "$RF_OUTPUT/build.log" 2>&1
 ctest --test-dir "$RF_OUTPUT/build" --output-on-failure \
  -R "^(native_restart_schedule|native_o2_configuration|phase6a1_direct_energy|phase6a1_thermal_ledger|passive_checkpoint_output_contract|phase6_breakpoint_contract|phase5d_component_tolerances|phase6_currentness_bytes)$" > "$RF_OUTPUT/contracts.log" 2>&1
 "$RF_OUTPUT/build/tests/native_restart_schedule" \
  "$RF_SOURCE/docs/validation/results-first-bnv-pilot/stage-a-recovery/evidence/main.accepted_states.tsv" > "$RF_OUTPUT/selector-roundtrip.log" 2>&1
 exec "$TC/python/bin/python" -B tools/linux/native_qualification/restart_preparation_snapshots.py "$RF_OUTPUT"
' > "$RF_OUTPUT/build-stage.log" 2>&1
