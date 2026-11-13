#!/bin/bash
set -euo pipefail
: "${COMP_SOURCE:?}" "${COMP_OUTPUT:?}"
TC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain
export PATH="$TC/bin:/opt/spack/linux-rocky9-zen2/gcc-11.4.1/cmake-3.29.4-gybrkgtp3fscoavg562fdufw33t7wnpu/bin:/opt/spack/linux-rocky9-zen2/gcc-11.4.1/gmake-4.4.1-cosp3uuqxxbzc6tea5rotv6umxpvom7q/bin:$PATH"
export TMPDIR="$COMP_OUTPUT/tmp" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$TMPDIR"
cd "$COMP_SOURCE"
scontrol show job "$SLURM_JOB_ID" > "$COMP_OUTPUT/build-allocation-$SLURM_JOB_ID.txt"
cmake -S "$COMP_SOURCE" -B "$COMP_OUTPUT/build" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS=-ffp-contract=off \
  -DGSL_ROOT_DIR="$TC/install/gsl-2.7.1" \
  -DZLIB_INCLUDE_DIR="$TC/sysroot/usr/include" -DZLIB_LIBRARY="$TC/sysroot/usr/lib64/libz.so" \
  -DPython3_EXECUTABLE="$TC/python/bin/python" \
  -DCOMPACTSTAR_LINUX_DEPENDENCY_LOCK="$TC/dependency-lock-linux.json" \
  -DCOMPACTSTAR_ZAKI_PREFIX="$TC/install/Release/Zaki" \
  -DCOMPACTSTAR_CONFIND_PREFIX="$TC/install/Release/CONFIND" \
  -DCOMPACTSTAR_QUALIFICATION_DEBUG_ZAKI_PREFIX="$TC/install/Debug/Zaki" \
  -DCOMPACTSTAR_QUALIFICATION_DEBUG_CONFIND_PREFIX="$TC/install/Debug/CONFIND" \
  -DCOMPACTSTAR_EOS_DATA_ROOT=/mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/compose \
  > "$COMP_OUTPUT/configure.log" 2>&1
cmake --build "$COMP_OUTPUT/build" --parallel "$SLURM_CPUS_PER_TASK" --target \
  bnv_competition_contract bnv_competition_short phase6a1_direct_energy phase6a1_thermal_ledger \
  passive_checkpoint_output_contract phase6_breakpoint_contract phase5d_component_tolerances \
  phase5d_response phase6_currentness_bytes > "$COMP_OUTPUT/build.log" 2>&1
ctest --test-dir "$COMP_OUTPUT/build" --output-on-failure --parallel "$SLURM_CPUS_PER_TASK" \
  -R '^(bnv_competition_contract|phase6a1_direct_energy|phase6a1_thermal_ledger|passive_checkpoint_output_contract|phase6_breakpoint_contract|phase5d_component_tolerances|phase5d_response|phase5d_independent_oracles|phase6_currentness_bytes)$' \
  > "$COMP_OUTPUT/contracts.log" 2>&1
