#!/bin/sh
# Reproduce the recovered original-stack Debug build in a NEW directory.
set -eu
: "${MAC_RECOVERY_BUILD:?Set an absolute fresh build directory}"
test ! -e "$MAC_RECOVERY_BUILD"
CS_SOURCE=$(git rev-parse --show-toplevel)
PREFIX=/Users/keeper/Documents/CompactStar/external/qualification/compactstar-migration/e263a6e-b0cbd510/Debug
cmake -S "$CS_SOURCE" -B "$MAC_RECOVERY_BUILD" -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DBUILD_TESTING=ON \
  -DPython3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3.12 -DGSL_ROOT_DIR=/opt/local \
  -DCOMPACTSTAR_ZAKI_PREFIX="$PREFIX/Zaki" -DCOMPACTSTAR_CONFIND_PREFIX="$PREFIX/CONFIND"
cmake --build "$MAC_RECOVERY_BUILD" --parallel 4 --target \
  chemical_production_contract particle_number_analytic passive_checkpoint_output_contract \
  phase6a1_direct_energy phase6a1_thermal_ledger heat_capacity_v1 cache_contract \
  cache_thermal_contract phase5d_component_tolerances phase6_breakpoint_contract phase6_currentness_bytes
ctest --test-dir "$MAC_RECOVERY_BUILD" --output-on-failure -j1 \
  --output-junit "$MAC_RECOVERY_BUILD/mac-tests.xml" \
  -R '^(chemical_production_contract|chemical_exact_oracles|chemical_curved_gc9|particle_number_analytic|passive_checkpoint_output_contract|phase6a1_direct_energy|phase6a1_thermal_ledger|heat_capacity_v1|cache_contract|cache_thermal_contract|phase5d_component_tolerances|phase6_breakpoint_contract|phase6_currentness_bytes)$'
