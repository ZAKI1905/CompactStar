#!/bin/bash
# Submit with --cpus-per-task=4 --mem=8G --time=01:00:00 after eku_environment.sh.
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
OUT="$TC/corrections/${SLURM_JOB_ID:?}"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
export COMPACTSTAR_UNIT1_CAPTURE_DIR="$OUT/captures"
mkdir -p "$COMPACTSTAR_UNIT1_CAPTURE_DIR"
readiness='chemical_production_contract|chemical_exact_oracles|chemical_curved_gc9|particle_number_analytic|passive_checkpoint_output_contract|phase6a1_direct_energy|phase6a1_thermal_ledger|heat_capacity_v1|cache_contract|cache_thermal_contract|phase5d_component_tolerances|phase6_breakpoint_contract|phase6_currentness_bytes'
for mode in Debug Release; do
 b="$TC/build/CompactStar-$mode"
 /usr/bin/time -v cmake --build "$b" --parallel 4 > "$OUT/build-$mode.log" 2>&1
 ctest --test-dir "$b" --show-only=json-v1 > "$OUT/inventory-$mode.json"
 set +e
 ctest --test-dir "$b" -R "^($readiness|rotochemical_trackr_pe|relativistic_unit_boundary|candidate_arithmetic_oracle|rotochemical_local_thermodynamics|rotochemical_trackr_freegas_local|rotochemical_trackr_npe|phase5d_protected_manifest|phase5d_coupled_oracles)$" --output-on-failure --output-junit "$OUT/tests-$mode.xml" -j1 > "$OUT/tests-$mode.log" 2>&1
 echo "$?" > "$OUT/tests-$mode.exit"
 set -e
 "$TC/python/bin/python" "$SRC/tests/linux/arithmetic_oracle.py" --executable "$b/tests/linux_arithmetic_probe" --output "$OUT/arithmetic-$mode.json" > "$OUT/arithmetic-$mode.log"
done
