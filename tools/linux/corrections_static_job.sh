#!/bin/bash
# Static characterization only. No evolution or campaign executable.
set -euo pipefail
SRC=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
OUT="$TC/corrections/static-${SLURM_JOB_ID:?}"
mkdir -p "$OUT/captures"
export PYTHONDONTWRITEBYTECODE=1 COMPACTSTAR_UNIT1_CAPTURE_DIR="$OUT/captures"
cd "$SRC"
for mode in Debug Release; do
 b="$TC/build/CompactStar-$mode"
 set +e
 /usr/bin/time -v ctest --test-dir "$b" -R '^(phase5b_(PB1|PB6|PB7|PB9-11|PB12|PB13|contracts)|phase5d_response|phase5d_independent_oracles|baryon_number_cmf|hartle_monopole_regression|hartle_monopole_physics_cmf)$' --output-on-failure --output-junit "$OUT/structure-$mode.xml" -j1 > "$OUT/structure-$mode.log" 2>&1
 echo "$?" > "$OUT/structure-$mode.exit"
 set -e
 "$TC/python/bin/python" "$SRC/tools/linux/arithmetic_environment.py" --build "$b" --output "$OUT/runtime-$mode.json"
done
b="$TC/build/CompactStar-Debug"
/usr/bin/time -v "$TC/python/bin/python" "$SRC/tests/analysis/chemical_trackr_budget.py" "$b/tests/chemical_trackr_fixture" "$OUT/characterization" > "$OUT/characterization.log" 2>&1
/usr/bin/time -v "$TC/python/bin/python" "$SRC/tests/analysis/chemical_production_validation.py" "$b/tests/chemical_production_fixture" "$OUT/characterization" "$OUT/production" --repeat > "$OUT/production.log" 2>&1
validation=$(find "$OUT/production" -name validation.json)
candidate=$(dirname "$validation")/t8192-r80000/candidate.json
"$TC/python/bin/python" "$SRC/tools/linux/compare_candidate_coefficients.py" --mac "$SRC/tests/baselines/phase5c_chemical_coefficients.json" --linux "$candidate" --validation "$validation" --output "$OUT/coefficient-agreement.json"
# Original Mac frozen certificate parser/currentness is exercised as its own input;
# this does not bind its sensitivities to a Linux star.
"$b/tests/phase6a1_frozen_monitor" /mnt/sdd/zaki/CompactStar/incoming/CompactStar-EKU-handoff-20261010/scientific/phase6/frozen.tsv > "$OUT/mac-certificate-parser.log" 2>&1
nm -C "$b/tests/phase6_breakpoint_qualification" | "$TC/python/bin/python" -c 'import sys,json; print(json.dumps([s.strip() for s in sys.stdin if any(t in s for t in ("TaskManager::","MicroBNVAna::","BNV_Analysis::"))]))' > "$OUT/legacy-symbols.json"
