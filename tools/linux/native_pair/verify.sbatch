#!/bin/bash
#SBATCH --job-name=native-p2-static
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_BUILD:?}" "${RF_COMMIT:?}" "${TC:?}"
cd "$RF_SOURCE"
test "$(git rev-parse HEAD)" = "$RF_COMMIT"
test -z "$(git status --porcelain --untracked-files=no)"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/verify-allocation.txt"
python3 tools/linux/native_pair/authenticate.py > "$RF_OUTPUT/verify-authentication.json"
native="$TC/../ba13-recovery/run1"
read -r archive_sha archive_name < "$RF_OUTPUT/capture/SHA256SUMS"
"$TC/python/bin/python" -I -B tools/linux/measure.py "$RF_BUILD/tests/native_p2_verify" \
 "$native/raw/chemical-profiles/target-0" "$native/raw/chemical-profiles/target-0/thermal" \
 "$RF_SOURCE/docs/validation/results-first-bnv-pilot/preparation-registry.json" \
 "$native/fit/certificate.tsv" "$native/raw/coefficients.tsv" \
 "$RF_OUTPUT/capture/sequence.archive" "$archive_sha" \
 "$RF_OUTPUT/verify-work" "$RF_OUTPUT/static-rhs.tsv" > "$RF_OUTPUT/verify.log" 2>&1
"$TC/python/bin/python" -I -B tools/linux/measure.py ctest --test-dir "$RF_BUILD" --output-on-failure \
 -R '^(phase6a1_direct_energy|phase6a1_thermal_ledger|passive_checkpoint_output_contract|phase6_breakpoint_contract|phase5d_component_tolerances|phase6_currentness_bytes)$' \
 > "$RF_OUTPUT/contracts.log" 2>&1
