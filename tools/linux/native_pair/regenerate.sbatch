#!/bin/bash
#SBATCH --job-name=native-p2-sequence
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_BUILD:?}" "${RF_COMMIT:?}" "${TC:?}"
cd "$RF_SOURCE"
test "$(git rev-parse HEAD)" = "$RF_COMMIT"
test -z "$(git status --porcelain --untracked-files=no)"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/sequence-allocation.txt"
python3 tools/linux/native_pair/authenticate.py > "$RF_OUTPUT/sequence-authentication.json"
# The executable itself refuses an existing capture directory. Never retry a producer run.
"$TC/python/bin/python" -I -B tools/linux/measure.py "$RF_BUILD/tests/native_p2_regenerate" \
 "$TC/../ba13-recovery/run1/raw/chemical-profiles/target-0" "$RF_OUTPUT/capture" \
 > "$RF_OUTPUT/sequence.log" 2>&1
