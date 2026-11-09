#!/bin/bash
#SBATCH --job-name=native-p2-retry-preflight
#SBATCH --partition=phyall
#SBATCH --nodelist=phyn004
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=4G
#SBATCH --time=00:15:00
# Non-evolution checks only. No automatic thermal submission.
set -euo pipefail
: "${RF_SOURCE:?}" "${RF_OUTPUT:?}" "${RF_SNAPSHOTS:?}" "${RF_RETAINED:?}" "${TC:?}"
cd "$RF_SOURCE"
scontrol show job "$SLURM_JOB_ID" > "$RF_OUTPUT/preflight-allocation.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TC/python/bin/python" -I -B tools/linux/native_qualification/runtime_audit.py \
  "$RF_SNAPSHOTS" > "$RF_OUTPUT/compute-runtime-audit.json"
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
"$TC/python/bin/python" -B tools/linux/native_pair/test_runtime.py
"$TC/python/bin/python" -B tools/linux/native_pair/test_snapshot.py "$RF_SNAPSHOTS"
"$TC/python/bin/python" -I -B tools/linux/measure.py \
  "$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py verify \
  "$RF_SNAPSHOTS" "$RF_OUTPUT/static-work" > "$RF_OUTPUT/static.log" 2>&1
cmp "$RF_OUTPUT/static-work/static-rhs.tsv" "$RF_RETAINED/static-rhs.tsv"
"$TC/python/bin/python" -I -B tools/linux/measure.py ctest --test-dir "$RF_RETAINED/build" --output-on-failure \
  -R '^(phase6a1_direct_energy|phase6a1_thermal_ledger|passive_checkpoint_output_contract|phase6_breakpoint_contract|phase5d_component_tolerances|phase6_currentness_bytes)$' \
  > "$RF_OUTPUT/contracts.log" 2>&1
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
printf 'PASS replacement compute-node bindings, static ownership/RHS, baseline identity, and retained contracts; zero thermal trajectories\n'
