# Source-only execution and read-only reproduction

The command below was submitted exactly once as job 19191, from authenticated
commit `37a9a21a85dfb03d9913d8e39739b294bb1acdd5`. It is an execution record,
not authorization to resubmit. Its exclusive output path now exists, so the
pre-submission guard refuses reuse. A new trajectory needs owner authorization.

```bash
set -euo pipefail
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_SNAPSHOTS="$RF_SOURCE/docs/validation/results-first-bnv-pilot/restart-preparation/snapshots"
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-restart-retry-20261109-r1/source-stage
cd "$RF_SOURCE"
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
test ! -e "$RF_OUTPUT"
mkdir -p "$RF_OUTPUT"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.log" \
  tools/linux/native_qualification/recovery_stage.sbatch
```

The source wrapper fixes one Julian year, source only, one logical execution
CPU, 4 GiB, 30 minutes and single-threaded libraries. No control dependency is
submitted. The retained replacement snapshots authenticate the unchanged stage
binary and inputs, including the approved restart selector.

Read-only assessment (no solver, new trajectory or reference integration):

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
python3 -B tools/linux/native_qualification/assess_source_retry.py \
 /mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-restart-retry-20261109-r1/source-stage \
 /mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-restart-retry-20261109-r1/reassessment
sha256sum -c docs/validation/results-first-bnv-pilot/stage-a-source-retry/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/stage-a-source-retry/evidence/external-SHA256SUMS
python3 -B tools/linux/native_pair/snapshot.py check \
 docs/validation/results-first-bnv-pilot/restart-preparation/snapshots
```

Use a fresh reassessment destination. The assessor authenticates prior trajectory
records against their retained external manifest before byte comparison. Its
110-digit Decimal arithmetic lifts original binary64 values exactly; its isolated
count-arithmetic bound uses 200 digits. This is postprocessing, not a substitute
for the prescribed high-precision chemical ODE reference calculation.

The full assessment and every produced raw file, including incomplete outputs,
remain under /mnt/sdd with external hashes. The compact summary and compressed
per-checkpoint inventory are committed. All old failure receipts remain unchanged.
