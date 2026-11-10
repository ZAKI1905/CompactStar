# Reproduce Stage A recovery evidence

The single thermal retry authorization was consumed by source job 19188 and
stopped on required restart coverage. These records do not authorize another
trajectory or a changed restart/cap rule.

Source implementation: c5d59b1 (two O2 literals, collector evidence/role boundary,
configuration test). Prelaunch setup fix: a909285 (publisher path typo).
Authenticated snapshot publication before execution: 509053d. The prior
4387e56 D1 resolution and all historical snapshots/failures remain preserved.

Actual job sequence, using the retained toolchain environment:

```sh
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-recovery-20261109-r1
export RF_COMMIT=c5d59b1e00fd731c08d77fe490ac187a743bf916
# 19186: build and seven tests passed; publisher-path typo stopped setup.
sbatch --export=ALL --output="$RF_OUTPUT/build-slurm-%j.log" tools/linux/native_qualification/recovery_build.sbatch
# 19187: after the committed typo fix, static verification and publication only.
sbatch --export=ALL --output="$RF_OUTPUT/static-slurm-%j.log" tools/linux/native_qualification/recovery_verify.sbatch
# 19188: authenticated source-only Stage A; no automatic control.
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-recovery-20261109-r1/source-stage
export RF_SNAPSHOTS="$RF_SOURCE/docs/validation/results-first-bnv-pilot/stage-a-recovery/snapshots"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.log" tools/linux/native_qualification/recovery_stage.sbatch
```

The existing output directories are evidence and must not be overwritten.
The build does not execute a stellar-regeneration target. All three allocations
were one logical execution CPU, 4 GiB, 30-minute ceilings on phyall/phyn004.

The source stopped by explicit `scancel 19188` after the read-only restart proof
had been persisted. No quarter-cap member was run: the proof establishes its
identical accepted history because every planned cap leaves all recorded GSL
proposals and targets unchanged. The collector's eventual restart refusal was
not falsely recorded as an observed native exception.

Read-only reproduction, from the repository root (no Slurm or solver needed):

```sh
run=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-recovery-20261109-r1/source-stage
python3 -B tools/linux/native_qualification/assess_recovery_restart.py \
  "$run/run/trajectory/source/tier-0-cap-0/main.accepted_states.tsv" \
  /tmp/stage-a-restart-recheck.json
python3 -B tools/linux/native_qualification/assess_stage_a_recovery.py \
  "$run" /tmp/stage-a-assessment-recheck
python3 -B tools/linux/native_qualification/summarize_stage_a_scientific_context.py \
  "$run" /tmp/stage-a-scientific-recheck
```

Compare the generated JSON/TSV with `evidence/{restart-gate,assessment,
scientific-context}.json` and `evidence/checkpoints.tsv`. The summary verifies
old reference bytes against the original external manifest and compares new
main/reference histories with their prior counterparts. Decimal is used only
for reporting exact lifts and arithmetic on retained records; no new IVP is
executed by these commands.

```sh
sha256sum -c docs/validation/results-first-bnv-pilot/stage-a-recovery/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/stage-a-recovery/evidence/external-SHA256SUMS
python3 -B tools/linux/native_pair/snapshot.py check \
  docs/validation/results-first-bnv-pilot/stage-a-recovery/snapshots
```

The original snapshot identities still identify the previous source tree and
binaries; they were not rewritten to authenticate the newly modified source.
New snapshots bind the approved source and new binaries. Physical input/archive
identities and runtime linkage are unchanged. No actual-output chemical,
thermal-roundoff, global-reference or full ledger uncertainty bound is implied
by small raw spreads or by the earlier 13-point D1 PASS.
