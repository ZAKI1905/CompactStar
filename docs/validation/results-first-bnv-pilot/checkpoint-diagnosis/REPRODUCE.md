# Reproduce the diagnosis without another trajectory

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE="$PWD"
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-checkpoint-diagnosis-20261109-r1
export RF_REPLAY_RETAINED=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r2
export RF_RETAINED=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-regeneration-20261108-r1
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check \
 docs/validation/results-first-bnv-pilot/native-runtime-retry/snapshots
sha256sum -c docs/validation/results-first-bnv-pilot/checkpoint-diagnosis/evidence/external-SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/checkpoint-diagnosis/SHA256SUMS
cat "$RF_OUTPUT/replay/output/checkpoint-1/failure.txt"
"$TC/python/bin/python" -B tools/linux/native_qualification/assess_checkpoint_failure.py \
 "$RF_OUTPUT" "$TC/../ba13-recovery/run1/raw/coefficients.tsv" \
 "$RF_OUTPUT/review-assessment.json" \
 "$RF_REPLAY_RETAINED/run/trajectory/source/reference-4-observations.tsv"
cmp "$RF_OUTPUT/assessment.json" "$RF_OUTPUT/review-assessment.json"
```

The assessment reads data and evaluates the already implemented quadratic energy
formula at already recorded states using 100-digit Decimal arithmetic. It calls
no numerical executable, RHS, stellar solve, checkpoint reconstruction or trajectory.
The finest-grid tolerance estimate is a necessary linearized budget check, not an
authorized numerical setting or additional replay.

The following are the **historical executed commands**, not permission to repeat
the consumed single trajectory authorization. Both scripts require a fresh build/
replay directory; the recorded command environment is above.

```bash
# Job 19175: compile two reporting executables and run synthetic regression tests.
sbatch --parsable --export=ALL --output="$RF_OUTPUT/build-slurm-%j.log" \
 tools/linux/native_qualification/diagnostic_build.sbatch
# After build/tests passed, job 19176: the one source-only diagnostic replay.
sbatch --parsable --export=ALL --output="$RF_OUTPUT/replay-slurm-%j.log" \
 tools/linux/native_qualification/diagnostic_replay.sbatch
```

The build manifest contains exact compile/link commands and hashes. It reuses
the original numerical static libraries, original flags (`-O3 -DNDEBUG`, C++17,
`-ffp-contract=off`) and original actual ELF linkage. `native_main` is extracted
verbatim from the original authenticated collector. No protected source file or
snapshot is replaced. Build output and large objects remain external.

Failure reporting is reusable via `tests/diagnostics/checkpoint_evidence.hpp`.
For this diagnosis it is used in the separate `checkpoint_replay.cpp` collector,
leaving the historical snapshot-bound collector intact. Every returned checkpoint
is serialized before acceptance is checked; no failed channel is dropped. The
main checkpoint candidate is O2 reconstructed from the retained main bracket,
not a nonexistent exact main sample at a strict-interior time. JSON explicitly
records that distinction. TSV rows carry all primary/reference values and budget
contributions. Human-readable failure text is flushed before termination.

The reporting test deliberately induces an Echem diagnostic failure and independent
state failures using the existing synthetic context; exact serialization and
unchanged refusal status are asserted. The saved test outputs are in
`evidence/regression/`. Main replay identity is independently checked against the
original accepted-state TSV; local RHS logs distinguish stage evaluations from
accepted steps. A recurrence returns exit 1 and never constitutes qualification.

Machine-readable nonfinite utilization in a generic failed state may be `inf` in
TSV, while JSON's unbounded main step cap is `null` (the main policy is uncapped).
The actual Echem failure record is finite. Raw Slurm receipts retain their original
whitespace. No original failure, certificate or conditional assumption is promoted.
