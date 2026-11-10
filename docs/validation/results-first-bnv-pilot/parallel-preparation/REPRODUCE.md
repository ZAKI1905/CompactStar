# 24-worker preparation and prepared source command

The implemented worker mode is `RF_STAGE_WORKERS=24`; omission or `1` retains
the serial path. Other counts are refused by the Stage A entry point. The
24-worker path requires 24 distinct physical cores in the process affinity and
permits only source Stage A. It does not authorize or run the matched control.

## Authorized source command — separate from preparation

This uses the authenticated replacement snapshots published after successful
preparation. Preserve the exclusive output path and stop on every failure.
The owner approved longer runs and the 24-worker implementation; the wrapper
requests 24 physical cores and 96 GiB. The command explicitly overrides its
two-hour default with the five-hour ceiling subsequently authorized by the owner.
The measured 6.06-fold speedup supports a rough 3–4 hour collection estimate.
This resource allowance does not change numerical
settings or convert collection completion into scientific qualification.

```bash
set -euo pipefail
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_SNAPSHOTS="$RF_SOURCE/docs/validation/results-first-bnv-pilot/parallel-preparation/snapshots"
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-parallel-20261110-r1/source-stage
cd "$RF_SOURCE"
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
test ! -e "$RF_OUTPUT"
mkdir -p "$RF_OUTPUT"
sbatch --time=05:00:00 --export=ALL --output="$RF_OUTPUT/slurm-%j.log" \
  tools/linux/native_qualification/parallel_stage.sbatch
```

The complete original matrix remains five references, 12 main trajectories,
12 endpoint-only repeats, four grids per main (240/480/960/1920 intervals),
and three quarter-cap restart tails. Only isolated checkpoint reconstruction
is parallel. The coordinator uses the unchanged ascending-order diagnostics
and `ComputeCheckpointR20`. No tree reduction, reordered sum, trajectory
segmentation, passive-checkpoint restart or altered acceptance condition occurs.

The source still needs actual-output high-precision chemical reference,
reference convergence, main-prefix and accumulated-error bounds, thermal
roundoff/enclosure assessment, full R18/R20 and restart accuracy. The collector
does not supply these proofs automatically. The prior 13-point D1 conditional
PASS remains local to those intervals. Control remains separately authorized
only after full source qualification; Stage B is not authorized.

## Preparation reproduction — no new main or stellar trajectory

Use the implementation commit recorded in `evidence/source-head.txt`, a clean
worktree with authenticated dependencies and a fresh output directory under
/mnt/sdd. This performs a build, synthetic process/transport tests, eight retained
contracts, restart selector/round-trip regression, retained-history reconstruction
and snapshot publication. It does not invoke the Stage A entry point.

```bash
set -euo pipefail
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/parallel-preparation-NEW
cd "$RF_SOURCE"
export RF_COMMIT=$(git rev-parse HEAD)
test ! -e "$RF_OUTPUT"
mkdir -p "$RF_OUTPUT"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.log" \
  tools/linux/native_qualification/parallel_preparation.sbatch
```

The preparation allocation is 24 physical cores/96 GiB/30 minutes, with four
compiler processes. Git checks are read-only and recorded before the build;
byte authentication still covers original dependencies and allows only the two
authorized previously bound changes, `tests/native/stage.cpp` and
`tests/CMakeLists.txt`. The exact transport schema is bound to unchanged
production-header hashes. Every new worker/helper/runner source is bound in the
replacement snapshots. Old snapshot files and failure records are preserved.

Read-only authentication:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/parallel-preparation/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/parallel-preparation/evidence/external-SHA256SUMS
python3 -B tools/linux/native_pair/snapshot.py check \
  docs/validation/results-first-bnv-pilot/parallel-preparation/snapshots
```

Native `packet.bin` files are internal same-build IPC evidence, not portable
restart checkpoints. Their complete binary64 fields are decoded only by the
matching authenticated executable. Scientific TSV/JSON evidence retains its
original format; timing and context-instance labels are bookkeeping, not altered
physical inputs or numerical allowances.
