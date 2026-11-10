# Restart amendment preparation and proposed source command

No Stage A trajectory was executed in this preparation. The command below is
for **separate owner authorization**, not a record of execution. It runs the
source only; it neither submits a control nor authorizes Stage B.

## Exact source-first retry command — not executed

Run from the EKU login environment only after explicit authorization against the
committed preparation receipt and snapshot identities:

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

The retained wrapper fixes `RF_STAGE_ROLE=source` and one Julian year; it requests
one logical execution CPU, 4 GiB and 30 minutes on phyall/phyn004. It uses fresh
output, the authenticated stage ELF and unchanged physical inputs. It does not
auto-submit a matched control (`tools/linux/native_qualification/recovery_stage.sbatch:20`).
The original aggregate Stage A ceiling and complete source qualification before
control remain in force; a later control authorization must account for elapsed
source/assessment time rather than resetting the aggregate allowance.

This is an **evidence-collection command**, not an automatic accuracy certificate.
Selector preparation does not discharge actual-output high-precision chemical
qualification, thermal enclosure/refinement assumptions, main-prefix/accumulated
error, thermal roundoff, global reference convergence, R18/R20, prefix energy
accounting or restart numerical accuracy. These remain required; unknown bounds
are UNQUALIFIED. Source failure stops the task; no control, retuning or resubmission
follows automatically. See `../restart-coverage/RETRY_PLAN.md:1` and
`../stage-a-recovery/REPORT.md:149` for the retained obligations.

## Preparation reproduction — no thermal trajectory

The build/static preparation uses the committed implementation at `f8a1837`.
From that source, with a **fresh** output directory under /mnt/sdd:

```bash
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/restart-preparation-NEW
export RF_COMMIT=$(git -C "$RF_SOURCE" rev-parse HEAD)
mkdir -p "$RF_OUTPUT"
cd "$RF_SOURCE"
sbatch --export=ALL --output="$RF_OUTPUT/build-slurm-%j.log" \
  tools/linux/native_qualification/restart_preparation.sbatch
```

The executed preparation directory is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/restart-preparation-20261109-r1`.
Do not overwrite it. The preparation wrapper builds `native_p2_stage`, the selector
and retained contracts, then runs the selector against the authenticated old
accepted history and its exact binary64 state round-trip. The publisher reuses
the unchanged old static verifier with fresh caches/work and the retained full
archive. It compares all twelve static RHS outputs byte-for-byte, and authenticates
the new ELF's runtime dependencies before publishing both matched snapshots.
No stellar-regeneration target or stage trajectory is invoked. A no-argument
stage smoke check deliberately refuses before input loading/integration.

Read-only checks on the final prepared tree:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/restart-preparation/SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/restart-preparation/evidence/external-SHA256SUMS
python3 -B tools/linux/native_pair/snapshot.py check \
  docs/validation/results-first-bnv-pilot/restart-preparation/snapshots
```

The old snapshots and manifests remain durable historical evidence. Their bound
old stage source is available at the recorded historical commit; they are not
rewritten to authenticate the amended working tree. New snapshot identities are
recorded in `snapshots/identities.json` and `evidence/snapshot-diff.json`.
