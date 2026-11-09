# Native snapshot verification and separately proposed thermal commands

Run on EKU, with source and all outputs on /mnt/sdd. No regeneration is required
for assembly or verification. The single authorized producer run was job19166;
its output is immutable input for all commands below. Do not resubmit
regenerate.sbatch under this delivery's already-consumed authorization.

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE="$PWD"
export RF_SNAPSHOTS="$RF_SOURCE/docs/validation/results-first-bnv-pilot/native-regeneration/snapshots"
export RF_RETAINED=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-regeneration-20261108-r1
python3 -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
python3 -B tools/linux/native_pair/test_snapshot.py "$RF_SNAPSHOTS"
```

For a fresh static reconstruction, use a fresh directory and Slurm (no evolution):

```bash
export RF_OUTPUT="$RF_RETAINED" RF_BUILD="$RF_RETAINED/build"
export RF_COMMIT=$(git rev-parse HEAD)
export RF_VERIFY_WORK="$RF_RETAINED/review-static-work"
export RF_VERIFY_TSV="$RF_RETAINED/review-static-rhs.tsv"
sbatch --export=ALL --output="$RF_RETAINED/review-static-%j.log" \
  tools/linux/native_pair/verify.sbatch
```

The script authenticates producer/input bytes, requires clean committed source,
loads capture/sequence.archive with its pinned SHA256 and never calls a solver.
The original producer authentication and numerical receipts are in evidence/.
To rebuild source use build.sbatch with RF_OUTPUT pointing at a fresh build
output root (4CPU/8GiB/15min). A newly linked binary has a new content identity;
retain the published binaries for exact delivered-snapshot replay. Do not
silently overwrite snapshot binary hashes after rebuilding.

The commands below are prepared, **NOT EXECUTED**, and require separate owner
thermal-trajectory authorization. They run the unchanged production evolution
and checkpoint machinery through tests/native/stage.cpp. The wrapper authenticates
both snapshots and their complete dependency manifests before numerical work.
It cannot select the full historical campaign: the only accepted bounds are1 and2048years.

Stage A, 1 Julian year, maximum30min aggregate, one CPU/4GiB:

```bash
export RF_YEARS=1
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r1
mkdir -p "$RF_OUTPUT"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.log" tools/linux/native_pair/stage.sbatch
```

After Stage A's actual-path assessment passes and owner authorization covers Stage B:

```bash
export RF_YEARS=2048
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-2048yr-r1
mkdir -p "$RF_OUTPUT"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.log" tools/linux/native_pair/stage.sbatch
```

No dependency submits Stage B automatically. RUN_PLAN.md and the preserved
native-assembly/NONLINEAR_PLAN.md specify the full maximum matrix, error norms,
uncertainty terms and stop rules. The collector's completion marker is not an
accuracy pass; actual-path assessment and unknown-error refusal remain required.
