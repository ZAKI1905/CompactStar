# Reproduction and authorization boundaries

All paths below are on EKU under `/mnt/sdd`. No numerical producer rebuild or
stellar regeneration is involved. This records the single authorized retry;
it is not permission for another retry after failure.

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE="$PWD"
export RF_RETAINED=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-regeneration-20261108-r1
export RF_SNAPSHOTS="$RF_SOURCE/docs/validation/results-first-bnv-pilot/native-runtime-retry/snapshots"
"$TC/python/bin/python" -B tools/linux/native_pair/test_runtime.py
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py check "$RF_SNAPSHOTS"
"$TC/python/bin/python" -B tools/linux/native_pair/test_snapshot.py "$RF_SNAPSHOTS"
```

Publication was performed at repair commit
`66419b813c15ebade0e2f1bc431666a3b4f07364` with the following command into the
then-absent destination; publication refuses an existing directory. The resulting
snapshots were committed in `ab01d7a`. Publication at a later HEAD changes the
recorded assembly commit and hence snapshot hashes; use the committed snapshots
for exact replay.

```bash
"$TC/python/bin/python" -B tools/linux/native_pair/snapshot.py publish \
  "$RF_RETAINED" "$RF_SNAPSHOTS"
```

Compute prerequisites, job 19173, completed successfully. The script authenticates
replacement manifests and dynamic linkage, tests refusals, performs a fresh static
assembly, compares its twelve RHS states exactly with the retained baseline, and
runs six retained contracts. It never calls thermal integration.

```bash
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-runtime-retry-20261109-r1
mkdir -p "$RF_OUTPUT"
"$TC/python/bin/python" -I -B tools/linux/native_qualification/runtime_audit.py \
  "$RF_SNAPSHOTS" > "$RF_OUTPUT/login-runtime-audit.json"
sbatch --parsable --export=ALL --output="$RF_OUTPUT/preflight-%j.log" \
  tools/linux/native_qualification/retry_preflight.sbatch
```

After all prerequisites passed, Stage A was submitted once, job 19174:

```bash
export RF_YEARS=1
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r2
mkdir -p "$RF_OUTPUT"
sbatch --parsable --export=ALL --output="$RF_OUTPUT/slurm-%j.log" \
  tools/linux/native_qualification/retry_stage.sbatch
```

The wrapper uses `srun --cpu-bind=threads`, checks that its affinity contains
exactly one logical CPU, then executes the unchanged historical stage script.
The snapshot loader rechecks byte hashes and actual runtime resolution before
invoking the unchanged stage ELF. No script submits Stage B automatically.
Stage B, if and only if every applicable Stage A gate passes, uses RF_YEARS=2048
and a fresh output directory, independently initialized from the same snapshot.
The full campaigns are excluded.

Stage A failed checkpoint diagnostic qualification. Stage B was not submitted.
The following is offline analysis of retained source data, with no solver, RHS,
interpolation or checkpoint replay. Use a fresh analysis output directory to avoid
overwriting the retained figures/assessment authenticated by the manifest.

```bash
export RF_STAGE_A=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r2
export MPLCONFIGDIR=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-plotting-env/matplotlib-config
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-plotting-env/bin/python -B \
  tools/linux/native_qualification/summarize_retry.py "$RF_STAGE_A" "$RF_STAGE_A/review-analysis"
sha256sum -c docs/validation/results-first-bnv-pilot/native-runtime-retry/evidence/external-SHA256SUMS
sha256sum -c docs/validation/results-first-bnv-pilot/native-runtime-retry/SHA256SUMS
```

The separate plotting environment uses matplotlib 3.9.4 and numpy 1.26.4; its full
package versions are `evidence/plotting-requirements.txt`. It was created under
`/mnt/sdd` with the existing Python's `-m venv` and installed without pip's cache.
It is not part of the numerical snapshot or scientific runtime. PNG/PDF plots
are diagnostic-only; timestamps in regenerated PDFs need not match original file
hashes. The retained source data and scientific summary are the reproducible
numerical evidence. Raw Slurm receipts preserve whitespace; other source/docs
pass Git whitespace checks.
