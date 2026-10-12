# Native BA13 static recovery recipe — candidate only

Use the authenticated Linux correction installation at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010`. No Mac numerical certificate
is a native input. Substantial work runs under Slurm; commands below have no
trajectory invocation. The source tree must remain unchanged during a stage.

```bash
cd /mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources/CompactStar
source tools/linux/eku_environment.sh
export PYTHONDONTWRITEBYTECODE=1
OUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-recovery
sbatch --partition=phyall --cpus-per-task=4 --mem=8G --time=01:00:00 \
  --export=ALL tools/linux/ba13_build.sh "$OUT/fresh-build-receipt"
# Wait for successful build, then for each distinct run label:
sbatch --partition=phyall --cpus-per-task=1 --mem=8G --time=06:00:00 \
  --export=ALL tools/linux/ba13_job.sh stars RUN_LABEL
# Wait for stars/COMPLETE.json, or use --dependency=afterok:STAR_JOB.
sbatch --partition=phyall --cpus-per-task=1 --mem=8G --time=01:00:00 \
  --export=ALL tools/linux/ba13_job.sh native-adapters RUN_LABEL
# Wait for native-adapters/COMPLETE.json; each array worker owns its files.
sbatch --partition=phyall --cpus-per-task=1 --mem=8G --time=02:00:00 \
  --array=0-20%8 --export=ALL tools/linux/ba13_job.sh coefficient RUN_LABEL
# Wait for every coefficient-0..20/COMPLETE.json before collecting.
"$TC/python/bin/python" -B tools/linux/ba13_collect.py \
  "$OUT/RUN_LABEL" "$OUT/RUN_LABEL/raw"
```

The characterizing input is the native current-source Phase-5C refinement run
`toolchain/corrections/static-14487/characterization/run-_kg8zl8c/result.json`,
SHA-256 `d70bc57a15fe2b1aaf242c9976d609ec9b7219f8f4377685744a55f79c62bb93`.
Its `t8192-r80000/freegas.tsv` SHA-256 is
`f6b9c4068a430cc356b54e39c21a0ca17e17b8c2b35d17b836638298d1112d7f`.
This label is the governed EOS-resolution parameter, not a claim that the final
onset-augmented EOS has only 8192 rows. The current table has 10375 data rows.

The original adapter (`ba13_job.sh adapters`) was actually run and refused the
inherited Debug B0 tail witness on the independently regenerated Release profile.
Its failure is retained. The native adapter recomputes the tail witness for every
actual profile through unchanged `chemical_trackr_budget.tail` before invoking
unchanged `chemical_production_evidence.certificates` and `transport`. It retains
all existing assertions, inherited center/background estimates, grids and goals.
It makes no claim that a B0 background estimate is independent per-star refinement.

The original C++ producer is compiled with zero source changes. The recovered
historical snapshots are under `historical/`; source identities are in
`recovered-identities.json`. The adapter uses the original `thermal_source`
function directly, so entropy mathematics is unchanged.

The final fitter is now recovered as the unchanged
`tests/bnv/compare_candidate.py` (see `fitter-identity.json`). Its historical
`merge-frozen` mode merges the original shards; the native collector performs the
same index-ordered concatenation with added completion/hash checks. The archived
final assembled tables and all shards are retained for verification.

Verify the original ZIP SHA-256
`990d3708a24c71e518d4bb59e6ad3c23ae47f337ba924007b083ca598e6b8e2b`,
extract into a fresh directory, and verify complete `SHA256SUMS` coverage before
using the package. Do not write fit outputs into the historical tree.

After collecting two independent native runs and building the validation targets:

```bash
# Run on the login node; the environment must be exported into Slurm.
source tools/linux/eku_environment.sh
ARCHIVE="$OUT/historical-package/BA13-original-evidence-20261012"
sbatch --partition=phyall --cpus-per-task=1 --mem=4G --time=00:30:00 \
  --export=ALL tools/linux/ba13_qualify.sh \
  "$OUT/run1/raw" "$OUT/repeat/raw" "$ARCHIVE" "$OUT/FRESH_QUALIFICATION"
```

The recipe runs the original fitter through atomic, hash-bound candidate output
publication; checks exact repeated certificate bytes; independently refits both
native and historical certificates using separate 80/112-digit normal equations;
compares all component errors; exercises the existing depletion/M20 monitor; and
binds the native candidate to its actual canonical native profile and exact B0.
It contains no trajectory command. The commands were executed individually in
jobs 14573–14577; the consolidated shell recipe preserves those stages. Job 14573
completed all fits, then its monitor build failed because modules were sourced on
the compute node. Job 14575 completed that build/control with the proper exported
login-node environment; both logs are retained.

A successful script exit reports completed static processing, **not** scientific
qualification. The platform report remains UNQUALIFIED for the historical error
assignments described in proposed ADR-0021. Do not run further certificate-dependent
pilots or promote these files. The separate provisional replacement draft was
never selected or executed.

Independent input-boundary controls use `phase6a1_ba13_profile_binding PROFILE
CERTIFICATE FRESH_WORK`: nominal native exact-profile equality, one-ULP profile
mutation refusal, and changed-certificate byte-currentness refusal through the
existing Phase-5 fixture and owners. These controls attach no sensitivity table
and invoke no trajectory.
