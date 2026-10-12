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

The final fit stage is **on hold for the owner's located Scriptorium archive**.
The drafts retained outside Git under
`ba13-recovery/provisional-fit-not-used/` are
provisional replacement drafts, not a selected reproduction of the historical
fit. Do not use their output as a qualified frozen certificate. Audit archived
raw tables, final JSON/covariance and fitting implementation before finalizing
stage 4 or running any certificate-dependent checks.

Independent input-boundary controls use `phase6a1_ba13_profile_binding PROFILE
CERTIFICATE FRESH_WORK`: nominal native exact-profile equality, one-ULP profile
mutation refusal, and changed-certificate byte-currentness refusal through the
existing Phase-5 fixture and owners. These controls attach no sensitivity table
and invoke no trajectory.
