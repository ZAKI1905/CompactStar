# Reproduce conditional endpoint accounting

Authority commit: `b7f1209`. Final numerical source:
`b0a080db9f334d6825412e1fb3b102f5bd0dba7b`. The delivery adds evidence and
reports afterward. All numerical work is EKU Slurm; no stellar or BNV evolution
executable is used. The old witness is authenticated input (`inputs.json:1`).

From a clean checkout of this candidate on EKU:

```bash
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-endpoint-reproduction
mkdir "$RF_OUTPUT"
RF_COMMIT=$(git -C "$RF_SOURCE" rev-parse HEAD)
export RF_COMMIT
sbatch --export="ALL,RF_COMMIT=$RF_COMMIT,RF_SOURCE=$RF_SOURCE,RF_OUTPUT=$RF_OUTPUT" \
  --output="$RF_OUTPUT/slurm-%j.out" --error="$RF_OUTPUT/slurm-%j.err" \
  "$RF_SOURCE/tools/linux/results_first/endpoint_acceptance.sbatch"
```

The output directory must be new and outside the repository. This runner checks
clean/exact source, pins one CPU, and enforces 120 CPU seconds/five minutes wall.
It copies materialized historical files into external storage and restores the
original ADR from the starting commit for the historical authority tests. The
old registries and their original checks are not edited. A broad Git archive is
not used because the supplied checkout intentionally omits unrelated Git blobs.

Expected execution: exit 0; 27 tests pass; historical audits 5697/5321 checks;
new independent audit 5804 checks with ten mutation refusals. Expected scientific
result: CONDITIONAL_STATIC_PASS, 420 conditional/21 exact/21 metric consumer-only
rows; all_restricted_pilot_gates_pass=false, nonlinear_trajectory_qualified=false,
actual_path_uncertainty=null. No operational certificate is issued.

The 13-test prerequisite suite can be reproduced with the same environment and
a different fresh RF_OUTPUT using `endpoint_prerequisites.sbatch`. It builds
only the named synthetic/analytic prerequisite executables; it does not invoke
`native_input_binding`, `BuildBoth`, stellar fixtures or any evolution runner.
Request: four CPUs, 8 GiB, 15 minutes wall. Its source/input/artifact identities
are recorded under `evidence/prerequisites/`. Exact binaries and the full build
remain externally at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-endpoint-prerequisites-20261109-r1/build`.

Final assessment output remains at
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-endpoint-acceptance-20261109-r4`.
`evidence/assessment` is copied verbatim; its JSON/gzip output is never hand-edited.
The preceding r1/r2/r3 directories are preserved externally, with compact logs,
receipts and any generated assessment retained in `attempts/`. The copied run
manifest names external historical-copy files as well; it is not a manifest
claim that all external files were duplicated into Git. The package SHA256SUMS
covers the actual delivered files, current source and consumed inputs.

From the repository root, verify delivery without numerical work:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/endpoint-acceptance/SHA256SUMS
```

Do not rerun an old manifest against the amended ADR and then rewrite its pins.
The historical decision is retained in Git and in the external historical test
copy; the new authority registry authenticates the accepted amendment separately.
There is no Mac artifact to stage. Any future substantial Mac-side artifact must
be placed on Scriptorium, as instructed.
