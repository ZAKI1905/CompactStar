# Reproduction and office-Mac review handoff

This package is a completed assessment with disposition E, not a completed
scientific certificate. Start with [REPORT.md](REPORT.md), [METHOD.md](METHOD.md),
[OWNER_DECISIONS.md](OWNER_DECISIONS.md) and [dispositions.json](evidence/dispositions.json).

The cluster source root is
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-closure,
branch candidate/linux-ba13-closure-20261019. External evidence is at
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-closure-20261019.
The historical inputs and original independent review remain at their original
locations. Exact source, input and external output identities are committed in
[source-hashes.json](evidence/source-hashes.json), [inputs.json](evidence/inputs.json)
and [external-manifest.json](evidence/external-manifest.json).

To inspect the full matrix without rerunning any calculation:

```sh
gzip -dc docs/validation/ba13-closure/evidence/matrix.tsv.gz
gzip -dc docs/validation/ba13-closure/evidence/matrix.json.gz
```

The matrix's component source pointers refer to the normalized P/R/X and
original physical values in external assessment-final/component-data.json.gz.
All original 13 thermal knots and 2049/8196 profile components are retained.
Every UNQUALIFIED row has proposed_total_utilization=null and named refusal
reasons. Do not interpret partial_utilization as a bound.

The numerical commands are entirely in [run.py](../../../tools/linux/ba13_closure/run.py).
Execution history:

| Slurm job | Action | Result |
|---|---|---|
| 18963 | srun visibility probe | Cancelled by this task before numerical work; subsequent sbatch route worked |
| 18964 | sbatch filesystem/toolchain probe | Completed |
| 18965 | Build and first N=512 replay | Schema gate stopped; retained decimals equal, four additive zero fields omitted by old schema |
| 18966 | Corrected schema gate and N=1024 benchmark | PASS; retired-estimator negative control refuses as required |
| 18967 | Remaining five fine central and six third-width runs | PASS execution |
| 18971 | First static M1–M5 assembly | Completed; original derived outputs retained at artifact root |
| 18977 | Final static assembly | Completed; surface-node localization, stricter order gate and sufficient no-switch conditions recorded under assessment-final |
| 18978 | All BA13 tests and initial independent audit | PASS: 30 tests, 111 checks |

Per-command CPU/wall/return-code/stdout/stderr hashes, process RSS, CPU affinity,
build commands, source hashes, Slurm allocations and final states are in the
external manifest. [resources.json](evidence/resources.json) totals all metered
new numerical phases, including the stopped replay and all static assemblies.
The initial historical authentication/fitter replay and login-node packaging
were not metered, so this is not an exact whole-session CPU claim. No 900-second
initial or 1800-second total numerical expansion was needed.

To run tests alone on the authenticated cluster:

```sh
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
cd "$ROOT/worktrees/ba13-closure"
"$ROOT/toolchain/python/bin/python" -B tests/bnv/test_ba13_headroom.py
"$ROOT/toolchain/python/bin/python" -B tests/bnv/test_ba13_recovery.py
"$ROOT/toolchain/python/bin/python" -B tests/bnv/test_ba13_closure.py
```

To reproduce numerical artifacts, first verify all input/output manifests.
The runner and assessor use write-once output files. Select a SEPARATE output
directory through BA13_CLOSURE_OUTPUT_ROOT; no input or source edit is needed.
Do not rerun into the retained package. For example, from the candidate checkout:

```sh
export BA13_CLOSURE_OUTPUT_ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-closure-independent-repeat
mkdir -p "$BA13_CLOSURE_OUTPUT_ROOT"
sbatch --wait --output="$BA13_CLOSURE_OUTPUT_ROOT/benchmark-%j.out" --error="$BA13_CLOSURE_OUTPUT_ROOT/benchmark-%j.err" tools/linux/ba13_closure/run.sbatch run benchmark
sbatch --wait --output="$BA13_CLOSURE_OUTPUT_ROOT/m1-%j.out" --error="$BA13_CLOSURE_OUTPUT_ROOT/m1-%j.err" tools/linux/ba13_closure/run.sbatch run m1
sbatch --wait --output="$BA13_CLOSURE_OUTPUT_ROOT/static-%j.out" --error="$BA13_CLOSURE_OUTPUT_ROOT/static-%j.err" tools/linux/ba13_closure/run.sbatch assess
sbatch --wait --output="$BA13_CLOSURE_OUTPUT_ROOT/audit-%j.out" --error="$BA13_CLOSURE_OUTPUT_ROOT/audit-%j.err" tools/linux/ba13_closure/run.sbatch audit
```

Each command must complete successfully before the next. This is a reproduction
recipe, not authorization to resume a stopped scientific extension. The retained
schema-stop correction used benchmark-v2; fresh reproductions use the corrected
benchmark gate directly. Execution receipts retain their exact earlier source
hashes. The later output-directory configuration changes I/O only; no equation,
scoring criterion, limit or numerical input changes. Publishing copies existing
evidence and verifies manifests; it does not run science.

No authenticated GitHub push was established; no push or merge was attempted.
A verified incremental bundle and compact review archive are provided at
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-closure-handoff-20261019.
The bundle is based on EXACT entry
286aafb63002093b3052e23a3972903011f36035. The handoff's HANDOFF.json records final
branch/SHA and SHA256SUMS authenticates the bundle and review archive. On the
office Mac, in the existing repository containing that starting commit:

```sh
shasum -a 256 -c SHA256SUMS
git bundle verify ba13-closure.bundle
git fetch ./ba13-closure.bundle candidate/linux-ba13-closure-20261019:refs/heads/review/ba13-closure-20261019
git show --stat review/ba13-closure-20261019
```

Review the candidate without merging canonical master, promoting a certificate,
or treating the candidate ADR interpretation as owner-ratified. Large raw
profiles/certificates stay outside Git; their byte identities are bound by the
input and external manifests. The original independent review remains read-only.
