# Reproduce and transfer the BA13 final decision

Entry: `84e1062db78b91e4fb0988074941b1aea6019917`.
Branch: `candidate/linux-ba13-final-decision-20261030`.
All commands below run from the isolated final-decision worktree unless stated
otherwise. Toolchain Python is the retained cluster Python 3.11.9. No new solver,
physical source, fitter or dependency build is needed.

## Inspect and authenticate

```sh
sha256sum -c docs/validation/ba13-final-decision/evidence/SHA256SUMS
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/python/bin/python -B tools/linux/ba13_final_decision/verify_delivery.py
gzip -dc docs/validation/ba13-final-decision/evidence/matrix.tsv.gz > /tmp/ba13-final-matrix.tsv
```

The delivery verifier requires the inherited external data at the manifest paths.
Every term retains its original source provenance; `source-manifest.json` binds
assessment inputs, run receipts bind executed code, and the final SHA256SUMS binds
the documentation and scripts. `authentication.json` is the entry snapshot. Its
`authenticate_entry.py` script describes the entry-time check; running it after
creating descendants is not an entry replay (its ancestor/clean-tree checks are
intentionally tied to entry). Use the delivery verifier on the finished branch.

The raw 21 MB matrix remains outside Git. Compact deterministic gzip matrix,
20-gap ledger, tangent reassessment, source manifest, tests, mutation failures,
Slurm outputs/allocation receipts, resource and external manifests are committed.
`publish.py` packages existing results and never solves a star.

## Optional fresh static replay

A new output directory is mandatory. No original artifact is overwritten.
These are static assessment/test jobs, not new stellar calculations. This directory
must be shared with the compute node. The existing launcher has the isolated
worktree path fixed explicitly; update that launcher path prospectively if moving
the checkout. Preserve entry/dependency hashes. Use the same single-thread settings.

```sh
export BA13_FINAL_OUTPUT_ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-final-decision-replay
mkdir "$BA13_FINAL_OUTPUT_ROOT"
sbatch --output="$BA13_FINAL_OUTPUT_ROOT/slurm-%j.out" --error="$BA13_FINAL_OUTPUT_ROOT/slurm-%j.err" tools/linux/ba13_final_decision/run.sbatch benchmark
# Inspect receipt-benchmark.json: PASS and <30 CPU s before expansion.
sbatch --output="$BA13_FINAL_OUTPUT_ROOT/slurm-%j.out" --error="$BA13_FINAL_OUTPUT_ROOT/slurm-%j.err" tools/linux/ba13_final_decision/run.sbatch validate
```

The final runner contains the strengthened missing-zero mutation. The historical
first control run and its later correction remain in the delivered receipts.
No need to rerun controls-v2 in a fresh replay: it was the bounded corrective test
run on the original evidence directory. The runner sums prior receipt CPU usage,
requires reserve and limits child CPU. Read process receipts; `sacct` reported
`Slurm accounting storage is disabled` for jobs 19138–19140.

## Office-Mac handoff

A verified incremental bundle and new-evidence archive are provided in
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-final-decision-handoff-20261030/`.
`HANDOFF.json` records the exact delivered head, base, branch, bundle/archive
SHA-256 and verification results. No GitHub publication was attempted; the bundle
is supplied without making network credentials a completion prerequisite.

Copy the handoff directory by the established office-Mac transport, then in a
repository containing the exact entry:

```sh
shasum -a 256 -c SHA256SUMS
git cat-file -t 84e1062db78b91e4fb0988074941b1aea6019917
git bundle verify ba13-final-decision.bundle
git fetch ba13-final-decision.bundle candidate/linux-ba13-final-decision-20261030
git switch -c review/ba13-final-decision FETCH_HEAD
```

If sparse checkout is active, materialize `docs/validation/ba13-final-decision/`
and `docs/validation/ba13-evidence/`. Review the report, methods, decision text and
matrix before any acceptance action. The archive contains only **new** external
evidence; the inherited X1/X2/BA13 data remain separate prerequisites, located and
hashed in their preserved manifests. Extract into a new directory and verify
manifest entries there; do not overwrite historical trees. Import is not a merge,
owner ratification, certificate promotion or permission to run BNV evolution.
