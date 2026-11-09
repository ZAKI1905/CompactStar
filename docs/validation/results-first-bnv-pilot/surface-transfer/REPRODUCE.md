# Reproduce this bounded local calculation

Numerical source commit: `30e24f488a6b0dbd5f1b7c24129256c67db8e184`.
The final delivery adds only documentation and copied evidence after that
commit. Use a clean EKU checkout of the source or delivery. All required
stellar states are already retained in Git; this calculation does not invoke
the stellar solver. `inputs.json` binds the consumed retained data, law,
accepted registry and prospective criteria. Their source provenance is in
`../cut-assessment/input-provenance.json`; no decimal inputs were regenerated.

From the clean candidate checkout, submit to EKU with a new external output
directory (the example name must not already exist):

```bash
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-surface-transfer-reproduction
mkdir "$RF_OUTPUT"
RF_COMMIT=$(git -C "$RF_SOURCE" rev-parse HEAD)
export RF_COMMIT
sbatch --export="ALL,RF_COMMIT=$RF_COMMIT,RF_SOURCE=$RF_SOURCE,RF_OUTPUT=$RF_OUTPUT" \
  --output="$RF_OUTPUT/slurm-%j.out" --error="$RF_OUTPUT/slurm-%j.err" \
  "$RF_SOURCE/tools/linux/results_first/surface_transfer.sbatch"
```

The pinned toolchain's Python supplies numpy and mpmath; the same environment
was used for the preceding cut assessment. The runner refuses a dirty/wrong
source, requires Slurm, pins one CPU, and enforces a 120-second aggregate CPU
ceiling. The allocation requests one CPU, 1 GiB and five minutes wall. It runs
22 tests, the 90/120-digit maps, six variations, endpoint diagnostics, and both
independent accounting/refusal audits. The historical auditor writes only to
an external copy of the historical ledger. No canonical old file is mutated.

Expected execution result: exit 0. Expected scientific result:
`FINITE_INTERIOR_TRANSFER_SUPPORTED_ENDPOINT_EXTENSION_UNSUPPORTED`,
`cut_closed=false`, with 399 approximate rows still unqualified. A passing
execution is not a passing C-CUT qualification. The unauthorized endpoint
scenario is labeled in every affected row. It must not be used as a runtime
certificate or an evolution authorization.

Final evidence is a verbatim copy of
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-surface-transfer-20261109-r3`.
The directory suffix uses the UTC date; Slurm's local job date is November 8.
The controller final record was added after the runner's manifest was written.
All assessment files are covered by the runner manifest; the delivery
SHA256SUMS additionally covers controller/log records and all delivered files.
JSON gzip output is deterministic (mtime=0); no numerical output is edited.

`attempts/19155` and `attempts/19156` preserve the failed/preliminary runs in
full. Their source commits are `3dd687003b1a6af8308843b4e79f2fe4664ebc6d` and
`67d8174179a6d34547e5e4636e473ca2edc90070`. Do not substitute the preliminary
ledger for the final ledger: its old numbers were serialized as strings, and
the final run includes the independent typed-record preservation check.

From the repository root, verify the complete delivery without running a
numerical calculation:

```bash
sha256sum -c docs/validation/results-first-bnv-pilot/surface-transfer/SHA256SUMS
```

The manifest excludes itself. It covers the new source/tests/documentation,
copied evidence and every input named in inputs.json. Historical package
manifests remain unchanged.
