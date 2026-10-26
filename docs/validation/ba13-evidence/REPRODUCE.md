# Reproduce the evidence candidate

The historical fitter, certificate, source profiles, prior review and X2 evidence
remain immutable inputs under `/mnt/sdd/zaki/CompactStar/linux-transition-20261010`.
The entry authentication binds their 884 manifest entries and historical source.
New numerical raw JSON, outputs, Slurm allocations, receipts, test failures and
successful checks are retained in `ba13-evidence-20261026/` and the compact
`evidence/` directory here. `matrix.json.gz` has all 462 rows and full provenance;
`matrix.tsv.gz` is the tabular readout. Unknown complete error/headroom fields
are null, not zero. `comparison.json.gz` retains normalized reference values and
signed errors. `external-manifest.json` binds the external evidence files.

Inspect without solving:

```sh
gzip -dc docs/validation/ba13-evidence/evidence/matrix.tsv.gz > /tmp/ba13-matrix.tsv
sha256sum -c docs/validation/ba13-evidence/evidence/SHA256SUMS
```

The SHA256SUMS check must be run from this worktree (paths are repository-relative).
Original tests and the 113-check audit remain reproducible from the original
`worktrees/ba13-closure` at 093d121. The new `run.py historical-audit` makes a
separate mirror of read-only inputs and writes new test receipts outside history.
The descendant's `validate.py` independently recomputes comparisons, verifies
historical identities, runs all 43 tests and executes six source mutations.

Fresh reproduction, only when authorized as a new numerical execution, uses an
empty output directory and the toolchain Python. It reuses authenticated retained
executables, never a third solver. Set BA13_EVIDENCE_OUTPUT_ROOT in the submitted
environment and override Slurm stdout/stderr to that directory. Run serially;
inspect each receipt before proceeding:

```sh
export BA13_EVIDENCE_OUTPUT_ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-evidence-replay
mkdir "$BA13_EVIDENCE_OUTPUT_ROOT"
sbatch --output="$BA13_EVIDENCE_OUTPUT_ROOT/slurm-%j.out" --error="$BA13_EVIDENCE_OUTPUT_ROOT/slurm-%j.err" tools/linux/ba13_evidence/run.sbatch run tangent
# Only after PASS_EXPERIMENT, submit the same wrapper with: run freeze
# Require the frozen-holdouts.json hash to match the committed frozen file.
# Then submit serially: run holdouts, run assemble, run historical-audit,
# run validate, run diagnostics. A scientific stop prohibits later stellar work.
```

The completed evidence uses jobs 19130–19137, with the disclosed parse-only failure
19135. A new reproduction does not erase that failure or any historical result.
No script overwrites numerical output. Numerical totals include subprocess CPU
time; build/hash-only packaging is separate. The runner requires Slurm, pins one
core, limits subprocess CPU and checks the cumulative 180 CPU-second reserve.
Slurm accounting storage was disabled; allocation and measured execution receipts
are the available accounting evidence.

The final runner adds an output-directory override and a source manifest to
support replay. The scientific assembly source is unchanged from its execution
commit 2cac8b0; the later reference-module change is the output-directory override
only. Receipt commits identify the exact versions used. All inherited historical
source/evidence blobs remain in Git; only the closure writer's latent nonzero
default is hardened prospectively.

Import the incremental bundle into a repository containing the exact base:

```sh
git bundle verify /path/to/ba13-evidence.bundle
git fetch /path/to/ba13-evidence.bundle candidate/linux-ba13-evidence-20261026
git switch -c review/ba13-evidence FETCH_HEAD
```

This imports an unratified candidate. It does not merge master, promote a
certificate or supply owner acceptance. The prior large evidence inputs are
still required for re-execution; their hashes and locations are in the inherited
manifests and entry-authentication record.
