# Reproduce the targeted retained-data assessment

The numerical execution source is
`c7a406803f57603b7c623404acfcd0eba7cc7115`; the final delivery changes only
documentation/evidence afterward. Use the clean candidate worktree on EKU.
All raw values required by the new calculation are in `retained-inputs.json.gz`.
`input-provenance.json` identifies and hashes 85 original source files; the
snapshot preserves decimal tokens. `inputs-manifest.json` authenticates the
snapshot and consumed existing Git-tracked evidence. Source constants are the
binary64 Zaki masses and conversion constants used by X2, explicitly widened
in `tools/linux/results_first/cut_response.py:31`.

Run with a fresh external output directory:

```bash
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-cut-reproduction
mkdir "$RF_OUTPUT"
RF_COMMIT=$(git -C "$RF_SOURCE" rev-parse HEAD)
export RF_COMMIT
sbatch --export="ALL,RF_COMMIT=$RF_COMMIT,RF_SOURCE=$RF_SOURCE,RF_OUTPUT=$RF_OUTPUT" \
  --output="$RF_OUTPUT/slurm-%j.out" --error="$RF_OUTPUT/slurm-%j.err" \
  "$RF_SOURCE/tools/linux/results_first/cut_response.sbatch"
```

The runner refuses dirty/wrong source, pins one CPU, runs the tests, evaluates
70/90-digit responses and checks a COPY of the existing ledger outside the
repository. Its child CPU allowance is at most the remaining part of 60 seconds.
The seven tests include independent momentum integrals, EOS derivatives, terminal
redshift, inverse-matrix response, an implicit-root path-tangent problem, and
refusals of invalid convergence estimates. The old ledger is not modified by its
auditor. Execution success means the assessment ran correctly; `cut_closed=false`
is still the scientific result. No matrix or operational certificate is generated.

`evidence/response.json.gz` is a deterministic gzip copy of the execution's
response.json, with no numerical editing. `evidence/summary.json` and receipts are
verbatim copies. The external manifest records original output hashes. Failed
and preliminary attempt logs remain in `evidence/prior/`; their source commits
are retained in Git. `SHA256SUMS` binds delivered source, inputs and evidence.
