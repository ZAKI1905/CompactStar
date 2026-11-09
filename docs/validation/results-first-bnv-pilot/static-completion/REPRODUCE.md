# Reproduce ADR-0022 static completion

Use the retained EKU dependency closure and the candidate worktree. The numerical
source commit is `6a4e109951df63e9d2a1002ad026fa99b6758331`; later delivery commits
add only documentation/evidence. Do not overwrite any retained input or output.
Choose a fresh RF_OUTPUT under /mnt/sdd. The source must be clean.

```bash
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-static-reproduction
export RF_COMMIT=$(git -C "$RF_SOURCE" rev-parse HEAD)
mkdir "$RF_OUTPUT"
sbatch --export=ALL --output="$RF_OUTPUT/slurm-%j.out" --error="$RF_OUTPUT/slurm-%j.err" \
  "$RF_SOURCE/tools/linux/results_first/static_completion.sbatch"
```

The job runs the three test suites, static assessment and independent serialized
ledger/mutation audit. The 60 CPU-second aggregate ceiling is enforced. No
stellar executable is invoked. Expected scientific exit is **2**, after retaining
`BLOCKED_UNSUPPORTED_CUT_RESPONSE`; all tests and independent audit must exit 0.
Check receipt.json rather than treating the scheduler's FAILED label as an
execution crash or treating successful tests as scientific qualification.

The source manifest authenticates 30 retained data files. The authority registry
authenticates the model, owner record, ADR and original approval/failure evidence.
The preserved ledger contains all 462 historical rows plus separate new records.
Inspect `assessment/matrix.tsv.gz`, `assessment/summary.json`,
`assessment/precision-witness.json`, and `assessment/independent-ledger-audit.json`.
The complete 353-row expectation list is
`assessment/amended-expectation-violations.json`.

Verify the delivered files from this directory with:

```bash
sha256sum -c SHA256SUMS
```

`evidence/external-final-manifest.json` binds both original external run directories,
including final controller receipts. Evidence is copied byte-for-byte; generated
artifacts are never hand-edited. In-memory mutation controls do not touch source,
archive or evidence files. Mac copies, if later requested, must go directly to a
verified writable Scriptorium mount; no Mac mount is claimed observable here.

The older `../PREPARATION_SHA256SUMS` binds the proposed ADR at the starting
commit, before acceptance. It is preserved unchanged. Its ADR entry intentionally
does not match today's accepted ADR; `evidence/historical-preservation.json`
verifies every old entry against the original Git blob and identifies that sole
authorized current-file difference. Use this directory's SHA256SUMS for this
delivery, not a rewritten historical manifest.
