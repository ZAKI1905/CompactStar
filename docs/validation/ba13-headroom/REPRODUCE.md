# Reproduction and handoff

Candidate only. The implemented assessment recommends B; it does not qualify BA13.
Source entry: `89b629dd9961c29e127badbd723455d00dc38553`.

The external artifact root is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-headroom`.
`evidence/external-manifest.json` authenticates every retained external output,
including failed diagnostics, binaries, command/resource receipts, nominal-profile
summaries and the uncompressed 462-row matrix. `COMPLETE.json` is a completion
marker for an assessment with explicit gaps, not a scientific qualification.

For the frozen original evidence use the existing X2 and BA13 recovery manifest
identities. No regeneration is needed. From the X2 worktree:

```sh
ROOT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
"$ROOT/toolchain/python/bin/python" -B tests/bnv/test_ba13_headroom.py
"$ROOT/toolchain/python/bin/python" -B tools/linux/ba13_headroom/audit.py
```

The audit rechecks 236 nominal manifest files and 257 X2 external artifacts,
all 462 historical thresholds/drifts, exact primary replay, the complete matrices
and direct-event checks. It is a completion/provenance check, not independent review.
Running the audit changes external `audit.json`/`tests.txt`; therefore verify the
manifest BEFORE regenerating those receipts. To preserve this retained package,
copy the external root first and adapt the explicit root in a separate checkout.

`run.py` and `run.sbatch` document all commands for the existing X2 route. Solves
are restricted to targets 0,10,20. `profile.py` and `supplement.py` read the existing
nominal profiles; they do not solve new stars. Numerical outputs are write-once.
Build phases are `build`, `build-v2`, `build-final`; the distinct executable names
retain each build. Do not rerun these phases into the retained artifact directory.
Use a separate artifact root and update the declared root in a separate checkout
for a full repeat. All solver inputs are in receipts, including exact Hc/B strings.

Read the committed full matrix without any recomputation:

```sh
gzip -dc docs/validation/ba13-headroom/evidence/matrix.tsv.gz
gzip -dc docs/validation/ba13-headroom/evidence/matrix.json.gz
```

The JSON preserves evidence types, assumptions, missing estimates, allocations,
normalization and correlation handling. A finite `partial_utilization` is not an
upper bound or completed total. `proposed_total_utilization` is null for every
UNQUALIFIED row. Historical fields have never been overwritten.

The office-Mac handoff directory is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/ba13-headroom-handoff-20261016`.
It contains an incremental bundle based on the full starting SHA, its SHA-256,
a verified bundle receipt and a compact artifact archive. On the existing Mac
repository with the starting commit, verify the SHA files and `git bundle verify`,
then fetch the candidate into a review branch. Do not merge master or promote a
certificate. External raw artifacts remain available through the authenticated
cluster manifest; the compact archive contains the matrices, pilot/profile/event
summaries and failure evidence for review. Next action: one independent review
of this construction and its explicit residual gaps.
