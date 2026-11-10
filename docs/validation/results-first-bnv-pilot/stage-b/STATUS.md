# Stage B handoff — 2026-11-10

The owner requested that monitoring stop and will notify the agent when the
Slurm jobs finish. The read-only monitor was stopped; scientific jobs were left
running. This is a checkpoint, not final Stage B qualification.

Repository: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
Branch: `candidate/results-first-bnv-pilot-20261108`.
Verification code: `6016c0a2790c438bb1ac1cc0ea710109636b8164`, descended from
the requested `c90b2205e553afcd102c7acd5f84b253e0a67f1b`.
The authenticated worktree was clean apart from untracked Python bytecode.
Earlier branch/worktree ancestry checks are retained in the campaign evidence.

## Completed and pending

Both original-t0 2,048-year thermal evolutions and all ten independent native
references completed. The source's 48 diagnostic grids and source audit completed.
The final monitor observation at 16:51:32 recorded 41/48 control grids complete:

| Job | Last observed state | Work |
|---|---|---|
| `19353_4` | RUNNING, phyn004 | Control REFINED diagnostics |
| `19353_5` | RUNNING, phyn004 | Control ULTRA diagnostics |
| `19352` | PENDING, dependency | Control audit |
| `19356` | PENDING, dependency | Combined numerical assessments |

Native tier/audit allocations use 48 logical CPUs, 24 physical workers, 48 GiB,
and unlimited partition walltime. Final job states must be inspected on resume;
the table is a historical snapshot, not a live status display.

The working [report](REPORT.md) records the preliminary matched endpoint
temperature difference, +15,325,803.20420946 K, and two actual source Cstar knot
crossings. The candidate isolated chemical-heating temperature contribution is
only +4.936 to +5.615 microkelvin. Endpoint reaction terms oppose about 12–15%
of the chemical driving, despite negligible thermal heating. These findings do
not replace completion of the original numerical gates or domain review.

## Resume after the owner reports completion

1. Inspect terminal states, exit codes, completion receipts, and scientific guards
   for the four jobs above. Preserve genuine failures and incomplete gates.
2. Inspect the results of `postprocess.sbatch` assessment job `19356`: residuals,
   observation residuals, physical assessment, chemical-response bounds, events,
   ledger, local references, accuracy, and coverage. Job completion alone does not
   establish that every scientific gate passes.
3. Run `matched.py`, then `qualification.py` against
   `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-20261110-r3`.
   These final aggregation steps are not included in job `19356`.
4. Generate and inspect final plots with `plot.py` using the existing
   `native-p2-plotting-env`; collect resources with `resource_report.py` across
   `stage-b-20261110-r1`, `stage-b-20261110-r2`, and `stage-b-20261110-r3`.
5. After all writers finish, run `seal.py` into a fresh seal directory, verify
   its two-pass manifest receipt, and run `collect_evidence.py`. Final plots,
   sealed manifests, and the complete numerical qualification are still pending.
6. Update the final report with the actual gates, limitations, resources, and
   exact commits; commit compact evidence and documentation. Do not rerun the
   original trajectories merely to produce the report.

Scripts are under `tools/linux/stage_b/` in the repository above. Heavy artifacts
remain under `/mnt/sdd`; no new horizon or stellar campaign is authorized.

## Checkpoint validation and authority

This checkpoint changes documentation only; no scientific output changes.
Authority remains GOVERNANCE §§1–5, ADR-0014 through ADR-0017 and Amendment 1,
and INV-15, as cited in the working report. Seven piecewise known-answer tests,
five Apply-residual tests, the audit-source equivalence test, and five standalone
report imports passed before this checkpoint. Historical failures and operational
refusals remain preserved. Full-matrix qualification and final manifest checks
are pending, not claimed as passed. New verification remains an unverified
scientific candidate under GOVERNANCE §5; Stage A owner acceptance does not
ratify Stage B. No new scientific ambiguity is resolved by this handoff.
