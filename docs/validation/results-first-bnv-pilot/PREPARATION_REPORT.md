# Restricted P2 preparation — 2026-11-08

**No paired trajectories ran. No thermal figure or source/control numerical result
is qualified.** The owner approval is recorded and remains effective, but the
conditional sidecar stops on a newly identified C-REF interpretation problem.
The single next scientific action is the adjudication in proposed ADR-0022.

## Source and storage

Approved delivery: `4509e13217ec442a2649ab0e2dcf950f3c0f5c97`.
Scientific candidate ancestor: `bb0f9c5af8ebaf86b5a59aece5551bb9dafbd264`.
Approval registry commit: `961ced1`. C-REF model/assay commit: `ac22c6c`.
Native verification source: `8492d6383fc5990add43bb9addd544a4abd378da`.
Review branch: `candidate/results-first-bnv-pilot-20261108`.
Repository: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
External outputs/builds:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-preparation-20261108`.

The initial branch/history audit found the approved delivery on this clean
worktree, descending from the reviewed Linux candidate and all 14 existing branch
tips. Master was not used as the source, merged or changed. Available storage was
248 GiB at preflight. All substantial data, builds, copied mutation fixtures and
temporary files remain under /mnt/sdd. No existing worktree/archive was moved or
deleted. This is a Linux session; the Mac's Scriptorium mount/writability is not
observable. No Mac copy was made. Direct transfer instructions are retained.

## Scientific gate result

| Gate | Result |
|---|---|
| Restricted preparation approval, S and seven conditions | Recorded separately; trajectory launch remains unauthorized |
| Metric consumer limitation | Accepted for restricted P2; no nodewise PASS and no historical certificate promotion |
| Conditional sidecar | BLOCKED: conservative pointwise C-REF completion fails 42 rows |
| Historical certificate/currentness behavior | All 462 historical rows preserved; runtime rows, monitor and refusal code unchanged |
| Experiment/accuracy definition | Existing 14,648.4375-year predeclaration unchanged; no 500 kyr card edits |
| Native ordinary-context identity | PASS: exact owning star, all four native hashes, 15 refusal cases; does not supply BA13 qualification |
| Source/control snapshots and paired runner | Stopped at the C-REF gate; not claimed ready |
| Actual R18/R20, nonlinear accuracy, signal resolution | Not evaluated on this pilot; require the future authorized trajectories |

C-REF bounds a shared residual by ten times observed production–X2 agreement but
does not specify its dependence on B. Taking the retained absolute pointwise
agreement gives the following **partial assigned error budgets**, before the
remaining positive terms. The unchanged utilization ceiling is 1.

| Quantity | max absolute P-X2 | C-REF allowance | Maximum partial utilization |
|---|---:|---:|---:|
| Envelope temperature | 2.549314522370949e-6 | 2.549314522370949e-5 | 1.8413247593801518 |
| Redshifted surface temperature | 2.550753074315182e-6 | 2.550753074315182e-5 | 1.8404747993416513 |

Both fail at every retained row: 42 failures. The review's <=1e-3 utilization-change
expectation fails at 355 rows. All other rows remain pending full completion;
passing a lower-budget check is not qualification. This does **not** measure the
star's true error or prove a sharper correlated estimate cannot pass.

Replacing absolute agreement by agreement of frozen-relative increments, or
assuming a B-independent common offset, would add a correlation statement not
explicit in the accepted text. It was not done. `CREF_COMPLETION_MODEL.md` freezes
the formula; `preparation-evidence/cref-summary.json` records the result. The full
diagnostic sidecar retains every historical field and status and lives externally,
with a committed SHA-256 manifest. It cannot be loaded as an operational
certificate. No missing term was silently set to zero or filled from headroom.

## Native identity and Bdot

The generic context now defaults to the same four historical hashes and offers a
separate, fixed native tuple for analytic controls. Caller-supplied replacement
hashes are not accepted. Exact owning-profile equality, lifetime source coverage,
token revocation and full byte-currentness remain mandatory. This is an engineering
identity adapter, not an accuracy certificate or a production-launch permission.

The binary64/80-digit witness reconciles the old absolute rate without a tolerance:

| Input | B0 (count) | Derived Bdot (count/s) |
|---|---:|---:|
| Historical Mac | 7.6169065187188905e56 | -2.4136520263641375e37 |
| Native Linux | 7.616906518718875e56 | -2.413652026364133e37 |

Both literals are reproduced exactly by `(-1e-12)*B0/31557600`. Native B0 differs
by -18 binary64 ulps (-2.058606263472188e-15 relative); Bdot differs by 10 ulps
(-1.956523322868249e-15 relative). The fractional drive and physical definition
are unchanged. The predicted depletion remains 1.46484375e-8. Historical cards
and assertions remain intact; the native value is separately checked against the
actual owning star and its coefficient table. This is no microscopic-rate input.

## Verification and resource record

Five new authority/algebra/refusal tests and all 12 retained final-decision tests
pass. The C-REF Slurm assay used one CPU on phyn004, a 60 CPU-second limit, and
1.155685 CPU seconds. It retained the failure and issued no certificate.

Native verification results and final Slurm receipt are recorded in the adjoining
preparation evidence. Job 19147 passes: four CPUs requested on phyn004, 8 GiB
memory requested, one-hour wall limit, total Slurm runtime 2m10s and exit 0;
the actual fixture used 121.286 seconds
wall, 120.122 CPU seconds and 316380 KiB maximum child RSS. Exact native B0/Bdot,
15 refusal cases and the separate byte-currentness check pass. All four retained
synthetic tests pass: direct energy, thermal ledger, passive checkpoints and
breakpoint contracts (CTest 0.06 seconds). These are not actual-pilot R18/R20
results. The successful native executable invoked zero actual BNV RHS calls.

Earlier attempts are preserved: job 19141 stalled in
`nfs_wait_bit_killable` on phyn001 and was cancelled after 6m44s; 19143 refused a
missing explicit GSL path; 19144 refused missing nested Debug package prefixes;
19145 found GNU time absent on the compute node; 19146 built successfully but
the retained fixture refused provider-source authentication from the wrong working
directory. These were launcher corrections, not scientific guard relaxations.
The final launcher uses the existing dependency lock, Release flags including
`-ffp-contract=off`, portable resource recorder and required repository directory.
Build 19146 took 74.254 seconds, 217.150 CPU seconds, maximum child RSS 432132 KiB.

Automatic approval review initially rejected submission because it read the
mutation test as modifying retained data. The test already copied its inputs;
explicit canonical-path and inode checks were added, and resubmission was approved.
All mutations occurred only in a fresh private fixture. No approval-review action
remains blocked. Slurm accounting storage is disabled; retained controller and
process-resource receipts are the available evidence.

## Deliverable limits and next action

There is no temperature difference, onset time, direct/reaction heating history,
integrated energy residual, actual-path accuracy measurement, or utilized frozen
certificate to report for a pilot that did not run. Synthetic contract checks do
not substitute for those results. The source/control comparison, figure/data and
plotting implementation remain gated; no result is presented as a physical BNV
prediction. Full-campaign authorization remains absent.

**Smallest next action:** adjudicate C-REF's common residual as a precise
B-dependent error model, including its increment, fit and normalization response,
against the retained arrays and original review. See proposed ADR-0022. Merely
choosing a smaller allowance or increasing a threshold is not a resolution. After
that one issue is settled, complete the other predeclared sidecar terms and the
bounded paired preparation, then request separate trajectory authorization.

The physical follow-on recommendation is unchanged: n -> chi + gamma, conditional
on chi's quantum numbers and validated escape/photon thermalization. It needs an
in-medium differential rate and consistent particle fates, with a source/partition
adapter; its first possible observable is surface-temperature/luminosity excess
versus age. No rate is selected. See PHYSICAL_BENCHMARK.md for the primary source
and energy bookkeeping; no new physical implementation was attempted.

Change classes: accepted scoped qualification policy, candidate validation
numerical-method, engineering identity/build, documentation and compact evidence.
Governing authority: GOVERNANCE §§1–5, ADR-0015, ADR-0017 and recorded owner
approval. No thermal, chemical, stellar or ledger equation changed. Scientific
impact is limited to making the already identified native input tuple explicitly
selectable and refusing an unsupported qualification interpretation. Candidate
status remains unpromoted. Files changed are enumerated by the branch diff and
the preparation manifest; historical package manifests remain valid.
