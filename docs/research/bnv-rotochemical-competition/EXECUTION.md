# Execution and verification record — 2026-11-12

Entry was the clean candidate branch at
81b808f5e85251c63622ed1a8bade2c206fdbf82. The original working directory was a
workspace container, so work used its existing results-first-bnv-pilot Git
worktree. All nine worktrees were enumerated and preserved. All local branch
tips were ancestors of the entry; nonancestral historical remote refs retain
their documented noncanonical status. The [entry receipt](evidence/entry-authentication.json)
records that distinction and exact Stage B review authentication.

The independent Stage B review's B — CONDITIONAL PASS was recorded as the
owner's restricted acceptance only after checking the reviewed
49253200a801ad9b336789dffb4e93d535a8e903 ancestry, imported review bytes, selected
evidence hashes and corrected report. C1–C7 remain unchanged. No Stage B solve,
certificate promotion, master merge or history rewrite was performed.

## EKU jobs and retained corrections

All scientific work below used Slurm on phyn004. The static map parallelizes
independent drive rows; native short runs share one expensive owner reconstruction.
Allocation receipts and output hashes remain in the external manifest.
Slurm accounting is disabled, so scontrol receipts and completion/log artifacts
are used instead of sacct.

| Job | Work / result |
|---|---|
| 19419 | r1 static map completed, 12 CPUs, about 29 s. Exact cancellation formed by subtracting large binary64 drives left spurious tiny channel imbalances. Retained as superseded output. |
| 19420 | r2 factored exact-null map, 12 CPUs, about 29 s. A stability assertion failed because ordinary eigenanalysis rounded a tiny negative mode to zero. No physical instability was inferred. |
| 19421 | Dependent verification cancelled after 19420 failed; no calculation result. |
| 19422 | r3 final map, 12 CPUs, 30.36 s. Near-degenerate eigenproblems use 55-digit eigenanalysis; all map residual/stability checks pass. |
| 19423 | Independent root checks computed; JSON serialization of a NumPy boolean failed. The log is retained. Explicit Boolean conversion fixed the output layer. |
| 19424 | Independent static verification and initial figure rendering pass. |
| 19425 | Build launcher failed before compilation: cmake missing from batch PATH. Its allocation receipt remains; the common configure.log was replaced by the successful configure. No failed compilation or scientific result is claimed. |
| 19426 | Explicit toolchain paths; fresh production build and 9/9 targeted contracts pass, 16 CPUs, 76 s. |
| 19427 | Short launcher failed before native invocation because compute-node PATH lacks git. The traceback is retained; source commit is now passed explicitly. |
| 19428 | Supplement calculated; JSON serialization of NumPy int64 failed. Retained log; converted count to Python int. |
| 19429 | Four native one-year matched cases, two tolerance tiers and passive schedule repeats pass; one CPU, 114.29 s including full owner reconstruction. |
| 19430 | Supplementary cancellation/timing calculations and thermal-root scan pass. |
| 19431 | Independent DOP853 short integrations and separately scaled thermal/chemical energy integrals pass. |
| 19432 | Final six PDF/PNG layouts rendered from existing data; no scientific map rerun. |
| 19433 | Added pointwise checks of each reaction/work power on its own scale; all pass with maximum relative disagreement 1.67e-15. No production or stationary-map rerun. |

These are bounded corrections to the new functionality. Prior BA12/BA12R and
Stage A/B failures/evidence remain unchanged. No long integration was launched.

## Proportionate checks

The [static verification](evidence/verification.json) tests complete polynomials,
three-variable independent roots, independent derivatives and full-system
stability, tightened tolerances, a second solver, both nulls, symmetry, zero drive,
dimensions and conventional RC limits. The 55-digit calculations cover 21
representative states; they are not a repeat of Stage B's certification campaign.
The [supplement](evidence/supplement.json) adds 183 near-null states, 32 period
examples and 15 thermal-manifold bracket scans.

The [CTest log](evidence/logs/contracts.log) records:
passive_checkpoint_output_contract, phase6a1_direct_energy,
phase6a1_thermal_ledger, phase5d_response, phase5d_independent_oracles,
phase5d_component_tolerances, phase6_breakpoint_contract,
phase6_currentness_bytes, and bnv_competition_contract — all pass.

The native S2 runner tests common spin ownership, exactly-once combined driving,
legacy P0/P2 partition reductions, zero-source reduction, correct signs,
chain-rule chemical work and exclusion of duplicate thermal work. It exercises
refusal of fast/unqualified spin, reaction-free combined mode, out-of-range
temperature and excessive depletion. Production tolerances are the existing
1e-9 and 1e-11 tiers with their component scales. Both give identical printed
endpoints and 12 accepted/zero rejected steps per case. Observation schedules
leave accepted states unchanged.

Maximum observed short depletion is 1.000026e-12 and maximum frozen utilization
is 3.47485e-6. Temperatures remain near 1e8 K and imbalances remain inside the
reviewed Stage B state box, while the **new driving functionality has its own
proportionate checks**, not a borrowed Stage B global-error certificate.
The [independent short check](evidence/short/independent-check.json) uses DOP853
at 2e-10 and 2e-12 relative tolerances with separately scaled energy integrals.

All six final PNG exports were visually inspected for axes, legends, labels,
limits and footers; vector PDFs carry the same Matplotlib figures. Reachability
color ranges deliberately saturate extreme channel-null modes, with explicit
annotations and explanation in the report. No physically valid regions are
inferred from a plot color alone.

## Artifact boundaries

The committed bundle contains source, scripts, original compact inputs, six
PDF/PNG figures, compact all-grid tables, full fiducial/representative tables,
short trajectories and selected verification evidence.
Full stationary outputs, superseded maps and native build files remain under
/mnt/sdd. [external-receipt.json](evidence/external-receipt.json) names the external
SHA-256 manifest. That manifest includes all retained scientific outputs/logs
in r1/r2/r3 and the new native binaries; build intermediates, dependency packages
and plotting caches are excluded explicitly. The retained complete stellar
archive has its own previously authenticated hash.

[final-source-sha256.json](evidence/final-source-sha256.json) binds the delivered
implementation and scripts. The static summary records the source bytes at
map time; later plotting-layout and output-serialization fixes have separate
final source/figure receipts. Native invocation records every production source
hash and the built executable hash before running it. The packaging checks
these bindings before commit. Committing and pushing candidate work does not
ratify the new scientific results.
