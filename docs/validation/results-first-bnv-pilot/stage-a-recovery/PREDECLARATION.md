# Authorized Stage A recovery — one bounded retry

Owner authorization starts from `4387e56c6b218411d755d27148f654aea5bc345a`
on `candidate/results-first-bnv-pilot-20261108`. The owner accepts the conditional
knot-aligned D1 receipt and authorizes precisely O2 chemical absolute tolerances
`(1e-24,1e-24) -> (5e-27,5e-27) MeV`, one source-first one-year Stage A retry,
and matched zero-source control only after the source passes every required gate.
No other solver setting, physical input or acceptance criterion changes.

Authority: current explicit owner instruction; GOVERNANCE §§2–3/5; accepted
ADR-0015/0016/0017 (especially ADR-0017 §§5–7); INV-11f and INV-15; the original
[RUN_PLAN](../native-regeneration/RUN_PLAN.md#L1) and
[NONLINEAR_PLAN](../native-assembly/NONLINEAR_PLAN.md#L1); review Section D and
the [conditional D1 resolution](../echem-refinement/d1-resolution/REPORT.md#L1).
Change classes: numerical-method (approved tolerances), engineering (execution
boundary and evidence), dependency/build, documentation/generated evidence.
This remains a scientific candidate pending owner assessment, not historical
certificate promotion. No pre-baseline exception or new physical equation.

## Preserved plan and stop rules

Both roles start independently from t=0, x=0, eta=(0,0); duration is exactly
31,557,600 seconds. Each role retains the original five references, twelve
candidate tier/cap runs, twelve endpoint-schedule repeats and three quarter-cap
restart tails. Preserve all four scientific grids (240/480/960/1920 intervals),
the original main/reference configurations, caps and original R20 trapezoids.
The source collector stops before control: collection completion alone never
authorizes the zero-source trajectory. The new wrapper cannot select Stage B.

Retain state/power norms, 0.05 ULTRA reference allocation, U<=0.20F, R18/gross
power<=1e-10, R20 total<=2e-4N, reconstruction/event<=5e-6N, convergence<=0.5
above independently bounded roundoff, restart/prefix accounting, exact-zero
control chemistry and matched signal U_delta<=0.10|delta_ref|. A 64-ulp floor
is not a bound on accumulated roundoff. Missing uncertainty is UNQUALIFIED, not
zero. Source, currentness, conditional frozen, temperature, depletion, finite
progress and rejection guards remain unchanged. At the first failed mandatory
criterion preserve the failed channel, state, allowance and partial output,
then stop. Timeout means INCOMPLETE, with no automatic resubmission.

The source-only selection and CheckpointEvidence serialization change execution
or reporting, not solver mathematics. Every actual checkpoint packet, including
a refusal, is written before the collector can throw. Previous collection
silently omitted the rejected packet; the new evidence makes any stop reviewable.
No checkpoint selection is based on the 13-point diagnostic study.

## Actual-output qualification

The original D1 FAIL and the later conditional PASS are both immutable evidence.
The 13 tested intervals do not qualify every new output. If native source
collection survives its mandatory gates, apply the knot-aligned 110-digit
chemical diagnostic to its actual strict-interior checkpoint intervals, with
the same fixed segment subdivisions, thermal offsets, RHS fidelity and C1–C3
criteria. Exact endpoints require main-prefix accounting rather than an
invented local reconstruction allowance. Deduplication is permissible only on
exact identity of the complete interval, initial state, thermal nodes and O2
output/allowances; any reused result must identify that identity.

Preserve the uniform continuous thermal-enclosure assumption and the assumption
of continued reference-refinement contraction. Sampled thermal agreement is
not a proof of uniform enclosure, nor are constant offsets rigorous extrema
over arbitrary time-dependent perturbations. Main-prefix error from t=0,
thermal RHS/accumulated roundoff, global reference and interpolation error,
event uncertainty, R18/R20 and restart consistency require their own actual-output
assessment. Report missing or unexecuted terms explicitly. No small raw spread
is substituted for a supported bound. A prior native mandatory failure stops
the retry before further local-reference or remaining matrix calculations.

Only if source qualification passes will a separately gated exactly matched
control be executed within the authorized Stage A resource ceiling, followed
by the paired state/power/energy/ledger and signal-resolution assessment. Any
unqualified source term prevents control. No automatic chained submission.

## Authentication, build and resource bounds

The worktree is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
It was clean at 4387e56, with all 15 local branch tips ancestors and nine
worktrees. Original snapshots and D1 manifests authenticated before edits;
see [preflight.json](preflight.json#L1).

Build into a fresh external tree; never overwrite original archives or binaries.
No stellar regeneration target is executed. Run the six retained contracts and
a compiled O1/O2 configuration contract, then rebuild owners from the unchanged
full archive and compare all twelve static RHS states byte-for-byte. Publish
new source/control identities only after these pass. The publisher checks that
the production source differs from 4387e56 by exactly the two approved literals,
and rejects any unlisted change to previously bound dependencies. Original
archive and physical configuration identities remain unchanged; record new ELF
and runtime bindings, source-file differences and the explicit O1/O2 settings.

Use Slurm phyall/phyn004, one logical execution CPU (verify affinity), 4 GiB,
30-minute ceiling, single-thread libraries. The build/static allocation is
separate from thermal execution. Stage A retains the original aggregate thermal
30-minute ceiling; no retries or extension. All large data/build/tmp outputs
stay under /mnt/sdd. Git tracks configuration, source, compact results and hashes.
No Mac-side work is needed.

## Delivery and exclusions

If both roles qualify, publish the first paired redshifted-temperature figure,
temperature difference, chemical imbalances, deposition, ledger consistency and
numerical uncertainties. State explicitly that this short mathematical P2
experiment primarily tests direct BNV deposition and does not establish a
measurable chemical-heating contribution. If a gate fails or the allocation
expires, report its exact obstacle and quantitative significance, preserve
partial results and stop. Never present unqualified source-only output as the
requested paired result.

No Stage B (2048 years), 14,648-year benchmark, 500 kyr campaign, stellar
regeneration, BA13 reopening, master merge, or historical certificate promotion.
