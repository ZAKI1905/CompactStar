# ADR-0022 static completion — scientific disposition

**BLOCKED_UNSUPPORTED_CUT_RESPONSE.** ADR-0022 is implemented and its amended
necessary C-REF check succeeds for all 42 previously failing temperature rows.
The full static qualification does not succeed: the required fixed-B response
to the accepted surface-pressure perturbation remains unsupported. No trajectory
authorization request is issued. Evidence: `evidence/summary.json`,
`evidence/ledger.json.gz`, `evidence/independent-ledger-audit.json`.

## Source and authority

Repository:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
Branch: `candidate/results-first-bnv-pilot-20261108`. Clean starting source:
`30f2ab39572abdf33f51f0ef9609d4fef9b90642`; all other local branch tips were its
ancestors. Other worktrees were inspected but not modified. This is the intended
source; master was not merged, used as the working source or changed.

Authority/model commit: `cd98e1c` (recorded before the assessor consumed it).
First implementation: `9d7f5a5e3cce3cce761b4e64b26d962d10bbf3ec`.
Final numerical execution: `6a4e109951df63e9d2a1002ad026fa99b6758331`.
The delivery commit adds only documentation and retained evidence after that
execution; exact delivery identity is returned with the handoff. The source
definitions are authenticated by `acceptance-registry.json` and the 30 retained
inputs by `inputs.json`.

Change classes: scoped scientific qualification semantics, validation
numerical-method, engineering, documentation and generated evidence. Authority:
GOVERNANCE §§1–5, ADR-0015/0017, accepted preflight §11, the existing restricted
approval and `OWNER_APPROVAL.md`. Implementation remains a scientific candidate;
acceptance of an assumption is not validation of absolute stellar accuracy.

## Scientific result

| Matrix disposition | Rows |
|---|---:|
| Exact P2 identity | 21 |
| Conditionally qualified species support | 21 |
| Metric consumer scope only; no nodewise qualification | 21 |
| Unqualified approximate quantities | 399 |

All 462 original records, including every historical field, status and term, are
preserved inside `evidence/ledger.json.gz`. The original pointwise assay's 42
failures and 355 expectation violations reproduce exactly. All 42 temperature
rows satisfy the amended necessary test, but remain unqualified in the full
matrix. `evidence/matrix.tsv.gz` is the compact row-by-row disposition.

The prior preparation manifest verifies against every original blob at the
starting commit. Its sole current-file difference is the expressly authorized
acceptance of ADR-0022; the old manifest is not rewritten to conceal that change.
The original proposed ADR is retained as `evidence/ADR-0022-proposed-at-start.md`.
See `evidence/historical-preservation.json` for the complete hash comparison.

| Group | Amended necessary partial maximum | Supported completion partial maximum | Full utilization |
|---|---:|---:|---|
| Envelope temperature | 0.5202795986735851 | 0.5202796022226576 | unknown |
| Redshifted surface temperature | 0.5186840956442711 | 0.5186840987286260 | unknown |
| N_mu (largest supported completion partial) | 0.7824646060132245 | 0.7824646077302390 | unknown |

No supported partial allowance exceeds the unchanged historical threshold.
This is not evidence that the missing allowance fits. No headroom supplies an
error estimate. `evidence/summary.json` reports every group's partial maximum.

**353 rows violate the earlier utilization-change expectation of <=1e-3.** Every
row identity, increase and amended partial utilization is listed in
`evidence/amended-expectation-violations.json`. The supported completion stage
also has 353 such violations, separately listed in
`evidence/completion-expectation-violations.json`. These are reported diagnostic
failures, not retuned acceptance limits.

## Uncertainty ledger

The machine-readable ledger contains a separate restricted entry for every row,
with every contribution, baseline allowance, evidence class, condition, source
identifier, partial total and nullable full total. Original ledger entries remain
in the same row unchanged. Its numerical scope is an inventory with explicit
unknowns, not a completed uncertainty enclosure.

| Required contribution | Completion/disposition |
|---|---|
| Production minus achieved-B X2, reference interpolation, original finite ledger | Retained exactly; original numerical and conditional evidence preserved |
| C-REF common amplitude/increment | Full ten-times allowance, full 2e increment, and baseline e retained; Ltilde_Me exception unchanged |
| Two temperature fit responses | Joint C-REF/C-CONT bound on total reference response; no shared-error B model |
| Other fit responses | Original all-row f*e allowance |
| Stage and arithmetic | Single finest-increment allowance, split only for bookkeeping and propagated through drift, fit and baseline |
| EOS/centre parameter remainder | No duplicate beyond the full estimates accepted under S/C-REF |
| C-TAN | Ratio-level GL4; second stencil tail; moving-cut/path residual; denominator separation and division correction; no independent endpoint cancellation claim |
| C-INPUT | Existing coordinate/floor responses; accepted second-order input term; fixed weights; full normalization kept separate |
| C-CONT/support | Original gap terms; endpoint known-error slope/sliver contribution; accepted topology condition |
| C-MUN | Equilibrium-reference identity only |
| C-CUT finite perturbation | **Unsupported: observable-specific response at fixed B is absent** |
| Baseline cut contribution | **Unsupported for relative rows** |
| Endpoint cut contribution | **Unsupported for endpoint rows** |
| Full normalization and total | **Unknown wherever their cut input is unknown** |

The known fit response does not erase the missing cut contribution to drift or
normalization. The accepted C-TAN path correction is not silently reused as a
bound on an independent cut-input perturbation. The accepted C-INPUT ranges
remain explicit conditions; this work does not invent larger or smaller ranges.

## Exact remaining scientific issue

C-CUT supplies `|delta ln Ps| <= 2.1e-12`, but not a bound on

    |Q(B, ln Ps + delta) - Q(B, ln Ps)| / S

at fixed B over that perturbation range. `continuation.json` supplies pressure-law
residual and geometry motion, `m1.json` supplies GL4/EOS/centre/stencil estimates,
`m3.json` supplies B-direction derivatives, and `diagnostics.json` supplies tangent
cross-geometry residuals. Their schemas and generating source were inspected.
The accepted one-dimensional K on the physical path does not bound this separate
cut response. A finite pressure residual is not an observable error allowance.

The single next scientific target is an authenticated fixed-B finite-cut response
bound, including the tangent and baseline consequences. A supported analytic
derivation from the retained terminal-state equations would be sufficient; a
local derivative without a finite-perturbation remainder would not. No such
bound is supplied by the retained package, and none is manufactured here. New
stellar calculations have **not** been demonstrated indispensable and were not
performed or authorized by this implementation. This is the owner's explicit
unsupported-allowance stop, not a request to reopen general solver development.

## Verification and resources

Final job **19149**, phyn004, executed at
`6a4e109951df63e9d2a1002ad026fa99b6758331`. All 26 tests passed: five original
C-REF tests, twelve final-decision tests and nine new authority/algebra/refusal
tests. The independent serialized-ledger audit passed **5,697 checks** and
refused all five in-memory mutations: dropped cut error, halved increment,
dropped denominator error, partial promoted to full, and changed historical
threshold. No retained file was mutated.

The 26 temperature-component fits were reconstructed independently at 70 and 90
digits; affine identities, all-row norms, deterministic discrepancy fits and
curvature bounds agree. Nine raw tangent ratios were reconstructed at both
precisions, including the correlated numerator/denominator identity and retained
order/width gates. Reported binary64 values agree exactly across precisions.
See `evidence/precision-witness.json` and `evidence/checks.json`.

Final aggregate CPU time: **10.640776448 s**, below the 60 CPU-second ceiling;
wall time **11.006109438836575 s**; maximum child RSS **283652 KiB**. The job
requested one CPU/task and constrained process affinity to one CPU. Slurm reports
two allocated CPU entries for the single-core allocation; its controller record
is retained without reinterpretation. Exit **2** is the deliberate scientific
refusal after writing evidence; tests and independent audit exit 0. The scheduler
therefore records FAILED/NonZeroExitCode. This is not a successful qualification.

Earlier job 19148 reproduced the same disposition and passed the 26 tests before
the independent serialized-ledger audit was added; CPU **5.535156524 s**. Its
receipt and summary are preserved in `evidence/prior-19148/`, and its complete
external outputs remain intact. Final outputs are under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/results-first-bnv-static-20261109-verified`.
All work and substantial files remain on EKU /mnt/sdd; no Mac-side work occurred.

Files changed: ADR-0022 acceptance; this new static-completion documentation,
registry, input manifest and evidence; `static_completion.py`,
`static_completion_run.py`, `static_validate.py`, `static_completion.sbatch`, and
`test_results_first_static.py`. Historical fitter, assessor, certificate, runtime
guards, original approval registry and prior evidence are unchanged. Scientific
impact: new conditional validation accounting and an explicit refusal, with no
change to any physical solver output. No new stellar calculation, BNV RHS call,
trajectory, certificate promotion or master merge occurred.
