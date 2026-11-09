# ADR-0022 restricted static completion — owner approval

Recorded 2026-11-09 UTC. Source of authority: the project owner's explicit
conversation instruction below. Starting branch:
`candidate/results-first-bnv-pilot-20261108`; starting commit:
`30f2ab39572abdf33f51f0ef9609d4fef9b90642`.
This extends the existing preparation approval; it does not replace the original
registry, certificate, failed assay, or preparation evidence.

## Exact owner instruction

> CompactStar — ADR-0022 Acceptance and Static Completion
>
> I approve ADR-0022 alternative 2, as amended by the independent C-REF adjudication, for the restricted RF-BNV-P2-20261108 preparation.
>
> Starting commit: "30f2ab39572abdf33f51f0ef9609d4fef9b90642"
>
> This approval is limited to the joint C-REF/C-CONT uncertainty propagation described in the adjudication.
>
> Requirements:
>
> 1. Retain the full ten-times pointwise C-REF allowance, the worst-case (2e) increment, baseline normalization, original historical thresholds, and unchanged fit.
> 2. Use the curvature-derived bound on the total reference fit response for the two surface-temperature groups. Do not introduce an unsupported B-dependent shared-error model.
> 3. Preserve C-TAN and C-INPUT separately, including their correlated tangent, coordinate, and nonlinear remainders.
> 4. Preserve the original 42 historical failures. Record that they now satisfy the amended partial C-REF test; do not automatically promote them to full PASS.
> 5. Explicitly report the 353 rows that violate the earlier utilization-change expectation.
>
> Proceed with the previously authorized full static sidecar completion using retained data.
>
> Implement the amendment, verify the fit identities and high-precision calculations, complete all remaining uncertainty terms, and rerun the full qualification matrix and relevant tests.
>
> Use the EKU cluster for numerical work, and Scriptorium for any new Mac-side worktrees or substantial files. Keep source and reproducibility records Git-tracked.
>
> Do not introduce new stellar calculations unless a specific, documented failure makes them indispensable. Do not reopen a general solver investigation.
>
> Stop and report if any required numerical allowance fails, remains unsupported, or depends on a new scientific assumption requiring owner approval.
>
> Deliver the updated scientific disposition, complete uncertainty ledger, exact commit, and evidence.
>
> No historical certificate promotion, master merge, or BNV trajectories.
>
> If the static qualification succeeds, prepare the exact authorization request for the first paired BNV-source versus no-BNV-control evolution.
>
> The objective is to finish this qualification step and proceed directly toward our first temperature-evolution figure.

## Scope

Change classes: scientific qualification semantics / validation numerical-method,
engineering, documentation and generated evidence. Authority: GOVERNANCE §§1–5,
ADR-0015/0017, accepted preflight §11, original OWNER_APPROVAL.md, this instruction
and ADR-0022. Production equations and runtime refusal behavior are unchanged.
Missing evidence must remain explicit, even if the partial utilization is below 1.
