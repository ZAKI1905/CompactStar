# ADR-0022 — C-REF common-residual interpretation for the restricted P2 pilot

Status: ACCEPTED FOR RESTRICTED STATIC PREPARATION; implementation qualification
pending. Original proposal date: 2026-11-08. Acceptance recorded: 2026-11-09 UTC.
Change class: scientific qualification semantics / validation numerical-method.
Authority: GOVERNANCE §§1–5; AGENTS §8; ADR-0015/0017; accepted restricted
RF-BNV-P2-20261108 preparation decision. ADR-0021 remains proposed outside that
scope. Blocks: restricted conditional sidecar and actual paired trajectories.

## Context

The owner accepted the review's Appendix A conditions, not arbitrary replacement
error bounds. C-REF says any shared residual is at most ten times observed
production–X2 agreement, with a separate Ltilde_Me exception. It does not specify
whether that agreement is pointwise absolute error or error in the frozen-relative
increment, nor a covariance or B-dependence for the shared residual.

The retained comparison distinguishes these objects. For envelope, max|P-X2| is
2.549314522370949e-6, whereas max|(P-X2)(x)-(P-X2)(0)| is
6.480784642270976e-8, in the historical normalized units. They cannot be substituted
for one another without a scientific statement. Error shared between two solvers
does not, by itself, mean error constant along the depletion sequence.

The precommitted static audit at `ac22c6c`, Slurm 19142, applies the conservative
pointwise interpretation to the existing joint inequality. Even before all other
missing positive contributions are added, envelope and Tsurface_inf exceed the
unchanged threshold at all 21 rows (42 failures). Their maximum assigned-budget
utilizations are 1.8413247593801518 and 1.8404747993416513. The review's <=1e-3
utilization-change expectation also fails (355 rows).

This is a failure of this error-budget completion, **not a measurement of physical
error** and not a proof that no sharper correlated completion can pass. All 462
historical rows are unchanged. The diagnostic sidecar is not an operational
certificate. See `../validation/results-first-bnv-pilot/CREF_COMPLETION_MODEL.md`
and its retained summary for the exact algebra and limitations.

## Decision

The owner selected alternative 2, specifically as amended by the independent
C-REF adjudication, starting at `30f2ab39572abdf33f51f0ef9609d4fef9b90642`.
Authority and the exact instruction are recorded in
`../validation/results-first-bnv-pilot/static-completion/OWNER_APPROVAL.md`.
The operative equations and limits are in the adjoining `MODEL.md`.

Retain the full ten-times pointwise allowance, the worst-case 2e increment,
baseline normalization, historical thresholds and unchanged fit. For envelope
and Tsurface_inf only, use the accepted C-CONT curvature and the unchanged
all-row fit to bound the TOTAL reference fit response. No B-dependent shared
error assumption is added. C-TAN and C-INPUT remain separate. The 42 original
failures are preserved; amended partial passage is not full qualification.
The earlier utilization-change expectation still fails on 353 rows at the
amended partial stage and must be reported. No trajectory, historical certificate
promotion, master merge or new stellar calculation is authorized by this ADR.

## Alternatives

1. Retain the pointwise absolute C-REF model and the published conservative joint
   propagation. The present completion fails; no pilot launch. No code or
   threshold changes. Existing historical reproduction remains unchanged.
2. Clarify C-REF with an explicit, scientifically justified **B-dependent common
   residual model**, bounding both the residual's frozen-relative increment and
   its response through the unchanged all-row fit and B0 normalization. This
   could distinguish a constant common offset from its varying remainder.
   The statement must specify the quantitative estimates and their provenance;
   replacing absolute discrepancy by increment discrepancy alone is insufficient
   justification. No bound is inferred from available headroom. Update only the
   separate owner registry and restricted assessor after acceptance, then repeat
   the bounded static completion. No historical certificate promotion follows.

## Validation and consequences

The smallest next action is a precise adjudication of that one C-REF statement,
using the already retained comparison arrays and the original review. This is not
authorization for another star, numerical solver or coefficient campaign. A
clarification must explain the failed pointwise check; it cannot silently change
the historical thresholds or mark unknown contributions zero. Surface-cut and
the remaining ledger terms still must be completed after the necessary C-REF
gate passes. Actual R18/R20, nonlinear accuracy and source/control signal
resolution remain pilot-time requirements, not results of this static audit.

## Provenance

The implementation remains an AI-authored unverified scientific candidate until
the authorized static checks establish its stated, conditional scope. The original owner's restricted preparation
approval remains valid; it has not been reclassified as pending. Only this newly
identified interpretation is unresolved. Separate trajectory authorization is
still required after all preparation gates pass.
