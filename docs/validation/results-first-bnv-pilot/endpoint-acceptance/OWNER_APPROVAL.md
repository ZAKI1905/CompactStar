# Restricted endpoint continuation — owner decision

Status: ACCEPTED CONDITIONAL SCIENTIFIC ASSUMPTION FOR STATIC PREPARATION ONLY.
Owner: project owner by explicit conversation instruction, session 2026-11-08.
Starting source: `afd4b06074a5f597ef4cbe38ce3c9cb04386f556`.
Authority: GOVERNANCE §§1–5; ADR-0022 accepted endpoint amendment;
`../OWNER_APPROVAL.md:1` and `../static-completion/OWNER_APPROVAL.md:1`.

Exact scientific decision:

> I approve, for the restricted RF-BNV-P2-20261108 experiment only, the conditional extension of the previously accepted equilibrium-family bands, curvature bounds, positive baryon derivative and fixed topology over the required endpoint extensions of at most (3.022\times10^{-26}) in the depletion coordinate.
>
> This is an explicit conditional scientific assumption for the C-CUT fixed-B surface-transfer assessment, not a newly established rigorous numerical bound.
>
> Do not extend the assumption to other stellar models, physical regimes or larger depletion intervals. Verify that no known endpoint transition contradicts it.

The owner also directs recording the decision, completing the C-CUT ledger from
the existing witness, recomputing the full restricted matrix, running relevant
numerical/audit/mutation tests, checking nonlinear and trajectory readiness, and
preparing an exact paired-evolution authorization request only if qualification
succeeds. The owner prohibits additional stellar calculations unless a new,
specific scientific failure demonstrates necessity, and prohibits master merge,
historical certificate promotion and BNV trajectories without separate authority.
EKU is the numerical host; substantial Mac artifacts would require Scriptorium.

Operational scope: keep the original physical interval [-1e-6,0] and use only
the retained eta from `../surface-transfer/evidence/assessment/bounds.json:1`
for remapping. The limit 3.022e-26 is a cap, not a replacement uncertainty.
Retain S, C-REF, C-TAN, C-CONT, C-INPUT, C-SUPPORT, C-CUT, C-MUN, the ADR-0022
thermal-fit amendment and the metric-consumer limitation. No shared-error
correlation, new stellar regime or larger depletion interval is accepted.

The historical 42 pointwise C-REF failures, 355 original/353 amended expectation
violations, six inconclusive endpoint grid checks and BA12/BA12R failures stay
historical failures/inconclusive results. Their status is not rewritten by this
decision. The new condition is named C-ENDPOINT in the successor sidecar solely
to keep it visible; that name creates no additional scope.
