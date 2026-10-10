# ADR-0017 Amendment 1: breakpoint-aware Phase-6 evolution

Status: **ACCEPTED / HUMAN-RATIFIED, 2026-10-10**, limited to the bounded
architecture at `0c3f9269c884161c9d89163660a47b141a93ace6`; see
[ratification](PHASE6_BREAKPOINT_RATIFICATION.md). The dated proposal text below
is retained as the exact accepted contract, with its scope limits unchanged.

Historical proposal status (2026-10-03): owner instruction authorizes
candidate implementation and bounded qualification before ratification. The full
campaign remains stopped. Historical BA12, BA12R and obs-32 FAIL records stand.
This replaces the prior reconstruction-only proposal; accepted ADR-0017 is unchanged.
Classes: numerical-method, architecture, engineering, documentation.

## Candidate contract

Amend ADR-0017 sections 4, 5.2-5.4, 7, 10 and 11 prospectively. Both Phase-6
main RKF45 and local rk8pd use the actual Cstar cache boundaries. An immutable
snapshot exposes exact temperatures, capacities and cache cell-selection semantics.
No table, interpolation law, source, ledger or Phase-5 solver changes are authorized.
The Phase-5 knot issue is a separate upstream numerical follow-up.

Each segment uses its active cell's smooth expression for all RK stages, including
excursions. Within that cell it retains the original RHS. Crossings are located on
each solve's own trajectory. Discarded trials restore time and all components and
reset both GSL evolve and step state. Step proposals carry across boundaries.
If a crossing bracket spans adjacent representable times, retain its integrated
right endpoint and the full time bracket; there is no representable interior root.
No state projection onto a knot is permitted. Event residuals, finite time spacing,
work and state uncertainty are recorded. Unresolvable direction or event identity
fails closed. A sampled monotonicity check is a numerical diagnostic, not a theorem
excluding every possible unsampled reversal; qualification is bounded to the tested
RHS and domains and requires independent step subdivision.

Scientific observations remain passive; exact accepted endpoints return stored
states. Boundary endpoints do not enter the scientific R20 quadrature grid.
Historical reproduction uses an explicitly selected historical path, never a fallback.

## Proposed uncertainty rule

Keep d_O, D_O1, F_O, F_i and U_O from ADR-0017. Require d_O <= D_O1 and
U_O + U_event <= 0.20 F_i for every component (x dimensionless, eta in MeV).
Do not adopt the review's scalar kappa/beta formula as a universal coupled bound.
Candidate numerical evidence instead compares separately located events and full
coupled endpoints under step subdivision and an independent RK method; reference
uncertainty includes further refinement and roundoff. A componentwise full-path
difference envelope is conservative relative to assigning all error only to x,
but is an empirical bounded-domain validation rule, not a rigorous global bound.
No identical one-step witness can qualify a crossing by itself.

Diagnostics must include the state envelope in all three components. R20 propagates
absolute thermal/outgoing/chemical endpoint uncertainty on its existing grid, without
cancellation between O1/O2. Main uncertainty belongs to the main accuracy assessment.
A missing justified envelope or failed reconstruction allocation is unresolved.

## Proposed successor accuracy assessment

See PHASE6_BREAKPOINT_VALIDATION_SPEC.md. Report original BA12R state, ledger and
solver-trend gates alongside the successor outcome. Never rename a failed original
gate PASS. Breakpoint-limited status is determined from accepted/event records before
scoring differences. Where steps really refine, report contraction; where they do
not, require independent accuracy evidence including reference uncertainty.
Equilibrium-tail rejection is assessed as stability/cost/progress evidence, with no
raised legacy threshold. Local windows do not establish campaign-global accuracy.

The independent review and owner acceptance are recorded in the ratification. No merge,
full campaign, physical rate/model, cluster, superfluidity, A18 or Regime-II work.

## Exact matched-control subspace

The first candidate's generic state box exposed a diagnostic refusal at exact
chemical equilibrium: the unchanged 2e-24 MeV chemical tolerance floor generates
an Echem box of about 9.2e5 MeV, against the unchanged zero-energy diagnostic
allocation of 2e-12 MeV. This is retained as a failed development result, not
rescored as success. No accuracy budget or reconstruction tolerance is changed.

The correction is an explicit mathematical constraint, not a zero inferred from
coincident numerical outputs. `OrdinaryMatterBnvHistory` may certify identically
zero sources over its whole history (default false). The frozen BNV owner also
requires `StaticZeroSpinHistory`; its constructor already fixes the governed Urca
processes. `UrcaImbalanceFunctions::HM(0)=HD(0)=0` exactly, and zero moving-reference
source and spin imply eta_dot=0 for every valid temperature. With both initial
etas exactly zero, the chemical solution and its numerical error are exactly zero.
Each of O1, O2, the subdivided rk8pd and independent RKF45 contexts must certify
the same identity, and every local RHS call checks zero chemical stage state and
derivative. The diagnostic owner must retain the certificate. Nonzero initial
etas, unknown histories and spin-on contexts receive no exemption.

Keep U_O, U_event and the combined allocation test unchanged. Only the propagated
chemical error components become zero under this checked identity; x keeps its
full envelope. Echem and all diagnostic/R20 budgets remain unchanged. Record the
constraint in checkpoint provenance. A false certificate and an uncertified
zero-output witness are regression refusals. This narrow constraint adds no
physical source model and changes neither a state value nor a main step.

## Bounded candidate result

The final implementation passed the declared bounded successor accuracy, coupled
uncertainty, reconstructed diagnostics, R18/R20, passive-schedule and regression
checks. Original BA12R-form contraction and rejection failures remain failures.
See [qualification evidence and limitations](PHASE6_BREAKPOINT_QUALIFICATION.md).
The equilibrium-tail stability probe limits extrapolation. This outcome makes the
candidate ready for one independent review; it does not ratify this amendment or
authorize a campaign retry. The [Phase-5 knot issue](PHASE5_CSTAR_KNOT_FOLLOWUP.md)
is a separate unimplemented follow-up.
