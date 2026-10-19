# ADR-0021: BA13 historical reproduction and unresolved uncertainty completion

Status: PROPOSED FOR OWNER REVIEW. Option A investigation authorized and piloted
2026-10-12; no new scientific method ratified and no certificate accepted.

## Context and authority

The owner authorized recovery of the original static BA13 producer and a Linux
candidate, preserving its fitting semantics. GOVERNANCE §§2–5, AGENTS §8,
INV-14/15, the accepted Phase-6A-1 preflight §11, ADR-0015/0016/0017 and proposed
ADR-0020 govern this work. This proposal records the remaining numerical-method
and scientific-semantic ambiguity; it does not modify accepted owners or policy.

The recovered fitter is Git blob `6790b8983d63fb908af43bf2c8e30e3a5fed7e94`,
SHA-256 `9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370`.
It passes the native 21-star data with maximum utilization 0.54075347985134414.
All 10303 component fits pass a separate 80/112-digit refit. The scientific
production owners, sample grid, fit window, thresholds and original source bytes
are unchanged. Reproduction success is distinct from independent qualification.

## Findings that prevent promotion

1. The historical fitter supplies zero errors for radius, envelope temperature,
   neutron chemical-potential profiles and metric profiles. The producer supplies
   `64*epsilon*value` for heat capacity, surface gravity and redshifted temperature.
   These do not establish structural, EOS or thermal refinement uncertainty.
2. Ltilde and direct-reference errors are absolute order-16/order-32 differences.
   Some differences are exactly zero. Agreement of those two rules alone does not
   bound common background, EOS, interpolation or quadrature error.
3. Both historical and native per-star chemical adapters inherit B0 background
   and center characterization. Native actual-cut tail witnesses were recomputed
   by unchanged Phase-5C mathematics. This does not establish independent
   background refinement bounds at all 21 target stars.
4. The historical envelope is cumulative maximum sampled drift plus the largest
   accepted residual. Numerical errors enter thresholds and fit covariance, but
   there is no separate numerical-error term in `drift_bound`. Preflight §11 says
   the full drift plus numerical and nonlinearity terms must remain within the
   threshold. Whether the existing threshold formulation discharges that phrase
   requires owner interpretation; this recovery does not silently choose one.
5. The historical weighted fit floor `64*epsilon*max(1,abs(value))` operates in
   each component's original units. For tiny dimensional Z/Ltilde values the floor
   dominates the supplied errors. It protects conditioning but is not a physical
   uncertainty estimate. Its covariance must not become a platform tolerance.

The measured platform differences are small, but several exceed the sum of the
declared errors. For example, heat capacity exceeds that sum at 214/273 knots;
N_mu at 15/21 stars; surface gravity and radius at all 21 stars. Exact-zero P2 and
support topology agree. See the componentwise report under
`docs/validation/ba13-recovery/evidence/platform-comparison-summary.json` and the
full external candidate archive. No global relative tolerance is substituted.

## Alternatives requiring a decision

A. Retain this output as a historically faithful static candidate only. Complete
independent per-quantity, per-target refinement/error propagation using current
owners and the existing Phase-5B/C/D validation methods. Keep nominal production
grids, 21 targets, 13 knots, thresholds and fit window fixed. Predeclare the error
construction and clarify the envelope interpretation before a successor fit.

B. Accept a narrower claim of historical-method reproduction, explicitly excluding
the missing numerical and platform-accuracy claims. This does not authorize a
qualified Linux certificate or certificate-dependent pilots under ADR-0020.

C. Amend the fitter's numerical treatment (including dimensional weighting or
envelope propagation) as a separately reviewed method change. Such output must
have a distinct identity and independent baseline; it cannot be described as the
original fitter's output. No such amendment is implemented here.

## Current disposition and next bounded step

The subsequent Option-A investigation is recorded in
[the uncertainty review](../validation/ba13-uncertainty/SCIENTIFIC_REVIEW.md).
It gives disposition **B: specific additional bounded evidence required**:
64 conditional radial holdout predictions pass; the EOS 32768 holdout refuses
on four rounded equal-value pairs near neutron onset; heat-capacity radial
convergence is irregular; Ltilde_Me changes about 6.35e-5 under radial refinement
although fixed-profile quadrature agrees near machine precision. The cache law
explains a shared absolute Cstar bias up to 6.56e-4, whose exact common factor
cancels in relative depletion ratios. Original covariance is unit dependent and
remains unsuitable as physical uncertainty. The 462-entry matrix contains only
21 exact-P2 PASS entries, with all other full-budget entries unresolved.

This evidence does not select Option C or amend the original fitter. Independent
review must resolve the joint error/envelope interpretation, then authorize any
required method change separately. Next numerical work is bounded onset/EOS
refinement and per-target background/tangent error completion, subject to the
predeclaration and its continuation estimate. Historical outputs and all 10,303
refits remain unchanged.


Leave platform agreement and scientific certificate qualification UNQUALIFIED;
retain passing static production, exact native repeats, independent refits and
runtime refusal controls as separate findings. Do not attach historical Mac
certificates to native stars. No promotion, campaign, extra certificate-dependent
pilot or nonlinear stability-ceiling claim follows from this recovery.

The next bounded step is independent owner review of the new uncertainty evidence
and envelope interpretation, followed by the specific bounded continuation above. Static
certificate qualification must precede the already scoped certificate-dependent
pilots; nonlinear stability-ceiling qualification remains a separate task.

## 2026-10-16 headroom assessment — proposal remains unratified

The user authorized candidate implementation of Option A and interpretation (2)
with all original right-hand thresholds frozen. That authorization does not
ratify either interpretation. The executed [headroom assessment](../validation/ba13-headroom/REPORT.md)
and [complete 21-by-22 matrix](../validation/ba13-headroom/MATRIX.md) supersede the
older *next numerical step* recommendations above, without altering any historical
disposition, fitter output or independent-review evidence.

The separate schema `compactstar.ba13.headroom-uncertainty-assessment.v1` preserves
all 462 historical drift/threshold fields. It implements headroom allocation,
deterministic paired corrections, B0 normalization, functional-error-aware common-B
propagation, changing inverse Jacobians and Neumann remainders, and a qualification
layer around the unchanged fitter. Three-target X2 extensions and all 21 retained
nominal profiles were evaluated. Exact replay and 14 focused fault/algebra tests
pass; the signed defect-transport correction fails its radius-direction check.
Evidence: `docs/validation/ba13-headroom/evidence/audit.json:1`, `pilots.json:1`,
`matrix.json.gz`, and `REPORT.md:1`.

Recommendation **B**: numerical and interval-coverage gaps prevent closure.
The matrix records 21 exact-P2 PASS and 441 UNQUALIFIED rows. The remaining work is
specifically quantity-dependent defect/roundoff propagation, remainders beyond
measured EOS/center/tangent perturbations, variation of inherited chemical
background uncertainty, uncertain achieved-B design/reference propagation, and
continuous paired-error coverage including production-node jumps. This is not
an added requirement for a formal proof of the entire stellar model. Nor does
observed small paired error or an enlarged numerical estimate justify increasing
any frozen threshold. A whole-cell L support allowance is too loose; that does
not establish a necessary production change.

The next action is **one independent review of the completed assessment and its
explicit residual gaps**, followed by owner adjudication. Do not reopen a general
oracle comparison, EOS/jitter sweep or tolerance ladder as the default next step.
No acceptance of BA13, ADR interpretation, narrower conditional certification,
production owner change, trajectory or certificate promotion follows automatically.

## 2026-10-19 closure study — proposed status unchanged

The authorized [closure study](../validation/ba13-closure/REPORT.md) corrects the
independent review's reference, evidence-class, tangent-coordinate and bookkeeping
defects and retires the inconsistent first-order defect transport. The updated
462-row sidecar retains 21 exact-P2 PASS, 0 FAIL and 441 UNQUALIFIED. The muon
sector remains unqualified: tangent reference, all-target transfer and continuous
curvature evidence are incomplete. Small measured errors do not establish closure.

[The owner memo](../validation/ba13-closure/OWNER_DECISIONS.md) records five pending
decisions: interval remainder, evidence standard, certified object, record format
and target-to-target magnitude transfer. No decision is ratified here. Option A
and interpretation (2), the original fitter, historical thresholds/certificate,
production owners and all prior failures remain unchanged. No successor fitter,
runtime certificate, trajectory or promotion is authorized by this assessment.
