# ADR-0021: BA13 historical reproduction and unresolved uncertainty completion

Status: PROPOSED FOR OWNER REVIEW; no option selected or implemented.

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

Leave platform agreement and scientific certificate qualification UNQUALIFIED;
retain passing static production, exact native repeats, independent refits and
runtime refusal controls as separate findings. Do not attach historical Mac
certificates to native stars. No promotion, campaign, extra certificate-dependent
pilot or nonlinear stability-ceiling claim follows from this recovery.

The next bounded step is owner review of the uncertainty construction/envelope
interpretation, followed by the predeclared static error-budget work. Static
certificate qualification must precede the already scoped certificate-dependent
pilots; nonlinear stability-ceiling qualification remains a separate task.
