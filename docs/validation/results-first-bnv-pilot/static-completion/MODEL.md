# Accepted joint propagation and retained-data completion model

Authority: OWNER_APPROVAL.md and ADR-0022. This document fixes the implementation
model before cluster assessment. No stellar integration or trajectory is part of
this calculation. Governing definitions: accepted preflight §11; historical
`tests/bnv/compare_candidate.py:131`; accepted review Appendix A at
`../evidence/INDEPENDENT_REVIEW.md:550`.

## Variables and the accepted correction

Q is the physical equilibrium output on the accepted cut law; I_X is the
piecewise linear X2 reference at achieved B, and S is the historical component
scale. q=Q/S, z=I_X/S, r=q-z. Q-I_X is an additive physical error; its magnitude
is absolute error; (Q-I_X)/I_X is fractional error; r is normalized additive
error. For signed relative drifts, d_q-d_z=(r_i-r_0-d_z*r_0)/(z_0+r_0).
An additive constant therefore does not disappear from baseline normalization.
Thermal outputs are actual temperatures at each fixed Tinf knot, after the
nonlinear envelope relation, not a temperature-state logarithm.

M=I-H is the historical weighted linear-fit residual operator: 11 fitted rows
(10 through 20), 21 evaluation rows, unchanged component weights. M1=Mx=0.
For the accepted |q''|<=K on a window of length L=1e-6, Taylor expansion about
the midpoint gives ||Mq||<=f*K*L^2/8, where f=||M||_infinity. Consequently

    ||M(q-z)|| <= F = f*K*L^2/8 + ||Mz||.

This bounds the TOTAL reference fit response, not the unseen common component
alone. It permits nonzero affine B dependence. It is analytic propagation
conditional on C-CONT, not a new empirical bound on shared error. C-CONT remains
conditional despite the retained holdouts.

Only envelope and Tsurface_inf use F. For all groups retain
e=10*max|P-I_X|/S, except the accepted Ltilde_Me allowance 10*9.1e-12.
The shared increment remains 2e at non-origin rows, zero at the origin.
Other groups retain f*e fit propagation. With the prior finite-ledger total A
and baseline error b, the amended necessary temperature check is

    [A*(1-b) + 2e*(i>0) + F] / (1-b-e).

The old finite fit entries remain conservatively charged. F is a bound on the
combined reference fit response; retaining those entries overcounts some known
fit uncertainty, safely. The original pointwise assay is replayed without edits;
its 42 failures and 355 expectation violations remain immutable. The amended
partial check must report 353 expectation violations separately. Expectation
violations are diagnostic, not replacements for historical thresholds.

## Completion of the ledger

Every historical row and term remains unchanged in the archived record. The
separate restricted ledger itemizes the following dispositions and provenance.

* Existing finite terms retain their values. Their standard/condition mappings
  are S and C-REF/C-TAN/C-CONT/C-INPUT/C-CUT/C-MUN as applicable. Source comparisons
  are measured differences, not bounds on common error. Old normalization is
  replaced only in the NEW ledger by propagation of the complete numerator.
* Stage plus arithmetic: one pointwise allowance of 15 times the finest GL4 tail
  under C-REF, split equally for bookkeeping. Propagate each through the paired
  difference and unchanged fit; include it in the baseline denominator. For
  tangents the tail is measured on the correlated ratio itself.
* EOS and centre: previously charged estimates are the full estimates under S
  and C-REF. The additional parameter remainder is zero because it would be a
  duplicate under those accepted models, not because missing error is zero.
* C-TAN: retain the finest ratio tail and the first narrow-width tail already
  charged; add the second narrow-width tail and the retained cross-geometry
  residual. Keep EOS/centre sensitivity on the ratio. Retain denominator
  separation at relative error <=1e-14 and conservatively include the division
  correction U*1e-14/(1-1e-14) for this accepted tangent sum U. Reconstruct the
  raw ratio and identity delta_t=(delta_A-t*delta_D)/D_new at high precision.
  No count discrepancy is divided by a small B spacing to invent tangent error.
* C-INPUT: retain the observed first-order coordinate and floor responses and
  baseline normalization. Add max_gap[(coordinate_estimate/achieved_spacing)^2
  * max|P-I_X|/S], as accepted. Historical weights remain fixed. No fit covariance
  is treated as physical uncertainty.
* C-CONT: preserve existing sample-gap remainders and topology hypotheses. For
  the endpoint sliver, bound the derivative by the last reference secant plus
  its endpoint bands and K*h/2; add slope*sliver+K*sliver^2/2. Do not treat a
  nominal endpoint as attained. Any missing endpoint input remains a dependency.
* C-SUPPORT gives zero numerical topology drift under its accepted no-crossing
  condition. P2_average remains the exact identity. C-MUN applies only to the
  equilibrium reference. Metric rows remain CONSUMER_SCOPE_ONLY, never nodewise
  PASS, with their complete historical record preserved.

## Moving-cut evidence stop

C-CUT supplies |delta ln Ps|<=2.1e-12. It requires propagation through an
observable-specific response at fixed B, not multiplication by the along-path
B slope. Neither a small pressure residual nor C-CONT's one-dimensional
curvature supplies that response. `continuation.json` records cut motion and
residuals; `m1.json` contains integration/EOS/centre/stencil estimates;
`m3.json` contains B-direction stencils, and `diagnostics.json` contains tangent
cross-geometry residuals. None contains a validated finite-cut sensitivity
envelope for each output at fixed B over the cut perturbation interval.

The assessor must not invent that input. It records
`moving_cut_input_remainder = null` with the missing input described explicitly,
and separately records the missing baseline and any endpoint dependence.
Deriving a new two-parameter sensitivity bound from terminal-state identities
would be a separate supported analytic calculation, not a number already present
in this ledger. Finite samples or an unbounded first derivative alone do not
justify a uniform finite-perturbation bound. No new assumption is accepted here.

For an approximate relative row, with unnormalized numerator N and total
baseline error b_total, normalization adds N*b_total/(1-b_total). If N or
b_total has an unsupported contribution, full normalization and full utilization
are null. Report the supported partial budget separately. Absolute t rows have
S=1 and no observable-denominator normalization; their missing cut/fit term still
prevents qualification. No partial sum is labeled a complete total.

## Verification and execution

Authenticate the accepted registry and all retained inputs; use fresh external
output storage. Reconstruct the temperature projection with independent
high-precision normal equations at 70 and 90 digits and verify affine identities,
norms and F. Recompute raw tangent ratios at both precisions. Preserve the
historical matrix and original 42 failures exactly. Refuse missing terms,
unseparated denominators, unauthorized assumptions or changed authority inputs.
Run the amended matrix and relevant algebra/refusal tests on one EKU Slurm CPU,
with a 60 CPU-second static-assessment ceiling. A blocked matrix is retained and
exits nonzero; it issues no operational certificate or trajectory permission.
