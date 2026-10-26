# BA13 evidence completion — prospective declaration, 2026-10-26

UNVERIFIED SCIENTIFIC CANDIDATE. Entry is
093d121135e35bc50ca2a807470906bcb28ef5c3, clean branch
candidate/linux-ba13-closure-20261019. New descendant worktree is
worktrees/ba13-evidence, branch candidate/linux-ba13-evidence-20261026.
All other local branches are ancestors of entry. The independent review's full
SHA256SUMS verifies; its authentication program passes 85/85, including 884
input entries. Historical FAIL, INCONCLUSIVE and UNQUALIFIED records are immutable.
In particular no historical M1 or M2 is rescored by the successor rule.

Authority: this task's explicit bounded-investigation authorization;
GOVERNANCE §§1–5; INV-02/03/04/06/14/15; accepted ADR-0009/0013/0015;
PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT §§8–9,11.
ADR-0021 Option A and interpretation (2) remain proposals. Classes: validation
numerical-method, engineering, documentation, generated evidence. No production
method, fitter, threshold, fit window, physical model or runtime record changes.
The Mac literature library is unavailable; the existing contracts, source
equations and authenticated review supply this investigation's derivations.

## Reference and error construction

P is the retained production sample, R its reconstruction, X the existing X2
method, Q the exact governed finite-pressure equilibrium response. At the
production achieved coordinate x=(B_P-B_P0)/B_P0 use

    P-Q = (P-I_X) + (I_X-I_Q) + (I_Q-Q).

I is PIECEWISE LINEAR interpolation at retained targets 0,10,20 in achieved B,
not target index or reconstructed B. Convex interpolation weights carry endpoint
X2 uncertainty without amplification. This deliberately uses the exact linear
interpolation remainder theorem, not a quadratic interpolant with an unexplained
linear-error allowance. Each target retains its own historical pressure cutoff.
For continuous interpretation propose log-linear interpolation of those cutoffs
in achieved x, retaining the measured deviations as an unresolved input term.

Record D_i=max(j<=i)|[(P_j-I_Xj)-(P_0-X_0)]| componentwise with historical scales.
Production stage, EOS-table, arithmetic, event, support-edge and envelope-node
selection errors at the samples are included in this measured comparison; they
are not replaced by X2 self-errors. R is diagnostic only. X2 errors attach ONLY
to X2, including root residual plus B-functional uncertainty and the positive
B_H hypothesis. Do not add reconstructed-B displacement c a second time to P-X.
Carry reference integration, EOS, center, stage/arithmetic finite sensitivities,
stencil error, endpoint alignment, reference interpolation, denominator error,
fit response, and continuous gaps separately. Unknown remainder values stay null.
No arbitrary unseen-bias allowance is inferred from the review's illustrative
1e-11. Cross-method discrepancies are observations, not bounds on common bias.

Frozen curvature proposal: K is the maximum of the PRIOR closure's K_hypothesis,
its wider_stencil_K, and ten times the prior review's target_curvature_max, per
quantity. This inherits an explicitly CONDITIONAL continuation proposal, not a
proven supremum or a newly tuned multiplier. Require positive sampled B_H,
fixed topology, and width/target curvature agreement within factor 2. The
interpolation envelope is K |x-a||x-b|/2 (maximum K h²/8); the sample-gap
envelope is K (gap+recorded_coordinate_uncertainty)²/8. These conditions do not
prove unsampled continuation. Any failure stops dependent numerical expansion.

For fit propagation use the actual historical weights and all 21 rows of I-H,
including the inner half. Charge deterministic (I-H)(P-I_X), operator-norm
propagation of reference/interpolation/alignment inputs, retained measured
design/floor responses, and explicitly missing nonlinear/design remainders.
No covariance is counted as physical uncertainty. Denominator excess uses the
FULL numerator including fit and inter-sample terms. Absolute tangents have S=1.

## Tangent experiment and quantitative gates

Run exactly twelve fixed-H evaluations: targets 0,10,20, signs ±, N=1024,2048,
H=retained X2 N4096 central H times (1±0.000025). Retain N512 endpoints and the
three historical widths 0.0001,0.00005,0.000025. Use 60-digit Decimal parsing.
T_N=(N_plus-N_minus)/(B_plus-B_minus), with correlated numerator AND denominator.
Require same-sign nonzero increments and p=log2(|d1/d2|) in [3.5,4.5] for each
tangent/target before calling |d2|/15 a validated numerical estimate. No floor,
absolute-value sign repair, added resolution or result-based multiplier.
N512 reference error, if needed, is |T2048-T512|+|d2|/15, not independent endpoint
errors. Width ratio must stay in [3,5], with narrow-width tail |delta_width|/3.
Report retained central cross-target secant agreement with extrapolated stencil
means and the curvature mean-value correction; this is a distinct consistency
test, not an independent physical solver. Retained X1 comparisons and anchor
residuals remain separate independent-method checks, with no invented tangent
accuracy from scalar agreement. Report common-B/root terms in absolute t units.
If any tangent order/width gate fails, stop new numerical work after the twelve
declared evaluations; retain finite increments as observations, not a tail bound.

Successor M1: discrepancy alone above the fixed 23% representation share is FAIL.
A justified complete total enclosure <= share passes the budget test; otherwise
refuse. Reference/discrepancy ratio is only a resolution diagnostic. A partial
sum below share is not a qualification and not evidence of a complete enclosure.
Other per-category planning shares remain exactly the historical FRACTIONS.
Estimated terms must be <=10% of their allocated share to stop refinement for
decision value; measured discrepancies are not refined to satisfy this rule.

## Prospective holdouts and stopping

Only after the tangent experiment gates pass, run target 5 then target 15,
each at its OWN retained B_solved and last-profile-node pressure. Use existing
X2 N2048 root mathematics, then the unchanged passive-functional executable at
the attained H. No new stellar solver, root ladder or N4096 evaluation. Each
holdout tests the 15 non-tangent scalar/vector X2 groups; it cannot validate the
tangent interpolator. Freeze a numeric criterion file before the first holdout:

    |X_holdout-I_X(x_holdout)| <= K|x-a||x-b|/2
       + convex endpoint numerical band + holdout numerical band
       + endpoint/holdout common-B alignment estimates.

For the predictive holdout band use the max retained N2048 integration+EOS+center
estimate at 0/10/20; record the all-target transfer hypothesis explicitly.
Arithmetic/stage ablations are finite sensitivities, not an exact-zero claim.
No common-bias floor is added. A PASS is a prospective interpolator consistency
result only. Any violated criterion stops reference transfer and target 15 if
not yet run. No fallback refit, new target or wider criterion.

Run on one pinned Slurm core, single-thread BLAS. Initial estimate 40–60 CPU s,
hard cumulative numerical ceiling 180 CPU s including analysis and tests;
build/hash-only work separately receipted. Check cumulative usage before each
subprocess, cap each at available reserve, and retain failures/partial outputs.
No automatic scope expansion. Stop at numerical gate failure, unjustifiable
reference uncertainty, failed continuation, or missing authority at qualification.
The user authorizes candidate experiments without ratifying their evidence class;
that is not authority to issue a scientific PASS. Finish static repairs, matrix,
tests, audit and reporting after a numerical stop.

## Evidence, profiles, and deliverables

Each nonzero contribution needs explicit evidence class, source hash/provenance,
and dependencies. No default analytic_bound. Numerical standards need recorded
acceptance; conditionals require itemized condition-to-authority acceptance.
The registry starts empty. Missing input or missing authority refuses qualification.
Analytic propagation does not upgrade uncertain inputs. Preserve old outcomes in
each successor row and include remaining headroom as partial whenever incomplete.

After the above, inspect mu_n_profile's actual field definition against the
equilibrium redshifted identity. Report its definition/authority separately from
the measured profile comparison. metric_structure requires a nodewise reference
for its full r,m,nu,lambda claim; P-R alone is insufficient. Do not implement an
observer while reference/continuation gates are unresolved. If prerequisites
close, a separately committed benchmark and coverage declaration must precede
observer runs: quadratic interior dense output does not inherit endpoint GL4
order. No observer development after a failed benchmark.

Re-run the original 30 tests and 113 audit checks without altering their recorded
historical results; add true mutations of the new assessor's reference binding,
alignment and evidence classification. Deliver 462 rows, explicit partial versus
justified totals, structured milestones, report, append-only ADR proposal, source
and artifact hashes, measured Slurm receipts, and verified incremental bundle
against entry if GitHub write authentication is absent. No master merge,
certificate promotion, BNV trajectory or full campaign.
