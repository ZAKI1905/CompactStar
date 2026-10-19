# BA13 closure equations and evidence rules

UNVERIFIED SCIENTIFIC CANDIDATE. This implements the authorized qualification
study of ADR-0021 Option A, interpretation (2); it confers no owner acceptance.
The immutable prospective limits are in [PREDECLARATION.md](PREDECLARATION.md).
Authority is GOVERNANCE §§1–5, INV-02/03/04/06/14/15 and the original Phase-6A-1
preflight §11. Original source models, cutoffs, fits and thresholds are retained.

## Objects, units and decomposition

Let P denote historical production, R a deterministic reconstruction on that
production profile, X the independent X2 reference, and Q the physical equilibrium
response. R is never an independent stellar reference. All differences below
use the original component scales S, including maximum-entry Z row scales.
Tangents use S=1 in absolute count/count component units. For Cstar, the same
known cache factor at each equal temperature knot cancels from the relative
drift; the comparison of heat-capacity coefficients removes precisely that
factor. All 13 original knots remain in the matrix and component data.

At a fixed coordinate, the identity is

    P-Q = (P-R) + (R-X) + (X-Q).

For reconstructed B error e_B=Br-Bp, define c=Q_x e_B/B0 and a=R-X-c.
Then a paired error decomposes into the signed paired representation correction,
aligned structural error, common-B contribution, and reference error:

    (P-Q)i-(P-Q)0 = Δ(P-R) + Δa + Δc + Δ(X-Q).

Tangents have no separate R: use a=P-X-c. The reported structural magnitude is
max over targets 0,10,20 of |a|. Twice that magnitude is an INDEPENDENT-ENDPOINT
transfer hypothesis for unsolved targets, never a validated supremum. At target
zero the pair is identically zero. Reference EOS/center/refinement terms describe
X-Q only; they no longer stand in for production structural error. This is the
D1 correction. See [the assessor](../../../tools/linux/ba13_closure/assess.py).

The signed common-B term is Δc, with the derivative evaluated separately for
each quantity, including all three t components. The common coordinate is
shared once. The production root residual does not replace e_B. X2 root
alignment separately uses

    |δH| <= (|computed B residual| + u_B)/m,

where m is a positive local B_H hypothesis; the study reports it as conditional.
See [m1.json](evidence/m1.json) and [m4.json](evidence/m4.json).

## Reference experiment and derivative estimates

Retained Hc and surface-pressure strings are fixed. At N=512,1024,2048 per
segment, define d1=X1024-X512 and d2=X2048-X1024. Compatible fourth order requires
same signed nonzero increments and p=log2(|d1/d2|) in [3.5,4.5]. Fine/coarse
Richardson tail ESTIMATES are |d2|/15 and 16|d1|/15. The unseen remainder is not
proved by this test. EOS96/192 and center-factor differences remain finite
sensitivities; even a zero spread would not imply a zero bound. Parsing retains
the 21-digit outputs in extended precision; the audit independently uses Decimal.

For a centered tangent, t=(N+−N−)/(B+−B−). Relative H widths 1e-4,5e-5,2.5e-5
test the h² stencil term. The narrow-stencil estimate is |t_narrow−t_middle|/3.
Without validated error correlation, the endpoint numerical contribution is

    u_t <= (u_N+ + u_N− + |t|(u_B+ + u_B−))/|B+−B−|.

Here endpoint magnitudes are CONDITIONALLY transferred from central refinement;
they are not measured offset-H refinement bounds. This intentionally retains
the denominator uncertainty. It dominates the t_mu paired discrepancy and
triggers the predeclared INCONCLUSIVE rule. Central N=4096 refinement cannot
close this offset-H/correlation question. No extra stellar evaluations follow.
Details: [m1.json](evidence/m1.json).

Both first-order signed and absolute defect transports are retired from closure.
The first-order intra-step propagator omits a leading O(h⁵) step contribution;
the review's second-order Magnus demonstration is preserved. The current
validation executable refuses the retired mode. Neither estimator supplies a
bound or an error allowance in the matrix. X2's physical solver is unchanged.

## Continuous response, production selection and normalization

With fixed topology and B_H>0, the chain rule gives

    Q_x = B0 Q_H/B_H,
    Q_xx = B0² (Q_HH B_H − Q_H B_HH)/B_H³.

For t=N_H/B_H the calculation also uses N_HHH and B_HHH from width differences.
The exact exponential manufactured test independently checks both chain rules.
The fixed-pressure H stencils do not establish the contribution from changing
the historical surface cut; that remainder is explicitly missing.

The predeclared K=10 max_three |Q_xx| is tested against narrower widths and each
of the 20 achieved-B gaps. The estimate is K*(dx+coordinate_estimate)²/8.
Agreement at three targets and multiple widths does not establish a global
supremum. No continuous bound is inferred from those finite differences.
See [m3.json](evidence/m3.json).

Sampled neutron cells and envelope indices are fixed. All 20258 interior radial
nodes are identical across 21 profiles; the final, moving surface node differs.
The sufficient no-switch condition L*max_gap/2 < minimum sampled cell margin
does not require monotonicity, but L itself must be justified. The matrix
records observed secants and the sufficient limits without equating them.
Physical equilibrium curvature receives no fictitious production-node jump.
An intermediate production-algorithm claim additionally needs adaptive step
control and terminal-event coverage. These are different certified objects.

Shared achieved-B normalization is propagated as

    x_true − x_prod = [e_B(i) − (1+x_prod)e_B(0)]/[B0+e_B(0)].

Given valid endpoint error inputs u0 and paired error up, the normalized drift
bound is (|ΔP|+up)/(S−u0). The assessor reports the denominator excess separately.
Its u0 input remains conditional; the formula does not validate that input.

## Fit and chemical inverse

The historical outer-half weights and their 64*binary64-epsilon value floor are
unchanged. The same signed representation error propagates through (I-H)e.
The study evaluates all 21 original residuals, inner-half excess, the change
from reconstructed achieved-B abscissae, and the value-dependent floor change.
No original covariance is added as physical uncertainty. Structural/reference
error in the fit inputs and unseen design remainders remain missing.
See [m4.json](evidence/m4.json).

For A=G^-1, E>=|δG|, K=|A|E, rho(K)<1, the full entrywise inverse majorant is

    M = (I−K)^−1 K |A|,
    M_linear=K|A|, M_nonlinear=M−M_linear.

At independent endpoints 0 and i, propagate |L|ᵀ(M0+Mi)|L|, with each endpoint's
own G, background, center and actual-cut tail matrices. BOTH nonlinear tails
are retained. This reproduces the review's ~5.6e-7 muon-row estimate with no
assumed signed-background correlation. An independent Neumann-series summation
checks the linear-solve implementation at every target. The per-target
background magnitude hypothesis remains unvalidated. The M5 construction is
an ALTERNATIVE to the direct P-R-X path and is not added to it, which would
count the same structural error twice. Full matrices: [m5.json](evidence/m5.json).

## Decision rules

The candidate classifier requires every required term, finite nonnegative
values, no unresolved remainder, and an explicit evidence class:

| Class | Eligibility |
|---|---|
| analytic_bound | Algebraic/established result with its actual input hypotheses |
| validated_numerical_estimate | Named standard plus externally accepted authority in the caller's registry |
| conditional_estimate | Named condition plus explicit external acceptance in that registry |
| unresolved | Never eligible |

No acceptance registry is populated in this study. A term cannot create
authority by writing an authority string. Definitional zeros carry a charged-to
pointer where applicable. Production binary64 arithmetic has no long-double
roundoff proxy. Exact identities require exact logical checks; they do not use
the scalar tolerance rule. Historical drift/threshold identity is checked
against authenticated original bytes. UNQUALIFIED rows have no total utilization.
See [the classifier](../../../tools/linux/ba13_headroom/enclosure.py),
[tests](evidence/test_ba13_closure.out) and [independent audit](evidence/audit.json).
