# BA13 headroom assessment — candidate method

UNVERIFIED SCIENTIFIC CANDIDATE. Entry is
`89b629dd9961c29e127badbd723455d00dc38553`. Numerical-method, engineering,
documentation and generated-evidence change; no production-method change.
Authority: GOVERNANCE §§2–5; INV-02/03/04/06/14/15; accepted ADR-0009,
ADR-0013 and ADR-0015; accepted Phase-6A-1 implementation preflight §11
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:655`).
ADR-0021 Option A, interpretation (2), is evaluated as PROPOSED, not ratified.

The question is whether historical frozen production quantities approximate the
depleted sequence over the entire specified interval within the original 462
right-hand thresholds. Those thresholds are immutable, including historically
error-dependent thresholds. The reference is the original native nominal data,
not an independent solver substituted into the certificate. The historical
fitter, its weights, covariance, values and output fields are unchanged.

Evidence vocabulary: **analytic_bound** is an algebraic upper bound subject to
explicit input hypotheses; **validated_numerical_estimate** is a conservative
numerical estimate with independent validation; **conditional_estimate** depends
on stated unvalidated continuation/regularity assumptions; **unresolved** means
no eligible estimate. An algebraic propagation theorem does not validate its
inputs. The contract allows numerical validation; it does not require a formal
computer-assisted proof of the entire stellar model. Conversely, sampled defects,
contraction and small root residuals alone do not establish error bounds.

## Budget before additional computation

For each approximate row let h = historical threshold minus historical envelope.
Allocate h by fractions: EOS/equilibrium 0.08; center 0.04; surface/events 0.06;
integration/propagation 0.10; implicit stages 0.02; floating point 0.02;
common-B alignment 0.10; production representation 0.23; differentiation/inverse
0.10; B0 normalization 0.05; nonlinear fit 0.10; inter-sample coverage 0.10.
These are engineering targets, **not claimed achieved errors**. Unused shares
remain reserve; no post-result enlargement of a right-hand threshold is allowed.
Exact P2 and topology have separate logical tests and no fractional-error budget.
The smallest t_mu margin is about 2.11e-7 in absolute component units.

Use X2 N=128 for an initial accuracy/cost benchmark; compare to already retained
fine solutions with their conditional tails, not exact truth. Add only passive
functionals to the same GL4 stages for particle numbers, both L channels, response
G and direct P0/P1 averages. Test at 0,10,20, then assess retained nominal profiles
at all 21 points. Local fixed-H offsets for derivatives fill a specific sensitivity
gap; they are not a tolerance ladder or a new equilibrium solver. Compare coarse
and fine stencil widths and record attained B coordinates. Retain any failure.
Use existing center/EOS/arithmetic ablations and independent response owners.

## Deterministic propagation

If X_i = y_i + e_i and X_0 = y_0 + e_0, represent shared deterministic errors
once. For S0=|y0|, |e0|<=u0<S0 and |ei-e0|<=upi,

    |Xi-X0|/|X0| <= (|yi-y0|+upi)/(S0-u0).

The same componentwise formula applies to the historical fixed profile scales;
for Z preserve the historical maximum-entry row scale (not a new row-sum scale).
Absolute tangents have S=1 and no scalar-ratio normalization. Common additive
errors cancel from differences only if their sensitivities agree. The thermal
cache multiplicative factor cancels exactly only at equal Tinf knots with the
same linear entropy law and the same cache grid. It remains in absolute Cstar.

For A=G^-1, E>=|delta G|, K=|A|E and rho(K)<1, the entrywise inverse majorant is
(I-K)^-1 K |A|. The linear part is K|A| and the remainder is their difference.
For paired inverses use separate A_i and A_0, including both remainders; subtracting
delta G before inversion is invalid. Z=L^T A L uses |L| on the majorant.

If Bhat(Hhat)-Btarget=r, the B functional error is at most uB and B'(H)>=m>0
on a containing root neighbourhood, |Hhat-Htrue|<=(|r|+uB)/m. Each quantity then
needs its own derivative bound. Secant sensitivity is numerical evidence, not a
proved lower bound on m. This applies also to production baryon quadrature.

For the unchanged historical fit projection H, source perturbations propagate
as (I-H)e. Source uncertainty is not added again as fit covariance. Achieved-B
uncertainty requires design/weight perturbation or a derivative enclosure.
All 21 residuals, including the inner half, must be evaluated. Between knots a
normalized second-derivative bound K yields K*dx^2/8 beyond the linear interpolant;
without K or another justified modulus the interval is UNQUALIFIED. A narrow
bump invisible at every sample is a deliberate counterexample. Production node
selection can jump; do not assert smoothness through it. Bound the jump separately.

## Execution and completion

New artifacts live under `ba13-headroom` outside Git; compact evidence under this
directory. Substantial work runs on one pinned CPU per Slurm array task, 1 GiB,
with process CPU/wall/RSS/return-code receipts. Initial numerical ceiling 1800
CPU seconds and per-case limit 180 seconds. Benchmark before array expansion.
No 21-star static rerun, EOS/jitter sweep, RKF45 ladder, trajectory or promotion.
The final matrix must have exactly 462 rows, with null for every missing term,
explicit reason and PASS/FAIL/UNQUALIFIED; finite partial totals are not passes.
One independent review of this construction is the next delivery action.
