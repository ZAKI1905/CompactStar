# Phase-6 campaign-readiness accuracy predeclaration

Status: PROPOSED numerical assessment; not a trajectory acceptance or campaign
launch authority. Written and committed before new numerical experiments on
2026-10-10. Authority: owner's readiness task, ADR-0017 Amendment 1 bounded
ratification (`PHASE6_BREAKPOINT_RATIFICATION.md:1`), ADR-0015/0016, and retained
bounded criteria (`PHASE6_BREAKPOINT_VALIDATION_SPEC.md:1`). No RHS, Cstar table,
interpolation, Phase-5D mathematics, physical source choice, or baseline changes.

## Fixed budgets and separate decisions

For y=(x,eta_e,eta_mu), x=ln(Tinf/Tref), eta in MeV, retain:

| Tier | rtol | atol x | atol eta_e, eta_mu |
|---|---:|---:|---:|
| BASELINE | 1e-7 | 1e-12 | 1e-18 |
| REFINED | 1e-9 | 1e-14 | 1e-20 |
| ULTRA | 1e-11 | 1e-16 | 1e-22 |

These are the existing scales (`tests/bnv/breakpoint_qualification.cpp:6`), now
proposed as global accuracy budgets, not interpreted as local-tolerance guarantees.
At each preselected scientific time and accepted endpoint define
D_ki=atol_ki+rtol_k max(|y_ki|,|y_ref_i|). Require
|y_ki-y_ref_i| + U_ref_i + U_output_ki <= D_ki independently for all components,
all tiers, source and control. At stored endpoints U_output=0; at observations it
includes reconstruction and event uncertainty. Do not count identical errors as
independent cancellation. Reference uncertainty must be <=0.05 D_ULTRA per
component; reconstruction retains its independent <=0.20 F_i gate. No larger
floor is introduced when roundoff dominates. Insufficient resolution is UNQUALIFIED.
Also evaluate temperature using the exact exp transform: relative uncertainty is
exp(U_x)-1, not a linearized approximation when that would underbound it.

For all 11 power channels in the existing ledger, use the same gross noncancelling
G_P definition (`tests/bnv/breakpoint_qualification.cpp:146`). Require error plus
reference and output uncertainty <=1e-7 G_P for BASELINE, <=1e-9 G_P for REFINED
and ULTRA (the latter two retain the accepted bounded scale). Reference allocation
<=5e-11 G_P. Report every channel and worst point, not just the maximum. The prior
control ratio 0.936051 and chemical ratios near 0.93 receive no waiver or rounding.

For matched source-minus-control differences at common physical times define
D_delta_ki=D_source_ki+D_control_ki. Require discrepancy from the independently
formed reference difference plus both reference/output uncertainties <=D_delta.
Do not use relative error divided by a vanishing signal. Separately require
U_delta <=0.10 |delta_ref| for any claimed resolved nonzero signal; otherwise
label that signal UNRESOLVED and publish its absolute uncertainty. Do not assume
correlated error cancels without evidence. Apply this rule to Tinf, both etas,
power channels, and finite-interval accumulated energies. Exact-zero control
chemistry is permitted only under the accepted owner-checked invariant subspace.

R18: |residual|/G_P <=1e-10, with the direct-energy identity evaluated from all
independent references as well as candidates. R18 is algebraic bookkeeping, not
an independent trajectory-error estimator. Retain all diagnostic propagated boxes,
including Echem, Cstar, actual/equilibrium powers and luminosities.

R20: retain the existing finite-interval normalization N_R20 and scientific-grid
composite trapezoid. Require (|residual|+U_main+U_reconstruction+U_event+U_quad)/N_R20
<=2e-4; separately U_reconstruction/N_R20<=5e-6 (event contributions included).
No event states may enter the quadrature grid. Quantify U_quad from independent
nested scientific grids h,h/2,h/4, including a tighter witness h/8. Use twice the
maximum adjacent discrepancy rather than an assumed order divisor. Require the
last discrepancy to contract <=0.5 unless all differences are below independently
computed rounding bounds; otherwise refine under a new declared extension or mark
UNQUALIFIED. Publish each term without sign cancellation, including absolute
endpoint thermal/chemical energy errors and trapezoid-weighted power errors.
Accumulated source-minus-control energies sum both uncertainty budgets.

## Independent reference and accumulated error

References must start from the same original physical initial state, owning their
own events and caches. Use rk8pd at rtol=1e-14, atol=(1e-19,1e-25,1e-25), then
rtol=1e-15, atol=(1e-20,1e-26,1e-26) at half and quarter maximum steps; independent
RKF45 at the tighter tolerances and quarter/eighth maximum steps. Initial
reference cap for each bounded window is duration/256, further limited by the
stability ceiling. Report the finest rk8pd; set U_ref componentwise to
2 max(all adjacent finest-refinement differences, independent finest-method
difference) + a conservative roundoff allowance >=64 ulps of the state.
The 64-ulp minimum is not an accumulated-roundoff proof: include operation-count
or high-precision witness evidence for the actual path, or retain UNQUALIFIED.
Require refinement discrepancy <=0.5 of its predecessor above the roundoff
allowance, and the 0.05 D_ULTRA allocation. Degenerate identical capped sequences
are not refinement. Record step identities/counts and explicit cap subdivision.
Endpoint observations must not reset the numerical integration.

Global assessment requires either continuous independent reference trajectories
from t=0 across the full authorized interval, or a verified error-transport bound.
For the latter, on each smooth cell propagate a componentwise enclosure
E_(n+1) <= |Phi_n| E_n + L_n + R_n + V_n, where Phi encloses the full coupled
variational flow over the state-error tube, L is a verified local defect bound,
R bounds roundoff, and V bounds independently located event error. Compose these
bounds through every prefix, starting with the initial-state uncertainty. Bound
Jacobian variation and nonlinear remainder over each tube; sampled eigenvalues,
embedded RK estimates, or endpoint finite differences alone do not bound Phi.
Include event saltation if f differs across a surface; for this continuous
piecewise RHS, prove continuity/transversality before using identity saltation.
Near tangency retain fail-closed behavior and unresolved bounds.

Gronwall/logarithmic-norm bounds may be conservative; a bound too large is not a
license to discard it. Multi-window warm starts from an archived candidate
state set E_start=0 only for local defect experiments; they cannot erase prefix
error. Local evidence can falsify global readiness, never prove it. R20 passing
cannot bound individual state errors through cancellation. A global assessment
must report accumulated uncertainty and maxima over the whole trajectory; claims
between sampled points additionally require a continuous enclosure or residual
bound. Otherwise restrict the claim to the declared output grid.

## Bounded falsification experiments and stopping rules

Keep historical BA12/BA12R FAIL and compute their original c1/c2/c3 and 10%
rejection gates separately. A successor PASS never renames a historical FAIL.
Use the four existing intervals from `PHASE6_BREAKPOINT_QUALIFICATION.md:126`,
plus matched source/control solves starting from common t=0 for the bounded
0..462269531250 s interval. Warm source/control/tail windows retain their local
qualification label. Each policy/tier must use the same physical start as its
references. No full 500 kyr solve is permitted by this task.

Compare unchanged RKF45 with a fixed ceiling independently derived below, then
half and quarter that ceiling, at all three tiers. Compare against independently
located rk8pd/RKF45 references. Keep all failed runs. Require the state/ledger
budgets above, rejection fraction <=10%, finite progress, exact endpoints,
bitwise schedule invariance at every tier, and event/reconstruction allocations.
Measure main, reconstruction, diagnostic and reference wall times separately.
No automatic tolerance retuning or budget relaxation. If explicit capped RKF45
passes, no stiff method is justified by this task; if it fails, preserve the
failure before considering any separately governed solver proposal.

## Stability-ceiling derivation fixed before experiments

Derive the stability polynomial from the actual GSL Fehlberg tableau, cross-check
it against the archived direct-GSL amplification matrices, and find the first
negative-real stability boundary (including the embedded estimator polynomial).
Use local full and chemical-response Jacobians, with finite-difference refinement,
at both archived tail endpoints as screening evidence. Candidate fixed ceiling:
0.5 times the smaller stability radius divided by an upper estimate of the largest
chemical decay rate, inflated by twice the measured Jacobian refinement difference
and a 1% screening margin. Compare explicitly with the review's 4.8e10 s scale.
No endpoint sample becomes a global Jacobian bound. Require subsequent nonlinear
window runs and refreshed Jacobians along every accepted state (with tube bounds
for a global claim); examine the scaled amplification norm as well as spectral
radius because nonnormal matrices can amplify perturbations transiently.

With unavailable qualified dependencies, authenticated archived-matrix analysis
and analytic linear experiments are diagnostic-only. They may nominate a cap but
cannot qualify nonlinear BNV stability, runtime or energy bookkeeping. Production
policy remains unchanged until qualified. Edge-case behavior is preserved unless
a narrow change has independent tests and separately recorded acceptance status.
