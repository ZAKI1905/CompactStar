# Isolating a chemical temperature contribution without a new trajectory

This comparison-bound candidate answers a narrower question than the paired
source/control separation. Define the contribution of LH by comparing the full
coupled model y with the same initial-value problem z whose thermal RHS omits
only LH/(T Cstar). Chemical kinetics and all other terms remain coupled. This is
a diagnostic mathematical definition, not a proposed alternative physical model.
The original candidate bound runs no such counterfactual trajectory and changes
no production equation. The later independent reviewer did integrate the
term-suppressed diagnostic equations; see the distinction below.
The net-beta contribution uses LH-DeltaLnu instead. The full source/control pair
continues to include both terms under ADR-0014/0015.

Write f=f0+(g,0,0), delta=y-z, and I(t)=integral g(y(s)) ds for nonnegative g.
Both vector fields have the same conservative scaled Lipschitz bound L from
PIECEWISE.md:1: omitting LH removes positive triangle-majorant terms; omitting
both LH and excess neutrino power also only removes those terms.
Gronwall gives ||delta|| <= exp(Lt) I(t). Consequently

    |delta_x - I(t)| <= (exp(Lt)-1) I(t),
    (2-exp(Lt)) I(t) <= delta_x <= exp(Lt) I(t).

Here exp(Lt)<2 over the specified horizon, so the lower bound stays positive.
For negative net-beta forcing, apply the same inequality to its magnitude and
reverse the sign. The polynomial coefficients verify DeltaLnu-LH>0 throughout
the chosen domain (the shared binary64 conversion constants are retained).
The temperature differences are T_y*(1-exp(-delta_x)) for positive heating and
-T_y*(exp(|delta_x|)-1) for net-beta cooling. Numerical reference uncertainty is
included in the temperature endpoints.

`tools/linux/stage_b/chemical_response.py:1` bounds I by interval rectangles:
each recorded scientific interval's entire true state is inside its left-node
error plus M_real*dt reach. Positive reaction monomials give lower/upper heating
powers; signed net-beta polynomial terms are bounded separately. Temperature and
Cstar extrema include every intersected cell boundary. Thus no smooth quadrature
assumption is made across a knot. The counterfactual's additional displacement
is checked below 1e-10 scaled units, far inside the stated domain margins.
The check uses exp(LH)*I_upper, which bounds the entire time-dependent tube,
not only the sharper endpoint displacement.
This is a conservative bound, not a leading-order energy/C estimate or a new
integration. The implementation retains the stated 110-digit arithmetic pads.

The primary-formula normalization is closed for both comparison flows. Let b
be the existing uniform anchored-versus-primary RHS perturbation bound. For
each full or term-suppressed flow, variation of constants gives
u_norm=b H exp(LH); the same triangle-majorant L bounds each anchored field.
The temperature response interval is expanded on each side by
2 Tmax (exp(u_norm)-1), using the declared domain's Tmax. This explicitly adds
both flow uncertainties without assuming their photon/direct offsets cancel.
It conservatively overlaps the normalization allowance already included in the
full trajectory's rectangle and temperature bounds.

Conditional sign resolution of this diagnostic does not make a microkelvin-scale
effect astrophysically significant, and does not discharge the full native Stage B
matrix's accuracy, event, reconstruction, or R20 obligations. The new comparison
argument and implementation have now received independent conditional review;
owner ratification under GOVERNANCE §5 remains pending.

A sharper bound uses the thermal damping sign. For both f0 definitions, the
net direct-minus-cooling power stays positive throughout the box, Cstar slopes
are positive, and cooling grows with x. Thus -Jxx <= partial f0_x/partial x <= 0.
Variation of constants gives a lower bound integral exp(-Jxx*(H-s))*g(s) ds
minus chemical feedback, and an upper bound I(H) plus chemical feedback.
The feedback magnitude is at most H*exp(LH)*I_upper*sum_j J_0j*s_j (j chemical),
using the already established full-state comparison bound. Weighted rectangle
integration is exact for the exponential weight. The report checks the original
10% signal-resolution criterion on the resulting interval's half-width and
midpoint. This tightening is an analytic sign/majorant argument, not an assumed
smallness or a measured response fit.

## Independent diagnostic values (documentation correction, 2026-11-12)

The complete channel table is in [REPORT.md](REPORT.md). The reviewer integrated
term-suppressed coupled three-ODE equations: full minus suppressed temperature
is +5.19 µK for LH, −2.59 µK for LH−DeltaLnu, −70.5 µK for the actual-potential
storage channel, and −73.1 µK combined. Only the first two have the analytic
intervals derived here; storage and total are independent convergence-supported
estimates. No new computation accompanies this citation. Source:
`review-evidence/REVIEW.md:281`, `review-evidence/outputs/trajectory-check.json:1`.
The diagnostic equations are not alternative energy-conserving physical models.
