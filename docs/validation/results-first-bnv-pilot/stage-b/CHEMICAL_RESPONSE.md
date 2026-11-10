# Isolating a chemical temperature contribution without a new trajectory

This comparison-bound candidate answers a narrower question than the paired
source/control separation. Define the contribution of LH by comparing the full
coupled model y with the same initial-value problem z whose thermal RHS omits
only LH/(T Cstar). Chemical kinetics and all other terms remain coupled. This is
a diagnostic mathematical definition, not a proposed alternative physical model.
No such counterfactual trajectory is run and no production equation is changed.
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
This is a conservative bound, not a leading-order energy/C estimate or a new
integration. The implementation retains the stated 110-digit arithmetic pads.

Conditional sign resolution of this diagnostic does not make a microkelvin-scale
effect astrophysically significant, and does not discharge the full native Stage B
matrix's accuracy, event, reconstruction, or R20 obligations. The new comparison
argument and implementation remain candidates for independent review under
GOVERNANCE §5.
