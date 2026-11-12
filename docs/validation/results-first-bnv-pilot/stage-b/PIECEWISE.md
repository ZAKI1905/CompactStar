# Piecewise verification candidate

Extends `../preliminary-pair/ENCLOSURE.md:1` and its accepted fixed-input scope.
The implementation is `tools/linux/stage_b/piecewise.py:1`. This is new verification
code, not a production integrator change or an already accepted Stage B enclosure.

The actual frozen Cstar table is lifted exactly from its binary64 values. At each
interior knot, both neighboring affine expressions equal the same stored capacity.
Thus Cstar and the vector field are continuous; their derivatives may jump.
The exact real knot uses ln(Tcell/(Tref*kB)); the native knot is the result of
the native binary64 calculation, not necessarily the correctly rounded real knot.
The reviewer inferred consistent native knots from retained events (the native
knot formula was not independently reproduced): offsets are +1.06e−16 in x at
knot 94 (about 15 ulp) and +1.8e−17 at knot 95, with timing shifts ≤5e−5 s.
These differences enter actual-path comparison/event error.

RKF45 event and event-trial c=1 stages in the two retained references overshoot
the knot by up to 1.3e−13 in x on the old cell's affine extension. Their RHS defect
is ≤82 ulp (~2.0e−26 s⁻¹), included in the actual-stage arithmetic audit. Main
histories have no retained stage traces: their wider strips are covered by the
measured state error against the independent reference (maximum about 0.3052 of
budget), not by that reference-stage strip measurement. Evidence:
`review-evidence/REVIEW.md:251`, `review-evidence/outputs/crossings-check.json:1`.

This event evidence covers only two transversal upward source crossings (94, 95).
No downward, near-tangential or control crossing is qualified. Endpoint bypass
and nextafter guards rest on authenticated guard code and successful task
receipts; individual neighbour values were not retained and the reviewer did not
rerun the guards (`review-evidence/REVIEW.md:495`).

Use the explicit real domain x in [-0.03,0.2], |eta_e|<4e-7 MeV,
|eta_mu|<8e-7 MeV, scaled by s=(1,1e-5,1e-5). Every reference step must check
its entire enclosure/reach against this domain; exhaustion stops the calculation.
For each intersected cell compute fresh positive polynomial Jacobian majorants,
real speed bound M_real, complex speed bound M and r=R/(2M), R=0.1.
In this HP enclosure, the cell's affine extension supplies the complex derivative estimate.
Cstar positivity is checked in that complex polydisc. No extrapolated crossing
rate is used to establish actual coverage. Runtime traces establish crossings.

Within a cell, the accepted Cauchy bound is R(h/r)^5/(1-h/r). Keep the entire
reference/true trajectory tube in that cell using distance > U+2 M_real h.
Use the maximum cellwise Jacobian norm L to propagate prefix error with exp(Lh).
Continuity of the piecewise field implies the same Lipschitz bound across knots.
The 110-digit/PI and conservative per-step 1e-80 arithmetic assumptions are retained.
Primary-formula photon/direct normalization remains a separately accumulated
forcing; the real RHS convention is unchanged from Stage A.

For a crossing bridge, do not use the analytic-cell remainder. For an autonomous
L-Lipschitz field bounded by M_real, the exact flow differs from its Euler step
by at most L M_real h²/2. Classical RK4 differs from the same Euler step by at
most L M_real h²/2, since its stage displacements are at most h M_real/2,
h M_real/2 and h M_real. Hence L M_real h² bounds their difference; the code
uses the larger L M_real h² exp(2Lh). All stage/reach tubes must stay inside
the real domain. A short bridge of at most 0.01 s confines this contribution.
The reference evaluates the actual piecewise field at all RK stages. Prefix
error is never reset at a crossing. This avoids requiring analytic smoothness
or an exactly represented event time in the high-precision reference.

Native event qualification remains the original contract: full event bracket,
residual, exact endpoint bypass, both nextafter neighbors, continuity and
transversality, and event uncertainty at the original thresholds. The offline
bridge bound does not automatically discharge those production obligations.
Continuous Hermite or quadrature bounds must likewise split at encountered knots;
a smooth fourth/second derivative bound cannot be applied over an unsplit knot.
A failure of a native checkpoint stops its affected task, even if an independent
reference can produce an approximate temperature.

The Stage B local check in `tools/linux/stage_b/local_reference.py:1` encloses
the coupled three-component solution from each exact retained native left state
and compares every scientific sample against its original stored radii. It does
not use the older auxiliary chemistry-only witness driven by a prescribed Hermite
temperature path (`tools/linux/preliminary_pair/high_precision.py:91`). Consequently
it does not assume, or claim to verify for Stage B, that witness's uniform 1e−12
thermal band. This is a different independent verification argument, not a change
to production RK8PD reconstruction, its component radii, or the original F_i
allocation. The original plan's state/power gates apply at scientific samples and
accepted endpoints; these checks do not certify arbitrary unsampled output curves.

## Actual-stage arithmetic at crossings

The one-cell Stage A trace parser cannot assume that a selected-cell block is a
single RK map at an event. The isolated audit class is generated from the exact
production header, changing only its class name and adding passive begin/end
callbacks around each GSL Apply. A source-equivalence test reverses those edits
and checks exact equality. Production binaries and their method are unchanged.

Each Apply record identifies its exact initial state, endpoint, RHS row interval
and rejected-attempt count. The residual audit selects the accepted final GSL
attempt and carries the bound from the recorded initial state's dependency.
Discarded event trials are not treated as successive physical evolution. Cached
first derivatives come from actual nonnegative-cell RHS records; a cell-selection
probe does not replace GSL's cached derivative. Every actual stage is checked
against the continuous piecewise real RHS and the global positive Jacobian.
Smooth-extension versus real-field defects in residual strips are included.
Unsupported trace shapes, missing dependencies or domain exhaustion stop the
audit, with the evidence retained. Known-answer tests cover ordinary maps,
discarded trials, rejected attempts, their combination and cached derivatives.

The reference budget retains twice the original refinement/method discrepancies,
the 64-ulp minimum, accumulated actual-stage/observation arithmetic, and the
independently enclosed total reference error. The last two deliberately overlap;
no positive term is dropped to obtain a favorable result. Restart comparisons
carry the full original-prefix error through the remaining horizon.
