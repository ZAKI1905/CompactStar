# Piecewise verification candidate

Extends `../preliminary-pair/ENCLOSURE.md:1` and its accepted fixed-input scope.
The implementation is `tools/linux/stage_b/piecewise.py:1`. This is new verification
code, not a production integrator change or an already accepted Stage B enclosure.

The actual frozen Cstar table is lifted exactly from its binary64 values. At each
interior knot, both neighboring affine expressions equal the same stored capacity.
Thus Cstar and the vector field are continuous; their derivatives may jump.
The exact real knot uses ln(Tcell/(Tref*kB)); the native knot is rounded binary64.
That rounding and native smooth-extension residual strips must enter actual-path
comparison/event error, rather than assuming native and real knots are identical.

Use the explicit real domain x in [-0.03,0.2], |eta_e|<4e-7 MeV,
|eta_mu|<8e-7 MeV, scaled by s=(1,1e-5,1e-5). Every reference step must check
its entire enclosure/reach against this domain; exhaustion stops the calculation.
For each intersected cell compute fresh positive polynomial Jacobian majorants,
real speed bound M_real, complex speed bound M and r=R/(2M), R=0.1.
The cell's affine extension is used only for its complex derivative estimate.
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
