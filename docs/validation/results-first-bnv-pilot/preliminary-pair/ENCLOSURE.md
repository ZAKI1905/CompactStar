# Conditional analytic enclosure for the retained one-cell thermal path

This is a new, unreviewed numerical-verification candidate. It changes no physical
equation, production step or criterion. It is conditional on the real-valued RHS
defined by the exported binary64 Cstar cell values, coefficients, and evaluated
photon/direct-power normalization anchors. An enclosure of that model is not
automatically an enclosure of unknown static normalization error. The separate
normalization calculation below now closes the evaluated anchors to the original
formula with fixed geometry operands. It does not requalify stellar parameters.

Scale the state by s=(1,1e-6,1e-6). On every real trajectory under consideration,
x is inside [-1e-4,2e-4], |eta_e|<3e-10 and |eta_mu|<8e-10. Verify this domain
including the error enclosure at every computed node. The whole real box is in
one Cstar cell. Analytically extend that cell's affine C(x) beyond its real cell
boundaries; this is solely for a complex-analytic derivative estimate, not a
change to the physical interpolation or to an integration that crosses a knot.

On a complex polydisc of radius R=0.1 in the scaled norm around any real-box
point, C is bounded away from zero by C(0)-|C'|(0.0002+R)>0. Bound exp(x) by
exp(0.0002+R), exp(-x) similarly, and each chemical polynomial by the sum of its
absolute coefficients at |eta|<=8e-10+R*1e-6. These give an explicit positive
bound M on the scaled vector field. `thermal_enclosure.py:constants` evaluates
these expressions with 110 digits and conservative padding. No observed
trajectory agreement enters M.

The analytic solution exists in the complex time disc r=R/(2M), staying within
R/2 of its initial point, by the bounded holomorphic Picard continuation/exit
argument. Its displacement is bounded by R/2 there. Every classical RK4 stage
also stays inside the state polydisc for |h|<=r: its stage displacement is at
most |h|M, and the final weighted displacement has the same bound. Both exact
and RK4 maps are analytic in h. Their coefficients agree through order four
(the classical RK4 order conditions). Cauchy's coefficient inequality therefore
gives the conservative one-step scaled error

    local_error <= R (h/r)^5 / (1-h/r),  0<h<r.

The real-box Jacobian majorant is computed from positive polynomial derivatives
and the affine Cstar denominator. Let L be its scaled infinity norm. Propagate
the local error with exp(L h). Include 1e-80 scaled units per high-precision step
for Decimal arithmetic: this is deliberately much larger than the 110-digit
operation-count error (fewer than 10,000 operations per step, bounded polynomial
arguments and denominators). Decimal elementary-function rounding and the stated
PI approximation remain implementation assumptions, cross-checked at 80/110
digits. Precision agreement alone is not the bound.

Cauchy's inequality also bounds the true thermal fourth time derivative by
24 R/r^4. For each Hermite segment of duration h, the exact-node cubic Hermite
remainder is bounded by (24 R/r^4) h^4/384. Recorded-node error adds at most the
largest endpoint x error plus h/4 times the largest endpoint derivative error.
The derivative error includes the independently evaluated real RHS and the
Jacobian-propagated RK4 enclosure. Thus the complete segment bound is

    max(endpoint x errors) + h/4 max(endpoint dx errors)
      + (24 R/r^4) h^4/384.

Check this uniform bound against the existing 1e-12 thermal-band assumption.
Separately compare all three reconstructed O2 state components and Echem against
the coupled high-precision endpoint and its analytic remainder, using the
original stored radii. This removes dependence on continued measured RK4
contraction for this local coupled check. It does not by itself certify the
production main/reference truncation, static anchor normalization, or event
coverage outside the actual one-cell trajectory. The following additional
calculations address the first two for this experiment.

## Normalization and global error

`normalization.cpp` reads R, M, surface redshift factors and the actual envelope
boundary redshift from the authenticated stellar archive. `normalization.py`
reevaluates the original envelope formula and -Bdot mu_B conversion at 110 digits,
using exact lifts of every frozen binary64 operand. The photon and direct-power
anchor discrepancies are respectively 1.66409e18 and 1.35463e18 erg/s, about
6e-16 and 4e-17 relative. They are numerical normalization discrepancies, not new
physical source definitions. A 1e-60 erg/s pad covers the stated 110-digit
arithmetic/96-digit PI implementation assumption.

Divide the sum of their uniform power discrepancies by Tmin Cmin to obtain an
additive scaled RHS perturbation b. Variation of constants gives b t exp(Lt)
as a conservative bound between the anchored and primary-formula solutions.
The one-year bounds in x are 7.306e-21 (source) and 4.028e-21 (control).
The same formula with the local interval length closes every local reconstruction
enclosure. No new trajectory is required. The anchor model's L suffices: write
the difference as its Lipschitz term plus the bounded forcing evaluated on the
perturbed solution; the real-box containment check closes the bootstrap.

`assess.py` encloses the original-t0 independent RK4 at every scientific sample
and retained accepted time using the analytic local remainder above, then adds
normalization perturbation. Error against this enclosed reference bounds total
main/reference error, including reference observation interpolation. The
predeclared twice-refinement differences and accumulated arithmetic term remain
in U_ref as well; this deliberately double-counts some error. No term is removed
to gain a PASS. `accepted.cpp` evaluates all five retained reference histories
at each main accepted endpoint using the original observation method. Its
observations receive the same independent reference and power assessment.

For control local reconstruction, compare output and retained left state to the
global enclosed solution. The triangle inequality plus exp(L span) times the
left-state error bounds the local solution starting from that retained state.
Zero chemical components use the existing owner-checked exact invariant. This
avoids unnecessary control witness integrations while checking every local box.

## Arithmetic residuals and refinement

`roundoff.py` lifts each observed RK stage, derivative and endpoint exactly,
parses the GSL rational RK tableau, and measures residuals of each stage equation
and the weighted endpoint equation against the independent real RHS. For stage j,
its state error is bounded by the incoming error, its stage-equation residual,
and h sum |a_jk| times (RHS evaluation defect + J times prior-stage error).
The endpoint uses the corresponding weighted sum. Retained-versus-replayed
endpoint differences are added. Every stage and endpoint enclosure is checked
inside the Jacobian box. This is an a posteriori arithmetic enclosure; RK
truncation is assessed separately with the global analytic reference.

This explicitly counts accumulation rather than assuming 64 ulps bounds a long
history. The 64-ulp minimum remains in the prescribed uncertainty. All observed
noncontracting refinement discrepancies fit within the independently evaluated
arithmetic allowances. Identical main cap histories only show inactive caps;
they are not evidence of cap convergence. Independent methods and the analytic
high-precision comparison supply the missing accuracy evidence.

`interpolation_roundoff.py` extends the same recurrence through every actual
reference observation: initialize from the retained reference's arithmetic
enclosure at the interpolation left endpoint, then audit its traced rk8pd stages.
Exact endpoints use the retained enclosure directly. The final reference budget
includes this observation arithmetic as well as the independently enclosed total
observation error. The 19,210-point audit required no new trajectory integration.

## Power and energy propagation

Positive complex-domain bounds for each power channel give gradient bounds by
Cauchy's inequality. A conservative 6 Pmax/R times scaled state uncertainty,
plus applicable normalization discrepancy, encloses each diagnostic power.
Constants and exact-zero channels are treated separately. The original 11-channel
budgets and noncancelling gross-power normalization are unchanged.

For a time-grid trapezoid, Cauchy gives |P''| <= 2 Pmax/r^2, hence the global
quadrature remainder H dt^2 |P''|/12. Thermal energy uses trapezoids in T instead:
on the occupied cell C(T)=a+b ln(T), so |C''(T)|<=|b|/Tmin^2 and the remainder
is at most sum |DeltaT|^3 |b|/(12 Tmin^2). Decimal reevaluation of the original
published-grid arithmetic bounds its binary64 accumulation error. Endpoint
thermal errors, chemical energy box propagation, integrated outgoing-power
errors and Cstar evaluation discrepancies contribute to U_main. The original
twice-nested-grid discrepancy remains in U_quad in addition to these bounds.

The independently established equilibrium-count rounding bound is added to the
source total without altering its recorded residual. No actual event occurs,
so U_event=0 has only this no-crossing scope. The original reconstruction and
combined reconstruction/event allocations are retained.

These are conditional mathematical enclosure arguments, not an assertion that
measured agreement proves bounds. The executable majorants, finite-precision
pads, fixed-input interpretation and derivation require independent domain review
under GOVERNANCE section 5. No outward-rounded interval package or proof assistant
has certified this implementation. Continuing observed RK4 contraction remains an
assumption of the older chemical-only diagnostic; the additional coupled analytic
enclosure avoids that assumption for the local state/energy and global tests.
