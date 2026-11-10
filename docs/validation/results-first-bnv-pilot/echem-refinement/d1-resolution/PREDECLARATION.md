# D1 convergence resolution — predeclared bounded calculation

2026-11-09. Scientific candidate; owner-authorized offline numerical-method test,
engineering, documentation and generated evidence. Starting commit
`871dbd77230f9e21b54ea0129473fe658236fd56`. Authority: owner's convergence-resolution
task, adjudication Section D, GOVERNANCE §§2–3/5, ADR-0017 §§5–7, INV-11f/INV-15.
No production source, coefficient, tolerance or thermal interpolant changes.
The original D1 FAIL receipt and its entire evidence directory remain intact.

## Diagnosis available before the new calculation

The original central refinement differences at points 1, 64 and 1440 are about
1e-48 to 5e-37 MeV, versus state radii about 2e-26 to 4e-23 MeV. They are many
orders above 110-digit arithmetic spacing. All three have nine Hermite segments:
eight nearly uniform segments and a tiny terminal segment. The original global
RK4 grids do not land exactly on the retained binary64 knots. The terminal
segments are 3.637978807091713e-12, 2.328306436538696e-10 and
3.725290298461914e-9 seconds respectively; only the final RK4 stage samples them.

Signed successive differences closely fit A/2^k+B/16^k (k=0,1,2). Independently,
extrapolating the preceding regular Hermite polynomial to the endpoint and
evaluating the endpoint RHS jump predicts A=-span*delta_f/(64*12), agreeing to
about 3e-16, 3e-11 and 4e-7 relatively at the three points. This supports an
unresolved endpoint-sampling contribution combined with fourth-order error.
Point 1440 has cancellation between those contributions, not Decimal roundoff.
The analysis script will preserve exact differences for both components and all
offsets, ratios to original radii, grid phases and this diagnostic model. The
model fit itself is not an error bound. Point 1919 also has a short terminal
segment and passes original C1; a short segment alone need not cause failure.

## Fixed reference method

One new retained-data calculation, on the original 13 points and offsets
`b=-1e-12,0,+1e-12`, at 110 Decimal digits. Preserve the exact source-bound RHS,
binary64 left states, all Hermite values/derivatives/knots, units (seconds, MeV,
x=ln(Tinf/1e8 K)) and local time origins. Reauthenticate every original input,
the original D1 receipt and implementation against the starting commit.
Recheck all local RHS rows at the original two-ulp limit and independently check
the double lifting using integer ratios. Keep rho <= 1e-13 and delta_x=1e-12.

Integrate each Hermite segment separately in normalized local coordinate s in
[0,1]: d eta/ds = segment_duration * rhs(x_segment(s)+b,eta). Use classical RK4
with exact rational coefficients and **K=8,16,32,64 steps per segment**. Restart
the grid, not the chemical state, at each knot. Every retained segment, however
short, receives K steps. Nominal labels N=64,128,256,512 mean 8K, not total
steps: actual totals are K times the recorded segment count. No step crosses a
knot. The central K=64 endpoint is eta_ref. This is a fresh integration of the
same local IVP, independent of the original global grid and production solver.

## Acceptance fixed before execution

1. Analytic tests pass before retained integrations: existing 11 tests plus
   segment-normalization, exact cubic-forcing integrals including a tiny terminal
   piece, fourth-order convergence and the sensitivity-bound algebra.
2. Authentication, exact initial-state conversion, RHS fidelity and thermal
   sample gates must pass. Preserve details before stopping on a failure.
3. **C1 unchanged:** for both components and all offsets, successive absolute
   differences contract by <=0.5 or the newer difference is below
   1e-60*abs(eta_K64). Handle zero differences explicitly. All 13 points must pass.
4. u_disc=2*abs(eta_K64(0)-eta_K32(0));
   u_offset=max_b abs(eta_K64(b)-eta_K64(0)), componentwise. Original C2/C3 use
   u_ref=u_disc+u_offset, unchanged recorded state radii and U_Echem, and four
   energy corners. No contribution suppressed; U<=0.20F is unchanged.
5. Additionally evaluate the conservative thermal sensitivity bound below and
   C2/C3 with u_supported_i=u_disc_i+max(u_offset_i,B_thermal). Both the original
   and the stronger inequalities must pass for a conditional qualification.
6. C4 remains report-only for O2 and six witnesses: error in MeV/ulp and ratio
   to enclosure_ctx+original u_ref. Ratio >1 invalidates that enclosure as a
   bound but is not an independent acceptance gate.
7. Any failed gate ends this authorized investigation after preserving its exact
   magnitude. No altered thresholds, additional levels, alternative method,
   repeated retained integration or point removal after results.

## Conservative sensitivity calculation and its scope

The source-bound rate is exactly a positive-coefficient polynomial
R_j=sum(k=0..3) a_jk eta_j^(2k+1) T^(6-2k), where
a_jk=L_j*c_k/(kBerg*11513*kB^(2k+1)*pi^(2k+2)),
c=(14680,7560,840,24). Use absolute L if needed for upper bounds.

A cubic Hermite segment is a cubic Bezier polynomial with control values
x_left, x_left+h*dx_left/3, x_right-h*dx_right/3, x_right.
An intentionally looser upper bound is
x_upper=max_node(x)+max_segment(h*max(abs(dx_left),abs(dx_right))/3).
Set T_upper=1e8*exp(x_upper+delta_x).

Let H be interval duration, y0=max_i abs(eta_left_i),
q=max_i sum_j abs(Z_ij)*abs(sigma_j) >= max_i abs(source_i),
Y=2*(y0+H*q). Evaluate positive polynomial bounds on |eta_i|<=Y,
T<=T_upper: beta=max_i sum_j |Z_ij| sum_k a_jk Y^(2k+1)T_upper^(6-2k).
Require y0+H*(q+beta)<Y. A first-exit argument then contains the central and
all thermally perturbed exact solutions in that domain. Set
L=max_i sum_j |Z_ij| sum_k (2k+1)a_jk Y^(2k)T_upper^(6-2k),
G=max_i sum_j |Z_ij| sum_k (6-2k)a_jk Y^(2k+1)T_upper^(6-2k).
The integral differential inequality gives
**B_thermal=delta_x*G*H*exp(L*H)** for either component, for any time-dependent
thermal perturbation bounded uniformly by delta_x. This is an algebraic bound,
not an inference that the constant-offset probes are extrema. Compute positive
upper bounds with Decimal ceiling arithmetic; denominators use floor arithmetic
and exp is advanced one Decimal representable value above its rounded result.
This bounds the specified real model with its source-bound pi constant.

The true thermal path being uniformly within the declared band remains an
assumption: node samples do not prove it. The discretization allowance remains
a convergence-supported estimate, conditional on continuation of the observed
contraction regime beyond the finest grid. For contraction <=0.5 thereafter,
the geometric remaining tail is <=the last difference, covered by twice that
difference. Four finite levels do not prove that continuation or constitute a
directed-rounding enclosure of the full IVP. Report this limitation explicitly.

## Execution and delivery

One EKU compute allocation, one logical execution CPU, 4 GiB, 30 minutes. All
large/new cluster artifacts stay under /mnt/sdd. Run tests and the fixed matrix
once; save results incrementally, exact failures, input hashes and execution
receipt. Commit this predeclaration and implementation before submission.

The original incremental bundle was authenticated (SHA-256
6cd9bc1a438df458f83d64335589889203babd781fc18b87ee44079955766553).
The office Mac/Scriptorium route was unavailable; the owner lacked laptop access.
The cluster's existing GitHub SSH identity successfully published the original
871dbd7 commit without force, and ls-remote confirmed the exact commit before
this new calculation. No original FAIL evidence was replaced.

Conditional success qualifies only these 13 local reconstructions under the
stated reference and thermal assumptions, pending domain-owner review. It does
not qualify every Stage A checkpoint, accumulated main-prefix error from t=0,
thermal RHS roundoff, R18/R20 or complete paired-run restart, chemical heating
in the one-year temperature observable, or Stage B readiness. Any subsequent
authorized Stage A run must diagnose its actual output checkpoints. The sole
possible future production change remains O2 chemical absolute tolerances
(1e-24,1e-24) -> (5e-27,5e-27) MeV; it is not applied here. No Stage A retry,
thermal retry, Stage B, master merge or certificate promotion.
