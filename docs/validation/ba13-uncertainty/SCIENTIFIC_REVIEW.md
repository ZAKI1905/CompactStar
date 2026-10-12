# BA13 static numerical uncertainty investigation

**Disposition B — specific additional bounded evidence required. No certificate
is qualified or promoted.** Option A was authorized as an investigation; the
owner has not ratified new error estimators, envelope semantics or a successor
method. The original fitter and certificates remain byte-identical. This study
establishes useful partial evidence and identifies missing errors; it does not
supply finite qualified bounds for all 21 targets.

## Identity and scope

Entry: CompactStar `8e6ac589d8d34e78447ee893e43e5673490eb728` on the recovered
Linux candidate. Separate branch `candidate/linux-ba13-uncertainty-20261012`,
worktree `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty`.
The predeclaration was committed as `b48c642` before any numerical calculation.
The first withheld predictions were committed in `1733092`. Original sources,
production scientific owners, dependencies and nominal 80000/8192 inputs were
not changed. Validation-only executables expose existing solver controls.

The historical archive SHA-256 is
`990d3708a24c71e518d4bb59e6ad3c23ae47f337ba924007b083ca598e6b8e2b`;
all 351 extracted manifest entries were verified again. The recovered fitter is
`9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370`.
Native run1/repeat raw and final manifests, exact profile/coefficient identities,
dependency lock and both retained 10,303-component 80/112-digit refits were
verified. Refit equality establishes arithmetic reproduction, not input accuracy.
See `authentication.json`, `evidence/reauthentication.json` and the retained
recovery records for exact source heads and lock bindings.

## Principal numerical findings

| Evidence | Finding | Disposition |
|---|---|---|
| Radial scalar/particle-number ladder | 64 eligible predictions pass withheld 320000 tests, including targets 5 and 15 | Conditional axis evidence only |
| EOS ladder | 32768 refuses at targets 0,10,20; four pairs of near-neutron-onset nodes round to equal pressure and/or energy | EOS holdout INCONCLUSIVE |
| Thermal radial ladder | Last increment grows by factors 1.236, 1.648, 2.730 at targets 0,10,20 | No geometric error bound |
| Independent thermal source/cache evaluation | Reconstructs all 546 native/Mac knot values within 6.52e-15 relative; common cache bias reaches 6.55536e-4 | Missing absolute interpolation error identified |
| Nominal thermal radial quadrature | Trapezoidal result differs from onset-split analytic-density integral by 7.12–7.16e-7 relative | Missing integration/representation error identified |
| Fixed-profile Ltilde_Me quadrature | Order/subdivision spread about 1.51e-16; nominal-to-320000 change 6.346–6.354e-5 | Common profile/support error is not bounded by order agreement |
| Direct P0/P1 averages | Fixed-profile independent agreement near one ulp; radial changes about 1.54e-9 / 2.14e-9 | Background error remains missing |
| Chemical response quadrature | Orders 16,32,64,96 and split32 Z spreads up to 1.38e-14 relative | Fixed-background arithmetic/quadrature diagnostic only |
| Fit units | Back-transformed covariance changes by many orders; Z row covariance ratios reach roughly 9.09e-94 | Historical floor is not physical uncertainty |

### Structure, common B, center, surface and support

There are 55 completed common-target structure/thermal cases and three refused
EOS cases. Every successful case uses the original absolute target B, including
B0; no refined run redefines B0. Maximum achieved residual is about 2.83e-14 of
target B; maximum bracket is 4.49e-11. Tightening tolerances to 5e-12 and 5e-13
retains identical selected stars in the three pilots. Endpoint R/M spreads and
central-density shifts are retained. The nominal-settings re-solve at target 15
is a distinct profile: achieved B changes by -2.74481e-15 relative and rho_c by
-18.125 g/cm^3 from the authenticated nominal certificate input. The absolute
B target is identical. Its partial holdout estimates apply to that study profile;
transfer back to the certificate input still requires an alignment enclosure.
Targets 0,5,10,20 reproduce original B and Cstar exactly. No certificate is
attached to the distinct target-15 profile. Neither small residuals nor endpoint spreads
prove a global target-sensitivity enclosure.

Common q=r/R mapping uses 2049 nodes, with 4097-node interpolation diagnostics and
withheld raw nodes. The report retains r,m,nu,lambda,nB,nn,np,ne,nmu and redshifted
neutron potential, linear/PCHIP spread, node counts, support edges and topology.
All tested topologies remain (3,2,1). These are full-radius diagnostics; they do
not replace the historical neutron-support coordinate or its exact certificate
fields. Mapping/interpolation spreads alone are not bounds. R uncertainty and
onset-location uncertainty remain part of a future joint enclosure.

The neutron-support edge moves by about 2.59e-5 in normalized radius between
nominal and finest radial cases. Ltilde_Me's last radial increment grows from
roughly 2.0e-5 to 4.3e-5 relative. This is direct evidence that a nearly exact
quadrature of a reconstructed profile can still have substantial profile error.
Ltilde_Mmu changes about 1.78e-7. Signed paired refinements show strong common
errors; they are retained in `paired-integrals.json`, not treated as independent
noise or assumed to cancel without a bound.

Center-radius changes (2,1,0.5,0.25 cm) produce maximum observed changes of about
1.79e-12 in R, 1.26e-13 in M and 4.45e-12 in N_mu; the small, irregular changes
are unresolved as global center bounds. Lower-density extensions to 1e-15/1e-16
retain identical selected R/M/N values. The existing central-pressure cutoff
still limits this test; equality does not prove convergence to vacuum. Existing
actual-cut analytic tail evidence is retained, while independent per-target
center/background remainders remain missing.

The EOS 32768 diagnostic finds equal-value rounding pairs around
7.356728903732e-9 fm^-3. The existing strict-ordering guard refuses safely. No
nodes were removed, sorted by value, smoothed or replaced to force a pass. This
is a refinement-tool precision limit requiring a separately reviewed response,
not evidence that nominal 8192 tables are unordered.

### Thermal errors and the 214/273 platform exceedances

The original adapter's entropy law is linear in temperature. Independently,

`Cstar/Tinf = kB^2 conversion * integral sum_s[pF_s sqrt(m_s^2+pF_s^2)]/(3 hc^3) exp(-nu) dV`.

The unchanged thermal owner caches C at 160 logarithmically spaced temperatures
and interpolates C linearly in log(T). For C proportional to T this gives the
explicit factor `F(T)=((1-w)Tlo+w Thi)/T`, with maximum F-1=6.55535864063e-4.
Independent nodal trapezoidal integration times F reproduces the producer to
6.51416e-15 relative. Twenty-seven 80/112-digit held-node Cv checks and 27
independent momentum-integral P1 checks pass fixed-input arithmetic tests.
These checks do not validate the star itself.

The 214/273 exceedances are reproduced. Maximum Mac/Linux Cstar difference is
2.36795e-13; independently summed nodal integrals differ by up to 2.34432e-13.
Thus the discrepancy is already present in nodal structure/composition inputs,
with small additional evaluation/accumulation differences. The original
64-epsilon assignment (about 1.42e-14 relative) omits both this input dependence
and much larger refinement/interpolation effects. A positive-sum gamma_N bound
on fixed nodal terms is about 4.50e-12, but it excludes input/libm errors and is
not adopted as a total or platform bound. No measured platform difference was
used to construct an acceptance tolerance.

Onset-split integration separately evaluates analytic composition and linearly
interpolated nodal composition. Nominal trapezoid versus independent integral
spread is about 7.1e-7. The analytic/nodal-composition distinction is about 4.5e-9
at nominal resolution. Refined independent integrals still change about 7.9e-7
from nominal to finest. This separates integration, composition reconstruction
and changing background contributions; none has yet received a global remainder
proof. Apparent contraction of the last two independent-integral increments is
insufficient for the predeclared three-increment plus holdout test.

F(T) is exactly common across targets for this linear source and fixed T, so it
cancels in ideal relative depletion ratios. Treating the absolute cache bias as
an independent error at every target would double count a known shared effect.
This cancellation does not apply automatically to structure, EOS, nodal
interpolation or rounding errors. Envelope and redshifted-temperature pilot
changes and all 13 knots are retained; their radial sequences also fail the
predeclared eligibility test. No surface-temperature accuracy claim follows.

### Chemical and event errors

All 21 nominal profiles receive a fresh leading regular-center estimate. Its
ratio to the inherited center matrix is 0.9999972913 to 1.0000000000000002 (the
upper excess is binary rounding). This checks the leading estimate only: it does
not establish the full center remainder or the inherited B0 background bound at
perturbed stars. Original center/background/tail matrices remain in the evidence.
Fresh response quadrature at targets 0,10,20 uses independent analytic response
mathematics, orders through 96 and cell bisection; Z spread is near arithmetic
resolution and is not a bound on shared profile error. Per-target new background
and tangent-stencil refinements remain missing, explicitly including t_n,t_e,t_mu
and both Z rows. Existing Phase-5B/C tail and tangent evidence is preserved.

P0=mu_B is checked by direct neutron-weighted integration; P1 uses an independent
momentum quadrature, retaining the historical P1 constants. No replacement of
those constants is smuggled into this comparison. P2 remains analytically exactly
zero at all targets. Neutron support and event measure are recorded. Exact-zero
P2 passes its own check; that cannot qualify the remaining event/profile errors.

### Fit conditioning and scientific uncertainty

All original 10,303 refits are unchanged. A separate diagnostic runs the unchanged
fitter after unit rescaling by inverse nominal scale times 2^-40,1,2^40, then
transforms coefficients and covariance back. Tiny dimensional Z/Ltilde quantities
are dominated by `64*epsilon*max(1,abs(value))` in their original units. Covariance
is strongly unit dependent; for Ltilde_Me the observed ratios range roughly
3.29e-64 to 3.97e-40. Fitted drift changes over the window remain small (largest
reported group maximum about 8.45e-14), but this does not rehabilitate covariance
as measurement error. No physical uncertainty is inferred from that floor, and
high-precision input arithmetic is not evidence of input accuracy.

## Frozen-envelope interpretation for owner adjudication

The historical implementation computes a cumulative maximum **observed** drift
plus the largest accepted outer-half residual. It records normalized numerical
uncertainty separately, uses it in covariance/residual testing and some right-hand
thresholds, and adds no separate input-error term to `drift_bound`. Historical
reference errors are combined by hypot; shared deterministic systematics do not
justify that independence assumption.

Preflight section 11 admits at least these readings: (1) bound physical drift,
reference/normalization uncertainty and nonlinear remainder on the left against a
physical limit; (2) retain the explicitly error-aware tabulated right threshold,
but still include the full appropriate numerical contribution on the left. A
third historical operational reading treats the right-hand 5u allowance as the
entire numerical accommodation. The phrase “full drift plus numerical and
nonlinearity terms” does not unambiguously establish that third reading. Reading
(1) also requires owner treatment of the explicit max/plus-5u threshold formulas.
We recommend explicit adjudication of (2) and its deterministic joint error
model, without relaxing particle-number limits or choosing whichever passes.

Proposed independent check: represent each shared structure/EOS/background error
once as a bounded variable theta, with point errors e_i(theta), achieved-x errors,
and uncertain reference scale S(theta)>0. Evaluate a supremum of the cumulative
normalized drift plus residual remainder over that same joint set. For a fixed
fit design/weights H, residual perturbations are `(I-H)e`, not independent fitted
noise. If x or scientific weights are uncertain, enclose H as well. A conservative
fallback bounds |yi-y0|+sup|ei-e0| over S0-U_S and bounds the residual linear forms
separately, then sums; report its conservatism. Do not add covariance a second time
after the same input errors have been propagated. Separately preserve any genuine
stochastic covariance model, if one is independently justified. At i=0, shared
e_i-e_0 cancels exactly. Keep nonlinear fit remainder distinct from source error.

Neither interpretation can presently be evaluated with finite complete bounds.
No historical `drift_bound`, `threshold` or `utilization` field was rewritten.
The matrix retains original utilization and leaves proposed utilization null
where required errors are missing. No threshold FAIL is inferred from unbounded
uncertainty, and no PASS is inferred from an observed small paired difference.

## Disposition matrix, cost and continuation

The machine-readable 21x22 matrix has 462 entries: 21 exact-P2 PASS, 121
INCONCLUSIVE and 320 UNBOUNDED. PASS is restricted to that exact check. The
remaining entries explicitly name missing bounds, original reported errors,
partial independent evidence, correlation treatment and utilization consequences.
“UNBOUNDED” means no established bound in this investigation, not infinite
physical error. No full-BA13 statement is made.

Resources and actual job identities are in `evidence/resources.json` and external
receipts. Measured numerical work is approximately 0.24 CPU-hours, plus about
0.054 build CPU-hours; final exact totals are in `evidence/summary.json`. Individual
star cases took at most about 22 seconds; the complete independent-integral pilot
took 242 seconds. Maximum measured RSS including build is 353916 KiB. Slurm
accounting storage is disabled; per-process resource measurements and scheduler
records substitute. Failed configure/submission attempts are documented in the
recipe and excluded from successful numerical CPU totals where unmeasured.

**Next step: independent scientific review, then targeted bounded refinement.**
A successor fitter is not yet proved necessary; accepting historical covariance
as physical error is unsupported. A continuation should first resolve the EOS
held-level precision limit and obtain three resolved increments plus a withheld
level for neutron-support reconstruction, Ltilde_Me and onset-aware Cv. Add
per-target response-background and tangent-step evidence only after these
estimators work. Then expand to all 21 targets with a joint error budget.

Cost estimate from this pilot: another three-target structure/onset block and
chemical-background checks should fit within 0.5–1 CPU-hour; 21-target expansion
at comparable costs is approximately 2–4 CPU-hours, with a new cap and stop rule
before higher resolutions. These are planning estimates, not authorization to
change a scientific owner or a claim that extra resolution will yield a bound.
No full BNV trajectory, long thermal regression, 500 kyr campaign, certificate
promotion, governed baseline change or nonlinear stability qualification occurred.
