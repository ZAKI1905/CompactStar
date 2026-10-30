# Methods, evidence scope and stop assessment

This is a validation-only candidate. Governing norms are the unchanged preflight
§11 (`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:691`).
The exact prospective claims, tests and stop were committed in
`PREDECLARATION.md:1` at d5b40b3. Existing source equations and authenticated local
records supply the derivations; the Mac literature library is unavailable here.
No independent author review of this new package is claimed.

## What an error estimate must cover

At each achieved production B, write P-Q=(P-I_X)+(I_X-I_Q)+(I_Q-Q).
P is the historical production number, X the established GL4 result, I the
piecewise-linear interpolator, and Q the physical equilibrium response with the
declared finite-pressure boundary law. Production arithmetic, quadrature, EOS
table, surface selection and support pinning at the samples are already in P-I_X.
They are not charged again individually. The missing terms are reference-side,
input/normalization, and between-sample physical response. Reconstruction R is
only a diagnostic for the metric group until an independent nodewise reference
exists. (`tools/linux/ba13_evidence/assemble.py:39`, `reference.py:63`.)

Each inherited term retains its source hashes and evidence class. The numerical
EOS/centre *sensitivity* entries retain that limited scope; their missing remainder
is not erased. A complete sum is null if any term or dependency is unresolved or
unaccepted. Required accuracy regions are a different data type from estimates
(`tools/linux/ba13_final_decision/model.py:19`); they cannot qualify a row.

## Reference uncertainty inventory

| Contribution | Established evidence | Classification and precise remainder |
|---|---|---|
| GL4 integration | N=32..4096, three central targets, signed fourth-order contraction; independent manufactured and off-stage checks | Validated numerical estimate for tested primary quantities under an asymptotic model; no arbitrary accuracy demand or guaranteed global tail |
| Implicit stages | Newton residual monitored; N1024/2048 binary64/long-double ablations, manufactured order and mass/B/metric closures | Residuals and ablations are observations; output error additionally needs a stable stage inverse and propagated response, including passive quantities |
| EOS evaluation/inversion | Independent momentum quadrature, thermodynamic/charge/beta checks, EOS96/192 sensitivity, bisection | Validated method checks and finite sensitivities; need output-level residual/conditioning or an explicitly accepted quantitative transfer/remainder model |
| Centre series | Analytic u/v expansions; omitted functional order O(x0^(5/2)), x0 proportional to N^-2; factor-100 ablations | Analytical order and numerical sensitivity; omitted coefficients and propagation are not measured by a zero difference |
| Surface | Hs is a prescribed segment endpoint; positive Hs avoids vacuum singularity | No event-search error for the prescribed Hs; pressure-to-H inversion and the physical moving-cut law still need uncertainty |
| Floating point/summation | Long double, plain/compensated and binary64 comparisons | Finite sensitivities; positive-sum gamma bounds cover fixed terms only, excluding trajectory/EOS/libm input and differentiated cancellation |
| Common-B alignment | Measured residual plus retained B-functional tail, positive sampled B_H | Conditional (residual+u_B)/m with uniform m>0 and a complete u_B; residual against computed B alone is insufficient |
| Tangent differentiation | Twelve successful offset-H evaluations, three grids and three widths | Correlated discretization estimate; separate joint derivative/stencil/surface remainder below |
| Common reference bias | X1 regularized/unregularized comparisons and model-anchor residuals | Independent consistency checks; neither differences nor agreement bound shared omitted bias without an error model |

Source scope: `docs/validation/x2-fixed-step/REPORT.md:42`,
`MATHEMATICS.md:68`, `evidence/ablations.json:1`,
`docs/validation/ba13-evidence/evidence/diagnostics.json:1`;
implementation: `tests/bnv/x2_fixed_step.cpp:44`, `:129`, `:186`.
The retained binary64/Newton N512 failure remains a failure. The retired defect
transport is not reused. No 1e-11 floor is adopted from the independent review's
explicitly illustrative ledger. A common bias constant in B may cancel a
particular difference; only a demonstrated constant/additive or multiplicative
identity can establish that, and normalization must still be treated correctly.

The missing quantity is a joint output-level residual estimate over this small
window, not an entirely new stellar solver. A validated sensitivity/conditioning
argument under named assumptions is allowed by the proposal. The retained
record does not yet state that quantitative argument for every required output.

## Tangents and correlation

For A=N+−N−, D=B+−B− and t=A/D, use the exact perturbation identity

    t_true-t = (a-t*d)/(D+d).

The refinement identity uses the coarser t and finer D. Our independent 60-digit
recomputation uses the twelve retained raw JSON outputs and N512 endpoints,
not merely the stored tangent report. All nine ratios reproduce the original
values to 1e-55 absolute, pass signed order [3.5,4.5] and width-ratio [3,5].
The retained Richardson tail |delta_grid|/15 is a discretization estimate under
its asymptotic model. Tests exercise nonzero correlation and exact cancellation
separately; they prohibit cancelling arbitrary joint errors.

A complete derivative uncertainty needs J>=|a-t*d| and uD>=|d|<|D|, giving
J/(|D|-uD), plus the stencil remainder and alignment. These are absolute count/count
units, never relative errors divided by t_mu. Three width ratios support a leading
h² term; a remainder after that model, including common error's derivative, is
not fixed by the grid tail. The central-secant checks are a different geometry
on the same method, not an independent solver derivative. Their 0–20 residuals
are about -1.20e-13, +1.23e-13 and +1.73e-16 (n,e,mu). Exceeding a finest grid
tail is not a physical-threshold failure; neither is it evidence for a zero tail.
(`docs/validation/ba13-evidence/evidence/diagnostics.json:1`.)

The minimal missing tangent component is thus J/uD plus a quantitative width and
surface-path remainder, with either independent derivative evidence OR a validated
remainder model. Requiring a second solver unconditionally would exceed the proposed
standard. No independent tangent solver or new stencil campaign is launched.

## Moving surface, continuous window and target transfer

Put y=ln P_s and regard Q and B as functions of Hc and y. For the declared cut law
y(x), the total equilibrium derivative obeys

    H'=(B0-B_y*y')/B_H,
    H''=-(B_HH*H'^2+2 B_Hy*H'*y'+B_yy*y'^2+B_y*y'')/B_H,
    Q''=Q_HH*H'^2+2 Q_Hy*H'*y'+Q_yy*y'^2+Q_H*H''+Q_y*y''.

At fixed cut, t_fixed=N_H/B_H. Along the moving cut,

    t_path=t_fixed+(N_y-t_fixed*B_y)*y'/B0.

These are chain-rule identities, not estimates of their inputs. They identify
precisely the cut derivative and mixed-response terms that fixed-cut H stencils
alone do not measure. Across-target curvature incorporates the actual cuts at
three locations; its agreement with fixed-cut curvature supports a hypothesis,
not an unsampled supremum. The analytic support boundaries Hn/Hmu and positive
sampled B_H support a smooth fixed-topology physical branch. Chemical response
also requires invertibility of its susceptibility; moving endpoints contribute
through these total derivatives. Discrete neutron-edge pinning and envelope-node
snapping belong to sampled representation error, not a physical discontinuity.

The established linear-interpolation theorem gives K h²/8 if |Q_xx|<=K. All 20
gaps are explicitly reevaluated, with their retained coordinate uncertainties,
for every reference-covered group. These remain conditional. The model tests a
C2 compact bump between data sites to show why finite observations cannot alone
prove a supremum; this is a logical falsifier, not a claim that the EOS has such
a bump. A numerical continuation estimate can be accepted with a named error
model and independently validated regime; an analytic supremum is not demanded.

The successful target-5 and target-15 predictions and criteria are unchanged.
They test 15 non-tangent groups using each target's own pressure and achieved B.
They do not validate tangent interpolation, arbitrary common bias or all interior
points. No further holdout is run: without a missing-error model to test, another
star has no declared falsification rule for the actual gap. The smallest later
static experiment must be chosen *after* that model is stated, with a numeric
prediction/uncertainty and stop rule committed before execution. There is no
blanket 21-star high-resolution plan.

The final achieved x is -9.999999990055996e-7: the nominal endpoint falls about
9.94e-16 beyond the observed range. This is far below the inherited coordinate
uncertainty scale and is part of C-INPUT/C-CONT, not a demand for tighter shooting.
An accepted derivative bound can cover the endpoint displacement; interpolation
must not silently extrapolate. Existing selected-grid margins are evidence of
sampled topology only; no adaptive-controller continuity claim is needed under
interpretation (2). Source: `docs/validation/ba13-closure/evidence/m3.json:1`,
`m4.json:1`, `docs/validation/ba13-evidence/evidence/continuation.json:1`.

## Chemical and thermal groups

Both Z rows use direct X2 full-inverse comparisons. Independently replaying the
review's 50-digit M5 calculation verifies, for A=G^-1, K=|A|E and rho(K)<1,

    M=(I-K)^-1 K|A|,
    paired Z majorant = |L|^T (M0+Mi) |L|.

Each endpoint and both Neumann remainders are retained. The alternative estimates
are approximately 1.9647e-7 (npe) and 5.5967e-7 (npmu), conditional on their
inherited background error sets. They are not added to P−X and do not replace
representation/alignment. This replay checks algebra and source inputs; it is
not new independent review or validation of the background assumption.
(`docs/validation/ba13-closure/METHOD.md:121`; new
`evidence/independent-reference-replay.json:1`.)

All N_n/e/mu, t_n/e/mu, Z rows and both modified-Urca coefficients retain individual
budgets. The large muon relative sensitivity gets no proxy or relaxed threshold.
Cstar, envelope and Tsurface_inf retain all 13 thermal knots. In this controlled
linear-temperature source C_cached(B,T)=F(T)C_linear(B,T), so F(T) cancels exactly
in same-knot relative depletion ratios. The absolute cache bias, up to 6.55536e-4,
remains a separate limitation on absolute thermal accuracy; it cannot be used to
claim failure of the relative frozen-background requirement or cancellation of
structure/EOS/roundoff. (`docs/validation/ba13-uncertainty/SCIENTIFIC_REVIEW.md:90`.)

The equilibrium redshifted field satisfies mu_n_inf=mu_B_inf at eta=0, with
H=log(mu_B/h0), nu=K-H in X2. This supports all 2049 mu_n_profile nodes under
the actual equilibrium/reference definition only. The actual eta-dependent
potential remains separate. mu_B, P0 and P1 retain their own source errors; P2
is the exact zero identity. Radius, gravity, envelope and surface temperature
retain the moving boundary and normalization dependencies.

Full metric_structure still requires r,m,nu,lambda at all original support
coordinates and the moving surface endpoint. No observer is implemented or run:
the complete reference/continuation model is missing, independently of owner
signature. If those prerequisites become sufficiently supported, the smallest
passive observer must leave stages/integration untouched, check interior
interpolation independently at withheld coordinates, record source-bound errors,
and retain the original support including the terminal surface. Quadratic
collocation interior values do not inherit fourth-order endpoint convergence.
Source: `docs/validation/ba13-evidence/REPORT.md:107` and X2 mathematics above.

## Headroom inversion, not invented uncertainty

Let F be the inherited finite added ledger including its normalization, D the
historical drift, T its unchanged threshold, b the baseline error estimate and
f=||I-H|| over all 21 fit rows. For additional uniform pointwise error e, extra
continuous gap remainder g and fit/input remainder z, a sufficient feasible region is

    (2+f+T)e + g + z <= (T-D-F)(1-b), b+e<1.

For absolute tangents the T and (1-b) normalization factors are absent.
The paired term costs 2e and the fixed-weight residual costs f e. To derive the
relative formula, write the old finite numerator as (D+F)(1-b), add
(2+f)e+g+z, divide by 1-b-e and require <=T. Nonlinear changes to the design or
weights belong in z; a constant/linear error is annihilated only for the fixed
projection. Shared B0 input is propagated through the exact coordinate identity
in `docs/validation/ba13-closure/METHOD.md:103`. No covariance is repurposed as
physical error. A cap is an axis intercept with other additional errors zero,
not a simultaneous budget allocation. Accepted baseline estimates and positive
lepton/support conditions remain prerequisites. These are feasibility results,
not a completed uncertainty construction or proof of unattainable accuracy.

The 12 inherited category shares are retained for traceability, explicitly as
planning diagnostics. No historical FAIL is rescored and no new substantive
scientific requirement is inferred from those shares.
