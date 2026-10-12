# BA13 uncertainty investigation predeclaration

Status: proposed Option-A investigation, explicitly authorized 2026-10-12;
not a ratified numerical method or certificate acceptance. Commit this document
before executing any new refinement. Entry is 8e6ac589d8d34e78447ee893e43e5673490eb728.
Worktree: /mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-uncertainty.
Branch: candidate/linux-ba13-uncertainty-20261012. Original worktree, historical
fitter and historical/native outputs are immutable inputs. Governing records:
GOVERNANCE §§2–5, preflight §11, INV-13/14/15, ADR-0015/0016/0017 and proposed
ADR-0020/0021. Changes are validation-tool numerical-method investigation,
engineering, build registration, documentation and generated evidence only.

## Scope and independent axes

Study all 22 historical quantity groups. Start physical refinements at targets
0, 10 and 20. Targets 5 and 15 are held out; do not tune estimators to them.
No full-BA13 statement until all 21 targets and all groups have bounded errors.
All 21 receive an explicit disposition matrix even when their bounds are missing.
Nominal 80000 radial / EOS 8192 data are never replaced.

Separate axes (one varied at a time unless explicitly a cross check):

* Radial sampling 20000, 40000, 80000, 160000; 320000 is the withheld level.
* EOS parameter 2048, 4096, 8192, 16384; 32768 is withheld. The existing
  onset-augmented generator and governed analytic Track-R provider remain owners.
* Center radius 2, 1, 0.5 cm; 0.25 cm held out through validation-only subclass
  access to the existing protected initial-radius parameter. No TOV RHS change.
* Surface: existing finite-pressure cut plus the independent Phase-5B/C pe-tail
  continuation/comparison inequalities. EOS lower-density extensions 1e-15 and
  1e-16 (nominal 1e-14) are diagnostics, not substitutes for a tail bound; preserve
  existing grid points when extending. If the central-pressure cutoff dominates,
  record this insensitivity rather than claiming convergence to vacuum.
* Onsets/interpolation: existing analytic provider versus tabulated/nodal
  reconstruction, with thresholds split at governed onsets. Never smooth a
  threshold or extrapolate through unavailable active charts.
* Achieved B: common absolute nominal B_target for every resolution, including
  target 0 (do not redefine B0 on each resolution). Pilot tolerance factors
  5e-11, 5e-12, 5e-13 times nominal B0; failure to close a bracket is retained.
* Integral order 16, 32, 64; held-out 96, plus cell subdivision and independent
  Simpson/adaptive integration where available. Equality of 16/32 is not a bound.
* Thermal: all 13 original Tinf knots; independently evaluate the analytic
  degenerate free-gas entropy/Cv law with existing governed masses and units,
  compare native thermal-source interpolation and radial quadrature, and
  propagate structural/EOS variations separately. Do not infer a bound from the
  214 Mac/Linux exceedances. The same-source exact law is an interpolation oracle,
  not a claim of realistic finite-temperature EOS physics.
* Tangents: keep production finite-difference recipe unchanged; validation can
  halve its central-density step ladder to distinguish truncation from shared
  structure/EOS error. Z: integrate the current analytic response on each target,
  characterize center/tail/background per target, compare with inherited B0
  components, and propagate through the matrix inverse with conditioning retained.

## Common target and coordinates

For each configuration bracket the same absolute B_target using canonical solves.
Record B_solved, bracket, residual, rho_c, R, M and all number errors. No comparison
at equal rho_c is represented as equal B. Residual sensitivity is measured from
bracket endpoint or separate nearby-target solves. Without a bounded local
sensitivity, B-alignment uncertainty remains missing even if the residual is small.
Keep direct common-rho diagnostics distinct if they are useful for attribution.

Profiles use normalized radius q=r/R over [0,1], at 2049 fixed knots including
endpoints, in addition to raw native nodes and the historical neutron-support
coordinate. Map center regularity explicitly and split at species thresholds.
Linear and shape-preserving cubic reconstructions provide interpolation spread;
that spread alone is not an error bound. Withheld raw nodes and doubled common
coordinate density test it. R uncertainty propagates into physical-coordinate
comparisons; extrema near onsets use absolute density scales, not small local
relative denominators. Support signatures, analytic P2=0 and prohibited channels
retain exact checks. A chart ambiguity is INCONCLUSIVE, not an interpolated pass.

## Error construction and convergence conditions

Separate truncation/quadrature, EOS/interpolation, center/tail, target alignment,
arithmetic and model/representation errors. Solver local rtol and repeat identity
are never global accuracy estimates. Differences below an independently estimated
arithmetic floor are unresolved, not zero scientific uncertainty.

For each scalar/component and axis, use at least three successive differences.
A provisional geometric-tail estimate is eligible only when differences contract
with both last ratios <=0.75, inferred orders are positive and agree within a
factor of two, topology is unchanged, and differences resolve arithmetic noise.
Let q be the larger last ratio. For nominal-to-finest difference D and finest
increment d, the conditional estimator is |D|+2|d|/(1-q), plus separately bounded
alignment/interpolation/arithmetic terms. It must predict the withheld refinement
without being enlarged afterwards. These fixed criteria were not derived from
Mac/Linux differences. They are empirical eligibility tests, not a theorem about
uncomputed levels. Without an analytic remainder/comparison justification or
independently validated estimator, the result remains INCONCLUSIVE.

Irregular/noncontracting sequences, unresolved subtraction, rank/conditioning
failure, onset ambiguity or a failed holdout are INCONCLUSIVE; a missing finite
bound is UNBOUNDED. An independently bounded excess over a governed threshold is
FAIL. PASS is restricted to the explicitly established bound/check, not inferred
from a small observed difference. A material owner defect is reported separately.

Independent error axes combine by a deterministic sum unless a demonstrated
correlation model supports a tighter enclosure. Never RSS unknown systematic
errors. Shared EOS/background errors across stars remain correlated. Measure
paired target-minus-B0 changes with the same refinement and retain signs; common
error cancellation is allowed only when justified analytically or by a validated
joint enclosure. No independent-noise interpretation of solver error is assumed.
Absolute scales: tangent as in historical fitter; P2 absolute 1; Z row maximum of
nominal B0 entries; other nonzero components use their nominal B0 magnitude.
Reference/normalization uncertainty is explicit: require S0-U_S>0, and use that
lower denominator when bounding a ratio. Near-zero components use a declared
physical absolute scale or remain UNBOUNDED; never choose a scale from observed
platform drift.

## Fit and envelope diagnostics

Preserve the 10303-component 80/112-digit refit and every original field. Diagnose
unit dependence by scaling each series by inverse nominal magnitude and powers
of two on either side of that scale, transforming outputs back; perform fits
with unchanged historical code on disposable inputs. Its 64*epsilon*max(1,|X|)
floor is conditioning protection, not an error source or platform allowance.
Compare fits using scientific input errors separately in high precision; unknown
correlations preclude a probabilistic covariance interpretation.

For independent qualification distinguish observed drift, error of each value,
B0/reference error, normalization error, fit remainder, monotone envelope and
threshold. Study deterministic interval/joint-error propagation: each underlying
systematic-error variable appears once; evaluate drift and fit-residual linear
forms with that shared variable rather than adding fitted covariance as another
independent error. A simple triangle bound is an explicitly conservative fallback,
not evidence of independence. Retain alternatives: error-expanded left envelope
against a fixed physical limit, versus the historical error-aware right threshold.
Report both where bounds exist; never choose the passing interpretation. Recommend
owner adjudication before any successor field or acceptance claim.

## Staging, validation and cost limits

Stage 1 is the highest-information three-target structure/EOS/thermal pilot plus
all-target serialized fit-conditioning/error-inheritance diagnostics. Reuse exact
nominal profiles and proven prior characterization only with hash bindings.
Stage 2 expands integral/chemical/tangent and holdout evidence only if pilot
integrity and convergence support it. Expansion to 21 refined targets requires
successful estimator holdouts and a cost estimate; it is not automatic.

Initial work cap: 4 measured CPU-hours of new pilot numerical processes; up to
8 single-core jobs, <=8 GiB each. Stop expansion if a representative refinement
exceeds 5 minutes or 2 GiB, or the projected next block exceeds that cap. Complete
running highest-information pilots, report actual costs and a costed continuation.
Build time is separately reported. There is no trajectory, long regression,
500 kyr run or dependency modernization. Every stage uses fresh directories,
retained failure logs, hash/source bindings and atomic completion publication.

Deliver measurements and a 21x22 matrix with original errors, proposed errors
(or null), correlation treatment, utilization consequence (or null), and explicit
PASS/FAIL/INCONCLUSIVE/UNBOUNDED. Final disposition A/B/C/D follows the user's
request. Missing evidence is never silently converted to a passing bound.
