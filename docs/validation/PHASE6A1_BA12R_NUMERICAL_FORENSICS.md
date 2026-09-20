# Phase-6A-1 BA12R numerical forensics

Date: 2026-09-20

Scope: analysis and documentation only, using existing source, authenticated
BASELINE/REFINED/ULTRA trajectory bytes, existing step logs, and existing
validation evidence.

## Disposition

**BA12R ROOT CAUSE NARROWED BUT NOT UNIQUE — MINIMAL DISCRIMINATING
EXPERIMENT READY FOR OWNER REVIEW.**

This is disposition **B**.

The present hierarchy is not an asymptotic tolerance-refinement family.  The
dominant architecture advances to every output checkpoint as a hard GSL `t1`
ceiling.  Consequently, 99.65--99.89% of the 8192 output intervals contain
exactly one accepted step, usually the entire fixed output spacing.  Tightening
the tolerance by 100 between adjacent tiers changes the total accepted-step
count by only 0--0.40%.

The few places where the controller does refine are not random.  They cluster
at crossings of the continuous but piecewise-smooth, log-temperature-linear
`Cstar` cache knots.  BASELINE, REFINED, and ULTRA resolve different subsets of
those derivative kinks.  The resulting offsets then propagate through long
one-step-per-output stretches.  This explains the opposite-signed B/R and R/U
errors and the long contraction plateaus near 0.91--0.93.  Existing evidence
does not yet uniquely separate the output-ceiling contribution from the
piecewise-smooth-knot contribution, so a bounded counterfactual is required.

Historical BA12 remains **FAIL** with maximum scaled difference
`1.7255917120989046`.  BA12R remains **FAIL** with its maximum
`d_RU/D_R=32.26041661869904`, maximum resolvable
`d_RU/d_BR=0.9329162771682704`, maximum ledger normalized difference
`7.611500166686796e-9`, source-minus-control stability utilization
`21.923678867554578`, and minimum ULTRA final-decade accepted step of 1 s for
both source and control.  All failure values are unchanged.  None is
reinterpreted as passing.

## Authentication and scope controls

- BA12R entry and final SHA:
  `971b2bac47c320bd262d8686841e679f703e5ccb`.
- R0-R3 recovery SHA:
  `d2822656aae8e264f4985f0d2f7228b90320a577`.
- Authenticated entry branch:
  `physics/phase6a1-bnv-recovery-ba12r`.
- Entry local, upstream, and live refs were clean and equal to the BA12R SHA.
- Canonical local, upstream, and live `master` remained clean at
  `bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`.
- Analysis branch:
  `analysis/phase6a1-ba12r-numerical-forensics`.
- Analysis worktree:
  `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-ba12r-numerical-forensics`.
- The branch was created from exact BA12R SHA `971b2bac...`.

The six trajectory hashes reauthenticated as:

| Mode | BASELINE | REFINED | ULTRA |
| --- | --- | --- | --- |
| source | `f331b5bbfd6bdfb0f8e93095052cf1edc712219180370d23383bed1c18b31093` | `1889336755a87377336e1b688e34032dab428304861bd2a63e142e27555b7a3a` | `584799310627541259cc1164c8a18f8950afb8dddd7f0f14a1f0b3e9dc557b52` |
| control | `be807f432c81b7f62594005b593bdd5bcfd301ddfd9ed2837178831fc2dbb9b8` | `be807f432c81b7f62594005b593bdd5bcfd301ddfd9ed2837178831fc2dbb9b8` | `2a742d0db6860e8606b94d3c7e0b9756d5656fdc3abe1a87452608c6ef4c55bb` |

No ODE trajectory, retry, tolerance level, test suite, test executable, or
candidate producer was run.  No production, test, CMake, baseline, trajectory,
or acceptance-threshold file was modified.  The scratch parser and its JSON
result are ignored build artifacts; their SHA-256 values are respectively
`14683d1810794f96c0a701fd45157ffa6bb2d951d4a988268a27833123d2b77c`
and `ef5170016224811571297fb20e48da6586a05df69d59b12d3c620d86f1167db2`.

## Exact execution architecture

The P2 call chain is:

1. `ba12r_ultra.cpp` selects the fixed P2 card and calls `Run` from
   `coupled_trajectory.cpp`.
2. `Run` constructs one `ControlledBnvSecularDriver`, one `EvolutionSystem`,
   and one `ScaledRKF45` object for the complete trajectory.
3. `Run` creates 8192 requested checkpoint times and divides them into nine
   batches: eight batches of 999 checkpoints and one batch of 200.
4. For each batch, `Run` calls `ScaledRKF45::Integrate` once.
5. Each `Integrate` call allocates new GSL RKF45 step, scaled-control, and
   evolve objects; sets `t=start`; and sets
   `h=min(1 s, 0.001*(first_checkpoint-start))`, which is exactly 1 s here.
6. Inside a batch the same GSL objects and updated `h` persist across all
   checkpoints.  They are not reset per output sample.
7. Every checkpoint is nevertheless passed to `gsl_odeiv2_evolve_apply` as
   `t1`.  GSL guarantees that `t1` is not exceeded and sets `t=t1` exactly on
   the final step.  Thus a final checkpoint-truncated step is forced at every
   requested output time.
8. At the next batch, the GSL objects are destroyed/reallocated and `h` is
   reset to 1 s.
9. The GSL RHS callback calls `EvolutionSystem::operator()`, which invokes the
   sole registered `ControlledBnvSecularDriver`.
10. `ControlledBnvSecularDriver` calls
    `FrozenControlledBnvRunContext::Evaluate`, which delegates its ordinary
    rotochemical term directly to `FrozenRotochemicalRunContext::Evaluate`.

`SecularEvolutionDriver` itself is **not** instantiated on this P2 path.  It is
another thin driver around the same `FrozenRotochemicalRunContext::Evaluate`
authority.  Therefore the relevant ordinary RHS is shared, but the P2 call
chain bypasses the `SecularEvolutionDriver` object.

Answers to the lifetime audit:

- `ScaledRKF45` lifetime: one object per complete trajectory.
- GSL step/control/evolve lifetime: one set per batch, nine sets per trajectory.
- `h` lifetime: persistent within one batch, discarded between batches.
- Initial `h`: exactly 1 s.
- Reset frequency: nine times, at output-interval indices
  `0,999,1998,2997,3996,4995,5994,6993,7992`.
- GSL objects reset every sample: **NO**.
- Adaptive history preserved across checkpoints: **YES within a batch; NO
  across the eight batch boundaries**.
- Output checkpoint forces a final truncated step: **YES at every output**.

## The 1-second steps

There are nine initialization-imposed 1 s accepted steps in each of the six
histories, one at the beginning of every batch.  The historical step schema did
not store every per-interval minimum, but the conclusion is exact: source code
assigns 1 s at all nine batch starts, and every batch-start interval has zero
rejections, so the initial 1 s attempt was accepted.  The ULTRA schema directly
confirms a 1 s interval minimum at exactly those nine indices and nowhere else.

| Evidence | accepted 1 s steps | location |
| --- | ---: | --- |
| BASELINE source | 9 | beginning of each batch |
| REFINED source | 9 | beginning of each batch |
| ULTRA source | 9 | beginning of each batch |
| BASELINE control | 9 | beginning of each batch |
| REFINED control | 9 | beginning of each batch |
| ULTRA control | 9 | beginning of each batch |

No serialized `last_step_s` equals 1 s; these steps occur at batch beginnings,
not output-interval endings.  They are **DRIVER/INITIALIZATION-IMPOSED**, not
solver-selected.  The final-decade 1 s BA12R trend failure therefore measures
batch construction, not late-time adaptive-step collapse.  That diagnosis does
not rescue BA12R or change its accepted gate result.

## Accepted-step scaling and output-grid dominance

| Mode/tier | accepted | rejected | 1-step intervals | 2-step intervals | 3+-step intervals |
| --- | ---: | ---: | ---: | ---: | ---: |
| source BASELINE | 8318 | 0 | 8183 (99.8901%) | 0 (0%) | 9 (0.1099%) |
| source REFINED | 8322 | 7 | 8180 (99.8535%) | 2 (0.0244%) | 10 (0.1221%) |
| source ULTRA | 8355 | 33 | 8164 (99.6582%) | 13 (0.1587%) | 15 (0.1831%) |
| control BASELINE | 8318 | 0 | 8183 (99.8901%) | 0 (0%) | 9 (0.1099%) |
| control REFINED | 8318 | 0 | 8183 (99.8901%) | 0 (0%) | 9 (0.1099%) |
| control ULTRA | 8329 | 12 | 8173 (99.7681%) | 9 (0.1099%) | 10 (0.1221%) |

Source ratios:

- `N_R/N_B = 1.0004808848280837`.
- `N_U/N_R = 1.0039653929343908`.
- `N_U/N_B = 1.004448184659774`.

Control ratios:

- `N_R/N_B = 1`.
- `N_U/N_R = N_U/N_B = 1.0013224332772301`.

For a smooth problem whose local-error tolerance controlled most RKF45 steps,
a 100-fold tolerance tightening would qualitatively produce a material step
increase (the usual embedded-method scaling suggests an order-unity factor,
not a required exact ratio).  Here the increases are 0--0.40%.  BASELINE's
8318 steps are exactly 8192 one-per-interval steps plus 14 extra startup steps
at each of nine batch restarts.  The output grid is therefore acting as an
integration grid, not merely an observation grid.  The sequences are
overwhelmingly checkpoint-limited, not tolerance-controlled.

## Signed state-error forensics

The following distributions include all 8193 checkpoints and therefore retain
the reported zero medians.  Percentiles are of absolute signed differences;
correlation and sign agreement use signed `e_BR` and `e_RU`.

| State/error | max absolute | median | p90 | p99 |
| --- | ---: | ---: | ---: | ---: |
| `x`, B/R | `1.35021969299709e-8` | 0 | `2.2033597169013355e-11` | `8.856804316614837e-9` |
| `x`, R/U | `5.964547356018812e-9` | 0 | `2.0012680401748673e-11` | `3.786425770702094e-9` |
| `eta_e`, B/R | `1.690127535411003e-14` | 0 | `1.669520732581746e-17` | `1.3891606337197277e-14` |
| `eta_e`, R/U | `9.307180607315274e-15` | 0 | `1.42169133301155e-17` | `6.640487021444414e-15` |
| `eta_mu`, B/R | `4.198143149608503e-14` | 0 | `3.237517684292367e-17` | `3.328098564355711e-14` |
| `eta_mu`, R/U | `2.075633431999896e-14` | 0 | `2.8342347101927535e-17` | `1.4012718762942866e-14` |

Maxima occur at:

- `x`: B/R index 31, `t=59709814453.125 s`; R/U index 163,
  `t=313958056640.625 s`.
- `eta_e`: B/R index 161, `t=310105810546.875 s`; R/U index 201,
  `t=387150732421.875 s`.
- `eta_mu`: B/R index 148, `t=285066210937.5 s`; R/U index 193,
  `t=371741748046.875 s`.

For every checkpoint where both signed differences are nonzero, their signs
are opposite: sign-agreement fraction is zero for all three states.  Pearson
correlations are `-0.6035442619082247` (`x`), `-0.8578018303874656`
(`eta_e`), and `-0.8232739632588751` (`eta_mu`).  ULTRA moves back across
REFINED toward BASELINE instead of continuing along a monotone refinement
sequence.

Fractions with `|e_RU| < |e_BR|` are 23.2149% (`x`), 23.9717% (`eta_e`), and
23.7276% (`eta_mu`).  Complementary `|e_RU| >= |e_BR|` fractions are 76.7851%,
76.0283%, and 76.2724%.  The latter include exact B/R/U equality at many late,
floor-limited checkpoints; the resolvable-only ratios below are the meaningful
contraction diagnosis.

## Resolvable contraction structure

| State | resolvable | median rho | p90 | p95 | p99 | maximum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `x` | 712 | `0.9082760028087328` | `0.9082788168450943` | `0.9329162594240223` | `0.9329162723377944` | `0.9329162771682704` |
| `eta_e` | 868 | `0.907646152761892` | `0.9082777979536288` | `0.908279362785158` | `0.9082922188056297` | `0.908293677374261` |
| `eta_mu` | 864 | `0.9081893863070942` | `0.9082800993461851` | `0.9082813739237146` | `0.9082935772615542` | `0.9082958835916675` |

By time decade, median `rho` evolves as follows:

| State | `1e10--1e11 s` | `1e11--1e12 s` | `1e12--1e13 s` |
| --- | ---: | ---: | ---: |
| `x` | `0.0242563` | `0.886557` | `0.908276` |
| `eta_e` | `0.123312` | `0.853132` | `0.908276` |
| `eta_mu` | `0.123290` | `0.873172` | `0.908277` |

Poor contraction is not a handful of isolated maxima.  Using the existing
accepted 0.10 criterion only as a diagnostic label, the principal contiguous
ranges are:

- `x`: indices 82--723, `1.5794208984375e11--1.392586962890625e12 s`
  (642 points), plus indices 12--30.
- `eta_e`: indices 119--879,
  `2.29208642578125e11--1.693062158203125e12 s` (761 points), plus 12--33.
- `eta_mu`: indices 119--875,
  `2.29208642578125e11--1.685357666015625e12 s` (757 points), plus 12--33.

Contraction degrades with age, enters a broad 0.91--0.93 plateau, and later
becomes floor-limited.  This is not an asymptotic convergence pattern.

## Control-only forensics

The BASELINE and REFINED control trajectory bytes and step bytes are
bit-identical.  Both chemical states remain exactly zero at ULTRA as well.
Nevertheless ULTRA changes the thermal state:

- Maximum control `x` drift from common B/R:
  `+3.7664721341812424e-10` at index 600,
  `t=1155673828125 s`.
- Nonzero control `x` differences: 8135 of 8193 checkpoints.
- Median absolute drift: `2.308113700166814e-10`.
- p90/p95/p99 absolute drift:
  `3.48967277297163e-10`, `3.5824526678496227e-10`,
  `3.705757811722776e-10`.
- Maximum control ledger drift: `Pnet`,
  `1.4904755921078965e24 erg/s`, normalized
  `9.483035685919317e-10`, at the same checkpoint.
- `Lgamma` supplies `1.4075130900186447e24 erg/s` of that difference;
  equilibrium neutrino differences are smaller.
- ULTRA uses 11 more accepted and 12 more rejected steps than the common B/R
  history; accepted-step increments differ in 10 intervals.

The control result proves that the ULTRA drift is not specifically caused by
BNV forcing or chemical feedback.  It is present in the source-free thermal
equation and shares the cache-knot/step-history pattern described below.

## Local failure windows

Seven-point windows (center +/- 3 output intervals) were inspected around all
three required locations.

### A. `eta_mu` stability maximum, index 216

- `t=416042578125 s`.
- `e_BR=+2.4715452183430588e-14`;
  `e_RU=-1.869212699188061e-14`.
- `x` contraction ratio in the same row: `0.9329162710278092`.
- REFINED `Tinf=160206314.97110257 K`.
- REFINED RHS:
  `x_dot=4.960354124113848e-13 s^-1`,
  `eta_e_dot=1.1839004772919047e-18 MeV/s`,
  `eta_mu_dot=2.5206623379868563e-18 MeV/s`.
- All three tiers take one accepted step equal to the full
  `1926123046.875 s` output spacing in every inspected interval, with zero
  rejections and no batch restart.
- Temperature, reaction rates, heat capacity, and powers vary smoothly;
  `sigma_e` and `sigma_mu` are constant; no serialized regime changes occur.

### B. `x` contraction maximum, index 179

- `t=344776025390.625 s`.
- `e_BR=+5.750012188610043e-9`;
  `e_RU=-5.36427996467026e-9`; `rho=0.9329162771682704`.
- REFINED `Tinf=153892714.50611147 K` and
  `Pnet=1.9693808600965022e34 erg/s`.
- All inspected tiers/intervals use one full-spacing accepted step, zero
  rejections, and no restart.  All diagnostics are smooth and monotone.

### C. `Pnet` ledger maximum, index 168

- `t=323588671875 s`.
- `Pnet_RU=+3.544875091136552e26 erg/s`.
- `x`: `e_BR=+6.185522871682281e-9`,
  `e_RU=-5.770574795782579e-9`, `rho=0.932916249037028`.
- REFINED `Tinf=151751878.35562435 K` and
  `Pnet=2.0593858651881065e34 erg/s`.
- The seven-point window again has one full-spacing step per tier/interval,
  zero rejections, no restart, constant sigma, and no regime transition.

None of the three maxima is itself a discontinuity, restart, rejection, or
step-count transition.  Their offsets were seeded earlier and propagated.

## RHS smoothness and cache-knot audit

| RHS feature | Classification | Finding |
| --- | --- | --- |
| linear controlled `B(t)` and fixed source | SMOOTH | affine in time; no state branch |
| moving-reference projection and actual potential | SMOOTH | linear algebra in source and eta |
| P2 full-retention partition | SMOOTH | fixed event mapping; no state branch |
| modified-Urca imbalance functions | SMOOTH | finite polynomials in `xi`; smooth through zero |
| fixed Me/Mmu process selection | SMOOTH | `L==0` branch fixed at construction, not crossed |
| iron envelope relation | SMOOTH | positive-domain power law; fixed model branch |
| other-neutrino contribution | SMOOTH | fixed disabled branch, identically zero |
| `Cstar(Tinf)` cache lookup | PIECEWISE_SMOOTH | continuous linear interpolation in `log(T)` with derivative kinks at 160 cache knots |
| frozen-validity interpolation | PIECEWISE_SMOOTH / FAIL-CLOSED ONLY | changes diagnostics and may refuse; does not change accepted RHS while valid |
| finite/domain/currentness checks | FAIL-CLOSED ONLY | no clipping or alternate accepted RHS |
| cache endpoint clamps | NONSMOOTH | not reached; trajectory remains inside `1e-5--1 MeV` cache domain |

The adaptive events correlate exactly with `Cstar` knot crossings:

- Source REFINED adds/rejects steps near cache knots 94, 95, and 96.
- Source ULTRA adds/rejects steps near knots 94--101, at output indices near
  11, 30, 53, 81, 116, 162, 227, and 354.
- Control ULTRA adds/rejects steps near knots 93, 92, 91, 90, 88, 87, and 84
  as its temperature falls.
- At each listed event the checkpoint temperature is within about 0.004--0.37%
  of the corresponding cache knot.

The required failure centers are 0.76--2.14% from their nearest cache knot, so
there is no branch transition at the center itself.  The most relevant source
ULTRA knot-99 refinement occurs at indices 162--163; its offset then propagates
through the Pnet and `x` failures at 168 and 179.  Knot 100 is crossed at
227--229, after the `eta_mu` stability maximum.  Thus the RHS is continuous,
but a piecewise-smooth thermal coefficient seeds tier-dependent adaptive
histories.

## Serialization and comparison audit

- `Capture` writes with `std::setprecision(17)`, equal to binary64
  `max_digits10`.
- Parsing those values as binary64 recovers the original serialized owner
  values exactly; declared `Q=0` is justified.
- `x_state` is written directly from the evolved thermal state.
- `eta_e` and `eta_mu` are read from the evolved chemical state and copied
  without unit conversion.
- State comparison performs no unit conversion.
- Derived temperature, RHS, and ledger quantities are recomputed by the same
  diagnostic owner at each accepted checkpoint, but those recomputations do
  not alter the three state values.

Even a conservative 64-ulp floor is far below the observed R/U differences.
The largest ratios `d_RU/(64 ulp)` are approximately 2.71 million (`x`),
3.99 million (`eta_e`), and 4.14 million (`eta_mu`).  Serialization and
binary64 roundoff are not the limiting mechanism.

## Local controller and global error

`ScaledRKF45` calls:

```text
gsl_odeiv2_control_scaled_new(
    eps_abs = 1,
    eps_rel = rtol,
    a_y = 1,
    a_dydt = 0,
    scale_abs = {atol_x, atol_eta_e, atol_eta_mu},
    dim = 3)
```

Therefore GSL's actual per-step requested scale is exactly:

```text
D_i(local) = atol_i + rtol |y_i|.
```

There is no derivative term, hidden component scale, or active `h|y'_i|`
term.  BA12R's formula has the same coefficients and structure, but substitutes
the cross-solution checkpoint magnitude
`M=max(|y_B|,|y_R|,|y_U|)` for the single trial step's `|y_i|`.  It is therefore
not the controller's recorded local error level, even though it is a
conservative checkpoint-scale analogue.

An accepted GSL step establishes only that the embedded local error estimate
for that proposed step satisfied the controller heuristic.  It does not bound
the accumulated global error, prove that two trajectories have entered an
asymptotic convergence regime, or make `D_R` a global-error denominator.
Consequently `d_RU/D_R >> 1` does not by itself prove an incorrect individual
step.  It proves failure of the separately accepted global cross-solution
BA12R criterion.  Local controller compliance and global hierarchy convergence
must remain distinct.

## Step-history correlation

At the three required maxima, all tiers have the same one-full-spacing-step
pattern, so there is no local correlation with rejection, restart, or step
count.  Whole-grid Pearson correlations between `|e_RU|` and ULTRA accepted
increments are 0.019 (`x`), 0.005 (`eta_e`), and 0.006 (`eta_mu`); correlations
with rejection increments are 0.064, 0.019, and 0.022.  `|Pnet_RU|`
correlations are 0.007 with accepted increments and 0.031 with rejections.
Restart correlations are between -0.006 and -0.008.

These small pointwise correlations are expected because a brief cache-knot
refinement seeds an offset that persists for hundreds of later intervals.
The causal pattern is temporal rather than coincident: the long poor-`rho`
ranges begin immediately after tier-specific refinement clusters (notably
ULTRA indices 81 and 116), and the largest broad plateau follows knot 99 at
162--163.  The nine 1 s restarts occur outside the early failure region except
at trajectory start and do not explain its maxima.

## Hypothesis classification

| Hypothesis | Classification | Evidence |
| --- | --- | --- |
| H1 — checkpoint-level reinitialization / `h` reset | SUPPORTED, LIMITED | exact nine batch resets explain all 1 s minima, but no reset occurs near the state/ledger maxima |
| H2 — output-grid segmentation dominates | SUPPORTED | 99.65--99.89% one-step intervals; 0--0.40% step-count growth under 100x tolerance tightening; every output is a hard `t1` ceiling |
| H3 — local tolerance is not a global estimator | SUPPORTED | controller scale is local; BA12R denominator is a cross-solution global diagnostic |
| H4 — binary64/roundoff plateau | DISFAVORED | differences are millions of 64-ulp floors and structured across long ranges |
| H5 — nonsmooth RHS / branch transitions | SUPPORTED AS PIECEWISE-SMOOTH KNOT EFFECT | adaptive events align with `Cstar` derivative kinks; no discontinuity or regime branch occurs at failure centers |
| H6 — stiffness or inadequate RKF45 | UNRESOLVED, NOT CURRENTLY SUPPORTED | low rejection fractions and no late step collapse; segmentation prevents a fair method assessment |
| H7 — serialization/comparison artifact | DISFAVORED | 17-digit round trip, `Q=0`, direct state serialization |
| H8 — actual implementation bug | SUPPORTED IN THE PHASE-6 VALIDATION ARCHITECTURE | output requests act as integration ceilings, defeating the intended tolerance-family interpretation; no physics-RHS bug demonstrated |
| H9 — other | SUPPORTED | tier-specific resolution of heat-capacity cache knots seeds opposite-signed offsets that then propagate on the checkpoint grid |

Another `rtol=1e-13` tier would not answer the identified question.  The
existing 100-fold tightenings barely change the step mesh because the output
spacing caps almost every step.  A still tighter tier would merely resolve
more cache knots or their neighborhoods while retaining the same confounding
segmentation.  No tighter tier should be authorized before the integration
ceiling is isolated.

## Minimal proposed diagnostic — not executed

The single best next experiment is a two-arm, bounded P2-source replay using
the same RKF45 method and the existing ULTRA tolerance.  It changes no physics,
tolerance, state definition, or acceptance threshold; it changes only whether
intermediate output times are passed to GSL as integration ceilings.

Predeclared configuration:

- Platform: the same authenticated local Mac toolchain as BA12R.
- Solver: GSL `gsl_odeiv2_step_rkf45` with the existing scaled controller.
- Tolerance: exact ULTRA `rtol=1e-11`,
  `atol=(1e-16,1e-22,1e-22)`.
- Source/card: exact `CPL-P2-LINEAR-QSS-v1` P2 source; all drive, boundary,
  frozen-context, and process settings unchanged.
- Initial time/state: `t0=0`, `x=0`, `eta_e=0`, `eta_mu=0`; this is the exact
  existing initial state, not a reconstructed mid-trajectory state.
- Final time: existing checkpoint index 240,
  `t1=462269531250 s`.
- Initial `h`: exactly 1 s once per arm.
- GSL step/control/evolve objects: persistent for the entire bounded arm.
- Arm S (segmented control): request the existing 240 checkpoint endpoints
  `t_k = k*1926123046.875 s`, `k=1,...,240`; record all 241 states including
  `t0`.  This must reproduce the archived ULTRA states and step history through
  index 240 exactly; otherwise the diagnostic is invalid.
- Arm E (endpoint-only counterfactual): give GSL only `t1` as the integration
  ceiling; record state and diagnostics only at `t0` and `t1`, while retaining
  every accepted/rejected internal-step statistic.  No intermediate output
  ceiling is permitted.
- Quantities compared at `t1`: exact binary64 `x`, `eta_e`, `eta_mu`; recomputed
  `Pnet`; accepted/rejected counts; minimum/maximum/internal step sequence; and
  the cache-knot crossing locations inferred from accepted steps.

Predeclared decision criterion, using the existing BA12R floor without a new
threshold:

1. Arm S must reproduce archived ULTRA evidence through index 240 exactly.
2. Define `F_i=max(D_U,i,64 ulp(M_i),Q_i)` exactly as BA12R, with `Q_i=0`.
3. If `|y_E-y_S| > 10 F_i` for any state at `t1`, checkpoint segmentation has
   a resolvable causal effect and H2 is confirmed as a recovery target.
4. If all three differences are `<=10 F_i`, H2 is falsified as the proximate
   cause on this window; the next separately authorized diagnostic should be a
   bounded `rk8pd` oracle over the same exact window and initial state to test
   H5/H6.  Do not authorize that second method in advance.

This first diagnostic comprises exactly two short source integrations and no
control, tolerance sweep, full trajectory, or candidate production.  It can be
implemented in Phase-6 test/validation space by calling the existing solver
with two checkpoint vectors.  It does **not** require changing
`ScaledRKF45`.  A production architecture change to decouple save cadence from
integration ceilings would require a later owner decision only if the
diagnostic confirms H2.

## Phase-5D and platform protection

Phase-5D baselines and governed bytes must remain unchanged.  The proposed
diagnostic belongs in Phase-6 test/validation space.  If subsequent evidence
shows the shared solver API itself cannot support passive sampling correctly,
that is a separate upstream numerical-governance issue; it must not be folded
into a Phase-6 candidate change.

The checkpoint ceiling, batch reset, and cache interpolation are
platform-independent algorithmic behaviors.  Exact last-bit outcomes may vary
slightly by platform, but roundoff is disfavored by several million-fold
margins.  The diagnostic should remain on the same local Mac platform as
BA12R.  Cluster qualification remains separate and is not implicated or
authorized.

## Explicit final questions

1. **Why are accepted-step counts almost equal?** Because almost every output
   interval is accepted as one full checkpoint-truncated step; tolerance only
   changes behavior near a few `Cstar` cache knots.
2. **Why is there a 1 s step in the final decade for both modes?** The ninth
   999-checkpoint batch begins inside the final decade and resets `h` to 1 s.
3. **Adaptive or imposed?** Imposed by `Integrate` initialization.
4. **Is the output grid controlling stepping?** Yes, overwhelmingly.
5. **Does `ScaledRKF45` preserve adaptive history?** Within each batch only;
   not across the nine `Integrate` calls.
6. **Does the GSL controller exactly use diagnostic `D_i`?** It exactly uses
   `atol_i+rtol|y_i|` locally with no derivative term.  BA12R uses the same
   formula with cross-solution `M`, so it is not literally the recorded local
   scale and is not a global-error estimator.
7. **Are B/R/U asymptotic?** No.
8. **Is roundoff limiting?** No; disfavored by multi-million 64-ulp margins.
9. **Is the RHS nonsmooth near failures?** No discontinuity at the failure
   centers.  Earlier continuous derivative kinks in the `Cstar` cache are
   strongly implicated.
10. **Does control exhibit ULTRA drift?** Yes: thermal `x` and its ledger drift;
    chemical states remain exactly zero.
11. **Is serialization relevant?** No.
12. **Is RKF45 unsuitable?** Not established.  The surrounding checkpoint
    architecture prevents a fair tolerance-convergence assessment; stiffness
    is not presently supported.
13. **Best next run?** The exact two-arm bounded segmented-versus-endpoint-only
    ULTRA P2-source experiment above.
14. **Does it require a new architectural decision?** The diagnostic requires
    only an owner-approved Phase-6 validation runner, not a shared-solver
    change.  A recovery architecture decision follows only if H2 is confirmed.
15. **Required authorization?** Explicit approval to add a Phase-6-only
    diagnostic runner and execute exactly those two bounded P2-source ULTRA
    integrations on the local Mac, with the stated initial state, endpoint,
    checkpoint vectors, comparison floor, and no candidate or merge.

## Recommended next action

Owner review should decide whether to authorize exactly the two-arm bounded
diagnostic above.  Do not retry ULTRA, introduce `rtol=1e-13`, change the
output grid of the governed campaign, change `ScaledRKF45`, use the cluster,
or begin candidate work under the present authorization.
