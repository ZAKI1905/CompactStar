# Phase-6A-1 checkpoint reconstruction preflight

Status: **PROPOSED — OWNER ACCEPTANCE REQUIRED BEFORE NUMERICAL ORACLE
EXECUTION**

Disposition: **CHECKPOINT RECONSTRUCTION PREFLIGHT COMPLETE — BOUNDED
CANDIDATE/ORACLE EXPERIMENT READY FOR OWNER REVIEW**.

Classification: **DOCUMENTATION-ONLY NUMERICAL-METHOD PREFLIGHT; NOT AN
IMPLEMENTATION; NOT A BNV CANDIDATE; NOT A GOVERNED BASELINE; NOT A PHYSICAL
RESULT**.

## 1. Authenticated entry, authority, and exclusions

This preflight was prepared on
`analysis/phase6a1-checkpoint-reconstruction-preflight` in
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-checkpoint-reconstruction-preflight`
from exact passive-observation SHA
`9af8912107ea77a2e2ea51517c1776f26a4b7b49`. At entry, the passive branch was
clean and its local, upstream, and live refs equaled that SHA. Canonical local
`master`, `origin/master`, and live `refs/heads/master` were clean and equal to
`bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`.

The passive experiment proved scheduling passivity only: its final state,
accepted/rejected counts, and all 232 accepted-step records exactly matched
Arm E, while all 241 requested observation indices were assigned to brackets
without becoming integration ceilings
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:275-305`). It
explicitly left reconstruction for a separately approved, independently
oracled task
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:313-328`).

Historical status is unchanged:

- BA12: **FAIL**.
- BA12R: **FAIL**.
- segmentation diagnostic: **CONFIRMED NUMERICAL EFFECT**;
- passive observation scheduling: **PASS**.

The previous BA12R state, ledger, matched-difference, and solver-trend failures
remain recorded without reinterpretation
(`docs/validation/PHASE6A1_CONTROLLED_BNV_BA12R.md:148-211`,
`docs/validation/PHASE6A1_CONTROLLED_BNV_BA12R.md:259-295`). No ODE
integration, reconstruction, oracle, interpolation implementation, test,
candidate production, or merge was performed for this preflight. Production
source, tests, CMake, `ScaledRKF45`, Cstar, Phase-5D, baselines, EOS/data, and
literature are unchanged.

## 2. Numerical object and ownership separation

The authoritative numerical evolution is one uninterrupted adaptive RKF45
initial-value solution represented by accepted endpoint records

```text
(t_n, y_n), (t_(n+1), y_(n+1)), ...
```

with no requested scientific output time supplied to the owner integrator as
an intermediate `t1`. For a requested time

```text
t_n < t_obs < t_(n+1),
```

checkpoint reconstruction is a separate post-hoc estimate of `y(t_obs)` from
the immutable accepted-step evidence. It is not an accepted RKF45 step, not a
new authoritative trajectory, and not permission to change the endpoint
history.

The boundary is strict:

```text
MAIN INTEGRATION
  owns: accepted endpoints, accepted/rejected adaptive history, final state
  may see: final integration ceiling only
  must not see: scientific observation schedule as t1 ceilings

POST-HOC RECONSTRUCTION
  consumes: immutable accepted endpoints and authenticated physics inputs
  owns: reconstructed states plus method/oracle evidence
  must not mutate: main state, main h, main GSL objects, main caches, or main
                   currentness/validity state
```

The passive probe already demonstrated the required main-path scheduling
shape: one persistent step/control/evolve/`h`, one final positive target, and
observer notification only after acceptance
(`tests/bnv/passive_observation_probe.cpp:215-280`). The future experiment may
copy accepted endpoint states after acceptance, but it may not call a
reconstruction method from inside the active GSL evolution.

## 3. Accepted-step information audit

The following table fixes what the future reconstruction runner may use.

| Information | Status | Required treatment |
| --- | --- | --- |
| `t_n`, `t_(n+1)` | already available | exact binary64 accepted endpoint times |
| `y_n`, `y_(n+1)` for `(x,eta_e,eta_mu)` | already available/capturable | initial state plus the accepted-state record; copy only after acceptance |
| accepted step `H_n=t_(n+1)-t_n` | already available | derive from exact endpoints; retained evidence also serializes it |
| suggested next `h` and cumulative rejections | already available | diagnostic only; never seed interpolation |
| `f_n`, `f_(n+1)` | not retained as values; computable post hoc | evaluate only in an isolated scratch context after the main run |
| Cstar cell before/after and number of crossed knots | already available for the bounded evidence | classify from `x` using the exact 160-point log-temperature grid |
| source/process/frozen-validity identities | diagnostic owners already provide them | re-evaluate from the reconstructed immutable state in the isolated context |
| other piecewise branch identities | fixed or owner-evaluable | serialize every material branch used by the diagnostic owner |
| GSL internal RKF45 stage vectors | not exposed by the current adapter/API | unavailable without shared-solver/GSL-internal work; not authorized |
| rejected trial states | not available | cumulative count only; not required by any candidate |

The retained passive evidence has 232 accepted-state rows containing all three
state components, while the historical internal projection has endpoint times,
thermal state, step size, and Cstar cells
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:291-309`). The
future experiment must authenticate those exact hashes before use:

```text
accepted states:
  7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459
observation brackets:
  0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d
internal accepted steps:
  fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8
schedule:
  43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e
```

No GSL-stage-dependent dense-output proposal is available under the present
architecture. Reconstructing or retaining GSL internal stages would require a
new shared-solver decision and is outside this plan.

## 4. RHS currentness and side-effect audit

Endpoint RHS evaluation is **not side-effect free on the live run objects**,
even though the returned RHS is deterministic for fixed authenticated inputs.

1. `EvolutionSystem::operator()` unpacks the supplied flat state into its
   referenced logical `StateVector` and clears/repopulates its shared
   `RHSAccumulator` (`CompactStar/Physics/Evolution/src/EvolutionSystem.cpp:93-146`).
2. The controlled driver invokes the same governed ordinary RHS and then the
   Phase-6 source/ledger owners
   (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:92-138`).
3. The ordinary RHS evaluates photon and other-neutrino owners through the
   shared `StarContext`
   (`CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:130-141`).
4. `StarContext::HeatCapacityStar_Tinf` may build its mutable 160-point Cstar
   cache and always updates the mutable interpolation hint `last_i` inside the
   cache domain
   (`CompactStar/Physics/Evolution/src/StarContext.cpp:753-801`).
5. The heat-capacity builder itself mutates the cached table
   (`CompactStar/Physics/Evolution/src/StarContext.cpp:807-887`). Other thermal
   diagnostic owners also maintain profile-keyed caches.
6. The reaction-free-only BNV path has an additional mutable stored reference,
   although the bounded P2 source uses coupled mode
   (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:116-125`).

Therefore no Hermite derivative, replay, oracle, or reconstructed diagnostic
may be evaluated through the context/EvolutionSystem later consumed by the
main trajectory. The future experiment must first finish and authenticate the
main endpoint archive, then build isolated scratch contexts from the same
authenticated inputs. Each numerical method/tier receives its own process or
fresh method-owned context/output root. Context construction, cache warming,
method order, and output serialization must be fixed and recorded. No method
may rely on cache state left by another method.

This separation makes post-hoc endpoint RHS evaluation acceptable: its mutable
effects are confined to disposable reconstruction contexts and cannot alter
the authoritative accepted-step history.

## 5. Candidate A — componentwise linear reconstruction

For

```text
H = t_(n+1)-t_n,
s = (t_obs-t_n)/H,
```

define

```text
y_A(t_obs) = (1-s)y_n + s y_(n+1).
```

With a twice-differentiable solution, the interpolation error is `O(H^2)`.
It is deterministic, requires no RHS evaluation, is cheap, and cannot mutate
the main path. Its disadvantages are a low formal order, no use of the known
ODE direction, and potentially material bias in long accepted steps. A Cstar
derivative kink does not introduce a state discontinuity, but it prevents any
claim of better-than-second-order behavior. Linear reconstruction is retained
as a basic/negative candidate; it is not presumed adequate.

## 6. Candidate B — cubic Hermite reconstruction

Let `f_n=f(t_n,y_n)`, `f_(n+1)=f(t_(n+1),y_(n+1))`, evaluated in the isolated
post-run context, and use

```text
h00(s) =  2s^3 - 3s^2 + 1
h10(s) =    s^3 - 2s^2 + s
h01(s) = -2s^3 + 3s^2
h11(s) =    s^3 -   s^2

y_B(t_obs) = h00 y_n + h10 H f_n + h01 y_(n+1) + h11 H f_(n+1).
```

With exact endpoint data and a sufficiently smooth solution this has
`O(H^4)` interpolation error. Endpoint RHS values must be evaluated once per
unique endpoint, post-run, in method-owned scratch contexts. They must not be
inserted into or called from the live main loop.

Cstar is continuous and linearly interpolated in `log(T)`, but its derivative
changes at the 160 grid knots
(`CompactStar/Physics/Evolution/src/StarContext.cpp:793-801`,
`CompactStar/Physics/Evolution/src/StarContext.cpp:836-854`). Across such a
knot the RHS is continuous but not generally differentiable. The solution is
expected to remain `C1`, while higher time derivatives can jump. A single
cubic spanning the kink therefore loses its smooth-interval fourth-order
guarantee and can degrade to approximately second-order local behavior. This
is why the knot-crossing stratum has a separate zero-failure requirement.

Hermite is cheap relative to replay after its endpoint derivatives exist, but
it is not selected merely for cost.

## 7. Candidate C — isolated local RKF45 replay

For every strict-interior `t_obs`, construct a new scratch integrator at exact
`(t_n,y_n)` and integrate only to `t_obs`. The sole replay ceiling is `t_obs`;
there is no intermediate output grid. Set the initial replay step exactly to

```text
h0 = t_obs - t_n.
```

The adaptive controller may reject/reduce that proposal. A replay stops if it
requires more than 100000 `gsl_odeiv2_evolve_apply` calls, encounters a GSL or
nonfinite failure, or trips any source/currentness/frozen-validity gate. An
observation exactly equal to `t_(n+1)` returns the stored accepted endpoint
bit-for-bit and performs no replay.

The fixed candidate hierarchy is:

| level | GSL method | `rtol` | `atol(x,eta_e,eta_mu)` | role |
| --- | --- | ---: | --- | --- |
| Replay-1 | RKF45 | `1e-11` | `(1e-16,1e-22,1e-22)` | candidate witness |
| Replay-2 | RKF45 | `1e-12` | `(1e-17,1e-23,1e-23)` | candidate value |

Both use the same scaled-controller convention as the existing adapter:
`D_i=atol_i+rtol*|y_i|`
(`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp:12-18`). Replay-2 must meet
the common oracle budget below; additionally require
`|y_Replay2-y_Replay1|<=F_i` at every component/observation. Failure
disqualifies local replay; its tolerances may not be changed after results.

Local replay is scientifically valid as a reconstruction candidate because it
solves the same IVP from an exact accepted left endpoint without affecting the
owner trajectory. It is not dense output of the original step and will not
reproduce its internal stages. Its cost and independent-oracle agreement must
therefore be evaluated rather than assumed.

## 8. Independent rk8pd oracle

The independent numerical oracle shares the physics RHS, exact left endpoint,
EOS, coefficients, source history, and validity rules, because it must solve
the same local IVP. It differs in the integration algorithm/order: GSL rk8pd
is the installed explicit embedded Prince-Dormand 8(9) method, whereas the
replay candidate uses RKF45 and the algebraic candidates use no integrator
(`/opt/local/include/gsl/gsl_odeiv.h:104-115`). This independence can expose
reconstruction/integration error; it cannot expose a shared physics-RHS,
coefficient, cache-mathematics, or accepted-left-endpoint defect.

For each strict-interior observation, both oracle levels start at exact
`(t_n,y_n)`, use only `t_obs` as the local ceiling, and set
`h0=t_obs-t_n`. The exact hierarchy is:

| level | GSL method | `rtol` | `atol(x,eta_e,eta_mu)` |
| --- | --- | ---: | --- |
| Oracle-1 | rk8pd | `1e-12` | `(1e-17,1e-23,1e-23)` |
| Oracle-2 | rk8pd | `1e-13` | `(1e-18,1e-24,1e-24)` |

For the single category-D observation that coincides bit-for-bit with an
accepted endpoint, both oracle values are defined to be that stored endpoint,
`d_O=U_O=0`, and no local integration is performed. This exact-return rule
also applies to every candidate method at that observation.

The retained endpoint-state magnitudes on the bounded interval reach
approximately `0.4924000883`, `4.8611e-7 MeV`, and `1.1121e-6 MeV`. At those
magnitudes the Oracle-2 requested scales are respectively about
`4.9241e-14`, `4.8612e-20 MeV`, and `1.1121e-19 MeV`: approximately 13.86,
7.17, and 8.21 times the corresponding 64-ulp floors. Oracle-2 is therefore
tighter than ULTRA while still meaningfully above the declared binary64
resolution floor. No `rtol=1e-14` tier is proposed because its thermal scale
would approach that floor.

## 9. Oracle self-qualification and resolution

For component `i` and observation `j`, define

```text
d_O(i,j) = |y_Oracle2(i,j)-y_Oracle1(i,j)|
M_O(i,j) = max(|y_Oracle1(i,j)|,|y_Oracle2(i,j)|)

D_O1(i,j) = atol_Oracle1,i + rtol_Oracle1 M_O(i,j)
D_O2(i,j) = atol_Oracle2,i + rtol_Oracle2 M_O(i,j)
F_O(i,j)  = max(D_O2(i,j), 64 ulp(M_O(i,j)), Q_O(i,j))
U_O(i,j)  = 2 max(d_O(i,j), F_O(i,j)).
```

`Q_O=0` when in-memory values or max-digits-10 round-tripped binary64 values
are compared; otherwise it is the summed serialization half-ulp uncertainty.
The factor two is a conservative oracle uncertainty rather than an assertion
that the two-level difference is an exact error estimator.

The common BA12R-scale `F_i` and reconstruction budget `B_i` are defined in
section 10. The oracle self-qualifies only if, for every component/observation,

```text
d_O <= D_O1
U_O <= 0.20 F_i.
```

If either condition fails anywhere, the oracle is unresolved and the complete
experiment stops. No candidate is labeled accurate or inaccurate at that
point, and no tolerance is altered.

## 10. Exact reconstruction error budget

The accepted BA12R state resolution scale is retained
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:275-319`). For
each candidate `C`, define a conservative magnitude that a bad candidate
cannot inflate:

```text
M_i = max(|y_n|,|y_(n+1)|,|y_Oracle1|,|y_Oracle2|)
D_U,i = atol_ULTRA,i + rtol_ULTRA M_i
F_i = max(D_U,i, 64 ulp(M_i), Q_i)
B_i = 0.25 F_i,
```

where `rtol_ULTRA=1e-11`,
`atol_ULTRA=(1e-16,1e-22,1e-22)`, and `Q_i=0` for in-memory or exact
max-digits-10 binary64 comparisons. Candidate `C` passes a state point only if

```text
E_C(i,j) = |y_C(i,j)-y_Oracle2(i,j)| + U_O(i,j) <= B_i(i,j).
```

This is the one predeclared state-reconstruction rule. It is stricter than
mechanically requiring `E_recon<=F_ULTRA`. Two reconstructed tiers can add at
most `0.5F` to an adjacent-tier state difference, only 5% of BA12R's
floor-limited `10F` allowance. A matched source-minus-control adjacent-tier
comparison can add at most four reconstruction errors, `F`, or 10% of `10F`.
The same bound is also far below the refined stability scale `D_R` because
`F` is normally governed by the 100-times-tighter ULTRA request. Thus
reconstruction cannot dominate refined/ULTRA stability, B/R/U contraction, or
matched differences.

For each power/ledger observable, retain BA12R's scale, evaluated from the
Oracle-2 governed diagnostic packet so a bad candidate cannot enlarge its own
budget:

```text
G_P = max(1 erg/s,
          |P_dir_actual|+|L_H|+|DeltaLnu|+|Lnu_eq|
          +|Lgamma|+|Lother|+|L_out_fluid|)
F_P = max(1e-11 G_P, 64 ulp(M_P), Q_P)
B_P = 0.25 F_P.
```

Here `M_P=max(|O_Oracle1|,|O_Oracle2|)` for the individual observable's ulp
floor. The oracle diagnostic uncertainty is formed by the same two-level rule, must
be `<=0.20F_P`, and every candidate diagnostic must satisfy
`|O_C-O_Oracle2|+U_O,diagnostic<=B_P`. This is applied separately to every
accepted BA12R ledger observable, not to a combined heating scalar
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:331-377`).

For integrated R20 evidence, additionally require the candidate-versus-
Oracle-2 composite-trapezoid luminosity integral difference and thermal/end
state propagation difference each to be `<=5e-6 N_R20`, one tenth of the
existing `5e-5` subsidiary budgets. The governed `|R20|/N_R20<=2e-4` gate is
unchanged.

## 11. Bounded validation set and Cstar stratification

The future experiment uses all 240 positive observation times in the
authenticated indices `1...240` schedule; index 0 supplies the exact initial
state. No result-selected subset is permitted. The schedule identity and size
are governed by the passive record
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:76-82`).

The authenticated passive brackets and internal-step records were joined by
accepted-step ordinal without evaluating a reconstruction. Categories A-C are
mutually exclusive. D and E are additional, potentially overlapping labels.

| Category | Exact definition | Positive observations |
| --- | --- | ---: |
| A — no Cstar knot | bracketing accepted step has `cstar_knots_crossed=0` | 237 |
| B — one Cstar knot | bracketing accepted step has `cstar_knots_crossed=1` | 3 |
| C — multiple Cstar knots | bracketing accepted step has `cstar_knots_crossed>1` | 0 |
| D — accepted endpoint | `t_obs==t_(n+1)` bit-for-bit | 1 |
| strict interior | `t_n<t_obs<t_(n+1)` | 239 |
| E — deep inside a long step | strict interior, `H>=2212334854.861145 s` and `0.25<=s<=0.75` | 81 |

The exact long-step threshold is the retained median accepted step from the
passive experiment
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:257-269`). The
240 observations occupy 192 distinct accepted-step brackets; 48 accepted
steps contain two observations and none contains more than two. The complete
main history crosses seven single Cstar knots, but only three of those accepted
steps bracket a requested observation. Every method must report its maximum
budget utilization and failures separately for A, B, D, E, and all points.
Absence of category C is reported; it is not fabricated or substituted.

## 12. Candidate acceptance and selection policy

The future comparison evaluates exactly:

1. Candidate A: componentwise linear interpolation;
2. Candidate B: endpoint-value/endpoint-RHS cubic Hermite;
3. Candidate C: Replay-2 RKF45, with Replay-1 as its fixed witness; and
4. a derived hybrid report using Hermite for A and Replay-2 for B/C, computed
   only from already-generated candidate values.

A candidate is eligible only if:

- the oracle self-qualifies at every point;
- it has zero state, diagnostic, R20, currentness, validity, and identity
  failures over all 240 observations;
- it has zero failures in every populated Cstar stratum;
- it is deterministic under a fresh-process byte-for-byte repeat on the same
  authenticated platform; and
- its evaluation never mutates the retained main evidence or another method's
  context.

For selection, define the candidate's worst utilization as the maximum of all
state, diagnostic, and integrated errors divided by their corresponding
budgets. A method with worst utilization `<=0.25` is **comfortably passing**.
If one or more single methods are comfortably passing, select the least costly
and simplest among them in the fixed order linear, Hermite, local replay. This
tie-break is allowed only after each tied method has a factor-four scientific
margin. If none is comfortably passing but one or more pass at `<=1`, select
the single method with the smallest worst utilization; exact binary64 ties are
resolved by the same cost/complexity order. If no single method passes, no
method is adopted.

There is **no pre-result scientific winner**. Linear is the expected basic
negative candidate; Hermite is the expected efficient smooth-interval
candidate; Replay-2 is the expected robust but expensive candidate. Those are
expectations, not dispositions.

The hybrid is not preferred. It is eligible only if Hermite passes all A
points, fails only B/C points, Replay-2 passes those B/C points, and the
combined policy passes every common gate. Its branch is deterministic from
the authenticated Cstar-cell transition count. Even then, prefer any passing
single method under the rule above. The hybrid exists to diagnose whether
knot localization can avoid full replay, not to introduce complexity by
default.

## 13. Observation schedule invariant

Future Phase-6 semantics must enforce:

```text
same initial state + physics + solver + tolerances + final time
    => identical accepted main-step history and final state
       for every requested observation schedule.
```

Scientific schedules are metadata processed against an immutable endpoint
archive after acceptance or after the complete trajectory. They never become
main GSL ceilings and never alter `h`. Structural inspection must show one
main positive target. Schedule variants (the exact 241-row list, final-only,
and the future 8193-row list) must produce identical bracket assignments when
applied post hoc to the same archive. If a future runner executes more than
one otherwise identical main trajectory to test this invariant, those runs
require separate explicit owner authorization; this preflight authorizes none.

## 14. Reconstructed diagnostics and R20

The reconstruction sequence is:

1. reconstruct only the state vector `y(t_obs)`;
2. load that immutable state into a method-owned scratch diagnostic context;
3. invoke the same governed `ControlledBnvSecularDriver::Evaluate` and
   diagnostic owners used by the trajectory; and
4. serialize state, source/currentness/validity identities, Cstar cell, all
   separated ledger values, and error budgets.

The BNV owner already derives the moving source, actual potential, direct
ledger, beta ledger, Cstar, surface temperature, and frozen metrics from time
and state (`CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp:100-184`).
No power, potential, rate, validity metric, or temperature is interpolated
independently.

R20 retains the existing composite trapezoid on the requested scientific
checkpoint grid. The previous accepted R20 contract and subsidiary budgets
remain the authority
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:356-377`).
The segmentation finding concerned using that grid as integration ceilings,
not using it as an independently owned output quadrature. Once states and
diagnostics are qualified at the requested times, there is no present evidence
requiring an R20 quadrature redesign. An adaptive-step quadrature would be a
new numerical owner and is explicitly outside this plan.

## 15. Consequence for a future clean BA12R campaign

The historical segmented BASELINE/REFINED/ULTRA trajectories remain immutable
failure evidence but cannot serve as a clean passive-observation convergence
hierarchy. Observation ceilings were shown to change the adaptive path and
bounded endpoint state (`docs/validation/PHASE6A1_BA12R_SEGMENTATION_DIAGNOSTIC.md:196-254`).

After reconstruction is independently qualified and separately authorized, a
scientifically valid BA12R retry requires six new uninterrupted full-duration
trajectories:

```text
source:  BASELINE, REFINED, ULTRA
control: BASELINE, REFINED, ULTRA.
```

Each trajectory must use one persistent main integrator from start to finish,
the exact existing tier tolerances and run card, and the 8193-point schedule as
post-hoc metadata only. All three controls are required because matched
source-minus-control convergence must use the same uninterrupted architecture
and tier; historical equality of segmented BASELINE/REFINED controls cannot
prove that property. Historical files remain useful only as provenance and
failure evidence. This preflight does not authorize any of the six runs.

## 16. Compute cost and process-level parallelism

For the bounded 240-observation dataset:

- linear is `O(240*3)` arithmetic and negligible;
- Hermite is `O(240*3)` after at most 233 unique endpoint RHS evaluations;
- each local RKF45 level performs an adaptive local IVP per 239 strict-interior
  points;
- each rk8pd oracle level performs the same 239 local IVPs with the higher-stage
  8(9) method; and
- diagnostic evaluation is required at every reconstructed state for every
  retained candidate/oracle result.

The passive owner run used 1985 RHS evaluations and 301.84 s including complete
fixture/context assembly
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:245-269`). If each
local interval accepts an initial full-interval proposal, the two RKF45 levels
and two rk8pd levels are already several thousand stage evaluations; kinks or
rejections increase that count. Exact wall time must be recorded rather than
predicted from the owner-run aggregate, because context construction and the
first Cstar-cache build are substantial fixed costs.

Independent observation replays are process-parallel only when each worker has
its own fully constructed context, GSL objects, cache state, and output root.
Use a fixed deterministic partition by observation index and serialize final
assembly/hash/comparison. Linear and Hermite algebra can be parallelized but do
not need to be. No scientific thread-level sharing is permitted.

An 8192-positive-time reconstruction has 34.1333 times as many observations as
the bounded 240-point experiment. Linear/Hermite remain linear and likely
secondary to context/diagnostic cost; four levels of local replay can dominate
the future campaign. Cost informs the tie-break only after scientific budgets
pass.

## 17. Eventual module ownership and ADR disposition

The first implementation of this experiment should remain in Phase-6 BNV
test/validation space and consume the retained endpoint artifacts. That
validation-only implementation requires no new ADR because it does not become
production numerical authority.

If a reconstruction method is adopted for Phase-6 candidate trajectories, its
first production owner should be a narrow Phase-6 BNV checkpoint-output
component, separate from `ScaledRKF45`. Adoption changes numerical-output
architecture and should be recorded in a new narrow ADR before production use.
Do not automatically create that ADR in the validation task. Promotion to a
generic/shared numerical utility, or any modification of the shared solver,
would cross Phase-5 ownership and requires separate governance and regression.

There is no Phase-5D implication in the proposed experiment or a Phase-6-only
owner. Phase-5D source, producer, comparator, baseline, and trajectory semantics
remain unchanged.

## 18. Exact future experiment order

After explicit owner acceptance, and not before:

1. authenticate this preflight SHA, Arm E/passive hashes, platform, production
   library, scientific inputs, and the exact 241-row schedule;
2. implement only a Phase-6 test/validation reconstruction harness and verifier;
3. ingest the existing accepted endpoint archive; do not run a new main
   trajectory;
4. reproduce the section 11 category counts exactly;
5. construct independent isolated contexts with fixed method/tier ownership;
6. evaluate linear and Hermite candidates at all 240 positive observations;
7. execute Replay-1 and Replay-2 for the 239 strict-interior observations;
8. execute Oracle-1 and Oracle-2 rk8pd for the same 239 observations;
9. self-qualify the oracle before adjudicating any candidate;
10. evaluate states and governed diagnostics against the exact section 10
    budgets, by stratum and globally;
11. form the derived hybrid report without additional integration;
12. perform one fresh-process deterministic repeat of every candidate method
    and require byte-identical normalized outputs; linear and Hermite repeats
    are algebraic, while both Replay-1 and Replay-2 are repeated for all 239
    strict-interior points; Oracle-1/Oracle-2 are not repeated beyond their
    nested self-qualification;
13. record wall/user/system time, RHS evaluations, accepted/rejected local
    steps, failures, per-stratum maxima, hashes, and method disposition; and
14. stop and return to the owner. Do not implement a production owner, create a
    candidate, rerun BA12R, or merge automatically.

No free tolerance, subset, category, budget, method, tie-break, or initial-step
choice remains for post-result selection.

## 19. Stop conditions

The future reconstruction validation stops if:

- any retained hash, platform, source, schedule, or currentness identity fails;
- any new main trajectory is required under the reconstruction authorization;
- the passive main history or final state changes;
- an observation time becomes a main GSL ceiling;
- the Oracle-1/Oracle-2 hierarchy fails section 9 anywhere;
- any tolerance, initial step, subset, category, or budget would need
  result-dependent adjustment;
- reconstruction or diagnostics mutate main state/cache/currentness later
  consumed by the main trajectory;
- isolated method ordering changes a normalized output;
- Cstar cell/knot classification is ambiguous or differs from section 11;
- a candidate crosses an unidentified material RHS branch;
- diagnostic-owner values cannot be evaluated consistently from a
  reconstructed immutable state;
- R20 would require an unapproved quadrature owner;
- satisfying the budget would require modifying `ScaledRKF45`, Cstar,
  Phase-5D, or any governed baseline;
- historical BA12 or BA12R would need reinterpretation; or
- a physical BNV rate/model or candidate conclusion would be required.

## 20. Explicit review answers

1. **What endpoint information is available?** Exact accepted endpoint times,
   all three state components, accepted step, suggested next `h`, cumulative
   rejection count, thermal state before/after, and Cstar-cell crossings.
   Endpoint RHS is post-hoc computable; GSL stages and rejected trial states
   are unavailable.
2. **Is endpoint RHS evaluation side-effect free?** **NO** on the live objects.
   It rewrites logical state/accumulator and mutates caches. It is acceptable
   only after the main run in isolated disposable contexts.
3. **Linear expected order?** `O(H^2)` for a smooth solution.
4. **Hermite expected order?** `O(H^4)` on a sufficiently smooth interval with
   exact endpoint values/derivatives.
5. **Cstar-kink effect on Hermite?** The continuous derivative kink removes the
   smooth-interval guarantee and can reduce the local behavior to about second
   order; knot intervals must pass independently.
6. **Is local RKF45 replay valid?** Yes as an isolated reconstruction candidate
   from exact `(t_n,y_n)`, not as part of or authority over the main trajectory.
7. **Exact oracle hierarchy?** rk8pd Oracle-1: `rtol=1e-12`,
   `atol=(1e-17,1e-23,1e-23)`; Oracle-2: `rtol=1e-13`,
   `atol=(1e-18,1e-24,1e-24)`.
8. **Oracle self-convergence?** `d_O<=D_O1` and conservative
   `U_O=2max(d_O,F_O)<=0.20F_i` at every component/point, before candidate
   adjudication.
9. **Exact reconstruction budget?** Conservative candidate error including
   oracle uncertainty must be `<=B_i=0.25F_i`; the analogous ledger budget is
   `0.25F_P`; R20 propagation/integral perturbations are `<=5e-6N_R20`.
10. **All observations or subset?** All 240 positive observations.
11. **Methods compared?** Linear, cubic Hermite, two-level local RKF45 replay,
    and a derived no-new-run Hermite/replay hybrid; two-level rk8pd is the
    oracle.
12. **Is a hybrid justified?** Only as a predeclared diagnostic/contingency if
    Hermite fails exclusively at knots. A passing single method is preferred.
13. **Preferred method before results?** No scientific winner. The
    predeclared selection rule chooses only after common budgets pass.
14. **Tie-break?** Among factor-four comfortably passing single methods, choose
    linear, then Hermite, then replay by cost/complexity; otherwise choose the
    smallest passing worst utilization, with exact ties in that same order.
15. **Must clean BA12R rerun uninterrupted B/R/U?** **YES**, all three source
    tiers.
16. **Must controls be rerun?** **YES**, all three uninterrupted matched-control
    tiers.
17. **Does R20 quadrature change?** **NO**. Retain composite trapezoid on
    qualified reconstructed scientific checkpoints and its existing budgets.
18. **Can reconstruction run process-parallel?** Yes, by observation partitions
    with completely independent contexts/caches/GSL/output roots; assembly is
    serialized.
19. **Runtime scaling?** Linear/Hermite are `O(N)`; local replay/oracle are
    `O(N * adaptive stages)`. Moving from 240 to 8192 observations multiplies
    the observation count by 34.1333 and may make replay dominant.
20. **New ADR needed?** Not for test/validation. **YES before adoption as
    Phase-6 production numerical-output architecture.**
21. **Any Phase-5D implication?** **NO.** Shared-utility or solver changes would
    require separate governance and are not proposed.
22. **Exact owner authorization required next?** Approval of this exact
    240-point, three-candidate/two-level-RKF45/two-level-rk8pd validation plan,
    its tolerances, initial-step rule, budgets, isolated-context semantics,
    category counts, deterministic repeat, stop conditions, and test/validation-
    only implementation. That authorization must explicitly permit the 1434
    local replay integrations (239 observations times six executions: the two
    RKF45 levels, their deterministic repeats, and the two rk8pd levels),
    while authorizing zero new main trajectories, zero BA12R runs, and no
    production adoption.

## 21. Owner decision requested

The owner is asked to accept or reject this exact bounded reconstruction
experiment. Acceptance would authorize a future Phase-6 test/validation-only
harness and exactly 1434 local one-accepted-interval replay/oracle integrations:
239 each for Replay-1, Replay-2, their two fresh-process deterministic repeats,
Oracle-1, and Oracle-2. Linear, Hermite, diagnostic evaluation, comparison,
and their deterministic repeats are not ODE integrations. The one accepted-
endpoint observation is returned exactly and requires no replay.

Acceptance would **not** authorize a new main source/control trajectory, a
full BA12R retry, `ScaledRKF45` or Cstar changes, Phase-5D changes, a production
reconstruction owner, an ADR, a BNV candidate, a physical BNV model/rate, or a
merge. The future task must return its evidence and method disposition to the
owner before any production adoption or clean six-trajectory BA12R campaign.
