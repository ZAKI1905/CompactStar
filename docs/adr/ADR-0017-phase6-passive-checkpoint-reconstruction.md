# ADR-0017: Phase-6 Passive Observation and Self-Qualified Local Checkpoint Reconstruction

## 1. Status, scope, and authority

**Status:** PROPOSED — OWNER RATIFICATION REQUIRED.

**Date proposed:** 2026-09-24.

This is a narrow Phase-6 numerical-output architecture proposal. It establishes no
production authority unless and until the human owner explicitly ratifies this exact ADR.
It does not implement checkpoint reconstruction, run or rerun a trajectory, create a BNV
candidate, qualify a production method, or merge any noncanonical validation ancestry.

The change class is **numerical-method + structural/architecture + documentation**. A
production reconstruction method would change numerical output and create a new Phase-6
checkpoint-output owner; governance therefore requires numerical rationale, an ADR, and a
same-change current-architecture description (`GOVERNANCE.md:36-57`). Accepted ADR-0015 and
ADR-0016 remain higher authority. This proposal changes neither one's physics or ownership
decision (`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:1-21`,
`docs/adr/ADR-0016-phase6-bnv-tangent-adapter-ownership.md:1-28`).

The canonical entry is
`bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`. The latest reconstruction-recovery evidence
source is the noncanonical branch head
`ae0d5607fe8889cec70604747dfe3bfea7e5fd97`. The governance branch was created directly
from the canonical entry; no failed or noncanonical implementation or validation head is an
ancestor of this proposal.

## 2. Clean-ancestry evidence import

Only durable Markdown validation records were imported. Each was copied from the exact
listed source commit; source and imported SHA-256 are equal. No executable source, test
harness, build file, solve matrix, scratch data, raw trajectory, generated artifact, or
noncanonical executable change was imported.

| Imported record | Source branch and SHA | Source SHA-256 | Imported SHA-256 |
| --- | --- | --- | --- |
| `docs/validation/PHASE6A1_BA12R_NUMERICAL_FORENSICS.md` | `analysis/phase6a1-ba12r-numerical-forensics` at `0aab1c2b5b748585326de35499d1830d0c788fae` | `daccc117237c6789c4f3fabc0d6cb0252fc88ce020d67eb3b8bba4eb33577cea` | `daccc117237c6789c4f3fabc0d6cb0252fc88ce020d67eb3b8bba4eb33577cea` |
| `docs/validation/PHASE6A1_BA12R_SEGMENTATION_DIAGNOSTIC.md` | `analysis/phase6a1-ba12r-segmentation-diagnostic` at `6057eb92339e5a0596baab6a652c6290d0930658` | `ab7f6b9f3a6d37a7c2c439f9b9c0a2a330fecc208bd801ae890fe89b0851e157` | `ab7f6b9f3a6d37a7c2c439f9b9c0a2a330fecc208bd801ae890fe89b0851e157` |
| `docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md` | `analysis/phase6a1-passive-observation-proof` at `9af8912107ea77a2e2ea51517c1776f26a4b7b49` | `684a0338b7a08591ac266d8c8bd7924d95ac86c3488ab115637b4f2b7d7f76ee` | `684a0338b7a08591ac266d8c8bd7924d95ac86c3488ab115637b4f2b7d7f76ee` |
| `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md` | `analysis/phase6a1-checkpoint-reconstruction-preflight` at `93e93c7f91a3cd8fced2f7a0961eda9c469c43fe` | `188861202879f7cb7b3bc7f106cbf28113881a2296779a20f58538bd80e72778` | `188861202879f7cb7b3bc7f106cbf28113881a2296779a20f58538bd80e72778` |
| `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_VALIDATION.md` | `analysis/phase6a1-checkpoint-reconstruction-validation` at `8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead` | `59093f30ecfc0543f73f60ba416bc23a80d2ad753a9051230fdfc5646ee552a2` | `59093f30ecfc0543f73f60ba416bc23a80d2ad753a9051230fdfc5646ee552a2` |
| `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_CANDIDATE_RECOVERY.md` | `analysis/phase6a1-checkpoint-reconstruction-candidate-recovery` at `ae0d5607fe8889cec70604747dfe3bfea7e5fd97` | `403090d6ada8c796b5b65c45d430a82fc55ce331f11bb1d9d4b34f607c4a05f7` | `403090d6ada8c796b5b65c45d430a82fc55ce331f11bb1d9d4b34f607c4a05f7` |
| `docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_IMPLEMENTATION.md` | `physics/phase6a1-bnv-recovery-r0-r3` at `d2822656aae8e264f4985f0d2f7228b90320a577` | `8ab91e364abc4eb843ecfdfbb03155d4ac2c18e20b4a4caf47d31dd01c4a2a3d` | `8ab91e364abc4eb843ecfdfbb03155d4ac2c18e20b4a4caf47d31dd01c4a2a3d` |
| `docs/validation/PHASE6A1_CONTROLLED_BNV_BA12R.md` | `physics/phase6a1-bnv-recovery-ba12r` at `971b2bac47c320bd262d8686841e679f703e5ccb` | `913ebc45fac54e60df640137f850b008ea612b88e30d153daa9f1f38da5a7bef` | `913ebc45fac54e60df640137f850b008ea612b88e30d153daa9f1f38da5a7bef` |

The R0-R3 record is retained because it authenticates the recovered Phase-6 RHS/context
provenance while demonstrating exact restoration of Phase-5D authority
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_IMPLEMENTATION.md:236-259`). The primary
BA12R record is retained because it is the durable source for the historical failure and its
non-candidate disposition (`docs/validation/PHASE6A1_CONTROLLED_BNV_BA12R.md:1-25`,
`docs/validation/PHASE6A1_CONTROLLED_BNV_BA12R.md:285-333`).

## 3. Evidence disposition

The original checkpoint-reconstruction comparison remains a failed campaign:

- Linear failed 239 of 240 observations.
- Hermite failed 226 of 240 observations, including 223 no-knot observations.
- Replay-2 failed knot observations 82 and 228. Its no-knot worst state-budget utilization
  was approximately `0.096283`; its Replay-1 witness failed at knot observation 228.
- Linear, Hermite, and Replay-2 each passed their R20 reconstruction budgets.
- The hybrid was ineligible because Hermite failed away from knots and Replay-2 failed at
  knot observations.
- No method was selected.

These are the unchanged adjudicated results
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_CANDIDATE_RECOVERY.md:335-356`). The
historical disposition remains **NO CHECKPOINT RECONSTRUCTION METHOD QUALIFIED**. Nothing in
this ADR converts that campaign to PASS
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_CANDIDATE_RECOVERY.md:250-266`).

The independent two-level local rk8pd oracle did self-qualify over every required
strict-interior observation. Its exact maxima were
`max(d_O/D_O1)=0.5794839113173227` and
`max(U_O/(0.20 F_i))=0.5793505315921852`, both at observation 117 for `x_state` in the
one-Cstar-knot category
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_VALIDATION.md:320-339`). The recovery
reused and reauthenticated all 478 oracle solves without rerunning them
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_CANDIDATE_RECOVERY.md:294-314`).

The 478 oracle solves accumulated approximately `1.6566 s` of local-solve wall time. For an
8192-positive-time schedule, meaning 8191 strict-interior local solves, the solve-only linear
extrapolations were approximately `27.7 s` for Oracle-1 and `29.1 s` for Oracle-2. These
measurements exclude context construction, serialization, scheduling, and diagnostic
overheads and are not acceptance criteria
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_VALIDATION.md:341-357`). They support a
uniform high-accuracy local method instead of a more complex knot-specific hybrid.

The rk8pd result is **SELF-QUALIFIED NUMERICAL ORACLE EVIDENCE**. It is not a qualified
production method. This ADR proposes it for future production reconstruction subject to
human ratification and a separately authorized implementation/validation task.

## 4. Proposed decision: main integration and passive scheduling

If ratified, a Phase-6 controlled BNV trajectory uses one uninterrupted adaptive RKF45 main
integration from its authorized initial state to its authorized terminal time. Requested
scientific output/checkpoint times are passive metadata. They must not:

- be passed to GSL as intermediate `t1` ceilings;
- reset `h`;
- reset or replace GSL stepper, control, or evolve state; or
- alter any accepted main step.

Segmentation was shown to change accepted-step history and the final endpoint even with the
same platform, method, tolerances, initial state, and final time
(`docs/validation/PHASE6A1_BA12R_SEGMENTATION_DIAGNOSTIC.md:192-222`,
`docs/validation/PHASE6A1_BA12R_SEGMENTATION_DIAGNOSTIC.md:265-288`). The passive experiment
then obtained exact final-state, accepted/rejected-count, accepted-history, validity, schedule,
and single-GSL-ceiling PASS results
(`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:275-305`).

The Phase-6 checkpoint-output layer owns an immutable requested observation schedule. After
each accepted main step it records every requested observation bracketed by immutable accepted
endpoints `(t_L,y_L)` and `(t_R,y_R)`. Scheduling remains observational only. For fixed
physics, initial state, solver configuration, and terminal time, changing the observation
schedule must not alter:

- the final main state;
- accepted-step count;
- rejected-step count; or
- accepted-step history.

This makes the already validated passive-observer invariant mandatory
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:457-474`).

## 5. Proposed decision: checkpoint reconstruction

### 5.1 Exact accepted endpoint

If a requested observation time equals an accepted main endpoint exactly, the checkpoint
layer returns that accepted main state exactly. It performs no replay, interpolation, or
recomputation (`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:270-282`).

### 5.2 Strict-interior checkpoint

For `t_L < t_obs < t_R`, each reconstruction level is an isolated local GSL rk8pd solve. It
starts from the exact accepted left endpoint `(t_L,y_L)`, uses `t_obs` as its only local
ceiling, sets the initial local step request to `t_obs-t_L`, and terminates at `t_obs`. It
uses the same physical RHS and immutable Phase-6 context semantics. It is not an accepted
main step and must not modify the authoritative main trajectory
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:49-90`,
`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:259-282`).

### 5.3 Two-level self-qualification

The proposed witness level is the already validated Oracle-1 configuration:

```text
GSL method = rk8pd
rtol       = 1e-12
atol       = (1e-17, 1e-23, 1e-23)
```

The proposed reported-state level is the already validated Oracle-2 configuration:

```text
GSL method = rk8pd
rtol       = 1e-13
atol       = (1e-18, 1e-24, 1e-24)
```

For component `i` at an interior observation, define the accepted resolution quantities:

```text
d_O(i)  = |y_O2(i)-y_O1(i)|
M_O(i)  = max(|y_O1(i)|, |y_O2(i)|)
D_O1(i) = atol_O1,i + rtol_O1 M_O(i)
D_O2(i) = atol_O2,i + rtol_O2 M_O(i)
F_O(i)  = max(D_O2(i), 64 ulp(M_O(i)), Q_O(i))
U_O(i)  = 2 max(d_O(i), F_O(i))

M_i     = max(|y_L(i)|, |y_R(i)|, |y_O1(i)|, |y_O2(i)|)
D_U,i   = atol_ULTRA,i + rtol_ULTRA M_i
F_i     = max(D_U,i, 64 ulp(M_i), Q_i)
```

Here `rtol_ULTRA=1e-11`, `atol_ULTRA=(1e-16,1e-22,1e-22)`, and `Q_O`/`Q_i`
retain the accepted serialization-resolution meanings. Each component must satisfy:

```text
d_O <= D_O1
U_O <= 0.20 F_i.
```

The definitions and factor-two uncertainty are unchanged from the accepted reconstruction
preflight (`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:293-350`). If any
component fails either inequality, that checkpoint is numerically unresolved and the
trajectory/candidate validation fails closed under the downstream gate. No local retuning,
retry, method switch, tolerance change, or result-dependent exception is permitted.

If self-qualification passes, the reported checkpoint state is exactly:

```text
y_checkpoint = y_rk8pd_Oracle2.
```

Oracle-1 is retained only as the numerical witness. The levels are not averaged, and the
reported level is not selected after seeing results.

### 5.4 Uniform rule across Cstar knots

The same two-level local rk8pd reconstruction applies to every strict-interior checkpoint.
There is no Cstar-knot method switch, Hermite fallback, RKF45 fallback, or knot-specific
tolerance. The oracle self-qualified in the one-Cstar-knot worst category, while the tested
RKF45 replay failed only at knot observations 82 and 228
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_VALIDATION.md:328-339`,
`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_CANDIDATE_RECOVERY.md:337-354`). This is
an architectural simplification, not a claim that the current Cstar piecewise-smooth RHS has
no kinks.

## 6. Isolation and diagnostics

Every local level uses fresh disposable context and GSL objects, the immutable accepted left
state, and the same authenticated source history, frozen coefficients, domain/star identity,
and validity/currentness rules. A replay may not mutate the main trajectory, main caches, or
another checkpoint replay. This isolation is mandatory because endpoint RHS/diagnostic
evaluation mutates live context/cache objects even when its returned values are deterministic
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:130-166`).

Independent local solves may run process-parallel only when each worker owns its context, GSL
objects, cache state, and output root. No scientific thread-level sharing or parallelism is
introduced (`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:527-559`).

At `t_obs`, the output layer must:

1. reconstruct and self-qualify the state;
2. create an isolated immutable diagnostic context at that state; and
3. evaluate `P_dir`, `LH`, `DeltaLnu`, `DeltaPbeta`, `Lnu_eq`, `Lnu_full`,
   `Lgamma`, `Lother`, `Pnet`, `mu_n_actual`, `sigma`, and frozen-validity diagnostics
   through the existing governed Phase-6 owners.

No luminosity, power, potential, rate, validity metric, or temperature is independently
interpolated (`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:476-491`).

## 7. R20 and uncertainty semantics

The scientific observation grid and existing composite-trapezoid R20 quadrature are
unchanged. R20 consumes diagnostics evaluated at qualified reconstructed checkpoint states;
it is not redesigned around adaptive RK stages. Reconstruction uncertainty remains a
separately propagated numerical contribution
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:493-501`).

In every future BA12R checkpoint and tier, the local two-level rk8pd self-qualification must
pass. Reconstruction numerical uncertainty remains distinct from main-integration
convergence, R20 quadrature, finite-temperature omission, and frozen-background uncertainty.
It must not be relabeled as main-integration error.

## 8. Consequence for a future clean BA12R

The historical segmented BASELINE, REFINED, and ULTRA source trajectories must not be reused
as the future clean BA12R hierarchy. After this ADR is ratified and a separate production
implementation is accepted and validated, a clean BA12R campaign requires six new
uninterrupted trajectories:

```text
source:  BASELINE, REFINED, ULTRA
control: BASELINE, REFINED, ULTRA.
```

Every tier uses the same passive observation schedule and accepted two-level local rk8pd
checkpoint reconstruction. The controls must be rerun because matched source-minus-control
convergence requires the same uninterrupted architecture and tier
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:503-525`).

Historical segmented BA12 and BA12R remain **FAIL** forever. They retain provenance and
failure-evidence value and are never rewritten as clean passive-observation results.

## 9. Explicit non-change boundaries

This proposal is Phase-6 numerical-output architecture only.

### 9.1 Phase-5D

No Phase-5D source, baseline, producer, comparator, provenance rule, or governed trajectory
byte changes. No governed Phase-5D trajectory byte is reinterpreted. If future evidence
demonstrates a shared Phase-5 integrator defect, that is a separate upstream governance issue
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:568-578`).

### 9.2 ScaledRKF45 and Phase-5 main driver

This ADR does not authorize modifying `ScaledRKF45`, its GSL controller, RKF45 coefficients,
or the Phase-5 main driver. Any ratified implementation must place the uninterrupted-run
wrapper/output architecture in a Phase-6-owned seam unless separately authorized. Promotion
to shared numerical machinery requires separate governance and regression
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:568-578`).

### 9.3 Cstar

This ADR does not authorize changing Cstar interpolation, smoothing, knot locations, or cache
resolution. The proposed local reconstruction must tolerate the current piecewise-smooth RHS.

### 9.4 BNV science

This ADR changes no BNV thermodynamic equation. ADR-0015, the moving-reference seam, `t`,
`sigma`, `Z`, actual-potential ledger, P0/P1/P2 semantics, thermal ledger, and
frozen-background scope remain unchanged. ADR-0016 tangent-adapter ownership remains
unchanged. No physical BNV model or rate is selected.

## 10. Rejected alternatives

- **Continue segmented main integration:** rejected because observation ceilings demonstrably
  change the adaptive trajectory.
- **Adopt Linear, Hermite, Replay-2, or their hybrid:** rejected because no original candidate
  method qualified under the predeclared complete budget.
- **Use a Cstar-knot special case:** rejected because the uniform rk8pd oracle self-qualified
  through the knot category and a method switch would add result-sensitive complexity.
- **Average or select between rk8pd levels:** rejected because Oracle-1 is the fixed witness and
  Oracle-2 is the fixed reported-state level.
- **Interpolate diagnostic powers:** rejected because governed diagnostic owners must evaluate
  the qualified reconstructed state.
- **Redesign R20 around adaptive stages:** rejected as unnecessary new numerical ownership.
- **Modify shared Phase-5 machinery:** rejected as outside Phase-6 ownership and this ADR.

## 11. Ratification and future implementation gates

This ADR is **PROPOSED — OWNER RATIFICATION REQUIRED**. It provides no immediate production
authority. Before any production change, the owner must explicitly ratify the exact decision.
A later separately authorized task must then implement only the Phase-6 owner, validate every
rule above, and return for acceptance before any BASELINE/REFINED/ULTRA rerun.

Until those gates are complete:

- do not implement production checkpoint reconstruction;
- do not run or rerun BASELINE, REFINED, or ULTRA;
- do not run BA12R;
- do not create a BNV candidate;
- do not modify Phase-5D, `ScaledRKF45`, Cstar, or BNV physics; and
- do not merge noncanonical validation ancestry.

The exact next action is human-owner review and explicit ratification or rejection of this
ADR. No implementation follows automatically.
