# Phase-6A-1 checkpoint reconstruction validation

Classification: **PHASE-6 CHECKPOINT-RECONSTRUCTION NUMERICAL VALIDATION
EVIDENCE; NOT A BNV CANDIDATE; NOT A GOVERNED BASELINE; NOT A PHYSICAL
RESULT**.

Historical status is unchanged:

- BA12: **FAIL**;
- BA12R: **FAIL**;
- segmentation diagnostic: **CONFIRMED NUMERICAL EFFECT**; and
- passive observation scheduling: **PASS**.

## Immutable pre-run declaration

Status: **PREDECLARED / NO LOCAL INTERVAL INTEGRATION EXECUTED**.

This section freezes the bounded validation experiment before implementation
or any local interval solve. It must not be amended after the predeclaration
commit. Execution evidence will be appended in a later commit.

### Authority and authenticated entry

The human owner accepted the checkpoint-reconstruction preflight at
`93e93c7f91a3cd8fced2f7a0961eda9c469c43fe` and authorized the exact bounded
candidate/oracle experiment declared there. The validation branch is
`analysis/phase6a1-checkpoint-reconstruction-validation` in worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-checkpoint-reconstruction-validation`,
created from that exact SHA. Its parent passive-observation evidence is
`9af8912107ea77a2e2ea51517c1776f26a4b7b49`.

At entry, the accepted preflight local/upstream/live refs were clean and equal
to `93e93c7f91a3cd8fced2f7a0961eda9c469c43fe`. Canonical local `master`,
`origin/master`, and live `refs/heads/master` were clean and equal to
`bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`.

This is a numerical-method validation change under `GOVERNANCE.md:43-57`.
The accepted preflight requires a Phase-6 test/validation-only implementation
and forbids a new main trajectory, production adoption, BA12R rerun, physical
model/rate, candidate, ADR, or merge
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:580-609`).

### Immutable passive inputs

No main ODE evolution is authorized. The existing uninterrupted passive
trajectory is immutable input. The following bytes reauthenticated before
this declaration:

| Artifact | SHA-256 |
| --- | --- |
| 241-row observation schedule | `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e` |
| 232 accepted endpoint states | `7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459` |
| 241 observation brackets | `0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d` |
| 232 internal accepted steps | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` |
| passive trajectory | `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c` |
| passive step summary | `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb` |
| authenticated production library | `b9b767dbc0114563e1d556e296b6d7fc9d680a9d90e8deae9b44357010dd6499` |

These are the accepted hashes in the design and passive records
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:110-128`,
`docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md:291-305`).
Before and after execution, every immutable passive input above must retain its
exact hash.

### Solve matrix and exact authority count

The complete machine-readable authority is
`docs/validation/phase6a1_checkpoint_reconstruction_solve_matrix.tsv`, SHA-256
`32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616`.
Its 240 rows record every positive observation identity, exact left/right
endpoint indices, times and three-component states, Cstar category, deep and
endpoint flags, and every applicable method/tier solve.

The matrix derives the accepted count rather than forcing it:

| Required local integration | Strict-interior observations | Count |
| --- | ---: | ---: |
| Oracle-1 rk8pd | 239 | 239 |
| Oracle-2 rk8pd | 239 | 239 |
| Replay-1 RKF45 | 239 | 239 |
| Replay-2 RKF45 | 239 | 239 |
| fresh-process Replay-1 deterministic repeat | 239 | 239 |
| fresh-process Replay-2 deterministic repeat | 239 | 239 |
| **TOTAL AUTHORIZED LOCAL INTEGRATIONS** |  | **1434** |

Linear, Hermite, their deterministic repeats, endpoint RHS calls, diagnostic
evaluation, comparison, and the derived hybrid are not ODE integrations. The
single exact-endpoint observation returns its stored state and consumes zero
local integrations. This is the decomposition accepted at
`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:580-605` and
`:699-708`.

If the tracked matrix hash, row count, category counts, or total differs at
execution, disposition is **SOLVE-MATRIX / AUTHORITY MISMATCH** and no local
integration may run.

### Frozen observation strata

The tracked matrix reproduces the exact accepted strata
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:381-409`):

| Stratum | Count |
| --- | ---: |
| positive observations | 240 |
| no Cstar knot (A) | 237 |
| one Cstar knot (B) | 3 |
| multiple knots (C) | 0 |
| exact accepted endpoint (D, overlapping) | 1 |
| strict interior | 239 |
| deep inside a long step (E, overlapping) | 81 |

Deep means strict interior, accepted interval
`H >= 2212334854.861145 s`, and `0.25 <= s <= 0.75`. No stratum may be
redefined after results.

### Candidate definitions

For an interval with `H=t_R-t_L` and `s=(t_obs-t_L)/H`:

```text
Linear:
  y_LIN = (1-s)y_L + s y_R.

Hermite:
  h00 =  2s^3 - 3s^2 + 1
  h10 =    s^3 - 2s^2 + s
  h01 = -2s^3 + 3s^2
  h11 =    s^3 -   s^2
  y_HER = h00 y_L + h10 H f_L + h01 y_R + h11 H f_R.
```

`f_L` and `f_R` are evaluated post-run in isolated disposable contexts.
Repeated evaluation of the same endpoint must be bit-identical. Raw Hermite
is used across Cstar knots without correction. Exact endpoints return stored
states for every method.

Local replay starts from exact `(t_L,y_L)`, advances only to `t_obs`, uses
`h0=t_obs-t_L`, and stops on 100000 GSL advance calls, any GSL/nonfinite
failure, or any currentness/validity refusal. Its fixed hierarchy is:

| Tier | GSL method | `rtol` | `atol(x,eta_e,eta_mu)` | Role |
| --- | --- | ---: | --- | --- |
| Replay-1 | RKF45 | `1e-11` | `(1e-16,1e-22,1e-22)` | witness |
| Replay-2 | RKF45 | `1e-12` | `(1e-17,1e-23,1e-23)` | candidate |

Replay additionally requires `abs(y_Replay2-y_Replay1)<=F_i` everywhere.
Definitions and fixed initial-step semantics are governed by
`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:223-257`.

### Independent oracle and self-qualification

Both independent local oracle levels use fresh GSL objects and disposable
contexts, start from exact `(t_L,y_L)`, use only `t_obs` as the ceiling, and
set `h0=t_obs-t_L`:

| Tier | GSL method | `rtol` | `atol(x,eta_e,eta_mu)` |
| --- | --- | ---: | --- |
| Oracle-1 | rk8pd | `1e-12` | `(1e-17,1e-23,1e-23)` |
| Oracle-2 | rk8pd | `1e-13` | `(1e-18,1e-24,1e-24)` |

For state component `i` and observation `j`, with exact max-digits-10
binary64 round trips so `Q_O=Q_i=0`:

```text
d_O = abs(y_Oracle2-y_Oracle1)
M_O = max(abs(y_Oracle1),abs(y_Oracle2))
D_O1 = atol_Oracle1 + rtol_Oracle1 M_O
D_O2 = atol_Oracle2 + rtol_Oracle2 M_O
F_O = max(D_O2,64 ulp(M_O),Q_O)
U_O = 2 max(d_O,F_O).
```

The oracle self-qualifies only if, without averaging or exception,

```text
d_O <= D_O1
U_O <= 0.20 F_i
```

for every component and observation. The exact endpoint has both oracle
values equal to the stored endpoint and `d_O=U_O=0`. Any failure stops the
experiment before pending replay solves and adjudicates no candidate
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:259-322`).

### State, ledger and R20 budgets

For every candidate state, define a candidate-independent magnitude and the
accepted BA12R/ULTRA scale:

```text
M_i = max(abs(y_L),abs(y_R),abs(y_Oracle1),abs(y_Oracle2))
D_U = atol_ULTRA + rtol_ULTRA M_i
F_i = max(D_U,64 ulp(M_i),Q_i)
B_i = 0.25 F_i
E_C = abs(y_C-y_Oracle2)+U_O.
```

Here `rtol_ULTRA=1e-11`,
`atol_ULTRA=(1e-16,1e-22,1e-22)`, and `Q_i=0`. A candidate passes a
state only if `E_C<=B_i`; all components and all 240 observations must pass.

State is reconstructed first. Every diagnostic is then evaluated through the
same governed Phase-6 owners in a disposable context; no power, rate,
potential, validity quantity, or temperature is interpolated independently.
For every accepted ledger observable:

```text
G_P = max(1 erg/s,
          abs(P_dir_actual)+abs(L_H)+abs(DeltaLnu)+abs(Lnu_eq)
          +abs(Lgamma)+abs(Lother)+abs(L_out_fluid))
F_P = max(1e-11 G_P,64 ulp(M_P),Q_P)
B_P = 0.25 F_P.
```

The diagnostic oracle uncertainty uses the same two-level construction,
must satisfy `U_O,diagnostic<=0.20F_P`, and the candidate must satisfy
`abs(O_C-O_Oracle2)+U_O,diagnostic<=B_P` for each separated BA12R ledger
observable. `Q_P=0` for exact max-digits-10 binary64 serialization.

R20 retains the requested checkpoint grid and composite trapezoid. For each
candidate, the candidate-versus-Oracle-2 luminosity-integral difference and
thermal/end-state propagation difference must each be
`<=5e-6 N_R20`. The governed `abs(R20)/N_R20<=2e-4` criterion is unchanged
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:324-379`).

### Hybrid and selection policy

The derived hybrid uses Hermite on A and Replay-2 on B/C, without additional
integration. It is eligible only if Hermite passes all A observations, fails
only B/C observations, Replay-2 passes every replacement point, and the
combined method passes every common state, ledger, R20, deterministic,
identity, currentness, and validity gate. Any passing single method is
preferred over the hybrid.

A method's worst utilization is the maximum of every state, ledger and R20
error divided by its corresponding budget. A single method with worst
utilization `<=0.25` is comfortably passing. If one or more singles are
comfortably passing, select the least costly in the frozen order linear,
Hermite, Replay-2. Otherwise select the passing single method with the
smallest worst utilization; exact binary64 ties use that same order. If no
single passes, consider only the strictly eligible hybrid. No threshold or
tie-break may change after results
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:412-455`).

### Isolation, determinism and parallel execution

No operation may use or mutate a live/main trajectory object. Every endpoint
RHS evaluation, local replay, oracle solve, and diagnostic evaluation runs in
an isolated disposable process/context with private mutable state,
accumulators, caches, GSL objects, and output root. Every solve starts from the
authenticated immutable left endpoint. The main passive inputs are read-only.

The fixed maximum concurrency is **2 local worker processes**. Jobs are
partitioned deterministically by solve ID. OpenMP/thread environment is not
changed from the authenticated platform. If process isolation or deterministic
serial assembly cannot be shown, execution falls back to one process without
changing scientific configuration. No scientific multithreading is added.

Oracle-1/Oracle-2 run first. Only after complete oracle self-qualification may
pending replay jobs execute. The atomic execution ledger records solve ID,
observation, method/tier, left endpoint, target, process ID, start/end time,
exit status, and output hash. Completed/attempted integrations may never exceed
1434. An external interruption retry counts against that cap.

### Stop conditions

Execution stops and returns to the owner if:

- an entry, passive artifact, production-library, schedule, or solve-matrix
  identity differs;
- row/category counts or the authoritative total are not exactly frozen;
- any new main trajectory is required or an observation becomes a main GSL
  ceiling;
- the passive main evidence changes;
- an oracle state or diagnostic fails self-qualification;
- a candidate would require retuning, a new method, or a changed budget;
- any reconstruction mutates another method's or the passive state/cache;
- repeated endpoint RHS or required replay output is nondeterministic;
- Cstar knot handling is ambiguous or differs from the matrix;
- governed diagnostic values cannot be produced from the reconstructed state;
- R20 requires a different quadrature owner;
- any production, `ScaledRKF45`, Cstar, Phase-5D, baseline, EOS/data, or
  literature change is required;
- historical BA12/BA12R would need reinterpretation; or
- a physical BNV rate/model, numerical candidate, production ADR, cluster job,
  full-suite run, or merge would be required.

No post-result repair or hidden retry is authorized.

### Pre-run scope statement

As of this declaration, local interval integrations executed: **0**. Main
trajectory integrations executed: **0**. No interpolation, endpoint RHS,
candidate diagnostic, or oracle value has been evaluated. No production,
test, CMake, baseline, EOS/data, or literature byte has yet changed on this
validation branch. The only new files are this declaration and its tracked
solve matrix.

## Execution result

**Disposition:** SIDE-EFFECT / CURRENTNESS / PROVENANCE FAILURE — RETURN TO
OWNER.

**Classification:** PHASE-6 CHECKPOINT-RECONSTRUCTION NUMERICAL VALIDATION
EVIDENCE — NOT BNV CANDIDATE — NOT GOVERNED BASELINE — NOT PHYSICAL RESULT.

The frozen pre-run declaration was committed as
`94788a2f0941eb9f9510ffb810c7f89b9e3ae26d`. The test-only validation harness
was then committed as `6aefb3833177b54d1165782f6fe7347d73ca4d1e` at:

- `tests/bnv/checkpoint_reconstruction_validate.cpp`;
- `tests/bnv/checkpoint_reconstruction_verify.py`.

The harness was compiled directly against the authenticated Phase-6 recovery
library; no CMake or production source was changed. Its executable SHA-256 was
`82eea47fba62f1440a1b87abec95dc4ac6287756ad44b990564794b9cf703044`.
The exact solve-matrix SHA-256 remained
`32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616`.

### Oracle execution and qualification

Oracle-1 and Oracle-2 ran first with the frozen configurations and maximum
process concurrency 2. Every one of the 239 strict-interior observations was
integrated independently at both tiers. The exact endpoint consumed no solve.
The two atomic solve ledgers therefore contain **478 unique successful local
integrations**. No main trajectory was run.

The two-level oracle passed every state, diagnostic, identity and currentness
self-qualification gate. The comparison artifact
`oracle-result.json` has SHA-256
`f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0`.
The exact maxima were:

- `max(d_O/D_O1) = 0.5794839113173227`, at observation 117,
  `x_state`, category B (one Cstar knot), not deep interior;
- `max(U_O/(0.20 F_i)) = 0.5793505315921852`, at observation 117,
  `x_state`, category B, not deep interior;
- maximum diagnostic oracle utilization `= 0.17336139714051704`, at
  observation 117 for `Pnet`, category B, not deep interior.

Oracle-1 accumulated 246 accepted and 13 rejected local steps, 3606 RHS
evaluations, `0.8088575 s` summed local-solve wall time,
`0.8085089999998445 s` user CPU and `0.00014300000000044832 s` system CPU.
Its median/p95/maximum solve wall times were respectively
`0.003145417 s`, `0.0031900920000000003 s`, and `0.026010459 s`.

Oracle-2 accumulated 249 accepted and 25 rejected local steps, 3801 RHS
evaluations, `0.847701083 s` summed local-solve wall time,
`0.84699500000022 s` user CPU and `0.0003089999999985604 s` system CPU.
Its median/p95/maximum solve wall times were respectively
`0.003123334 s`, `0.003195542 s`, and `0.045704834 s`.

A solve-only linear extrapolation to 8191 strict-interior points is
`27.721137165271966 s` for Oracle-1 and `29.052383141644352 s` for Oracle-2
before context-construction, serialization, scheduling, and diagnostic
overheads. This extrapolation is operational only and is not scientific
evidence.

### Candidate-phase stop

After oracle qualification, an initial launcher authentication attempt used a
nonqualifying historical evidence path. It stopped before `RunPhase`, created
no candidate output/work root, and executed zero integrations. The intended
immutable historical PRETRAJECTORY record was then authenticated at failed
implementation SHA `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c`, with SHA-256
`1d4e77780fde474f24782879d6ac43fb5bd6d01f3919c359437c0c5f01d327a6`.

The candidate launcher subsequently passed that gate, created isolated
`linear` and `linear-repeat` scratch roots, and began their context
construction. It then failed before producing any observation result because
the supplied profile argument was a trajectory TSV rather than the required
EOS/profile directory. Both child processes failed during EOS import; no
candidate solve ledger or result row exists, and **zero candidate local
integrations** were executed. The launcher reported method failure and the
committed harness makes such a method-batch failure non-retriable. Consistent
with the frozen no-repair/no-hidden-retry rule, the candidate phase was not
restarted and no path, code, tolerance, budget, or method was changed after
the failure.

Consequently linear, Hermite, Replay-1, Replay-2 and the conditional hybrid
were **not adjudicated**. There are no candidate state utilizations, ledger
utilizations, R20 reconstruction utilizations, deterministic endpoint-RHS
result, or selected method. The remaining **956** authorized replay
integrations were skipped. Total executed local integrations are **478**,
within the hard maximum 1434.

### Immutability and status

The passive artifacts remained byte-identical after execution:

- trajectory: `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c`;
- accepted endpoints: `7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459`;
- observation brackets: `0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d`;
- internal steps: `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8`;
- step summary: `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb`.

No production, `ScaledRKF45`, Cstar, Phase-5D, baseline, EOS/data, literature,
or passive-main byte changed. Main trajectory integrations remained zero.
Historical BA12 remains **FAIL**, BA12R remains **FAIL**, and passive
scheduling remains **PASS**. No BNV candidate, physical rate/model,
production ADR, full-suite run, cluster job, or merge was created or performed.
