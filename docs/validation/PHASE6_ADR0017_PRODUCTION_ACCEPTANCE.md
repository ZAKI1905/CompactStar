# Phase-6 ADR-0017 production implementation owner acceptance

## Status and authority

**Status:** OWNER-ACCEPTED / CANONICAL INTEGRATION AUTHORIZED.

**Clean BA12R campaign:** NOT YET AUTHORIZED.

On 2026-09-25, the human owner explicitly accepted the already-qualified
ADR-0017 production implementation at
`47a24317c604d73a7a8eece720823231036c765f`.  This record consumes the bounded
qualification evidence already produced.  It performs no ODE integration,
rk8pd reconstruction, BA12/BA12R rerun, source/control trajectory, numerical
candidate creation, or cluster work.

| Authority | SHA or identity |
|---|---|
| canonical entry | `ffd597e00162fa4efc477fe0d23878d7ad05ff6c` |
| qualified branch | `physics/phase6-adr0017-production-implementation` |
| production commit | `dee330df5b8231ff55beb4adb75824272056db21` |
| focused test/qualification commit | `a0679f01e2ffe259a4a03081721a892b2cf8bee4` |
| qualified implementation and validation commit | `47a24317c604d73a7a8eece720823231036c765f` |
| acceptance commit | the commit containing this record, with subject `docs: accept adr0017 production implementation` |

The canonical entry is the exact merge base and exact ancestor of the
qualified implementation.  The qualified history is a single-parent chain
from canonical entry through predeclaration, production, focused
test/qualification, and qualification-result commits.  It contains no merge
commit and no ancestry from the noncanonical checkpoint-reconstruction,
passive-observation, BA12R, failed controlled-BNV, or recovery branches.

## Accepted production architecture

The owner accepts the following architecture without modification:

1. `UninterruptedBnvTrajectory` is the authoritative uninterrupted adaptive
   RKF45 main evolution.
2. `PassiveObservationSchedule` is metadata-only and has no integration
   authority.
3. `Rk8pdCheckpointReconstructor` owns strict-interior checkpoint states.
4. An observation coincident with an accepted main endpoint uses that exact
   accepted state without replay, interpolation, or recomputation.
5. Strict-interior reconstruction uses isolated two-level local GSL rk8pd.
   The witness uses `rtol=1e-12` and
   `atol=(1e-17,1e-23,1e-23)`.  The reported state uses `rtol=1e-13` and
   `atol=(1e-18,1e-24,1e-24)`.  Componentwise qualification requires
   `d_O=|y_O2-y_O1|`, `d_O<=D_O1`, and
   `U_O=2 max(d_O,F_O)<=0.20 F_i`.  Failure is fail-closed.  No fallback,
   local retuning, averaging, or knot-specific reconstruction is authorized.
6. The same method applies on smooth and Cstar-knot intervals.
7. Diagnostics are evaluated from the qualified reconstructed state through
   governed Phase-6 diagnostic owners.  Powers and luminosities are not
   separately interpolated.
8. R20 retains the existing scientific observation grid and composite-
   trapezoid semantics, with reconstruction uncertainty propagated
   separately.

## Accepted qualification evidence

The owner accepts the bounded qualification record in
`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md`, including:

| Evidence | Accepted result |
|---|---|
| final Arm-E state | `x=0.49240008824076903`, `eta_e=-2.5123474256442210e-7`, `eta_mu=-4.7906773046561003e-7`; exact |
| main accepted / rejected steps | `232 / 60`; exact |
| internal accepted-step SHA-256 | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8`; exact |
| positive GSL `t1` targets | one: `462269531250 s` |
| strict-interior checkpoints | `239/239` self-qualified; zero unresolved |
| Cstar-knot checkpoints | all three self-qualified under the same rule |
| production O1/O2 values | all `478` bit-identical to validated numerical authority |
| reconstructed diagnostic packets | all `241` bit-identical |
| R20 | residual `5.7186274768377777e39 erg`; normalizer `1.0940924194731047e46 erg`; normalized `5.2268230499136139e-7` |
| reconstruction uncertainty | `7.0988433612780846e31 erg`, separate from R20 |

The retained authentication and result records report PASS.  The retained
checkpoint and performance artifacts have SHA-256
`a182f07affea4ff38aff827ac7000d3ae831a3f8dcb4f827fffc51b07c815171`
and `1274ad3a2cdb868b47611d9333ffb249a1f1da082bb051abd42f0dbcb6a46fa5`,
respectively.  This acceptance does not rerun the qualification.

## Owner-authorized path correction

The owner authorizes exactly one correction inside the frozen predeclaration
text:

```text
docs/CURRENT_ARCHITECTURE.md
```

becomes:

```text
docs/architecture/CURRENT_ARCHITECTURE.md
```

The correction occurs after qualification and only in this new descendant.
It is explicitly human-owner-authorized, typographical, path-only, and
documentation-only.  It changes no equation, numerical literal, tolerance,
qualification gate, run configuration, expected hash, scope exclusion,
implementation behavior, evidence, or result.  The historical predeclaration
commit `0c25df48a871403e716f1926390d2fec41c15c05` remains unchanged.

## Protection and non-change result

Read-only reauthentication before acceptance recorded:

- `11/11` governed baselines unchanged;
- `33/33` protected Phase-5D paths unchanged;
- `15/15` tracked Phase-5B/C/D, ADR, and imported-evidence authorities
  unchanged;
- ADR-0015 SHA-256
  `795c8bc851644de00c8c2ae37ea40c13c5a061bee691ef855d78fcd2096e1da5`
  unchanged;
- ADR-0016 SHA-256
  `b2b0bba07ebb9d926473e0123e3b28a777656368a4d642f8f5172c85711412d5`
  unchanged;
- `CompactStar/Physics/Rotochemical/ScaledRKF45.hpp` SHA-256
  `27bd510b99cfe67370be4c2f354840b81eb94b5cf773627d1ed4f4cc328ef3d0`
  unchanged;
- no tracked EOS/data, baseline, or literature path changed from canonical
  entry through the qualified implementation;
- all `22/22` entries in `literature/SHA256SUMS.txt` pass; and
- the acceptance delta from the qualified implementation is documentation-
  only.  Production and test bytes remain exactly those at
  `47a24317c604d73a7a8eece720823231036c765f`.

Phase-5D, `ScaledRKF45`, Cstar, ADR-0015, ADR-0016, governed Phase-5
baselines, EOS/data inputs, literature, production source, test source, CMake
files, and numerical evidence/results are unchanged by acceptance.

## Scientific status and exclusions

- ADR-0015: **ACCEPTED**.
- ADR-0016: **ACCEPTED**.
- ADR-0017: **ACCEPTED / HUMAN-RATIFIED**.
- ADR-0017 production implementation: **QUALIFIED / OWNER-ACCEPTED /
  CANONICAL INTEGRATION AUTHORIZED**.
- Phase-6 main integration: **uninterrupted RKF45**.
- Phase-6 observation: **passive**.
- Phase-6 strict-interior reconstruction: **self-qualified two-level isolated
  rk8pd**.
- Historical BA12: **FAIL**, permanently.
- Historical BA12R: **FAIL**, permanently.
- Clean BA12R campaign: **NOT YET AUTHORIZED**.
- EKU cluster numerical platform: **NOT YET QUALIFIED**.
- Canonical Phase-6 numerical BNV candidate: **NONE**.
- Governed BNV numerical baseline: **NONE**.

Explicit exclusions remain: no clean BA12R rerun; no six fresh
BASELINE/REFINED/ULTRA source/control trajectories; no numerical BNV
candidate; no physical BNV rate/model; no A18; no superfluidity; no
Regime-II/MixedStar; no Phase-5D, `ScaledRKF45`, or Cstar modification; no
change to ADR-0015 or ADR-0016; and no governed Phase-5 baseline modification.
No fallback, local retuning, averaging, or knot-specific reconstruction is
authorized.

## Acceptance-task execution accounting

| Operation | Count or result |
|---|---:|
| ODE trajectories run | `0` |
| rk8pd reconstructions run | `0` |
| BA12 run | `NO` |
| BA12R run | `NO` |
| cluster jobs submitted | `0` |
| six clean trajectories | `0` |
| numerical candidate created | `NO` |

## Integration gate and next action

This record authorizes only a non-force push of the acceptance commit followed
by fast-forward-only integration to canonical `master`.  No merge commit,
squash, cherry-pick, rebase, or force push is authorized.  Canonical entry must
remain the exact ancestor, both worktrees must remain clean at their gates, and
local/upstream/live parity must be demonstrated after each push.

After successful canonical integration, the exact next action is to open a
fresh **DOCUMENTATION / PLANNING** task from canonical master for **EKU CLUSTER
NUMERICAL-PLATFORM QUALIFICATION**.  That planning task must define repository
and bootstrap strategy, exact compiler/CMake/GSL/toolchain capture, Slurm
resources, process-level parallelism, isolated scratch/output roots,
same-code/source authentication, minimum Phase-5 regression qualification,
bounded ADR-0017 cross-platform qualification, comparison tolerances,
reproducibility/provenance, a future job-array strategy, stop conditions, and
the criteria required before the six-run clean BA12R campaign may be
authorized.  It must not submit cluster jobs automatically.
