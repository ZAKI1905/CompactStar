# Phase-6A-1 passive-observation experiment

Status before execution: **PREDECLARED / NOT YET EXECUTED**

Classification: **PHASE-6 NUMERICAL-ARCHITECTURE VALIDATION EVIDENCE;
NOT BNV CANDIDATE; NOT GOVERNED BASELINE; NOT PHYSICAL RESULT**

## Immutable pre-run declaration

This section freezes the experiment before any new ODE integration. It must
not be edited after execution. The result record will be appended separately.

### Authority and authenticated entry

The human owner accepted the segmentation diagnostic at
`6057eb92339e5a0596baab6a652c6290d0930658` and authorized exactly one
bounded Phase-6 P2 source integration to test a scheduling-only passive
observer. The authorization excludes checkpoint-state reconstruction and does
not reopen BNV thermodynamics.

The experiment branch is
`analysis/phase6a1-passive-observation-proof` in worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-passive-observation-proof`,
created from exact SHA `6057eb92339e5a0596baab6a652c6290d0930658`.
At entry, the segmentation branch local/upstream/live refs were clean and
equal to that SHA. Canonical local `master`, `origin/master`, and live
`refs/heads/master` were clean and equal to
`bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`. Failed implementation
`e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` was not an ancestor of either
the experiment entry or canonical master.

This is a numerical-method validation experiment under `GOVERNANCE.md` section
2. It adds only Phase-6 test/validation code and this record. ADR-0016's
Phase-6 ownership boundary and the recovered Phase-5B tangent authority remain
unchanged. `ScaledRKF45`, all Phase-5 production code, Phase-5D, the Cstar
cache, governed baselines, EOS/data, and literature are protected.

Historical status is immutable:

- BA12: **FAIL**.
- BA12R: **FAIL**.
- Segmentation diagnostic: **CONFIRMED NUMERICAL EFFECT**.

### Authenticated Arm E reference

The ignored diagnostic artifacts were authenticated in place against
`docs/validation/PHASE6A1_BA12R_SEGMENTATION_DIAGNOSTIC.md`:

| Artifact | Required and authenticated SHA-256 |
| --- | --- |
| Arm E trajectory | `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c` |
| Arm E checkpoint steps | `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb` |
| Arm E internal steps | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` |

Arm E ended at `t = 462269531250 s` with exact max-digits-10 values:

| Component | Reference value |
| --- | ---: |
| `x_state` | `0.49240008824076903` |
| `eta_npe` / `eta_e` (MeV) | `-2.5123474256442210e-7` |
| `eta_npmu` / `eta_mu` (MeV) | `-4.7906773046561003e-7` |
| accepted steps | `232` |
| rejected steps | `60` |

The Arm E normalized internal-step representation contains 232 accepted-step
rows with ordinal, final ceiling, previous and accepted times, accepted step,
suggested next `h`, cumulative rejections, thermal state before/after, and
Cstar-cell crossings. Exact byte identity to that representation is the
predeclared accepted-history gate. The archived representation does not carry
the two chemical components at every internal step; their final endpoint
values are independently gated bit-for-bit by the Arm E trajectory.

Arm E used maximum frozen utilization `0.007921712322551693` and maximum
`|DeltaB|/B0 = 1.46484375124617e-08` on the bounded interval.

### Authenticated observation schedule

The exact Arm S endpoint list has 241 rows, indices `0 ... 240`, with 240
positive requested observations. Its SHA-256 is
`43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e`.
Index 0 is `t = 0`; index 240 is `t = 462269531250 s`. No schedule value may
be regenerated, reordered, or numerically altered.

### Exact run configuration

| Property | Predeclared value |
| --- | --- |
| source/card | source ON; `CPL-P2-LINEAR-QSS-v1` |
| source drive | `Bdot/B0 = -1.0e-12 yr^-1`; `Bdot = -2.4136520263641375e37 count/s` |
| partition | P2, full retention |
| spin | OFF |
| processes | Me/Mmu enabled; De/Dmu disabled |
| initial state | `(x_state,eta_npe,eta_npmu) = (0,0,0)`; `Tinf = 1e8 K` |
| start/final time | `0 s`; `462269531250 s` |
| solver | GSL RKF45 |
| relative tolerance | `1e-11` |
| absolute tolerances | `(1e-16,1e-22,1e-22)` |
| initial `h` | exactly `1 s`, assigned once |
| GSL state lifetime | one persistent stepper/control/evolve/`h` from start to finish |
| positive-time GSL `t1` values | exactly one distinct value: `462269531250 s` |
| observer input | immutable requested schedule plus accepted-step ordinal and bracketing times only |
| new source integrations | exactly one |

The executable will use the same local Mac platform, Apple clang toolchain,
Debug flags, GSL 2.7.1 libraries, authenticated Phase-6 input bytes, recovered
production static library, and linked dependency libraries used by the Arm E
diagnostic. A pre-run machine-readable proof must authenticate these identities
and the reference/schedule hashes. Any difference other than attaching the
passive schedule is a stop condition before integration.

### Passive observer semantics

The observation schedule is metadata. After each accepted adaptive step, and
only then, the validation adapter may notify the observer with immutable
`(t_previous, t_new, accepted_step_index)`. For every requested positive time
satisfying `t_previous < t_obs <= t_new`, the observer records only:

- observation index;
- requested `t_obs`;
- previous and new accepted-step indices;
- `t_previous`; and
- `t_new`.

Index 0 is recorded as the initial observation at `t = 0`. The observer must
not receive a state vector. It must not interpolate, evaluate the RHS, replay
the RHS, reintegrate, use dense output, assign a nearest state, or reconstruct
`y(t_obs)` in any way.

The integration loop may pass only the final endpoint to
`gsl_odeiv2_evolve_apply` as `t1`. It will audit every advance-call target,
report both the call count and the set of distinct positive `t1` values, and
refuse any intermediate schedule value used as a target.

### Frozen qualification criteria

PASS requires every gate P1-P6. There is no tolerance and no post-result
repair, retuning, code change, or retry.

1. **P1 — final state exact.** At `462269531250 s`, all three binary64 state
   components equal Arm E exactly; unequal component count is zero.
2. **P2 — step counts exact.** Accepted/rejected counts are exactly `232/60`.
3. **P3 — accepted-step history exact.** All 232 archived internal-step rows
   compare exactly in every deterministic recorded field, including binary64
   fields. The preferred whole-file hash is exactly
   `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8`.
4. **P4 — source/currentness/validity exact.** The trajectory's source,
   moving-reference currentness, tangent currentness, frozen validity,
   process selection, domain/star/partition identities, maximum frozen
   utilization, and maximum `|DeltaB|/B0` equal Arm E exactly. Exact whole-file
   trajectory identity is preferred.
5. **P5 — schedule coverage exact.** All 241 indices are present and ordered;
   every positive observation is assigned exactly once; missing, duplicate,
   and out-of-order counts are zero; final index 240 may equal the final
   accepted endpoint.
6. **P6 — one GSL ceiling.** There is exactly one distinct positive GSL `t1`
   value, `462269531250 s`; none of the 239 intermediate positive requested
   observations appears as a target. Source inspection plus runtime target
   audit is required.

The performance record will include wall/user/system CPU time, accepted and
rejected counts, observer callbacks, requested observations, positive-time
bracket assignments, and minimum/median/maximum accepted step.

### Failure and interpretation boundary

If any P1-P6 gate fails, execution stops without repair or rerun and the result
returns to the owner. If all pass, the only permitted conclusion is:

> **PASSIVE OBSERVATION ARCHITECTURE QUALIFIED FOR SCHEDULING PASSIVITY ONLY.**

This experiment does not qualify checkpoint-state reconstruction. Linear or
Hermite interpolation, RK dense output, RHS-based interpolation, reintegration,
and nearest-step assignment are explicitly excluded. A future owner-approved
task must compare reconstruction candidates against an independent numerical
oracle before adoption.

## Execution result

Not yet executed. The predeclaration commit SHA, implementation, machine-readable
single-variable proof, one-run evidence, qualification result, and final
disposition will be appended after execution without changing the immutable
pre-run declaration above.
