# Stage A restart coverage — protocol change proposed, not applied

**OWNER APPROVAL REQUIRED for the exact scheduling rule below. No new trajectory
is authorized or executed by this proposal.** The active runner, runtime snapshots,
production solver and previous failure records remain unchanged.

## Decision and authority

The smallest correction is to select an **existing accepted endpoint before
midpoint**, using the existing state-restart path. Do not force a midpoint step,
change a cap, add an artificial Cstar knot, or restart from a passive packet.
The prepared [patch](proposed-selector.patch) changes only harness selection and
adds a pure selector test; it has been applied only in a disposable test directory.

This is compatible with ADR-0017, but **not already authorized by the fixed
Stage A run protocol**. That protocol explicitly says “first accepted endpoint
at or beyond half-time” in both
`../native-assembly/NONLINEAR_PLAN.md:65` and
`../native-regeneration/RUN_PLAN.md:35`. Replacing that rule changes the tested
restart initial state and interval. It is therefore a substantive experimental
protocol amendment requiring owner approval, even though its implementation is
a small test-harness edit. The owner's current instruction permits preparation
of this proposal, not unilateral substitution of a different fixed rule.

ADR-0017 §4 requires uninterrupted main integration and passive observations;
§5.1 permits exact stored accepted states; §5.2 explicitly says reconstructed
strict-interior checkpoints are not accepted main steps
(`../../../adr/ADR-0017-phase6-passive-checkpoint-reconstruction.md:108`, `:140`,
`:146`). Amendment 1 authorizes independently located **actual Cstar** boundaries,
not arbitrary clock-time breakpoints
(`../../PHASE6A1_PIECEWISE_RECONSTRUCTION_PROPOSAL.md:16`, `:33`;
`../../PHASE6_BREAKPOINT_RATIFICATION.md:25`). The proposal needs **no amendment
to ADR-0017**, physical model, accuracy threshold or integrator.

Classification: numerical-method/test-protocol (changed restart interval),
engineering (immutable-history selector and tests), documentation and generated
evidence. Main scientific values are unchanged by selection; future restart
values will be evaluated over a different interval. This remains an unratified
scientific candidate under `GOVERNANCE.md:171`; that file's §§1–3 govern authority
and classification. No pre-baseline correctness exception is invoked.

## Exact proposed rule, Stage A only

For **each role and each tolerance tier's quarter-cap main**, after obtaining its
own complete accepted history on [0,T], T=31,557,600 s:

1. Select the latest accepted endpoint satisfying **0 < t_restart <= T/2**.
2. Require a contiguous complete history and such an endpoint; otherwise stop.
   Do not substitute t=0, a reconstructed checkpoint, another tier or another role.
3. Round-trip that endpoint's exact binary64 time and three-component state using
   the existing restart transport, retaining its actual accepted ordinal, role,
   authenticated dependency identities and full original prefix evidence.
4. Use the existing fresh state/GSL restart path from t_restart to **T**, preserving
   absolute time for B(t), all tolerances, caps, physical inputs and guards.
5. Retain the original comparison against the uninterrupted reference, all prefix
   state errors and cumulative energies. Do not zero the prefix uncertainty or
   claim bitwise adaptive continuation. Missing justified uncertainty still fails.

This guarantees a tail at least **half the original horizon** whenever coverage is
accepted. With the original initial proposal of 1 s, a successful complete Stage A
main has an early accepted endpoint well before midpoint; nevertheless the
selector fails closed if the supplied history does not. The generic helper also
rejects empty, gapped, nonprogressing, truncated and out-of-horizon histories.
The bound is a coverage guarantee, **not a proof of restart accuracy**. Half-year
subsequent evolution is a nonvacuous test window; sufficiency is decided by the
unchanged numerical comparisons, not duration alone.

This replaces only the Stage A half-time selection sentence in the two run plans
cited above. Stage B retains its existing selection rule and remains unauthorized.
The patch gates the new policy on the existing one-year `Year` horizon, in the
experiment harness, and uses the original rule otherwise. It changes no generic
solver infrastructure (`tests/native/stage.cpp:107`). No extra matrix member is
introduced; the original three tails per role remain three tails.

## Root cause and actual evidence

The successful adaptive main takes a final step across half-time. The original
selector consequently chooses the terminal endpoint and its existing strict
interior guard refuses. The original caps were retained at hundreds of years,
so they do not guarantee a midyear endpoint on a one-year experiment. This is a
test-coverage mismatch, not a demonstrated inaccurate solution
(`../stage-a-recovery/REPORT.md:16`, `tests/native/stage.cpp:108`,
`CompactStar/Physics/BNV/BreakpointStepper.hpp:61`).

| Actual source BASELINE quantity | Value |
| --- | ---: |
| Accepted steps | 12 |
| First accepted endpoint | 1 s |
| Midpoint | 15,778,800 s |
| Original selected endpoint, ordinal 12 | 31,557,600 s |
| Original restart tail | **0 s — FAIL, preserved** |
| Proposed selected endpoint, ordinal 11 | 10,025,181.919443566 s |
| Proposed subsequent evolution | **21,532,418.080556434 s** |
| Tail / full horizon | **0.682321154985057** |
| Smallest planned cap | 9,680,951,939.25 s |

These values are reproduced from authenticated actual accepted records in
[result.json](evidence/result.json); that file also records exact Decimal lifts of
the binary64 restart state and time. No new RHS evaluation or integration is used.
The actual half-time observation, index 120, is strictly inside step 12, between
ordinal 11 and year-end. Its own provenance identifies it as reconstructed:
[midpoint-checkpoint.json](evidence/midpoint-checkpoint.json). Native checkpoint
qualification does not change that provenance or promote it to a main endpoint.

The quarter-cap member was **not executed**. The preserved
[inactive-cap proof](../stage-a-recovery/evidence/restart-gate.json) checks every
actual proposal and target and establishes path equivalence for all three finite
caps on this observed BASELINE history, with no rejected trials or events. The
smallest cap exceeds both the horizon and all proposals. Thus the proposed rule
would select ordinal 11 for that quarter-cap case too; this is an inference from
the authenticated proof, not a measured quarter-cap result. Other tiers and the
control must select from their own future histories; no universal equivalence is
assumed (`../stage-a-recovery/REPORT.md:44`).

## Alternatives considered

| Alternative | Consequence |
| --- | --- |
| Latest accepted endpoint at/before midpoint — proposed | No main-step change; at least half-year tail; only restart protocol changes |
| Last nonterminal endpoint with no midpoint bound | Can leave an arbitrarily short tail; does not guarantee useful temporal coverage |
| Treat passive midpoint O2 packet as accepted | Incorrect provenance; adds reconstructed-state uncertainty and changes the restart test premise |
| Force an intermediate main target or split the main at midpoint | Changes accepted main history and violates the passive schedule contract without further authorization |
| Insert artificial Cstar knot / reuse a knot event time | Not an actual cache boundary; event ownership and physical meaning are wrong |
| Reduce quarter cap to a fraction of one year | Existing cap machinery could create coverage, but changes numerical policy and main results; larger scope than selecting an existing endpoint |

These distinctions follow ADR-0017 §§4–6 and Amendment 1, and the current API
(`CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp:137`,
`CompactStar/Physics/BNV/BreakpointStepper.hpp:16`). Existing breakpoint machinery
does not expose an accepted arbitrary wall-clock restart event. It should not be
extended when the immutable accepted history already supplies a useful point.

## What restart coverage does and does not discharge

The current harness constructs fresh `RunState`, driver, EvolutionSystem and GSL
objects for a restart, retaining the immutable physical owner. It transports time,
state, role and prefix ordinals; its comparison file currently contains the tail
and uninterrupted main final states (`tests/native/stage.cpp:9`, `:110–118`;
`tests/rotochemical/coupled.hpp:27`). State restart is not process-state recovery
of GSL's adaptive history. The run plan explicitly distinguishes these contracts
(`../native-assembly/NONLINEAR_PLAN.md:31`).

Selecting a valid endpoint fixes only coverage. Dependency authentication,
currentness/cache isolation, reference comparison, cumulative prefix energies,
thermal roundoff and accumulated error remain required in the complete assessor.
The comparison file by itself does not implement these proofs or carry numeric
prefix budgets. Existing static scalar/RHS round-trip checks likewise do not
qualify a restarted trajectory (`tests/native/verify.cpp:39`). A future run must
preserve full prefix evidence and report justified error propagation, rather than
assuming it vanishes because the stored state is exact binary64.

The earlier Echem improvement was native self-qualification for 129 packets,
not independent qualification of every future checkpoint. The actual-output
100-digit reference, reference convergence, thermal enclosure assumptions,
main-prefix/thermal-roundoff accounting, R18/R20 and matched control remain
unresolved as recorded in `../stage-a-recovery/REPORT.md:149`. Nothing here
supersedes the original D1 FAIL, later 13-interval conditional PASS, or Stage A FAIL.

## Validation and provenance

Entry was clean at `d0e3c12fe0a7fe8e1cea4aa4144dec23e14e4236`, branch
`candidate/results-first-bnv-pilot-20261108`, worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All 15 local branch tips are ancestors; nine worktrees exist. All 58 compact and
1,894 external recovery manifest entries authenticated before analysis
([preflight.json](evidence/preflight.json)).

Slurm **19189** completed the standalone selector checks on phyn004 in **2 s**:
one logical execution CPU (affinity [0]), 4 GiB allocation, five-minute ceiling.
Slurm reserves two hardware threads; only one is used. GCC 14.1.0 compiled with
C++17, `-O2 -ffp-contract=off -Wall -Wextra -Werror`. Twelve synthetic boundary/
refusal tests and the actual-history selection/exact state round-trip passed.
No thermal, stellar, RHS or reference integration occurred
([tests.log](evidence/tests.log), [completion.txt](evidence/completion.txt)).
Slurm accounting storage is disabled; no measured peak RSS is claimed.

The patch applies cleanly to the starting harness. Only its standalone selector
test was compiled; the patched full stage executable was **not built or run**.
Active `tests/native/stage.cpp`, CMake, production code and snapshots remain at
their prior bytes. The additions are this proposal, inert patch, reproducibility
tools and generated evidence, with hashes in SHA256SUMS. No scientific acceptance
is inferred from these selector tests. Reproduction and the bounded future plan
are in [RETRY_PLAN.md](RETRY_PLAN.md).

**Requested decision:** approve replacing the Stage A restart selection with the
latest accepted endpoint at or before half-time, for each quarter-cap tier and
role, using the exact prepared patch. This approves a protocol repair only;
another Stage A trajectory still requires separate explicit authorization.
