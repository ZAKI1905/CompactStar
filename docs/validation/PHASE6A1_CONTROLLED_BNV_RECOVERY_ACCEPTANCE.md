# Phase-6A-1 controlled-BNV recovery-plan owner acceptance

> **OWNER-ACCEPTED RECOVERY PLAN**
>
> **ADR-0016 HUMAN-RATIFIED**
>
> **R0 STRUCTURAL RECOVERY AUTHORIZED AFTER CANONICAL INTEGRATION**
>
> **RECOVERY IMPLEMENTATION NOT YET RESUMED**

**Date:** 2026-09-19

**Change class:** documentation/governance only

**Canonical authority:** `961dfa0de6f76df71df4cb98edc8e1b35a5c21b1`

**Failed implementation evidence:** `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c`

**Owner-accepted recovery preflight:** `987d78f17906851dbd91342a2008ba65cf02bba1`

**Accepted preflight document SHA-256:**
`28717431c5c1521768cd1afe58780e06044d09d9bb71eeb31529b5af3b07ea68`

The accepted preflight document is imported byte-identically from its failed-
implementation ancestry onto a fresh branch descended only from the canonical authority.
Its internal `PROPOSED` label is retained as immutable pre-acceptance history; this record
is the later human-owner decision. The imported document identifies both authenticated
lineages and makes clear that the failed implementation is not rejected as scientific work
but remains noncanonical (`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:11-42`).

This record changes no scientific equation, production source, test, build file, governed
baseline, EOS/data byte, literature byte, run card, tolerance, or numerical result. It does
not resume implementation, run a trajectory, move the tangent adapter, or create a BNV
candidate or baseline.

## 1. Human-owner decision

The human-owner decision is recorded faithfully, with symbols normalized to repository
notation:

> I accept the Phase-6A-1 controlled BNV recovery plan at SHA
> `987d78f17906851dbd91342a2008ba65cf02bba1`.
>
> 1. Historical BA12 remains **FAIL**.
> 2. New BA12R three-level BASELINE/REFINED/ULTRA requalification is accepted exactly
>    as predeclared.
> 3. If BA12R passes, REFINED is nominal and ULTRA is witness only.
> 4. Only `CPL-P2-LINEAR-QSS-v1` source and its exact matched no-BNV control may
>    receive new ULTRA runs.
> 5. Relocation of the typed Phase-6 equilibrium-baryon tangent adapter from
>    `CompactStar/Analysis` to the Phase-6 BNV production module is authorized,
>    preserving Phase-5B derivative authority,
>    `t = (partial N_y^eq / partial B)_Omega`, `B_B = B_n + B_p`, numerical budgets,
>    currentness, and BA2-BA5 semantics.
> 6. The Phase-5D producer, comparator, provenance rules, governed baseline, and
>    existing Phase-5 source bytes must not change.
> 7. The relocation must be recorded by narrow ADR-0016 before production relocation
>    occurs.
> 8. Existing successful evidence may be reused exactly as specified.
> 9. No physical BNV rate/model, `n -> chi + gamma`, A18, superfluidity, Regime-II,
>    MixedStar, sliding background, variable `Z`, or widened depletion domain is
>    authorized.
> 10. Every accepted stop condition remains binding; no post-result retuning is
>     authorized.

This acceptance authorizes the bounded recovery plan. The owner has now separately ratified
ADR-0016; R0 structural recovery is authorized only after canonical integration of that
ratification. This remains neither acceptance of the failed implementation nor acceptance of
any numerical candidate, and recovery implementation has not resumed.

### 1.1 ADR-0016 human ratification

The human owner explicitly ratified ADR-0016, **Phase-6 BNV Tangent Adapter Ownership**, at
governance SHA `69b999a062636fb0c03212eca83c77d124ac6f86`.

The ratified owner and paths are exactly:

```text
CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp
CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp
CompactStar::Physics::BNV::EquilibriumBaryonTangent
```

The ratification preserves the adapter solely as a typed Phase-6 consumer of the governed
Phase-5B `Analysis::EquilibriumSequenceNumberDerivative`, not an independent structural-
response authority. It preserves
`t = (partial N_y^eq / partial B)_Omega`, zero-spin
`t_i = B_i / (B_n + B_p)`, the existing tangent values and propagated errors,
closure/currentness semantics, BA2-BA5 behavior, and M2/M17/M18 falsifiers.

It authorizes no change to the Phase-5B derivative authority, existing Phase-5 source bytes,
Phase-5D producer, comparator, provenance rules, or governed baseline; no `k` substitution;
no raw-`G_y S_y` changing-baryon response; and no independent recomputation of `t`. Fresh
Phase-5D provenance and governed-artifact identity remain mandatory before any new BNV
trajectory. ADR-0015 physics and the accepted recovery plan, including BA12R, are unchanged.

## 2. Failures remain historical failures

BA12 remains **FAIL** for `CPL-P2-LINEAR-QSS-v1` at
`t = 23113476562.5 s`: the baseline/refined `x_state` difference is
`8.8250672047873735e-9`, the accepted scaled difference is
`1.7255917120989046`, and the historical limit is `1`
(`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:44-92`). No old threshold is relaxed
or reinterpreted.

BA15 remains **FAIL** because the failed implementation added two files below the
recursively authenticated `CompactStar/Analysis` tree, increasing fresh Phase-5D source
provenance from 91 to 93 entries while changing no existing governed source byte
(`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:94-124`). This is an architectural
ownership/provenance conflict, not a changed Phase-5 scientific result.

## 3. Accepted BA12R numerical requalification

The owner accepts the exact hierarchy and equations in the recovery preflight, without
altering the historical BASELINE or REFINED runs:

| tier | `rtol` | `atol` for `(x_state, eta_e, eta_mu)` | accepted role |
|---|---:|---|---|
| BASELINE | `1e-7` | `(1e-12, 1e-18, 1e-18)` | immutable historical run |
| REFINED | `1e-9` | `(1e-14, 1e-20, 1e-20)` | nominal candidate if BA12R passes |
| ULTRA | `1e-11` | `(1e-16, 1e-22, 1e-22)` | convergence witness only |

For each component/checkpoint, BA12R retains the predeclared definitions of `d_BR`,
`d_RU`, `D_R`, `D_U`, parsing uncertainty `Q`, and
`F = max(D_U, 64 ulp(M), Q)`. It requires refined-to-ULTRA stability
`d_RU / D_R <= 1`; contraction `d_RU / d_BR <= 0.10` when `d_BR > 10 F`; and
floor-limited `d_RU <= 10 F` otherwise
(`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:275-329`). The independently
normalized ledger, endpoint-energy, R18/R20, frozen-validity, matched-control, and
solver-trend requirements remain exactly those predeclared
(`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:331-397`).

Only the `CPL-P2-LINEAR-QSS-v1` source and its exact matched no-BNV control may receive
new ULTRA runs. REFINED is fixed as nominal and ULTRA as witness if every BA12R criterion
passes (`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:399-419`). No P0, P1,
reaction-free, or other new ULTRA run is authorized.

## 4. Accepted ownership recovery and immutability

The owner authorizes a narrow ADR-0016 proposal for relocating the typed adapter to:

```text
CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp
CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp
CompactStar::Physics::BNV::EquilibriumBaryonTangent
```

The adapter remains only a typed consumer/view of the governed Phase-5B
`Analysis::EquilibriumSequenceNumberDerivative`. It may introduce no structural solve,
independent derivative authority, `k` substitution, or raw-`G_y` source route. The complete-
star/domain convention, `B_B = B_n + B_p`, component values and propagated errors,
`b^T t` closure, and sequence/star/domain/currentness refusals remain unchanged
(`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:126-205`).

The recovery must not edit the Phase-5D provenance producer, add ignore paths, change the
comparator, regenerate or replace the governed baseline, change Phase-5B code, change any
existing Phase-5 source byte, or weaken provenance. Fresh empty-scratch Phase-5D provenance
and the governed artifact must match exactly before any recovered trajectory.

## 5. Evidence reuse and required reruns

The reusable historical evidence is exactly the set in the accepted preflight: BA1-BA10b,
M1-M21, P0/P1/P2 event oracles, R18, R-a/R-b/R-c, the static BA13 certificate,
BASELINE/REFINED raw trajectories and hashes, matched controls, protected hashes, and the
BA12/BA15 failure evidence (`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:421-438`).
Reuse does not prove a modified binary.

After a separately ratified ADR-0016 and R0, the dependency-linked reruns are the focused
BA2-BA5 build/tests, M2/M17/M18, nonzero-eta R18, focused BA10a/BA10b, retained
certificate parsing/currentness, fresh Phase-5D producer/comparator regression, and focused
Phase-6 regression (`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:440-460`). Only after
those pass may the two ULTRA trajectories and BA12R run; the unfinished gates and full suite
remain later stages (`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:462-478`).

## 6. Scope exclusions and stop conditions

This acceptance does not authorize a physical BNV rate or model, `n -> chi + gamma`, A18,
superfluidity, Regime-II, MixedStar, sliding stellar background, variable-`Z` evolution,
or depletion beyond `abs(DeltaB)/B0 <= 1e-6`. ADR-0015 physics and every Phase-5
authority remain unchanged.

Every stop condition in the accepted recovery preflight remains binding, including any
need to modify governed Phase-5 bytes or provenance machinery; unexpected BA2-BA5,
M2/M17/M18, R18, BA10a/BA10b, certificate, or Phase-5D identity change; any ULTRA or
BA12R failure; any post-result adjustment; depletion-domain violation; or need for excluded
physics (`PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:506-529`). No tolerance,
card, drive, duration, grid, boundary, or budget may be changed in response to a result.

## 7. Governed status and next action

- Recovery plan: **OWNER-ACCEPTED**.
- Historical Phase-6A-1 implementation: **FAILED / NOT CANONICAL**.
- Historical BA12 and BA15: **FAIL**.
- ADR-0016: **ACCEPTED / HUMAN-RATIFIED**.
- R0 structural recovery: **AUTHORIZED AFTER CANONICAL INTEGRATION**.
- Recovery implementation: **AUTHORIZED BUT NOT YET RESUMED**.
- Canonical BNV candidate: **NONE**.
- Governed BNV baseline: **NONE**.
- Physical BNV rate/model: **NONE**.

After canonical integration, the exact next action is to create a fresh bounded recovery
implementation branch/worktree with explicit read-only access to the failed implementation as
source evidence, without importing its ancestry. Do not begin R0 automatically.
