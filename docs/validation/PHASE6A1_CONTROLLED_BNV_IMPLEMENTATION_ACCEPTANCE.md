# Phase-6A-1A — controlled-BNV implementation-plan owner acceptance

> **OWNER-ACCEPTED / READY FOR BOUNDED CONTROLLED IMPLEMENTATION / NOT YET
> IMPLEMENTED / NOT A GOVERNED NUMERICAL BNV BASELINE.**

**Date:** 2026-09-18
**Change class:** documentation/governance only
**Accepted plan:** `docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md`
**Governing contract:** `docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md`

Authority: `GOVERNANCE.md:54`; `AGENTS.md:11`; accepted ADR-0015 at
`docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md:1`; and the owner-accepted plan
at `docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:1`.

This record changes no scientific equation, production source, test, CMake file,
governed baseline, EOS/data file, literature file, numerical method, computed result,
or physical interpretation. It creates no BNV artifact and begins no implementation.

## 1. Authenticated identity

| Identity | SHA / value |
|---|---|
| Canonical `master` at entry | `15a224c804605877103cf8464d48a89781097f74` |
| Owner-accepted implementation plan | `76ee0277f412be412f856e5fc0c83c58b63c86e0` |
| Accepted branch | `analysis/phase6a1-controlled-bnv-implementation-preflight` |
| Accepted worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-bnv-implementation-preflight` |

At entry, candidate local/upstream/live refs all equalled the accepted-plan SHA and
canonical local/origin/live refs all equalled the canonical-entry SHA. Both worktrees
were clean. Canonical `master` was the exact merge base and an ancestor of the
accepted branch. The commit carrying this record defines `PHASE6A1_ACCEPTANCE_SHA`;
Git reports that SHA after commit and push rather than self-referencing it here.

## 2. Review closure

The final bounded one-hunk confirmation of
`842fe8d38cbccb4ea3e277e55be0ed28db2bfa5b..76ee0277f412be412f856e5fc0c83c58b63c86e0`
returned **0 BLOCKING / 0 MATERIAL / 0 NONBLOCKING / 0 NOTE**. It confirmed the
R20 identity unchanged, the corrected residual normalizer, narrowed M11/M12/M21
detector claims, predeclared reaction-rate resolution, local/infinity event-energy
mapping, current/proposed code citations, target-`B` tolerance, and absence of
load-bearing semantic regression. The accepted plan remains the complete authority
for the validation ladder and falsifier suite.

## 3. Human-owner decision

The human-owner decision is recorded faithfully, with mathematical notation
normalized to the repository's notation:

> I accept the Phase-6A-1 controlled BNV implementation plan at SHA
> `76ee0277f412be412f856e5fc0c83c58b63c86e0`.
>
> I authorize implementation of the first controlled abstract Regime-I neutron-sink
> experiment according to that preflight, including the governed moving-reference
> seam `{Bdot, sigma}`, Phase-5B sequence tangent `t`, actual-potential direct-energy
> ledger, P0/P1/P2 mathematical fixtures, frozen-background validity budget, matched
> no-BNV controls, BA1-BA17 validation ladder, and M1-M21 falsifier suite.
>
> This acceptance does not authorize a physical BNV rate, `n -> chi + gamma` physics,
> realistic A18, superfluidity, Regime-II/MixedStar physics, variable-`Z`/sliding-
> background evolution, or any widening of the declared
> `|DeltaB|/B0 <= 1e-6` frozen-background validity domain.
>
> The first numerical BNV artifact remains a candidate until independently validated
> and separately accepted.

This acceptance is implementation authority for that bounded plan. It is not
implementation validation and does not promote any future numerical result.

## 4. Explicit exclusions and retained limits

This acceptance does **not** authorize:

- a physical BNV rate or physical BNV model;
- `n -> chi + gamma` physics;
- realistic A18;
- superfluidity;
- Regime-II or MixedStar physics;
- variable-`Z` or sliding-background implementation; or
- any widening beyond `|DeltaB|/B0 <= 1e-6`.

The first numerical BNV result remains **CANDIDATE / NOT GOVERNED** until it is
independently validated and separately accepted. No BNV numerical baseline exists.

## 5. Governed status after acceptance

- Phase-6A-1 implementation plan: **OWNER-ACCEPTED**.
- Controlled BNV implementation: **AUTHORIZED / NOT YET BEGUN**.
- ADR-0015: **UNCHANGED / GOVERNED**.
- BNV numerical baseline: **NONE**.
- Physical BNV model/rate: **NONE**.

## 6. Exact next action

Create a fresh implementation branch/worktree from `PHASE6A1_ACCEPTANCE_SHA` and
execute the owner-accepted Phase-6A-1 controlled BNV implementation plan. Do not
begin implementation automatically in this acceptance task.
