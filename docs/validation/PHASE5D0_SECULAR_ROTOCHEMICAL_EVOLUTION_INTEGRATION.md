# Phase-5D-0INT — secular rotochemical evolution contract integration

> **PHASE-5D SECULAR ROTOCHEMICAL EVOLUTION CONTRACT — HUMAN-RATIFIED AND CANONICALLY
> INTEGRATED / PRODUCTION IMPLEMENTATION NOT YET BEGUN.**

This is a documentation/status integration record. It reconciles the already human-ratified
ADR-0014/preflight history with the intervening canonical Phase-5C integration. It changes no
production source, test, governed baseline, EOS/data, literature, or accepted scientific contract.

## 1. Authority and identity

| Item | Authenticated value |
|---|---|
| Canonical entry | `f0106c2bbaff4e7c10750c3bd6f618f447c1720e` |
| Phase-5D ratification | `3fca17ebec1d70c57dc85fbf01f920c5387217b4` |
| Common ancestor / exact merge base | `27727016856a6a25a46e447c70e380722ea8ddbf` |
| Integration branch | `integration/phase5d-secular-rotochemical-contract` |
| Integration worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-integration` |
| Merge strategy | `git merge --no-ff --no-commit 3fca17ebec1d70c57dc85fbf01f920c5387217b4` followed by one governed merge commit |
| Intended parent 1 | `f0106c2bbaff4e7c10750c3bd6f618f447c1720e` |
| Intended parent 2 | `3fca17ebec1d70c57dc85fbf01f920c5387217b4` |

The ratified Phase-5D linear history is preserved unchanged:

```text
27727016856a6a25a46e447c70e380722ea8ddbf
  -> 080c5bcbb7c10242b6146da3d9fbee961b3d82e6
  -> 109ebfbbb9543c5f8984f85b097a137f6cce8754
  -> 5f04b5ef7cefc7ceb0d73fb0b3927bfbb28508be
  -> 3fca17ebec1d70c57dc85fbf01f920c5387217b4
```

The commit carrying this record defines `PHASE5D0_INTEGRATION_SHA`; its identity is reported after
commit and push rather than embedded self-referentially here.

## 2. Changed paths and conflicts

The Phase-5D branch contributed these six documentation paths:

- `docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md`
- `docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_PREFLIGHT.md`
- `docs/validation/PHASE5D0_SECULAR_ROTOCHEMICAL_EVOLUTION_RATIFICATION.md`
- `docs/SCIENTIFIC_INVARIANTS.md`
- `docs/MODERNIZATION_ROADMAP.md`
- `docs/architecture/CURRENT_ARCHITECTURE.md`

This integration adds this record as the seventh documentation path. The merge encountered three
content conflicts, all in status documentation:

| Conflict path | Semantic resolution |
|---|---|
| `docs/SCIENTIFIC_INVARIANTS.md` | Retained the canonical Phase-5C integration/portability status and added the accepted ADR-0014 contract plus the ratified INV-11a-f sub-statuses. |
| `docs/MODERNIZATION_ROADMAP.md` | Retained Phase-5C closure and made the now-satisfied Phase-5C dependency historical; added Phase-5D contract integration and implementation-ready/not-implemented status. |
| `docs/architecture/CURRENT_ARCHITECTURE.md` | Retained the canonical Phase-5B/Phase-5C live architecture and recorded ADR-0014's accepted future ownership separation without claiming an implemented secular pipeline. |

There were no source, test, baseline, CMake, EOS/data, or literature conflicts. The net tree change
from canonical entry is documentation-only.

## 3. Preserved canonical state

Phase-5C remains **IMPLEMENTED / VALIDATED / INDEPENDENTLY REVIEWED / HUMAN-RATIFIED /
CANONICALLY INTEGRATED / GOVERNED-REGRESSION PROTECTED / CLOSED FOR GENERIC/FREE-GAS
COEFFICIENT SCOPE**. ADR-0013 remains **ACCEPTED / IMPLEMENTED / INTEGRATED**. GC1-GC12, GC9b,
and GC14 remain PASS under their ratified evidence classifications; GC13 remains
**SOURCE-LIMITED / BLOCKED**.

Phase-5B compiler portability is preserved exactly: only `provenance.build.compiler` is portable
execution provenance. Phase-5C compiler portability is preserved exactly: only
`provenance.toolchain.compiler` is portable execution provenance. Both fields remain mandatory and
truthful; no scientific or other provenance field is exempt.

The governed baseline count remains 10. The entry hashes are:

| Artifact | SHA-256 |
|---|---|
| Phase-5B governed baseline | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` |
| Phase-5C governed baseline | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| Phase-5C reviewed candidate | `a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b` |

No baseline is regenerated or rewritten by this integration.

## 4. Accepted Phase-5D contract, unchanged

ADR-0014 is **ACCEPTED / HUMAN-RATIFIED / CANONICALLY INTEGRATED**. The scientific preflight is
**COMPLETE / INDEPENDENTLY REVIEWED / HUMAN-RATIFIED / NOW CANONICALLY PRESENT**. The controlled
non-superfluid free-gas contract is ready for a separately governed implementation task; production
implementation has **NOT BEGUN**.

The merge does not alter any ratified physics: the evolved state is
`(eta_npe^infinity, eta_npmu^infinity)`; `dot eta = -Z R + 2 W Omega Omega_dot`;
`eta^infinity=e^nu eta_local`; reaction-rate integration uses `e^(+nu)` and neutrino luminosity
uses `e^(2nu)`; `Ltilde` has units `erg s^-1 K^-q`; chemical power is formed in MeV/s and converted
exactly once at the thermal boundary; the benchmark enables `{Me,Mmu}` and disables `{De,Dmu}`;
the corrected final `H_M` denominator is `pi^8` without claiming a published erratum;
`D_a subseteq D`; microphysics normalization and global stellar integration remain separate;
`Z/W` are frozen v1 inputs; and spin history remains externally supplied.

## 5. Retained review caveats

All owner-ratified caveats remain controlling:

- the benchmark enabled-process set requires explicit owner authority;
- the below-guard direct-Urca sliver remains a constructible negative control;
- RE10b numerical tolerance must be predeclared and use nontrivial `nu`/`lambda`;
- its wrong-domain fixture must have nonzero `S_a` outside `D_a`;
- the `G_y`/domain-consistency mutation relationship remains explicit;
- upstream INV-11 authority remains cited;
- status banners are evidence labels, not independent validation;
- the historical M10/M11 mislabel is retained as one algebraic mutation, not duplicate coverage;
- ordering-dependent integrations must assert innermost-first profile ordering;
- disconnected or outer-region direct-Urca support requires a production negative control; and
- the controlled benchmark does not validate realistic `Ltilde` construction or normalization.

Realistic FR2005 reproduction remains source-limited and blocked on authenticated A18/APR,
phase/crust construction, YKGH2001 normalization and `alpha_n`, effective masses, support authority,
and benchmark arrays or governed digitization. No realistic A18 implementation is present.

## 6. INV-11 status

| Subpart | Canonical status after reconciliation |
|---|---|
| INV-11a | Coefficient-object redshift semantics **PARTIALLY RESOLVED upstream** by ADR-0013 and extended/clarified by accepted ADR-0014. |
| INV-11b | Evolved eta-state ownership **CONTRACT RESOLVED / IMPLEMENTATION PENDING**. |
| INV-11c | Reaction sign/index convention **CONTRACT RESOLVED / IMPLEMENTATION PENDING**. |
| INV-11d | Thermal energy ledger/no-double-counting semantics **CONTRACT RESOLVED / IMPLEMENTATION PENDING**. |
| INV-11e | Frozen coefficient lifetime/update policy **CONTRACT RESOLVED / IMPLEMENTATION PENDING**. |
| INV-11f | ODE/source coupling **UNRESOLVED / IMPLEMENTATION + VALIDATION PENDING**. |

Global INV-11 remains **UNRESOLVED**.

## 7. Scope closure

No Phase-5D production code or test is introduced. Eta evolution, Urca rates, nonequilibrium
neutrino correction, chemical heating, thermal coupling, A18, and BNV have **NOT BEGUN**. The
reconciliation authorizes only the next separately governed controlled free-gas production task.
