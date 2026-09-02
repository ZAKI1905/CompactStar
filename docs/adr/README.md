# Architecture and Scientific Decision Records

An ADR records a decision that is **expensive to reverse and hard to reconstruct later** —
typically a scientific convention, a numerical method, or an ownership boundary.

## When to write one

Write an ADR when a change is **scientific-semantic** or **structural/architecture** under
`GOVERNANCE.md` §2, or when a fail-closed condition (§3) is hit and the resolution requires
human judgement.

Do **not** write one for engineering refactors, documentation edits, or anything a reader could
reconstruct from the diff alone.

## Lifecycle

```
PROPOSED  →  ACCEPTED   →  (later)  SUPERSEDED by ADR-NNNN
          ↘  REJECTED
```

- **PROPOSED** — written, not decided. Carries no authority.
- **ACCEPTED** — ratified by the owner. Becomes **normative**, ranking directly below
  `GOVERNANCE.md` in the authority hierarchy.
- **REJECTED** — considered and declined. Kept, because the reasoning has value.
- **SUPERSEDED** — replaced. Never deleted; the successor is named in the header.

**Only the project owner may move an ADR to ACCEPTED.** An AI agent may draft, may present
alternatives, and may recommend — but must not mark one accepted.

## Numbering

`ADR-NNNN-short-kebab-title.md`, allocated sequentially and never reused.

## Index

| ID | Title | Status | Decided |
|---|---|---|---|
| [ADR-0001](ADR-0001-species-profile-semantics.md) | Species profile semantics: densities or fractions? | **ACCEPTED** | 2026-08-31 — species columns are dimensionless fractions `Y_i = n_i/n_B` |
| [ADR-0002](ADR-0002-thermal-heat-capacity-ownership.md) | Heat-capacity ownership for the evolved thermal degree of freedom | **ACCEPTED** | 2026-08-31 — one canonical `C_⋆(T∞)`, the GR-integrated EOS-based stellar heat capacity |
| [ADR-0003](ADR-0003-profile-cache-provenance-and-invalidation.md) | Profile-derived cache provenance and dependency-complete invalidation | **ACCEPTED** | 2026-09-01 — provenance is `(StarProfile identity, Version())`; Q1 = S1, Q2 = Option A |
| [ADR-0004](ADR-0004-proper-volume-and-metric-measure-ownership.md) | Proper-volume measure and metric-factor ownership | **ACCEPTED** | 2026-09-01 — Q1 = Option B (dependency-neutral primitive owns the mathematics, `GeometryCache` owns the cached arrays); Q2 = `MixedStar` governed now, migration deferred; Q3 = hybrid physical-domain contract (regular-center limit at `r=m=0`, fail closed otherwise, no `1e-15` clamp) |
| [ADR-0005](ADR-0005-canonical-tov-numerical-primitive.md) | Canonical TOV numerical primitive and sequence workflow interface | **ACCEPTED** | 2026-09-02 — `SingleStarSolveToTOVPoints` is the canonical numerical primitive; Q1 = retain `Solve()` as a subordinate workflow orchestrator, Q2 = preserve the `_Sequence.tsv` contract, Q3 = P3 staged migration, Q4 = preserve the `Analysis`/export hooks |

## Anticipated

Identified during Phase-0 reconnaissance or deferred by an accepted ADR; not yet drafted:

| Subject | Blocking |
|---|---|
| Hartle first-order normalization and physical Ω (INV-07) | Phase 4 |
| Chemical-imbalance definition, ordering, redshift frame (INV-11) | Phase 5 |
| Thermal-balance architecture — where the division by `C_⋆(T∞)` occurs (deferred by ADR-0002 §6) | Phase 3 |
| MixedStar modernization scope — including the six inline proper-volume sites ADR-0004 §15 declines to migrate | Phase 3 |
| Generated artifacts in version control | Phase 1 |
| Dependency ownership — Zaki / Confind | Phase 1 |
| Direct-Urca muon channel omission (INV-16) | Phase 5 |
