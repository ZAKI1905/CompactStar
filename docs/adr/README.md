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

## Anticipated

Identified during Phase-0 reconnaissance or deferred by an accepted ADR; not yet drafted:

| Subject | Blocking |
|---|---|
| Hartle first-order normalization and physical Ω (INV-07) | Phase 4 |
| Chemical-imbalance definition, ordering, redshift frame (INV-11) | Phase 5 |
| Thermal-balance architecture — where the division by `C_⋆(T∞)` occurs (deferred by ADR-0002 §6) | Phase 3 |
| Canonical TOV integration path | Phase 3 |
| MixedStar modernization scope | Phase 3 |
| Generated artifacts in version control | Phase 1 |
| Dependency ownership — Zaki / Confind | Phase 1 |
| Direct-Urca muon channel omission (INV-16) | Phase 5 |
