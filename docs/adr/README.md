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
| [ADR-0006](ADR-0006-hartle-first-order-physical-normalization.md) | Physical normalization and unit contract for Hartle first-order rotation | **ACCEPTED** | 2026-09-02 — Q1 = A + D (public spin input is physical `Ω` in rad s⁻¹, carried by an explicit typed quantity), Q2 = A (seed strictly internal), Q3 = A (one canonical geometric representation + named physical accessors), Q4 = A (seed-free normalized response exposed through `NStar`); binding clarification: **no implicit physical spin on `NStar` construction** |
| [ADR-0007](ADR-0007-hartle-second-order-monopole-response.md) | Hartle O(Ω²) monopole structural response at fixed central energy density | **ACCEPTED** | 2026-09-02 — variable `p₀*`; fixed-`ε_c` family with no surface condition; seed-free coefficients per `Ω_geom²` materialized only at an explicit `AngularVelocity`; `dε/dp` owned by the EOS/TOV layer on the star's own `ε(p)` interpolant; surface `ADEQUATE AS-IS` with explicit `R_*` semantics; `l = 0` is the Phase-4 deliverable that unlocks Phase 5 (`l = 2` a separate future extension, **not** validated); atomic candidate replacement in 4C-I1; `GOVERNANCE.md` §3.1 record accepted. **Modified at acceptance:** the homogeneous sequence-derivative response is **not** public in 4C. **Implemented:** 4C-I0 (EOS `dε/dp` authority, 2026-09-02) and **4C-I1 (2026-09-03 — §3.1 correction EXECUTED: candidate deleted, governed monopole response in place; independent physical validation still pending, no baseline)**; **4D (2026-09-03): `HARTLE MONOPOLE VALIDATION FAILED` on stepped crusts — implementation independently verified (analytic `9.7e-9`, DS(CMF)-1 `≤ 3.8e-7`, Chandrasekhar–Miller 1974 and Hartle–Thorne 1968 second-order tables reproduced), but the smooth `dε/dp` authority omits Hartle's internal delta-function shells at the crust's density steps: ≈ `4.6 %` of `δM̂` on DS(CMF)-1 (§7 item 11 `1.04e-3` vs `1e-3`). Not repaired; no baseline; **amendment required** — drafted as ADR-0008 (PROPOSED 2026-09-03). **4D-RV (2026-09-03): corrected revalidation on the ADR-0009 backgrounds — CHARACTERIZED; ADR-0008 D monotonicity not met. 4D-DA (2026-09-04): D adjudicated — `HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED`; first baseline established in Phase 4D-BL (`docs/validation/PHASE4D_MONOPOLE_BASELINE.md:1`)** (`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md` §21) **Phase 4E (2026-09-04): COMPLETE; existing normalized structural interface ratified; Phase 5 not begun** (`docs/validation/PHASE4_CLOSEOUT.md:1`). |
| [ADR-0008](ADR-0008-measure-complete-eos-energy-density-source.md) | Measure-complete EOS energy-density source for the Hartle monopole response | **ACCEPTED** | 2026-09-03 — amends/supersedes ADR-0007 P2/P5/P6 only where EOS energy-density variation is represented by a measure rather than a pointwise radial derivative: `dm̂₀|_EOS = −4πr²ξ̂₀ dε` (H67 eq. 93), realized as a per-segment secant source on the profile partition (optional EOS-knot refinement), an exact jump operator `4πr_t²(ε⁻−ε⁺)ξ̂₀(r_t)` at declared discontinuities (certified `2e-10…4e-9` on a two-layer star), and the surface shell as the terminal atom of the same measure; `dε/dp` retained for the centre series and diagnostics; Phase-5 `δN̂_i` must be measure-complete; no baseline before re-validation. Scratch: DS(CMF)-1 `δM̂` converged to `4e-5` (vs 1.6 % erratic today), +4.8 %. Evidence: `docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md`. **Owner adjudication recorded 2026-09-03:** Q1 Option C (measure); Q2 profile `ε` nodes mandatory, knot snapshot optional/future; Q3 per-segment source with profile nodes as mandatory integration boundaries, no threshold; Q4 no steepness rule; Q5 exact jump, no smoothing; Q6 EOS owns `p_t, ε^±`, profile owns `r_t`, no detector yet; Q7 one unified measure, terminal atom once; Q8 `dε/dp` retained for the centre series and diagnostics; Q9 provenance unchanged; Q10 point-constructed stars unchanged; Q11 Phase-5 must be measure-complete; Q12 no baseline before independent revalidation. **Implemented 2026-09-03 (Phase 4D-RI):** per-segment measure source with profile nodes as hard integration boundaries, terminal atom once, `dε/dp` retained for the centre series; constant-density output bitwise unchanged, smooth-EOS equivalence `≤ 1.3e-5`, same-partition accounting `≤ 3.8e-7`, DS(CMF)-1 `δM̂` `+3.2…+6.5 %`, detectors D1–D4 fire. `CORRECTION IMPLEMENTED — INDEPENDENT REVALIDATION REQUIRED — NO BASELINE`. Evidence: `docs/validation/PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md`. **Revalidated 2026-09-03 (Phase 4D-RV):** independent Stieltjes/EOS-knot measure oracle `≤ 5.3e-5`, A/B/C/E/F/G/J met, M1–M10 fire; **Validation D monotonicity NOT MET** at the node-placement floor (spread `3.7e-5`) — `HARTLE O(OMEGA^2) MONOPOLE RESPONSE CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE`. **4D-DA (2026-09-04): post-validation clarification of D** (the monotonicity heuristic is not a property of the accepted `O(h)` profile-partition measure; D′1–D′5 met) — `HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED`; first baseline established in Phase 4D-BL (`docs/validation/PHASE4D_MONOPOLE_BASELINE.md:1`) **Phase 4E (2026-09-04): COMPLETE; existing normalized structural interface ratified; Phase 5 not begun** (`docs/validation/PHASE4_CLOSEOUT.md:1`). |
| [ADR-0009](ADR-0009-tov-surface-event-and-termination.md) | Ordinary-star TOV surface event and integration-termination contract | **ACCEPTED — SOURCE CONFORMED / NUMERICALLY VALIDATED / ARTIFACT MIGRATION COMPLETE** | 2026-09-03 — accepted Q1–Q14 unchanged; V7a fixed-εc impact and V7b target-mass contract clarified separately. Exact candidate restored; both-EOS sweeps, independent locator, derivatives, downstream checks and mutation detectors pass. Prior stops and owner migration-envelope adjudication retained. Seven artifacts promoted with exact producer reproduction; 41/41 and 20/20 serial tests pass. Corrected Phase-4D independent revalidation is ready, not run; first monopole baseline still requires it. Current hashes and evidence: `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:889`, `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:960`. **Phase 4E (2026-09-04): COMPLETE; existing normalized structural interface ratified; Phase 5 not begun** (`docs/validation/PHASE4_CLOSEOUT.md:1`). |

## Current Phase-4 closeout and baseline chronology

GOVERNANCE §3.1 condition 7, deferred by ADR-0008, is discharged by Phase 4D-BL after independent VERIFIED authority.
Eight current regression artifacts and final complete suites: `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:250`, `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:224`.
Historical status clauses in the index retain their dated scope. **Phase 4E is COMPLETE**:
existing normalized response objects are the ratified Phase-5 structural interface;
INV-08 is CLOSED/VERIFIED for its ordinary-NStar fixed-ε_c l=0 scope. INV-09 remains
unresolved. **PHASE 5 NOT YET BEGUN.** Current contract, eight unchanged hashes and
final serial suites: `docs/validation/PHASE4_CLOSEOUT.md:64`, `docs/validation/PHASE4_CLOSEOUT.md:158`,
`docs/validation/PHASE4_CLOSEOUT.md:183`, `docs/validation/PHASE4_CLOSEOUT.md:219`. Normal reviewed integration is next.

## Anticipated

Identified during Phase-0 reconnaissance or deferred by an accepted ADR; not yet drafted:

| Subject | Blocking |
|---|---|
| ~~Hartle first-order normalization and physical Ω (INV-07)~~ — **ADR-0006 ACCEPTED 2026-09-02** | — |
| ~~Hartle O(Ω²) monopole — perturbation variable, l = 0 equations, fixed-`ε_c` condition, exterior `δM`, EOS `dε/dp` authority, `Ω²`-normalized exposure (INV-08, INV-09)~~ — **ADR-0007 ACCEPTED 2026-09-02** | — |
| Ordinary-star TOV surface event / termination contract — **ADR-0009 ACCEPTED, SOURCE CONFORMED, NUMERICALLY VALIDATED, ARTIFACT MIGRATION COMPLETE** | 4D-RV + 4D-DA independently VERIFIED; 4D-BL first baseline established. Phase 4E complete, Phase 5 not begun; `docs/validation/PHASE4_CLOSEOUT.md:201` |
| ~~ADR-0008 Validation D monotonicity clause vs the measured O(h) node-placement floor~~ — **clarified by owner adjudication, Phase 4D-DA 2026-09-04** (`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md` §21; EOS-knot partition stays optional) | first baseline established (4D-BL); Phase 4E complete, Phase 5 not begun; `docs/validation/PHASE4_CLOSEOUT.md:201` |
| Hartle O(Ω²) **`l = 2`** sector — shape, ellipticity, quadrupole moment `Q`. A separate future rotation extension; ADR-0007's acceptance explicitly does **not** validate it and it blocks neither Phase-4 completion, Phase 5, nor the BNV thermal program | future rotation extension |
| Public ownership of `B_i` / `∂N_i/∂ε_c` (the sequence-derivative response) — deferred out of Phase 4C at ADR-0007's acceptance | Phase 5 |
| Chemical-imbalance definition, ordering, redshift frame (INV-11) | Phase 5 |
| Thermal-balance architecture — where the division by `C_⋆(T∞)` occurs (deferred by ADR-0002 §6) | Phase 3 |
| MixedStar modernization scope — including the six inline proper-volume sites ADR-0004 §15 declines to migrate, and the two-fluid rotation path (`FindMixedMomInertia`, `Solve_Mixed`, stub second order) | post-Phase-3 (deferred by ADR-0004 §0-Q2; needs focused coverage first) |
| Solar-mass-in-km authority (`Zaki::Physics::SUN_M_KM` vs `GSL_CONST_CGSM_SOLAR_MASS`, `6.2e-5` apart) | unadjudicated; deferred out of Phase 3 by the roadmap |
| Core-library vs project-specific-extension boundary (BNV / Decay / DarkCore modules; `Analysis` as the extension seam) | future repository-organisation task (owner clarification, Phase 3F) |
| Generated artifacts in version control | Phase 1 |
| Dependency ownership — Zaki / Confind | Phase 1 |
| Direct-Urca muon channel omission (INV-16) | Phase 5 |
