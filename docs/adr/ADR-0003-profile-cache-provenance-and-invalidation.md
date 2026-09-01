# ADR-0003 — Profile-derived cache provenance and dependency-complete invalidation

| | |
|---|---|
| **Status** | **ACCEPTED** |
| **Date accepted** | 2026-09-01 — by owner adjudication |
| **Date drafted** | 2026-09-01 |
| **Change class** | **Structural / architecture** (`GOVERNANCE.md:51`) |
| **Drafted at** | `3cd78fb140c005e9c38d165406b0a8e3add4e546` |
| **Roadmap increment** | Phase 3B (`docs/architecture/PHASE3_CONSOLIDATION_PLAN.md`) |
| **Governing invariant** | INV-12 — remains **unresolved** until implementation lands and validates |

---

## 1. Context

`GOVERNANCE.md` fail-closed condition 4 (`:72`) — *"Uncertain cache validity. It cannot be shown
when a cached quantity is invalidated"* — is **active**. Phase 2B-3 measured five distinct
hazards and demonstrated each under controlled reproduction
(`docs/validation/CACHE_CORRECTNESS.md`). The Phase-3 entry audit classified the repair as
**structural**, not engineering, so an accepted ADR is required before any production change.

This ADR proposes the semantic contract. It writes no code.

## 2. Measured failure evidence

Re-authenticated against source at the commit above.

| # | Hazard | Source evidence | Measured error |
|---|---|---|---|
| **H1** | `GeometryCache` carries no provenance | public surface is the constructor, `Size()` and column accessors only — no version, no source, no `Invalidate()` (`GeometryCache.hpp:86-136`) | 51.6 % geometry divergence, **undetectable by any caller** |
| **H2** | `C_⋆` cache key omits geometry | `StarContext.cpp` keys on `m_cv_cache.prof_version` and `m_cv_cache.thermo_tag` only, while `BuildHeatCapacityCache_` genuinely consumes the supplied `GeometryCache` | **50 %** |
| **H3** | `ProfileVersionedCache` keys on the numeric version alone | `ProfileCache.hpp` `Get()`: `const uint64_t v = ProfileVersion(sc); if (!m_built \|\| v != m_version)` | **85.7 %** |
| **H4** | `NeutrinoCooling` payload omits star identity and geometry | `NeutrinoCooling.hpp:265` `mutable NeutrinoCoolingProfileCache cache_`, keyed via H3 | **80 %** |
| **H5** | `StarContext` binds column pointers once | `BindColumnsOrThrow_()` is called only from the constructor (`StarContext.cpp:130,133`); `RefreshDerivedCachesIfNeeded_` never re-binds | column addresses move when the column count changes |

## 3. Cache/provenance actor map, and what this ADR governs

| Actor | Classification | Governed here? |
|---|---|---|
| `StarProfile`, `StarProfile::Version()` | provenance source | **yes** — supplies the token |
| `StarContext` bound column pointers | view into a profile | **yes** (H5) |
| `m_rho_gcm3`, `m_Yq_cache`, DU mask/boundary | DERIVED FROM STAR PROFILE | **yes** |
| `m_cv_cache` (`C_⋆`) | DERIVED FROM STAR + EXTERNAL DEPENDENCY (thermo, geometry) | **yes** (H2) |
| `GeometryCache` | DERIVED FROM STAR PROFILE (snapshot) | **yes** (H1) |
| `ProfileVersionedCache<T>` | generic, DERIVED FROM STAR PROFILE | **yes** (H3) |
| `NeutrinoCooling` payload | DERIVED FROM STAR + EXTERNAL DEPENDENCY (geometry) | **yes** (H4) |
| `RotationSolver::fast_k_`, `stored_omega_bar_` | ALGORITHM-LOCAL ACCELERATION | **no** — cannot outlive its own solve; no provenance hazard demonstrated |
| `TOVSolver` EOS splines / `gsl_interp_accel` | keyed to the imported **EOS table**, rebuilt by `ImportEOS` | **no** — not profile-derived |
| `TimeSeriesObserver` name→pointer maps | LOOKUP/BOOKKEEPING | **no** — not scientific state |

**This ADR governs only the class that actually shares the provenance problem.** It does not
propose an application-wide dependency framework.

## 4. Decision question

> What provenance must a profile-derived cache carry, what must its validity condition include,
> and who is responsible for detecting and responding to a mismatch?

## 5. Semantic contracts proposed

- **C1 — identity ≠ version.** `Version()` is the *revision counter of one particular profile*,
  not a global star identity. Two independently built profiles routinely both report `1`
  (measured, `cache_contract` C6a). No reusable derived cache may treat the numeric version alone
  as sufficient provenance **across profiles**.
- **C2 — every profile-derived snapshot knows its source**: which profile, and which version.
- **C3 — every cache key is dependency-complete**: it must include every dependency whose change
  can change the payload — and **nothing else**. (See §8: driver options are measurably *not* a
  dependency of the neutrino payload, and must not be added.)
- **C4 — stale state fails closed.** Production must not silently reuse a payload built against
  different provenance. §9 assigns the responsible owner and response per object.
- **C5 — view validity follows structural mutation.** If `StarContext` retains pointers into a
  profile, a sanctioned mutation that may reallocate columns must result in re-binding before
  use, or the architecture must prohibit such mutation. §7 chooses.
- **C6 — same-star repeated access stays cheap.** Correctness must not be bought by rebuilding
  on every call. The Phase-2B-verified same-profile behavior is preserved exactly.

## 6. Alternatives considered — provenance mechanism

| | Mechanism | Correctness | Surface | Lifetime | Verdict |
|---|---|---|---|---|---|
| **A** | `{const StarProfile* source; uint64_t version;}` | closes H1–H4 | minimal; no `StarProfile` data-model change | pointer identity is lifetime-scoped; address reuse conceivable after destruction | **RECOMMENDED** |
| **B** | stable generated `profile_id` + version | also closes H1–H4; immune to address reuse; better diagnostics | changes `StarProfile` construction/copy/move/import semantics | none | viable, larger; **owner question Q2** |
| **C** | derived caches permanently bound to one profile at construction | closes H3/H4 by API | forces drivers to be reconstructed per star | none | **rejected** — see below |
| **D** | anything smaller supported by source | — | — | — | none found |

**Why A.** `StarContext` is constructed **only by callers** — `main/Test/spin_therm_evol_2_main.cpp:138`
and the test harnesses; **no library code constructs one**. The caller therefore already owns
both the profile and everything derived from it, and already outlives them. Option A's lifetime
requirement is thus a *statement* of existing practice rather than a new burden.

**The residual risk, stated plainly.** Under A, a destroyed profile whose address is reused by a
new profile at the same version would be indistinguishable. This ADR closes that **by contract,
not by mechanism**: *a profile-derived cache must not outlive the `StarProfile` it was built
from.* If that contract is violated the cache is already invalid for reasons no token can fix.
If the owner wants the hole closed mechanically, Option B is the answer — **Q2**.

**Why not C.** Permanent binding would forbid reusing a `ProfileVersionedCache` across profiles.
The canonical passive-cooling run constructs fresh drivers per star, so nothing today would
break — but future sequence work (many stars in one process) would then require reconstructing
every driver per star. C is strictly more restrictive than A for no additional correctness.

## 7. Proposed decision

**Adopt Option A: a typed provenance token `(source profile, version)`, carried by every
profile-derived cache, with an explicit lifetime contract.**

### 7.1 `GeometryCache`

- Remains an **immutable snapshot**. It is not refreshed in place; a changed profile means
  constructing a new one.
- **Binds to exactly one (profile, version)** and records it at construction.
- Exposes provenance conceptually as `SourceProfile()`, `SourceVersion()` and a
  `Matches(const StarContext&)` predicate — enough for a holder to ask *"am I stale?"*, which is
  precisely what it cannot do today.
- **No `Refresh()`.** **The caller owns rebuilding**, consistent with callers already owning
  construction.

### 7.2 `StarContext` — Contract **S1** (mutable-profile-aware)

`RefreshDerivedCachesIfNeeded_()` re-runs `BindColumnsOrThrow_()` on a version change, before
rebuilding derived payloads.

**Why S1 and not S2 (immutable binding).** The Phase-2B-verified contracts *require* in-place
mutation with a live context: `cache_contract.cpp` constructs a context at `:186` then mutates
the profile at `:250`; `cache_thermal_contract.cpp` does the same at `:235`/`:251`. S2 would
invalidate tests that Phase 2B-3 verified and that this ADR is obliged to preserve (C6).
Production never mutates a bound profile — `NStar` and `Pulsar` mutate only during construction,
before any context exists — so S1 costs production nothing.

`BindColumnsOrThrow_` already throws when required columns are absent, so a structural change
that invalidates the layout **fails closed** rather than silently rebinding to the wrong column.

### 7.3 `ProfileVersionedCache<T>` — Contract **P1**

Key on **(profile identity, profile version)**. Not P2 (construction-time binding), for the
Option-C reason above.

**External dependencies stay at the consumer**, typed: the generic template knows *profile
provenance only*; a consumer that also depends on geometry or a thermo table extends its own
validity condition. **No `void*` dependency list and no dependency graph.**

## 8. Dependency-complete keys (C3) — derived, not invented

| Cache | Complete dependency set | Change vs today |
|---|---|---|
| `m_rho_gcm3` | profile identity + version | + identity |
| `m_Yq_cache` | profile identity + version (species content participates in the revision) | + identity |
| DU mask / boundary | profile identity + version | + identity |
| `m_cv_cache` (`C_⋆`) | profile identity + version, **thermo identity**, **geometry provenance** | **+ geometry** (H2) |
| `NeutrinoCooling` payload | profile identity + version, **geometry provenance** | **+ identity + geometry** (H3/H4) |

Two negative results, both measured rather than assumed:

- **Thermo needs no version concept.** `CompOSE_Thermo` exposes no setter and no reload — only
  `IsLoaded()`. It is immutable after construction, so **pointer identity is sufficient** and the
  current `thermo_tag` is already correct.
- **Driver options are NOT a cache dependency.** `NeutrinoCooling::SetOptions` (`:207`) can
  change options after construction, but `BuildNeutrinoCoolingCache` reads **no** option — the
  `include_*` flags and `global_scale` are applied at evaluation time in `ComputeDerived`, not
  baked into the payload. Adding options to the key would be inventing a dependency C3 forbids.

## 9. Mismatch response — who enforces what

| Object | Detects | Proposed response |
|---|---|---|
| `GeometryCache` | nothing itself — it is an immutable snapshot | exposes provenance so **consumers** can decide |
| `StarContext` derived caches (`rho`, `Y_q`, DU) | profile revision | **automatic rebuild** (already correct today) |
| `StarContext` column views | version change | **re-bind**; `BindColumnsOrThrow_` **throws** if the layout no longer satisfies the contract |
| `m_cv_cache` | thermo identity or geometry provenance mismatch | **automatic rebuild** |
| `HeatCapacityStar_Tinf(..., geo)` when the *caller supplies* a geometry whose provenance does not match this context's profile | argument check | **fail closed — throw.** This is unambiguously a caller error; no current caller does it, so nothing regresses |
| `NeutrinoCooling` payload | star identity or geometry provenance mismatch | **automatic rebuild** |

**Rebuild is the default; throwing is reserved for caller errors.** Rebuilding a stale payload
restores correctness silently on paths that today return wrong numbers, whereas throwing there
would convert latent staleness into new crashes in code that currently "works". A *supplied*
mismatched geometry is different in kind: silently rebuilding would discard what the caller
explicitly asked for.

## 10. Consequences

**Positive.** H1–H5 all close. INV-12 becomes resolvable. `GOVERNANCE.md` fail-closed condition 4
is discharged for the governed class. Provenance becomes queryable, so future sequence and spin
work can reuse contexts safely.

**Costs.** `GeometryCache` grows two members and a predicate. `ProfileVersionedCache` grows one
member. `StarContext` re-binds on version change (a cost paid only when the profile actually
changes — never on the canonical path). One new argument check may throw where nothing threw
before.

**Explicitly preserved (C6).** Repeated same-profile access stays cheap; sanctioned same-profile
mutation still rebuilds density, `Y_q`, DU, `C_⋆` and the neutrino payload; thermo identity still
invalidates `C_⋆`; the canonical passive-cooling path pays **no** additional rebuild. **Caching
must not be weakened or disabled to achieve correctness.**

## 11. Migration plan (for the accepted ADR, not authorized yet)

1. Introduce the provenance token type and have `StarProfile`/`StarContext` supply it.
2. `GeometryCache` records and exposes it. *(No behavior change yet.)*
3. `StarContext::RefreshDerivedCachesIfNeeded_` re-binds columns.
4. `m_cv_cache` key gains geometry provenance; the supplied-geometry mismatch check is added.
5. `ProfileVersionedCache` key gains profile identity.
6. `NeutrinoCooling` extends its validity condition with geometry provenance.
7. Convert each `--audit-known-hazards` reproduction into a passing correctness assertion.

Steps 1–3 are separately mergeable from 4–6.

## 12. Validation requirements for implementation

**Must stay green:** `cache_contract`, `cache_thermal_contract`, `heat_capacity_v1`,
`heat_capacity_real_star`, `photon_cooling_conformance`, `passive_cooling_regression`, and the
full **13/13** plus **8/8**.

**Must invert.** The five reproductions currently in `--audit-known-hazards` become **passing
correctness assertions**, and the audit mode must stop reproducing them:

| Today (audit demonstrates the bug) | After implementation (asserted) |
|---|---|
| H1 stale geometry undetectable | a holder can determine staleness |
| H2 second call with different geometry reuses the table | rebuilds, or fails closed on a supplied mismatch |
| H3 equal-version different profile reuses payload | **MUST NOT** reuse |
| H4 reused driver serves star A's payload for star B | **MUST NOT** reuse |
| H5 pointers not re-bound | re-bound, or fails closed |

**Known-bug output must not be retained as expected behavior.**

**Golden artifacts.** The five protected artifacts
(`passive_cooling_cmf_1p6_debug.tsv`, `tov_dscmf1_reference.tsv`,
`grid_convergence_cmf_1p6_{debug,trajectory}.tsv`, `hartle_I_dscmf1_debug.tsv`) **must remain
byte-identical**. The canonical run provably never enters a hazard state
(`PHASE2B_CLOSURE.md` §6), so eliminating stale-state paths it never reaches cannot move a
number. **Byte-identity is the default expectation, not a hoped-for outcome.**

## 13. Explicit non-scope

This ADR does **not** decide: proper-volume mathematical ownership (Phase 3D / INV-04 ADR); the
canonical TOV path (Phase 3E ADR); thermal Pattern A vs B (ADR-0002 §6, deferred); physical
constants or unit authority; `k_B`; the solar-mass conversion; Hartle normalization (INV-07);
O(Ω²); rotochemical semantics; BNV; general observer caching; or any application-wide dependency
graph.

## 14. Owner adjudication — DECIDED

Both questions were adjudicated by the owner on 2026-09-01. The alternatives analysis in §6 and
§7 is retained deliberately: it records *why* these were the choices, and what was rejected.

| Question | Decision |
|---|---|
| **Q1 — `StarContext` mutation model** | **ACCEPT S1.** A `StarContext` remains valid across a sanctioned in-place `StarProfile` mutation. On a `Version()` change, before any cached view or payload is used, it must revalidate/re-bind the column views, invalidate dependent payloads, and rebuild lazily. If the changed profile no longer satisfies the expected structural schema it must **FAIL CLOSED** rather than retain stale pointers. This does **not** imply arbitrary mutation is always safe — it means mutations performed through the sanctioned profile-edit/version mechanism are part of the supported contract. |
| **Q2 — provenance identity mechanism** | **ACCEPT Option A.** Runtime provenance is `(source StarProfile object identity, Version())`. Binding lifetime rule: *a profile-derived context/cache/snapshot MUST NOT outlive the StarProfile from which it was built*; pointer identity is therefore valid for the lifetime of the source. **Do NOT add** generated UUIDs, global counters, serialized profile IDs, persistent provenance IDs, or copy/move identity semantics. Persistent or cross-process profile identity would be a separate architectural requirement and decision. |

### The questions as originally posed (retained for provenance)

**Q1 — `StarContext` mutation model.** Should a `StarContext` remain valid across sanctioned
in-place `StarProfile` mutations, automatically re-binding its views and rebuilding derived
caches (**S1**)? Or should a changed profile require constructing a new context, making a context
an immutable snapshot binding (**S2**)?

> **Recommendation: S1.** Evidence is not merely stylistic: the Phase-2B-verified tests mutate in
> place with a live context and expect rebuilds, so S2 would invalidate verified contracts, while
> production never mutates a bound profile and so pays nothing for S1. **Asked anyway** because
> it encodes intended *future* usage — if profiles are meant to be immutable once a context
> exists, S2 is simpler and safer, and those tests would need re-scoping.

**Q2 — provenance identity mechanism.** Is `(profile address, version)` under an explicit
"cache must not outlive its profile" contract sufficient (**Option A**)? Or is a stable generated
profile ID wanted (**Option B**)?

> **Recommendation: A**, since callers already own and outlive both objects. **Asked anyway**
> because it has material future consequences: Option B is immune to address reuse and gives
> better provenance diagnostics, and if profiles will later be copied, serialized, or held in
> long-lived multi-star sequence containers, retrofitting B is more disruptive than adopting it
> now.

No other question is raised: every remaining choice in this ADR is determined by the measured
correctness requirements.


---

## 15. Implementation record (Phase 3B)

**The accepted semantic decision was not modified during implementation.** Nothing in the
evidence contradicted it, so no return to owner adjudication was required.

**Accepted decisions implemented:** Q1 = **S1**, Q2 = **Option A**.

**Components changed**

| Component | Change |
|---|---|
| `Physics/Evolution/ProfileProvenance.hpp` | **new** — typed `{const StarProfile*, uint64_t}`, `IsSet()`, `==`/`!=`. Includes only `<cstdint>` and forward-declares `StarProfile`. No UUID, registry, allocation or serialization. |
| `Physics/Evolution/GeometryCache.{hpp,cpp}` | records provenance at construction; adds `Provenance()`, `SourceProfile()`, `SourceVersion()`, `Matches(ctx)`. **No geometry array changed.** |
| `Physics/Evolution/StarContext.{hpp,cpp}` | `Provenance()`; column views become `mutable` and `BindColumnsOrThrow_()` `const`; S1 re-bind ordering; `HeatCapacityCache::geo_prov`; supplied-geometry mismatch throws |
| `Physics/Evolution/ProfileCache.{hpp,cpp}` | `ProfileVersionedCache` keys on `(identity, version)`; `ProfileProvenanceOf()`; `BuiltAgainst()` |
| `Physics/Driver/Thermal/NeutrinoCooling.hpp`, `src/NeutrinoCooling_Details.cpp` | `cached_geo_prov_`; geometry-change invalidation; fail-closed on a `geo`/`star` mismatch via `ok=false` |

**Three defects were found by the new contracts and fixed** — none was among the original five,
and each is recorded because it shows the tests doing real work:

1. `BuildHeatCapacityCache_` moves a locally built `HeatCapacityCache` into `m_cv_cache` at the
   end, so an earlier `m_cv_cache.geo_prov = …` was silently discarded. The validity condition
   then mismatched on every call and rebuilt the 160-point table each time — the passive-cooling
   regression went from 11 s to a **1500 s timeout**. `geo_prov` is now set on the object that is
   moved in. Caught by runtime, not by any assertion.
2. `RefreshDerivedCachesIfNeeded_` early-returned on `!IsValid()` *before* the revision check, so
   a failed re-bind left the context silently degraded rather than failing closed. The revision
   check now precedes that guard. Caught by `cache_contract` **P4b**.
3. `MassDensity_gcm3`, `DirectUrcaMask`, `DirectUrcaLastAllowedIndex` and
   `DirectUrcaBoundaryRadius_km` tested `IsValid()` *before* refreshing, short-circuiting the
   retry — the same wrong order as the H5 hazard itself. They now refresh first. The constructor
   already throws when `r`/`m` are absent, so `IsValid()` can only become false via a failed
   re-bind; the reorder therefore cannot change behavior for any validly constructed context.
   Caught by `cache_contract` **P4c**.

**Validation**

- **13/13** authenticated and **8/8** self-contained CTests pass.
- All five golden artifacts **byte-identical**, and — a stronger check than hashes alone — the
  artifacts **re-emitted after the change** are byte-identical to the frozen pre-change files:
  `passive_cooling_cmf_1p6_debug.tsv`, `grid_convergence_cmf_1p6_{debug,trajectory}.tsv`,
  `hartle_I_dscmf1_debug.tsv`. `heat_capacity_real_star`'s 247-line deterministic output is
  unchanged.
- Four controlled detector mutations were applied and reverted byte-identically: removing
  identity from the generic key (P2b–P2d fail), removing geometry from the `C_⋆` condition
  (T3b, U7.d.2 fail), suppressing the column re-bind (P3c–P3d, P4a–P4b fail), and removing the
  neutrino provenance guards (the test aborts on an exception escaping `ComputeDerived`, because
  the `StarContext` guard still fires — two independent guards).

**INV-12 disposition:** moved to **RESOLVED for profile-derived caches**, with the scope limit
stated explicitly. Algorithm-local, EOS-keyed and bookkeeping caches were outside this contract
and are unchanged.
