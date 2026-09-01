# Cache correctness — Phase 2B-3 measurement, Phase 3B repair (INV-12)

> **STATUS: REPAIRED under ADR-0003 (ACCEPTED).**
>
> This document has two layers, deliberately kept separate:
>
> - **Phase 2B-3 (2026-08-31) — the measurement.** Five hazards were reproduced under
>   controlled conditions and quantified. Those numbers are **historical provenance** and are
>   preserved verbatim below. They describe the code **before** Phase 3B.
> - **Phase 3B (2026-09-01) — the repair.** ADR-0003 introduced a `(profile identity, version)`
>   provenance contract. All five hazards are now **enforced contracts under CTest**, not
>   reproductions. The `--audit-known-hazards` mode is **gone**: known-bug output is no longer
>   emitted or treated as expected behavior.
>
> **Read every "CONFIRMED"/"KNOWN HAZARD" verdict below as SUPERSEDED CURRENT BEHAVIOR** — an
> accurate record of what was measured, not a description of the code today.

## 0. Before → after

| # | Hazard (Phase 2B-3) | Measured error, before | Now enforced by | After |
|---|---|---|---|---|
| **H1** | `GeometryCache` had no provenance; staleness unaskable | 51.6 % geometry divergence, undetectable | `cache_contract` **P1a–P1d** | `Matches()` reports staleness; provenance preserved; caller rebuilds |
| **H2** | `C_⋆` key omitted the geometry | **50 %** | `heat_capacity_v1` **U7.d.1–3**, `cache_thermal_contract` **T3a–T3c** | foreign geometry **fails closed**; equivalent provenance still reuses |
| **H3** | `ProfileVersionedCache` keyed on version alone | **85.7 %** | `cache_contract` **P2a–P2d** | keyed on `(identity, version)`; builder re-runs per star |
| **H4** | `NeutrinoCooling` reused across equal-version stars | **80 %** | `cache_thermal_contract` **T4a–T4d** | rebuilds per star; mismatched geometry fails closed |
| **H5** | `StarContext` column views never re-bound | addresses moved on a column-count change | `cache_contract` **P3a–P3d, P4a–P4c** | re-binds before use; invalid schema throws; recovers when repaired |

**Two implementation defects were found by these very tests and fixed** (neither was in the
original five):

- `RefreshDerivedCachesIfNeeded_` early-returned on `!IsValid()` *before* the revision check, so
  a failed re-bind left the context silently degraded rather than failing closed. The revision
  check now precedes that guard. (Caught by **P4b**.)
- Four accessors — `MassDensity_gcm3`, `DirectUrcaMask`, `DirectUrcaLastAllowedIndex`,
  `DirectUrcaBoundaryRadius_km` — tested `IsValid()` *before* refreshing, short-circuiting the
  retry. `IsValid()` is a post-binding predicate, so that was the same wrong order as H5 itself.
  The constructor already throws when `r`/`m` are absent, so `IsValid()` can only become false
  via a failed re-bind; the reorder therefore cannot change behavior for any validly constructed
  context. (Caught by **P4c**.)

---

## 1. Test organization, and why the split matters

Two executables, each with two modes:

**Phase 2B-3 layout (historical).** Two modes: the registered CTest asserted only
currently-correct contracts, while `--audit-known-hazards` reproduced the defects for this
document without ever asserting that stale behavior was correct.

**Phase 3B layout (current).** The audit mode is **removed**. Every case is an ordinary
pass/fail contract in the registered CTest, because ADR-0003 makes each one a requirement.
`WILL_FAIL` is not used anywhere, and no known-bug output is retained as expected behavior.

---

## 2. Dependency / key matrix

`PV` = `StarProfile::Version()`. "Identity" means the profile *object*, not its version.

| Cache | Owner | Actually depends on | Current key | Same-version repeat | Profile mutation | Cross-profile identity | Geometry dependency | Canonical-run exposure | Status |
|---|---|---|---|---|---|---|---|---|---|
| `m_rho_gcm3` (mass density) | `StarContext` | `eps` column | `PV` (shared gate) | **stable** (C2a) | **rebuilds** (C1c/C1d) | n/a — owned per context | none | not reached | **VERIFIED FOR CURRENT CONTRACT** |
| `m_Yq_cache` (charge fraction) | `StarContext` | strong-sector species fractions | `PV` (shared gate) + bare null-check | **stable** (T1b) | **rebuilds** (T1d) | n/a — owned per context | none | not reached | **VERIFIED FOR CURRENT CONTRACT** |
| `m_durca_mask` / boundary | `StarContext` | `n_B`, `Y_n`, `Y_p`, `Y_e`, `r` | `PV` (shared gate) | **stable** (C4a) | **rebuilds** both ways and tracks a partial region (C3b/C3c/C3d) | n/a — owned per context | none | not reached | **VERIFIED FOR CURRENT CONTRACT** |
| `m_cv_cache` (`C_star(T_inf)`) | `StarContext` | `PV`, `thermo`, **`GeometryCache`**, `n_B`, `Y_q`, `nu` | `(PV, &thermo)` | **stable** (U7.a, T1b) | **rebuilds** (U7.b, T1d) | n/a — owned per context | **OMITTED FROM KEY** | not reached (one geometry per run) | **KNOWN HAZARD — geometry omitted** |
| `GeometryCache` | free-standing snapshot | `r`, `m`, `nu`, `Lambda` at construction | **no key at all** | n/a — immutable | **never rebuilds**; no staleness query exists | carries **no** source identity or version | is the geometry | not reached (constructed once, never outlived) | **KNOWN HAZARD — no provenance** |
| `ProfileVersionedCache<T>` | each consumer | whatever the builder reads | `PV` **only** | **stable** (C7b) | **rebuilds** (C7d) | **COLLIDES** (audit C) | not in key | not reached (one context per run) | **KNOWN HAZARD — identity omitted** |
| `NeutrinoCooling` payload | the driver instance | `StarContext`, **`GeometryCache`**, `rho`, DU boundary | `PV` only (via the generic cache) | **stable** (T2b) | **rebuilds** (T2c/T2d) | **COLLIDES** (audit D) | **OMITTED FROM KEY** | not reached (one driver, one context) | **KNOWN HAZARD — identity + geometry omitted** |
| `StarContext` column pointers | `StarContext` | the `DataSet` column objects | bound **once** in the constructor | n/a | **never re-bound** | n/a | n/a | not reached (no structural mutation during evolution) | **KNOWN POINTER-REBIND HAZARD** |

`RefreshDerivedCachesIfNeeded_()` was confirmed against source to reset `m_rho_gcm3`,
`m_durca_mask` / `m_durca_last_allowed` / `m_durca_boundary_r_km`, `m_Yq_cache` and
`m_cv_cache`, and to advance `m_cached_version` — and to **not** call
`BindColumnsOrThrow_()`.

---

## 3. Version is not identity

`StarProfile::m_version` starts at `0` and increments per `Touch()` / per closed
`EditScope`. Two independently constructed profiles that received the same number of edits
therefore carry the **same** version while describing different stars. Measured (C6a): two
profiles built by the same helper both report `Version() == 1` with mean densities
`1.35e14` and `6.73e14 g/cm³`.

**`profile version != profile identity`.** A version-only cache owned by an object that may
be reused across `StarContext`s can collide even when every profile's own mutation tracking
is working perfectly. This is the root cause of hazards C and D.

---

## 4. Verified contracts (the registered regression criterion)

`cache_contract` — 20 checks:

- **C1** mass density builds, and rebuilds by *exactly* the mutation factor (`×2` then
  `×0.25`, relative deviation `0.0`); the version increments on each sanctioned edit.
- **C2** five repeated queries at a fixed version return the identical payload object.
- **C3** the Direct-Urca mask rebuilds across the kinematic threshold **in both directions**,
  and tracks a partial allowed region (boundary index `31` of `64`, exactly the construction).
  With charge neutrality `Y_e = Y_p`, INV-16's `kFn <= kFp + kFe` reduces to `Y_n <= 8 Y_p`;
  the fixtures sit far from that boundary on purpose. The unresolved muon channel is not
  touched.
- **C4** five repeated DU queries return identical mask sum, boundary index and radius.
- **C5** `GeometryCache` reproduces its source geometry at construction (`r`, `exp(nu)`, area).
- **C6** version-vs-identity, above.
- **C7** `ProfileVersionedCache` builds once, does **not** re-run at a fixed version, records
  the version it built against, re-runs after a sanctioned mutation with the correct new
  payload, and honors explicit `Invalidate()`.
- **C8** in-place value edits leave `DataColumn` addresses stable, so pointers bound in the
  `StarContext` constructor stay valid — the boundary of what is currently safe.

`cache_thermal_contract` — 10 checks:

- **T1** `C_star` tracks `(1 + Y_q)` **exactly** (`Y_p: 0.10 -> 0.30`, relative deviation
  `0.0`) and returns when reverted. The fixture is `Q2 = slope·T·(1+Y_q)`, so
  `c_V = T·n_B·slope·(1+Y_q)` analytically. This is the only probe of the `Y_q` cache
  available, since `StarContext::ChargeFractionYq()` is private — a deliberate behavioral
  test rather than a new public accessor.
- **T2** the `NeutrinoCooling` payload is deterministic on repeat, rebuilds after an
  in-place mutation with `L_nu` scaling **exactly** with `rho` (`×3`, relative deviation
  `0.0`, since with `nu = Lambda = 0` the coefficients are `K ∝ ∫ rho · 4πr² dr`), and
  rebuilds again when a composition change closes the DU channel (`L_nu_DU -> 0`, MU
  unchanged). This is §14's same-star semantics and must not be confused with the
  equal-version cross-star hazard in §5.

Deliberately **not** duplicated: `heat_capacity_v1` U7 already covers `C_star` repeat
stability, its rebuild on profile-version change, and its rebuild on a different
`CompOSE_Thermo` object. **Thermo-object dependency: tracked — YES**, re-authenticated by
running that test (green).

---

## 5. Known hazards — reproduced and quantified

> **HISTORICAL — Phase 2B-3 measurement.** Everything in this section describes the code
> **before** Phase 3B and is retained as provenance. Each hazard is now an enforced contract
> under CTest (§0); the reproductions and the `--audit-known-hazards` mode no longer exist.

### HAZARD A — `GeometryCache` has no provenance · KNOWN INV-12 HAZARD

Construct `G_old`, mutate the profile through the sanctioned API (`r ×1.30`,
`Lambda = 0.20`), construct `G_new`. `StarContext`'s own caches rebuild correctly.
`G_old` keeps its construction-time values: `WV[mid]` `3.574e2` vs `7.378e2`, a **51.6 %**
divergence.

An immutable snapshot returning its construction-time values is not itself wrong. **The
defect is that `GeometryCache` exposes no source profile identity, no source version and no
`Invalidate()`, so a holder of `G_old` cannot ask whether it is stale.**

### HAZARD B — heat-capacity key omits the geometry · CONFIRMED

One `StarContext`, one profile version, one thermo; two semantically different
`GeometryCache` objects of identical size (`e^Lambda = 1` vs `2`, so proper volume doubles).

| | `C_star` (erg/K) |
|---|---|
| truth, flat geometry | `1.909267e35` |
| truth, fat geometry | `3.818534e35` |
| 1st call with `geo_flat` | `1.909267e35` |
| 2nd call with `geo_fat` | `1.909267e35` ← **stale** |

**Relative error served on the second call: 50 %.** The truths were computed through fresh
`StarContext`s so no cache is shared. The supplied `GeometryCache` is a genuine input to
`BuildHeatCapacityCache_` but is absent from the key `(prof_version, &thermo)`.

### HAZARD C — `ProfileVersionedCache` omits profile identity · CONFIRMED

Two profiles, `A.Version() == B.Version() == 1`, energy densities differing `7×`. One cache
object, `Get(A)` then `Get(B)`: the builder ran **once**. B was served A's payload —
`1.346591e14` instead of `9.426136e14`, a **85.7 % relative error**, with no diagnostic.

### HAZARD D — `NeutrinoCooling` reused across equal-version stars · CONFIRMED

The concrete consequence of C on the thermal path. Two stars, equal version, `5×` energy
density apart; one reused driver instance.

| | `L_nu` (erg/s) |
|---|---|
| reused driver, star A | `5.640605e38` |
| reused driver, star B | `5.640605e38` ← **stale** |
| **fresh** driver, star B | `2.820302e39` ← correct |

**CONCRETE SILENT CROSS-STAR CACHE COLLISION**, relative error **80 %**, no error, no
warning, no diagnostic flag.

### HAZARD E — `StarContext` column pointers are never re-bound · KNOWN POINTER-REBIND HAZARD

Assessed from source plus **safe** container-identity observation on a profile with **no
`StarContext` bound to it**. No dangling pointer is ever formed or dereferenced; no
undefined behavior is invoked anywhere in this audit.

| mutation | column address | consequence |
|---|---|---|
| in-place value edit | unchanged | bound pointers stay valid (this is C8) |
| clear + refill, same column count | unchanged | bound pointers survive this path |
| **column count `11 -> 14`** (e.g. re-import with more species) | **CHANGED** | **any pointer bound earlier is dangling** |

`BindColumnsOrThrow_()` runs only in the constructor; `RefreshDerivedCachesIfNeeded_()`
invalidates payloads without re-binding. A structural mutation therefore bumps the version —
so payloads rebuild — while `m_r`/`m_m`/`m_nu`/`m_lam`/`m_nb`/`m_pre`/`m_eps` may no longer
refer to live storage. **The version gate gives false confidence here.**

---

## 6. Canonical passive-cooling exposure

The central question: *can any of these hazards silently fire during the established
passive-cooling baseline?*

The canonical procedure is

```
StarProfile (finalized)  ->  ONE StarContext  ->  ONE GeometryCache
                         ->  drivers constructed fresh in scope
                         ->  one continuous Integrate(t0, t1), no structural mutation
```

`passive_cooling_regression` now asserts this on **every** observer callback — not only at
the checkpoints it keeps, so a mid-run swap could not slip through. Measured over **602
driver-context observations** in the canonical run:

| assertion | observed |
|---|---|
| profile version unchanged across the integration | start `5`, finish `5` |
| profile version never changed mid-run | **1** distinct version |
| exactly one `GeometryCache`, and it is this run's | **1**, identity matches |
| exactly one `StarContext` — so no driver is reused across stars | **1**, identity matches |
| exactly one `CompOSE_Thermo` — so the `C_star` key is stable | **1** |
| ADR-0002 Pattern A: photon and neutrino share one `C_star` | retained, green |

**KNOWN CACHE HAZARDS ARE NOT REACHED BY THE CANONICAL PASSIVE-COOLING BASELINE.**

Every hazard above requires either a structural mutation during evolution (A, E) or a
cache/driver reused across two contexts (B, C, D). The canonical run does neither, and now
proves it.

**Scope limit, stated explicitly: this is a claim about one procedure, not about the API.**
Any code that re-solves a profile behind a live `GeometryCache`, or reuses a driver or a
`ProfileVersionedCache` across stars, remains exposed. The general API is **not** safe merely
because this baseline is.

---

## 7. Regression detector proof

Three controlled regressions were applied to production source, measured, and reverted.
`git diff` confirms the tree is byte-identical afterwards, and the full suite is green again.
No scientific formula was perturbed.

| # | Controlled regression | Result | Discrimination |
|---|---|---|---|
| 1 | `StarProfile::Version()` pinned to `1` (version bump suppressed) | `cache_contract` **FAIL** (8), `cache_thermal_contract` **FAIL** (4) | broad, as expected |
| 2 | `RefreshDerivedCachesIfNeeded_` invalidation block bypassed, version still advanced | **FAIL**: C1c, C1d, C3b, C3d, C7e, T1d, T2c, T2d | surgical — every other check stayed green |
| 3 | only the DU mask/boundary reset removed | **FAIL**: C3b, C3d, T2d | exactly the DU-dependent checks, nothing else |

Detector 2 is the sharpest evidence: it leaves `Version()` correct and breaks only
invalidation, and precisely the invalidation-dependent assertions fire.

**Honest limit.** Detector 1 does *not* break the §6 canonical-path assertions — a frozen
version still looks immutable. Those assertions guard against a *changing* profile and
against multiple contexts/geometries; they are not, and are not claimed to be, a detector for
a suppressed version bump. Detectors 2 and 3 cover the invalidation path itself.

---

## 8. Phase-3 repair requirements

Distilled from the evidence above. **Not implemented here.** Several designs satisfy each;
the requirement is stated, not the implementation.

1. **Cache identity must be `(profile identity, profile version)`, not version alone.**
   Evidence: hazards C and D, at 85.7 % and 80 % silent error. A stable per-profile identity
   token (address is *not* sufficient — storage can be reused) plus the version.
2. **`GeometryCache` must carry enough provenance to answer "am I stale?"** — source profile
   identity and source version, plus either rejection or rebuild on mismatch. Evidence:
   hazard A. Today no caller *can* ask.
3. **Cache keys must be dependency-complete.** Any external object that materially changes a
   cached payload must appear in that payload's key. Concretely: `GeometryCache` in the
   `C_star` key (hazard B, 50 % error) and in the `NeutrinoCooling` payload key (hazard D).
   `CompOSE_Thermo` is already keyed and correct — keep it.
4. **`StarContext` must re-bind its column pointers when the profile's structure changes**,
   not merely invalidate derived payloads. Evidence: hazard E. Rebinding inside
   `RefreshDerivedCachesIfNeeded_` is the obvious candidate but not the only one; a structural
   version distinct from the value version would also work.
5. **Consolidate the five coexisting invalidation rules** enumerated in INV-12. The audit
   found no case where the *shared* gate (rule 2) is wrong; the defects are all in what the
   keys omit, not in when the gate fires.

**No new ADR is required.** The audit uncovered no scientific-semantic ownership ambiguity —
every finding is an engineering/architectural defect with an unambiguous owner. ADR-0002's
`C_star` ownership is unaffected.

---

## 9. Reproduction

```bash
ctest --test-dir build -L cache --output-on-failure
```

There is no longer a separate hazard-audit mode; every contract above runs in the
registered tests.
