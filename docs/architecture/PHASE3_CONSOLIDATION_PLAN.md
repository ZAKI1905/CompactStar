# Phase 3 entry — behavior-preserving consolidation scoping audit

> **STATUS: `PHASE-3 ENTRY/SCOPING AUDIT COMPLETE`.** Planning only — **no production source,
> test, baseline, ADR or invariant was modified.**
>
> **Headline finding, and it changes the plan:** the roadmap states *"Every item is engineering
> class"* (`MODERNIZATION_ROADMAP.md:298`), but **three of the five Phase-3 items are
> STRUCTURAL** under `GOVERNANCE.md:51`, and **two sit on active fail-closed conditions**
> (`GOVERNANCE.md:70,72`). `GOVERNANCE.md` (rank 1) and `SCIENTIFIC_INVARIANTS.md` (rank 3)
> both outrank the roadmap (rank 4), so **ADRs are required** — they are not optional
> paperwork, and Phase 3 cannot be executed as pure refactoring.

| Field | Value |
|---|---|
| **Starting `master`** | `d2040d899d10a3d0da54a5a2facd7d9cbec34850` |
| **Phase-3 branch** | `modernization/behavior-preserving-consolidation` |
| **Worktree** | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-consolidation` |
| **Branch HEAD = upstream** | `d2040d89…` |
| **Phase-3 prerequisite** | "Phase 2B baselines exist" (`:322`) — **SATISFIED** (`PHASE2B_CLOSURE.md`) |
| **Change class of this audit** | documentation |

---

## 1. Baseline at entry

| Configuration | Result |
|---|---|
| Full, authenticated `COMPACTSTAR_EOS_DATA_ROOT` | **13/13 PASS** |
| Self-contained (no data root) | **8/8 PASS** (5 external-data tests excluded by the CMake guard) |

Protected baseline artifacts — any Phase-3 increment must leave these **byte-identical** unless
it predeclares a tolerance:

| Artifact | SHA-256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `edaa5e3bc5984e2d8cf5acee6664816083ecf83415a9355582614cc3baf42d6e` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `grid_convergence_cmf_1p6_debug.tsv` | `3be2005f8cfdae2798637e4d51674461c9f56dc36ea48d79ad9459109dcc3c88` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `e04f748536d27331bbd383ce4aa11547d5a4f12f927b5df9192d36522a986194` |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |

## 2. Governing-class correction (rank 1 and 3 over rank 4)

`GOVERNANCE.md:51` — **Structural / architecture**: *"Moves ownership, changes boundaries,
promotes or retires a code path"* → **ADR · `CURRENT_ARCHITECTURE.md` updated in the same
change.**

`GOVERNANCE.md` fail-closed conditions:
- **`:70` #3 — "Uncertain authoritative code path. Two or more implementations are live and no
  document names the authority."** → **ACTIVE** for TOV (§3).
- **`:72` #4 — "Uncertain cache validity. It cannot be shown when a cached quantity is
  invalidated."** → **ACTIVE** for the five INV-12 hazards.

`SCIENTIFIC_INVARIANTS.md:147` — INV-04: *"Single owner. **Structural change; requires ADR.**"*

| Roadmap Phase-3 item | Roadmap says | Actual governing class | ADR? |
|---|---|---|---|
| Canonical TOV path | engineering | **STRUCTURAL** + fail-closed #3 | **YES** |
| Single owner for unit conversions | engineering | **mixed** — engineering for exact duplicates, **scientific/authority** for the conflicting solar mass | **partly** |
| Single owner for proper volume | engineering | **STRUCTURAL** (INV-04 says so explicitly) | **YES** |
| One cache-invalidation rule | engineering | **STRUCTURAL** + fail-closed #4 | **YES** |
| Classify dead/unreachable code | engineering | engineering (classification only) | no |

**This is the single most important output of the audit.** Attempting Phase 3 as
"behavior-preserving refactoring" would violate rank-1 and rank-3 authority.

---

## 3. Target A — TOV path ownership map

**Both paths are genuinely LIVE**, exercised by different built programs.

| | Path 1 — sequence scan | Path 2 — profile |
|---|---|---|
| entry | `TOVSolver::Solve(Axis, dir, file)` (`TOVSolver.cpp:1648`) | `NStar::SolveTOV_Profile` (`NStar.cpp:1117`) |
| inner loop | `RadiusLoop` (`:1741`), called `:1714` | `SingleStarSolveToTOVPoints` (`:2500`) |
| target-mass search | none — scans an `ec` axis | `SolveToProfile` (`:2673`), 25 coarse + bisection, `mass_tol = 1e-4 M☉` |
| star construction | `n_star.Append(...)` (`:1867`) + `FinalizeSurface()` | `BuildFromTOV` via public `NStar(points, labels)` |
| `m` conversion | `y[0]/GSL_CONST_CGSM_SOLAR_MASS` (`:1867`) | same, then `× SUN_M_KM` in `BuildFromTOV` |
| exercised by | `spin_therm_evol_main.cpp:67`, `tov_debug_main.cpp:165`, `main/Examples/sig_omega.cpp:64` | `spin_therm_evol_2_main.cpp:124`, `ns_build_main.cpp:81` |
| **validation coverage** | **NONE** | **all seven Phase-2 harnesses** |
| status | LIVE, **unvalidated** | LIVE, **VALIDATED** (2B-2, 2B-4A, 2B-4B) |

`TOVSolver.cpp:2574` states outright: *"Radius loop — copy of RadiusLoop, but pushing
TOVPoint."* The duplication is **both** orchestration *and* the stepping/step-scaling machinery;
genuinely shared below are `ODE`, `GetEDens`, `p_of_e`, `PressureCutoff` and the EOS splines.

**Which harness uses which path** — every Phase-2 harness uses Path 2:

| Harness | Entry used |
|---|---|
| `passive_cooling_regression`, `heat_capacity_real_star` | `SolveTOV_Profile` |
| `tov_reference_cmf`, `grid_convergence_cmf`, `hartle_moment_inertia_cmf` | `SolveToProfile` + `SingleStarSolveToTOVPoints` + `NStar(pts)` |
| `tov_reference_analytic` | `TOVSolver::ODE` directly |
| `hartle_moment_inertia_analytic` | `NStar(pts)` |

**Are they numerically identical? Unknown — and that is the finding.** No evidence compares
them. Path 1 has zero regression coverage, so subordinating it **cannot be protected by the
current suite**. Any retirement must first add an equivalence harness, or accept that three
built programs change behavior unobserved.

**`BuildFromTOV → Find_MomInertia` side effect.** Confirmed (`NStar.cpp:319`): first-order
Hartle runs on **every** profile construction. `hartle_I_dscmf1_debug.tsv` and the two Hartle
CTests now depend on it. **Do not relocate or remove it in Phase 3** — it is load-bearing for
Phase-2 evidence and its removal is a separate, measured change.

**Verdict: designating a canonical TOV path requires an ADR** (structural + fail-closed #3).
Do **not** pick Path 2 merely because it is newer or better validated — that is the *evidence*
for the ADR, not a substitute for it.

## 4. Target B — unit-conversion ownership map

| Conversion | Expression | Locations | Numerically identical? | Validated consumer? | Class |
|---|---|---|---|---|---|
| `k_B` MeV/K | `8.617333262145e-11` | `PhotonCooling_Details.cpp:347`, `NeutrinoCooling_Details.cpp:214` | — | yes (passive cooling, V1) | — |
| `k_B` MeV/K | `8.617333262e-11` | `CompOSE_Thermo.cpp:566,723` | **NO — differs at 1.7e-11 relative** | yes (heat capacity) | **ENGINEERING** (below every tolerance, but *not* bitwise) |
| km³→cm³ | `1.0e15` | `NeutrinoCooling_Details.cpp:94`, `StarContext.cpp:761` | **YES, exact** | yes | **ENGINEERING CONSOLIDATION** |
| MeV fm⁻³→erg cm⁻³ | `1.602176634e33` | `CompOSE_Thermo.cpp:569,726` | **YES, exact** | yes | **ENGINEERING CONSOLIDATION** |
| `M☉`→km | `Zaki::Physics::SUN_M_KM = 1.476625038050` | `NStar.cpp` (`BuildFromTOV`), `NStar.hpp:289`, `TOVSolver.cpp:2040` | — | yes | — |
| g→`M☉` | `GSL_CONST_CGSM_SOLAR_MASS = 1.98892e33` | `TOVSolver.cpp:1867` and ~6 print sites | **CONFLICTS with the above at 6.2e-5** | yes | **SCIENTIFIC/UNIT AUTHORITY QUESTION** |
| `G`, `c` | `GSL_CONST_CGSM_*` | `TOVSolver::ODE` and variants | consistent within TOV | yes (`3.5e-16`) | leave alone |

### The solar-mass conflict — classification

`SUN_M_KM = 1.476625038050 km` is **exactly the IAU nominal `GM☉/c²`**. GSL's own pair
(`G = 6.673e-8`, `M☉ = 1.98892e33`) gives `1.476716 km`. The composite `g → M☉ → km` round trip
therefore carries a **6.2e-5** systematic offset into the geometric `m [km]` column and hence
into `λ`, the `GeometryCache`, and every proper-volume integral.

**Classification: `SCIENTIFIC/UNIT AUTHORITY QUESTION` → `OWNER/ADR ADJUDICATION REQUIRED`.**

Choosing either constant **changes numbers** and is therefore *not* behavior-preserving.
`PHASE2B_CLOSURE.md` §13 already classified it as Phase-3 debt; this audit refines that: it must
**not** be folded into any engineering consolidation commit. The two Hartle constants are
*more* accurate, but "more accurate" is a scientific decision, not a refactor.

### Partition (§10)

| Group | Content | Class | Expected numeric impact | Before/after TOV? |
|---|---|---|---|---|
| **U-1** | km³→cm³, MeV fm⁻³→erg cm⁻³ — exact duplicates | ENGINEERING | **bit-identical** | either; safest first |
| **U-2** | the two `k_B` precisions | ENGINEERING (documented tolerance) | `≤1.7e-11` relative — below every tolerance but **not bitwise** | after U-1 |
| **U-3** | geometric `G`/`c`/km↔cm inside TOV | leave alone | — | not scheduled |
| **U-4** | solar-mass authority | **SCIENTIFIC — ADR** | `6.2e-5` on `λ`, geometry, all volume integrals | **deferred out of Phase 3** |

## 5. Target C — proper-volume ownership map (INV-04)

| Location | Quantity | Measure | Live? | Validated? |
|---|---|---|---|---|
| `GeometryCache.cpp:230` | canonical weights | `w_V = 4πr² e^Λ`; `w_V e^ν`; `w_V e^{2ν}` | **LIVE** | yes (V1, passive cooling, Hartle) |
| `GeometryCache.cpp:194` | λ fallback | `DeriveLambdaFromMR_(r, m, 1e-15)` — re-derives λ when absent | LIVE | indirectly |
| `NStar.cpp:279` / `:563` | baryon-number integrand | `4π r² n_B` (**coordinate**, then × `e^Λ` elsewhere) | LIVE | `seq.b` only |
| `NStar.cpp:1067` | baryon-number integrand | `4π r² n_B / √(1−2m/r)` — inline `e^Λ` | LIVE | weakly |
| `MixedStar.cpp:163,179` | mixed-star volume | inline equivalent | compiled, **unexercised** | no |
| `DarkCore_Analysis.cpp:97` | dark-core integrand | `4π r²` — **coordinate volume, correctly not proper** | compiled, unexercised | no |

**Not proper-volume — must not be swept in:** `RotationSolver.cpp:309,329,349` are the Hartle
ODE *coefficients* `4πr(ε+p)/(1−2m/r)`; `RotationSolver.cpp:1054-1088` are O(Ω²) source terms;
`TOVSolver::ODE` `4πr²ρ` is the mass ODE. Sharing the token `4πr²` does not make them the same
measure. Distinguish: **proper volume** (`e^Λ`), **coordinate volume** (none), **redshift
weighting** (`e^ν`, `e^{2ν}` — separately applied, correctly), **baryon-number current**
(INV-14), **Hartle radial measure**.

**Canonical owner candidate: `GeometryCache`** — it already owns the canonical form and the
redshift variants, and every validated consumer already reads it.

**Floating-point warning (§12).** `NStar.cpp:1067` computes `1/√(1−2m/r)` inline, while
`GeometryCache` computes `exp(Λ)` from a tabulated/derived `Λ`. These are algebraically equal
but **not bitwise equal** — different operation order and a `1e-15` clamp. Centralizing the
baryon-number integrand on `GeometryCache` therefore **cannot promise bit identity**; it needs a
predeclared tolerance. `seq.b` is not covered by any golden baseline, which lowers the risk but
also means the change would be **unobserved** — an argument for adding coverage first.

## 6. Target D — cache dependency/provenance map

Beyond the five measured hazards, the full inventory:

| Cache | Owner | Payload | Real dependencies | Current key | Canonical-run exposure |
|---|---|---|---|---|---|
| `m_rho_gcm3` | `StarContext` | ρ(r) | `eps` | profile version | none |
| `m_Yq_cache` | `StarContext` | `Y_q(r)` | strong-sector species | version + null-check | none |
| `m_durca_mask` | `StarContext` | mask + boundary | `n_B`, `Y_{n,p,e}`, `r` | profile version | none |
| `m_cv_cache` | `StarContext` | `C_⋆(T∞)` | version, **thermo**, **GeometryCache**, `n_B`, `Y_q`, `ν` | `(version, &thermo)` — **geometry omitted** | none |
| `GeometryCache` | free-standing | metric weights | `r, m, ν, Λ` at construction | **none at all** | none |
| `ProfileVersionedCache<T>` | each consumer | arbitrary | builder-dependent | **numeric version only** | none |
| `NeutrinoCooling::cache_` | driver instance | `K_DU`, `K_MU` | context, **geometry**, ρ, DU boundary | numeric version only | none |
| `StarContext` column pointers | `StarContext` | 7 raw pointers | `DataSet` column objects | **bound once in ctor** | none |
| **`TOVSolver` EOS splines + accels** | `TOVSolver` | `ε(p)`, `ρ(p)`, `ρ_i(p)` | the imported EOS table | **rebuilt only by `ImportEOS`** | none (one import per solver) |
| `RotationSolver::fast_k_`, `stored_omega_bar_` | `RotationSolver` | bracket index, `ω̄(r)` | profile arrays via `SetFastProfilePtrs_` | **none** | recomputed each `FindNMomInertia` |
| `TimeSeriesObserver` lookup maps | observer | pointers by name | registered drivers | string key | non-scientific |

**Smallest coherent rule the evidence supports** — not a framework, three additive pieces:

1. **A provenance token** `(profile identity, profile version)` that `StarProfile` can hand out,
   where *identity* is a stable per-object id (an address is **not** sufficient — storage is
   reused; that is exactly how hazard C reproduces).
2. **`GeometryCache` carries that token** and exposes it, so a holder can ask "am I stale?".
   That alone closes hazard A and makes B and D expressible.
3. **Consumers extend their own key** with the dependency identities they actually read —
   `C_⋆` adds the geometry token, `NeutrinoCooling` adds geometry + identity. No generic
   dependency-tracking machinery.

Plus, separately: **`RefreshDerivedCachesIfNeeded_` re-binds columns** (hazard E), or structural
mutation is prohibited on a bound `StarContext`. The audit prefers **re-binding**, because
prohibiting it would silently break the `NStar(points)` constructor path.

`TOVSolver`'s EOS splines are **out of scope** — they are EOS-keyed, not profile-keyed, and no
hazard was demonstrated. Do not fold them into this rule.

## 7. Target E — dead / unreachable classification (nothing to delete in Phase 3)

| Item | Classification |
|---|---|
| `TOVSolver::Solve` / `RadiusLoop` | **DUPLICATE OF CANONICAL PATH — but LIVE and unvalidated.** Three built programs use it. Not dead. |
| O(Ω²) second-order Hartle (`RotationSolver.hpp:227,268,275`) | **UNREACHABLE CANDIDATE — protected by `GOVERNANCE.md` §5.** INV-08 unverified. **Do not schedule deletion.** |
| `RotochemicalCache` | **UNREACHABLE CANDIDATE — not compiled**, protected §5. Do not delete. |
| `MixedStar` / `Solve_Mixed` / `DarkCore_Analysis` | **LEGACY BUT LIVE** (compiled; `tov_debug_main` references MixedStar). Unvalidated. |
| `main/Test/*` demo programs | **LEGACY BUT LIVE** — 8 built executables; three exercise the Path-1 TOV. |
| Commented rotation fragments (`RotationSolver.cpp:196-209`, `894-895`) | **DEAD** — safe to remove, but zero benefit and nonzero diff noise. |
| `build_xcode/` (123 tracked files) | **GENERATED ARTIFACT** — pre-existing debt; `GOVERNANCE.md` requires an ADR on what may be tracked. Not Phase-3 consolidation. |
| Doxygen output under `docs/doxygen_output/` | **GENERATED ARTIFACT** — same disposition. |

**Nothing is recommended for deletion in Phase 3.** The classification itself is the deliverable.

## 8. Thermal-balance ADR (Pattern A vs B) — **DEFER, no Phase-3 blocker**

Authenticated from source: `PhotonCooling_Details.cpp:356,364` and
`NeutrinoCooling_Details.cpp:898,971` each obtain the **same** canonical `C_⋆(T∞)` from
`ctx.star->HeatCapacityStar_Tinf(..., ctx.geo)` and divide their own luminosity by it — textbook
Pattern A. `passive_cooling_regression` asserts the two denominators are identical to `1e-12`.
**No double-counting and no ownership ambiguity exists.** ADR-0002 §6 already records that
Pattern A *"does not foreclose Pattern B later"*.

Neither cache nor proper-volume work becomes materially simpler under Pattern B. **Explicitly
kept deferred. Do not open this ADR in Phase 3.**

## 9. Cross-target dependency graph (derived from source)

```
   [U-4 solar-mass authority]  ── deferred out of Phase 3 (changes numbers) ──┐
                                                                             │
   TOVSolver::Solve ──┐                                                      ▼
   (LIVE, unvalidated)│                                            m[km] ── λ ── GeometryCache
                      ├─→ StarProfile ──→ StarContext ──→ GeometryCache ──→ w_V, w_V e^ν, w_V e^2ν
   SolveToProfile ────┘        │              │                    │
   (LIVE, validated)           │              │                    ├─→ C_⋆ cache  (key omits geometry)
                               │              │                    ├─→ NeutrinoCooling cache
                               │              │                    └─→ [proper-volume owner]
                               │              └─→ rho / Yq / DU caches (version-gated, correct)
                               └─→ BuildFromTOV → Find_MomInertia → SeqPoint::I  (Hartle side effect)

   U-1/U-2 (exact + k_B duplicates) ── touch thermal + EOS only, no structural coupling
```

**Ordering consequences read off the graph:**

- **Cache repair does *not* depend on TOV consolidation.** Caches key on `StarProfile` *version
  and identity*, which both TOV paths produce identically; the hazards are in the keys, not in
  who built the profile. Repairing caches first will **not** need redoing after TOV work.
- **Proper-volume consolidation *does* interact with caches** — `C_⋆` reads `GeometryCache::WV`
  (`StarContext.cpp:758`). But the direction matters: giving `GeometryCache` a provenance token
  (cache work) does **not** change `w_V`'s value, whereas moving `NStar`'s inline `1/√(1−2m/r)`
  onto `GeometryCache` **does** change floating-point evaluation order. **Cache provenance
  first, proper-volume second** — the reverse order would force the cache keys to be restated
  after the ownership move.
- **Unit U-1 is independent of everything** — pure constant deduplication, bit-identical.
- **TOV consolidation is the highest-risk item and depends on nothing else**, but it requires an
  ADR *and* new equivalence coverage for Path 1 before anything can be retired.

## 10. Baseline-protection and preservation matrix

| Increment | Required regression evidence | Preservation standard |
|---|---|---|
| **3A** unit U-1 (exact duplicates) | full 13/13; all five artifact hashes byte-identical | **BIT-IDENTICAL REQUIRED** |
| **3B** cache provenance + dependency-complete keys | `cache_contract`, `cache_thermal_contract`, `passive_cooling_regression`, `heat_capacity_v1`, `heat_capacity_real_star`; **the `--audit-known-hazards` reproductions must stop reproducing**; full 13/13 | **BIT-IDENTICAL REQUIRED** on the canonical run (it provably never enters a hazard state — `PHASE2B_CLOSURE.md` §6) |
| **3C** unit U-2 (`k_B`) | `heat_capacity_v1`, `heat_capacity_real_star`, `passive_cooling_regression`, full 13/13 | **NEW TOLERANCE MUST BE PREDECLARED** (`≤1.7e-11` relative, derived from the constants themselves *before* running) |
| **3D** proper-volume ownership | `heat_capacity_v1` (V1), `passive_cooling_regression`, neutrino diagnostics, `hartle_moment_inertia_*`, baryon number; full 13/13 | **NUMERICALLY IDENTICAL TO ROUNDOFF** for `w_V` consumers; **NEW TOLERANCE PREDECLARED** for the `NStar` baryon-number path |
| **3E** TOV canonicalization | *new* Path-1↔Path-2 equivalence harness **first**; then `tov_reference_analytic`, `tov_reference_cmf`, `grid_convergence_cmf`, `passive_cooling_regression`, `hartle_moment_inertia_cmf`; full 13/13 | **BIT-IDENTICAL REQUIRED** for Path 2 (it must not move at all); Path 1 behavior change must be **measured and documented** |
| **3F** dead-code classification | compile + link + full 13/13; reference search | **NOT APPLICABLE** (no deletion) |

**No tolerance may be selected after observing new output.** 3C's and 3D's tolerances are
derivable from the constants and the algebra in advance; that is why they are stated here.

## 11. Recommended Phase-3 sequence

| # | Title | Scope | Class | Prereq | ADR | Independently mergeable |
|---|---|---|---|---|---|---|
| **3A** | Centralize exactly-duplicated unit constants | `KM3_TO_CM3` (2 sites), `MeVfm3_to_ergcm3` (2 sites) → one owner | engineering | none | no | **yes** |
| **3B** | Cache provenance and dependency-complete keys | profile identity+version token; `GeometryCache` carries it; `C_⋆` and `NeutrinoCooling` keys extended; `StarContext` re-binds columns | **STRUCTURAL** (fail-closed #4) | 3A optional | **YES** | **yes** |
| **3C** | Unify `k_B` precision | one `k_B` across thermal + EOS | engineering w/ tolerance | 3A | no | **yes** |
| **3D** | Single owner for the proper-volume measure | `GeometryCache` canonical; retire `NStar`/`MixedStar` inline forms | **STRUCTURAL** (INV-04) | **3B** | **YES** | yes |
| **3E** | Canonical TOV path | designate Path 2; subordinate Path 1 — **after** building equivalence coverage | **STRUCTURAL** (fail-closed #3) | 3B, 3D | **YES** | yes |
| **3F** | Dead/unreachable classification record | document only; delete nothing | documentation | none | no | **yes** |

**Deliberately *not* the roadmap bullet order.** The roadmap leads with TOV; the dependency
graph says TOV is the *riskiest* item, is the only one with an unvalidated live path, and
benefits from cache/geometry ownership being settled first.

**Deferred out of Phase 3 entirely:**
- **U-4 solar-mass authority** — changes numbers; needs an owner/ADR decision, not a refactor.
- **Pattern B thermal balance** — no blocker (§8).
- **O(Ω²) and rotochemical candidates** — protected under `GOVERNANCE.md` §5.
- **`build_xcode/` and Doxygen artifact debt** — needs its own tracked-artifact ADR.
- **CI** — outstanding *Phase-2B* infrastructure debt, parallel to Phase 3, not a Phase-3 item
  and not a blocker for it.

## 12. The first implementation task

**`3A` — centralize the exactly-duplicated unit constants.**

Chosen on dependency and risk, not convenience. It is the only Phase-3 item that is
simultaneously: **provably bit-identical** (the constants are literally equal — `1.0e15` and
`1.602176634e33` at four sites), **fully covered** by the existing 13/13 suite plus five golden
artifacts, **prerequisite to nothing and blocked by nothing**, and **impossible to need
redoing** after any later increment. It reduces real architectural ambiguity — the roadmap's
"~15 local re-derivations" — with essentially zero scientific risk.

It deliberately excludes `k_B` (3C, not bitwise) and the solar mass (deferred, changes numbers).

**Stop conditions for 3A:** any artifact hash changes; any test fails; any consolidation would
merge two numerically *different* constants; the change grows beyond constant ownership.

## 13. Stop conditions for Phase 3 generally

- Any increment that would alter a golden artifact without a **predeclared** tolerance.
- Any structural increment attempted **without** its ADR.
- Any attempt to fold the solar-mass authority decision into an engineering commit.
- Any deletion of O(Ω²) or rotochemical candidate code.
- Removal or relocation of the `BuildFromTOV → Find_MomInertia` side effect without separate
  measurement.
- Retiring TOV Path 1 before equivalence coverage exists.
