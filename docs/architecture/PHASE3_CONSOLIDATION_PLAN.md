# Phase 3 entry — behavior-preserving consolidation scoping audit

> **STATUS: `3A COMPLETE` · `3B COMPLETE` · `3C COMPLETE` · `3D GOVERNANCE / ADR PROPOSED — AWAITING OWNER ADJUDICATION`.**
> The Phase-3 entry/scoping audit that produced this document is complete and its findings
> stand except where §5 is corrected below. Increments 3A, 3B, 3C and 3D have landed and are
> validated, and 3E-0 has measured the two TOV paths. **ADR-0004 is ACCEPTED (2026-09-01) and
> implemented in 3D**; INV-04 is now
> `GOVERNED (ACCEPTED) — CANONICAL VALIDATED PATH CONFORMED; LEGACY MIGRATIONS DEFERRED`.
> **3E is COMPLETE**: [`ADR-0005`](../adr/ADR-0005-canonical-tov-numerical-primitive.md) is
> **ACCEPTED (2026-09-02)** and implemented through I1, I2 and I4; `SingleStarSolveToTOVPoints`
> is the sole ordinary-star radial implementation and `RadiusLoop` is deleted. **Phase 3 is
> closed** — see [`PHASE3_CLOSEOUT.md`](../validation/PHASE3_CLOSEOUT.md).
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
| `DarkCore_Analysis.cpp:97-100` | dark-core integrand | ~~`4π r²` — coordinate volume~~ **CORRECTED by ADR-0004 §4.3: line `:98` divides by `√(1−2M_tot/r)`, so this IS a proper-volume weight, and `:100` builds the `e^{ν}` variant** | compiled, unexercised | no |

**Not proper-volume — must not be swept in:** `RotationSolver.cpp:309,329,349` are the Hartle
ODE *coefficients* `4πr(ε+p)/(1−2m/r)`; `RotationSolver.cpp:1054-1088` are O(Ω²) source terms;
`TOVSolver::ODE` `4πr²ρ` is the mass ODE. Sharing the token `4πr²` does not make them the same
measure. Distinguish: **proper volume** (`e^Λ`), **coordinate volume** (none), **redshift
weighting** (`e^ν`, `e^{2ν}` — separately applied, correctly), **baryon-number current**
(INV-14), **Hartle radial measure**.

**Canonical owner candidate: `GeometryCache`** — it already owns the canonical form and the
redshift variants, and every validated consumer already reads it.

> **Superseded in part by ADR-0004 (ACCEPTED 2026-09-01).** The 3D audit found this wording — *"`GeometryCache`
> canonical; retire `NStar`/`MixedStar` inline forms"* — **not implementable as written**. `MixedStar`
> has no `StarProfile` (`MixedStar.hpp:269,272`) so it cannot construct a `GeometryCache` at all;
> `NStar` builds its integrand inside an open `EditScope` (`NStar.cpp:90,277-281`), where a
> provenance-bearing `GeometryCache` would conflict with **ADR-0003**; and `GeometryCache` throws
> without a `ν` column (`GeometryCache.cpp:178-179`), which the measure does not need. ADR-0004
> proposes splitting the **mathematical** owner (a dependency-neutral primitive) from the
> **cached-representation** owner (`GeometryCache`, unchanged). See ADR-0004 §8, §10, §17.
>
> The audit also raised the inventory from the six rows above to **twenty-two** live or compiled
> proper-volume / redshifted-proper-volume sites (ADR-0004 §4), and recorded **six** mutually
> inconsistent degenerate-input behaviors where INV-03 records one clamp (ADR-0004 §11).

**Floating-point warning (§12).** `NStar.cpp:1067` computes `1/√(1−2m/r)` inline, while
`GeometryCache` computes `exp(Λ)` from a tabulated/derived `Λ`. These are algebraically equal
but **not bitwise equal** — different operation order and a `1e-15` clamp. Centralizing the
baryon-number integrand on `GeometryCache` therefore **cannot promise bit identity**; it needs a
predeclared tolerance. `seq.b` is not covered by any golden baseline, which lowers the risk but
also means the change would be **unobserved** — an argument for adding coverage first.

**Now measured (ADR-0004 §6, §7).** Pointwise the canonical `w_V` and the inline product agree to
**3 ULP** (max relative 6.4e-16). The *integral* is far less sensitive: across 1.0/1.4/1.6/2.0 M☉,
`seq.b` moves by at most **1 ULP — 1.64e-16 relative** — and is bitwise unchanged on three of the
four stars. The predeclared bound, derived from the operation counts **before** any implementation
exists, is **`|ΔB|/B ≤ 1.0e-15`** (ADR-0004 §7.1). Separately measured: `GeometryCache`'s
stored-Λ and derived-Λ routes are **bit-identical** on every node of all four stars (ADR-0004 §9).

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
| **3C** unit U-2 (`k_B`) | `heat_capacity_v1`, `heat_capacity_real_star`, `passive_cooling_regression`, full 13/13 | **PREDECLARED `≤1.7e-11` — HELD FOR THE PHYSICS, NOT MET FOR THE EVOLVED TRAJECTORY.** The bound is **left exactly as predeclared** — it was not widened after the fact. Resolution: the bound *is* satisfied wherever the quantity is a physical observable. Fixed-state response is `1.68e-11`; the static `grid_convergence_cmf_1p6_debug.tsv` moves only its two heat-capacity columns, at `1.70e-11` and `1.71e-11`; `M`, `R`, and all luminosities at fixed `T` move **exactly zero**. The `1.3e-6` seen on the *evolved* passive-cooling trajectory is adaptive step placement, not displacement: with `k_B` held **fixed** and only `rtol` moved `1e-6 → 1e-8`, the same trajectory moves `1.146e-6`, and the old-vs-new difference shrinks 3.8× when the integrator is tightened. Step placement is not a physical observable, so the predeclared bound was applied to the wrong quantity there. Owner adjudicated: adopt and re-baseline. Evidence: [`PHASE3C_BOLTZMANN_AUTHORITY.md`](../validation/PHASE3C_BOLTZMANN_AUTHORITY.md) §14. |
| **3D** proper-volume ownership | `proper_volume_contract`, `geometry_cache_measure_contract`, `baryon_number_cmf`, `heat_capacity_v1`, `passive_cooling_regression`, `hartle_moment_inertia_*`; full suite | **MET.** Every `w_V` consumer **bit-identical**: all five protected artifacts byte-identical, `GeometryCache`'s `ExpLambda`/`WV`/`WVExpNu`/`WVExp2Nu` bitwise on both Λ routes. **`\|ΔB\|/B` predeclared `≤1.0e-15` in ADR-0004 §7.1 before implementation; measured `1.368e-16`** (one ULP, on the 1.0 M☉ star only — 1.4/1.6/2.0 bitwise), matching the ADR's blind §7.2 prediction exactly. Structure (`M`, `R`, `ec`) bitwise on all four stars. **NOT MIGRATED, as planned**: the zero-caller scalar accessor (separate INV-14 defect), TOV Path 1, all six `MixedStar` sites, and every candidate site. |
| **3E** TOV canonicalization | *new* Path-1↔Path-2 equivalence harness **first**; then `tov_reference_analytic`, `tov_reference_cmf`, `grid_convergence_cmf`, `passive_cooling_regression`, `hartle_moment_inertia_cmf`; full 13/13 | **BIT-IDENTICAL REQUIRED** for Path 2 (it must not move at all); Path 1 behavior change must be **measured and documented** |
| **3F** dead-code classification | compile + link + full suite (**executed: 19/19, 10/10**); reference search | **NOT APPLICABLE** (no deletion; documentation-only) |

**No tolerance may be selected after observing new output.** 3C's and 3D's tolerances are
derivable from the constants and the algebra in advance; that is why they are stated here.

## 11. Recommended Phase-3 sequence

| # | Title | Scope | Class | Prereq | ADR | Independently mergeable |
|---|---|---|---|---|---|---|
| **3A** | ✅ **COMPLETE** — Centralize exactly-duplicated unit constants | `KM3_TO_CM3` (2 sites), `MeVfm3_to_ergcm3` (2 sites) → `CompactStar/Units.hpp`. **Bit-identical**: 352 lines of deterministic output unchanged, all five golden hashes unchanged, 13/13 + 8/8. Evidence: [`PHASE3A_UNIT_DUPLICATES.md`](../validation/PHASE3A_UNIT_DUPLICATES.md) | engineering | none | no | **yes** |
| **3B** | ✅ **COMPLETE.** Cache provenance and dependency-complete keys | profile identity+version token; `GeometryCache` carries it; `C_⋆` and `NeutrinoCooling` keys extended; `StarContext` re-binds columns. Contract: [`ADR-0003`](../adr/ADR-0003-profile-cache-provenance-and-invalidation.md) (**ACCEPTED**). All five hazards are now enforced CTests; INV-12 RESOLVED for profile-derived caches; goldens byte-identical. | **STRUCTURAL** (fail-closed #4) | 3A ✅ | **YES — accepted** | **yes** |
| **3C** | ✅ **COMPLETE** — Adopt the authoritative ZakiLib Boltzmann constant | All four production consumers now use `Zaki::Physics::K_BOLTZ_EV * 1.0e-6`; **no `k_B` literal and no GSL Boltzmann constant remain on any production path**. No ZakiLib change, no archive rebuild. Test oracles derive `k_B` independently in `long double` from the SI-exact defining constants and agree with production to **1 ULP**. `tov_dscmf1_reference.tsv` and `hartle_I_dscmf1_debug.tsv` **byte-identical** as predicted (no thermal constant enters them); the three thermal artifacts re-baselined with recorded hashes. 13/13 + 8/8. Evidence: [`PHASE3C_BOLTZMANN_AUTHORITY.md`](../validation/PHASE3C_BOLTZMANN_AUTHORITY.md) | **numerical-method / constant-authority** | 3A ✅ | no | **yes** |
| **3D** | ✅ **COMPLETE IN VALIDATED SCOPE — LEGACY MIGRATIONS DEFERRED.** Canonical owner for the proper-volume measure | Contract: [`ADR-0004`](../adr/ADR-0004-proper-volume-and-metric-measure-ownership.md) (**ACCEPTED 2026-09-01**). Q1 = **Option B**: `CompactStar/Geometry.hpp` owns the mathematics (`f`, `Λ`, `e^Λ`, `w_V`) **and the domain semantics**; `GeometryCache` stays the cached owner; consumers keep their own physics factors and unit conversions. **No new dependency edge** — `Core` does not include `Physics/Evolution`. Q2 = `MixedStar` governed, migration deferred. Q3 = **hybrid physical-domain contract** — exact regular-centre limit, fail closed otherwise, **no `1e-15` clamp** — replacing six mutually inconsistent legacy behaviors. Conformed: `NStar::BuildFromTOV` Λ + baryon integrand, `GeometryCache::DeriveLambdaFromMR_`. **INV-04 is NOT resolved**: TOV Path 1, the scalar accessor, `MixedStar` and candidate code remain governed-but-nonconformant. Evidence: [`PHASE3D_PROPER_VOLUME.md`](../validation/PHASE3D_PROPER_VOLUME.md) | **STRUCTURAL + SCIENTIFIC-SEMANTIC** (INV-04, INV-03) | 3B ✅ | **YES — accepted** | yes |
| **3E** | ✅ **CANONICAL TOV OWNERSHIP COMPLETE.** 3E-0, I1, I2 and I4 done; I3 optional and not taken | ADR-0005 **ACCEPTED**. `SingleStarSolveToTOVPoints` is the canonical numerical primitive and, after **I4**, the **only** ordinary-star radial implementation: `Solve(Axis)`, `SolveToProfile` and `GenTestSequence` all delegate to it, and **`RadiusLoop` is deleted**. `GenTestSequence` was covered **before** migration (16/16 bitwise) and its output is **byte-identical** after. **I2** had already conformed Path-1 geometry to `CompactStar::Geometry`. Both `_Sequence.tsv` and `_TestSequence.tsv` contracts preserved and guarded; no caller changed. **Fail-closed condition #3 CLOSED for ordinary visible-sector TOV radial ownership** (not `MixedStar`, not dark-sector). Evidence: [`PHASE3E_I4_RADIUSLOOP_RETIREMENT.md`](../validation/PHASE3E_I4_RADIUSLOOP_RETIREMENT.md) | **STRUCTURAL** (fail-closed #3) | 3B, 3D ✅ | **YES — accepted** | yes |
| **3F** | ✅ **COMPLETE** — Dead/unreachable/legacy classification and Phase-3 closeout | documentation only; nothing deleted; core-vs-project boundary recorded per owner clarification. Evidence: [`PHASE3_CLOSEOUT.md`](../validation/PHASE3_CLOSEOUT.md) | documentation | none | no | **yes** |

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

**`3A` — centralize the exactly-duplicated unit constants. ✅ COMPLETE** — see
[`PHASE3A_UNIT_DUPLICATES.md`](../validation/PHASE3A_UNIT_DUPLICATES.md). **`3B` ✅ COMPLETE**
(ADR-0003 accepted). **`3C` ✅ COMPLETE** — see
[`PHASE3C_BOLTZMANN_AUTHORITY.md`](../validation/PHASE3C_BOLTZMANN_AUTHORITY.md).
**`3D` ✅ COMPLETE in validated scope** (ADR-0004 accepted and implemented) — see
[`PHASE3D_PROPER_VOLUME.md`](../validation/PHASE3D_PROPER_VOLUME.md).
**`3E-0` ✅ COMPLETE** and **`3E-I1` ✅ COMPLETE** (ADR-0005 accepted and the radial numerics
unified) — see
[`PHASE3E_I1_CANONICAL_TOV.md`](../validation/PHASE3E_I1_CANONICAL_TOV.md). **`3E-I2` ✅ COMPLETE** and **`3E-I4` ✅ COMPLETE** — Path-1 geometry conformed, then the
duplicate radial loop retired; see
[`PHASE3E_I2_PATH1_GEOMETRY.md`](../validation/PHASE3E_I2_PATH1_GEOMETRY.md) and
[`PHASE3E_I4_RADIUSLOOP_RETIREMENT.md`](../validation/PHASE3E_I4_RADIUSLOOP_RETIREMENT.md).
**Phase 3E's roadmap requirement — establish the canonical TOV path and retire or clearly
subordinate the duplicate — is met**, so 3E is complete. **`3E-I3` remains OPTIONAL** (converging
`Append`+`FinalizeSurface` onto `BuildFromTOV`); it is postprocessing consolidation, it does not
block 3E, and it should not be converted into a completion requirement after the fact.
The 3E-0 harness record:
[`TOV_PATH_EQUIVALENCE.md`](../validation/TOV_PATH_EQUIVALENCE.md). It also corrected the caller
count: `TOVSolver::Solve` has **six** live callers, not three (`spin_therm_evol_main`,
`tov_debug_main`, `sig_omega`, `Table_5-8_Glenn`, `coulatt`, `polytrope`).
*(Historical: the next increment at that point was 3E-G.)* The ADR that section anticipated
became **ADR-0005, ACCEPTED 2026-09-02 and fully implemented** (I1, I2, I4).
**Superseded — retained as history:**
[`ADR-0005`](../adr/ADR-0005-canonical-tov-numerical-primitive.md), status **PROPOSED**. It
carries no authority under `docs/adr/README.md` until the owner accepts it, and **no 3E source
change is authorized**. The 3E-G governance task changed documentation only.

> *(This paragraph previously pointed at ADR-0004 and described 3D as governance-only. That was
> leftover 3D-G-era prose: ADR-0004 is the **accepted** proper-volume ADR and 3D is implemented.
> Corrected in the 3E-G documentation commit.)*

Chosen on dependency and risk, not convenience. It is the only Phase-3 item that is
simultaneously: **provably bit-identical** (the constants are literally equal — `1.0e15` and
`1.602176634e33` at four sites), **fully covered** by the existing 13/13 suite plus five golden
artifacts, **prerequisite to nothing and blocked by nothing**, and **impossible to need
redoing** after any later increment. It reduces real architectural ambiguity — the roadmap's
"~15 local re-derivations" — with essentially zero scientific risk.

It deliberately excludes `k_B` (3C, not bitwise) and the solar mass (deferred, changes numbers).

**Stop conditions for 3A:** any artifact hash changes; any test fails; any consolidation would
merge two numerically *different* constants; the change grows beyond constant ownership.

## 12b. Owner design intent — recorded, NOT implemented

> Direction supplied by the owner during Phase 3B-G. This section records **intent for future
> work**. Nothing here describes current state, and **nothing in it was implemented or moved by
> this or any preceding task.**

### ZakiLib / CompactStar units boundary

- **ZakiLib is the reusable general-purpose physics utility library.** Broadly reusable physical
  constants and ordinary conversions belong there when appropriate: the speed of light, Newton's
  constant, the Boltzmann constant, ordinary metric conversions, and similarly general
  physics-textbook quantities.
- **ZakiLib's values were deliberately researched and sourced at modern precision.** Do **not**
  assume GSL supersedes them merely because CompactStar already links GSL. (Phase 2B-4B measured
  a case in point: `Zaki::Physics::SUN_M_KM` is exactly the IAU nominal `GM☉/c²`, while GSL's own
  `G`/`M☉` pair differs from it by `6.2e-5`.)
- **CompactStar should own the specialized compact-object / nuclear-astrophysics conversions and
  domain conventions** — natural/geometric-unit density conversions, fm-based compact-object
  conversions, cgs ↔ geometric pressure and energy-density conversions, and other conventions
  whose meaning is specific to neutron-star work. `CompactStar/Units.hpp` (added in 3A) is
  consistent with this split.
- **`M☉ ↔ km` is a boundary case and is NOT adjudicated by this statement**, particularly because
  the current choices differ numerically. It remains the deferred authority question.

**No constants were moved by this task, and ZakiLib was not modified.**

### CONFIND / ROOT — future dependency direction

- **CONFIND should fundamentally be a numerical library, not a visualization framework.** Its
  desired output is numerical data/files suitable for post-processing.
- **Visualization belongs downstream, preferably in Python.**
- **If** a future audit confirms ROOT is used only for plotting/visualization — and not for
  numerical algorithms, serialization contracts, or required public interfaces — **the intended
  direction is to remove ROOT from CONFIND**, and thereby remove the unnecessary transitive ROOT
  dependence from CompactStar.
- This is **future dependency work, not Phase 3B cache scope.** CONFIND was not inspected in
  depth and its dependency configuration was not altered.

---

## 13. Stop conditions for Phase 3 generally

- Any increment that would alter a golden artifact without a **predeclared** tolerance.
- Any structural increment attempted **without** its ADR.
- Any attempt to fold the solar-mass authority decision into an engineering commit.
- Any deletion of O(Ω²) or rotochemical candidate code.
- Removal or relocation of the `BuildFromTOV → Find_MomInertia` side effect without separate
  measurement.
- Retiring TOV Path 1 before equivalence coverage exists.
- **Repairing `NStar::BaryonNumIntegrand`'s missing `1e54` (or `MixedStar`'s commented-out
  equivalent) inside a proper-volume change.** That is an INV-14 scientific-semantic defect,
  deliberately held separate from 3D — see ADR-0004 §14.
