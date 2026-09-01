# CompactStar — Current Architecture

> **STATUS: DESCRIPTIVE.** Authoritative for component boundaries and ownership
> (`GOVERNANCE.md` authority rank 6). Describes **only** behavior that is compiled and reachable.
>
> **Current snapshot:** `11ffe45` — the head of roadmap Phase 1 (increments 1A–1D).
>
> **Evidence scope is deliberately mixed, and the layers must not be conflated:**
>
> | Claim area | Authenticated at | How |
> |---|---|---|
> | Build, configuration, warning policy, test plumbing (§5) | **`11ffe45`** | Re-authenticated during this Phase-1 sync against the actual CMake files and a clean configure/build/test |
> | `RotationSolver` merge integrity and buildability | **`57334d8`** | Re-authenticated during Phase 1B against the three-way merge `9f70f14` = `3639d71` + `e60e656` |
> | Everything else — EOS, TOV, Evolution, Drivers, Microphysics, caches | `9f70f14`, from the `d91c31b` audit | Phase-0 reconnaissance; see `docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md` |
>
> **The repository was not re-audited in full at `11ffe45`.** Phase 1 was a build/test phase and
> touched no scientific source. Unchanged-subsystem claims still carry their original Phase-0
> evidence scope and its caveats.
>
> **Historical base commit:** `9f70f14` (`master` merge of owner work `3639d71` + April candidate
> `675b4a9`) — retained because most component-level evidence below still derives from it.
>
> For intended design, see [`TARGET_ARCHITECTURE.md`](TARGET_ARCHITECTURE.md). **Do not read
> intent from this document.**

## Labels used here

| Label | Meaning |
|---|---|
| **LIVE** | Compiled, reachable, and exercised by a program that runs |
| **COMPILED, UNEXERCISED** | Built into the library; no program calls it |
| **UNREACHABLE SCAFFOLDING** | Compiled but cannot be called — no public accessor |
| **NOT COMPILED** | Source exists; absent from every CMake source list |
| **EMPTY** | Zero-byte or body-less file |
| **CANDIDATE** | Present but scientifically unverified (`GOVERNANCE.md` §5) |

---

## 1. What actually executes

One program exercises the full pipeline end to end: `main/Test/spin_therm_evol_2_main.cpp`.

```
CompOSE EOS table
   ↓
TOVSolver::SolveToProfile → SingleStarSolveToTOVPoints     LIVE
   │   integrates in CGS with explicit G, c  (INV-02)
   ↓   unit conversion happens once, at NStar::Append
NStar::prof_ : StarProfile        LIVE — canonical star representation
   │   geometric units: km, km⁻²  ·  versioned (m_version + EditScope)
   ↓
RotationSolver::FindNMomInertia → SeqPoint::I              LIVE, O(Ω) only
   ↓
StarContext(prof_)                LIVE
   ├─→ GeometryCache               LIVE — deep copy, NO version gate (INV-12)
   ├─→ DirectUrcaMask              LIVE — hand-rolled gate
   ├─→ HeatCapacity C(T∞)          LIVE — third, different gate
   └─→ MassDensity g/cm³           LIVE
   ↓
DriverContext {star, geo, envelope, thermo, cfg}    5 raw non-owning pointers
   ↓
StateLayout{Spin:1, Thermal:1}    ← n_eta = 0; the Chem block is never configured
   ↓
EvolutionSystem                   LIVE — plain registration-order loop
   ├─ MagneticDipole      LIVE      dΩ/dt = −K·sign(Ω)·|Ω|ⁿ
   ├─ NeutrinoCooling     LIVE      L_ν ÷ C_⋆(T∞)  ← governed denominator (ADR-0002)
   │                                placeholder Q₀ normalizations
   └─ PhotonCooling       LIVE      L_γ ÷ C_⋆(T∞)  ← same denominator
                                    CONFORMING with ADR-0002 (Pattern A)
   ↓
GSLIntegrator (RKF45, rtol 1e-6 / atol 1e-10)       LIVE
   │   default was MSBDF and unusable — implicit, no Jacobian (Phase 2B-1R)
   │   MSBDF now REJECTED with a diagnostic until Jacobian support exists
   ↓
TimeSeriesObserver + DiagnosticsObserver            LIVE
```

**Everything not on this diagram is scaffolding.**

---

## 2. Component status

### Core

| Component | Status | Note |
|---|---|---|
| `StarProfile` | **LIVE** | Canonical for `NStar`. Column enum + `m_version` + `EditScope` RAII |
| `StarProfileView` | **COMPILED, UNEXERCISED** | `NStar::View()` has zero callers |
| `NStar` | **LIVE** | Migration to `prof_` complete; legacy `ds` is commented out, not dual-live |
| `MixedStar` | **COMPILED, UNEXERCISED** | No surviving `main/` uses it. Master-grid totals added by `3639d71` |
| `TOVSolver` | **LIVE** | Two live integration paths — see §3 |
| `TOVSolver_Thread` | **COMPILED, UNEXERCISED** | Bookkeeping subclass, 124 lines |
| `RotationSolver` — O(Ω) | **LIVE** | Runs on every star build; feeds `SeqPoint::I`. Its profile-backed interpolation path was restored in Phase 1B — see below. **`I = J/Ω` VALIDATED as a scale-free observable (Phase 2B-4B)**; the absolute first-order normalization remains **unresolved** (INV-07) |
| `RotationSolver` — O(Ω²) | **UNREACHABLE SCAFFOLDING · CANDIDATE** | `rot_solver` private, no accessor. Equations untouched by Phase 1B; still unratified under `GOVERNANCE.md` §5 |
| `StarBuilder` | **LIVE** | On the file-reading path only |
| `SeqPoint`, `Prog` | **LIVE** | |

### Evolution

**LIVE:** `State`, `StateVector`, `StateLayout`, `StatePacking`, `RHSAccumulator`, `StarContext`,
`GeometryCache`, `DriverContext`, `EvolutionSystem`, `EvolutionConfig`, `GSLIntegrator`,
`IObserver`, `TimeSeriesObserver`, `DiagnosticsObserver`, `DiagnosticCatalog`,
`DiagnosticPacket`, `RunPaths`, `RunBuilder`, `RunObservers`.

| Component | Status | Note |
|---|---|---|
| `ProfileVersionedCache` | **LIVE — one client** | Only `NeutrinoCooling`. Not used by DU mask, heat capacity, or rotochemical |
| `UnitContract` / `UnitVocabulary` | **COMPILED, UNEXERCISED** | Plumbing real; every producer returns an empty contract |
| `RotochemicalCache` | **NOT COMPILED · CANDIDATE** | Absent from `Physics/Evolution/CMakeLists.txt`; `Build()` has zero callers. **Known nonconformant with ADR-0001** — see below |
| `CheckpointObserver` | **EMPTY** | 0 bytes; not in any CMakeLists |

### Drivers

| Driver | Status | Note |
|---|---|---|
| `MagneticDipole` | **LIVE** | `use_moment_of_inertia` is a logged no-op |
| `NeutrinoCooling` | **LIVE** | Divides `L_ν,∞` by the governed `C_⋆(T∞)` — **conformant** with ADR-0002 (`NeutrinoCooling_Details.cpp:889`, `:968`). Emissivity normalizations self-labeled placeholders; `K_PBF = 0.0`. Carries a null-deref ordering defect — see §4 |
| `PhotonCooling` | **LIVE** | Divides `L_γ,∞` by the canonical `C_⋆(T∞)` from `StarContext::HeatCapacityStar_Tinf` — **conformant** with ADR-0002 (Phase 2A-3). `Options::C_eff` removed. Requires `ctx.star` and `ctx.thermo` when actively cooling; a disabled driver still needs neither |
| `Rotochemical` | **NOT COMPILED · CANDIDATE** | `Driver/Chem/CMakeLists.txt` sources list is empty |
| `HeatingFromChem` | **NOT COMPILED** | Header only; no `.cpp`; commented out of CMakeLists |
| `AccretionTorque`, `BNVSpinTorque`, `BNVSource`, `WeakRestoration`, `Coupling` | **EMPTY** | 0–1 byte files; five appear in `install(FILES …)` rules |
| Envelope: `IEnvelope`, `SurfaceGravity` | **LIVE** | `SurfaceGravity` is the most complete file in the layer |
| `EnvelopePotekhin1997` | **LIVE** | Iron fit real; Accreted variant is a self-declared placeholder |
| `EnvelopePotekhin2003` | **LIVE, MISLEADING** | Forwards verbatim to the 1997 fit |

### Microphysics

| Component | Status |
|---|---|
| `Rates/Urca.hpp` | **NOT COMPILED** — mislabeled copy of `RateChannels.hpp`; no includer; `# add_subdirectory(Rates)` commented out |
| BNV channels / analysis | **COMPILED, UNEXERCISED** — no evolution driver connects them |

---

## 3. Ambiguous authority — unresolved

These are live conflicts. Under `GOVERNANCE.md` §3 they are fail-closed until adjudicated.

1. **Two live TOV integration paths.** `RadiusLoop` (sequence scans) and
   `SingleStarSolveToTOVPoints` (modern profile path) are an acknowledged copy-paste —
   `TOVSolver.cpp:2574` says *"copy of RadiusLoop."* Both are live. No document names a canonical one.
   The `SingleStarSolveToTOVPoints` path — and `SolveToProfile`/`NStar::SolveTOV_Profile` above it —
   is **VALIDATED** as of Phase 2B-2 against the exact Schwarzschild interior solution (to `3.5e-16`)
   and the official CompOSE `eos.mr` (`M_max` to `2.8e-4`, radii to 0.20–0.35 %); see
   `docs/validation/TOV_REFERENCE.md`. That document also records two deliberately unrepaired
   characteristics: the surface is the EOS table floor rather than vacuum, and the default
   `r_max = 70 km` with `radial_res = 10000` leaves ~80 % of the radial grid outside the star.
2. **Two `NStar` profile-construction blocks** — `BuildFromTOV` and
   `InitFromTOVSolver`+`Append`+`FinalizeSurface`, with duplicated hardcoded column layouts.
3. **Proper volume defined in three places** (INV-04).

**Resolved since the Phase-0 audit:**

- **Species semantics** are no longer ambiguous. Per **ADR-0001 (ACCEPTED 2026-08-31)**,
  `StarProfile::BaryonDensity` stores `n_B` in fm⁻³ and species columns store dimensionless
  fractions `Y_i = n_i/n_B`, with `n_i = Y_i n_B` derived at the point of use. `TOVSolver` and
  `NStar` preserve the EOS-supplied values without normalization. See INV-01.
- **Heat-capacity ownership** is no longer ambiguous. Per **ADR-0002 (ACCEPTED 2026-08-31)**,
  the thermal degree of freedom has exactly one physical heat capacity — `C_⋆(T∞)`, the
  GR-integrated EOS/CompOSE-based stellar heat capacity, designated to
  `StarContext::HeatCapacityStar_Tinf` — and every energy channel divides by it. See INV-15.
  **The ownership question is decided; the code does not yet conform** — see the next section.

---

### `PhotonCooling` — ADR-0002 conformance (Phase 2A-3)

**Corrected.** Previously the driver divided `L_γ,∞` by a driver-local constant
`Options::C_eff = 1e40`, hand-set at the call site, while `NeutrinoCooling` used the canonical
`C_⋆(T∞)` — two different heat capacities summed into one state element.

Now:

- `PhotonCooling` obtains `C_⋆(T∞)` from `StarContext::HeatCapacityStar_Tinf`, the same call
  `NeutrinoCooling` makes, so the live equation is
  `C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞`.
- `PhotonCooling::Options::C_eff` **no longer exists**; it was removed rather than deprecated,
  because ADR-0002 permits a constant heat capacity only as an explicit *whole-balance*
  approximation, which a per-channel option can never be.
- The driver requires `ctx.star` and `ctx.thermo` when actively cooling and **fails closed** with a
  diagnostic if either is absent. A deliberately disabled driver
  (`radiating_fraction <= 0` or `global_scale <= 0`) still returns zero and needs neither.
- A new `C_star_erg_K` diagnostic exposes the denominator so later regressions can audit it.
- **ADR-0002 Pattern A is preserved.** No central thermal-energy manager, power-only driver
  interface, or new RHS ownership was introduced; Pattern B remains the deferred Phase-3 question.

Verified by `tests/thermal/photon_cooling_conformance.cpp`, which asserts
`C_star_erg_K == HeatCapacityStar_Tinf(...)`, `dT∞/dt == −L_γ/C_⋆`, and
`d ln T∞/dt == (dT∞/dt)/T∞`, and which **fails against the old `1e40` denominator** (confirmed by
a controlled temporary regression).

**Consequence for historical outputs:** every passive-cooling product generated before this change
is superseded as a validation reference. They are retained, not deleted, and must not be used as a
regression baseline.

**`main/Test/spin_therm_evol_main.cpp` is no longer compatible** with the governed contract: it
builds an empty `StarContext` and sets no `ctx.thermo`, so its photon cooling now fails closed and
contributes zero. That program is an infrastructure demo that relied on the constant precisely
because it has no thermodynamics; wiring a real EOS into it is outside this correction and was not
attempted. `spin_therm_evol_2_main.cpp` supplies both and is unaffected.

### `RotochemicalCache` — ADR-0001 nonconformance

Recorded for Phase 5. **Not repaired here.**

- Status is unchanged: **NOT COMPILED · CANDIDATE**. It is absent from every CMake source list,
  `Build()` has zero callers, and it has never produced output.
- Under the accepted ADR-0001 semantics, its per-species integration path uses the **wrong
  primitive**: `RotochemicalCache.cpp:147` passes the raw stored `Y_i` column into
  `ComputeEnclosedNumber` (`:25-44`) and `ComputeStructuralDerivative` (`:47-104`), both of which
  document and treat the argument as `n_i` in fm⁻³ (`RotochemicalCache.hpp:116,133`). No `× n_B`
  is applied, so `N_i`, `A_i`, and `B_i` are all computed from fractions where densities are
  required.
- This is **pre-existing unvalidated candidate code** from `675b4a9`, not a regression introduced
  by ADR-0001. The ADR made an existing latent defect explicit.

## 4. Known structural hazards

> **Cache provenance is governed by ADR-0003 (ACCEPTED) and was implemented in Phase 3B.**
> The five hazards that Phase 2B-3 measured are **repaired and now enforced by CTest**; the
> `--audit-known-hazards` mode no longer exists. Historical measurements are preserved, labelled
> superseded, in `docs/validation/CACHE_CORRECTNESS.md`.

**Current cache contract.**

- **Runtime provenance is `(StarProfile identity, StarProfile::Version())`**, carried by the
  typed `ProfileProvenance` in `CompactStar/Physics/Evolution/ProfileProvenance.hpp`. There is no
  generated UUID, no global registry and no persistent/serialized profile ID.
- **Lifetime rule:** a profile-derived context, cache or snapshot **must not outlive the
  `StarProfile` it was built from**. Pointer identity is meaningful for exactly that long. In
  practice `StarContext` is constructed only by callers — `main/Test/spin_therm_evol_2_main.cpp`
  and the test harnesses; no library code constructs one — so the caller already owns and
  outlives both.
- **`GeometryCache` is an immutable provenance-carrying snapshot.** It records
  `(profile, version)` at construction and exposes `Provenance()`, `SourceProfile()`,
  `SourceVersion()` and `Matches(ctx)`. There is deliberately **no `Refresh()`**: a changed
  profile means the caller constructs a new one. Its geometry arrays (`R`, `Mass`, `Area`,
  `ExpNu`, `ExpLambda`, `WV`, `WVExpNu`, `WVExp2Nu`, …) are **numerically unchanged** by this
  work — proper-volume ownership remains Phase 3D.
- **`StarContext` follows contract S1** — it stays valid across a sanctioned in-place profile
  mutation. On a revision change `RefreshDerivedCachesIfNeeded_` **re-binds the column views
  first**, then invalidates derived payloads, then advances the cached revision **last**. A
  profile that no longer satisfies the required schema makes `BindColumnsOrThrow_` throw; because
  the revision is advanced only on success, the context is never left falsely marked current, and
  a later repair restores it. The revision check deliberately precedes the `IsValid()` guard, and
  the four derived accessors refresh before testing `IsValid()`, so a failed re-bind fails closed
  rather than degrading to silent nulls.
- **`ProfileVersionedCache<T>` keys on `(identity, version)`**, not the version alone. It knows
  profile provenance **only**; consumer-specific dependencies stay typed at the consumer, and
  there is no generic `void*` dependency list or dependency graph.
- **Heat capacity `C_⋆(T∞)` is keyed on `(profile version, thermo identity, geometry
  provenance)`.** A caller-supplied `GeometryCache` whose provenance does not match the context's
  current profile **throws** — silently substituting a locally built geometry would discard what
  the caller asked for. An *equivalent* geometry rebuilt for the same `(profile, version)` is
  interchangeable and does not force a rebuild: the key is provenance, not object address.
  `CompOSE_Thermo` needs no version because it exposes no setter or reload — pointer identity
  suffices.
- **The `NeutrinoCooling` payload depends on `(profile identity, version, geometry
  provenance)`.** Driver `Options` are **not** in the key: `SetOptions` exists, but the payload
  builder reads no option — `include_*` and `global_scale` are applied at evaluation time. A
  `DriverContext` whose `geo` does not belong to its `star` fails closed through the driver's
  ordinary `ok=false` diagnostic, without terminating.

**Out of the ADR-0003 contract, unchanged:** `RotationSolver`'s bracket/`omega_bar` acceleration
caches (algorithm-local, cannot outlive their own solve), `TOVSolver`'s EOS splines and
`gsl_interp_accel` (keyed to the imported EOS table, not to a profile), and the
`TimeSeriesObserver` name→pointer maps (bookkeeping, not scientific state). No provenance hazard
was demonstrated for any of them.

- **No driver dependency graph.** `IDriver::DependsOn()` / `Updates()` are declared, overridden
  by every driver, and **never called**. `EvolutionSystem.cpp:125-134` iterates in registration
  order. The comment at `IDriver.hpp:102-103` referring to an "Evolution graph" is stale.
- **Declared-but-undefined symbols reachable from compiled code** —
  `Physics::Spin::DipoleFieldEstimate` and `CharacteristicAge` are called at `Pulsar.cpp:204,212`
  with no definition anywhere.
- **Null-deref ordering — ✅ RESOLVED (Phase 2A-3).** `NeutrinoCooling_Details.cpp` used to
  dereference `ctx.star` in its `HeatCapacityStar_Tinf` call twelve lines before the
  `if (!ctx.star)` guard, which therefore could never fire. The guard now precedes first use.
  No emissivity, rate, cache, option, or numerical constant changed.
- **Stiff early transient, explicit stepper.** `T∞` falls from `1e9 K` to `1.2e7 K` between 100
  and ~300 yr under the placeholder neutrino normalizations, so the passive system is effectively
  stiff there. The default `RKF45` is adequate on measured evidence, but one output cadence
  (`samples_per_decade = 300`) produces a non-finite state in that window while 50–250 and 400 all
  complete and agree to 4.3e-6. Cause undiagnosed; recorded in
  `docs/validation/PASSIVE_COOLING_BASELINE.md`. A stiff method with a real Jacobian is the likely
  long-term answer, and `MSBDF` remains unavailable until one exists.
- **Heat-capacity cache key — ✅ REPAIRED (Phase 3B, ADR-0003).** The key was
  `(profile version, thermo pointer)` while the supplied `GeometryCache` was a genuine input, so a
  later call at the same profile version with a different geometry silently reused the earlier
  table — measured at **50 %** error in Phase 2B-3. The key now includes the geometry's
  provenance, and a supplied geometry belonging to a different profile or revision **fails
  closed**. See the cache contract above.
- **Canonical-baseline exposure — NONE, and now additionally protected by the contract above.** `passive_cooling_regression` now asserts, on every one
  of 602 driver-context observations in the canonical run, that exactly one profile version, one
  `GeometryCache`, one `StarContext` and one `CompOSE_Thermo` are in play from start to finish.
  Every hazard above needs either a structural mutation during evolution or a cache/driver reused
  across contexts; the canonical run does neither.

---

## 5. Build reality

Re-authenticated at **`11ffe45`** after roadmap Phase 1. Full evidence and commands:
[`docs/build/MACOS_BUILD.md`](../build/MACOS_BUILD.md).

### Configuration and build

- **A clean macOS checkout configures out of source, and the library builds and links.**
  The seven historical/private program directories (`BNV_2022`, `EOS`, `SpinDown_2022`,
  `SpinDown_2023`, `BNV_2023`, `BNV_2024`, `LightDM_2024`) are gitignored, absent from a clean
  checkout, and are now **optional**: each is added only if it carries a `CMakeLists.txt`,
  otherwise a `STATUS` message is emitted (`main/CMakeLists.txt:7-23`). Tracked
  `main/Examples` and `main/Test` remain **required** and unconditional (`:28-29`), so a
  misspelled required directory still fails configure.
- **Python discovery still needs help on the reference Mac.**
  `find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)`
  (`CMakeLists.txt:82`) fails when the first Python CMake finds lacks NumPy. The canonical
  configure command therefore passes `-DPython3_EXECUTABLE`, or a project `.venv` supplies one.
- **Builds are genuinely out of source.** `configure_file` writes the generated
  `CompactStarConfig.h` into the **binary tree** at
  `${CMAKE_BINARY_DIR}/generated/include/CompactStar/Core/CompactStarConfig.h`
  (`CMakeLists.txt:110-116`), and that generated root precedes the source root on the include
  path (`:137`). A fresh configure leaves tracked source state unchanged.
  **Artifact debt remains:** a stale generated copy is still tracked at
  `CompactStar/Core/CompactStarConfig.h`, alongside an editor duplicate
  `CompactStar/Core/CompactStarConfig 2.h`. Both are inert for the build — the binary-tree copy
  wins — and neither was removed; retiring them is a generated-artifact decision.

### Platform support — unchanged

- **macOS is still the only authenticated development platform.** This repository ships Zaki and
  Confind as prebuilt static archives for `Darwin/{arm64,x86_64}` only
  (`dependencies/lib/{Zaki,Confind}/Darwin/…`), and `CMakeLists.txt:95-101` hard-fails when the
  archive for the host platform is absent. Configuration therefore still fails on other platforms.
- **Dependency source exists externally but is not integrated here.** Per project-owner authority,
  ZakiLib (private) and CONFIND (public) exist as external version-controlled repositories. That
  is **not** the same as authenticated binary equivalence or cross-platform support: their current
  source revisions have **not** been shown to reproduce the archives CompactStar consumes, and
  their build contracts have not been authenticated. This repository contains no dependency
  source, no submodule, no `FetchContent`, and no package-manager integration.
  **Cross-platform status is unchanged.**

### Build configuration, warnings, diagnostics

- **A default development build type now exists.** `Debug` is selected when the user supplies
  none, and only for single-configuration generators (`CMakeLists.txt:35-43`); a user-supplied
  `-DCMAKE_BUILD_TYPE` is never overridden and multi-config generators are untouched.
  `Debug` is a **development** default, not a production or scientific configuration, and
  **`Debug` and `Release` are not claimed to be numerically equivalent** — no comparison exists.
- **A warning policy exists for CompactStar-owned code.** `-Wall -Wextra`, applied **`PRIVATE`**
  to the `CompactStar` target for AppleClang/Clang/GNU (`CMakeLists.txt:161-166`), so the flags
  reach neither consumers of the installed library, nor `main/`, nor `tests/`.
- **Dependency headers are `SYSTEM`-classified** (`CMakeLists.txt:144-149`), so GSL and the
  vendored Zaki/Confind headers compile with `-isystem` and their diagnostics are excluded from
  CompactStar's inventory. This suppresses reporting; it does not fix those dependencies.
- **`-Werror` is not enabled.** A clean `Debug` build succeeds **with warnings**; they were
  inventoried during Phase 1D and **none was repaired**. Several classes touch scientific/EOS
  code — hidden overloaded virtuals in the EOS particle hierarchy, VLAs in `SigmaOmegaRho*`, and
  set-but-unused variables in `TOVSolver` — and require classified review before any change.
  See [`docs/build/MACOS_BUILD.md`](../build/MACOS_BUILD.md) for the full table.
- **No sanitizer policy and no assertion policy exist.** Neither has been established. Production
  source still contains no assertions; the checks introduced in Phase 2A–2B live in `tests/`, not
  in the library.

### Automated tests

- **CTest infrastructure exists.** `include(CTest)` at top level provides standard `BUILD_TESTING`
  and `enable_testing()`; `tests/` is added only when testing is on (`CMakeLists.txt:222-226`).
- **`tests/` is the canonical automated-test root**, holding **13 registered CTests** as of the
  Phase-2B closure audit (`tests/CMakeLists.txt`). Eight are self-contained; five are registered
  only when `COMPACTSTAR_EOS_DATA_ROOT` supplies the authenticated CompOSE tables
  (`tests/CMakeLists.txt:71` — the guard *excludes* them, it does not skip them):
  `heat_capacity_real_star`, `passive_cooling_regression`, `tov_reference_cmf`,
  `grid_convergence_cmf`, `hartle_moment_inertia_cmf`. Measured: **13/13 pass** with the data
  root, **8/8 pass** without it.
- **`main/Test/` remains manual demo/debug programs**, not a suite. None of its eight executables
  is registered with CTest.
- **The smoke test is infrastructure validation, not a scientific baseline.** It establishes that
  CompactStar headers compile, that a test executable links a genuine out-of-line symbol from
  `libCompactStar.a` (`CompactStar::Core::Prog::GetName()`, verified present as a defined `T`
  symbol rather than an unresolved import), that transitive dependencies link, and that the binary
  runs and exits zero under CTest. **It asserts no scientific value of any kind.**
- **A scientific test suite now exists**, built across Phase 2A–2B: heat capacity
  (`tests/thermal/heat_capacity_v1.cpp`, `heat_capacity_real_star.cpp`), photon conformance
  (`photon_cooling_conformance.cpp`), the stepper contract
  (`evolution_stepper_contract.cpp`), the passive-cooling regression baseline
  (`passive_cooling_regression.cpp`), TOV references (`tests/core/tov_reference_*.cpp`), cache
  contracts (`tests/evolution/cache_contract.cpp`, `tests/thermal/cache_thermal_contract.cpp`),
  grid convergence (`tests/thermal/grid_convergence_cmf.cpp`) and the scale-free Hartle moment
  of inertia (`tests/rotation/hartle_moment_inertia_*.cpp`). Each has a durable evidence record
  under `docs/validation/`. **Scope limits are recorded in those records and are not superseded
  by the existence of a green suite.**
- **No CI exists.** No `.github/`, `.gitlab-ci.yml`, or equivalent is present.

### Tracked-artifact debt — unchanged

- **Generated artifacts still dominate the tracked tree** — Doxygen output, a committed Xcode
  build directory, and ~400 MiB of committed run results including three ~90 MB profiler traces.
  The **~84%** figure comes from the Phase-0 reconnaissance and is retained with that evidence
  scope; it was not recomputed here. Phase 1 added **3** tracked files to a tree of ~2,650, so it
  cannot have moved the figure materially.

---

## 6. What this document does **not** claim

- It does **not** claim the rotochemical pipeline is operational. **It is not compiled.**
- It does **not** claim second-order Hartle is validated. It is unreachable and unverified.
- It claims of the O(Ω) solver **only** that its scale-free observable `I = J/Ω` is validated
  (Phase 2B-4B: equation match against published Hartle, analytic and numerical cancellation
  of the arbitrary normalization, agreement with an independent solver to 9.5e-9 analytic /
  2.1e-5 on the CMF sequence, and the correct Newtonian limit). It does **not** claim the
  absolute first-order normalization is correct: `ω̄`, `Ω`, `J`, the `init_omega_bar` seed and
  the `Ω [s⁻¹]` annotation remain **unresolved** (INV-07).
- It does **not** claim placeholder emissivities represent real microphysics.
- It does **not** claim the passive-cooling regression validates the physics. A regression
  baseline now exists (`tests/baselines/passive_cooling_cmf_1p6_debug.tsv`, CTest
  `passive_cooling_regression`, requiring the authenticated external CMF data), and
  `GOVERNANCE.md` §3.1 condition 7 is satisfied — but the neutrino emissivity normalizations it
  freezes are self-labelled placeholders, so it detects change, not correctness.
- It does **not** claim CI exists. The regression must be run by hand.
- It does **not** claim the corrected cooling *trajectory* has been validated. Only the denominator
  identity and the verified `C_⋆(T∞)` are established; the neutrino emissivity normalizations
  remain self-labelled placeholders.
- It does **not** claim the Phase-1B `RotationSolver` repair validated any physics. The library
  compiles; INV-07 and INV-08 are exactly as unresolved as before.
- It does **not** claim a scientific test suite exists. One infrastructure smoke test exists, and
  it asserts no scientific value. Scientific validation is roadmap Phase 2A/2B.
- It does **not** claim `Debug` is a production or scientifically validated configuration, nor
  that `Debug` and `Release` agree numerically.
- It does **not** claim the external ZakiLib/CONFIND repositories constitute cross-platform
  support. Their revisions have not been shown to reproduce the archives consumed here.
