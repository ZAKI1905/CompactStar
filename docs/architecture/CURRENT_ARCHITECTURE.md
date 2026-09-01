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
   └─ PhotonCooling       LIVE      L_γ ÷ constant C_eff = 1e40
                                    NONCONFORMING with ADR-0002 (INV-15)
   ↓
GSLIntegrator (MSBDF, rtol 1e-6 / atol 1e-10)       LIVE
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
| `RotationSolver` — O(Ω) | **LIVE** | Runs on every star build; feeds `SeqPoint::I`. Its profile-backed interpolation path was restored in Phase 1B — see below. **Buildable, not validated** (INV-07) |
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
| `PhotonCooling` | **LIVE · NONCONFORMANT** | Divides `L_γ,∞` by the driver-local constant `C_eff = 1e40` (`PhotonCooling_Details.cpp:320`; `PhotonCooling.hpp:229`) instead of the governed `C_⋆(T∞)`. **Violates ADR-0002** (INV-15). No correction has landed — see below |
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

### `PhotonCooling` — ADR-0002 nonconformance

Recorded for Phase 2A. **Not repaired here.**

- Status is unchanged: **LIVE**. This driver runs in the one program that exercises the full
  pipeline (`main/Test/spin_therm_evol_2_main.cpp`), so the nonconformance is **on the live
  thermal path and affects every result produced to date** — unlike the ADR-0001 nonconformance
  below, which is in code that has never been compiled.
- Under ADR-0002 the sole physical denominator is `C_⋆(T∞)`. `PhotonCooling_Details.cpp:320`
  instead divides by `drv.GetOptions().C_eff`, a driver-local constant defaulting to `1.0e40`
  erg K⁻¹ (`PhotonCooling.hpp:229`) and hand-set at the call site with the comment
  `// Change this!` (`spin_therm_evol_2_main.cpp:245`; also `spin_therm_evol_main.cpp:178`).
- There is **no coupling from `PhotonCooling` to the stellar heat-capacity path at all**:
  `HeatCapacityStar_Tinf` appears nowhere in the driver's four files. `ctx.star` is used only for
  the envelope `Tb` mapping and surface gravity (`PhotonCooling_Details.cpp:153-177`).
- The live equation is therefore `dT∞/dt = −L_γ/1e40 − L_ν/C_⋆(T∞)`: **two different heat
  capacities summed into one state element**, with nothing in the code or the diagnostic output
  signalling the mismatch.
- **This is not resolved in code.** ADR-0002 settles which quantity is authoritative; it
  authorizes no source change. The correction is roadmap **Phase 2A**, is scientific-semantic
  class, and will change numbers the code produces by an amount that must be measured.
- The driver's Doxygen (`PhotonCooling.hpp:55-62`, `:120-123`, `:214-229`;
  `PhotonCooling.cpp:27,36`) documents the constant-`C_eff` equation as the driver's physics and
  becomes wrong on the day the source is corrected.

### `RotationSolver` — Phase-1B merge repair

Recorded because it changes what "LIVE" can mean in this document.

- **Until Phase 1B the library did not compile at all.** Merge `9f70f14` combined owner commit
  `3639d71` with the candidate lineage `e60e656` and resolved three hunks toward the candidate,
  leaving `RotationSolver` internally inconsistent: the header lost the owner's ten profile-backed
  interpolation members and four method declarations while the `.cpp` kept using them;
  `FindNMomInertia` lost the declarations of `R` and `r`, the start-radius scan, the `P`/`E`/`M`
  fetches and its `SetFastProfilePtrs_` call; and seven active lines spliced over a dead
  commented-out function left the file at brace depth 1. Strictly, **no component in this document
  could have been "compiled and reachable" while that held.**
- **Phase 1B (`57334d8`) restored the intended union of both parents.** The owner's first-order
  O(Ω) and MixedStar master-grid machinery from `3639d71` is authoritative and active again; the
  library compiles and links. The first-order numerical path was verified equivalent to
  `3639d71` — its 663-line kernel is byte-identical and the repaired `FindNMomInertia` differs
  from the owner's by additions only, all of them write-only storage.
- **Nothing was validated by this.** The repair was engineering class. The second-order candidate
  region is byte-identical to before — no equation, `j²`, `δM`, `dε/dp`, or boundary condition was
  touched — it remains unreachable (`rot_solver` is private, `SolveHartle2_N` and
  `GetHartleResult` have no external callers, INV-08 unchanged), and `init_omega_bar` was not
  touched, so **INV-07 remains unresolved**.

**Buildable and reachable is not numerically validated.** The `LIVE` labels in this document now
mean what they say; they still make no claim about correctness.

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

- **`GeometryCache` has no version gate.** Deep-copies at construction, never rebuilt. If a
  profile is re-solved afterwards, downstream integrals silently use stale geometry (INV-12).
- **`StarContext` binds raw column pointers once** and never re-binds. A reallocating profile
  mutation bumps the version — payloads rebuild — while the seven cached pointers dangle.
- **No driver dependency graph.** `IDriver::DependsOn()` / `Updates()` are declared, overridden
  by every driver, and **never called**. `EvolutionSystem.cpp:125-134` iterates in registration
  order. The comment at `IDriver.hpp:102-103` referring to an "Evolution graph" is stale.
- **Declared-but-undefined symbols reachable from compiled code** —
  `Physics::Spin::DipoleFieldEstimate` and `CharacteristicAge` are called at `Pulsar.cpp:204,212`
  with no definition anywhere.
- **Null-deref ordering** — `NeutrinoCooling_Details.cpp:889` dereferences `ctx.star` (the
  `HeatCapacityStar_Tinf` call); the `if (!ctx.star)` guard is twelve lines later at `:901` and
  cannot fire. **Confirmed present at `ba49e10`.** Engineering-class defect, tracked separately
  from the INV-15 heat-capacity decision, but scoped into Phase 2A because routing `PhotonCooling`
  through the same context path exercises the same unguarded pattern.
- **Heat-capacity cache key omits the geometry** — `StarContext::HeatCapacityStar_Tinf` accepts an
  optional `GeometryCache` and falls back to a locally constructed one
  (`StarContext.cpp:754-755`), but keys its cache only on `(profile version, thermo pointer)`
  (`:712-714`). A later call at the same profile version with a different `GeometryCache` silently
  reuses the earlier table (INV-12).

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
- **No sanitizer policy and no assertion policy exist.** Neither has been established; the
  repository still contains no assertions.

### Automated tests

- **CTest infrastructure exists.** `include(CTest)` at top level provides standard `BUILD_TESTING`
  and `enable_testing()`; `tests/` is added only when testing is on (`CMakeLists.txt:222-226`).
- **`tests/` is the canonical automated-test root**, holding exactly **one** test:
  `compactstar_library_smoke` (`tests/CMakeLists.txt:19`;
  `tests/smoke/compactstar_library_smoke.cpp`).
- **`main/Test/` remains manual demo/debug programs**, not a suite. None of its eight executables
  is registered with CTest.
- **The smoke test is infrastructure validation, not a scientific baseline.** It establishes that
  CompactStar headers compile, that a test executable links a genuine out-of-line symbol from
  `libCompactStar.a` (`CompactStar::Core::Prog::GetName()`, verified present as a defined `T`
  symbol rather than an unresolved import), that transitive dependencies link, and that the binary
  runs and exits zero under CTest. **It asserts no scientific value of any kind.**
  There is no scientific test suite.
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
- It does **not** claim the O(Ω) solver is numerically correct. It is live and untested, and
  its normalization is unresolved (INV-07).
- It does **not** claim placeholder emissivities represent real microphysics.
- It does **not** claim the heat-capacity inconsistency is resolved in code. ADR-0002 decides the
  physical owner; `PhotonCooling` still divides by a constant, and no source correction has landed.
- It does **not** claim `StarContext::HeatCapacityStar_Tinf` is numerically validated. It is the
  designated implementation of the governed quantity, and it is untested.
- It does **not** claim the Phase-1B `RotationSolver` repair validated any physics. The library
  compiles; INV-07 and INV-08 are exactly as unresolved as before.
- It does **not** claim a scientific test suite exists. One infrastructure smoke test exists, and
  it asserts no scientific value. Scientific validation is roadmap Phase 2A/2B.
- It does **not** claim `Debug` is a production or scientifically validated configuration, nor
  that `Debug` and `Release` agree numerically.
- It does **not** claim the external ZakiLib/CONFIND repositories constitute cross-platform
  support. Their revisions have not been shown to reproduce the archives consumed here.
