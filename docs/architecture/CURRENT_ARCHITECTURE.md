# CompactStar — Current Architecture

> **STATUS: DESCRIPTIVE.** Authoritative for component boundaries and ownership
> (`GOVERNANCE.md` authority rank 6). Describes **only** behavior that is compiled and reachable.
>
> **Base commit:** `9f70f14` (`master` merge of owner work `3639d71` + April candidate `675b4a9`).
> Evidence for non-`RotationSolver`/`MixedStar` components derives from the `d91c31b` audit; see
> `docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md` §N.
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
| `RotationSolver` — O(Ω) | **LIVE** | Runs on every star build; feeds `SeqPoint::I` |
| `RotationSolver` — O(Ω²) | **UNREACHABLE SCAFFOLDING · CANDIDATE** | `rot_solver` private, no accessor |
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

- **A clean clone cannot configure.** `main/CMakeLists.txt` requires seven gitignored,
  absent directories.
- **Non-macOS builds are impossible.** Zaki and Confind ship only as prebuilt
  `Darwin/{arm64,x86_64}` static archives with no source; configuration hard-fails elsewhere.
- **`configure_file` writes into the source tree**, so no build is truly out-of-source.
- **No warning policy, no default optimization, no sanitizers, no assertions.**
- **No automated tests, no CI.** `main/Test/` is manual demo programs, not a suite.
- **~84% of tracked files are generated artifacts** — Doxygen output, an Xcode build directory,
  and ~400 MiB of committed run results including three ~90 MB profiler traces.

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
