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
   ├─ NeutrinoCooling     LIVE      placeholder Q₀ normalizations
   └─ PhotonCooling       LIVE      L_γ ÷ hardcoded C_eff = 1e40   (INV-15)
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
| `RotochemicalCache` | **NOT COMPILED** | Absent from `Physics/Evolution/CMakeLists.txt`; `Build()` has zero callers |
| `CheckpointObserver` | **EMPTY** | 0 bytes; not in any CMakeLists |

### Drivers

| Driver | Status | Note |
|---|---|---|
| `MagneticDipole` | **LIVE** | `use_moment_of_inertia` is a logged no-op |
| `NeutrinoCooling` | **LIVE** | Emissivity normalizations self-labeled placeholders; `K_PBF = 0.0` |
| `PhotonCooling` | **LIVE** | Hardcoded `C_eff = 1e40` (INV-15) |
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
3. **Heat capacity has two owners** (INV-15) — the two contributions are summed into the same
   state slot using different `C`.
4. **Proper volume defined in three places** (INV-04).
5. **Species semantics undefined** (INV-01, ADR-0001).

---

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
- **Null-deref ordering** — `NeutrinoCooling_Details.cpp:889` dereferences `ctx.star`; the guard
  is at `:901`.

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
