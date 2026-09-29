# Phase-6 external-dependency (ZakiLib 2.0 / CONFIND 2.0) and TaskManager stellar-equivalence migration — reconnaissance and predeclaration

**Status:** PREDECLARED — RECONNAISSANCE COMPLETE — OWNER REVIEW REQUIRED — NO IMPLEMENTATION AUTHORIZED

**Disposition (§32):** **C — ZAKI / TASKMANAGER NUMERICAL DRIFT REQUIRES A SEPARATE PATCH BEFORE
COMPACTSTAR MIGRATION.**

**Date:** 2026-09-28. **Platform:** local Mac only. No Linux host or EKU cluster was accessed.

**Change class of this record:** documentation. The migration it predeclares is dependency/build +
structural architecture, contains an engineering-class plotting deletion, and (per §26) requires a
dependency-side numerical-preservation change in ZakiLib before CompactStar may consume it.

**Execution authority:** none. This file is the only repository change. No CMake, production
source, test, baseline, dependency, EOS/data, or literature byte is changed. Zaki and CONFIND are
unchanged. All builds, probes and prototype runs described below were performed in a
session-scratch directory outside every repository; they are reconnaissance observations, not
governed qualification (§31).

Placement: `docs/governance/` does not exist in this repository. Dependency/cluster-track records
live in `docs/validation/PHASE6_*` (for example `docs/validation/PHASE6_ADR0018_ACCEPTANCE.md`,
`docs/validation/PHASE6_LINUX_CLUSTER_BOOTSTRAP_PREFLIGHT.md`), and predeclarations use the
`_PREDECLARATION` suffix (`docs/validation/PHASE5D1_CONTROLLED_EVOLUTION_PREDECLARATION.md`). This
file follows that convention.

---

## 0. Findings that change the migration plan

1. **Canonical Zaki `c8c6813` is not arithmetic-equivalent to the historical vendored `libZaki.a`,
   in Debug as well as Release.** The vendored arm64 archive contains **zero** fused
   multiply-add instructions; canonical Zaki built by its own CMake with AppleClang 21 defaults
   (`-ffp-contract=on`) contains 78 (Debug) / 120 (Release): `Segment::GetIntersection` (all three
   overloads) and `Curve2D::GetIdx` in both modes, plus `Curve2D::Intersection`,
   `Curve2D::Bisect` and `Segment::Intersects` in Release. A differential
   probe over 162,619 values found 321/326 intersection x-coordinates and 319/326 y-coordinates
   different (≤ 531 ULP), and 774–1694 of 2,400 `IntegralTable` `I_1/I_2/I_3` values different,
   in **both** Debug and Release (§11).
2. **The drift reaches TaskManager's final scientific outputs on a real stellar grid.** On a real
   20×20 DS(CMF)-1 + Fermi-gas mixed-star grid (400 TOV solves, then `Precision_Task` and
   `FindLimits`), the unpatched canonical stack reproduced the dark EOS, the 400-row sequence, the
   critical curve and the M_tot intersection byte-for-byte, but changed the `Bisect` cut of B_tot
   contour 4 (one point fewer in `B_tot_2.01_4.tsv`). That shifted the stride-10 contour points
   `Precision_Task` solves, so **6 of 38 text outputs differ**, including all three final
   `BNV_tau_*` limit files, by up to **5.2 % relative** (neutron), 2.4 % (Λ) and 0.39 % (Σ⁻) (§15).
3. **FMA contraction is the sole Debug-mode cause, and it is removable without source change.**
   Canonical Zaki compiled with `-ffp-contract=off` reproduced the vendored oracle bit-for-bit on
   all 162,619 probe values (whole-output SHA-256 identical) and reproduced every text output of
   the real stellar prototype byte-for-byte (§11, §15).
4. **Release adds libm-call rewriting that contraction control does not fix.** Canonical Zaki
   Release imports `___exp10` in `Axis::operator[]` and (inlined) `GridVals_2D::Interpolate`
   (87/46,002 random log-axis nodes differ by 1 ULP; 0 on every historical TaskManager grid and
   thread partition tested), and rewrites `pow(x,2)` into `x*x` in `IntegralTable` (2/2,400 `I_2`
   values, 1 ULP) even under `-fno-builtin` (§11).
5. **A clean preservation design exists and was demonstrated in scratch.** A scratch-only
   prototype of the proposed ZakiLib 2.0.1 change (§26: private `-ffp-contract=off` plus a
   noinline `-fno-builtin -fno-lto` historical `pow` boundary for `Axis::operator[]`, the
   `Quantity` printer and `IntegralTable`) reproduced the vendored oracle bit-for-bit in **Debug
   and Release** (all 162,619 probe values) and reproduced **all 38 T2 text outputs** in both
   modes; ZM-1's 12 tests still pass.
6. **CONFIND 2.0 introduces no TaskManager-visible numerical drift in the evidence gathered**, and
   the historical arm64 `ContourFinder::Plot` is proven a stack-only no-op; `SetPlotConnected`
   writes two fields that no historical code reads (§8). Its N-3 long-path change does alter
   TaskManager output **file names** under long working roots (demonstrated, §15).
7. **Coord3D:** in every existing CompactStar executable the comparators are not linked at all; in
   a TaskManager-bearing OLD link they resolve to vendored `libConfind.a(Cont2D.cpp.o)`; in the
   prospective NEW link to canonical `libCONFIND.a(Cont2D.cpp.o)` (Debug) or nowhere (Release,
   inlined). A default-flags consumer TU using `std::set<Coord3D>` was shown to capture all three
   symbols and link an FMA-fused `XYDist2` (§13).
8. **Pre-existing TaskManager defects bound the fixture design:** its multi-thread path races on a
   static sequence guarded by a per-object mutex, and recomputes per-thread axis endpoints so that
   stellar outputs depend on thread count (demonstrated: 4 threads changed one stellar mass).
   Governed TaskManager authority must be serial (§22).
9. Compiled CompactStar C++ uses no Python/NumPy directly; the linkage exists only for the
   vendored Zaki plotting bridge, and its `-undefined dynamic_lookup` currently masks an unlinked
   zlib (`_compress` is a flat-namespace lazy bind). CTest still needs a Python interpreter with
   numpy/scipy/mpmath (§9).

---

## 1. Authenticated identities

### 1.1 Repositories (all tracked trees clean)

| Repository | Path | HEAD = local `master` = `origin/master` = live `master` |
|---|---|---|
| CompactStar | `/Users/keeper/Documents/CompactStar/repo/CompactStar` | `812463ac9ed374f64ac9cadd500066ab723d3a6c` |
| ZakiLib | `/Users/keeper/Documents/CompactStar/external/Zaki` | `c8c68131b04e9d216673725075bef38df81e6041` |
| CONFIND | `/Users/keeper/Documents/CompactStar/external/CONFIND` | `b0cbd510fd3fd0c772fa50499cd749287cb39e7b` |

Live values were read with `git ls-remote origin refs/heads/master`. Migration branch
`physics/external-deps-taskmanager-migration`, worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration`,
created fresh from `812463a`. The unmerged ADR-0018 experiment branch
`physics/adr0018-mac-source-equivalence` (`abf2b62`, four-path override CMake at `fc55386`) stopped
at its candidate-build gate and is superseded by this plan; it must not be merged.

Ancestry audit (local and `origin` branches ahead of `master`): no branch touches
`TaskManager.{hpp,cpp}`, `TOVSolver_Thread.{hpp,cpp}`, `dependencies/` or `main/CMakeLists.txt`;
only the superseded ADR-0018 experiment branch touches root `CMakeLists.txt`; eleven non-canonical
Phase-6A-1 branches (the historical 77-test topology,
`PHASE6_LINUX_CLUSTER_BOOTSTRAP_PREFLIGHT.md:454-464`) touch `tests/CMakeLists.txt`. The M1 test
registration must be reconciled with any of those only if the owner later integrates them.

### 1.2 Historical vendored artifacts (CompactStar `812463a`)

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `dependencies/lib/Zaki/Darwin/arm64/libZaki.a` | 2,478,624 | `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef` |
| `dependencies/lib/Confind/Darwin/arm64/libConfind.a` | 703,208 | `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2` |
| `dependencies/lib/Zaki/Darwin/x86_64/libZaki.a` | 1,079,232 | `45e2d9c8d784a601252438c8021866b996626289aa30a30f368e2a34786ab09f` |
| `dependencies/lib/Confind/Darwin/x86_64/libConfind.a` | 903,904 | `869fd9232ed1e1c24d227a2282392d612d7c6a0a119263275d887fd229794520` |
| `dependencies/include/Zaki/**/*.hpp` (31 files) | — | manifest `02e48f0cd309148c3680f77d7fe17435997aac8c665ddf9d99bd90c254e60263` |
| `dependencies/include/Confind/*` primary (7 files; the six `* 2.hpp` are byte-identical duplicates) | — | manifest `c7852fcb49edc43d0245c27ee3445311360b1c76b5351ae6381ea35bba87cce8` |
| `dependencies/include/matplotlibcpp.hpp` | — | `60fd88a9b631be394a5b5cbe415105240f7a0bcb5887b7d1ae90cc9d71c2b3d6` |

Manifest construction: `LC_ALL=C`-sorted relative paths, one `shasum -a 256 <path>` line each,
SHA-256 of that text. The arm64 hashes equal ADR-0018 §2 and `PHASE6_ADR0018_ACCEPTANCE.md:70-71`.
All 31 vendored Zaki headers are byte-identical to Zaki `b9ddeba`; seven differ from `c8c6813`
(`Func2D`, `GSLMultiFWrapper`, `Math_Core`, `MemFuncWrapper`, `Newton`, `Logger`, `DataSet`; §6).

Binary facts established here:

- arm64 `libZaki.a`: `LC_BUILD_VERSION` minOS/SDK 13.0, no DWARF, `-O0` (every value spilled),
  asserts live (`___assert_rtn` in `Axis::operator[]`), **0 FMA-class instructions** in all 25
  members, `Axis::operator[]` Log branch calls libc++ `pow<int,double>` → libm `_pow`; no Coord3D
  comparator definitions.
- arm64 `libConfind.a`: 0 FMA; `Cont2D.cpp.o` defines weak external `Coord3D::operator<`,
  `operator==`, `XYDist2`; `ContourFinder::Plot` is 8 instructions storing its arguments to its own
  stack frame; `SetPlotConnected` stores to object offsets `0x188` and `0x47`, which only the
  constructor, copy constructor and copy assignment otherwise touch.
- x86_64 archives are **not** equivalent Mac authority: optimized `libZaki.a` (`___exp10` in
  `Axis::operator[]` and constant initializers), different source revision, and an exported
  `DataColumn::pow(double const&)` that mismatches the vendored header; x86_64 `Plot` calls
  `ConvertToDataSet` → `SortNew`, re-sorting contours in place before `GetContourSet`. No governed
  evidence used them (all governed evidence is arm64).

### 1.3 Reference toolchain (this Mac)

macOS 26.6.2 (25G83), Darwin 25.6.0, arm64; Apple clang 21.0.0 (`clang-2100.3.34.2`) at
`/usr/bin/clang++`; SDK 27.0; CMake 4.2.1 (`/opt/homebrew/bin/cmake`); GSL 2.7.1 at `/opt/local`;
libomp `/opt/local/lib/libomp/libomp.dylib`; Python 3.12.10 (`/Users/keeper/miniforge3/bin/python3`)
with NumPy 2.3.1 and matplotlib 3.10.3. This equals the ADR-0017 qualification platform
(`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:29-45`).

Discovery hazard recorded for G0: CompactStar's configure resolves GSL headers/libraries at
`/opt/local` but `GSL_CONFIG_EXECUTABLE=/Users/keeper/miniforge3/bin/gsl-config`, so CMake reports
"found version 2.7" for a 2.7.1 installation. Canonical Zaki's own configure (with
`-DGSL_ROOT_DIR=/opt/local`) reports 2.7.1.

---

## 2. Governing authority read

`GOVERNANCE.md`; `AGENTS.md`; ADR-0018 in full; `docs/validation/PHASE6_ADR0018_ACCEPTANCE.md`;
`docs/validation/PHASE6_LINUX_CLUSTER_BOOTSTRAP_PREFLIGHT.md`; the unmerged
`PHASE6_ADR0018_MAC_SOURCE_EQUIVALENCE.md` result (stopped at build gate);
`docs/build/MACOS_BUILD.md`; `docs/SCIENTIFIC_INVARIANTS.md` header and INV-13;
`docs/validation/PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md` (platform, protected bytes);
`docs/validation/PHASE5B_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md` and
`PHASE5C2_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md`; `docs/adr/README.md`. ZakiLib:
`docs/modernization/ZAKILIB_ZM1_PYTHON_FREE_CORE.md`, `ZAKILIB_ZM1_ACCEPTANCE.md`,
`ZAKILIB_PLOTTING_CONSUMER_MIGRATION.md`. CONFIND: `CONFIND_2_0_ACCEPTANCE.md`,
`COMPACTSTAR_MIGRATION.md`, `CONFIND_2_0_MODERNIZATION.md`, `CONFIND_2_0_REPORT.md`,
`CONFIND_2_0_M1_CORRECTION.md`, `evidence/m1/INDEPENDENT_REVIEW.md`, and the relevant sections of
`CONFIND_RECONNAISSANCE_PREFLIGHT.md` (§§1–4, 8–18, 22–39).

Governing constraints carried forward:

- ADR-0018 §3 clause 1 keeps the vendored Darwin archives as the default Mac authority, clause 11
  places replacing them outside ADR-0018, and §5 fixes the equivalence rule: deterministic
  control/treatment scientific fields at **0 ULP**, no post-result tolerance.
- The Phase-5B/5C portability ratifications admit exactly one non-equality field (the compiler
  string) and state that a toolchain change never excuses scientific drift.
- CONFIND acceptance mandates, for this migration, N-1 Coord3D link qualification, log-axis
  `pow`/`exp10` characterization outside CONFIND, TaskManager-native `pow` characterization,
  preservation of raw contour order, `mass_curve[0]`, intersection order and `Bisect(...).first`,
  and real-stellar evidence (`CONFIND_2_0_ACCEPTANCE.md` "Mandatory future CompactStar
  qualification" 1–7; `COMPACTSTAR_MIGRATION.md` 1–12).
- ZM-1 froze Zaki constants, linear interpolation, `DataColumn` semantics, `Coord3D` ordering,
  `DataSet::Solve`, default export precision and import coercion
  (`ZAKILIB_ZM1_PYTHON_FREE_CORE.md:199-213`).

---

## 3. Scope and non-scope

In scope (later, after owner acceptance): replace vendored Zaki/CONFIND headers and archives in
CompactStar's active build by the canonical installed packages `Zaki::Zaki` and `CONFIND::CONFIND`;
delete CompactStar's plotting call sites; remove Python/NumPy from the compiled build; qualify
same-Mac equivalence of every governed and TaskManager result.

Not in scope: changing any CompactStar equation, tolerance, baseline, EOS input or test inventory
beyond the predeclared T1/T2 additions; changing CONFIND semantics (Coord3D/radius
de-duplication, `RMDuplicates`, `SortNew`, chaining, saddle, NaN, exact equality, grid indexing,
curve order); changing Zaki semantics; repairing TaskManager's pre-existing defects (§22); CompactStar
Release-vs-Debug authority (§21); Linux and EKU cluster qualification (§30).

---

## 4. Frozen authorities

| Authority | Frozen object | Source |
|---|---|---|
| Mac dependency authority | vendored arm64 `libZaki.a` `3dd4789a…`, `libConfind.a` `09ed1a7c…`, vendored headers (§1.2) | ADR-0018 §3.1; `PHASE6_ADR0018_ACCEPTANCE.md:64-74` |
| Governed baselines | the 11 files and SHA-256 values of `PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:79-93` (re-verified by the test reconnaissance) | ADR-0017 record |
| CONFIND numerical semantics | historical serial CONREC, Coord3D radius de-duplication, `RMDuplicates`, `SortNew`, chaining, saddle, NaN, exact equality, indexing, order | `CONFIND_2_0_ACCEPTANCE.md` |
| Zaki semantics | ZM-1 frozen list | `ZAKILIB_ZM1_PYTHON_FREE_CORE.md:199-213` |
| TaskManager | **no baseline exists.** Its authority is the behaviour of the OLD build (§18) on the T1/T2 fixtures (§16), captured before any migration change | this record |
| Build authority | Debug, AppleClang 21 `clang-2100.3.34.2`, effective flags `-g -std=c++17 -arch arm64 -fPIC -pthread -Wall -Wextra -Xclang -fopenmp` (library) | `PHASE6_ADR0017_PRODUCTION_IMPLEMENTATION.md:35-44` |

---

## 5. Current dependency consumption (CompactStar `812463a`)

Link wiring: `CMakeLists.txt:76` GSL, `:79` OpenMP, `:82`
`find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)`, `:92-101` archive paths
from `CMAKE_SYSTEM_NAME`/`CMAKE_HOST_SYSTEM_PROCESSOR` with FATAL_ERROR if absent, `:146-150`
SYSTEM includes (`${DEP_INC_DIR}` = `dependencies/include`), `:168-181` PUBLIC link order
`libZaki.a → libConfind.a → GSL::gsl → GSL::gslcblas → Python3::Python → Python3::Module →
Python3::NumPy → OpenMP::OpenMP_CXX`. `main/Examples` and `main/Test` also add `${DEP_INC_DIR}`.
Only `Darwin/{arm64,x86_64}` archives exist; the *host* processor selects them; no Linux path
exists.

Observed OLD executable link line (tests): `... -Xlinker -undefined -Xlinker dynamic_lookup <obj>
../libCompactStar.a <vendored>/libZaki.a <vendored>/libConfind.a /opt/local/lib/libgsl.dylib
/opt/local/lib/libgslcblas.dylib /Users/keeper/miniforge3/lib/libpython3.12.dylib
/opt/local/lib/libomp.dylib`. `Python3::Module` injects `dynamic_lookup`. No zlib is linked;
`_compress` (imported by vendored `DataSet.cpp.o`, `Math_Core.cpp.o`, `Cont2D.cpp.o` through the
inline `VecSaver`) remains a flat-namespace lazy bind in every executable.

| Dependency | Classification | Evidence |
|---|---|---|
| Zaki headers/symbols | DIRECT_COMPACTSTAR_USE + SCIENTIFIC_NUMERICAL; 201 out-of-line symbols referenced by libCompactStar objects (55 `DataSet`, 27 `DataColumn`, 16 free `Vector` ops, 13 `Curve2D`, 4 `GridVals_2D`, 1 `Axis`, 2 `Segment`, 3 `I_1/I_2/I_3`, 43 `Physics` constants, String/File/Util) | `nm -u` of the 83 OLD objects intersected with vendored `libZaki.a` definitions |
| CONFIND | DIRECT_COMPACTSTAR_USE (TaskManager only) + SCIENTIFIC_NUMERICAL; 14 symbols | only `TaskManager.cpp.o` references CONFIND |
| Python / NumPy | TRANSITIVE_ZAKI_USE + PLOTTING_ONLY + BUILD_LEGACY (link); CTest Interpreter use is test tooling | §9 |
| OpenMP | BUILD_LEGACY for CompactStar (no pragma, no `omp_*`, no `_OPENMP` test in CompactStar code); TRANSITIVE_CONFIND_USE (vendored `libConfind.a` Parallel/Ludicrous modes, unused by TaskManager) | threading audit; `-Xclang -fopenmp` is FP-arithmetic-neutral at `-O0` across 32,041 functions (§20) |
| zlib | TRANSITIVE_ZAKI_USE / TRANSITIVE_CONFIND_USE (inline `VecSaver` compression); **not linked explicitly today**; masked by `dynamic_lookup` | `nm -u`, `dyld_info -fixups` |
| GSL | DIRECT_COMPACTSTAR_USE + SCIENTIFIC_NUMERICAL (TOV splines, integrators) and TRANSITIVE_ZAKI_USE (`gsl_spline`, `gsl_spline2d` in DataSet/GridVals) | `CMakeLists.txt:76,173-174` |

No existing executable links TaskManager or CONFIND: every `main/` program and every test was
checked with `nm` (0 `CONFIND::`, 0 `TaskManager::` symbols).

---

## 6. Zaki API surface versus canonical 2.0

| API family consumed by CompactStar | Classification | Notes |
|---|---|---|
| `Physics` constants (43 referenced; all 132 probed) | UNCHANGED (bitwise, Debug and Release) | probe CONST 0/132 |
| `DataSet` construction, Import/Export, `Interpolate`, `Evaluate`, `Derivative`, `Integrate`, `MakeSmooth`, `Col`, `AddColumn`, etc. | UNCHANGED on the probed valid path; ADDITIVE_COMPATIBLE lifetime fixes (accelerator/resource release, negative-index normalization); **ABI_CHANGED layout** (the `PlotParam plt_par` member is removed, `DataSet.hpp` diff) | full recompile against canonical headers is mandatory; no object may mix header generations |
| `DataSet` plotting (`Plot`, `LogLogPlot`, `SemiLogXPlot`, `SemiLogYPlot`, `PlotParam`, `SetPlotPars`, `ResetPlotPars`) | PLOTTING_REMOVED | §8 |
| `DataColumn` arithmetic, `exp`, `log10`, `pow` | UNCHANGED (probe DC_* 0/2,382) | |
| `Axis::operator[]` | Linear branch UNCHANGED; Log branch **NUMERICALLY_SENSITIVE** in Release (`___exp10`, 87/46,002 nodes); BEHAVIORALLY_CHANGED_ERROR_PATH in Release (`assert(i <= res)` compiled out under `NDEBUG`) | canonical `Zaki/Math/src/Math_Core.cpp:1504-1513`; historical `b9ddeba` `:1433-1443` |
| `GridVals_2D` ctor/`Interpolate`/`Evaluate` | NUMERICALLY_SENSITIVE (Release `___exp10` inlined into `Interpolate`); ADDITIVE_COMPATIBLE (deep copy/move, `delete[]`, re-`Interpolate` cleanup); bilinear GSL `gsl_spline2d` unchanged | canonical `:1397`, `:1437`; CompactStar never copies a `GridVals_2D` |
| `Segment::GetIntersection`, `Curve2D::Intersection`, `Curve2D::GetIdx`, `Curve2D::Bisect`, `Segment::Intersects` | **NUMERICALLY_SENSITIVE** (FMA in Debug for `GetIntersection`/`GetIdx`; additionally `Intersection`/`Bisect`/`Intersects` in Release) | canonical `:789-848`, `:1026-1100` |
| `Segment::P`, `Curve2D::MakeSmooth`, `Curve2D::Export/Import`, `Coord2D` | UNCHANGED (probe) | export via `VecSaver::Export1D<Coord2D>` whose only FMA sizes a zlib buffer |
| `Curve2D::Plot` (both overloads) | PLOTTING_REMOVED | |
| `Coord3D` (`operator<`, `operator==`, `XYDist2`) | UNCHANGED source (canonical `Zaki/Physics/Coordinate.hpp:61-80`); link-resolution sensitive (§13) | |
| `IntegralTable` `I_1/I_2/I_3` | **NUMERICALLY_SENSITIVE** (FMA both modes; Release `pow(x,2)`→`x*x`) | used by `CompactStar/EOS/src/Baryon.cpp:214,245,252,315,333`, `Particle.cpp:160,167`; linked only into `main/Examples` programs, not into any test |
| GSL wrappers (`GSLFuncWrapper`, `GSLMultiF*`), `Newton`, `Func2D`/`MemFuncWrapper` | UNCHANGED (header include hygiene only; `Newton` `std::abs`) | |
| `Directory`, String utilities, `CSVIterator`/`CSVRow`, `VecSaver` | UNCHANGED on valid paths; `Directory::Create` now `std::filesystem::create_directories` (terminates on relative paths) | |
| `Logger`, `Instrumentor` | UNCHANGED semantics (member-initializer reorder only); historical ODR layout mismatch between vendored `libConfind.a` and `libZaki.a` disappears | CONFIND reconnaissance §18 item 9 |
| `Range::LenAbs` (`std::abs`) | UNCHANGED for CompactStar (0 uses) | |
| MISSING | **none** — all 74 non-plotting CompactStar TUs and a plotting-free TaskManager compile against canonical headers with 0 errors | §15 |

---

## 7. CONFIND API surface

The only caller is `CompactStar/Core/src/TaskManager.cpp`. Symbols referenced at `812463a`:
`ContourFinder()`, `~ContourFinder()`, `SetGrid(Grid2D const&)`, `Base::SetWrkDir(Directory const&)`,
`SetContVal(vector<double> const&)`, `SetGridVals(GridVals_2D*)`, `GetContourSet() const`,
`ExportContour(Directory const&, FileMode const&)`, `Cont2D::size() const`, `Cont2D::GetVal() const`,
`Cont2D::operator[](unsigned long const&) const`, `Cont2D::ConvertToCurve2D()` — all **UNCHANGED**
and defined by canonical `libCONFIND.a` (Debug and Release, verified with `nm`) — plus
`SetPlotConnected(bool)` and `Plot(Directory const&, char const*, char const*, char const*)` —
**PLOTTING_REMOVED**. No THREADING_API use (TaskManager passes precomputed `GridVals_2D`; the new
`Evaluate(EvaluatorFactory, ThreadingOptions)` is not used). BUILD_ONLY: `find_package(CONFIND 2.0
CONFIG)`, `Threads::Threads`. Minimal post-migration API: the 12 unchanged symbols above.

Error-path changes relevant to TaskManager (CONFIND N-3, documented, non-numerical): empty found
contour → empty `Curve2D` (historically `SortNew` read `pts[0]`: UB); `SetWrkDir` recursive, may
throw; `SetGrid` rejects zero resolution/unknown scale; export paths longer than the historical
149-character payload are no longer truncated (§15).

---

## 8. Plotting inventory (complete; `812463a`)

Independent inventory over tracked `CompactStar/`, `main/`, `tests/`; the prior ZakiLib ledger's
224 lines match HEAD exactly, but its "active" rule (not line-commented) misses `#if 0`.

| File | Active plot calls | Inactive | Active setup (`PlotParam`/`SetPlotPars`/`ResetPlotPars`) | Commented setup | `SetPlotConnected` active/commented |
|---|---:|---:|---:|---:|---|
| `Core/src/TaskManager.cpp` | 5 | 9 (incl. `:455` in `#if 0`) | 0 | 0 | 3 / 1 |
| `EOS/src/CompOSE_EOS.cpp` | 9 | 1 | 16 | 1 | — |
| `Extensions/MixedStar/src/DarkCore_Analysis.cpp` | 18 | 0 | 0 | 0 | — |
| `Microphysics/BNV/Analysis/src/BNV_Analysis.cpp` | 3 | 18 | 4 | 0 | — |
| `Microphysics/BNV/Analysis/src/BNV_Sequence.cpp` | 8 | 5 | 12 | 5 | — |
| `Microphysics/BNV/Analysis/src/Decay_Analysis.cpp` | 6 | 8 | 0 | 4 | — |
| `Microphysics/BNV/Channels/src/BNV_B_Chi_Photon.cpp` | 10 | 3 | 20 | 4 | — |
| `Microphysics/BNV/Channels/src/BNV_B_Chi_Transition.cpp` | 2 | 1 | 4 | 1 | — |
| `Microphysics/BNV/Internal/src/BNV_Chi.cpp` | 14 | 5 | 26 | 2 | — |
| **Total** | **75** | **50** | **82** | **17** | **3 / 1** |

By API (active): 70 `DataSet`, 2 `Curve2D::Plot` (`TaskManager.cpp:415` member, `:623` static),
3 CONFIND `Plot` (`:347`, `:519`, `:628`). 238 active and 107 commented `plt_par.*` member calls go
away with their sites. `main/Examples/Table_5-8_Glenn.cpp:31` and `rotating_ns.cpp:32` include
`matplotlibcpp.hpp` directly but are not built and use nonexistent CompactStar types.

Every active call sits in a TU compiled into libCompactStar and **none is reachable** from any
tracked executable or registered test. Historical behaviour (`b9ddeba` proxy + vendored binaries):
the DataSet plot overloads and `PlotParam::Use` read columns through const accessors and write no
`DataSet` data, spline, accelerator, precision or working-directory state; `SetPlotPars` mutates
only `plt_par`, which only plot methods read; `Curve2D::Plot` reads copies from const
`GetXVals/GetYVals`. Side effects are limited to file creation, the embedded Python interpreter,
and throw/abort paths. A scratch OLD-linked run of both `Curve2D::Plot` overloads returned
normally, wrote PDFs, and left the curve's bit hash, FPCR (`0`) and rounding mode unchanged.

Recommendations: **DELETE 43** (data already exported by the same path, or input file), **EXPORT_DATA
31** (no current export of the plotted derived data; add a numerical export or accept loss of the
figure — owner-visible, non-numerical), **OWNER_REVIEW 1** (`TaskManager.cpp:347`: arm64 `Plot` is a
no-op, so deletion reproduces the Mac authority; x86_64 `Plot` would re-sort before
`GetContourSet`, so deletion is *not* x86_64-historical behaviour — the review is to confirm the
arm64 authority). `TaskManager.cpp:519`, `:628` and the three `SetPlotConnected` calls are provably
DELETE. Items adjacent to plot blocks that must be kept (they feed exports):
`BNV_Sequence.cpp:440,1491-1494,1530` (`MakeSmooth(5)` changes exported data),
`DarkCore_Analysis.cpp:370,392,412`, `Decay_Analysis.cpp:549,557`, `BNV_Analysis.cpp:265,316,361`,
`BNV_B_Chi_Photon.cpp:787-789,1119-1121,1256,1769-1771,1905`, `BNV_Chi.cpp:946,1302,1381,1545-1547`,
`CompOSE_EOS.cpp:708,714,720`.

Removal staging: see §20. No plotting API is recreated in Zaki or CONFIND.

---

## 9. Python / NumPy

- Compiled CompactStar C++: no `Python.h`, `PyObject`, `Py_*`, NumPy or matplotlibcpp in any built
  TU (`CompOSE_EOS.cpp:11` include is commented; the two direct-matplotlibcpp examples are unbuilt).
- The link exists only because vendored `libZaki.a` embeds matplotlibcpp (`DataSet.cpp.o` 29,
  `Math_Core.cpp.o` 25 undefined `_Py*` symbols; `libConfind.a` none). `Python3::Module` adds `-undefined dynamic_lookup`,
  which hides unresolved symbols (the unlinked zlib today).
- CTest: 18 of 76 registrations run Python scripts through `${Python3_EXECUTABLE}`
  (`tests/CMakeLists.txt:26,32,41,55,77,88,599,610,616,620,639`); packages: numpy, scipy, mpmath
  (recorded stack Python 3.12.10, NumPy 2.3.1, SciPy 1.16.0, mpmath 1.4.1; no lock file).

**Conclusion: YES — after plotting removal and package consumption, the compiled CompactStar
targets can be Python/NumPy-free.** `find_package(Python3 COMPONENTS Interpreter)` remains for the
Python-driven tests only. Python scripts under `tests/`, `docs/`, and
`main/Test/results/spin_therm_evol_2/plot.py` remain valid external utilities and are not deleted.

---

## 10. Build authority

Governed evidence (ADR-0017 production qualification, Phase-5B/5C/5D integrations on
AppleClang 21) is **Debug** CompactStar (`-g`, no `-O`, so `-O0`) linked to the **`-O0`,
asserts-live, FMA-free vendored archives** — a mixed configuration, since CompactStar's own code
is compiled by AppleClang 21 with its default `-ffp-contract=on` (e.g. 4–5 `fmadd` in
`TaskManager.cpp.o` at `-O0`). Baselines that record configuration: B7 `"Debug"`, B8
`"assertions-enabled"`; `*_debug.tsv` names. No governed baseline records Zaki, CONFIND or Python
identity, so the dependency swap needs no new provenance allowlist.

"Same-Mac equivalence" therefore means: same Mac, same AppleClang, same CompactStar compile and link
flags, same GSL/zlib runtime, same inputs; only the dependency identity differs; every
deterministic scientific field and output byte identical, per build mode, OLD versus NEW.

---

## 11. Zaki numerical differential (pow/exp10 and FMA)

Probe (scratch): one C++17 TU compiled identically (`-O0 -ffp-contract=off`) against (a) vendored
headers + vendored arm64 `libZaki.a` (OLD), (b) canonical c8c6813 Debug, (c) canonical Release,
plus contraction-off variants. All probe-side inputs use explicit libm calls in the probe TU, so
only out-of-line Zaki code differs. Output: one record per value, `%a` and IEEE bits. OLD output is
run-to-run identical (SHA-256 `376bc3fd…e1d0`, 162,619 records).

| Section (records) | NEW Debug | NEW Release | Debug `-ffp-contract=off` | Release `-ffp-contract=off` | Release `-ffp-contract=off -fno-builtin` |
|---|---|---|---|---|---|
| Axis nodes, historical TaskManager grids (1,575) | 0 | 0 | 0 | 0 | 0 |
| Axis nodes, TaskManager thread partitions T=1..10 (1,057) | 0 | 0 | 0 | 0 | 0 |
| Axis nodes, 300 random log axes (46,002) | 0 | **87 (1 ULP)** | 0 | **87** | 0 |
| GridVals_2D bilinear, Log 99×150 / Lin 99×99 / Log 49×49 at nodes, midpoints, 4,000 random points each (66,652) | 0 | 0 | 0 | 0 | 0 |
| `Curve2D::Intersection` x / y (326 each) | **321 / 319** (≤176 / ≤531 ULP) | **321 / 319** | 0 | 0 | 0 |
| `Segment::GetIntersection` x / y (480 each) | **479 / 479** (≤3,564 / ≤318,976 ULP) | **479 / 479** | 0 | 0 | 0 |
| intersection count, `Intersects`, `GetIdx`, `Bisect` sizes and points, `Segment::P`, `MakeSmooth` (26,806) | 0 | 0 | 0 | 0 | 0 |
| `I_1` / `I_2` / `I_3` (2,400 each) | **1,625 / 774 / 1,694** | **1,625 / 774 / 1,694** | 0 | `I_2` **2 (1 ULP)** | `I_2` **2 (1 ULP)** |
| Constants (132) | 0 | 0 | 0 | 0 | 0 |
| DataSet on DS(CMF)-1 EOS: Evaluate p, ρ, Derivative, Integrate; DataColumn ops (11,583) | 0 | 0 | 0 | 0 | 0 |
| Whole output SHA-256 | `3fbf1285…` | `602a0cc9…` | **`376bc3fd…` = OLD** | `189e2718…` | `eabbb137…` |

Codegen facts (scratch builds of `c8c6813`): FMA-class instructions Debug 78 / Release 120 (vendored
0); canonical CONFIND 0 in both modes; Release `Math_Core.cpp.o` imports `___exp10` from
`GridVals_2D::Interpolate` (2 sites, inlined `Axis::operator[]`), `Axis::operator[]` (1) and the
`Quantity` printer (1); Debug calls libm `_pow` through libc++ `std::__math::pow<int>`. The residual
Release `I_2` difference is LLVM rewriting `pow(x,2)` to `x*x` on the `llvm.pow` intrinsic that
libc++ emits: 8 historical `pow` calls in `I_2` become 3, and `-fno-builtin` does not block it.

**ZM-1's own characterization suite (12 tests) passes in every variant above**, including FMA-on
and `exp10` builds: it does not discriminate either effect, just as CONFIND's C22 did not
discriminate M-1.

Scratch prototype of the §26 design (not a Zaki change): Debug and Release archives contain 0
FMA-class instructions and no `___exp10`; `I_2` keeps its 8 historical `pow` calls (through the
noinline boundary); probe output SHA-256 equals OLD (`376bc3fd…e1d0`) in **both** modes, i.e. 0 of
162,619 records differ; ZM-1 12/12 pass.

Downstream contour-input differences: none on the historical TaskManager grids (nodes and bilinear
values identical); the FMA drift enters TaskManager through `Curve2D::Intersection` →
`m_tot_range[1]` → `B_tot_grid.Evaluate` → B_tot contour levels, and through `Bisect`/`GetIdx` (§15).

---

## 12. TaskManager-native `pow(10, …)`

Sites on the TaskManager path: `TaskManager.cpp:162` (×2, per-thread dark-axis endpoints) and
`CompactStar/EOS/src/Model.cpp:36` (dark-EOS ρ grid in `FindEOS`, reached by `FindDarkEOS`).
Governed sites reached by tests (threading/pow audit): `Core/src/TOVSolver.cpp:2669`
(`SolveToProfile` coarse bracket), `Physics/Evolution/src/EvolutionConfig.cpp:74`, and test-side
`tests/eos/structure1/table.hpp:71`, `tests/rotochemical/trajectory.hpp:20`,
`tests/core/tov_reference_cmf.cpp:211-212`, `tests/core/tov_surface_contract.cpp:132-134`,
`tests/rotation/hartle_thorne_1968_hw_eos.hpp:164-165`,
`tests/eos/rotochemical_trackr_freegas_structure.cpp:153`.

`TaskManager.cpp` compiled with its OLD compile command at `-O0` imports `_pow` (via
`std::__math::pow<int>`); at `-O3 -DNDEBUG` `Task` imports `___exp10`. Both contain the same 4–5
`fmadd` (AppleClang default contraction). Verbatim arithmetic of `:144-170` compiled both ways:
0/770 range endpoints differ for the historical grids × T=1..10; **89/40,000** differ on random log
grids.

Consequences: (a) OLD-vs-NEW at the same build mode is unaffected (CompactStar code is compiled
identically on both sides); (b) CompactStar **Release** differs from the **Debug** governed authority
wherever these sites see a discriminating argument — a pre-existing CompactStar property, not a
migration effect. No CompactStar patch is required for migration equivalence; a CompactStar
historical-pow helper is required only if the owner later demands Release reproduction of Debug
authority (§26, owner decision OD6).

---

## 13. Coord3D link resolution

| Link | Provider of `Coord3D::operator<`, `operator==`, `XYDist2` |
|---|---|
| Every current `main/` and test executable (OLD) | **not linked** (no executable pulls `TaskManager.cpp.o`, so no CONFIND member is loaded) |
| OLD TaskManager-bearing scratch driver (full OLD link line, `-Wl,-map`) | vendored `libConfind.a(Cont2D.cpp.o)` — unfused, `-O0` |
| NEW Debug scratch driver (plotting-free TaskManager, canonical Zaki + CONFIND Debug) | canonical `libCONFIND.a(Cont2D.cpp.o)` — `-ffp-contract=off`, `-O0` |
| NEW Release scratch driver | **no out-of-line definition** (inlined; immune) |
| NEW Debug + a default-flags consumer TU using `std::set<Coord3D>` placed before the archives | **the consumer TU** for all three symbols; linked `XYDist2` contains `fmadd d0, d0, d1, d2` |

Canonical Zaki (Debug/Release) defines none of the three symbols; no CompactStar object at `812463a`
defines or references them (only `Segment(Coord3D, Coord3D)` is referenced by `TaskManager.cpp.o`).
Consumer flags therefore can alter the comparator only by adding an instantiating TU — including
future T1/T2 harnesses. The OLD executable also exports the three comparators as weak
definitions (`dyld_info -exports`: `[weak-def]`; the link map shows `.stub`/`.got` entries), so the
gate checks both the static link map and that no loaded image exports a competing weak definition.

Mandatory qualification (G4, §24): link map + `nm -m` + disassembly on every CONFIND-bearing
executable, plus a generalized weak-symbol audit of FP-bearing Zaki/CONFIND inline symbols. A
reconnaissance preview on the full T2 executables (OLD link map versus the NEW Debug link with the
Zaki 2.0.1 prototype) found 1,379 Zaki/CONFIND-namespace symbols in both images; 3 changed
provider class — two `Directory(Directory&&)` constructors (Zaki archive → CONFIND archive) and one
`std::vector<Coord2D>` buffer relocation (Zaki archive → CompactStar TU) — none with FP arithmetic.

---

## 14. TaskManager numerical flow and one-ULP sensitivity

```
inputs: visible .eos (path relative to wrk, TOVSolver.cpp:745-753) ; FindDarkEOS(m): Fermi_Gas,
        rho 1e-5..50 fm^-3, FindEOS(2000) [log grid, pow(10,.) Model.cpp:36] -> EOS/Fermi_Gas_<m>mn.eos (%.8e)
Work(T) -> Task(t): partition d_ax (TaskManager.cpp:144-170; pow(10,.) :162) -> TOVSolver::Solve_Mixed
        (TOVSolver.cpp:1858-2131): d loop p_of_e_dark(Axis_t[d]); v loop p_of_e(v_ax[v]); RK8PD
        -> MixedStar::SurfaceIsReached -> static sequence (TOVSolver_Thread.cpp:27-44)
        -> <wrk>/NStar/Dark_Core/<m>/<m>_<v>x<d>_Sequence.tsv (14 cols, %.8e = 9 digits)
ImportSequence (DataSet::Import; wrk must be absolute) -> sequence_grid
FindCriticalCurve (:237-484): GridVals_2D B_vis(col 6), M_tot(col 3+9) by index cols 0,1
        -> CONFIND 92 levels (0.70..0.955 step 0.005; 0.960..0.999 step 0.001) x max(col 6)
        -> [SetPlotConnected :342, Plot :347] -> GetContourSet (:352, raw cell-scan order)
        -> per level argmax of bilinear M_tot over raw points p = 0..n-2 (:366-398; strict >)
        -> MakeSmooth(10) (:409) -> Append visible-max point (col 2, col 8 at MaxIdx(col 3)) (:412-413)
        -> [Curve2D::Plot :415] -> Export Crit_curve_smooth.tsv (:417, %-15.10e = 11 digits)
FindMtotContour(M) (:487-546): CONFIND {M} -> ConvertToCurve2D (SortNew) -> mass_curve (:506)
        -> critical_curve.Import(TSV) (:512) -> [Plot :518-522] -> ExportContour (:523)
        -> Intersection(critical, mass) (:529) ; |I| != 1 -> error, stop ; else Intersection_<M>.tsv
        -> FindBtotContour(M, {mass_curve[0], I[0]}) (:537, :543)
FindBtotContour (:549-630): B_tot(col 6+12).Interpolate (:563); L_min = Evaluate(mass_curve[0]),
        L_max = Evaluate(I[0]) (:578-579); 10 levels (:584-587) -> ConvertToCurve2D ->
        Intersection (:603) ; |J| != 1 -> skip index ; Bisect(J[0]).first (:612) -> label %.9e ->
        B_conts/<M>/B_tot_<M>_<i>.tsv (:618) -> [Curve2D::Plot :623, Plot :627-629]
Precision_Task (:633-672): Contour.Import(B_tot_<M>_<i>.tsv) (exit on missing) -> TOVSolver
        (radial res default 10000) Solve_Mixed(Contour) with DarkCore_Analysis -> <M>_<i>_Sequence.tsv,
        BNV_rates_<M>_<i>.tsv ; FindLimits -> BNV_tau_*.tsv
```

Branches where a 1-ULP change can alter a discrete result: per-thread axis construction (NaN when
`(res+1)/T == 1`, underflow when `res+1 < T`); `p_of_e` clamps and Brent termination; RK8PD
trial-stage cutoff events (core→mantle, surface); 9-digit sequence quantization; CONREC vertex
status and case selection; CONREC log back-transform; `RMDuplicates` radius classes and order;
`SortNew` start point (`mass_curve[0]`) and nearest-neighbour ties; the critical argmax (strict `>`,
all-but-last point, `size()-1` underflow on an empty contour); `MakeSmooth`; `MaxIdx` first-max tie
in file row order; 11-digit critical-curve round trip; `Intersects` (parallel, closed-x-domain,
shared-vertex double count) and therefore the `|I| == 1` / `|J| == 1` gates and which intersection
is `[0]`; bilinear `Evaluate` levels; `GetIdx`/`Bisect` strict `<` cut index; `%.9e` label and
11-digit coordinate re-import; `Solve_Mixed(Contour)` trial thresholds. The one demonstrated drift
(§15) is the `Bisect` cut.

---

## 15. Real-stellar prototype (scratch reconnaissance, not governed)

Driver: scratch program calling the real TaskManager API — `TaskManager(T)`, absolute `SetWrkDir`,
`SetVisEOSDir("EOS/DS(CMF)-1_with_crust.eos")` (SHA-256 `5747dd73…47ae5dd`), `FindDarkEOS(0.8)`,
`SetChiMass(0.8)`, `SetGrid({{8e14,5e15},19,"Log"},{{5e13,5e16},19,"Log"})` (historical ranges of
`TaskManager.cpp:270-271` at 20×20), `Work()`, `ImportSequence`, `FindCriticalCurve()`,
`FindMtotContour(2.01)` (which runs `FindBtotContour`), then `Precision_Task(2.01)` and
`FindLimits(2.01)`.

NEW builds used the 74 non-plotting CompactStar TUs recompiled with their OLD compile commands
against canonical headers (0 errors) plus a copy of `TaskManager.cpp` with exactly its eight
plotting statements (`:342`, `:347-348`, `:415-416`, `:518-522`, `:623-625`, `:627-629`) commented
out; for the full chain, a copy of `DarkCore_Analysis.cpp` with its 18 active plot lines commented
was added. Scratch copies only; no repository source was edited.

| Run | Work wall time | Result versus OLD serial |
|---|---:|---|
| OLD serial, twice | 39.1 s / 39.0 s (400 stars, ≈0.1 s/star) | identical to each other: 15 text outputs (hashes below) |
| OLD T=2 | 22.3 s | all 400 rows present; rows identical after sorting by (d_idx, v_idx); file row order differs |
| OLD T=4 | 11.4 s | all rows present; **1 row differs** (d_idx 13, column 3 = M) |
| NEW Debug (unpatched c8c6813 + CONFIND b0cbd51) | 39.2 s | dark EOS, sequence, `Crit_curve_smooth.tsv`, `Intersection_2.01.tsv`, 9 of 10 B_tot files identical; **`B_tot_2.01_4.tsv` differs** (OLD line 29 `1.8238189670e+15 5.3978659336e+14` absent in NEW) |
| NEW Release deps (unpatched) | 33.7 s | same single difference |
| NEW Debug, Zaki `-ffp-contract=off` | 39.2 s | **all 15 text outputs byte-identical** |
| NEW Release deps, Zaki `-ffp-contract=off -fno-builtin` (contour stage on the OLD sequence) | — | **all contour-stage outputs byte-identical** |
| NEW Debug / Release deps, scratch Zaki 2.0.1 prototype (§26) | 39.2 s / 33.9 s | **all 15 text outputs byte-identical** |

Full chain (Work + contours + `Precision_Task(2.01)` + `FindLimits(2.01)`), 38 text outputs, with a
plotting-free scratch copy of `DarkCore_Analysis.cpp` (its 18 active plot lines commented) added to
the NEW library; these NEW executables link **without** `-undefined dynamic_lookup` and without
Python:

| Run | `Precision_Task` / `FindLimits` | Result versus OLD |
|---|---:|---|
| OLD, twice | 26.5 s / 1.3–1.4 s (the latter includes matplotlib PDFs) | identical to each other (0/38 differ); whole T2 ≈ 70 s serial |
| NEW Debug, unpatched `c8c6813` | 26.2 s / 0.002 s | **6/38 differ**: `B_tot_2.01_4.tsv`, `2.01_4_Sequence.tsv` (contour point 30 is a different star), `BNV_rates_2.01_4.tsv` (≤ 19 % relative), `BNV_tau_neutron.tsv` (≤ 5.2 %), `BNV_tau_lambda.tsv` (≤ 2.4 %), `BNV_tau_sigmam.tsv` (≤ 0.39 %) |
| NEW Debug, Zaki 2.0.1 prototype | 26.4 s / 0.002 s | **0/38 differ** |
| NEW Release deps, Zaki 2.0.1 prototype | 22.2 s / 0.002 s | **0/38 differ** |

OLD serial output SHA-256: `EOS/Fermi_Gas_0.8mn.eos` `0d3e2103…c93ddc`;
`0.8_19x19_Sequence.tsv` `c58088f2…31cddc`; `Crit_curve_smooth.tsv` `79ec12d1…5b880f`;
`Intersection_2.01.tsv` `19459377…10d7d55a`; `B_tot_2.01_4.tsv` `99383931…7869ec` (NEW unpatched
`5e524562…a19323`).

Additional observations:

- The M_tot contour export **content** is identical OLD/NEW, but its **name** is not: historical
  CONFIND truncates `wrk + "/" + "NStar/Dark_Core/0.8/M_tot_2.01" + "_2.01e+00.tsv"` at 149
  characters, producing `M_tot` for a 123-character root and `M_` for a 126-character root; CONFIND
  2.0 writes the full name. Full names need a root of at most 105 characters.
- `ImportSequence` → `DataSet::Import` resolves `"/" + (wrk + "/" + rel)`, so TaskManager requires an
  absolute working directory.
- The OLD embedded-Python `Curve2D::Plot` path ran and wrote PDFs; PDFs are not deterministic and are
  not observables.
- **Conclusion: historical TaskManager behaviour is capturable end to end and serial runs are
  deterministic; the drift of finding 1 reaches the final BNV limit outputs at the percent level;
  the contraction-off Zaki build (Debug) and the 2.0.1 prototype (Debug and Release) remove it.**

---

## 16. TaskManager fixtures and frozen observables

Both fixtures are added in stage M1 as **test-only** code on the OLD stack (no production change),
then carried unchanged through M2–M6. All observables are serialized as defined in §23.

### 16.1 T1 — deterministic algorithmic fixture (CTest, data-free, ≈1 s)

Input: the OLD serial 20×20 sequence file from T2 (400 rows, 14 columns), committed as a fixture
with its SHA-256, plus `v_ax={{8e14,5e15},19,"Log"}`, `d_ax={{5e13,5e16},19,"Log"}`, `M=2.01`,
`chi=0.8`; and a synthetic discriminating curve set (intersection/`GetIdx`/`Bisect` near-ties at
CompactStar scale, from the §11 probe generator) so that T1 discriminates arithmetic drift even
where the realistic data happens not to.

- **T1a (real API):** `ImportSequence` → `FindCriticalCurve` → `FindMtotContour(2.01)`. Frozen: bytes
  of `Crit_curve_smooth.tsv`, `Intersection_2.01.tsv`, `B_tot_2.01_{0..9}.tsv` (and the set of
  indices present), the M_tot export content, exit status, count of "There should be exactly one
  intersection!" messages.
- **T1b (white-box replica of `:237-630` using the same public calls):** records per level: level
  hex, found flag, raw contour size and full raw `(x,y,z)` sequence from `GetContourSet`; in OLD
  only, the raw sequence before and after `SetPlotConnected`+`Plot` (must be identical); per level
  argmax value, point and index `p*`; critical curve before `MakeSmooth`, after, and after the
  visible-maximum append; the re-imported critical curve; `mass_curve` full sequence and
  `mass_curve[0]`; `M_tot` intersection vector; `L_min`, `L_max`, the 10 levels; per level the
  converted curve, intersection vector, `GetIdx`, `Bisect` first/second sizes and first-half points.
  **Validity rule:** T1b's serialized exports must equal T1a's files byte-for-byte in the same build.
- Acceptance: every T1a byte and every T1b field bitwise equal OLD versus NEW in each declared build
  mode; OLD Plot-neutrality equality holds.

### 16.2 T2 — real stellar fixture (guarded, ≈70 s serial measured on the OLD stack)

Inputs: `DS(CMF)-1_with_crust.eos` (`5747dd73…`, copied into `<root>/EOS/`), dark Fermi gas
`m_χ = 0.8 m_n` generated by `FindDarkEOS(0.8)`, grid as T1, **`TaskManager(1)`** (serial), `M = 2.01`.
Calls: `FindDarkEOS(0.8)`, `SetChiMass(0.8)`, `SetGrid`, `Work()`, `ImportSequence`,
`FindCriticalCurve()`, `FindMtotContour(2.01)`, `Precision_Task(2.01)`, `FindLimits(2.01)`.

Preconditions (fail closed if violated): absolute root of at most 100 characters (asserted:
every CONFIND export path < 150 characters); fresh empty root; EOS hash verified; registered only
when `COMPACTSTAR_EOS_DATA_ROOT` and a short root variable (e.g. `COMPACTSTAR_T2_ROOT`) are set.
Frozen outputs: bytes of `EOS/Fermi_Gas_0.8mn.eos`, `0.8_19x19_Sequence.tsv`,
`Crit_curve_smooth.tsv`, `M_tot_2.01_2.01e+00.tsv`, `Intersection_2.01.tsv`,
`B_tot_2.01_{0..9}.tsv`, `2.01_{i}_Sequence.tsv`, `BNV_rates_2.01_{i}.tsv`,
`BNV_tau/BNV_tau_{neutron,lambda,sigmam}.tsv`; the set of files; exit status; message counts.
**T2b hex shadow:** a test-side replica of `Task` for T=1 (same `TOVSolver_Thread` settings,
`SetRadialRes(3.0e4)`, `SetMinIdxOffset(0)`) with an `Analysis` hook that records every
`MixedStar` sequence point in binary64; valid only if its exported sequence equals T2's file
byte-for-byte. Multi-thread runs (T=2,4) are recorded as **diagnostic only**, compared OLD-vs-NEW at
the same T after sorting rows, and are not an authority (race; §22).

Acceptance: all frozen T2 bytes and T2b binary64 fields identical OLD versus NEW in each declared
build mode. PDFs excluded.

Fixture fixing rule: M1 first runs T2 on the OLD stack only. If the 20×20 chain fails a
precondition (not exactly one M_tot intersection, a missing B_tot index, an empty-contour crash,
`Precision_Task` exit), the owner re-declares the grid (for example 30×30 on the same ranges)
before any NEW build exists; the grid, EOS, `m_χ`, `M` and thread count are never changed after a
NEW result has been produced.

### 16.3 Frozen observables by function (answers §14–§16 of the task)

- **FindCriticalCurve:** for all 92 levels: level hex, found flag, raw size, raw point sequence
  `(x,y,z)` hex in `GetContourSet` order, argmax value/point/index; smoothed curve hex; appended
  visible-maximum point; `Crit_curve_smooth.tsv` bytes; OLD raw-before-Plot == raw-after-Plot.
- **FindMtotContour:** raw contour, converted `Curve2D` sequence, `mass_curve[0]` hex, re-imported
  critical curve hex, intersection vector (size and hex), `m_tot_range`, `Intersection_2.01.tsv`
  bytes, M_tot export content.
- **FindBtotContour:** `L_min`, `L_max`, 10 level hex, per level raw contour, converted curve,
  intersection vector, `GetIdx`, `Bisect(...).first/second` sizes and points, label string, exported
  file bytes and index set; downstream `Precision_Task` inputs (`Contour.Import` curve and `val`),
  `Solve_Mixed(Contour)` sequences and `DarkCore_Analysis` outputs (T2).

---

## 17. Existing tests and baselines

76 registered tests at `812463a` (52 always + 24 under `COMPACTSTAR_EOS_DATA_ROOT`; 58 C++, 18
Python; unchanged since `232565a`). `adr0017_production_qualification` is built but not registered
(run through `tests/bnv/adr0017_production_verify.py`). No test exercises TaskManager, CONFIND,
MixedStar, `GridVals_2D`, `Curve2D`, `Grid2D`, `IntegrateTrapz`, `Newton`, or any plot API; every
test uses Linear axes only. Zaki paths exercised by governed tests: constants, `DataSet`/`DataColumn`
(linear GSL interpolation, `Integrate`), `Axis` Linear, `GSLFuncWrapper`, File/String I/O — all
bitwise unchanged in §11.

A fresh OLD Debug build of `812463a` with the data root registered 76 tests and passed 76/76 in
this reconnaissance (§31).

Governed baselines and comparator classes (from the test reconnaissance): B1
`baryon_number_dscmf1_reference.tsv` relative ≤ 1e-15; B5 `hartle_monopole_dscmf1_debug.tsv` raw
bytes; B6 `passive_cooling_cmf_1p6_debug.tsv` relative 1e-5 (states), 1e-4 (luminosities and — in
code — `dLnTinf_dt`, although the document says 1e-5); B7 `phase5b_structural_response.json` raw
bytes / 0-ULP except compiler string; B8 `phase5c_chemical_coefficients.json` exact recursive
equality except compiler string; B9 `phase5d1_controlled_evolution.json` empty difference inventory;
B2, B3, B4, B10, B11 are protected **only** by SHA pins (`phase5d_protected_manifest`) and are not
regenerated by CTest (their producers use `--emit`).

---

## 18. OLD and NEW builds

### 18.1 OLD (reference) build

CompactStar `812463a` (for M1: `812463a` + test-only T1/T2 commit), vendored arm64 archives and
headers of §1.2, reference toolchain of §1.3, fresh out-of-tree build, environment without
`CONDA_PREFIX`:

```
cmake -S <src> -B <build-old-{debug,release}> -DCMAKE_BUILD_TYPE={Debug,Release} \
  -DPython3_EXECUTABLE=/Users/keeper/miniforge3/bin/python3 \
  -DCOMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose
```

Recorded: `compile_commands.json`, `CMakeCache.txt` (GSL/Python/OpenMP resolutions), archive and
header hashes, `otool -L`, link maps. OLD-Release keeps the `-O0` vendored archives (that is the
historical composition).

### 18.2 NEW (treatment) build

Same CompactStar source plus only migration changes (M2 plotting deletion, M3 CMake/package
consumption). Dependencies are built outside every repository from clean clones at exact SHAs:

```
env -u CONDA_PREFIX cmake -S <zaki@SHA> -B <b> -DCMAKE_BUILD_TYPE={Debug,Release} \
  -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -DGSL_ROOT_DIR=/opt/local -DZAKI_BUILD_TESTS=ON -DZAKI_BUILD_EXAMPLES=OFF \
  -DCMAKE_INSTALL_PREFIX=<zaki-prefix>          # ctest must pass 12/12, then install
env -u CONDA_PREFIX cmake -S <confind@b0cbd51> -B <b> -DCMAKE_BUILD_TYPE={Debug,Release} \
  -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DZaki_DIR=<zaki-prefix>/lib/cmake/Zaki \
  -DGSL_ROOT_DIR=/opt/local -DCONFIND_BUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=<confind-prefix>
                                                # ctest must pass 29/29, then install
```

Prefixes: versioned, read-only after install, outside Git, keyed by dependency, source SHA, build
mode and toolchain (for example
`/Users/keeper/Documents/CompactStar/external/packages/zaki-<sha12>-<mode>-applclang2100.3.34.2/`),
each with a provenance manifest (source SHA and dirty flag, CMake cache excerpt, compiler and
flags, GSL/zlib resolutions, archive SHA-256, installed-header manifest SHA-256, member list, FMA
and `exp10` counts, test log hash). The Zaki SHA is the owner-accepted numerical-preservation
release (§26) — **not** `c8c6813` as built by default, which fails G3/G6 on present evidence.

Reference values from the scratch builds of `c8c6813`/`b0cbd51` (evidence only, not identities):
installed Zaki header manifest `8bd68eb12291731c6a86d04a1ce8ec08040fe2fe307ec9e45f3bc53b8370e910`
(34 files; `Version.hpp` reports 2.0.0, SHA `c8c6813…`, dirty 0); installed CONFIND header
manifest `0985d1168d7912e72b7dbf4f12c45ff6d7f61831ea381e4dd2251d89de21c5f7` (7 files); Zaki
12/12 and CONFIND 29/29 tests passed in Debug and Release.

### 18.3 Package authentication (fail closed)

- New cache variables `COMPACTSTAR_ZAKI_PREFIX` and `COMPACTSTAR_CONFIND_PREFIX` (absolute, required
  on every platform; no vendored default after M3);
  `find_package(Zaki <accepted version> CONFIG REQUIRED NO_DEFAULT_PATH PATHS ${COMPACTSTAR_ZAKI_PREFIX})`
  and the same for CONFIND (no `CMAKE_PREFIX_PATH`, environment, registry, `/usr/local`,
  `/opt/homebrew`, `/opt/local` or miniforge resolution of Zaki/CONFIND).
- Verify `Zaki_DIR`, `CONFIND_DIR` and every imported location lie inside the declared prefixes.
- Parse the installed `Zaki/Version.hpp`: `ZAKI_SOURCE_GIT_SHA` equals the accepted SHA,
  `ZAKI_SOURCE_DIRTY_KNOWN == 1`, `ZAKI_SOURCE_DIRTY == 0`, exact version string.
- CONFIND embeds no SHA: require the prefix provenance manifest and verify archive and header
  manifest hashes against it and against a predeclared expected-identity record (a generated
  `Version.hpp` for CONFIND is a possible later 2.0.x hygiene item, not required).
- Pin transitive discovery: `GSL_ROOT_DIR=/opt/local` with the version read from
  `gsl/gsl_version.h` (2.7.1), not from a `gsl-config` found on `PATH`; zlib resolved explicitly
  and recorded (SDK `libz.tbd` in the scratch builds).
- Print every resolved path and hash at configure time; qualification tooling cross-checks them.

---

## 19. ADR-0018 status and vendored-artifact disposition

Migration conflicts with ADR-0018 §3 clauses 1 and 3 (vendored archives as default Mac
authority), clause 2's four raw-path overrides (imported targets carry include directories,
compile features and dependency closure), clause 5 as worded (CONFIG-mode `find_package`), and
clause 11 (replacement is outside ADR-0018); it carries forward clauses 6–10 and the §5
equivalence rule, and generalizes clause 4's fail-closed explicit-path principle to every
platform. **Recommended smallest action: B — a short follow-up ADR-0019 ("Canonical external
Zaki/CONFIND package consumption") that supersedes ADR-0018 §3 clauses 1–3, 5 and 11 for
Zaki/CONFIND resolution, retains the rest, and names the accepted dependency identities.** ADR-0018
is not rewritten; on ADR-0019 acceptance only its status header gains the successor name, per
`docs/adr/README.md` lifecycle.

Vendored artifacts: **move** `dependencies/lib/**` and the vendored headers to an explicit
historical/reference tree (for example `historical/dependencies/`) in M7, excluded from every
active build and install rule (`CMakeLists.txt:188-191` currently installs them), with their hashes
recorded in ADR-0019. They remain the oracle for T1/T2/probe re-capture and must not be deleted.

---

## 20. Staging of plotting, Python, and link ownership

- Plotting (M2): a separate preparatory branch step on the **vendored** stack deletes the 75 active
  calls, 82 active setup lines and 3 `SetPlotConnected` calls, keeps every adjacent export, and
  deletes commented plot debris; for the 31 EXPORT_DATA sites a data export is added only where the
  owner wants the figure's data retained. Acceptance: all artifacts of M1 (76 tests' outputs,
  regenerated SHA-pinned baselines, T1, T2) byte-identical. Mechanical deletion is justified by §8;
  the test proves it.
- Python (M3): remove `Development` and `NumPy` components and `Python3::Python/Module/NumPy`
  linkage; keep `Interpreter` for tests. Required evidence: no `-undefined dynamic_lookup` in any
  link command; `otool -L` without libpython; `nm -u` of every executable resolves at link time;
  `_compress` bound to an explicitly linked zlib.
- Link ownership after migration: CompactStar direct = GSL (own TOV/integrator code), Threads
  (TaskManager/TOVSolver_Thread use `std::thread`), `Zaki::Zaki`, `CONFIND::CONFIND`; transitive =
  zlib (through `Zaki::Zaki`), GSL (also through Zaki), Threads (through CONFIND). OpenMP: CompactStar
  has no OpenMP use; the only OpenMP consumer was vendored `libConfind.a`. Static evidence: the
  `-Xclang -fopenmp` flag changes no FP-arithmetic sequence in any of 32,041 functions of the 83
  CompactStar TUs (differences are EH/control-flow layout, string-label numbering and destructor
  aliasing). Keep `OpenMP::OpenMP_CXX` through M3 to hold compile lines constant; remove it in a
  separate step (M3b) under the same exact gates.

---

## 21. Debug and Release authority

Recommended matrix (owner decision OD2): **Debug OLD vs Debug NEW** and **Release OLD vs Release
NEW**, both mandatory, each exact. OLD-Release is CompactStar Release + `-O0` vendored archives;
NEW-Release is CompactStar Release + Release dependencies (Zaki preservation release + CONFIND
2.0 Release). NEW-Release must therefore reproduce `-O0` dependency arithmetic, which requires the
Zaki patch of §26 (CONFIND already does, M-1). Release is compared only OLD-vs-NEW: governed
baselines are Debug by construction (B8 records `assertions-enabled`; B7 `"Debug"`), and CompactStar's
own Release-vs-Debug differences (§12) are outside this migration. Debug-only qualification is not
recommended: the owner's CONFIND M-1 decision rejected it for expensive stellar evaluators.
Release+LTO is not part of the matrix (no CompactStar LTO exists).

---

## 22. Threading and performance

Existing TaskManager threading (`TaskManager.cpp:39-56,106-127,130-183`): `std::thread` fan-out,
worker count `min(requested, hardware_concurrency())` (host-dependent), static contiguous
partition of `d_ax`, per-thread `TOVSolver_Thread`. Defects (pre-existing, not repaired here):
`TOVSolver_Thread.hpp:61-65` guards the static `mixed_seq_static` with a per-instance `std::mutex`,
and the destructor decrements `tov_counter` before a separate zero check
(`TOVSolver_Thread.cpp:27-44`): concurrent `Combine` is a data race and the last-finishing thread
can export before others combine; per-thread endpoint recomputation makes node values depend on T
(demonstrated); `(res+1)/T == 1` yields NaN nodes; `res+1 < T` underflows; `TaskManager.hpp:58`
uses libc++-internal `std::__1::chrono`. Measured OLD speed-up: 39.1 s (T=1) → 22.3 s (T=2) →
11.4 s (T=4).

CONFIND's threaded `Evaluate(EvaluatorFactory, ThreadingOptions)` could parallelize expensive
stellar sampling only if each worker owns an independent solver/EOS state; TaskManager computes all
stars before contouring, so nesting TaskManager threads with CONFIND workers would oversubscribe.
**Migration keeps TaskManager's existing threading unchanged, uses CONFIND with worker count 1
(the default), and qualifies serial TaskManager only.** A later performance phase (separately
predeclared) may: repair the static-sequence race; replace endpoint recomputation with global-index
node evaluation; then either keep TaskManager's partition or move per-star evaluation behind a
single deterministic executor (one level of parallelism, never both).

Performance sanity measured here: T2 Work 39.0 s OLD, 39.2 s NEW Debug, 33.7 s with Release
dependencies (CompactStar Debug).

---

## 23. Characterization artifact format

TSV, UTF-8, `\n`, fixed column order, explicit record ordering (as produced, never re-sorted except
where §16 declares sorting). Every floating field carries both `%a` hex-float and the 16-hex-digit
IEEE-754 bit pattern; integers in decimal; strings verbatim. Header block records: artifact schema
version; CompactStar SHA; dependency identities (Zaki SHA and dirty flag, CONFIND SHA, prefixes,
archive SHA-256, installed-header manifest SHA-256) or the vendored hashes; compiler string; build
type and effective compile/link flags; GSL version and path; EOS SHA-256; grid (ranges, res,
scales); level list; thread count; working-root length. A SHA-256 ledger covers every artifact;
comparison excludes only the declared provenance lines (dependency identity/path/hash, build
directory, timestamps).

---

## 24. Gates G0–G9 (pass/fail fixed now)

| Gate | Pass criterion (all must hold; any failure = stop, return to owner) |
|---|---|
| **G0** identity | CompactStar SHA as declared; OLD uses exactly the §1.2 hashes; NEW dependency prefixes contain the accepted Zaki/CONFIND SHAs (clean), predeclared archive and header-manifest hashes, Zaki 12/12 and CONFIND 29/29 tests passed in the same mode; configure resolved only the declared prefixes, `/opt/local` GSL 2.7.1 and the recorded zlib; compiler `clang-2100.3.34.2` |
| **G1** existing suite | OLD and NEW each pass the identical 76-test set (plus declared T1/T2 additions) in Debug with the data root; no test removed, renamed or relaxed; governed comparators unchanged |
| **G2** plotting neutrality | (a) static: vendored `Plot` stack-only, `SetPlotConnected` fields unread (disassembly hashes recorded); (b) OLD T1b raw-before-Plot == raw-after-Plot for all 92 levels; (c) the M2 plotting-deletion build on the vendored stack reproduces every M1 artifact byte-for-byte (76 tests' outputs, regenerated B2–B4/B10/B11, T1, T2) |
| **G3** Zaki differential | the §11 probe, formalized, plus all T1/T2 grids: **0 differing records** OLD vs NEW in every section, Debug and Release; every NEW Zaki and CONFIND archive member contains 0 FMA-class instructions and no `___exp10` import, and no libm call of a Zaki function's `-O0` inventory is rewritten into a different function or into inline arithmetic (elimination of an identical duplicate call excepted) |
| **G4** Coord3D/weak symbols | for every CONFIND-bearing executable: link map shows the three comparators provided by `libCONFIND.a(Cont2D.cpp.o)` (Debug) or absent (Release); no CompactStar/test object defines them (`nm -m`); linked `XYDist2` has no FMA-class instruction; weak-definition exports checked (no loaded image exports a competing definition); generalized audit over every Zaki/CONFIND-namespace symbol in the final image: provider class (CompactStar TU / Zaki archive / CONFIND archive) equals OLD's for every symbol containing FP arithmetic, class changes are permitted only for symbols with no FP arithmetic, and no chosen definition contains an FMA-class instruction where OLD's did not |
| **G5** T1 | all T1a bytes and T1b fields bitwise equal OLD vs NEW, Debug and Release; T1b validity rule holds |
| **G6** T2 | all T2 output bytes and the file set equal OLD vs NEW (serial), Debug and Release; T2b binary64 fields equal; preconditions held |
| **G7** governed rechecks | NEW Debug: the 6 content-compared baselines pass their own comparators; B2, B3, B4, B10, B11 (SHA-pin only) regenerated with `--emit` under OLD and NEW and byte-identical OLD vs NEW (M1 records OLD-vs-committed; any OLD-vs-committed mismatch is a pre-existing reconstruction finding and a stop before M3); Phase-5B/5C/5D regressions and the bounded ADR-0017 qualification (`adr0017_production_qualification` + verifier) pass; **plus** OLD-vs-NEW byte identity of every artifact the suite writes (tolerance comparators cannot mask drift); Release: OLD-vs-NEW byte identity of the same artifacts |
| **G8** package/link/symbol | no Python/NumPy/matplotlib in any compile or link line; no `-undefined dynamic_lookup`; no libpython/libomp (after M3b) in `otool -L`; all executables link with no unresolved symbol; zlib explicit; no vendored path in any command; install rules exclude historical artifacts |
| **G9** performance | T2 serial Work wall time NEW ≤ 1.25 × OLD per mode (median of 3, same machine, idle); CONFIND worker count 1; declared as a sanity bound, not a numerical tolerance |

---

## 25. No-post-hoc-tolerance rule

Every OLD-vs-NEW comparison is exact: byte identity for text artifacts, bitwise identity for
binary64 fields, identical sets and orders. Pre-existing tolerance comparators (B1 1e-15, B6
1e-5/1e-4, internal test bounds listed by the test reconnaissance) keep governing their own
baseline checks but do not qualify OLD-vs-NEW drift. The only admissible differences are the
declared provenance fields of §23, PDF files, and stdout/stderr text other than the declared
message counts and exit status. No tolerance, allowlist entry or reclassification may be added
after a result is seen; a failing candidate is rejected and returned to the owner with a new
predeclaration required.

---

## 26. Staging, and the Zaki and CompactStar patch decisions

**Central decision (task §36): B — Zaki requires a bounded numerical-preservation patch before
CompactStar migration.** Proposed ZakiLib **2.0.1** (patch release: no API change; restores
historical arithmetic on the shipped API; ZM-1 artifacts stay byte-identical, as observed):

1. PRIVATE `-ffp-contract=off` on every Zaki library TU for AppleClang/Clang/GNU (not exported to
   consumers), mirroring CONFIND.
2. A private `HistoricalMath` TU compiled `-fno-builtin -fno-lto`, noinline wrappers returning libm
   `pow(x, y)`, routing: `Axis::operator[]` Log branch (`Math_Core.cpp:1509`) and therefore
   `GridVals_2D::Interpolate`; `IntegralTable` `I_1/I_2/I_3` (`IntegralTable.cpp:14-80`) integer-power
   calls; plus every other `pow` site found by a complete Debug-vs-Release libm call-inventory
   audit of the library, including the formatting-only `Quantity` stream printer, so that the
   archive-level G3 check (no `___exp10`, no rewritten `pow`) is absolute rather than a list of
   exceptions.
3. Discriminating tests generated **from the vendored arm64 oracle** (not from source builds):
   the §11 sections, random log axes, intersection/`GetIdx`/`Bisect` near-tie fixtures, `I_x` with
   cancellation cases; required: starting Debug and Release FAIL, corrected Debug, Release and
   Release+LTO PASS; FMA and `exp10` absent from corrected archives.
4. Separate predeclaration, independent review, owner acceptance and tag decision, as for CONFIND M-1.

Feasibility (scratch prototype of items 1–2 on a copy of `c8c6813`; one new private TU, 13 routed
`pow` calls, one CMake block): bitwise equal to the vendored oracle on all 162,619 probe records
and on all 38 T2 text outputs, in Debug and Release (§11, §15). This removes the stop condition
"no clean preservation strategy can be predeclared"; it is evidence for Z1, not a Zaki change.

Alternatives for the owner: (A) consume `c8c6813` as-is — **rejected by evidence** (G3/G6 fail);
(C) avoid affected APIs in CompactStar — not feasible without copying Zaki code; (D) consume
`c8c6813` built with a consumer-imposed "historical arithmetic" recipe (`-O0` or
`-ffp-contract=off` plus `-fno-builtin`) — no source change, exact in scratch Debug, but fragile
(a default rebuild silently reintroduces drift), leaves Release `I_2` inexact, and puts numerical
authority in a build recipe.

**CompactStar-native `pow` (task §37): no patch required for migration** (§12). If the owner later
requires Release to reproduce Debug authority, a separate CompactStar `HistoricalPow10` helper for
`TaskManager.cpp:162`, `Model.cpp:36`, `TOVSolver.cpp:2669`, `EvolutionConfig.cpp:74` (and
test-side sites) plus a CompactStar contraction policy would be needed, with its own
predeclaration.

Stages derived from the evidence:

| Stage | Content | Exit gate |
|---|---|---|
| Z1 (ZakiLib repo) | Zaki 2.0.1 numerical-preservation patch as above | owner acceptance of 2.0.1 |
| M0 | owner review of this record; ADR-0019 drafted (PROPOSED) and accepted | ADR-0019 ACCEPTED |
| M1 | test-only T1/T2 harness on `812463a`; OLD Debug and Release reference capture (all suite artifacts, `--emit` regeneration of B2–B4/B10/B11, T1, T2, probe) | captures reproducible twice |
| M2 | plotting-only deletion on the vendored stack | G1, G2 |
| M3 | package consumption (ADR-0019 mechanism), vendored paths removed from the active build, Python linkage removed, OpenMP kept | G0, G1, G3–G8 |
| M3b | OpenMP linkage removal | G1, G5–G8 unchanged |
| M4 | (only if OD6 = yes) CompactStar Release-authority program — **not part of migration** | separate |
| M5 | Coord3D and weak-symbol final-link qualification | G4 |
| M6 | full Mac matrix, Debug and Release | G0–G9 |
| M7 | owner review; fast-forward integration; vendored artifacts moved to historical tree; ADR-0018 status header, invariant header and `CURRENT_ARCHITECTURE.md` updated | owner acceptance |
| M8 | Linux qualification (separate predeclaration) | — |
| M9 | EKU cluster qualification (separate predeclaration) | — |

---

## 27. Top five risks (by scientific consequence)

1. **FMA-contraction drift in canonical Zaki** (`GetIntersection`, `GetIdx`, `Intersection`,
   `Bisect`, `IntegralTable`): demonstrated, in Debug and Release, to change a `Bisect` cut and
   thereby the final BNV lifetime limits by up to 5.2 %.
2. **TaskManager's discrete, untested selection logic** (raw order, argmax, `|I| == 1` gates,
   `mass_curve[0]`, `Bisect` index, 9/11-digit round trips, empty-contour UB) turns ULP differences
   into different stars; no baseline existed before this record.
3. **Optimized libm rewriting** (`pow(10,x)`→`__exp10` in `Axis`/`GridVals_2D`; `pow(x,2)`→`x*x`):
   0.19 % of general log nodes; latent for grids other than those tested.
4. **Coord3D weak-symbol interposition:** inactive today, demonstrated mechanism; any future
   consumer TU instantiating the comparator can silently change CONFIND de-duplication in Debug.
5. **Link/package identity drift:** implicit GSL discovery (`gsl-config` from miniforge),
   `dynamic_lookup` masking unresolved symbols, `DataSet` layout change forbidding mixed header
   generations, non-equivalent x86_64 vendored archives.

Plot-call side effects rank below these: proven neutral on arm64.

---

## 28. Owner decisions required before implementation

- **OD1** Zaki path: 2.0.1 numerical-preservation patch (recommended) vs build-recipe preservation
  (D) — option A is excluded by evidence.
- **OD2** Qualify both Debug and Release OLD-vs-NEW (recommended) vs Debug only.
- **OD3** Governance form: ADR-0019 superseding ADR-0018 §3 clauses 1–3, 5, 11 (recommended).
- **OD4** Vendored artifacts moved to a historical tree and retained (recommended) vs kept in place
  vs deleted.
- **OD5** Accept the T2 fixture (real 20×20 DS(CMF)-1, `m_χ = 0.8`, `M = 2.01`, serial, ≈70 s
  measured, short-root precondition) as a guarded governed test.
- **OD6** Whether CompactStar Release must ever reproduce Debug authority (a separate future
  program; not required for migration).
- **OD7** Declare multi-threaded TaskManager outside governed authority until its race and
  thread-count dependence are repaired in a separate change (recommended).
- **OD8** For the 31 EXPORT_DATA plotting sites, whether their derived data must gain numerical
  exports or the figures are simply dropped.
- **OD9** Confirm arm64 behaviour as the TaskManager authority at `TaskManager.cpp:347` (x86_64
  `Plot` re-sorted contours).

---

## 29. Stop conditions (implementation)

Stop and return to the owner if: any identity in §1 or a declared dependency identity differs; the
OLD reference cannot be reproduced twice; historical TaskManager outputs cannot be captured; a
plotting deletion changes any artifact; any G3/G5/G6/G7 record differs; the Coord3D provider or a
weak FP symbol cannot be authenticated; a fixture is nondeterministic at T=1; any fix would change
frozen CONFIND or Zaki semantics; any tolerance, allowlist entry or reclassification would be
needed after a result; a Python/zlib/GSL resolution differs from the declared one; the Zaki
preservation release is not owner-accepted.

---

## 30. Central questions

- **Q1** Consume canonical Zaki 2.0 as-is? **No** (G3 and G6 fail on present evidence).
- **Q2** Does optimized `pow`→`exp10` alter results? It alters 0.19 % of general Release log-axis
  nodes (1 ULP) but none on the historical TaskManager grids tested; the drift that reaches
  TaskManager output is FMA contraction, present in Debug too.
- **Q3** Zaki 2.0.x patch first? **Yes** (§26).
- **Q4** CompactStar historical-pow patch? **Not for migration**; only for a future Release-authority
  program.
- **Q5** Remove all plotting without changing science? **Yes on arm64 authority** (static and
  dynamic evidence); prove per G2.
- **Q6** Python/NumPy out of the compiled build? **Yes**; the interpreter stays for tests.
- **Q7** Remaining APIs: Zaki families of §6 minus plotting; CONFIND's 12 symbols of §7.
- **Q8** Final Coord3D provider: OLD vendored `libConfind.a(Cont2D.cpp.o)`; NEW Debug canonical
  `libCONFIND.a(Cont2D.cpp.o)`; NEW Release none (inlined).
- **Q9** Proof: G4 link map, `nm -m`, disassembly, weak-definition export check and generalized
  weak-symbol audit, plus T1/T2 exactness.
- **Q10** T1: §16.1. **Q11** T2: §16.2. **Q12–Q14** Frozen outputs: §16.3.
- **Q15** Rerun: the complete 76-test suite with data root, `--emit` regeneration of B2–B4/B10/B11,
  the Phase-5D coupled/regression pipeline, and the bounded ADR-0017 qualification; none is
  "provably unaffected" until run, although every Zaki path they use was bitwise unchanged in §11.
- **Q16** OLD build: §18.1. **Q17** NEW build: §18.2–18.3.
- **Q18** Debug and Release both: **yes** (§21).
- **Q19** Keep vendored artifacts: **yes, moved to a historical tree** (§19).
- **Q20** Threading conflict: yes if nested; migration keeps TaskManager threads and CONFIND at one
  worker (§22).
- **Q21** Linux: §30.1. **Q22** EKU: §30.2.

### 30.1 Deferred Linux qualification (requirements)

Accepted dependency releases built on Linux from the same SHAs; GCC contraction policy honoured
(GNU dialect defaults to `-ffp-contract=fast`; Zaki sets `CMAKE_CXX_EXTENSIONS ON`); glibc libm
`pow`/`log10`/`exp10` differ from Apple libm, so Mac bitwise equality is **not** expected and Linux
needs its own T1/T2/governed baselines under a separate predeclaration; `TaskManager.hpp:58`
`std::__1::chrono`, `Tags.hpp:72` `using enum`, GNU-ld link order (target-based linking), GSL
version, Python stack, and the shared `/tmp/compactstar-tov-surf-ir/` path are prerequisites.

### 30.2 Deferred EKU cluster qualification

Only after Linux qualification and owner source-authority acceptance: the existing LB0–LB10 and
CQ0–CQ7 gates of `PHASE6_LINUX_CLUSTER_BOOTSTRAP_PREFLIGHT.md` §15 apply, with ADR-0019 package
prefixes replacing the four-path overrides, serial T2 on a compute node under Slurm, and no
cross-platform bitwise claim.

---

## 31. Reconnaissance evidence provenance and controls

All builds and runs used a session-scratch directory outside every repository: clones of Zaki
`c8c6813`, CONFIND `b0cbd51` and CompactStar `812463a`; scratch Zaki builds (Debug, Release, and
`-ffp-contract=off` / `-fno-builtin` variants) and CONFIND builds (Debug, Release); an OLD
CompactStar Debug build of `812463a` (83 library TUs, 18 s parallel build, 76 tests registered);
probe `zaki_diff.cpp` (SHA-256 `880cb715…30f268`), comparator `cmp_sections.py` (`e8d3f224…c10b5`);
link-map drivers; the real-stellar prototype (OLD, unpatched NEW, contraction-off NEW, and the
scratch Zaki 2.0.1 prototype, with scratch plotting-free copies of `TaskManager.cpp` and
`DarkCore_Analysis.cpp`). These are reconnaissance observations pending governed re-capture in
M1/M6; they are not qualification evidence and are not preserved in Git.

OLD full-suite run (scratch clone of `812463a`, Debug, vendored archives, data root set,
`ctest -j6`, AppleClang `clang-2100.3.34.2`): **76/76 passed** in 5,712.6 s wall, including every
content-compared governed regression (`baryon_number_cmf`, `hartle_monopole_regression`,
`passive_cooling_regression`, `phase5b_structural_response_regression`,
`phase5c_chemical_coefficient_regression` 358.5 s, `phase5d_coupled_oracles` 1,724.3 s,
`phase5d1_controlled_evolution_regression` 2,458.5 s) and the SHA-pin check
`phase5d_protected_manifest`. The test reconnaissance found no earlier record of a full-suite run
under `clang-2100.3.34.2`; the current governed baseline is reconstructible on this Mac with the OLD
stack.

Controls: canonical CompactStar, Zaki and CONFIND tracked trees unchanged and `master` at the §1.1
SHAs at exit; vendored archives unchanged; no CMake, source, test, baseline, EOS/data or literature
byte changed; no cluster access; no push.

---

## 32. Disposition and next action

**C — ZAKI / TASKMANAGER NUMERICAL DRIFT REQUIRES A SEPARATE PATCH BEFORE COMPACTSTAR MIGRATION.**

The migration architecture (canonical packages, fail-closed identity, plotting deletion,
Python removal, T1/T2 authority, exact OLD-vs-NEW gates) is sound, but canonical ZakiLib
`c8c6813` changes TaskManager's stellar outputs (final BNV limits by up to 5.2 %) through FMA
contraction in Debug and Release, and its Release build additionally rewrites historical libm
calls. A bounded preservation design is demonstrated feasible in scratch. **Exact next action:**
return this record
to the owner for OD1–OD9; on OD1 = patch, open a separately predeclared ZakiLib 2.0.1
numerical-preservation task (Z1) whose oracle is the vendored arm64 `libZaki.a`; CompactStar M0/M1
may proceed in parallel only after owner acceptance of this record.


---

## 33. Owner implementation authorization — 2026-09-29

The owner explicitly authorizes this implementation pass on the existing
`physics/external-deps-taskmanager-migration` branch, from predeclaration commit
`6ab1783b8003295224dcf2d720d500a174f9ef9e`. The original body above is preserved.
This section records the owner's current decisions; historical prerequisites
and recommendations above retain their original chronological meaning.

| Decision | Authorized implementation contract |
|---|---|
| OD1 | Consume exact clean Zaki 2.0.1 candidate `e263a6e180c5c417198e7778bd21fc9c0a32dc33`, provisionally authorized for this migration candidate. Do not substitute master, another SHA, or a system package. |
| OD2 | Require OLD Debug == NEW Debug and OLD Release == NEW Release, exactly. |
| OD3 | Write ADR-0019 as the owner-authorized successor for active dependency resolution; preserve ADR-0018 historically. Candidate Mac authority only until later acceptance/integration. |
| OD4 | Preserve historical vendored archives and headers permanently as immutable oracles, at their current paths; remove them from active include/link/install resolution. |
| OD5 | Accept DS(CMF)-1, dark mass 0.8 neutron masses, 20 by 20 logarithmic grid, M=2.01, TaskManager(1), absolute working root <=100 characters, through Precision_Task and FindLimits. |
| OD6 | Release need not equal Debug. Same-mode OLD-versus-NEW is the authority. |
| OD7 | Multi-threaded TaskManager is outside governed migration authority; preserve its implementation and use one thread. Record deferred work only. |
| OD8 | Delete visualization when numerical data already exists/exported. Preserve useful otherwise unavailable numerical data through deterministic export, without one export per plot or a new plotting layer. |
| OD9 | Historical arm64 semantics govern TaskManager plotting removal; historical arm64 CONFIND Plot is a proven no-op. |

One combined independent review occurs AFTER the complete implementation and
qualification, covering exact Zaki candidate e263a6e and the final CompactStar
candidate together. No intermediate or separate Zaki independent review is
required. Neither repository is merged or tagged here. Canonical CompactStar
remains `812463ac9ed374f64ac9cadd500066ab723d3a6c`; canonical CONFIND remains
`b0cbd510fd3fd0c772fa50499cd749287cb39e7b`. Dependency sources are read-only.
Local Mac only; no Linux or cluster authority or access.

No post-hoc tolerance, baseline rewrite, numerical-semantic repair, or
TaskManager threading redesign is allowed. Existing G0–G9 exact numerical,
provider, provenance and baseline gates remain controlling. OD1 authorizes the
actual bounded Zaki candidate, including its documented formatting-only
Quantity exp10 site; the historical proposed archive-wide no-exp10 design is
not a requirement to alter that candidate. Frozen numerical pow paths must
remain exact, and no other provider/arithmetic exception is introduced.

Entry authentication: canonical CompactStar and CONFIND local/origin/live
master identities match the SHAs above; all relevant trees are clean. The
migration branch contains only the original documentation predeclaration over
canonical master. Zaki's registered candidate worktree is
`/Users/keeper/Documents/CompactStar/worktrees/ZakiLib-2.0.1-fp-preservation`,
clean at the exact authorized SHA. Historical arm64 SHA-256 values match:
Zaki `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef`;
CONFIND `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2`.
AppleClang is 21.0.0 (`clang-2100.3.34.2`), macOS 26.6.2 build 25G83,
CMake 4.2.1, historical GSL `/opt/local` 2.7.1. Build-time network is forbidden;
read-only live-ref authentication and the final non-force migration push are
separate owner-authorized repository operations.
