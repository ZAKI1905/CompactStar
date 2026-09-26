# Phase-6 Linux / cluster bootstrap governance preflight

**Status:** PROPOSED — OWNER ACCEPTANCE REQUIRED

**Date:** 2026-09-26

**Canonical entry:** `232565a32303a4953f3f516d1d5286b6663f8f99`

**Imported discovery evidence:**
`docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md`, SHA-256
`c7e71653d44fb22606f6dafebdae24cb814a5629292573f0f2c15cfb7ff10e2e`

**Change class:** dependency/build + structural architecture + documentation

**Execution authority:** none. No build, test, install, cluster access, source transfer, or
qualification is authorized by this document.

## 1. Outcome and hard boundary

The EKU environment discovery is complete, but cluster qualification execution remains
**BLOCKED PENDING LINUX DEPENDENCY BOOTSTRAP**. CQ0-CQ7 remain valid in concept only downstream
of the LB0-LB9 gates in section 15. CQ0 is not authorized now.

Model B—separately built, authenticated static Zaki/CONFIND libraries supplied to CompactStar
by explicit paths—is the recommended architecture. The current authenticated Darwin archives
remain unchanged and remain the default Mac numerical authority. No implicit system search,
network fetch, or vendoring of new Linux binaries into CompactStar is proposed.

The audit found a material source-authority blocker: the current local Zaki candidate is
header-compatible with CompactStar, but the available `ZAKI1905/CONFIND` HEAD is materially
older than the vendored CONFIND interface. It cannot be treated as the source of the current
archive or used as a Linux dependency candidate. A matching CONFIND source snapshot must be
recovered and authenticated before LB1. The disposition is therefore **C**:

> ZAKI / CONFIND HAS MATERIAL SCIENTIFIC DIVERGENCE FROM CURRENT COMPACTSTAR DEPENDENCY
> CONTRACT — RETURN TO OWNER.

This is a dependency-contract divergence, not evidence that the vendored Darwin archive is
wrong. No scientific equation, Phase-5 baseline, or production source has changed.

## 2. Local-only evidence chain

This task did not access the EKU cluster. The local Mac cannot SSH directly to the cluster;
cluster access occurs from a separate work/campus environment. The sole cluster-side input was
the manually transferred ZIP:

| Item | Authenticated value |
|---|---|
| ZIP | `/Users/keeper/Downloads/PHASE6_EKU_CLUSTER_PREFLIGHT_TRANSFER.zip` |
| ZIP SHA-256 | `49718055eeacc2b81cb25a29a69a782985ab52436702760553ee336cc04760e6` |
| Exact members | `PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md`; `PHASE6_EKU_CLUSTER_PREFLIGHT_TRANSFER_MANIFEST.txt` |
| Extracted/imported preflight SHA-256 | `c7e71653d44fb22606f6dafebdae24cb814a5629292573f0f2c15cfb7ff10e2e` |
| Byte-identical import | YES |
| Import commit | `3338d05705af691321d6aa192f9f66d2945dcf38` |

The imported file is immutable discovery evidence. Its historical cluster values are not
rewritten or refreshed from Mac inference. UNKNOWN values remain UNKNOWN. In particular, this
document does not claim that a transferred historical observation is current cluster state.

## 3. Repository authentication

### 3.1 CompactStar

The canonical checkout was clean on `master`; HEAD, local `master`, `origin/master`, and live
`refs/heads/master` all equalled the canonical entry SHA. The documentation worktree is
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6-linux-cluster-bootstrap-preflight`
on branch `docs/phase6-linux-cluster-bootstrap-preflight`.

### 3.2 Zaki

| Field | Authenticated value |
|---|---|
| Path | `/Users/keeper/Documents/CompactStar/external/Zaki` |
| Status | clean tracked state; ignored build/output/documentation files present |
| Branch | `master` |
| HEAD / local origin / live origin | `b9ddebaded24962468954846f47238aec2726fd4` |
| Ahead / behind | `0 / 0` |
| Remote | `git@github.com:ZAKI1905/ZakiLib.git` |
| Origin HEAD | `master` |
| Tags | none |
| Latest commit | `Initial import of Zaki library`, 2026-08-31T22:12:39-04:00 |
| License | MIT (`License.txt`) |
| Build / target | CMake; static target `Zaki`; archive `libZaki.a` |
| Declared standard | C++17 |
| Candidate SHA | `b9ddebaded24962468954846f47238aec2726fd4` — candidate only |

Public headers are the `Zaki/{File,Math,Physics,String,Util,Vector}` trees. The current CMake
unconditionally requires GSL, zlib, Python development, and NumPy and links Python/NumPy into
the one static target. `DataSet` plotting uses `matplotlibcpp`. ROOT setup is commented out.
Examples are added unconditionally and contain historical `/opt/local/include/root6` paths.
Install rules install the target, generated config, and license but not the public header tree.

An ignored `build_keeper/libZaki.a` is not provenance: it predates the repository import and
has SHA-256 `4311773a127e8244ea0c29a9614dee9c073843f6cf17886e9751c43da730f20f`, which differs
from CompactStar's vendored archive. Current HEAD is selected only as a behavioral-equivalence
candidate, not claimed as the historical archive source.

### 3.3 CONFIND

| Field | Authenticated value |
|---|---|
| Path | `/Users/keeper/Documents/CompactStar/external/CONFIND` |
| Status | clean |
| Branch | `master` |
| HEAD / local origin / live origin | `89c5d9b731534e4289d9f686549d9f0ac178e567` |
| Ahead / behind | `0 / 0` |
| Remote | `https://github.com/ZAKI1905/CONFIND.git` |
| Origin HEAD | `master` |
| Tags | none |
| Latest commit | `Update README.md`, 2021-09-27T11:27:46-04:00 |
| License | GPLv3 (`LICENSE`) |
| Build / target | legacy Makefile; `genlib` produces `libconfind.a` |
| Declared standard | C++14 |
| Candidate SHA | **NONE / UNRESOLVED** |

The Makefile selects `g++` on Linux and `clang++-mp-9.0` elsewhere, uses `-fopenmp`, and
hard-codes ROOT under `/home/zaki/HEP_Tools/root_dir` or `/opt/local/lib/root6`. It expects a
lowercase `zaki` dependency path and library name. The public API uses `.h` headers, lowercase
`<zaki/...>` includes, and live ROOT types. No CMake target or install interface exists.

## 4. Authenticated Darwin artifacts and provenance

CompactStar selects both archives from
`dependencies/lib/<dependency>/${CMAKE_SYSTEM_NAME}/${CMAKE_HOST_SYSTEM_PROCESSOR}` and fails
configure when either file is absent (`CMakeLists.txt:89-101`). The library is linked to both
archives (`CMakeLists.txt:170-180`).

| Dependency | SHA-256 | Size | Archive members | Source provenance |
|---|---|---:|---|---|
| Zaki | `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef` | 2,478,624 B | `Instrumentor`, `MemoryManager`, `Profile_Timer`, `Logger`, `ObjObserver`, `Simple_Timer`, `Banner`, `Directory`, `String_Basic`, `TextBox`, `CSVIterator`, `CSVRow`, `DataColumn`, `DataSet`, `IntegrateTrapz`, `TempDC`, `Func2D`, `Math_Core`, `NDimContLevel`, `IntegralTable`, `Newton`, `Constants`, `Coordinate`, `DateTime`, `Sun` objects plus symbol table | **PARTIAL** — consumed headers match candidate HEAD, but no immutable record maps the archive to a source SHA/tag and the retained ignored build differs |
| CONFIND | `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2` | 703,208 B | `Base`, `Cell`, `Cont2D`, `Common`, `ContourFinder` objects plus symbol table | **UNKNOWN** — no SHA/tag mapping and available repository HEAD is materially incompatible |

The Zaki archive has unresolved Python C API, GSL spline/interpolation, and zlib symbols. The
CONFIND archive exposes the evolved contour API described below and has OpenMP and zlib
references but no ROOT references. These inventories inform the future link closure; they do
not establish source identity.

## 5. Consumed headers and ABI consistency

### 5.1 Zaki

Every directly consumed Zaki header is byte-identical to candidate
`b9ddebaded24962468954846f47238aec2726fd4`:

| Header (under `Zaki/`) | SHA-256 | Result |
|---|---|---|
| `File/CSVIterator.hpp` | `924b2bac9789c6a5978b2b20cf6f9726fae13f10f78364ffb0709ff6c4dbf819` | BYTE_IDENTICAL |
| `File/VecSaver.hpp` | `8f6ef4e30716917ce043746dd3b9968ca6d76adc0674d7b7907d76ea193f3461` | BYTE_IDENTICAL |
| `Math/GSLFuncWrapper.hpp` | `3fde46023598cec451a015b4532a1111743786753644a969d05e2989972b3160` | BYTE_IDENTICAL |
| `Math/GSLMultiFWrapper.hpp` | `4de7609b1d235d70dc299a4cd724c2a4a106d3b43b3585df2f61f67b5d6a9b44` | BYTE_IDENTICAL |
| `Math/GSLMultiFdfWrapper.hpp` | `90e18fa5221fb287d2db98aab1b120d20ae765567768e3ec491e302c96fb395f` | BYTE_IDENTICAL |
| `Math/IntegralTable.hpp` | `2a7a191cde90f06b751132dad36bab10c36956fac1bc20f9fa2f83857257d1b5` | BYTE_IDENTICAL |
| `Math/Math_Core.hpp` | `d3bdd7da10abd15fc66449964421e9b08288fe66f14cce33231efcec3367362c` | BYTE_IDENTICAL |
| `Math/Newton.hpp` | `aa5b36736ea022ed29be6d1511def0305859c1d0d686fa3cbaef54c78b68e34e` | BYTE_IDENTICAL |
| `Physics/Constants.hpp` | `c7a3e1a7a32fd91808f6a1ed130cad0c3caae4643bb00135529060cfeeba719b` | BYTE_IDENTICAL |
| `String/Banner.hpp` | `46dce1269a80ce4938901e97f278f309a4276f5d972d5493fea9246ebf6461b9` | BYTE_IDENTICAL |
| `String/Directory.hpp` | `5a7a561cf8b3db62aa827ff1f550c3d5758d0b98f349312cbfa79c377f251fce` | BYTE_IDENTICAL |
| `Util/Instrumentor.hpp` | `2c73af797fc5e0298b8d117c866441c986d583aae06e18df56b2754ce16987bc` | BYTE_IDENTICAL |
| `Util/Logger.hpp` | `ef52fc252bd9853860bdeb27b39ac93fffb2c75c093ae5b9475602df095bb5d5` | BYTE_IDENTICAL |
| `Vector/DataColumn.hpp` | `d99f79c6c1e8c66559491abaf9b71126ca27e4b3f05d86145ae12966289a6624` | BYTE_IDENTICAL |
| `Vector/DataSet.hpp` | `66283b2703b3621badbc434725f93870d798613af2ca4c8ff3668a56b07f9647` | BYTE_IDENTICAL |
| `Vector/IntegrateTrapz.hpp` | `a9668875517bae8037657db427ae086e9793ac75b6e29494dceb777c9a5a9018` | BYTE_IDENTICAL |
| `Vector/Vector_Basic.hpp` | `df90bbc5583a4dfe8f6e931aa1b040cf4da146c1046938eac9f09d6d8296c134` | BYTE_IDENTICAL |

This supports testing the Zaki HEAD; it does not prove that HEAD built the archive.

### 5.2 CONFIND — material divergence

CompactStar consumes `dependencies/include/Confind/ContourFinder.hpp`, SHA-256
`6b95b0a30c5e9f03aee022ed4e1994a3aed34b470f31472c90a550e6f20a09dd`. The logical
external counterpart is `include/ContourFinder.h`, SHA-256
`dd1c0ada1f04239bd089de63e0c2fee58f9fb6302621c79ae74ebe97d5a9f069`.
The result is **MATERIAL_DIVERGENCE**, not a spelling-only difference. Supporting comparisons:

| Vendored header / SHA-256 | External header / SHA-256 | Result |
|---|---|---|
| `Base.hpp` / `3e089b6d9429edd9dc859956fe1ce63e632dde1dfef1cacac273df4b3f9d2e7b` | `Base.h` / `5daf213205a8ab17073173a3be21a2573910d0ab61e2374aec89f6cb478885bf` | MATERIAL_DIVERGENCE |
| `Cell.hpp` / `b4b7c155fdf3439dfb4a835954cb3a9a7ca40a1628954776d4c26c73611ffea4` | `Cell.h` / `48328daca803cbfa8f409377a17f0fd92035a01719b40ad64e5cc78a9989a38d` | MATERIAL_DIVERGENCE |
| `Common.hpp` / `aaeb12aed6f620b10034b10d79faccd04313bda48cf77ba8247c2d47edf22bb8` | `Common.h` / `408f5acfd4a7e5d4fd012ab2f37cf33640fbbe1e15a666084783fd878864302d` | MATERIAL_DIVERGENCE |

The vendored interface adds the `Zaki::Math::GridVals_2D*` input, `GetContourSet()`, and
`Cont2D::ConvertToCurve2D`; the public HEAD lacks those CompactStar-required APIs and exposes
ROOT graph/legend types. Header layout, include spelling, data model, and exported ABI differ.
This is a material dependency-version risk and blocks nomination of CONFIND HEAD.

## 6. Actual CompactStar dependency use

### 6.1 Zaki API inventory

The complete current API families and callers are:

| Header/API family | Symbols | CompactStar callers | Category / purpose |
|---|---|---|---|
| `Physics/Constants.hpp` | particle and quark masses; `MEV_2_INV_FM`; `MEV_FM3_2_{Dyn_CM2,G_CM3,INV_KM2}`; `INV_FM4_2_{Dyn_CM2,G_CM3,INV_KM2}`; `K_BOLTZ_EV`; `LIGHT_C_{M_S,KM_S}`; `Q_E`; `SUN_M_KM`; `YR_2_SEC`; `GEV_2_S`; `Element` | `AngularVelocity.hpp`, `Units.hpp`; `Core/src/{MixedStar,RotationSolver,TOVSolver}.cpp`; `EOS/**`; `Extensions/**`; `Microphysics/BNV/**`; thermal/evolution sources | GENERAL_PHYSICS_CONSTANT_OR_CONVERSION, with compact-star ownership risk for solar-mass/geometric and EOS-density conversions |
| `Math/GSL*`, `Newton`, `IntegralTable`, `Math_Core` | `GSLFuncWrapper`, `GSLMultiFWrapper`, `GSLMultiFdfWrapper`, `Newton`, `I_2`, `I_3` | solver/EOS/BNV sources including `MixedStar.cpp`, `TOVSolver.cpp`, `SigmaOmega*.cpp`, `Particle.cpp`, `Baryon.cpp`, BNV channels | GENERIC_NUMERICAL_UTILITY |
| `Math` geometry/data types | `Axis`, `Range`, `Coord2D`, `Segment`, `Curve2D`, `Grid2D`, `GridVals_2D`, `Quantity`, `Cond_Polygon` | `Core/{TOVSolver,TaskManager,RotationSolver,Pulsar}*`, `LightDM*`, `BNV*`, `SpinState.hpp` | GENERIC_NUMERICAL_UTILITY; plotting methods are OTHER |
| `Vector` | `DataColumn`, `DataSet`, `DataSet::PlotParam`, `Exists`, trapezoid/vector helpers | broad `Core/**`, `EOS/**`, `Extensions/**`, `Microphysics/**`, `Physics/Evolution/**`, `State.hpp`, `Geometry.hpp` | GENERIC_NUMERICAL_UTILITY plus data output/plotting |
| `File` | `CSVIterator`, `VecSaver`, `FileMode::Write` | `Core/src/{TOVSolver,RotationSolver,TaskManager}.cpp`, `EOS/src/{Common,Model}.cpp`, BNV/extension analyses | GENERIC_NUMERICAL_UTILITY / data output |
| `String` | `Directory`, `Banner`, `EndsWith`, `Multiply`, `Pars`, `Strip` | broad core/EOS/analysis/observer/run-path code | GENERIC_NUMERICAL_UTILITY |
| `Util` | logger levels/manager, instrumentation macros | `TaskManager.cpp` and instrumented production sources | GENERIC_NUMERICAL_UTILITY / OTHER |

The boundary is sufficient to qualify an external library, but it is broad. Zaki remains the
reusable general-purpose library. Compact-star/nuclear-astrophysics-specific conversions
should ultimately be CompactStar-owned; that later cleanup is not a prerequisite to this
source-equivalence gate and is not authorized here.

### 6.2 CONFIND API inventory

`CompactStar/Core/src/TaskManager.cpp` is the sole caller. It constructs `ContourFinder`, then
uses `SetGrid`, `SetWrkDir`, `SetContVal`, `SetGridVals`, `SetPlotConnected`, `Plot`, and
`GetContourSet`; it consumes `Cont2D::size`, indexed values, `GetVal`, `ConvertToCurve2D`, and
`ExportContour`. The path performs numerical contour generation, data export, and plot output.

The authenticated current archive has no unresolved ROOT symbol and the vendored headers do
not require ROOT. Therefore ROOT is **not required by CompactStar's current consumed CONFIND
ABI**. Plotting is still present through Zaki's plotting layer. The minimum Linux dependency
must reproduce the entire current caller-visible interface above; the older public repository
cannot be trimmed into that target without first recovering the intervening source authority.

## 7. Portability findings

| Dependency | Finding | Classification |
|---|---|---|
| Zaki | AppleClang is forced only under `APPLE`; Linux may use GCC | MECHANICAL_BUILD_FIX (intrusive Mac behavior, not Linux blocker) |
| Zaki | Python development + NumPy + embedded matplotlib linkage is unconditional for the single library target | MECHANICAL_BUILD_FIX; ABI/link-closure qualification required |
| Zaki | GSL and zlib are unconditional legitimate dependencies | NONE once authenticated paths exist |
| Zaki | examples are unconditional and contain MacPorts ROOT paths | MECHANICAL_BUILD_FIX |
| Zaki | install rules omit public headers | MECHANICAL_BUILD_FIX |
| Zaki | ROOT and OpenMP logic is commented, not an active library dependency | NONE for current target |
| CONFIND | Make-only build, lowercase naming/path assumptions, no install interface | MECHANICAL_BUILD_FIX |
| CONFIND | absolute Linux/macOS ROOT paths and live ROOT public types | API/ABI_CHANGE if removed; MECHANICAL_BUILD_FIX only for path discovery |
| CONFIND | hard-coded compiler split and `clang++-mp-9.0` | MECHANICAL_BUILD_FIX |
| CONFIND | unconditional OpenMP and old lowercase Zaki dependency | MECHANICAL_BUILD_FIX plus ABI/link qualification |
| CONFIND | public HEAD lacks CompactStar-required methods/types | **API/ABI_CHANGE / MATERIAL dependency risk** |
| CONFIND | scientific equivalence of an unrecovered newer source | **UNKNOWN** |

Linux support for the authenticated Zaki candidate appears mechanically attainable without a
CompactStar-visible semantic change. Linux support cannot yet be assessed for the actual
CONFIND contract because its matching source is missing.

## 8. Language-standard audit

CompactStar declares strict C++17 with extensions disabled (`CMakeLists.txt:49-51`,
`:109-116`), but `CompactStar/Physics/State/Tags.hpp:72` contains `using enum StateTag;`, a
C++20 language feature. Existing Mac documentation records AppleClang accepting it as an
extension (`docs/build/MACOS_BUILD.md:455-477`). No other active `using enum` was found in
CompactStar, Zaki, or CONFIND.

Verdict: **COMPILER_EXTENSION_CURRENTLY_USED**. The present source is not formally C++17
portable. GCC commonly diagnoses later-standard constructs under pedantic modes while
accepting some as extensions; no compile was authorized here. A future CQ1 must record GCC
14.1.0's actual configure/build result. Do not silently raise the project standard: owner must
choose either a narrow C++17 spelling repair or a separately governed C++20 change.

## 9. Recommended dependency resolution model

**Recommendation: MODEL B — authenticated external static archives.**

- Model A (`FetchContent`, submodule, `add_subdirectory`) would couple three build systems,
  invite build-time network access, and make CompactStar own external source orchestration.
- Model B preserves dependency boundaries and Darwin authority, makes source/toolchain hashes
  explicit, supports build-once/read-only consumption, and isolates license/build closure.
- Model C would add platform binaries to CompactStar Git and blur source/toolchain provenance.
- No Model D is superior for the current contract.

The narrow future CompactStar interface should use cache entries consistent with the existing
`COMPACTSTAR_*` cache-variable namespace:

| Variable | Type | Semantics |
|---|---|---|
| `COMPACTSTAR_ZAKI_LIBRARY` | `FILEPATH` | exact `libZaki.a` |
| `COMPACTSTAR_ZAKI_INCLUDE_DIR` | `PATH` | root containing `Zaki/...` |
| `COMPACTSTAR_CONFIND_LIBRARY` | `FILEPATH` | exact `libConfind.a` |
| `COMPACTSTAR_CONFIND_INCLUDE_DIR` | `PATH` | root containing `Confind/...` |

On Darwin, an unset variable defaults to the exact current vendored path, including the shared
`dependencies/include` root. On non-Darwin, all four variables are mandatory, absolute, and
validated for file/directory existence plus sentinel headers. Configuration must fail closed.
No `find_library`, `find_path`, filesystem probing, environment fallback, `FetchContent`, or
network is allowed. Configure output records resolved paths; qualification tooling records
library/header hashes and rejects path/hash drift. This is design only; CMake is unchanged.

## 10. Candidate source authority

| Dependency | Candidate for Mac equivalence | Historical identity claim | Gate |
|---|---|---|---|
| Zaki | `b9ddebaded24962468954846f47238aec2726fd4` | NO; consumed-header match only | owner may approve it for LB1 testing after ADR decision |
| CONFIND | **NONE / UNRESOLVED** | NO | recover source matching vendored `.hpp` API and archive ABI; authenticate repository/commit/license |

No tag, version macro, archive string, build record, or local provenance establishes an exact
vendored source SHA. Archive timestamps and commit dates are circumstantial only and are not
source identity. The CONFIND mismatch prevents an acceptable joint source-authority proposal.

## 11. Predeclared Mac source-build equivalence gate

This experiment is required before either source SHA may become Linux cluster authority. It
is not run by this task.

1. Owner ratifies or revises ADR-0018 and explicitly approves one exact candidate SHA per
   dependency for testing. CONFIND first requires recovery of a matching source candidate.
2. Create clean disposable macOS arm64 checkouts/worktrees at those SHAs. Record source tree
   manifests, licenses, submodules, and clean status.
3. Build separate static archives in scratch with a declared AppleClang/CMake/build-mode,
   SDK/deployment target, C++ standard, GSL, zlib, Python/NumPy, OpenMP, flags, and linker/ar
   configuration. Never replace or rewrite a vendored archive.
4. Record source/header/object-member/exported-symbol/archive hashes and the full link closure.
5. After separately authorized CMake override implementation, configure two disposable builds
   of canonical CompactStar: control with vendored Darwin inputs and treatment with explicit
   source-built archive/include paths. All non-dependency inputs and toolchain settings match.
6. Run focused dependency contracts: every consumed header compiles; all referenced symbols
   link; Zaki constants/conversions, numerical wrappers, CSV/data/string utilities, and
   CONFIND contour/grid/export behavior match the authenticated contract. ROOT absence on the
   consumed CONFIND path is an explicit detector.
7. Run all relevant Phase-5B, Phase-5C, and Phase-5D governed regressions and the bounded
   ADR-0017 Phase-6 qualification against control and treatment.
8. Produce a signed-off comparison manifest. Any mismatch outside the predeclared fields
   below rejects the candidate; no post-result tolerance selection is permitted.

### 11.1 Acceptance classes

**BYTE_IDENTICAL_REQUIRED**

- every consumed installed header against the approved source tree;
- all authenticated scientific input and governed baseline bytes;
- deterministic textual/JSON/TSV scientific payloads after excluding only explicitly declared
  provenance fields;
- Phase-5D and ADR-0017 accepted-step/checkpoint schedules, counts, validity flags, and
  deterministic serialized binary64 fields under the same compiler/toolchain.

Archive bytes are not required identical because `ar` headers, object paths, and build metadata
can differ. Their ordered member set and normalized exported ABI must be exactly equal unless
the owner pre-accepts a documented non-consumed difference.

**EXACT_SEMANTIC_IDENTITY**

- dependency-focused API results and error/fail-closed behavior;
- test inventory and pass/fail outcome;
- link closure for the consumed path, including no ROOT requirement for current CONFIND use;
- provenance/currentness decisions.

**NUMERICAL_TOLERANCE**

No new cross-build tolerance is authorized: treatment minus control is **0 ULP / exact
binary64 identity** for deterministic same-Mac scientific fields. Each governed regression
also retains its already-declared internal analytic/reference tolerance; those tolerances do
not permit control/treatment drift. If exactness fails because the source build demonstrably
changes compiler/archive mechanics, the candidate is rejected and the owner must predeclare a
new experiment—results from this gate cannot be rescued by choosing a tolerance afterward.

Allowed byte differences are limited to paths, timestamps, build IDs, archive container
metadata, and the explicit dependency source/path/hash provenance fields. They are not
scientific numerical tolerance fields.

## 12. Future cluster toolchains

### 12.1 Compiler and external libraries

Use the transferred discovery's Spack GCC 14.1.0 and CMake 3.29.4 unless CQ1 exposes an
incompatibility. After Mac equivalence and owner source acceptance, build static Zaki and
CONFIND archives from the same accepted SHAs under a key such as:

`/mnt/sdd/zaki/CompactStar/toolchains/deps/<dependency-key>/{include,lib,provenance}`.

The key/manifests capture source SHA and tree manifest, compiler, standard, flags, CMake/make,
GSL/zlib/Python/OpenMP link closure, headers, member/symbol inventories, archive hashes, and
licenses. Builds occur under Slurm on a compute node if nontrivial, never as heavy login-node
work.

### 12.2 GSL

The transferred evidence identifies Spack GSL 2.6, while the Mac numerical authority uses
2.7.1. Bootstrap a user-local GSL 2.7.1 under versioned toolchain storage, statically linked
where practical. Record upstream source version and checksum/signature, configure flags,
compiler, static/shared policy, installed header hashes, archive hashes, and license/source
obligations. Matching 2.7.1 removes an unnecessary RKF45/rk8pd qualification variable.

### 12.3 Python

Local source audit finds exactly three non-stdlib runtime packages for CQ producers,
comparators, and validators: NumPy, SciPy, and mpmath. NumPy is broad across analysis and
rotochemical scripts; SciPy is used by `tests/relativity/revalidate_background.py` and
`tests/analysis/chemical_trackr_budget.py`; mpmath is used by chemical reference/budget/
evidence scripts and rotochemical oracles. `tests/bnv/adr0017_production_verify.py` is stdlib.

For the transferred Python 3.11.9 authority, propose a versioned venv with exactly
`numpy==2.3.5`, `scipy==1.17.1`, and `mpmath==1.4.1`. These releases declare Python 3.11
support and SciPy 1.17.1 accepts NumPy 1.26.4 through below 2.7. LB5 must freeze wheel/source
filenames and SHA-256 hashes in a hash-locked requirements manifest; no Conda is assumed.
Developer/doc extras are excluded.

## 13. Manual transfer protocol

Every Mac/work-machine/cluster handoff uses one `<descriptive-name>.zip` containing
`TRANSFER_MANIFEST.txt`. The producer manifest states purpose, source environment, target
environment, source Git SHA, authorized file list, per-file SHA-256, scope, and credential
exclusion. The producer reports the ZIP SHA-256 externally because a ZIP cannot reliably
contain its own final hash.

The receiver manually copies the ZIP, verifies the external ZIP hash, extracts to scratch,
verifies every member hash, and imports only authorized files. Filenames never establish
source state. SSH keys, Git credentials, tokens, cookies, and passwords are forbidden.

Normal cluster-side Git access is preferred for CompactStar and public/authenticated sources.
If private ZakiLib cannot be fetched with an already configured user-owned cluster credential,
the fallback is a Mac-created, manually transferred ZIP containing an exact `git bundle` plus
manifest at the accepted commit. An exact source tree with manifest is second choice. No
credential is bundled and no repository is made public.

## 14. Proposed persistent layout

Carry forward the historically authenticated root `/mnt/sdd/zaki/CompactStar`; do not create
anything from this Mac task.

| Directory | Purpose |
|---|---|
| `repo/` | primary authenticated CompactStar clone |
| `worktrees/` | governed Git worktrees |
| `toolchains/` | versioned GSL, Python venvs, and external dependency builds |
| `builds/` | disposable/keyed CMake build trees |
| `releases/` | immutable qualified executable/library bundles |
| `inputs/` | authenticated immutable copied scientific inputs when required |
| `qualification/` | compact durable CQ evidence/manifests |
| `campaigns/` | durable campaign manifests and results |
| `scratch/` | large transient per-job data |

The login node is limited to Git/network acquisition, small file management, job submission,
and evidence collation. Compilation, tests, numerical runs, and nontrivial dependency builds
run in Slurm allocations on compute nodes.

## 15. Ordered bootstrap gates

No step may skip its predecessor.

| Gate | Required result |
|---|---|
| LB0 | Owner ratifies or revises ADR-0018. |
| LB1 | Both exact candidate source SHAs pass the Mac source-build equivalence gate in section 11. |
| LB2 | Owner explicitly accepts the two source SHAs as cluster dependency authority. |
| TRANSFER GATE A | Mac packages the accepted dependency-authority manifest and bootstrap metadata in a governed ZIP; a human transfers and the receiving environment verifies it. |
| LB3 | Authorized cluster task creates the internal directory structure. |
| LB4 | Cluster task builds/authenticates user-local GSL 2.7.1. |
| LB5 | Cluster task creates the pinned, hash-locked Python environment. |
| LB6 | Cluster task obtains exact canonical CompactStar source. |
| LB7 | Cluster task obtains exact accepted Zaki and CONFIND source revisions. |
| LB8 | Compute-node jobs build and authenticate Linux static archives. |
| LB9 | Compute-node job configures/builds CompactStar with explicit dependency overrides and captures all resolved hashes. |
| LB10 | Only after LB0-LB9 evidence is accepted may CQ0-CQ7 begin. |

The exact CQ0 gate is: accepted ADR-0018; Mac-equivalent and owner-accepted Zaki and CONFIND
source SHAs; verified manual-transfer authority; authenticated GSL/Python/dependency artifacts;
and a successful fail-closed Linux configure/build at the canonical CompactStar SHA with all
dependency paths and hashes recorded. This preflight reaches none of those execution gates.

## 16. Test inventory: 76 versus 77

Static registration audit at the canonical SHA finds 52 always registered tests plus 24 under
the valid `COMPACTSTAR_EOS_DATA_ROOT` guard: **76 currently registered tests**. The imported
preflight independently records that count (`PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md:763-790`).

The historical 77 count belongs to a different inventory/branch: the noncanonical Phase-6A-1
recovery topology removed one canonical checkpoint-output registration and added two controlled
BNV registrations, a net increase of one. Historical 77-identity campaign records are not
rewritten. Future cluster documents must say “76 currently registered tests at canonical SHA
`232565a32303a4953f3f516d1d5286b6663f8f99`” and identify any later inventory SHA explicitly.

## 17. ADR and status synchronization

Cross-platform dependency resolution changes the build/dependency architecture and therefore
requires ADR-0018. The proposal is
`docs/adr/ADR-0018-cross-platform-external-dependency-resolution.md`, status **PROPOSED — OWNER
RATIFICATION REQUIRED**. It carries no implementation authority.

No current-architecture, roadmap, or scientific-invariant status is changed by an unratified
proposal. The durable proposed state is recorded here instead:

- cluster environment: DISCOVERED / NOT QUALIFIED;
- CompactStar Linux build: BLOCKED PENDING CROSS-PLATFORM DEPENDENCY SUPPORT;
- ADR-0018: PROPOSED;
- Zaki SHA: CANDIDATE ONLY; CONFIND SHA: UNRESOLVED;
- Mac source equivalence: NOT RUN;
- cluster bootstrap: NOT AUTHORIZED;
- Slurm jobs: 0.

## 18. Blockers and owner decisions

**Blocking facts**

1. No source matching the vendored/current CONFIND header and ABI has been authenticated.
2. Current source contains a C++20 `using enum` extension under the declared C++17 contract.
3. ADR-0018 is proposed, not ratified.
4. No Mac source-equivalence experiment has run.

**Owner decisions required**

1. Ratify or revise ADR-0018 and its Model B/fail-closed override contract.
2. Authorize focused recovery of the CONFIND source revision that matches the vendored API;
   do not approve `89c5d9b...` as the candidate.
3. Approve Zaki `b9ddebad...` only as an LB1 equivalence candidate.
4. After a matching CONFIND candidate is found, explicitly approve both SHAs for the Mac-only
   source-equivalence experiment.
5. Separately adjudicate the `using enum` C++17 conformance repair versus a C++20 policy.

## 19. Final controls

This preflight accessed no cluster, created no cluster directory, ran no remote job or build,
and performed no remote source mutation. It modified neither Zaki nor CONFIND. CompactStar
production source, CMake, tests, baselines, EOS/data, and literature are unchanged. The
imported discovery file remains byte-identical to its authenticated transfer.

**Exact recommended next action:** return ADR-0018 and this blocked bootstrap/equivalence plan
to the human owner. The next task, only after explicit approval, should be Mac-only ADR-0018
ratification plus recovery/authentication of a CONFIND source snapshot matching the vendored
contract and preparation—not execution—of the source-equivalence campaign.
