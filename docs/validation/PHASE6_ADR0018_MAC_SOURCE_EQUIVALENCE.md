# Phase-6 ADR-0018 Mac source-build equivalence result

**Status:** STOPPED AT EXTERNAL CANDIDATE BUILD GATE — SOURCE AUTHORITY NOT QUALIFIED

**Disposition B:** **EXTERNAL CANDIDATE BUILD BLOCKED WITHOUT SOURCE
MODIFICATION — SOURCE AUTHORITY NOT QUALIFIED.**

**Date:** 2026-09-26

## 1. Authority and immutable sequencing

| Item | Identity |
|---|---|
| original canonical entry | `232565a32303a4953f3f516d1d5286b6663f8f99` |
| ADR-0018 acceptance/canonical | `812463ac9ed374f64ac9cadd500066ab723d3a6c` |
| predeclaration | `e0879af0cbbc01bdf4094ad597073fefdf8ee225` |
| override implementation | `fc55386ae7549deda4665923d0623e02cc17878a` |
| experiment branch | `physics/adr0018-mac-source-equivalence` |
| experiment worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-adr0018-mac-source-equivalence` |

The predeclaration was committed and pushed before any candidate configure or
build.  The override implementation was then committed and pushed separately.
No candidate-specific CompactStar configuration, treatment result, test, or ODE
existed before either commit.

ADR-0018 remains **ACCEPTED / HUMAN-RATIFIED**.  Its build mechanism is present
only on the unmerged experiment branch.  Neither candidate has been accepted
as a governed dependency source authority.

## 2. Mac platform

| Field | Authenticated value |
|---|---|
| macOS | `26.6.2` build `25G83` |
| Darwin | `25.6.0` |
| architecture | `arm64` |
| compiler | `/usr/bin/clang++`; Apple clang `21.0.0 (clang-2100.3.34.2)` |
| CMake | `/opt/homebrew/bin/cmake`; `4.2.1` |
| GSL | `/opt/local`; `2.7.1` |
| Python selected for candidate builds | `/Library/Frameworks/Python.framework/Versions/3.12/bin/python3`; `3.12.4`; NumPy `2.4.0` |
| build type | `Debug` |
| C++ standard | C++17; unchanged |

The OS, Darwin, architecture, compiler, CMake, GSL, and Debug configuration
match the qualified ADR-0017 Mac reference.  Candidate commands explicitly
removed the ambient Conda and Homebrew compiler selectors and selected Apple
clang plus the declared Python.  No unsafe floating-point option was used.

## 3. Minimal override implementation

The branch adds exactly these cache variables:

- `COMPACTSTAR_ZAKI_LIBRARY`;
- `COMPACTSTAR_ZAKI_INCLUDE_DIR`;
- `COMPACTSTAR_CONFIND_LIBRARY`;
- `COMPACTSTAR_CONFIND_INCLUDE_DIR`.

Static review established:

- no-variable Darwin mode resolves to the same historical
  `dependencies/lib/{Zaki,Confind}/Darwin/arm64/*.a` and
  `dependencies/include` paths;
- any explicit override requires all four values;
- non-Darwin mode requires all four;
- archive paths must be absolute existing regular `.a` files accepted by
  `${CMAKE_AR} -t`;
- include paths must be absolute existing directories containing the required
  Zaki and CONFIND sentinel headers;
- there is no `find_library`, `find_path`, system-path probing, FetchContent,
  network access, or fallback;
- resolved mode and paths are printed;
- link order, compile definitions, floating-point flags, and C++17 are
  unchanged.

The stop occurred before C0/C1, so scientific neutrality of the mechanism has
not been experimentally established.  The implementation is not merged.

## 4. Candidate identities reauthenticated

### 4.1 Zaki

| Field | Value |
|---|---|
| repository | `/Users/keeper/Documents/CompactStar/external/Zaki` |
| Git SHA | `b9ddebaded24962468954846f47238aec2726fd4` |
| tracked state | clean before and after |
| tracked paths | 121 |
| source-manifest SHA-256 | `a8ed2907812354ab2d46bc941806f7bec1e681eaaba2b7d785dc3de60c0e842b` |
| source modified | **NO** |

### 4.2 CONFIND

| Field | Value |
|---|---|
| preserved root | `/Users/keeper/Documents/CompactStar/external/recovered/CONFIND/ed76163c22e0a1f8/source` |
| exact source manifest | `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083` |
| transfer ZIP SHA-256 | `be41ee1c71b88627f0b11703c9d732abd7378ebbc162ee10649af771dec9ae0f` |
| preserved source modified | **NO** |
| divergent current Git repository HEAD | `89c5d9b731534e4289d9f686549d9f0ac178e567` |
| current Git repository tracked state | clean and unchanged |

## 5. Zaki candidate build

Zaki was configured from the exact Git worktree into external disposable build
space.  The source tree was not edited.  Its repository CMake target `Zaki` was
built in Debug with Apple clang and GSL 2.7.1.  The target compiled these 25
members:

`Instrumentor`, `MemoryManager`, `Profile_Timer`, `Logger`, `ObjObserver`,
`Simple_Timer`, `Banner`, `Directory`, `String_Basic`, `TextBox`,
`CSVIterator`, `CSVRow`, `DataColumn`, `DataSet`, `IntegrateTrapz`, `TempDC`,
`Func2D`, `Math_Core`, `NDimContLevel`, `IntegralTable`, `Newton`, `Constants`,
`Coordinate`, `DateTime`, and `Sun`.

| Field | Result |
|---|---|
| build | PASS |
| archive | `/Users/keeper/Documents/CompactStar/external/equivalence/ADR0018/20260926-812463a-a8ed2907-ed76163c/dependencies/Zaki/build/libZaki.a` |
| size | 10,174,504 bytes |
| SHA-256 | `5d80e83a00cb74fa97aae650df61a47848c304605c49262b3f4dce81da2e5832` |
| source modification | none |
| warnings | three Python 3.12 deprecation warnings from `matplotlibcpp.hpp`; no error |

The archive has a `__.SYMDEF` table followed by the 25 object members named
above.  Full compile commands and source/object evidence remain under the
evidence root.  The archive was not copied to a vendored path and is not a
governed dependency artifact.

## 6. CONFIND candidate build and hard stop

The external wrapper requested exactly the five recovered implementation files
predeclared for the consumed contract:

- `source/Base.cpp`;
- `source/Cell.cpp`;
- `source/Cont2D.cpp`;
- `source/Common.cpp`;
- `source/ContourFinder.cpp`.

It used the recovered include tree, candidate Zaki headers and archive, Apple
clang, strict C++17, Debug, GSL 2.7.1, zlib, the declared Python/NumPy, and
OpenMP.  It did not request or link ROOT.

The exact source pair does not compile.  The blocking diagnostics are:

| Recovered source | Lines | Candidate Zaki conflict |
|---|---:|---|
| `Cont2D.cpp` | 109, 110 | direct access to `DataColumn::label`, private in candidate `DataColumn.hpp:211` |
| `Cont2D.cpp` | 114, 115 | direct access to `DataColumn::vals`, private in candidate `DataColumn.hpp:214` |
| `ContourFinder.cpp` | 1381, 1386 | direct access to `DataColumn::label`, private in candidate `DataColumn.hpp:211` |

`Base.cpp`, `Cell.cpp`, and `Common.cpp` compiled before the parallel build
stopped.  `Cont2D.cpp` and `ContourFinder.cpp` failed.  No
`libConfind.a` was created, so no CONFIND archive SHA/member list exists.

This is the concrete load-bearing effect anticipated by the historical Zaki
header analysis: `DataSet.hpp`/`DataColumn.hpp` are not merely plotting or
logging surfaces.  The recovered CONFIND implementation depends on the older
public data-container convention.  The accepted current Zaki candidate exposes
a materially different source API to this exact CONFIND snapshot.

Resolving the error requires at least one unauthorized action: editing current
Zaki candidate headers, editing recovered CONFIND source, substituting an older
Zaki header/source snapshot, or writing a compatibility implementation.  The
predeclaration says such a need is a hard stop.  None was attempted.

## 7. Qualification matrix

| Gate | C0 | C1 | T | Reason |
|---|---|---|---|---|
| configure/build | NOT RUN | NOT RUN | NOT RUN | stopped before CompactStar configuration |
| override dependency identity | N/A | N/A | N/A | no CMake cache created |
| focused dependency contract | NOT RUN | NOT RUN | NOT RUN | exact candidate pair cannot link |
| Phase-5B | NOT RUN | NOT RUN | NOT RUN | downstream of candidate-build gate |
| Phase-5C | NOT RUN | NOT RUN | NOT RUN | downstream of candidate-build gate |
| Phase-5D | NOT RUN | NOT RUN | NOT RUN | downstream of candidate-build gate |
| bounded ADR-0017 | NOT RUN | NOT RUN | NOT RUN | downstream of candidate-build gate |

Consequently there are no C0/C1 or C1/T scientific comparison results, final
states, step counts, internal-history hashes, reconstruction results, or R20
values from this experiment.  The accepted historical ADR-0017 values were
predeclared but not regenerated.  BA12, BA12R, clean BA12R, CQ0-CQ7, and all
cluster work remained untouched.

## 8. Authority protection

| Protected authority | Before | After / result |
|---|---|---|
| vendored Zaki archive | `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef` | identical |
| vendored CONFIND archive | `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2` | identical |
| dependency headers | canonical entry bytes | no changed path |
| 11 governed baselines | canonical entry bytes | 11/11 unchanged by branch diff |
| 33 protected Phase-5D paths | canonical entry bytes | 33/33 unchanged by branch diff |
| Phase-5 scientific sources | canonical entry bytes | unchanged |
| Phase-6 scientific sources | canonical entry bytes | unchanged |
| `ScaledRKF45` | canonical entry bytes | unchanged |
| `Cstar` | canonical entry bytes | unchanged |
| EOS/data and literature | canonical entry bytes | unchanged |

The only permanent branch changes before this result record are the committed
predeclaration and `CMakeLists.txt` override mechanism.  No external source
repository or preserved source was modified.  No cluster was accessed, no
file was transferred to a cluster, and zero Slurm jobs were submitted.

## 9. Evidence

Evidence root:

`/Users/keeper/Documents/CompactStar/external/equivalence/ADR0018/20260926-812463a-a8ed2907-ed76163c/`

It contains the Zaki candidate archive, Zaki source and compile-command
manifests, the external CONFIND wrapper, CONFIND compile commands, partial
object evidence, and `BUILD_FAILURE.txt`.  C0/C1/T build and output roots were
created as empty containers but were never configured or used.

Observed setup/build timings are operational provenance only: Zaki configure
about 2.7 s, Zaki build a few seconds, CONFIND configure about 1.7 s, and the
parallel CONFIND failure within one second.  They are not scientific outputs.

## 10. Conclusion and owner decision

The authorized pair

- Zaki `b9ddebaded24962468954846f47238aec2726fd4`, and
- CONFIND manifest
  `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`

is **not qualified** for source authority.  There is no evidence of a
dependency-induced numerical difference because execution correctly stopped
earlier at source/API compatibility.

The exact next action is an owner architecture decision, not another build.
The owner must choose whether to recover a Zaki source snapshot matching the
recovered CONFIND public-member convention, authorize a bounded compatibility
implementation/source adaptation with a new predeclaration, or replace/remove
CompactStar's CONFIND dependency through a separate architecture decision.
ADR-0018 itself need not be revoked, but this candidate pair cannot advance to
Linux authority or cluster bootstrap.
