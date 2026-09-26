# Phase-6 CONFIND source-authority recovery

## 1. Status and authority

**Status:** SOURCE-FORENSICS COMPLETE — OWNER DECISION REQUIRED.

**Date:** 2026-09-26.

**CompactStar authority examined:**
`232565a32303a4953f3f516d1d5286b6663f8f99`.

**Bootstrap-preflight parent:**
`939f278e806da6ece7db62b8fe6d91a3cb5e38f3`.

**Recovery classification:** **B — EXACT CONFIND SOURCE SNAPSHOT RECOVERED WITHOUT GIT
REVISION.**

**Final disposition:**

> **EXACT CONFIND SOURCE SNAPSHOT RECOVERED — READY FOR OWNER DECISION ON SOURCE
> PRESERVATION AND EQUIVALENCE TEST.**

This record authenticates a non-Git source snapshot whose exact public headers and structurally
matched implementation expose every CONFIND ABI/API entry currently consumed by canonical
CompactStar. It does **not** establish that the recovered implementation bytes were the exact
inputs used to create the authenticated arm64 archive. It does not ratify the snapshot, ADR-0018,
or either external dependency. No build, test, ODE, production edit, cluster access, or checkout
of an historical external-repository revision occurred.

## 2. Entry authentication and boundaries

The canonical repository at `/Users/keeper/Documents/CompactStar/repo/CompactStar` was clean on
`master`. `HEAD`, local `master`, local `origin/master`, and live `origin/master` all resolved to
`232565a32303a4953f3f516d1d5286b6663f8f99`.

The existing worktree and branch were clean and exactly at the authorized preflight SHA, so they
were reused:

- branch: `docs/phase6-linux-cluster-bootstrap-preflight`;
- worktree:
  `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6-linux-cluster-bootstrap-preflight`;
- entry SHA: `939f278e806da6ece7db62b8fe6d91a3cb5e38f3`.

The complete bootstrap preflight, proposed ADR-0018, imported EKU discovery, CompactStar CMake,
vendored headers, active caller, external CONFIND history, and recovered local material were read.
ADR-0018 remains **PROPOSED — OWNER RATIFICATION REQUIRED**. This task did not access EKU or any
other cluster and did not touch the independent C++17/`using enum` blocker.

## 3. Authenticated vendored contract

### 3.1 Archive identity

| Field | Authenticated value |
|---|---|
| Path | `dependencies/lib/Confind/Darwin/arm64/libConfind.a` |
| SHA-256 | `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2` |
| Size | 703,208 bytes |
| Format | Darwin arm64 static archive containing Mach-O arm64 objects |
| Public contract version header | `CONFIND_VERSION_STR "1.0"`; release date `09, 21, 2023` |

The repository working-copy timestamps are 2026 checkout metadata and are not provenance. A
separate retained Google Drive copy of this archive has the same SHA-256 and preserves mtime
`2023-09-16T11:34:14-0400`; that timestamp is a clue, not source authority.

### 3.2 Header contract

Canonical CompactStar directly includes `Confind/ContourFinder.hpp`. The following six public
headers and generated config form its transitive CONFIND contract and the library's matching
public/build interface. Files with a ` 2` duplicate name in the vendored directory are
byte-identical to the corresponding canonical name and add no second identity.

| Vendored path | SHA-256 | Relevant public content |
|---|---|---|
| `dependencies/include/Confind/Base.hpp` | `3e089b6d9429edd9dc859956fe1ce63e632dde1dfef1cacac273df4b3f9d2e7b` | `CONFIND::Base`; work-directory API |
| `dependencies/include/Confind/Bundle.hpp` | `8a85cdd53267da5882c158bc49d46d33b27f82915f8f21c9f4aa6dbd9cc0b13c` | thread-work bundle; inline constructors and Zaki grid/function types |
| `dependencies/include/Confind/Cell.hpp` | `b4b7c155fdf3439dfb4a835954cb3a9a7ca40a1628954776d4c26c73611ffea4` | contour cell/vertex/triangle internals used by `ContourFinder` |
| `dependencies/include/Confind/Common.hpp` | `aaeb12aed6f620b10034b10d79faccd04313bda48cf77ba8247c2d47edf22bb8` | `CONFIND::Color` |
| `dependencies/include/Confind/Cont2D.hpp` | `b22cb1e6ebe62d3c91d3cd17b3de627e549f5142e290ec18b06553895cd15a7f` | `CONFIND::Cont2D`, point/value/curve/data conversion |
| `dependencies/include/Confind/ContourFinder.hpp` | `6b95b0a30c5e9f03aee022ed4e1994a3aed34b470f31472c90a550e6f20a09dd` | `CONFIND::ContourFinder`, grid/contour/plot/export API |
| `dependencies/include/Confind/ConfindConfig.h` | `da6209f01d0fc7c0b7dbf620696885fd0a19a156b2c9c85103692c6111e6f114` | version `1.0`, release date `09, 21, 2023` |

### 3.3 CompactStar-consumed signatures

The only active CompactStar caller is `CompactStar/Core/src/TaskManager.cpp`. Comments were not
counted as calls. Its authenticated contract is:

| Class/header | Exact callable signature represented by the header/archive | Active purpose |
|---|---|---|
| `ContourFinder` / `ContourFinder.hpp` | `ContourFinder()` and `~ContourFinder()` | construct/destroy finder |
| inherited `Base` / `Base.hpp` | `void SetWrkDir(const Zaki::String::Directory&)` | select output root |
| `ContourFinder` | `void SetGrid(const Zaki::Math::Grid2D&)` | configure two-dimensional grid |
| `ContourFinder` | `void SetContVal(const std::vector<double>&)` | configure contour levels |
| `ContourFinder` | `void SetGridVals(Zaki::Math::GridVals_2D*)` | provide precomputed grid values |
| `ContourFinder` | `void SetPlotConnected(const bool = true)` | request connected contour plotting |
| `ContourFinder` | `void Plot(const Zaki::String::Directory&, const char* const = nullptr, const char* const = nullptr, const char* const = nullptr)` | generate contour plot through Zaki `DataSet` plotting |
| `ContourFinder` | `std::vector<CONFIND::Cont2D> GetContourSet() const` | retrieve contour collection |
| `ContourFinder` | `void ExportContour(const Zaki::String::Directory&, const Zaki::File::FileMode&)` | export contour data |
| `Cont2D` / `Cont2D.hpp` | `size_t size() const` | point count |
| `Cont2D` | `double GetVal() const` | contour level value |
| `Cont2D` | `Zaki::Physics::Coord3D operator[](const size_t&) const` | indexed point access |
| `Cont2D` | `Zaki::Math::Curve2D ConvertToCurve2D()` | convert for intersections and downstream analysis |

The authenticated arm64 archive exports the corresponding demangled symbols with these types,
including libc++ `std::__1` container spellings. `SetGridVals`, `GetContourSet`, and
`ConvertToCurve2D` are therefore facts from both the header and binary, not inferred APIs.

## 4. Search inventory

Read-only searches covered reasonable user-owned locations under `/Users/keeper/Documents`,
`Downloads`, `Desktop`, relevant `/Users/keeper/Library/CloudStorage` roots, CompactStar-related
directories, Git repositories, transferred backups, ZIP/TAR member names, old archives, and
source/build filenames. Search keys included case variants of CONFIND, `ContourFinder`,
`Cont2D`, `SetGridVals`, `GetContourSet`, and `ConvertToCurve2D`. System/private OS areas and
shell credential material were excluded.

| Candidate | Result |
|---|---|
| `/Users/keeper/Documents/CompactStar/external/CONFIND` | 2019-2021 Git history; old `.h`, lowercase-Zaki, ROOT-facing contract; **MATERIAL_DIVERGENCE** |
| `/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/Confind` | five-file 2023 implementation plus build system and GPL-3.0 license; selected implementation snapshot |
| `/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/dependencies/include/Confind` | all selected public headers byte-identical to canonical; selected header snapshot |
| same retained CompactStar tree, `dependencies/lib/Confind/Darwin/arm64/libConfind.a` | byte-identical canonical archive; strongest adjacency/provenance corroboration |
| same retained CompactStar tree, `dependencies/lib/Confind/Darwin/x86_64/libConfind.a` | 2023-09-21 archive from the later source/header epoch; corroborating, not Darwin-arm64 authority |
| retained `Research/Tools/Coding/DMSS` and 2019-2021 DMSolarSignal trees | older header/archive variants; only partial matches |
| retained `Confind/Backup/2020` Git object store and dated ZIPs | old reachable/public design; no authenticated 2023 contract |
| manually copied Downloads Google Drive tree | incomplete copy/placeholders; no additional candidate bytes |
| canonical CompactStar and its worktrees | authenticated consumer headers/archive only; no CONFIND implementation source |

Some Google Drive entries were dataless placeholders. Reading them caused the local sync provider
to hydrate content without changing file bytes; resulting 2026 ctimes are access/materialization
metadata and are explicitly not provenance. No content was written to the snapshot.

## 5. Local and remote Git archaeology

The external repository remained at clean `master`, SHA
`89c5d9b731534e4289d9f686549d9f0ac178e567`, tree
`1f5a154835e73772d383a03dc18daa7781dc1763`.

- All 64 reachable commits, all local/remote refs, tags, and reflogs were inspected.
- The reachable range is `f6f34a983b4498734ae1500d967d3dc2e8e94fa6` (2019-12-26) through
  current HEAD (2021-09-27).
- `git fsck --full --no-reflogs --unreachable` reported no unreachable objects. No
  `--lost-found` write was performed.
- Historical `SetGridVals` declarations exist in the old interface, but no reachable tree has
  `GetContourSet` or `ConvertToCurve2D`; none has the exact `.hpp` header set.
- No checkout, reset, fetch into the working repository, or ref mutation occurred.

Read-only `git ls-remote --symref` queries to both the canonical URL
`https://github.com/ZAKI1905/CONFIND.git` and the historical
`https://github.com/ZAKI1905/ContourFinder.git` advertised only `master` at
`89c5d9b731534e4289d9f686549d9f0ac178e567`; no remote tags or additional branches were
advertised. Local reachable history therefore contains every currently advertised remote ref.
The required 2023 state is not remotely reachable. Whether it was never pushed or remote history
was later rewritten is **UNKNOWN**.

## 6. CompactStar history

The current headers and both Darwin architecture archives entered CompactStar together in:

| SHA | Tree | Date | Subject | Relevance |
|---|---|---|---|---|
| `7fd01327a919f736e71b2385f1e3970110810b0f` | `3eceb6cb1bc7d6c5614e6dbef5474c978ed98e44` | 2025-04-17 | `uploaded on Git` | first introduction of current CONFIND headers and arm64/x86_64 archives |
| `dfb44433dd56cc820889ed07de079e320db18032` | `1301b4a3ee475daa2e762f9a74463b69588628ca` | 2025-11-27 | `Second major release of CompactStar...` | reorganized project; removed dependency `.DS_Store`; did not replace CONFIND authority bytes |

The introducing commit contains no CONFIND implementation source and no repository SHA/tag or
build record. No older CompactStar commit contains the selected 2023 implementation. Git history
therefore establishes the artifact introduction point but not an external source revision.

## 7. Archive forensics

### 7.1 Members and timestamps

| Member | Size | Preserved archive timestamp |
|---|---:|---|
| `__.SYMDEF` | 43,560 | 2023-09-16 11:34 |
| `Base.cpp.o` | 109,312 | 2023-09-16 11:33 |
| `Cell.cpp.o` | 119,512 | 2023-09-16 11:33 |
| `Cont2D.cpp.o` | 196,960 | 2023-09-16 11:33 |
| `Common.cpp.o` | 10,920 | 2023-09-16 11:33 |
| `ContourFinder.cpp.o` | 222,448 | 2023-09-16 11:33 |

Members retain UID/GID `501/20`. Every object is Mach-O arm64 and declares macOS platform 1,
minimum OS 13.0, SDK 13.0. libc++ `std::__1` mangling is present. The objects contain no DWARF
producer record, compiler-version string, Git SHA, or CONFIND source-version string, so the exact
compiler release is **UNKNOWN**.

### 7.2 Source-path and link clues

Embedded strings disclose historical build paths including:

- `//Users/keeper/Library/CloudStorage/GoogleDrive-mzake001@ucr.edu/My Drive/Work/Coding/Confind/src/Base.cpp`;
- the analogous `Cell.cpp`, `Cont2D.cpp`, and `ContourFinder.cpp` paths;
- `.../Confind/include/Confind/Bundle.hpp`;
- `.../Confind/dependencies/include/Zaki/File/VecSaver.hpp`.

The exact former UCR Drive path is no longer present locally. The selected EKU Drive snapshot has
the same project structure and file names. This is strong lineage evidence, not byte proof that
the later implementation files were the arm64 compiler inputs.

Undefined-symbol inspection shows Zaki APIs, OpenMP runtime calls, and zlib `compress`. It shows
no ROOT graphics/runtime symbol. The active archive therefore structurally excludes ROOT from
the consumed link path. No absolute Git identity or build command is embedded.

## 8. Recovered source snapshot

### 8.1 Composition

The source root retains all five implementation files and the CMake build description, but its
`include/Confind` directory currently retains only generated config files. Exact public headers
are retained in the adjacent historical CompactStar consumer tree and are byte-identical to the
canonical contract. Accordingly, the candidate is an explicit **composite immutable snapshot**:

- implementation/build root:
  `/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/Confind`;
- exact header root:
  `/Users/keeper/Library/CloudStorage/GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/dependencies/include/Confind`.

The implementation provides every library-owned `.cpp`; the header root provides every
library-owned public/inline header; the source root provides CMake lists, config template,
README, and GPL-3.0 license. Zaki, GSL, zlib, OpenMP, and Python/NumPy are external dependencies
and are not incorporated into CONFIND source identity.

### 8.2 Complete candidate manifest

The following is the complete manifest for the candidate's CONFIND-owned source, public headers,
build declarations, README, and license. Each line is `SHA-256`, byte size, and logical path.
The SHA-256 of this newline-terminated, logical-path-sorted manifest is
`ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`.

```text
fd5219e5010182cc74cd75c6564f3895119ba65c7848c0f0179bbdcf05951acb  880  README.md
9e1b1e95422e175400e8c5e51df54100024533a5d030268072a6c16f55567c56  7001  build/CMakeLists.txt
d6ac63e3a1f7e1060b2badd420514764cd71d6bded79965749400d18feee028f  346  build/ConfindConfig.h.in
2862cc18f330c9f401c496fb9a56b0ce5bc0a01c5cc86b1cc5ed29e33d8c4e96  74  build/include-CMakeLists.txt
963226b1c056908776469da7f8ec8fc769ef613a28384c58c3eaa9e19a3f69fd  138  build/src-CMakeLists.txt
3e089b6d9429edd9dc859956fe1ce63e632dde1dfef1cacac273df4b3f9d2e7b  2633  include/Confind/Base.hpp
8a85cdd53267da5882c158bc49d46d33b27f82915f8f21c9f4aa6dbd9cc0b13c  3856  include/Confind/Bundle.hpp
b4b7c155fdf3439dfb4a835954cb3a9a7ca40a1628954776d4c26c73611ffea4  7153  include/Confind/Cell.hpp
aaeb12aed6f620b10034b10d79faccd04313bda48cf77ba8247c2d47edf22bb8  692  include/Confind/Common.hpp
da6209f01d0fc7c0b7dbf620696885fd0a19a156b2c9c85103692c6111e6f114  244  include/Confind/ConfindConfig.h
b22cb1e6ebe62d3c91d3cd17b3de627e549f5142e290ec18b06553895cd15a7f  2475  include/Confind/Cont2D.hpp
6b95b0a30c5e9f03aee022ed4e1994a3aed34b470f31472c90a550e6f20a09dd  7368  include/Confind/ContourFinder.hpp
3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986  35149  legal/LICENSE
6e74812f8eab42063f1cac6fe4418c2a80c0328ae3115741e26f97b1f03f6ba1  5450  source/Base.cpp
bc58cd5d51c8852d7d016424fade34841b015eeddbcad9fa52669f10f391f670  15456  source/Cell.cpp
942a3b6d7b4f7aa62acde01f73b9788d938f2056ce30a69e8b0430ab338ec257  1762  source/Common.cpp
fca14d98fa4fc9646883eebcf762f8b4b09c21d02808ec3e913f2d5f9ce9497e  10318  source/Cont2D.cpp
8eaf44570eca16b1a7387e06416d3501fd8ddb674f4321e02c228847e31d8b6d  47834  source/ContourFinder.cpp
```

No Git SHA or Git tree can honestly be assigned to this composite snapshot.

### 8.3 Metadata limits

The two most relevant implementation files have preserved mtimes on 2023-09-21, five days after
the arm64 member timestamps: `Cont2D.cpp` at 15:11 and `ContourFinder.cpp` at 16:19. The exact
headers are also timestamped 2023-09-21. A retained x86_64 archive was built at 16:19-16:21 that
day and corroborates that epoch. These facts support snapshot cohesion but prevent claiming that
the later bytes were the exact 2023-09-16 arm64 compiler inputs.

## 9. Match scoring

### 9.1 Headers

Every candidate CONFIND header in section 8 is byte-identical to the authenticated canonical
header in section 3. Classification: **EXACT_HEADER_MATCH** for all six public headers and config.

The current Git repository and every reachable remote revision remain
**MATERIAL_DIVERGENCE**. Older retained DMSS/DMSolarSignal snapshots are **PARTIAL_MATCH**: some
base headers match, but `Common.hpp`, `Cont2D.hpp`, `ContourFinder.hpp`, and config differ.

### 9.2 Implementation and symbols

The recovered sources define every active CompactStar-consumed signature in section 3 with the
same CONFIND/Zaki types, and every one has a corresponding exported arm64 symbol. Classification
for the **consumed** interface: **COMPLETE_SYMBOL_CONTRACT_MATCH**.

This is not whole-archive identity. The 2023-09-21 source defines
`Cont2D::ConvertToDataSet()` and the same-day x86_64 archive exports it, while the authenticated
2023-09-16 arm64 archive does not export that symbol even though the public header declares it.
CompactStar does not call `ConvertToDataSet` directly; recovered `ContourFinder::Plot` does.
That delta reinforces the requirement for build/link and behavior equivalence before authority
can be accepted.

## 10. Zaki convention and dependency fit

The recovered candidate uses the authenticated convention required by CompactStar:

- uppercase `Zaki/...` include paths;
- namespace spellings `Zaki::Math`, `Zaki::Physics`, `Zaki::String`, `Zaki::File`,
  `Zaki::Vector`, and `Zaki::Util`;
- public types `Grid2D`, `GridVals_2D`, `Curve2D`, `Coord3D`, `Directory`, `FileMode`, and
  `DataSet`.

Of the 15 Zaki headers directly included by the recovered CONFIND headers/sources, 12 retained
snapshot copies are byte-identical to CompactStar's authenticated Zaki headers. Three differ:
`Util/Instrumentor.hpp`, `Util/Logger.hpp`, and `Vector/DataSet.hpp`. `DataSet` is relevant to
`Plot`; the other two affect instrumentation/logging integration. The accepted future experiment
must therefore build the recovered CONFIND snapshot against the separately governed Zaki
candidate `b9ddebaded24962468954846f47238aec2726fd4`, not silently use the snapshot's bundled Zaki
archive/headers. Compatibility and link closure remain experiment gates.

Classification: **CONSUMED ZAKI NAMING/TYPE CONVENTION MATCH; DEPENDENCY-VERSION EQUIVALENCE
REQUIRED**. Zaki's status remains **CANDIDATE FOR SOURCE-EQUIVALENCE TESTING ONLY**; it was not
modified or ratified here.

## 11. ROOT and plotting separation

ROOT includes, objects, accessors, and drawing code are commented out in the exact headers and
recovered implementation. ROOT discovery/linking is commented out in the recovered CMake file,
and the authenticated archive has no undefined ROOT symbol. **ROOT is not required by the
CompactStar-consumed path.**

The active `ContourFinder::Plot` implementation instead converts contours to
`Zaki::Vector::DataSet` and calls its PDF plotting methods. The recovered CMake file requires
Python development components and NumPy “for matplotlib.” Because canonical CompactStar actively
calls `Plot`, Python/NumPy and the accepted Zaki `DataSet` plotting implementation are part of the
candidate link/runtime qualification unless a separate, explicitly authorized architecture
change removes that call. This task made no such change.

## 12. Candidate decision and limitations

**Selected status:** exact non-Git source snapshot, manifest
`ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`.

**Header confidence:** exact.

**Consumed structural ABI confidence:** high; complete symbol-contract match.

**Historical arm64 source identity confidence:** moderate, not exact. The archive embeds the same
project structure but predates the two latest source/header mtimes by five days.

**Scientific/behavioral confidence:** unknown until the predeclared same-Mac source-build
equivalence experiment runs. Source recovery is not behavior qualification.

The source set is complete enough to attempt a reproducible library build after preservation:
all CONFIND-owned implementation, inline/public headers, config template, build lists, and license
are present. It must first be copied byte-for-byte into an owner-approved durable recovery
repository/branch or immutable source bundle with this manifest. No such history or bundle was
created here.

## 13. ADR-0018 consequence

The Model-B, fail-closed dependency-resolution decision does not need architectural expansion.
However, the proposed ADR currently requires a human-ratified **Git source SHA**. The recovered
CONFIND authority has no defensible Git revision. Factual consistency therefore requires a
narrow wording adjustment, without ratification:

> An external dependency's exact source identity is a human-ratified Git commit SHA when one is
> available; otherwise it is a human-ratified immutable source-manifest SHA-256 covering every
> source, public/inline header, build declaration, and license file. Either identity is a
> qualification key, and any member-byte change invalidates qualification.

ADR-0018 remains **PROPOSED — OWNER RATIFICATION REQUIRED**. The adjustment does not select or
accept this snapshot and does not authorize LB1.

## 14. Exact next action

Return this manifest and ADR-0018 to the human owner. The owner must first decide whether to
preserve the composite snapshot byte-for-byte in a dedicated CONFIND recovery branch/repository
or immutable source bundle, retaining the manifest identity. Only after preservation should the
owner explicitly authorize a **Mac-only** source-build equivalence experiment for both:

- Zaki candidate Git SHA `b9ddebaded24962468954846f47238aec2726fd4`;
- CONFIND candidate source-manifest SHA
  `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`.

That future experiment must use the unchanged vendored Darwin build as control, build candidate
archives separately, supply explicit override paths, verify link closure including Plot/Python
behavior and no ROOT requirement, and apply the already predeclared Phase-5B/C/D and bounded
ADR-0017 equivalence criteria. It may not overwrite the vendored archives. This record itself
authorizes none of those actions.

## 15. Negative attestation

- Cluster/SSH access: **0**.
- Builds: **0**.
- Tests/ODEs: **0**.
- CompactStar production/CMake changes: **none**.
- Test, baseline, EOS/data, and literature changes: **none**.
- External CONFIND Git repository changes: **none**.
- External Zaki Git repository changes: **none**.
- ADR-0018 ratification: **no**.
- C++ standard work: **not touched**.
