# CompactStar external dependency migration qualification

Status: complete local Mac migration candidate, pending one combined
independent review and subsequent owner acceptance/integration.
This record establishes no canonical merge or release.

Supporting records: [plotting ledger](COMPACTSTAR_PLOTTING_REMOVAL_LEDGER.md),
[ADR-0019](../adr/ADR-0019-canonical-source-built-dependency-authority.md),
[append-only predeclaration](PHASE6_EXTERNAL_DEPS_TASKMANAGER_MIGRATION_PREDECLARATION.md),
and [evidence index](evidence/external-deps-migration/resume/README.md).

## Authority and scope

The owner authorized the exact provisionally qualified Zaki candidate before
one combined independent review, and subsequently replaced object/provider
identity as an automatic stop with scientific and functional equivalence.
The original predeclaration and the earlier G4 stop remain historical records.
No new numerical tolerance, physics equation, scientific constant, EOS method,
TOV algorithm, BNV formula, thermal/rotochemical method, or TaskManager thread
implementation was selected or changed.

CompactStar canonical entry is `812463ac9ed374f64ac9cadd500066ab723d3a6c`.
Worktree: `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration`.
Branch: `physics/external-deps-taskmanager-migration`.
Change classes: dependency/build, migration-specific structural cleanup,
validation and documentation. Mac arm64 only. No merge, tag, Linux or cluster work.

## Active dependency architecture

CMake requires absolute `COMPACTSTAR_ZAKI_PREFIX` and
`COMPACTSTAR_CONFIND_PREFIX`. It authenticates source/build provenance,
archive bytes, every installed header, and package-config files against
`cmake/dependency-lock.json`, then uses exact CONFIG-mode discovery with
`NO_DEFAULT_PATH`. Explicit prefixes override stale package cache entries;
imported target archive and include paths must resolve inside the authenticated
prefix. Seven positive/negative controls cover missing prefixes, mode mixing,
stale cache entries, modified headers, modified archives, wrong source identity,
and an already-created target pointing at an unauthenticated archive.

Active targets are `Zaki::Zaki`, `CONFIND::CONFIND`, direct `GSL::gsl` and
`GSL::gslcblas`, and `Threads::Threads`. Zaki supplies zlib transitively.
Python is an interpreter for configuration/validation, never a C++ link input.
There are no active OpenMP source calls/pragmas, OpenMP link flags, vendored
Zaki/CONFIND include paths or archives, Python/NumPy link targets, or
`-undefined dynamic_lookup` in the final build. Historical artifacts stay
at their original paths and are also excluded from installation.

Qualified environment: AppleClang 21 (`clang-2100.3.34.2`), CMake 4.2.1,
macOS 26.6.2 (25G83), arm64; actual linked GSL is `/opt/local` 2.7.1,
not the miniforge `gsl-config` 2.7 installation. SDK and directly queried runtime zlib are 1.2.12; the GSL runtime also reports 2.7.1.
Python 3.12.10 with NumPy/SciPy/mpmath runs Python-driven tests.
The candidate lock is intentionally specific to these Mac package artifacts;
rebuilds require recorded provenance and renewed qualification.

Package root:
`/Users/keeper/Documents/CompactStar/external/qualification/compactstar-migration/e263a6e-b0cbd510`.
Each mode has separate `Debug/{Zaki,CONFIND}` or `Release/{Zaki,CONFIND}` prefixes.

| Dependency/mode | Source SHA | Archive SHA-256 |
|---|---|---|
| Zaki Debug | `e263a6e180c5c417198e7778bd21fc9c0a32dc33` | `f30e68b261504e1ca3b9af905dc98a0ce70fa25f04c3c19ab81e3b2629fb074f` |
| Zaki Release | same | `3b7f375ff789cd2c0edf58f97274c06fa047e4f477bf327214f9406461baad9d` |
| CONFIND Debug | `b0cbd510fd3fd0c772fa50499cd749287cb39e7b` | `1f36fc6ef6964d4342a025539a5a26c6c55faf04f70d1fe9319c0a4e70cc32c7` |
| CONFIND Release | same | `4cee9b8361fec75023d74b8a68bf7efae54f69cf1464e478173bcfb5d38b524d` |

Installed-header manifest SHA-256 (same in each mode): Zaki
`51931f9eac6823cc21609741c7f60485d0dbfab8d815d0a541e897cbcc2528a2`;
CONFIND `0e02d55747cfa7bce70fdc9c660a4b95ec540dd4ce0c9bd94fbb68bccd5fd99b`.
Manifest format is sorted SHA-256, two spaces, prefix-relative filename and LF.
Package-config hashes and complete compiler/build provenance are in the lock
and evidence manifests. Sources were clean and neither external repository
was modified.

Historical arm64 archive hashes remain:

- Zaki: `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef`.
- CONFIND: `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2`.

## Plotting and link cleanup

The ledger accounts for 75 active plotting calls in nine compiled translation
units: TaskManager 5, DarkCore_Analysis 18, CompOSE_EOS 9, BNV_Analysis 3,
BNV_Sequence 8, Decay_Analysis 6, BNV_B_Chi_Photon 10,
BNV_B_Chi_Transition 2 and BNV_Chi 14. Zero removed-API calls remain.
Two unbuilt legacy demos also lose 12 matplotlibcpp calls. Plot parameter,
axis and legend setup is removed; five CompOSE legend-selection branches
were deleted after their only consumer disappeared. Numerical operations,
including existing smoothing and exports, remain in their original order.

Twelve export sites preserve otherwise unexported derived tables (four in
BNV_Sequence, six in BNV_Chi, two in BNV_B_Chi_Transition). PlotFermiE delegates
to its existing ExportFermiE companion. No export was added for every removed
plot; already available data remains the external visualization input.
The full per-call classification and retention paths are in the plotting ledger.
No plot operation was found to mutate subsequent numerical state. The direct
three-stack experiment (OLD original, OLD with plotting-only deletions, NEW)
reproduces all 13 T1 and 38 T2 files in both modes. Remaining legacy analysis
plot sites have source-level data-flow assessment and build coverage; this is
not a claim that every legacy analysis entry point was exercised at runtime.

Removing permissive undefined-symbol linking exposed two historically
undefined spin diagnostics: CharacteristicAge and DipoleFieldEstimate.
They now throw explicit logic_error exceptions; a test verifies that behavior.
No unqualified formula or dipole normalization was invented. These diagnostics
remain unavailable for scientific use. The root install now owns the public-header tree and generated configuration
header, bypassing stale per-directory install recipes while preserving their
frozen source hashes. The license install path is also corrected. Fresh installations retain all
127 public headers, with exact source/generated-header contents. All 85
static-archive member payloads match the build; CMake install invokes ranlib
and changes only the symbol-index timestamp. The dependency-package archive
locks remain strict whole-file hashes.

## TaskManager and arithmetic qualification

T1 uses the actual FindCriticalCurve, FindMtotContour and FindBtotContour
methods on the frozen real 20x20 stellar sequence. All 13 exported files match.
A separate test-only generated TU reads 53,484 exact hex/bit records without
replacing arithmetic or control conditions: all 92 critical levels, found
flags, raw contour counts/order/x/y/z, every argmax update/index/selected point,
smoothed/final/reimported critical curves, mass_curve[0], intersections,
baryon bounds/levels, GetIdx and ordered Bisect points. Both historical
observer runs also match the 13 ordinary OLD outputs. The production library
has no observer hooks.

The 1,510-record replay freezes real near-tie geometry and selection behavior.
The known unpatched Zaki 2.0.0 executable differs in 22 records in each mode;
the exact authorized candidate differs in zero. This is a discriminating
control, not a golden file regenerated from NEW. Fresh Zaki characterization against the installed packages matches all 162,619 fixed and 251,515 stress records in each
mode (zero differing records).

T2 generates the real DS(CMF)-1 grid with dark mass 0.8 neutron masses,
20x20 logarithmic nodes and target mass 2.01, then runs contours,
Precision_Task and FindLimits. Every run uses one TaskManager thread and a
working root shorter than 100 characters. All 38 deterministic outputs match
OLD byte-for-byte independently in Debug and Release. This includes Sequence,
BNV_rates, BNV_tau and the final neutron, Lambda and Sigma- lifetime limits:
zero differing files and zero differing bytes. Debug is not compared to Release.

Final lifetime-table SHA-256 values (OLD and NEW identical):

| Mode | Species | SHA-256 |
|---|---|---|
| Debug | neutron | `6705e3f78022ec608dcfc344624659d82f263b5deb67c25ac34e9aeabe9315d0` |
| Debug | lambda | `4c5ac872ee4609cc35256dd4dbd4d11a5a0b3941b5343e10587d89fbc571cce2` |
| Debug | sigmam | `1d59f55be6e9a5cbf85976885b54536212a2c0fceb0f8e2ed12719e9c000dac6` |
| Release | neutron | `553fe165b3a2d02602447bc967d12c3093e3148c8301836e619dfc1ef88eebf5` |
| Release | lambda | `f1862cad2957e63a224fc77e6c376e981544f4bb4e56a4605a2a888ee1ac5aa9` |
| Release | sigmam | `484b594d46bfb32b33c7bcdd0a1ae5b64030467d740a86858290a2207e5ac948` |

## Final provider assessment

The final TaskManager-bearing images are taskmanager_migration and the
separate test observer executable. Neither libCompactStar nor another
consumer TU defines a Coord3D comparator. Debug resolves weak XYDist2,
operator< and operator== from the pinned CONFIND Cont2D object. Distance
arithmetic remains subtract/multiply/add without contraction. Release inlines
the comparator inside CONFIND; no consumer override is present. Exact raw
contour/selection records and T2 outputs independently confirm the functional
result.

The full final-image comparison finds four changed provider classes among
1,581 common Debug symbols and five among 784 Release symbols. Directory
moves, vector storage helpers, destructors, assignments and a switch table
have no scientific arithmetic. The one FP-bearing provider change is the
already reported Release VecSaver::Export1D<Coord2D> compression-buffer
allocation calculation. Its FMA lies only in WriteBinary; governed contour
exports use Text. It therefore changes neither scientific values nor buffer
allocation on the qualified path. No general qualification of every legacy
binary-export size is claimed. Full link maps, symbols and disassembly are
retained; provider identity alone is not a gate under the resume authority.

## Scientific regression matrix

OLD Debug passes 76/76 (6,070.05 s); NEW Debug passes 81/81 (6,076.96 s).
Both include the complete long Phase-5D coupled-oracle and Phase-5D1 tests.
OLD/NEW Phase-5D1 receipts in both modes reproduce the governed artifact
byte-for-byte and pass ten comparator controls per invocation. The repaired
NEW Release Phase-5D1 rerun passes in 2,599.13 s.
The final rebuilt source additionally passes 39/39 affected tests in each
mode. All 66 executables per mode are inspected again for link closure.
Relative to the full-suite image inventory, only the repaired heat-capacity
test changes in Release; Debug additionally rebuilds two spin-diagnostic
consumers after header documentation changes. All other executable hashes,
including TaskManager and the Phase-5/6 trajectory executables, are unchanged.
Final Release passes 77/77 in 275.42 s, with the separately completed fresh
Phase-5D1 rerun providing the 78th applicable result. OLD Release has 73/73
applicable results across its full initial attempt and documented repairs.
The three Debug-authority tests remain fully exercised in both Debug suites.
Python-driven registration is 18 OLD tests per mode, 22 NEW Debug and 19 NEW
Release (the latter excludes those three Debug-authority tests).

The first complete Release attempts expose five cross-mode reference failures
on both OLD and NEW. Three have specifically Debug-only authority: Phase-5B
structural-response, Phase-5C coefficient regression, and the coupled-oracle
certificate handshake. CMake now registers those three only in Debug. OLD/NEW
Release generated Phase-5B and Phase-5C artifacts are byte-identical.
The Hartle-monopole and baryon-number comparisons support same-mode Release
references; they remain registered using separately frozen canonical OLD
Release outputs under the migration fixtures. Both NEW emissions match those
files byte-for-byte, the comparison rules are unchanged, and existing Debug
baselines and explicit ADR-0012 overrides are untouched.

The initial NEW Release Phase-5D1 attempt completed numerical generation but
rejected exactly two provenance fields: the Diagnostics and Observers install
CMake hashes. Every scientific field matched. Both CMake files were restored
byte-for-byte, installation moved to the root rule, and a complete fresh
Phase-5D1 rerun passed with exact baseline bytes and all ten controls.
No provenance comparison or baseline was weakened.

The initial OLD Release heat_capacity_v1 aborted because concurrent suites
shared and deleted the same temporary fixture directory. A serial OLD rerun
passes. The candidate test now allocates its own directory atomically; its
numerical assertions remain unchanged. Initial logs and repaired results are
retained separately. No slow test is omitted.

Five hash-only baselines reproduce their committed Debug hashes. All five
also compare exactly OLD versus NEW Release under the same-mode authority.

The bounded ADR-0017 qualification retains its accepted scientific meaning:
one uninterrupted RKF45 main trajectory and 239 strict-interior O1/O2
reconstructions, with existing endpoint and observation qualification.
It does not rerun BA12/BA12R or turn their historical FAIL into a pass, does
not create a BNV baseline, and confers no cluster authority. NEW passes: all
478 O1/O2 states and diagnostics are bit-identical, all 239 strict interiors
qualify with zero failures, and R20 and its reconstruction uncertainty remain
exact. Main accepted/rejected steps remain 232/60. The whole executable takes
8,623.44 s including full provenance checks; its reported checkpoint component
time alone is 1.8617 s. The fresh OLD executable also passes in 8,181.16 s.
OLD, NEW and the accepted historical result are identical in every scientific
result field. All 241 checkpoint records match after excluding only the twelve
explicit wall/CPU timing columns; eight non-timing performance fields also match.
Every main output file is byte-identical OLD versus NEW.

## Performance

Back-to-back measurements (seconds):

| Fixture | OLD Debug | NEW Debug | OLD Release | NEW Release |
|---|---:|---:|---:|---:|
| T1 | 15.740 | 0.057 | 16.242 | 0.019 |
| T2 | 87.729 | 69.759 | 55.609 | 30.130 |
| Phase-5D coupled oracles | 1834.18 | 1825.75 | Debug authority | Debug authority |

Every paired T1/T2 run also rechecks exact scientific outputs. An earlier NEW
Debug T2 timing of 115.3 seconds overlapped the busy parallel part of the suite;
repeat back-to-back runs resolve that contention rather than attributing it
to dependency performance. The representative paired measurements show no
material regression and are below the earlier 1.25x criterion. Long-run
whole-process ADR-0017 elapsed time is 8,181.16 s OLD and 8,623.44 s NEW
(1.054x). Timing data is separated
from deterministic scientific records. Runs overlap other qualification work;
they are not an idle-machine benchmark. The owner resume replaces the earlier
arbitrary 1.25x stop with investigation of material regressions.

ADR-0017's component timing excludes its expensive full provenance rereads.
Whole-process elapsed time is recorded separately; a fresh OLD run uses the
same canonical source, inputs and Mac to distinguish inherited overhead from
a migration regression. No provenance guard is bypassed for speed.

## Gates and disposition

| Gate | Result | Evidence / scope |
|---|---|---|
| G0 package identity | PASS | Exact mode-specific archives, headers/configs/source provenance; seven discovery controls |
| G1 existing tests | PASS | OLD Debug 76/76, NEW Debug 81/81; applicable OLD Release 73/73, NEW Release 78/78; all long tests included |
| G2 plot removal | PASS within governed scope | 75 removed-API ledger rows; OLD original = OLD plot-free = NEW for T1/T2; source/build assessment for other legacy analysis paths |
| G3 Zaki arithmetic | PASS | 162,619 fixed and 251,515 stress records exact per mode; 1,510 replay records exact, unpatched negative control differs in 22 |
| G4 Coord3D/providers | PASS under resume authority | Expected contour arithmetic/provider; no consumer comparator interposition; only benign diagnostic provider changes |
| G5 TaskManager T1 | PASS | 13/13 files and 53,484 observer records exact in each mode |
| G6 TaskManager T2 | PASS | 38/38 files exact in each mode; all three BNV lifetime tables exact |
| G7 governed science | PASS | Phase-5B/5C/5D/coupled/5D1 and five hash baselines pass; fresh OLD and NEW bounded ADR-0017 match accepted historical authority exactly |
| G8 link closure | PASS | 67 link commands and 66 final executables inspected per mode; no vendored/Python/OpenMP inputs or dynamic_lookup; explicit GSL/zlib/Threads closure |
| G9 performance | PASS | Representative paired T1/T2 and long coupled-oracle timing show no material regression; whole-process ADR timing recorded separately |

The historical strict provider gate and hard performance stop are superseded
only by the owner's appended resume authorization. Numerical comparisons and
existing governed tolerances remain unchanged. The initial Release failures,
their classification/repairs and original logs are retained; no failed attempt
is presented as an original all-green run.

Nonblocking findings are the still-unimplemented spin diagnostics, source/build
rather than universal runtime coverage of legacy analysis export paths, the
Mac-specific package lock, and the separately deferred TaskManager concurrency
redesign. The inherited full-provenance cost in ADR-0017 remains visible.
No unresolved scientific regression or build/dependency blocker remains.

Disposition **B — MIGRATION COMPLETE WITH NONBLOCKING CAVEATS — READY FOR ONE
INDEPENDENT REVIEW**. The next and only next step is one combined independent
review of exact Zaki e263a6e, canonical CONFIND b0cbd510 and the complete
CompactStar migration candidate. Do not merge, separately integrate Zaki,
tag, start Linux/cluster qualification, or begin a new scientific campaign here.

## Commit and repository accounting

| Change | SHA |
|---|---|
| authorization | `8be1b6a083e1adac1dba7d64c75aa437533ab69a` |
| resume authority | `a7c4934fc6cea7e2b58b7d3437dcd910ab2e873c` |
| historical references | `d8d4604009c93a51b8f6657284d56b5df9e13659` |
| plotting removal | `02db887979ae333798543e9ebde2c75a935c41d1` |
| package python openmp migration | `a940e09aba3831803ad154cf77539be1c9cce54f` |
| explicit unavailable spin | `fa359772abd83e2b0ef676885aaed485bb4bbbc7` |
| taskmanager regressions | `d26e29b8e78b2df0101034103e1cb4b381af5a18` |
| ADR0019 | `164e214413d5a272b1792f815ab9d7d18aa62bf8` |
| protected install repair | `d75e8fffff7e6fbabdba0238272ce99d8c0449b6` |
| Release same mode references | `865c96eb0a55c5d5e4d424681bc076e2a336d549` |
| heat capacity fixture isolation | `cb2aa04b1ef3daf92963bb61f7fabe50880ed17c` |


Qualification source/code tip is `cb2aa04b1ef3daf92963bb61f7fabe50880ed17c`;
subsequent commits record evidence and documentation only. The exact final
candidate SHA is the documentation commit containing this report and is pinned
in the owner handoff. The final evidence commit is `293a396470f343138bf61f289df9ebebf0028898`.

CompactStar canonical master remains `812463ac9ed374f64ac9cadd500066ab723d3a6c`;
CONFIND remains `b0cbd510fd3fd0c772fa50499cd749287cb39e7b`;
Zaki's authorized candidate remains `e263a6e180c5c417198e7778bd21fc9c0a32dc33`.
External source trees are unchanged and clean. Historical dependencies,
ADR-0018, all existing governed baseline bytes, the 33 protected Phase-5D paths
and all baseline-recorded scientific production hashes are unchanged. The
original predeclaration body is unchanged; implementation results are appended.
No physics change beyond the described migration, no tolerance change,
no master update, no release tag, no cluster access and no Linux qualification.

Compact evidence and reproducible local recipes are under
`docs/validation/evidence/external-deps-migration/resume/`, authenticated by its
SHA256SUMS. Large build trees and raw characterization tables remain in the
external execution root. The migration branch is to be left clean and pushed
non-force; final owner handoff records the checked local/upstream/live tip.
