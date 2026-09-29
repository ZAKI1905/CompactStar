# External dependency migration — link-authority stop

**Disposition D — PACKAGE / LINK / IDENTITY AUTHORITY FAILED — RETURN TO OWNER.**

The identities pass. The Release provider/arithmetic clause of G4 fails in a
fresh TaskManager-bearing diagnostic link. There is no complete migration
candidate, no T1/T2 numerical comparison from this pass, and no demonstrated
stellar numerical difference. This is an implementation stop, not review-ready.

The owner authorized exact Zaki candidate e263a6e provisionally and one combined
review after implementation. Lack of a separate Zaki review is **not** the
blocker. The blocker is the preserved CompactStar G4 criterion, quoted below.

## Authenticated scope and completed work

Entry authentication and owner authorization are in the append-only
[predeclaration](PHASE6_EXTERNAL_DEPS_TASKMANAGER_MIGRATION_PREDECLARATION.md:1063).
`MIGRATION_AUTHORIZATION_SHA` is
`8be1b6a083e1adac1dba7d64c75aa437533ab69a`.
Canonical CompactStar local/origin/live master remain
`812463ac9ed374f64ac9cadd500066ab723d3a6c`; CONFIND local/origin/live master remain
`b0cbd510fd3fd0c772fa50499cd749287cb39e7b`. Exact Zaki candidate remains clean at
`e263a6e180c5c417198e7778bd21fc9c0a32dc33`.

Fresh canonical OLD CompactStar Debug and Release builds succeeded. The five
hash-only OLD Debug baseline artifacts were regenerated with their existing
producers and all five reproduce the committed bytes. Separate package builds
and installs used the exact clean dependency SHAs, AppleClang 21.0.0
(clang-2100.3.34.2), CMake 4.2.1 and historical `/opt/local` GSL 2.7.1.
Zaki passed its **12 registered tests** in each mode; CONFIND passed **29/29**
in each mode. Zaki's prior expanded 20/20 suite was not enabled in these builds
and is not a fresh result of this pass. Prefixes and archive/header/config hashes
are recorded in [packages.json](evidence/external-deps-migration/packages.json).

Registered branch ancestry matches the predeclaration's known divergence:
the superseded ADR-0018 experiment touches root CMake, and historical Phase-6
branches touch test registration. None was integrated or substituted.

Change classes: repository documentation/evidence only; external build and
read-only artifact diagnosis. No production, CMake, maintained test, baseline,
EOS/data, dependency source, or historical oracle bytes changed in any repository.
The sole disposable source patch removes the previously qualified 26 plotting
statements in TaskManager and DarkCore_Analysis, with exact patch SHA-256
`3e82a276d5425f68529ef7648c5a7900337fff84edeea685eab2fe5105884b52`.
Its 23 plot calls and three setup statements are separately manifested.

## Reproduced G4 failure

The original predeclaration G4 requires:

> provider class (CompactStar TU / Zaki archive / CONFIND archive) equals OLD's
> for every symbol containing FP arithmetic, class changes are permitted only
> for symbols with no FP arithmetic, and no chosen definition contains an
> FMA-class instruction where OLD's did not

Source: [predeclaration line 839](PHASE6_EXTERNAL_DEPS_TASKMANAGER_MIGRATION_PREDECLARATION.md:839).
The user also requires returning to the owner on an FP-bearing provider failure,
without post-result reclassification or tolerance.

The exact function is `Zaki::File::VecSaver::Export1D<Zaki::Math::Coord2D>`.

| Release evidence | OLD | NEW diagnostic |
|---|---|---|
| Selected provider | historical `libZaki.a(Math_Core.cpp.o)` | consumer `libCompactStar.a(TaskManager.cpp.o)` |
| Linkage from `nm -m` | weak external | non-external, was private external |
| Floating arithmetic | `ucvtf; fmul; fadd; fcvtzu` | `ucvtf; fmadd; fcvtzu` |
| Relevant address | multiply `0x100092368`, add `0x100092370` | fused operation `0x100020808` |

The arithmetic sizes a compressed-export buffer:
`unsigned long sizeDataCompressed = (sizeDataOriginal * 1.1) + 12`.
The historical and candidate public `VecSaver.hpp` are byte-identical, SHA-256
`8f6ef4e30716917ce043746dd3b9968ca6d76adc0674d7b7907d76ea193f3461`.
This is not a changed source expression; different provider selection chooses
differently compiled arithmetic. It is outside the text-export branch used by
T2. No compressed-output or stellar-output difference is claimed.

Zaki's prior qualification already disclosed this limited-closure provider
change as nonblocking for that dependency task. That record does not amend
CompactStar's stricter G4. Exact-candidate authorization does not supply an
explicit exemption from the simultaneously retained FP-provider criterion.
No exception has been inserted for this result.

The focused diagnostic built the actual TaskManager driver against canonical
OLD libCompactStar, then independently rebuilt the 12 CompactStar archive
members selected by that link against the fresh pinned packages. It preserved
consumer arithmetic flags, original archive member ordering, and the governed
one-thread driver. The only disposable source edits were the authenticated
two-file plotting deletion. These are real final executable link maps and
machine code, but **not a complete final migrated CompactStar build**.

Provider inventory: Debug has four changed provider records among 1,581 common
namespace-containing symbol records, none with FP arithmetic. Release has five
changed records among 784, including the one FP-bearing VecSaver function.
The counts include template instantiations, stubs, GOT entries and switch tables;
they are not directly interchangeable with the reconnaissance's 1,379 count.
The reproduction does not recover the anticipated three/no-FP result.

Focused Coord3D selection agrees with the expected pattern: OLD selects
historical `libConfind.a(Cont2D.cpp.o)`; NEW Debug selects pinned
`libCONFIND.a(Cont2D.cpp.o)`; NEW Release has no out-of-line comparators. This
does not substitute for the unperformed all-final-target/interposition audit.

Evidence: [provider report](evidence/external-deps-migration/provider-probe-audit.json),
[Release OLD full link map](evidence/external-deps-migration/Release-old.map.gz),
[Release NEW full link map](evidence/external-deps-migration/Release-new.map.gz),
[OLD disassembly](evidence/external-deps-migration/Release-old-selected-disassembly.txt),
[NEW disassembly](evidence/external-deps-migration/Release-new-selected-disassembly.txt),
and [bundle inventory/reproduction](evidence/external-deps-migration/README.md).

## Stop accounting and next action

Departure from stage ordering: the known provider risk was checked early,
after fresh package/OLD builds and before completing M1 or broader plotting
edits. It failed; no later implementation or scientific gate was attempted.
Full CompactStar suites, T1/T2 execution, G2/G3 integrated characterization,
Phase-5/6 migration regression, ADR-0019, package discovery, Python/OpenMP cleanup
and performance qualification therefore remain **NOT RUN / NOT IMPLEMENTED**.
No nonzero numerical difference is invented from this artifact-level failure.

Both historical archives remain preserved at their specified hashes. All three
canonical/candidate repository identities remain unchanged and clean; only this
migration branch carries documentation and compact text/compressed evidence.
No merge, tag, Linux qualification or cluster access occurred. No process or
numerical campaign remains running from this pass.

**Exact next action:** owner adjudication of the G4 scope for the specifically
identified compressed-buffer sizing function. If the owner chooses a narrow
exemption, record it explicitly before resuming; otherwise a separately authorized
provider-preservation solution is required. Do not silently treat the prior Zaki
caveat as a CompactStar waiver. After that decision, resume the incomplete
migration and its original exact scientific gates. The single independent review
of both candidates remains after complete implementation; it has not begun.

## Complete requested 101-field accounting

NOT RUN and NOT MEASURED never mean PASS or zero difference. All results below
are from this pass unless explicitly identified as historical.

| # | Required field | Result |
|---:|---|---|
| 1 | CompactStar canonical entry SHA | 812463ac9ed374f64ac9cadd500066ab723d3a6c |
| 2 | Migration branch/worktree | physics/external-deps-taskmanager-migration; /Users/keeper/Documents/CompactStar/worktrees/CompactStar-external-deps-taskmanager-migration |
| 3 | Migration authorization SHA | 8be1b6a083e1adac1dba7d64c75aa437533ab69a |
| 4 | OLD-reference commit SHA | NOT CREATED; fresh OLD builds use canonical 812463a; five baseline reproductions retained |
| 5 | Plotting-removal SHA | NOT CREATED; qualified two-file patch used only in disposable diagnostic export |
| 6 | External-package migration SHA | NOT CREATED |
| 7 | Python/OpenMP cleanup SHA | NOT CREATED |
| 8 | T1/T2 qualification SHA | NOT CREATED |
| 9 | ADR-0019 SHA | NOT CREATED |
| 10 | Final documentation SHA | STOP_RECORD_SHA: commit containing this report; exact resolved SHA supplied in final handoff |
| 11 | Final migration candidate SHA | NONE; stop-record branch tip is not a completed migration candidate |
| 12 | Zaki candidate SHA | e263a6e180c5c417198e7778bd21fc9c0a32dc33 |
| 13 | CONFIND canonical SHA | b0cbd510fd3fd0c772fa50499cd749287cb39e7b |
| 14 | OLD Zaki archive SHA-256 | 3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef |
| 15 | OLD CONFIND archive SHA-256 | 09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2 |
| 16 | Zaki Debug package archive SHA-256 | f30e68b261504e1ca3b9af905dc98a0ce70fa25f04c3c19ab81e3b2629fb074f |
| 17 | Zaki Release package archive SHA-256 | 3b7f375ff789cd2c0edf58f97274c06fa047e4f477bf327214f9406461baad9d |
| 18 | CONFIND Debug package archive SHA-256 | 1f36fc6ef6964d4342a025539a5a26c6c55faf04f70d1fe9319c0a4e70cc32c7 |
| 19 | CONFIND Release package archive SHA-256 | 4cee9b8361fec75023d74b8a68bf7efae54f69cf1464e478173bcfb5d38b524d |
| 20 | Package header manifest SHA-256 | Zaki both modes: 51931f9eac6823cc21609741c7f60485d0dbfab8d815d0a541e897cbcc2528a2; CONFIND both modes: 0e02d55747cfa7bce70fdc9c660a4b95ec540dd4ce0c9bd94fbb68bccd5fd99b |
| 21 | Fail-closed package discovery | NOT IMPLEMENTED in CompactStar; explicit source/build/prefix identity recorded for diagnostic packages |
| 22 | Active vendored Zaki usage | YES, repository CMake remains unchanged; migration incomplete |
| 23 | Active vendored CONFIND usage | YES, repository CMake remains unchanged; migration incomplete |
| 24 | Historical artifacts preserved | YES; authenticated hashes unchanged |
| 25 | Total active plot calls before | 75 per frozen predeclaration inventory; source untouched |
| 26 | Total active plot calls after | 75 in repository; source untouched |
| 27 | Plotting removal ledger count | Complete ledger NOT CREATED; external two-file diagnostic manifest has 23 calls + 3 setup statements |
| 28 | Plot deletions with numerical side effects | 0 repository deletions; diagnostic patch retains prior 0-side-effect classification; no new G2 runtime proof |
| 29 | Replacement numerical exports added | 0 |
| 30 | Direct Python C++ uses remaining | 0 in compiled sources; two unbuilt matplotlibcpp examples remain |
| 31 | Direct NumPy C++ uses remaining | 0 in compiled sources |
| 32 | Python linked into numerical executables | YES in unchanged OLD/repository build; absent from diagnostic NEW link only |
| 33 | Python Interpreter retained for tests | YES |
| 34 | dynamic_lookup present | YES in unchanged OLD/repository link; absent from diagnostic NEW link only |
| 35 | zlib closure correct | Diagnostic NEW links explicit -lz and packages export ZLIB::ZLIB; final migration G8 NOT QUALIFIED |
| 36 | GSL identity | /opt/local GSL 2.7.1; libgsl.27.dylib SHA256 9667d3f90d26c5cfb5229472adf455c0bf6406f946764ad165ec5d4b2121543b; gslcblas hash in stop-state.json |
| 37 | OpenMP active source use | 0 found in CompactStar/main/tests source scan |
| 38 | OpenMP linked after migration | Migration incomplete; OpenMP retained in unchanged OLD and staged diagnostic |
| 39 | TaskManager threading changed | NO |
| 40 | TaskManager governed thread count | 1 in preserved test driver; no T1/T2 run in this pass |
| 41 | OLD Debug full suite result | NOT RUN in this pass; 76/76 is prior reconnaissance, not fresh authority |
| 42 | NEW Debug full suite result | NOT RUN |
| 43 | OLD Release relevant-suite result | NOT RUN; build succeeded |
| 44 | NEW Release relevant-suite result | NOT RUN |
| 45 | T1 Debug differences | NOT MEASURED |
| 46 | T1 Release differences | NOT MEASURED |
| 47 | FindCriticalCurve raw-order differences | NOT MEASURED |
| 48 | FindCriticalCurve argmax differences | NOT MEASURED |
| 49 | mass_curve[0] Debug difference | NOT MEASURED |
| 50 | mass_curve[0] Release difference | NOT MEASURED |
| 51 | Intersection differences | NOT MEASURED |
| 52 | GetIdx differences | NOT MEASURED |
| 53 | Bisect differences | NOT MEASURED |
| 54 | T2 output count | 38 required; 0 captured in this pass |
| 55 | T2 Debug differing files | NOT MEASURED |
| 56 | T2 Release differing files | NOT MEASURED |
| 57 | Neutron BNV lifetime difference | NOT MEASURED |
| 58 | Lambda BNV lifetime difference | NOT MEASURED |
| 59 | Sigma- BNV lifetime difference | NOT MEASURED |
| 60 | Coord3D final provider | Focused NEW Debug diagnostic: libCONFIND.a(Cont2D.cpp.o); Release no out-of-line definitions; full final-build audit NOT RUN |
| 61 | Coord3D interposition gate | Focused provider pattern matches; all-target/loaded-image gate NOT COMPLETED |
| 62 | External symbol provider differences | Focused Debug 4/1581 common records; Release 5/784 common records |
| 63 | FP-bearing provider differences | Debug 0; Release 1, VecSaver::Export1D<Coord2D>; G4 FAIL |
| 64 | G0 result | PARTIAL: entry/package identities authenticated; fail-closed consumer discovery NOT IMPLEMENTED |
| 65 | G1 result | NOT RUN |
| 66 | G2 result | NOT RUN; diagnostic mechanical patch recorded |
| 67 | G3 result | Integrated differential NOT RUN; fresh Zaki 12/12 per mode, prior expanded 20/20 not rerun |
| 68 | G4 result | FAIL at early TaskManager Release provider gate |
| 69 | G5 result | NOT RUN |
| 70 | G6 result | NOT RUN |
| 71 | G7 result | PARTIAL: five OLD Debug hash-only artifacts exact; NEW and Phase-5/6 matrix NOT RUN |
| 72 | G8 result | NOT COMPLETED; diagnostic explicit links retained |
| 73 | G9 result | NOT RUN |
| 74 | Phase-5B result | Migration regression NOT RUN |
| 75 | Phase-5C result | Migration regression NOT RUN |
| 76 | Phase-5D result | Migration regression NOT RUN |
| 77 | Phase-5D coupled-oracle result | NOT RUN |
| 78 | Phase-5D1 result | NOT RUN |
| 79 | Phase-6 / ADR-0017 result | NOT RUN; accepted bounded historical status unchanged |
| 80 | Five hash-only baselines result | OLD Debug 5/5 exact committed bytes; NEW NOT RUN |
| 81 | T1 OLD wall time | NOT MEASURED |
| 82 | T1 NEW wall time | NOT MEASURED |
| 83 | T2 OLD Debug wall time | NOT MEASURED |
| 84 | T2 NEW Debug wall time | NOT MEASURED |
| 85 | T2 OLD Release wall time | NOT MEASURED |
| 86 | T2 NEW Release wall time | NOT MEASURED |
| 87 | Performance <=1.25x | NOT QUALIFIED |
| 88 | ADR-0019 status | NOT CREATED |
| 89 | No-post-hoc tolerance respected | YES; no numerical tolerance, baseline or FP-provider exemption added for this failure |
| 90 | CompactStar physics changed beyond migration | NO; all repository production bytes unchanged |
| 91 | Zaki source modified in this task | NO |
| 92 | CONFIND source modified | NO |
| 93 | Canonical CompactStar master modified | NO |
| 94 | Linux qualified | NO |
| 95 | Cluster accessed | NO |
| 96 | Release tag created | NO |
| 97 | git diff --check | PASS before report commit; cached/unstaged checks repeated at finalization |
| 98 | Blockers | Release FP-bearing VecSaver provider/arithmetic change violates G4; implementation stopped |
| 99 | Nonblocking findings | Exact source identities and both package builds pass; OLD Debug five baselines exact; no stellar drift claimed; prior Zaki expanded suite differs from fresh 12-test inventory |
| 100 | Disposition | D — PACKAGE / LINK / IDENTITY AUTHORITY FAILED — RETURN TO OWNER |
| 101 | Exact recommended next action | Owner adjudication of the named compressed-buffer sizing G4 scope; no automatic exemption or repair. Then resume implementation under explicit authority before the single combined review. |
