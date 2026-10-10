# Phase-6 cluster transfer preparation — Mac, 2026-10-10

**HANDOFF PREPARED; MAC REGRESSION GAP RECOVERED; LINUX UNQUALIFIED.**
This is an additive operational report. It does not alter the historical readiness
report, numerical acceptance limits, or failed BA12/BA12R results. The full BNV
campaign remains unauthorized (`PHASE6_BNV_CAMPAIGN_READINESS.md:1`,
`PHASE6_CAMPAIGN_ACCURACY_PREDECLARATION.md:1`).

## Scope, authority, source identities

Change classes: documentation and engineering transfer tooling, plus generated
validation evidence; no scientific-semantic, numerical-method, or build-resolution
implementation changes. Authority is the owner's Mac transfer request,
`GOVERNANCE.md:36`, accepted ADR-0018/0019, and
`COMPACTSTAR_EXTERNAL_DEPENDENCY_MIGRATION_ACCEPTANCE.md:3`.
All additions are under `docs/validation` on `docs/phase6-cluster-preparation`,
created in `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6-cluster-preparation`
from the named readiness candidate. The source audit is
`cluster-handoff/evidence/source-build-tests.json:1`.

| Identity | Exact commit | Role |
|---|---|---|
| CompactStar master | `9a5c1eca40758e4c39f9f80db7324e12ad23e505` | Canonical, unchanged |
| Breakpoint candidate | `0c3f9269c884161c9d89163660a47b141a93ace6` | Ratified bounded numerical architecture |
| Readiness candidate | `5433fc2da71ba380264f490b87926d850ff9f375` | Preparation base; scientific gaps remain |
| Zaki production | `e263a6e180c5c417198e7778bd21fc9c0a32dc33` | Exact source to build on Linux |
| Zaki accepted master | `7aca56cfdf3875cf2a8c0acea4acdac44f69c9f5` | One acceptance-document commit after production |
| CONFIND production/master | `b0cbd510fd3fd0c772fa50499cd749287cb39e7b` | Exact source to build on Linux |

Canonical, readiness, Zaki and CONFIND checkouts were clean. Their registered
worktrees, branch names, HEADs and readiness ancestry are recorded in the audit.
Master is an ancestor of the readiness candidate; ancestry does not ratify all
intervening commits. The earlier 24-commit integration boundary remains as recorded
in `PHASE6_BREAKPOINT_RATIFICATION.md:1`. No canonical ref was changed.

Production, dependency/CMake resolution, all tests, Cstar inputs and baseline bytes
remain equal to readiness. The production numerical diff from the accepted
breakpoint SHA is empty. All 33 protected paths match accepted bytes: 30 retain
historical Phase-5D pins; three are already accepted amendment deltas. All 11
baseline files are hashed separately (`cluster-handoff/evidence/source-build-tests.json:1`).

## Recovered Mac environment and regression result

Scriptorium is mounted as APFS volume UUID
`99EE5D09-F272-401A-95CB-0D19755A3C24`. Original symlink paths resolve into the
previous offload tree. Independent hashing and the unchanged production verifier
both authenticated **four of four** original package prefixes: Zaki/CONFIND,
Debug/Release. Archive, header-manifest and CMake-config digests match the existing
lock; source SHA, mode and clean provenance match. No replacement package was
built (`cluster-handoff/evidence/mac-recovery.json:1`,
`original-verifier-Debug.json:1`, `original-verifier-Release.json:1` in that directory).

Mac lock SHA-256 remains
`e1307a042180034ad028f2ff07dbc7fd80752494ff4b2ff61cc8ca0eec2089c5`.
A fresh CompactStar Debug build used the original Debug prefixes, AppleClang
21.0.0 (`clang-2100.3.34.2`), CMake 4.2.1, Python 3.12.10 and GSL 2.7.1 from
`/opt/local`. Compiler/CMake strings match the original package provenance; GSL
dylib and ODE-header hashes match the prior readiness audit. New binary, cache,
compile-command, linkage and Python-module identities are recorded explicitly;
no claim of a byte-identical CompactStar rebuild is made.
Evidence: `cluster-handoff/evidence/source-build-tests.json:1`,
`mac-configure.log:1`, `mac-build.log:1`; reproduction recipe:
`cluster-handoff/run_mac_tests.sh:1`.

**13/13 affected regressions PASS, zero failures/skips, serial execution, 15.93 s.**
This is a Debug result. Release packages were authenticated; a Release test rerun
was not part of this 13-test recovery. Optional external-EOS tests were not
registered in this targeted build. No numerical window or campaign was run.
The exact test list and raw CTest/JUnit evidence are in
`cluster-handoff/evidence/mac-tests.log:1`, `mac-tests.xml:1`, and
`source-build-tests.json:1`.

## Scientific input recovery and transfer contents

All six original Phase-6 input identities now PASS historical hashes. The formerly
missing frozen BNV certificate is
`9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab`;
coefficients are
`3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6`.
The certificate, coefficient and Cstar/thermal bytes were copied unchanged.
Recovery closes availability/authentication gaps, not the outstanding global
accuracy or tail-stability qualification (`cluster-handoff/evidence/mac-recovery.json:1`).

The transfer also includes both CompOSE sets (14 files), authenticated against
`PASSIVE_COOLING_BASELINE.md:268` and `TOV_REFERENCE.md:124`, and three obs32 input
files after verifying their original review manifest. The transfer's
`scientific/input-map.json:1` records 23 logical identities and relocation paths;
`cluster-handoff/copy_inputs.py:1` is the reproducible copier. The larger CompOSE
thermal table is relevant to later data-dependent regressions; unrelated stored
results and disposable builds are excluded.

GSL 2.7.1 source was downloaded from GNU over HTTPS and independently matched to
the original MacPorts source tarball and installed Portfile checksum:
`dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b`
(7,509,209 bytes). The original Portfile is included. No GPG signature verification
or complete reproduction of historical MacPorts patches/build flags is claimed
(`cluster-handoff/evidence/gsl-source.json:1`).

Each source bundle retains the exact original commit graph and tree objects,
selected refs, and source blobs; no history rewriting occurred. The v3 filtered
bundle design excludes Mac binary objects, generated Xcode/Doxygen output, old
results, credentials and duplicate archives. All active source/CMake/test files
are supplied. Omitted blobs and tip paths are explicit in bundle manifests.
These are source handoffs, not full binary-history backups; ordinary `git clone`
is not the import path. Use the supplied verified importer. A trial import checked
refs, full commit sets, supplied file SHA-256, clean sparse status and `git fsck`.
Implementation: `cluster-handoff/package_sources.py:1`, `import_sources.py:1`;
format authority: <https://git-scm.com/docs/gitformat-bundle>.

The final ZIP's `SHA256SUMS`, source manifests, unavailable-input manifest and
`verify_transfer.py` establish transport integrity. Its adjacent `.sha256` receipt
records the ZIP digest (outside the archive to avoid self-reference). The final
payload is made from committed source. See the external transfer receipt for its
exact preparation commit and final archive hash.

## Linux disposition and next action

The current independent requirements inventory and explicit Linux bootstrap
instructions are `cluster-handoff/README.md:1`. Remaining work is:

- Native Linux compiler/CMake/GSL/zlib/Python/toolchain authentication and build.
- A separately reviewable, fail-closed Linux lock selector; keep the Mac lock
  unchanged. Current verifier hardcodes it (`cmake/verify_dependencies.py:29`).
- Resolve the C++20 `using enum` construct under the C++17 contract
  (`CompactStar/Physics/State/Tags.hpp:72`).
- Relocate Mac-bound harness paths/library suffixes and separate platform metadata
  comparisons without weakening existing Mac oracles.
- Reconcile the historical cluster plan with ADR-0019, breakpoint Amendment 1 and
  the campaign-global accuracy predeclaration; predeclare tests and stop rules.

The older preflight describes a superseded dependency graph in places. Current
accepted targets require GSL, zlib and Threads; C++ Python/NumPy, ROOT and OpenMP
are not active requirements. Historical cluster modules remain unverified until
a campus-side inventory. Linux wheels/native dependencies and toolchain binaries
are deliberately not substituted with Mac packages; their status is explicit in
`unavailable-required-inputs.json` in the transfer.

| Readiness item | Status | Scope |
|---|---|---|
| Named source identities and ancestry | PASS | Local Mac refs and original objects |
| Original four Mac packages | PASS | Archive/header/config/provenance authentication |
| Frozen Phase-6 inputs | PASS | All six identities; exact recovered bytes |
| Targeted Mac regressions | PASS | Debug 13/13; Release rerun DEFERRED |
| Manual source/input handoff | PASS after final transport verification | Filtered history, source and scientific input payload |
| Linux dependency artifacts/platform lock | UNQUALIFIED | No Linux build or qualification performed |
| Linux C++17/harness portability | UNQUALIFIED | Specific blockers preserved above |
| Source/control BASELINE, REFINED, ULTRA global accuracy | UNQUALIFIED | No status promotion from artifact recovery |
| Tail stability policy | UNQUALIFIED | Readiness finding remains open |
| Historical BA12/BA12R | FAIL | Unchanged historical outcomes |
| Full BNV campaign / canonical integration | DEFERRED | Not authorized; no launch or merge |

**Recommended next action:** transfer and verify the ZIP on the campus laptop,
import the source bundles, then prepare the actual Linux host/toolchain inventory
and narrowly scoped platform/harness candidate work. Numerical platform
qualification requires its own authorization and evidence. Do not advance to a
full campaign or infer Linux scientific equivalence from successful compilation.
