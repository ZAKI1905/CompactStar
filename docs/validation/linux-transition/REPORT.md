# CompactStar EKU Linux qualification — 2026-10-10

**Candidate only; numerical equivalence remains unqualified.** No BNV
source/control trajectory or campaign ran. Canonical master and Mac baselines are
unchanged. See [detailed results](DETAILS.md), [predeclared criteria](PREDECLARATION.md),
[authentication](authentication.json), and [evidence index](evidence/index.json).

| Milestone | Result |
|---|---|
| Linux build successful | **YES**, Debug and Release |
| Linux functional tests passed | **PARTIAL**: readiness 13/13 in each mode; broader results below |
| Linux numerical equivalence qualified | **NO** |
| Linux scientific campaign authorized | **NO** |

## Import and implementation

The supplied ZIP matches the owner's SHA-256 and external receipt:
`da21bd7d745ca03a145176e632da15b2a4a6ad6b4ce8b95f888de138be65e737`.
The exact supplied verifier/importer passed 55 payload files, 23 scientific
identities, three bundles, original refs/history sets, source hashes and fsck.
All six frozen Phase-6 inputs retain their authenticated bytes. The isolated
workspace is `/mnt/sdd/zaki/CompactStar/linux-transition-20261010`.

The candidate descends from cluster preparation `d8e7208b71e01168394e41f8f8f32b027b5d0a63`
and readiness `5433fc2da71ba380264f490b87926d850ff9f375`. Zaki production remains
`e263a6e180c5c417198e7778bd21fc9c0a32dc33`; CONFIND remains
`b0cbd510fd3fd0c772fa50499cd749287cb39e7b`. Master remains `9a5c1eca40758e4c39f9f80db7324e12ad23e505`.

A separate explicit Linux lock authenticates native Debug/Release packages,
compiler/flags and selected dependency/runtime closures before loading package
configs. The accepted Mac lock is byte-identical. Six authentication unit tests
and four real negative CMake configurations verify refusal behavior. Portability
changes qualify C++17 enum names, add required standard headers, replace a libc++
internal clock type, and forward native compiler/library paths in test harnesses.
Scientific equations, constants, tolerances and historical failure records remain
unchanged. PWNSpall and unrelated workspaces were not modified.

## Builds and tests

Native toolchain: GCC 14.1.0, CMake 3.29.4, C++17, `-ffp-contract=off`, isolated
glibc 2.34-231 / zlib 1.2.11 sysroot, authenticated static GSL 2.7.1. Historical
`-fno-builtin -fno-lto` protections remain. Python 3.11.9 uses NumPy 2.2.6,
SciPy 1.15.3 and mpmath 1.3.0. Archive/install/wheel hashes and flags are retained.

| Test group | Outcome |
|---|---|
| GSL | 58/59 initial PASS, one packaging SKIP; packaging rerun PASS |
| Zaki | 12/12 PASS in each mode |
| CONFIND | 4/29 effective PASS in each mode; 25 exact numerical failures remain |
| CompactStar readiness | 13/13 PASS in Debug and Release |
| Broad Debug / Release | Debug **72 PASS / 10 FAIL**; Release **72 PASS / 7 FAIL** in corrected bounded scope |
| Authenticated frozen monitor | PASS |
| Bounded source/control/events | All refused before evolution at exact profile identity guard |

Scope correction: two standard Debug scalar thermal regressions with a 1 Myr
horizon completed before their horizon was identified and cancellation took
effect. These exceeded the conservative bounded-only plan. Their results remain
recorded; further long-horizon thermal runs and the entire 10 Gyr
controlled-evolution regression were excluded. The interrupted static Hartle
check was rerun. The first Python suite attempt and packaging/build failures are preserved
alongside repaired reruns. No failed numerical test was relabeled as passing.

## Numerical assessment and remaining blockers

The passed Debug coefficient production result is byte-repeatable on Linux.
Its G/Q/Z/I/W differences from Mac fit inside the predeclared sum of supplied
propagated uncertainties; the largest ratio is 1.62e-4. The unchanged exact Mac
coefficient regression still fails. Zaki and CONFIND historical tables also have
nonzero differences without an independent cross-platform acceptance authority.

The proton-electron neutron-onset check also fails: relative error 7.95e-16
exceeds its unchanged 3e-16 threshold. Independent diagnostic execution isolates
this to PE-V6; the other twelve PE checks pass.

U11 also selects a different synthetic density layer at a near-roundoff boundary.
The CMF baryon-number regression reaches 5.00e-12 relative drift against its
unchanged 1e-15 bound. The supplied Mac 13-test record does not establish that
every older regression passes at the readiness revision.

The regenerated star differs from the frozen Mac profile by at most 1.92e-13 km
in radius and 1.17e-15 in metric nu. These small differences do not authorize
certificate rebinding. Protected-source pin mismatches also remain. Therefore
bounded source/control trajectories, event uncertainty, and trajectory R18/R20
behavior are not qualified; ledger and breakpoint unit passes are narrower evidence.

## Execution planning and next action

All substantial work used Slurm. The final four-CPU, 12-GiB build allocation took
2m58s, including incremental work; this is not a clean-build benchmark. Failed
one-second evolution startups do not measure ODE performance. Accounting did not
provide reliable peak memory. phyn001 was drained; other nodes were usable.

Use one physical core per serial ODE. After resolving certificate authority,
benchmark at most two concurrent bounded pilots, then size six separate
SOURCE/CONTROL × BASELINE/REFINED/ULTRA jobs from measured wall time, memory and
I/O. Reference/reconstruction jobs need separate authenticated dependencies and
limits. Campaign reservations are not yet justified.

Branch: **`candidate/linux-transition-20261010`**. This report's exact candidate
SHA is returned in the task response and obtained with `git rev-parse HEAD`.
The candidate is committed locally; **push is pending explicit destination
approval** for `https://github.com/ZAKI1905/CompactStar.git`. Automatic approval
review rejected the push before execution because potentially private repository
contents would be exported to an unconfirmed destination.

Recommended next action: approve that destination for candidate review, then
obtain one independent scientific/numerical review of dependency drift,
profile/certificate identity, the PE-V6 onset threshold, and protected provenance before further bounded
qualification. No canonical integration or long campaign is authorized.
