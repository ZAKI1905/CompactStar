# EKU Linux transition qualification — 2026-10-10

This is an unqualified Linux candidate. Debug and Release build successfully and
both pass the 13 affected Mac-readiness tests. Scientific qualification is partial;
the frozen Mac certificate cannot currently attach to the Linux-generated star.
No BNV source/control trajectory or campaign ran. The thermal-regression scope
correction below is explicitly retained.

| Milestone | Result |
|---|---|
| Linux build successful | YES, Debug and Release |
| Linux functional tests passed | PARTIAL: readiness 13/13 in each mode; broader results below |
| Linux numerical equivalence qualified | NO |
| Linux scientific campaign authorized | NO |

## Authenticated import

The transferred outer ZIP contained the original handoff ZIP. Its SHA-256 is
`da21bd7d745ca03a145176e632da15b2a4a6ad6b4ce8b95f888de138be65e737`, matching
both the owner's independently supplied digest and the adjacent external receipt.
The supplied verifier passed 55 payload files, 23 scientific identities, and all
three bundle digests. The exact supplied importer restored each filtered bundle;
refs, commit sets, supplied source hashes, omitted-file absence, sparse clean
status, and `git fsck --full` passed. No bundles were manually reconstructed.

Sources are isolated under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/sources`.
The original input trees and import evidence remain under `incoming` and
`bootstrap-20261010`. See `authentication.json` and archived verifier/import logs.

| Identity | Commit |
|---|---|
| CompactStar readiness | `5433fc2da71ba380264f490b87926d850ff9f375` |
| Candidate base / cluster preparation | `d8e7208b71e01168394e41f8f8f32b027b5d0a63` |
| Preserved canonical master | `9a5c1eca40758e4c39f9f80db7324e12ad23e505` |
| Zaki 2.0.1 production | `e263a6e180c5c417198e7778bd21fc9c0a32dc33` |
| CONFIND 2.0 production | `b0cbd510fd3fd0c772fa50499cd749287cb39e7b` |

All six frozen Phase-6 inputs retain authenticated bytes, including embedded Mac
paths. A separate adapter maps logical input keys to relocated paths. Original
canonical and historical refs remain intact. PWNSpall and unrelated workspaces
were not modified. The importer intentionally retains an offline promisor remote
for omitted historical binary blobs; this is not a complete binary-history clone.

## Native dependency and toolchain evidence

The login node is Rocky 9.7; compute node phyn002 is Rocky 9.4. Compute nodes lack
development headers and Tcl for environment modules. Modules were loaded on the
login node and exported through Slurm. All substantial builds and calculations
ran in Slurm allocations. Local administration and hashing used the login node.

The candidate uses shared GCC 14.1.0, CMake 3.29.4, Make 4.4.1, binutils 2.42,
C++17, `-ffp-contract=off`, Debug `-g`, and Release `-O3 -DNDEBUG`.
HistoricalMath translation units retain `-fno-builtin -fno-lto`. No fast-math,
native tuning, or LTO was enabled. Zaki uses GNU C++17 extensions; CONFIND and
CompactStar use their explicit C++17 configuration.

The isolated native sysroot stages installed-RPM-content-verified glibc
2.34-231.el9_7.10, its development headers, kernel headers
5.14.0-611.34.1.el9_7, and zlib 1.2.11-40.el9. GCC wrappers explicitly select this
sysroot and its dynamic loader/runtime. This runtime differs from the compute
node's system glibc 2.34-100.el9_4.2. This is an installed-RPM integrity check,
not a newly obtained independent package-signature attestation. RPM ownership
differences were retained in logs; content mismatches were not accepted.

Existing GSL 2.6 was unsuitable. The authenticated handoff GSL 2.7.1 archive
(`dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b`) was
built static with `-O2 -g -ffp-contract=off`. Its 59 upstream suites initially
reported 58 PASS, one pkg-config SKIP, zero FAIL. After making pkgconf available,
the packaging check passed separately; the initial skip remains recorded.

Zaki passed 12/12 tests in both modes. CONFIND initially passed 3/29 in each mode;
the installed-consumer test failed because its subprocess lacked GSL/zlib search
paths. With explicit prefixes it passed, giving 4/29 effective passes per mode.
The other 25 exact historical numerical comparisons per mode remain failures.
Within Linux, stellar serial and worker counts 1, 2, 3, 4, 6, 8 produced identical
tables in each build mode. This does not establish cross-platform equivalence.

Python 3.11.9 uses NumPy 2.2.6, SciPy 1.15.3, and mpmath 1.3.0. The initial
NumPy 1.26.4 environment lacked `numpy.trapezoid`; the first broad-suite attempt
was stopped, its failure log retained, and versions updated before the retry.
Both sets of wheel hashes and actual package versions are recorded.

Installed archives, headers, CMake configs, install manifests, compile commands,
tool executables, GSL/runtime/compiler closures, source identities, and original
test exit codes are recorded in the Linux lock and associated evidence. Native
artifacts remain outside Git under the isolated `toolchain` directory.

## Portability and lock implementation

Darwin retains `cmake/dependency-lock.json` byte-for-byte, SHA-256
`e1307a042180034ad028f2ff07dbc7fd80752494ff4b2ff61cc8ca0eec2089c5`.
Linux requires an explicit absolute `COMPACTSTAR_LINUX_DEPENDENCY_LOCK`, with
separate Debug/Release artifacts, exact source identities, compiler/flags,
GSL/zlib selection and hashed closure. Authentication runs before package CONFIG
files execute. Missing or mismatched selections fail closed. Six authentication unit tests and four real negative CMake configurations
(missing lock, wrong-mode package prefix, changed compiler flags, Mac lock on
Linux) passed. Native macOS execution was not repeated on this cluster.
Neither a lock nor
its `CANDIDATE_UNQUALIFIED` status confers scientific acceptance. Absolute prefixes
bind this installation; relocation/rebuilding requires reviewed reauthentication.

The candidate replaces C++20 `using enum` with C++17-qualified enum names, adds
direct standard-library includes, and replaces a libc++ internal clock type with
the public type corresponding to the existing clock calls. Linux harnesses receive
explicit compiler/dependency paths and static GSL selection. No scientific
equation, constant, historical baseline, or scientific threshold was changed.
The migration harness also uses the native temporary directory on Linux while
retaining `/private/tmp` on Darwin. Four original path failures are retained;
T1/replay pass after repair, while T2 and observation comparisons still fail.
The original Mac paths remain in the Darwin harness branches. ADR-0020 is proposed
for review, not ratified.

## Scientific inventory and outcomes

| Inventory | PASS | FAIL | Excluded from bounded scope |
|---|---:|---:|---:|
| Debug corrected bounded scope | 72 | 10 | 4 |
| Release corrected bounded scope | 72 | 7 | 4 |

Debug additionally completed two now-excluded scalar thermal tests (both PASS)
before the scope correction. Thus its observed inventory contains 74 PASS,
10 FAIL and two unexecuted deferred entries out of 86 registered. Release has
83 registered tests: 79 executed and four deferred. Initial packaging/Python/
temporary-directory failures remain in separate logs; repaired reruns are
explicitly attributed in the outcomes JSON. See `TEST_INVENTORY.tsv` for every
test and evidence source. The failure inventory is:

| Test | Debug | Release |
|---|---|---|
| `baryon_number_cmf` | FAIL | FAIL |
| `hartle_monopole_regression` | FAIL | FAIL |
| `phase5b_structural_response_regression` | FAIL | not registered |
| `phase5c_chemical_coefficient_regression` | FAIL | not registered |
| `phase5d_coupled_oracles` | FAIL | not registered |
| `phase5d_protected_manifest` | FAIL | FAIL |
| `relativistic_unit_boundary` | FAIL | FAIL |
| `rotochemical_trackr_pe` | FAIL | FAIL |
| `taskmanager_T2` | FAIL | FAIL |
| `taskmanager_observations` | FAIL | FAIL |

Independent Phase-5D response, independent oracles, component tolerances and
harness controls pass in both modes; the protected manifest and coupled-oracle
guards remain distinct failures. The breakpoint contract exercises both RKF45
and rk8pd. Static/analytic and CMF EOS/TOV/rotation/surface checks pass except
for the explicitly listed regressions and near-boundary checks.


Scope correction: the initial broad inventory included standard scalar thermal
regressions with a 1 Myr physical horizon. `passive_cooling_regression` and
`grid_convergence_cmf` completed before cancellation took effect. This exceeded
the conservative bounded-only plan; these results are retained, not hidden.
Further runs of those tests and `tov_surface_grid_cmf` were excluded. The static
Hartle test interrupted by cancellation was rerun in the resumed allocation.
No BNV source/control trajectory was executed.

The inventory intentionally defers `phase5d1_controlled_evolution_regression`,
which generates the entire 10 Gyr mathematical trajectory, beyond the declared
bounded plan. Its test and baseline remain unchanged. The 13 readiness tests cover
chemical production/exact/curved oracles, particle-number analytic behavior,
passive checkpoint output, direct BNV energy, thermal ledger, heat capacity,
cache contracts, component tolerances, breakpoint integration, and currentness.

Additional authenticated frozen-monitor execution passed and detected its M20
runtime-monitor condition. The relocated source, control, and event attempts all
refused before evolution with `UNRESOLVED characterization profile differs from
owning canonical star`. The guard compares the Linux-generated star's r, m, nu,
and nb against the frozen Mac profile exactly. No guard was weakened, certificate
rebound, or substitute input fabricated. Consequently, bounded source/control
trajectories, event uncertainty and trajectory R18/R20 behavior remain unqualified.
Passing ledger and breakpoint unit contracts is narrower evidence.

The relativistic unit-boundary test fails U11 in both modes: the synthetic
threshold selects layer 0 rather than expected layer 1. Its row-1 density
reconstructs as 9,999,999,999.9999981 g/cm^3 at a requested 1e10 threshold.
`FindTbIndex` compares stored geometric density against the forward-converted
threshold and selects the outermost inner-side layer. This is consistent with
a near-boundary rounding sensitivity in the extended-precision fixture; the
rule and exact expectation remain unchanged. U1/U2/U5 had passed before refusal.

The CMF baryon-number historical regression also fails its unchanged 1e-15
relative bound: observed sequence discrepancies reach 4.999e-12, with associated
structure-baseline differences. The baryon-number magnitude check passes.
The reported solver refusals at low central densities occur during the mass
search; the final four target stars were constructed. These errors are retained.

The Mac readiness evidence covers the 13 affected tests, not a fresh run of every
old regression at this source revision. Therefore the broader failures cannot
all be attributed solely to the operating system without matched Mac reruns.

The proton-electron EOS test fails PE-V6: the production neutron-onset density
is 7.356728903732837544e-9 fm^-3 versus the independent long-double expression
7.356728903732831693e-9 fm^-3, absolute difference 5.85e-24 fm^-3 and relative
7.95e-16, exceeding its unchanged 3e-16 threshold. The independent expression
passes its separate reviewed high-precision reference check. The x86_64 long
double significand has 64 bits; no assumption of Mac long-double behavior was
substituted. A separate diagnostic calls every unchanged PE-V1 through PE-V13
check individually: only V6 fails. This remains a scientific tolerance/arithmetic
review item, not an authorized production formula or threshold repair.

Historical protected-source pins also differ: three pre-existing accepted
amendment paths and this candidate's mechanical particle-number compiler adapter.
The old pin manifest remains untouched. Exact names/digests and test refusals are
archived; passing independent response tests does not override these guards.

## Numerical differences and assessment

`PREDECLARATION.md` fixed the criteria before numerical execution. Existing
analytic residual bounds remain authoritative for their named tests. Cross-OS
trajectory assessment requires both runs' justified uncertainties within the
existing component budgets; R18 <= 1e-10, R20 <= 2e-4 and reconstruction <= 5e-6
remain unchanged. No blanket bitwise equality or arbitrary epsilon was imposed
as a cross-platform accuracy criterion. Exact provenance identities still apply.

| Diagnostic | Observed Mac/Linux difference | Assessment |
|---|---|---|
| Zaki fixed historical table | 1,595 changed floats in 162,619 rows; max relative 9.15e-7 | Near-zero integral_3: absolute 5.55e-17 at about 6.067e-11; exact oracle fails |
| Zaki stress table | 1,886 changed floats; max relative 2.15e-13, max 1,184 ULP | No independent cross-platform acceptance bound supplied |
| CONFIND stellar table | 418 changed floats in 81,457 rows; max relative 4.14e-15, max 32 ULP | Max absolute 8 occurs at roughly 1e15 coordinates; exact oracle fails |
| Static profile r | max absolute 1.92e-13 km; 118 changed rows | Diagnostic only |
| Static profile m | max absolute 6.66e-16 geometric km; 4,398 changed rows | Diagnostic only |
| Static profile nu | max absolute 1.17e-15; all 20,259 rows changed | Exact certificate attachment fails |
| Static profile nb | max absolute 1.11e-16 fm^-3; max relative 1.09e-10 near zero | No certificate rebinding authority |

The freshly passed Debug coefficient production result was also compared
componentwise with the governed Mac artifact under the predeclared combined
uncertainty diagnostic. G/Q/Z maximum relative differences are approximately
2.66e-14 / 1.15e-14 / 1.19e-14; I/W are 5.79e-11 / 5.78e-11. Every component
fits inside the sum of the two artifacts' supplied propagated uncertainties;
the largest budget fraction is 1.62e-4 (I). This is diagnostic consistency, not
new governed-baseline acceptance, a continuum proof, or certificate rebinding.

The Debug migration observation table differs in 37 of 53,484 float records,
each by one ULP; its 13 final T1 artifacts and 1,510-record replay match exactly.
T2 differs in 22 of 38 recorded output hashes, including static BNV lifetime
outputs; no new tolerance was applied. These are static calculations, not source/
control evolution. The observation and T2 comparisons remain failed.

The Zaki Debug/Release historical producer tables match each other exactly, and
the dependency comparisons found no row-count or nonnumeric-field differences.
Two independently regenerated Linux Debug coefficient artifacts were also
byte-identical (SHA-256 `ab87361ba9663d31b2aeac982fb90d637cbeeb534d34605338a432906ed59920`).
These observations characterize drift; they do not establish numerical or
scientific equivalence. The static diagnostic reused the existing EOS/solve
implementation but consumed no frozen certificate. Raw results and per-column
absolute/relative/ULP summaries are retained separately from pass/fail tests.

## Slurm performance and execution planning

The available nodes have AMD EPYC 7352 CPUs, 48 physical cores / 96 logical CPUs
and 125,000 MB reported RAM. phyn001 was drained (`Kill task failed`); phyn002-004
were usable. `select/cons_tres` uses CR_CORE rather than memory-aware scheduling.
One requested CPU may allocate two logical CPUs belonging to one physical core.
Slurm accounting storage is disabled, and sstat returned unusable resource fields;
do not infer measured peak memory from requested reservations.

The final build allocation (13209, four CPUs, 12 GiB requested) completed in
2m58s, including incremental Debug work, Release build and readiness checks.
This is not a clean-build benchmark or a proof of byte-identical rebuilds.
Track-R characterization took 137/138 seconds (Debug/Release); production
chemical-response validation took 472/218 seconds. These are measured static/local
workloads, not estimates of ODE throughput. Library and breakpoint runtime smoke
checks also passed on phyn003 and phyn004, using the intended staged loader and
libraries. The bounded source/control/event refusals
took approximately 1.1 seconds and 20 MiB each; those are startup costs, not
successful ODE performance measurements.

For future authorized pilots, use one physical core per serial ODE and set
OMP/OpenBLAS/MKL threads to one. Preserve distinct outputs/manifests for SOURCE
and CONTROL at BASELINE, REFINED and ULTRA. A six-element array or six explicit
jobs is appropriate; reference/reconstruction jobs should depend on successful,
authenticated inputs and have separate limits. Begin with at most two concurrent
pilot roles after certificate qualification, measure wall time and peak RSS, then
choose concurrency from measured memory and shared-I/O use. Four CPUs and 8–12 GiB
are reasonable provisional build reservations based on successful allocations;
successful evolution measurements are missing, so campaign RAM/wall-time values
cannot yet be justified. No six-role campaign was submitted.

## Blockers, review and next action

1. Review the Linux candidate and its arithmetic/runtime provenance independently.
   Exact dependency oracles have not qualified Linux numerical equivalence;
   PE-V6 also exceeds its declared onset threshold.
2. Establish a scientifically justified Linux profile/certificate authority or a
   reviewed cross-platform identity rule with propagated uncertainty. The small
   observed differences alone cannot authorize consuming the Mac certificate.
3. Resolve protected-source provenance through its owning authority, retaining
   historical manifests and failures. Rerun the affected bounded tests afterward.
4. Obtain successful bounded source/control, event and reference/reconstruction
   evidence before recommending campaign resources or requesting a campaign.

Tested final source/harness commit: `de45f22ea811370052d1a8c690b2b7d1780dd7d6`.
Candidate branch: `candidate/linux-transition-20261010`. The final report commit
is identified by `git rev-parse candidate/linux-transition-20261010`; the exact SHA
is also returned with this report. Canonical master remains unchanged.

The candidate is committed locally. Automatic approval review rejected the
candidate-only push to `https://github.com/ZAKI1905/CompactStar.git` before it
executed, because this destination was not explicitly confirmed for potentially
private repository contents. The remote matches the repository documentation and
published canonical identities, but destination approval is still required before
retrying. No push or canonical integration is claimed.

Recommended next action: approve that exact candidate destination for review,
then conduct one independent scientific/numerical review focused on the Linux
profile identity and dependency-oracle drift. Do not authorize a long campaign.

## Evidence and reproduction

`evidence/index.json` maps retained evidence files to original paths and SHA-256
digests. Compressed files are byte-preserving gzip copies of textual logs/JSON/XML;
The old bootstrap `qualification-status.txt` predates the transferred archive
and is retained as historical evidence; this report supersedes its status.
Raw logs, native installations, generated scientific profiles and bounded-attempt
directories remain in `toolchain`. Failed build attempts, the first Python suite
attempt, historical numerical failures, and successful repairs are all retained.
The scripts and cluster-specific reproduction notes are under `tools/linux`.
