# Phase-6A-1 controlled-BNV recovery implementation (R0-R3)

## 1. Status, authority, and scope

**Status (2026-09-19): R0-R3 RECOVERY PASS / RECOVERY BRANCH ONLY / NOT
MERGED / NO PHASE-6 CANDIDATE ARTIFACT.**

The canonical recovery entry is
`bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`. The historical failed
implementation `e56e6e50040dbcd9dcbee1acecf58843f3dddf1c` was used only as a
tree-level source of evidence; it is not an ancestor of this recovery branch.
The recovery branch is `physics/phase6a1-bnv-recovery-r0-r3` in
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-bnv-recovery-r0-r3`.

ADR-0016 is accepted and human-ratified. It assigns the adapter to the Phase-6
BNV module while retaining `Analysis::EquilibriumSequenceNumberDerivative` as
the Phase-5B derivative authority (`docs/adr/ADR-0016-phase6-bnv-tangent-adapter-ownership.md:22-42`).
The governed mathematics remains
`t = (partial N_y^eq / partial B)_Omega`, with `B_B = B_n + B_p` and
`t_i = B_i/B_B` at zero spin; the reduced `B_n+B_e+B_mu` is only a closure
check (`docs/adr/ADR-0016-phase6-bnv-tangent-adapter-ownership.md:84-118`).
No Phase-5D producer, comparator, provenance rule, baseline, or existing
Phase-5 source byte was changed (`docs/adr/ADR-0016-phase6-bnv-tangent-adapter-ownership.md:120-136`).

This record stops after R3. It does not claim BA12R, ULTRA, trajectory BA17,
BA11, BA16, candidate production, a full-suite result, a physical BNV model,
or canonical integration. That ordering is required by the accepted recovery
plan (`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:462-478`).

## 2. R0 ownership result

The recovered owner is exactly:

```text
header: CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp
source: CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp
type:   CompactStar::Physics::BNV::EquilibriumBaryonTangent
```

The forbidden historical paths
`CompactStar/Analysis/EquilibriumBaryonTangent.hpp` and
`CompactStar/Analysis/src/EquilibriumBaryonTangent.cpp` are absent. Repository
search found no stale production/test include or type reference. The adapter
performs no solve and computes no independent derivative; it consumes the
unchanged Phase-5B sequence derivative as required by ADR-0016
(`docs/adr/ADR-0016-phase6-bnv-tangent-adapter-ownership.md:65-78`).

The failed implementation was not merged, cherry-picked, rebased, or made an
ancestor. Files were reconstructed directly from `git show e56e6e5:<path>`.
All non-identical bytes below are limited to the ADR-0016 filesystem path,
namespace/type ownership, dependent include/type spellings, and the BNV-local
build registration. Diff inspection found no other production or test change.

## 3. Recovered production/build manifest

`historical path` names the source-evidence blob at the failed implementation
SHA. The tangent rows deliberately map the new BNV path to the old Analysis
path. There are 18 production files, plus two build-registration files.

| recovered path | historical path | historical SHA-256 | recovered SHA-256 | classification |
|---|---|---|---|---|
| `CompactStar/Physics/BNV/BnvDiagnostics.hpp` | same | `3fdda08e04e51bf1edb4a1b8c976138a9039262bf95aec2e9d74ab929d0808b8` | `3fdda08e04e51bf1edb4a1b8c976138a9039262bf95aec2e9d74ab929d0808b8` | IDENTICAL |
| `CompactStar/Physics/BNV/CMakeLists.txt` | same | `74a26a7c51424331c2d50443f4bbc30f08bc03b2462ba06ac2242d866be20a46` | `fe593d62348743f7715325a752ec8ee8214630235d32c60556685cd9e033a069` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp` | same | `b668b90029f80a47563ccabda01aa0dd891a449db7784923e0072aa99ca7d32c` | `b668b90029f80a47563ccabda01aa0dd891a449db7784923e0072aa99ca7d32c` | IDENTICAL |
| `CompactStar/Physics/BNV/DirectEnergyLedger.hpp` | same | `356e8cd00f9862afc43855237147484e8fc8dc4d3febaf10ff4bf0b143b7bc61` | `0a20dd65c0a15d5c614e0efa903abce6782143dbb8c16f439fc6a2293590d0cc` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/EquilibriumBaryonTangent.hpp` | `CompactStar/Analysis/EquilibriumBaryonTangent.hpp` | `7218b1b337037ed11697dca7c321f595be201679fa8683a21c03bcd5c2931771` | `202ca840e48782398424eda76aea6c563a7dc5d119e465d00c3a88e86d00eced` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/FrozenBnvValidity.hpp` | same | `d74d1b6e2537077adeebc01a35f0fd7cf22ac92ad7069b7f65da270398bda0c4` | `d74d1b6e2537077adeebc01a35f0fd7cf22ac92ad7069b7f65da270398bda0c4` | IDENTICAL |
| `CompactStar/Physics/BNV/FrozenControlledBnvRunContext.hpp` | same | `ab730e5567af2d4a3e047ce53fdb8c1a51db1f39c0887f8975e4a6995d6edeb7` | `d9644edcc315597c00cf0923745ae08eabfcf6c520f2e693bb6f6dd910be151c` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/MovingReferenceSource.hpp` | same | `1871988896dce8564d3d9606bea729b18c7a945bc5c3338b946ace5778630a46` | `ae754809660cd82aec256d7d249c55e31ce29d42566045f499aa1aebcce52ef0` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/OrdinaryMatterSource.hpp` | same | `58a4dc71e59aac30e30a3dea7ef8eccd0f544c9b27aadf517a8c8361ecd6cb1e` | `58a4dc71e59aac30e30a3dea7ef8eccd0f544c9b27aadf517a8c8361ecd6cb1e` | IDENTICAL |
| `CompactStar/Physics/BNV/ProductFate.hpp` | same | `20fe2e12e225764c4ef03ac855925dca1b80ddcffcfa9594af8de8940331875a` | `20fe2e12e225764c4ef03ac855925dca1b80ddcffcfa9594af8de8940331875a` | IDENTICAL |
| `CompactStar/Physics/BNV/StaticZeroSpinHistory.hpp` | same | `4dddab4cd66df82476fbd3a64ac78d6f626448c6bf51f739f1d31d28f0b55f2d` | `4dddab4cd66df82476fbd3a64ac78d6f626448c6bf51f739f1d31d28f0b55f2d` | IDENTICAL |
| `CompactStar/Physics/BNV/src/ControlledBnvSecularDriver.cpp` | same | `2abe4754badbcebf33201008dc2a364dfd3e1bd8f2ba5eac9890ace92c637abc` | `2abe4754badbcebf33201008dc2a364dfd3e1bd8f2ba5eac9890ace92c637abc` | IDENTICAL |
| `CompactStar/Physics/BNV/src/DirectEnergyLedger.cpp` | same | `390a30f2219f54b86b2f142b4d00de5843d576b73f19deb0e059c0a841aa77d7` | `4edbc2c7d692e2f73b2bf55c7c95f8de061d2c98ba4643ab90330918171734af` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/src/EquilibriumBaryonTangent.cpp` | `CompactStar/Analysis/src/EquilibriumBaryonTangent.cpp` | `a26b6e45f6c3e42a49e42d2d84293b9cf3bbc5657ab114ee274940772aa2d265` | `388b7576617312704ffb9cd312129c95cfed1b4d9512a2dbf90e72ab25de8f1b` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/src/FrozenBnvValidity.cpp` | same | `cd31381eea59d9cc6bedf5e099d547850f30647a756b20d68d444acf559f6654` | `cd31381eea59d9cc6bedf5e099d547850f30647a756b20d68d444acf559f6654` | IDENTICAL |
| `CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp` | same | `70bad4ac67cdf29c11777a94aa160aefbe0ebb73ed530b44be958727b9b2e658` | `56b63e9d8e8e3b0349191404843edf2d82856532e87838d09b1939e24b4c5468` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/src/MovingReferenceSource.cpp` | same | `445a1c0bf3655e547bfac3b8494ec11643c2d3ab61888e8120c1d771bd8c002a` | `1751e03f0171fbdafefe57caa17d1903aaebce8f4fd8d11bd401f7ba6cc84f25` | INTENTIONALLY_CHANGED |
| `CompactStar/Physics/BNV/src/StaticZeroSpinHistory.cpp` | same | `9b96e95cf17e88b595975e23d1c4b10c9ecffec98992b6e3e105eb3d0249dba6` | `9b96e95cf17e88b595975e23d1c4b10c9ecffec98992b6e3e105eb3d0249dba6` | IDENTICAL |
| `CompactStar/Physics/CMakeLists.txt` | same | `3a83f6722593451164919ce52dd34de6dddd07bc46df2f6ed45a90f00a001e33` | `3a83f6722593451164919ce52dd34de6dddd07bc46df2f6ed45a90f00a001e33` | IDENTICAL |
| `tests/CMakeLists.txt` | same | `1cfc5c002e02efdaa3eac9b1e668f73ab45f1d039daa4d7a242af6ee85ff9029` | `1cfc5c002e02efdaa3eac9b1e668f73ab45f1d039daa4d7a242af6ee85ff9029` | IDENTICAL |

The BNV CMake file changed only to install/compile the adapter from its new BNV
paths. The other eight intentionally changed production files changed only
include paths, namespace qualification, or dependent tangent parameter/member
types. Phase-5 build/source ownership is untouched.

## 4. Recovered test manifest

There are 15 recovered test/support files. Five changed only by the mechanical
include/type relocation. Historical `tests/bnv/produce_candidate.py`, SHA-256
`4a82eab9380d5deec70ee74ae26447c0ea139f1bcab8f971ee66ec59a3cb1190`,
is **NOT_RECOVERED** because R0-R3 neither produces nor serializes a Phase-6
candidate.

| recovered path | historical SHA-256 | recovered SHA-256 | classification |
|---|---|---|---|
| `tests/bnv/campaign_fixture.hpp` | `e1222d20572299a61bbec02c3e1b4c565dcd00ba1797db1d25fa7ebb51e7c08c` | `a1e5be5b831c4f384a8dd511152919c1b241691dabc068c471668f6eef258250` | INTENTIONALLY_CHANGED |
| `tests/bnv/compare_candidate.py` | `2f010d6ec3a35696567dd30d7c59764a4ea84fbdeb57e90bd8343a70e3e9a843` | `2f010d6ec3a35696567dd30d7c59764a4ea84fbdeb57e90bd8343a70e3e9a843` | IDENTICAL |
| `tests/bnv/controlled_neutron_sink_fixture.hpp` | `ebf302a634d71323dc7b8b3d9647e0e1a19097065948d709efb462b2a146a642` | `ebf302a634d71323dc7b8b3d9647e0e1a19097065948d709efb462b2a146a642` | IDENTICAL |
| `tests/bnv/coupled_trajectory.cpp` | `52d7d219c3f5e60e823649b3873610bc51213474cda269b266cb6df224af5cfd` | `52d7d219c3f5e60e823649b3873610bc51213474cda269b266cb6df224af5cfd` | IDENTICAL; NOT EXECUTED |
| `tests/bnv/direct_energy.cpp` | `6c4b54482952d6021e396a2cef6e07d4fca290ba5f3ded1959fdf3f13ff9a481` | `6c4b54482952d6021e396a2cef6e07d4fca290ba5f3ded1959fdf3f13ff9a481` | IDENTICAL |
| `tests/bnv/energy_partition_fixtures.hpp` | `9579c53a33f2ad60be9579e613a9e402c6e3a3f42e37999d20832962e8822afe` | `9579c53a33f2ad60be9579e613a9e402c6e3a3f42e37999d20832962e8822afe` | IDENTICAL |
| `tests/bnv/frozen_certificate_fixture.hpp` | `fbbab21644bebd77c41f1f256a92f9b2072ed6bda65fd88b20751b0e5494212d` | `96aa9189b59da73000c8f5608588b2e0667f6d874ef8f3c8730278fc6afacbc2` | INTENTIONALLY_CHANGED |
| `tests/bnv/frozen_monitor.cpp` | `cef2b184969cfc69a21c608f80f89372d68eb2200833ae57435e946ae92e84fe` | `cef2b184969cfc69a21c608f80f89372d68eb2200833ae57435e946ae92e84fe` | IDENTICAL |
| `tests/bnv/frozen_sensitivity.cpp` | `efc7af412f7cc99e7b8eda452c90b2e2ea643e46d1113505cd6d3060d9fd4416` | `30d59d0026902246a2e96778914b6a79f0aaf4a10ac76d519390d13bba314e80` | INTENTIONALLY_CHANGED; NOT EXECUTED |
| `tests/bnv/matched_control.cpp` | `6c2ef89c4eb0f97e2540276bc5a13441bb8a5aae432c12e4337615f6d8975efe` | `6d1b5b85cc7131a003b441d14c51784e14ccdfc66812298d486fe765425227df` | INTENTIONALLY_CHANGED |
| `tests/bnv/ode_refinement.cpp` | `b0214a594751e6962aba55ce47883ed1a0c7f1845c77cd10046580d1f1a7a702` | `b0214a594751e6962aba55ce47883ed1a0c7f1845c77cd10046580d1f1a7a702` | IDENTICAL; NOT EXECUTED |
| `tests/bnv/produce_frozen_certificates.py` | `fb609dbd549763800cc47e962e4f12333835916ea299058388dcb6c1c592ad79` | `fb609dbd549763800cc47e962e4f12333835916ea299058388dcb6c1c592ad79` | IDENTICAL; NOT EXECUTED |
| `tests/bnv/qss_beta_bounds.cpp` | `5814222fa29eaf5dce349908c0f2753f8d6e173df8fb600b6b64b59a11a342f1` | `5814222fa29eaf5dce349908c0f2753f8d6e173df8fb600b6b64b59a11a342f1` | IDENTICAL; NOT EXECUTED |
| `tests/bnv/source_projection.cpp` | `d8d2e7ff70aa4486da42dc592cbdd28c8a321b949fdce39d48af2c4b459d7522` | `be932fdcce5d74bfaf7b5dc42e4eae4a14e800e4a178650851235d82f16ff3c7` | INTENTIONALLY_CHANGED |
| `tests/bnv/thermal_ledger.cpp` | `e3ecc6363bc8ada634cf1336dcc8ab02678d7e2a74303d06cc6480bf2e5f9e62` | `e3ecc6363bc8ada634cf1336dcc8ab02678d7e2a74303d06cc6480bf2e5f9e62` | IDENTICAL |

## 5. R1 semantic-neutrality result

The clean Debug build used
`build/phase6a1-recovery-r0-r3-debug`. Focused tests used independent roots.
The required accepted fixture remains:

```text
raw t = (0.9657700849496014,
         0.030852171225661786,
         0.0033777438247248118)
errors = (1.4428870413124226e-7,
          4.3530476723815904e-9,
          2.1785265203189970e-9)
tau_t = 1.5082028542937021e-7
canonical B_B = 1.6831408136820063e59
canonical B_B error = 1.2652133015974696e52
```

The closed-accessor output is
`(0.96577008494961347, 0.030852171225661786,
0.0033777438247248118)`. The fixture values/errors and required behavior are
the predeclared recovery authorities
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:174-189`).

| focused obligation | result |
|---|---|
| BA2 tangent fixture, `B_B`, closure, identities/currentness | PASS |
| BA3 moving-reference projection | PASS; worst baryon residual `4.2064129956997931e-12`, worst lift residual `4.0927261579781771e-12` |
| BA4 neutron-sink drive | PASS; `(-1.4302859054377826e-55, -3.6280967205205911e-55)` |
| BA5 sliding null / negative `t`, `k`, raw-`G_y` controls | PASS |
| M2 `k` substitution | DETECTED |
| M17 stale tangent | DETECTED |
| M18 source/t domain mismatch | DETECTED |
| nonzero-eta R18 | PASS; residual `-20548825889.525021` |
| generic R18 | PASS; residual `-1.4108007606527171`, budget `1639.3180435981362` |
| BA10a spin-on construction/RHS | PASS; bit identity |
| BA10b zero-spin construction/RHS | PASS |
| direct-energy ledger | PASS; M13-M16 detected |
| retained 21-star certificate refit | PASS; max utilization `0.54075347894390757`, index 20, `N_mu` |
| retained runtime currentness | PASS; just-inside `-9.9999997155751408e-7`, just-outside `-1.0000000283979203e-6`; M20 detected |

The retained certificate output hashes are JSON
`6643054d8cba36035af15a15615baf9920d9a56adf3f5211f13ea9f2745f4409`
and TSV
`9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab`,
identical to historical evidence. The 21 expensive structural star solves were
not rerun, as required when retained evidence authenticates
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:440-455`).
Every comparison-bearing result matched the historical recorded result; no
tolerance was changed and there was no unexpected semantic value change.

## 6. R2 Phase-5D governed regression and BA15 closure

R2 ran serially from the previously absent recovery-specific root
`build/phase6a1-recovery-r0-r3-debug/r2-phase5d`. The fresh producer/comparator
returned zero and all positive/negative controls passed.

| R2 evidence | result |
|---|---|
| governed scientific provenance entries | `91` |
| fresh scientific provenance entries | `91` |
| exact map/set difference | empty |
| old Analysis tangent header in fresh provenance | absent |
| old Analysis tangent source in fresh provenance | absent |
| governed artifact SHA-256 | `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd` |
| fresh artifact SHA-256 | `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd` |
| baseline/producer/comparator edit | none |
| R2 result | PASS |

This is the exact proof required by the recovery plan
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:191-205`). It
closes the historical BA15 root cause: both formerly offending Analysis paths
are absent, while the governed 91-entry scientific source map and artifact bytes
are restored exactly.

## 7. R3 focused regression

After R2 passed, R3 reran the affected source/projection/direct/wrapper and
retained-certificate paths. Source projection reproduced BA2-BA8, M2, M17,
M18, and R18 exactly. Direct energy passed with M13-M16 detected. BA10a and
BA10b passed. The retained certificate refit/currentness and its two hashes
were exact. The focused thermal-ledger executable passed BA8, BA9, and its
static energy-partition BA17 fixture with M8-M12 detected. The static fixture
does not run `phase6a1_coupled_trajectory`; no BNV trajectory executable was
run. R3 result: **PASS**.

## 8. Timing and resource evidence

Peak RSS was not available non-invasively on this host: macOS `/usr/bin/time -l`
was denied access to the required `sysctl`. CPU and wall times below are from
`/usr/bin/time -p` where available. The initial configure itself succeeded, but
its `time -l` wrapper returned nonzero after configuration solely because of
that denied RSS query.

| stage/job | output or scratch root | wall s | user s | sys s | classification |
|---|---|---:|---:|---:|---|
| configure Debug | `build/phase6a1-recovery-r0-r3-debug` | 4.11 | 1.52 | 1.73 | SERIAL_SHARED_ROOT |
| minimum target build | same build root | 15.52 | 51.76 | 9.06 | SERIAL_SHARED_ROOT |
| R1 source projection | `r1-source-projection.z1z81k/work` | 192.51 | 186.49 | 2.47 | SAFE_PROCESS_PARALLEL |
| R1 direct energy | no output files | 0.28 | 0.00 | 0.00 | SAFE_PROCESS_PARALLEL |
| R1 matched control | `r1-matched-control.fXEqDh/work` | 388.78 | 378.10 | 5.07 | SAFE_PROCESS_PARALLEL |
| R1 certificate refit | `r1-retained-certificate.OiTiiX` | 1.20 | 1.08 | 0.07 | SERIAL_SHARED_ROOT |
| R1 frozen monitor | retained certificate TSV | 0.42 | 0.01 | 0.00 | SAFE_PROCESS_PARALLEL |
| R2 fresh Phase-5D regression | `r2-phase5d/regression-savjodiw/fresh-governed` | 2572.40 | 2492.26 | 57.56 | SERIAL_PROVENANCE_AUTHORITY |
| R3 source projection | `r3-source-projection/work` | 191.88 | 189.80 | 1.63 | SAFE_PROCESS_PARALLEL |
| R3 matched control, first execution | `r3-matched-control/work` | approximately 390; terminal stream truncated | unavailable | unavailable | SAFE_PROCESS_PARALLEL |
| R3 direct energy | no output files | 0.01 | 0.00 | 0.00 | SAFE_PROCESS_PARALLEL |
| R3 certificate refit | `r3-retained-certificate` | 1.18 | 1.10 | 0.06 | SERIAL_SHARED_ROOT |
| R3 frozen monitor | retained certificate TSV | 0.01 | 0.00 | 0.00 | SAFE_PROCESS_PARALLEL |
| R3 matched control confirmation | `r3-matched-control-confirm/work` | 388.50 | 385.52 | 2.22 | SAFE_PROCESS_PARALLEL |
| R3 thermal-ledger build | shared Debug build root | 1.48 | 0.87 | 0.39 | SERIAL_SHARED_ROOT |
| R3 thermal ledger | no output files | 0.30 | 0.01 | 0.00 | SAFE_PROCESS_PARALLEL |
| final protected manifest check | `final-protected-manifest.json` | 0.30 | 0.11 | 0.11 | SERIAL_PROVENANCE_AUTHORITY |

The first R3 matched-control process exited, but its orchestration
stream was truncated before the terminal lines and timing could be retained.
It was therefore repeated once, without a source/fixture/budget change, in a
new output root; the confirmation supplies the reviewable BA10a/BA10b PASS and
timing above. Source projection and the first matched-control process ran in
parallel with independent roots as permitted by the recovery plan
(`docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_PREFLIGHT.md:485-504`).

## 9. Protected-authority and commit evidence

The repository helper regenerated and checked the exact 33-path Phase-5D
protected manifest: 33 protected paths, ten entry-era baselines, and three
special artifacts passed. All 11 current files under `tests/baselines` are
byte-identical to canonical entry. Tree comparison also found no change to
Phase-5B, Phase-5C, the Phase-5D producer/comparator/schema, Rotochemical
production, EOS/data, or literature paths. The separately authenticated
literature manifest passed 22/22 files.

Reviewable commits before this validation record are:

1. `1ded3000dc67aa70f5929ed350a372f829c08ff6` —
   `feat: recover phase6 bnv implementation ownership`
2. `465652f1bb9a3eef01884c0fef12ef34610362c6` —
   `test: verify phase6 bnv recovery neutrality`

The validation record and architecture/status updates are a separate commit.
No commit was amended or squashed after evidence referred to it.

## 10. Stop/pass disposition

**A. R0-R3 RECOVERY PASS — ADR-0016 TANGENT OWNERSHIP SEMANTICALLY NEUTRAL —
PHASE-5D PROVENANCE EXACTLY RESTORED — READY FOR SEPARATE BA12R ULTRA TASK.**

No blocker remains within R0-R3. The only nonblocking observation is that peak
RSS was unavailable non-invasively and the first R3 matched-control terminal
stream was truncated; the exact confirmation rerun closes the latter evidence
gap.

The exact next action is a separate, newly authorized BA12R execution task. It
may run only the `CPL-P2-LINEAR-QSS-v1` ULTRA source and exact matched no-BNV
ULTRA control, followed by the predeclared three-level BA12R comparison. This
recovery task does not start it.
