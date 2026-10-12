# CompactStar Linux transition — independent scientific and numerical review

Reviewer: Claude Fable 5.1 (independent, read-only). Date: 2026-10-11.
Host: ekucn1.eku.edu login node (AMD EPYC 7232P, Zen 2, FMA/AVX2, no AVX-512; same feature
class as compute nodes phyn002-004, EPYC 7352). All bounded experiments ran on the login node
against the staged sysroot libm (glibc 2.34-231.el9_7.10), which is byte-identical to the login
node's `/lib64/libm.so.6` and is the libm every candidate binary loads through its rpath
(`ldd toolchain/build/CompactStar-Debug/tests/rotochemical_trackr_pe`). No repository, baseline,
lock, or certificate was modified. No BNV trajectory, thermal regression, or campaign was run.
Scratch work, scripts and outputs: `experiments/` beside this file; checksums in `SHA256SUMS`.

Path convention below: `W` = `/mnt/sdd/zaki/CompactStar/linux-transition-20261010`;
`CS` = `W/sources/CompactStar`.

---

## 0. Disposition

**C — Material corrections or independent numerical authority needed.**

- No Linux-produced number was found outside an independently justified physical or numerical
  bound. Where an independent high-precision oracle exists, the Linux arithmetic is equal to or
  better than the Mac arithmetic, with one exception: the Track-R neutron-onset density, where
  Linux is 7 ULP from the exact value and Mac happens to be correctly rounded (§5.1). That
  difference is 8e-16 relative and physically irrelevant, but it is a real precision loss on x86.
- No accepted contract permits reclassifying any exact-oracle failure. The CONFIND acceptance,
  the Zaki floating-point contract and both Phase-5 portability ratifications declare Mac-only
  authority or forbid excusing a changed number. Every failed test therefore remains FAIL.
- Qualification is blocked by missing authority, not by a Linux defect: there is no governed
  cross-platform acceptance contract, no Linux certificate authority, and four protected-source
  pins are stale (three from Mac-side commits, one from this candidate).

---

## 1. Authenticated candidate identity

| Item | Value | Verified by |
|---|---|---|
| CompactStar branch / HEAD | `candidate/linux-transition-20261010` @ `f4e7c18f568045febc5bf65c6e95200752196dad`, working tree clean | `git rev-parse HEAD`, `git status` |
| Readiness ancestor / cluster prep / master | `5433fc2` / `d8e7208` / `9a5c1ec` (master unchanged) | `git log`, receipt |
| Zaki / CONFIND checkouts | `e263a6e180c5c417198e7778bd21fc9c0a32dc33` / `b0cbd510fd3fd0c772fa50499cd749287cb39e7b` (detached, clean) | `git rev-parse HEAD` |
| REPORT.md | SHA-256 `8c0735f5…8953` matches `candidate-receipt.json` | `sha256sum` |
| Mac lock | `cmake/dependency-lock.json` SHA-256 `e1307a04…89c5`, unchanged | `sha256sum`; `git diff 5433fc2..f4e7c18` |
| Linux lock | `cmake/dependency-lock-linux.json` byte-identical to `W/toolchain/dependency-lock-linux.json`; status `CANDIDATE_UNQUALIFIED` | `cmp` |
| Six frozen Phase-6 inputs | independently hashed: frozen.tsv `9217e278…`, coefficients.tsv `3efa060d…`, certificate.txt `7fc892b5…`, entry.json `5cfbf4b1…`, profile_tree/profile.tsv `e9cd03b0…`, profile_tree/freegas.tsv `7cd44c92…`; the tree-level hash `23311486…` is the implementation agent's recorded value in `execution.json`, not independently recomputed | `sha256sum` vs readiness doc and `authentication.json` |
| Baselines / fixtures / governed JSON | no change between `5433fc2` and `f4e7c18` | `git diff --stat -- tests/baselines tests/dependency_migration/fixtures docs/validation/*.json` (empty) |
| Push state | not on GitHub; incremental bundle `2befc706…` verified, requires `5433fc2` on remote | `bundle-verification.txt`, `handoff.json` |

Change set `5433fc2..f4e7c18`: 363 files. Source-class content: ten C++ files (standard-header
includes, `using enum` → qualified names in `Tags.hpp`, libc++ clock alias in `TaskManager.hpp`),
three CMake/Python authentication files, four test scripts, nineteen `tools/linux` files. The
remainder is documentation and gzip evidence. No equation, constant, tolerance, baseline or
fixture changed. One changed test script is a frozen-context protected path (§7.4).

---

## 2. Independent test and failure classification

Registered inventory: Debug 86 (72 PASS in scope + 2 out-of-scope 1 Myr PASS + 10 FAIL + 2
deferred), Release 83 (72 PASS + 7 FAIL + 4 deferred). My classification of every failure:

| Test (modes) | Class | Mechanism | Numerical defect? |
|---|---|---|---|
| CONFIND `historical_serial`, `historical_samples_workers_{1,2,3,4,6,8}`, `scheduler_stress_120_runs` (8/mode) | **Text-identity artifact** | Single differing line (fixture `C13_explicit_radius_exact_ulp`): `%a` rendering of the subnormal 2^-1074 (`0x0.0000000000001p-1022` on glibc vs `0x1p-1074` on Apple); the printed 64-bit pattern `0000000000000001` is identical | No. Not numerical at all |
| CONFIND `m1_historical_small`, `m1_small_workers_{…}`, `m1_scheduler_stress_120_runs` (8/mode) | **Mac-libm rounding encoded in oracle** | pow(10,x) at a deliberately constructed near-tie; glibc is correctly rounded, Apple is not (§3.2); exact-equality de-duplication turns one ULP into 4-vs-6 curve points | No. Linux is the mathematically correct side |
| CONFIND `m1_historical_stellar`, `m1_stellar_workers_{…}` (7/mode) | libm ULP drift | 418/81457 rows, ≤32 ULP, max rel 4.1e-15 through log-axis pow/log10/exp | No |
| CONFIND `m1_randomized_{log,linear}_600` (2/mode) | libm ULP drift | order-sensitive bit hashes of all points; any ULP change alters the hash | No |
| CONFIND `installed_package_consumer` | environment | missing GSL/zlib search paths; passes with explicit prefixes | No |
| `phase5b_structural_response_regression` (Debug) | exact cross-compiler identity gate | values: A,N ~1e-15; B,K,I_phys ~5e-11; error-estimate fields 1e-3; residuals up to 86% (roundoff-dominated); regenerated EOS table hash differs (§4.3) | No |
| `phase5c_chemical_coefficient_regression` (Debug) | exact cross-compiler identity gate | G/Q/Z ≤2.7e-14, I/W ≤5.8e-11; provenance constants string embeds platform-dependent onset values | No |
| `phase5d_protected_manifest`, `phase5d_coupled_oracles` | **pre-existing stale source pins** | `ChemicalResponse.cpp` hash at `5433fc2` already ≠ pin; fails on any platform after 2026-10-03 (§7.4) | Not Linux-related |
| `rotochemical_trackr_pe` PE-V6 (both) | **latent production precision loss exposed by a wider oracle** | production onset +7 ULP on x86; Mac correctly rounded only through FMA contraction (§5.1) | Real but negligible (8e-16) |
| `relativistic_unit_boundary` U11 (both) | **fixture construction is platform-arithmetic dependent** | fixture value 1 ULP below threshold on x86, 1 ULP above on arm64 (§5.2) | No |
| `baryon_number_cmf` (both) | same-grid bound applied cross-grid | dR ≤1.3e-10, dB/B ≤5e-12 Debug; within-Mac Debug/Release already 1.1e-10 / 6e-13 (§5.3) | No |
| `hartle_monopole_regression` (both) | byte identity | worst 6.5e-10 Debug / 1.2e-10 Release; within-Mac Debug/Release 5.7e-10 (§4.4) | No |
| `taskmanager_T2`, `taskmanager_observations` (both) | exact hash identity of legacy pipeline | T1: 37/53484 records at 1 ULP. T2: lifetimes differ up to 2.2/2.5/4.3 % (§6.3) | Not a Linux defect; **material conditioning finding** (cross-platform drift ≈4× the within-Mac build-mode drift) |

Positive evidence the failures concentrate in exact-identity gates: every tolerance-gated
analytic, EOS, Phase-5B/5C/5D component, breakpoint, ledger and readiness contract passes on
Linux in both modes, including the two 1 Myr scalar thermal regressions that completed before
the scope correction, and the Linux coefficient artifact is byte-repeatable
(`phase5c-within-linux-repeat.json`). Within Linux, CONFIND worker counts 1–8 and Debug/Release
produce identical tables (`cmp`, this review).

---

## 3. CONFIND root causes

### 3.1 Governing authority
`CONFIND_2_0_ACCEPTANCE.md` establishes "library-level authority on the qualified local Mac
configuration only" and explicitly does **not** establish Linux qualification or cross-platform
bitwise authority (requirement 7: "Keep Linux and EKU cluster qualification later and separate").
The 25 exact failures per mode are therefore the expected outcome of running Mac-only oracles on
Linux. They remain FAIL under the existing contract; nothing in them is a defect of the port.

### 3.2 The m1-small "pow split" fixture (decisive)
The fixture (`CONFIND/tests/m1_characterize.cpp`) evaluates `HistoricalPow10(x)` = libm
`pow(10,x)` at interpolated log-coordinates. I reconstructed the two inputs from `Cell.cpp`
case36/case5 and evaluated them with the sysroot libm (ctypes) and a 60-digit mpmath oracle
(`experiments/outputs/confind_pow.txt`):

| input x | exact 10^x relative to rounding midpoint | correctly rounded | glibc pow | glibc exp10 | Mac table |
|---|---|---|---|---|---|
| `0x1.22fad6cb53501p-10` (the level, RAW 4/5) | −0.0001 ULP (near-perfect tie) | `…b13f` | `…b13f` | `…b13f` | `…b140` |
| `0x1.22fad6cb53600p-10` = 1−(1−L) (RAW 0/2) | +0.575 ULP | `…b140` | `…b140` | `…b140` | `…b140` |

Apple's pow returns the upper neighbour at a point where the true value is one ten-thousandth of
an ULP below the midpoint; glibc's value at L equals Apple's exp10 value, which is why the
fixture comment reads "pow=…140, exp10=…13f" on Apple. The frozen reference encodes that rounding; glibc is correctly
rounded. The acceptance document already records that the "unchanged radius-based
de-duplication amplifies these to four versus six curve points"; the amplification is a
discontinuity of the historical algorithm, not a Linux behaviour.

### 3.3 The serial family
`diff toolchain/build/CONFIND-Debug/tests/serial.tsv CONFIND/tests/reference/serial.tsv` is a
single line pair: RAW 4 of fixture `C13_explicit_radius_exact_ulp`, whose explicit point uses
`nextafter(0.0, 1.0)` = 2^-1074, differing only in the `%a` spelling of that subnormal; both carry
bit pattern `0000000000000001`. Eight tests fail on text, not numbers.

### 3.4 Symbols and flags
`nm -u` on all four installed archives: CONFIND imports `pow`(1) and `log10`(2); Zaki imports
`pow`, `log10`, `exp`, `log`, `asinh`, `sqrt`; **no `exp10` anywhere**. HistoricalMath TUs are
compiled `-ffp-contract=off -fno-builtin -fno-lto` in both modes
(`floating-point-build-audit.json`). The CONFIND/Zaki `pow`→`exp10` hazard addressed by their
M-1 corrections is absent on Linux. The remaining differences are pure libm rounding.

---

## 4. Zaki and CompactStar numerical differences

### 4.1 Zaki historical tables (Mac-only authority by CMake gate, `Zaki/tests/CMakeLists.txt:59`)
I reconstructed every log-axis node of the fixed table from the harness formula and reproduced
the Linux table exactly (46002/46002 random-axis nodes, 251/251 stellar-range nodes), which
validates the attribution below (`experiments/outputs/zaki_axis_attrib.txt`):

| random-axis differing nodes (278 total) | count |
|---|---|
| glibc `log10` returned a non-correctly-rounded input; glibc pow correctly rounded at its own input; Mac value equals correctly rounded 10^(correctly rounded input) | 207 |
| identical input; glibc pow correctly rounded, Apple pow not | 53 |
| identical input; Apple pow correctly rounded, glibc pow not | 17 |
| glibc input not correctly rounded and Mac value also not | 1 |

Survey: glibc pow(10,x) is not correctly rounded for 11 of 20000 random inputs in the stellar
exponent range; the reconstructed axis input x_i differed from the correctly rounded input for
289 of 46002 nodes (each input depends on log10 at two axis endpoints, so no per-call rate is
claimed). Neither libm is uniformly correctly rounded; all differences are last-bit. Section maxima: axes/grids ≤1.1e-15 relative;
`integral_1/2/3` up to 9.2e-7 relative at near-zero values (absolute ≤2.8e-14) from the
asinh-based closed form in `IntegralTable.cpp:48-57` (difference of nearly equal terms);
`stress_grid` up to 1184 ULP (2.1e-13 relative) from interpolation-weight cancellation.

### 4.2 Phase-5C coefficients G, Q, Z, I, W
From `phase5c-coefficient-differences.json` (Linux Debug production vs governed Mac artifact,
Structure-1 midpoint, ρc = 1.10e15):

| component | max relative | max ratio to supplied combined `numerical_error` |
|---|---|---|
| G | 2.66e-14 | 7.9e-8 |
| Q | 1.15e-14 | 5.8e-8 |
| Z | 1.19e-14 | 5.7e-8 |
| I | 5.79e-11 | 1.6e-4 |
| W | 5.78e-11 | 9.8e-5 |

Two classes, each with an independently justified origin: closed-form quantities (G, Q, Z;
Phase-5B A, N) drift at the libm/long-double-width level ~1e-14; structure-derived quantities
(I = K/c², W = Z·I; Phase-5B B, K, I_phys; radii, Hartle integrals) drift at ~5e-11, which is
the rk8pd driver tolerance (`TOVSolver.cpp:2501,2557`: eps_abs = eps_rel = 1e-10) propagated
through the solve. The ratios show the cross-platform differences are 4 to 8 orders of magnitude
inside the artifacts' own declared numerical uncertainties. Under the ratified portability rule
("a changed compiler string never excuses a changed number") the exact regression nonetheless
remains FAIL.

### 4.3 Provenance strings make exact identity unreachable even if numbers agreed
The Phase-5B/5C artifacts embed the provider's onset constants in
`particle_constants_and_provider_source` and in the regenerated EOS table's authority string:
neutron onset `7.3567289037328318e-09` (Mac) vs `7.3567289037328375e-09` (Linux), muon onset
`0.45698480541241987` vs `0.45698480541241826`, hence different `particle_constants_source_sha256`
and `EOS_table_sha256` (`6505db38…` vs `7cd44c92…`). The regenerated free-gas table itself drifts
(per column, Linux vs handoff): eps 9.6e-13, pressure 4.1e-12, n_B 3.6e-15, species fractions
≤5.7e-13, one fraction of order 1e-14 at 5.8e-7 relative (absolute 4e-17).

### 4.4 Hartle monopole and CMF structure (regenerated with capture, `experiments/hartle-captures/`)
| pairing | worst column relative drift |
|---|---|
| Linux Debug vs Mac Debug | 6.5e-10 (`surface_shell_mass_over_Omega2_km3`), R 1.3e-10 |
| Linux Release vs Mac Release | 1.2e-10 |
| control: Mac Debug vs Mac Release | 5.7e-10 |
| control: Linux Debug vs Linux Release | 1.2e-10 |

Cross-platform drift equals the within-platform build-mode drift already accepted by the
project through separate Debug/Release references.

---

## 5. PE-V6, U11 and CMF verdicts

### 5.1 PE-V6 — verdict: legitimate detection of a latent production precision loss; test premise invalid on arm64
Reproduced bit-exactly (`experiments/outputs/cxx-experiments.txt`):

| quantity | value | ULP from exact |
|---|---|---|
| exact onset from the binary constants (mpmath, 60 digits) | 7.35672890373283169196e-9 | — |
| correctly rounded double | `0x1.f98cf18cad246p-28` = 7.3567289037328317540e-9 | 0 |
| Linux production (`TrackRFreeGasThermodynamics.cpp:219`, `LocalThermodynamics.cpp:198`) | `0x1.f98cf18cad24dp-28` | **+7** |
| Linux 80-bit independent expression | 7.3567289037328316926e-9 | 8.5e-20 relative |
| Mac production (governed artifact constant `7.3567289037328318e-09`) | `0x1.f98cf18cad246p-28` | 0 |
| Mac independent (Phase-5A5 review log `7.356728903732833408335e-09`) | reproduced only with FMA contraction on | — |

Cause: production rounds μ_e and ħc to double before its long-double section; on x86 the
subsequent 80-bit arithmetic faithfully propagates those input roundings (≈3× each through the
cube), giving +7 ULP. On arm64 `long double` is double, AppleClang contracts
`(mn−mp)(mn+mp)+me²` into one fused operation (no `-ffp-contract` flag exists in the CompactStar
build), and the result is correctly rounded by accident. The Phase-5A5 independent review
(P-2) warned that the oracle "cannot exceed double precision" and that "a future platform change
[could] quietly invalidate" it. The 3e-16 gate sits below the formula's own double-rounding
error bound (~8e-16), so it is unattainable by the current formula against an exact reference.

Checkable fixes (ADR required, numerical-method class; the test is ratified):

| variant | x86 (executed) | arm64-equivalent (simulated: `long double` = double, no contraction) | note |
|---|---|---|---|
| current | +7 | +5 | Mac's actual 0 is a contraction accident |
| μ_e and ħc computed in long double | 0 | +5 | platform-dependent |
| **double throughout with explicit `std::fma` for μ_e, no long double** | **0** | **0 predicted** | bit-identical to the governed Mac constant on x86; predicted on arm64 because basic operations, sqrt and fma are correctly rounded everywhere and `acos(-1.0)` is assumed to return the correctly rounded π; to be confirmed on the Mac before ratification |

### 5.2 U11 — verdict: portability-sensitive test fixture; production correct
`relativistic_unit_boundary.cpp` builds the synthetic column as
`(1.0001 − i·0.0001)·1e10L·G/c/c·1e10L` with `long double` G and c. Since `1.0001 − 0.0001 == 1.0`
exactly, layer 1 is meant to equal the threshold, but the fixture reaches it through a
different operation order and width than the production conversion
(`RelativityUnits.hpp:18,25`, pure double, compile-time constants). Reproduced
(`experiments/scripts/u11.cpp`): x86 80-bit gives `0x1.fe38f2694b4c7p-28`, one ULP **below** the
threshold `…4c8` (inverse 9999999999.9999981, matching the log), arm64-equivalent gives `…4c9`,
one ULP **above**. The expected index 1 depends on the accidental sign of a one-ULP mismatch.
`FindTbIndex` (`TbDefinition.cpp:62`) behaves as contracted. Fix: construct the layer through
the production forward conversion of exactly 1e10 (equality by construction) or use a value
strictly above the threshold; rule and expectation unchanged.

### 5.3 CMF baryon number — verdict: same-grid bound applied across grids; no error
ADR-0004 §7.1 derives the 1e-15 bound for "the same quadrature on the same grid" with common-mode
summation. Cross-platform, the grid itself comes from the ODE. Evidence: central density
bitwise identical for all four stars (dec = 0), so the mass bisection
(`SolveToProfile`, mass_tol 1e-4 M☉) took identical discrete decisions; all drift is continuous
ODE drift at or below rtol 1e-10:

| pairing | worst dR | worst dB/B |
|---|---|---|
| Linux Debug vs Mac Debug | 1.26e-10 | 5.0e-12 |
| Linux Release vs Mac Release | 2.5e-11 | 2.6e-13 (two stars bitwise identical) |
| control: Mac Debug vs Mac Release references | 1.1e-10 | 6.0e-13 |

Debug differing more than Release is consistent with the TOV right-hand side
(`TOVSolver.cpp:1407-1535`) calling libm `pow(r,2.)`, `pow(r,3.)`, `pow(c,2.)` at `-O0`, whereas
both optimizers replace `pow(x,2.)` and constant powers by exact multiplications; I did not
verify the generated code.

---

## 6. Frozen profile, certificate and the legacy T2 pipeline

### 6.1 Frozen profile differences are localized and explained
From `profile-diagnostic/profile.tsv` (20259 rows; Linux star regenerated from the authenticated
freegas.tsv, no certificate consumed):

| column | rows differing | where | max | explanation |
|---|---|---|---|---|
| r | 118 | rows 1–142 (≤1 ULP) and row 20258 (108 ULP = 1.92e-13 km) | 1.5e-14 rel | central grid rounding; surface from the terminal pressure-coordinate step `LocateSurface` |
| m | 4398 | throughout | 7 ULP, 6.7e-16 abs | ODE drift |
| ν | all | all | 84 ULP, 1.2e-15 abs | surface boundary condition ν(R) offsets every row by 1.4e-16 |
| n_B | 2947 | 6812–10522, surface | 1.1e-16 abs (1.1e-10 rel near zero) | EOS interpolation at tiny values |

### 6.2 Certificate verdict
The guard is exact double equality on every row (`chemical_production_fixture.cpp:68`,
`rotochemical/fixture.hpp:73`); refusal is correct fail-closed behaviour. Numerically, the
profile perturbation is ≤1.5e-14 relative in r and ν and ≤1e-15 in m; the certificate's
quantity-level drift bounds are ≥6e-15 and its thresholds ≥3.4e-7 (frozen.tsv: 3.4e-7 for t_mu,
3.1e-6 for t_e, 1e-4 for the rest; maximum utilization 0.54 for N_mu). The perturbation is below
every threshold by seven orders of magnitude, but these are different objects: I did not
propagate the perturbation through the certificate's sensitivities, so this bounds
explainability, not rebinding. The differences are explainable and numerically controlled; they
do not authorize rebinding, because the certificate's sensitivities were derived from the Mac
star. A Linux certificate can in principle be generated and qualified
natively — the Linux star solves and is byte-repeatable — but the 21-star producer's identity
and cost are not in the handoff (`authentication.json` records only the Mac build-directory
origin `phase6a1-frozen-certificate-final`), and the Mac readiness record itself marks
certificate currentness UNQUALIFIED. That producer must be identified and run on Linux under a
predeclared contract before any Linux certificate exists.

### 6.3 T2 legacy lifetime pipeline — material conditioning finding (not a Linux defect)
The implementation report records "22 of 38 hashes differ" without magnitudes. Quantified
(`experiments/outputs/t2_compare.txt`, `t2_detail.txt`, `t2_seq.txt`, `t2_macctl.txt`):

| output | Linux Debug vs Mac Debug | control: Mac Debug vs Mac Release |
|---|---|---|
| Fermi-gas EOS table | 97/6006 tokens, ≤2.5e-7 (low-density cancellation) | identical |
| 19×19 TOV sequence | 4 tokens ≤4.9e-9 | 11 tokens, up to **2.4 %** (row 388, extreme dark-core corner); Linux Debug equals Mac Debug there |
| critical curve | 24 tokens ≤1.4e-9 | 24 tokens ≤1.4e-9 |
| B-contour sequences | one point per contour (k=0–5,8) moves, up to 1.6e-3 in ε_c (k=5) | up to 3.7e-4 |
| BNV lifetimes n / Λ / Σ⁻ | **2.2 % / 2.5 % / 4.3 %** | **0.59 % / 0.60 % / 0.66 %** |

Mechanism (code: `TaskManager.cpp:581-590`): each constant-B contour is cut at its intersection
with the critical curve; near the maximum-mass turning point the two curves are nearly tangent,
so a ~1e-9 shift in the critical curve moves the intersection ("death" point) by 1.6e-3 in
ε_c, which the lifetime integration turns into percent-level changes. The Linux Debug value of
the worst entry equals the Mac Release value. Cross-platform lifetime drift (2.2–4.3 %) is roughly four times the within-Mac
build-mode drift (0.6 %); both arise from the same mechanism, but the magnitudes are not the
same. The migration record states "Debug is not
compared to Release" (`COMPACTSTAR_EXTERNAL_DEPENDENCY_MIGRATION_QUALIFICATION.md:137`), so this
intrinsic sensitivity was never characterized. The pipeline's analysis module (`Microphysics/BNV/Analysis/src/BNV_Analysis.cpp`) is
compiled into libCompactStar but is not referenced by `Physics/BNV`, `tests/bnv`, `Core`, or
any Phase-6 document I grepped. If its lifetimes are used
scientifically, they carry an undisclosed numerical uncertainty at the percent level on every
platform; this maps to CONFIND acceptance requirement 5.

---

## 7. Scientific-source and provenance audit

### 7.1 Locks and platform assumptions
- Mac lock byte-preserved; Darwin path refuses any Linux lock; Linux path requires an absolute
  lock and authenticates compiler bytes/version, flags per mode, GSL source digest, closures
  (GSL, sysroot incl. `usr/lib64/libm.so.6` `d4d4aca3…`, compiler tree, tools), link inputs,
  per-mode package archives/headers/configs and provenance (`cmake/linux_dependencies.py`). Six
  unit tests and four negative configurations exercise refusal.
- Compute nodes carry system glibc 2.34-100; the binaries load the staged 2.34-231 libm via
  rpath, so the arithmetic is pinned independently of the node (verified by `ldd`).
- Gap: the lock pins machine type, not CPU feature class. glibc selects FMA/AVX2 libm variants at
  run time (IFUNC). Login and compute nodes are both Zen 2 and my login-node ctypes calls matched
  compute-generated tables bitwise, so it held here; the acceptance contract should pin it.
- Gap: the Mac CompactStar build sets no `-ffp-contract`; AppleClang contracts by default.
  Proof: PE-V6 Mac values are reproduced only with contraction on (§5.1). Mac baselines therefore
  depend on AppleClang's contraction decisions, which no Mac record pins. Linux: all 162
  CompactStar TUs `-ffp-contract=off`, no `-march`, zero FMA instructions in either library
  (`objdump`).
- Gap: `long double` in three production files (`Analysis/src/ChemicalResponse.cpp`,
  `EOS/src/LocalThermodynamics.cpp`, `Physics/Rotochemical/src/GlobalUrcaChannelCoefficient.cpp`;
  imports `cbrtl`, `hypotl`, `expl`, `sqrtl`) makes results platform-width-dependent by
  construction (64-bit on arm64, 80-bit on x86_64). The Mac authority was established at 64 bits.
- ADR-0020 (PROPOSED) never names libm identity, long-double width or contraction as the known
  sources of cross-platform non-identity.

### 7.2 Changed C++ sources (all ten reviewed)
Includes added (`<iomanip>`, `<atomic>`, `<chrono>`, `<filesystem>`, `<algorithm>`, `<cmath>`);
`Tags.hpp` enum qualification; `TaskManager.hpp` `std::__1::chrono::steady_clock::time_point`
→ `std::chrono::high_resolution_clock::time_point` (the only consumers already call
`high_resolution_clock::now()`, `TaskManager.cpp:45,60`; timing only; on libstdc++ this
aliases `system_clock`, not monotonic — nonblocking). No scientific-semantic change found.

### 7.3 Scientific inputs and baselines
All six frozen inputs byte-authenticated (§1). `tests/baselines`, migration fixtures and
governed JSON unchanged. Adapter maps logical input keys to relocated paths without altering
bytes (`execution.json`).

### 7.4 Protected-source pins (frozen-context plan `f7116c1`, 2026-09-11, 33 paths; pins match file states from 2026-09-06/07, e.g. `ChemicalResponse.cpp` at `4d78bf4`)
| path | changed by | date | in master? | nature |
|---|---|---|---|---|
| `CompactStar/Analysis/src/ChemicalResponse.cpp` | `d222b60` | 2026-10-03 | no (readiness ancestry) | Mac-side amendment, pin not updated |
| `tests/analysis/chemical_production_contract.cpp` | `d222b60` | 2026-10-03 | no | same |
| `CompactStar/Analysis/src/ParticleNumberResponse.cpp` | `a8a43e9` | 2026-10-03 | no | same |
| `tests/analysis/produce_particle_number_reference.py` | **`bd07325` (this candidate)** | 2026-10-10 | no | **candidate edit to a protected path** |

The last Mac record of `phase5d_protected_manifest` passing is dated 2026-09-18
(`PHASE6A0_BNV_THERMAL_FIRST_LAW_PREFLIGHT.md:1210`), before the amendments; the guard must fail
on every platform since. The candidate's own edit is mechanical (compiler path selection) and
disclosed in DETAILS.md, but it is an edit to a frozen-context protected file without the owning
authority's re-pin, and must be adjudicated as such.

### 7.5 Unauthorized scientific semantic changes
None found in C++ sources, CMake, baselines, fixtures, locks or inputs. The one governance
violation is the protected-path edit above (change class: engineering/test harness, but the file
is provenance-protected).

---

## 8. Proposed corrective actions

1. **Re-pin the four protected paths through the frozen-context plan's amendment process**
   (owner authority, documentation class), explicitly adjudicating `bd07325`. Unblocks
   `phase5d_protected_manifest` and `phase5d_coupled_oracles` in both modes on every platform.
2. **PE-V6**: ADR for a numerical-method change computing the onset in double with an explicit
   `std::fma` for μ_e and no `long double`; verified bit-identical to the governed Mac constant
   on x86 and predicted identical on arm64 (confirm on the Mac before ratification). Alternatively rederive the threshold from the
   operation count (≥8e-16). Do not weaken the gate without the derivation.
3. **U11**: rebuild the synthetic layer through the production forward conversion (equality by
   construction). Test change to a ratified ADR-0012 test; expectation unchanged.
4. **CONFIND serial family**: compare the 64-bit pattern column, not the `%a` text (or normalize
   subnormal spelling). CONFIND is frozen; this is a CONFIND-side change for its owner.
5. **Pin arithmetic semantics on both platforms**: add `-ffp-contract=off` to the CompactStar
   build (or record the Mac contraction state in the Mac lock), record CPU feature class in the
   Linux lock, and amend ADR-0020 to name libm identity, long-double width and contraction as the
   governed non-identity sources.
6. **Audit `long double` in the three production files**; decide per quantity whether width
   dependence is acceptable or should be replaced by double plus compensated/fma arithmetic.
7. **T2 lifetimes**: characterize the near-tangent intersection sensitivity (refine the 19×19
   grid or the intersection solver) and publish a numerical uncertainty for `BNV_tau` outputs;
   convert the hash gate to a tolerance gate derived from that characterization.
8. **Capture regression tables** (`COMPACTSTAR_UNIT1_CAPTURE_DIR`) in every qualification run
   so Hartle/CMF drift can be quantified without regeneration.
9. **Identify and record the 21-star certificate producer** and its cost before planning any
   Linux certificate.

---

## 9. Proposed cross-platform acceptance methodology

Separate five questions; a test answers exactly one.

1. **Same-platform reproducibility (bitwise).** Within one platform, build mode and lock,
   regenerate twice → identical bytes. Already demonstrated on Linux (coefficient artifact,
   worker invariance, Debug/Release CONFIND tables). Keep exact gates here, platform-tagged.
2. **Cross-platform numerical agreement (tolerance by quantity class).** Derive bounds from
   independent error analysis, not from observed drift: (a) closed-form/libm class — ≤ a few ULP
   for direct function outputs, ~1e-14 relative for short chains, with explicit cancellation
   handling (absolute bounds for near-zero quantities such as `integral_3`, species fractions,
   residuals); (b) structure/ODE class — relative drift ≤ ~10× the integrator rtol (1e-9 here)
   for R, M, B, I, K, W and Hartle integrals, with dec bitwise to prove identical discrete
   bisection decisions; (c) discrete outputs (contour point counts, selected layer indices, curve
   topology, hash records) are **never** compared across platforms by identity — compare raw
   points with tolerance plus a topology check. Report ULP, absolute and relative drift for every
   quantity; missing uncertainty means UNQUALIFIED.
3. **Physical accuracy (platform-independent truth).** Where an analytic or high-precision oracle
   exists (onset densities, unit conversions, Schwarzschild interior, free-gas identities),
   gate each platform against the oracle, not against the other platform. PE-V6 should become
   such a test with a threshold derived from the formula's operation count.
4. **Frozen provenance.** Mac certificates stay Mac-bound. A Linux certificate is produced
   natively by the same producer from the Linux star, bound by the same exact-equality guard,
   and its quantities are compared to the Mac certificate under class 2 bounds. Pin libm bytes,
   CPU feature class, contraction state and long-double width in the lock of each platform.
5. **Governed baseline acceptance.** Linux-native golden baselines (Debug and Release separately)
   are ratified by the owner only after classes 1–4 pass; exact Mac oracles receive an explicit
   platform-authority tag and are reported as NOT-APPLICABLE on other platforms rather than as
   failures, once a ratified contract says so. Until then every current failure stays FAIL.

---

## 10. Blockers, material findings, nonblocking findings

**Blockers (for qualification, not defects of the port)**
- B1. No governed cross-platform acceptance contract; all exact oracles are Mac-only by their
  own accepted documents.
- B2. No Linux certificate authority; producer identity and cost unrecorded.
- B3. Stale protected pins (three Mac-side, one candidate-side) fail two tests on every platform.

**Material**
- M1. PE-V6: production onset is +7 ULP on x86; Mac correctness is a contraction accident; the
  gate is below the formula's error bound (§5.1).
- M2. T2 lifetime outputs change by 2–4 % cross-platform and 0.6 % between Mac build modes from a
  near-tangent intersection; uncertainty undisclosed (§6.3).
- M3. Mac contraction state unpinned; Mac baselines depend on it (§7.1).
- M4. `long double` production arithmetic is width-dependent by construction (§7.1).
- M5. Two ratified tests encode platform-dependent constructions (U11 fixture; CONFIND `%a`
  text identity).
- M6. CPU feature class (libm IFUNC dispatch) not pinned.

**Nonblocking**
- N1. `high_resolution_clock` aliases `system_clock` on libstdc++ (timing only).
- N2. ADR-0020 silent on arithmetic-semantics sources of non-identity.
- N3. Regression tables not captured during the qualification run (regenerated here).
- N4. DETAILS.md reports T2 as hash differences without magnitudes.
- N5. `installed_package_consumer` initial failure was environmental and repaired.
- N6. Zaki `integral_3` relative drift 9.2e-7 is near-zero cancellation (absolute 2.8e-14).

---

## 11. Smallest practical next implementation task

Re-pin the four protected paths through the frozen-context amendment process, with the
candidate's `bd07325` edit explicitly adjudicated (accept and re-pin, or revert and route the
compiler selection outside the protected file). It is platform-independent, documentation class,
weakens no gate, unblocks two tests in both modes, and is required before any Mac rerun as well.
The smallest code tasks that follow, each needing an ADR because the tests are ratified: the U11
fixture construction fix, then the PE-V6 double-plus-fma onset.

---

## 12. What this review did not do
No BNV source/control trajectory, no 1 Myr or 10 Gyr regression, no full suite rerun, no Slurm
job (every experiment ran in under a minute on the login node, same CPU class and identical libm
bytes). The Mac side was not executed; Mac arithmetic was inferred from committed artifacts,
review logs and exact reproduction of their values. The 21-star certificate producer was not
located or run.
