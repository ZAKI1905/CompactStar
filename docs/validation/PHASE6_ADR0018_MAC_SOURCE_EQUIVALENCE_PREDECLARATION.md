# Phase-6 ADR-0018 Mac source-build equivalence predeclaration

**Status:** PREDECLARED / COMMITTED BEFORE BUILD — SOURCE-EQUIVALENCE EXECUTION AUTHORIZED, NOT YET RUN

**Date:** 2026-09-26

**Canonical input:** `812463ac9ed374f64ac9cadd500066ab723d3a6c`

**Experiment branch:** `physics/adr0018-mac-source-equivalence`

## 1. Scope and authority

This record fixes the experiment, comparisons, equality classes, and stop
conditions before candidate compilation, candidate-specific configuration, or
observation of a treatment result.  It implements the experiment authorized by
`PHASE6_ADR0018_ACCEPTANCE.md`; it does not grant dependency source authority,
replace either vendored Darwin archive, authorize cluster work, or alter a
scientific model.

The accepted source-equivalence candidates are only:

- ZakiLib Git commit
  `b9ddebaded24962468954846f47238aec2726fd4`;
- CONFIND exact non-Git source manifest
  `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`.

The authenticated control archives remain:

- Zaki `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef`;
- CONFIND `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2`.

## 2. Pre-run platform gate

The local platform was authenticated before this declaration:

| Field | Required and observed value |
|---|---|
| macOS | `26.6.2` build `25G83` |
| Darwin | `25.6.0` |
| architecture | `arm64` |
| compiler | `/usr/bin/clang++`; Apple clang `21.0.0 (clang-2100.3.34.2)` |
| CMake | `/opt/homebrew/bin/cmake`; `4.2.1` |
| GSL | `/opt/local`; `2.7.1` |
| build type | `Debug` |
| Python | CPython `3.12.4`; one identical resolved interpreter for C0/C1/T |

The configure commands will explicitly bind `/usr/bin/clang` and
`/usr/bin/clang++`.  No unsafe floating-point option (`-ffast-math`, `-Ofast`,
associative math, reciprocal math, or finite-math-only) is permitted.  C0, C1,
and T must have identical non-dependency configure values, Python, GSL, SDK,
OpenMP runtime, environment, inputs, and source revision except for the
predeclared override implementation and dependency identities.

## 3. Configurations

| ID | CompactStar source | dependency selection |
|---|---|---|
| C0 | exact canonical input `812463ac...` | implicit authenticated Darwin defaults |
| C1 | experiment source containing only the minimal ADR-0018 override mechanism | all four explicit overrides point to the same vendored headers and archives used by C0 |
| T | byte-identical CompactStar source and settings as C1 | all four overrides point to separately source-built candidate headers and archives |

C0 therefore isolates the contemporaneous canonical authority, C0 versus C1
isolates the override mechanism, and C1 versus T isolates the external
dependency artifacts.  Separate empty CMake build directories and separate
output roots are mandatory.  No cache is shared.

The planned evidence root is outside CompactStar Git and outside vendored
dependency paths:

`/Users/keeper/Documents/CompactStar/external/equivalence/ADR0018/20260926-812463a-a8ed2907-ed76163c/`

It will contain dependency recipes, commands, source/header/object/archive
manifests, the candidate archives, C0/C1/T configure manifests, outputs, and
comparison reports.  Candidate source repositories and the preserved CONFIND
tree remain read-only inputs.

## 4. Source and input identities

The Zaki tracked-file manifest contains 121 paths and has manifest SHA-256
`a8ed2907812354ab2d46bc941806f7bec1e681eaaba2b7d785dc3de60c0e842b`.
The CONFIND reconstructed source manifest and transfer ZIP must remain,
respectively:

- `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`;
- `be41ee1c71b88627f0b11703c9d732abd7378ebbc162ee10649af771dec9ae0f`.

The bounded ADR-0017 input bindings are fixed as follows; every identity is
authenticated before each run:

| Input | Bound path | SHA-256 / tree identity |
|---|---|---|
| solve matrix | `.../CompactStar-phase6a1-checkpoint-reconstruction-validation/docs/validation/phase6a1_checkpoint_reconstruction_solve_matrix.tsv` | `32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616` |
| profile root | `.../CompactStar-phase6a1-controlled-bnv-implementation/build/phase6a1-controlled-bnv-debug/phase6a1-fresh-qualification/chemical-characterization/run-dhp6aivf/t8192-r80000` | tree `233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72` |
| qualification certificate | sibling `certificate-t8192-r80000.txt` in `phase6a1-fresh-qualification` | `7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453` |
| thermal source | `.../phase6a1-controlled-bnv-debug/phase6a1-shared-inputs/thermal` | tree `1cdb98507958b50fbd1e1eff3eb3bfa4f48c790fef86a8740855884fc43b9d94` |
| entry manifest | `.../CompactStar-phase6a1-controlled-bnv-implementation/docs/validation/phase6a1_controlled_bnv_entry_hashes.json` | `5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d` |
| frozen certificate | `.../phase6a1-controlled-bnv-debug/phase6a1-frozen-certificate-final/certificate.tsv` | `9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab` |
| frozen coefficients | same directory, `coefficients.tsv` | `3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6` |
| oracle root | `.../CompactStar-phase6a1-checkpoint-reconstruction-validation/build/phase6a1-checkpoint-reconstruction-validation/run/oracle` | tree `21f1ff9eaf23b078f86aa2b10ec45be55f84a98b9aedc4b94f4a7c4bf999f866` |
| oracle result | sibling build root, `oracle-result.json` | `f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0` |

Ellipses above expand only to
`/Users/keeper/Documents/CompactStar/worktrees`; execution manifests must store
the full absolute paths.

## 5. Acceptance classes fixed before results

### 5.1 BYTE_IDENTITY

- CompactStar scientific source, governed inputs, EOS/data, literature, eleven
  baselines, and the 33 protected Phase-5D paths.
- Vendored archive and governed dependency-header bytes before and after.
- C0/C1/T input bytes and all non-dependency configuration values.
- Every consumed candidate header for which the source-recovery record claims
  an exact header match.
- Normalized deterministic scientific output files after removal only of the
  explicitly listed provenance fields below.
- IEEE-754 binary64 encodings of all pairwise equality-bearing scalar fields;
  this is 0 ULP, not a numeric tolerance.

### 5.2 EXACT_SEMANTIC_IDENTITY

- pass/fail classifications, fixture identities, channel masks, phase/spin
  switches, units, row counts, labels, source/domain/revision identities,
  accepted/rejected step counts, observation source classes, and checkpoint
  qualification states;
- focused dependency-contract return values, container dimensions and values,
  contour counts/indexing/values, and exported numeric content;
- C0/C1 and C1/T scientific JSON/TSV payloads when serializers carry permitted
  non-scientific provenance differences.

### 5.3 GOVERNED_NUMERICAL_TOLERANCE

Only existing analytic/reference comparisons may use their already committed
test tolerances.  Those tests must pass unmodified.  A governed comparator's
tolerance cannot excuse any nonzero C0/C1 or C1/T difference in a field listed
above as exact.  No new or relaxed tolerance is introduced by this experiment.

### 5.4 PROVENANCE_DIFFERENCE_ALLOWED

The following may differ and must be reported rather than normalized away
silently: absolute source/build/output paths; dependency mode; dependency
archive/header/source identities; object/archive/executable hashes; compile and
link command paths; wall/CPU timing; temporary directory names; and the already
governed portable compiler-provenance exclusion in Phase-5C.  No physical,
numerical, solver, input, state, diagnostic, or scientific-output field is in
this class.

## 6. Exact build and test inventory

All CTest invocations use `--output-on-failure`, run serially for evidence
stability, and execute in each configuration's own build tree.

### 6.1 Dependency contracts

1. Configure-time fail-closed cases: no overrides on Darwin (default mode), all
   four valid overrides, each of four single/partial overrides, nonexistent
   paths, non-directory includes, non-regular archives, invalid archive
   content, and missing sentinel headers.  The negative cases must fail and
   logs must show no implicit discovery or fallback.
2. Existing `compactstar_library_smoke` and all tests below exercise compiled
   and linked Zaki APIs.
3. One directly necessary focused executable is fixed now:
   `tests/dependencies/adr0018_confind_contract.cpp`, CMake target and CTest
   name `adr0018_confind_contract`.  It uses a deterministic in-memory grid
   and a configuration-local temporary output directory and covers the
   CompactStar-consumed CONFIND contract:
   `ContourFinder`, grid/work-directory/contour configuration,
   `SetGridVals` overloads, `GetContourSet`, `ConvertToCurve2D`, `Cont2D`
   indexing/value access, `Plot`, and contour export.  It must also cover the
   Zaki types crossed by that API.  ROOT must not be linked unless an actually
   consumed call proves it unavoidable.  The exact same executable and inputs
   run under C1 and T; C0 runs an equivalent executable built from the same
   committed test source.

### 6.2 Phase-5B

Exact CTest names:

`phase5b_PB1`, `phase5b_PB6`, `phase5b_PB7`, `phase5b_PB9-11`,
`phase5b_PB12`, `phase5b_PB13`, `phase5b_contracts`, and
`phase5b_structural_response_regression`.

Exact driver command:

`ctest --test-dir <build> --output-on-failure -j 1 -R '^(phase5b_PB1|phase5b_PB6|phase5b_PB7|phase5b_PB9-11|phase5b_PB12|phase5b_PB13|phase5b_contracts|phase5b_structural_response_regression)$'`

The registered commands are those in `tests/CMakeLists.txt`: the staged tests
invoke `tests/analysis/run_phase5b_validation.py` with
`phase5b_freegas_validation`, the named stage, and the configuration-local
`tests/phase5b-evidence`; the regression invokes
`particle_number_response_regression.py` with
`produce_particle_number_reference.py`, `phase5b_freegas_validation`, the
configuration source root, governed baseline
`tests/baselines/phase5b_structural_response.json`, and a fresh local output
root.  Each configuration must pass the governed baseline independently; the
newly generated scientific payloads are then compared C0/C1 and C1/T under
the exact classes above.

### 6.3 Phase-5C

Exact CTest names:

`chemical_production_contract`, `chemical_exact_oracles`,
`chemical_curved_gc9`, `chemical_trackr_budget`,
`chemical_production_validation`, and
`phase5c_chemical_coefficient_regression`.

Exact driver command:

`ctest --test-dir <build> --output-on-failure -j 1 -R '^(chemical_production_contract|chemical_exact_oracles|chemical_curved_gc9|chemical_trackr_budget|chemical_production_validation|phase5c_chemical_coefficient_regression)$'`

The regression executes `chemical_coefficient_regression.py` with the
committed producer, the characterization and production executables, the
configuration source root, governed baseline
`tests/baselines/phase5c_chemical_coefficients.json`, reviewed candidate
`docs/validation/phase5c_chemical_coefficients_candidate.json`, and a fresh
configuration-local output root.  Governed qualification and C0/C1/T pairwise
identity are separate requirements.

### 6.4 Phase-5D

Exact CTest names:

`phase5d_response`, `phase5d_independent_oracles`,
`phase5d_component_tolerances`, `phase5d_protected_manifest`,
`phase5d_harness_controls`, `phase5d_coupled_oracles`, and
`phase5d1_controlled_evolution_regression`.

Exact driver command:

`ctest --test-dir <build> --output-on-failure -j 1 -R '^(phase5d_response|phase5d_independent_oracles|phase5d_component_tolerances|phase5d_protected_manifest|phase5d_harness_controls|phase5d_coupled_oracles|phase5d1_controlled_evolution_regression)$'`

Each configure supplies the same authenticated
`COMPACTSTAR_EOS_DATA_ROOT=/Users/keeper/Documents/CompactStar/data/compose`,
so the last test is registered.  The expensive regression uses the exact
registered producer `tests/rotochemical/fresh_context.py`, source root, EOS
root, governed baseline
`tests/baselines/phase5d1_controlled_evolution.json`, and a fresh local output
root.  No run card, baseline, or tolerance may change.  Governed qualification
and pairwise scientific identity are both required.

### 6.5 Bounded ADR-0017 Phase-6

First run CTest `passive_checkpoint_output_contract`.  Then run exactly one
`adr0017_production_qualification` main integration per configuration with the
authenticated matrix, profile, certificate, thermal source, entry manifest,
frozen certificate, coefficients, a fresh work root, fresh output root, and an
oracle-authentication flag produced by the unmodified verifier's
`authenticate` command.  The verifier `verify` command consumes that output
and the retained oracle root/result.

The exact contract driver is
`ctest --test-dir <build> --output-on-failure -j 1 -R '^passive_checkpoint_output_contract$'`.
The bounded harness positional command is
`<build>/tests/adr0017_production_qualification <matrix> <profile> <certificate> <thermal> <fresh-work> <entry> <frozen> <coefficients> <fresh-output> <oracle-qualified-flag>`.
The flag is created only after
`python3 tests/bnv/adr0017_production_verify.py authenticate --oracle-root <oracle-root> --oracle-result <oracle-result> --matrix <matrix> --profile <profile> --certificate <certificate> --thermal <thermal> --frozen <frozen> --coefficients <coefficients> --entry <entry> --output <authentication-json>`
passes; its exact contents are `ADR0017_ORACLE_AUTHENTICATED\n`.  Final
verification is
`python3 tests/bnv/adr0017_production_verify.py verify --output <fresh-output> --oracle-root <oracle-root> --oracle-result <oracle-result> --result <result-json>`.

The fixed fixture is `CPL-P2-LINEAR-QSS-v1`, `0` to `462269531250 s`,
`x0=eta_e0=eta_mu0=0`, P2, spin OFF, Me/Mmu ON, De/Dmu OFF; uninterrupted main
RKF45; `rtol=1e-11`; component `atol=(1e-16,1e-22,1e-22)`; 241-point passive
schedule; two isolated rk8pd reconstructions for every strict-interior point.

Each run independently must satisfy all self-qualification checks, every Cstar
knot observation, R18, and R20.  The accepted Mac authority is:

| Field | Exact requirement |
|---|---|
| final x | binary64 of `0.49240008824076903` |
| final eta_e | binary64 of `-2.512347425644221e-7` |
| final eta_mu | binary64 of `-4.7906773046561003e-7` |
| accepted / rejected | `232 / 60` |
| internal-step-history SHA-256 | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` |
| normalized R20 | binary64 of `5.2268230499136139e-7` (the serializer may print the equivalent last-digit spelling already governed by the verifier) |

The schedule, accepted states, brackets, internal history, trajectory, step
summary, every checkpoint state/diagnostic packet, reconstruction
classification, R18, and R20 are equality-bearing.  C0 must reproduce the
accepted authority.  C0/C1 and C1/T must be bit-identical in every scientific
state and diagnostic field and byte-identical in the normalized evidence.
Only timing, paths, dependency provenance, and binary identities may differ.
The retained 478-solve historical oracle is read-only comparison authority;
the production run performs the separately authorized bounded reconstruction
work and does not rerun BA12 or BA12R.

## 7. Three differing historical Zaki headers

The recovered CONFIND historical context exposed differing copies of
`Zaki/Util/Instrumentor.hpp`, `Zaki/Util/Logger.hpp`, and
`Zaki/Vector/DataSet.hpp`.  Their roles are not inferred from their names:

- `Instrumentor.hpp` controls compiled profiling/instrumentation macros and can
  affect generated calls and timing/provenance; no scientific difference is
  permitted.
- `Logger.hpp` controls emitted diagnostics and error paths and is runtime code,
  even when ordinary successful execution uses it only for messages; it cannot
  alter a scientific state or pass/fail decision.
- `DataSet.hpp` defines a container used in CompactStar scientific input,
  interpolation, integration, and serialized output paths and in the recovered
  CONFIND ABI.  It is scientifically load-bearing and is covered by dependency
  contracts plus Phase-5B/C/D and ADR-0017 exact treatment comparisons.

Current Zaki candidate headers consumed directly by canonical CompactStar are
authenticated separately.  The classifications above do not waive any ABI,
header, output, or numerical comparison.  Treatment qualification is the
authority.

## 8. Candidate build boundary

Zaki is built from an authenticated scratch copy or an external wrapper so its
Git worktree is not modified by generated configuration headers.  The complete
compiled source list, exact compiler/flags/includes, generated non-source
configuration, object hashes, archive members, archive hash, and link
dependencies are captured.  No consumed implementation may be omitted.

CONFIND is built by an external wrapper from the preserved read-only snapshot.
The required implementation set is predeclared as `Base.cpp`, `Cell.cpp`,
`Cont2D.cpp`, `Common.cpp`, and `ContourFinder.cpp`, subject to a pre-build
definition/reference audit that may add an existing exact source file but may
not edit one.  The build uses candidate Zaki headers/archive.  ROOT is excluded
unless the consumed structural contract proves a link requirement before the
treatment build.  Any need to change a scientific source or header is a stop,
not a patch opportunity.

## 9. Stop conditions

Execution stops immediately if a candidate, preserved source, transfer ZIP,
vendored archive/header, governed input, baseline, or protected-path identity
differs; the Mac platform materially differs; C0 fails current authority; C1
fails C0 equivalence; any override is partial, searches, or silently falls
back; a candidate needs scientific-source modification; archive provenance is
ambiguous; a focused contract or governed regression fails; any predeclared
exact C0/C1 or C1/T field differs; ADR-0017 self-qualification, R18, or R20
fails; or interpreting a result would require a post-result tolerance or class
change.

No result can ratify either candidate.  A complete pass establishes only:

**SOURCE-BUILD EQUIVALENCE QUALIFIED — READY FOR OWNER ACCEPTANCE AS LINUX
DEPENDENCY SOURCE AUTHORITY.**
