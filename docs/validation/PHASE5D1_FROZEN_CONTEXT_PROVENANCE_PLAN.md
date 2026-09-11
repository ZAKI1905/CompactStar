# Phase-5D-1P frozen-context provenance-preservation plan

**Disposition A: PHASE-5D PROVENANCE-PRESERVING FROZEN-CONTEXT ARCHITECTURE
AUTHORIZED — PHASE-5B/5C REGRESSION GATES RESTORED — READY TO REIMPLEMENT
COUPLING DOWNSTREAM.**

**Status:** ARCHITECTURE / PROVENANCE AUTHORIZATION; PRE-IMPLEMENTATION;
PRE-ODE; PRE-TRAJECTORY. This documentation-only record authorizes a later
downstream reimplementation. It is not an evolution or scientific candidate.

**Change class:** documentation. The future work described here will be a
bounded structural/architecture and numerical-adapter implementation under
accepted ADR-0014; this task changes no production or test source and changes
no scientific result. `GOVERNANCE.md:35-57` defines the authority order and
`GOVERNANCE.md:64-78` defines the relevant change classes.

## 1. Canonical SHA

Canonical local, `origin`, and live `master` were authenticated at
`d019ae390be4f5e3daba05039903485cb497e397`. The archived entry-ref snapshot
records the same local/origin canonical identity
(`build/phase5d-audit/provenance-gate-stop/identity/refs.json:1-8`).

## 2. Branch SHA and worktree

The authorized branch is
`physics/phase5d-controlled-rotochemical-evolution`, in
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5d-controlled-evolution`.
The committed entry branch SHA was
`74ddc16839fbea3675d5a9177eecf59f3bf3cb11`; local, upstream, and live branch
refs all matched before any restoration
(`build/phase5d-audit/provenance-gate-stop/identity/refs.json:2-4`).

## 3. Failure summary

The uncommitted coupled draft changed two equality-bearing governed Analysis
sources while attempting to reduce repeated full-file currentness reads during
RHS evaluation. The Phase-5B and Phase-5C regressions correctly rejected the
changed scientific source identities. Separate still-dirty reruns captured raw
rc `8` for each regression. The artifact comparison found no scientific-field
drift; that equality does not waive source identity and does not make the
failures false positives
(`build/phase5d-audit/provenance-gate-stop/manifest.json:14-16`).

The failed draft remained uncommitted, never entered the authorized trajectory,
and created neither an evolution SHA nor a candidate SHA. Its passing draft
oracles are evidence about useful behavior only, not candidate authority. The
accepted contract requires a changed dependency to refuse before scientific
access and requires result provenance to retain the semantic response,
coefficient, star/geometry, spin, ordering, solver, and initial-condition
identities (`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:422-429`).

## 4. Complete governed hash closures

The exact closure is the set of repository paths emitted as keys of each
producer's `provenance.source_sha256` and compared as equality-bearing artifact
provenance. All paths in the two sets below are class **A — equality-bearing
governed source identity**. The only class **B — portable execution
provenance** items are JSON fields, not repository paths:
`provenance.build.compiler` for Phase-5B and
`provenance.toolchain.compiler` for Phase-5C. No other portability exception is
created.

`PHASE5B_GOVERNED_SOURCE_HASH_SET` (19):

```text
CompactStar/Analysis/ParticleNumberResponse.hpp
CompactStar/Analysis/src/ParticleNumberResponse.cpp
CompactStar/AngularVelocity.hpp
CompactStar/Core/NStar.hpp
CompactStar/Core/RotationSolver.hpp
CompactStar/Core/StarProfile.hpp
CompactStar/Core/src/NStar.cpp
CompactStar/Core/src/RotationSolver.cpp
CompactStar/Core/src/TOVSolver.cpp
CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp
CompactStar/Geometry.hpp
CompactStar/RelativityUnits.hpp
CompactStar/Units.hpp
tests/analysis/particle_number_homogeneous.hpp
tests/analysis/phase5b_freegas_validation.cpp
tests/analysis/produce_particle_number_reference.py
tests/eos/structure1/local_oracle.hpp
tests/eos/structure1/table.hpp
tests/eos/structure1/tov_oracle.hpp
```

`PHASE5C_GOVERNED_SOURCE_HASH_SET` (17):

```text
CompactStar/Analysis/ChemicalResponse.hpp
CompactStar/Analysis/ParticleNumberResponse.hpp
CompactStar/Analysis/src/ChemicalResponse.cpp
CompactStar/Analysis/src/ParticleNumberResponse.cpp
CompactStar/EOS/LocalThermodynamics.hpp
CompactStar/EOS/TrackRFreeGasThermodynamics.hpp
CompactStar/EOS/src/LocalThermodynamics.cpp
CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp
docs/validation/phase5c2_preproduction_evidence.json
tests/analysis/chemical_curved_gc9.py
tests/analysis/chemical_exact_oracles.py
tests/analysis/chemical_production_contract.cpp
tests/analysis/chemical_production_evidence.py
tests/analysis/chemical_production_fixture.cpp
tests/analysis/chemical_production_validation.py
tests/analysis/chemical_trackr_budget.py
tests/analysis/chemical_trackr_fixture.cpp
```

The overlap is exactly
`CompactStar/Analysis/ParticleNumberResponse.hpp`,
`CompactStar/Analysis/src/ParticleNumberResponse.cpp`, and
`CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp`. Therefore
`UNION_GOVERNED_SOURCE_HASH_SET` contains 33 paths.

The following inspected paths are class **C — not repository inputs to these
`source_sha256` maps**:

```text
tests/analysis/particle_number_response_regression.py
tests/analysis/produce_chemical_coefficient_reference.py
tests/analysis/chemical_coefficient_regression.py
tests/baselines/phase5b_structural_response.json
tests/baselines/phase5c_chemical_coefficients.json
docs/validation/phase5c_chemical_coefficients_candidate.json
docs/validation/PHASE5B_INV09_GLOBAL_RESPONSE_INTEGRATION.md
docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_INTEGRATION.md
docs/validation/PHASE5B_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md
docs/validation/PHASE5C2_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md
docs/adr/ADR-0011-particle-number-structural-response.md
docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md
tests/CMakeLists.txt
```

These paths remain comparator/governance/artifact authorities and are not
implementation-editable merely because they are outside the source-hash set.
Generated EOS/profile/partition/particle-constant hashes are equality-bearing
artifact provenance but are generated values, not additional repository paths.

## 5. Exact protected-source paths and branch-entry hashes

The later fail-fast guard uses this complete sorted union and these exact
branch-entry bytes:

| SHA-256 | Protected repository path |
|---|---|
| `56111729f309e3066371fb6b70eb1803ee1fe8c729bebf2aeae3298a9102ddbe` | `CompactStar/Analysis/ChemicalResponse.hpp` |
| `df5b55a1ee679428052c0ee5797c1cdc7a61a4ad35270ea13d470c722957378e` | `CompactStar/Analysis/ParticleNumberResponse.hpp` |
| `985910361aa58a525abea9d6e1726138b856ca4743966bb08a5e0ee3ec4f3a49` | `CompactStar/Analysis/src/ChemicalResponse.cpp` |
| `0b03159c07c81ec35ada39f3fd5eabb410fe4e7c0c4cc46ab8133eb674c0f5fd` | `CompactStar/Analysis/src/ParticleNumberResponse.cpp` |
| `aa32cf06f172ca36ed5061c19db66f72dafbea9056807b3334c52d8ad178db79` | `CompactStar/AngularVelocity.hpp` |
| `7aa71d70496f359c48f6a26b8dc17424a98b4e30e6c5f5c0627010d340504e2b` | `CompactStar/Core/NStar.hpp` |
| `7d3dbfab42734c47a55c98877b855a37d14d37df71450966a8ecf514c557421d` | `CompactStar/Core/RotationSolver.hpp` |
| `07c2fdb783d03b938ba0b7de9e3fb2cd6ca62a6ad6b9c430b4101e29d25f8237` | `CompactStar/Core/StarProfile.hpp` |
| `69a4b403370feeb6c1ba7d469082c2a6f5d824bf604abbae268f505204feb08e` | `CompactStar/Core/src/NStar.cpp` |
| `b3e4b8db36952f5d1cb60412947a071b30a22de9a755f77571a580f3f8856982` | `CompactStar/Core/src/RotationSolver.cpp` |
| `d5da96e512453ae9ce18d7e6b3334b82bf17dd1def3230dd24ad267bc1fc7ad8` | `CompactStar/Core/src/TOVSolver.cpp` |
| `cd28ac9ad7ef6efd9f910d78991d242de71b0e1a03992718e7dca260a083c601` | `CompactStar/EOS/LocalThermodynamics.hpp` |
| `57f62a1e636cf3cfa23db9289a42302a484581c8c8f413fba9057cca6217c710` | `CompactStar/EOS/TrackRFreeGasThermodynamics.hpp` |
| `aa6fdb64094be2f53c1bfb52c91d5f102f9c17d4f37d08343ba096303a1fab6d` | `CompactStar/EOS/src/LocalThermodynamics.cpp` |
| `ae6c1c570e076b6462e0558592de39d8be922eeb5118cb40b0b296baba0800e8` | `CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp` |
| `cad68a4773f899006c037ee00116dbc2d01c4b239c06b8f6fac9d55f00e12167` | `CompactStar/Geometry.hpp` |
| `4277eb3caff0623f2615e8435a6632f954cf505a4be9aa9a7a1509b8488549fd` | `CompactStar/RelativityUnits.hpp` |
| `fbeb722e03d40ce7bf23680987e11a38da11566df292f8123d3aa16b0145c5b6` | `CompactStar/Units.hpp` |
| `883f9879103fde113b49ed311b38b5f4c4876d84ffd2aa7987270c5c08f4e481` | `docs/validation/phase5c2_preproduction_evidence.json` |
| `877d1e71800179f578dce3a36e113af49785c1e10204555fe67deb2b108595a4` | `tests/analysis/chemical_curved_gc9.py` |
| `c2d6fb699e1a322b9e4e93ce1fcbb76c60f0d5efb8fc04348988d8d43e50902e` | `tests/analysis/chemical_exact_oracles.py` |
| `1c9cdeff75e76fbe7a93c54e42764ac4878c7f4c48d13a9c05a5a1b9ebeb56c2` | `tests/analysis/chemical_production_contract.cpp` |
| `9124324d82fd4bc877fb99c35e7fe62a0d72f4fea9f03b8bf4fed61dde46e0ee` | `tests/analysis/chemical_production_evidence.py` |
| `c16e204cb48edc30e06332e6ae729156ed459fe66ce2bccd60cb460c676fb0e0` | `tests/analysis/chemical_production_fixture.cpp` |
| `6e4a7726744d1e601faf42ff9b3270710ff327972b901b1381f3cf03bbb90185` | `tests/analysis/chemical_production_validation.py` |
| `397923d762d6ff6687b9b2d0bb4d7bbee5abd97265426fc38d225db6aeb8e56f` | `tests/analysis/chemical_trackr_budget.py` |
| `6d725cce39fc05e0275673dbb36cffd27eeb0776b8745c7d62d5675778ec6978` | `tests/analysis/chemical_trackr_fixture.cpp` |
| `333258fcabd537a003fba7ccabfe7909139da850106c399431e0c3fc0fd2a243` | `tests/analysis/particle_number_homogeneous.hpp` |
| `2bb980f1916f0498e626fd8601ba023e6f211f93d5b8655b40aa62283c613fa4` | `tests/analysis/phase5b_freegas_validation.cpp` |
| `194f2000bdac9260c55a35135a38ce6a9864fb115a3e309f6887d69e56c207c9` | `tests/analysis/produce_particle_number_reference.py` |
| `d46d9ac1ef23df8754d32329e0df4497624d7da87d3f327704e0e0c81c2008ca` | `tests/eos/structure1/local_oracle.hpp` |
| `cf164e22626cda08050461776a1e87195538b6913c24f6ec47fe41e57ef4e160` | `tests/eos/structure1/table.hpp` |
| `d36df7a0f818b73227920bb50e6f8e5c113679328921d0c560de549c04b62686` | `tests/eos/structure1/tov_oracle.hpp` |

## 6. Exact failed-draft source mismatches

| Path | Failed-draft SHA-256 | Required entry SHA-256 | Closure impact |
|---|---|---|---|
| `CompactStar/Analysis/src/ParticleNumberResponse.cpp` | `af723087a773a1913021342e0336918ae98f0ccf7cf7cc112ded35f741aeccc6` | `0b03159c07c81ec35ada39f3fd5eabb410fe4e7c0c4cc46ab8133eb674c0f5fd` | Phase-5B and Phase-5C |
| `CompactStar/Analysis/src/ChemicalResponse.cpp` | `8d1ca9338a67e5261c0e348c192b44ca22f142220f07ae00e6bc73d1ae559ab8` | `985910361aa58a525abea9d6e1726138b856ca4743966bb08a5e0ee3ec4f3a49` | Phase-5C |

The archive independently records these exact pairs
(`build/phase5d-audit/provenance-gate-stop/evidence/exact-source-hash-mismatches.json:1-10`).
The other 31 union paths matched entry and governed artifact bytes.

## 7. Scientific payload equality observation

The dirty Phase-5B artifact differed from its baseline at the allowed compiler
field and the forbidden `ParticleNumberResponse.cpp` source hash. The dirty
Phase-5C artifact differed at the two forbidden source hashes. Removing only
those enumerated provenance fields produced exact parsed equality; scientific
payload drift was therefore **NO**. The baseline/candidate hash authorities
remain `7588f0e9...5fa`, `7027aa61...fe7`, and `a6b430b2...33b`, as recorded by
the prior qualification (`docs/validation/PHASE5D1_STRUCTURAL_RESOLUTION_QUALIFICATION.md:248-268`).

## 8. Why the provenance failure remains binding

Scientific-value equality does not make a changed producer source the ratified
producer. The regression contract intentionally binds scientific source
identity and permits only the already-ratified compiler field to vary. No
allowlist expansion, re-baseline, new compiler exception, or “science
unchanged” waiver is authorized. The regressions were correct to fail, and any
future need to edit a union path stops Phase-5D for a separate versioned
upstream extension.

## 9. Archive location and manifest

The complete failed draft is preserved at
`build/phase5d-audit/provenance-gate-stop/`. Its reconstruction authority is
`manifest.json`, authenticated by `manifest.json.sha256`. The archive contains
the full binary-capable tracked patch, the empty staged patch, exact snapshots
of every tracked and untracked draft path, separate dirty-tree Phase-5B and
Phase-5C logs/rc files, generated comparison evidence, coupled-oracle logs,
compiler/build identity, exact source mismatch records, status snapshots, and
branch/canonical refs. It records 12 tracked changes, zero staged paths, 16
untracked paths, tracked-patch SHA-256
`76eae7e209cb7567f45171520d456b40dc809ea35d8adb8c0bf15ee0f8316368`,
and empty staged-patch SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
(`build/phase5d-audit/provenance-gate-stop/manifest.json:2-16`).

The manifest SHA-256 is
`4815dc19157875de0b58de851824e4fadf72a17b59555f2de8056e971f91263c`.
All 112 manifest-listed archive files and all 91 archived-original/evidence
entries were rehashed successfully before restoration. Compiler identity was
Apple clang 21.0.0, target `arm64-apple-darwin25.6.0`
(`build/phase5d-audit/provenance-gate-stop/identity/cxx-version.txt:1-4`).
The scratch archive is ignored and is not committed scientific authority.

## 10. Failed-draft change classification

The two Analysis changes replaced character/stream-buffer reads with 64 KiB
bulk reads. Both are class 2, performance optimization, with secondary class 5,
implementation convenience. Neither is required scientific capability, a
lifetime repair, nor a provenance bypass. The committed implementations already
perform exact source-byte checks
(`CompactStar/Analysis/src/ParticleNumberResponse.cpp:72-79` and
`CompactStar/Analysis/src/ChemicalResponse.cpp:257-264`).

Every non-governed draft path was archived, assigned one of the required three
classifications, and removed from the active tree:

| Exact draft path | Classification | Required handling |
|---|---|---|
| `CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp` | REIMPLEMENT | Move controlled-source composition downstream; do not rewrite generic cooling. |
| `CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Cache.hpp` | REIMPLEMENT | Replace the injected generic cache payload with run-context ownership. |
| `CompactStar/Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp` | REIMPLEMENT | Compose same-Ltilde beta power inside the controlled context. |
| `CompactStar/Physics/Evolution/DriverContext.hpp` | REIMPLEMENT | Replace the raw spin pointer with downstream strong ownership; no generic edit. |
| `CompactStar/Physics/Evolution/EvolutionConfig.hpp` | REIMPLEMENT | Put the exact three-component tolerance contract in the downstream adapter. |
| `CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp` | REIMPLEMENT | Put telemetry/component control in the bounded downstream adapter. |
| `CompactStar/Physics/Evolution/Integrator/src/GSLIntegrator.cpp` | REIMPLEMENT | Use the downstream scaled-RKF45 adapter; do not alter the generic solver. |
| `CompactStar/Physics/Evolution/src/EvolutionConfig.cpp` | REIMPLEMENT | Log/validate the downstream adapter's frozen configuration instead. |
| `docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_IMPLEMENTATION.md` | REJECT | Reject the dirty addendum as authority; this plan is the sole addendum. |
| `tests/CMakeLists.txt` | REIMPLEMENT | Re-add only audited downstream tests in the resumed task. |
| `CompactStar/Physics/Driver/Thermal/ImmutableThermalSource.hpp` | REIMPLEMENT | Recreate as a Rotochemical-owned source-qualified immutable thermal object. |
| `CompactStar/Physics/Evolution/ISpinHistory.hpp` | REIMPLEMENT | Recreate the interface/implementation under Rotochemical with strong ownership. |
| `CompactStar/Physics/Rotochemical/FrozenChemicalCoefficients.hpp` | REIMPLEMENT | Retain the concept but qualify once and avoid per-RHS upstream accessors. |
| `CompactStar/Physics/Rotochemical/SecularEvolutionDriver.hpp` | REIMPLEMENT | Build the composite driver around the sealed context and pre-accumulation validation. |
| `docs/validation/PHASE5D1C_POST_REQUALIFICATION_STOP_REPORT.md` | REUSE CANDIDATE | Evidence content only; consume from archive and do not restore this draft file. |
| `docs/validation/phase5d1c_stop_evidence.json` | REUSE CANDIDATE | Evidence content only; consume from archive and do not restore this draft file. |
| `docs/validation/phase5d1c_upstream_regression_failure.txt` | REUSE CANDIDATE | Raw failure content only; consume from archive and do not restore this draft file. |
| `tests/rotochemical/component_tolerances.cpp` | REIMPLEMENT | Test the downstream scaled adapter and exact flat ordering. |
| `tests/rotochemical/coupled.hpp` | REUSE CANDIDATE | Audit and rerun the useful numerical design before credit. |
| `tests/rotochemical/coupled_oracles.hpp` | REUSE CANDIDATE | Audit and rerun; newest disconnected-support control has no pass credit. |
| `tests/rotochemical/evolution.cpp` | REIMPLEMENT | Rebuild around the new frozen context; do not restore the draft executable. |
| `tests/rotochemical/fixture.hpp` | REIMPLEMENT | Construct exclusively through public semantic APIs and retained owners. |
| `tests/rotochemical/injected_faults.hpp` | REUSE CANDIDATE | Audit and rerun each fault against the eight guards. |
| `tests/rotochemical/qualified_suite.py` | REUSE CANDIDATE | Audit path/identity assumptions and rerun against the new executable. |
| `tests/rotochemical/trajectory.hpp` | REIMPLEMENT | Rebuild only in the later task; no trajectory is authorized here. |
| `tests/rotochemical/validate_trajectory.py` | REUSE CANDIDATE | Retain validation design only; every trajectory check is still unexecuted. |

The newest disconnected-support control compiled but did not run, and
trajectory validation never ran.

Nothing from the archive may be copied blindly. Its intended behavior must be
reimplemented inside the boundary in section 26.

## 11. Public-API feasibility result

**YES; no essential upstream public API is missing.** Existing public APIs can
construct, validate, retain, and read the qualified semantic objects without a
union-path edit:

- `FixedBaryonNumberResponse::Compute`, inherited `RequireCurrent`,
  `WholeStarIPhysical`, and `WholeStarEquilibriumNumberRate`
  (`CompactStar/Analysis/ParticleNumberResponse.hpp:84-96` and
  `CompactStar/Analysis/ParticleNumberResponse.hpp:135-151`);
- `GlobalChemicalNumberResponse::Compute`, `Values`, `NumericalError`,
  `Support`, `Partition`, `Diagnostics`, `Lifetime`, and `RequireCurrent`
  (`CompactStar/Analysis/ChemicalResponse.hpp:130-185`);
- `ChemicalImbalanceResponse::Compute`, `Values`, `Channels`, `PaperZ`,
  `Global`, and `RequireCurrent`
  (`CompactStar/Analysis/ChemicalResponse.hpp:194-253`);
- `RotochemicalSpinDrive::Compute`, `Values`, `IPhysical`, `Evaluate`, and
  `RequireCurrent` (`CompactStar/Analysis/ChemicalResponse.hpp:255-310`);
- `GlobalUrcaChannelCoefficient::Compute`, `RequireCurrent`,
  `LuminosityCoefficient`, `Entry`, `Selection`, `ChemicalDomain`,
  `DomainIdentity`, `MetricIdentity`, `PartitionKm`, and `Order`
  (`CompactStar/Physics/Rotochemical/GlobalUrcaChannelCoefficient.hpp:25-51`);
- `RotochemicalReactionResponse` and `RotochemicalThermalPower::From`, which
  already keep the Ltilde/F/H/rate/heating ledger on one coefficient authority
  (`CompactStar/Physics/Rotochemical/RotochemicalReactionResponse.hpp:7-55`).
- thermal construction/reference surfaces:
  `EOS::CompOSE_Thermo(directory,Options)`, `IsLoaded`, const
  `TGrid_MeV`/`NbGrid_fm3`/`YqGrid`, and
  `CvDensity_cgs_ForCooling`
  (`CompactStar/EOS/CompOSE_Thermo.hpp:78-182`);
  `Evolution::StarContext::HeatCapacityStar_Tinf`
  (`CompactStar/Physics/Evolution/StarContext.hpp:135-154`);
  `Driver::Thermal::PhotonCooling(Options)` and
  `Detail::ComputeDerived`
  (`CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp:273-321` and
  `CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp:94-103`);
  and `Driver::Thermal::NeutrinoCooling(Options)` with
  `NeutrinoCooling_Details::ComputeDerived`
  (`CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp:141-205` and
  `CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Details.hpp:104-170`).
  The resumed implementation snapshots/adapts these surfaces downstream; it
  does not use mutable generic caches as the controlled same-Ltilde authority.

`ChemicalLifetime` owns the stars, provider, revision, and source paths needed
for currentness (`CompactStar/Analysis/ChemicalResponse.hpp:47-64`). Global
construction snapshots those identities/bytes, and its currentness method
checks revision, owner coverage, profile pointer/version, and source bytes
(`CompactStar/Analysis/src/ChemicalResponse.cpp:434-480` and
`CompactStar/Analysis/src/ChemicalResponse.cpp:662-674`). W construction checks
that every structural source is covered by that lifetime and retains both Z and
the structural response (`CompactStar/Analysis/src/ChemicalResponse.cpp:744-772`).

No scratch proof program is necessary: the public ownership chain is explicit,
and the archived draft already constructed accepted radial-80000 Z, W, and
Ltilde values before failing the source-identity regression. This conclusion
does not validate the archived coupled implementation.

## 12. Proposed frozen-context ownership

The resumed implementation shall introduce
`CompactStar::Physics::Rotochemical::FrozenRotochemicalRunContext` (name may be
adjusted only to repository convention). It is immutable after successful
construction and privately owns or retains:

- the qualified central/background identity and all required star owners;
- `shared_ptr<const GlobalChemicalNumberResponse>`;
- `shared_ptr<const ChemicalImbalanceResponse>` and a typed private Z snapshot;
- `shared_ptr<const FixedBaryonNumberResponse>`;
- `RotochemicalSpinDrive` and typed private W/I snapshots;
- `shared_ptr<const GlobalUrcaChannelCoefficient>`, its controlled selection,
  support, domain/metric/partition identities, and private Ltilde snapshots;
- a frozen downstream reaction evaluator using those same Ltilde snapshots;
- private `StarContext`/`GeometryCache` authorities;
- an immutable source-qualified thermal table, heat-capacity source, photon
  envelope source, and separately represented non-controlled neutrino source;
- exactly one `shared_ptr<const PrescribedSpinHistory>` for this benchmark;
- normalization classification and exact source/path/hash identities;
- the radial-80000/EOS-8192/rho-c qualification certificate and metadata; and
- the frozen flat state/tolerance/order contract.

Private snapshots are not bare provenance-free arrays: the same sealed object
retains the accepted semantic authorities, owners, validation metadata, and
hashes that authorized every copied value. No mutable semantic handle is
exposed.

## 13. Construction-time dependency flow

Before production context construction, the candidate-validation harness shall
verify the 33-path entry manifest, the unchanged baseline hashes, and both
governed regressions. It then supplies the authenticated qualification
certificate/source-manifest identities to the context; the production context
does not open a Phase-5B/5C baseline. The expensive context
construction/qualification phase shall:

1. record and bind the authenticated qualification certificate and source
   manifest identities supplied by that harness;
2. construct/retain the radial-80000 star, profile, EOS/provider, geometry,
   fixed-baryon structural response, and complete chemical lifetime;
3. construct `GlobalChemicalNumberResponse`, then
   `ChemicalImbalanceResponse`, then `RotochemicalSpinDrive` through their
   public factories and unchanged goals;
4. construct one `GlobalUrcaChannelCoefficient` from the same chemical-domain
   owner and the predeclared controlled normalizations;
5. construct and validate the heat-capacity, photon, non-controlled-neutrino,
   and immutable thermal-table authorities;
6. construct the sole prescribed spin-history owner;
7. call `RequireCurrent()` on Global, Z, fixed structural response, W, and
   Ltilde, and validate every source/lifetime/version token;
8. require exact `{Npe,NpMu}` semantic order, 2-by-2 Z, two W/I entries,
   `{Me,Mmu}` enabled, `{De,Dmu}` disabled, matching domain/partition/metric,
   declared units, and shared owners;
9. require radial 80000, EOS 8192, `rho_c=1.10e15 g cm^-3`, and exact
   profile/model/certificate/source hashes; and
10. copy the already-accepted typed Z/W/I/Ltilde values into private immutable
    fields and seal the context.

The revised predeclaration fixes radial 80000 and freezes every other benchmark
input (`docs/validation/PHASE5D1_CONTROLLED_EVOLUTION_PREDECLARATION_REVISION.md:9-50`).

## 14. Runtime dependency flow

Each cheap RHS evaluation shall:

1. validate the cheap live tokens described in section 15 before scientific
   accumulation;
2. decode exactly `(ln(Tinf/1e8 K), eta_npe, eta_npmu)` and validate finite
   state/time;
3. sample `Omega` and `OmegaDot` once from the sole spin owner;
4. evaluate F/H, rates, equilibrium/incremental neutrino power, and chemical
   heating from the context's same-Ltilde frozen authority;
5. evaluate photon, non-controlled neutrino, and heat-capacity terms from the
   owned immutable thermal authorities;
6. form both chemical rows with the full semantic Z matrix and frozen W;
7. validate the complete local result, signs, units, and finiteness; and only
8. then commit one complete update to the RHS accumulator.

RHS evaluation must not rebuild radial-80000 G/Z/W/Ltilde, reread baseline JSON,
reread Analysis source files, reconstruct upstream internals, or mutate any
owned authority. The accepted dependency direction keeps rotochemical physics
downstream of Analysis and keeps the generic evolution core independent of
rotochemical physics (`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:396-420`).

## 15. Currentness policy

Full upstream currentness is mandatory at context construction and again
immediately before run start: Global, Z, fixed response, W, and Ltilde all call
their existing `RequireCurrent()` contracts. During evolution, the sealed
context uses immutable run-owned snapshots and exposes no mutator. Before each
RHS it performs cheap checks only on dependencies still capable of external
mutation: central/profile pointer and `StarProfile::Version`, geometry
provenance, `ChemicalRevision` liveness/serialization, provider metadata,
lifetime-owner identities, and the spin-history generation/currentness token.
Any mismatch refuses before a local result reaches the accumulator.

The policy does not weaken or bypass existing factories. It closes one v1 run
over values accepted while current, retains the source semantic objects, and
prevents later mutation from changing the values being evolved. The upstream
contracts already make profile version and exact source identity explicit
(`CompactStar/Analysis/src/ParticleNumberResponse.cpp:247-257` and
`CompactStar/Analysis/src/ChemicalResponse.cpp:810-834`).

ADR-0014 literally says that a read-only spin history is supplied through
`DriverContext` (`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:383-391`).
This plan satisfies the ownership and registration-order invariant downstream
by having the one composite controlled driver own the frozen context and its
one spin history. If owner review later requires a literal member on generic
`Evolution::DriverContext`, that one-file seam requires separate explicit
authorization; it is not authorized here.

## 16. Lifetime policy

The context holds strong shared ownership of every object required by the
semantic lifetime chain. `ChemicalLifetime` keeps all contributing stars,
provider, revision token, and source paths alive. Z retains Global; W retains Z
and the fixed response; Ltilde retains the same Global domain. The context also
holds the central star, thermal data, geometry, photon/heat-capacity sources,
and spin history directly. Destruction of caller handles cannot expire a run
dependency. Foreign owners, changed coverage, changed profile versions, or an
expired revision/spin token refuse.

## 17. No-baseline-runtime policy

Runtime reading of `tests/baselines/phase5b_structural_response.json`,
`tests/baselines/phase5c_chemical_coefficients.json`, or the reviewed candidate
is forbidden. Baselines remain regression comparators only. Production Z and W
come from their unchanged public semantic factories, and Ltilde comes from
`GlobalUrcaChannelCoefficient::Compute`. The prior qualification likewise used
the governed baseline only after primary construction passed
(`docs/validation/PHASE5D1_STRUCTURAL_RESOLUTION_QUALIFICATION.md:248-254`).

## 18. No-upstream-source-edit policy

Every path in `UNION_GOVERNED_SOURCE_HASH_SET` is immutable for Phase-5D v1.
No source buffering, allowlist expansion, baseline rewrite, or upstream
convenience edit is permitted. A scientifically necessary change to any union
path is a hard stop and requires a separately versioned/revalidated upstream
task.

## 19. Thermal source-byte ownership

A Phase-5D-owned `FrozenThermalSource` under `Physics/Rotochemical/` shall own
the exact paths, bytes, SHA-256 values, parsed/interpolated identity, and
lifetime of `eos.t`, `eos.nb`, `eos.yq`, and `eos.thermo`. Construction reads
and hashes the files, constructs an internal
`shared_ptr<const EOS::CompOSE_Thermo>`, then rereads/compares the source bytes
to reject construction races. The sealed parsed table and bytes have no mutable
alias. Runtime uses that immutable representation and never an anonymous
callback over ephemeral buffers.

This preserves the already-declared requirement to own a const thermal table
and input bytes (`docs/validation/PHASE5D1_CONTROLLED_EVOLUTION_PREDECLARATION.md:28-34`).
No file-buffering change in `ParticleNumberResponse.cpp` or
`ChemicalResponse.cpp` is authorized.

## 20. Spin owner

The context owns one immutable prescribed-dipole spin-history object. Both
`Omega(t)` and `OmegaDot(t)` are sampled once from it per RHS and the same
sample feeds every chemical row and any controlled spin consumer. The
rotochemical module does not rederive a torque. Foreign equal-valued instances,
expired callback captures, token changes, or registration-order dependence
refuse. The unchanged spin law is recorded in the predeclaration revision
(`docs/validation/PHASE5D1_CONTROLLED_EVOLUTION_PREDECLARATION_REVISION.md:33-37`).

## 21. Semantic channel validation

Construction uses `ChemicalImbalanceResponse::Channels()` and
`PaperZ(row,column)`, not positional copying. It requires exact semantic order
`{Npe,NpMu}`, a 2-by-2 matrix, both cross terms, two W/I entries in the same
order, and common Global/domain ownership. Runtime accesses typed channels and
must not infer them from a flat index. The contract fixes this state ordering
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:69-91`).

## 22. Empty-channel behavior

Heat capacity, photon cooling, and non-controlled neutrino sources are owned
and validated independently of controlled process count. An empty controlled
selection returns exactly zero controlled reaction rate, neutrino increment,
and chemical heating while still evaluating finite positive `Cstar`, photon
cooling, and any explicitly enabled non-controlled source. There is no division
by a channel count or channel luminosity. Invalid/missing heat capacity refuses;
an empty controlled selection does not.

## 23. Pre-RHS dependency checks

The controlled driver builds a complete local evaluation before invoking any
accumulator update. Stale spin, star/profile, geometry, chemical domain,
provider/revision, owned source-identity/certificate token, thermal ownership,
channel identity, or a nonfinite result throws with the accumulator unchanged.
Full Analysis/source-byte currentness runs at construction and immediately
before run start; the RHS does not reread source files. No partial thermal or
chemical contribution may survive a dependency failure.

## 24. Radial qualification provenance

The context records and requires radial resolution 80000, EOS resolution 8192,
central density `1.10e15 g cm^-3`, exact profile/model/certificate hashes,
domain/partition/metric identities, units, channel order, and normalization
classification. The qualified r80000 W factory accepted without bypass, in
Npe/NpMu order and MeV s-squared units
(`docs/validation/PHASE5D1_STRUCTURAL_RESOLUTION_QUALIFICATION.md:204-227`).
Runtime diagnostics emit these frozen identities and never rebuild them.

Radial 80000 remains qualified and frozen. The failed coupled draft did not
invalidate the revised predeclaration; it invalidated only the approach that
modified governed upstream source bytes. No new radial predeclaration is
created.

## 25. Eight architecture guards

| Guard | Owner | Construction-time invariant | Runtime invariant | Failure mode | Required test |
|---|---|---|---|---|---|
| 1. Callback currentness | Frozen context and its `shared_ptr<const SpinHistory>` | Nonempty identity, live lifetime/generation token, immediate `RequireCurrent` | Validate the cheap token before sampling and before accumulation | Throw with untouched accumulator | Expire/mutate the token and test ephemeral-capture refusal |
| 2. Shared spin owner | Frozen context | Every controlled consumer receives the pointer-identical one owner | Sample once; share the exact Omega/OmegaDot sample independent of driver order | Reject foreign/stale owner | Foreign equal-valued owner, stale same owner, and registration-order permutation |
| 3. Semantic Z channels | Context retaining Z/W plus typed snapshots | Exact `{Npe,NpMu}`, 2-by-2 Z, two W/I entries, both cross terms, common Global/domain owner | Typed access only; both cross terms used | Refuse wrong order/shape/owner | Constructor-refusal tests for wrong semantic order, wrong shape, and foreign Z/W owner through a test builder; separately swap flat layout and drop each cross-Z term as arithmetic mutations |
| 4. Empty-channel thermal denominator | Context thermal ledger | Heat capacity/photon/non-controlled sources validate independently of process selection | Empty controlled set gives zero controlled terms but retains cooling and positive finite Cstar | Bad thermal/Cstar refuses; emptiness does not | Empty controlled set: Cstar positive, heating zero, cooling derivative negative |
| 5. Dependencies before RHS | Context plus composite controlled driver | All dependencies validated before driver registration | Check cheap live tokens, then complete and validate the local result before the first `AddTo`; no per-RHS source-file read | Throw with bit-unchanged sentinel accumulator | One fault per stale spin/domain/source-identity token, foreign star/geometry, and thermal owner |
| 6. Radial/provenance diagnostics | Context qualification record | Exact radial-80000/EOS-8192/rho-c/profile/model/certificate/domain/order/unit identities | Emit frozen identities; no G/Z/W/Ltilde rebuild | Refuse construction or pre-RHS check | Alter radial count, profile/certificate hash, partition endpoint, and foreign matched geometry |
| 7. Thermal source-byte boundary | `FrozenThermalSource` | Own paths, exact bytes/hashes, const parsed table; pre/post parse equality | Read only the owned immutable representation; no disk/source callback | Missing/read/hash/parse/race refusal | Mutate during construction; mutate after construction/before run and require pre-run refusal; mutate after run start and prove sealed values unchanged without per-RHS reread; destroy caller handles and prove retained lifetime |
| 8. Upstream governed-source immutability | Entry hash manifest and candidate-validation harness | Hash all 33 paths exactly; Phase-5B/5C PASS before candidate | Repeat before first trajectory; no union writer | Immediate STOP; no candidate/trajectory | Negative scratch-manifest perturbation and positive exact-restoration regression run |

## 26. Intended editable-file boundary

The next implementation may edit only:

```text
CompactStar/Physics/Rotochemical/**
CompactStar/Physics/CMakeLists.txt
tests/rotochemical/**
tests/CMakeLists.txt
docs/validation/PHASE5D1_*.md          (new implementation/validation/status records only)
docs/validation/phase5d1_*.json        (new implementation evidence only)
```

New production context, driver, spin-history, frozen-thermal, and scaled-RKF45
adapter files must live under `CompactStar/Physics/Rotochemical/`. Every path in
the 33-path union is forbidden. Existing Phase-5B/5C baselines/candidates and
benchmark inputs are also read-only. Existing Phase-5D predeclaration and
qualification records may be cited but not revised by the implementation task.

No edit to generic `EvolutionSystem`, generic `GSLIntegrator`, generic
`EvolutionConfig`, generic `DriverContext`, or the thermal core is authorized
by this record. If downstream implementation proves that a narrow generic seam
is unavoidable, it must stop and obtain separate explicit authorization with
dedicated tests; a broad solver or thermal redesign is not an admissible seam.

## 27. Component-tolerance plan

The frozen flat state remains
`(ln(Tinf/1e8 K), eta_npe, eta_npmu)`. A new Phase-5D-owned scaled-RKF45 adapter
under `Physics/Rotochemical/` shall use GSL's scaled control with baseline
`rtol=1e-7`, absolute scale `(1e-12,1e-18,1e-18)`, and refinement
`rtol=1e-9`, absolute scale `(1e-14,1e-20,1e-20)`. It forwards to public
`EvolutionSystem::operator()` and notification methods, validates exact
three-component layout/positive finite tolerances, transports callback
failures, and owns its stiffness/step telemetry.

The exact GSL scaled-control convention is `eps_abs=1`, `eps_rel=rtol`,
`a_y=1`, `a_dydt=0`, and `scale_abs=component_atol`. Thus each listed absolute
scale is the actual per-component absolute tolerance, not an additional
multiplier with an unspecified convention.

Current generic configuration exposes only scalar `atol`, and the generic
integrator allocates the scalar-control GSL driver
(`CompactStar/Physics/Evolution/EvolutionConfig.hpp:178-190` and
`CompactStar/Physics/Evolution/Integrator/src/GSLIntegrator.cpp:475-483`).
Therefore a downstream adapter is required, but a generic integrator edit is
**not** required. A later optional generic `component_atol` feature would be a
separate authorized numerical-method change with shape/finite/positive and
legacy-parity tests. No solver redesign is authorized. The tolerances and their
refinement are unchanged from the revised predeclaration
(`docs/validation/PHASE5D1_CONTROLLED_EVOLUTION_PREDECLARATION_REVISION.md:38-44`).

## 28. Same-Ltilde plan

One context-owned `GlobalUrcaChannelCoefficient` is the authority for each
controlled channel. At construction its exact Me and Mmu Ltilde values and
support/identity are copied into the sealed downstream evaluator. Those same
values feed equilibrium beta cooling, F correction, H reaction rate, and
chemical heating. The controlled driver computes the entire beta ledger from
that authority; it does not inject into or modify generic `NeutrinoCooling`.

For the controlled run only, historical placeholder DU/MU cooling is not
registered as a separate source. Photon cooling and an explicitly represented
non-controlled neutrino source remain separate; PBF is disabled for this v1
fixture. This prevents double counting without an upstream cooling rewrite.
The same-Ltilde requirement is normative
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:440-454`).

## 29. Previous coupled-oracle evidence

The archived dirty draft retained the following useful, non-candidate evidence:

| Draft check | Preserved result |
|---|---|
| spin-only analytic maximum relative error | `1.1160340242261544e-16` |
| reaction-only maximum relative error | `2.8436935598444193e-08` |
| coupled RE9 same-Ltilde ledger | PASS |
| active/dead-channel Lyapunov checks | PASS |
| dependency guards exercised by that run | PASS |
| qualified coupled-oracle raw rc | `0` |

The oracle log, executable/input/CMake hashes, and exact build identity are
preserved in the archive manifest around
`build/phase5d-audit/provenance-gate-stop/manifest.json:286-371`. The future
implementation must rerun all of these after provenance-preserving
reimplementation. The newest disconnected-support test was compiled but not
run; no trajectory evidence exists.

## 30. Regressions restored to PASS

After the archive was verified, every draft path was explicitly restored or
removed, and all 33 union paths were hashed equal to branch-entry/governed
bytes. The relevant producer executables were rebuilt from the restored tree.

| Restored gate | Result | Raw rc | Generated-versus-baseline behavior |
|---|---|---:|---|
| `phase5b_structural_response_regression` | 1/1 PASS | 0 | Scientific/equality fields equal; only ratified compiler path differs. Generated SHA `05fe5e87689df6a84312db1a15f1b768e613373a9c18bf4d31da9579da96bee6`; baseline SHA `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`. |
| `phase5c_chemical_coefficient_regression` | 1/1 PASS | 0 | Generated artifact is byte-identical to baseline, SHA `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`. |

The raw PASS records are
`build/phase5d1p-audit/phase5b-post-restore.log:1-14` and
`build/phase5d1p-audit/phase5c-post-restore.log:1-14`. The comparison evidence
records the Phase-5B compiler-only portable path and controls
(`build/tests/phase5b-regression-evidence/regression-ko10vadi/comparison-evidence.json:1-24`)
and exact Phase-5C equality
(`build/tests/phase5c-regression-evidence/regression-lufue5dc/comparison-evidence.json:1-45`).
No baseline or reviewed candidate changed.

## 31. Exact resumed-implementation gates

The next task must pass these gates in order:

1. authenticate this plan SHA, branch/worktree, clean entry, canonical master,
   predeclaration, and exact benchmark inputs;
2. hash all 33 union paths against section 5 before any implementation and
   after every source edit batch;
3. restrict writes to section 26 and stop on any required upstream/generic
   exception;
4. construct only through unchanged public semantic APIs; call all upstream
   `RequireCurrent()` contracts and close all eight guards;
5. prove typed Z/W/I/Ltilde, lifetime, source bytes, channel order, units,
   radial-80000, same-Ltilde, spin-owner, empty-channel, and pre-accumulation
   dependency behavior;
6. rerun the retained spin-only, reaction-only, RE9, Lyapunov, dependency,
   disconnected-support, component-tolerance, source-byte, and mutation
   oracles;
7. hash the 33 paths again and require both governed regressions PASS before a
   candidate commit;
8. immediately before the first trajectory, repeat all 33 hashes and both
   governed regressions; and only then
9. execute the unchanged predeclared trajectory and its stiffness,
   convergence, thermal-sign, and quasi-steady diagnostics.

Any mismatch stops without trajectory or candidate. The hash guard is an
earlier fail-fast requirement and does not replace either regression.

## 32. No trajectory in this task

No ODE was entered, no benchmark trajectory was generated, no coupling was
restored after the documentation commit, and no `PHASE5D1_EVOLUTION_SHA` or
`PHASE5D1_CANDIDATE_SHA` exists. This task authorizes architecture only.

## 33. Global INV-11

Global INV-11 remains **UNRESOLVED**. Contract-level INV-11b-e status and open
INV-11f implementation/validation status are unchanged
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:532-545`). This
plan, its future implementation, passing regressions, or a later candidate does
not by itself constitute global invariant closure or scientific ratification.

## 34. A18

A18 was not implemented or begun. Realistic FR2005/A18 normalization remains
outside this controlled mathematical architecture benchmark. The accepted ADR
excludes A18 implementation (`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:481-486`).

## 35. BNV

BNV was not begun. No BNV source, rate, heating, code, or claim is authorized;
the accepted contract explicitly forbids beginning BNV before standard
rotochemical evolution is validated
(`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:467-477` and
`docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md:588`).

## Remaining caveats

- The downstream frozen-context design is authorized, not implemented or
  independently reviewed.
- The literal ADR-0014 `Evolution::DriverContext` transport wording is satisfied
  semantically by a context-owned composite controlled driver; a literal
  generic member, if later required by owner review, is outside this boundary.
- Exact component tolerances require the authorized downstream scaled-RKF45
  adapter; no generic integrator edit is authorized.
- Archived draft oracle evidence must be rerun; its newest disconnected-support
  check and every trajectory check remain unexecuted.
- No realistic source normalization or global INV-11 closure follows from this
  authorization.

## One recommended next action

Reimplement the Phase-5D coupled secular evolution from
PHASE5D1_FROZEN_CONTEXT_PLAN_SHA strictly within the authorized downstream
file boundary. Construct a provenance-preserving immutable frozen run context
from unchanged governed Phase-5B/5C public semantic APIs, rerun all previously
passing spin-only/reaction-only/RE9/Lyapunov/dependency oracles, verify the
Phase-5B/5C governed regressions before the first trajectory, then execute the
unchanged radial-80000 predeclared controlled trajectory with stiffness,
convergence, thermal-sign and quasi-steady diagnostics. Do not modify any
governed upstream source-hash path, do not change benchmark inputs, do not
claim realistic FR2005 normalization, do not implement A18, and do not begin
BNV.

That next action is not begun by this record.
