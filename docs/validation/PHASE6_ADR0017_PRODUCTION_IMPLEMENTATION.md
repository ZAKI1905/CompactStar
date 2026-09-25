# Phase-6 ADR-0017 production implementation and bounded qualification

## Status and authority

**Status:** PRE-RUN DECLARATION — PRODUCTION IMPLEMENTATION CANDIDATE.

This record implements and boundedly qualifies the accepted, human-ratified,
canonical ADR-0017 architecture.  It does not authorize a Phase-6 numerical
candidate, a physical BNV rate/model, a BA12/BA12R rerun, or canonical
integration.

| Authority | Identity |
|---|---|
| canonical entry, local `master`, `origin/master`, live `master` | `ffd597e00162fa4efc477fe0d23878d7ad05ff6c` |
| historical passive-observation authority | `9af8912107ea77a2e2ea51517c1776f26a4b7b49` |
| segmentation diagnostic authority | `6057eb92339e5a0596baab6a652c6290d0930658` |
| reconstruction evidence source | `ae0d5607fe8889cec70604747dfe3bfea7e5fd97` |
| implementation branch | `physics/phase6-adr0017-production-implementation` |
| implementation worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6-adr0017-production-implementation` |

The branch and worktree were absent at entry and were created from the exact
canonical entry.  Canonical `master` must remain at the entry SHA throughout.

## Authenticated local platform

The qualification platform is the same local numerical platform as the
accepted Arm-E and reconstruction evidence:

| Item | Authenticated value |
|---|---|
| machine | local Mac only; no cluster |
| macOS | `26.6.2` build `25G83` |
| kernel | Darwin `25.6.0` |
| architecture | `arm64` |
| compiler | Apple clang `21.0.0` (`clang-2100.3.34.2`) |
| compiler target | `arm64-apple-darwin25.6.0` |
| CMake | `4.2.1` |
| GSL | `2.7.1` |
| build type | `Debug` |
| configured C++ flags | empty; Debug flags `-g` |
| effective C++ flags | `-g -std=c++17 -arch arm64 -fPIC -pthread -Wall -Wextra -Xclang -fopenmp` |
| floating-point flags | no `-ffast-math`, finite-math-only, unsafe-math, explicit FMA, or contraction-changing flag |

A material platform or effective-floating-point-flag mismatch is a stop before
qualification integration.  The cluster may not be substituted.

## Exact production implementation map

The recovered Phase-6 BNV source architecture is source evidence, not
canonical ancestry.  Production code remains Phase-6-owned.

| Responsibility | Exact production owner |
|---|---|
| Phase-6 physics/RHS and diagnostic context | `CompactStar/Physics/BNV/FrozenControlledBnvRunContext.hpp`, `CompactStar/Physics/BNV/src/FrozenControlledBnvRunContext.cpp` |
| controlled BNV secular RHS | `CompactStar/Physics/BNV/ControlledBnvSecularDriver.hpp`, `CompactStar/Physics/BNV/src/ControlledBnvSecularDriver.cpp` |
| immutable accepted-step record and uninterrupted main RKF45 integration | `CompactStar/Physics/BNV/PassiveCheckpointOutput.hpp`: `AcceptedStepRecord`, `UninterruptedBnvTrajectory`; `CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp` |
| passive observation schedule and bracketing only | same files: `PassiveObservationSchedule`, `ObservationBracket` |
| isolated strict-interior rk8pd O1/O2 and exact-endpoint shortcut | same files: `Rk8pdCheckpointReconstructor` |
| reconstruction self-qualification and fail-closed status | same files: `ReconstructionQualification`, `CheckpointStatus` |
| diagnostic evaluation from the selected checkpoint state | same files: `CheckpointOutput`; callback into the existing `FrozenControlledBnvRunContext` diagnostic owner |
| reconstruction provenance and uncertainty | same files: `CheckpointProvenance`, per-component `d_O`, `U_O`, and `F_i` |
| unchanged scientific-grid composite-trapezoid R20 | same files: `ComputeCheckpointR20` |
| production build registration | `CompactStar/Physics/CMakeLists.txt`, `CompactStar/Physics/BNV/CMakeLists.txt` |
| architecture documentation | `docs/architecture/CURRENT_ARCHITECTURE.md` |
| P1--P10 contract tests | `tests/bnv/passive_checkpoint_output_contract.cpp` |
| bounded actual-fixture qualification harness and verifier | `tests/bnv/adr0017_production_qualification.cpp`, `tests/bnv/adr0017_production_verify.py` |

The scheduler owns requested times and bracket metadata only.  It does not call
the RHS, reconstruct state, evaluate diagnostics, or mutate the main
trajectory.  Each strict-interior level receives a new disposable evaluation
context and new GSL step/control/evolve objects; immutable fixture/source
inputs may be shared, mutable contexts and caches may not.  The exact accepted
main endpoint bypasses all reconstruction and RHS replay.  No change to
`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp`, any Cstar interpolation or
cache, or any governed Phase-5 path is permitted.

## Entry protected-byte manifest

### Eleven governed baselines

| Repository path | Entry SHA-256 |
|---|---|
| `tests/baselines/baryon_number_dscmf1_reference.tsv` | `90d607519cbdf3c4a0bf6ef50cc8fd22a8526b5db0354dc319e96854da29041d` |
| `tests/baselines/grid_convergence_cmf_1p6_debug.tsv` | `b48519c3e948e9979a385d19facee2777d15955eeb8711b4bdd46b81fef74741` |
| `tests/baselines/grid_convergence_cmf_1p6_trajectory.tsv` | `d5b753932c0523e67a7f25b460c7494bec1a006a8d01c9e43124cb2e78f0720f` |
| `tests/baselines/hartle_I_dscmf1_debug.tsv` | `034ecddbd9bd847650429d7dc87d0331ec9e87aca3862ff87594e4bff5b707dd` |
| `tests/baselines/hartle_monopole_dscmf1_debug.tsv` | `caaa0ac0d3219cda0a9fb518b27688afc23c6cdad1ec76a2bcd7359614a8d4e8` |
| `tests/baselines/passive_cooling_cmf_1p6_debug.tsv` | `8fef2314673fceb939f859612f4befe94117115d6d6b3ad0dcc59d1faa68c9f9` |
| `tests/baselines/phase5b_structural_response.json` | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` |
| `tests/baselines/phase5c_chemical_coefficients.json` | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| `tests/baselines/phase5d1_controlled_evolution.json` | `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd` |
| `tests/baselines/tov_dscmf1_reference.tsv` | `3d9af9129a6a4ffde9e0f8c5507a160f968a861c0cf9f3b089cceecab86b701a` |
| `tests/baselines/tov_path_equivalence_dscmf1.tsv` | `5c0f4b3bdb70921f8f2a869af10edc4d8f5ae3963a9d150e11ef859d21e1c678` |

### Thirty-three protected Phase-5D paths

| Entry SHA-256 | Protected path |
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

### Authorities, imported evidence, and scientific inputs

| Artifact | Entry SHA-256 |
|---|---|
| `docs/adr/ADR-0011-particle-number-structural-response.md` | `e6e8b27c3a4212a824e96f64155d039bbbd888336d0ad48983c255dbea3f0ae3` |
| `docs/adr/ADR-0013-corrected-rotochemical-chemical-coefficients.md` | `f9ece5ddd8dc14d267efde95b47d9113f5296e620e567fb953c9ed6257e8fd70` |
| `docs/adr/ADR-0014-secular-rotochemical-evolution-contract.md` | `bc38a49cd16b87ae72eeeb943dad20b192e3d291e4127cc0c3c68aeef6e067a6` |
| `docs/adr/ADR-0015-bnv-open-system-thermal-ledger.md` | `795c8bc851644de00c8c2ae37ea40c13c5a061bee691ef855d78fcd2096e1da5` |
| `docs/adr/ADR-0016-phase6-bnv-tangent-adapter-ownership.md` | `b2b0bba07ebb9d926473e0123e3b28a777656368a4d642f8f5172c85711412d5` |
| `docs/adr/ADR-0017-phase6-passive-checkpoint-reconstruction.md` | `9c05fac09dec7a025d28550415fdbd6a95c0cfb6605b10cff1c033b4ac2c27a9` |
| `docs/validation/PHASE5D1_CONTROLLED_ROTOCHEMICAL_EVOLUTION_INTEGRATION.md` | `69b32ba4b7d1e7de1cb134cb35be95f5f8b2c0d40235245577dbd594cd038291` |
| `docs/validation/PHASE5D1_GOVERNED_ARTIFACT_PREPARATION_PLAN.md` | `a272b3678521c89554ae85ba7e01c2ba0828b71102c658ee00b04f8623c3cac2` |
| `docs/validation/PHASE6A1_CONTROLLED_BNV_RECOVERY_IMPLEMENTATION.md` | `8ab91e364abc4eb843ecfdfbb03155d4ac2c18e20b4a4caf47d31dd01c4a2a3d` |
| `docs/validation/PHASE6A1_BA12R_NUMERICAL_FORENSICS.md` | `daccc117237c6789c4f3fabc0d6cb0252fc88ce020d67eb3b8bba4eb33577cea` |
| `docs/validation/PHASE6A1_BA12R_SEGMENTATION_DIAGNOSTIC.md` | `ab7f6b9f3a6d37a7c2c439f9b9c0a2a330fecc208bd801ae890fe89b0851e157` |
| `docs/validation/PHASE6A1_PASSIVE_OBSERVATION_EXPERIMENT.md` | `684a0338b7a08591ac266d8c8bd7924d95ac86c3488ab115637b4f2b7d7f76ee` |
| `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md` | `188861202879f7cb7b3bc7f106cbf28113881a2296779a20f58538bd80e72778` |
| `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_VALIDATION.md` | `59093f30ecfc0543f73f60ba416bc23a80d2ad753a9051230fdfc5646ee552a2` |
| `docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_CANDIDATE_RECOVERY.md` | `403090d6ada8c796b5b65c45d430a82fc55ce331f11bb1d9d4b34f607c4a05f7` |
| qualification certificate | `7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453` |
| thermal source tree | `1cdb98507958b50fbd1e1eff3eb3bfa4f48c790fef86a8740855884fc43b9d94` |
| frozen certificate | `9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab` |
| frozen coefficients | `3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6` |
| profile root tree | `233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72` |
| profile `profile.tsv` | `e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890` |
| profile `model.txt` | `3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4` |
| profile `freegas.tsv` | `7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a` |
| literature manifest | all entries PASS under `shasum -a 256 -c literature/SHA256SUMS.txt` |

All listed bytes are immutable.  Final protection rehashes them and refuses
any mismatch.

## Historical Arm-E and oracle authority

| Artifact | Authenticated SHA-256 |
|---|---|
| exact 241-point passive schedule | `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e` |
| 232 accepted endpoint states | `7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459` |
| 241 observation brackets | `0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d` |
| 232 internal accepted steps | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` |
| passive trajectory | `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c` |
| passive step summary | `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb` |
| solve matrix | `32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616` |
| 478-solve oracle root tree | `21f1ff9eaf23b078f86aa2b10ec45be55f84a98b9aedc4b94f4a7c4bf999f866` |
| Oracle-1 batch tree / ledger | `89bc3cb0bde7e58f60e1e44fdee9332e920ae163870213fe79d7bb81ba65ecf9` / `0fec6717e153ce19ecafda0fe0ed5c01d1be3934342ac9776f05ad2da769fbcc` |
| Oracle-2 batch tree / ledger | `f0bf0fb4e1e0acef58da1d770c0a71d87c4021264017f75cbfa8aa9da1a0af2d` / `b1bc3feb6bce96b7e4f1feecf2660b1f88781941d66d25ae832e4103a3248f74` |
| historical oracle result | `f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0` |
| historical harness executable | `82eea47fba62f1440a1b87abec95dc4ac6287756ad44b990564794b9cf703044` |

The historical 478 rk8pd solves are read-only numerical authority and are not
rerun to create comparison evidence.

## Bounded qualification fixture and exact expected outputs

Exactly one main integration is authorized:

- `CPL-P2-LINEAR-QSS-v1`, source ON, P2, spin OFF, Me/Mmu ON, De/Dmu OFF;
- `t=0` to `462269531250 s`;
- initial `(x,eta_e,eta_mu)=(0,0,0)`;
- one uninterrupted adaptive RKF45 state with `rtol=1e-11`,
  `atol=(1e-16,1e-22,1e-22)`, initial `h=1 s`;
- exactly one distinct positive GSL `t1`, `462269531250 s`;
- passive schedule hash `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e`.

The main path must reproduce Arm E exactly:

```text
x      = 0.49240008824076903
eta_e  = -2.5123474256442210e-7
eta_mu = -4.7906773046561003e-7
accepted steps = 232
rejected steps = 60
internal accepted-step SHA-256 = fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8
```

Any difference stops before reconstruction.  Schedule-invariance is tested by
feeding a second schedule to the retained accepted-step bytes only; it may not
run a second main ODE trajectory.

The one exact accepted-endpoint observation must return the exact stored main
state with zero rk8pd calls, zero RHS replay, and no numerical tolerance.  The
remaining 239 strict-interior observations are the complete reconstruction
set.  The three predeclared one-knot observations are the solve-matrix category
B rows; no special method or tolerance is allowed there.

## Frozen reconstruction hierarchy and self-qualification

Each strict-interior observation starts both levels at the exact `(t_L,y_L)`
and stops at the exact `t_obs`:

| Level | Method | `rtol` | `atol(x,eta_e,eta_mu)` |
|---|---|---:|---|
| O1 witness | GSL `rk8pd` | `1e-12` | `(1e-17,1e-23,1e-23)` |
| O2 reported | GSL `rk8pd` | `1e-13` | `(1e-18,1e-24,1e-24)` |

With exact max-digits-10 binary64 round trips (`Q_O=Q_i=0`), for every state
component:

```text
d_O  = abs(y_O2-y_O1)
M_O  = max(abs(y_O1),abs(y_O2))
D_O1 = atol_O1 + rtol_O1 M_O
D_O2 = atol_O2 + rtol_O2 M_O
F_O  = max(D_O2,64 ulp(M_O),Q_O)
M_i  = max(abs(y_L),abs(y_R),abs(y_O1),abs(y_O2))
D_U  = atol_ULTRA + rtol_ULTRA M_i
F_i  = max(D_U,64 ulp(M_i),Q_i)
U_O  = 2 max(d_O,F_O)
```

The checkpoint is qualified only when `d_O<=D_O1` and
`U_O<=0.20 F_i` for every component.  Otherwise its status is
`NUMERICALLY_UNRESOLVED`, O2 is not reported as a qualified checkpoint, and
downstream diagnostic/R20 qualification fails closed.  There is no fallback,
retry, tolerance change, or knot exception.

## Predeclared comparison budgets and gates

Because production uses the same platform, GSL algorithm, exact left state,
initial-step rule, isolated Phase-6 context semantics, and max-digits-10
binary64 data, the predeclared production-versus-historical oracle budget is
**bit-identical binary64 state and diagnostic values** for every O1 and O2
result.  File serialization, field ordering, timing, and newly added production
provenance fields are excluded from byte identity, but parsed scientific
binary64 values, identities, currentness, validity, and knot classifications
must be exact.  No post-result relaxation is allowed.

R20 retains the exact 241-point scientific grid and existing composite
trapezoid definition.  Production O2 checkpoint diagnostics, end-state terms,
`R20`, and `N_R20` must be bit-identical to the historical Oracle-2 evaluation.
Reconstruction uncertainty is additionally computed and reported separately;
it is not folded into or used to redesign R20.  The inherited subsidiary
comparison bound remains `5e-6 N_R20` and the governed closure gate remains
`abs(R20)/N_R20<=2e-4`.

Qualification gates, in order:

1. P1--P10 focused unit/contract and mapped negative/mutation tests pass before
   any actual-fixture ODE qualification run.
2. The single bounded main trajectory exactly reproduces Arm E; passive
   scheduling uses only the terminal `t1` and leaves main history unchanged.
3. Historical oracle bytes and manifests authenticate before reconstruction.
4. All 239 O1/O2 pairs self-qualify; all three one-knot observations pass the
   same uniform rule; production scientific values match the oracle exactly.
5. Reconstructed diagnostics match exactly and R20 preserves its existing
   semantics and budgets.
6. Final protected-byte, platform, branch, remote, and canonical checks pass.

P1--P10 cover passive authority, endpoint bypass, exactly O1/O2, pass,
fail-closed synthetic failure, no fallback/retuning, context isolation,
diagnostic-state ownership, R20 checkpoint consumption, and schedule
ordering/duplicates/missing refusal.  Any failed gate stops the later stage.

## Production output schema

The retained output exposes main integration provenance, passive schedule
identity, `MAIN_ENDPOINT` versus `RK8PD_RECONSTRUCTED`, left/right accepted-step
identities, O1/O2 configurations, per-component `d_O`, `U_O`, `F_i`, status,
self-qualified flag, selected reconstructed state, the complete existing
Phase-6 diagnostic packet, and a separate reconstruction-error contribution.
Main and reconstructed state semantics may not be mixed silently.

## Explicit exclusions and run accounting

This task authorizes one bounded source main trajectory, 239 O1 local solves,
and 239 O2 local solves after their gates.  The local reconstructions may use
conservative process parallelism with independent contexts, GSL objects, and
output roots; a scientific solve itself is not multithreaded.

It explicitly excludes fresh BASELINE/REFINED/ULTRA full P2 trajectories,
matched controls, BA12, BA12R, the future six-trajectory source/control
campaign, the full 77-test campaign, Phase-5D or shared-solver changes, Cstar
changes, governed baseline changes, any physical BNV rate/model selection, a
BNV candidate, A18, superfluidity, Regime-II/MixedStar, cluster use, automatic
merge, or any automatic next phase.

## Results append-only boundary

The pre-run content above is frozen by the commit named
`docs: predeclare adr0017 production qualification`.  Results below this
boundary are appended in later commits; the pre-run declaration is not
amended.

---

## Qualification results

**Final branch status:** PRODUCTION IMPLEMENTATION CANDIDATE — bounded
qualification PASS; awaiting explicit owner review and acceptance.  This is
not a governed/canonical implementation and does not create a Phase-6
numerical candidate or BNV baseline.

### Immutable implementation commits

| Purpose | Commit |
|---|---|
| pre-run declaration | `0c25df48a871403e716f1926390d2fec41c15c05` |
| Phase-6 production implementation | `dee330df5b8231ff55beb4adb75824272056db21` |
| focused tests and bounded qualification machinery | `a0679f01e2ffe259a4a03081721a892b2cf8bee4` |

After qualification, the human owner explicitly authorized correction of the
implementation-map path from `docs/CURRENT_ARCHITECTURE.md` to
`docs/architecture/CURRENT_ARCHITECTURE.md`.  The correction is typographical
and path-only.  It changes no equation, numerical requirement, implementation
behavior, run card, tolerance, qualification gate, expected hash, scope
exclusion, evidence, or result.  The historical predeclaration commit remains
unchanged; the correction occurs only in its owner-authorized acceptance
descendant.

### Stop-gate execution

| Gate | Result |
|---|---|
| exact local platform/toolchain | PASS |
| P1--P10 focused contract suite before actual-fixture integration | PASS, 10/10 |
| historical oracle authentication | PASS: 478 existing results and both ledgers authenticated; zero oracle reruns |
| bounded main Arm-E equivalence | PASS before reconstruction |
| passive scheduling invariance | PASS by two schedules replayed over the same retained history; no second main ODE |
| exact positive endpoint shortcut | PASS |
| 239 strict-interior O1/O2 pairs | PASS, zero unresolved |
| reconstructed diagnostics and R20 | PASS |
| final protected-byte checks | PASS |

Two verifier-preflight attempts stopped before numerical integration: the
first resolved the solve matrix to its authenticated historical worktree path;
the second bound the hash audit to the preserved ledger's actual tier-specific
file naming.  Both used fresh output roots.  They performed no ODE or rk8pd
solve.  The final authentication root passed before the sole bounded main run.

### Production owners and isolation audit

`UninterruptedBnvTrajectory` owns one persistent RKF45 step/control/evolve
state, persistent `h`, immutable accepted-step records, and exactly one
positive terminal `t1`.  `PassiveObservationSchedule` owns only ordered
requested times and accepted-step brackets.  It neither calls the RHS nor
reconstructs/evaluates a checkpoint.  `Rk8pdCheckpointReconstructor` owns the
exact-endpoint shortcut, two fresh strict-interior rk8pd solves, the ratified
self-qualification rule, fail-closed output status, post-selection diagnostic
evaluation, and separate reconstruction uncertainty.

The shared frozen Phase-6 physics fixture is immutable in coupled mode.  Every
O1, O2, and diagnostic evaluation receives a fresh wrapper with its own
`RunState`, `EvolutionSystem`, controlled driver, and GSL objects.  The only
mutable reference-only cache in `FrozenControlledBnvRunContext` is not used by
the coupled path.  The main trajectory, its integrator, and all other
checkpoint wrappers remain inaccessible to a reconstruction worker.  P7 also
mutates an isolated synthetic context and confirms no cross-context or main
state effect.  Actual-fixture accounting created 958 distinct wrappers.

No shared Phase-5 solver or governed Phase-5 code was changed.  In particular,
`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp` and all Cstar interpolation
and cache owners are byte-identical to entry.

### Main integration and passive scheduling

Exactly one new bounded source main integration ran.  It used
`CPL-P2-LINEAR-QSS-v1`, source ON, P2, spin OFF, Me/Mmu ON, De/Dmu OFF, and the
predeclared ULTRA RKF45 settings from `0` to `462269531250 s`.

| Evidence | Production result | Arm-E comparison |
|---|---:|---|
| final `x` | `0.49240008824076903` | exact binary64 |
| final `eta_e` | `-2.5123474256442210e-7` | exact binary64 |
| final `eta_mu` | `-4.7906773046561003e-7` | exact binary64 |
| accepted / rejected | `232 / 60` | exact |
| distinct positive GSL `t1` targets | `1` (`462269531250 s`) | exact |
| schedule SHA-256 | `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e` | exact |
| accepted-state SHA-256 | `7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459` | exact |
| bracket SHA-256 | `0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d` | exact |
| internal-step SHA-256 | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` | exact |
| trajectory SHA-256 | `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c` | exact |
| step-summary SHA-256 | `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb` | exact |

The second schedule was pure postprocessing of those retained accepted-step
bytes.  It did not invoke the main RHS or GSL integration and therefore could
not change main history.

### Endpoint, reconstruction, diagnostics, and knots

The initial observation uses the stored initial main state.  The sole positive
observation coincident with an accepted endpoint is observation 240 at
`462269531250 s`; it returned the exact final main state above with zero rk8pd
invocations, zero RHS replay, and no numerical tolerance.

All 239 strict-interior observations launched exactly O1 and O2.  All passed
the uniform componentwise requirements `d_O<=D_O1` and
`2 max(d_O,F_O)<=0.20 F_i`; there were zero unresolved checkpoints, no retry,
no fallback, and no tolerance change.  The maxima were both at observation
117, component `x`:

- maximum `d_O/D_O1 = 0.5794839113173227`;
- maximum `U_O/(0.20 F_i) = 0.5793505315921852`.

The three one-knot observations 82, 117, and 228 all self-qualified under the
same path and tolerances.  No knot-special code or tolerance exists.

Parsed production O1 and O2 binary64 values are bit-identical to all 478
authenticated historical oracle results.  Reconstructed diagnostic values are
also bit-identical at all 241 checkpoints for `P_dir_eq`, `P_dir_actual`,
`L_H`, `DeltaLnu`, `DeltaPbeta`, `Lnu_eq`, `Lnu_full`, `Lgamma`, `Lother`,
`Pnet`, `mu_B`, `mu_n_actual`, both `sigma` components, `Echem`, Cstar,
temperature, and baryon count.  Source/domain/revision/partition/product-fate
identities are exact.  Diagnostics were evaluated from O1/O2 checkpoint states
through the existing Phase-6 owner; no diagnostic was independently
interpolated.

### R20 and reconstruction uncertainty

The production result retained the exact 241-point grid, diagnostic-state
inputs, composite-trapezoid quadrature, normalizer, and existing thresholds:

| Quantity | Result |
|---|---:|
| `R20_residual_erg` | `5.7186274768377777e39` |
| `N_R20_erg` | `1.0940924194731047e46` |
| `R20_normalized` | `5.2268230499136139e-7` |
| separate reconstruction uncertainty | `7.0988433612780846e31 erg` |

R20, its normalizer, and its normalized value are bit-identical to the
Oracle-2 evaluation.  The uncertainty is reported separately and is below the
predeclared `5e-6 N_R20` subsidiary bound.  R20 semantics and thresholds were
not changed; adaptive-step quadrature was not introduced.

### Performance evidence

Actual reconstruction concurrency was one process.  No scientific local solve
was multithreaded.

| Stage | Wall (s) | CPU (s) |
|---|---:|---:|
| shared frozen physics-fixture construction | `264.32115375000001` | `264.22510299999999` |
| one bounded main RKF45 integration | `23.151261792` | `23.147572` |
| 239 O1 solves | `0.79458454400000011` | `0.79248000000000041` |
| 239 O2 solves | `0.83425795600000019` | `0.83415400000000051` |
| 958 disposable wrapper constructions | `0.03927129900000001` | `0.039194000000000034` |
| reconstructed diagnostic evaluations | `0.14473092599999998` | `0.14456300000000008` |
| total checkpoint-output overhead | `1.8128447250000002` | `1.8103910000000010` |

Linear checkpoint-only extrapolation to 8192 observations is
`62.129753734205025 s` serial.  Simple ideal two-process scaling is
`31.064876867102512 s`; it was not executed and is not claimed as a measured
or guaranteed parallel runtime.  The one-time fixture construction and main
integration are separate from those checkpoint-only estimates.

Retained untracked build evidence has SHA-256
`a182f07affea4ff38aff827ac7000d3ae831a3f8dcb4f827fffc51b07c815171`
for `checkpoints.tsv` and
`1274ad3a2cdb868b47611d9333ffb249a1f1da082bb051abd42f0dbcb6a46fa5`
for `performance.tsv`.

### Final protection and scope accounting

Final rehashing passed for all 11 governed baselines, all 33 protected
Phase-5D paths, 15 tracked Phase-5B/C/D/ADR/imported-evidence authorities, the
authenticated EOS/profile/thermal/frozen inputs, both oracle trees and
ledgers, and every literature-manifest entry.  Changed-path inspection found
no Phase-5D, `ScaledRKF45`, Cstar, EOS/data, baseline, or literature change.

The focused ADR-0017 contract test passed after implementation.  The full
77-test campaign was not run.  BA12 and BA12R were not rerun; no full
BASELINE/REFINED/ULTRA trajectory or matched control ran; the future six-run
campaign did not run.  No BNV candidate or baseline was created, no physical
BNV rate/model was selected, and A18, superfluidity, Regime-II/MixedStar, and
cluster work were not begun.

### Disposition

**ADR-0017 PRODUCTION IMPLEMENTATION QUALIFICATION PASS — UNINTERRUPTED MAIN
INTEGRATION, PASSIVE SCHEDULING AND SELF-QUALIFIED rk8pd CHECKPOINT OUTPUT MATCH
VALIDATED NUMERICAL AUTHORITY — READY FOR OWNER ACCEPTANCE.**

Recommended next action: return this noncanonical branch to the owner for
explicit acceptance.  Do not run the six clean source/control trajectories.
Only after owner acceptance and canonical integration should a separate
EKU-cluster production-qualification task be considered; only after that
separate qualification may the owner be asked to authorize the six-run clean
BA12R campaign.

## Postqualification owner acceptance

**Status:** QUALIFICATION PASS / OWNER-ACCEPTED / CANONICAL INTEGRATION
AUTHORIZED.

On 2026-09-25, the human owner explicitly accepted the qualified ADR-0017
production implementation at
`47a24317c604d73a7a8eece720823231036c765f`, including production commit
`dee330df5b8231ff55beb4adb75824272056db21` and test/validation commits
`a0679f01e2ffe259a4a03081721a892b2cf8bee4` and
`47a24317c604d73a7a8eece720823231036c765f`.  The detailed acceptance record is
`docs/validation/PHASE6_ADR0017_PRODUCTION_ACCEPTANCE.md`.

This acceptance authorizes fast-forward-only canonical integration of the
qualified implementation plus the documentation-only acceptance descendant.
It does not create a Phase-6 numerical BNV candidate or governed BNV baseline,
authorize a clean BA12R campaign, qualify the EKU cluster, select a physical
BNV rate/model, or authorize any additional numerical execution.  Historical
BA12 and BA12R remain permanently **FAIL**.  The exact next task after
successful canonical integration is a fresh documentation/planning task for
EKU-cluster numerical-platform qualification; it may not submit cluster jobs.
