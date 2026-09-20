# Phase-6A-1 checkpoint reconstruction candidate recovery

Classification: **PHASE-6 CHECKPOINT-RECONSTRUCTION RECOVERY VALIDATION
EVIDENCE; NOT BNV CANDIDATE; NOT GOVERNED BASELINE; NOT PHYSICAL RESULT**.

Historical status is immutable and is not rewritten by this recovery:

- historical validation attempt at
  `8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead`: **Disposition E — SIDE-EFFECT /
  CURRENTNESS / PROVENANCE FAILURE**;
- BA12: **FAIL**;
- BA12R: **FAIL**;
- passive observation scheduling: **PASS**;
- oracle self-qualification: **PASS**.

## Immutable pre-run recovery declaration

Status: **PREDECLARED / NO NEW LOCAL INTEGRATION EXECUTED**. This section is
frozen before any of the 956 authorized replay integrations. It must not be
amended after the predeclaration commit; execution evidence is appended in a
later commit.

### Authority and authenticated entry

The human owner explicitly accepted the checkpoint-reconstruction validation
failure at `8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead` and authorized a fresh
bounded candidate-only recovery that reuses the 478 authenticated rk8pd oracle
integrations, corrects only the candidate context binding, and executes exactly
the remaining 956 replay integrations. The accepted reconstruction preflight
`93e93c7f91a3cd8fced2f7a0961eda9c469c43fe` still controls the solve matrix,
candidate definitions, oracle hierarchy, budgets, knot categories, tie-break,
hybrid eligibility and stop conditions
(`docs/validation/PHASE6A1_CHECKPOINT_RECONSTRUCTION_PREFLIGHT.md:223-455`,
`:580-609`). None of those is changed here.

The recovery branch is
`analysis/phase6a1-checkpoint-reconstruction-candidate-recovery` in worktree
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-checkpoint-reconstruction-candidate-recovery`,
created from exact failed-validation SHA
`8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead`. At entry the failed branch was
clean with local/upstream/live refs equal to that SHA; canonical local
`master`, `origin/master` and live `refs/heads/master` were clean and equal to
`bd697ffdc474863d7a39f42e17ad8e8dbf105e5d`. The failed validation history and
its evidence are preserved; nothing is squashed or rewritten.

The owner does not authorize, and this recovery will not perform: any new
oracle integration; any main trajectory; a BA12/BA12R rerun; a tolerance,
observation, method, budget, category or tie-break change; any `ScaledRKF45`,
Cstar, Phase-5D, production, baseline, EOS/data or literature change; a
physical BNV rate/model; production adoption; a candidate artifact; cluster
execution; or a merge.

### Exact root cause of the historical failure

The historical candidate phase failed because the harness received a
trajectory TSV where the EOS/profile directory was required. The traced data
flow at `8b783dbe73cc504b8aa00a5e9aeb48677bcf1ead` is:

| Step | Location | Finding |
| --- | --- | --- |
| launcher | shell invocation of the candidate phase (not archived) | 12 unlabeled positional arguments; slot 3 received a trajectory TSV path instead of the profile/EOS directory |
| parameter | `tests/bnv/checkpoint_reconstruction_validate.cpp:315-319` (`ParseInputs`) | `argv[3]` bound positionally to `Inputs::profile` with no type or content validation |
| consumer 1 | `tests/bnv/checkpoint_reconstruction_validate.cpp:215` (`BuildContext`) | `Fixture(inputs.profile, ...)` |
| consumer 2 | `tests/bnv/checkpoint_reconstruction_validate.cpp:219` | `Campaign::Tangent(fixture, inputs.profile, ...)` |
| consumer 3 | `tests/bnv/checkpoint_reconstruction_validate.cpp:222` | `Campaign::Qualification(inputs.profile, ...)` |
| expected semantic type | `tests/rotochemical/fixture.hpp:44-70`; `tests/bnv/campaign_fixture.hpp:19-40` | a directory containing `freegas.tsv`, `model.txt`, `profile.tsv` |
| failure point | `tests/rotochemical/fixture.hpp:49-50` | `source(profile_dir/"freegas.tsv")` then `solve(table_path, ...)` imported `<trajectory.tsv>/freegas.tsv`; the EOS import failed inside the child process before any observation result, leaving empty `owning-star/central` scratch directories |
| authoritative bytes | `CompactStar/Physics/Rotochemical/FrozenRotochemicalRunContext.hpp:89-92` | production qualification requires `profile.tsv`, `model.txt`, `freegas.tsv` with SHA-256 `e9cd03b0…`, `3ea70de7…`, `7cd44c92…` |

The oracle phase of the same attempt used the correct directory (its contexts
constructed and qualified), so the defect was confined to the candidate launch
binding. The historical `.error` file was not written because the EOS import
terminated the child outside the harness exception path; the harness reported
the method-batch failure as non-retriable, and no candidate integration ran.

### Unique authenticated EOS/profile binding

The corrected binding is derived from existing authenticated provenance, not
from directory names:

| Authority | Recorded identity | Resolved path / identity | Hash / provenance |
| --- | --- | --- | --- |
| production qualification literals | `profile.tsv`, `model.txt`, `freegas.tsv` byte hashes | `FrozenRotochemicalRunContext.hpp:89-92` | `e9cd03b0b8449806f6c9883d75de1d3dff0cf1481675efc45a56519655d40890`, `3ea70de79e15b70c5a6d68f48335d18047ff80e60b55a9acdb78084e9be4d6d4`, `7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a` |
| passive-trajectory pre-run proof | `scientific_input_hashes` entry for the profile directory | `build/phase6a1-passive-observation-proof/preflight.json` (passive worktree) | tree SHA-256 `233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72` |
| governed Phase-5D candidate/provenance record | the same three file hashes for another instance of `t8192-r80000` | `docs/validation/phase5d1_controlled_evolution_candidate.json:10832-10843` | identical file hashes |
| Phase-6 entry manifest | entry-manifest identity consumed by `Qualification` | `.../docs/validation/phase6a1_controlled_bnv_entry_hashes.json` | `5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d` |
| reused oracle contexts | contexts built from the same directory qualified with these literals | failed-attempt `run/oracle/oracle1`, `oracle2` (478 successful solves) | identity `domain_identity=0:0:0:whole star through P=0; ...`, `revision_identity=phase6a1-uniform-proper-neutron-sink-v1` |

Resolved unique binding (directory, tree SHA-256
`233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72`):

```text
/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-controlled-bnv-implementation/build/phase6a1-controlled-bnv-debug/phase6a1-fresh-qualification/chemical-characterization/run-dhp6aivf/t8192-r80000
```

Its `profile.tsv`, `model.txt` and `freegas.tsv` hashes equal the production
qualification literals exactly. No other plausible binding remains.

### Harness binding correction (validation-harness code only)

`tests/bnv/checkpoint_reconstruction_validate.cpp` replaces the 12-positional
interface with named arguments. `--profile-root` is validated before any
context construction: it must be a directory (a regular file such as a
trajectory TSV is refused by name), must contain the three files, and each file
must reproduce the production qualification hash. There is no fallback search,
no working-directory discovery, no "first existing path" and no environment
variable. The harness never consumes a trajectory file; the solve matrix is
its only endpoint authority. Additional harness-only changes: single-method
`stage` execution in one fresh child process (concurrency 1, no retry), a
no-solve `dry-run` mode, an explicit `--oracle-qualified-flag` path, and an
atomic global execution ledger keyed by the frozen 956 solve IDs. `Integrate`,
`EvaluateMethod`, `BuildContext`, `ReadMatrix`, `WriteMeta`, the diagnostic
`Save` path, the replay/oracle configurations, the initial-step rule and the
apply-call cap are byte-identical to the failed-validation harness.
`tests/bnv/checkpoint_reconstruction_verify.py` is unchanged; its `oracle` and
`final` subcommands remain the sole scientific adjudicators. Production source,
`ScaledRKF45`, Cstar, Phase-5D, baselines, EOS/data and literature are not
changed.

### Solve matrix and complete accounting

The tracked solve matrix
`docs/validation/phase6a1_checkpoint_reconstruction_solve_matrix.tsv` has
SHA-256 `32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616`
(reauthenticated; 240 rows; A=237, B=3, C=0, exact endpoint 1, strict interior
239, deep 81). The decomposition is unchanged:

| Local integration | Count | Recovery status |
| --- | ---: | --- |
| Oracle-1 rk8pd | 239 | REUSED, already executed |
| Oracle-2 rk8pd | 239 | REUSED, already executed |
| Replay-1 RKF45 (`rtol=1e-11`, `atol=(1e-16,1e-22,1e-22)`) | 239 | NEW — Stage A |
| Replay-2 RKF45 (`rtol=1e-12`, `atol=(1e-17,1e-23,1e-23)`) | 239 | NEW — Stage B |
| Replay-1 fresh-process repeat | 239 | NEW — Stage C |
| Replay-2 fresh-process repeat | 239 | NEW — Stage D |
| **TOTAL AUTHORIZED LOCAL INTEGRATIONS** | **1434** | 478 reused + 956 new |

Exactly 956 new integrations are authorized. Their solve IDs are the frozen
list `replay1-obs-NNN`, `replay2-obs-NNN`, `replay1-repeat-obs-NNN`,
`replay2-repeat-obs-NNN` over the 239 strict-interior observations, SHA-256
of the newline-joined list
`300022427e5f9efae258e830b1381ea2cd39ddebae72c6f503e995afacfd3298`. Retries
count against the 956; the harness refuses a solve ID that is unauthorized or
already present in the execution ledger and refuses any 957th integration.
The initial-`h` policy (`h0=t_obs-t_left`), local semantics, 100000 apply-call
cap, thermal-domain guard and exact-endpoint rule are those frozen in the
accepted preflight. The exact endpoint observation consumes no integration.

### No-oracle-rerun rule and reused oracle evidence

No rk8pd integration is authorized. The reused evidence is the failed
attempt's untracked, read-only oracle root
`/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase6a1-checkpoint-reconstruction-validation/build/phase6a1-checkpoint-reconstruction-validation/run/oracle`,
authenticated from stored bytes only before this declaration
(`build/phase6a1-checkpoint-reconstruction-recovery/assembly/pre-run-authentication.json`):

| Reused item | SHA-256 |
| --- | --- |
| oracle root tree | `21f1ff9eaf23b078f86aa2b10ec45be55f84a98b9aedc4b94f4a7c4bf999f866` |
| `oracle1` batch tree | `89bc3cb0bde7e58f60e1e44fdee9332e920ae163870213fe79d7bb81ba65ecf9` |
| `oracle2` batch tree | `f0bf0fb4e1e0acef58da1d770c0a71d87c4021264017f75cbfa8aa9da1a0af2d` |
| `oracle1/solve_ledger.tsv` (239 rows) | `0fec6717e153ce19ecafda0fe0ed5c01d1be3934342ac9776f05ad2da769fbcc` |
| `oracle2/solve_ledger.tsv` (239 rows) | `b1bc3feb6bce96b7e4f1feecf2660b1f88781941d66d25ae832e4103a3248f74` |
| historical `oracle-result.json` | `f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0` |
| historical harness executable | `82eea47fba62f1440a1b87abec95dc4ac6287756ad44b990564794b9cf703044` |

Every stored Oracle-1/Oracle-2 result was authenticated: 239 unique solve IDs
per tier equal to the strict-interior matrix rows; ledger left endpoints and
targets bit-identical to the matrix; every result and meta file hash equal to
its ledger entry; 241 result plus 241 meta files per tier with no missing or
duplicate observation; identical identity/currentness fields in all rows; the
exact endpoint and initial observation consumed no integration. The oracle
self-qualification was recomputed from stored bytes only through the unchanged
verifier `oracle` subcommand: the recomputed result JSON reproduced SHA-256
`f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0` exactly,
with `max(d_O/D_O1)=0.5794839113173227` and
`max(U_O/(0.20F))=0.5793505315921852` at observation 117, `x_state`, category
B, not deep, and worst diagnostic utilization `0.17336139714051704` at
observation 117, `Pnet`. **PASS, reproduced.**

### Immutable inputs

| Artifact | SHA-256 |
| --- | --- |
| 241-row observation schedule | `43ec23ada72bfa59c4e89672ae2987e9927164db765f05afb7b45c924cc0c67e` |
| 232 accepted endpoint states | `7f968f6c2fcbf43f442285b48d3b825fa9cf241604f9dbbd88c52c25b96ff459` |
| 241 observation brackets | `0428c176a2d9233add816621a53cc795063480462458d195f3af88d965998b3d` |
| 232 internal accepted steps | `fe9afd8c1ddfc7abefca3f1e57a76c62345be1026d7042f41692650f531742c8` |
| passive trajectory | `8c9531c87b53d8bd189e50802f3d1dd526db6f9b0d256b120c4d635e1e31a98c` |
| passive step summary | `912b0e6300c745f02d733da91fb1072e927aadcd4616b944cc8d02c8a26b2ecb` |
| authenticated production library | `b9b767dbc0114563e1d556e296b6d7fc9d680a9d90e8deae9b44357010dd6499` |
| qualification certificate | `7fc892b50bfbcdd2c5d963e0ff866777f3628f40fc282e1ce5a0b0364200a453` |
| thermal source tree | `1cdb98507958b50fbd1e1eff3eb3bfa4f48c790fef86a8740855884fc43b9d94` |
| frozen certificate | `9217e278eb0a4b06b193a2a02f4d59f285d74e4033427b87ccd51181c3f7beab` |
| frozen coefficients | `3efa060d99a7a92266679aadb31885abbed9c940404247f1d3009cafa94832d6` |
| entry manifest | `5cfbf4b1d623b58fa2a94101631dce20c602950eec8eb514dc1fcd3c7a04620d` |
| historical pretrajectory authority | `1d4e77780fde474f24782879d6ac43fb5bd6d01f3919c359437c0c5f01d327a6` |

These and the oracle evidence are hashed before and after execution and must
be identical. The production library, platform (macOS 26.6.2 / 25G83, Darwin
25.6.0, arm64, Apple clang 21.0.0 `clang-2100.3.34.2`, GSL 2.7.1) and the
compile/link recipe (`-g -std=c++17 -arch arm64 -pthread -Xclang -fopenmp`,
the BA12R Debug generated include, `libCompactStar.a` `b9b767db…`, Zaki,
Confind, GSL, Python 3.12, OpenMP) are those of the failed attempt; production
sources are byte-identical between the library build SHA and this branch.

### Fresh output and scratch roots

All new outputs live under
`build/phase6a1-checkpoint-reconstruction-recovery/` with distinct subroots:
`candidate/{linear,linear-repeat,hermite,hermite-repeat,replay1,replay2,replay1-repeat,replay2-repeat}`,
`candidate/dryrun`, `candidate-work/<method>`, `assembly/` (authentication,
execution ledger, stage checks, determinism, final adjudication). No
failed-attempt candidate output or work root is reused or modified; the oracle
root is read-only.

### Inherited budgets, tolerances and policies (unchanged)

State: `M_i=max(|y_L|,|y_R|,|y_O1|,|y_O2|)`, `F_i=max(atol_ULTRA+1e-11 M_i, 64 ulp(M_i))`,
`B_i=0.25F_i`, `U_O=2max(d_O,F_O)`, pass iff `|y_C-y_O2|+U_O<=B_i`. Ledger:
`G_P`, `F_P=max(1e-11G_P,64 ulp)`, `B_P=0.25F_P` per separated observable.
R20: composite trapezoid on the scientific checkpoint grid; each candidate-
versus-Oracle-2 contribution `<=5e-6 N_R20`; governed `|R20|/N_R20<=2e-4`.
Replay witness `|y_R2-y_R1|<=F_i`. Determinism: byte identity of fresh-process
repeats. Strata A=237, B=3, C=0, D=1, strict 239, E=81. Hybrid eligible only if
Hermite fails exclusively at B/C and Replay-2 passes those points. Selection:
factor-four comfortably passing singles in the order linear, Hermite, Replay-2;
otherwise smallest worst utilization; a passing single outranks hybrid; no
method preselected. None of these values may change after results.

### Execution plan and stop conditions

Order: (1) no-solve dry-run context construction for a no-knot interior
observation, each one-knot observation, the deepest-interior observation and
the exact endpoint; (2) algebraic linear and linear-repeat; (3) Hermite and
Hermite-repeat with isolated endpoint RHS, then endpoint-RHS determinism
check; (4) Stage A Replay-1; Stage B Replay-2; Stage C Replay-1 repeat; Stage D
Replay-2 repeat, each in one fresh child process (process concurrency 1) with a
stage check after each; (5) replay determinism; (6) unchanged `final`
adjudication; (7) before/after immutability hashes. Execution stops and returns
to the owner, without repair, hidden retry or reinterpretation, if: any input,
library, matrix, platform or oracle identity differs; the profile root fails
its byte authentication; the dry run fails; endpoint RHS is nondeterministic; a
stage fails, undercounts, overcounts or reports an error; any solve ID is
unauthorized or repeated; the 956 cap would be exceeded; a repeat differs; any
identity, currentness or validity field differs; R20 needs another quadrature;
or a tolerance, method, subset or budget would need adjustment. Main trajectory
integrations: 0. New integrations executed at this declaration: 0.
