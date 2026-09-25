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

## Execution result

Appended after execution. The immutable pre-run declaration above is the
version committed at `35292e5f69c535a138fa1d7bcf5ccf430f109814` and was not
amended. The harness binding correction was committed separately at
`a8ab1eeeb13035c54fec5a345f0453a5cb23d653` before any candidate integration.
Classification is unchanged: **PHASE-6 CHECKPOINT-RECONSTRUCTION RECOVERY
VALIDATION EVIDENCE; NOT BNV CANDIDATE; NOT GOVERNED BASELINE; NOT PHYSICAL
RESULT**. Historical BA12 FAIL, BA12R FAIL, passive scheduling PASS and oracle
self-qualification PASS are not rewritten.

### Disposition

**C — NO CHECKPOINT RECONSTRUCTION METHOD QUALIFIED — RETURN TO OWNER.**
No single method and no hybrid satisfied the complete predeclared budget. No
tolerance, method, subset, budget or observation set was adjusted; no further
numerical execution is authorized by this record.

### Binding correction and authorities

- Root cause confirmed exactly as declared: the historical harness bound the
  positional `argv[3]` (passed as the passive trajectory TSV) to the
  EOS/profile directory consumed by `Fixture`, `Campaign::Tangent` and
  `Campaign::Qualification`; the child terminated inside EOS import of
  `<trajectory.tsv>/freegas.tsv` before any exception path.
- Correction confined to `tests/bnv/checkpoint_reconstruction_validate.cpp`
  (source SHA `adc2ce2bf2143809506622b6f66de83044490fa5cc5aed7dd7b5da6ff155f8de`):
  named arguments; `ValidateProfileRoot` refusing a regular file and requiring
  the three authenticated byte hashes; per-solve authorization against the
  frozen 956 solve IDs (`300022427e5f9efae258e830b1381ea2cd39ddebae72c6f503e995afacfd3298`)
  with an append-only execution ledger; one method per fresh child process;
  oracle stages refused. `Integrate`, `EvaluateMethod`, `BuildContext`,
  `ReadMatrix`, `WriteMeta` and `Callback` unchanged.
- `tests/bnv/checkpoint_reconstruction_verify.py` unchanged
  (`0a76c1bb48ee75d55082bd4a4f877c2e7c004104647c7fad5a3d9f244ba56a42`); its
  `oracle` and `final` subcommands were the sole scientific adjudicators.
- Production sources (`CompactStar/`, `EOS/`, `data/`, `literature/`,
  `tests/baselines/`, `CMakeLists.txt`) untouched: forbidden-path diff against
  the preflight base `93e93c7f91a3cd8fced2f7a0961eda9c469c43fe` is empty.
- Harness binary `aae834374615815246087befdd5d1dbb00708cede7471b14aeb82724540887a4`,
  linked against the byte-identical `libCompactStar.a`
  `b9b767dbc0114563e1d556e296b6d7fc9d680a9d90e8deae9b44357010dd6499` used by
  the failed harness and the reused oracle.

### Solve matrix, oracle reuse and accounting

- Matrix SHA `32f3277cdf3984318fe2323da825de1d5e337778bb05106725fb0fcfa0529616`
  authenticated before and after; strata A=237, B=3 (82, 117, 228), C=0, D=1
  (240), strict interior 239, deep interior 81.
- Reused oracle evidence: root tree
  `21f1ff9eaf23b078f86aa2b10ec45be55f84a98b9aedc4b94f4a7c4bf999f866`
  (oracle1 `89bc3cb0bde7e58f60e1e44fdee9332e920ae163870213fe79d7bb81ba65ecf9`,
  oracle2 `f0bf0fb4e1e0acef58da1d770c0a71d87c4021264017f75cbfa8aa9da1a0af2d`;
  ledgers 239/239). Oracle self-qualification recomputed from stored bytes by
  the unchanged verifier: result SHA
  `f27304bd6f85b6b7be20537978ef37777104ff10c4baf25c9ab007e3675d99a0`
  reproduced byte-for-byte; maxima 0.5794839113173227 (state) and
  0.5793505315921852 (obs 117, x_state, stratum B); Pnet diagnostic
  0.17336139714051704 (obs 117). The oracle was not rerun.
- New integrations executed: **956** = Replay-1 239 + Replay-2 239 +
  Replay-1 repeat 239 + Replay-2 repeat 239, exactly the authorized set
  (956 unique solve IDs, set equality with the frozen list, no unauthorized or
  repeated ID, no retry). Verifier campaign count 1434 = 478 reused + 956
  new; `solve_count_pass` true. Main trajectory integrations: 0. Oracle
  integrations: 0. BA12/BA12R integrations: 0.
- Execution ledger (`51125b384eef5cb64ec9ed5f38af96b63736087a253146595b27bf0c32885efd`)
  shows one pid per replay stage (836, 2559, 4527, 6297), disjoint stage time
  intervals in the declared order, i.e. process concurrency 1. The unchanged
  verifier's `final` output carries a hard-coded label
  `process_concurrency: 2` inherited from the historical two-process design;
  it is a verifier literal, not a measurement, and is superseded by the ledger.

### Dry run and determinism

- Dry run (no ODE solve) at observations 1, 82, 117, 228, 37 and 240: PASS;
  all identity, currentness and validity fields equal to the reused oracle
  context; the exact endpoint (obs 240) reproduces the Arm E final state.
- Endpoint RHS determinism (isolated disposable contexts): 478 evaluations,
  203 unique endpoints, 0 mismatches — PASS.
- Algebraic repeats: linear and Hermite byte-identical 241/241; their solve
  ledgers are empty (no integration).
- Replay repeats: Replay-1 and Replay-2 byte-identical 241/241, 0 unequal
  results, max state difference 0, identical step counts and ledger result
  hashes — PASS.

### Adjudication (unchanged `final`)

| Method | Failures | Max state util | Max ledger util | Worst util | R20 normalized | Verdict |
|---|---|---|---|---|---|---|
| Linear | 239 (236 A + 3 B) | 2.0648e8 | 1.149e6 (Pnet, obs 204) | 2.0648e8 | 0.010367 (pass) | FAIL |
| Hermite | 226 (223 A + 3 B) | 2131.3 (p95 380, median 33.2) | 111.85 (Pnet, obs 228) | 2131.3 | 0.002613 (pass) | FAIL |
| Replay-2 | 2 (B: obs 82, 228) | 6.7056 (obs 82, x_state) | 3.4813 (Pnet, obs 228) | 20.24 (replay witness) | 0.002613 (pass) | FAIL |

- Replay-2 stratum detail: A (no knot) max state 0.096283, max ledger
  0.08013, p95 state 0.0857, median 0.0799; E (deep interior) max 0.096283;
  D (exact endpoint) utilization 0; B: obs 82 state 6.7056 / ledger 1.1546,
  obs 117 state 0.6023 / ledger 0.1803 (pass), obs 228 state 3.9564 /
  ledger 3.4813. Replay witness `|R2−R1|/F` maximum 20.239777365061528 at
  obs 228, x_state (stratum B) — witness FAIL.
- R20: all three methods pass; endpoint-propagation, luminosity-difference
  and thermal-difference utilizations 5.1e-5 / 0.519 / 5.1e-5 (linear),
  5e-9 / 2.1e-6 / 5e-9 (Hermite), 9e-11 / 1.1e-8 / 9e-11 (Replay-2). No
  additional quadrature was required.
- Hybrid: ineligible — Hermite fails at 223 stratum-A observations, not
  exclusively at knots; Replay-2 itself fails at knots 82 and 228.
- Selection: none (`selected_method: null`, "no method satisfied the complete
  predeclared budget"). No factor-four comfortable pass exists.

### Performance (measurement only; no scheduling claim)

956 replay solves: wall sum 1.8254 s, CPU 1.8222 s, median 0.0015 s, p95
0.0041 s, max 0.0121 s, 523.7 solves per wall-second. Per-tier 8191-point
solve-only extrapolation: Replay-1 12.6 s, Replay-2 18.6 s. Stage wall
times (context construction dominated, concurrency 1): linear 289 s,
linear-repeat 288 s, Hermite 287 s, Hermite-repeat 287 s, Replay-1 301 s,
Replay-2 306 s, Replay-1 repeat 304 s, Replay-2 repeat 305 s; dry run 292 s.

### Immutability before and after

Pre-run (`pre-run-authentication.json`) and post-run
(`post-run-authentication.json`, `019550cfe506d4d49aac1d76095d0aebd6e686497d7a0c1f9df2cb917fca9db6`)
records agree on every hash field: passive trajectory, steps, internal steps,
accepted states, brackets, schedule, library, matrix, profile root tree
`233114862a2ab6826151519114e4bcb74a01f6a8e739bac0502ee62884688d72`,
certificate, thermal, frozen certificate, coefficients, entry manifest,
pretrajectory, oracle root tree, oracle result SHA and platform. Only the
recorded repository `head` differs (`8b783dbe…` → `a8ab1eee…`). The reused
oracle root was not written.

### Artifacts (untracked build root, not committed)

`build/phase6a1-checkpoint-reconstruction-recovery/`: harness binary,
`oracle-qualified.flag`, `candidate/{dryrun,linear,linear-repeat,hermite,
hermite-repeat,replay1,replay2,replay1-repeat,replay2-repeat}`,
`assembly/` (authorized solve IDs, execution ledger, authentication records,
oracle requalification `f27304bd…`, dry-run and stage logs, stage checks,
stage timeline, determinism records, `final-result.json`
`95a370f248c55e075d6e29389cadb0e68ca0b819a08767cedf114a3cb3838899`,
`recovery-summary.json`
`66c5dc33b827f222c0901a86b53a932218c05645f46ca5203b42b31d76d5c314`).

### Boundaries preserved

No oracle rerun; no main trajectory; no BA12/BA12R rerun; no tolerance,
observation, method or budget change; `ScaledRKF45`, `Cstar`, Phase-5D and all
production sources unmodified; no BNV candidate; no production ADR; no
cluster use; not merged. Recommended next action: return to owner; no further
numerical execution under this declaration.
