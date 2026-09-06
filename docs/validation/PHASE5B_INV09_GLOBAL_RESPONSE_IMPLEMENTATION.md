# Phase-5B INV-09 global particle-number structural response

**PHASE-5B GLOBAL PARTICLE-NUMBER RESPONSE IMPLEMENTED AND CANDIDATE VALIDATED — INDEPENDENT REVIEW REQUIRED**

ADR-0011 structural contract: **IMPLEMENTED IN CANDIDATE**. PB1–PB14: **CANDIDATE PASS**. Independent review and human ratification: **PENDING**. INV-09 remains **INTENDED BUT UNVERIFIED**; INV-11 remains **UNRESOLVED**. This record and its reference artifact confer no human ratification, governed baseline installation, or canonical integration.

## Authentication, authority, and preserved history

Canonical master, origin/master, live master, candidate HEAD, and merge-base were authenticated at `a43d02227bf53c3242d3212f81dd71963804f3aa`. The canonical checkout was clean. The authorized existing candidate was `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5b-global-response`, branch `physics/phase5b-global-particle-response`, 0 ahead / 0 behind, without an upstream or live branch at entry. No new branch/worktree, cleanup, reset, restore, stash, rebase, or discarded work was used.

Entry tracked diff SHA-256: `47f95960fb0d08fe458b88949cb99cde906470ac4ee4186be04ae761d1e0bbc7`. All five modified tracked files and all eight untracked-file hashes matched the forensic snapshot. The six source/test untracked hashes also matched the historical stop record. The two remaining document hashes were checked against `/private/tmp/compactstar-phase5b-forensic.rVQXC8/entry-untracked-sha256.txt`. Evidence and exact hashes: `docs/validation/phase5b_resume_evidence.json:1`.

One agent, one production writer. Change classes: numerical representation, additive scientific-semantic structural driver, API/provenance, test/build registration, documentation and candidate artifact. Authority: `GOVERNANCE.md:39`, accepted ADR-0011 Q1–Q3 and PN1–PN8 (`docs/adr/ADR-0011-particle-number-structural-response.md:31`), ADR-0001 nodal species semantics, ADR-0003 profile identity/version lifetime semantics, ADR-0008 signed-measure semantics, ADR-0009 declared finite-cut surface, and the explicit Phase-5B-R owner instructions. ADR-0012 is complete for ordinary-NStar A1 and closes the Q4 unit prerequisite (`docs/adr/ADR-0012-relativistic-unit-boundary.md:3`). No pre-baseline correctness exception is invoked; this additive structural result has no pre-existing governed response baseline.

The original PB6 execution stopped before emitting a radial response row. It did not run PB7 or later gates. Its PB1–PB5 claims were partial and its factor-16 tail enclosure was unaccepted. The immutable historical `phase5b_stop_evidence.json` remains byte-identical, SHA-256 `cde603c6576babb50cd7622c10a9f7868cc78eaf766f07709c8d0bb7dea97d62`. Its commit/push/status fields describe that historical stop, not the later candidate disposition. The old implementation document was additionally preserved at `/private/tmp/compactstar-phase5br.5Txqjk/original-implementation.md` before this update.

## PB6 diagnosis and representation correction

The first edit added refusal diagnostics only. It preserved the comparison, segment values, measure, tolerances, jump behavior and terminal term. The necessary target was rebuilt, and the same 8192-interval / 10375-row table and first radial-resolution-10000 midpoint star were regenerated in fresh output. The table hash reproduced exactly: `7cd44c92e1e7206e0e68e3fed7e3f0ca68e79ab4517d02b96ff78b9be23d3f1a`.

Original refusal: `STOP unrepresented gap/jump in number measure`. Fresh diagnostics identified neutron species index 0, label `10`, segment index **2517**, expected shared original profile node **2516**, an ordinary continuous node with no declared jump, not a terminal or explicit-domain boundary.

| Endpoint | Radius [km] | Density [fm^-3] | Density hexfloat | Density bits |
|---|---:|---:|---|---|
| Previous segment right | 12.642323193955173 | 1.7560842199532513e-07 | `0x1.791dc03ee3ep-23` | `0x3e8791dc03ee3e00` |
| Current segment left / canonical source | 12.642323193955173 | 1.7560842199532527e-07 | `0x1.791dc03ee3e05p-23` | `0x3e8791dc03ee3e05` |

Both radius values are `0x1.948de95eeffddp+3`, bits `0x402948de95eeffdd`: **0 ULP** apart. Density is **5 ULP** apart. The independently evaluated `a + 1*(b-a)` endpoint lost five ULP relative to the stored node; the next segment used `t=0` and recovered the original node. All five root-cause conditions A–E were established before repair. Diagnostic-only source and log hashes are retained in `phase5b_resume_evidence.json:1`.

The correction builds a canonical endpoint vector, copies original profile nodes exactly, interpolates only an explicit clipped domain endpoint once, and constructs adjacent segments from the same stored endpoint. No arbitrary radius/density snapping, epsilon continuity, invented atom, threshold detector, or TOV/Hartle change was made. The exact equality refusal remains (`CompactStar/Analysis/src/ParticleNumberResponse.cpp:180`, `:258`). The optional EOS-knot validation adapter likewise constructs each refinement endpoint once and shares it between segments (`tests/analysis/phase5b_freegas_validation.cpp:79`).

Endpoint micro-falsifiers: **8/8 PASS**. Exact shared continuity passes; independent 1-ULP radius and density perturbations refuse; an undeclared finite jump refuses; a declared jump contributes exactly one signed atom; continuous onset creates no atom; the finite terminal contributes once; duplicate and omitted terminal mutations both fail the independent boundary oracle (`tests/analysis/particle_number_analytic.cpp:225`).

## Public API, ownership, domain, and provenance

One generic `CompactStar::Analysis` module supplies `ParticleNumbers`, `FixedCentralEnergyNumberResponse`, `EquilibriumSequenceNumberDerivative`, and `FixedBaryonNumberResponse`. The existing `Core::Analysis` extension seam remains unchanged. The structural module is independently callable, is not a thermal/BNV driver, and creates no persistent cache hierarchy (`CompactStar/Analysis/ParticleNumberResponse.hpp:97`).

`NumberDomain` explicitly distinguishes `WholeStar` and `ExplicitFixedIsobar`, with a named outer pressure and optional moving inner isobar. No core is inferred from saturation, crust, Urca, species onset, arbitrary pressure, or an FR2005 figure. `ReduceDomain` always uses the whole-star baryon constraint. For a core/shell, the reported enclosed baryon response is boundary flux rather than a zero-baryon residual (`CompactStar/Analysis/src/ParticleNumberResponse.cpp:474`).

Results are eager values with validity/refusal reasons. Numerical access checks the owning profile address/version under the existing lifetime rule, both Hartle sources, EOS identity/revision/domain and exact table contents. Sequence ownership retains all **15 separately solved stars**, each own profile/grid/metric/composition/surface; reductions retain central plus sequence dependencies. Metadata includes species charges/order, units, central-state definition, actual achieved states, steps, branch policy, radial resolution and node counts, and every realized tail correction/bound plus callback identity/revision. The callback contract is pure under an immutable declared policy; no lazy callback is retained. Refinement knots are absent in production and explicitly identified as a table-backed validation partition. Runtime addresses are not persistent object IDs and are not serialized as such (`CompactStar/Analysis/ParticleNumberResponse.hpp:33`, `:52`, `:64`, `:109`).

The existing NStar/RotationSolver provenance plumbing remains byte-identical to the authenticated entry state. No Hartle/TOV equation, coefficient, normalization, integrator, tolerance, step policy, center condition or surface matching changed. The new count-result tail dependency now also retains both owning Hartle sources. Stale profile, first-order, monopole, EOS revision/domain/table bytes and changed sequence contributors all refuse; near-zero B_B refuses without fallback (`tests/analysis/phase5b_freegas_validation.cpp:385`).

## PN1–PN8 and source traceability

| Contract / source | Public symbol / operation | Sign, units, domain |
|---|---|---|
| PN1; FR2005 (20),(21) | `ParticleNumbers::Compute` | `1e54 integral 4pi r^2 n_i/sqrt(1-2m/r) dr`; n_i=Y_i*n_B reconstructed at nodes; count; no lapse; declared D |
| PN2; FR2005 (24), Hartle (109)–(114) | `FixedCentralEnergyNumberResponse::Compute` | density measure `-integral w*xihat dn_i`, positive metric/velocity volume terms, outer-minus-inner boundary; count km^2 per q=Omega_geom^2; no l=2 |
| PN3; FR2005 (26),(27) | `EquilibriumSequenceNumberDerivative::Compute` | complete canonical stars; derivative with respect to geometric epsilon_c; count km^2; fixed EOS/domain policy |
| PN4–PN6; FR2005 (18) | `FixedBaryonNumberResponse::Compute` | depsilon_c/dq=-A_B/B_B; K=A+B*shift; raw baryon/charge constraints; whole-star total baryon number |
| PN7–PN8; FR2005 (30),(31) | `ReduceDomain`, `FixedIsobarIGeometric` | I_geom=K-Y_outer Q_outer+Y_inner Q_inner; count km^2; explicitly declared core/shell |
| Whole-star PN7/PN8 limit | `WholeStarIPhysical` | I_geom=K only at whole-star zero total baryon flux; I_phys=K/c^2 [count s^2] |
| FR2005 (30), structural term in R2006 (7) | `WholeStarEquilibriumNumberRate` | +2 I_phys Omega_phys Omega_dot_phys [count/s] for equilibrium numbers. R2006 (7) writes the excess-number equation, where this equilibrium driving is subtracted. No excess-number evolution is implemented. |

Literature remained read-only. Direct FR2005 formula text and R2006 excess/equilibrium sign were checked against the catalogued sources; hashes freshly authenticated: FR2005 `f184d7d1d7030b61a021eb5c7ac14b1f1b30c7ea69e9d53473d153cfb069ea88`, R2006 `a286f15e083e52becd95b3000cbb5ec3ed97148681cf10a43f1a1cc5c4d23ae8`, Hartle 1967 `ed263946e9bc13842399b5c9e9c2eae31823e7323bc81b456fb5174697cefc35`. The accepted preflight source ledger and derivation remain authority (`docs/validation/PHASE5B0_INV09_GLOBAL_RESPONSE_PREFLIGHT.md:75`, `:191`, `:392`). No numerical free-gas core I_Omega or source-qualified Mmax reproduction is claimed.

Every positive-radius profile node is a measure integration boundary; nodal density products are interpolated linearly within a segment. The regular missing-center interval uses the existing candidate center extension. Exact density increments and compensated accumulation are used. Continuous neutron/muon onsets have no atoms; explicitly declared finite jumps have signed Stieltjes atoms; the finite terminal atom is the boundary term exactly once. Shells subtract the inner moving boundary (`CompactStar/Analysis/src/ParticleNumberResponse.cpp:131`, `:258`).

Production B uses x=ln(rho_c/rho_ref), steps 0.001/0.0005/0.00025, achieved-abscissa three- and five-point estimates, and materializes dN/d epsilon_c. No inverse-mass fit or central-star geometry reuse is used. The homogeneous/sensitivity solution remains private test-only code and was executed fresh at two numerical settings (`CompactStar/Analysis/src/ParticleNumberResponse.cpp:349`; `tests/analysis/particle_number_homogeneous.hpp:1`).

## Primary fixture and raw structural results

The selected fixture is the human-ratified Structure-1 midpoint, rho_c=1.10e15 g/cm^3, table 8192 intervals / 10375 rows, radial setting 80000 / 20259 actual profile nodes. It is source-compatible with the printed FR2005 common static state, not a published I_Omega numerical benchmark. Species order `[n,p,e,mu]` maps to profile labels `[10,11,0,1]`, baryon charges `[1,1,0,0]`, electric charges `[0,1,-1,-1]`. No Mmax tuning is used.

| Species | N [count] | A [count km^2] | B [count km^2] | K=I_geom [count km^2] | I_phys [count s^2] |
|---|---:|---:|---:|---:|---:|
| n | 7.5683115419394785e+56 | 3.9647718619373394e+59 | 1.6255270466118125e+59 | 1.1793897334337809e+58 | 1.3122480530141585e+47 |
| p | 4.8594976779411839e+54 | 1.8404826326625264e+57 | 5.7613767070193853e+57 | -1.179389733433783e+58 | -1.3122480530141607e+47 |
| e | 4.8119204549793175e+54 | 1.8292534582866583e+57 | 5.192854858061696e+57 | -1.0459711462551772e+58 | -1.1637998545112904e+47 |
| mu | 4.7577222961866916e+52 | 1.1229174375868009e+55 | 5.6852184895566917e+56 | -1.3341858717812756e+57 | -1.484481985023382e+46 |

| Aggregate | Value |
|---|---:|
| A_B | 3.9831766882639647e+59 |
| B_B | 1.6831408136820063e+59 |
| B_B_error | 1.2652133015974696e+52 |
| baryon_budget | 6.0150855100049548e+52 |
| baryon_residual | -2.0906948623622459e+43 |
| central_epsilon_km_minus2 | 0.00081671852064503763 |
| charge_budget | 5.2889072695244279e+51 |
| charge_residual | -4.7821160485099105e+45 |
| conditioning | 1 |
| depsilon_dq | -2.3665142309456835 |
| depsilon_dq_error | 1.7868631534267983e-07 |

These are raw values: no dependent species was overwritten to manufacture baryon or charge closure. Independent neutral reconstruction is a separate validation result. The candidate JSON contains full ingredient errors, stencil estimates, every achieved state and per-star tail inputs (`docs/validation/phase5b_structural_response_candidate.json:1`).

## PB1–PB14 validation and uncertainty

| Gate | Kind and actual result |
|---|---|
| PB1 | PASS: analytic proper-volume oracle relative 1.2112533198660458e-13; independent 32-point nodal count <=1e-8; independent direct phase-space / enthalpy TOV counts lie within the predeclared sum of radial, EOS, oracle-setting and tail budgets; raw/reconstructed baryon and charge identities. |
| PB2 | PASS: analytic lapse cancellation, 1/3 angular average, l=2 orthogonality and independent integration by parts (5.5511151231257827e-16 relative). |
| PB3 | PASS: signed jump/ramp measure (4.1422421048764591e-13 same-representation relative); shrinking-ramp weight envelopes decrease from 3.7037e-4 to 3.7037e-10; exact endpoint falsifiers preserved. |
| PB4 | PASS: all four independent analytic contributions; radial ladder 256–16384; finest relative 3.0267983675003052e-9; each omission fires. |
| PB5 | PASS: independent nonlinear toy current at q=[1e-3,5e-4,2.5e-4,1.25e-4]; quotient errors [0.0030768930592284249,0.0015382563567367669,0.00076908051549473555,0.00038452799411281546] approach linearly; exact physical/geometric conversion and q/seed/c^2 mutations. |
| PB6 | PASS after confirmed representation repair: radial settings [10000,20000,40000,80000], all nonzero species, independent immutable EOS knots and both onset cells; unchanged 1e-3 per-species / 1e-8 sum|A| floor. Finest knot relative differences [9.5908665557846801e-9,3.9219973058174648e-11,1.3312970348933061e-10,2.8115286899092991e-8]. |
| PB7 | PASS: fresh canonical A1 sequence versus independently derived homogeneous oracle; relative differences n/p/e/mu=[9.9431884303946561e-7,8.4305491032665714e-7,8.0256827628133465e-7,1.2128547399736789e-6], each <=2e-4. Oracle setting differences <=4.21e-9. This qualitatively reproduces a small convention-coherent discrepancy; historical Unit-1 numbers were never a target. |
| PB8 | PASS: 15 complete independent stars, three step scales, three-/five-point estimates and actual centers recorded; implemented adjacent-step and conditioned error gates <=2e-4. Tables are in the candidate JSON and resume evidence. |
| PB9 | PASS: independent long-double PN4 expansion, B_B conditioning, raw baryon residual within propagated budget; sign and baryon-map mutations fire. |
| PB10 | PASS: independent finite-current extrapolation plus independently quadratured neighboring counts, neutral reconstruction and raw charge residual within propagated budget. Charge/order mutation fires. |
| PB11 | PASS: new complete stars at each linearly shifted central energy; independent finite-q current, full Jacobian and angular Lorentz factor. Baryon q-halving ratios [0.25086803070455282,0.24941095600815485,0.25057825160373287]; all species quotient-error ratios [0.4982,0.5050]. Suppressing the shift fires. |
| PB12 | PASS: +0.410% dEdP-only mutation on 133 noncentral nodes gives exactly zero in owning Hartle fields, N, A and K. EOS settings [2048,4096,8192], maximum K difference 2.1665905514112183e-7 <=1e-3. Four complete central configurations at 0.99/1.01 of neutron and muon onset; lower-density stellar radius search explicitly uses existing SetMaxRadius(50000 km); no solver equation/tolerance changes. |
| PB13 | PASS: independently integrated direct-vacuum current/tail, comparison inequalities, +/-0.0005 and +/-0.00025 neighboring-center tail derivatives, and all four PB11 shifted-center tails. Species absent before the surface have exact zero tails. The provisional factor-16 envelope is replaced only after the comparison test passed. |
| PB14 | PASS: executable analytic whole/core/shell PN7/PN8 mapping, flux sign, count/c^2/physical-driver units, domain refusals and source traceability. No source-authenticated free-gas core domain is invented. |

PB6 full per-resolution/per-partition N and A, density measure, metric, velocity, outer boundary, quadrature error, onset-cell contributions and tail estimates are preserved as machine-readable evidence lines in `phase5b_resume_evidence.json:1`. Onset diagnostics explicitly describe the cells crossing support termination, not an arbitrary wider threshold region. These cells contain no atoms.

PB13 midpoint: R_cut=12.766174760730493 km; direct-vacuum continuation R=12.7681549034222 km. Each charged-species tail has N_tail=1.0959016269095933e43 count, bounded by 2.7545531560385803e43. Interior A tail=6.6962407445770081e49 count km^2; subtract the already included cutoff atom 6.6935294316024069e49, giving net 2.7113129746011584e46. The triangle/comparison replacement bound is 1.3396213421634935e50, below 1e-6 of the response scale. dN_tail/dx at two steps is 4.1787221246134659e42 and 4.1787223016945989e42; neutron/muon tails are zero. Direct setting differences are explicitly recorded. No finite-cut radius is silently relabelled P=0.

The test-side tail enclosure follows the preflight section 9 inequalities: monotone pe density/energy, a Schwarzschild lower-gravity radius enclosure, bounded mass and proper volume, positive first-order comparison, upper monopole/displacement bounds and an explicit positive lower bound for phat. All comparison-domain conditions are checked before accepting the bound. This is a fixture-specific validation adapter, not generic EOS physics (`tests/analysis/phase5b_freegas_validation.cpp:303`; `tests/analysis/particle_number_homogeneous.hpp:1`).

Reported `Errors()` budgets include in-segment quadrature/roundoff, declared tail bounds and (for B/K) propagated stencil/ingredient estimates. They are not claimed to be a universal rigorous bound on EOS and profile interpolation. Those errors are separately measured by PB1/PB6/PB7/PB12 and retained in the validation artifact. Exact raw conservation is not enforced by editing a species. No accepted tolerance was widened.

## Mutation and provenance falsifier matrix

| No. | Controlled mutation | Executed firing evidence |
|---:|---|---|
| 1 | n_i -> Y_i | FIRED: PB1-PB5-analytic / `MUTATION FIRED n_i -> Y_i` |
| 2 | omit density measure | FIRED: PB1-PB5-analytic / `MUTATION FIRED omit density-measure term` |
| 3 | omit metric | FIRED: PB1-PB5-analytic / `MUTATION FIRED omit metric term` |
| 4 | omit velocity | FIRED: PB1-PB5-analytic / `MUTATION FIRED omit velocity term` |
| 5 | omit outer boundary | FIRED: PB1-PB5-analytic / `MUTATION FIRED omit outer-boundary term` |
| 6 | double terminal atom | FIRED: endpoint-micro / `ENDPOINT 7 PASS` |
| 7 | omit terminal atom | FIRED: endpoint-micro / `ENDPOINT 8 PASS` |
| 8 | invent continuous-onset atom | FIRED: PB1-PB5-analytic / `MUTATION FIRED invent continuous-onset atom` |
| 9 | wrong angular factor | FIRED: PB1-PB5-analytic / `MUTATION FIRED wrong angular 1/3` |
| 10 | wrong q normalization | FIRED: PB1-PB5-analytic / `MUTATION FIRED wrong q normalization` |
| 11 | wrong c squared conversion | FIRED: PB1-PB5-analytic / `MUTATION FIRED wrong c squared conversion` |
| 12 | reuse central metric | FIRED: contracts / `MUTATION FIRED reuse central-star metric` |
| 13 | flip fixed-baryon sign | FIRED: PB9-11 / `MUTATION FIRED flip fixed-baryon` |
| 14 | wrong baryon map | FIRED: PB9-11 / `MUTATION FIRED wrong baryon species map` |
| 15 | wrong charge map/order | FIRED: PB9-11 / `MUTATION FIRED wrong charge map/order` |
| 16 | stale profile | FIRED: contracts / `MUTATION FIRED stale profile provenance` |
| 17 | stale first order | FIRED: contracts / `MUTATION FIRED stale/mismatched first-order` |
| 18 | stale monopole | FIRED: contracts / `MUTATION FIRED stale/mismatched monopole` |
| 19 | near-zero B_B | FIRED: contracts / `MUTATION FIRED near-zero B_B` |
| 20 | suppress central shift | FIRED: PB9-11 / `MUTATION FIRED suppress PB11 central shift` |
| 21 | nodal dn_i/dp source | FIRED: PB1-PB5-analytic / `MUTATION FIRED nodal dn_i/dp` |

All **21/21** required mutations fired. The **8/8** additional endpoint micro-falsifiers passed. Extra controls cover lapse, count units/metric/domain, l=2 contamination, signed measure, EOS revision/domain/exact bytes, sequence contributors and core/whole confusion. Mutants exist only as controlled test-side inputs/calculations or disposable star state; no mutation remains in production. The test EOS byte mutation is restored before subsequent checks. No scientific failure was patched through or tolerance loosened.

## Timing, repeatability, suites and protected scope

Representative serial timings: ParticleNumbers 0.278014917 s; A 0.233659167 s; full three-scale B sequence 19.135139583 s; K assembly with repeated provenance validation 5.487549042 s; PB11 four-star finite-spin ladder 12.934065750 s. Timing is measured operational evidence, not a guaranteed performance contract. No persistent optimization cache was added.

Canonical producer: `tests/analysis/produce_particle_number_reference.py:1`, invoking the `emit` mode of `phase5b_freegas_validation`. Both fresh generations solved all contributors and wrote isolated outputs; neither read historical scratch coefficients. Source hashes, exact table authority, schema, units, profile versions/node counts, actual centers, stencil and tail policy are retained. Runtime profile addresses are deliberately not serialized or converted into persistent IDs.

Candidate artifact: `docs/validation/phase5b_structural_response_candidate.json`.

Generation 1 SHA-256: `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`.
Generation 2 SHA-256: `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`.
Exact byte-repeatability: **PASS**. This is a candidate reference, **not an installed governed baseline**.

| Validation selection | Inventory / outcome |
|---|---|
| Focused | 8/8 PASS; affected analytic executable additionally rechecked 1/1 after final linking |
| Complete data-free, serial | 38 / PASS |
| Complete external-data, serial after data-free | 61 / PASS |

The initial focused launch overlapped final linking of the additive physical-rate method; its analytic check was rerun successfully after linking. Complete suites use completed builds. All test invocations receive fresh isolated directories through `run_phase5b_validation.py`; full suites are serial. Independent build/reference preparation used separate output directories without shared generated files. A duplicate test-harness main declaration caused one build error and was corrected before any affected test ran; no scientific gate failed.

Historical baseline diff: empty. TOVSolver, Geometry, RelativityUnits and Track-R thermodynamic physics are byte-identical to HEAD. RotationSolver/NStar changes remain the authenticated provenance-only entry edits. Literature is unchanged and no paper entered Git. Production Analysis has no Track-R masses, thresholds, Sigma ceiling, fixture density, inferred core cutoff, BNV or chemical dependency. Legacy RotochemicalCache is not activated (`phase5b_resume_evidence.json:1`).

NO BTILDE, CHEMICAL Z/W, EVOLVED ETA, WEAK-RATE EVOLUTION, OR BNV IS IMPLEMENTED.

## Remaining disposition and exact next action

All required candidate gates, mutations, repeatability and regression suites passed. The commit containing this record is PHASE5B_SHA; its exact hash and post-push local/upstream/live equality are reported in the final handoff. No merge or human ratification is performed.

Independent scientific/numerical review, human ratification and any governed baseline/integration remain separate. Source-qualified free-gas core I_Omega and Mmax are not claimed. Generic EOS adapters remain responsible for their declared domain/tail authority; the pe tail comparison here is specific to the validation adapter. INV-09 remains INTENDED BUT UNVERIFIED; INV-11 remains UNRESOLVED.

Exactly one recommended next action after successful candidate completion: Run an independent Claude Opus 5 XHIGH scientific/numerical review of PHASE5B_SHA, including the confirmed PB6 representation repair, PN1-PN8, all PB1-PB14 evidence, finite-spin closure, sequence-derivative independence, measure-complete dn_i implementation, provenance/refusal behavior, mutation ladder, domain/source mapping, and the deterministic candidate structural artifact before any human ratification or INV-09 closure. Do not begin that action automatically.
