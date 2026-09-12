# Phase-5D-1G0 governed-artifact preparation plan

## Status and authority

This record is the pre-result plan for repairing the Phase-5D-1 reproducibility
harness and preparing a promotion candidate. It is not a scientific review, a
governed-baseline installation, or a canonical integration record. No result
from the repaired producer may be used until this plan is committed.

The controlling identities are:

| Item | Identity |
|---|---|
| Canonical `master` | `d019ae390be4f5e3daba05039903485cb497e397` |
| Controlled-evolution implementation | `d3670f6d4e021def0483909b6d2fdeed1c6973a4` |
| Historical candidate commit | `3486b972f71f57e8351fa8320c1ffb250fcd5c42` |
| Human-ratification entry | `77fda93107b474e2e5420dde1bd84ec55762d727` |
| Frozen-context plan | `f7116c1408c06f976527f86d4397ad6d4540dedf` |
| Historical candidate | `docs/validation/phase5d1_controlled_evolution_candidate.json` |
| Historical candidate SHA-256 | `47da15fa9e32be095a78d1b78ed17c3ffa14e5ef9be5528db3780a76afd0c079` |

The owner-supplied independent Phase-5D-1R review reported 0 BLOCKING, 0
MATERIAL, 7 NONBLOCKING, and 12 NOTES. It regenerated the scientific fields and
trajectories successfully. The controlled mathematical/architecture frozen-v1
science is therefore human-ratified and is immutable in this task.

## Reproducibility findings and dependency inventory

The existing qualification/oracle/trajectory/candidate chain has 29 logical
hidden-state nodes. Seven are sufficient to break the tracked harness in a
fresh checkout: the hard-coded main `build/` CTest tree, its `CMakeCache.txt`, a
pre-existing `entry-hashes.json`, two saved `-entry.rc` files, and the saved
oracle raw return code and run provenance. Four additional hidden nodes serve
qualification, ten are ignored orchestration scripts, and eight are ignored
candidate/final-gate sidecars.

There are 15 literal `ROOT / "build"` or `root / "build"` sites (nine tracked
and six ignored), five additional current-directory-relative hard-coded
`build` sites, no `../build` sites, two trusted saved entry return-code files at
one tracked code site, one hard-coded candidate-worktree path in producer code,
65 absolute worktree-path occurrences in the historical candidate, and one
sibling scratch-build dependency. The historical candidate contains 108
absolute-path occurrences across 38 concrete selectors and 40 repo-relative
`build/...` occurrences.

Every dependency in the repaired pipeline shall be classified as exactly one
of:

1. tracked repository source;
2. authenticated governed external scientific/provenance source;
3. fresh-generated qualification artifact;
4. run-local execution scratch; or
5. forbidden hidden prior-build state.

The fifth class is prohibited. In particular, the repaired producer shall not
inspect arbitrary sibling build directories, pre-existing executables,
pre-existing qualification/thermal/trajectory files, saved return-code text,
or historical candidate output as numerical authority.

## Fresh-context pipeline

The tracked entry point will be `tests/rotochemical/fresh_context.py`. It will
take an explicit `--source-root` and `--scratch-root`; all generated state will
be created below the scratch root. Repository discovery used by subordinate
tools will be based on tracked file location or explicit input, never on a
developer worktree name. A nonempty output location will be refused unless it
is the driver's own newly-created run directory.

For each complete run the driver shall perform this directed sequence:

1. Authenticate a clean source snapshot at the requested source identity and
   build Debug targets into `<scratch-root>/build` with the repository's
   existing compiler-portability rules. Record the actual compiler,
   configuration, and freshly built executable hashes in an execution sidecar.
2. Generate a run-start manifest directly from the tracked frozen-context
   definition. Hash the current source paths for the ratified 33 protected
   paths, the ten governed baselines, and the three specifically protected
   artifacts. No pre-existing manifest is an input. Rehash against this
   run-start manifest at every release gate.
3. Generate the radial/EOS qualification under
   `<scratch-root>/qualification`: run `chemical_trackr_budget.py` with the
   fresh `chemical_trackr_fixture`; run `chemical_structural_envelope`; run
   `phase5b_freegas_validation` for PB6, PB7, PB9-11, PB12, and PB13 in fresh
   directories; recompute the M1/M2/PB structural formulas, `V_K`, and `V_I`
   from those fresh results; read only the predeclared goals from
   `phase5c2_preproduction_evidence.json`; write the transport and certificate
   with the tracked Phase-5C certificate machinery; and invoke the freshly
   built, unchanged `chemical_production_fixture` to create the EOS-8192,
   radial-80000 profile, free-gas table, model, `Z`, `I`, `W`, and Ltilde
   context. Previously generated numerical files may be validation targets but
   never seeds.
4. Construct `eos.t`, `eos.nb`, `eos.yq`, and `eos.thermo` from the fresh
   profile under `<scratch-root>/thermal`. Validate their freshly computed
   hashes against the ratified deterministic byte targets, then retain them as
   run-local inputs. No historical thermal file is read.
5. Execute the Phase-5B and Phase-5C governed regressions in the current build.
   Capture the actual process return codes and require both to be zero before
   trajectory release. Any `.rc` file is output evidence only.
6. Regenerate spin-only, reaction-only, Lyapunov, same-Ltilde coupled RE9,
   architecture-guard, component-tolerance, and mutation-control oracles from
   the current executable and context. A trajectory gate may consume only the
   oracle result produced in the same invocation and bound to the same
   executable hash.
7. Rehash protected and thermal inputs, then run the unchanged P0 = 1 ms
   baseline, refined, initial-T1e7, initial-T1e9, initial-xi1, and initial-xi20
   trajectories. Generate each TSV, Jacobian diagnostic, and complete step
   table afresh.
8. Run the producer-authoritative validator over these outputs. It shall emit
   checkpoints, final/range/convergence results, quasi-steady results,
   checkpoint-bracket/log-interpolation crossing diagnostics, Jacobian and
   stiffness diagnostics, and the step statistics described below.
9. Emit the scientific promotion-candidate JSON and a separate execution
   sidecar. Only after generation may an independent comparator read the
   historical candidate or an existing promotion candidate.

This is the definition of **FRESH-CONTEXT REPRODUCIBLE**: every generated
computational dependency arises from the current invocation plus authenticated
tracked or governed external source authority. A suite may be labelled
`self-contained` only when that definition is true; the inaccurate current
labels will otherwise be removed.

The core controlled producer is data-free. The complete authenticated suite
may additionally consume the existing externally authenticated CompOSE/data
and literature authorities, but only after their governed hashes have been
verified. Missing authentication is a stop, not a skip.

## Producer-authoritative schema

The promotion candidate will be
`docs/validation/phase5d1_controlled_evolution_promotion_candidate.json`.
Its closed schema rejects unknown members and has these exact top-level
sections:

| Section | Equality treatment | Contents |
|---|---|---|
| `schema` | exact | schema name/version and producer contract version |
| `classification` | governed transition | candidate/baseline state and controlled mathematical scope |
| `scientific_contract` | exact | normalization, enabled reactions, frozen definitions, ledger, equations, and unchanged caveat identifiers |
| `fixture` | exact | radial/EOS resolution, central density, initial state, and fixed thresholds |
| `source_provenance` | exact | canonical, frozen-plan, science-implementation identity, and logical protected-source hashes |
| `coefficient_authority` | exact | `Z`, `W`, `I`, Ltilde and support metadata derived by the fresh qualification |
| `thermal_authority` | exact | logical thermal-source hashes and the implemented envelope provenance label |
| `spin_driver` | exact | B, P0, PPdot, spin-law identity, run interval, and mathematical-only classification facts |
| `solver` | exact | baseline/refined tolerances, output grid, limits, and directly emitted step statistics |
| `oracles` | exact | producer-generated scientific oracle results and mutation-family count |
| `trajectory_checkpoints` | exact | the ratified checkpoint payload and ordering |
| `trajectory_summary` | exact | final states, ranges, ledgers, and absolute-power convergence |
| `convergence` | exact | ODE and initial-condition convergence |
| `quasi_steady` | exact | fixed criterion and producer-generated result |
| `thermal_sign_crossings` | exact | fixed roots and producer-generated checkpoint-bracket/log-interpolation diagnostics |
| `numerical_diagnostics` | exact | radial comparison, Jacobian/stiffness, cancellation, and step-budget diagnostics |
| `protected_upstream` | exact | logical 33-path and governed-baseline manifest identity and final unchanged status |
| `suite_results` | exact for named gate/status; live receipts in sidecar | required suite names and zero-failure/zero-unexplained-skip status |
| `producer_provenance` | exact | tracked producer path/version, source identity, and stable logical artifact identifiers |

`science_implementation_sha` remains
`d3670f6d4e021def0483909b6d2fdeed1c6973a4`; the harness run HEAD is execution
information and must not masquerade as a scientific implementation change.
No free-form prose is permitted in the governed JSON.

The historical audit found 434 logical leaf selectors: 31
`SCIENCE_EQUALITY`, 67 `SOURCE_PROVENANCE_EQUALITY`, 25 `FIXTURE_EQUALITY`, 32
`NUMERICAL_METHOD_EQUALITY`, 58 `EXECUTION_PROVENANCE_INFORMATIONAL`, 20
`GOVERNANCE_CLASSIFICATION`, 132 `DERIVED_DIAGNOSTIC_EQUALITY`, and 69
`NARRATIVE_DOCUMENTATION_ONLY`. The first, second, third, fourth, and seventh
classes are equality-bearing. Governance fields are governed by an explicit
state transition. Execution fields are typed informational evidence. Narrative
fields are excluded from the artifact and remain in validation records.

### Historical-to-promotion mapping

The comparator will implement this closed mapping; it will not expose a
generic recursive ignore-path option:

| Historical selector | Promotion section |
|---|---|
| `canonical`, `frozen_context_plan_sha`, `implementation_sha`, `branch_lineage`, equality-bearing `source_hashes`, `protected_gate`, `protected_entry` | `source_provenance`, `protected_upstream` |
| `fixture`, `normalizations`, `initial`, `normalization_classification`, equality-bearing `state_layout` and fixed diagnostic criteria | `scientific_contract`, `fixture` |
| `spin` | `spin_driver` |
| `solver`, `stiffness.step_statistics` | `solver`, `numerical_diagnostics` |
| `Ltilde`, `semantic_W_I`, `semantic_Z`, `Ltilde_support_km` | `coefficient_authority` |
| `thermal_source_hashes` and envelope formula/provenance identity | `thermal_authority` |
| `checkpoints` | `trajectory_checkpoints` |
| `final`, `ranges`, `ledger_residuals`, `power_convergence_absolute` | `trajectory_summary` |
| `ODE_convergence`, `initial_condition_convergence` | `convergence` |
| `quasi_steady` | `quasi_steady` |
| `incremental_root_crossings`, `full_root_crossings`, and aggregate crossing numerics | `thermal_sign_crossings` |
| equality-bearing radial/Jacobian/stiffness/cancellation numerics | `numerical_diagnostics` |
| equality-bearing oracle measurements and suite science gates | `oracles`, `suite_results` |
| `classification`, `candidate_only`, `governed_baseline`, `protected_baseline_installed`, and other typed state facts | `classification` |

Array order and exact key sets are significant. Numbers receive no newly
invented tolerance: exact scientific values must match except where the
ratified comparison already specifies a numerical comparison. A schema
migration may move or type a value but may not conceal a changed trajectory
number.

### Direct step-statistics decision

The two previously postprocessed key names are retained as direct
`NUMERICAL_METHOD_EQUALITY` output:

- `maximum_accepted_steps_per_checkpoint`;
- `last10_accepted_steps_per_checkpoint`.

For each baseline and refined complete `.steps` table the validator will compute
integer per-checkpoint accepted counts as successive differences of cumulative
accepted counts, using zero before the first checkpoint. It will emit the
maximum and the final ten values directly. No postprocessor may patch them.

Crossing times remain deterministic checkpoint-bracket/log-interpolation
diagnostics derived from the generated checkpoint table. This task adds no
event localization and changes no crossing semantics.

## Fields excluded from governed equality

The following exact categories are excluded from the scientific artifact:

- human interpretation, conclusions, recommendations, review dispositions,
  team-review reports, reproduction prose, future-work prose, and free-form
  caveats;
- commit prose, development notes, figure descriptions, archived-draft notes,
  and explanatory definitions already fixed in the governing documentation;
- absolute source/scratch/executable/log paths, temporary-directory names,
  shell command strings, CMake-cache hashes, compiler path, host identity,
  wall-clock times/durations, live branch refs, and run-specific executable
  locations;
- raw return-code filenames, log filenames, evidence-file sizes, and other
  host-specific execution receipts.

The 31 historical top-level additions not emitted by the original validator
are: `validation_status`, `implementation_sha`, `branch_lineage`,
`execution_source_parity`, `semantic_authorities`, `thermal_authorities`,
`Ltilde_support_km`, `state_layout`, `protected_entry`,
`archived_draft_audit`, `oracle_results`, `architecture_guards`,
`all_six_trajectories_raw_rc`, `trajectory_validator_raw_rc`,
`evidence_files`, `caveats`, `future_work`,
`aggregate_incremental_power_crossings`,
`net_power_cancellation_diagnostic`, `implementation_team_reviews`,
`trajectory_figure`, `development_execution_notes`,
`reproduction_prerequisites`, `aggregate_full_power_crossings`,
`pretrajectory_oracle_results`, `preflight_mutation_inventory_coverage`,
`final_test_only_supplement`, `complete_validation`, `final_gate`,
`INV11_status`, and `final_disposition`. Their equality-bearing machine values
are mapped as declared above; their prose and execution receipts are not
copied into the governed JSON.

## Execution provenance and governance treatment

Each fresh run will emit a separate typed execution sidecar. It records the
absolute scratch and source roots, executable locations and hashes, compiler,
configuration, actual commands, actual raw return codes, log identities,
wall-clock details, host-specific information, and run HEAD. These values are
reported and validated for type/completeness but are not scientific equality
fields. Paths will be recorded honestly; no fabricated normalized absolute
path will replace them. Stable logical identifiers and content hashes remain in
the scientific artifact.

The promotion candidate has typed classification:

```text
classification = promotion_candidate
candidate_only = true
governed_baseline = false
benchmark_scope = CONTROLLED_MATHEMATICAL_ARCHITECTURE
physical_spin_interpretation = false
```

It also records P0 = 0.001 s, the fresh fixture mass/radius values, and an
explicit `super_kepler = true` caveat identifier. These are machine-readable
facts, not a claim of physical pulsar realism. The P0 = 1 ms trajectory is not
replaced or altered.

A later promotion-candidate-to-governed-baseline comparison may allow only:

```text
classification: promotion_candidate -> governed
candidate_only: true -> false
governed_baseline: false -> true
```

No scientific, source-provenance, fixture, numerical-method, diagnostic, or
producer-provenance difference is allowed. This task does not perform that
transition.

## Sole metadata correction

The historical label `iron Potekhin1997 envelope` is retained unchanged in the
historical candidate. The promotion mapping makes exactly this value-level
correction:

```text
historical thermal_authorities.photon envelope description
  -> thermal_authority.envelope_provenance =
     "FR2005 Eq. (49) / PCY97 fully accreted"
```

This corrects provenance description only. The formula, implementation,
thermal-source bytes, checkpoints, and every trajectory number must be
unchanged. No other post-hoc metadata correction is authorized.

## Determinism and negative controls

After the repair implementation is committed, the complete producer will run
from two independent clean source snapshots with independent source, build,
qualification, thermal, oracle, and trajectory trees. They will share no
mutable generated state. Producer-authoritative JSON bytes and the
equality-bearing trajectory/Jacobian/step outputs must be identical under the
same toolchain. Execution sidecar differences are reported only under their
informational semantics. A difference in any equality-bearing field stops the
task; exclusions will not be widened after observation. A third run may
diagnose but cannot vote.

The predeclared negative controls are:

1. a fake sibling build containing expected-looking files is ignored;
2. a forged saved `-entry.rc = 0` cannot authorize an actually failing command;
3. a missing qualification producer or required output fails rather than
   skips;
4. a protected-source mutation after entry-manifest creation fails;
5. a thermal-source mutation after construction fails;
6. absence of the historical candidate does not impede generation and affects
   only the later comparison;
7. absence of a pre-existing promotion candidate does not impede generation;
8. changing the absolute scratch root leaves scientific artifact bytes
   identical.

Additional fail-closed checks will reject a nonfresh output tree, missing,
duplicate, or extra manifest membership, wrong radial-80000/EOS-8192 identity,
oracle/executable mismatch, and stale candidate or trajectory seeds.

## Acceptance and stop conditions

Preparation succeeds only when both fresh runs are byte-identical for the
producer-authoritative artifact; the historical-to-promotion projection has no
scientific difference; the only substantive metadata correction is the
predeclared envelope label; all negative controls pass; live Phase-5B and
Phase-5C regressions pass with raw return code zero; retained response and
coupled suites, trajectory validation, complete data-free suite, and complete
authenticated suite have zero failures and zero unexplained skips; and all ten
governed baselines, the historical candidate, Phase-5C candidate, 33 protected
paths, EOS/data, literature, and scientific production code remain unchanged.

The producer construction phase shall not read the historical Phase-5D
candidate, the promotion candidate, or any future Phase-5D baseline. Dedicated
comparison occurs only after producer output is closed. Tests will prove that
both candidate files can be absent during generation.

The task stops if qualification or thermal inputs cannot be regenerated from
governed authority; hidden build state or saved return-code text remains an
input; scientific production code requires modification; two fresh runs
differ; historical scientific equality fails; a governed upstream regression
fails; or any protected scientific artifact changes. The existing unrelated
caveats remain unchanged: the Lother/PBF-only interface, RKF45 stability-bound
behavior, heat-capacity cache kinks, Sommerfeld/free-gas limitations, dot(Z),
DU, superfluidity, super-Kepler physical inadmissibility, realistic A18 work,
and BNV. Global INV-11 remains unresolved.

No governed Phase-5D baseline is installed, no canonical merge is performed,
and no realistic A18 or BNV work begins in this task.
