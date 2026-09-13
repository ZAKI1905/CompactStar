# Phase-5D-1 controlled rotochemical evolution — governed integration record

## Disposition

> **PHASE-5D CONTROLLED NON-SUPERFLUID ROTOCHEMICAL EVOLUTION —
> GOVERNED / COMPLETE VALIDATION PASS / READY FOR CANONICAL INTEGRATION
> FOR CONTROLLED FROZEN-V1 SCOPE.**

This record closes the branch-side validation needed to fast-forward the
controlled Phase-5D-1 implementation into canonical history. It does not widen
the scientific contract, change the ratified trajectory, regenerate or replace
the governed baseline, begin realistic A18/FR2005 work, or begin BNV. Canonical
integration is not claimed until the documented branch has been pushed and
`master` has been advanced by fast-forward only.

## Authenticated authority

| Authority | Authenticated value |
|---|---|
| Canonical entry | `d019ae390be4f5e3daba05039903485cb497e397` |
| Candidate branch | `physics/phase5d-controlled-rotochemical-evolution` |
| Evolution implementation | `d3670f6d4e021def0483909b6d2fdeed1c6973a4` |
| Candidate evidence | `3486b972f71f57e8351fa8320c1ffb250fcd5c42` |
| Human ratification | `77fda93107b474e2e5420dde1bd84ec55762d727` |
| Governed-artifact plan | `b53f26afeaa29e2e43744e27edda2cb5a94fa7c6` |
| Fresh-context repair | `da1789ea152e2f0f3101f1ffecd54144e9d093d8` |
| Artifact-preparation record | `4a0d4fbc087dd683e97897edfd3eaa944a7bddd4` |
| Governed baseline commit | `92c7b4ef384284be8870f49c98c385fb4e230084` |
| Python isolation commit | `fa61467b4bbfdebd048b153b603c9f93bd94e9c8` |

The historical ratified artifact remains
`docs/validation/phase5d1_controlled_evolution_candidate.json`, SHA-256
`47da15fa9e32be095a78d1b78ed17c3ffa14e5ef9be5528db3780a76afd0c079`.
The fresh-context promotion candidate remains
`docs/validation/phase5d1_controlled_evolution_promotion_candidate.json`,
SHA-256
`0fc183202f9ca5651d0c600546cd0ca3893629b4db076c19d31fb34db84c21fa`.
The governed artifact is
`tests/baselines/phase5d1_controlled_evolution.json`, SHA-256
`2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd`.
Its installation increased the governed baseline count from ten to eleven.

The Phase-5B baseline remains SHA-256
`7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`;
the Phase-5C baseline remains
`7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7`;
and the reviewed Phase-5C candidate remains
`a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b`.

## Governed promotion evidence

The governed producer is `tests/rotochemical/fresh_context.py`; the schema
comparator is `tests/rotochemical/compare_artifacts.py`; and the canonical
expensive regression is
`tests/rotochemical/phase5d1_controlled_evolution_regression.py`. The producer
starts from tracked current source and a supplied nonexistent scratch root. It
does not read the historical candidate, promotion candidate, or governed
baseline to construct its result. Execution-path and command provenance remain
in `execution-sidecar.json`, outside scientific byte equality.

Two independent empty-scratch governed runs completed before baseline
installation:

| Run | Scratch root | Governed artifact SHA-256 | Result |
|---|---|---|---|
| 1 | `/private/tmp/phase5d1int.hSkoSH/run-1` | `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd` | PASS, all producer raw return codes 0 |
| 2 | `/private/tmp/phase5d1int.hSkoSH/run-2` | `2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd` | PASS, all producer raw return codes 0 |

The two scientific artifacts were byte-identical; no third-run vote was used.
A post-governed-commit regeneration from
`/private/tmp/phase5d1int.hSkoSH/post-governed` produced the same SHA-256.

The promotion-candidate-to-baseline comparator recorded exactly these three
differences and no others:

```text
classification.candidate_only: true -> false
classification.classification: promotion_candidate -> governed
classification.governed_baseline: false -> true
```

There were no unexpected scientific, source-provenance, fixture,
numerical-method, diagnostic, or metadata differences. The baseline is direct
producer output and was not edited after generation. The producer directly
emits authoritative step statistics, checkpoint-bracket/log-interpolated root
crossings, and quasi-steady diagnostics; no event-localization or hand
postprocessing change was made.

The non-self-comparison controls all passed: generation completed with the
baseline absent and the comparison then refused; mutations of a numerical
trajectory value, Z or W, `Ltilde`, fixture metadata, solver tolerance, source
provenance, envelope label, and physical-spin interpretation all failed; and
only the dedicated promotion-to-governed classification transition was
accepted. The actual baseline was never mutated.

The retained fresh-context negative controls also passed: fake sibling build
ignored; forged saved rc unable to override live failure; missing qualification
refused; protected-source mutation refused; thermal-source mutation refused;
historical candidate absence did not block production; promotion candidate
absence did not block production; and alternate scratch roots preserved
scientific bytes. Before and after complete producer runs, the 33 upstream
protected paths were unchanged. Live Phase-5B and Phase-5C governed regressions
returned rc 0.

## Original complete-suite stop and root cause

The first baseline-promotion validation attempt correctly stopped during the
complete authenticated suite. `chemical_trackr_budget` passed in 91.14 seconds,
but its direct source-tree import of `chemical_production_evidence.py` allowed
CPython 3.12 to create:

```text
tests/analysis/__pycache__/chemical_production_evidence.cpython-312.pyc
```

The next test, `phase5d1_controlled_evolution_regression`, invoked the fresh
producer. Its clean-source gate correctly refused the untracked bytecode in
0.25 seconds. The preserved audit copy is
`/private/tmp/phase5d1int.hSkoSH/interrupted-full-suite-generated-pycache/chemical_production_evidence.cpython-312.pyc`,
SHA-256
`fb482c5b454904dd52c5c48a510a9156e4901eec06510d9737069fee6cdcdf15`.
No scientific artifact, governed baseline, or protected source had drifted.

The failure was not repaired by ignoring, deleting, or cleaning contamination,
and the Phase-5D clean-worktree gate was not weakened. The protected script
`tests/analysis/chemical_trackr_budget.py` remains byte-identical at SHA-256
`397923d762d6ff6687b9b2d0bb4d7bbee5abd97265426fc38d225db6aeb8e56f`.

## Python CTest isolation repair

The CTest registration audit found 18 Python validation tests. They run the
configured `Python3_EXECUTABLE` with build-tree working directories while
loading scripts and sibling imports from source directories. There was no
common registration helper and none had a bytecode-isolation environment
property. The smallest systematic launch-level repair is in
`tests/CMakeLists.txt`: the 18 repository Python validation registrations receive
`PYTHONDONTWRITEBYTECODE=1` through one CTest property list.

This policy prevents the root cause rather than hiding it. Bytecode caching has
no scientific value in validation, arbitrary user Python behavior is unchanged,
no `.gitignore` rule was added, and no after-test deletion or pre-regression
cleanup occurs.

Focused and ordering-dependent evidence from build root
`/private/tmp/phase5d1int-r1.ePbGT5/build` is:

| Gate | Result |
|---|---|
| Focused `chemical_trackr_budget` | PASS, 1/1, rc 0; source clean; no source `__pycache__` or `.pyc` |
| Forward order: chemical then Phase-5D regression | PASS/PASS, rc 0/0; no cleanup; source clean |
| Reverse order: Phase-5D regression then chemical | PASS/PASS, rc 0/0; no cleanup; source clean |
| Registered Python contamination sweep | PASS, 18/18 registrations; cleanliness checked after every test; no additional contaminator |

The forward Phase-5D artifact, reverse Phase-5D artifact, authenticated-suite
Phase-5D artifact, and post-full-suite Phase-5D artifact each had SHA-256
`2606916915b2da5c051b1a06637c0a371a74751c6e77b2775bc629a63bc9f6dd`.

## Complete validation after isolation repair

All validation below was run from the committed isolation state. The complete
authenticated suite was restarted from the beginning rather than resumed.

| Gate | Actual result |
|---|---|
| Phase-5B governed regression | PASS, 1/1, raw rc 0 |
| Phase-5C governed regression | PASS, 1/1, raw rc 0 |
| Phase-5D governed regression | PASS, 1/1, raw rc 0; artifact SHA matched baseline |
| Complete data-free suite | PASS, 51/51, raw rc 0; zero failures, zero unexplained skips |
| Complete authenticated suite | PASS, 75/75, raw rc 0; zero failures, zero unexplained skips; 4716.50 s |
| Immediate post-full-suite Phase-5D regression | PASS, 1/1, raw rc 0; 2042.66 s; artifact SHA matched baseline |

The full suite ran `chemical_trackr_budget` first and the governed Phase-5D
regression immediately afterward. The chemical test passed in 95.96 seconds;
the fresh governed regression then passed in 2011.26 seconds instead of
refusing a dirty source tree. After the data-free suite, after the authenticated
suite, and after the immediate post-suite regression, `git status --porcelain`
was empty and explicit source-directory searches found no new `__pycache__`
directory or `.pyc` file.

All ten pre-existing baselines, the historical and promotion Phase-5D
candidates, the Phase-5C candidate, EOS/data, literature, 123 production
scientific paths, and all 33 protected paths remain byte-identical to entry.
In particular, `tests/analysis/chemical_trackr_budget.py`,
`CompactStar/Analysis/src/ParticleNumberResponse.cpp`, and
`CompactStar/Analysis/src/ChemicalResponse.cpp` are unchanged. The only
non-documentation change after governed baseline promotion is the CTest launch
isolation in `tests/CMakeLists.txt`.

## Frozen scientific classification and retained caveats

The governed baseline retains:

```text
benchmark_scope = CONTROLLED_MATHEMATICAL_ARCHITECTURE
classification = governed
candidate_only = false
governed_baseline = true
P0_s = 0.001
physical_spin_interpretation = false
super_kepler = true
envelope_provenance = FR2005 Eq. (49) / PCY97 fully accreted
```

The envelope correction is metadata-only and has zero numerical effect. The
historical ratified candidate remains untouched. The P0=1 ms history is
super-Kepler / physically inadmissible for the free-gas fixture during early
evolution and is governed only as the ratified mathematical/architecture
control, never as realistic pulsar evolution.

The remaining controlled-model caveats are unchanged: the Lother/PBF-only
interface; cancellation-sensitive net thermal power; scaled-RKF45
stability-bound late evolution and no solver claim beyond this benchmark;
heat-capacity cache kinks; Sommerfeld/free-gas limitations; frozen coefficients
with omitted dot(Z); no Direct Urca; no superfluidity; no crust/envelope chemical
dependence; prescribed rather than state-coupled spin; no physically admissible
spin-history claim; and no realistic A18/FR2005 normalization.

## INV-11 and future boundary

Upon canonical fast-forward, the controlled frozen-v1 dispositions are:

| Subpart | Scoped disposition |
|---|---|
| INV-11b | **RESOLVED FOR CONTROLLED FROZEN-V1 SCOPE** |
| INV-11c | **RESOLVED FOR CONTROLLED FROZEN-V1 SCOPE** |
| INV-11d | **RESOLVED FOR CONTROLLED FROZEN-V1 SCOPE** |
| INV-11e | **RESOLVED FOR FROZEN-V1 COEFFICIENT SEMANTICS** |
| INV-11f | **RESOLVED FOR CONTROLLED FROZEN-V1 ODE/SOURCE COUPLING** |

Global INV-11 remains **UNRESOLVED** for realistic A18/FR2005 normalization,
dot(Z), Direct Urca, superfluidity, crust/envelope chemical dependence,
state-coupled spin, physically admissible spin histories, and BNV. Realistic
FR2005/A18 is **SOURCE-LIMITED / BLOCKED**. No A18 implementation was performed,
and BNV is **NOT BEGUN**.
