# Phase-5C corrected chemical coefficients — canonical integration record

**Date:** 2026-09-10

**Status:** **PHASE-5C CORRECTED CHEMICAL COEFFICIENTS — CANONICAL INTEGRATION
VALIDATED / GOVERNED REGRESSION INSTALLED / READY TO CLOSE PHASE-5C COEFFICIENT
SCOPE.**

This record closes only the accepted generic/free-gas coefficient scope. It does not implement
secular evolution or expand the scientific authority ratified for Phase-5C.

## 1. Authenticated history and isolation

| Item | Authenticated value |
|---|---|
| Canonical entry `master` | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Human-ratified Phase-5C head | `27727016856a6a25a46e447c70e380722ea8ddbf` |
| Phase-5C portable-provenance policy | `2ec126d27920689de9488ead2a3cafafce6dbbb7` |
| Phase-5C governed-regression commit | `dbc6422e4ebe7347940e9d948e7dcbdd838ce514` |
| Phase-5B portable-provenance repair | `bcd6e0e498ab1acfee1077a056c7f4ebee22e122` |
| Integration branch | `integration/phase5c-corrected-chemical-coefficients` |
| Integration worktree | `/Users/keeper/Documents/CompactStar/worktrees/CompactStar-phase5c-integration` |

The canonical entry is an ancestor of every listed integration commit. The integration worktree
was created at the human-ratified Phase-5C head and remained isolated from canonical `master`
through all pre-master gates. The separately ratified Phase-5D branch and ADR-0014 history were
not imported.

## 2. Reviewed candidate and governed artifact

| Artifact | Path | SHA-256 |
|---|---|---|
| Reviewed candidate | `docs/validation/phase5c_chemical_coefficients_candidate.json` | `a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b` |
| Governed baseline | `tests/baselines/phase5c_chemical_coefficients.json` | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |

The reviewed candidate remained byte-identical. The governed artifact is separate and truthfully
classified with `candidate_only = false`, `governed_baseline = true`, and the governed
classification string. Its compiler is the actual producing compiler, Apple LLVM 21.0.0
(`clang-2100.1.1.101`); it does not copy the reviewed candidate's Apple LLVM 17.0.0
(`clang-1700.6.4.2`) string.

The exact reviewed-candidate to governed-baseline difference allowlist is:

1. `candidate_only`
2. `governed_baseline`
3. `classification`
4. `provenance.toolchain.compiler`

The exact fresh-governed-run to governed-baseline portability allowlist is only
`provenance.toolchain.compiler`. No scientific field is excluded and no numerical tolerance is
used. Deep comparison found the candidate and governed baseline identical everywhere outside the
four listed paths. Thus the independently produced LLVM17 and LLVM21 artifacts have identical
equality-bearing scientific and provenance payloads. This is positive portability evidence, not
proof of compiler independence.

## 3. Fresh producer and independence

The canonical producer is `tests/analysis/produce_chemical_coefficient_reference.py`; the
regression is `tests/analysis/chemical_coefficient_regression.py`. The producer creates the
Structure-1 fixture in fresh scratch, invokes the actual production chemical-coefficient
executable, and obtains `G_y`, `Q`, `Z`, `I_phys`, `W`, numerical errors, validation envelopes,
certificates, support/rank information, and onset/refusal/tail evidence from production code and
governed configuration/source authority. It accepts a caller-supplied output directory.

Source inspection and the registered regression establish that production executes before the
baseline is loaded for comparison. The producer does not read the governed baseline, the reviewed
candidate as numerical authority, or any historical scratch coefficient result. Expected result
bytes cannot seed production. The regression's scientific-field mutation must fail, which also
demonstrates that the comparison is not self-comparison.

Three isolated Apple-LLVM-21 generations were completely byte-identical:

| Generation | SHA-256 |
|---|---|
| 1 | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| 2 | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| 3 | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| Post-governed-commit regeneration | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |
| Final pre-master regeneration | `7027aa6179fe9111d76586d467d007c4223ec9d26ece0fee773327973c77bfe7` |

Every listed governed generation is byte-identical to the committed baseline. Comparator controls
pass for compiler-only difference and fail as required for missing compiler, empty compiler,
altered scientific value, altered scientific source/provenance, and any other toolchain field.
Compiler provenance is retained and reported before comparison.

## 4. Two fail-closed portability stops

The first integration attempt correctly stopped when the reviewed LLVM17 Phase-5C candidate and
fresh LLVM21 artifact differed at `provenance.toolchain.compiler` under the then-controlling
classification-only allowlist. No other field differed. The human-ratified narrow policy is
recorded in `docs/validation/PHASE5C2_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md` and precedes
the governed-baseline commit.

The resumed data-free suite then correctly stopped at 44/45, raw rc 8, because the historical
LLVM17 Phase-5B baseline and fresh LLVM21 output differed only at
`provenance.build.compiler`. The exact one-field Phase-5B portability decision is recorded in
`docs/validation/PHASE5B_GOVERNED_REGRESSION_PORTABILITY_RATIFICATION.md` and implemented at
`bcd6e0e498ab1acfee1077a056c7f4ebee22e122`. The historical Phase-5B baseline remained unchanged
at SHA-256 `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa`.

Two separate fresh current-compiler Phase-5B generations both had SHA-256
`05fe5e87689df6a84312db1a15f1b768e613373a9c18bf4d31da9579da96bee6` and were raw-byte
identical. Comparison with the historical baseline passed after excluding only
`provenance.build.compiler`; every equality-bearing Phase-5B field matched. The Phase-5B controls
proved same-compiler raw-byte enforcement and rejection of missing/empty compiler, changed
coefficient, numerical error, EOS SHA, architecture, build configuration, other provenance, and
compiler-plus-science mutations. Both stops remain visible as valid fail-closed outcomes.

## 5. Fresh validation inventory

All tests used the authenticated external EOS/provider data where required.

| Gate | Inventory/result | Raw rc | Failures/skips |
|---|---:|---:|---|
| Focused Phase-5C, including governed regression | 6/6 PASS | 0 | none |
| Complete data-free (`-LE external-data`) | 45/45 PASS | 0 | none |
| Complete authenticated repository suite | 68/68 PASS | 0 | none |
| Focused Phase-5B governed regression after repair | 1/1 PASS | 0 | none |

No failure was absorbed and no unexplained skip occurred. The current CTest inventory, rather
than a hard-coded historical count, supplied these totals.

## 6. Immutability audit

The baseline count was 9 on canonical entry and is 10 with the new Phase-5C governed artifact.
All original nine baseline files are byte-identical to entry. The Phase-5B historical baseline,
EOS/data files, literature files, and reviewed Phase-5C candidate are byte-identical to entry.
There is no Phase-5C production-source or central-value change beyond the already accepted
candidate history. `G_y`, `Q`, `Z`, `I_phys`, `W`, both UQ tracks, support/rank semantics,
onset/refusal/tail machinery, lifetime/staleness semantics, and all predeclared acceptance goals
remain frozen. Python bytecode was confined to confirmed scratch; no bytecode or housekeeping
change is committed.

## 7. Scientific status and retained caveats

ADR-0013 is **ACCEPTED / IMPLEMENTED / INTEGRATED** for the governed generic/free-gas coefficient
scope. GC1-GC12, GC9b, and GC14 are **PASS** under their ratified evidence classifications. GC13
is **SOURCE-LIMITED / BLOCKED**. INV-09 remains **VERIFIED / RESOLVED** for its governed structural
scope. Global INV-11 remains **UNRESOLVED**.

The following caveats remain controlling:

- NB-1: the current single-enum `PaperZ` axis representation is accepted only because the
  symmetric two-channel map is one-to-one; asymmetric or expanded input/output spaces require
  semantically distinct axis typing.
- NB-2: `ChargeNeutralNumberSusceptibility::NumericalError()` is local numerical
  solve/congruence uncertainty, not total provider/model uncertainty; provider/background effects
  remain accounted globally.
- NB-3: Python bytecode and `__pycache__` creation is housekeeping only.
- The validation envelope is conservative and PB11-constraint dominated; it is not achieved
  accuracy. The PB7 micro-difference is immaterial and the frozen envelope is not post-hoc refit.
- GC9 M20 is an equivalent transformed-input mutant, not a literal mutation of the production
  weight expression.
- Current refusal-edge extrema reasoning relies on the cellwise-linear background representation.
- The old-M tail mutant has a thin but independently positive separation; outward-rounding
  discipline remains required.
- GC13/A18 remains source-limited and blocked; global INV-11 remains unresolved.

No eta evolution, weak rates, heating/cooling or thermal rotochemical coupling, realistic A18
closure, superfluidity, Phase-5D implementation, or BNV is implemented or begun. The separately
human-ratified Phase-5D ADR-0014/preflight branch remains out of this history and must be reconciled
in a separate governed task.
