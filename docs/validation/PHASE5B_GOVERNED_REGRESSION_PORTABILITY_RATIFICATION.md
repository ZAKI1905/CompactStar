# Phase-5B governed-regression compiler-portability ratification

**Date:** 2026-09-10

**Decision:** **HUMAN-RATIFIED — `provenance.build.compiler` IS RETAINED
EXECUTION-ENVIRONMENT PROVENANCE, NOT A CROSS-COMPILER SCIENTIFIC-IDENTITY
EQUALITY FIELD.**

This decision is restricted to the exact JSON path `provenance.build.compiler`. It changes no
Phase-5B artifact, scientific value, uncertainty, tolerance, source identity, or other
provenance requirement.

## 1. Authenticated authority and integration context

| Item | Authenticated value |
|---|---|
| Canonical `master` | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Governed Phase-5B baseline | `tests/baselines/phase5b_structural_response.json` |
| Governed baseline SHA-256 | `7588f0e9cd62f5b6be48bb725e4d0ba6e47b64f50d843117e70c17a78beeb5fa` |
| Current integration branch | `integration/phase5c-corrected-chemical-coefficients` |
| Phase-5C integration entry | `dbc6422e4ebe7347940e9d948e7dcbdd838ce514` |
| Interrupted data-free result | `44/45 PASS`, raw CTest rc `8` |
| Sole failed test | `phase5b_structural_response_regression` |
| Baseline compiler | `Apple clang version 17.0.0 (clang-1700.6.4.2)` |
| Current compiler | `Apple clang version 21.0.0 (clang-2100.1.1.101)` |
| Exact differing path | `provenance.build.compiler` |

The stopped current-toolchain artifact had SHA-256
`05fe5e87689df6a84312db1a15f1b768e613373a9c18bf4d31da9579da96bee6`.
Deep comparison against the historical governed baseline found exactly one difference: the
truthfully recorded compiler string above. Every coefficient (`A`, `B`, `K`, `I_phys`), count,
numerical error, summary value, contributor, tail value, unit, fixture/domain field, EOS identity,
EOS revision and bytes, source hash, architecture, build configuration, and other provenance field
was exactly equal.

The old regression correctly failed under its then-controlling raw-byte rule. The accepted
ADR-0011 and Phase-5B implementation, ratification, integration, and invariant records require
deterministic artifacts and retained provenance; none explicitly declares compiler-version
equality to be scientific identity across compilers. Mandatory provenance retention and mandatory
cross-environment value equality are distinct requirements.

## 2. Ratified provenance classification

`provenance.build.compiler` is execution-environment provenance. It must be present, be a nonempty
string, record the actual compiler used, and be retained and reported before comparison. It may
not be normalized, fabricated, removed, or copied from a historical artifact.

Every other field is equality-bearing. This includes `provenance.build.architecture`,
`provenance.build.configuration`, EOS identity/revision/table bytes, source identity, Hartle
inputs, fixture, domain, species, contributors, tails, units, every central value, and every
numerical error. The complete and exclusive portability allowlist is:

```text
provenance.build.compiler
```

A compiler-only difference passes while both exact compiler strings are reported. A compiler
difference accompanied by any other difference fails. No numerical tolerance or fuzzy comparison
is introduced, and a compiler change never excuses scientific drift.

## 3. Governed comparator contract

Both artifacts must parse as JSON and carry present, nonempty compiler strings. If those strings
are identical, complete raw-byte identity remains mandatory. If they differ, the comparator
retains both values and requires exact parsed equality after replacing only the two compiler values
in comparison copies. No command-line allowlist and no generic provenance exclusion are permitted.

Negative controls must reject missing/empty compiler provenance; numeric-coefficient,
numerical-error, EOS-table-hash, architecture, configuration, other-provenance, and
compiler-plus-science mutations. Same-current-compiler independent generations remain subject to
complete raw-byte identity.

The observed Apple-LLVM-17 to Apple-LLVM-21 equality-bearing payload equality is positive
portability evidence, not proof of compiler independence.

## 4. Scientific and governance disposition

The historical governed Phase-5B baseline remains byte-identical and authoritative. This decision
does not alter `A`, `B`, `K`, `I_phys`, any error, any baseline byte, any Phase-5B physics, or any
accepted evidentiary caveat. INV-09 remains **VERIFIED / RESOLVED** for its governed structural
scope; global INV-11 remains **UNRESOLVED**.
