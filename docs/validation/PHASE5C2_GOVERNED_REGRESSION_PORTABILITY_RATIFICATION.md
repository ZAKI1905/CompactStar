# Phase-5C-2INT-R governed-regression portability ratification

**Date:** 2026-09-09

**Decision:** **HUMAN-RATIFIED — THE COMPILER-VERSION STRING IS RETAINED EXECUTION
PROVENANCE, NOT A CROSS-TOOLCHAIN SCIENTIFIC-IDENTITY EQUALITY FIELD.**

This decision is deliberately narrow. It permits no numerical tolerance, no compiler-string
normalization, and no exemption for any field other than the exact JSON path
`provenance.toolchain.compiler`.

## 1. Authenticated authority and previous stop

| Item | Authenticated value |
|---|---|
| Canonical `master` | `49ab2b8c2881b6ef7b9309307d18cea51d557f72` |
| Human-ratified Phase-5C head | `27727016856a6a25a46e447c70e380722ea8ddbf` |
| Reviewed candidate | `docs/validation/phase5c_chemical_coefficients_candidate.json` |
| Reviewed candidate SHA-256 | `a6b430b2356987b54a7c82c87d0b00f3e1d9698f6f299c32c80efcc40c63833b` |
| Reviewed compiler | `Apple LLVM 17.0.0 (clang-1700.6.4.2)` |
| Fresh current compiler | `Apple LLVM 21.0.0 (clang-2100.1.1.101)` |
| Previous disposition | `PHASE-5C GOVERNED REGRESSION / BASELINE VALIDATION FAILED — CORRECTION REQUIRED` |

The previous integration attempt stopped before baseline installation, commit, push, or master
movement. That stop was correct under its then-controlling rule: only governance/classification
metadata was allowed to differ, while the fresh compiler string differed at the non-classification
path `provenance.toolchain.compiler`.

The preserved artifacts were reauthenticated before this decision. Deep comparison found the
reviewed candidate and fresh current-toolchain output identical in every equality-bearing field:
`G`, `Q`, `Z`, `I`, `W`, every numerical-error and validation-envelope entry, certificates,
support/rank data, source and physics metadata, goals, constants, methods, domains, partitions,
and refusal/onset/tail policies. The only additional non-classification difference was exactly
`provenance.toolchain.compiler`. The candidate remained byte-identical at its reviewed SHA.

## 2. Governing-text finding

`GOVERNANCE.md` requires recorded toolchain versions and reproducible build instructions for a
dependency/build change; it does not declare compiler-version equality to be scientific identity
(`GOVERNANCE.md:53`). The implementation record requires the actual compiler string to be retained
and says that the candidate records the toolchain reported by its executable
(`docs/validation/PHASE5C2_CORRECTED_CHEMICAL_COEFFICIENT_IMPLEMENTATION.md:147`,
`:215-216`). The predeclaration required byte-identical candidate reproduction in the reviewed
release environment; it did not define compiler-version equality for all future governed runs
(`docs/validation/PHASE5C2_PRODUCTION_ACCEPTANCE_PREDECLARATION.md:162-165`). No ACCEPTED rule
explicitly makes the compiler-version string a cross-toolchain scientific-identity field.

Recording provenance and requiring all provenance fields to compare equal are distinct operations.
The former remains mandatory. This decision governs the latter narrowly for one field.

## 3. Ratified provenance classification

### 3.1 Scientific/computational identity provenance

Equality-bearing identity includes governed source-code/version identity; EOS/provider bytes;
physical constants; basis and units; metric convention; domain and partition policy; quadrature
method/order; support, onset, refusal and tail policies; accuracy goals; structural-source
identity; predeclaration identity; linked-library identity; language level; architecture;
platform; configuration; and every scientific value, error, envelope and certificate. Any
unexplained difference in these fields is a regression failure.

### 3.2 Execution-environment provenance

For this decision the sole portable execution-provenance field is
`provenance.toolchain.compiler`. It must exist, be nonempty, and contain the compiler identity and
version truthfully emitted by the build/runtime provenance mechanism. Its value is preserved and
reported before comparison. It may differ across toolchains without failing the artifact
comparison only when every equality-bearing field is exactly equal.

The exact portable-field allowlist is therefore:

```text
provenance.toolchain.compiler
```

This decision does not generalize to Python, operating system, architecture, CMake, GSL or other
library versions, build configuration, C++ language level, or any other toolchain/environment
field. A difference in any such field stops for separate adjudication.

## 4. Comparator and failure contract

The comparator must retain and report both compiler strings, compare exact structured values with
no floating tolerance, and exclude only the single path above. A changed compiler string never
excuses a changed number. A compiler difference plus any scientific, numerical, source,
scientific-provenance, or other environment difference is a regression failure.

Required negative controls change a scientific number, a scientific source/provenance field, an
environment field other than the compiler, remove the compiler, and empty the compiler; every one
must fail. A compiler-only mutation must pass while reporting the difference. The exclusion list
must not be widened after observing a failure.

The observed LLVM17-to-LLVM21 equality-bearing payload identity is positive portability evidence.
It is not proof of compiler independence and creates no numerical tolerance.

## 5. Explicit prohibitions and scope

- Do not normalize, fabricate, remove, or copy an earlier compiler string into a fresh artifact.
- Do not install, recover, or pin Apple LLVM 17 merely to reproduce provenance text.
- Do not widen scientific tolerances or artifact-comparison exclusions.
- Do not change `G`, `Q`, `Z`, `I`, `W`, error/envelope semantics, support/rank semantics,
  onset/refusal/tail machinery, lifetime/staleness semantics, or acceptance goals.
- Do not change the reviewed candidate artifact.
- Do not implement eta evolution, weak rates, heating/cooling, A18, Phase-5D, or BNV.

This ratification changes no Phase-5C physics. At this policy commit no governed Phase-5C baseline
has yet been installed. Baseline generation and canonical integration remain subsequent gated
steps.
