# ADR-0019 — Canonical Source-Built Zaki and CONFIND Dependency Authority

| Field | Value |
|---|---|
| **Status** | ACCEPTED for the owner-authorized Mac migration candidate; canonical integration awaits one combined independent review and owner acceptance |
| **Date** | 2026-09-29 |
| **Change class** | dependency/build + structural |
| **Governing authority** | Owner OD1–OD9 and resume authorization; ADR-0017; successor to the active-dependency portions of ADR-0018 |
| **Affected invariants** | Existing numerical contracts and governed baselines unchanged |
| **Blocks** | Canonical integration of this migration candidate |

## Context

Historical CompactStar consumed bundled Darwin archives and plotting headers.
The plotting bridge introduced Python/NumPy linkage and permissive undefined
symbol resolution. The owner authorized exact source-built successors and
same-mode local Mac scientific equivalence as the acceptance criterion.
ADR-0018 remains unchanged and valid for its historical scope and evidence.

## Decision

Active builds consume `Zaki::Zaki` from source
`e263a6e180c5c417198e7778bd21fc9c0a32dc33` (2.0.1) and `CONFIND::CONFIND`
from source `b0cbd510fd3fd0c772fa50499cd749287cb39e7b` (2.0.0).
Both absolute package prefixes are explicit and required. CONFIG discovery
has no default search path; archive, installed-header, package-config,
source identity and build mode are checked before loading package targets.
The candidate lock records the qualified Mac artifacts. Rebuilding packages
requires explicit provenance and requalification, not silent hash updates.

Historical headers and archives stay at their existing paths as immutable
oracles. They are excluded from active include/link/install resolution:

- Zaki arm64 SHA-256: `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef`.
- CONFIND arm64 SHA-256: `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2`.

C++ plotting and Python/NumPy linkage are removed. Python remains an
interpreter for validation. GSL remains a direct CompactStar dependency;
zlib comes from Zaki and portable threading uses `Threads::Threads`.
Unused OpenMP linkage is removed. TaskManager's implementation is unchanged;
one thread governs this migration. Useful derived plot-only data is exported
as numerical tables; visualization is external.

## Alternatives

Retaining active vendored archives would preserve the obsolete plotting/link
chain. Unpinned package search could substitute scientifically different
arithmetic. Requiring identical provider ownership or machine code would
reject benign compiler/template changes without establishing scientific
behavior. The chosen gate is exact deterministic numerical/output equality,
with symbol/disassembly inspection as diagnostic evidence.

## Consequences and validation

OLD Debug versus NEW Debug and OLD Release versus NEW Release are separate
authorities. T1 ordering, selection and replay records and all 38 T2 outputs,
including three species' BNV lifetime tables, must match exactly. Existing
Phase-5/6 baseline rules remain intact, including Debug-only governed
certificates. No new tolerance is authorized. Provider differences are
acceptable only when they do not change scientific values, branching,
ordering, solver inputs, serialization or memory safety in the qualified paths.
The migration report records actual results and remaining limitations.

This is a local Mac dependency candidate. It establishes no Linux, cluster,
Windows, cross-platform bitwise, or multithreaded TaskManager authority.
TaskManager concurrency redesign and future CONFIND evaluator parallelism
remain separate work. No merge or release tag is authorized here.

## Provenance

The human owner selected the exact sources, authorized this successor ADR
(OD3), and authorized implementation before a single combined independent
review. The agent implemented and documented the decision. The next step is
that combined review of Zaki 2.0.1, CONFIND 2.0 and CompactStar; integration
requires subsequent owner acceptance.
