# ADR-0020: Explicit native Linux dependency candidate

Status: PROPOSED FOR INDEPENDENT REVIEW; implementation authorized by the owner's
2026-10-10 Linux-transition task, not scientifically ratified or integrated.

Dependency/build and structural changes extend ADR-0019's explicit external
package authority to native Linux. Darwin retains the byte-identical Mac lock.
Linux must select an absolute `COMPACTSTAR_LINUX_DEPENDENCY_LOCK`; absent locks,
wrong platform/mode/source/compiler/flags, changed package bytes or changed
GSL/zlib/runtime/toolchain closure are fatal before package CONFIG files execute.
Cross compilation and implicit package fallback are unsupported. Exact Zaki
2.0.1 and CONFIND 2.0.0 source identities remain those accepted by ADR-0019.

A lock is an integrity record for a candidate build, not a signature, numerical
acceptance, or campaign authorization. The versioned isolated prefixes are
reviewed with this candidate; changing a prefix or its content requires new
provenance and requalification. Absolute paths intentionally bind this cluster
installation. A relocatable design would need explicit reauthentication.

Alternatives: overwriting the Mac lock would destroy accepted authority; accepting
any package with a matching version would omit source and arithmetic provenance;
vendoring Linux binaries would blur external dependency ownership. Retain the
separate externally built packages and authenticate before loading CMake code.

`Tags.hpp` only qualifies existing enum names to honor C++17. Numerical equations,
constants, baselines and failure records are unchanged. Linux harness adapters
retain the historical Mac verification paths; they must not counterfeit Mac
platform evidence. Cross-platform comparison follows the predeclared criteria in
`docs/validation/linux-transition/PREDECLARATION.md:1`. A combined independent
scientific/numerical review precedes any canonical integration or long campaign.
