# Phase-6 ADR-0018 owner acceptance

## Status and authority

**Status:**

```text
ADR-0018 HUMAN-RATIFIED
MAC SOURCE-BUILD EQUIVALENCE AUTHORIZED
DEPENDENCY SOURCE AUTHORITY NOT YET GRANTED
CLUSTER BOOTSTRAP NOT AUTHORIZED
```

On 2026-09-26 the human owner explicitly ratified ADR-0018 in its proposed form at documentation
branch SHA `a5926e5c5fc59a01240b4bd60765a6f640f19dd7`. The canonical entry is
`232565a32303a4953f3f516d1d5286b6663f8f99`; it is the exact ancestor of the documentation
branch. The acceptance commit is the commit containing this record with subject
`docs: ratify cross-platform dependency resolution`.

This is a documentation/governance acceptance step. It performs no build, test, ODE, dependency
promotion, archive replacement, cluster access, source transfer, Slurm action, CQ stage, clean
BA12R, scientific-model change, or C++ standard change.

## Accepted dependency-resolution decision

The owner accepts the following as one narrow architecture decision:

1. Existing authenticated Darwin Zaki and CONFIND archives remain unchanged and remain the
   dependency authority for existing governed Mac evidence.
2. Cross-platform builds use separately authenticated external static dependencies supplied
   through explicit fail-closed include/library overrides.
3. No implicit system-library or header search is permitted.
4. No build-time network fetch is permitted.
5. An exact external Git SHA, or a complete immutable source-manifest identity where no
   defensible Git revision exists, together with consumed-header hashes, archive hashes,
   toolchain, flags, and build configuration enters qualification provenance.
6. Any dependency identity change invalidates the applicable platform qualification.
7. Compilation alone does not establish scientific authority.
8. A candidate source must pass same-Mac source-build equivalence before it may be considered as
   Linux dependency source authority.
9. ADR-0018 does not replace or re-ratify the current vendored Darwin archives.

The accepted override names are `COMPACTSTAR_ZAKI_LIBRARY`,
`COMPACTSTAR_ZAKI_INCLUDE_DIR`, `COMPACTSTAR_CONFIND_LIBRARY`, and
`COMPACTSTAR_CONFIND_INCLUDE_DIR`. Their implementation is authorized only in the separately
isolated experiment branch after canonical integration of this record.

## Accepted equivalence candidates

The owner accepts the following identities strictly as **SOURCE-EQUIVALENCE CANDIDATE ONLY**:

| Dependency | Candidate identity | Current authority status |
|---|---|---|
| ZakiLib | Git SHA `b9ddebaded24962468954846f47238aec2726fd4` | candidate only; not governed Linux authority |
| CONFIND | complete source-manifest SHA-256 `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083` | candidate only; not governed Linux authority |

The recovered CONFIND candidate is an exact non-Git composite snapshot preserved at
`/Users/keeper/Documents/CompactStar/external/recovered/CONFIND/ed76163c22e0a1f8/`. Its transfer
ZIP `CONFIND_RECOVERED_SOURCE_ed76163c22e0a1f8.zip` has SHA-256
`be41ee1c71b88627f0b11703c9d732abd7378ebbc162ee10649af771dec9ae0f`. Recovery and preservation
evidence are in `PHASE6_CONFIND_SOURCE_RECOVERY.md` and
`PHASE6_CONFIND_SOURCE_PRESERVATION.md`. No historical CONFIND Git SHA is assigned.

## Existing Darwin authority

The authenticated vendored artifacts remain unchanged:

| Artifact | SHA-256 |
|---|---|
| `dependencies/lib/Zaki/Darwin/arm64/libZaki.a` | `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef` |
| `dependencies/lib/Confind/Darwin/arm64/libConfind.a` | `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2` |

They are the control dependencies for the authorized Mac experiment and remain the dependency
authority for current governed Mac evidence regardless of that experiment's outcome.

## Authorized Mac-only experiment

After fast-forward-only canonical integration, the owner authorizes a fresh noncanonical branch
to:

- predeclare the complete C0/C1/T experiment and equality classes before any build;
- add the minimum CMake-only fail-closed dependency override mechanism;
- build the exact external candidates outside their source repositories and CompactStar's
  vendored dependency paths without changing source bytes;
- run focused dependency contracts, required Phase-5B/C/D regressions, and the bounded ADR-0017
  qualification; and
- compare C0 versus C1 and C1 versus T under predeclared same-platform rules.

C0 is canonical CompactStar with implicit vendored dependencies. C1 is the experiment source
with explicit overrides pointing to the same vendored dependencies. T differs from C1 only by
using source-built candidate archives and candidate headers. Deterministic equality-bearing
scientific fields require exact/0-ULP equality; existing analytic/reference tests retain their
own governed acceptance tolerances but cannot hide a treatment/control difference.

## Status synchronization and exclusions

- ADR-0015: **ACCEPTED**.
- ADR-0016: **ACCEPTED**.
- ADR-0017: **ACCEPTED / HUMAN-RATIFIED**.
- ADR-0018: **ACCEPTED / HUMAN-RATIFIED**.
- Zaki source SHA: **EQUIVALENCE CANDIDATE ONLY**.
- CONFIND source manifest: **EQUIVALENCE CANDIDATE ONLY**.
- Mac source equivalence: **AUTHORIZED / NOT YET RUN**.
- EKU cluster: **DISCOVERED / NOT QUALIFIED**.
- CQ0: **BLOCKED**.
- Historical BA12: **FAIL**.
- Historical BA12R: **FAIL**.
- Clean BA12R: **NOT AUTHORIZED**.

The owner does not authorize candidate promotion, vendored-archive replacement, a C++ standard
change, cluster build or access, manual cluster transfer, cluster bootstrap, CQ0-CQ7, Slurm,
clean BA12R, scientific-model change, baseline change, or unrelated refactoring. The separate
C++17 `using enum` issue remains unresolved and outside this decision.

## Integration gate

This acceptance delta must remain documentation-only. After a non-force branch push and parity
check, canonical `master` may advance from
`232565a32303a4953f3f516d1d5286b6663f8f99` to the acceptance commit by fast-forward only. No
merge commit, squash, cherry-pick, rebase, or force push is authorized.

Successful canonical integration authorizes only creation of the fresh Mac experiment branch and
execution of the bounded experiment above. Even a passing result must return the candidate
identities to the owner for a separate Linux/cluster source-authority decision.
