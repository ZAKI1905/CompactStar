# ADR-0018: Cross-platform external dependency resolution

## 1. Status, scope, and authority

**Status:** ACCEPTED / HUMAN-RATIFIED.

**Date proposed:** 2026-09-26.

**Date accepted:** 2026-09-26.

**Change class:** dependency/build + structural architecture.

**Blocks:** Linux dependency bootstrap LB1-LB10 and cluster qualification CQ0-CQ7.

The human owner ratified this decision and separately authorized the bounded Mac-only
source-build equivalence experiment recorded in
`docs/validation/PHASE6_ADR0018_ACCEPTANCE.md`. Acceptance does not itself change CMake, build a
dependency, select a Linux artifact, replace a Darwin archive, grant dependency source authority,
or authorize cluster work.

The decision is deliberately limited to resolution and provenance of ZakiLib and CONFIND.
It changes no scientific equation, numerical method, Phase-5 baseline, Phase-5 ownership,
EOS/data input, or ADR-0017 contract. The evidence and complete audit are in
`docs/validation/PHASE6_LINUX_CLUSTER_BOOTSTRAP_PREFLIGHT.md`.

## 2. Context

CompactStar currently derives both dependency archive paths from
`CMAKE_SYSTEM_NAME` and `CMAKE_HOST_SYSTEM_PROCESSOR`, stores them in ordinary variables, and
fails when the expected files do not exist (`CMakeLists.txt:89-101`). Both archives and one
shared dependency include root are linked/exposed by the `CompactStar` target
(`CMakeLists.txt:140-180`). This preserves a deterministic Darwin build but makes canonical
Linux configuration impossible because no Linux archives exist at the hard-coded locations
and command-line values cannot override the ordinary `set()` assignments.

The authenticated Darwin artifacts are governed inputs:

- `libZaki.a` SHA-256
  `3dd4789a20c35064b3133bb863c54c4f64e7df31c94b83201d68f5463902dfef`;
- `libConfind.a` SHA-256
  `09ed1a7c43a83b42f64ee8e0bda3b879af970126ee75159a802179a4d0a49eb2`.

No durable record maps either archive to an exact source SHA. The current Zaki source candidate
has byte-identical consumed headers. The current public/local CONFIND repository does not: its
`.h`/ROOT-facing interface lacks methods and types required by CompactStar's vendored `.hpp`
contract. Consequently ADR ratification alone cannot nominate CONFIND source or begin LB1.

The transferred EKU discovery found GCC 14.1.0 and CMake 3.29.4 but only GSL 2.6, and confirmed
that the present hard-coded dependency layout blocks a Linux configure. That evidence remains
historical and immutable in
`docs/validation/PHASE6_EKU_CLUSTER_QUALIFICATION_PREFLIGHT.md:142-190`.

## 3. Accepted decision

The owner accepts all of the following as one narrow decision.

1. The existing authenticated Darwin `libZaki.a` and `libConfind.a` remain unchanged and remain
   the default Mac dependency artifacts. No source-built replacement is implied.
2. CompactStar will support explicit, fail-closed dependency-path overrides named:
   `COMPACTSTAR_ZAKI_LIBRARY`, `COMPACTSTAR_ZAKI_INCLUDE_DIR`,
   `COMPACTSTAR_CONFIND_LIBRARY`, and `COMPACTSTAR_CONFIND_INCLUDE_DIR`.
3. On Darwin, each unset override defaults to the exact existing vendored archive/include
   location. Thus an ordinary Mac configure retains present behavior.
4. On non-Darwin platforms, all four paths must be supplied explicitly. They must be absolute,
   exist, be the expected file/directory kind, and contain expected sentinel headers. A missing
   or invalid path is a configure-time fatal error.
5. CompactStar must not use `find_library`, `find_path`, uncontrolled prefix/environment
   probing, or implicit system fallback for Zaki/CONFIND.
6. CompactStar and its scientific build must perform no `FetchContent`, submodule update,
   package download, or other build-time network access for these dependencies.
7. External dependency binaries are built separately, outside CompactStar Git, from exact
   human-ratified source identities: a Git commit SHA when available, or an immutable complete
   source-manifest SHA-256 when no defensible Git revision exists. Cluster artifacts live under
   versioned toolchain storage and are consumed read-only.
8. Qualification provenance includes the dependency Git SHA or complete source-manifest
   SHA-256, every
   consumed/installed header hash, archive hash and member/symbol inventory, compiler and
   version, language standard, flags, build system/version, link closure, and resolved absolute
   paths.
9. Any dependency source identity or member byte, consumed-header byte, archive byte,
   compiler/toolchain key,
   or declared link-closure change invalidates the applicable platform qualification until
   requalified.
10. Before a source identity can become Linux cluster authority, it must pass the separately
    predeclared same-Mac source-build equivalence gate against CompactStar using the current
    authenticated vendored Darwin archives.
11. Replacing, rebuilding, or re-ratifying the current Darwin vendored archives is explicitly
    outside ADR-0018.

Resolved paths must be printed during configuration. Hashing remains a qualification-tooling
responsibility rather than an implicit CMake filesystem search; the provenance manifest must
bind the printed paths to authenticated hashes before a build may qualify.

## 4. Source-authority gate

ADR-0018 governs how an accepted dependency is resolved. It does not decide which source
identity is scientifically acceptable.

- Zaki `b9ddebaded24962468954846f47238aec2726fd4` is a **candidate only** for Mac
  equivalence because all consumed headers are byte-identical. It is not asserted to be the
  historical archive source.
- CONFIND `89c5d9b731534e4289d9f686549d9f0ac178e567` is **not a candidate**. Its public interface
  materially diverges from the current CompactStar dependency contract. Source for the matching
  interface has since been recovered as a composite non-Git snapshot identified by complete
  source-manifest SHA-256
  `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`; recovery evidence and
  limitations are in `docs/validation/PHASE6_CONFIND_SOURCE_RECOVERY.md`. It remains a candidate
  only and requires owner preservation and authorization before equivalence testing.

The owner must separately approve each exact candidate source identity for the Mac experiment
and, only after it passes, separately accept it as cluster authority. Source existence, build
success, a manifest, or current repository HEAD is never ratification.

## 5. Same-Mac equivalence requirement

The mandatory control is canonical CompactStar built with the unchanged vendored Darwin
archives. The treatment is the same CompactStar source, compiler, SDK, configuration, inputs,
and environment, but with separately built candidate archives/headers supplied by the explicit
overrides. Vendored files are never overwritten.

Acceptance requires:

- byte-identical consumed headers;
- exact caller-visible API and normalized exported ABI;
- focused Zaki constants/conversions, numerical, data, string, and file behavior;
- focused CONFIND grid/contour/data-export behavior and no ROOT requirement on the consumed
  path;
- unaffected governed Phase-5B, Phase-5C, and Phase-5D regressions;
- unaffected bounded ADR-0017 Phase-6 qualification;
- exact deterministic scientific payload, binary64 fields, accepted-step histories,
  observation/checkpoint schedules, counts, and validity/provenance decisions between control
  and treatment.

Archive-container timestamps, object paths, build IDs, and explicitly declared provenance
fields may differ. No new numerical cross-build tolerance is authorized: deterministic
control/treatment scientific fields require 0 ULP. Existing analytic/reference test tolerances
remain their own acceptance rules but cannot mask control/treatment drift. A candidate that
fails is rejected; tolerance cannot be selected after observing results.

## 6. Alternatives considered

### Alternative A — integrate dependency source into CompactStar

Submodules, `FetchContent`, or `add_subdirectory` would let one configure build everything.
They would also couple unrelated build systems, blur source ownership, complicate license and
toolchain provenance, and risk build-time network behavior. Rejected as the recommendation.

### Alternative B — authenticated external static archives

Build exact accepted sources separately, store them under versioned platform/toolchain keys,
and pass explicit paths. This preserves current dependency boundaries and Darwin authority,
supports read-only reuse, and makes invalidation inputs inspectable. **Recommended.**

### Alternative C — commit prebuilt Linux archives to CompactStar

This would imitate the current Darwin layout but grow platform binaries in Git and make source
and toolchain authority easier to lose. Rejected.

### Alternative D — implicit package/system discovery

This reduces configuration arguments but permits accidental ABI/version fallback and weakens
reproducibility. Rejected and explicitly forbidden by the proposed decision.

## 7. Consequences

- A later bounded implementation may change only the dependency-resolution portion of CMake
  and its build documentation/tests; this ADR itself implements nothing.
- Current Mac defaults and governed results remain valid because their input artifacts remain
  byte-identical and selected by default.
- Linux builds become possible only with explicit authenticated paths and complete provenance.
- Dependency builds, GSL/Python environments, source acquisition, and network behavior remain
  outside the CompactStar scientific configure/build.
- A source or toolchain change is a new qualification key, not an in-place update.
- Model B becomes the prescribed boundary, but owner preservation and acceptance of the
  recovered CONFIND snapshot remain prerequisites.
- Cluster CQ0-CQ7 remain blocked until LB0-LB9 are completed and their evidence is accepted.

## 8. Validation required before implementation and use

Before a CMake implementation can be accepted:

1. owner ratifies ADR-0018;
2. tests demonstrate unchanged Darwin default resolution and fail-closed non-Darwin behavior;
3. tests demonstrate no implicit search or network path;
4. configure evidence prints the four resolved paths;
5. protected production, test, baseline, EOS/data, and literature bytes outside the authorized
   implementation remain unchanged.

Before Linux artifacts can qualify, the source-authority and same-Mac equivalence gates in
sections 4-5 must pass, then the cluster bootstrap must record exact GSL, Python, compiler,
header, archive, and source provenance. No CQ stage may compensate for a skipped predecessor.

## 9. Non-scope

ADR-0018 does not:

- replace the Darwin archives or choose new Mac authority;
- ratify any Zaki or CONFIND source identity;
- redesign Zaki/CONFIND APIs or move visualization to Python;
- change the C++ standard or repair `using enum`;
- install GSL/Python or build any dependency;
- change an equation, tolerance, baseline, test inventory, or scientific input;
- authorize direct Mac-to-cluster SSH, source transfer, Slurm, CQ0-CQ7, or a merge.

## 10. Ratification record and next gate

This proposal was drafted from the authenticated canonical source, two read-only local external
repositories, the authenticated Darwin artifacts, and the byte-identically imported cluster
preflight. A later local-only forensic task recovered a matching non-Git CONFIND snapshot and
recorded it in `docs/validation/PHASE6_CONFIND_SOURCE_RECOVERY.md`; the snapshot is not accepted
source authority unless the owner preserves and authorizes it.

On 2026-09-26 the human owner explicitly accepted the four cache-variable names and all eleven
decision clauses in section 3. The owner also accepted, strictly as **SOURCE-EQUIVALENCE
CANDIDATE ONLY**:

- Zaki Git SHA `b9ddebaded24962468954846f47238aec2726fd4`; and
- CONFIND immutable source-manifest SHA-256
  `ed76163c22e0a1f8ba3f71f62f5a56527528850650bbf7127da9c0c908d14083`.

The complete acceptance record is `docs/validation/PHASE6_ADR0018_ACCEPTANCE.md`. The owner
separately authorized a bounded Mac-only source-build equivalence experiment after canonical
integration of this acceptance. Neither candidate is thereby accepted as a governed Linux or
cluster dependency source authority. Compilation alone is insufficient; the same-Mac gate in
section 5 must pass and a later explicit owner decision must grant source authority. Cluster
bootstrap, transfer, CQ0-CQ7, Slurm, clean BA12R, vendored-archive replacement, C++ standard
change, and scientific-model change remain unauthorized.
