# Linux transition qualification predeclaration

Status: candidate, 2026-10-10; no campaign or canonical integration authority.
The owner explicitly authorized cluster bootstrap, native builds, portability
repairs, and bounded numerical qualification in this task. Change classes:
dependency/build, engineering, documentation and generated evidence. No scientific
semantics or numerical algorithm changes are planned. Authority: `GOVERNANCE.md:36`,
ADR-0018 and its ADR-0019 successor; exact dependency identities are those in
`cmake/dependency-lock.json:1`. The Mac lock and all historical baselines remain
immutable. The Linux lock authenticates candidate artifacts, not scientific results.

## Source and scope

Base: cluster preparation d8e7208b71e01168394e41f8f8f32b027b5d0a63, clean imported
sparse checkout; one worktree. Readiness 5433fc2 is an ancestor and has identical
production/CMake/test bytes. Canonical master 9a5c1ec is preserved. Exact importer
passed all three repositories' refs, commit sets, source files and fsck. Transfer
verifier passed 55 payload files, 23 scientific identities and three bundles.
Raw evidence: `/mnt/sdd/zaki/CompactStar/bootstrap-20261010/import.log:1` and
`/mnt/sdd/zaki/CompactStar/bootstrap-20261010/transfer-verification.log:1`.

## Criteria fixed before numerical runs

1. Dependency functional tests and scientific regressions retain every existing
   threshold and baseline. Report original exact-hash failures as failures, even
   if a separate floating-point diagnostic explains them. Never rewrite or relax
   Mac oracles. Run Debug and Release separately and serially within each suite.
2. Algebraic/analytic static, EOS and Phase-5B/5C/5D checks use their already
   declared residual/error bounds. Passing establishes only the named contracts.
3. Cross-platform bounded source/control assessment uses the componentwise scales
   and uncertainty rules in `docs/validation/PHASE6_CAMPAIGN_ACCURACY_PREDECLARATION.md:14`:
   BASELINE (rtol=1e-7, atol=(1e-12,1e-18,1e-18)); REFINED
   (1e-9,(1e-14,1e-20,1e-20)); ULTRA (1e-11,(1e-16,1e-22,1e-22)).
   For a Mac/Linux discrepancy, both runs' justified uncertainties must fit in the
   combined component budgets; missing uncertainty means UNQUALIFIED, not PASS.
   These existing component budgets reflect the requested scientific accuracy;
   machine epsilon or bitwise equality is not substituted for trajectory accuracy.
4. Retain R18 <=1e-10, R20 <=2e-4, reconstruction <=5e-6 and all event and
   reference allocations from that predeclaration. Retain historical BA12/BA12R
   FAIL. Local windows cannot establish global trajectory accuracy.
5. Discrete input identities, certificate decisions, schemas, enum mappings and
   schedule contracts remain exact. Across platforms, report step-count and event
   differences; do not require identical adaptive step histories as an unsupported
   scientific criterion. Within-platform deterministic schedule invariance remains
   governed by its existing tests.
6. Historical dependency floating-point tables have Mac-only authority. Exercise
   available Linux producers and compare diagnostic absolute/relative/ULP drift;
   do not label nonzero drift equivalent without a relevant accuracy authority.
7. First run dependency/data-free contracts including all 13 readiness tests,
   then applicable authenticated EOS/response and bounded Phase-6 tests. Every
   failure is retained and investigated. Test timeouts remain explicit. No 500 kyr
   run, no full campaign, no automatic baseline or tolerance regeneration.

## Toolchain candidate

GCC 14.1.0, CMake 3.29.4 and GNU Make 4.4.1 via shared Spack; native x86_64.
Compute nodes lack development headers and Tcl module support. Load modules on
login before Slurm export. Stage RPM-content-verified login glibc 2.34-231.el9_7.10,
headers and zlib 1.2.11-40.el9 in an isolated sysroot and explicitly select its
loader/runtime. Ownership-only RPM differences arise in the sandbox; no content
mismatches were accepted. Hash the staged files and effective runtime closure.
This is a distinct runtime candidate from the node's system glibc 2.34-100.el9_4.2.
Use supplied GSL 2.7.1 source, static, -O2 -g -ffp-contract=off; no fast-math,
LTO or native tuning. Preserve Zaki/CONFIND per-source HistoricalMath protections.
Python wheel versions and hashes will be recorded before their scientific use.
A successful toolchain smoke is not CompactStar build or scientific qualification.
