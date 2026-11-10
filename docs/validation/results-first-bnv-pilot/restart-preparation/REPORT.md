# Stage A restart amendment — preparation PASS, trajectory not run

**PREPARATION PASS.** The exact approved four-file patch from `3557da9` is
implemented at `f8a1837bd582eb13d0a8d4783eb9501ea79a25c2`. The complete amended
stage harness builds; all eight selected contracts, the retained-history selector
and exact binary64 round-trip pass. Twelve static RHS states reproduce the prior
reference byte-for-byte. Both replacement runtime snapshots authenticate.
**No new Stage A trajectory was submitted or executed.** The exact source-first
command is in [REPRODUCE.md](REPRODUCE.md:6), for separate owner authorization.

## Implemented scope

For each Stage A role/tier's quarter-cap main, the harness selects the latest
accepted noninitial endpoint at or before half-time. It requires a complete,
contiguous, progressing history and a suitable endpoint, otherwise throws an
explicit refusal. The selected tail is nonzero and at least half the horizon.
No initial-state or passive-output fallback exists
(`tests/native/restart_schedule.hpp:9`, `tests/native/stage.cpp:108`).

The test horizon remains one year. The recorded BASELINE history selects ordinal
11 at 10,025,181.919443566 s, leaving 21,532,418.080556434 s. This selection is a
regression on retained data, not a new restart trajectory or a universal accuracy
claim ([selector-roundtrip.log](evidence/selector-roundtrip.log:1)). The previous
inactive-cap proof still explains the old quarter-cap coverage failure; no new
quarter-cap solve occurred (`../stage-a-recovery/REPORT.md:44`).

The exact approved patch adds the selector/header test, includes it in the stage
harness, and registers its CTest target. Source hashes match all four previously
tested proposal hashes. Stage B retains its original rule and remains unauthorized.
The two fixed run plans have an explicit accepted Stage A amendment notice while
retaining their original text as historical provenance
(`../native-assembly/NONLINEAR_PLAN.md:3`, `../native-regeneration/RUN_PLAN.md:3`).

The production `CompactStar` tree is unchanged relative to `3557da9`. Integration
algorithms, step controls, caps, initial step, all tolerances including O2
(5e-27,5e-27) MeV, physical inputs, and every numerical acceptance threshold are
unchanged. All 86 production-library archive member payloads are byte-identical
to the previous build; the whole archive hash differs only because member
timestamp headers changed ([library-regression.json](evidence/library-regression.json:1)).

Authority is the explicit owner approval of the `3557da9` proposal, preserved in
[OWNER_APPROVAL.md](OWNER_APPROVAL.md:1). ADR-0017 §§4–6 and Amendment 1 continue
to govern passive observations, accepted endpoints and actual Cstar events.
Classification under GOVERNANCE §2: numerical-method/test-protocol, engineering,
documentation and generated evidence. No scientific-semantic or solver-method
change occurs. The protocol is owner-approved; trajectory accuracy remains
unqualified, and scientific candidate status under GOVERNANCE §5 is retained.

## Preparation results

Slurm job **19190** completed on phyn004 with exit 0 in **4 min 53 s**, one logical
execution CPU (affinity [0]), a 4 GiB allocation and a 30-minute ceiling. Slurm
reserved two hardware threads; execution was pinned to one. This is a build/static
preparation job, not use of another Stage A trajectory authorization
([completion](evidence/completion.txt:1), [affinity](evidence/build-affinity.txt:1)).

| Check | Result |
| --- | --- |
| Full `native_p2_stage` build, retained toolchain/Release/`-ffp-contract=off` | PASS |
| Eight CTest contracts, including O1/O2 configuration and restart selector | 8/8 PASS |
| Selector synthetic coverage/refusals | 12 cases PASS |
| Retained accepted-history endpoint selection | ordinal 11 PASS |
| Exact binary64 time and three-component state round-trip | PASS |
| Existing static ownership/currentness/source-control checks | PASS |
| Twelve retained static RHS states | byte-identical PASS |
| New stage ELF/runtime linkage | PASS |
| Stage invoked without arguments | expected argument refusal before input loading/integration |
| Matched source/control snapshot authentication | PASS |
| New thermal/stellar trajectories | **0** |

The eight CTest contracts are passive checkpoint output, direct energy, thermal
ledger, component tolerances, breakpoint handling, currentness bytes, O2
configuration and restart selection. Their logs and the separate real-history
round-trip are retained ([contracts](evidence/contracts.log:1),
[static verifier](evidence/verify.log:1),
[smoke refusal](evidence/stage-linkage-smoke.log:1)). Existing negative tests
intentionally refuse invalid ownership/archive inputs; they are passing contract
tests, not suppressed scientific failures. Compiler warnings remain in the raw
build log; no new warning cleanup was included.

## Snapshot publication was necessary and bounded

Both roles bind the same stage source/ELF. Changing the harness therefore requires
new source **and** control snapshot identities even though the production solver
and physical inputs do not change. The publisher permits only the approved stage
source difference among previously bound inputs, verifies the exact approved patch,
binds the new selector/header/test and preparation provenance, and replaces only
the stage executable path. The unchanged static verifier ELF is retained.
Physical/common numerical fields and runtime linkage are compared exactly;
only common restart-test metadata is added
([snapshot-diff.json](evidence/snapshot-diff.json:1),
`tools/linux/native_qualification/restart_preparation_snapshots.py:16`).

| Snapshot | SHA-256 |
| --- | --- |
| New source | `4fa012bda09d5022df7dce6ad7278978b7f638c3d15786d84750a96da256c4d5` |
| New control | `584251e18fb131cb8f141e7972de6fcfaecf52a53f818d5f5f6a3f3c8d19fb27` |
| Preserved previous source | `7ed42a3ef02d989caa06246940de2cfe75fe7bd554b2650f6314f778fb77c175` |
| Preserved previous control | `9226ea1f2ed0542b53b25012fc4cc9d42221efa4093a20174db3b297fd5daf5a` |
| Unchanged full-owner archive | `26dd2b060ad8ed10d26ec46b67a9b957e3517fe6718708ba91d09c5f204158cf` |

The new stage binary is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/restart-preparation-20261109-r1/build/tests/native_p2_stage`.
Both new snapshots name implementation commit `f8a1837bd582eb13d0a8d4783eb9501ea79a25c2`;
the later receipt commit records the resulting artifacts without changing those
bound source bytes. Old snapshots/ELFs/failure receipts remain historical and are
not rewritten to claim authentication of the amended source tree.

## Provenance and remaining scope

Entry was clean at `3557da9861274033f642ca08f690b4fb59105797`, branch
`candidate/results-first-bnv-pilot-20261108`, worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All 15 local branch tips are ancestors; nine worktrees exist. Proposal and recovery
manifests (15+12 and 58+1894 entries), and the previous runtime pair, authenticated
before the final patch application ([preflight.json](evidence/preflight.json:1)).
An initial preflight cleanliness assertion caught two Python-generated bytecode
files; the exact patch was reversed, those generated files removed, and clean
authentication repeated with `-B` before reapplication. No numerical execution
occurred during that setup correction; the preflight receipt preserves it.

Large build/static outputs remain under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/restart-preparation-20261109-r1`.
Git contains the implementation, approval/plan annotations, maintained preparation
tools, new snapshots, compact receipts and manifests. Original D1 FAIL and later
conditional PASS, Stage A Echem/restart failures, original plans' historical text
and the `3557da9` proposal remain preserved. No master merge, BA13 reopening,
stellar solve, Stage B or certificate promotion occurred.

This receipt establishes **preparation only**. A future source run must still pass
actual-output chemical/Echem reference checks, reference convergence, thermal
enclosure/refinement assumptions, prefix and accumulated-error accounting,
thermal roundoff, R18/R20 and restart accuracy. Missing bounds are UNQUALIFIED.
The earlier 13-point conditional PASS is not universal checkpoint qualification.
Control remains contingent on complete source qualification. The provided source
command collects evidence; it does not automatically certify or advance a run.
