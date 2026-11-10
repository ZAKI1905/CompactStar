# 24-worker preparation receipt — PASS, scientific candidate

The independent checkpoint scheduling implementation passes preparation. The
original source trajectory is not yet scientifically qualified. EKU job 19194
completed successfully in 13m31s, using 24 physical execution cores, 96 GiB and a
30-minute ceiling (`evidence/completion.txt:1`). No new main integration,
reference trajectory or stellar solve occurred (`evidence/regression.json:24`).

Authority and class: owner's explicit 24-worker approval, GOVERNANCE §2
engineering/build/documentation/generated evidence, and ADR-0017 §§4–6 plus
Amendment 1 for immutable accepted histories and isolated reconstruction
(`OWNER_APPROVAL.md:1`). The original entry authentication, all-branch ancestry
and worktree inventory are retained in `evidence/preflight.json:1`. Implementation
HEAD is recorded in `evidence/source-head.txt:1`; preparation started clean
(`evidence/source-status.txt`, empty). Source of truth remains the candidate
branch in the absolute worktree recorded in OWNER_APPROVAL.md.

## Implementation and numerical impact

`tests/native/stage.cpp:31` selects one or 24 workers; `CollectGrid` dispatches
isolated checkpoint reconstructions and aggregates in original ascending order.
The new process pool (`tests/native/process_pool.hpp:1`) pins persistent child
processes to distinct physical cores and stops dispatch on failure. Complete
packets cross a same-build binary64 transport (`tests/native/checkpoint_transport.hpp:1`),
with strict schema/size checks and atomic publication. Internal binary packets
are not portable restart files. Context instance identifiers are deterministic
bookkeeping labels; mathematical inputs and allowances are unchanged.

The sequential main integrator, references, repeats, restart selection and
restart tails are unchanged. No production CompactStar source changed. All 84
production archive object payloads match the previous build byte for byte
(`evidence/library-regression.json:1`). The six new native helper/test files,
stage harness, tests/CMakeLists.txt and three parallel preparation/execution tools
are the maintained implementation; reproduction is in `REPRODUCE.md:1`.

## Validation and performance

All eight retained contracts passed (`evidence/contracts.log:1`). The restart
selector and actual ordinal-11 binary64 round trip passed
(`evidence/selector-roundtrip.log:1`). Pool tests cover 24 distinct pinned cores,
73 tasks, exception/crash/SIGTERM cleanup and bounded dispatch
(`evidence/pool-tests.log:1`). Transport tests cover signed zero, subnormals,
round trip and malformed input refusal (`evidence/transport-tests.log:1`).

On the authenticated retained 240-interval history, both serial and 24-worker
reconstruction reproduced 967 original files/results each, including the exact
R20 row. Serial collection took 396.438780075 s; 24 workers took 65.428782841 s:
6.0590884143-fold speedup (`evidence/regression.json:1`). Static verification
reproduced all 12 RHS rows exactly and runtime linkage is identical
(`evidence/snapshot-diff.json:1`). These are scheduling regressions, not an
independent accuracy bound. Extrapolation to 43,248 checkpoint packets gives
about 3.26 hours, hence a practical 3–4 hour collection estimate; this excludes
subsequent independent high-precision scientific adjudication and is not a
wall-time guarantee (`REPRODUCE.md:1`).

## Authentication and preserved attempts

Replacement source identity:
`cbe2cd7f4042b7586149b28428d6ea86bfe245997a520c5ca82829734b5b8a9f`.
Replacement control identity:
`880dee60847576fe829dc0c8fd854862c664768a101453d7430f3c1f4bf126e6`.
The control snapshot was published but not executed. Old identities, binaries
and snapshots remain intact. The stellar owner archive remains
`26dd2b060ad8ed10d26ec46b67a9b957e3517fe6718708ba91d09c5f204158cf`
(`evidence/snapshot-diff.json:1`). Snapshot checking authenticates executable,
source bindings, physical inputs and runtime dependencies before execution.

Preparation 19192 failed compilation because the regression included a second
main entry point; commit f55879d corrected the test entry guard
(`attempts/19192/build-tail.log:1`). Preparation 19193 failed its prebuild check
before compilation/numerics. Its original wrapper did not retain status output,
so the exact cause remains unknown; the last observed operation was Git status
on the shared filesystem. Commit 7ebbae0 records read-only Git checks without
optional index writes (`attempts/19193/build-stage.log:1`). No scientific gate
was changed to resolve either engineering failure. Both attempts and successful
19194 outputs are authenticated in `evidence/external-SHA256SUMS:1`; raw logs
retain original whitespace. The earlier serial TIMEOUT and original D1 FAIL
receipts are unchanged (`../stage-a-source-retry/REPORT.md:1`).

## Authorized next execution and remaining obligations

The owner's subsequent instruction explicitly authorizes submitting the source
job with a five-hour ceiling (`SOURCE_RUN_AUTHORIZATION.md:1`). The exact
source-only command is in `REPRODUCE.md:1`; it uses a fresh exclusive output path
and overrides the unchanged authenticated wrapper's two-hour default. It does
not submit the control. No Stage B, stellar regeneration, BA13 reopening,
master merge or certificate promotion is included.

Collection completion is not full qualification. Actual-output high-precision
chemical-state/Echem reconstruction, reference convergence, main-prefix
accuracy, accumulated and thermal roundoff, R18/R20 and restart consistency
remain required. The prior 13-interval knot-aligned conditional PASS does not
qualify every actual checkpoint or accumulated trajectory error. Uniform
continuous thermal enclosure and refinement assumptions must remain explicit;
sampled thermal agreement and two constant-offset runs alone do not establish
uniform extrema. No chemical-heating observability claim follows from this
short direct-deposition experiment. Candidate status remains under GOVERNANCE
§5 pending scientific adjudication. No unresolved source-of-truth conflict was
found; the historical 19193 prebuild failure remains unexplained.
