# Native runtime repair and Stage A retry

**Runtime repair passed; Stage A ran and failed checkpoint diagnostic qualification.**
Stage B was not submitted. The one authorized retry was Slurm job **19174** on
phyn004. Five one-year source references, the first BASELINE source main, and its
endpoint-schedule repeat completed. The first nonzero passive checkpoint's
diagnostic uncertainty failed the existing acceptance gate. Control and restart
runs were never reached. No parameter, tolerance, algorithm, physical input, or
certificate was changed to obtain acceptance (`evidence/assessment.json:1`,
`evidence/stage.log:42`).

## Loader mismatch and minimal repair

Original job 19171 pinned `/usr/lib64/ld-linux-x86-64.so.2` through host `ldd`
output. Login ekucn1 runs Rocky 9.7 and that loader has SHA-256
`ad96200950148efb34358c6e980e06b15a807644ce0a003e1d8a71f9f0820fcd`.
Compute phyn004 runs Rocky 9.4; its host loader has SHA-256
`17a1a80f01feb4a94c5e357f20b62db00c9c7746e32ac7edfad63b1ebc9a5a09`.
This was the sole old-manifest mismatch. Original snapshot files and failure
evidence remain unchanged (`../native-qualification/REPORT.md:1`).

The actual ELF PT_INTERP for both frozen native executables is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/sysroot/usr/lib64/ld-linux-x86-64.so.2`.
Its hash is `ad96200950148efb34358c6e980e06b15a807644ce0a003e1d8a71f9f0820fcd`
on both nodes. That interpreter's own `--list` resolves the shared toolchain's
libz, libm, libc and GCC's libstdc++/libgcc; the unused host loader is absent.
The replacement publisher records these actual runtime files, lexical/resolved
paths, and PT_INTERP. Loading rechecks resolution and authenticates every file;
missing/unrecognized dependencies, changed resolution, unbound interpreters,
and undeclared preloads/auditors refuse. No genuine dependency is suppressed
(`evidence/runtime-audit-summary.json:1`, `evidence/snapshot-diff.json:1`,
`tools/linux/native_pair/runtime.py:1` from repository root).

This is an engineering/dependency-publication repair, with documentation and
generated evidence. Governing authority is the owner's explicit retry request,
GOVERNANCE sections 2–3, accepted ADR-0015/0016/0017 and INV-11f/INV-15; candidate
owner transport remains as described in ADR-0024. No platform expansion is
claimed. Both static preflight and thermal retry remain pinned to phyn004.
Login/compute variation was checked; other compute nodes require their own
authentication and are not qualified by this result (`PREDECLARATION.md:1`).

## Replacement identities and unchanged science

| Snapshot | SHA-256 |
|---|---|
| P2 source | `c6a29044e8330f079a1593be166238e11eafb66745186c5322c04c272358ec18` |
| Exact-zero control | `0588ddf4518b3dec0f1cf103dcb8a3b356cf4631f2394df41b6a2d6baadf4ceb` |

The exact machine-readable diff against each original is
`evidence/snapshot-diff.json:1`. Only `assembly_commit`, `files`, and the added
`runtime_bindings` metadata differ. File inventory: 508 → 509; remove the unused
host loader, update `snapshot.py`, add `runtime.py` and `test_runtime.py`.
Every physical field and all other file hashes remain identical. This includes
the full stellar archive, source definitions, coefficient tables, initial state,
integration tolerances, boundary law, conditional sidecar, and scientific sources.
Both replacement roles have identical common content; only their source object
differs (`snapshots/identities.json:1`, `evidence/snapshot-diff.json:1`).

Unchanged full-owner archive SHA-256:
`26dd2b060ad8ed10d26ec46b67a9b957e3517fe6718708ba91d09c5f204158cf`.
Unchanged stage ELF SHA-256:
`1e90cf96d2a39b96a56b7b8f8b7432f6c8c3fec1a92fb7217e55bd907693bfcd`.
Unchanged verifier ELF SHA-256:
`712e8c5b3ebd56ac93871324a1dd97a02b15cdb3192fd5497b5169cbebaa6923`
(`snapshots/source.json:1`). No native binary was rebuilt.

Repair commit: `66419b813c15ebade0e2f1bc431666a3b4f07364`.
Snapshot publication commit: `ab01d7a`.
The execution wrapper/predeclaration was committed in `af5ae3f`.
Starting branch was clean at the requested `46e0b7c`; ancestry and worktree
preflight are recorded in `evidence/preflight.json:1`.

## Prerequisites passed before evolution

Job **19173** checked 509/509 dependency hashes on phyn004, as on ekucn1, with
matching actual ELF linkage. All four focused runtime-discovery/refusal tests
and snapshot corruption/missing/mismatched/invalid-state tests passed. Static
assembly exercised original ownership/currentness and certificate refusals,
source/control isolation, R18 and initial-state round trips. Twelve RHS states
matched the retained baseline byte-for-byte, SHA-256
`ed01bc8da5e509a0807e2700c52689920ad3d6f11d06c9e37d4df504697972e4`.
All six retained contracts passed: direct energy, thermal ledger/R20, passive
checkpoint, breakpoint, component tolerances, and currentness. These static
successes did not establish nonlinear checkpoint accuracy
(`evidence/preflight-19173.log:1`, `evidence/static.log:1`,
`evidence/static-rhs.tsv:1`, `evidence/contracts.log:1`).

## Stage A failure and completed evidence

The exact stop was `STOP native checkpoint diagnostics unresolved`.
It occurred for the **source, BASELINE, uncapped** main, at checkpoint index 1
on the 240-interval scientific grid: **131,490 s = 1.521875 days**.
Only the initial row was saved in `uncertainties-240.tsv`. State reconstruction
passed before `EvaluateDiagnostics` changed the checkpoint status to unresolved.
At least one diagnostic failed **U <= 0.20 F**. The unchanged collector did not
serialize the failing diagnostic name, its U/F terms, or the rejected packet;
those specific values cannot be recovered from the saved output alone
(`tests/native/stage.cpp:82`,
`CompactStar/Physics/BNV/src/PassiveCheckpointOutput.cpp:623` from repository root;
`evidence/assessment.json:1`). No diagnostic replay or further trajectory was run.

Completed raw evidence consists of five source references (each 1,921 observation
rows spanning exactly one Julian year), one 12-step source main with zero rejected
steps, and its repeat with only endpoint observations. Accepted state histories
are byte-identical between main and repeat. The main has no Cstar-knot event.
References reached the exact endpoint with finite values and temperatures in the
declared range. The collector passed each reference's rejection ceiling before
writing its observations. These facts do not resolve the failed diagnostic gate
(`evidence/assessment.json:1`, `tests/native/stage.cpp:34` from repository root).

The original R20 grid file contains its header only; no R20 result was computed.
The buffered diagnostic capture was never published, so no dynamic checkpoint R18
packet is available. Restart tails, the other 11 source policies/tiers, and all
control trajectories were not reached. Independent accumulated-roundoff and
reference-interpolation bounds remain unestablished. Small raw reference spreads
cannot be substituted for those bounds. No convergence PASS, energy-conservation
PASS, actual-event qualification, or source-minus-control result is claimed
(`evidence/assessment.json:1`, `../native-regeneration/RUN_PLAN.md:1`).

## Preliminary source-only outputs — diagnostic, not qualified

The following endpoint is from completed reference 4 (independent-method RKF45,
rtol 1e-15, atol `(1e-20,1e-26,1e-26)`, cap 15,408.984375 s). It is a descriptive
choice, not a certified truth solution. Full time series and figures remain
external; the initial/final rows are in `evidence/source-endpoints.tsv:1`.

| Quantity at one year | Diagnostic value |
|---|---:|
| Redshifted temperature | 100,007,404.51388617 K |
| Increase from source initial temperature | 7,404.51388617 K |
| Electron imbalance | -1.0894196614137575e-10 MeV |
| Muon imbalance | -2.7634376967724932e-10 MeV |
| Neutrino luminosity | 1.9125473233321554e32 erg/s |
| Photon luminosity | 2.7754493815996762e33 erg/s |
| Actual BNV direct deposition | 3.3583239602973686e34 erg/s |
| Reaction heating | 5.054307851959247e15 erg/s |
| Excess neutrino power | 7.581461777938871e15 erg/s |

All powers are at infinity. This remains the controlled mathematical neutron-sink
benchmark, not a microscopic decay-rate prediction. The temperature increase is
relative to the source initial state, **not** a source-minus-control difference.
Across all five reference grids, maximum raw spread is 3.2526065174565133e-19 in
x, 2.714186384939622e-25 MeV in eta_e, and 5.945360652724886e-25 MeV in eta_mu.
Reported T-infinity values agree at their stored binary64 precision; this is not
a zero uncertainty claim (`evidence/assessment.json:1`).

Figures are explicitly marked DIAGNOSTIC ONLY:

- [Source thermal, chemical and power curves](/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r2/analysis/source-diagnostic-curves.png)
- [Raw reference spread, not an error bound](/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r2/analysis/reference-spread.png)

Matching PDFs are beside the PNGs. The plotting script reads retained data only;
it calls no solver, RHS or checkpoint reconstruction. Its isolated plotting
environment does not modify the frozen numerical toolchain
(`tools/linux/native_qualification/summarize_retry.py:1` from repository root).

## Resources, preservation, and remaining blocker

| Work | Slurm job | Wall measurement | Peak child RSS | Outcome |
|---|---:|---:|---:|---|
| Static replacement verification | 19173 | 112.751 s for verification | 535,120 KiB | PASS; six contracts also passed |
| Single Stage A retry, including assembly | 19174 | 167.087 s | 434,912 KiB | FAIL checkpoint diagnostics |

Stage A requested one task/CPU, 4 GiB and 30 minutes. Slurm's whole-core allocation
reserved two SMT threads, while the `srun` execution mask was exactly **CPU [0]**;
the wrapper refuses if affinity is not exactly one logical CPU. Single-thread
libraries remained enforced. Slurm reports FAILED / exit 1 after 2m47s, with no
requeue/restart. Its resources stayed within the declared execution/memory/time
limits (`evidence/cpu-affinity.txt:1`, `evidence/allocation.txt:1`,
`evidence/stage-A-final.txt:1`, `evidence/stage.log:54`).

Runtime publication is repaired. The remaining immediate blocker is diagnostic
self-qualification at the first nonzero checkpoint; the absent failing-channel
receipt prevents a more specific diagnosis without separately bounded diagnostic
work. No repair or retuning of that numerical gate was attempted. Stage B and the
full restricted comparison are **not qualified to proceed**. The one retry
authorization has been consumed (`evidence/assessment.json:1`).

Source changes are limited to publication/runtime helpers and tests; added
qualification scripts provide preflight, CPU affinity, and read-only summary/
plots. No stellar regeneration, BA13 reopening, historical certificate promotion,
master merge, 14,648-year run or 500 kyr campaign occurred. BA12/BA12R stay FAIL;
conditional BA13 retains its original assumptions and unchanged bytes. This
delivery remains a scientific candidate, not human ratification.

All raw results are retained under `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/`:
`native-p2-runtime-retry-20261109-r1` and `native-p2-qualification-1yr-r2`.
`evidence/external-SHA256SUMS` authenticates both trees, including full reference
histories and figures. Git contains source, replacement snapshots, exact manifest
differences, compact receipts, this report and reproduction commands. Raw Slurm
receipts preserve their original whitespace. See `REPRODUCE.md:1` for commands;
none grants another thermal retry.
