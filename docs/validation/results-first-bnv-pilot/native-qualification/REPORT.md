# Native P2 thermal qualification — startup refusal

Stage A did not start numerical evolution. Slurm job **19171** failed at snapshot
dependency authentication on phyn004, before invoking the native C++ runner.
Stage B was **not submitted** because Stage A did not pass. No thermal retry,
stellar regeneration, tolerance change, or certificate promotion was performed.
This supersedes the earlier delivery's claim of no remaining runtime blocker;
its static assembly results remain historical evidence, not a thermal result
(`../native-regeneration/REPORT.md:1`, `evidence/stage.log:1`).

## Execution and identities

Owner authorization covered the fixed 1-year matrix and, conditionally on its
acceptance, an independent 2048-year matrix at starting commit
`3436063dec533d0e7be7b9f01b74e95c64fd0da2`. Repository:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`,
branch `candidate/results-first-bnv-pilot-20261108`. Initial worktree was clean;
all 15 local branch tips were ancestors of HEAD; nine worktrees were present
(`evidence/preflight.json:1`). These exact authorized snapshot identities remain:

| Role | SHA-256 |
|---|---|
| Mathematical P2 source | `93f1dab5269eb36259598839426f33ce28132430ea92602c4892317d832e48d8` |
| Exact-zero control | `fe692927032cc1387e9bbc404d12615d4c170f07c4abd79fdf036ece6a8a94cc` |

Both roles bind the same 508-file inventory. Login-node authentication passed
before submission and again afterward; snapshot tampering, missing/mismatched
dependency, and invalid initial-state refusal tests passed. The submitted wrapper
refused on the compute node. No `run/` directory exists: authentication precedes
its creation and the C++ invocation (`tools/linux/native_pair/snapshot.py:93` from repository root,
`evidence/assessment.json:1`, `evidence/post-run-snapshot-check.log:1`,
`evidence/snapshot-refusal-tests.log:1`).

## Exact blocker

The exception was:

```text
ValueError: missing/mismatched authenticated dependency: /usr/lib64/ld-linux-x86-64.so.2
```

| Location | Host OS | Host loader SHA-256 |
|---|---|---|
| ekucn1 login / expected in snapshots | Rocky Linux 9.7 | `ad96200950148efb34358c6e980e06b15a807644ce0a003e1d8a71f9f0820fcd` |
| phyn004 compute / actual | Rocky Linux 9.4 | `17a1a80f01feb4a94c5e357f20b62db00c9c7746e32ac7edfad63b1ebc9a5a09` |

A separate **read-only** Slurm audit, job **19172**, checked all dependencies:
508/508 match on login, 507/508 on compute, with precisely the host loader above
differing. No solver, reconstruction, RHS, or trajectory ran in that audit.
Direct SSH was unavailable; the audit used Slurm instead
(`evidence/assessment.json:1`, `evidence/runtime-audit-final.txt:1`).

The ELF interpreter is the shared toolchain's
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/sysroot/usr/lib64/ld-linux-x86-64.so.2`.
Its SHA-256 is `ad96200950148efb34358c6e980e06b15a807644ce0a003e1d8a71f9f0820fcd`
on both nodes. Its `--list` output resolves the same library paths on both nodes;
all those files match their snapshot hashes. The host loader is absent from that
resolved list. Publication currently parses **all** absolute file tokens in host
`ldd` output, including its host-loader alias, and resolves/pins them. This explains
the unnecessary host-specific binding; it does not authorize bypassing it
(`tools/linux/native_pair/snapshot.py:46` from repository root;
`evidence/assessment.json:1`).

Smallest proposed repair: make publication authenticate the ELF's actual
interpreter and its resolved dependencies, with explicit missing-library refusal;
test publication/authentication on the target compute node; publish new source
and control identities retaining the same stellar archive, coefficients, binary,
physical settings, and numerical settings. Do not merely ignore a mismatched
dependency. No implementation of this proposal or snapshot replacement was made.
The exact-hash authorization and stop rule require an explicit decision before a
new thermal attempt with replacement snapshots. No new stellar calculation is
needed to address this observed startup blocker. Successful startup would still
leave every nonlinear gate to be established (`../native-regeneration/RUN_PLAN.md:1`).

## Scientific and resource disposition

| Required output or gate | Result |
|---|---|
| Source/control trajectories, T-infinity, chemical imbalances | Not produced; zero samples |
| Neutrino/photon luminosities, BNV deposition, reaction heating | No trajectory diagnostics |
| R18 / R20, nonlinear convergence, accumulated roundoff | Not evaluated |
| Checkpoint/restart and actual breakpoint agreement | Not evaluated |
| Source-minus-control temperature / signal uncertainty | Not available |
| Preliminary figures | None: no completed or partial trajectory exists |
| Stage B / readiness for full restricted comparison | Not run / **not qualified** |

These are absent measurements, not zero-valued results or scientific accuracy
failures. The failure is the required startup authentication gate. Static RHS
values from the earlier package must not be plotted as evolution
(`evidence/assessment.json:1`).

Stage A requested one task, one CPU per task, 4 GiB, 30 minutes, with single-thread
libraries and `--hint=nomultithread`. Slurm recorded FAILED / exit 1 after one
second. Wrapper measurement: **0.665205 s**, **19,220 KiB** maximum child RSS,
0.238702 user CPU s and 0.119861 system CPU s. No thermal resource estimate can be
inferred from this startup-only measurement (`evidence/allocation.txt:1`,
`evidence/slurm-final.txt:1`, `evidence/stage.log:1`).

Slurm reserved/billed **two logical CPUs**, despite CPUs/Task=1 and one thread per
core requested. The node reports two threads/core; this is consistent with a
whole-core reservation. It must not be reported as a literal one-logical-CPU
allocation. Execution was a single authentication process, with no thermal solver
started. The read-only audit requested one CPU/4 GiB/2 minutes and completed within
Slurm's one-second timestamp resolution (`evidence/slurm-node.txt:1`,
`evidence/runtime-audit-allocation.txt:1`, `evidence/runtime-audit-final.txt:1`).

## Scope, preservation, and review status

Changes: this report, reproduction instructions, compact evidence/manifests, and
`tools/linux/native_qualification/{runtime_audit.py,runtime_audit.sbatch,collect_refusal.py}`.
Classes: documentation, generated evidence, and engineering for read-only evidence
collection. No numerical or physical source changed. Governing authorities:
GOVERNANCE sections 2–3, accepted ADR-0015 (energy meaning), ADR-0016 (owner
currentness), ADR-0017 (checkpoint contract), INV-11f/INV-15, the candidate owner
transport ADR-0024, and the owner-authorized bounded plan. Outputs remain an
unverified scientific candidate; no new accuracy baseline or qualification exists.

BA12/BA12R failures remain FAIL. Conditional BA13 stays at 420 conditional rows,
21 exact identities, 21 metric-consumer-only rows, maximum utilization
0.782464607730239, with all prior assumptions intact. Its bytes are part of the
authenticated snapshot inventory, not recalculated here
(`../native-regeneration/REPORT.md:75`, `evidence/assessment.json:1`).

External evidence remains under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r1`.
`evidence/external-SHA256SUMS` authenticates every retained external file, including
the full per-node inventories; compact receipts are Git-tracked. No master merge,
full 14,648.4375-year benchmark, or 500 kyr campaign occurred. Reproduction commands
are in `REPRODUCE.md:1`; they distinguish evidence inspection from prohibited
automatic resubmission.
