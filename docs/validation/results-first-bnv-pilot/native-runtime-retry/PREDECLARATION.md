# Runtime-only repair and single conditional retry

Owner authorization at starting commit `46e0b7cd1389048a5b6c06843358105868a39903`
permits correcting runtime publication, replacing snapshot identities, and one
conditional Stage A retry. Stage B remains conditional on all applicable Stage A
gates. Original snapshots, job 19171, and its failure report are preserved.

The mismatch was an unused host loader included by parsing host `ldd` output.
Publication now reads ELF PT_INTERP and queries that interpreter with `--list`.
It records lexical/resolved paths and hashes every actual loader/library, refuses
missing/unrecognized linkage and undeclared preloads/auditors, and rechecks the
resolution at each snapshot load. The supported execution target remains
phyall/phyn004. Login/compute OS variation is tested explicitly; no other compute
node is qualified by inference (`../native-qualification/REPORT.md:1`,
`tools/linux/native_pair/runtime.py:1` from repository root).

Change classes: dependency/build and engineering, with documentation/generated
evidence. Governing authority: owner authorization; GOVERNANCE sections 2–3;
ADR-0015 energy definitions, ADR-0016 owner/currentness, ADR-0017 checkpoint
contracts; INV-11f and INV-15. No numerical source, executable, stellar input,
coefficient, source definition, initial condition, or integration setting changes.
Candidate status and all historical BA12/BA12R failures and conditional BA13
assumptions remain unchanged. No new equation or scientific source derivation.

Before retry: authenticate both replacement snapshots on phyn004; run focused
publication/refusal tests; reconstruct all owners from the same archive; compare
all twelve static RHS states byte-for-byte to retained output; exercise existing
ownership, missing/mismatched certificate and source/control guards; pass the six
retained direct-energy, R20, checkpoint, breakpoint, tolerance, and currentness
contracts. Any failure stops before thermal evolution. The preflight is job 19173,
one task/CPU, 4 GiB, 15 minutes, with no thermal integration calls.

Thermal execution reuses the exact frozen `native_p2_stage` ELF and original
`tools/linux/native_pair/stage.sbatch`. A small Slurm wrapper binds its one task
to one logical CPU and verifies affinity before starting. Slurm may reserve both
SMT siblings of a physical core; the explicit execution affinity distinguishes
that reservation from actual CPU access. Memory limit remains 4 GiB and aggregate
stage wall limit 30 minutes. OMP/BLAS/MKL remain single-threaded.

Stage A end = 31,557,600 s; Stage B end = 64,629,964,800 s, independently
initialized. The full unchanged matrix, state/power norms, R18/R20 requirements,
restart checks, accumulated-roundoff requirement and first-failure/timeout stop
rules are `../native-regeneration/RUN_PLAN.md:1` and
`../native-assembly/NONLINEAR_PLAN.md:1`. Collection completion is not accuracy
acceptance. Unknown error terms remain UNQUALIFIED; no automatic extension or
retuning. Source runs before control in the frozen collector, so an early source
failure can prevent any control output. No full campaign is permitted.

External preflight data:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-runtime-retry-20261109-r1`.
Single Stage A retry output:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/native-p2-qualification-1yr-r2`.
Snapshots and compact evidence are Git-tracked; substantial data remain external
with SHA-256 manifests. No automatic Stage B submission is installed.
