# Source-only Stage A launch receipt — job 19195

Submitted once from preparation receipt commit fac8080 after expected-commit,
branch, clean-worktree, snapshot authentication and exclusive-output checks
(`evidence/preflight.txt:1`, `evidence/submission-commit.txt:1`). The owner explicitly
authorized this submission and five-hour ceiling
(`../parallel-preparation/SOURCE_RUN_AUTHORIZATION.md:1`). Classification:
engineering execution/documentation/generated evidence under GOVERNANCE §2;
ADR-0017 §§4–6 and Amendment 1 continue to govern checkpoint reconstruction.
No further implementation changes were made for this launch.

Slurm job 19195 started 2026-11-10 00:07:52 EST on phyn004. Allocation is one
node, one coordinator task, 24 physical worker cores, 96 GiB, five hours;
Slurm reserves 48 logical CPUs because SMT siblings belong to the allocated
cores. Workers select one execution CPU per physical core. The enforced end
ceiling is 05:07:52 EST (`evidence/allocation.txt:1`,
`evidence/coordinator-affinity.txt:1`). Submission used the exact command in
`../parallel-preparation/REPRODUCE.md:1`, including `--time=05:00:00` overriding
the unchanged authenticated wrapper's default.

Output remains under
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-a-parallel-20261110-r1/source-stage`.
The source snapshot identity is
`cbe2cd7f4042b7586149b28428d6ea86bfe245997a520c5ca82829734b5b8a9f`;
the matched, unexecuted control identity is
`880dee60847576fe829dc0c8fd854862c664768a101453d7430f3c1f4bf126e6`
(`../parallel-preparation/evidence/snapshot-diff.json:1`). Both snapshots passed
fresh checking immediately before submission. Original physical inputs,
tolerances, sequential integrator, matrix, accepted-endpoint restart selector,
accuracy gates and guards are preserved. Only checkpoint scheduling uses 24
workers, with byte-identical retained regression results
(`../parallel-preparation/evidence/regression.json:1`).

Status at this receipt: RUNNING, NOT SCIENTIFICALLY QUALIFIED. The 3–4 hour
estimate concerns checkpoint collection; actual-output high-precision reference,
main-prefix/accumulated and thermal roundoff, reference convergence, R18/R20,
restart consistency and conditional uniform thermal enclosure/refinement
obligations remain. The 13-point study does not discharge actual-output gates.
No matched control, automatic retry, Stage B, stellar solve, BA13 reopening,
master merge or certificate promotion has been launched.

Only immutable launch evidence is hashed here. Live logs/output must be
inventoried and authenticated after the job stops; this receipt does not claim
an immutable digest of changing files or a completed numerical result. The
original failures are preserved. Preparation code and evidence are committed
locally; public GitHub publication is pending explicit confirmation after
automatic review rejected public release of cluster evidence
(`evidence/publication-status.json:1`). This does not affect execution.
