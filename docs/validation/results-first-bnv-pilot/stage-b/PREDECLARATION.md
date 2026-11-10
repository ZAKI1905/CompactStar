# Stage B: 2,048-year paired evolution

Owner authority: 2026-11-10 request, recorded in `../preliminary-pair/OWNER_ACCEPTANCE.md:1`.
Starting source: `c90b2205e553afcd102c7acd5f84b253e0a67f1b`.
Worktree: `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`;
branch `candidate/results-first-bnv-pilot-20261108`. All fifteen local branch tips
are ancestors of the requested source; nine worktrees exist (`evidence/preflight.json`).
Only pre-existing Python bytecode was untracked at entry. No source-of-truth conflict.

The original matrix, equations, initial state, accepted stellar inputs, thresholds,
stops and Stage B restart rule remain `../native-assembly/NONLINEAR_PLAN.md:49–134`.
Both roles start at original t=0 and end at exactly 64,629,964,800 s. Source is the
mathematical full-retention P2 benchmark; control is exact zero. No microscopic
rate, realistic decay prediction or astrophysical constraint is inferred.

`tools/linux/stage_b/driver.cpp:1` uses a generated copy of the retained collector
with only reference/tier scheduling filters. There are ten independent references
and six independent tier tasks, each retaining four caps, endpoint-only schedule
repeat, four nested scientific grids, event neighbors, and quarter-cap restart.
A scientific guard failure stops that affected task; partial evidence is retained.
No retry or extension of the matrix is authorized by this document. New numerical
verification remains a candidate under GOVERNANCE §5 until assessed/reviewed.

Resources supersede the obsolete one-CPU/30-minute proposal. Build separately;
run independent references concurrently on one core each, and six tier tasks
with 24 physical reconstruction workers each. Initial tier memory: 48 GiB each,
based on the completed Stage A allocations (`../preliminary-pair/REPORT.md:180`).
Use Slurm phyall, no node pin, unlimited partition walltime; further independent
high-precision/reconstruction assessments use at least 24 workers where beneficial.
The cluster currently has four idle 96-CPU/125,000-MB nodes. Heavy data remains
under `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/stage-b-20261110-r1`.

Before launch: authenticate pair snapshots/runtime bindings and frozen model files;
run retained contracts; verify Cstar continuity and fresh per-cell majorants;
regress the generalized Jacobian against the accepted Stage A function. During
and after execution: check domain containment, actual crossings, localization,
transversality, both nextafter neighbors, event error, all original accuracy and
ledger gates. Passing preflight does not promote Stage B to numerical qualification.

Classes: documentation, engineering (scheduling), numerical-method verification,
and generated evidence. Authority: GOVERNANCE §§1–5; ADR-0014 controlled chemical
state/signs, ADR-0015 energy ledger, ADR-0016 ownership, ADR-0017 §§4–7 and Amendment 1;
INV-15 heat-capacity ownership. No production equation, interpolation, tolerance,
solver, threshold or frozen input changes. No new stellar solve, master merge,
14,648-year extension, or historical 500 kyr run.
