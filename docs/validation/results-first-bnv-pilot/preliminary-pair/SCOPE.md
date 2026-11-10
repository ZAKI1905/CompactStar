> Status update, 2026-11-10: the owner accepted the restricted Stage A
> CONDITIONAL PASS after independent review; see `OWNER_ACCEPTANCE.md:1`.
> Stage B is now authorized subject to runtime and piecewise numerical preflight.
> The original execution scope below is retained as historical provenance.

# Exploratory matched control and source evidence recovery

Owner authorization: the 2026-11-10 request explicitly permits a one-year zero-BNV
control before source qualification, isolated source reconstruction, additional
high-precision calculations, and efficient Slurm parallelism. This exception is
exploratory only. Formal Stage A criteria remain `../native-assembly/NONLINEAR_PLAN.md:95`.
Stage B and the historical campaign remain unauthorized.

Entry: clean branch `candidate/results-first-bnv-pilot-20261108`, HEAD
`191fac300ec3770bc79f505e3afd606075fc65c3`, in
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
All other local branch tips were ancestors; nine worktrees were present.

Classes: engineering (evidence and scheduling harnesses), numerical-method
(independent offline high-precision assessments), documentation and generated
evidence. Production equations, libraries, solver tolerances and scientific
thresholds are unchanged. Authority: GOVERNANCE sections 1–5, ADR-0014 controlled
chemical conventions, ADR-0015 energy accounting, ADR-0016 tangent ownership,
ADR-0017 sections 4–6 and Amendment 1. All new work is an unverified scientific
candidate, not human ratification.

The control uses the authenticated existing pair snapshot. Five original-t0
references run independently; three main tiers each retain all four cap policies,
four passive grids, endpoint-only repeats, and the quarter-cap restart. Parallel
tasks own distinct output files and mutable solver histories. Only scheduling
loop filters are applied to a generated copy of the retained stage harness.
The original stage source and libraries are not edited (`tools/linux/preliminary_pair/run.py`).

Source recovery first tests 16 checkpoints spanning initial/final endpoints,
early chemical-energy sensitivity, and late thermal sensitivity. Four packet
files must reproduce exactly. Full coverage then visits the remaining unique
1920-grid points, with an explicit mapping to every retained packet. Each local
witness uses the existing breakpoint stepper and prior /4,/8,/16 rk8pd/RKF45
settings (`tools/linux/preliminary_pair/recover.cpp`).

High-precision chemical assessments reuse the existing 110-digit model and
knot-aligned RK4 at 8,16,32,64 steps per Hermite segment. Uniform thermal-band
sensitivity uses the existing directed-rounding Gronwall calculation. Agreement
at sampled nodes does not prove the thermal band. Continuing refinement
contraction is an explicit assumption whenever twice a difference estimates
remaining discretization error. Original-t0 coupled high-precision comparisons
will separately evaluate reference/main/restart discrepancies with the retained
binary64 Cstar cells and physical coefficients. Evaluated photon and deposition
anchors require explicit normalization uncertainty; they do not independently
revalidate the stellar background or the static certificate.

All figures and reported physical comparisons remain
**PRELIMINARY — NOT NUMERICALLY QUALIFIED** until all required terms are justified.
Measured differences, conditional estimates and proved arithmetic statements
must be reported separately. Unknown positive uncertainty terms are never zero.
