# Minimum bounded follow-up — not an execution authorization

This plan accompanies [PROPOSAL.md](PROPOSAL.md). Its only proposed protocol
change is selection of the latest positive accepted endpoint at/before half-time
for Stage A quarter-cap restart tests. The patch remains unapplied to the active
runner. Prior thermal retry authorization was consumed; no new trajectory has
been submitted. Original failures remain unchanged.

## Before any trajectory

1. Obtain explicit approval of the precise protocol change. Apply
   `proposed-selector.patch`, annotate the two fixed run plans with the accepted
   Stage A exception and retain their historical rule/status. Do not rewrite any
   historical failure or predeclaration. Commit the approved runner/configuration.
2. Build the stage harness and registered selector test with the retained EKU
   toolchain. Run the selector tests and the existing configuration/static
   contracts. The production library, O2 tolerances (5e-27,5e-27) MeV, other
   numerical settings, physics and scientific criteria must remain unchanged.
3. Republish/authenticate source and control runtime snapshots because the stage
   source/ELF identity changes. Bind the new helper and runner, preserve old
   snapshots/ELFs, verify exact common configuration and physical/archive hashes.
   This needs no stellar regeneration and no trajectory.
4. Before requesting execution, make the actual-output assessor reviewable:
   accepted restart identity/time/state and fresh-cache checks, original-prefix
   error/energy accounting, reference convergence/interpolation and roundoff,
   actual-checkpoint high-precision chemistry/Echem, thermal uncertainty and
   R18/R20. The current collector is not a complete accuracy certificate. Missing
   justified bounds are UNQUALIFIED; neither selector success nor small raw
   disagreements substitutes for them.
5. Return the exact committed runner/snapshot identities for **one separately
   authorized Stage A attempt**, with source-first ordering and the guards below.
   Approval of the protocol alone does not start it.

These preparation obligations are already present in
`../native-regeneration/RUN_PLAN.md:9`, `:39`, `:49` and
`../native-assembly/NONLINEAR_PLAN.md:118`; they are not new accuracy requirements.

## Smallest retry retaining the fixed qualification matrix

No extra pilot solve or forced-breakpoint experiment is needed. Once separately
authorized and ready, use one Stage A attempt with the original [0,31,557,600 s]
horizon and unchanged matrix. Per role: five independent original-t0 references,
12 mains (three tiers, four caps), 12 passive-schedule repeats and three quarter-cap
restart tails. Select each tail from its own accepted history using the proposed
rule. Do not add a new cap, tolerance level, restart family or horizon.

Run the source first. Stop on any scientific, runtime, missing-bound or restart
coverage failure, preserving partial evidence. Only if the source passes every
required gate may the exactly matched zero-source control run and the paired
comparison/figure be produced. No automatic retries or Stage B.

Keep the original one logical execution CPU, 4 GiB and **30-minute aggregate
Stage A wall ceiling**, single-thread libraries and /mnt/sdd storage
(`../native-regeneration/RUN_PLAN.md:57`). Budget source, its required diagnostics,
then control within that ceiling; timeout is INCOMPLETE, with no automatic
extension. The prior partial run does not prove this full workload fits.
Build/static preparation is separate from the thermal allocation.

A single tail from the archived equivalent BASELINE history could test only
that local restart. It would not complete the original three-tier/two-role matrix,
the actual quarter-cap collection, or the outstanding accuracy accounting.
Likewise archived references may support regression, but their substitution for
fresh required members is **not authorized here**. Avoid an extra coverage-only
thermal run that consumes resources without qualifying Stage A.

Use original state/power/R18/R20 budgets without relaxation, retain every prefix
contribution, and keep source absolute time B(t). Apply the independently
converged chemical-reference diagnostic to actual outputs with its conditional
continuous-thermal-enclosure and refinement assumptions. No prior 13-point PASS
is universal qualification. Identical inactive caps still do not demonstrate
cap convergence (`../native-assembly/NONLINEAR_PLAN.md:74`, `:78–112`).

## Reproduce this delivery without a trajectory

From the candidate worktree, verify `SHA256SUMS`. The following **selector-only**
job creates a fresh scratch copy of the proposed harness patch, compiles the
standalone tests and reads retained data; it cannot call the thermal runner:

```sh
source /mnt/sdd/zaki/CompactStar/linux-transition-20261010/toolchain/environment.sh
export RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
export RF_CHECK_OUTPUT=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/restart-coverage-NEW
mkdir -p "$RF_CHECK_OUTPUT"
sbatch --export=ALL --output="$RF_CHECK_OUTPUT/slurm-%j.log" \
  tools/linux/native_qualification/restart_selection_check.sbatch
```

The executed output was
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/restart-coverage-20261109-r1`;
do not overwrite it. `check/result.json` records input/candidate hashes and exact
state/time values; `tests.log` records the test results. The patch was applied
only under that scratch directory. `git apply --check` on the active worktree
validates applicability without making changes. No commit, runtime snapshot or
source trajectory is promoted by this proposed patch.
