# Reproducing the preliminary pair

All substantial data and builds are on EKU under `/mnt/sdd`; no Mac/laptop build
or data copy was made. Start from commit `191fac300ec3770bc79f505e3afd606075fc65c3`
plus the candidate commit containing this directory. Production files and
retained numerical archives are unchanged.

```bash
RF_SOURCE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot
RF_BASE=/mnt/sdd/zaki/CompactStar/linux-transition-20261010
RF_OUTPUT="$RF_BASE/preliminary-pair-20261110-r1"
RF_PY="$RF_BASE/toolchain/python/bin/python"
RF_PLOT_PY="$RF_BASE/native-p2-plotting-env/bin/python"
cd "$RF_SOURCE"
```

The exact native compile/link commands, libraries, runtime bindings and hashes
are in the external `build/build.json`, `recovery-build/build.json`, and the
other `*-build/build.json` or completion receipts. Native tools reuse retained
Release archives with `-ffp-contract=off`; they do not rebuild the solver.
Snapshot validation authenticates both roles, the original stellar archive,
coefficients, frozen-domain inputs and executable dependencies. Source input
paths are intentionally bound to the completed retained source campaign.

## Read-only evidence verification

Check the manifest SHA against `evidence/external-receipt.json`, then verify all
external files. The receipt gives the absolute manifest path. The compressed
manifest is outside the run directory to avoid self-reference.

```bash
RF_EVIDENCE="$RF_BASE/preliminary-pair-evidence-20261110-r1"
sha256sum "$RF_EVIDENCE/external-SHA256SUMS.gz"
gzip -cd "$RF_EVIDENCE/external-SHA256SUMS.gz" | sha256sum --check --quiet
cd "$RF_SOURCE/docs/validation/results-first-bnv-pilot/preliminary-pair/evidence"
sha256sum --check --quiet SHA256SUMS
```

The original 217,761-file source manifest is checked before reconstruction and
again after all work. `coverage.json` maps every source packet to the recovered
unique checkpoint; `control-inventory.json` separately audits the control.

## Calculation order on a fresh output root

The native builders and collectors refuse existing destinations. Preserve the
recorded root; use a new `/mnt/sdd` directory for another execution. Use Slurm
throughout. Native pools require actual physical-core affinity: for a 24-worker
task allocate `--cpus-per-task=24 --hint=nomultithread` and invoke
`srun --ntasks=1 --cpus-per-task=24 --cpu-bind=cores`. An allocation alone does
not give a default one-CPU srun step 24 cores. Five/seven-hour time limits were
available; no shortened scientific horizon or resource ceiling was imposed.

1. `run.py authenticate OUT` and `run.py build OUT` authenticate the retained
   inputs and build the scheduling-only control harness. `run.py reference OUT I`
   for I=0..4 runs independent native references. `run.py tier OUT I` for I=0..2
   runs each complete four-cap main/repeat/restart matrix with 24 physical workers.
   These independent jobs run concurrently. `job.sbatch` supplies the wrapper;
   its RF_MODE is `reference` or `tier` and the array index supplies I.
2. `recovery.py build OUT` compiles the isolated reconstruction tool.
   `recovery.py sample OUT` checks the predeclared 16-point sample; only after
   its PASS, `recovery.py full OUT` processes the remaining 1,905 points.
   `recovery.sbatch` uses RF_SOURCE, RF_OUTPUT and RF_MODE. Both calculations
   use 24 physical workers. `coverage.py OUT` verifies all 43,248 source packets.
3. `export.py OUT` exports the fixed heat-capacity cells and diagnostic anchors.
   `high_precision.py global OUT --role ROLE --n N --precision P` runs both roles
   at N=1,2,4,8 with P=110 and N=8 with P=80. These ten jobs are independent.
   `high_precision.py sample OUT --workers 24` precedes
   `high_precision.py full OUT --workers 48`; allocate corresponding CPUs.
   The original chemical model and knot-aligned integrator are reused.
4. `history_trace.py OUT build`, then `history_trace.py OUT source` and `control`
   generate six independent retained-history arithmetic traces per role. Each
   role needs six distinct physical cores. `roundoff.py OUT ROLE NAME` assesses
   NAME=`main` or `reference-0`..`reference-4`; run these twelve tasks independently.
   The final script explicitly checks every stage/endpoint enclosure's domain.
5. `thermal_enclosure.py OUT 24` checks all local continuous thermal paths and
   coupled endpoint errors using the analytic remainder in `ENCLOSURE.md`.
   `normalization.py OUT` exports fixed geometry operands and evaluates the
   primary photon/direct-power formulas at high precision.
6. `interpolation.py OUT build`, then `interpolation.py OUT source` and `control`
   reproduce all 9,605 observations per role using 24 physical workers each.
   `accepted.py OUT` observes the five retained reference histories at every
   main accepted endpoint. These are bounded local calculations, not source
   original-t0 reruns. `interpolation_roundoff.py OUT 48` then applies the RK
   arithmetic recursion to all 19,210 observations with their actual accumulated
   reference prefixes, using 48 physical workers. It must finish before the
   final assessment; observed 64-ulp floors are never substituted for this proof.
7. After all control tiers finish, `final.sbatch` runs `inventory.py`, `assess.py`
   and `supplement.py`. Run `summarize.py OUT` after these pass. These commands
   write assessments but do not integrate production trajectories. `--partial`
   exists only for interim assessments and cannot produce a final PASS receipt.
8. `compare.py OUT DESTINATION` under RF_PLOT_PY creates the temperature CSV,
   PNG/PDF figures and measured comparison. Use the final supplement so its
   conditional temperature uncertainty appears in the figure caption.
9. After every writer has exited, `seal.py OUT EXTERNAL_EVIDENCE_DIRECTORY`
   authenticates the source again, hashes all run files with 24 workers, then
   independently rereads and verifies every file. Keep its Slurm log outside OUT.
   Commit the compact receipt, assessments, figures and reproduction sources;
   retain large traces, packets, builds and the compressed manifest externally.

For this execution, source/control native references used original tolerances,
caps, physical state and independent histories. `control_stage.inc` is generated
from `tests/native/stage.cpp` with only reference/tier loop filters. Its build
receipt records every substitution. The reference and roundoff sbatch examples
contain this run's explicit EKU paths; adjust those paths when choosing a fresh
root. Do not overwrite completed evidence.

## Preserved unsuccessful attempts

Snapshot reads on phyn001/002 stalled on filesystem metadata. Those preparations
were canceled and their directories/logs retained; successful native assemblies
ran on phyn003/004. Initial pool launches failed the existing physical-affinity
guard because srun inherited one core; the corrected wrapper requests its full
allocation explicitly. One recovery compile failed on a string/`const char*`
argument and was corrected before numerical execution. A shared history step
ended with its parent allocation; a subsequent six-logical-CPU launch exposed
the distinction from six physical cores. Both attempts remain external. Two
I/O probes failed because `/usr/bin/time` was absent. These are infrastructure
or evidence-harness failures, not hidden scientific retries. No completed source
trajectory was repeated, no prior failure was deleted, and no threshold changed.

Supplemental job 19280 stopped at an exact-zero assertion: two algebraically
equivalent 110-digit equilibrium-luminosity sums left a tiny subtraction residue
in the independent DeltaLnu diagnostic. The evaluator now applies the exact
polynomial identity DeltaLnu=0 at eta=(0,0); the native output never changed.
Job 19281 passed the strengthened contraction test using computed arithmetic
enclosures alone, without a 64-ulp allowance as proof. Final assessments before
the added observation arithmetic are preserved in `before-interpolation-arithmetic/`.
The final read-only assessment on phyn003 later stalled in filesystem access;
job 19283 and its unstarted dependent seal were canceled, with their logs retained,
and the assessment/seal were moved to phyn004. This repeated no native trajectory.
Before final sealing, the figure's displayed uncertainty was rounded outward
from 2.8e-8 to 2.9e-8 K. The earlier figure, script and summary are preserved in
`plot-before-outward-rounding/`; the interrupted artifact-only attempt is retained
as `seal-attempt1.log` and `post-authentication-attempt1/`. No numerical assessment
or native trajectory was repeated for this display correction.

Formal status remains **PRELIMINARY — NOT NUMERICALLY QUALIFIED**. Reproducing
these results does not ratify the new enclosure argument under GOVERNANCE §5.
