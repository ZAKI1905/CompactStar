# Stage B evidence and reproduction

All commands refer to the candidate worktree on EKU:
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
The starting source is `c90b2205e553afcd102c7acd5f84b253e0a67f1b`.
Do not rerun the scientific matrix merely to reproduce a report. The archived
histories and diagnostics are inputs to the read-only verification scripts.

## Immutable model and production configuration

`tools/linux/preliminary_pair/run.py` authenticates the original source/control
snapshots in `linux-transition-20261010/parallel-preparation/snapshots`, archive,
certificates, coefficient/tangent files, frozen capacities, runtime libraries,
and compiled input bindings. `stage_b/preflight.py` additionally authenticates
the accepted Stage A thermal/normalization model and computes fresh piecewise
constants. The six retained contract tests and piecewise regression tests must
pass. The original matrix and thresholds remain in
`../native-assembly/NONLINEAR_PLAN.md`; the owner superseded only its obsolete
CPU/time ceilings and accepted the explicitly scoped Stage A conditions.

The high-precision reference lifts actual stored binary64 inputs, uses 110-digit
Decimal and the retained arithmetic padding, and follows the actual piecewise
capacity field from t=0. `reference.py`, `piecewise.py`, `hp.sbatch` retain the
scientific sample grid, all crossings, and propagated bounds. The source/reference
models and native binaries are separate implementations. No solver or equation
change is used to obtain a favorable comparison.

## Retained execution directories

`linux-transition-20261010/stage-b-20261110-r1` contains original references,
uncapped native histories, high-precision references and early operational refusals.
`stage-b-20261110-r2` continues diagnostic collection with sufficient physical cores.
`stage-b-20261110-r3` preserves immutable trajectories from interrupted NFS-bound
workers and supplies the final unified matrix. Completed source BASELINE/REFINED
configuration directories are linked to their original r2 locations, preserving
bytes and provenance. Incomplete diagnostic packets remain in r2. Each migration
has an exact retained-input manifest; no incomplete grid is used for an R20 sum.

Every build records commands, library/source hashes and runtime binding. Every
native dispatch records arguments, snapshot/build hashes and its Slurm job.
`continuation.py` rehydrates completed accepted records, including event metadata;
it does not reconstruct an adaptive process from only the terminal state.
Original-t0 main/reference evolutions are not silently repeated. Operational
pre-dispatch refusals, resource-affinity refusals, and interrupted I/O attempts
remain in their original logs. A scientific refusal is not converted into an
operational retry. The retained `tests/native/pair.hpp` deliberately runs
`MustRefuse` probes for missing certificates and mismatched frozen owners. Their
expected REFUSAL messages remain in successful assembly logs; actual task
completion/STOP records distinguish them from a failed scientific calculation.

## Read-only verification order

Run in Slurm on the healthy node with the pinned interpreter
`linux-transition-20261010/toolchain/python/bin/python`. Each independent source
or control audit uses 48 allocated hardware threads, providing 24 distinct physical
cores, 48 GiB and unlimited partition walltime. Set `RF_STAGE_WORKERS=24`,
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` and include `toolchain/bin`
on PATH. Source the toolchain module setup only on the login node, never inside
the compute-node shell. Job scripts use the explicit interpreter path.

1. `supplement.py OUT source|control` traces every retained reference RK map and
   observes all five references at every distinct main/restart accepted endpoint.
2. `observation_trace.py OUT source|control` traces the actual rk8pd observation
   operation at every scientific and accepted endpoint, asserting exact replay
   agreement. `native_points.py` evaluates native powers at all retained main and
   restart endpoints without evolving a trajectory.
3. `residuals.py OUT 24` and `observation_residuals.py OUT 24` propagate actual-stage
   arithmetic, including the inherited reference prefix. `AuditStepper.hpp` is
   generated as an isolated class with passive Apply begin/end callbacks. The
   source-equivalence test restores those hooks/naming and compares it byte for
   byte to production. Event trials and rejected attempts are recorded explicitly;
   only the selected path contributes as successive evolution.
4. `assess.py OUT 24`, `chemical_response.py OUT`, and `events.py OUT` produce the
   physical curves, analytic response bounds and actual crossing comparisons.
5. `ledger.py OUT 24` assesses original R18/R20 values, all nested-grid terms and
   independent piecewise quadrature/count bounds. `local_reference.py OUT 24`
   checks reconstructed coupled states from their retained exact left endpoints.
6. `accuracy.py OUT 24` combines original reference differences, accumulated
   arithmetic, total high-precision comparison bounds, reconstruction uncertainty,
   and restart-prefix inheritance. Missing terms remain explicitly incomplete.
   `coverage.py OUT` verifies all 24 policies, 24 endpoint schedule repeats,
   96 grids and six restarts, including equality at common nested-grid times.
   `matched.py OUT` adds source/control uncertainties and budgets, reports endpoint
   resolved signals, and records the triangle-inequality proof on the shared
   scientific grid. Each role checks its own adaptive endpoints separately.
   `qualification.py OUT` collects those completed gates, binding each input by
   SHA-256. It does not confer owner acceptance. The earlier `assessment.json`
   remains the physical-curve generation receipt with its then-pending obligations.
7. `plot.py OUT` uses the retained plotting environment
   `linux-transition-20261010/native-p2-plotting-env/bin/python`. The computation
   interpreter itself does not include matplotlib.
8. `seal.py SEAL_DIRECTORY R1 R2 R3` performs two complete SHA-256 byte passes and
   an unchanged inventory check. Full compressed manifests stay with the cluster
   evidence; compact receipts and ordered directory roots are committed to Git.

Run `test_piecewise.py`, `test_apply_residuals.py` and `test_audit_source.py` with
the same interpreter. New verification code and numerical bounds remain candidates
until assessed under GOVERNANCE §5; reproducing them is not owner acceptance.
