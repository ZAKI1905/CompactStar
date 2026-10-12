# Bounded qualification execution scope

Pre-scoring contract commits: bbe45ad (ADR-0020), d5125fd (width audit).
Engineering candidate: 44d51e5, with subsequent validation/handoff commits.

Slurm jobs initially submitted 2026-10-11:
* 14463: incremental Debug/Release builds; 13 readiness tests; PE/U11;
  arbitrary-precision oracle; local EOS; original manifest/coupled failures.
* 14464: CONFIND exact payload serial and native repeat/worker controls;
  Zaki original suite (Mac-gated, not numerical coverage); explicit independent
  Linux log-axis/libm attribution checks and held-out scalar tests.
* 14465 afterok 14463: static Phase-5B, Hartle/CMF and response independent
  oracles; fresh EOS/profile refinement, coefficients and exact repeat;
  conditional componentwise uncertainty comparison. Empirical bounds only.
* 14466–14468 afterok 14463: actual arithmetic/runtime receipts on phyn002,
  phyn003 and phyn004. phyn001 is drained, no inference about it.

All test commands are enumerated in tools/linux/corrections_*_job.sh. Structural
and coefficient jobs use static star solves only. Test horizons were inspected:
no passive_cooling_regression, grid_convergence_cmf, full controlled evolution,
1 Myr, 10 Gyr or 500 kyr campaign command is included. Historical protected
manifest/coupled_oracles fail before any trajectory on the unchanged historical
pins; their original exact policy is intentionally preserved.

Certificate-dependent source/control/event pilots are not launched: the original
21-star producer is missing. The initial candidate already recorded correct
exact-profile rejection; rebinding Mac sensitivities would invalidate the science.
R18/R20, nonlinear stability, event/diagnostic uncertainty and source-minus-control
qualification therefore remain pending a coherent native fixture. Linear Jacobian
diagnostics do not confer campaign or nonlinear stability authority.

Dependency overlay is on a separate CONFIND candidate branch. The accepted
b0cbd510 build/install artifacts and lock remain byte-unchanged; the overlay runs
them without replacing dependency master or silently updating the lock.

Job 14469 adds the separately committed Zaki integral overlay: 36 values in each
build mode, independent quadrature and cancellation-aware absolute bounds.
The runtime receipt includes actual loaded libm/libc/loader/libstdc++ bytes,
CPU dispatch features, rounding mode and MXCSR controls, dispatch environment,
FP widths, compiler commands and full locked GSL/zlib closure. The comparison
requires three observed nodes with exact arithmetic/payload agreement; CPU model
names alone are never the criterion. Receipts are a proposed runtime class,
not a silent mutation of either dependency lock.

GitHub starting branch was independently observed at the exact f4e7c18 hash.
The existing HTTPS push dry-run failed: no available username/credential and
terminal prompts disabled. No credentials were changed. An incremental bundle
handoff for the office environment is required, including dependency branches.

Job 14470 adds 16 held-out CONFIND stellar-field cells in Debug/Release: exact
raw intersection counts, independent scalar libm checks and rational-coordinate
error propagation. This supplements the near-tie case without changing the
contour algorithm or accepting arbitrary global topology drift.
