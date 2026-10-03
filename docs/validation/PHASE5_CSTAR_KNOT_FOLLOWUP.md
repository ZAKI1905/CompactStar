# Separate follow-up: Phase-5 Cstar-cell integration

Status: **RECORDED / NOT IMPLEMENTED / NO NEW PHASE-5 AUTHORITY** (2026-10-03).

The Phase-6 obs-32 investigation identifies a numerical issue upstream of the
checkpoint layer: integrating the existing piecewise-smooth Cstar interpolation
across cache knots can contaminate an otherwise high-order RK step. Phase-5 uses
the same thermal cache, but this task leaves
`CompactStar/Physics/Rotochemical/ScaledRKF45.hpp`, the Phase-5 numerical paths,
EOS values, interpolation and governed baselines unchanged.

A separate owner-authorized task should first reproduce and quantify the issue
on Phase-5's own trajectories and baselines. It should then decide whether to
adopt cache-cell segmentation, with its own event solves, uncertainty propagation,
regression qualification and acceptance. The Phase-6 bounded tests do not qualify
Phase-5 trajectories or justify silently replacing their accepted baselines.
No cross-project migration, external dependency work or Phase-5 rerun is authorized
by this follow-up record.

The shared Analysis file-reader changes in the Phase-6 candidate are engineering
only: every source byte is still read and compared on every full-currentness call.
They neither repair this knot issue nor change Phase-5 numerical mathematics.

Related records: [candidate amendment](PHASE6A1_PIECEWISE_RECONSTRUCTION_PROPOSAL.md),
[successor validation specification](PHASE6_BREAKPOINT_VALIDATION_SPEC.md), and
[historical obs-32 forensics](PHASE6A1_OBS32_FORENSICS.md).
