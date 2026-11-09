# Targeted C-CUT assessment

UNVERIFIED SCIENTIFIC CANDIDATE. Starting source:
`878229af791fc787123cc121f3a346c8a29ba470`, branch
`candidate/results-first-bnv-pilot-20261108`, worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`.
The worktree was clean; every other local branch tip is an ancestor. Other
worktrees exist and are unchanged. This is the owner's intended source.

Authority: the owner's Final Targeted C-CUT Assessment instruction; GOVERNANCE
sections 1–5; the accepted restricted OWNER_DECISION.md sections A–D; its
authenticated independent review Appendix A (S, C-TAN, C-CONT, C-CUT); ADR-0022.
The older static-completion records and all thresholds remain immutable.
Change classes: validation numerical-method, engineering, documentation, generated
evidence. No change to production equations, outputs, policy, or certificates.

The Mac literature catalog is unavailable on EKU. This calculation differentiates
the already governed historical definitions and the retained X2 equations; it
does not adopt a new EOS, atmosphere, or scientific source. Relevant source:
`tests/bnv/x2_fixed_step.cpp:106`, `:183`; `x2_uncertainty.cpp:55`;
`tools/linux/ba13_closure/assess.py:58`; `tests/bnv/compare_candidate.py:131`.

## Before numerical execution

Use only retained central and offset outputs at targets 0, 10, 20. Evaluate the
terminal EOS analytically, differentiate the terminal TOV identities, and apply
the fixed-B central adjustment. Retain separate fixed-central and fixed-B
responses. Reconstruct the path-tangent mixed response, with the accepted
log-linear cut-law slope. All postprocessing runs in an EKU Slurm allocation;
no stellar integration or trajectory is authorized by this predeclaration.

Use 70 and 90 decimal digits to check postprocessing stability, without claiming
to increase the precision of the retained 64-bit-mantissa stellar outputs.
Evaluate three central-width levels at N=512 and three grid levels at the narrow
width. Apply S's signed-increment/order diagnostic: 2 +/- 0.5 for central widths,
4 +/- 0.5 for GL4 grids. A rounding-limited or wrong-order result is inconclusive,
not a zero error estimate. A tiny response or ample threshold headroom is not
validation. Exact algebraic identities are separately tested.

For each output report the actual response per unit ln(Ps), normalized using its
historical component scale. Also report the first-order product with 2.1e-12,
explicitly as a diagnostic, never as a finite-cut allowance without its remainder.
Pressure matching at the endpoints is not a measurement of output accuracy.

Completion requires either an analytic finite-transport bound from the governing
equations, or an S-qualified sensitivity/convergence estimate including the
remainder for fixed-B matching and tangent response. It must apply on the full
accepted B interval, including baseline and endpoint. C-CONT's on-law curvature
alone must not be asserted to bound an independent pressure direction. Failure
to establish a needed contribution leaves the existing matrix blocked. A next
calculation, if needed, must name the missing quantity and prospective tests;
no broad stellar refinement campaign follows automatically.
