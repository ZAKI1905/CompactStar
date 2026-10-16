# X1 independent oracle-anchor result: INCONCLUSIVE

All 18 prescribed solves completed, but the frozen oracle-qualification gates
failed. **No absolute-anchor PASS is awarded.** No unconditional numerical
uncertainty bound or physical-accuracy certificate has been established.

This is a validation-only numerical-method candidate under GOVERNANCE §§2–5
and INV-01/02/03/04/06/14/15, explicitly authorized by the X1 instruction.
Production owners, fitter, historical evidence and certificate status are
unchanged. The 0.5 CPU-hour EOS/jitter continuation remains NO-GO.

## A. Regularization mathematics

**The coordinate transformation is correct on the specified segments.** With
t=-sqrt(H-H_on), dx/dt=-2t; every structure and integral RHS receives this
factor, and its onset endpoint measure vanishes. Neutron and muon onsets are
located analytically, the envelope and surface cuts remain explicit, and the
common-B solve uses its own secant. Neither executable reads production tables,
profiles or Steffen splines. Both use the independent 96-point direct
momentum-integral EOS. The equation-by-equation audit, units, endpoint orientation
and limitations are in [MATHEMATICS.md](MATHEMATICS.md).

The independent checks passed: analytic transformed moments <=2.23e-16
relative error; transformed-versus-original RHS residual zero at the fixed
test states; zero endpoint measure; correct neutron/muon support on both sides;
and compensated-summation stress residual zero. A fixed finite-difference
dP/dH check agreed with epsilon+P to 2.17e-8. Across all regularized shots,
epsilon+P=n_B mu_B residual <=1.34e-15, charge <=1.57e-18, beta equilibrium
<=4.45e-16, proper-volume chain rule <=4.45e-16, and cancellation-normalized
mass-chain identity <=3.33e-16. Evidence:
[mathematical-checks.out](evidence/mathematical-checks.out:1),
[runs.json](evidence/runs.json:1).

The redshifted chemical-potential identity agrees to 2.23e-16 in the fixed
checks. Its metric is constructed from H, so this is an EOS/equilibrium
consistency test, not an independently integrated metric test. The accepted
functional increments for B, L, C and N_mu use Kahan accumulation. RK stages,
EOS quadrature and structure variables remain binary64 without a roundoff
enclosure. The center expansion and zero initial C contribution are inherited
shared approximations. These limitations survive the correct transformation.

The usual Hermite O(h^4) assertion is **not justified uniformly across onset**:
a massive species contributes density proportional to (delta H)^(3/2), leading
to a possible (delta r)^(7/2) term in nu and an unbounded fourth derivative.
Exact endpoint derivatives remove the leading linear-interpolation O(h²) error,
but do not supply the missing smoothness theorem. No frozen band was adjusted
in response to this finding.

## B. Empirical contraction predictions

A1 uses two successive increment ratios from the four prescribed regularized
tolerances. Both must be <=1/3.3, and the last must also be <=0.3 for the frozen
U formula. All five quantities are monotone across this ladder.

| Quantity | First ratio, targets 0 / 10 / 20 | Last ratio, targets 0 / 10 / 20 | A1 |
|---|---|---|---|
| Ltilde_Me | .107956 / .107955 / .107955 | .092281 / .092272 / .092280 | PASS all |
| Cstar coefficient | .569179 / .569082 / .569172 | .139150 / .139051 / .139114 | **FAIL all** |
| R | .098899 / .098937 / .098911 | .059410 / .059543 / .059376 | PASS all |
| M | .086802 / .086849 / .086275 | .070577 / .069292 / .065277 | PASS all |
| neutron edge | .098919 / .098963 / .098934 | .059401 / .059598 / .059257 | PASS all |

Thus 12/15 target-quantity contraction predictions pass; C fails at all three
targets. The first C ratio corresponds to the 2e-11→2e-12 increment divided by
the 2e-10→2e-11 increment. It cannot be discarded merely because the final
ratio is smaller. This is a decisive failure of that frozen empirical forecast.
It does not prove that the solver diverges.

The unregularized four-level ladder combines the retained 2e-10/2e-12 values
with the new 2e-11/2e-13 solves. Its increments, ratios and signs, plus those
for mu_inf and achieved B, are all retained in
[scores.json](evidence/scores.json:1). B ratios are irregular at the shooting
residual scale; the alignment gate, not contraction of B residuals, governs B.
Neither formulation's observed ratios are a general theorem of convergence.

## C. Absolute-anchor predictions and formulation agreement

The table gives final regularized residuals X/X_inf-1. Units are **1e-12**.
Each entry shows residual / original full-precision band, rounded for display.

| Quantity | Target 0 | Target 10 | Target 20 | A2 |
|---|---|---|---|---|
| L | +9.789 / 13.163 | +9.510 / 13.392 | +4.703 / 16.380 | INCONCLUSIVE |
| C | +2.334 / 3.720 | +2.050 / 4.041 | -3.269 / 7.251 | INCONCLUSIVE |
| R | +2.634 / 6.859 | +3.156 / 6.875 | +1.614 / 8.139 | INCONCLUSIVE |
| M | -0.018 / 2.396 | -0.047 / 2.377 | -0.119 / 2.507 | INCONCLUSIVE |
| neutron edge | +0.260 / 6.696 | +0.165 / 6.782 | -1.773 / 8.084 | INCONCLUSIVE |

All residuals happen to lie inside the original bands even without U_emp.
This is descriptive consistency, **not a scored A2 PASS**: the global A1 gate
failed, so the frozen rule stops A2 scoring. All 15 A2 outcomes remain
INCONCLUSIVE; none is scored PASS or FAIL. The model's onset-crossing allowance
remains explicitly **heuristic and post hoc**: 10*1e-10*0.045 km=4.5e-11 km,
sized with knowledge of the observed target-20 pressure offset. Agreement here
does not retrospectively make that allowance an a-priori theorem.

There is a second failed gate. At 2e-13, unregularized-minus-regularized
relative differences have magnitudes about 8.06e-12 for L, 4.98e-12 for C,
1.86–3.20e-14 for R, 1.56–2.82e-14 for the edge, and **3.20–3.47e-13 for M**.
The frozen 10*U_emp(M) limits are only 5.83–6.88e-14. Mass fails this required
agreement check at every target. L,C,R,edge meet that check descriptively.
This test uses the regularized estimate alone; the failure does not identify
which formulation supplies the discrepancy or prove a missing physical term.
No new allowance for the unregularized method was introduced.

U_emp=|last increment|rho/(1-rho) is reported only as a **conditional empirical
tail estimate**. For C it is about 1.24e-12, numerically below B_C, but A1's
failure prevents treating that number as an eligible uncertainty estimate.
Mass illustrates why even small apparent tail estimates cannot stand alone.
For mu_inf and B, no frozen 10*U band exists; their disagreements are reported
without inventing one.

## D. Cstar decisiveness and frozen directions

C decreases from the retained unregularized 2e-12 value by
5.01388e-11 / 5.01196e-11 / 5.01381e-11 at targets 0 / 10 / 20.
**The C direction prediction passes at all three targets**, with a change near
the proposed 4–5e-11. L and M also decrease and the neutron edge increases:
all 12 frozen A2s sign tests pass. R had no sign prediction. There is no
wrong-direction prediction failure in this experiment.

C is decisive as a direction observation and as a failed A1 forecast, but
**not decisive as an absolute-anchor test**. Its final residual falls within
the original band; the required uncertainty qualification does not pass.
Remaining absolute uncertainty is not rigorously bounded. A numerical sign
test PASS is not a proof of the sign of the exact physical error.

## E. Paired depletion responses

All new solves use the original physical B targets, including the individual
retained surface pressures. Maximum achieved |B/B_target-1| is 7.87e-14,
below 5e-13. Both variants' finest common-B disagreement is <=1.31e-14.
These are residuals of the oracle B functional; they do not bound its common
quadrature/center error. Original production B targets and achieved B are
preserved in [paired-production.json](frozen/paired-production.json:1).

All **72** variant/tolerance/target/quantity comparisons to the retained oracle
pass the unchanged 1e-12 A3 limit; the maximum difference is 3.52e-14.
All **432** corresponding comparisons to the six Hermite-production phases
pass the unchanged 1e-11 limit; maximum difference is 5.13e-12.
These are separate paired diagnostics despite A1 failure, as predeclared.

| Final regularized relative drift | Target 10 | Target 20 |
|---|---|---|
| L | +9.046047777e-7 | +1.809209155e-6 |
| C | +7.350082272e-7 | +1.470015425e-6 |
| R | +6.415312044e-7 | +1.283061643e-6 |
| M | -4.753428654e-7 | -9.506857648e-7 |
| neutron edge | +6.287134494e-7 | +1.257426111e-6 |
| mu_inf | +9.413198726e-8 | +1.882637831e-7 |

Absolute discrepancies are measured against a model anchor, not a known exact
solution. Near-common offsets cancel deterministically in ratios. For example,
the finest formulation disagreement in L is about 8e-12 at every target, while
paired solver differences are much smaller. `anchor_error_decomposition` in
the score record separates the target-zero common offset from target-dependent
residuals. Agreement of paired responses therefore does not establish absolute
accuracy. Small residual targeting differences are measured, but no additional
derivative bound was constructed to rigorously apportion their effects. The
remaining discrepancies may combine solver, shared EOS, center and anchor-model
errors; this experiment does not uniquely attribute them. No 21-star enclosure
follows from these three targets.

## F. Bounds actually established

**No unconditional global numerical uncertainty bounds.** The evidence consists
of algebraic identities, finite floating-point checks, observed convergence,
prospective sign outcomes and deterministic paired comparisons. There is no
validated integration remainder, global stability/roundoff enclosure, certified
center remainder, bounded EOS quadrature/inversion error or certified common-B
root sensitivity. Empirical geometric contraction cannot substitute for these.
The regularization is mathematically correct but the regularized solver remains
an uncertified absolute reference. Historical BA13 qualification is unchanged.

## G. Smallest next step

Predeclare a **separately implemented fixed-step, known-order solver in s**,
for example Gauss-Legendre collocation, with a fixed step-halving schedule at
these same B references and explicit center-series checks. C's nonuniform
contraction and the mass formulation discrepancy are the immediate targets.
Keep this experiment's failures and bands frozen. No such new solve was run.
For unconditional certification, additionally justify a global integration
remainder and propagation through EOS, floating point, center/surface and the
common-B root, or use validated integration. Step-halving agreement alone would
still be empirical evidence.

## Authentication, resource receipts and reproducibility

* Entry: `dad05b17056f3ccb1e3871e6e4bb24ebe5db5df4`; independent review SHA256SUMS
  verified before and after, including all calculation artifacts and anchor.json.
  All 650 retained pilot artifacts authenticate. Review files were only read.
* Prediction freeze: `91e8254`. Solver/scorer gate: `d67a10e`. First numerical
  execution: `d2796ac20b24acee681dfc301926e7a90cbb10d3`, after a pre-result
  launcher correction. No solver or scoring change followed observed results.
* Slurm 17021 stopped before any numerical work because compute-node PATH lacked
  Git. Its receipt is retained. Slurm **17022 completed**, exit 0, node phyn004,
  all 18 prescribed solves, all mathematical checks and scoring.
* The independent 60-digit Decimal audit passed under Slurm **17023**, checking
  all A1, A2s, formulation-agreement and A3 classifications against raw outputs,
  the unchanged A2 gates/bands, and source authentication. It added no stellar
  solve. Receipt: [independent-audit.json](evidence/independent-audit.json:1).
* Measured numerical run **3.04996 CPU-s**, including controller and children;
  independent audit **0.087927 CPU-s**; combined **3.137887 CPU-s**. Maximum
  per-process measurement **0.267799 CPU-s**; peak RSS upper bound **37,732 KiB**.
  Slurm accounting storage is disabled; getrusage receipts and Slurm allocation/
  final-state records are retained. Slurm rounded allocation to two hardware
  threads, but every numerical process was pinned to logical CPU **4**, inherited
  by all children, with single-thread environment variables and a 1 GiB limit.
  Thus execution used one CPU; the two-thread allocation is disclosed, not
  reported as a one-thread Slurm allocation. All ceilings were respected.
* Distinct executable SHA256 identities: unregularized
  `ce88435ad795afa9097da11500d3fa6d652b9bb20ea880c6c8ee24c470b584d3`;
  regularized `5dcc971ff9d20627a6052d297b444c0e99ec3f15c147d5366ea366ddf4a36dc1`.
  Exact source, dependency and input hashes are in
  [execution-gate.json](execution-gate.json:1). Numerical executables remain
  outside Git; the full external artifact manifest is committed.
* Changed files are only the new validation oracle, X1 tools, and X1 documents/
  evidence. Original source, fitter and P5/P7 failure evidence are unchanged.
  No production solve, certificate promotion, master merge or campaign occurred.

Reproduction and Mac transfer instructions: [REPRODUCE.md](REPRODUCE.md).
The code remains an UNVERIFIED SCIENTIFIC CANDIDATE pending human review;
neither this mathematical audit nor a commit constitutes ratification.
