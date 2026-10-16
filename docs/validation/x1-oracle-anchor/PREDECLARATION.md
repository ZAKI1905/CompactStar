# X1 independent oracle-anchor qualification — prospective gate

UNVERIFIED SCIENTIFIC CANDIDATE. Validation-only numerical-method change,
with engineering, documentation and generated evidence. Authority: explicit
2026-10-16 X1 user instruction; GOVERNANCE §§2–5; INV-01/02/03/04/06/14/15
(density, units, metric, proper volume, finite-pressure surface, baryon number,
thermal meaning). ADR-0021 remains PROPOSED. No production owner is changed.

Source: `dad05b17056f3ccb1e3871e6e4bb24ebe5db5df4`, clean onset-pilot
worktree. All existing local branches are ancestors. Separate descendant branch
`candidate/linux-x1-oracle-anchor-20261016`, worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/x1-oracle-anchor`.
The independent review is read-only throughout. Its complete SHA256SUMS verifies;
byte-identical copies, including anchor.json and its generating scripts, are in
`frozen/`. Retained oracle metrics and resource input strings match Git evidence.
`authentication.json` records source hashes; `targets.json` freezes exact physical
B, initial rho guesses and individual finite surface pressures.

## Frozen experiment

Unchanged unregularized source: `tests/bnv/onset_structure_oracle.cpp` and its
two independent EOS/TOV headers. New source: `tests/bnv/x1_regularized_oracle.cpp`.
No production table, spline or stellar profile may enter either solver.
Targets exactly 0, 10, 20. New unregularized tolerances exactly 2e-11, 2e-13;
regularized exactly 2e-10, 2e-11, 2e-12, 2e-13. Retained unregularized 2e-10
and 2e-12 are used without rerunning. No additional stellar tolerance or target.
Run sequentially on one Slurm CPU, memory <=1024 MiB. Numerical work, including
checks and scoring, <=120 measured user+system CPU seconds. Each solve has a
60 CPU-second hard limit; halt at first failure, excess, authentication mismatch
or |B/target-1| >5e-13. Reserve time before starting any subsequent solve.
Build/hashing bookkeeping is separately receipted. Hash executed sources,
executables and scripts in a second committed gate before submission. No edits
to solver or scorer based on numerical outcomes; any defect invalidates affected
evidence and is reported, without an unbudgeted replacement experiment.

## Original predictions and classifications, unchanged

The normative numerical anchor is the full-precision
`frozen/anchor.json:targets`, fields `predicted_analytic_limit` and
`band_total_excluding_oracle`; rounded review tables are explanatory only.
The copied review §9 preserves original words, including its overclaim of
"certification", which this experiment explicitly does not adopt.

* A1: each target and each L,C,R,M,r_e has three successive increments from
  four regularized tolerances. Both final ratios of absolute increments must
  be <=1/3.3, and the last ratio rho must be <=0.3 as separately specified by
  the original U formula. Record signs and monotonicity separately. Differences
  use a common final-value denominator; zero denominator in a contraction
  ratio is undefined and cannot PASS. Record FAIL for a violated prediction;
  the experiment and A2 are then INCONCLUSIVE. A1 is an empirical test only.
* U_emp = |last increment| rho/(1-rho), relative to final value. This is a
  conditional geometric-tail estimate, not a mathematical global bound.
* A2: globally gated on A1 and completion/alignment/independence checks. For
  each quantity, U_emp > B_X means INCONCLUSIVE. Otherwise the original test
  |X_final/X_inf-1| <= B_X+U_emp is empirical PASS, else FAIL. No bands grow.
  If A1 fails, show residuals descriptively but do not score A2.
* A2s is separately scored even if A1 is inconclusive (explicit user request):
  from retained unregularized 2e-12 to regularized 2e-13, C,L,M must decrease,
  r_e must increase. Wrong sign is recorded FAIL regardless of size. Zero
  motion is INCONCLUSIVE. The expected C decrease is about 4–5e-11 and M
  decrease about 1e-12; those approximate magnitudes are descriptions, not
  new acceptance bands. R has no frozen sign prediction.
* A3 is separately diagnostic if A1 fails: for each variant/tolerance and
  target 10 or 20, drift = X_i/X_0-1, X=L,C,r_e,R,M,mu_inf. Difference from
  retained 2e-12 drift <=1e-12; from each of the six retained Hermite-phase
  production drifts <=1e-11. Missing original phase evidence is INCONCLUSIVE,
  never silently omitted. No inference for the other 18 targets.
* Independence: regularized/unregularized at 2e-13 must differ by <=10 U_emp
  for the five A1 quantities. No new U or band is invented for mu_inf or B.
  Failure makes overall qualification INCONCLUSIVE; retain raw disagreement.

PASS means prospective empirical consistency of the candidate error model.
FAIL means a frozen prediction failed; INCONCLUSIVE means insufficient oracle
uncertainty control or incomplete run. Neither PASS nor contraction establishes
unconditional physical accuracy or a historical certificate.

## Heuristic provenance and mathematical limits

The original production band contains 10*rtol*H_P = 4.5e-11 km, with
rtol=1e-10 and H_P=0.045 km, mapped to L,r_e,R. The factor ten was chosen
with knowledge of the target-20 observed 4.57e-10 ln(P) crossing offset.
It is post hoc and heuristic, not a crossing-error theorem. No component
is enlarged or refitted here. The original P5 and P7 failures remain failures.

Regularization uses s=sqrt(H-H_on) on the high-density side, with H=H_on+s^2,
dx=-2s ds for x=Hc-H, and zero transformed endpoint measure at s=0. Both
neutron and muon boundaries must be explicit. Compensate accepted increments
of B,L,C (and N_mu); preserve the RKF45 tableau and the original error scales.
Independent algebra/analytic-moment tests and EOS thermodynamic identities
must check the transformation, orientation, support and measure. These are
not a second stellar solver. Check mu_B exp(nu)=constant and proper-volume
identities. Keep center initialization and EOS quadrature limitations visible.

A bounded global error needs certified regularity and derivative/remainder
bounds, stable error propagation, EOS inversion/quadrature and floating-point
enclosures, center/surface treatment and common-B inversion sensitivity.
An independent fixed-step known-order solver with predeclared step halving
would be a useful next cross-check, but its empirical ratios alone also do not
prove an unconditional bound. Classical Hermite O(h^4) requires appropriate
fourth-derivative control; investigate onset smoothness rather than assume it.

No 0.5 CPU-hour EOS/jitter continuation, production solve, owner change,
fitter modification, certificate promotion, 21-star expansion, BNV trajectory,
campaign, or master merge. Commit compact evidence; large artifacts outside
Git with hashes. If GitHub authentication is absent, produce and verify an
incremental bundle for the office Mac.
