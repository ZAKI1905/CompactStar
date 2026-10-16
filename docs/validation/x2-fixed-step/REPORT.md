# X2 independent fixed-step oracle: empirical PASS, reference qualification INCONCLUSIVE

The independently implemented GL4 solver passes the declared formulation and
step-halving tests at targets 0,10,20. The two X1 failures now have reproducible
numerical explanations at target 0: an onset endpoint support error changes
Cstar's adaptive contraction, and a discrete mass-chain defect dominates the
regularized/unregularized mass discrepancy. **No rigorous global reference bound
or BA13 qualification is awarded. All original X1 failures remain unchanged.**

This is a validation-only numerical-method candidate under GOVERNANCE §§2–5,
INV-02/03/04/06/14/15, ADR-0004 and ADR-0009, explicitly authorized by the X2
instruction. It is an UNVERIFIED SCIENTIFIC CANDIDATE pending human review.
Equations, assumptions and residual construction: [MATHEMATICS.md](MATHEMATICS.md).

## Authentication and prospective boundary

Entry worktree `/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/x1-oracle-anchor`
was clean on `candidate/linux-x1-oracle-anchor-20261016`, HEAD
`6cc7e20713e3c0210faa9036316650816f0290d5`, the requested source of truth. All ten
existing local branch tips were ancestors of that commit; no competing modification
to an X2 path existed. The separate descendant is `worktrees/x2-fixed-step`, branch
`candidate/linux-x2-fixed-step-20261016`.

All 48 X1 manifest entries, 14 frozen review entries, 41 execution inputs,
15 independent-review inputs and 62 external artifacts authenticated. The source,
raw outputs, audit, frozen predictions and failure dispositions were checked again
at delivery. [authentication.json](authentication.json) and the independent audit
retain those identities. X1's three C A1 failures, three mass-agreement failures,
and fifteen A2 INCONCLUSIVE outcomes remain exactly as recorded.

Predeclaration `75d8c02` preceded every X2 numerical run; implementation/execution
freeze `d313cc3` preceded pilot 17024. Prior X1 observations were explicitly known.
Withheld 2048/4096 grids were specified initially and released at `bc3b851` because
the base ladder could not yet resolve the X1 mass scale. The shell mechanism tests
were a separately labeled post-trace hypothesis, frozen at `bd8c279`/`2137846`
before their results. They introduce no new tolerance or resolution and do not
retrospectively change an initial prediction. See the committed predeclarations.

## Independent solver correctness and convergence

No old adaptive driver, LocalOracle implementation, production EOS table or
interpolator enters the new solver. Only authenticated physical constants are
shared. Fresh long-double momentum quadrature, equilibrium inversion, implicit
Gauss stages, center expansion, common-B secant and off-stage checks are implemented
in [x2_fixed_step.cpp](../../../tests/bnv/x2_fixed_step.cpp).

The smooth manufactured problem measures order 4.00025. Center residual orders
are 2.00004 and 0.999982 as predicted. Thermodynamic/charge/beta residuals are
below 4.7e-17, the finite-difference dP/dH check is 4.47e-10, transformed moments
are within 2.2e-19, and both species support checks pass. The analytic onsets and
positive target-specific surface are explicit segment boundaries. Center z,
onset s, and shell sqrt(H) maps apply their Jacobians to ALL equations.

N is the number of fixed steps in EACH of five segments. The complete ladder
32,64,128,256,512,1024,2048,4096 ran for exactly three targets. The final observed
orders and conditional relative Richardson tails are:

| Quantity | Final p, targets 0 / 10 / 20 | Largest conditional relative tail |
|---|---|---:|
| M | 3.99982 / 3.99949 / 3.99990 | 9.79e-16 |
| R | 3.99993 / 4.00014 / 4.00002 | 1.49e-15 |
| neutron edge | 4.00000 / 4.00009 / 3.99996 | 1.46e-15 |
| Cstar coefficient | 3.99996 / 3.99999 / 3.99990 | 2.31e-15 |
| Ltilde_Me | 4.00000 / 4.00007 / 3.99989 | 2.46e-15 |

All 24 scalar/Hc order tests pass the predeclared two-consecutive-order and sign
criteria (includes rb and mu_inf). Achieved B is a root-alignment diagnostic,
**INCONCLUSIVE for order**, not a quantity expected to contract. Its maximum
relative alignment residual is 5.08e-17, below 2e-16. This is alignment to the
computed B functional, not a bound on its exact value.

At N=4096 the independent off-stage mass closure is <=9.51e-17, baryon closure
<=9.30e-16, and independently integrated metric log-mu drift <=1.30e-18. They
contract from the pilot. The separate stage-mass closure is <=1.23e-15. These
PASS the declared residual tests; they are not equivalent to global error bounds.
Full signed increments and all values: [scores.json](evidence/scores.json),
[runs.json](evidence/runs.json), and exact decimal output strings in
[raw-decimal-text.json](evidence/raw-decimal-text.json).

## Cstar: explain the failed initial contraction without erasing it

At target 0, the X1 regularized original C ratios are 0.569179 and 0.139150.
The first remains a frozen FAIL. Holding Hc fixed changes them to 0.309746 and
0.129110: common-B shooting affects the observed cancellation, but does not by
itself identify the responsible integration component.

The decisive trace is the sub-neutron shell. Neutron-only JL must be CONSTANT
there, yet X1 regularized adds 1.00113944e-4, 1.0625582e-5, 9.84712e-7 and
1.04239e-7 at its four historical tolerances. At the retained Hc,
`Hc-(Hc-Hn)-Hn = +7.914675859e-18`: the first lower-side RHS evaluates `H>Hn`
and enables JL. The EOS neutron density can still round to zero, so the support
predicate and EOS support disagree. The high-side s Jacobian does not regularize
this first point of the subsequent unregularized segment.

Separate diagnostic copies, all at the SAME historical tolerances, isolate it:

| Target-0 diagnostic | First / last common-B C ratio | Finest C relative to X2 N=4096 |
|---|---|---:|
| Original X1 regularized | 0.569179 / 0.139150 | +1.3600e-12 |
| Segment-topology JL support | 0.090032 / 0.203207 | +1.3527e-12 |
| Also remove rounded lower-side neutron density | identical to previous | +1.3527e-12 |
| Also use sqrt(H) in the outer shell | 0.107757 / 0.069953 | -2.9744e-13 |

All 24 support-corrected shot/fixed-H traces have EXACT zero sub-neutron JL
increment. Removing a rounded neutron EOS contribution has no further effect;
phantom-neutron C is therefore not supported as the explanation. Changing JL's
incorrect endpoint activation changes the adaptive mesh and removes the irregular
initial C contraction at target 0. This is an integration/controller effect
triggered by a boundary support predicate, not evidence of divergent C physics.
The remaining original C error is concentrated in the shell; its coordinate
change reduces that error and leaves an interior/other-component residual.
The shell-s mass result does not uniformly improve; no corrected X1 variant is
selected as a replacement oracle or declared a passing historical certificate.

Center factors and corrected initial C change original finest C only at about
1e-15–1e-14 for the nominal/smaller center; a deliberately larger uncorrected
center produces the expected missing-center contribution (~7.8e-13). In X2,
100x center shifts, 96->192 EOS quadrature and plain/compensated accumulation
change C by at most ~6.2e-18 at the tested resolutions. Successful binary64
comparisons are ~1e-15. These observations do not account for the original
1e-12-scale shell residual. The raw-H fixed-step ablation has C offsets
3.31e-7,3.48e-8,2.83e-9 at N=512,1024,2048: onset mapping materially matters.
Its single measured order is 3.214, **not** the anticipated eventual 3/2;
that asymptotic prediction remains INCONCLUSIVE on the allowed raw-H grids.

Mechanism attribution is PASS for the target-0 support trigger and shell error
localization. Transfer of that mechanism to targets 10/20 is supported by their
nearly identical X1/X2 patterns, but remains an inference: no target-10/20
instrumented mechanism replays were run.

## Mass: discrete integration defect, not an arbitrary roundoff allowance

The passive X1 mass integral uses accepted RKF45 stage values and the radial
mass equation, without adding a controlled channel. All EIGHT target-0 baseline
replays reproduce all NINE original reported outputs bit-for-bit. Thus the
instrumentation has not changed their accepted numerical result.

At common B and the retained 2e-13 tolerance:

| Target-0 X1 variant | Reported M (km) | Passive radial Mint (km) | (Mint/M)-1 |
|---|---:|---:|---:|
| regularized | 0.92093274313331197 | 0.92093274313328632 | -2.785e-14 |
| unregularized | 0.92093274313360662 | 0.92093274313328691 | -3.473e-13 |

The two passive integrals differ by only 6.41e-16 relative. The difference of
their discrete mass-chain defects accounts for **99.8%** of the reported mass
disagreement. In the unregularized trace the defect grows mainly between muon
onset and Hb, and decreases roughly tenfold per tolerance decade. This is the
expected signature of integration truncation in the nonlinear derived mass
`m=v*u^(3/2)`, whose continuous identity need not be preserved by RKF45.
The pointwise algebraic chain identity tested in X1 could not detect it.

At identical retained Hc, X1-minus-X2 mass differences are -7.42e-14
(regularized) and +1.87e-13 (unregularized), while B differences are -1.04e-13
and -1.78e-13. Thus common-B shooting compensates differently for each B
functional's integration error; it does not remove the mass-chain defect.
Center variants alter the mass difference by only ~1e-15–2.1e-14; changing
regularized accumulation changes M by 1.6e-15. No roundoff allowance is added
to the original gate. Mass mechanism attribution is PASS at target 0, with
shared EOS/finite-arithmetic limitations explicitly remaining.

## Cross-method, anchor and paired agreement

X2 N=4096 minus X1 finest, in relative units:

| Quantity | vs regularized, range over three targets | vs unregularized |
|---|---:|---:|
| M | -8.27e-14 to -7.39e-14 | -4.21e-13 to -3.96e-13 |
| R | +1.69e-13 to +1.76e-13 | +1.95e-13 to +2.02e-13 |
| neutron edge | +1.66e-13 to +1.74e-13 | +1.90e-13 to +1.96e-13 |
| L | -8.86e-12 to -8.85e-12 | -8.00e-13 to -7.86e-13 |
| C | -1.37e-12 to -1.36e-12 | -6.35e-12 to -6.34e-12 |

These are small physical differences but exceed the new discretization tail
estimates. The predeclared cross-method uncertainty classification is
**INCONCLUSIVE**, because X1 supplies no eligible global reference uncertainty.
A small difference is not a certificate. In particular, the false JL shell term
explains most of the regularized L offset, so treating all X1 offsets as an
independent random error or an arithmetic floor would be wrong.

X2 residuals against the unchanged model anchor, in units of 1e-12:

| Quantity | Target 0 | Target 10 | Target 20 |
|---|---:|---:|---:|
| C | +0.974 | +0.678 | -4.635 |
| L | +0.942 | +0.648 | -4.152 |
| M | -0.094 | -0.121 | -0.201 |
| R | +2.810 | +3.325 | +1.788 |
| neutron edge | +0.434 | +0.331 | -1.601 |

All lie inside their original bands descriptively. C's remaining anchor residual
is far larger than the observed X2 discretization tail, and target dependent.
Center/EOS/summation ablations do not explain it. The model-anchor extrapolation,
common-B model alignment and its documented heuristic allowances remain candidates;
there is no exact physical anchor to uniquely separate them. This absolute-anchor
attribution remains INCONCLUSIVE. X1 A2 is never rescored.

X2 target-20 relative drifts are C=1.470015418e-6, L=1.809209147e-6,
M=-9.506857716e-7, R=1.283061640e-6, edge=1.257426109e-6 and
mu_inf=1.882637835e-7. Paired responses remain much closer than absolute values
because shared offsets cancel. They do not establish absolute accuracy or an
all-21-target bound. Original physical B strings and each surface cutoff are retained.

## Physical headroom versus reference certification

All 462 historical component/sample headrooms were calculated BEFORE X2 solves,
from the original thresholds and drift envelopes, not from a refitted certificate.
The table below shows the minimum over the retained 21 samples; it uses no new
21-target numerical run. These margins are conditional on the historical envelope
and its unresolved numerical-error interpretation. Thresholds themselves sometimes
contain unqualified historical errors. They are not independently certified budgets.

| Historical component | Minimum T - drift | Norm |
|---|---:|---|
| Cstar | 9.8540992e-05 | relative |
| Ltilde_Me | 9.877332e-05 | relative |
| Ltilde_Mmu | 7.4318907e-05 | relative |
| N_e | 9.5116345e-05 | relative |
| N_mu | 4.5924652e-05 | relative |
| N_n | 9.9028029e-05 | relative |
| P0_average | 9.9811736e-05 | relative |
| P1_average | 9.9738543e-05 | relative |
| P2_average | 0.0001 | exact-zero identity required |
| Tsurface_inf | 9.9212939e-05 | relative |
| Z_row_npe | 9.7836982e-05 | historical scaled profile/row norm |
| Z_row_npmu | 6.4907789e-05 | historical scaled profile/row norm |
| envelope | 9.9024632e-05 | relative |
| metric_structure | 9.6698758e-05 | historical scaled profile/row norm |
| mu_B | 9.9811736e-05 | relative |
| mu_n_profile | 9.9811735e-05 | historical scaled profile/row norm |
| radius | 9.8716936e-05 | relative |
| species_support | 1 | exact topology required |
| surface_gravity | 9.6294885e-05 | relative |
| t_e | 2.8915626e-06 | absolute tangent component |
| t_mu | 2.106701e-07 | absolute tangent component |
| t_n | 9.6256315e-05 | absolute tangent component |

For nonzero scalar ratios, if each endpoint's actual relative error is bounded
by e, the ratio error is bounded by 2e*|ratio|/(1-e). Solving for e using the
headroom allows about 4.93e-5 per endpoint for C and 4.94e-5 for L. The narrowest
tangent margin is 2.1067e-7 in **absolute t_mu units**, allowing half that per
endpoint if all other errors were absent. It is not a 1e-7 relative demand.
Profiles/matrices need their historical scaled norms; support and exact P2 cannot
be certified by scalar allowances. [headroom.json](evidence/headroom.json) gives
every component, sample, original threshold and norm distinction.

Thus three different questions must remain separate:

* Assessing the frozen-relative approximation needs justified uncertainty below
  its remaining headroom, often around 1e-5–1e-4 in relative units. It does not
  intrinsically need 1e-13 accuracy.
* Absolute stellar/thermal accuracy requires an independently stated application
  tolerance and the full physical/model error budget. None was supplied here;
  C's linear controlled-gas coefficient does not qualify production thermal knots,
  cache interpolation, realistic EOS microphysics or a trajectory.
* Certifying an independent reference requires a justified enclosure for EOS,
  inversion, center/surface, integration stability, arithmetic and the B root.
  Observed 1e-15 discretization tails alone do not supply that enclosure.

## Bounds actually established and next useful experiment

**No unconditional global numerical error bound is established.** The GL4
Richardson estimates are conditional on continuing the observed asymptotic
expansion without a hidden common bias. Center/EOS/arithmetic sensitivities are
finite comparisons, not upper bounds. Exact ratio-error propagation above is
conditional on valid endpoint bounds, which this experiment has not established.
Off-stage closures and sampled defects test consistency, not a norm-bounded global
solution error; an approximate trajectory can satisfy them while sharing a bias.

A useful next experiment is a bounded residual/defect enclosure for this now
understood GL4 path: interval EOS quadrature/inversion, a center-series remainder,
interval or adjoint stability propagation of the defect, arithmetic enclosure,
and a nonzero bounded dB/dH to propagate through the common-B solve. The
regular-singular center needs separate treatment; a sampled maximum defect and
an assumed Lipschitz constant would not suffice. The target enclosure should be
chosen from the physical headroom before tightening reference accuracy. Another
blind RKF45 tolerance ladder is not useful. A separately authorized validation of
the discovered endpoint-support correction could be useful, without modifying X1.

BA13 remains unqualified for tangents, Z, full particle numbers/support, both
Urca channels, production Cstar knots/cache, profiles/metric, envelope/surface
quantities, source averages and common-B/global uncertainty propagation. The
historical matrix's exact-P2 entries are unchanged; none of its other components
is promoted by these three stars. No fitter/covariance or envelope interpretation
was amended. No certificate, baseline, master branch, BNV trajectory or owner changed.

## Failures, resources and scientific disposition

| Item | Disposition |
|---|---|
| X1 source/results/failure authentication | PASS |
| X2 formulation checks and three-target GL4 empirical order | PASS |
| Target-0 C support-trigger/shell localization | PASS, with stated scope |
| Target-0 mass-chain attribution | PASS, 99.8% of observed disagreement |
| Raw-H ablation's eventual 3/2 order | INCONCLUSIVE at prescribed resolutions |
| First binary64 N=512 diagnostic | FAIL: stage residual |
| Corrected binary64 N=512 diagnostic | FAIL: Newton convergence |
| Corrected binary64 N=1024 and 2048 diagnostics | PASS; no double ladder qualification |
| X1 historical C contraction and mass gates | FAIL, preserved |
| Cross-method reference qualification / rigorous global bound / BA13 | INCONCLUSIVE |

The initial binary64 loop stopped at 32 epsilon, looser than its unchanged 1e-15
stage gate. Its separate corrected copy uses a tighter stop; N=512 still fails
Newton convergence and is retained. No primary extended-precision solve failed.
See [DIAGNOSTIC_CORRECTION.md](DIAGNOSTIC_CORRECTION.md).

All substantial calculations ran on Slurm, one pinned logical CPU (affinity [4]),
1 GiB memory. Slurm allocated two hardware threads despite cpus-per-task=1;
only one was used for numerics. Jobs 17024–17031 consumed 148.018673 measured
numerical/controller-loop CPU seconds. Independent audit 17034 PASS consumed
0.180683 CPU seconds; combined measured 148.199356. Failed audit launchers 17032
(missing Git PATH) and 17033 (missing /usr/bin/time, before Python) each lasted
one second; charging their entire two-thread allocations adds four CPU-seconds.
Measured numerical cost is safely below 900. Interpreter startup/source hashing
outside the controller timers is not included in that precise loop total; no
claim of exact whole-allocation CPU measurement is made. Compiler costs are
separate build receipts. All decisions had ample remaining budget; no expansion
to 21 targets or new resolution occurred.

Final controller records for 17024–17031 expired before collection. Their saved
allocations, per-command CPU/return-code/hash receipts, raw outputs and stderr are
preserved; no missing Slurm final state is fabricated. Final failed audit states
and successful 17034 state are retained. The independent audit verified eight
bit-identical X1 replays, 24 GL4 order tests, 24 zero-shell-JL tests, 462 headrooms,
all execution gates and all authenticated external artifacts.

Changed files are ONLY the new X2 test solver, X2 tools and X2 documentation/evidence.
The source-bound external manifest contains raw traces, diagnostic generated
sources/executables and receipts. No large profiles were needed. Small exact
outputs and derived summaries are committed. GitHub CLI/credentials are absent;
a verified incremental bundle against 6cc7e20 is supplied for the office Mac.
See [REPRODUCE.md](REPRODUCE.md) for reproduction and transfer.
