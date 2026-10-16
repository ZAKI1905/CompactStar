# BA13 headroom-scaled uncertainty assessment

**Recommendation B: identified numerical and interval-coverage gaps prevent
closure.** This is an executed candidate assessment with a complete 21-by-22
inventory, not a completed qualified uncertainty construction. The matrix has
**21 exact-P2 PASS, 441 UNQUALIFIED, zero demonstrated threshold FAIL**. Missing
terms remain null or are explicitly named missing remainders of partial estimates.
No incomplete total is called a pass. See [matrix](MATRIX.md),
[evidence/summary.json](evidence/summary.json:1) and the complete numerical rows
in `evidence/matrix.json.gz` / `evidence/matrix.tsv.gz`.

The observed errors do **not** demonstrate that production changes are necessary
to satisfy the frozen-relative thresholds. Nor do the results support A or a
complete conditional claim D: some quantity-specific propagation remains absent,
not merely awaiting acceptance of one interpretation. The next action is **one
independent review of this assessment, propagation code, failure evidence and
specified residual gaps**, followed by owner adjudication. Another general
oracle comparison is not recommended.

## Identity, contract and preservation

Entry: `89b629dd9961c29e127badbd723455d00dc38553`, clean branch
`candidate/linux-x2-fixed-step-20261016`, worktree
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/x2-fixed-step`.
All other local branch tips were ancestors at entry; there was no competing file
owner. Authentication is in [authentication.json](authentication.json:1).
Changes are numerical-method validation, engineering, documentation and generated
artifacts. No production owner, fitter, baseline or certificate is modified.
Optional passive observers extend X2; the final null-observer executable exactly
reproduces the entire retained JSON, including decimal strings and shooting
statistics, at targets 0,10,20. See [audit](evidence/audit.json:1).

The acceptance question follows the accepted Phase-6A-1 preflight §11
(`docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:655`),
GOVERNANCE §§2–5 and INV-02/03/04/06/14/15. It requires full drift, numerical error
and nonlinearity within the threshold. It neither permits replacing that requirement
with small empirical residuals nor demands a formal computer-assisted proof of the
entire stellar model. ADR-0021 Option A and interpretation (2) are the **proposed**
construction evaluated; owner ratification is not implied.

The 462 threshold/drift pairs are exactly X2's frozen headroom input, authenticated
against the original retained native nominal certificate `run1/fit/certificate.json`
(SHA-256 `4212cad05bb4244ab75e80b8dfea023e6d0fcee69086c9993bd48d26a40e41fb`).
This identity is explicit: it is not the Mac replay certificate, a refined star
sequence, or an X2 substitution. The native historical fitter SHA-256 remains
`9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370`.
The audit verifies 236 original nominal manifest files, 257 original X2 external
artifacts and 387 protected evidence/review files. Historical X1, P5, P7 and BA13
scientific dispositions are unchanged ([audit](evidence/audit.json:1)).

X2 establishes independent formulation checks and empirical fourth-order
convergence at three targets. Its Richardson tails are conditional estimates,
not unconditional global bounds. The support and mass mechanisms are directly
demonstrated at target 0; transfer to 10/20 remains inference. See the preserved
`docs/validation/x2-fixed-step/REPORT.md:1` and `MATHEMATICS.md:1`.

## What was implemented and measured

[Method](METHOD.md:1) allocates every approximate row's headroom across twelve
contributions; fractions sum to one. The original units survive: scalar relative,
absolute tangent, historical scaled maximum-component Z/profile norms, and exact
identities. The implementation does not replace the historical Z row scale with
a row-sum norm. Allocation is an accuracy target, never an error estimate.

The code in `tools/linux/ba13_headroom/enclosure.py:1` implements deterministic
paired errors, uncertain normalization, functional-error-aware root propagation,
entrywise matrix inverse majorants with nonlinear remainder, the changing paired
inverse Jacobian, fixed-design residual perturbations `(I-H)e`, and a conditional
inter-sample curvature bound. The test suite includes all six requested faults,
a between-knot bump invisible at the samples, and refusal of a partial estimate
with an unresolved remainder. **14 tests pass**; these validate propagation rules
and counterexamples, not their unestablished physical inputs
([test output](evidence/tests.txt:1)).

`tests/bnv/x2_uncertainty.cpp:1` reuses the existing GL4 stages and transformed
center/onset equations. It adds passive integrals for all four particle numbers,
both enabled Urca channels, all independent G components, and direct P0/P1 averages.
Z is formed with the existing chemical map. No extra state drives structure or
changes accepted steps. N=128/512 compare the new functionals; retained N=4096
values validate the already available structural quantities. Fine solutions are
never defined as exact truth. Additional fixed-H offsets have one purpose:
number tangents and dQ/dB, using their actual computed B differences. Both widths
and center/EOS perturbations of those stencils are retained
([pilots](evidence/pilots.json:1)).

The two derivative widths give dB/dH around `1.91804e57` counts, with relative
width sensitivity about `3.3e-9`. These are positive local numerical secants,
not a validated lower derivative bound over an uncertain root enclosure. The
nominal-profile independent B reconstruction differs from the production B
functional by about `-7.50e-9` relative. Its paired change from target 0 to 20
is about `-2.33e-12`. The computed shooting residual cannot replace this functional
uncertainty. The matrix carries measured paired alignment effects separately
from residual, functional uncertainty and derivative-enclosure gaps
(`external profile-v2-{0,10,20}.json`, bound by [manifest](evidence/external-manifest.json:1)).

### The defect propagation test gives a specific failure

The observer computes the implicit-stage tangent transition matrix and transports
an off-stage collocation defect with four/eight-node comparison; it assumes no
arbitrary global stability constant. The center remains a separate analytic
startup. The intra-step transport is a first-order Jacobian approximation, and
neither its remainder nor a nonlinear derivative tube is enclosed. At target 0,
the N=128 signed radius correction is `-5.94e-11`, while retained fine-minus-coarse
radius is `+1.20e-10`. The signed correction therefore fails even as a directional
predictor. It is not used to assert an error bound. The mass prediction is also
inaccurate in magnitude ([pilots](evidence/pilots.json:1)).

The absolute numerical transport at N=512 is approximately `8.07e-7` relative
for R and `5.54e-6` for M at all three targets, much larger than Richardson tails
but relevant to the headroom. It encloses these finite observed differences;
that does not validate the uncomputed transport/quadrature/nonlinear remainder.
Only R has a directly evaluated stage-error transport in the matrix; other
functional stage errors are explicitly missing rather than borrowed from R.
This is a **specific missing estimator remainder and quantity sensitivity**, not
a demand for proof of the entire stellar model ([pilots](evidence/pilots.json:1)).

EOS96/192 and center-factor perturbations are recorded, including derivatives
of their effects on tangents. Some differences disappear at the binary64 summary
precision; exact 21-digit solver output is retained externally. A zero measured
sensitivity is not entered as a zero upper bound. EOS quadrature/inversion and
center-series/initial-integral remainders remain named gaps. Floating-point
operation-count allowances cover accumulation only; input/elementary-function
propagation remains explicit. The matrix does not use these missing terms to
claim a completed numerical estimate (`matrix.json.gz`).

## Historical production, chemical and sample coverage

All 21 **existing nominal profiles** were read and authenticated. No blanket
21-star solve was run. Three-target profile checks preceded application to the
remaining 18 profiles through a Slurm array. These computations provide:

* Signed neutron-edge and production spatial-C discrepancies; independent
  density reconstructions for all particle numbers; both Urca channels.
* Full independent fixed-profile G, reconstructed-background G, Z, current
  center leading terms, and the retained per-target actual-cut tail certificates.
* All thermal knots and deterministic cache factors, physical versus selected
  envelope-base nodes, direct P0/P1 averages, chemical-potential profiles and
  remapped metric profiles. Original production fields remain inputs, untouched.

Sources: `profile.py:1`, `supplement.py:1`, and all 21 profile/supplement entries
in [external-manifest.json](evidence/external-manifest.json:1).

The cache law reproduces production to below `7e-15` relative at the pilots.
Its shared factor can cancel at equal Tinf knots for this linear entropy law and
unchanged cache grid. Its absolute bias (previously measured up to `6.56e-4`)
is retained as an absolute-C issue, not added independently at both endpoints
or cancelled from unrelated quantities. Production L is about `9.29e-5` high
at target 0 against the onset-aware same-profile reconstruction, but the target-20
paired correction is only about `3.04e-6`. A simple whole-cell support allowance
is `1.68e-4` relative for L at a single endpoint, already larger than its frozen
threshold. It cannot serve as a successful independent-endpoint budget.
This is failure of that loose allowance, **not** proof that production must change
(`profile-v2-0.json`, `matrix.json.gz`; earlier cache law:
`docs/validation/ba13-uncertainty/SCIENTIFIC_REVIEW.md:1`).

The inverse calculation retains the inherited background as one shared matrix
perturbation, separate current center/tail terms, distinct inverse Jacobians and
both Neumann remainders. The nonlinear remainder is nonzero (target-0 absolute
Z entries about `4.5e-67` to `8.4e-66`). This completes the propagation formula,
not validation of the inherited B0 background at depleted targets. The change
of that background error still needs an estimate. No common G cancellation is
assumed after nonlinear inversion (`profile-v2-0.json`, `enclosure.py:1`).

The finite-pressure surface is the defined stellar quantity under INV-06 and
ADR-0009, not a vacuum target. The accepted chemical-tail term is retained in the
chemical inverse error set only. It is not double counted as a radius error.
Event-location and production-node errors remain separate
(`docs/SCIENTIFIC_INVARIANTS.md:406`, `matrix.json.gz`).

## Frozen envelope and residual gaps

The unchanged outer-half fit is evaluated at all 21 samples. Extra inner-half
residuals and the residual of the same signed input correction under `(I-H)` are
reported separately. Historical covariance is retained historically and is never
added a second time as physical uncertainty. Uncertain achieved-B design and
unseen nonlinear residuals remain missing. The complete interval includes all
20 gaps between achieved sample coordinates. The reported sampled curvature is
only a diagnostic; the separate sufficient-curvature column derives from the
allocated budget and is **not** an asserted bound
([inner-half.json](evidence/inner-half.json:1), `matrix.json.gz`).

For example, the t_mu inter-sample share would allow a normalized second-derivative
bound around `6.74e7` in its absolute tangent norm per squared depletion coordinate;
its sampled finite-difference curvature is about `351`. That large separation
suggests practical numerical coverage may be achievable, but cannot exclude a
between-sample feature by itself. Production support/envelope node selection also
needs a jump allowance. Simply declaring all production quantities C2 would be
false. The deliberate bump test verifies this failure mode
(`matrix.json.gz`, [tests](evidence/tests.txt:1)).

Every UNQUALIFIED row lists its specific missing estimates. The important shared
and family-specific gaps are:

| Scope | Remaining estimate, not a new scientific convention |
|---|---|
| All approximate groups | EOS/inversion and center remainder; quantity-dependent integration/stage/arithmetic propagation; common-B functional/sensitivity uncertainty; shared B0 normalization; uncertain-design/nonlinear fit remainder; inter-sample coverage |
| t_n, t_e, t_mu | Remainder beyond measured center/EOS derivative perturbations, stencil remainder and achieved-B coordinate propagation; preserve absolute units |
| Both Z rows | Depleted-target variation of inherited background error, in addition to the now implemented changing-Jacobian and inverse-remainder propagation |
| Metric profiles | Full nodewise structural-error transport; the computed support remapping is only one production contribution |
| P0/P1 and chemical-potential profiles | Remaining event-measure/background and nodewise transport after fixed-profile integration/identity checks |
| Exact support | Exact topology passes at all sampled profiles; continuous no-crossing coverage is not established |

At target 20, incomplete partial utilizations are roughly 0.043 for L_Me,
0.015 for C, 0.351 for the muon Z row, 0.378 for t_mu, and 0.541 for N_mu.
These figures include only recorded finite terms and are explicitly **not total
qualified utilizations**. Missing terms cannot be set to zero. No evidence here
establishes recommendation C (`matrix.json.gz`).

## Validation, resources and delivery

The audit passes the 462-row identity/budget checks, three exact decimal X2 replays,
all nominal/profile/tail checks, and 14 fault/algebra tests. The first fixed-H
replay comparison failed because shooting-history diagnostic maxima differed;
replaying the original shooting history passes. Three profile jobs failed after
calculation when NumPy int64 failed JSON serialization; their partial outputs and
errors remain, and corrected outputs validate serialization before publication.
The signed-defect failure remains a scientific limitation
([audit](evidence/audit.json:1), [resources](evidence/resources.json:1)).

Retained numerical measurements plus a conservative charge for the three failed
profile allocations total **237.21 CPU seconds**; builds total **9.99 CPU seconds**.
Peak recorded RSS is **149476 KiB**. This excludes earlier short local assembly/
audit invocations, hashing and interpreter startup outside timers, so it is not
an exact whole-session CPU claim. Slurm final records for early jobs expired;
per-command receipts and allocations survive. Failed-profile 13/12/12-second
runtime observations are explicitly transcribed, not represented as saved raw
scheduler receipts. The 1800-second numerical planning ceiling was not approached
([resources](evidence/resources.json:1)).

Large outputs stay in the isolated external directory with authenticated manifests.
Compact compressed matrices and small summaries are committed. GitHub CLI is absent;
[REPRODUCE.md](REPRODUCE.md:1) describes the incremental bundle and existing
Mac workflow. No master merge, dependency change, campaign, trajectory, physical
rate selection, certificate promotion or ADR ratification occurred.
