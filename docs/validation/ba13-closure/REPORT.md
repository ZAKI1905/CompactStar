# BA13 frozen-background closure assessment — disposition E

**Full BA13 closure was not achieved.** The corrected 21x22 matrix has
**21 PASS, 0 FAIL, 441 UNQUALIFIED**. All PASS rows are the exact historical
P2 identity. No approximate muon-sector row is qualified; the four numerical
groups have 84 UNQUALIFIED rows, and species support has 21 more. No demonstrated
threshold violation requires a production-method change. This is an executed
UNVERIFIED SCIENTIFIC CANDIDATE, not owner ratification or
CANDIDATE_READY_FOR_INDEPENDENT_REVIEW for a complete construction.
Sources: [summary.json](evidence/summary.json), [complete matrix](MATRIX.md),
[machine dispositions](evidence/dispositions.json).

The smallest remaining scientific question is whether a quantity-specific,
accepted continuation construction can supply uniform bounds/validated
estimates for aligned structural errors and curvature on the tiny fixed-topology
H interval. Three targets and small numerical terms alone do not supply them.
The tangent reference also needs offset-H endpoint refinement or a justified
correlated derivative-error estimate. These are bounded questions; another
general stellar-solver investigation is unwarranted.

## Authentication, scope and governance

Entry is 286aafb63002093b3052e23a3972903011f36035 on the clean
candidate/linux-x2-fixed-step-20261016 worktree. The isolated descendant is
candidate/linux-ba13-closure-20261019 under worktrees/ba13-closure. Predeclaration
fdffaaf was committed before new scored closure calculations. All local commits
touching the changed validation paths were ancestors of entry; no competing
owner was found. The original independent review is unmodified: all 40 manifest
entries authenticate. Its copied, unchanged reproduction script passes all 31
checks; 391 headroom, 257 X2 and 236 nominal-input files authenticate. Both
historical certificate files reproduce byte for byte. See
[PREDECLARATION.md](PREDECLARATION.md), [authentication](evidence/authentication.json)
and [audit](evidence/audit.json).

The certificate SHA-256 remains
4212cad05bb4244ab75e80b8dfea023e6d0fcee69086c9993bd48d26a40e41fb.
The fitter SHA-256 remains
9b3ef1ce17df7e1f8d2d7c40e111ee0844abfddc2cc54b8f93cb3f2ce1322370.
GOVERNANCE §§1–5, INV-02/03/04/06/14/15 and accepted Phase-6A-1 preflight §11
govern this validation numerical-method/engineering/documentation/artifact change.
ADR-0021 Option A and interpretation (2) remain PROPOSED. Production physical
owners, dependency locks, accepted baselines, historical source/certificate
identities, canonical master and X1/X2/P5/P7/BA13 failures are preserved. The
changed-file list and exact hashes are in [source-hashes.json](evidence/source-hashes.json).

## Confirmed defect corrections

| Review item | Implemented result |
|---|---|
| D1 wrong reference | Explicit P-R, aligned R-X, common-B and X-Q decomposition. All production structural magnitudes are separate from X2 reference error; historical scales and units survive. |
| D2 defective estimator | Both signed first-order and absolute transports retired from closure; requesting that mode refuses. Historical code and failed records remain. No fitted factor or rigorous-bound claim is introduced. |
| D3 evidence acceptance | Every conditional term requires a named external acceptance; numerical estimates require an accepted standard. Missing flags are not the sole guard. Forged authority, missing terms, nonfinite values and threshold mutations refuse. |
| D4 tangent common-B | All three tangent components carry measured B-functional effects in absolute units: maxima t_n=1.0925e-12, t_e=6.6358e-13, t_mu=4.2896e-13. Functional/sensitivity remainders remain explicit. |
| D5 interval interpretation | Physical response and intermediate production output are distinct. Sampled selection identity is confirmed; monotonicity is not extrapolated without a derivative condition. |
| D6 bookkeeping | Root propagation is wired for the X2 reference; no production binary64 arithmetic allowance uses 2^-64. Missing arithmetic/stage remainders remain unresolved. Definitional zeros have attribution. |

These changes are in [enclosure.py](../../../tools/linux/ba13_headroom/enclosure.py),
[the corrected assessor](../../../tools/linux/ba13_closure/assess.py) and
[x2_uncertainty.cpp](../../../tests/bnv/x2_uncertainty.cpp). The old assessment
entry delegates to the isolated closure output; it cannot overwrite historical
headroom artifacts. Equations and input hypotheses are in [METHOD.md](METHOD.md).

## Muon-sector calculations

All six prescribed N=1024/2048 fixed-Hc functional evaluations and six third-width
N=512 stencil evaluations completed, without root shooting or production solves.
The benchmark plus remaining numerical evaluations used 13.526425 CPU seconds;
the earlier schema-gate stop is charged separately. Third-width t_mu increment
ratios are approximately 4.0000 at all three targets. All assessed functional
fourth-order tests meet the frozen order criterion.
Sources: [m1.json](evidence/m1.json), [resources.json](evidence/resources.json).

At target 20, in the ORIGINAL historical scale:

| Quantity | Measured paired production-X2 difference | Partial reference estimate | M1 result |
|---|---:|---:|---|
| t_mu, absolute | 2.22591e-10 | 4.17816e-9 | INCONCLUSIVE: reference estimate exceeds discrepancy |
| N_mu, relative | 1.15734e-10 | 2.16004e-15 | Estimator test passes; incomplete reference band |
| Z_row_npmu, historical row norm | 8.28165e-11 | 3.14897e-15 | Estimator test passes; incomplete reference band |
| Ltilde_Mmu, relative | 7.19856e-10 | 5.46604e-15 | Estimator test passes; incomplete reference band |

These are measured differences and conditional estimates, not qualified totals.
Some differ slightly from the review's diagnostics because this calculation
uses the same original scale for both references and carries denominator error
separately, whereas the review's diagnostic ratios used each reference's own B0
quantity. No norm or threshold has been changed. The independently added tangent
endpoint denominator errors are conservative but conditional; center-refining
to N=4096 would not establish offset-H derivative correlation. The permitted
extra resolution was therefore not spent.

M5 independently reproduces the full-endpoint inverse majorant:

| Row | Conditional majorant | Minimum headroom | Fraction of headroom |
|---|---:|---:|---:|
| Z_row_npe | 1.96465176e-7 | 9.7836982e-5 | 0.2008% |
| Z_row_npmu | 5.59672494e-7 | 6.4907789e-5 | 0.8623% |

Each endpoint has its own G and full background+center+actual-cut tail error set;
both nonlinear Neumann remainders and the full Z matrix are retained. No perfect
background correlation is assumed. The remaining input is the per-target
background magnitude hypothesis. This is an ALTERNATIVE structural-error route,
not an extra contribution added on top of the direct X2 comparison.
See [m5.json](evidence/m5.json) and [independent audit](evidence/audit.json).

## M2–M4 and continuous coverage

All 21 profiles have 20259 radial nodes. Their first 20258 nodes are identical;
only the finite-pressure surface node differs. Thus the review's literal
complete-grid identity prerequisite is false, while its fixed interior support
grid premise is supported. This is an INCONCLUSIVE transfer gate, not a physical
threshold FAIL. A prospective construction should isolate the moving terminal
event; grid identity still cannot prove a uniform structural-error magnitude.
No three-target magnitude was promoted across the other 18 samples.
Source: [m2.json](evidence/m2.json).

Neutron support cell k=20118, its radial endpoints, and envelope node 20105 agree
across all 21 profiles. Sampled species topologies agree. The neutron edge stays
at least 3.90856e-4 km from a cell boundary while moving 1.58981e-5 km. The
envelope margin is 1.68429e-6 in nu against sampled motion 1.07126e-7. Endpoint
threshold ordering has Hc-Hmu≈0.01972 and Hn-Hs≈0.000820. These are strong
diagnostics, not a proof of between-sample monotonicity.

A sufficient Lipschitz condition would allow edge motion up to 15634 km per
unit x and envelope displacement up to 67.37 per unit x; observed maximum
secants are 15.90 and 0.1071. The Lipschitz suprema are not established. The
physical-response construction carries no fictitious node-jump allowance;
production-output certification additionally requires adaptive/event behavior.
Source: [m3.json](evidence/m3.json).

All 20 intervals were evaluated. K=10 max_three |Q_xx| gives maximum conditional
interval remainders 9.08e-15 (absolute t_mu), 6.37e-12 (N_mu), 5.70e-12 (muon Z
row), and 6.58e-13 (Ltilde_Mmu). Narrower-width checks agree with the wider-width
prediction; the sampled curvature ratios across the three targets are below
1.00007. None establishes the required global supremum or moving-cut derivative
remainder. Those terms remain conditional. See [m3.json](evidence/m3.json).

M4 propagates the measured B functional, includes the common B0 denominator in
the abscissa correction, evaluates all original outer/inner fit residuals and
perturbs the WLS design and value-dependent weight floor. Same-source signed
corrections use (I-H)e; no fit covariance is counted as physical uncertainty.
Unbounded structural/reference residual-input and design remainders remain
explicit. See [m4.json](evidence/m4.json).

## Validation, stops and remaining work

All 17 existing BA13 tests and 13 new controls pass. The independent audit passes
112 checks, including Decimal N_mu/t_mu calculations, all 20 changing-Jacobian
Neumann-series comparisons, exact inventory/threshold preservation and all 884
input-manifest entries. This verifies the implementation and honest refusals;
it is not scientific acceptance. See [audit.json](evidence/audit.json) and the
three committed test receipts. The wider production/trajectory suite was not
run: no production owner changed and trajectories are outside this task.

The first replay gate stopped because an older schema omitted four zero-valued
transport fields. Every retained numerical decimal matched. The corrected gate
explicitly checks those original fields and the additional zeros, then passes.
The stopped run, traceback, receipts and Slurm FAILED state are retained. There
was no scientific replay discrepancy and no change to a numerical tolerance.
Source: [resources.json](evidence/resources.json) and external manifest.

M1 t_mu and M2 transfer remain INCONCLUSIVE; M3 is conditional. Under the
predeclared stop rules, no qualified full-inventory extension or nodewise observer
was launched. All 22 groups still appear in the corrected matrix, using retained
profiles and passive functional evidence where available. Metric_structure and
mu_n_profile have no independent nodewise reference. A small passive step
observer is technically possible: record the existing collocation polynomial,
evaluate it on the original support coordinates, and verify dense-output order
independently. Quadratic interior dense output cannot inherit endpoint GL4 order
without a test. No physical solver change is needed, but this would not resolve
the present transfer/continuous-window gate.

Minimal next work, after choosing an evidence standard:

1. Establish a uniform aligned structural-error and curvature continuation
   condition over Hc in [0.13154851072782064,0.13154890784581055], treating the
   moving final node separately. A named conditional owner acceptance is a
   different scope from a validated enclosure.
2. Resolve the tangent reference with twelve offset-H N=1024/2048 evaluations
   at the already declared narrow width, or a justified correlated derivative
   estimate. Based on the benchmark, the twelve evaluations cost about 20 CPU
   seconds. This is a future calculation, not executed after the transfer stop.
3. Supply quantity-specific EOS/center/stage/arithmetic and fit-input remainder
   evidence under that standard. Their measured sensitivities are already saved;
   no arbitrary floor or tolerance multiplier may fill the gap.
4. Once M1–M6 genuinely pass, test a passive node observer at the same three
   targets and N=512/1024/2048. Solver cost extrapolates to about 12 CPU seconds,
   plus observer/postprocessing cost that must first be benchmarked. This is
   specifically for the two named profile groups, not a new 21-star campaign.

The five owner decisions are in [OWNER_DECISIONS.md](OWNER_DECISIONS.md).
Full execution, numerical CPU/RSS/wall receipts, exact hashes, manifests and
office-Mac bundle instructions are in [REPRODUCE.md](REPRODUCE.md).
