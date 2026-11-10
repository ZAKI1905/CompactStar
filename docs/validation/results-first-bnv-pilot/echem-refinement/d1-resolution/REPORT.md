# D1 reference convergence resolved — conditional local qualification

2026-11-09. **PASS_CONDITIONAL for the existing 13 local chemical
reconstructions.** The knot-aligned reference passes unchanged C1, C2 and C3
at every point. It also passes the predeclared stronger state/energy comparisons
using a conservative bound for arbitrary time-dependent perturbations within
the assumed thermal band. The original global-grid **D1 FAIL remains FAIL** and
is preserved byte-for-byte at [the original receipt](../d1/REPORT.md#L1).

This is a convergence-supported local qualification under the explicitly stated
thermal and discretization assumptions below, not an unconditional IVP enclosure
or domain-owner ratification. The evidence is sufficient for the review's local
accuracy requirements under those assumptions; no further reference precision is
needed for this task. [Machine-readable receipt](evidence/summary.json#L1).

## What failed in the original reference

For the central thermal path, the signed adjacent differences are below. All
numbers are MeV; the columns mean eta_fine minus eta_coarse. Their absolute
values are obtained directly by removing the sign. All three offsets, the exact
110-digit differences and each absolute difference divided by its accepted
chemical-state radius are in [original-differences.tsv](evidence/original-differences.tsv#L1).
[Generated tables](evidence/TABLES.md#L1) give the normalized comparisons directly.

| Point | Component | 64→128 | 128→256 | 256→512 | Original radius |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | e | +1.003645971445e-47 | +5.018247885140e-48 | +2.509125069315e-48 | 2.177742570168e-26 |
| 1 | mu | +3.166027501807e-47 | +1.583021899246e-47 | +7.915114588944e-48 | 4.040191348982e-26 |
| 64 | e | +1.992748735308e-42 | +1.009859717149e-42 | +5.057726929100e-43 | 7.637549006414e-25 |
| 64 | mu | +6.269928859276e-42 | +3.184618005966e-42 | +1.595412351471e-42 | 1.958952551619e-24 |
| 1440 | e | -1.245167598036e-37 | -9.296011462462e-40 | +3.368247284566e-39 | 1.720438397395e-23 |
| 1440 | mu | -4.676207499014e-37 | -7.609301076985e-39 | +1.033291333377e-38 | 4.316781172688e-23 |

The absolute differences / original radius are respectively:

| Point | Component | 64→128 | 128→256 | 256→512 |
| ---: | --- | ---: | ---: | ---: |
| 1 | e | 4.608653e-22 | 2.304335e-22 | 1.152168e-22 |
| 1 | mu | 7.836331e-22 | 3.918185e-22 | 1.959094e-22 |
| 64 | e | 2.609147e-18 | 1.322230e-18 | 6.622186e-19 |
| 64 | mu | 3.200654e-18 | 1.625674e-18 | 8.144211e-19 |
| 1440 | e | 7.237502e-15 | 5.403281e-17 | 1.957784e-16 |
| 1440 | mu | 1.083263e-14 | 1.762726e-16 | 2.393662e-16 |

The original ratios identify a reference-grid defect, not a demonstrated loss
of physical accuracy. These differences are about 1e75–1e83 times a 110-digit
endpoint spacing, so they are not near Decimal arithmetic cancellation.
Subtraction still resolves many digits. At point 1440 the small middle
difference is cancellation **between discretization-error contributions**.

Each failed path has eight nearly uniform Hermite pieces and a ninth, very short
terminal piece. Their durations are 3.637978807091713e-12, 2.328306436538696e-10
and 3.725290298461914e-9 seconds at points 1, 64 and 1440. None of the original
global-grid endpoints lands exactly on an interior retained knot. Their offsets
from nominal endpoints vary in sign and normalized size with refinement; thus
stages straddle different pieces. No interior RK4 stage resolves the short
terminal piece: it is sampled only at the final k4 stage. Exact knot offsets,
step phases and segment widths are preserved in
[original-failure-analysis.json](evidence/job-19185/diagnostic/original-failure-analysis.json#L1).

There is a quantitative explanation, beyond observing short segments. Fit the
signed differences to d_k=A/2^k+B/16^k, k=0,1,2. The first term represents the
unresolved endpoint contribution and the second fourth-order error. The fitted
first two differences predict the third to relative residuals below 3e-16,
3e-11 and 4e-7 at the three points. Independently extrapolating the preceding
regular Hermite polynomial to the endpoint predicts the first-order coefficient
as A=-H*delta_f/(64*12), from the final RK4 weight h/6. It agrees with the fit to
the same orders. A is positive and B negative. At points 1 and 64 the first term
dominates, producing ratios slightly above 0.5; at 1440 the contributions cancel
and change sign. This is a diagnostic explanation, not an error enclosure.
Point 1919 also has a short terminal piece but passed original C1; the mere
presence of such a piece is not sufficient to cause a failed ratio.

## Fixed resolution and observed convergence

The [predeclaration](PREDECLARATION.md#L1), implementation and six additional
analytic tests were committed and published as `7842938` **before** compute job
19185. The original source-bound RHS, physical constants, coefficients, exact
binary64 initial states and complete piecewise-Hermite thermal path are
unchanged. Every original input was reauthenticated. All 14,796 RHS rows / 29,592
chemical components pass the unchanged 2-ulp fidelity test; maxima remain
0.7404581109 and 0.8022895880 ulp. An independent integer-ratio conversion check
passes for both left chemical components at every point. New and old left states
and Hermite nodes also compare exactly. No model-definition discrepancy was found.

RK4 now advances separately on each Hermite piece with normalized local
coordinate s in [0,1]. It carries eta continuously through knots and uses
K=8,16,32,64 steps **per piece**, including the tiny terminal pieces. Nominal
labels N=64,128,256,512 mean 8K; the failed points have actual totals
72,144,288,576. There is no step across a knot and no thermal-node deletion or
replacement. All three offsets and all 13 intervals were calculated once at
110 digits: 156 integrations. This is an independent local reference grid,
not a production trajectory or a change to the thermal solver.

Every component and offset passes unchanged C1. The largest ratio anywhere is
**0.06250000005471215**, consistent with fourth-order contraction, versus the
unchanged limit 0.5. No floor-based exception or relaxed threshold was needed.
The finest central differences at the originally failed points are:

| Point | e: abs(eta_K64-eta_K32), MeV | mu: abs(eta_K64-eta_K32), MeV |
| ---: | ---: | ---: |
| 1 | 1.609635343806e-55 | 7.275305613169e-55 |
| 64 | 1.204049044049e-46 | 4.433354957500e-46 |
| 1440 | 5.475782033242e-40 | 2.019651734525e-39 |

The largest change from the original to the aligned reference, normalized by
the corresponding state radius, is only **8.8502921542e-16** (point 1919, mu).
That comparison supports the scientific irrelevance of the original failed
ratios; acceptance rests on the new convergence and uncertainty tests, not on
small disagreement alone. Full results and ratios are in
[results.json](evidence/job-19185/diagnostic/results.json#L1), with per-point JSON
files alongside it.

## Accuracy budgets and supported uncertainty

The original review prescription is unchanged: u_disc=2*abs(K64-K32) and
u_offset=max_b abs(eta_K64(b)-eta_K64(0)). It yields maximum C2 ratio
**0.001370086143058** (point 16, mu) and C3 ratio **0.000498235835324** (point 4).
All original allowances, including U<=0.20F and all recorded contributions,
remain intact. The largest u_disc is 3.9406946686e-38 MeV, at most
6.7176173879e-16 of the corresponding state radius. The largest u_offset is
2.6757870102e-26 MeV.

The additional sensitivity test bounds the exact model's response to any
time-dependent |delta_x(t)|<=1e-12, rather than identifying constant offsets as
extrema. The positive-polynomial derivative and first-exit/Gronwall derivation
is in the predeclaration. Every chemical-domain check passes; the largest
first-exit upper bound / domain radius is 0.500025836 and largest L*H is
2.583562e-5. Conservative Decimal rounding is used for that algebraic sensitivity
bound. Its largest thermal contribution is **8.5633122384e-26 MeV**. This does
not prove that the actual thermal solution lies in the assumed band.

For the stronger check,
u_supported_i=u_disc_i+max(u_offset_i,B_thermal). All components and points still
pass. The worst supported C2 ratio is **0.004303188584056** (point 1919, e);
the worst supported C3 ratio is **0.001660771997586** (point 1919). Thus the
combined discrepancy and supported reference allowance use at most 0.431% of
the original chemical-state budget and 0.167% of the original energy budget.

At that most limiting point, the distinct quantities are:

| Quantity, point 1919 | Actual disagreement | Supported reference allowance | Original accepted budget |
| --- | ---: | ---: | ---: |
| eta_e, MeV | 1.17913911e-26 | 8.56331224e-26 | 2.26400753e-23 |
| eta_mu, MeV | 4.63410489e-26 | 8.56331224e-26 | 5.86620887e-23 |
| Echem, MeV | 9.77980205e16 | 2.23148289e18 | 1.40252901e21 |

The energy allowance uses exactly four-corner propagation of the corresponding
chemical reference allowance. The largest O2 discrepancy in ulps across all
points is 1.9680116806 (point 1, e). Actual errors, discretization allowances,
offset probes, uniform thermal bounds, original budgets, and both versions of
C2/C3 are reported separately for **every** component and point in
[chemical-states.tsv](evidence/chemical-states.tsv#L1),
[chemical-energy.tsv](evidence/chemical-energy.tsv#L1) and
[the generated tables](evidence/TABLES.md#L1).

C4 remains **report-only**. All seven states at every point are reported against
the new reference in MeV and ulp, with ratios to the original enclosure_ctx plus
the review's u_ref. The maximum ratio is 0.9865931756011 (point 128, RKF45/4, e).
No audit enclosure is falsified by these comparisons. A ratio above one would
invalidate that audit enclosure as a bound without independently failing D1;
ratios below one do not prove the audit enclosure.

Two conditions are material to the qualification. First, sampled agreement
(maximum rho=3.1101367210e-20, fixed limit 1e-13) does not prove a uniform
continuous thermal enclosure with delta_x=1e-12. Second, twice the finest
refinement difference is a **convergence-supported estimate**, conditional on
continuation of the contracting regime beyond the finest tested grid. If
subsequent differences contract by <=0.5, their remaining geometric tail is at
most the last difference and is covered by the factor two. Four tested levels
do not themselves prove that premise. This task has not constructed a rigorous
directed-rounding enclosure of the full IVP. It has established the review's
conditional local reference accuracy with ample budget margin.

## Authentication, execution and preservation

The actual worktree is
`/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/results-first-bnv-pilot`,
branch `candidate/results-first-bnv-pilot-20261108`, starting at the requested
`871dbd77230f9e21b54ea0129473fe658236fd56`. It was clean before agent-created
bytecode; that bytecode was removed. All 15 local branch tips were ancestors;
nine worktrees existed. [Preflight receipt](preflight.json#L1).

The incremental bundle passed its recorded SHA-256, bundle verification and
preserved manifest checks before use. Its hash is
`6cd9bc1a438df458f83d64335589889203babd781fc18b87ee44079955766553`.
The inherited sparse-history connectivity limitation in the prior handoff
remains recorded: no new full-tree connectivity claim is made. The office Mac
was unavailable to the owner. The existing cluster GitHub SSH identity worked,
so the requested candidate branch was pushed without force and verified at
exactly `871dbd7` before predeclaration and numerical execution. The original
HTTPS-authentication failure receipt was not rewritten.

Job 19185 rechecked all 354 external files, source-bound audit definitions,
conditional-roundoff evidence, physical tables and owner-supplied adjudication:
361 input identities, all passing before and after. The adjudication provenance
remains the original owner-supplied artifact, without claiming an independently
signed reviewer digest. The original D1 evidence and three original maintained
files also compare byte-for-byte against commit 871dbd7 before and after the
calculation. Original inputs were read only.

All **17 analytic tests passed** before retained integration: the original 11,
plus normalized-segment and tiny-piece exact integrals, unequal-piece
fourth-order convergence, polynomial equivalence and sensitivity checks. The
single numerical allocation completed with exit code 0 on phyn004, one logical
execution CPU (affinity [0]), 4 GiB, 30-minute ceiling; elapsed Python time
10.323 seconds, max RSS 44,232 KiB. Slurm reserved two hardware threads but the
process executed on one. No repeated retained-data calculation occurred.
See [allocation](evidence/job-19185/allocation.txt#L1),
[completion](evidence/job-19185/completion.txt#L1),
[test log](evidence/job-19185/tests-aligned.log#L1),
[execution log](evidence/job-19185/execution.log#L1) and
[capture hashes](evidence/capture.json#L1).

Files added are the maintained aligned-reference, test, launch and reporting
scripts under `tools/linux/native_qualification`, this sibling resolution
directory, and an addendum to the refinement report. Change classes:
numerical-method **test candidate**, engineering, documentation, generated
evidence. Governing authority: owner's bounded resolution task; review Section
D; GOVERNANCE §§2–3/5; ADR-0017 §§5–7; INV-11f/INV-15. No production refactor or
pre-baseline exception is involved. Scientific impact: new offline reference
numbers only; no production solver output, model or setting is changed.
AI-authored work remains a scientific candidate pending domain-owner review.

## Scope and owner decision

The local chemical reconstruction is now qualified **conditionally for these
13 intervals**, under the stated thermal/reference assumptions. The main-prefix
error from t=0 and thermal RHS roundoff obligations remain open. This does not
qualify every Stage A checkpoint, complete paired-run R18/R20 or restart,
chemical heating in the one-year temperature observable, or Stage B readiness.
Any subsequently authorized Stage A attempt must diagnose its **actual output
checkpoints**, not inherit universal qualification from this study.

The sole proposed future production change remains O2 chemical absolute
tolerances **(1e-24,1e-24) -> (5e-27,5e-27) MeV**, with all other settings
unchanged. This receipt is returned for the owner's decision. The change was
not implemented; no stellar solve, Stage A retry, thermal retry, Stage B run,
master merge or certificate promotion was performed. This bounded task stops
here, with no additional numerical investigation.
