# BA13 frozen-background closure predeclaration — 2026-10-19

UNVERIFIED SCIENTIFIC CANDIDATE. This declaration precedes all newly scored
closure calculations. Entry: 286aafb63002093b3052e23a3972903011f36035,
clean candidate/linux-x2-fixed-step-20261016. Descendant branch:
candidate/linux-ba13-closure-20261019. Isolated worktree:
/mnt/sdd/zaki/CompactStar/linux-transition-20261010/worktrees/ba13-closure.
The other five worktrees are retained. All local commits touching the proposed
validation files belong to the entry's ancestry (c5a0ff4, 232a7e7, 97b2d13).
There is no competing implementation of these files. The independent review's
40 SHA256SUMS entries verify. Its unchanged reproduction script, copied outside
the review, passes all 31 authentication/reproduction checks; all 884 input
manifest entries verify. The fitter replay reproduces both certificate files
byte for byte. Authentication is recorded in evidence/authentication.json.

Authority: GOVERNANCE.md §§1–5, AGENTS.md §§1–11; accepted Phase-6A-1 preflight
§11 (docs/validation/PHASE6A1_CONTROLLED_BNV_IMPLEMENTATION_PREFLIGHT.md:655);
INV-02/03/04/06/14/15; accepted ADR-0009/0013/0015. ADR-0021 Option A and
interpretation (2) remain PROPOSED. Changes are validation numerical-method,
engineering, documentation, and generated evidence. No production physical
owner or historical result is changed. Numerical/validation policy slots are
unwritten; they supply no additional authority. The Mac literature directory
in AGENTS.md is unavailable here; this task uses the authenticated review's
explicit derivations and the existing source equations, with no new physical
source model.

## Corrections D1–D6

D1: explicitly distinguish historical production P, fixed-production-profile
reconstruction R, independent X2 X, and physical solution Q. At common coordinates,
P-Q=(P-R)+(R-X)+(X-Q). Use signed paired differences for deterministic P-R;
use independently bounded endpoint structural errors unless shared sensitivities
are justified. X2 self-error belongs only to X-Q. Retain actual historical units,
all component scales, all 13 thermal knots, and Z's maximum-entry row norm.

D2: retire BOTH first-order signed and absolute defect transports from closure
qualification. Their existing code and failed records remain historical evidence.
The closure executable refuses to enable that estimator. Use the prescribed GL4
refinement experiment instead; Richardson agreement is a numerical estimate,
not an analytic remainder bound. No expensive estimator repair is planned.

D3: qualification requires a finite nonnegative value and appropriate evidence
for EVERY required term. Classes: analytic_bound (with input hypotheses),
validated_numerical_estimate (named accepted standard/authority),
conditional_estimate (named unresolved hypothesis; requires explicit acceptance),
and unresolved. There is no default acceptance of a conditional estimate and no
invented owner approval. A propagation theorem does not upgrade its inputs.

D4: evaluate dt_n/dx, dt_e/dx, dt_mu/dx, and propagate the measured B-functional
coordinate discrepancy in ABSOLUTE tangent units. Record uncertainty of the
sensitivity and of the functional itself. Root residuals do not replace B errors.

D5: distinguish continuous physical equilibrium response from intermediate
production-algorithm output. Verify 21 sampled topologies, complete radial grids,
neutron cell and envelope node, threshold distances, slopes, and margins. Equal
sampled indices alone prove neither monotonicity nor absence of a hidden switch.
Do not add a fictitious production-node jump to the physical response. Continuous
selection claims need a derivative bound or named conditional Lipschitz estimate.

D6: use charged-to pointers for definitional zero terms, wire root_error for X2
root alignment only, and never apply long-double roundoff to production binary64.
Unenclosed arithmetic/stage/elementary-function effects remain named gaps.
Preserve original source and certificate identities and all X1/X2/P5/P7 failures.

## M1–M6: frozen calculations and refusal criteria

M1: at fixed retained Hc strings for targets 0,10,20, evaluate passive functionals
at N=1024,2048, using retained N=512. Verify two signed increments and fourth-order
compatibility (p in [3.5,4.5]); outside this range report INCONCLUSIVE rather than
assert a Richardson bound. Report |last increment|/15 as the fine-grid tail
estimate; for N=512 use |Q1024-Q512|*16/15. No arbitrary floor/multiplier.
Include existing EOS96/192 and center-small/large observations, keeping zero
spreads and unseen remainders explicit. Add the third relative stencil width
2.5e-5 at N=512 (six solves), test h² ratio near 4 (range [3,5]), and retain endpoint
integration uncertainty. Decimal strings, not binary64 output parsing, govern
small differences. M1 target gate: R=|delta_P-delta_X|+B_X below the unchanged
23% representation share; discrepancy alone above share is FAIL; B_X >= measured
discrepancy is INCONCLUSIVE. One additional N=4096 passive-functional resolution
at the same three Hc values is allowed ONLY if M1 remains inconclusive and can
reasonably become decisive within budget. Unresolved EOS, stage, center, tangent
or arithmetic remainders prevent an unconditional band regardless of magnitude.

M2: compare complete radial grids at all 21 retained profiles before magnitude
transfer. Compute deterministic signed representation correction separately at
each sample. Candidate structural magnitude is max_three |R-X| plus the X2 band,
with independent endpoint addition, not assumed signed correlation. Grid identity
is necessary for the proposed transfer, never sufficient. Transfer requires an
independent sensitivity/enclosure or prospective validation over the full target
window; otherwise mark it conditional/unqualified, stop dependent numerical
expansion, and name the missing condition. Do not launch 21 structural solves.

M3: evaluate every one of the 20 achieved-B gaps. In fixed topology,
Q_xx = B0² (Q_HH B_H - Q_H B_HH)/B_H³. Obtain derivatives from the existing
centered stencils and test width dependence. For tangents, differentiate N_H/B_H;
third derivatives require additional justified evidence. Test K=10*max_three
|Q_xx|, prospectively checking narrower stencils where available; three-target
curvature variation >2 is INCONCLUSIVE. K*dx²/8 above the fixed 10% share is FAIL.
Factor 10 and sampled agreement are hypotheses, not a global theorem: without a
validated supremum the term remains conditional. Expand uncertain achieved-B
gaps explicitly. Include all 21 residuals of the unchanged outer-half fit.

M4: propagate both reconstruction and X2-versus-production B-functional errors,
including their common B0 coordinate and denominator. Propagate same signed
corrections through historical (I-H), perturb achieved-B design and weight floor,
and avoid adding the same source again as fit covariance. Test each M4 term
against its original planning share; an unbounded remainder stays unqualified.
Absolute tangent normalization is exactly S=1, never |t0|. Wire the reference
root_error with a positive local B slope hypothesis, not a production-root proxy.

M5: recompute full 3x3 inverse and 2x2 Z independent-endpoint majorants for all
21 samples using each G, full background+center+actual-cut tail set, distinct
Jacobians and BOTH nonlinear Neumann remainders. Compare with the review's
~5.6e-7 muon-row result. No perfect B0 correlation. Require rho<1. Per-target
background magnitude transfer remains named conditional until justified.
Do not double charge the same background through two alternative comparisons.

M6: run all 14 existing enclosure tests and relevant existing BA13 tests; add
negative/mutation controls for unaccepted conditional evidence, unaccepted
numerical standards, absent terms, invalid/NaN values, exact checks, changed
thresholds, wrong reference, common-B omission and inverse nonlinear tails.
Missing or conditional totals cannot PASS. Exact P2 and support identities retain
logical rules. No historical fitter edit, certificate promotion or ratification.

## Extension, resources and stops

Only after M1–M6 meet their criteria extend the qualified construction to other
groups. Reuse all retained profiles for static diagnostics even if qualification
stops. A possible validation-only node observer would use the existing X2
collocation path on the ORIGINAL source-support coordinates with independent
resolution checks; no physical solver change. Do not implement/run it if the
muon transfer or interval prerequisite remains unsupported. Record the exact
small remaining nodewise evidence for metric_structure and mu_n_profile.

Initial numerical ceiling 900 CPU seconds; total task numerical ceiling 1800
CPU seconds including any extension and validation. One pinned core per Slurm
job, 1 GiB, single-thread BLAS; per numerical process 90 CPU seconds, 120 seconds
wall, controller hard check before each case. Benchmark N=1024 before remaining
runs. Predeclared initial work: six N=1024/2048 runs (~120 CPU s), six third-width
runs (~30 CPU s), static M2–M5/authentication/tests (<60 CPU s). Builds and hashing
receipts are separate. Before any extension publish an itemized CPU estimate.
A source/runtime replay must match retained nominal X2 decimal output.

Stop numerical expansion at first predeclared FAIL, invalid estimator,
unautheticable evidence, unjustified transfer, or demonstrably unclosable gap.
Complete independent code repairs, static evidence, matrix, report and owner memo
when such a stop applies. Do not adjust thresholds, shares, source model, fit
window, historical certificate bytes or uncertainty budgets to obtain PASS.
No physical BNV or thermal trajectory, broad EOS sweep, blanket star sequence,
canonical merge, push without authenticated access, or certificate promotion.

Deliver a complete reproducible 21x22 matrix with explicit terms, units,
correlations, evidence classes/authority and reasons; proposed_total=null when
unqualified. Disposition A–E is scientific, separate from implementation tests.
Keep large artifacts outside Git with verified manifests. Supply an incremental
bundle against the exact starting SHA if GitHub authentication is unavailable.
