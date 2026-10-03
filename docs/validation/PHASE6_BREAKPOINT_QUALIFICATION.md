# Phase-6 breakpoint correction — bounded candidate qualification

**BOUNDEDLY QUALIFIED NUMERICAL CANDIDATE / PROPOSED AMENDMENT / NOT RATIFIED.**
The full campaign remains stopped. Ready for **one independent review**, followed
by owner acceptance; no merge or campaign restart is authorized by these results.

## Authority, identity and scientific impact

Branch: `physics/phase6-breakpoint-aware-evolution`, descending from forensics
`669898c6741883ecacdc7e2fa2599cfa3df92918` and canonical
`9a5c1eca40758e4c39f9f80db7324e12ad23e505`. The final numerical implementation is
`fa01369654f93a3ae473e6f51fe2a480d0bdfc26`; the subsequent crossing harness only
adds explicit reference-uncertainty and roundoff records. The evidence manifest
records exact executable and source hashes for each run.

Classes: numerical-method, architecture, engineering, documentation and generated
evidence. Authority is the owner's 2026-10-03 bounded implementation instruction,
accepted ADR-0017's passive-output contract, and its explicitly **proposed**
[Amendment 1](PHASE6A1_PIECEWISE_RECONSTRUCTION_PROPOSAL.md). The successor
[validation specification](PHASE6_BREAKPOINT_VALIDATION_SPEC.md) was committed
before new scoring. ADR-0015/0016 source/ledger semantics, Cstar table values and
interpolation, and Phase-5 numerical mathematics are unchanged. Numerical main
states and accepted histories intentionally change where the old method crossed
cache cells. Historical BA12, BA12R and obs-32 records remain byte-identical.

## Preserved review and permitted historical use

The actual transferred files were found in `/Users/keeper/Downloads`, authenticated,
and imported separately in commit `0441b75813c775bd33e984d3f8b9a7f12ae25b49`. The original archive is
853543 bytes and was not repacked. All 76 internal SHA256SUMS entries verified;
the standalone and bundled reports are byte-identical.

- Archive SHA-256: `5ba87000f94f4b91743c91f62ac024f459bb22922c962aa47d899571c8763b4e`.
- Review SHA-256: `85403e6b13ef8cbfc8e71a9a40eba662c6e2d456270a785bd069a171fd72c0f3`.
- Reviewed SHA: `669898c6741883ecacdc7e2fa2599cfa3df92918`.

The [import README and validity index](evidence/obs32-independent-review/README.md)
retain the exact member identities. Used: the readable review/README, archived
source/control states and crossing inventory, the reviewed fixture assembly,
GSL cache-reset demonstration, forensic cases, and run4's **complete** observation
files. The sign-correct run2 `analysis.json` is the authority for historical
control predictions; the new comparisons do not use metadata's wrong sign.
Valid run2 reproductions, fixed-step/event perturbation studies, replays and
unsplit results support the diagnosis and remain available without unnecessary
expensive repetition. Run3's in-memory gate/error summaries are retained as valid
summaries over 1277 observations per the README, not as substitute rows.

Excluded: run1 (CWD/context failure); every run2 split M-* row/column (stale GSL
state); the **entire** run3 M-control-observations.tsv (corruption from data row
599, including repeated rows); metadata's incorrect control-prediction sign.
No historical file was repaired, truncated into a substitute, or retrospectively
promoted to valid evidence. Preservation alone confers no numerical validity.

## Implemented architecture and rollback

`StarContext::ReadHeatCapacityCells` provides a typed, provenance-bearing copy of
actual cache temperatures/capacities and existing cell-selection semantics.
`CstarSegments` uses those knots, never an analytic grid approximation, and keeps
the existing RHS bitwise within the active cell. Temporary RK stage excursions
use that cell's smooth expression. Cache values and the interpolation law are
unchanged; the historical analytic classifier remains only historical metadata.

The reusable three-state `BreakpointStepper` drives both main RKF45 and local
rk8pd. Each solve locates its own crossings, carries all three integrated state
components continuously, and retains a useful step proposal across boundaries.
A discarded trial restores t/y/h and resets **both** GSL evolve and step state.
Root trials also start with fresh caches. If no representable time lies inside a
one-ULP crossing bracket, it retains the integrated right state and the full
bracket. There is no state projection. Stepper and accepted-step data structures expose cell/knot/direction, event
residual, bracket and time ULP, local error, numerical rejection, event work,
discarded trials, resets, stage excursions and initialization separately.

Exact scientific endpoints return stored main states. Interior reconstruction is
breakpoint-aware. Scientific times never enter the main stepper or event solve;
boundary states are not appended to the R20 quadrature grid. The qualification is
serial within each process; no concurrent shared-cache use is qualified.

## Uncertainty rule and the exact-control repair

Keep the original d_O, D_O1, F_i and U_O allocations. Require d_O <= D_O1 and
U_O + U_event <= 0.20 F_i componentwise. The candidate does **not** adopt the
review's scalar propagation formula as a universal coupled bound. Instead,
U_event_i = 2 max(|O2_i-H_i|, |H_i-I_i|) + 64 ulp(|O2_i|), with tighter half-step
rk8pd H and quarter-step independent RKF45 I, each owning its events. Units are
one for x and MeV for the two etas. Endpoint propagation factor is one because
full coupled propagation has already been integrated to the observation. This
is an empirical bound in the tested domain, not a global theorem.

Further eighth-step rk8pd and sixteenth-step RKF45 solves from the same saved
left state test the envelope, including reference disagreement. All diagnostic
state-box corners are evaluated with a factor-two margin. Absolute thermal,
luminosity and chemical-energy contributions enter R20 uncertainty without
cancellation. The reported envelope checks cover all three state components and
18 diagnostic quantities. The original obs-32 candidate error/envelope ratio is
0.00387239; smooth cases 26, 36 and 18 retain bitwise old O1/O2 states.

Development run05 correctly refused 189 control diagnostics: propagating the
unchanged chemical atol floor into Echem at exact eta=0 produced about 9.22e5 MeV
against a 2e-12 MeV allocation. The repair preserves an owner-certified exact
subspace: identically zero source, static zero spin, governed H(0)=0, and both
initial etas exactly zero imply eta(t)=0. Every local stage/derivative is checked,
all four contexts must agree on the certificate, and the diagnostic owner must
retain it. False certificates are rejected; nonzero starts retain full uncertainty. Allocation
fields remain unchanged; only the exact chemical error is zero. The control
history, observation states, reference errors and allocation TSV are **bitwise
unchanged** after this repair. No budget or tolerance was relaxed.

## Bounded accuracy and original gates

Each state ratio below is (absolute error + qualified reference uncertainty)/D_i;
one is the limit. Ledger ratios use the retained 1e-9 G_P accuracy scale for
REFINED/ULTRA. Envelope ratios include independent reference disagreement. R20
limits remain 2e-4, with reconstruction contribution/N <= 5e-6. Maximum normalized
R18 residual across these windows is 1.31872e-16.

| Window | State ratio | Ledger ratio | Error/envelope | max abs(R20/N) | max U_recon/N | Successor |
|---|---:|---:|---:|---:|---:|---|
| source | 0.929636 | 0.0034748 | 0.0349417 | 1.2873e-05 | 4.88392e-11 | PASS |
| control | 0.276599 | 0.936051 | 0.0120531 | 1.24901e-05 | 2.16404e-10 | PASS |
| bounded | 0.911509 | 0.00398206 | 0.0443038 | 5.25076e-07 | 1.60816e-10 | PASS |
| tail | 0.345048 | 2.11035e-13 | 0.0920124 | 8.74729e-10 | 4.11186e-11 | PASS |

The windows are source 56152132796.889648..1e12 s; control
6540338562008.4297..9e12 s; the prior 0..462269531250 s case; and tail
13008696323235.514..13508696323235.514 s. Warm windows start from the same archived
physical states as their references. They do not remove earlier prefix error or
prove full-campaign global accuracy. The prior case starts from x=eta_e=eta_mu=0.

Original BA12R-form failures on the **new bounded trajectories** remain failures:

| Window | State c1/c2/c3 failures | Ledger c1/c2/c3 failures | B/R/U rejection fractions (old limit 10%) |
|---|---|---|---|
| source | 0 / 0 / 0 | 0 / 0 / 0 | 5.41%, 3.55%, 2.21% |
| control | 0 / 23 / 12 | 0 / 98 / 2 | 0.00%, 16.67%, 12.50% |
| bounded | 0 / 2 / 0 | 0 / 0 / 0 | 5.66%, 4.08%, 1.93% |
| tail | 0 / 2 / 0 | 0 / 0 / 0 | 14.29%, 22.22%, 15.38% |

These are not replacements for historical campaign verdicts, which remain FAIL.
The successor establishes accuracy independently of an assumed fixed 0.10
contraction ratio. It reports ratios as convergence evidence, including failures;
it does not label every failing point boundary-limited. A predeclared record-only
classification finds one common single-segment boundary-limited interval in source
and one in control, and none in the bounded/tail cases. Other adaptive step sequences
really differ. Independent-method agreement, explicit within-cell refinement and
reference uncertainty supply the positive accuracy evidence for the successor.

## Coverage, stability and cost

All 37 archived crossing intervals (8 source, 29 control), three tiers and three
components pass: 333 endpoint/event rows. Worst endpoint ratio is 0.130807;
worst event-state ratio is 0.0725445. Maximum endpoint/event reference uncertainty
is 0.00218015/0.00172813 of D_ULTRA, below the predeclared 0.05 allocation.
Run4 comparisons use all 490 source and 1277 control times; maximum tier-scaled
errors plus reference uncertainty are 0.934930 and 0.277794. Archived-reference
agreement is within 0.002063/0.003571 D_ULTRA; no uncertainty was invented for the
historical reference itself.

Old Arm E reproduces exactly (232 accepted, 60 rejected, exact final state).
Old obs-32 reproduces its failure exactly: d/D=3.2368362848836014 and
U/(0.20F)=3.1220859973768409. Corrected main histories are qualified separately.
Analytic coupled heating/cooling, knot starts, direction changes, stage overshoot,
independent event/endpoint solutions, missed/wrong-cell crossings, stale rollback,
shared-event cancellation, artificially identical O1/O2, exact endpoints and both
nextafter sides are covered. Main accepted times/states/proposals/rejection counts
are bitwise invariant under alternate passive schedules in each window. R20 uses
only the declared 65/65/241/33-point scientific grids.

| Window | B/R/U accepted | rejected | boundaries | event trials | RHS calls |
|---|---|---|---|---|---|
| source | 70, 163, 399 | 4, 6, 9 | 7, 7, 7 | 23, 21, 22 | 776, 1516, 3430 |
| control | 6, 10, 14 | 0, 2, 2 | 5, 5, 5 | 21, 16, 15 | 222, 226, 250 |
| bounded | 50, 94, 203 | 3, 4, 4 | 7, 7, 7 | 26, 24, 23 | 634, 976, 1840 |
| tail | 6, 7, 11 | 1, 2, 2 | 0, 0, 0 | 0, 0, 0 | 55, 69, 101 |

The tail's 14.3–22.2% rejection fractions fail the original gate. Its bounded state,
ledger, R20, progress and cost checks pass, but the direct GSL linear stability
probe confirms a limitation: local eigenvalues are approximately -4.77e-12,
-3.86e-11 and -3.15e-11 s^-1; at the largest accepted ULTRA step (~1.466e11 s),
the amplification spectral radius is 8.704, versus 0.705 at half that step.
Jacobian step refinement changes the scaled matrix by about 4e-10 relatively.
Thus a near-equilibrium error-controlled step can be locally unstable even while
this finite window meets accuracy. This is not dismissed or converted into a
higher rejection threshold; longer-term stability remains a campaign limitation.

## Engineering, regressions and remaining limits

Checked writers detect open/write/flush/close failure and validate bytes, schema,
row count, ordering and completeness before atomic publication. Controlled write
failure, reordered/missing rows and schema corruption are tested. Calculation
completion is distinct from successful artifact publication and numerical
qualification. Isolated outputs remain small; no historical evidence or unrelated
files were deleted. The final package includes no executables or build trees.

All 13 affected regression tests pass, including cache/thermal contracts,
Phase-5 component tolerances, chemical exact/curved oracles, particle numbers,
BNV energy/ledger, passive output and the new breakpoint/currentness tests.
The two shared Analysis changes are checked full-byte block reads only; they keep
every currentness check and all floating-point mathematics. Both source mutation
and missing-file refusal are tested. External library/header hashes, protected
baseline trees, Rotochemical/EOS code and historical failure documents are unchanged.
The upstream Phase-5 knot issue is recorded as a
[separate follow-up](PHASE5_CSTAR_KNOT_FOLLOWUP.md), not implemented here.

The empirical envelope is limited to these fixtures, methods and windows. Sampled
stage monotonicity is not a theorem excluding every unsampled reversal. No full
500 kyr run, global stability proof, physical BNV model/rate, superfluidity, A18,
Regime-II/MixedStar or cluster work was performed. The recommendation is one
independent review of this working bounded correction and amendment, then owner
acceptance. Full-campaign readiness is **not established**.

Inspect the [evidence manifest and assessments](evidence/phase6-breakpoint-candidate/README.md)
for exact commands, hashes, raw compact results, failure records and provenance.
