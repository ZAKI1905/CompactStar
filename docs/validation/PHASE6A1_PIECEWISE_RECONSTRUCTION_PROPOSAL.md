# PROPOSED ADR-0017 amendment: event-aware local reconstruction

Status: **PROPOSED — NOT ACCEPTED; NO IMPLEMENTATION OR CAMPAIGN AUTHORITY**.
Date: 2026-09-30. Classes: numerical-method and reconstruction architecture.
This document does not edit or supersede accepted ADR-0017. It records a concrete
option for owner decision after one independent scientific/numerical review.

## Problem and evidence

The [bounded observation-32 report](PHASE6A1_OBS32_FORENSICS.md) and diagnostic
commit `86c8978` preserve exact reproduction of the archived unresolved O1/O2
checkpoint. Cstar is continuous and linear in log(T); its dC/dT jumps by 7.5094%
at the crossed cache knot. The complete thermal RHS is continuous but its
x-derivative has the predicted kink. There is no demonstrated context, lookup or
budget defect. Tight unsplit rk8pd converges toward the split result. Three nearby
smooth controls pass, and split RKF45 O3 agrees within 3.19e-16 in x.

The prototype's split O1/O2/O3 equality uses one common O3 event time. That is good
local evidence but **cannot bound common event-location error**. Its root bracket
is a floating-point numerical-IVP bracket, not a certified exact-IVP bracket.
No production acceptance can be inferred from this prototype alone.

## Exact proposed semantic change

Replace the uniform *unsplit* strict-interior rk8pd replay with a uniform
*event-partitioned* strict-interior rk8pd replay. Apply it prospectively to every
checkpoint, regardless of whether the old rule would fail. There is no failure
trigger, retry, knot-specific alternate solver, result selection, or tolerance
retuning. A no-crossing interval has one segment and uses the same rule.

1. Known interpolation interfaces are numerical metadata supplied by the frozen
   RHS owner. Cstar interfaces are the actual authenticated cache temperatures,
   not a separately regenerated approximate logarithmic index. Expose immutable
   metadata through an explicit interface; the diagnostic private-member accessor
   is not a proposed production API. Do not alter tables or interpolation values.
2. For each disposable O1 and O2 context, start at the exact saved accepted left
   state and solve to the requested observation. Use the unchanged governed
   rk8pd tolerance vectors on each smooth segment. Retain the accepted
   main-state bypass when an observation exactly matches an accepted endpoint.
3. Find each encountered interface from the condition `Tinf(t)=T_knot` using a
   bounded event solve. Preserve the full physical state continuously; do not
   project x, eta or energy onto the interface. Reset local numerical solver
   state between segments. Internal boundaries affect only the disposable replay,
   never the persistent authoritative main solve or scientific schedule.
4. Event discovery must account for direction, repeated crossings and endpoints.
   Endpoint temperatures alone do not prove the absence of an interior crossing.
   A validated monotonic segment certificate or equivalent safeguarded event
   search must exclude missed crossings. If that cannot be established within
   fixed work limits, return unresolved. No extrapolation or silent continuation.
5. Localization error is part of reconstruction uncertainty. Preserve
   `d=abs(O2-O1)`, `D_O1`, `F_O`, and `F_i` exactly as currently defined. Let
   `U_RK=2 max(d,F_O)`. Add a componentwise conservative propagated event-error
   bound `U_event`, covering all crossings and endpoint effects, and require
   **both** `d<=D_O1` and `U_RK+U_event<=0.20 F_i`. Report O2 only on success.
   Thus event handling consumes the existing budget; it does not enlarge it.
   This additive accounting is a proposed change, not the current ADR rule.
6. A shared root alone must not qualify `U_event`. Independently resolve event
   location and endpoint propagation using tighter numerical evidence plus a
   justified bound on localization and sensitivity, including binary64 floors.
   Propagate bracket uncertainty through both adjoining segments. A small root
   residual or identical tier results alone is insufficient. The concrete bound
   must be predeclared and independently reviewed before implementation is called
   qualified; the present experiments do not supply that proof.
7. Preserve full/cheap currentness contracts, immutable owners and isolated
   evolving/GSL state. A cache hint may be shared only where lookup invariance is
   established; it cannot become event identity or provenance authority.
8. Diagnostics are evaluated from the qualified reconstructed state. Their error
   qualification must include event uncertainty. R18, R20's unchanged scientific
   grid/composite trapezoid and uncertainty allocation, source/ledger equations,
   frozen ceiling and candidate gate chain all remain binding. Internal knot
   endpoints are not extra R20 grid points or extra thermal heat terms.
9. Export event metadata, localization brackets, residuals, directions, solver
   states/counts and componentwise uncertainty contributions. Any missing event,
   ambiguous ownership, unbounded uncertainty, GSL pathology or currentness
   failure is unresolved and halts the campaign under the same fail-closed rule.

This is an architectural contract for a proposed implementation, not finished
production code. In particular, item6 is a required numerical design/validation
obligation; acceptance of a direction must not be represented as proof of it.

## Alternatives requiring an owner choice

| Alternative | Evidence / tradeoff |
|---|---|
| Retain ADR-0017 unchanged and keep the campaign stopped | Correct current disposition; provides no new eligible checkpoint. |
| Uniformly tighten unsplit reconstruction | O3/O4 support this bracket's local value; O4 is near binary64 limits. Must qualify all required smooth/knot cases and costs, rather than infer general sufficiency from one bracket. |
| Uniform event-partitioned rk8pd (proposed option above) | Strong bounded split evidence; requires event discovery/error accounting and explicit ADR amendment. |
| Another uniformly qualified local method | Possible, but independent unsplit RKF45 here is nonmonotonic; no replacement method is selected or qualified by this investigation. |

Changing Cstar interpolation/physics, relaxing the 0.20 budget, special-casing
observation32 or choosing a passing result after failure are excluded alternatives.

## Required evidence before production adoption

Predeclare a separate bounded qualification covering smooth intervals, exact
accepted endpoints, exact/adjacent knot temperatures in both directions, multiple
interfaces, rejected stage crossings, endpoint knots, nonmonotone/unsupported
searches and fail-closed work limits. Verify independent event-error accounting,
method comparison, diagnostic uncertainty and invariance to valid cache histories.
Use analytic piecewise-smooth test IVPs and the preserved production bracket.
Reauthenticate unchanged source/physics inputs and demonstrate main-history
passivity under changed observation schedules. None of these new ODE experiments
is authorized by this proposal alone.

The saved BASELINE main knot-spanning step differs at its right endpoint from two
agreeing local O3 methods by 2.38e-7 in x. The review must explicitly disposition
this finding and the distinct roles of main hierarchy error and reconstruction
error; successful local reconstruction cannot certify a coarse main trajectory.
No main-method change is proposed or authorized here. The deferred fresh
BASELINE/REFINED/ULTRA and matched-control hierarchy remains a later, separately
authorized campaign, after architecture acceptance and implementation qualification.

## Decision boundary

Recommended next action: one independent scientific/numerical review of the exact
reproduction, root/cell/RHS evidence, proposed event-error contract and main-step
finding; then an explicit owner decision among the alternatives. Do not implement
this amendment or resume the campaign automatically. Historical BA12, BA12R and
the observation32 STOP remain immutable FAIL. No candidate, physical BNV model or
rate, A18, superfluidity, Regime-II/MixedStar or variable-Z evolution is authorized.
