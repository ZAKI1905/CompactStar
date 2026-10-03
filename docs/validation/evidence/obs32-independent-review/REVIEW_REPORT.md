# Phase-6A-1 obs-32 independent review — disposition **B**

Reviewed: `analysis/phase6a1-obs32-reconstruction-forensics` @ `669898c6741883ecacdc7e2fa2599cfa3df92918`
(entry `7d2777561a17b56b466fb024ef13104da6c0e86a`, canonical master `9a5c1eca40758e4c39f9f80db7324e12ad23e505`).
Review date: 2026-10-01. Read-only: no repository, worktree, ADR or production file was modified.
Evidence: this bundle (see `README.md`).

**Disposition B — ROOT CAUSE CONFIRMED; MAIN-TRAJECTORY IMPLICATION REQUIRES MORE BOUNDED NUMERICAL WORK.**

The forensics' root cause is right, and it is now a quantitative model. But the proposal leaves out
the main trajectory, and the evidence shows the main trajectory is where most of the damage is. The
authoritative RKF45 under-resolves **every** Cstar knot crossing:

- BASELINE local errors range from 0.07× to 162× its own tolerance across 37 archived crossings.
- In bounded windows, errors at REFINED and ULTRA are 26–95× their tolerances.

Making knots integration boundaries in the main fixes that: every tier's true error drops to ≤0.93×
its tolerance. Even so, the predeclared BA12R contraction rules still fail in the control's
late-cooling regime. So the amendment needs widening and one more bounded qualification before
ratification.

## What was run

Everything was bounded and local; no trajectory of more than ~70 main steps was run.

- **Build:** a scratch executable linked to the existing `build/obs32` library, whose hashes match
  `packages.json`. All six inputs re-authenticated. It ran from a scratch working directory using
  byte-identical copies of the four source files the fixture reads.
- **Reproductions:**
  - the production reconstructor on obs 32;
  - bitwise replays of 60 archived main steps from their saved states and saved step-size proposals;
  - 67-step and 6-step continuations of the archived main.
- **Comparisons:**
  - a tight reference: rk8pd at rtol 1e-15, split at every knot;
  - fixed-step order studies;
  - event-time perturbation;
  - three-tier windows (unsplit, split, split with a per-segment "frozen" Cstar cell) scored with
    the BA12R state-gate logic;
  - small direct C tests of GSL's own behaviour.
- **Harness bug, caught and excluded:** the first split-main windows were invalid. GSL's evolve
  object reuses the last step's end derivative, so undoing a step without `gsl_odeiv2_evolve_reset`
  corrupts the next one. Confirmed directly against GSL, fixed, rerun; only corrected runs are used.

## Core mechanism (items 1–3, 6)

- **The kink:**
  - Cstar is interpolated linearly in log T (`StarContext.cpp:793`).
  - The tabulated Cstar is exactly proportional to T in this range (local power 1.000000 for knots 62–102).
  - So **every** knot has the same slope jump: r − 1 = **7.50942544%**, with r = 10^(5/159).
  - In ẋ this is a jump in ∂ẋ/∂x of −κ·ẋ, with κ = 0.07244. η̇ has no Cstar dependence.
- **The error model:** an explicit RK step that straddles the kink at fraction θ of the step has
  local error [y'']·h²·K(θ), with K(θ) = Σbᵢ(cᵢ−θ)₊ − (1−θ)²/2. GSL's error estimate follows the
  same form with the weight differences in place of bᵢ.
- **Verified against GSL:** the kernels match GSL's rkf45 to 6 digits. GSL propagates the 5th-order
  solution.
- **How badly RKF45 is fooled:** the ratio of true error to estimated error is A(θ). Its median over
  θ is 11, it is 43 as θ→0 and 9 as θ→1, and it is unbounded near θ = 0.413 and 0.756.
- **Fit to the archive:** with exact θ, the model predicts all 37 knot-crossing local errors
  (8 source, 29 control) to within **0.2–6.2%, signs included**.
- **Step 21 specifically:** the model also matches GSL's estimate for each of step 21's five trials.
  The estimate-to-tolerance ratio fell 2.08, 1.70, 1.43, 1.23, then 1.083 — accepted, because GSL
  only rejects above 1.1 — while the true error was 18.4× tolerance.

## Main-trajectory windows (items 7, 21–25)

**Source window: t_L → 1e12 s, 490 observations**

| | Unsplit (M1 = archived, bitwise) | Knot boundaries (M2) | Boundaries + frozen cell |
|---|---|---|---|
| Max true error ÷ tolerance (B / R / U) | 18 / 71 / 66 | ≤0.85 / ≤0.74 / ≤0.83 | ≤0.87 / ≤0.90 / ≤0.93 |
| Worst REFINED-vs-ULTRA difference ÷ REFINED tolerance | **70.7** (1442 of 1470 fail) | 0.73 | 0.89 |
| Worst tier-to-tier contraction ratio | 1.80 (37 fail) | **0.036** (prediction 0.025) | 0.99 (1 fail) |

**Control window: late cooling** — worst REFINED-vs-ULTRA ratio: 27.3 unsplit, 0.023 split, 0.275
frozen. However, the contraction and floor rules fail 624 times (split) and 433 + 260 times (frozen).
Every failing difference is ≤0.28× the REFINED tolerance.

- **Why the control still fails:** knot spacing (~5.5e11 s) caps the step length for BASELINE and
  REFINED alike. The tiers end up sharing step sequences, so the hierarchy is no longer asymptotic
  even though every tier is accurate.
- **Persistence:** the knot-95 error decays with τ ≈ 7 kyr. It is still 1.6e-7 at +3 kyr and 3.9e-8
  at +16 kyr, and reaches roundoff by about +160 kyr. Source errors are fully damped by 500 kyr.
- **Control accumulation:** the control damps only weakly (propagation factors 0.40–0.99) and its
  errors are biased (21 of 29 positive). A linearized estimate gives ≈ +4.6e-6 in x at 500 kyr,
  about 22× its tolerance.
- **Context:** the historical BA12R maximum of 32.3 is consistent with this mechanism. Segmentation
  explained only 4–6% of that failure, and the BA12R forensics had already flagged knots (H5/H9).

## Event location (item 4)

- The root is unique: all 2409 logged right-hand-side samples in the bracket have
  ẋ ∈ [1.751, 1.793]e-12 s⁻¹.
- A smooth-side Newton solve converges in 3 IVPs to **zero** residual in the exact knot test. It
  lands 8.4e-5 s from the forensics' root.
- Because the right-hand side is continuous, misplacing the event only produces a second-order state
  error. Measured slopes match the model: 2.60e-17 per second of offset (model 2.63e-17) for an
  early boundary, 7.2e-18 (model 7.35e-18) for a late one.
- **Proposed bound:** U_event = Σ_k P_k·κ_k·δx_k·(β_M·|Δx_adj| + δx_k/2), where:
  - β_M = 0.0418 for rk8pd and 0.18 for RKF45;
  - δx_k = the event residual plus the segment's own RK uncertainty U_RK;
  - Δx_adj = the x-change over the step next to the boundary;
  - P_k = a predeclared bound on how the error propagates to t_obs.
  For obs 32 this is ≈ 3e-5 of F_O.
- **Combination rule:** d_O ≤ D_O1 **and** U_O + U_event ≤ 0.20·F_i, each level locating its own
  event (no shared root).

## Findings

**BLOCKING (2)**
1. The main under-resolves every knot crossing, and the proposal explicitly excludes main changes.
2. With the main unsplit, the BA12R hierarchy cannot pass: the worst REFINED-vs-ULTRA ratio reaches
   27–71 in bounded windows.

**MATERIAL (7)**
1. Split O1/O2/O3 agreement is a degenerate witness: every leg is one identical rk8pd step, so d = 0
   by construction. The real convergence evidence comes from independent split RKF45 and fixed-step
   order studies, and those do support it.
2. How the right-hand side is evaluated within a segment is unspecified. With the production lookup,
   the c=1 stage of the landing step overshoots the knot, and split RKF45 drops to order 4.0. Freezing
   the segment's cell restores full order: 1.0e-11 → 1.8e-15 at one step.
3. U_event is undefined in the proposal (derived and validated above).
4. The GSL evolve derivative cache must be invalidated after a discarded trial step. Keeping the
   step-size proposal across a boundary is safe: GSL leaves it unchanged after a final step clipped
   at a ceiling (verified).
5. BA12R contraction and floor rules fail in the late-cooling regime even with knot boundaries; this
   needs an owner decision before any run.
6. Knot errors in the control are biased and only weakly damped, which matters for
   source-minus-control comparisons.
7. A predicted BA12R §6.2 solver-trend refusal: in the final decade the source BASELINE has a 0.149
   rejection fraction from a stability-limited cycle at equilibrium, unrelated to knots.

**NONBLOCKING (4)**
1. The forensics' bisection bracket width is a floating-point width, not an accuracy statement
   (already disclosed).
2. The analytic `CstarCell` classifier must not be used as event identity; it did agree with the
   actual cache for all 37 archived crossings.
3. The original ADR-0017 qualification used ULTRA bracket lengths (~2.2e9 s). BASELINE brackets are
   ~8.1e9 s, and the error scales with length squared.
4. Phase-5D/ScaledRKF45 trajectories cross the same knots and need a separate upstream review
   (ADR-0017 §9.1).

**NOTE (5)**
1. Interpolating log C against log T would be exact and kink-free in this range, and would remove a
   mid-cell bias of up to 6.6e-4. That is a cross-phase physics change, not a fix for this campaign.
2. The background inventory found Cstar to be the only active piecewise-smooth dependence in the
   right-hand side.
3. The GSL kernel verification and the RKF45 propagated order are as stated above.
4. This Mac's data volume had ~534 MiB free of 460 GiB during the runs. It silently corrupted one
   scratch output file; evidence writers should fail closed on write errors.
5. The scratch evidence is not committed.

## Recommended amendment wording (ADR-0017 Amendment 1, not a new ADR)

- **A1 Scope:** amends §§4, 5.2–5.4, 7, 10, 11. Cstar values and knots, tolerances, ADR-0015/0016,
  R18/R20, the frozen domain, ScaledRKF45 and Phase-5D are all unchanged.
- **A2 Registered breakpoints:** the right-hand-side owner exports immutable, provenance-bound
  breakpoints: the exact cache knots T_k, the exact cell-selection test, and κ_k. Any non-smooth
  dependence that is not registered fails closed.
- **A3 Uniform rule:** every integration — the main at every tier and every reconstruction level —
  integrates breakpoint-free segments only, and no accepted step may straddle a crossing. Within a
  segment the right-hand side is the active cell's smooth piece, bitwise equal to production while
  the state is inside that cell. The rule is applied a priori and is never triggered by a failure.
- **A4 Main integration:**
  - observation times are never ceilings;
  - a trial step that crosses a breakpoint is discarded, the crossing is located on the integrator's
    own trajectory, and the integrator lands exactly there;
  - boundary states are accepted main states;
  - the step-size proposal carries over, and the derivative cache is invalidated;
  - the accepted history must not depend on the observation schedule.
- **A5 Reconstruction:** exact-endpoint semantics are unchanged; O1 and O2 each locate their own events.
- **A6 Uncertainty:** the U_event formula and combination rule above. d = 0 counts as a witness only
  on breakpoint-free segments.
- **A7 Fail-closed:** a missed or ambiguous crossing, a non-monotone approach, an exhausted work
  limit, a residual beyond tolerance, a GSL failure or a currentness failure all halt the run.
  No retry.
- **A8 Unchanged:** R20 grid and quadrature (boundary states are not grid points). Main-integration
  event error stays part of the main hierarchy, not reconstruction uncertainty. BA12, BA12R and the
  obs-32 campaign remain FAIL.

## Validation ladder; must pass before any BASELINE retry

- **V0** — ratification of the amendment, plus the breakpoint interface replacing the private-member accessor.
- **V1** — obs 32 still fails, bit-for-bit, under the old rule; the new rule resolves it without tuning. *(required)*
- **V2** — analytic piecewise IVPs and convergence of the event solve and split legs. *(required)*
- **V3** — independent-method agreement.
- **V4** — brackets without knots are bitwise unchanged. *(required)*
- **V5** — all 37 archived crossings: authenticated replay, and the new main stays within tolerance at all three tiers. *(required)*
- **V6** — bounded three-tier windows scored with the full BA12R state and ledger logic: early source
  heating, late control cooling, and the equilibrium tail, including §6.2 statistics. *(required,
  after the owner has decided BA12R applicability for windows where step lengths are capped by breakpoints)*
- **V7** — the accepted history is invariant to the observation schedule. *(required)*
- **V8** — exact and adjacent endpoints. *(required)*
- **V9** — U_event propagated into diagnostics and the R20 budget.
- **V10** — the old path reproduces Arm-E bit for bit, and the new path qualifies all 240 observations. *(required)*

## 45-item report

1. 669898c6741883ecacdc7e2fa2599cfa3df92918 — entry 7d27775 and master 9a5c1ec authenticated
2. Yes, read-only
3. Yes, bitwise — d/D = 3.2368362848836014, U/(0.2F) = 3.1220859973768409; chemical components pass
4. Unsplit O1–O4 reproduce bitwise. They converge to the split value (−3.9e-13, +1.4e-14, −6.4e-16,
   −4.2e-17), but not asymptotically: signs alternate and fixed-step "orders" range from −3.8 to +8
5. Split value 0.12472490938613233 confirmed, but O1/O2/O3 equality is degenerate (single identical
   steps). Convergence is supported independently
6. Split RKF45 O3 is +3.2e-16 from split rk8pd; O2→O3→O4 decrease monotonically
7. Observations 26, 36 and 18 qualify (d = 0, U/0.2F ≤ 0.0964)
8. Continuous, linear in log T
9. 7.50942544% (= r − 1), identical at every knot in range
10. Lipschitz and piecewise C^∞; ∂ẋ/∂x jumps by −κẋ; the solution is C¹ with a jump in y''. A step
    straddling the kink has O(h²) local error with a θ-dependent constant
11. Unique (ẋ > 0 throughout); every knot is crossed once
12. Very well conditioned in state, at binary64 resolution in time; the error it induces is second order
13. U_event = Σ P_k·κ_k·δx_k·(β_M|Δx_adj| + δx_k/2), validated to ~2%
14. d ≤ D_O1 and U_O + U_event ≤ 0.20 F_i, with events located per level
15. A uniform rule if it applies to every integration a priori; a forbidden fallback only if
    triggered by failure
16. Yes — same t_R, same y_L, bitwise main replay
17. δx = +2.3800513596072825e-7
18. δT/T = 2.38005e-7
19. δT = 27.09 K
20. 18.4× the BASELINE tolerance scale; 17.0× its own estimate; 1.84e5× F_x; δPnet/G_P ≈ 8.2e-8;
    R20 effect ≈ 7.5e-9
21. M1: bitwise to the archive; errors up to 18× tolerance (BASELINE) and 26–71× (REFINED/ULTRA)
22. M2: ≤0.85× tolerance (x component 0.043×) at negligible extra cost
23. M3: rk8pd rtol 1e-15 split at every knot, self-consistent to ~1e-16
24. Yes — true errors ≤1× tolerance at every tier; contraction 0.036 in the source window, but the
    late-cooling contraction rule fails
25. Persists ~τ ≈ 7 kyr; gone by ~160 kyr
26. 8 in the source (knots 94→101) and 29 in the control (93→65)
27. High for the hierarchy gates (structural, biased); low for R20 and the source final state
28. **Yes**
29. **Yes** (same rule)
30. **B**, with frozen-cell segments; C noted as a separate future item
31. A minimal generic interface — Cstar is the only thing that populates it
32. Not ready as written
33. Amendment to ADR-0017 (re-ratified); Phase-5D implication routed separately
34. V0–V10 above
35. V0, V1, V2, V4–V8, V10, plus the owner's BA12R-applicability decision
36. BA12: FAIL preserved
37. BA12R: FAIL preserved
38. Current campaign: FAIL preserved
39. No physical BNV model or rate selected
40. BLOCKING: 2
41. MATERIAL: 7
42. NONBLOCKING: 4
43. NOTE: 5
44. **Disposition B**
45. Next action: owner declines to ratify the proposal as written; authorizes one bounded,
    predeclared qualification task that revises the amendment to A1–A8 and implements it in the
    Phase-6 integration layer on a branch; records the BA12R-applicability decision first; runs
    V1–V10; then returns for ratification. The campaign stays stopped.
