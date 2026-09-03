# ADR-0009 — Ordinary-star TOV surface event and integration-termination contract

| Field | Value |
|---|---|
| **Status** | **ACCEPTED — 2026-09-03** (owner adjudication in TOV-SURF-I); **SOURCE CONFORMED / NUMERICALLY VALIDATED** (TOV-SURF-I-R); artifact migration required |
| **Date** | 2026-09-03 |
| **Change class** | scientific-semantic (definition of `R_*`; INV-06) **and** numerical-method (event location; right-hand-side domain contract) **and** structural (completion status; fail-closed publication) — the strictest class governs |
| **Governing authority** | ADR-0005 §7.1 (the primitive owns "the termination test"; §7.3 froze `SolveToProfile`'s fallback for Phase 3E only); ADR-0003; ADR-0004 (`EvaluateNu` surface BC); ADR-0007 P7 / ADR-0008 Q7 (surface = EOS-floor node, terminal atom); INV-06, INV-13; TOV-RR-01 |
| **Affected invariants** | INV-06 (locator amended, value unchanged); INV-13 (partition is sampling only); INV-08/INV-09 (dependency only) |
| **Blocks** | regeneration of every durable artifact carrying `R_*` until implementation and validation; the corrected Phase-4D revalidation and the first monopole baseline (ADR-0008 Q12) |
| **Evidence record** | `docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md`; defect: `docs/validation/TOV_RADIAL_RES_2500_AUDIT.md` |

> **Relation to prior ADRs.** Amends ADR-0005 only in what "the termination test" of the canonical
> primitive means and in reopening `SolveToProfile`'s fallback policy (frozen by §7.3 for Phase 3E
> only); amends INV-06's *locator* while keeping its *value*; leaves ADR-0007/ADR-0008's surface
> semantics (`R_*` = EOS-floor surface, terminal atom) intact. No accepted Decision is rewritten.

## Context

**What the code does** (`CompactStar/Core/src/TOVSolver.cpp`, authenticated at `e2fe0ad`):
`TOVSolver::ODE` returns `GSL_EBADFUNC` when a **trial** state has `y[1] < PressureCutoff()`;
`gsl_odeiv2_driver_apply` treats that as a fatal user-function error; the loop in
`SingleStarSolveToTOVPoints` `break`s **without storing** the failing target, so the accepted-state
test `if (y[1] <= p_cut) break;` never fires (0/87 audited solves), `R_*` is the last output target
the driver happened to reach (one step short, scattering 33 m / 0.24 % with `radial_res`), and at a
near-discontinuous EOS feature the same guard ends the star 0.56 km early for scattered
`(radial_res, ε_c)` pairs — 4 of 29 resolutions at 1.6 M☉, and ≈3 % of central densities near 2 M☉
at the **default 10000**. `SolveToProfile` bisects on `tmp.back().m` of any returned profile and
falls back to the closest coarse sample on bracket failure, so a truncated star is silently
compensated by a shifted `ε_c` (`grid_convergence_cmf` `B_fixed_mass`, +0.22 %). Two earlier records
observed the symptom (`TOV_REFERENCE.md` §5; `GRID_CONVERGENCE.md` §7's "where the driver first
fails") and deferred it.

**What was measured** (evidence record §5–§8, production physics without the fatal guard):

| question | measurement |
|---|---|
| bracketed event locators (C bisection, D pressure-coordinate landing) across 8 output partitions | `R_* = 13.473358` km at **every** partition; spread `5.6e-11` / `7.8e-11`; `M` spread `1.1e-11` |
| today's node-quantized locator (A), made non-fatal | spread `1.6e-3` (21 m), always overshooting |
| interpolation (B) | `1.4e-4` (linear), `6.8e-4` (log) |
| driver reset per segment vs inherited `h`, under the corrected contract | `≤ 4e-11` |
| trial continuation `ε(max(p,p_cut))` vs `ε(max(p,p_floor))` at a raised cutoff | identical to 10 digits |
| the five `radial_res = 10000` densities that truncated | all reach the event; `R` +0.30–0.34 km, `M` +`6e-4` |
| `dM/dp_c` over partitions / step sizes | `3e-6` / `5e-7` (Richardson pair) |
| shift of the four canonical stars' `R_*` | +1.2 / +1.0 / +5.0 / +3.7 m; `M` `≤ 7e-10`; `z` `1e-4…5e-4` |

## Decision

The owner accepted Q1–Q14 on 2026-09-03 in the TOV-SURF-I implementation request.
This acceptance precedes source mutation. The following decisions are authoritative;
the questions and scratch measurements below retain their historical evidentiary role.

| Question | Accepted decision |
|---|---|
| Q1 | `R_*` is the unique crossing of the accepted continuous TOV solution, `p(R_*) = p_cut`. |
| Q2 | Keep `PressureCutoff() = max(1e-15 p_c, eos_tab.pre[0])`; the cutoff value is not reopened. |
| Q3 | A trial pressure at or below the cutoff is not a fatal RHS condition. The surface is an accepted-solution integration event. |
| Q4 | For internal RK trial states only, evaluate `ε(max(p, p_floor))` through the existing EOS clamp while retaining trial `p` in the pressure terms. No broader EOS continuation or sub-cutoff publication. |
| Q5 | Canonical locator: pressure-coordinate terminal integration from the last accepted state above the cutoff to exactly `p_cut`. Independent test-side cross-check: bracketed radial re-integration/root refinement. |
| Q6 | The radial target partition is sampling only and does not define the star. |
| Q7 | No canonical driver reset per output segment. If post-fix evidence requires it, stop for owner review. |
| Q8 | Fail closed before `SURFACE_REACHED`; no partial star is authoritative. |
| Q9 | Explicit reset-per-solve completion status: `SURFACE_REACHED`, `GSL_FAILURE`, `EOS_DOMAIN_FAILURE`, `R_MAX_EXHAUSTED`, `INVALID_INITIAL_STATE`, `PARTITION_INVALID`, or semantic equivalents. Integer success is positive only for `SURFACE_REACHED`; failed output is cleared. |
| Q10 | Append exactly one final node at `R_*`, with `p=p_cut`, event mass, EOS energy/baryon/species densities and authoritative derivative at `p_cut`, TOV `ν′`, and initial `ν=0` for later reconstruction. Nothing below the cutoff is stored. |
| Q11 | Preserve all seven durable artifacts as historical evidence. Migration is a separate later task after successful validation. |
| Q12 | Apply evidence §12 V1–V12 and its predeclared scientific bounds; do not substitute artifact regression tolerances. |
| Q13 | Amend INV-06's locator to the accepted-solution crossing stored as the final node; retain its cutoff value and EOS-floor semantics. |
| Q14 | ADR-0008 physics is unchanged. Corrected Phase-4D revalidation and the first monopole baseline remain blocked until TOV implementation, validation, and artifact migration. |

This is the normal scientific-semantic migration path, **not** the GOVERNANCE §3.1
pre-baseline exception: acceptance → source correction → independent validation →
separate artifact migration. Acceptance alone does not claim source conformance.

## Historical implementation disposition — TOV-SURF-I (2026-09-03)

Acceptance remains in force. The candidate compiled, but the 2.0-M☉ target-mass
workflow changed its returned mass by `4.8353e-5` relative, exceeding V7's `1e-9`
bound. Work stopped under the owner's impact-envelope rule. Candidate source was
preserved outside the checkout and restored to the authenticated starting source;
no source conformance or completed V1–V12 validation is claimed. All seven artifacts
remain unchanged. INV-06/INV-13 validation, artifact migration, corrected Phase-4D
revalidation and the first monopole baseline remain blocked. Evidence:
`docs/validation/TOV_SURFACE_IMPLEMENTATION.md:96` and `:181`.

## Decision questions

- **Q1 — Scientific definition of `R_*`.** (A) last requested output node before integration
  failure [today]; (B) first accepted output node with `p ≤ p_cut`; (C) **the unique radius at
  which the accepted continuous TOV solution satisfies `p(R_*) = p_cut`**; (D) other.
- **Q2 — `PressureCutoff()` value.** Keep `p_cut = max(1e-15 p_c, eos_tab.pre[0])` unchanged
  (authenticated: floor-dominated on DS(CMF)-1, `p_c`-dominated on HT68 HW)? The value and the
  floor-to-vacuum question are **not** reopened.
- **Q3 — May the right-hand side treat a trial `p ≤ p_cut` as fatal?** (A) yes [today]; (B) no —
  the surface is an integration **event** on accepted states, never a user-function error.
- **Q4 — Trial-state continuation below `p_cut`.** Which bounded continuation may the right-hand
  side use for internal stages: `ε(max(p, p_floor))` (the existing `SafeSplineEval` clamp) with `p`
  as is; `ε(max(p, p_cut))`; a frozen state; other — under the requirements that it is trial-only,
  that no sub-cutoff state is published, and that `R_*` is insensitive to it within a predeclared
  bound.
- **Q5 — Event-locating algorithm.** (A) first target below the cutoff; (B) interpolation in the
  bracketing interval; (C) bracketed re-integration / root refinement; (D) pressure-coordinate
  integration of the terminal interval landing on `p_cut` exactly.
- **Q6 — Is the target-radius partition scientific or sampling-only?** (A) part of the definition
  of the solution; (B) an output/sampling device — `M`, `R_*`, `z_surf`, the central-density root
  and the sequence derivative must be invariant under it to the governed tolerance.
- **Q7 — Driver adaptive-history / reset policy.** (A) reset per output segment as canonical;
  (B) optional stabilization; (C) unnecessary after correct event handling; (D) forbidden as a
  contract element because it couples integration to sampling.
- **Q8 — GSL failure / partial-profile policy.** (A) publish the points before failure [today];
  (B) publish with a warning/status; (C) **fail closed** — no authoritative profile unless the
  governed event was reached — propagated through `SingleStarSolveToTOVPoints`, `SolveToProfile`
  (bisection and coarse fallback), `Solve`, sequence derivatives and test builders.
- **Q9 — Completion status.** An explicit internal solve status
  `{SURFACE_REACHED, GSL_FAILURE, EOS_DOMAIN_FAILURE, R_MAX_EXHAUSTED, INVALID_INITIAL_STATE,
  PARTITION_INVALID}` that callers can query, so a valid complete star is distinguishable from a
  non-empty partial vector; the `int` return of the primitive stays `> 0` only for
  `SURFACE_REACHED`. No large public redesign.
- **Q10 — Surface-point fields.** Store a final `TOVPoint` **at** `R_*` with `p = p_cut` exactly,
  `ε = ε_EOS(p_cut)` (the floor row on floor-dominated tables), `m = m(R_*)`, `ν'` from the TOV
  relation at that state, `ρ`, species and `dε/dp` from the EOS authority evaluated at `p_cut`
  (inside the domain, so the derivative is defined); preceding nodes unchanged; nothing below the
  cutoff stored. INV-06's "surface quantities are the last grid point" then holds by construction.
- **Q11 — Artifact migration / rebaseline policy.** Old artifacts remain historical evidence;
  after validation, class-A/B/D artifacts regenerated once in a dedicated commit with every
  old→new delta recorded against the impact map; `grid_convergence_cmf`'s `2500` rows replaced by
  the corrected star and the resolution set extended (evidence §13).
- **Q12 — Validation sweep and bounds.** Evidence §12: V1–V12 with the predeclared bounds
  (`p(R_*)/p_cut − 1 ≤ 1e-8`; `M` partition invariance `≤ 1e-9`; `R_*` `≤ 1e-8`; `dM/dp_c` `≤ 1e-5`);
  sweep `3e14 … 1.6e15 g cm⁻³` (≥150 densities plus the five that failed) × `{2500, 5000, 10000,
  20000, 40000}` × three uniform grids × both governed EOSs; regression tolerances set *after*
  regeneration from same-build repeatability, never from the scientific bounds.
- **Q13 — INV-06.** Amend the invariant's **locator** ("the last grid point the driver reached" →
  "the crossing `p(R_*) = p_cut` of the accepted solution, stored as the final node") while keeping
  its **value** and the EOS-floor semantics of ADR-0007 P7 / ADR-0008 Q7.
- **Q14 — Phase 4D / ADR-0008.** ADR-0008's physics is not reopened; the corrected Phase-4D
  revalidation and the first monopole baseline **wait** for this ADR's implementation and
  validation (evidence §14), since `R_*`, `I²/R_*³` and the terminal atom's `ε_*` move.

## Alternatives

### Alternative A — keep the fatal trial-state guard, add a driver reset per output segment (audit R2)

- **Statement** — leave `TOVSolver::ODE` as is; reallocate the driver at each output target.
- **Required code changes** — one line in the loop.
- **Migration risk** — none measurable on normal stars (bit-identical in the audit), which is
  precisely the problem: it masks the guard without correcting the surface, `R_*` stays
  node-quantized and one step short, and a differently placed stiff feature can still kill a
  trial stage.
- **Validation needed** — would fail V3 and V4 by construction.
- **Implications** — rejected as the contract; admissible later only as an optional stabilization
  under Q7 if invariance tests demand it.

### Alternative B — non-fatal right-hand side, first accepted node below the cutoff as `R_*` (locator A)

- **Statement** — remove the guard; keep today's loop semantics.
- **Required code changes** — the guard, the store-then-test order.
- **Migration risk** — every star's `R_*` moves *up* by up to one output step and still depends
  on the partition (`1.6e-3`); artifacts move without the partition dependence being cured.
- **Validation needed** — fails V4.
- **Implications** — rejected: the partition would still define the star.

### Alternative C — non-fatal right-hand side + deterministic bracketed event locator + fail-closed status (recommended)

- **Statement** — Q1(C), Q3(B), Q4 narrowest continuation, Q5(D) with (C) as cross-check, Q6(B),
  Q7(C), Q8(C), Q9, Q10.
- **Required code changes** (for the later repair increment, not this ADR) — `TOVSolver::ODE`
  loses the fatal guard (EOS evaluated at `max(p, p_floor)`, which `SafeSplineEval` already does);
  `SingleStarSolveToTOVPoints` detects the event on accepted states, locates the crossing by
  pressure-coordinate integration of the terminal interval (or bracketed re-integration), appends
  the surface point, and returns an explicit status; `SolveToProfile`, `Solve` and the test builders
  refuse any profile without `SURFACE_REACHED`; `TOVSolver.hpp` gains the status type.
- **Migration risk** — every durable artifact carrying `R_*` moves by the impact map (+1–5 m at the
  canonical masses; `z` `1e-4…5e-4`; `M` `≤ 7e-10`); the two `grid_convergence_cmf` artifacts change
  materially in their `2500` rows; `hartle_I_dscmf1_debug.tsv` is no longer bitwise.
- **Validation needed** — V1–V12 (evidence §12).
- **Implications for existing outputs** — all seven artifacts are regenerated once, after
  validation, with recorded deltas; historical files kept as evidence.

### Alternative D — event-aware integration that never evaluates beyond the domain (Q5 option C only, no continuation)

- **Statement** — bracket the crossing purely by re-integration from the last accepted state with
  an ever-shrinking target, never letting a stage cross `p_cut`.
- **Required code changes** — as C, plus a bracketing loop that rejects any stage below the cutoff.
- **Migration risk** — same outputs as C (measured agreement `< 1e-10`), higher cost (≈25 short
  integrations), and the final probe still evaluates at the boundary.
- **Validation needed** — as C.
- **Implications** — admissible equivalent; recommended only as the cross-check of C's locator.

## Consequences

- `R_*` becomes a property of the star (`p(R_*) = p_cut` on the accepted solution), not of the
  output grid; INV-06's locator is amended; INV-13 gains "the radial output partition is a sampling
  device".
- The right-hand side may never signal the surface as a fatal error; the surface is an event.
- No orchestrator may publish, bisect on, or differentiate a profile that did not reach the
  event; `SolveToProfile`'s coarse fallback is reopened and constrained.
- Every durable artifact carrying `R_*` is scheduled for one governed regeneration; the
  `grid_convergence_cmf` refinement set is re-decided on the corrected solver.
- The corrected Phase-4D revalidation and the first monopole baseline are sequenced after this
  ADR's implementation and validation.

## Validation

Evidence record §12, V1–V12, with the predeclared bounds, the sweep, and the separation of
scientific from regression tolerances. Baseline migration per evidence §11: regenerate once, after
V1–V8 pass, in a dedicated commit with recorded deltas; any delta outside the impact map is a STOP.

## Historical recommendations (accepted above)

| # | recommendation | basis |
|---|---|---|
| Q1 | **C** — `R_*` is the unique crossing `p(R_*) = p_cut` of the accepted solution | locators C/D agree to `8e-11` across 8 partitions; today's spread `1.6e-3` |
| Q2 | **unchanged** `p_cut = max(1e-15 p_c, p_floor)` | authenticated on both EOSs; not the defect |
| Q3 | **B — no fatal guard**; the surface is an event on accepted states | the guard is the root cause (audit §17); with it removed every tested star reaches the event |
| Q4 | `ε(max(p, p_floor))`, `p` as is, trial-only, never published | identical to `ε(max(p, p_cut))` to 10 digits; it is the existing `SafeSplineEval` clamp |
| Q5 | **D** (pressure-coordinate landing) canonical, **C** (bracketed re-integration) as cross-check; A and B rejected | D/C `< 1e-10` apart and partition-invariant; A `1.6e-3`, B `1.4e-4` |
| Q6 | **B — sampling only**; `M`, `R_*`, `z_surf`, `ε_c` root, `dM/dp_c` invariant to the governed tolerance | measured `1.1e-11`, `7.8e-11`, `3e-6` |
| Q7 | **C — unnecessary**; not a contract element; optional stabilization only if a post-fix invariance test demands it | reset vs inherited `≤ 4e-11` under the corrected contract |
| Q8 | **C — fail closed** through every orchestrator; no `ε_c` compensation | the `B_fixed_mass` +0.22 % `ε_c`; the factor-4 derivative |
| Q9 | explicit status enum, queried by callers; `int` return `> 0` only for `SURFACE_REACHED` | minimal, no public redesign |
| Q10 | store the surface point at `R_*` with `p = p_cut`, `ε(p_cut)`, `m(R_*)`, `ν'`, `ρ`, species, `dε/dp(p_cut)` | INV-06's "last grid point" then holds by construction |
| Q11 | regenerate affected artifacts once, after validation, with recorded deltas; old files kept | evidence §11 |
| Q12 | the §12 sweep and bounds; regression tolerances from repeatability after regeneration | measured margins of ≥ 2 decades |
| Q13 | amend INV-06's locator, keep its value and the floor semantics | ADR-0007 P7 / ADR-0008 Q7 untouched |
| Q14 | ADR-0008 not reopened; corrected Phase-4D revalidation and first monopole baseline wait for this ADR | `R_*` +5 m at 1.6 M☉; shell `ε_*` changes |

## Provenance

Drafted by the AI agent (TOV-SURF-G, 2026-09-03) from the TOV-RR-01 audit and the scratch
measurements recorded in `docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md`. The draft changed no production, test, CMake or baseline file. The owner accepted Q1–Q14
on 2026-09-03 in TOV-SURF-I; this acceptance commit changes only this ADR and precedes
implementation. No artifact has been regenerated.

## Post-acceptance owner clarification — V7 (TOV-SURF-I-R, 2026-09-03)

This validation clarification leaves the accepted scientific Decision Q1–Q14 above
unchanged. The owner splits V7 into two experiments:

- **V7a — fixed central density:** compare old and corrected solves at identical
  inherited εc. The existing `|ΔM/M| <= 1e-9` and radius impact map apply here.
  The four inherited εc values are `454550405078491.75`, `616488270506054.5`,
  `731253342677476.12`, and `1298349261929558.8 g/cm³`.
- **V7b — target-mass contract:** independently solve targets 1.0, 1.4, 1.6, 2.0 M☉.
  Only complete `SURFACE_REACHED` stars may participate; retain stable-branch
  selection, exclude failures, use no nearest-sample fallback, reach `p=p_cut`,
  and reproduce identical results on repeated calls. Require the unchanged
  absolute residual `|M-returned - M-target| < 1e-4 M☉`. εc may move.

The `1e-9` bound belongs only to V7a; it is not widened or applied to two independently
solved target-mass results. Coarse `N=24`, stable-branch criterion, first-acceptable
bisection stopping, and production `mass_tol=1e-4 M☉` remain unchanged. Commit
`96c1425` retains the valid historical STOP caused by applying the fixed-density
impact bound to independent target-mass workflows. The resume diagnosis and exact
traces are recorded in `docs/validation/TOV_SURFACE_IMPLEMENTATION.md:262`.

## Implementation validation outcome — TOV-SURF-I-R (2026-09-03)

The exact restored candidate conforms to Q1–Q14. V7a and V7b pass under the owner
clarification; both-EOS surface/partition sweeps, sequence derivatives, downstream
surface checks and V11/V12 detectors pass. No production reset; no PressureCutoff,
mass tolerance, coarse N or stable-branch change. The historical stop remains in Git
and the evidence record. All seven original artifact hashes remain unchanged.
**TOV SURFACE EVENT IMPLEMENTED AND VALIDATED — ARTIFACT MIGRATION REQUIRED.**

Evidence: `docs/validation/TOV_SURFACE_IMPLEMENTATION.md:339`. Corrected Phase-4D revalidation and the first monopole
baseline must use the migrated artifacts and remain subsequent tasks. This outcome
does not reopen ADR-0008 physics or authorize Phase 5.
