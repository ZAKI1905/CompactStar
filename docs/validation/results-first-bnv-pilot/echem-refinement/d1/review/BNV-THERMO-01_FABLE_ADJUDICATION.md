# BNV-THERMO-01 — Coupled-roundoff adjudication

Reviewer: Fable 5.2. Date: 2026-11-09. Mode: read-only scientific review.
No agents, cluster jobs, trajectories, tolerance changes or repository edits were made.

## 0. Authentication

| Item | Result |
|---|---|
| Branch | `candidate/results-first-bnv-pilot-20261108` on `origin` (pushed during this review) |
| HEAD | `ad4e91aee3837b8b07e7ed723afdd40b765011fe` — matches expected |
| Review copy | fresh single-branch clone, `/Volumes/Scriptorium/Research/CompactStar/worktrees/bnv-thermo-01-review-20261109`, clean tree, connectivity check clean |
| Committed bundle | `echem-refinement/SHA256SUMS`: 47 of 47 files OK |
| Production source | `git diff 56dd195 HEAD -- CompactStar/` is empty; O2 default is still `1e-13 / (1e-18, 1e-24, 1e-24)` at `PassiveCheckpointOutput.cpp:375` |
| Not verified | The 354 external files under `/mnt/sdd/.../echem-refinement-20261109-r1` (cluster unreachable from this Mac). Everything below marked "recomputed" uses committed files only. |

## A. Scientific verdict

**Conditionally resolvable (outcome B). I do not certify the coupled-roundoff qualification on the present evidence.**

The 73.31 % figure is a budget of *declared tolerance constants*, not of measured error. It shows that the Echem gate passes **if** the true error of the reconstructed chemical state lies inside the declared state box. The evidence that it does is strong (agreement at the few-ulp level against a box of 1600+ ulp) but it is measured agreement plus one same-author, never-validated conditional audit. Under the study's own predeclaration, that is UNQUALIFIED.

One bounded offline diagnostic (section D) converts the premise into a direct measurement. No new solver, no relaxed criterion and no physics correction is needed.

## B. Evidence assessment

### B1. Astra's six status statements

| Statement | Verdict | How checked |
|---|---|---|
| Tested tolerance 5e-27 MeV | Confirmed | It is `O2.absolute[1..2]` of the local rk8pd reconstruction, injected test-only by `const_cast` in `tests/diagnostics/checkpoint_refinement.cpp` |
| Maximum utilization 73.31 % | Confirmed, recomputed | Rebuilt from the code formula and committed states: 0.7327 linearized vs 0.73308 committed (difference is binary64 corner subtraction) |
| Identical reconstructed states | Confirmed, recomputed | `candidates.tsv`: one distinct `(x, eta_e, eta_mu)` per checkpoint across all six atols, 13 of 13 points |
| No production tolerance change | Confirmed | Empty library diff; default unchanged |
| No Stage A retry | Consistent with records | Jobs 19177–19181 are two builds, one failed harness, one study; cannot be independently confirmed without cluster access |
| Coupled roundoff unresolved | Confirmed | See B3 |

### B2. What the tolerance and the budget actually are

The tested number enters the computation in two places:

```text
(a) GSL step acceptance:   |yerr_i| <= atol_i + rtol*|y_i|        (a_y=1, a_dydt=0)
(b) declared error floor:  F_O,i = max(atol_i + rtol*M_i, 64 ulp)
                           r_i   = 2*max(|O2_i-O1_i|, F_O,i) + U_event,i
                           P     = 2*max over 8 corners |Echem(O2 +/- r) - Echem(O2)|
                           U     = 2*max(|Echem(O2)-Echem(O1)|, 64 ulp) + P   <= 0.20 F,  F = 1e-11*|Echem|
```

Role (a) is inactive (B4). Only role (b) changes the result. Decomposition of the 73.31 % at the worst point (index 1, t = 16,436.25 s), recomputed:

| Contribution to U / (0.20 F) | Share | Nature |
|---|---:|---|
| Relative floor, `2*rtol*|eta|` with rtol = 1e-13 | 40.00 % | declared constant |
| Absolute floor, `2*atol` with atol = 5e-27 | 30.61 % | declared constant (the tested quantity) |
| `U_event` (witness discrepancy + its 64 ulp) | 1.67 % | mostly a fixed ulp floor; measured part ~2 ulp |
| Explicit `2*64 ulp(Echem)` | 0.99 % | fixed floor |
| Total | 73.3 % | |

At every point the utilization is affine in atol (recomputed, 13 of 13):

```text
utilization(t) ~= 0.427 + 0.306 * (atol / 5e-27) * (16436 s / t)
```

So the maximum sits at the earliest checkpoint and falls as 1/t. At index 1 the gate would pass up to atol ~ 9.4e-27, and the study's 0.80 margin up to ~6.1e-27. The value 5e-27 was selected by the predeclared rule (largest listed value meeting 0.80), not after the fact.

### B3. Major numerical claims

| Claim | Source | What it establishes | What it does not |
|---|---|---|---|
| 73.31 % maximum utilization | `candidates.tsv`; recomputed | The declared box at 5e-27 propagates to an Echem uncertainty inside the gate at 13 sampled points | That the true error is inside the box |
| O1/O2 Echem disagreement = 0 | `state.tsv`, audit records | Nothing about accuracy: O1 and O2 are the same single rk8pd step (B4) | Two-level convergence |
| Witness agreement <= 2.0e-15 relative Echem | `witnesses.tsv`; recomputed: eta differs by <= 10 ulp among O2 and six witnesses, x by <= 4 ulp | Two methods and three step counts agree at the arithmetic floor; truncation differences are not visible | A bound; all share the left state, the binary64 RHS and GSL |
| Embedded RK error estimates | `witness-steps.tsv` (744 rows); recomputed: <= 1.6 ulp(eta), <= 4.2 ulp(x), zero rejections | The estimators return rounding noise; step control never binds | A truncation bound |
| 100-digit Echem algebra, <= 2.9e-16 relative | `assessment.json` | The quadratic form is well conditioned; no cancellation in the energy algebra | Anything about the state |
| Conditional audit, O2 enclosure <= 5.9e-16 relative Echem (<= 2.7 ulp in eta) | `conditional-roundoff.json`, `audit_chemical_roundoff.py` | A running a-posteriori bound, in exact rational arithmetic, on chemical RK arithmetic plus chemical RHS evaluation. I reviewed the recursion, tableau extraction and context classification and found them structurally correct. | It excludes thermal coupling, truncation and the main prefix; it was never validated on a known case; I could not re-execute it (inputs are external) |
| Audit falsification test (new, this review) | committed files only | 546 pairwise comparisons (13 points x 21 pairs x 2 components): every observed difference <= sum of the two enclosures; worst ratio 0.954 | The enclosures are unfalsified and non-vacuous, not proven |
| Five global references agree to 1.1e-15 relative Echem | Astra report / prior diagnosis | Main-prefix error is small as measured | Certified reference uncertainty; I did not recompute it |

### B4. The identical-state result

Committed audit records show `accepted_steps = 1` for all 156 O1/O2 contexts, with zero rejections. Each reconstruction is one rk8pd step across the whole local span. GSL accepts it at every atol, including 1e-29 where the threshold is essentially `rtol*|eta|` (>= 5.7e-27 MeV) alone; the recorded witness estimates are ~1e-29 MeV. (The O2 step's own estimate is not in the committed files.)

| Candidate explanation | Status | Evidence |
|---|---|---|
| Stopping criterion no longer affects the returned state | **Established** | One accepted attempt and zero rejections at every atol; witness estimates 2–3 orders below the threshold |
| Finite-precision floor of the error estimator | **Established** | Estimates <= 1.6 ulp |
| Dependence between nominally separate checks | **Established** | O1 and O2 are the identical computation, so `d_O = 0` by construction |
| Output rounding concealing differences | Excluded | 17 significant digits round-trip binary64 exactly |
| Genuine convergence under refinement | Not shown by this comparison | No refinement occurred; the tolerance sweep varied nothing in the arithmetic |

The identical states therefore justify exactly one conclusion: the returned state does not depend on atol over [1e-29, 1e-24]. They carry no information about its error.

## C. Mathematical explanation

### C1. Origin of the difficulty

Over the pilot year the chemical state is a driven, almost undamped quadrature:

```text
d(eta_i)/dt = s_i - sum_j Z_ij R_j(T, eta_j),     s = -Z sigma  (constant)
|beta_i| / |s_i| <= 3.6e-5 over the year   =>   eta_i(t) = s_i t (1 - O(1e-5))
Echem = (1/2) eta^T Z^-1 eta ~= (1/2) t * sum_i sigma_i |eta_i|   (quadratic, grows as t^2)
grad Echem = Z^-1 eta ~= -sigma t
```

The state allocation for eta has an absolute floor (`F_i = 1e-22 MeV + 1e-11*|eta|`), but the Echem allocation is purely relative (`1e-11*|Echem|`). Propagating a fixed absolute state floor `a` through the quadratic gives a relative energy uncertainty of order `8a/|eta|`, which diverges as eta -> 0. With a = 1e-24 MeV and |eta| ~ 6e-14 MeV this is 61 times the allowance. That is the original failure, and it is a mismatch between two declared allocations, not an integration error. Astra's classification E is correct.

### C2. Actual error structure

- **Truncation.** The non-polynomial content of eta(t) is the 1e-5 relaxation term with a coefficient varying as exp(6x), and 6*Delta x <= 4.4e-4 over the year. For an order-8 step this is tens of orders below one ulp. Consistent with every recorded estimate.
- **Roundoff.** One RK step from an exact binary64 start: a few ulp. The audit says <= 2.7 ulp for O2; witnesses scatter by <= 10 ulp.
- **Thermal-to-chemical coupling** (the term the audit excludes). Each reaction term scales as T^(6-2k), so `|dR/dx| <= 6|R|`:

```text
|delta eta_i| / |eta_i| <= 6 * delta_x * beta_i t / (2|eta_i|) <= 1.1e-4 * delta_x     (worst case, year end)
```

  With the x gate radius (<= 1.7e-17) this is <= 2e-21 relative. A thermal error would have to reach ~2e-10, seven orders above its own gate, to consume 10 % of the smallest chemical box. The "conditional on recorded thermal stages" caveat is therefore closed analytically for the chemical channel.
- **Chemical-to-thermal coupling.** `DeltaPbeta / Pnet` is 2e-26 at the first checkpoint and ~8e-20 at one year, below half an ulp of the thermal power sum. In Stage A the thermal trajectory is insensitive to the chemical state to the last bit.

Against this, the box at 5e-27 is never smaller than `2e-13*|eta|`, about 1600–3400 ulp.

### C3. Why a tighter tolerance did not, and could not, improve accuracy

Lowering atol changed neither the step sequence nor any arithmetic. It only shrank the declared box. That is legitimate if and only if the true error is inside the smaller box. The binary64 refinement test in the predeclaration (contraction <= 0.5) cannot decide this, because every difference is already at the rounding floor and does not contract. Contraction is a proxy for error when the truth is unknown. In exact arithmetic the proxy is unnecessary: a converged high-precision reference measures the error directly.

### C4. Physical cross-check of the chemical-heating ledger

No inconsistency found. From `FrozenControlledBnvRunContext.cpp`:

```text
C T dx/dt        = P_dir_actual + LH - Lnu_full - Lgamma - Lother
P_dir_actual     = P_dir_eq + eta.sigma           (R18 identity)
dEchem/dt        = -eta.R - eta.sigma,   LH = eta.R
=> d/dt (E_eq + Echem + U_th) = -(Lout + Lnu_full + Lgamma + Lother)      (R20), identically
```

The chemical work `eta.sigma` and the reaction heating `eta.R` each appear exactly once with each sign, so there is no double counting. Numerical checks on committed data:

| Check | Result |
|---|---|
| `mu_B * |Bdot|` vs `P_dir` | equal to 16 digits |
| `DeltaLnu / LH`, `DeltaPbeta / LH` | 1.5000 and -0.5000 (correct small-xi modified-Urca limit) |
| `Echem` vs source work `-(1/2) t eta.sigma` | relative deficit equals `-(1/2) beta/s` at all 13 points, 7.9e-9 to 1.5e-5; no constant offset, so no visible Z asymmetry |

Context for the owner: at one year Echem is 2.7e-15 of the thermal energy change, and chemical heating is ~1.5e-19 of the direct BNV power. The disputed gate protects ledger diagnostics, not the Stage A temperature observable. This is not a reason to relax it.

## D. Minimal next action — Codex task BNV-THERMO-01-D1

**Purpose.** Measure the total numerical error of the O2 chemical state against a converged 100-digit reference on the 13 retained local intervals, with the thermal path treated as an explicit band.

**Class.** Offline analysis of retained artifacts, same class as `assess_refinement.py`. Python standard library `decimal`, precision >= 100 digits. No solver binary, no trajectory, no tolerance change, no production edit.

**Inputs (read-only; verify against `evidence/external-SHA256SUMS` first).**
- `$study/study/output/`: `candidates.tsv`, `witnesses.tsv`, `witness-steps.tsv`, `local-rhs.tsv`, and `point-*/atol-5e-27/{checkpoint.json,state.tsv,diagnostics.tsv}` for all 13 points.
- `$inputs/coefficients.tsv`, `$inputs/tangent.tsv`.
- Committed `evidence/conditional-roundoff.json`.
- The real-valued RHS definition (`rhs`, `reaction`, constants) in `audit_chemical_roundoff.py`. Reuse the definition, not the audit recursion.

**Procedure.**
1. *RHS fidelity.* At every row of `local-rhs.tsv`, require `|f_recorded - f_100digit| <= 2 ulp(f_recorded)` for `de` and `dmu`.
2. *Thermal band.* Per point, build `xhat(tau)` on `[t_left, t_obs]` by cubic Hermite interpolation through the accepted rk8pd/8 witness nodes (values from `witness-steps.tsv`, derivatives from `local-rhs.tsv`). Let `rho` be the largest `|x_node - xhat|` over the accepted rk8pd/16 and RKF45/16 nodes and the O2 endpoint. Fix `delta_x = 1e-12` and require `rho <= 1e-13`.
3. *Reference.* Integrate `d eta/d tau = rhs(xhat(tau) + b, eta)` from the exact binary64 left state to `t_obs` for `b in {-delta_x, 0, +delta_x}`, using classical RK4 with rational coefficients (not a GSL tableau) at N = 64, 128, 256, 512 steps.
4. Define `eta_ref = eta_512(b=0)` and `u_ref,i = 2*|eta_512 - eta_256|_i + max_b |eta_512(b) - eta_512(0)|_i`.

**Acceptance (existing criteria only; no new thresholds).**
- C1, reference convergence: successive differences contract by <= 0.5, or are below `1e-60*|eta|`.
- C2, state premise: at every point and for i in {e, mu}: `|eta_O2,i - eta_ref,i| + u_ref,i <= r_i`, with `r_i` the recorded `propagated_state_radius` at atol 5e-27.
- C3, energy: `|Echem_O2 - Echem(eta_ref)| + u_E <= U_Echem` recorded at 5e-27, with `u_E` the four-corner propagation of `u_ref`.
- C4, report only: for O2 and the six witnesses at each point, `|eta_ctx - eta_ref|` in ulp and its ratio to `enclosure_ctx + u_ref`. A ratio above 1 does not fail D1 but voids the audit as a bound and must be stated.

**Required output.** The script under `tools/linux/native_qualification/`; one JSON with per-point `eta_ref` (40 digits), contraction ratios, both parts of `u_ref`, `rho`, `|O2 - ref|` in MeV and ulp, `r_i`, the C2/C3/C4 ratios and input hashes; a short addendum to the echem-refinement report; updated SHA256SUMS.

**Failure conditions.** Any hash mismatch; step 1 fails (offline RHS is not the production RHS — treat as a continuing blocker and investigate before anything else); `rho > 1e-13`; C1, C2 or C3 fails at any point; any change to `delta_x`, N, method, point list or thresholds after seeing results. On failure: stop, report the measured size, component and time. No retry with altered parameters.

**If D1 passes.** The missing accumulated-roundoff clause of `PREDECLARATION.md` is satisfied for the local reconstruction at O2 chemical atol 5e-27 MeV. The verification receipt may be published, and the owner's already-declared single conditional Stage A retry becomes eligible for the owner's go-ahead: source first, then matched zero-source control, original `RUN_PLAN.md` matrix and criteria, with the sole change O2 chemical atol 1e-24 -> 5e-27. Nothing else is authorized: not Stage B, not another tolerance, not a second retry. Stage A's on-output assessments (U_main, accumulated roundoff on actual outputs, reference interpolation, events, R18/R20, restart, matched signal) all remain; the D1 script should be re-applied to the actual Stage A checkpoints to supply the chemical accumulated-roundoff term.

### Necessary and sufficient evidence per outcome

| Outcome | Requires |
|---|---|
| A. Full qualification | An independent bound showing the O2 error is inside `r_i` at every tested point. Not available today. |
| B. Conditional | D1 passing C1–C3. Sufficient for the reconstruction floor. |
| C. Blocker | D1 failing C2/C3 (a real precision obstacle; report its size), or failing step 1 (model-definition discrepancy). |

## E. Final readiness statement

- `COUPLED_ROUNDOFF_QUALIFIED: NO`
- `CHEMICAL_HEATING_PHYSICS_BLOCKER: NO`
- `TESTED_TOLERANCE_SUFFICIENTLY_CERTIFIED: NO` (pending D1; the value itself is accuracy-neutral and was selected by the predeclared rule)
- `STAGE_A_READY: NO` (becomes eligible for the owner's go-ahead when D1 passes and its receipt is published)
- `STAGE_B_READY: NO`
- `MINIMUM_REQUIRED_NEXT_ACTION:` Codex task BNV-THERMO-01-D1 (section D): 100-digit chemical reference with a thermal band on the 13 retained local intervals, judged by C1–C3.
- `REMAINING_UNCERTAINTIES:`
  1. External traces were not accessible; the audit numbers were reviewed and cross-tested but not re-executed by me.
  2. Scope is 13 of 1920 checkpoints. The 1/t law puts the maximum at the earliest checkpoint only while no checkpoint precedes 16,436 s and reconstructions stay single-step; Stage A must confirm on actual outputs.
  3. Main-prefix (from t = 0) error is outside this gate. It is measured at ~1e-15 by five references but certified by nobody; it belongs to Stage A's U_main assessment.
  4. Thermal RHS evaluation roundoff has no high-precision audit. It is immaterial to Echem and to the temperature signal, but it is still owed to Stage A's strict x-channel accounting.
  5. The O1/O2 two-level check is degenerate in this regime and yields no error information; only the half-step and RKF45 witnesses do.
  6. The absolute-floor state allocation and the relative-only Echem allocation are incommensurate at small |eta|. The tension will recur if eta ever passes near zero; that would be an owner policy question, not a tolerance question.
  7. Stage A does not exercise chemical heating in the temperature channel (contribution below one ulp), so it cannot validate that coupling; the Stage A source-minus-control temperature difference is direct deposition only.
