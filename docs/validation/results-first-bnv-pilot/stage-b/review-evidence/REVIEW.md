# CompactStar — Independent scientific and numerical review of Stage B (2,048-year P2 benchmark)

| | |
|---|---|
| Target | `49253200a801ad9b336789dffb4e93d535a8e903` on branch `candidate/results-first-bnv-pilot-20261108`. Its parent is `a4894de`; it is 24 linear commits after the accepted Stage A commit `c90b220`. |
| Documents | `stage-b/REPORT.md`, `PIECEWISE.md`, `CHEMICAL_RESPONSE.md`, `PREDECLARATION.md`, `REPRODUCE.md`, the committed evidence, and the external roots r1/r2/r3 with their sealed manifest |
| Reviewed | 2026-11-11 on the EKU cluster, from a clean export of the commit. Slurm work ran on phyn004. phyn001 and phyn002 stalled in NFS during Python start-up; those jobs were cancelled and their logs are kept. |
| Reviewer | Claude (AI). Under GOVERNANCE §5 this is verification evidence for the owner, not ratification. |
| Constraints kept | No source change, commit, merge, certificate promotion, native solver run or stellar evolution. Nothing outside this directory was written (`outputs/no-write-check.txt`). |

**Benchmark scope.** The reviewed experiment is the controlled mathematical P2 benchmark: full retention, a fixed Ḃ, fixed accepted stellar inputs, and an exact-zero control. It is not a physical BNV decay model or an astrophysical constraint.

**Note on the review brief.** The R20 values quoted in the brief (source 9.034e-5 N, control 6.99e-9 N) are the **Stage A** one-year totals. The Stage B values are 4.794e-5 N and 4.036e-5 N. Section 6 reviews the Stage B ledger.

---

## 1. Disposition: **B — CONDITIONAL PASS**

Stage B is suitable for explicit owner acceptance under the conditions in section 9. **No numerical correction is required.** I found no load-bearing error, no unsupported assumption that carries weight, and no underestimated uncertainty.

Every central number was reproduced from retained data or by an independent calculation:
- the temperature separation;
- the crossing times;
- the native accuracy;
- the chemical-heating interval;
- R18 and R20.

The candidate's full assessment chain also re-executes with identical results.

I stop short of A for three reasons:
1. The event qualification covers only two transversal upward crossings.
2. The report's isolated "≈5 µK chemical heating" figure answers a narrower question than the net chemical thermal response. That response is a cooling of −73 µK, dominated by actual-potential storage. The report gives this channel's power (3.223e23 erg/s) and the stored E_chem, and says the isolated figure excludes it, but it does not give its temperature effect.
3. GOVERNANCE §5 ratification remains the owner's act.

### Answers to the six questions

**A. Is the +15.33 million K source/control difference numerically reliable?**
Yes. Three independent determinations agree:

| Determination | Endpoint ΔT∞ [K] |
|---|---:|
| Candidate's 110-digit Decimal piecewise reference (HP) | 15,325,803.204209456 |
| Reviewer's mpmath piecewise RK4, independent code and step structure | 15,325,803.204209456 (−2.0e-17 K from HP) |
| Five native rk8pd/RKF45 references | HP + 6.2e-9 to 3.6e-8 K, the binary64 quantum of T (1.5e-8 K) |

The HP's own enclosure, 7.47e-9 K, holds. It is very conservative: the HP's actual error is about 6e-25 in scaled state. The signal exceeds the enclosure by a factor of at least 2.05e15 at every one of the 1,920 nonzero times.

This is a fixed-input numerical statement. The declared omitted finite-temperature source term alone has a scale of 6.26 K, nine orders of magnitude above the numerical uncertainty.

**B. Is the chemical-heating contribution resolved, and is the net chemical response heating or cooling?**
The isolated reaction-heating contribution is resolved:
- My direct counterfactual integration gives **+5.1887 µK**. This lies inside the candidate's analytic interval [+4.936, +5.615] µK, at 37% of its width.
- The net reaction term LH − ΔL_ν is **−2.5943 µK**, inside [−2.811, −2.464] µK. The reactions cool on balance, because ΔL_ν/LH = 1.500 at ξ ≈ 5e-5.
- The actual-potential channel (P_dir,act − P_dir,eq, which feeds chemical storage) gives **−70.53 µK**. It is not covered by the candidate's isolated bounds.
- **The net chemical thermal response is cooling: −73.12 µK in total.** The channels superpose linearly to 1.4e-13.

All of these are microkelvin effects on a 1.1e8 K star. They are mathematically resolved and physically negligible within this benchmark.

**C. Does the calculation demonstrate a physical BNV effect?**
No. It demonstrates the correct thermal and chemical response of the declared model to a prescribed mathematical sink:
- fixed Ḃ;
- full-retention deposition μ_B|Ḃ|;
- fixed tangent coefficients t_j;
- frozen background.

No microscopic rate, decay channel, retention physics or astrophysical constraint is tested.

**D. Which numerical assumptions remain conditional?**
See section 9:
- fixed inputs;
- the 110-digit arithmetic with a 1e-80 per-step pad (corroborated independently, not machine-proved);
- horizon- and box-specific constants (exp(LH) = 1.558, which would reach 2 at about 3,200 yr);
- event evidence limited to two transversal upward crossings;
- an R20 gate that cannot test the chemical ledger.

**E. What is needed before the results can appear as a validated controlled benchmark in a paper?**
No new trajectory or campaign is needed:
1. owner ratification;
2. a complete chemical-channel table, including actual-potential storage and the total;
3. separation of numerical, model-form and fixed-input uncertainty;
4. a statement of what R20 can and cannot test;
5. the documentation clarifications D1–D4 in section 9.

**F. What can Stage B establish about longer-term BNV-driven chemical disequilibrium?**
The state at 2,048 yr:
- The imbalances are still in their growth phase, at only 12–15% of their quasi-stationary values.
- The linear relaxation times are 13,000–16,000 yr at T_end.

Structural results from the benchmark's own equations:
- While the imbalance grows or is stationary, the chemical thermal channel is a net heat sink.
- At stationarity it equals −ΔL_ν, about 1.1e-10 of deposition.
- Net reaction heating needs ξ > 4.91, an imbalance of about 48 keV at T_end, which is 9×10⁴ times the 2,048-yr value.
- Chemical heating can rival full-retention deposition only if Σ t_j|η_j| ~ μ_B, which would need an imbalance of about 25 GeV.

So **under P2 full retention, chemical heating can never compete with direct deposition** (section 11).

---

## 2. Authentication and provenance

**35 of 35 checks pass** (`outputs/authentication.json`).

**Commit and scope**
- History is linear from c90b220 to 4925320.
- Every change is a Stage B addition, or one of the two disclosed Stage A documentation edits (REPORT and SCOPE) plus the added OWNER_ACCEPTANCE.
- Nothing changed under `CompactStar/`, `tests/`, `tools/linux/preliminary_pair` or `native_qualification`.
- The Stage A evidence directory and the historical BA12R documents are byte-unchanged.
- The worktree is at the target commit with no tracked change. Its 19 untracked `.pyc` files predate this review.

**Sealed external evidence**
- My own full pass (`scripts/verify_manifest.py`, job 19415) re-hashed all 663,232 manifest paths (518,068 unique files, 3.23 GB counted with links): **0 mismatches**.
- The walked inventory and symlink list equal the manifest and receipt.
- A negative control was detected.
- The manifest SHA-256 is `4eae8a01…25e0`, as reported.

**Committed evidence**
- Every Git display copy equals the declared 17-significant-digit transform of its external original. I re-implemented the transform and compared structurally; there are 30 bindings.
- `manifest-receipt.json` is byte-identical to the sealed receipt.
- The deliverable zip matches its `.sha256` file.

**Executed code**
- The HP references record the hashes of `piecewise.py`, `reference.py`, `capacity.tsv` and `anchors.tsv`; all four equal the committed files.
- `piecewise.py`, `reference.py`, `driver.cpp` and `run.py` are unchanged since the launch commit ceb2049.
- All finalisation script hashes equal the committed files.
- The GSL `rkf45.c` and `rk8pd.c` sources are now hash-bound in preflight, closing Stage A finding F5, and are unchanged since the Stage A review.

**Inputs and run cards**
- `coefficients.tsv` and `tangent.tsv` equal their committed hashes. Z is exactly symmetric (Z01 = Z10 in binary64).
- Main RKF45 tiers use 1e-7/1e-9/1e-11 with atol (1e-12, 1e-18, 1e-18) scaled by 1e-2 per tier.
- References use rk8pd at 1e-14 with cap H, rk8pd at 1e-15 with caps H/2 and H/4, and RKF45 at 1e-15 with caps H/4 and H/8.
- Main caps are ∞, 38,723,807,757 s, 19,361,903,878.5 s and 9,680,951,939.25 s.
- All of these match `NONLINEAR_PLAN.md:72–79` (`tests/bnv/breakpoint_qualification.cpp:6–8`, `tests/native/stage.cpp:107–131`).

**Initial states and source/control identity**
- Both roles start at t = 0 with x = 0, η = (0, 0), T∞ = 1e8 K and B0 = 7.6169065187188748e56. This holds in every native file and in both HP records.
- The source has Ḃ = −2.4136520263641328e37 s⁻¹; the control has Ḃ = 0.
- Every time grid is exactly END·i/1920; the step is the integer 33,661,440 s.
- The roles ran as separate Slurm jobs and snapshots. The control chemistry is exactly zero everywhere.

**Reproducibility of the 58 checks.** I re-executed the committed assessment chain read-only from the clean export: residuals, observation residuals, assess, chemical_response, events, ledger, local_reference, accuracy, coverage, matched, qualification and final_summary. Inputs came through a shadow `OUT` directory that contains links only to input files.
- **58 of 58 checks are true.**
- Of the 83 outputs, 45 are byte-identical, including accuracy, ledger, matched, events, coverage, chemical-response and all CSVs.
- The other 38 differ only in `elapsed_s`, output paths, and the input hashes that follow from those (`outputs/rerun-comparison.json`).

**Historical failures.** BA12 is still FAIL at 1.7255917120989046. BA12R c1/c2/c3 are still all FAIL: state 32.26 / 0.933 / 43.55, ledger 7.6e-9 / 0.933 / 94.83, and matched 21.92 / 4.10 / 41.03. The values are identical to the c90b220 record, and the governing BA12R documents are byte-unchanged.

## 3. Numerical accuracy of both trajectories

**Independent reference.** `scripts/model.py` and `scripts/integrate.py` implement my own model.
- It is written from the frozen binary64 inputs and the declared formulas, in mpmath at 50 digits.
- It uses fixed-step RK4 with n = 16, 32 and 64 substeps per scientific interval.
- No step straddles a knot. A step that would change cell is redone with x as the independent variable, so the solution lands exactly on the real knot.

Results:
- Convergence is cleanly fourth order (observed order 4.0). The n = 32 and n = 64 runs differ by 6.8e-25 (source) and 1.7e-30 (control); the Richardson estimate is 4.5e-26.
- Against the candidate's HP at all 1,921 points the maximum scaled difference is 5.8e-25 (source) and 1.4e-30 (control). That is 4.0e-8 and 2.1e-9 of the HP's own enclosure U.
- The 5.8e-25 equals the HP's expected RK4 truncation at its 1e6 s step. The HP is therefore accurate far beyond its stated enclosure.

**Native main histories, measured** (`scripts/accuracy_check.py`). The norm is |native − truth|/D with the plan's D, against my Richardson trajectory.

| | Maximum measured ratio | Candidate's gated maximum (adds U terms) |
|---|---:|---:|
| 1,921-point scientific grid, all 24 configurations | 0.3052 (source ULTRA, η_μ, endpoint) | 0.3249 |
| Accepted main endpoints | 0.3052 | 0.3123 |
| Restart-tail endpoints | 0.2932 | **0.7051** |

The gated 0.705 sits at t = 34,756,364,997 s, component η_e, second step of the ULTRA quarter-cap tail. It decomposes exactly (`outputs/closeout-check.json`):
- 0.168 is the measured tail error at that record;
- 0.533 is the inherited prefix bound, the restart state's scaled error of 6.2571e-14 (an x error) carried into η by the max-norm;
- about 0.004 is the U terms.

The scientific-grid values are themselves the native rk8pd reconstructions, so this independently checks reconstruction at all 46,104 samples. Native power channels on the grid are within 5.3e-4 of their G_P budgets (`outputs/powers-check.json`).

**Adaptive control is genuinely exercised in Stage B**, unlike Stage A:
- The source uses 21 / 26 / 40 accepted steps with no rejections.
- Endpoint ΔT errors fall 8.0e-3 → 1.7e-4 → 2.3e-6 K per 100× tolerance step, as RKF45 tolerance proportionality predicts (100^0.8 ≈ 40).
- For the source, caps 1–2 never bind, and the quarter cap changes only the BASELINE history (22 steps). REFINED and ULTRA have one history each.
- The control's capped histories differ from its uncapped ones at most tier/cap combinations (17–22 steps).
- Identical cap subsets are, correctly, not claimed as cap convergence.

## 4. Source-minus-control temperature (`outputs/temperature-check.json`)

**Time alignment and representation**
- Source and control times are identical strings at all 1,921 points and equal END·i/1920 exactly, in the HP, the 5 native references and all 12 native pairs.
- The endpoint is exactly 64,629,964,800 s.
- Temperatures are in K. T = 1e8·e^x is evaluated in 120-digit Decimal from the stored x strings.
- The reviewer does no interpolation: every comparison is at stored common times. The native grid values are themselves rk8pd reconstructions, verified in §3.

**HP endpoint enclosure split.** U_Δ = 7.467e-9 K:
- 78% (5.82e-9 K) is the anchored-versus-primary-formula normalization allowance for the two roles;
- 22% (1.64e-9 K) is truncation and prefix (U_x = 1.44e-17 source, 8.8e-22 control).

**Native pairs about truth**

| Pair | Endpoint error [K] | Reported bound [K] | Measured share of bound |
|---|---:|---:|---:|
| BASELINE | 8.021e-3 | 8.022e-3 | 99.99% |
| REFINED | 1.70e-4 | 1.71e-4 | 99.3% |
| ULTRA | 2.27e-6 | 3.36e-6 | 66% |

Every pair lies inside its bound at all times. The largest grid deviation is 8.24e-3 K, which is 5.2e-10 relative.

As in Stage A, the BASELINE and REFINED bounds are essentially the measured deviation from the finest reference. They are valid bounds, not a-priori error estimates.

**Uncertainty classes**

| Class | Size |
|---|---|
| (i) Numerical, fixed inputs | ≤ 7.5e-9 K (HP); ≤ 8e-3 K for the coarsest native pair |
| (ii) Declared omitted model term (finite-temperature source weighting) | ~6.26 K scale |
| (iii) Fixed accepted inputs | Not requalified: C* table, Urca coefficients, tangent, anchors, frozen background with 2.05e-9 depletion. This is scope, not error. |

## 5. C\* crossing qualification (`outputs/crossings-check.json`)

**Crossing times**

| Knot | My time | Years | Position in HP 0.01 s bracket |
|---|---:|---:|---|
| 94 | 21,351,230,838.2063 s | 676.5796777 | Inside, 1.2 ms after its left edge |
| 95 | 58,866,837,228.7729 s | 1,865.3775074 | Inside, 8.3 ms after its left edge |

The times converge to all printed digits between n = 16 and n = 64. The control never leaves cell 93; its lowest x, −0.0146, is well above the lower knot at −0.0250.

**All 20 source histories pass**: 12 mains, 3 restart tails and 5 references.
- No accepted step straddles a knot.
- Every crossing step ends within 1 ulp of the native binary64 knot. I inferred that knot from the ULTRA events (x − residual) because the native knot formula was not reproduced; the inferred value is consistent across all 20 histories. It lies +1.06e-16 (knot 94) and +1.8e-17 (knot 95) from the real knot, which shifts crossing times by ≤ 5e-5 s.
- The native cell changes only at the event, and the next step starts in the new cell.
- Events match the crossing steps.
- Restart tails start after knot 94 and cross only knot 95.
- All 12 endpoint-schedule repeats are byte-identical in accepted states and events. Observation is passive.

**Native event-time errors against truth.** All are early, consistent with the small positive x error.

| History | Knot 94 | Knot 95 |
|---|---:|---:|
| BASELINE | −2.50 s | −40.69 s |
| BASELINE quarter cap | −1.13 s | −4.68 s |
| REFINED | −0.135 s | −0.82 s |
| ULTRA | −0.0018 s | −0.0098 s |
| Restart tails | — | −4.24 / −0.65 / −0.0093 s |
| References | ≤ 7e-5 s | ≤ 7e-5 s |

The reported 40.697 s full-trajectory bound covers the true 40.695 s.

**Native event location, restoration and neighbour guards** (`outputs/closeout-check.json`)
- Native event brackets are 0 to 2 ulp of t wide (≤ 7.6e-6 s), with x residuals of 0 to 2 ulp. The local root is therefore resolved to binary64. The timing error is entirely trajectory error (table above).
- **Stepper and evolve restoration.** In all 27 events (12 mains and 3 tails), the step after the event equals the step suggestion the controller held before the event truncated the step, to ≤ 0.5 ulp of t. The step was resumed, not restarted from a small trial. There are two expected exceptions:
  - the BASELINE quarter cap, where the cap binds;
  - BASELINE and REFINED after knot 95, where the next step is the final step truncated at the horizon end.
- **Exact-endpoint bypass and both nextafter neighbours** are enforced as stop-on-failure requirements in the retained collector (`tests/native/stage.cpp:140–150`), for the main histories. All six tier tasks completed with exit 0 and "Stage B task complete". The individual neighbour reconstructions are not retained as values, so this evidence rests on the guard code and the authenticated completion receipts (§13).

**Active-cell RHS.** I evaluated all 50,000-odd recorded stage RHS values of the five references with my model, using the cell each stage actually used:
- The thermal derivative agrees to ≤ 10.2 ulp. This is the known systematic photon-anchor offset, Stage A F3.
- Stages within 2 ulp of a knot sometimes use the adjacent cell. The effect is ≤ 0.04 ulp.
- In the two RKF45 references, the c = 1 stage of event and event-trial steps overshoots the knot by ≤ 1.3e-13 in x and is evaluated on the old cell's affine extension. The RHS effect is ≤ 82 ulp (2.0e-26 s⁻¹). The candidate's arithmetic audit already counts it: the audited maximum RHS defect is 2.0e-26 for references 3 and 4.

This is the "smooth-extension strip" that PIECEWISE.md requires to be measured. It is bounded and included.

The 12 main histories are RKF45 as well, and their event steps have the same kind of stage overshoot. Their larger tolerances make the strip wider, but no stage traces are retained for the mains. For them the strip is covered only by the measured state error against truth: at most 0.305 of budget at every grid, accepted and restart point (§3), with no change in accuracy across either crossing.

**Continuity**
- The native knots differ from the real knots by +1.06e-16 (knot 94, about 15 ulp) and +1.8e-17 (knot 95). The effect is a timing shift of 4.6e-5 s or less.
- Across both crossings, the native grid accuracy ratio is unchanged (BASELINE 0.0011; ULTRA 0.0084 and 0.0145).
- Native C\* equals my piecewise C(x) to ≤ 8.7e-16 relative at all 1,921 points.
- dT/dt is continuous: 2.20800e-4 → 2.20779e-4 K/s at knot 94.

There is no artificial thermal discontinuity. Stage B genuinely exercises the breakpoint-aware architecture across tiers, caps, restarts and references.

## 6. Chemical heating and cooling (`outputs/trajectory-check.json`)

**Endpoint reconstruction with my model** (all values agree with the candidate's diagnostics to the printed digits):

| Quantity | Value |
|---|---|
| η_e, η_μ | −2.1199137e-7, −5.3112668e-7 MeV |
| ξ_e, ξ_μ | −2.16e-5, −5.41e-5 |
| L_H | 4.1450912065e22 erg/s |
| ΔL_ν | 6.2176368096e22 erg/s |
| ΔL_ν/L_H | 1.5000 (small-ξ limit 22020/14680) |
| L_H − ΔL_ν | −2.0725e22 erg/s |
| P_dir,eq − P_dir,act | 3.2230e23 erg/s (9.597e-12 of P_dir,eq) |
| L_H/P_dir,eq | 1.234e-12 |
| E_chem | 9.8696e33 erg |
| Reaction / source in η̇ | 12.26% (e), 15.06% (μ) |
| Reaction loss / storage input | 12.86% |

**How the candidate isolated the chemical-heating effect.** The candidate defines it as the difference between the full flow y and a flow z whose thermal right-hand side omits only L_H/(TC), with chemistry kept coupled. The difference is bounded analytically: Gronwall gives δ ≤ e^{LH}·I, and variation of constants with damping in [−J_xx, 0] gives the sharper interval.

I checked the argument and its premises:
- The damping sign holds: −0.86·J_xx ≤ ∂f₀ₓ/∂x ≤ −3.8e-12 at 400 sampled points per cell, for both suppressed flows.
- The net direct-minus-cooling power is positive.
- ΔL_ν ≥ L_H everywhere in the box.
- The cross majorant holds.
- exp(LH) = 1.558 < 2, which the lower bound needs.

Then I computed z directly, as a term-suppressed 3-ODE flow integrated by my independent code (a reviewer diagnostic, not a new stellar evolution):

| Thermal channel removed | ΔT at 2,048 yr | Candidate's analytic interval |
|---|---:|---:|
| L_H only | **+5.18865 µK** | [+4.93637, +5.61461] µK, inside |
| L_H − ΔL_ν (net reaction) | **−2.59433 µK** | [−2.81136, −2.46412] µK, inside |
| Actual-potential correction (P_act → P_eq) | **−70.528 µK** | not covered |
| All three chemical channels | **−73.123 µK** | not covered |

- The discretisation error is about 1e-17 K.
- The flow actually damps at ∂f_x/∂x ∈ [−5.35, −4.27]e-12 s⁻¹, against the candidate's majorant J_xx = 6.86e-12 s⁻¹. This explains why the true L_H value lies in the lower half of the interval.

At 1,024 yr the four values are +0.543, −0.272, −19.24 and −19.51 µK. The effects grow roughly quadratically, because η grows linearly.

**Conclusions on chemistry**
1. The candidate's L_H and net-reaction intervals are correct and conservative. Their signs and magnitudes are independently confirmed. The resolution is in the fixed-input mathematical sense.
2. The reaction channel cools on balance: ΔL_ν = 1.5·L_H for ξ ≪ 1.
3. The dominant chemical thermal effect is not reaction heating. It is the actual-potential deficit that the source diverts into chemical storage: ∫(P_eq − P_act) dt = 1.065e34 erg, which is 13.7 × ∫L_H.
4. **Net, chemistry cools the star by 73 µK relative to an equilibrium-deposition flow.**

The report gives the actual-potential power and E_chem and says the isolated bounds exclude that channel, but it does not give its temperature effect. That is condition C6 and correction D1.

**Chemical work is counted exactly once.** The ADR-0015 actual-potential ledger is:
- thermal: C·dT/dt = P_act + L_H − L_ν,full − L_γ;
- chemical storage: dE_chem/dt = (P_eq − P_act) − L_H;
- sum: d(U + E_chem)/dt = P_eq − L_ν,full − L_γ.

L_H appears once as a transfer and ΔL_ν once as a loss. The identity dE_chem/dt = −erg·η·σ − L_H needs Z symmetric, which holds exactly.

Pointwise at all 1,921 points (BASELINE and ULTRA), the native diagnostics satisfy:
- Ė_source·erg = P_eq − P_act to ≤ 2.0e-16 G;
- Ė_reaction·erg = −L_H to 1e-28 G;
- ΔP_β = L_H − ΔL_ν exactly;
- L_ν,full = L_ν,eq + ΔL_ν;
- P_net matches its channels to ≤ 1.7e-16 G;
- P_eq = erg·μ_B·|Ḃ| to 4e-17.

Integrated, ΔE_chem = ∫(P_eq − P_act) − ∫L_H to 8.5e-6 relative. That residual is the binary64 quantum of P_act (4.7e18 erg/s) against a 3.2e23 erg/s deficit, not a ledger error.

## 7. R18, R20 and energy conservation (`outputs/ledger-check.json`)

**Full retention.** L_out,fluid, L_esc,star, J_X, L_other and the three escape-energy rates are exactly zero at every point of all 96 grids, in both roles. No deposited energy escapes, so P_dir,eq = erg·μ_B·|Ḃ| is the whole source term of the ledger.

**R18.** R18/G ≤ 1.96e-16 in every configuration. That is the binary64 resolution of P_act: one ulp of 3.36e34 erg/s is 4.7e18 erg/s, about 1.3e-16 G. Recomputing P_eq − P_act − erg·Ė_source exactly from the stored values gives the same level, ≤ 3.9e-16 G. The direct-energy identity holds to rounding.

**R20, recomputed exactly.** I recomputed R20 in exact 120-digit arithmetic from the native diagnostics on all four grids of all 24 configurations. The original native sums are reproduced to ≤ 4.7e-15 N.

Decomposition of the source residual:

| Term | Size |
|---|---|
| Recorded residual, 1,920 grid | 1.135e-8 N |
| of which baryon-count quantisation | 1.046e-8 N (92%). The recorded B_end is 0.167 ulp off the exact line B0 + Ḃt, as in Stage A F4. |
| Residual with the exact count | 9.0e-10 N (BASELINE), 2.7e-10 N (ULTRA) |
| after removing the h² trapezoid term | 6.6e-10 N (BASELINE), 3.9e-11 N (ULTRA) |

The BASELINE remainder equals C·ΔT_err = 1.48e38 erg/K × 8.0e-3 K = 6.1e-10 N, which is the native temperature error. The control residual is 6.0e-11 N and is pure trapezoid error: nested ratios 0.2500, Richardson remainder −3.5e-13 N.

**Allowance**

| Allowance term | Source | Control |
|---|---:|---:|
| Total / N | 4.794e-5 | 4.036e-5 |
| U_quad / N | 4.786e-5 (99.8%) | 4.036e-5 |
| of which outgoing-power rectangle bound | 4.783e-5 | 4.036e-5 |
| U_count / N | 6.26e-8 (0.13%) | 0 |
| U_main / N | 7.1e-10 | 4.6e-13 |
| U_reconstruction / N | 8.9e-10 | 8.3e-10 |
| U_event / N | 2.6e-11 | 2.6e-11 |
| Limit | 2e-4 | 2e-4 |

**The brief's premise no longer holds for Stage B.** The equilibrium-count term is 0.13% of the source allowance, not a substantial part. U_count is one ulp of B0 times μ_B = 1.21e38 erg. It is added on top of a residual that already contains the count error, which is a conservative double count, not an unsupported cancellation.

The rectangle bound is valid, because trapezoid error ≤ h·(sup − inf) per interval. It is first-order, however:
- 2.4×10⁵ times the actual outgoing trapezoid error for the source;
- 7.2×10⁵ times for the control.

It alone explains the control's jump from 6.99e-9 N (Stage A) to 4.04e-5 N.

**What R20 can and cannot test.** The ledger is conservative, has no double counting that matters and uses no cancellation. It is also uninformative for chemistry:
- ΔE_chem = 5.1e-12 N and ∫L_H = 4.0e-13 N;
- both are below even the true residual (about 1e-9 N), let alone the allowance.

R20 confirms total-energy closure at about 1e-9 N. The chemical ledger rests on the pointwise identities above and on the native–HP–reviewer agreement, not on R20. The report says so; condition C6 keeps it explicit.

## 8. Precision and accumulated-error assessment (ENCLOSURE / PIECEWISE)

**Derivations**
- **Per-cell Cauchy bound** R(h/r)⁵/(1−h/r) on the cell's affine extension. It is valid because every reference step keeps the true tube and all RK stages inside the active cell (distance > U + 2M_real·h).
- **Crossing bridge.** The exact-flow-versus-Euler error is LM_real h²/2 and the RK4-versus-Euler error is LM_real h²/2, so their sum LM_real h² holds. The code's extra factor e^{2Lh} is a further overestimate.
- **Lipschitz continuity across knots** follows from the continuity of C\*. The preflight confirms all knot jumps < 1e-100.
- **Prefix propagation** uses e^{Lh}. The normalization forcing uses b·t·e^{Lt}.

I re-derived all four; they are correct.

**Majorants** (`outputs/constants-check.json`)
- My independently derived per-cell J, L, M_real, complex M and r equal the preflight values to ≤ 1.1e-40 relative in all five cells for both roles.
- At 400 points per cell, including all corners, the actual Jacobian never exceeds its majorant. The chemical rows are attained at the corners, so the majorant is a supremum, as in Stage A, and the 1e-80 pad makes it strict.
- Real speed is ≤ 0.99 of M_real. Complex-polydisc speed is ≤ 0.85 of M.
- exp(L·H) = 1.558 for the source and 1.456 for the control. It would reach 2 at 3,202 yr and 3,776 yr respectively.

**Classification of the numerical claims**

| Class | Claims |
|---|---|
| 1. Established mathematical bounds, given their stated premises | Piecewise HP enclosure (Cauchy + bridge + Gronwall); normalization closure; chemical comparison intervals; R20 rectangle bound; event sign-bracket certification (true crossing inside the 0.01 s bracket, transversal speed > 0) |
| 2. Independently validated numerical estimates | HP accuracy (reviewer 5.8e-25); ΔT to 2e-17 K; native state errors ≤ 0.305 of budget; crossing times; chemical channel values; R18/R20 decomposition; actual-stage arithmetic audits (reproduced by re-execution; the RHS-defect magnitudes agree with my independent evaluation) |
| 3. Explicitly conditional | 110-digit arithmetic with 1e-80 per-step and PI-literal pads (not machine-proved; corroborated by mpmath to 1e-25); fixed inputs; anchored RHS (closed by b; I recompute b = 3.3010e-28 / 2.1982e-28 s⁻¹ for source / control from the anchor-versus-exact-formula differences in `normalization.json` to 4.5e-26 relative, and the exact direct normalization equals erg·μ_B·|Ḃ| exactly); domain box and horizon |
| 4. Unsupported | None found |

**Shared errors between "independent" references.**
- The five native references share the production RHS and its anchors. Their mutual contraction tests discretisation only, not implementation.
- The two genuinely independent implementations are the candidate's Decimal HP and my mpmath model. Both agree with the native code at the ulp level (systematic 9–10 ulp anchor offset).
- All three share the declared model definition by design, and that is what the fixed-input scope means.

**Precision adequacy.** Every scientific conclusion needs only:
- about 1e-3 K for ΔT;
- about 10% of 5 µK (5e-7 K, or 4e-15 in x) for the chemical diagnostics.

The available accuracy (1e-17 K reviewer, 7.5e-9 K HP enclosure) exceeds this by large margins. No further precision work is warranted.

## 9. Conditions for acceptance and documentation corrections

| | Condition |
|---|---|
| C1 | **Fixed inputs.** The C\* table, Urca coefficients, tangent row, photon/direct anchors and frozen background are accepted as given (2.05e-9 depletion, frozen utilisation ≤ 1.1e-3). BA13 and static uncertainty are not requalified. |
| C2 | **Model scope.** P2 full retention is a prescribed mathematical sink. The omitted finite-temperature weighting has a ~6.26 K scale. Numerical precision claims apply to the declared real RHS. |
| C3 | **Arithmetic premise.** The HP enclosure assumes 110-digit evaluation errors below the 1e-80 per-step pad. This is corroborated by an independent mpmath implementation, not formally proved. |
| C4 | **Horizon and box specific.** x ∈ [−0.03, 0.2], \|η_e\| < 4e-7, \|η_μ\| < 8e-7 MeV, exp(LH) < 2. The constants must be regenerated for longer horizons; η_μ already uses 66% of its box. |
| C5 | **Event evidence.** Two transversal upward crossings (knots 94, 95) in the source only. There are no downward, near-tangential or control crossings. |
| C6 | **Chemistry.** The microkelvin results are fixed-input mathematical diagnostics. The isolated L_H and net-reaction bounds exclude actual-potential storage. R20 cannot test the chemical ledger. |
| C7 | **Ratification.** Formal acceptance is the owner's act under GOVERNANCE §5. |

**Documentation corrections.** None changes a number or blocks acceptance.
- **D1.** Add the full chemical-channel table to REPORT.md: L_H, net reaction, actual-potential storage, total. Either extend `chemical_response.py` to the storage forcing (g = erg·η·σ/(TC) has a fixed sign, so the same comparison argument applies), or cite this review's diagnostic values. State the net chemical response as cooling, −73 µK.
- **D2.** PIECEWISE.md calls the native knot "rounded binary64". It sits 15 ulp above the real knot at knot 94. RKF45 event-step stages also overshoot the knot by up to 1.3e-13 on the old cell's extension. Both are bounded and audited, but the text should describe them.
- **D3.** State that R20's allowance is 99.8% a first-order rectangle bound, about 10⁵–10⁶ times the actual trapezoid error. Give the actual residual decomposition: count quantisation 1.05e-8 N; trajectory error ≤ 1e-9 N. Optionally replace the rectangle bound with a curvature bound so the gate becomes informative.
- **D4.** Present the three uncertainty classes of section 4 separately next to the headline ΔT.

## 10. Independent checks performed

| Script | Check |
|---|---|
| `verify_manifest.py` | Full re-hash of 663,232 sealed paths; inventory, symlinks, negative control |
| `authenticate.py` | 35 provenance, scope, binding, code, input and run-card checks |
| `compare_rerun.py` (+ `rerun.sbatch`) | Read-only re-execution of the whole committed assessment chain; 58/58 |
| `model.py`, `integrate.py`, `trajectory_check.py` | Independent piecewise model and integrator; HP comparison; ΔT; crossings; five chemical counterfactual flows |
| `temperature_check.py` | ΔT from HP, 5 references and 12 native pairs; alignment, initial states, enclosure split |
| `accuracy_check.py`, `powers_check.py` | Measured native state and power errors at every grid, accepted and restart point |
| `crossings_check.py` | 20 histories: knots, events, strips, stage-level cell selection, grid continuity, schedule repeats |
| `ledger_check.py` | Exact R20 recomputation and decomposition; R18 and ADR-0015 identities |
| `constants_check.py` | Per-cell majorants, sampling, comparison-bound premises, exp(LH), net-heating threshold |
| `longterm.py` | Relaxation times, quasi-stationary imbalance, stationary energetics, structural ratio |
| `closeout_check.py` | Step restoration after events, native event brackets, exact decomposition of the 0.705 gate, full-retention zeros, recomputed normalization bound b |

## 11. What Stage B establishes about longer-term chemical disequilibrium; next physics milestone

**Established within the benchmark** (`outputs/longterm.json`)
1. At 2,048 yr the imbalances (η_μ = −5.3e-7 MeV, ξ_μ = −5.4e-5) are still in near-linear growth. They are 12.2% and 15.2% of their linear quasi-stationary values η\* = (−1.73, −3.50)e-6 MeV.
2. The relaxation times are 16,170 and 13,196 yr at T_end. They have shortened 2.2× as T rose, since the rates scale as T⁶.
3. The chemistry's thermal effect is (P_act − P_eq) + L_H − ΔL_ν = −dE_chem/dt − ΔL_ν. That is ≤ 0 whenever the imbalance grows or is stationary. **Chemistry is a net heat sink throughout this regime.**
4. At stationarity, L_H equals the actual-potential deficit (7.5e-11 of P_dir) and the net thermal effect is −ΔL_ν (−1.1e-10 of P_dir).
5. Reaction heating exceeds the excess neutrino loss only for ξ > 4.910, about 48 keV at T_end. That is 1.4×10⁴ times η\*.
6. The exact benchmark identity (P_eq − P_act)/P_eq = −Σ t_j η_j/μ_B (t_e = 0.0309, t_μ = 0.0034, μ_B = 868 MeV) bounds every chemical channel by about Σ t_j|η_j|/μ_B. To compete with **full-retention** deposition, the imbalance would need to be about 25 GeV.

Stage B therefore answers the priority question negatively *for P2*: under full retention, chemical heating or cooling never competes with direct deposition. The ratio is ~1e-11 now and ≤ ~1e-10 at stationarity.

**Orientation from the linearised benchmark equations** (not a prediction). For ξ ≪ 1:
- η\* ∝ Ḃ/k(T), with k ∝ T⁶;
- the stationary chemical heating therefore scales as L_H\* ∝ Ḃ²/k(T), while deposition scales as Ḃ.

At T_end, L_H\* is 7.5e-11 of deposition and about 6e-10 of the passive cooling L_γ + L_ν. Parity with cooling would need the benchmark Ḃ scale to be about 4×10⁴ times larger. That would put ξ\* ≈ 15, beyond the net-heating threshold and outside the linear regime the scaling assumes. Wherever chemistry becomes thermally relevant, the calculation must therefore handle the nonlinear (ξ ≳ 5) Urca polynomials and a time-dependent T.

**Next physics milestone (recommendation, not authorised here).** The scientifically interesting regime needs a benchmark in which the deposited thermal energy per destroyed baryon is decoupled from μ_B while the lepton-number driving σ_j is kept. One example is a retained-energy fraction ε, with direct heating ε·μ_B|Ḃ|, including ε → 0 as an invisible-channel control.

The competition threshold is then ε ~ Σ t_j|η_j\*|/μ_B (about 1e-10 at T_end in this background). The net-heating regime additionally needs ξ ≳ 5.

Required steps:
1. An owner-approved ADR defining that parameter. It is a scientific-semantic change, and still not a microscopic BNV rate.
2. An analytic quasi-stationary predictor, using the linear relaxation matrix above, to choose ε and horizons before any run.
3. Regenerated piecewise constants on a larger η box. The ~5×10⁴ yr needed to approach stationarity exceeds both the 3,200-yr exp(LH) limit and the present η box.

The existing Stage B machinery and this review's independent model can be reused for both the enclosure and the cross-check.

## 12. Findings index

| | Finding | Severity |
|---|---|---|
| F1 | ΔT∞ = +15,325,803.204209456 K verified three ways; enclosure 7.47e-9 K valid; native pairs inside bounds | Confirmation |
| F2 | Piecewise enclosure derivation and majorants correct; HP accurate to 6e-25 (4e-8 of U) | Confirmation |
| F3 | Crossings at 676.5796777 / 1,865.3775074 yr; 20 histories clean; strips bounded and audited | Confirmation, plus doc D2 |
| F4 | L_H +5.189 µK and net reaction −2.594 µK inside the candidate's intervals; storage −70.5 µK; total −73.1 µK (cooling) not reported | Moderate (documentation, D1) |
| F5 | R20 conservative; count term 0.13%, not substantial; rectangle bound 10⁵–10⁶× loose; R20 blind to chemistry | Minor (D3) |
| F6 | The 0.705 restart gate is 0.168 measured + 0.533 inherited prefix + about 0.004 U; true native errors ≤ 0.305 everywhere | Confirmation |
| F7 | Assessment chain re-executes: 58/58; 45 byte-identical outputs, 38 metadata-only | Confirmation |
| F8 | The brief's R20 numbers are Stage A's | Note |

No MATERIAL or BLOCKING finding.

## 13. Limits of this review

- The physics of the declared RHS and the accepted stellar inputs were taken as given (C1, C2). Agreement between three implementations shows they compute the same function, not that the function is the right physics.
- My counterfactual flows are a reviewer diagnostic: 3-ODE integrations of the declared model with one thermal channel removed. They are not stellar evolutions or native solver runs.
- Not independently re-derived:
  - the native rk8pd reconstruction radii F_i and the local-enclosure computation (reproduced by re-execution, and covered by my direct measurement at every reconstructed sample);
  - the actual-stage arithmetic recursion for each reference (reproduced; spot-checked through the RHS-defect magnitudes).
- Majorant sampling is a sample of 400 points per cell plus corners; the majorants themselves were re-derived analytically.
- Not re-executed; authenticated by manifest, completion receipts and committed logs only:
  - native trace generation (`supplement.py`, `observation_trace.py`, `native_points.py`), which needs the native binaries;
  - the preflight with its six ctest contracts;
  - the candidate's 13 unit tests (7 piecewise, 5 Apply-residual, 1 audit-source);
  - the native event-neighbour guards, evidenced by guard code plus exit-0 completion.

  My conclusions do not depend on them. The state, power, event and ΔT results come from direct comparison of stored native outputs with my independent trajectory.
- Compute nodes phyn001 and phyn002 stalled during Python start-up (in-kernel NFS waits). Those jobs were cancelled before computing anything, and their logs are in `logs/cancelled-*`. All results come from phyn004 or the login node.

## 14. Files

All scripts run as `toolchain/python/bin/python -I -B scripts/<name>.py <BASE> <REVIEW_DIR>`. The Slurm wrappers are `scripts/*.sbatch`.

`clean-tree/` is the export of 4925320 and is regenerated by `scripts/export_clean_tree.py`. `rerun/` holds the re-executed candidate outputs; its input links point read-only into r3.

`SHA256SUMS` covers every regular file in this directory except `clean-tree/`. That includes the logs of the cancelled jobs and the `rerun/` outputs, but not `rerun/`'s input links. It is verified in `SHA256SUMS.check.txt`.
