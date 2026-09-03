# TOV-RR-01 — audit of the `radial_res = 2500` TOV anomaly

> **FORMAL STATUS: `TOV RADIAL_RES ANOMALY AUDIT COMPLETE — GOVERNANCE DECISION REQUIRED`**
>
> The anomaly reproduces, and it is **not** specific to `radial_res = 2500`. On DS(CMF)-1 a
> scattered set of resolutions — and, at the **production default 10000**, a scattered set of
> central densities — produce a star whose integration **dies at the crust–core transition**,
> six orders of magnitude in pressure above the surface, losing the entire outer crust
> (`ΔM ≈ 2.3e-3 M☉`, `ΔR ≈ 0.56 km`). One durable artifact
> (`grid_convergence_cmf_1p6_debug.tsv` / `_trajectory.tsv`, `radial_res = 2500` rows) was
> computed on such a star.
>
> **Root cause:** `TOVSolver::ODE` returns `GSL_EBADFUNC` whenever a **trial** pressure falls
> below `PressureCutoff()`. A user-function error is **fatal and non-retryable** for
> `gsl_odeiv2_driver_apply`, so one wild internal Runge–Kutta stage at a near-discontinuous EOS
> feature permanently ends the star. **Amplifying factor:** the single driver's step size `h` is
> inherited across target segments, so which `(radial_res, ε_c)` pairs are hit is a fine-grained
> function of the target partition.
>
> **No production change was made and none is proposed here.** Both candidate repairs are
> scientific/numerical-method changes requiring an ADR.

| Field | Value |
|---|---|
| **Starting commit** | `e2fe0adf53d7975d3a5c57e84a8d3481ffd2ce41` (Phase 4D-RI), branch `analysis/tov-radial-res-2500-audit`, worktree `worktrees/CompactStar-tov-rr2500` |
| **master** | `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Change class** | **documentation** — production diff **NONE**; every experiment is scratch, outside the tree |
| **Suites at entry** | full **31/31** (241.72 s), self-contained **17/17** (23.28 s) |
| **Governing authority** | ADR-0005 (`SingleStarSolveToTOVPoints` is the sole ordinary-star radial primitive); INV-06 (surface convention); INV-13; `GOVERNANCE.md` §2 change classes |
| **Scratch** | `scratchpad/rr/targets.py`, `rr_sweep.cpp`, `rr_part.cpp` — not tracked |

---

## 1. Starting commit and authentication

Fresh worktree from `e2fe0ad`, clean, seven durable artifact hashes as `PHASE3_CLOSEOUT.md` §1
(re-verified unchanged at exit, §18).

## 2. Exact reproducer

```
EOS   data/compose/DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos
ec    7.312533e14 g/cm^3
call  TOVSolver::SetRadialRes(N); TOVSolver::SingleStarSolveToTOVPoints(ec, pts)
```

Everything below is **TOV-only** — `M`, `R`, `p(r)`, `ε(r)`, node layout, sequence derivative. No
Hartle quantity was used to establish the anomaly.

## 3. Resolution sweep

| N | n | R [km] | M [M☉] | p_last [dyn/cm²] | termination | spans transition | dM/dp_c [km³] | verdict |
|---|---|---|---|---|---|---|---|---|
| 1000 | 266 | 13.472558 | 1.599975772 | 3.60e25 | GSL-break | yes | 9671.25 | normal |
| 1250 | 331 | 13.442528 | 1.599975757 | 2.71e26 | GSL-break | yes | 9671.25 | normal |
| 1500 | 398 | 13.471508 | 1.599975771 | 3.94e25 | GSL-break | yes | 9671.25 | normal |
| 1750 | 461 | 13.440008 | 1.599975755 | 3.10e26 | GSL-break | yes | 9671.25 | normal |
| **2000** | 529 | 13.462933 | 1.599975769 | 7.72e25 | GSL-break | yes | **20361.81** | **derivative anomalous** |
| 2250 | 594 | 13.447319 | 1.599975761 | 2.09e26 | GSL-break | yes | 9671.25 | normal |
| 2350 | 620 | 13.452221 | 1.599975764 | 1.57e26 | GSL-break | yes | 9671.25 | normal |
| 2400 | 634 | 13.459258 | 1.599975768 | 9.98e25 | GSL-break | yes | 9671.25 | normal |
| 2450 | 646 | 13.457437 | 1.599975767 | 1.13e26 | GSL-break | yes | 9671.25 | normal |
| 2475 | 654 | 13.450331 | 1.599975763 | 1.75e26 | GSL-break | yes | 9671.25 | normal |
| 2490 | 658 | 13.460671 | 1.599975768 | 9.06e25 | GSL-break | yes | 9671.25 | normal |
| 2495 | 659 | 13.461752 | 1.599975769 | 8.40e25 | GSL-break | yes | 9671.25 | normal |
| 2499 | 661 | 13.469616 | 1.599975771 | 4.62e25 | GSL-break | yes | 9671.25 | normal |
| **2500** | 641 | **12.904228** | **1.597681092** | **4.34e31** | GSL-break | **no** | 9734.39 | **TRUNCATED** |
| 2501 | 661 | 13.458845 | 1.599975768 | 1.03e26 | GSL-break | yes | 9671.25 | normal |
| 2505 | 662 | 13.465297 | 1.599975770 | 6.49e25 | GSL-break | yes | 9671.25 | normal |
| **2510** | 643 | **12.908594** | **1.597999359** | **3.74e31** | GSL-break | **no** | **18511.48** | **TRUNCATED** |
| **2525** | 646 | **12.894285** | **1.596917316** | **5.80e31** | GSL-break | **no** | 9736.90 | **TRUNCATED** |
| 2550 | 674 | 13.469243 | 1.599975771 | 4.76e25 | GSL-break | yes | 9671.25 | normal |
| 2600 | 687 | 13.468008 | 1.599975771 | 5.26e25 | GSL-break | yes | 9671.25 | normal |
| 2750 | 726 | 13.451844 | 1.599975764 | 1.60e26 | GSL-break | yes | 9671.25 | normal |
| **3000** | 768 | **12.904391** | **1.597693209** | **4.32e31** | GSL-break | **no** | **19876.92** | **TRUNCATED** |
| 3500 | 922 | 13.460008 | 1.599975768 | 9.49e25 | GSL-break | yes | 9671.25 | normal |
| **4000** | 1055 | 13.461621 | 1.599975769 | 8.48e25 | GSL-break | yes | **20802.71** | **derivative anomalous** |
| 5000 | 1319 | 13.463458 | 1.599975769 | 7.43e25 | GSL-break | yes | 9671.25 | normal |
| 7500 | 1977 | 13.470201 | 1.599975771 | 4.40e25 | GSL-break | (node inside) | 9671.25 | normal |
| 10000 | 2635 | 13.468323 | 1.599975771 | 5.13e25 | GSL-break | yes | 9671.25 | normal |
| 20000 | 5268 | 13.471018 | 1.599975772 | 4.11e25 | GSL-break | yes | 9671.25 | normal |
| 40000 | 10535 | 13.472681 | 1.599975772 | 3.56e25 | GSL-break | yes | 9671.24 | normal |

**Answers to the first question.** `2500` is **not** unique and **not** a narrow band: the
truncated set sampled here is `{2500, 2510, 2525, 3000}` — scattered, with `2499`, `2501` and
`2505` all normal. A second, distinct symptom (a correct central star but a corrupt sequence
derivative, because only one of the two perturbed stars truncates) appears at `{2000, 4000, 2510,
3000}`. The `7500` "no" is a false flag of the strict bracket test — a node landed *inside* the
transition window, so no single interval spans both edges; the star is normal.

**Three observations that hold at every resolution.**

1. **No resolution ever terminates through the governed cutoff test.** Every run ends with
   `gsl_odeiv2_driver_apply` returning status **9 = `GSL_EBADFUNC`** (87/87 solves). The loop's
   own `if (y[1] <= p_cut) break;` never fires on this EOS.
2. `p_cut = max(1e-15 p_c, eos_tab.pre[0]) = 3.351885e25` dyn/cm² — the EOS table floor, not the
   `1e-15 p_c` term (`p_c = 9.31e34`, so `1e-15 p_c = 9.3e19`).
3. Normal runs stop at `p_last ≈ 1.0–8.1 × p_cut`; truncated runs stop at `p_last ≈ 1.1e6 × p_cut`.

## 4. Deterministic repeats

| Test | Result |
|---|---|
| two fresh solvers, same process, N = 2500 | **bitwise identical** (`n = 641`, `R = 12.904228156540`) |
| two fresh solvers, same process, N = 10000 | **bitwise identical** (`n = 2635`, `R = 13.468323075955`) |
| same solver object solved twice, N = 2500 / 10000 | identical |
| fresh process | identical |

The anomaly is fully deterministic and reproducible.

## 5. `SetRadialRes` audit

`TOVSolver::SetRadialRes(const size_t&)` assigns `radial_res` and nothing else. Type `size_t`,
default `10000` (`TOVSolver.hpp:490`). **No clamping, no validation, no derived quantity
recomputed, no static or shared state.** `radial_res` is read exactly once per solve, inside
`SingleStarSolveToTOVPoints`, as `step = (r_max - r_min)/radial_res`.

**Ordering test:**

| Sequence | Result |
|---|---|
| 2500 in a fresh solver | `n = 641`, `R = 12.904228157`, `M = 1.597681092159` |
| 2500 after 10000 in the same solver | **identical** |
| 10000 in a fresh solver | `n = 2635`, `R = 13.468323076`, `M = 1.599975770860` |
| 10000 after 2500 in the same solver | **identical** |

`p_cut` and `init_press` identical across all of them. **No state-reset defect (category D) and no
EOS/interpolator state defect (category E).**

## 6. Target-radius generation

Authenticated loop (`TOVSolver.cpp`, `SingleStarSolveToTOVPoints`):

```cpp
double step = (max_log_r - min_log_r) / radial_res;     // r_min = 1 cm, r_max = 70e5 cm
double step_scale = 1.0;
for (double log_r_i = min_log_r; log_r_i <= max_log_r; log_r_i += step * step_scale)
{
    double ri = log_r_i;
    if (driver_apply(driver, &r, ri, y) != GSL_SUCCESS) break;   // <-- no point stored
    if      (ri <   100.0) step_scale = 0.005;
    else if (ri <  1000.0) step_scale = 0.025;
    else if (ri < 10000.0) step_scale = 0.05;
    else if (ri < 100000.0) step_scale = 0.25;
    else                    step_scale = 1.0;
    /* store TOVPoint */
    if (y[1] <= p_cut) break;
}
```

The variables named `log_r_*` are **ordinary radii in cm**; the spacing is piecewise uniform, not
logarithmic. `step_scale` is updated *after* the apply and therefore governs the *next* increment
— which is the intended reading and introduces no off-by-one.

The sequence was materialized exactly in IEEE-754 double arithmetic (`scratchpad/rr/targets.py`)
for `N ∈ {1250 … 40000}`. For `N = 2500`: `step = 2799.9996 cm`, effective spacing
`14 / 70 / 140 / 700 / 2800 cm` in the five bands, 450 targets beyond 1 km, 665 targets up to
13.5 km, 4 targets inside `[12.90, 12.99] km`. **No skipped interval, no duplicate, no backwards
target, no unusual gap, no early loop termination.** The generation at 2500 is entirely regular
and unremarkable beside its neighbours. **Category A is ruled out.**

## 7. `step_scale` threshold analysis

Distance from each threshold to the last target below it, for `N` around the anomaly:

| N | 100 cm | 1000 cm | 10000 cm | 100000 cm |
|---|---|---|---|---|
| 2450 | 13.286 | 41.857 | 113.287 | 684.729 |
| 2475 | **0.010** | 37.384 | 57.587 | 118.206 |
| 2490 | 0.606 | 43.177 | 117.475 | 17.087 |
| 2499 | 0.961 | 46.619 | 13.007 | 237.110 |
| **2500** | 1.000 | 47.000 | 17.001 | 277.014 |
| 2501 | 1.039 | 47.381 | 20.993 | 316.887 |
| 2550 | 2.922 | 65.667 | 75.472 | 36.269 |

`N = 2500` is unexceptional at every threshold; the closest approach in the whole sample is
`N = 2475` (0.010 cm, `7.1e11` ulp from 100 cm), which is **normal**. **No floating-point
loop-phase or borderline-comparison effect is involved.**

## 8. Floating-point loop phase

Measured, not inferred: the ulp distances above are `1e11–1e15` ulp — i.e. macroscopic, far from
any representable-boundary effect. The `<` comparisons at 100/1000/10000/100000 cm are never
close to a tie. **Ruled out.**

## 9. GSL target-partition experiment

Four orchestrations of the **same production RHS** (`TOVSolver::ODE`), varying only the driver
loop (scratch, `rr_part.cpp`):

| N | **A** production targets, one driver | **B** same targets, fresh driver per segment | **C** same targets, clamping RHS | **D** common 7 m grid, one driver |
|---|---|---|---|---|
| 2450 | R 13.45744 M 1.5999758 break | R 13.45744 M 1.5999758 break | R 13.48601 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| **2500** | **R 12.90423 M 1.5976811** break | R 13.46423 M 1.5999758 break | R 13.49223 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| **2510** | **R 12.90859 M 1.5979994** break | R 13.46636 M 1.5999758 break | R 13.49425 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| **2525** | **R 12.89429 M 1.5969173** break | R 13.44874 M 1.5999758 break | R 13.47646 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| 2550 | R 13.46924 M 1.5999758 break | R 13.46924 M 1.5999758 break | R 13.49669 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| **3000** | **R 12.90439 M 1.5976932** break | R 13.46439 M 1.5999758 break | R 13.48772 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| 5000 | R 13.46346 M 1.5999758 break | R 13.46346 M 1.5999758 break | R 13.47746 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |
| 10000 | R 13.46832 M 1.5999758 break | R 13.46832 M 1.5999758 break | R 13.47532 M 1.5999758 **cutoff** | R 13.46801 M 1.5999758 break |

- **B — resetting the driver (and hence the step size `h`) at every target removes the truncation
  completely.** Every resolution reaches the surface with `M = 1.5999758`.
- **C — making the RHS clamp instead of abort removes the truncation completely AND restores the
  governed `cutoff` termination.** `R` becomes 13.476–13.497 km, i.e. production's `R` is
  systematically **one target step short** even in the normal case.
- **D — a common grid gives every resolution the identical answer**, confirming that the target
  partition is the only channel through which `radial_res` reaches the physics.

## 10. The step size that decides it

Inherited `h` carried into each segment near the transition (production RHS, scratch loop):

| N = 2450 (normal) | | N = 2500 (truncated) | | N = 2550 (normal) | |
|---|---|---|---|---|---|
| target [km] | h_before [cm] | target [km] | h_before [cm] | target [km] | h_before [cm] |
| 12.886008 | 2344.9 → ok, p 7.08e31 | 12.876228 | 2499.4 → ok, p 8.65e31 | 12.892773 | 2266.7 → ok, p 6.03e31 |
| **12.914580** | **1381.9 → ok, p 3.10e31, ε 4.74e13 (crossed)** | **12.904228** | **1040.0 → ok, p 4.34e31, ε 7.06e13 (not yet crossed)** | **12.920224** | **1352.6 → ok, p 2.58e31, ε 4.47e13 (crossed)** |
| 12.943151 | 71.1 → ok (h collapsed) | **12.932228** | **1118.8 → FAIL 9, state not advanced** | 12.947675 | 43.6 → ok (h collapsed) |
| … to the surface | | *star ends here* | | … to the surface | |

The decisive difference is **not** that `h` is larger at 2500 — it is *smaller* (1119 cm) than at
2450 (1382 cm). What differs is the **phase**: whether the ≈0.44 m transition layer falls inside
the first trial step of a segment. When it does, an internal Runge–Kutta stage of that rejected
trial lands at `p < p_cut`, `TOVSolver::ODE` returns `GSL_EBADFUNC`, and
`gsl_odeiv2_driver_apply` — for which a user-function error is **fatal, not a step rejection** —
aborts without advancing the state (`p_after` unchanged, `ε` NaN). In the surviving runs the same
controller simply rejects the trial, collapses `h` to 40–70 cm, and walks through the feature.

That phase is a chaotic function of the whole preceding integration, which is why the bad set is
scattered rather than a band.

## 11. Surface / cutoff analysis

Last stored nodes, `N = 2450` (normal): `p` falls `5.53e29 → 1.13e26` over the final 13 targets,
`m` flat at `1.5999757`; the run then fails at target 13.486008 km with `p = 3.61e25 = 1.08 p_cut`
and **that point is never stored** (the `break` precedes the store). Hence:

- the governed surface test `p ≤ p_cut` is **unreachable in practice** on this EOS — the RHS guard
  always fires first;
- `R_*` is systematically **one target step short** of where the cutoff would place it
  (experiment C: 13.486 vs 13.457 at N = 2450, 13.475 vs 13.468 at N = 10000);
- across normal resolutions `R_*` scatters over **13.4400–13.4727 km (33 m, 0.24 %)** with `M`
  constant to `1.7e-8` — the residual `R` jitter seen in the Phase-4D convergence study.

Surface quantization explains the `R` jitter of the **normal** runs. It does **not** explain the
truncated runs, which end 0.56 km early and lose real mass.

## 12. Crust–core sampling

DS(CMF)-1's crust–core feature: `ε 4.93e13 → 6.71e13 g/cm³` over `p 3.519e31 → 3.570e31`
(`Δε/ε = 36 %` over `Δp/p = 1.5 %`), ≈0.44 m thick in the star at `r ≈ 12.91 km`.

Distinguishing the categories the task lists: it is **(A) the target grid does not resolve it —
but the ODE state is supposed to traverse it, and in the normal runs it does.** The truncated runs
are **(B) the star terminates before it**: the last stored node sits at `p = 4.3e31`, `ε = 7.06e13`
— just on the high-pressure side — and the integration never gets across. It is **not** a wrong
TOV branch (identical `ε_c`, `p_c`, `p_cut`) and **not** dropped output nodes (the store follows a
successful apply).

## 13. State-reset audit

`init_press` is reassigned every solve; the driver is allocated and freed inside the function; the
GSL interpolation accelerators are EOS-level and produced bitwise-identical results in every
ordering; `p_of_e_prec`, `r_min`, `r_max`, `central_eps_floor_factor` are not touched by the loop.
Measured: **no ordering dependence, bitwise determinism, fresh-process reproducibility.** No
state-reset defect was found. *(Not repaired, nothing to repair.)*

## 14. Sequence derivative, independently

`M(ε_c)` at `ε_c(1 ± rel)` with three step sizes, `dM/dp_c` from the solver's own central
pressures:

| N | rel = 5e-4 | rel = 1e-3 | rel = 2e-3 | interpretation |
|---|---|---|---|---|
| 2450 | 9671.2422 | 9671.2485 | 9671.2691 | stable, correct |
| **2500** | **9734.3873** | **9734.3905** | **9734.4043** | stable but **0.65 % biased** — all three stars truncated (common-mode) |
| **2550** | **38041.49** | 9671.2470 | 9671.2679 | **one-sided truncation** at 5e-4 (`R⁻ = 12.893`, `R⁺ = 13.469`) → derivative wrong by **4×** |
| 5000 | 9671.2424 | 9671.2474 | 9671.2686 | stable, correct |
| 10000 | 9671.2415 | 9671.2468 | 9671.2680 | stable, correct |

So the original `9734` observation is a **true difference in the solved `M(p_c)` curve** of the
truncated branch, not differentiation noise, not a central-pressure inversion artifact and not a
state-history artifact. The far more dangerous case is the *one-sided* one: a resolution whose
central star is perfectly normal can still yield a sequence derivative wrong by a factor of four.

## 15. Affected resolutions and central densities

| N | class | primary symptom | cause category |
|---|---|---|---|
| 1000–2499 (sampled), 2501, 2505, 2550–2750, 3500–40000 | NORMAL | `R` jitter ≤ 33 m | surface quantization |
| **2000, 4000** | **TRANSITIONAL** | central star normal, `dM/dp_c` off by 2× | one perturbed star truncated |
| **2500, 2510, 2525, 3000** | **ANOMALOUS** | star truncated at the transition, `ΔM = −2.3e-3 M☉`, `ΔR = −0.56 km` | fatal RHS guard × partition-dependent `h` |

**Central-density fragility at fixed resolution** (150 densities, `3.0e14 … 1.6e15 g/cm³`):

| radial_res | truncated | where |
|---|---|---|
| 5000 | **6 / 60 (10.0 %)** | 6.27e14, 7.03e14, 7.23e14, 8.57e14, 1.11e15, 1.35e15 |
| **10000 (production default)** | **5 / 150 (3.3 %)** | 1.322e15, 1.352e15, 1.414e15, 1.430e15, 1.530e15 — all `M ≈ 2.00–2.04 M☉` |
| 20000 | 0 / 60 in this sample | — |

**The production default is affected.** Not at the `ε_c` of the DS(CMF)-1 1.6 M☉ star, but at
≈3 % of central densities, concentrated in the 2 M☉ region that matters most observationally.

## 16. Scientific impact on prior work

**Directly contaminated — one durable artifact:**

`tests/baselines/grid_convergence_cmf_1p6_debug.tsv` and `grid_convergence_cmf_1p6_trajectory.tsv`
are produced by `tests/thermal/grid_convergence_cmf.cpp`, whose refinement sequence is
`kResolutions = {2500, 5000, 10000, 20000, 40000}` (`:107`). Reproduced exactly:

| row | `N_profile` | `R_km` | `achieved_M` | `ec_gcm3` |
|---|---|---|---|---|
| `A_fixed_ec 2500` (artifact) | 641 | **12.90422815654** | **1.597681155885** | 7.312533e14 |
| `A_fixed_ec 2500` (this audit) | 641 | **12.904228157** | **1.597681092** | 7.312533e14 |
| `A_fixed_ec 5000` | 1319 | 13.463458 | 1.599975833 | 7.312533e14 |
| `B_fixed_mass 2500` (artifact) | 641 | **12.90422815654** | 1.600076364 | **7.328604e14** |

Both `radial_res = 2500` rows are the truncated star. In `B_fixed_mass` the mass root-finder
**raised `ε_c` by 0.22 %** to recover the 0.0023 M☉ the missing crust took with it, so that row
describes a *different star* than the label implies. The derived thermal columns move accordingly
(`z_surf` 0.7965 vs 0.8056, `compactness` 0.3657 vs 0.3510, `B` 2.1211e57 vs 2.1246e57).

The convergence conclusion of that study — quantities approaching a limit as `dr_eff → 0` — was
therefore drawn with its **coarsest point off the solution branch**. The remaining four
resolutions are clean and monotone, so the conclusion is recoverable; the artifact's first row is
not a coarse-grid value of the same star.

**Not contaminated (verified by direct reproduction):**

| configuration | verdict |
|---|---|
| `SolveToProfile` 1.0 / 1.4 / 1.6 / 2.0 M☉ at res 10000 (the `I` golden, passive cooling, baryon number, all Phase-4 rotation work) | **normal** (`R` 13.426 / 13.545 / 13.468 / 12.712) |
| `tov_reference_cmf` configs (70 km, 10000), (70, 40000), (20, 10000), (20, 80000) | **normal** |
| `TaskManager` 30000, `BNV_Sequence` 90000, `LightDM`/`BNV_Chi` 10000 & 100000 | **normal** at this `ε_c` |

**Phase-4D and 4D-RG/RI conclusions.** Both excluded 2500 explicitly and used
5000/10000/20000/40000, all of which are normal at the 1.6 M☉ `ε_c`. **Excluding 2500 was
sufficient for those records** — with one qualification now measurable: the residual non-monotone
`δM̂` convergence recorded in `PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md` §15 (`R_*` moving
13.463458 → 13.471018 km over 5000→20000) is exactly the surface-quantization jitter of §11, i.e.
a TOV-background effect as that record already stated. Nothing in Phase 4D-RI's monopole
conclusions changes. **ADR-0008's implementation is not blocked**; its revalidation increment
should, however, treat `R_*` jitter as a known background artifact rather than a monopole
convergence property.

## 17. Root-cause classification

**Primary: `C` — SURFACE-DETECTION / OUTPUT-QUANTIZATION DEFECT.**
`TOVSolver::ODE` applies the surface test `y[1] < PressureCutoff()` to **trial** states inside the
right-hand side and signals it with `GSL_EBADFUNC`, which `gsl_odeiv2_driver_apply` treats as a
fatal user-function error rather than a step rejection. The loop already performs the correct test
on *accepted* states (`if (y[1] <= p_cut) break;`), which the guard pre-empts and makes
unreachable. Experiment C: removing only the fatal guard removes the anomaly at every resolution
and restores the governed termination.

**Amplifying factor: `B` — GSL TARGET-PARTITION COUPLING.** The single driver's step size `h` is
inherited across target segments, so the trial-step phase at the stiff feature — and therefore
which `(radial_res, ε_c)` pairs die — is set by the target grid. Experiment B: removing only the
inheritance also removes the anomaly.

**Physical trigger:** the DS(CMF)-1 crust–core near-discontinuity, ≈0.44 m thick, thinner than
every production target spacing.

**Downstream symptoms:** missing outer crust (`ΔM = −2.3e-3 M☉`, `ΔR = −0.56 km`); `ε_c`
displacement of 0.22 % when a mass root-finder compensates; sequence derivatives wrong by 0.65 %
(common-mode) or 4× (one-sided); a contaminated row in one durable artifact.

**Explicitly ruled out:** A (target generation — regular and unremarkable at 2500), D (state
reset — no ordering dependence), E (EOS/interpolator state — bitwise deterministic), G (the
anomaly reproduces), H (the cause is resolved).

## 18. Repair options — described, NOT implemented

| # | Option | Effect measured here | Class |
|---|---|---|---|
| R1 | Make the RHS non-fatal at low pressure (clamp `p` to `p_cut`, or return a finite RHS) and let the existing accepted-state test locate the surface | anomaly gone at every N; termination becomes `cutoff`; `R` grows by one target step (4–52 m) | **scientific/numerical-method — ADR required**: it changes `R_*`, hence `INV-06`'s surface semantics and every downstream artifact |
| R2 | Reset the driver (or `h`) at each target segment | anomaly gone at every N; `R`/`M` otherwise unchanged in the normal cases (`13.45744`, `13.46832` reproduced exactly) | **numerical-method — ADR required**: it changes the integration partition, though measurably *not* the normal-case answer |
| R3 | Locate the surface by interpolating between the last accepted node and the first sub-cutoff state | removes the systematic one-step `R` deficit as well | **scientific — ADR required** (INV-06) |
| R4 | Refuse to publish a profile that terminated by GSL failure above some pressure (fail closed) | would have caught every truncated star in this audit | behaviour-preserving **guard** for valid stars; still needs governance because it turns silent corruption into a hard failure |

R2 is the only option that left every *normal* result bit-for-bit unchanged in the measurements
above, which makes it the least invasive; it does **not**, however, fix the unreachable cutoff test
or the one-step `R` deficit. **No option was implemented and none is recommended here as a
decision** — that is the owner's.

## 19. Governance recommendation

An ADR is required before any repair. It should adjudicate, together: the surface-detection
contract (INV-06 — where `R_*` is and how it is found), whether the ODE right-hand side may carry
a fatal domain guard at all, whether the integration partition is part of the scientific contract,
and whether a GSL-terminated profile may be published. Repairing any one of these in isolation
would change `R_*` (and therefore every durable artifact) without the others being settled.
Until then: **`radial_res` values should not be varied casually, sequence sweeps at any resolution
can silently contain truncated stars, and `grid_convergence_cmf`'s `2500` row should not be read
as a coarse-grid value of the 1.6 M☉ star.**

## 20. Non-scope

No production file was modified (`git diff -- CompactStar/` = 0 lines at every stage). No repair,
no EOS change, no Hartle change, no baseline created or regenerated, no `radial_res = 2500`
investigation of its *thermal* consequences, no merge. The `grid_convergence_cmf` test and its
artifacts were **read, not modified** — re-deciding them belongs to the ADR above.
