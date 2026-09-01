# Grid convergence — Phase 2B-4A (nonrotating TOV → cooling slice)

> **STATUS: `TOV→COOLING CONVERGENCE CHARACTERIZED`.**
>
> Refinement **does** contract, on every observable and on the full cooling trajectory. The
> thermal-integrator error is **2.7×10³ times smaller** than the radial-refinement error, so
> radial convergence is cleanly isolated. The target-mass search contributes **exactly zero**
> differential error above `radial_res = 5000`.
>
> But the **measured order is ≈ 0.70, not 2** — and it drifts downward with refinement
> (0.85 → 0.70), so no convincing asymptotic regime is reached. Combined with the absence of
> any governed continuum-accuracy budget (§9), that is `CHARACTERIZED`, not `VERIFIED`.
>
> **This is the NONROTATING SLICE ONLY. The full ratified `TOV → Hartle → cooling`
> convergence item remains INCOMPLETE**, pending first-order Hartle validation and INV-07.

| Field | Value |
|---|---|
| **Change class** | numerical validation / verification |
| **Governing invariant** | INV-13 |
| **Starting commit** | `e1d369f8915e0825dd2ed80d45bbd1a9cd09909b` |
| **Harness** | `tests/thermal/grid_convergence_cmf.cpp` (CTest `grid_convergence_cmf`, labels `thermal;scientific;external-data;convergence`) |
| **Artifacts** | `tests/baselines/grid_convergence_cmf_1p6_debug.tsv`, `..._trajectory.tsv` |
| **Runtime** | ~106 s |
| **Build / toolchain** | `Debug`, AppleClang 17.0.0.17000604, CMake 4.2.1, GSL 2.7.1, macOS 15.6.1 arm64 |

---

## 1. What INV-13 asked for, and what this answers

INV-13: interpolation is **linear**, `DataSet::Integrate` is the **trapezoid rule**, nominally
`O(Δr²)` for smooth integrands — followed by the explicit warning that *"nominal order is not
observed order … The convergence order of the complete calculation must be measured, not
inferred from the interpolation scheme."*

**Measured. The answer is ≈ 0.70, not 2.** §7 explains why, and the explanation is structural
rather than accidental.

## 2. Scope boundary (binding)

Validated here: `TOV structure → StarProfile/GeometryCache → C_star, L_nu, L_γ → passive
thermal evolution → T_inf(t)`.

**Not** validated, claimed, or touched: the Hartle observable, `init_omega_bar`, INV-07, spin
evolution, chemical imbalance, heating.

**Disclosure.** `NStar::BuildFromTOV` calls `Find_MomInertia()` → `RotationSolver::
FindNMomInertia()` as a side effect of constructing *any* `StarProfile`. That is pre-existing
production behavior on the canonical path — the established passive-cooling baseline already
incurs it — and cannot be avoided without deviating from the very construction this study
measures. It is **inert** for the observable studied here (`couple_spin = false`, no spin
driver registered, `I` never enters `T_inf(t)`), and `I` is **never reported as a convergence
observable**. No statement about its accuracy is made or implied.

## 3. Data authority

Re-authenticated against the durable records; all hashes match.

| File | SHA-256 | Recorded in |
|---|---|---|
| `DS(CMF)-1_with_crust.eos` | `5747dd73…d47ae5dd` | `TOV_REFERENCE.md` |
| `eos.nb` | `d9c8e78c…f74edef8` | `TOV_REFERENCE.md` |
| `eos.thermo` (cold) | `41644499…c13e900ce3` | `TOV_REFERENCE.md` |
| `DNS-CMF-Hadronic-with-electrons/eos.thermo` | `a456fb85…0aa0724b` | `HEAT_CAPACITY_V1.md`, `PASSIVE_COOLING_BASELINE.md` |

No new dataset was introduced; no EOS data is committed.

## 4. Numerical controls, re-authenticated from source

**TOV.** `radial_res = 10000`, `r_max = 70 km`, `r_min = 1 cm`. Step
`= (r_max − r_min)/radial_res`, multiplied by a `step_scale` that refines toward the centre
(`0.005` below 1 m, rising to `1.0` beyond 1 km) — so the grid is **strongly non-uniform**
(§6). `SetRadialRes` and `SetMaxRadius` are public. Surface termination: the integration stops
when the driver first fails with `GSL_EBADFUNC`, raised by `ODE` when `p < PressureCutoff()`,
where `PressureCutoff() = max(1e-15 · p_c, eos_tab.pre[0])`; for this EOS the **table floor**
dominates. Central-density search in `SolveToProfile`: 25 coarse log-spaced samples, then
**bisection in linear `ε_c`** on the stable-branch bracket, `mass_tol = 1e-4 M☉`,
`max_iter = 40`.

**Thermal** (copied verbatim from `passive_cooling_regression.cpp`, so this is the same cooling
problem): `M = 1.6 M☉`, `T_inf(0) = 1e9 K`, `100 yr → 1 Myr`, RKF45, `rtol = 1e-6`,
`atol = 1e-10`, `samples_per_decade = 150`, `SaveCadence::LogTime`, Potekhin-1997 iron
envelope, PhotonCooling (`EnvelopeTbTs`, `radiating_fraction = 1`) + NeutrinoCooling
(DU + MU, no PBF), `couple_spin = false`, `n_eta = 0`. Nine epochs
`1e2 … 1e6 yr`.

## 5. Production equivalence of the variable-resolution builder — BIT-EXACT

`NStar::SolveTOV_Profile` owns its `TOVSolver` privately, so the harness constructs the solver
itself in order to call `SetRadialRes`. Everything downstream is the identical production
code: `SolveToProfile`, then `BuildFromTOV` via the **public** `NStar(points, labels)`
constructor. **No production API was added.**

At `radial_res = 10000`, all twelve checks return relative difference **exactly 0**:

| check | result |
|---|---|
| radial point count | 2635 = 2635 |
| achieved mass, radius, central density, central pressure, baryon number | rel `0` |
| surface redshift factor `e^ν(R)` = 0.805709 | rel `0` |
| `C_star`, `L_ν`, `L_γ` at `1e8 K` | rel `0` |
| **full 9-epoch cooling trajectory** | `max |Δln T_inf| = 0` |

The study is therefore performed on the production construction, not on a lookalike.

## 6. The resolution matrix and the ACTUAL radial spacings

`1/radial_res` is **not** the grid spacing — the grid is non-uniform by construction, so the
real spacings were measured:

| `radial_res` | `N_profile` | `Δr_eff = R/(N−1)` [km] | `Δr_min` [km] | `Δr_max` [km] | ratio max/min |
|---|---|---|---|---|---|
| 2 500 | 641 | `2.0163e-2` | `1.400e-4` | `2.800e-2` | 200 |
| 5 000 | 1 319 | `1.0215e-2` | `7.000e-5` | `1.400e-2` | 200 |
| 10 000 | 2 635 | `5.1133e-3` | `3.500e-5` | `7.000e-3` | 200 |
| 20 000 | 5 268 | `2.5576e-3` | `1.750e-5` | `3.500e-3` | 200 |
| 40 000 | 10 535 | `1.2790e-3` | `8.750e-6` | `1.750e-3` | 200 |

`Δr_max` is exactly `r_max/radial_res`, confirming the outer grid is uniform at that step, with
a 200× refined region inside 1 km. Successive `Δr_eff` ratios are `1.974, 1.998, 1.999, 2.000`
— close enough to 2 that order fitting uses the **actual measured ratio**, not an assumed one.

**`radial_res = 2500` is outside the asymptotic regime and is excluded from order fitting.**
It is retained, reported, and used as the detector (§10) — never silently dropped. Its radius
is wrong by 4.2 % and its 1 Myr temperature by 2 %.

## 7. Experiment A — fixed `ε_c*`, isolating radial discretization

`ε_c* = 7.312533426775e14 g/cm³`, the central density of the canonical 1.6 M☉ passive-cooling
star at the default resolution. Held fixed at every refinement, solved via
`SingleStarSolveToTOVPoints`.

| `radial_res` | `M` [M☉] | `R` [km] | `e^ν(R)` | `C_star(1e8 K)` | `L_ν(1e8 K)` |
|---|---|---|---|---|---|
| 2 500 | 1.59768116 | 12.9042282 | 0.79646477 | `2.4111008e38` | `4.4389016e39` |
| 5 000 | 1.59997583 | 13.4634581 | 0.80563049 | `2.4283250e38` | `4.4450566e39` |
| 10 000 | 1.59997583 | 13.4683231 | 0.80570917 | `2.4280622e38` | `4.4526531e39` |
| 20 000 | 1.59997583 | 13.4710181 | 0.80575272 | `2.4280435e38` | `4.4568634e39` |
| 40 000 | 1.59997584 | 13.4726806 | 0.80577958 | `2.4280646e38` | `4.4594614e39` |

### Observed orders (three finest, actual spacing ratio)

| quantity | successive |differences| (5k→10k, 10k→20k, 20k→40k) | observed order |
|---|---|---|
| `M` | `1.47e-9`, `6.49e-10`, `3.49e-10` | **at the roundoff floor** (rel `2.2e-10`) — no order reported |
| `R` | `4.865e-3`, `2.695e-3`, `1.662e-3` | **0.697** |
| `e^ν(R)` | `7.87e-5`, `4.36e-5`, `2.69e-5` | **0.698** |
| `L_ν(1e8 K)` | `7.60e36`, `4.21e36`, `2.60e36` | **0.697** |
| `d ln T_inf/dt` | `3.33e-10`, `1.75e-10`, `1.05e-10` | **0.730** |
| `L_γ(1e8 K)` | `4.28e30`, `8.18e29`, `3.28e29` | 1.320 |
| `C_star(1e8 K)` | `2.63e34`, `1.87e33`, `2.11e33` | **not reliably measurable** — sign flip at rel `8.7e-6` |
| baryon number `B` | `4.22e51`, `6.62e52`, `2.19e52` | **not reliably measurable** — non-monotone |

No order is reported where the differences are zero, sign-flipped, or at the quantity's own
roundoff floor. Producing a number there would be meaningless, so none is produced.

### Why ≈ 0.7 and not 2

`R`, `e^ν(R)` and `L_ν` all land on the **same** order, 0.697–0.698. That is not coincidence:
all three are controlled by the stellar radius, and **the radius is not a quadrature — it is a
termination event.** The DS(CMF)-1 table stops at `n_B = 1e-7 fm⁻³`, still inside the outer
crust (see `TOV_REFERENCE.md` §5), and the integration ends when the adaptive GSL driver first
*fails*. The solver's own verbose output shows the mechanism directly: the driver aborts a step
whose internal trial evaluation dips below the cutoff even though the recorded grid point is
still above it, so the truncation point depends on the requested step size. A quantity fixed by
"where the driver first fails" inherits neither the trapezoid's `O(Δr²)` nor any clean order,
and every integral over the star inherits its error.

`C_star` is the exception that confirms this: it is dominated by the dense core rather than the
surface, and it has already converged to `~1e-6` relative by `radial_res = 10000` — where the
differences hit a floor and flip sign.

**The observed order also drifts downward** with refinement (`R`: 0.852 on 5k/10k/20k → 0.697
on 10k/20k/40k). That is not a settled asymptotic regime, and it is the main reason this
increment is `CHARACTERIZED` rather than `VERIFIED`.

## 8. Experiment B — fixed target mass, the production workflow

| `radial_res` | achieved `M` | `|M − 1.6|` | selected `ε_c` [g/cm³] | `R` [km] | `T_inf(1 Myr)` [K] |
|---|---|---|---|---|---|
| 2 500 | 1.60007636 | `7.64e-5` | `7.32860432e14` | 12.9042282 | 457 529.85 |
| 5 000 | 1.59997583 | `2.42e-5` | `7.31253343e14` | 13.4634581 | 466 896.48 |
| 10 000 | 1.59997583 | `2.42e-5` | `7.31253343e14` | 13.4683231 | 466 966.06 |
| 20 000 | 1.59997583 | `2.42e-5` | `7.31253343e14` | 13.4710181 | 466 853.99 |
| 40 000 | 1.59997584 | `2.42e-5` | `7.31253343e14` | 13.4726806 | 466 802.31 |

**Target-mass-search contribution: zero above `radial_res = 5000`.** The bisection selects a
*bit-identical* `ε_c` at 5 000, 10 000, 20 000 and 40 000 — successive differences are exactly
`0.000e+00`. Experiment B is therefore the *same physical star* as Experiment A at those
resolutions, and the entire trajectory difference is radial discretization, uncontaminated by
the root finder. The mass search stays inside its own `1e-4 M☉` tolerance everywhere
(worst `7.64e-5`, at the coarse point). No central density was adjusted after the fact.

### Trajectory convergence

`D(h₁,h₂) = max_k |ln[T_inf(h₁,t_k)/T_inf(h₂,t_k)]|` over the nine canonical epochs.

| pair | `D_max` | RMS | final epoch | worst epoch |
|---|---|---|---|---|
| 2 500 → 5 000 | `2.0265e-2` | `8.899e-3` | `2.027e-2` | 1 Myr |
| 5 000 → 10 000 | `4.536e-4` | `3.661e-4` | `1.490e-4` | 300 yr |
| 10 000 → 20 000 | `2.400e-4` | `2.250e-4` | `2.400e-4` | 1 Myr |
| 20 000 → 40 000 | `1.435e-4` | `1.296e-4` | `1.107e-4` | 1 000 yr |

Contraction order of the norm over the two finest pairs: **0.742** — consistent with the
structural order, as it must be. `T_inf(1 Myr)` alone gives 1.117.

The per-epoch error between the default and finest grids is remarkably **flat**
(`3.82e-4` at 300 yr, `3.51e-4` at 1 Myr): the radial discretization sets the star's structure
once, and the resulting offset is carried through the whole trajectory rather than accumulating.

## 9. Thermal-integrator floor — demonstrably subdominant

Test configuration only; production defaults untouched. Both tolerances tightened 100×
(`rtol 1e-6 → 1e-8`, `atol 1e-10 → 1e-12`) with the physical problem unchanged:

| grid | `D_max` from tightening the time integrator | radial-refinement `D_max` at the same grid |
|---|---|---|
| `radial_res = 10000` | `1.432e-7` | `3.820e-4` (vs finest) |
| `radial_res = 40000` | `7.422e-8` | — |

**Conclusion: time-integration error is demonstrably subdominant — smaller by a factor of
2.7×10³.** Radial convergence is cleanly isolated; the stepper floor is nowhere near limiting.

## 10. Default-grid error, and the coarse-grid detector

| comparison | `R` | `C_star(1e8 K)` | `T_inf` trajectory `D_max` |
|---|---|---|---|
| `radial_res = 10000` vs finest computed (40 000) | `3.234e-4` | `9.929e-7` | `3.820e-4` |
| `radial_res = 10000` vs Richardson continuum | **`5.220e-4`** | n/a (order unmeasurable) | **`4.456e-4`** (from `T_inf(1 Myr)`) |
| `radial_res = 2500` vs finest | `4.219e-2` | — | `2.006e-2` |

**Detector proof.** The coarse member of the real refinement sequence — no artificial bug, no
modified source — is detected at **130× the default-grid radius error** and **53× its
trajectory error**, and refinement then reduces the discrepancy monotonically
(`2.03e-2 → 4.54e-4 → 2.40e-4 → 1.44e-4`). The harness can distinguish a materially
under-resolved star from the converged one.

### Richardson estimates, where the measured order supports them

| quantity | finest computed | extrapolated | `|diff|/finest` |
|---|---|---|---|
| `R` [km] | `13.4726806` | `13.4753568` | `1.986e-4` |
| `T_inf(1 Myr)` [K] | `466 802.31` | `466 758.07` | `9.476e-5` |
| `C_star`, `M` | — | **not justified** (order not reliably measurable) | — |

No production or golden value is replaced by an extrapolated one.

## 11. Accuracy adequacy — UNRESOLVED, deliberately

Per the task's own rule, the acceptance criterion was not invented after seeing the numbers.
The governed validation records were searched for an existing continuum-accuracy budget:

- `PASSIVE_COOLING_BASELINE.md`: `kTolState = 1e-5`, `kTolLumin = 1e-4`. These are explicitly
  **regression-detector** tolerances derived from repeat-run, tolerance and cadence
  measurements. They are **not** a continuum-accuracy requirement and are not reinterpreted as
  one here.
- `TOV_REFERENCE.md`: 0.5 % on `R`, 1 % on `M_max`. These are **external-comparison** budgets
  against `eos.mr`, explicitly stated to be looser than the observed agreement.
- `HEAT_CAPACITY_V1.md`, `CACHE_CORRECTNESS.md`: no downstream discretization budget.

**No authoritative downstream discretization budget exists.** The result is therefore reported
as a number, not a verdict:

```
default-grid (radial_res = 10000) continuum error
    R                 ~ 5.2e-4
    T_inf trajectory  ~ 4.5e-4
    C_star            ~ 1e-6
```

and accuracy adequacy is classified **unresolved**.

One observation worth recording without overstating it: the default grid's *continuum* error in
`T_inf` (`~4.5e-4`) is about **45× larger** than the passive-cooling regression's state
tolerance (`1e-5`). These measure different things — the baseline reproduces bit-identically
and remains a valid regression reference — but `1e-5` is a reproducibility tolerance and must
not be read as physical accuracy. Whether `~5e-4` is adequate is a question for whichever
downstream result needs it, and that budget does not yet exist.

## 12. Classification

**`TOV→COOLING CONVERGENCE CHARACTERIZED`.**

Satisfied: refinement contracts on every observable and on the trajectory norm; the
integrator floor is subdominant by 2.7×10³; the mass search contributes zero differential
error above 5 000; the default-grid error is quantitatively bounded; the coarse grid is
detected; no numerical pathology appeared.

Not satisfied for `VERIFIED`: the observed order is sub-linear (~0.70) **and drifting
downward**, `C_star` and `B` orders are not reliably measurable, and no governed accuracy
budget exists against which "sufficiently accurate" could be judged.

## 13. What remains

**Full `TOV → Hartle → cooling` convergence remains INCOMPLETE.** This increment covers the
nonrotating slice only. The Hartle leg requires the first-order moment-of-inertia validation,
which is itself blocked on the unresolved INV-07 normalization.

A second, separate observation for a future numerical-method task — **not acted on here,
because production is frozen**: the ~0.70 order traces to the surface-termination convention
(§7), not to the interpolation scheme. If a higher-order stellar structure is ever wanted, the
outer boundary is the thing to fix, not the quadrature. That is a numerical-method change under
`GOVERNANCE.md` and would move every existing baseline.

## 14. Reproduction

```bash
ctest --test-dir build -L convergence --output-on-failure
```

To regenerate the artifacts:

```bash
./build/tests/grid_convergence_cmf "$COMPACTSTAR_EOS_DATA_ROOT" --emit-dir tests/baselines
```

The registered CTest deliberately does **not** pass `--emit-dir`, so running the suite never
rewrites the recorded matrix.
