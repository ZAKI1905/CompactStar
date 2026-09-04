# First-order Hartle moment of inertia — Phase 2B-4B

> **STATUS: `HARTLE-I VERIFIED`** — as a **scale-free observable only**.
>
> The production first-order equation is derived, term by term, to be the published Hartle
> equation. The arbitrary central normalization cancels from `I = J/Omega` **analytically**
> and **numerically to 3.3e-15 over six decades** of central amplitude. Production agrees
> with a genuinely independent solver to **9.5e-9** on an exact analytic background and to
> **1.3e-5 – 2.1e-5** across the real CMF sequence, reproduces the Newtonian limit
> `I/(M R^2) -> 2/5`, and lands inside the published universal-relation bands.
>
> **INV-07 remains UNRESOLVED.** Absolute `omegabar`, physical `Omega`, physical `J`, the
> meaning of `init_omega_bar = 5e-3`, and the `Omega [s^-1]` unit annotation are all
> untouched by this work. **No claim is made about O(Omega^2) Hartle physics.**
>
> **No production source was modified.** The detector mutations in §21 were applied,
> measured and reverted; `git` confirms a byte-identical tree.

| Field | Value |
|---|---|
| **Change class** | scientific / numerical validation |
| **Starting commit** | `a24fe95f05574c7021f00abdf7c1a08cc8547aae` |
| **Harnesses** | `tests/rotation/hartle_moment_inertia_analytic.cpp` (self-contained), `tests/rotation/hartle_moment_inertia_cmf.cpp` (external data), `tests/rotation/hartle_reference.hpp` (independent solver) |
| **Artifact** | `tests/baselines/hartle_I_dscmf1_debug.tsv` |
| **Build / toolchain** | `Debug`, AppleClang 17.0.0.17000604, CMake 4.2.1, GSL 2.7.1, macOS 15.6.1 arm64 |

---

## 1. Literature authorities

- **J. B. Hartle (1967)**, *Slowly Rotating Relativistic Stars. I. Equations of Structure*,
  ApJ **150**, 1005 — the frame-dragging equation, its self-adjoint form, and the
  angular-momentum integral.
- **Breu & Rezzolla (2016)**, MNRAS **459**, 646 — the `Ibar(C)` universal relation, used as
  a **secondary sanity check only**.
- **Lattimer & Schutz (2005)**, ApJ **629**, 979 — the empirical `I`–`M`–`R` relation, again
  **sanity only**.
- Schwarzschild (1916) constant-density interior, as in MTW §23.7 / Shapiro & Teukolsky §5.5.

## 2. Metric, units, and the production call path

Convention `g_tt = -exp(2 nu)`, `g_rr = exp(2 lambda)`. `StarProfile` stores geometric units:
`r, m [km]`, `p, eps [km^-2]`, `nu, lambda` dimensionless. Hence `omegabar [km^-1]`,
`J [km^2]`, `Omega [km^-1]`, `I [km^3]`.

```
NStar::BuildFromTOV(...)                  <- runs for ANY profile construction
  └─ NStar::Find_MomInertia()             (NStar.cpp: seq.I = Find_MomInertia())
       └─ RotationSolver::FindNMomInertia()
            ├─ init_omega_bar = 5e-3      hard-coded (RotationSolver.cpp:701)
            ├─ y[0] = init_omega_bar, y[1] = 0     at r0 = first grid radius > 0
            │    (else kR_EPS_KM = 1e-6 km)
            ├─ gsl_odeiv2_step_rk8pd, h0 = 1e-1, abstol = 1e-10, reltol = 1e-10
            ├─ ODE_N_Fast + GetHartleOmegaCoeff_N_Fast / GetHartleDOmegaCoeff_N_Fast
            └─ J = R^4 y[1]/6 ; Omega = y[0] + R y[1]/3 ; I = J/Omega
```

Confirmed against source, not assumed.

### 2.1 Unit annotation audit — a real mislabel, recorded not repaired

| field | annotated | actually stored | verdict |
|---|---|---|---|
| `HartleResult::Omega` | `[s^-1]` | **`km^-1`** (geometric; the file multiplies by `LIGHT_C_KM_S` elsewhere when `s^-1` is wanted) | **WRONG** |
| `HartleResult::J` | `[km^2]` | `km^4 * km^-2 = km^2` | correct |
| `HartleResult::I` | `[km^3]` | `km^2 / km^-1 = km^3` | correct |
| `omega_bar`, `domega_bar` | `[km^-1]`, `[km^-2]` | as annotated | correct |

The mislabel is INV-07's recorded secondary defect. **It does not invalidate `I`**: the ratio
`J/Omega` is formed from two quantities in the same (geometric) system, so `I [km^3]` is
correct regardless of how `Omega` is labelled. Not repaired here — production is frozen.

## 3. The canonical Hartle equation, and that it is linear homogeneous

With `omegabar = Omega - omega(r)` and `j = exp[-(nu + lambda)]`, Hartle (1967) gives

```
(1/r^4) d/dr [ r^4 j d(omegabar)/dr ] + (4/r) (dj/dr) omegabar = 0 ,     omegabar'(0) = 0.
```

Every term is first degree in `omegabar` or its derivatives and there is no source term, so
the equation is **linear and homogeneous**.

### 3.1 Normalization cancellation — the proof

Under `omegabar -> A omegabar` (any constant `A != 0`): `omegabar' -> A omegabar'` and
`omegabar'' -> A omegabar''`, so every term scales by `A` and the equation is satisfied
identically. The centre condition `omegabar'(0) = 0` is itself scale-invariant. Therefore

```
J     = R^4 omegabar'(R)/6        ->  A J
Omega = omegabar(R) + R omegabar'(R)/3  ->  A Omega
I     = J/Omega                   ->  UNCHANGED.
```

`I` is the scale-free content of the first-order solution. `J` and `Omega` are not, which is
exactly INV-07.

### 3.2 The surface relations, derived independently

Outside the star `eps = p = 0` and `j = 1`, so `(r^4 omegabar')' = 0`, hence
`r^4 omegabar' = const`. Matching to the asymptotically flat exterior `omega = 2J/r^3`
(so `omegabar = Omega - 2J/r^3`) gives `omegabar' = 6J/r^4`. Evaluating both at `r = R`:

```
J     = R^4 omegabar'(R)/6                                  (matches production)
Omega = omegabar(R) + 2J/R^3 = omegabar(R) + R omegabar'(R)/3   (matches production)
```

## 4. Production-equation equivalence — `EQUATION MATCH`

Expanding the canonical form:

```
j omegabar'' + (4j/r + j') omegabar' + (4 j'/r) omegabar = 0
=>  omegabar'' = -4 omegabar'/r - (j'/j) omegabar' - (4/r)(j'/j) omegabar
```

The required identity, derived from `j = exp[-(nu+lambda)]` plus the TOV relations
`lambda' = (4 pi r eps - m/r^2)/(1-2m/r)` and `nu' = (m + 4 pi r^3 p)/(r^2 (1-2m/r))`:

```
j'/j = -(nu' + lambda')
     = -[ (m + 4 pi r^3 p)/r^2 + 4 pi r eps - m/r^2 ] / (1 - 2m/r)
     = -4 pi r (eps + p) / (1 - 2m/r)          <-- the m/r^2 terms cancel exactly
```

Substituting:

```
omegabar'' = -4 omegabar'/r
           + [4 pi r (eps+p)/(1-2m/r)] omegabar'
           + [16 pi (eps+p)/(1-2m/r)] omegabar
```

Production, `RotationSolver.cpp` `ODE_N_Fast` with its two coefficient helpers:

```
f[1] = -4 y[1]/r + [16 pi (p+e)/(1-2m/r)] y[0] + [4 pi (p+e) r/(1-2m/r)] y[1]
```

**Term by term identical. Classification: `EQUATION MATCH`.**

## 5. The independent reference solver

`tests/rotation/hartle_reference.hpp`. It does **not** call `ODE_N_Fast`,
`GetHartleOmegaCoeff_N_Fast` or `GetHartleDOmegaCoeff_N_Fast`, and does not restate their
expanded form. It integrates the **conservative** system in different state variables:

```
state:  (omegabar, q),   q = r^4 j omegabar'          [km^2]
d(omegabar)/dr = q / (r^4 j)              j = exp[-(nu+lambda)]  from the metric columns
dq/dr          = 16 pi r^4 (eps+p) exp(lambda-nu) omegabar
```

The RHS is built from `nu`, `lambda` and `(eps+p)` directly — never from production's
`1/(1-2m/r)` helpers. Linear interpolation of the tabulated background (INV-13; nothing
smoother is justified). GSL `rk8pd`, tolerances tightened well below production's.

**Centre treatment (independent).** `omegabar'(0) = 0` implies `q = r^4 j omegabar' = O(r^5)`,
so `q(r0) = 0` is the regular start with truncation error `O(r0^5)`. Production's
`kR_EPS_KM` prescription was deliberately not copied.

## 6. Centre sensitivity

Analytic background, seed 1.0, four start radii over four decades:

| `r0` [km] | `I` [km^3] |
|---|---|
| `3.250e-3` | `1.5713286903e2` |
| `1.0e-3` | `1.5713286903e2` |
| `1.0e-4` | `1.5713286903e2` |
| `1.0e-5` | `1.5713286903e2` |

Relative spread **`2.9e-15`** — far below the production/reference gap of `9.5e-9`, so the
reference is admissible as the authority.

## 7. Multi-seed invariance (mandatory experiment)

Analytic background, `M/R = 0.154`:

| seed | `omegabar(R)` | `J` | `Omega` | `I = J/Omega` |
|---|---|---|---|---|
| `1.000e-3` | `1.253e-3` | `2.298e-1` | `1.462e-3` | `1.5713286903e2` |
| `5.000e-3` | `6.265e-3` | `1.149e0` | `7.311e-3` | `1.5713286903e2` |
| `1.000e0` | `1.253e0` | `2.298e2` | `1.462e0` | `1.5713286903e2` |
| `1.000e3` | `1.253e3` | `2.298e5` | `1.462e3` | `1.5713286903e2` |

- `J` proportional to the seed to **`9.4e-15`**
- `Omega` proportional to the seed to **`1.0e-14`**
- **`I` invariant to `3.3e-15` across six decades**

The analytic proof of §3.1 is confirmed numerically to machine precision.

## 8. Surface vs volume extraction

The volume form, derived by integrating the conservative equation
(`R^4 j(R) omegabar'(R) = -4 int r^3 j' omegabar dr` with `j(R) = 1` and
`j' = -4 pi r (eps+p) e^{lambda-nu}`):

```
I = (8 pi/3) int_0^R r^4 (eps + p) exp(lambda - nu) [omegabar(r)/Omega] dr
```

evaluated by quadrature over the **interior solution**, not by reading the ODE state at `R`.

| case | `I_surface` [km^3] | `I_volume` [km^3] | relative |
|---|---|---|---|
| analytic, `M/R = 0.154` | `1.5713286903e2` | `1.5713288849e2` | `1.24e-7` |
| CMF 1.0 M☉ | — | — | `6.77e-7` |
| CMF 1.4 M☉ | — | — | `6.62e-7` |
| CMF 1.6 M☉ | — | — | `6.76e-7` |
| CMF 2.0 M☉ | — | — | `8.05e-7` |

This detects a corrupted **`J`** extraction (§21 D1 fires here). It cannot detect a corrupted
**`Omega`**, since both forms divide by the same `Omega` — that failure mode is covered by the
weak-field limit instead (§21 D2). The two tests are complementary by construction.

## 9. Weak-field limit — `I/(M R^2) -> 2/5`

Exact constant-density interiors at decreasing `M/R`, `R = 13 km`:

| `M/R` | reference `I/(M R^2)` | dev. from 0.4 | **production** `I/(M R^2)` | dev. from 0.4 |
|---|---|---|---|---|
| 0.150 | 0.462892464 | `1.57e-1` | 0.462892468 | `1.57e-1` |
| 0.100 | 0.438966166 | `9.74e-2` | 0.438966168 | `9.74e-2` |
| 0.050 | 0.418226223 | `4.56e-2` | 0.418226225 | `4.56e-2` |
| 0.020 | 0.407023177 | `1.76e-2` | 0.407023177 | `1.76e-2` |
| 0.010 | 0.403469508 | `8.67e-3` | 0.403469508 | `8.67e-3` |
| 0.005 | 0.401724450 | `4.31e-3` | 0.401724450 | `4.31e-3` |
| 0.002 | 0.400687334 | `1.72e-3` | 0.400687334 | `1.72e-3` |

Both approach `2/5` **monotonically**, and the residual falls essentially linearly in `M/R` —
the expected leading relativistic correction, not numerical error. This is a sequence, not a
hand-picked point.

## 10. Production vs reference on the analytic star

Production first reconstructs the analytic metric from `nu'` to `3.4e-10`
(`max |nu_prod - nu_exact|`), so the two solvers genuinely see the same background. Then:

```
I_production = 1.5713287051e2 km^3
I_reference  = 1.5713286903e2 km^3
relative difference = 9.455e-09
```

## 11. Real CMF sequence — DS(CMF)-1_with_crust

Data identity is the already-authenticated set (`TOV_REFERENCE.md` §3); hashes re-checked and
matching.

| `M` [M☉] | `R` [km] | `C = M/R` | `I_prod` [km^3] | `I_ref` [km^3] | prod/ref | surf/vol |
|---|---|---|---|---|---|---|
| 1.000044 | 13.426323 | 0.109985 | `8.699576e1` | `8.699435e1` | `1.61e-5` | `6.77e-7` |
| 1.400022 | 13.545323 | 0.152621 | `1.356161e2` | `1.356143e2` | `1.32e-5` | `6.62e-7` |
| 1.599976 | 13.468323 | 0.175416 | `1.595871e2` | `1.595837e2` | `2.15e-5` | `6.76e-7` |
| 2.000095 | 12.712323 | 0.232325 | `1.937231e2` | `1.937193e2` | `1.97e-5` | `8.05e-7` |

### 11.1 Physical magnitude and universal relations (SANITY ONLY)

`I[g cm^2] = I[km^3] x 1e15 x c^2/G`, with `c^2/G = 1.346590922e28 g/cm`, i.e. the factor
`1.346590922e43`.

> **Recorded, not silently absorbed — figure corrected by the Phase-2B closure audit.**
> `Zaki::Physics::SUN_M_KM = 1.476625038050 km` is **exactly the IAU nominal `GM_sun/c^2`**
> (`1.32712440018e26 / c^2 = 1.476625038 km`), i.e. it is the accurate, standard constant.
> GSL's own pair `G = 6.673e-8`, `M_sun = 1.98892e33 g` gives `G M_sun/c^2 = 1.476716 km`.
> The internally consistent discrepancy is therefore **6.2e-5**, not the `2.8e-4` originally
> stated here — that larger figure came from comparing against CODATA `G = 6.67430e-8` rather
> than the `G` the build actually uses. The gap is the well-known consequence of `G` being far
> less precisely known than `GM_sun`. It affects **only** the cgs number quoted for reference;
> every validation comparison in this document is performed in consistent geometric units and
> is untouched by it. Classified in `PHASE2B_CLOSURE.md` as Phase-3 unit-consolidation debt.

| `M` [M☉] | `I` [g cm^2] | `I/(M R^2)` | `Ibar = I/M^3` | Breu–Rezzolla | ratio | Lattimer–Schutz ratio |
|---|---|---|---|---|---|---|
| 1.0 | `1.1715e45` | 0.32681 | 27.017 | 25.294 | 1.068 | 1.048 |
| 1.4 | `1.8262e45` | 0.35754 | 15.350 | 14.742 | 1.041 | 1.045 |
| 1.6 | `2.1490e45` | 0.37238 | 12.102 | 11.766 | 1.029 | 1.036 |
| 2.0 | `2.6087e45` | 0.40589 | 7.520 | 7.554 | 0.996 | 0.998 |

Magnitudes are squarely in the neutron-star range. Both relations agree within **7 %**, well
inside their own quoted EOS scatter, with the deviation shrinking toward high mass. **These are
fits, not truth**; they are consistent with, and do not substitute for, the independent solver.

## 12. Publication provenance for DS(CMF)-1 `I` — none found

The official CompOSE DS(CMF)-1_with_crust distribution ships `eos.mr` with exactly **two**
columns, `R [km]` and `M [M☉]` — it contains **no moment of inertia**. No EOS-specific
published `I` sequence for this model was located in the distribution or in project history.

Classification of what exists:
- CompOSE `eos.mr` — **INDEPENDENT**, but carries no `I`; unusable for this purpose.
- Any prior CompactStar-generated `I` plot or sequence — **SAME-LINEAGE / NOT INDEPENDENT**;
  it cannot validate the solver that produced it and none is used as truth here.
- Breu–Rezzolla / Lattimer–Schutz — **INDEPENDENT but EOS-averaged fits**; sanity only.

The validation authority is therefore the independent solver of §5, exactly as intended.

## 13. Radial convergence of `I`

Fixed `ec* = 7.312533426775e14 g/cm^3` (the canonical 1.6 M☉ central density), on the 2B-4A
resolution matrix:

| `radial_res` | `N` | `dr_eff` [km] | `I_prod` [km^3] | `I_ref` [km^3] | prod/ref |
|---|---|---|---|---|---|
| 2 500 | 641 | `2.016e-2` | `1.59084426e2` | `1.59082407e2` | `1.27e-5` |
| 5 000 | 1 319 | `1.022e-2` | `1.59588111e2` | `1.59584580e2` | `2.21e-5` |
| 10 000 | 2 635 | `5.113e-3` | `1.59587141e2` | `1.59583716e2` | `2.15e-5` |
| 20 000 | 5 268 | `2.558e-3` | `1.59577623e2` | `1.59575294e2` | `1.46e-5` |
| 40 000 | 10 535 | `1.279e-3` | `1.59580720e2` | `1.59578035e2` | `1.68e-5` |

Successive `|dI|`: `9.69e-4`, `9.52e-3`, `3.10e-3` — **non-monotone, with a sign change**, so
**no convergence order is reported**; producing one would be meaningless.

**This answers §15's question directly.** `I` inherits the surface-location jitter that
Phase 2B-4A traced to the step-dependent termination event at the EOS table floor — but the
`r^4` weighting **does not amplify it**:

```
relative spread over radial_res 5000-40000:   I = 6.57e-5     R = 6.85e-4
```

`I` is an order of magnitude **less** grid-sensitive than the radius itself. `radial_res =
2500` is outside the asymptotic regime (2B-4A) and is excluded from the spread — reported,
not dropped.

Decisively, the **production/reference agreement is resolution-independent**, staying in
`[1.27e-5, 2.21e-5]` across a 16x range of grid spacing. The residual grid sensitivity
therefore belongs to the TOV background, **not** to the Hartle solver.

## 14. Reference numerical floor

| case | reference floor | production/reference gap | ratio |
|---|---|---|---|
| analytic star | `2.17e-15` | `9.46e-9` | `4.4e6` |
| CMF 1.4 M☉ | `1.26e-15` | `1.32e-5` | `1.1e10` |

**Reference error is demonstrably subdominant.** The production/reference comparison is
numerically resolved.

## 15. Detector proof

Four controlled production mutations, each applied, measured and reverted; `git` confirms a
byte-identical tree after every one, and the full suite is green afterwards.

| # | Mutation | Result |
|---|---|---|
| D1 | `J = R^4 y1/6` → `/5` | **FIRES** — 2 checks fail; `I` off by `2.000e-1`; weak-field limit breaks (`2.02e-1` from 2/5) |
| D2 | `Omega = y0 + R y1/3` → `/2` | **FIRES** — production/reference gap `6.68e-2` |
| D3 | ODE sign `-4 y1/r` → `+4 y1/r` | **FIRES** — 3 checks fail; gap `3.605e0`; weak-field limit diverges to `3.9e2` |
| D4 | centre `omegabar'(r0) = 0` → `1e-3`, `1e20`, `1e24` | **does not fire — and provably cannot** |

### 15.1 Why D4 cannot be a detector (a result, not a gap)

Instrumenting `FindNMomInertia` directly (temporary, reverted) with `omegabar'(r0) = 1e24`
versus `0`:

```
unmutated :  y0(R) = 6.22143542e-3   y1(R) = 2.31526129e-4
mutated   :  y0(R) = 4.14762362e18   y1(R) = 1.54350753e17
ratio     :  6.666667e20             6.666667e20      <-- IDENTICAL
```

Both components scale by the *same* factor, so the perturbation excited the **same regular
solution at a larger amplitude**, not a different solution. The mechanism: the seeded flux is
`q(r0) = r0^4 j omegabar'(r0) = 1e-20 x 1e24 = 1e4`; the irregular mode `omegabar ~ r^-3`
decays to nothing across the star while injecting the constant `q/(3 r0^3) = 3.3e18` — which
matches the observed `y0(R) = 4.1e18` — and a constant is precisely the regular solution.
Since `I = J/Omega` is invariant to amplitude (§3.1), `I` is unchanged.

**Conclusion: the task's suggested fourth detector is ineffective for this observable by
construction**, because the centre condition can only perturb `I` through a mode that both
decays as `r^-3` and renormalizes amplitude. This is a robustness property of the centre
treatment, not a blind spot in the validation. Three independent detectors do fire.

### 15.2 A maintenance hazard noted in passing

`RotationSolver.cpp` contains the byte-identical extraction block

```
ang_mom_J = pow(r_surface, 4) * y[1] / 6.;
ang_vel_Omega = y[0] + r_surface * y[1] / 3.;
```

at **two** live sites — line 488 in `Solve_Mixed` (the mixed dark/visible star path) and line
737 in `FindNMomInertia`. They are different functions serving different models, so this is
not a defect. It is recorded because a text-targeted edit silently hits the wrong one; any
future repair must address both deliberately.

## 16. The precise validation boundary

**VALIDATED**

- the production first-order frame-dragging equation, as the published Hartle equation;
- the surface extractions `J = R^4 omegabar'(R)/6` and `Omega = omegabar(R) + R omegabar'(R)/3`;
- the derived observable `I = J/Omega`, in geometric `km^3`, including its Newtonian limit,
  its agreement with an independent solver, and its physical magnitude.

**NOT VALIDATED — INV-07 REMAINS UNRESOLVED**

- the absolute value of `omegabar(r)` and `omega(r)`;
- physical `Omega`, physical `J`, and any requested physical spin;
- the meaning or correctness of `init_omega_bar = 5e-3`;
- the `HartleResult::Omega [s^-1]` annotation (measured wrong in §2.1, not repaired);
- **anything at O(Omega^2)** — the second-order candidate is untouched and unexamined.

Any quantity quadratic in `omegabar` still inherits the square of the arbitrary factor. Phase 4
remains blocked.

## 17. Reproduction

```bash
ctest --test-dir build -L hartle --output-on-failure
```

## TOV-SURF-MR artifact supersession note — 2026-09-03

The preceding record is historical. Current I values use the corrected TOV background. B3a now tests I against its own finest-grid value; event-radius spread no longer bounds a sampled integral. Independent first-order physics bounds remain unchanged.
Current artifact hashes, exact producer reproduction, regression evidence and the
remaining Phase-4D dependency gate: `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:889`, `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:960`,
`docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:1003`. No production scientific source changed during migration.
