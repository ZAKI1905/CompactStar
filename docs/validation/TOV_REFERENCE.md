# TOV reference validation — Phase 2B-2

> **STATUS: TOV REFERENCE VERIFIED.**
>
> The production TOV path reproduces the **exact Schwarzschild constant-density interior
> solution to machine precision** (max relative deviation `3.5e-16`), and reproduces the
> **official CompOSE `eos.mr`** mass–radius sequence for DS(CMF)-1_with_crust to
> **0.03 % in maximum mass** and **0.20–0.35 % in radius**. The radius residual is
> systematic, single-signed, and **fully attributed** to a surface-boundary convention
> (§5) — not to a defect in the equations, units, or integration.
>
> **No production scientific source was modified in this phase.** Three pre-existing
> characteristics are documented here as findings, deliberately **not repaired**: the
> converter's split nucleon-mass convention across the crust/core transition (§3.2), the
> outer-boundary convention (§5), and the default radial resolution (§6).

| Field | Value |
|---|---|
| **Change class** | verification / numerical validation |
| **Object under test** | `TOVSolver::ODE`, `TOVSolver::SingleStarSolveToTOVPoints`, `TOVSolver::SolveToProfile`, `NStar::SolveTOV_Profile` |
| **References** | exact Schwarzschild interior solution; official CompOSE `eos.mr`; published CMF anchors |
| **Harnesses** | `tests/core/tov_reference_analytic.cpp` (self-contained), `tests/core/tov_reference_cmf.cpp` (external data) |
| **Artifact** | `tests/baselines/tov_dscmf1_reference.tsv` |
| **Build / toolchain** | `Debug`, AppleClang 17.0.0.17000604, CMake 4.2.1, GSL 2.7.1, macOS 15.6.1 arm64 |

---

## 1. What was verified, and what was not

**Verified.** That the production code integrates the standard relativistic
Tolman–Oppenheimer–Volkoff system in CGS with correct constants, signs, GR correction
factors, EOS ingestion, center condition, and central-density search — to the precision
stated above, against references that are *independent of CompactStar*.

**Not verified here.** The DS(CMF)-1 equation of state itself; rotation; anything thermal.
No CompactStar-generated output is used as a reference anywhere in this phase.

---

## 1b. The live production path, its units, and its equations

```
NStar::SolveTOV_Profile(eos_file, target_M, out_dir)      <- what production calls
  └─ TOVSolver::ImportEOS(file)                           <- Hidden_ImportEOS_Vis
  └─ TOVSolver::SolveToProfile(target_M_solar, …)         <- ec bisection on the stable branch
       └─ TOVSolver::SingleStarSolveToTOVPoints(ec, …)    <- gsl_odeiv2 rk8pd driver
            └─ TOVSolver::ODE(r, y, f, params)            <- the right-hand side
  └─ NStar::BuildFromTOV(…)
```

Both live callers use exactly this path: `main/.../spin_therm_evol_2_main.cpp` and
`tests/thermal/passive_cooling_regression.cpp`. (`TOVSolver::RadiusLoop`, the second live TOV
path used for sequence scans, is not exercised by production profile construction and is not
covered by this document.)

**Units inside the integrator are CGS throughout**: `r [cm]`, `y[0] = m [g]`,
`y[1] = p [dyne/cm²]`, `GetEDens(p)` returns mass density `rho [g/cm³]`. Only at the point a
`TOVPoint` is stored is `r` converted to km and `m` to solar masses
(`GSL_CONST_CGSM_SOLAR_MASS = 1.98892e33 g`). The integrator uses
`gsl_odeiv2_step_rk8pd` with absolute and relative tolerance `1e-10`.

**The equations, as implemented** (`TOVSolver.cpp`, `TOVSolver::ODE`):

```
f[0] = dm/dr = 4 pi r^2 rho(p)

f[1] = dp/dr = -(G/r^2) (rho(p) + p/c^2) (m + 4 pi r^3 p / c^2)
               / (1 - 2 G m / (c^2 r))
```

with `G = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT` and
`c = GSL_CONST_CGSM_SPEED_OF_LIGHT`. This is the standard CGS Tolman–Oppenheimer–Volkoff
system. `rho(p)` is a `gsl_interp_steffen` (monotonicity-preserving cubic) spline of the
tabulated `eps` against `p`, built in linear coordinates over the 1191-point table.
§2 verifies this right-hand side directly rather than by inference.

---

## 2. Tier A — the equation itself, against an exact solution

`tests/core/tov_reference_analytic.cpp` builds a synthetic **constant-density** EOS in a
temporary directory, loads it through the production `TOVSolver::ImportEOS` path, and then
calls the production `TOVSolver::ODE` at exact interior-solution values.

For constant mass density `rho0`, with

```
m(r) = (4pi/3) rho0 r^3 ,  M = (4pi/3) rho0 R^3 ,  y = 2GM/(Rc^2) ,
A(r) = sqrt(1 - y r^2/R^2) ,  B = sqrt(1 - y) ,
```

the Schwarzschild 1916 interior solution and its closed-form derivative are

```
p(r)  = rho0 c^2 (A - B) / (3B - A)
dp/dr = -2 rho0 c^2 y r B / ( R^2 A (3B - A)^2 )
```

(the derivative obtained by differentiating `p(r)`, so it is an independent expression, not
a restatement of the TOV equation). See Misner, Thorne & Wheeler, *Gravitation* §23.7, or
Shapiro & Teukolsky §5.5.

### Result

| compactness `2GM/Rc²` | R | M | max rel. dev. `dm/dr` | max rel. dev. `dp/dr` |
|---|---|---|---|---|
| 0.30 | 9.821 km | 0.998 M☉ | `0.0` | `3.5e-16` |
| 0.50 | 12.679 km | 2.147 M☉ | `0.0` | `2.8e-16` |

Evaluated at `r/R = 0.02, 0.10, 0.25, 0.50, 0.75, 0.90, 0.98`. Both compactness values are
genuinely relativistic and safely below the Buchdahl bound `8/9`. Two further guards
confirm `dp/dr < 0` and that the GR correction factors are actually active — at
`y = 0.5, r = R/2` the production result is **2.15×** the Newtonian value, so a
Newtonian-limit implementation could not pass.

**The production right-hand side is correct to roundoff.** This is the strongest statement
in this document, and it disposes of the equations, the CGS constants `G` and `c`, the
signs, the `(rho + p/c^2)` and `(m + 4 pi r^3 p/c^2)` source terms, and the
`(1 - 2Gm/(c^2 r))` denominator.

The fixture is synthetic and incompressible; it asserts no neutron-star property.

---

## 3. Data authentication

The local data root is **byte-identical** to the official CompOSE distribution held in the
owner's cold archive:

| File | SHA-256 |
|---|---|
| `DS(CMF)-1_with_crust.eos` | `5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd` |
| `eos.nb` | `d9c8e78c2fcf37fe770fecfc2d3a211d840a28299821a56c77e66f9ff74edef8` |
| `eos.thermo` | `416444999ccac569e2c9b34808888949c36d759f30cce25dab0d42c13e900ce3` |
| `eos.mr` (**the reference**) | `98e38994b4191b1024f5f5edbada0e6c1abe9e0cb6d81de94e47be08e3e8efa8` |

`eos.mr` is 179 rows of `R [km]`, `M [M☉]`, produced by the **CompOSE project's own TOV
solver** from `eos.nb`/`eos.thermo`. Provenance per the distribution `README`: CMF #1 with
crust, from <https://compose.obspm.fr>, citing Dexheimer & Schramm, *ApJ* **683** (2008) 943;
*PRC* **81** (2010) 045201; Dexheimer, *PASA* **34** (2017); Dexheimer, Gomes, Klähn, Han &
Salinas, *PRC* **103** (2021) 025808.

None of these files is committed; they are large third-party data reached through
`COMPACTSTAR_EOS_DATA_ROOT`.

### 3.1 EOS dual representation — a measured discrepancy

Production TOV integrates the tabulated `DS(CMF)-1_with_crust.eos`; the reference `eos.mr`
was computed from `eos.nb`/`eos.thermo`. These are two representations of the same EOS on
the same 1191-point grid. Comparing them point by point
(`e = n_B m_n (Q7+1)`, `p = Q1 n_B`, `1 MeV/fm³ = 1.602176634e33 erg/cm³`):

| quantity | agreement |
|---|---|
| `n_B` | exact (`0.0`) |
| `p` | `<= 4.3e-9` relative, everywhere |
| `eps`, above `n_B = 4.0e-2 fm^-3` | `<= 5e-9` relative |
| `eps`, below `n_B = 4.0e-2 fm^-3` | **uniformly larger by `6.8866e-4`** |

The jump is perfectly sharp at `n_B = 4.0e-2 fm^-3` — the crust/core transition — and the
offset is reproduced exactly by

```
939.5653 / 938.9187125 - 1 = 6.886512e-4
```

### 3.2 Converter provenance — KNOWN, and it is CompactStar's own code

The origin is **not** an inconsistency in the third-party distribution. It is a deliberate,
explicit branch in CompactStar's own CompOSE-to-`.eos` converter,
`CompactStar/EOS/src/CompOSE_EOS.cpp` (`.eos` written at line 430):

| column | line | rule |
|---|---|---|
| pressure | 283–291 | `p = Q1 * n_B * MEV_FM3_2_Dyn_CM2` — no mass enters, which is why `p` matches `eos.thermo` to `4.3e-9` everywhere |
| energy density | 294–311 | `e = (Q7 + 1) * n_B * m * MEV_FM3_2_G_CM3`, with **`m = 939.5653` hardcoded** when `line_num-1 < crust_core_x_idx`, and `m = m_n` (parsed from the `eos.thermo` header, `938.9187125`) otherwise |

The source comments the branch as *"Checking if we are in the crust"*. `crust_core_x_idx` is
the first index at or above the constructor's `transition_density` argument
(`CompOSE_EOS.cpp:169`); for the shipped table that argument was `4.0e-2 fm^-3`, which is
precisely where the measured jump sits.

So production and CompOSE genuinely disagree about the crust energy density, by construction:
CompactStar normalizes the crust to the free neutron mass and the core to the CMF reference
mass, while `eos.thermo` — and therefore the `eos.mr` reference — expresses the entire table
against the single header mass.

**Which convention is physically correct is a scientific-semantic question and is explicitly
out of scope here**; this document neither adjudicates it nor changes it. What matters for
this validation is quantitative, and it was measured: rescaling the crust to the reference
convention and re-solving moves `R(1.4)` by **less than the 7 m radial grid step**, against a
residual of 30 m. The convention difference is therefore **real, understood, attributed to a
specific line of CompactStar source — and not the cause of the radius residual.** The cause is
§5.

---

## 4. Tier B — the production path against the official sequence

`tests/core/tov_reference_cmf.cpp`. Reference `M_max = 2.068869 M☉ @ R = 11.877 km`,
`R(1.4) = 13.5749 km`, `R(1.6) = 13.4948 km`.

### 4.1 Target-mass path — `NStar::SolveTOV_Profile`

This is the entry point production actually uses (`spin_therm_evol_2_main.cpp`,
`passive_cooling_regression.cpp`). Radii are compared at the **achieved** mass.

| target | achieved M | `R` production | `R` official | rel. dev. |
|---|---|---|---|---|
| 1.0 | 1.000044 | 13.426323 km | 13.472825 km | `3.45e-3` |
| 1.4 | 1.400022 | 13.545323 km | 13.574875 km | `2.18e-3` |
| 1.6 | 1.599976 | 13.468323 km | 13.494859 km | `1.97e-3` |

The central-density bisection hits every target mass to better than `1e-4` relative.

### 4.2 Fixed-central-density path and maximum mass

A 61-point sweep of `TOVSolver::SingleStarSolveToTOVPoints` over
`ec = 3e14 … 5e15 g/cm³`:

| quantity | production | official | rel. dev. |
|---|---|---|---|
| `M_max` | 2.069444 M☉ | 2.068869 M☉ | **`2.78e-4`** |
| `R(M_max)` | 11.8233 km | 11.8770 km | `4.52e-3` |
| `R(1.4)` via this path | 13.5428 km | 13.5749 km | `2.36e-3` |

`R(M_max)` is the loosest number in this document; near the turning point `R` varies steeply
with `ec` and the maximum is located only to within one sweep interval. That is a property of
the test's sampling, not of production.

The two production entry points — fixed central density and target mass — agree with each
other at 1.4 M☉.

### 4.3 Published sanity anchor

Independent of `eos.mr`: `M_max = 2.069 M☉` and `R(1.4) = 13.55 km` fall inside the published
CMF brackets, and `2GM/Rc² = 0.517` at maximum mass is comfortably sub-Buchdahl. These
brackets are deliberately wide; their purpose is to catch a gross scaling error that a
self-consistent comparison against `eos.mr` could not.

---

## 5. FINDING — the outer-boundary convention (documented, not repaired)

The radius residual in §4.1 is **systematic, single-signed** (production is always smaller)
and **grows toward low mass**. That is the signature of a surface boundary, not of a solver
error.

`PressureCutoff() = max(1e-15 * p_c, eos_tab.pre[0])`. For this EOS the table floor
dominates, so integration stops at `p = 3.351885e25 dyne/cm²`, `rho = 1.658808e8 g/cm³`
(`n_B = 1e-7 fm^-3`) — still **inside the outer crust**. The DS(CMF)-1 table simply does not
extend to vacuum. The layer below that floor carries negligible *mass* but tens of metres of
*radius*.

Estimating that omitted height hydrostatically — near the surface `p << rho c^2`, so
`dp/dr = -rho g_eff` with `g_eff = GM / (r^2 sqrt(1 - 2GM/(rc^2)))`, and with a local
polytrope `rho = rho0 (p/p0)^(1/Gamma)` fitted to the two lowest tabulated rows
(`Gamma = 1.3474`, i.e. relativistic degenerate electrons):

```
delta_r = (1/g_eff) * (p0/rho0) / (1 - 1/Gamma)
```

| M | residual `R_ref - R_prod` | omitted layer `delta_r` | fraction |
|---|---|---|---|
| 1.000 | 0.04650 km | 0.09399 km | **0.4948** |
| 1.400 | 0.02955 km | 0.06449 km | **0.4583** |
| 1.600 | 0.02654 km | 0.05392 km | **0.4921** |

The residual is a **constant fraction (0.46–0.49) of the omitted layer at every mass** — an
8 % spread across a 1.6× range in surface gravity and a 1.75× range in the residual itself.
A unit error, a wrong constant, or a defective integration would show no such scaling.

**Conclusion.** The entire mass-dependence of the radius residual is explained by where the
integration is terminated. CompOSE carries its sequence to a lower outer boundary than the
production cutoff. This is a **boundary-convention difference between two solvers on a table
that does not reach vacuum**, and it bounds the achievable radius agreement at roughly
0.2–0.35 % for this EOS.

Not repaired: the production scientific source is frozen for this phase, and changing the
surface convention is a numerical-method change that would move every existing result,
including the Phase-2B-1 passive-cooling baseline.

---

## 6. FINDING — default radial resolution (documented, not repaired)

`TOVSolver` defaults to `r_max = 70 km` with `radial_res = 10000`, i.e. a uniform **7 m** step
beyond 1 km. `NStar::SolveTOV_Profile` never calls `SetMaxRadius`, so production always runs
at the 70 km default — for a ~13.5 km star, ~80 % of the grid is outside the star. The header
comment on `SetMaxRadius` already anticipates exactly this ("if a maximum is known, say 15 km
for NS, it would increase the radial resolution").

Measured at `ec = 9e14 g/cm³`:

| `r_max` | `radial_res` | step | M | R |
|---|---|---|---|---|
| 70 km | 10 000 | 7.00 m | 1.79814532 | 13.26532310 km |
| 70 km | 40 000 | 1.75 m | 1.79814532 | 13.26793060 km |
| 20 km | 10 000 | 2.00 m | 1.79814532 | 13.26800337 km |
| 20 km | 80 000 | 0.25 m | 1.79814532 | 13.26799462 km |

The mass is converged to **9 significant figures** and does not move at all. The radius moves
by `2.0e-4` relative (2.7 m) between the production default and an 80×-finer grid, and is
converged by `r_max = 20 km, res = 10000`. So the default costs ~2.7 m of radius accuracy —
an order of magnitude below the boundary-convention effect of §5, and well inside the
tolerance budget. It is a cheap improvement that is nevertheless a numerical-method change
and therefore out of scope here.

---

## 7. Center and surface condition audits

**Center.** `r_min = 1 cm`; `p(r_min) = p_c`; `m(r_min) = (4/3) pi r_min^3 rho(p_c)`. This is
the correct leading-order regular series, and it is reproduced to `<1e-6` relative
(the residual is pure floating-point rounding: `3769911184307730.0` vs
`3769911184307730.5` g). `r_min / R = 7.5e-7`, so the central offset is negligible.

**Surface.** Integration terminates at the EOS table floor (§5), `p_surf/p_c = 3.2e-10`, and
the mass has stopped changing to `<1e-6` relative over the final step. The last retained
point sits within 1.36× of the table floor, i.e. within one grid step of it — confirming the
surface is set by the table and by nothing else.

---

## 8. Consistency with the thermal baseline

The 1.6 M☉ star underlying `PASSIVE_COOLING_BASELINE.md` (`M = 1.599976`, `R = 13.468323 km`)
is reproduced exactly by this harness and sits `1.97e-3` from the official reference —
inside the same budget as every other mass. The Phase-2B-1 thermal baseline therefore rests
on a TOV solution that is independently corroborated. This check compares the *cooling
star's own parameters* against the **external** reference; it does not use CompactStar output
as a reference.

---

## 9. Tolerance policy

Tolerances were fixed **a priori** from the numerics, not tuned to observed results:

| source | magnitude |
|---|---|
| radial grid quantization of `R` (§6) | `~5.2e-4` |
| crust energy-density convention (§3.1) | `~6.9e-4` |
| reference `eos.mr` interpolation & 0.001 km quantization | `~1e-4` (measured: linear vs local cubic differ by `<=0.001 km`) |
| target-mass bisection | `<1e-4` |

Adopted with margin: **0.5 % on radius**, **1 % on maximum mass**, `2x` on `R(M_max)` for the
sampling reason in §4.2. Every observed deviation is inside budget. Tier A carries no budget —
it is checked at `1e-8` and delivers `1e-16`.

The budget is *looser* than the observed agreement because Tier A already certifies the
differential equation to roundoff. Tier B's job is to catch gross defects in EOS ingestion,
unit handling, the center/surface treatment and the central-density search — and the
attribution in §5 shows the residual it does see is understood rather than merely tolerated.

---

## 10. Reproduction

```bash
cmake -S . -B build -DPython3_EXECUTABLE="$(command -v python3)" \
      -DCOMPACTSTAR_EOS_DATA_ROOT=/path/to/compose
cmake --build build -j8
ctest --test-dir build -L tov --output-on-failure
```

`tov_reference_analytic` needs no external data. `tov_reference_cmf` reports CTest **SKIP**
(exit 77) if `eos.mr` is absent from the distribution, and fails loudly on any real
disagreement.
