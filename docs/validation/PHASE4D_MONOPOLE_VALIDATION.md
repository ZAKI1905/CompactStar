# Phase 4D — independent physical validation of the governed Hartle O(Ω²) monopole response

> **FORMAL STATUS: `HARTLE MONOPOLE VALIDATION FAILED`** — precisely scoped:
>
> * The **implementation** of ADR-0007 (`RotationSolver::ODE_HartleMonopole_`,
>   `NStar::ComputeHartleMonopoleResponse()`, `HartleMonopoleResponse`, materialization) is
>   **independently verified**: two test-only solvers in different variables that never call
>   production's ODE, the continuum limit, the regular-centre series, the Newtonian limits, the
>   surface shell, the exact homogeneous-family derivative, and the published second-order tables
>   of Chandrasekhar & Miller (1974) and Hartle & Thorne (1968) all agree with production within
>   their predeclared bounds. Detectors M1–M9 all fire.
> * The **accepted contract** is physically incomplete on a tabulated EOS with density
>   discontinuities: ADR-0007 P2 term 1 with the smooth P5 derivative authority cannot carry
>   Hartle's `dE/dP` delta functions at *internal* density steps. On DS(CMF)-1 the dominant step is
>   the table's crust–core transition (`Δε/ε = 36 %` over `Δp/p = 1.5 %`); the omitted internal
>   shells are worth **≈ 4.6 % of `δM̂`** (≈ `1e-3` of the homogeneous `δM̂`), established by an
>   independent route (§14). ADR-0007 §7 item 11 is **NOT MET** (`1.04e-3` vs `1e-3`).
> * Per the 4D mandate the discrepancy is **reported, not repaired**. **No monopole baseline was
>   created.** `GOVERNANCE.md` §3.1: **CORRECTION EXECUTED — IMPLEMENTATION INDEPENDENTLY VERIFIED —
>   PHYSICAL VALIDATION FAILED ON STEPPED CRUSTS — NO BASELINE YET.**
>
> **SLOW-ROTATION DISCLAIMER.** Everything here establishes the O(Ω²) *coefficients*. Nothing here
> bears on the accuracy of truncating a rapidly rotating star at O(Ω²).

| Field | Value |
|---|---|
| **Starting HEAD** | `377bc4a7f0233c6b44ca0a66393d26045761fcd3` (4C-I1), branch `physics/rotation-correctness`, upstream equal, clean, **8 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Change class** | **test-only + documentation** — production diff **NONE** (`git diff -- CompactStar/` = 0 lines at every stage; detector mutations reverted byte-identically, SHA-256 verified) |
| **Governing authority** | ADR-0007 §7 (predeclared validation requirements) and §8 (detectors); ADR-0006 (first order); ADR-0003 (provenance); ADR-0005 (`I` and the goldens bitwise); INV-08, INV-09, INV-13 |
| **Primary sources** | Hartle (1967) ApJ 150, 1005 — eqs. (87)–(90), (97)–(100), (105)–(108); Hartle & Thorne (1968) ApJ 153, 807 — Tables 1, 3, 5; Chandrasekhar & Miller (1974) MNRAS 167, 63 — eqs. (18)–(19), (32), Table I, misprint list |
| **Pre-task baseline** | full **27/27** (212.65 s), self-contained **14/14** (13.70 s); seven artifact hashes as `PHASE3_CLOSEOUT.md` §1 |
| **Post-task suites** | full **30/30** (223.54 s), self-contained **16/16** (18.41 s); seven hashes unchanged (§20) |

---

## 1. Claim under test and what is not claimed

**Claim (ADR-0007).** For an ordinary `NStar` the fixed-central-energy-density `l = 0` response —
`m̂₀ = m₀/Ω²`, `p̂₀* = p₀*/Ω²`, `δp̂₀`, `ξ̂₀`, `δM̂` — computed by `RotationSolver::ODE_HartleMonopole_`
from the Phase-4B-verified `s = ω̄/Ω`, `s' = ω̄'/Ω` and the Phase-4C-I0 `dε/dp` authority is a
correct numerical implementation of Hartle's primary equations on CompactStar's governed background
conventions (`ν_H = 2ν`, `j² = e^{−2ν}(1 − 2m/r)`).

**Not claimed:** `l = 2`, the quadrupole, the slow-rotation truncation error, baryon-conserving
sequences, rotochemical coefficients. Not used as oracles anywhere: the retired candidate, any
`HartleResult` value, any `PHASE4_ROTATION_ENTRY.md` diagnostic, any 4C-I1 CMF diagnostic value.

## 2. Read-only versioning hazard (recorded, not fixed)

`NStar::GetSequence()`'s **non-const** overload calls `SeqMutable()` → `Touch()`, so reading the
sequence through a mutable star bumps `Version()` and invalidates the version-keyed monopole
response (ADR-0003). Every 4D harness reads through `const NStar &`. Classification:
**PRE-EXISTING CONSERVATIVE INVALIDATION / API-HYGIENE DEBT** — the response is never stale, it is
recomputed; no production change in 4D. No supposedly-const operation bumped `Version()` in any
run.

## 3. The independent reference (`tests/rotation/hartle_monopole_reference.hpp`)

Production integrates Hartle's `(m₀, p₀*)` pair, (97)+(100). The reference integrates the
**`(m₀, h₀)` pair**, (97)+(98), with `p₀*` an *algebraic by-product* of the first integral (90):

```
dm̂₀/dr = 4πr²(ε+p)(dε/dp) p̂₀* + (1/12) r⁴ e^{−2ν}(1−2m/r) s'² + (8π/3) r⁴ (ε+p) e^{−2ν} s²
dĥ₀/dr = m̂₀ (1+8πr²p)/(r−2m)² + 4π(ε+p) r² p̂₀*/(r−2m) − (1/12) r³ e^{−2ν} s'²
p̂₀*    = γ̂ − ĥ₀ + (1/3) r² e^{−2ν} s²
```

Different state vector, different pressure-side equation, its own binary-search interpolation, its
own centre initialisation (the leading power of each right-hand side used numerically, not the
closed-form series transcribed), its own exterior arithmetic, tolerances `1e-13/1e-16`
(production `1e-10`). It never calls `ODE_HartleMonopole_`, `ComputeHartleMonopoleResponse`,
`HartleMonopoleResponse::At`, `Geometry::MetricDenominator` or any production helper.
**Re-authentication:** Chandrasekhar & Miller (1974) eqs. (18)–(19) state exactly this `(m₀, h₀)`
system in CompactStar's metric convention; their misprint list confirms the Hartle (117)
`(a − M) → (a − 2M)` misprint found in 4C-G. The equivalence of (100) to (98)+(90) rests only on the
tabulated `ν'` being the derivative of the tabulated `ν` — `O(h²)` on a linearly interpolated
background — and 4D **measures** that floor (§9) instead of assuming it away.

A second, stronger oracle is the **continuum solver** on the closed-form Schwarzschild interior
(`SolveUniformContinuum`): no tabulation at all, first and second order carried together, `p₀*`
obtained **both** from (90) and by integrating (100), so Hartle's first integral is tested at the
continuum level.

**Two chains** are reported everywhere: **SECOND-ORDER-ISOLATED** (production's verified `s`, `s'`
+ the independent second order) and **FULLY INDEPENDENT** (`hartle_reference.hpp`'s own first
order + the independent second order).

## 4. Test organisation

| Test | Data | Content |
|---|---|---|
| `hartle_monopole_physics_analytic` | self-contained | Experiments A, B, C, D, E, J, R, S on the exact constant-density interior |
| `hartle_monopole_physics_cmf` | external DS(CMF)-1 | Experiments F, G, H, I, J, R on 1.0/1.4/1.6/2.0 M☉ |
| `hartle_monopole_published` | self-contained (fixture embedded) | Experiment K: Chandrasekhar & Miller 1974 Table I; Hartle & Thorne 1968 Tables 3/5 on the printed HW EOS |

Fixture: `tests/rotation/hartle_thorne_1968_hw_eos.hpp` (§17).

## 5. Predeclared bounds and the two labelled re-scopings

All bounds are ADR-0007 §7 as accepted; **none was widened**. Two implementation choices in the
harness were corrected after a first run and are recorded verbatim in the source:

* **Item 7 (Newtonian limits).** The ADR says *intercept*. The first run asserted the
  weakest-field value instead (stricter than the ADR) and measured `4.5e-3` (`δM̂/R³`) and
  `5.3e-3` (`3Mξ̂/R⁴`) at `M/R = 0.002`. The intercept is the accepted criterion; the leading
  deviation is first order in `M/R` by the post-Newtonian expansion of the closed-form interior
  (measured ratio between `M/R = 0.002` and `0.001`: `2.000`, `2.000`). `M/R = 0.001` was added
  afterwards and is labelled.
* **Item 6 (first integral).** The first run asserted the *tabulated* form at `N = 16001` and
  measured `1.918e-9` against `1e-9`, with exact `h²` scaling (ratios `4.00` at every doubling,
  `4.8e-10` at `N = 32001`, added afterwards). That is the `O(h²)` floor of a piecewise-linear `ν`
  against a tabulated `ν'` — an input property, not a solver property. The ADR criterion is taken
  at the continuum level (`6.1e-15`); the tabulated value is reported.

Two further harness defects corrected after a first CMF run, both recorded in the source: the
near-vacuum window `(ε+p)r² < 1e-6` also admitted the central nodes (where `I²/r³ ≈ 1e19`) and was
replaced by `ε < 1e9 g/cm³, r > R_*/2`; and the FD-sensitivity line was asserted against the ADR
bound although it had been declared diagnostic-grade — the assertion was withdrawn, the number kept.

## 6. Experiment R — reference admissibility

| Background | Reference self-movement (rtol `1e-11…1e-15`, `r₀ = 1e-5…1e-7 km`) | Disagreement with production | Ratio |
|---|---|---|---|
| analytic, N = 4001 | `0.0` (bitwise; with rk4 at rtol `1e-7` it moves `6.6e-15`, so the metric is live) | `9.7e-9` | `0` |
| DS(CMF)-1, 1.6 M☉ | `9.96e-15` | `1.96e-7` | `5.1e-8` |

Bound `0.1` — met. The reference is converged to double precision on both backgrounds; every
disagreement below is production-side or a formulation floor.

## 7. Experiment A — regular-centre series (bound `1e-8`, first ten nodes)

Fine-centre fixture: first twelve intervals at `1e-4 km` from `r₀ = 1e-5 km`, then 2000 uniform
nodes; the `O(r²)` series correction over the first ten nodes is `< 1e-9`, so the bound is meaningful.

| Coefficient | Worst relative error (nodes 0–9) |
|---|---|
| `p̂₀*/r² → (1/3) j_c² s_c²` (`a₂ = 0.2786`) | `1.26e-10` |
| `m̂₀/r⁵ → (4π/15)(ε_c+p_c)(dε/dp+2) j_c² s_c²` (`b₅ = 3.385e-4`) | `5.02e-10` |
| `ξ̂₀/r → j_c² s_c²/[4π(ε_c+3p_c)]` (`c₁ = 229.0`) | `1.24e-10` |

**Met.**

## 8. Experiment B — full profile against the `(m₀, h₀)` solver (bound `1e-7` at N = 4001)

`M = 2.0 km, R = 13.0 km` exact interior; `dε/dp = 0` supplied through the governed mechanism.

| N | ISO `m̂₀` | ISO `p̂₀*` | ISO `ξ̂₀` | ISO `δM̂` | FULL `m̂₀` | FULL `p̂₀*` | FULL `ξ̂₀` | FULL `δM̂` | (90) residual |
|---|---|---|---|---|---|---|---|---|---|
| 2001 | `0` | `3.89e-8` | `3.89e-8` | `3.48e-8` | `1.93e-8` | `2.76e-8` | `2.76e-8` | `2.41e-8` | `1.23e-7` |
| 4001 | `0` | `9.72e-9` | `9.72e-9` | `8.71e-9` | `4.82e-9` | `6.89e-9` | `6.89e-9` | `6.03e-9` | `3.07e-8` |
| 8001 | `0` | `2.43e-9` | `2.43e-9` | `2.18e-9` | `1.21e-9` | `1.72e-9` | `1.72e-9` | `1.51e-9` | `7.67e-9` |

Exact `h²` scaling; worst node at the surface (`p̂₀*`) and at `r = 6.5e-3 km` (`m̂₀`, FULL). On a
`dε/dp = 0` star the isolated `m̂₀` equation decouples and both integrators reproduce it to
`2.5e-32` RMS — so term 1 is exercised only on DS(CMF)-1 (§12). **Met** in both chains.

## 9. Experiment B-cont and E — continuum solver and the first integral

| N | vs continuum `m̂₀` | `p̂₀*` | `ξ̂₀` | `δM̂` | continuum (100)-vs-(90) residual |
|---|---|---|---|---|---|
| 1001 | `8.18e-8` | `7.93e-8` | `8.01e-8` | `4.22e-8` | `1.1e-15` |
| 2001 | `2.05e-8` | `1.98e-8` | `2.00e-8` | `1.06e-8` | `2.9e-15` |
| 4001 | `5.12e-9` | `4.96e-9` | `5.01e-9` | `2.64e-9` | `6.1e-15` |
| 8001 | `1.28e-9` | `1.24e-9` | `1.25e-9` | `6.60e-10` | `6.0e-15` |

Production converges to the continuum solution as `h²`. **Continuum first integral** (bound `1e-9`,
explicit absolute scale = max `|p̂₀*|`): `6.1e-15` — **met**. Tabulated first integral on production's
own fields: `1.23e-7 → 3.07e-8 → 7.67e-9 → 1.92e-9 → 4.80e-10` for `N = 2001…32001`, ratios `4.00`
(§5).

## 10. Experiments C and D — the surface shell and the Newtonian limits

**C.** `δM̂ = m̂₀(R⁻) + shell + I²/R³` exact (`0.0`); the shell agrees with `4πR²ε ξ̂₀(R)` from the
independent solver to `9.7e-9`; on the constant-density star **shell/δM̂ = 0.896**, `m̂₀(R)/δM̂ =
0.097`, `(I²/R³)/δM̂ = 0.0077`. The shell is load-bearing.

**D.** `R = 13 km`, `M/R = 0.15 … 0.001`:

| M/R | `δM̂/R³` | `3Mξ̂/R⁴` | without shell `(δM̂−shell)/R³` |
|---|---|---|---|
| 0.15 | 0.675828 | 0.608084 | `6.77e-2` |
| 0.05 | 0.888872 | 0.868116 | `2.08e-2` |
| 0.01 | 0.977553 | 0.973525 | `4.03e-3` |
| 0.002 | 0.995502 | 0.994701 | `8.01e-4` |
| 0.001 | 0.997751 | 0.997350 | `4.00e-4` |

Both monotone; production and the independent solver identical to the digits shown. Linear-in-`M/R`
intercepts: `|δM̂/R³ − 1| = 1.07e-6`, `|3Mξ̂/R⁴ − 1| = 4.9e-7` (bound `5e-3`) — **met**. Omitting
the shell gives `4e-4` instead of `1` (`→ 0.4 M/R`, i.e. `IΩ²/M`): detector M4's target.

## 11. Experiment J on the analytic star — the homogeneous family (bound `1e-3`)

Test-side only (ADR-0007 P11 as modified). Sources off, `p̂₀*_c = 1`, `δp_c = (ε₀+p_c)`; the exact
constant-density sequence has `dM/dp_c = (3/4) R y_R (3y_R − 1)²/ρ₀`. Result: `δM̂_hom = 20.20` vs
`(dM/dp_c) δp_c = 20.20`, **rel `3.05e-9`**; with `dε/dp = 0` the homogeneous `m̂₀` vanishes
identically (mass moves only through the shell). **Met.** This certifies the reference's
(98)+(90) machinery against exact analytic knowledge.

## 12. Experiment G — DS(CMF)-1, four stars (bound `1e-4`)

| M☉ | `R_*` km | `(dε/dp)_c` | `m̂₀(R_*)` | `p̂₀*(R_*)` | `ξ̂₀(R_*)` | shell | `I²/R_*³` | `δM̂` |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 13.42632 | 5.2825 | 930.99 | 35.338 | 3364.7 | `9.8e-4` | 3.127 | 934.11 |
| 1.4 | 13.54532 | 3.7908 | 885.53 | 34.252 | 2111.8 | `6.3e-4` | 7.400 | 892.93 |
| 1.6 | 13.46832 | 3.2786 | 815.55 | 32.446 | 1617.1 | `6.2e-4` | 10.424 | 825.97 |
| 2.0 | 12.71232 | 2.5624 | 537.88 | 24.053 | 704.52 | `2.6e-4` | 18.268 | 556.15 |

| M☉ | ISO worst `m̂₀ / p̂₀* / δp̂₀ / ξ̂₀ / δM̂` | FULL worst `m̂₀ / p̂₀* / ξ̂₀ / δM̂` |
|---|---|---|
| 1.0 | `1.5e-8 / 1.0e-7 / 1.0e-7 / 1.0e-7 / 1.5e-8` | `6.7e-6 / 6.8e-6 / 6.8e-6 / 3.0e-6` |
| 1.4 | `2.2e-8 / 1.6e-7 / 1.6e-7 / 1.6e-7 / 2.2e-8` | `1.1e-5 / 1.2e-5 / 1.2e-5 / 3.0e-6` |
| 1.6 | `2.7e-8 / 2.0e-7 / 2.0e-7 / 2.0e-7 / 2.7e-8` | `1.8e-5 / 1.9e-5 / 1.9e-5 / 5.7e-6` |
| 2.0 | `4.8e-8 / 3.8e-7 / 3.8e-7 / 3.8e-7 / 4.7e-8` | `3.6e-5 / 3.7e-5 / 3.7e-5 / 3.9e-6` |

Worst nodes at `r ≈ 12.4–13.0 km` (crust); the fully independent chain is dominated by the Phase-4B
first-order gap (`~1e-5`). **Met** in both chains on every star. These values are diagnostics of this
run, not expected answers for anything.

## 13. Experiment F — near-vacuum identity (bound `1e-6`)

No profile node lies beyond `R_*`, so the vacuum statement is examined on the outermost nodes
(`ε < 1e9 g/cm³`, `r > R_*/2`): matter-source-corrected `m̂₀ + S + I²/r³` spread over δM̂ =
`7.2e-9 / 1.2e-8 / 1.2e-8 / 1.7e-8` (10/7/5/3 nodes); raw spread `4.3e-6 / 3.1e-6 / 2.1e-6 /
1.0e-6` — the crust still sources `m̂₀` at `ρ ~ 1e9`. Matching arithmetic exact; `s'(R_*) = 6I/R_*⁴`
to `1.6e-16`. **Met as testable**; the constancy beyond the surface is established analytically by
Experiments C/D.

## 14. Experiment J on DS(CMF)-1 — the finding (bound `1e-3`, **NOT MET**)

Homogeneous solve (sources off, `p̂₀*_c = 1`) on the 1.6 M☉ background; sequence derivative from
`M(ε_c ± 1e-3 ε_c)` with the solver's own central `p_c, ε_c` read back, at radial resolution
10000 / 20000 / 40000.

| res | `dM/dp_c` km³ | `(dε/dp)_c` sequence vs authority | `δM̂_hom` (ODE) | `(dM/dp_c) δp_c` | rel (`p_c` form) | rel (`ε_c` form) |
|---|---|---|---|---|---|---|
| 10000 | 9671.245 | `9.9e-7` | 6.003196 | 5.996177 | `1.171e-3` | `1.172e-3` |
| 20000 | 9671.245 | `9.9e-7` | 6.002267 | 5.996176 | `1.016e-3` | `1.017e-3` |
| 40000 | 9671.243 | `9.9e-7` | 6.002438 | 5.996175 | `1.044e-3` | `1.045e-3` |

Resolution-independent; the central-condition inversion is ruled out (`9.9e-7`). **Diagnosis,
proven test-side:**

1. *Integrability of the nodal `dε/dp` column.* `Σ ½(dε/dp)_{i,i+1} Δp` over the profile vs the
   profile's own `ε_* − ε_c`: whole star `2.4 % / 2.3 % / 2.3 %`, **crust (`ε < 1e14`) `17.9 % /
   16.8 % / 17.2 %`** at 10000 / 20000 / 40000. A sampled derivative that misses 17 % of the crust's
   density change at every resolution is missing *sub-node* density steps.
2. *The steps are in the table.* The DS(CMF)-1 crust table (1191 rows) has no vertical intervals,
   but its crust–core transition is a near-jump: `ε = 4.93e13 → 6.71e13 g/cm³` (**`Δε/ε = 36 %`**)
   over `p = 3.519e31 → 3.570e31` (**`Δp/p = 1.5 %`**), followed by `→ 8.39e13` over a 191 %
   pressure rise. In the star that jump is a layer `Δr ≈ Δp/[(ε+p)ν'] ≈ 0.6 m` thick — thinner
   than the node spacing (5 / 2.5 / 1.3 m), so the Steffen derivative sampled at the flanking nodes
   sees only the smooth parts. It is ≈ 21 % of the crust's `Δε`: the 17 % deficit.
3. *Stieltjes evaluation.* With `dp/dr = −(ε+p)ν'` the source `4πr²(ε+p)(dε/dp)p̂₀* dr` equals
   `4πr² p̂₀* (−dε)/ν'`. Summing it against the **column** reproduces the ODE's `m̂₀_hom(R_*)` to
   `2.2e-5 / 8.9e-6 / 1.4e-6` (the discretization is sound); summing it against the **profile's
   own `ε` steps** gives `5.99555 / 5.99582 / 5.99573`, which matches the independent sequence
   derivative to **`1.0e-4 / 5.9e-5 / 7.4e-5`** — inside the bound by an order of magnitude.
4. *Size on the sourced response.* The same evaluation applied to production's sourced `p̂₀*`
   gives the contribution production omits: **`4.78 % / 4.46 % / 4.58 %` of `δM̂`** (positive:
   production underestimates `δM̂`). The sign flip relative to the homogeneous case (`−1e-3`) is
   real: `p̂₀*_hom(R_*) < 0` (the star shrinks with `p_c`), `p̂₀*(R_*) = +32` for the sourced solution.
   This is the same ≈ 5 % that the retired-FD substitution (Experiment I) produced — the FD's
   "crust noise" of 4C-I0 was the step.

**Physics.** Hartle (1967) eq. (97) contains `dE/dP`; at a density discontinuity it is a delta
function, and its integral is a mass shell `4πr_i² Δε_i ξ₀(r_i)` — exactly the mechanism ADR-0007
P4/P6 adopted for the *surface* and did not adjudicate for *internal* discontinuities. On
DS(CMF)-1 the internal shell at the crust–core transition (and the smaller ones) is ≈ 4.6 % of
`δM̂`, ≈ 1 % of `ξ̂₀(R_*)`, negligible on `p̂₀*` in the core (the step is in the outermost crust).
Whether the transition is a genuine first-order transition or a matching artifact of the
"with crust" table is a decision for the ADR amendment; either way the tabulated EOS has it and
production's `δM̂` does not. **This is a substantive physical discrepancy of the accepted contract:
reported here, not repaired in 4D, ADR-0007 amendment required.** It also bounds Experiment H:
`δM̂`'s erratic order (2.35) is node placement relative to the same layer.

## 15. Experiment H — radial convergence at fixed `ε_c = 7.3125e14 g/cm³` (1.6 M☉)

| res | nodes | `R_*` km | `M` | `δM̂` | `ξ̂₀(R_*)` | `m̂₀(R_*)` | `p̂₀*(R/2)` |
|---|---|---|---|---|---|---|---|
| 5000 | 1319 | 13.463458 | 1.59997583 | 839.27834 | 1610.2622 | 828.84159 | 10.954706 |
| 10000 | 2635 | 13.468323 | 1.59997583 | 825.97015 | 1617.0515 | 815.54503 | 10.969607 |
| 20000 | 5268 | 13.471018 | 1.59997583 | 828.57691 | 1617.2796 | 818.15938 | 10.978080 |

Observed orders: `δM̂` **2.35**, `ξ̂₀(R_*)` 4.9, `m̂₀(R_*)` 2.35, `p̂₀*(R/2)` 0.81; Richardson
residual on `δM̂` at 20000 `7.7e-4`. **Reported** (INV-13: no pass criterion invented). At the
default resolution `δM̂` carries ≈ `3e-3` grid dependence, dominated by the crust and the surface
node; §14 explains the erratic order.

## 16. Experiment I — EOS-derivative sensitivity (§7 item 10: **not independently testable**)

| Source of `dε/dp` | `δM̂` | `ξ̂₀(R_*)` |
|---|---|---|
| A — Steffen authority (production) | 825.970148 | 1617.05153 |
| B — retired profile FD through the governed explicit-supply path (diagnostic) | 867.208633 | 1601.06671 |
| C — tabulated `c_s²` | **CONDITIONAL CHECK UNAVAILABLE** — the governed import has no such column; the raw `eos.thermo`'s three extra columns are undocumented (one negative, one zero, one equal to `n_B`) |

B−A: `δM̂` `5.0e-2`, `ξ̂₀` `9.9e-3`, profile `m̂₀` `5.8e-2` (worst at `R_*`). The ADR bound (`1e-3`)
applies to an independent source, which does not exist; the 5 % is the crust–core step (§14).

## 17. Experiment K — published results (bound `2e-2`)

### 17.1 Chandrasekhar & Miller (1974) Table I — homogeneous configurations

Transcribed from the 400-dpi journal scan (p. 73), all 20 rows, 4–5 significant figures; units
from the table footnote (`R_S = 2M`, `I` in `R_S³`, `ϖ₁` in `J/R_S³`, `ξ₀` in `J²/R_S³`, `δM/M` in
`J²/R_S⁴`). `R/R_S = 1.125` is the Buchdahl limit (`p_c` infinite) and is skipped; 19 rows compared
on exact interiors (`R = 12 km`, N = 4001). Conversions: `s(R) = ϖ₁ I_tab`; `ξ̂₀(R)/R_S³ = ξ₀ I_tab²`;
`δM/M = δM̂ R_S⁴/(M I²)`.

| Quantity | Worst relative disagreement (19 rows) |
|---|---|
| `I/MR²` | `3.7e-4` |
| `s(R) = ϖ₁ I` (first order) | `3.7e-4` |
| `ξ₀(R)` (second order, including its **sign change** between `R/R_S = 1.2` and `1.3`) | `7.3e-4` |
| `δM/M` in C&M's convention (§17.2) | `2.7e-4` |

**Met**, 19/19, at the tables' own precision (4th digit).

### 17.2 C&M omit the surface shell — a literature finding

C&M eq. (32) imposes that `m₀` "join continuously" with the exterior at the boundary. For a density
discontinuity that drops the shell `4πR²εξ₀(R)`. Their tabulated `δM/M` therefore equals
`[m₀(R⁻) + J²/R³]/M`, whose Newtonian limit is `IΩ²/M` — not the physical fixed-central-density
`Ω²R³/M` that Experiment D verifies (`δM̂/R³ → 1`). Compared like-for-like (shell excluded) the
agreement is `2.7e-4` at every compactness; the shell-inclusive value differs by a factor that grows
as `(5/2)(R/R_S)` (`0.68` at 1.15, `3.8` at 2.0, `493` at 100). The interior `m₀` integration is
thereby validated by a published relativistic calculation; the shell by C, D and J.

### 17.3 Hartle & Thorne (1968) Tables 3 and 5 — Harrison–Wheeler EOS

Fixture `tests/rotation/hartle_thorne_1968_hw_eos.hpp`: HT68 Table 1, all 72 rows `(P, E, ε)` in
`cm⁻²`, transcribed twice from the ADS scan at 300 dpi. **Transcription uncertainty:** one smudged
entry (the `ε` of the `P = 5.19E-29` row; leading digit illegible; log-log interpolation gives
`6.8E-29`, recorded as `6.82E-29`) — `ε` feeds only the baryon-density column, which neither the
TOV nor the Hartle equations consume, so it cannot affect any compared quantity; every `P, E` entry
is three significant figures. Conversion with HT68's own `Gc⁻² = 0.742×10⁻²⁸ cm/g` (their model; the
3-digit constant is `< 1e-3` on any compared quantity). HT68 interpolated "logarithmically between
table entries"; the fixture is densified **log-log linear to 40 points per decade** before import
so that CompactStar's Steffen spline follows *their* interpolant rather than inventing a different
EOS on a 2-per-decade table — their interpolant, not an improvement of it. CompactStar's cutoff
`max(1e-15 p_c, table floor)` omits an outer layer of metres on these models (shell/δM̂
`≈ 1e-8`).

| `E_c` | `R` km (tab) | `M/M☉` (tab) | `R_g/R` | `Ω` s⁻¹ | `ω_s/Ω` | `ω_c/Ω` | `δR/R` | `δM/M` |
|---|---|---|---|---|---|---|---|---|
| 1e14 | 35.882 (36.0) | 0.2649 (0.266) | 0.2559 (0.255) | 872 (867) | 1.43e-3 (1.41e-3) | 0.0621 (0.062) | 0.2452 (0.246) | 0.0528 (0.0522) |
| 3e14 | 20.797 (20.8) | 0.4035 (0.405) | 0.3553 (0.355) | 2440 (2430) | 7.23e-3 (7.19e-3) | 0.1183 (0.118) | 0.2005 (0.201) | 0.1288 (0.128) |
| 1e15 | 14.133 (14.2) | 0.5524 (0.554) | 0.4053 (0.404) | 5096 (5070) | 1.90e-2 (1.88e-2) | 0.2166 (0.217) | 0.1801 (0.181) | 0.1645 (0.163) |
| 3e15 | 10.171 (10.2) | 0.6594 (0.661) | 0.4341 (0.433) | 9119 (9100) | 3.61e-2 (3.59e-2) | 0.3498 (0.350) | 0.1630 (0.163) | 0.1626 (0.162) |
| 6e15 | 8.386 (8.41) | 0.6816 (0.684) | 0.4452 (0.444) | 12380 (12300) | 4.76e-2 (4.73e-2) | 0.4432 (0.443) | 0.1538 (0.154) | 0.1497 (0.149) |
| 1e16 | 7.469 (7.48) | 0.6661 (0.668) | 0.4451 (0.444) | 14560 (14500) | 5.22e-2 (5.20e-2) | 0.5026 (0.503) | 0.1520 (0.152) | 0.1371 (0.137) |
| 3e16 | 5.957 (5.96) | 0.5754 (0.577) | 0.4254 (0.425) | 19000 (19000) | 5.16e-2 (5.16e-2) | 0.6160 (0.616) | 0.1608 (0.161) | 0.1074 (0.107) |
| 1e17 | 5.138 (5.15) | 0.4609 (0.462) | 0.3810 (0.380) | 21230 (21200) | 3.85e-2 (3.81e-2) | 0.7110 (0.711) | 0.1874 (0.188) | 0.0824 (0.0818) |

Worst relative disagreements: `R` `4.7e-3`, `M` `4.1e-3`, `R_g/R` `3.6e-3`, `Ω` `6.8e-3`, `ω_s/Ω`
`1.3e-2` (a 3-s.f. quantity of order `1e-3`), `ω_c/Ω` `2.9e-3`, **`δR/R` `4.9e-3`**, **`δM/M`
`1.1e-2`**. **Met**, 8/8, within HT68's stated "about 1 per cent or better". The HW table as
printed has no internal density steps, so this comparison is blind to §14 by construction.

## 18. Experiment S — materialization (secondary)

`+100` and `−100 rad/s` bit-identical; `Q(2π·100) = 4 Q(π·100)` to `0.0`; zero spin exact zeros.
The 4C-I1 contract test (`hartle_monopole_contract`) re-runs green in both suites.

## 19. Detectors M1–M9

Each mutation applied at one production site of `RotationSolver.cpp` in a separate build tree,
the three 4D harnesses rebuilt and run (CMF in quick mode: G + F), the source restored and
verified byte-identical (SHA-256 `db7d213d…` before and after). Production diff after the sweep:
**NONE**.

| Detector | Fires on | Signature |
|---|---|---|
| M1 drop `e^{−2ν}` from term 3 | analytic B/A/C; published K1c/K1d/K2h; CMF G | analytic `4.2e-1`; C&M `ξ₀` `5.7`; HT68 `δR/R` `5.8e-2` |
| M2 flip the sign of term 3 | analytic; published K1c/K1d/K2h/K2i; CMF G | analytic `2.0`; HT68 `δM/M` `8.3e-2` |
| M3 omit `I²/R_*³` | analytic Bc/Ca; CMF matching arithmetic | `7.7e-3` (published blind: `< 2e-2`) |
| M4 omit the surface shell | analytic Bc/Cb/Cc/Da/Dc **only** (CMF and published blind — the shell is `1e-6`–`1e-8` there) | `δM̂/R³` intercept `1.0` instead of `0` |
| M5 profile-FD `dε/dp` | CMF G (analytic and HW blind: `dε/dp = 0` / smooth table) | `δM̂` `5–7e-2` |
| M6 impose `δp₀(R_*) = 0` by shooting | everything: analytic `1.5e4`, published, CMF | centre series `2e12` |
| M7 literal-zero start | analytic A only (node 0: the initial value itself) | downstream displacement at `r₀ = 1e-5 km` is `≈ 6e-13`, below every other bound — **as predicted by ADR-0007 §8**, this detector is weak at production `r₀`; it fires on the initial-value check, not on physics |
| M8 seed leak (raw `ω̄` for `s`) | analytic (13 lines), published, CMF (12 lines); **and** `hartle_monopole_contract` C5a/C5b/Sa (`4.0e10` at seed `1e3`) | §7 items 1–2 as the ADR expects |
| M9 drop `(1 + 8πr²p)` | analytic B/Cb/Eb; published K1c (`1.2`); CMF G (`2.5e-3`–`1.3e-2`) | — |

All nine fire on at least one 4D line; M8 additionally on the 4C-I1 contract (§7 items 1–2).

## 20. First-order protection

`git diff -- CompactStar/` = 0 lines. `tests/baselines/*.tsv` SHA-256 unchanged before and after,
in particular `hartle_I_dscmf1_debug.tsv = ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`;
`hartle_moment_inertia_cmf` and both Hartle contract tests green.

## 21. Suites

Run after every source file reached its committed state (CTest's own summary lines):

| Configuration | Result | Time | New tests |
|---|---|---|---|
| full (`COMPACTSTAR_EOS_DATA_ROOT` set) | **30/30 passed** (was 27/27) | 223.54 s (was 212.65 s) | `hartle_monopole_physics_analytic` 1.81 s, `hartle_monopole_published` 0.89 s, `hartle_monopole_physics_cmf` 4.22 s |
| self-contained | **16/16 passed** (was 14/14) | 18.41 s (was 13.70 s) | `hartle_monopole_physics_analytic` 1.84 s, `hartle_monopole_published` 0.88 s |

`git diff --check` clean; no scratch file in the tree; the detector build tree and every log live in
the session scratchpad.

## 22. Baseline decision and `GOVERNANCE.md` §3.1

**No `tests/baselines/hartle_monopole_dscmf1_debug.tsv` was created** and no regression
comparison registered: the status is not VERIFIED, and freezing a `δM̂` that omits a ≈ 4.6 % physical
term would enshrine it (§3.1 condition 2 in reverse). §3.1 state: **CORRECTION EXECUTED —
IMPLEMENTATION INDEPENDENTLY VERIFIED — PHYSICAL VALIDATION FAILED ON STEPPED CRUSTS — NO BASELINE
YET.** Condition 4's independent evidence now exists for the implementation (this record) and will
carry over to the amended contract; condition 7 remains deferred until that amendment is validated.

## 23. Status, invariants, Phase-4E disposition, next action

**Status: `HARTLE MONOPOLE VALIDATION FAILED`** — scoped exactly as in the banner. Every
implementation-level line of ADR-0007 §7 passed (items 3, 4, 5, 6-continuum, 7, 8, 12; 9 reported;
2 re-run); item 10 is not independently testable; **item 11 is NOT MET** on the tabulated crust and
the cause is a physical omission of ≈ 4.6 % in `δM̂`.

**INV-08:** `GOVERNED (ADR-0007 ACCEPTED) — O(Ω²) MONOPOLE SOURCE CONFORMED; IMPLEMENTATION
INDEPENDENTLY VERIFIED; PHYSICAL VALIDATION FAILED ON TABULATED CRUSTS WITH DENSITY STEPS (δM̂ ≈ −4.6 %
ON DS(CMF)-1) — ADR-0007 AMENDMENT REQUIRED`. `l = 2` remains not implemented, not validated; the
owner's decision that it blocks nothing stands. **No O(Ω²) `δM̂` from a stepped tabulated EOS may be
cited as a result.**

**INV-09:** not resolved; it inherits the same internal-step physics on the `δN_i` integrals.

**Phase 4E disposition.** The Phase-5 structural fields (`m₀/Ω²`, `δp₀/Ω²`, `ξ₀/Ω²`, `δM/Ω²`,
`ω̄/Ω`, `I`) are already carried by `HartleMonopoleResponse`/`PhysicalHartleMonopole` (4C-I1) and
are verified at the implementation level here; 4E is *not* satisfied, because the contract they
implement is physically incomplete on the project's EOS.

**Recommended smallest next action:** an **ADR-0007 amendment** adjudicating internal density
discontinuities of a tabulated EOS — either explicit internal shells `4πr_i²Δε_i ξ₀(r_i)` at
detected steps, or term 1 integrated against `dε` on the EOS table (the Stieltjes form of §14, which
needs no derivative at all) — followed by its implementation, a re-run of Experiments J, H and I on
DS(CMF)-1 plus the sourced Stieltjes check, and **then** the first monopole baseline. **Phase 5 must
not begin before that.**
