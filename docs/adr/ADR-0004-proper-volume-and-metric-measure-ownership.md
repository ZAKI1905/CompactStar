# ADR-0004 — Proper-volume measure and metric-factor ownership

| | |
|---|---|
| **Status** | **ACCEPTED — 2026-09-01** |
| **Date drafted** | 2026-09-01 |
| **Change class** | **Structural / architecture** (`GOVERNANCE.md:51`) |
| **Drafted at** | `0dd11a80ccaffedb093095aaea903e6e689e40c3` |
| **Roadmap increment** | Phase 3D (`docs/architecture/PHASE3_CONSOLIDATION_PLAN.md` §11) |
| **Governing invariants** | **INV-04** (proper-volume measure) — remains `VERIFIED CURRENT BEHAVIOR / LEGACY split`; **INV-03** (metric convention) |
| **Prerequisite** | Phase 3B — ADR-0003 **ACCEPTED**; provenance contract in place |
| **Accepted at** | `5eb6bdd314b34b1acccb2df59a92eb7eb813c145` (acceptance commit; implementation follows separately) |

> **ACCEPTED.** The owner adjudicated Q1, Q2 and Q3 (§0 below). This ADR is now **normative**,
> ranking directly below `GOVERNANCE.md`.
>
> The audit that produced the draft wrote no code: no production source, test, baseline or
> invariant status was changed by it. Implementation is Phase 3D-I and lands in a **separate
> commit after** this acceptance.
>
> **INV-04 is NOT marked fully resolved by this acceptance.** The canonical validated path is
> conformed; TOV Path 1, the scalar `NStar` accessor, `MixedStar` and the candidate scientific
> code remain governed-but-nonconformant. See §24 of the implementation record.

---

## 0. Owner decisions (binding)

Adjudicated 2026-09-01. The alternatives (§17), the measurements (§6, §7, §9) and the open
questions as originally posed (§22) are retained unchanged below, as the lifecycle requires.

### Q1 — ownership boundary: **OPTION B — a dependency-neutral mathematical primitive**

Three distinct owners, per §3:

| Role | Owner |
|---|---|
| **A. Mathematical** | a dependency-neutral CompactStar primitive owning `f(r,m)`, `Λ(r,m)`, `e^{Λ}(r,m)`, `w_V(r,m) = 4π r² e^{Λ}` **including domain/failure semantics** |
| **B. Cached representation** | `Physics::Evolution::GeometryCache` — retains `ExpLambda`, `WV`, `WVExpNu`, `WVExp2Nu` |
| **C. Consumer integrands** | `NStar`, the thermal drivers and others retain `n_B·w_V`, `c_V·w_V`, `Q_ν·w_V·e^{2ν}` **and their own unit conversions** |

**`Core` must NOT be made to depend on `Physics/Evolution` merely to obtain `w_V`.**

### Q2 — `MixedStar`: **GOVERNED NOW, SOURCE MIGRATION DEFERRED**

This ADR governs `MixedStar`'s future mathematical contract — including the requirement that its
metric factor use **total enclosed mass** (§15). `MixedStar` source is **not** modified in 3D. A
later task must establish focused coverage **first**; superficial tests must not be manufactured
to unlock the migration.

### Q3 — degenerate inputs: **HYBRID PHYSICAL-DOMAIN CONTRACT**

For finite `r`, `m`:

| Case | Behavior |
|---|---|
| `r == 0` **and** `m == 0` | **regular-center limit**: `f = 1`, `Λ = 0`, `e^{Λ} = 1`, `w_V = 0` |
| `r < 0` | **fail closed** |
| `r == 0` **and** `m != 0` | **fail closed** |
| `r > 0` **and** `f = 1 − 2m/r ≤ 0` | **fail closed** — includes the horizon/Schwarzschild-radius condition and any configuration outside the domain of the static stellar metric |
| `r > 0` **and** `f > 0` | evaluate normally, **no artificial clamp** |

Non-finite `r` or `m`: **fail closed**.

There is **no `1e-15` clamp** in the canonical primitive, and **no tolerance band around
`r = 0`** — the regular-center case is the exact `r == 0, m == 0` case. This task introduces no
new negative-mass policy; finite `m` is otherwise governed by the `f > 0` domain condition.

This is option **(d)** of §22-Q3 for the invalid cases, combined with the exact analytic limit at
the regular center. Per §22's supplementary note, adopting it inside `GeometryCache` is a
**scientific-semantic** change and not merely structural — which is why §4 of the implementing
brief classifies 3D as spanning both classes. The ordinary physical domain is untouched: measured
`max 2m/r = 0.481` across the four authenticated stars (§11), so **no validated production datum
reaches any fail-closed branch.**

---

## 1. Context

`SCIENTIFIC_INVARIANTS.md:151` (INV-04) records the canonical proper-volume measure as
`w_V(r) = 4πr² e^{Λ}` and states its proposed action outright at `:162`: *"Single owner. Structural
change; requires ADR."* `GOVERNANCE.md:73` fail-closed condition 5 — *"Uncertain ownership. Two
components compute or own the same quantity and no document says which is authoritative"* — is
therefore **active** for this quantity.

The Phase-3 entry audit (`PHASE3_CONSOLIDATION_PLAN.md` §5) proposed `GeometryCache` as the
canonical owner and phrased increment 3D as *"`GeometryCache` canonical; retire `NStar`/`MixedStar`
inline forms."* **This ADR audits that wording and finds it not architecturally sound as stated**,
for reasons recorded in §7 and §8. It proposes a different resolution that satisfies INV-04 more
completely and at lower behavior-preservation risk.

`ARCHITECTURE.md` calls `WV()` *"the universal radial integration measure."* It is not universal:
§4 below inventories **twenty-two** live or compiled sites that compute a proper-volume or
redshifted-proper-volume factor, of which exactly **three** read `GeometryCache`.

### 1.1 Baseline at drafting

Re-authenticated at `0dd11a8`, clean tree, upstream equal, 7 ahead of `master` (`d2040d8`), 0 behind.

| Configuration | Result |
|---|---|
| Full, authenticated `COMPACTSTAR_EOS_DATA_ROOT` | **13/13 PASS** |
| Self-contained (no data root) | **8/8** (5 external-data tests excluded by the CMake guard) |

All five protected artifacts hash-match the post-3C authenticated record
(`docs/validation/PHASE3C_BOLTZMANN_AUTHORITY.md` §14, rows at `:341-343` for the three
re-baselined thermal artifacts; `tov_dscmf1_reference.tsv` and `hartle_I_dscmf1_debug.tsv`
unchanged since 3A).

---

## 2. Relationship to INV-03 and INV-04

INV-03 (`SCIENTIFIC_INVARIANTS.md:135-145`) fixes the **metric convention**:
`ds² = −e^{2ν}dt² + e^{2Λ}dr² + r²dΩ²`, with `Λ = −½ ln(1 − 2m/r)` and the degenerate argument
clamped `denom ≤ 0 → 1e-15`.

INV-04 (`:149-162`) fixes the **measure built from it**: `w_V = 4πr² e^{Λ}`, with redshifted
variants `w_V e^{ν}` and `w_V e^{2ν}`.

These are two statements about one object. INV-03 says what `Λ` *is*; INV-04 says what is *built*
from it. This ADR proposes an owner for the pair and **does not alter either invariant's
statement**. It also does not alter INV-03's clamp — §11 records that the clamp is not in fact
uniform across implementations and routes the choice to the owner.

---

## 3. Three distinct meanings of "owner"

The Phase-3 wording collapses three questions that have different answers. This ADR keeps them
apart, and that separation is the substance of the proposal.

| | Question | Concerns |
|---|---|---|
| **A. Mathematical owner** | Who owns the *definition* `e^{Λ} = (1 − 2m/r)^{−1/2}` and `w_V = 4πr² e^{Λ}`, including domain, clamp and failure semantics? | one formula, one clamp, one degenerate-input rule |
| **B. Cached-representation owner** | Who owns the precomputed radial arrays `WV()`, `WVExpNu()`, `WVExp2Nu()`, `ExpLambda()`? | lifetime, provenance (ADR-0003), rebuild cost |
| **C. Consumer-specific integrand owner** | Who owns `n_B w_V`, `c_V w_V`, `Q_ν w_V e^{2ν}`, `ṅ w_V e^{ν}`? | the *physics* factor and its unit conversion — **not** the measure |

**A ≠ B.** A single mathematical owner does not require every caller to hold a `GeometryCache`.
**C is not the measure.** A baryon-number integrand contains `w_V`; it is not `w_V`, and its
`10^{54}` unit conversion is governed by INV-14, not by this ADR (§14).

---

## 4. Current implementation map

Complete, re-authenticated at `0dd11a8`. Classification vocabulary is the one required by the
Phase-3D audit brief.

### 4.1 Canonical measure and its producers

| # | Site | Expression | Class | Reachability | Coverage |
|---|---|---|---|---|---|
| 1 | `GeometryCache.cpp:211,223,241` | `m_area=(4π)r²` · `m_expLam=exp(Λ)` · `m_wV=m_area*m_expLam` | **PROPER VOLUME — canonical** | **LIVE** | `heat_capacity_v1`, `heat_capacity_real_star`, `passive_cooling_regression`, `grid_convergence_cmf`, `cache_contract` |
| 2 | `GeometryCache.cpp:244` | `m_wVExpNu = m_wV * m_expNu` | **REDSHIFTED PROPER VOLUME** | **COMPILED, UNEXERCISED** — no production consumer | none |
| 3 | `GeometryCache.cpp:247` | `m_wVExp2Nu = m_wV * m_exp2Nu` | **REDSHIFTED PROPER VOLUME** | **LIVE** — `NeutrinoCooling_Details.cpp:68` | `passive_cooling_regression` |
| 4 | `GeometryCache.cpp:119-147` `DeriveLambdaFromMR_(r,m,1e-15)` | `Λ = −½ ln(max(1−2m/r, 1e-15))`, `r ≤ 0 → denom = 1` | **MATHEMATICAL DEFINITION — fallback route** | LIVE code path, **never taken** on a profile built by `NStar` (§9) | indirect |
| 5 | `NStar.cpp:211-227` (`BuildFromTOV`) | same formula, same `1e-15` clamp, writes the `MetricLambda` column | **MATHEMATICAL DEFINITION — producer** | **LIVE** (TOV Path 2) | all Phase-2 harnesses |
| 6 | `NStar.cpp:702-721` (`Append`) | byte-for-byte the same block | **MATHEMATICAL DEFINITION — producer** | **LIVE** (TOV Path 1) | **none** |

### 4.2 Consumer-specific integrands carrying an inline metric factor

| # | Site | Expression | Class | Reachability | Coverage |
|---|---|---|---|---|---|
| 7 | `NStar.cpp:277-281` (`BuildFromTOV`) | `r².pow(2)` → `*= 4π n_B` → `/= (1−2m/r).sqrt()` → `*= FM3_TO_KM3` (`:23` = `1e54`) | **CONSUMER-SPECIFIC INTEGRAND** with inline proper volume | **LIVE** — produces `seq.b` (`:314-316`) | star built by every harness; **no golden artifact asserts `B`** |
| 8 | `NStar.cpp:562-569` (`FinalizeSurface`) | **identical algebra, identical operand order** | same | **LIVE** (TOV Path 1, via `TOVSolver::SurfaceIsReached()` `TOVSolver.cpp:2928-2931`) | **none** |
| 9 | `NStar.cpp:1053-1068` `BaryonNumIntegrand(double)` (decl. `NStar.hpp:393`) | `4π·r·r·n_B/√f`, **no `1e54`**; `r ≤ 0 → 0`, `n_B ≤ 0 → 0`, `f ≤ 0 → 0` | **CONSUMER-SPECIFIC INTEGRAND** + **separate INV-14 defect** | **COMPILED, UNEXERCISED** — zero callers | none |
| 10 | `MixedStar.cpp:238-248` (`SurfaceIsReached`) | visible `B`: `r².pow(2)` → `*= 4π ρ_vis` → `/= (1 − 2 m_tot/r).sqrt()` → `*= pow(10,54)`; `m_tot` = `mass_tot_dc` or its `GetSubSet` | **CONSUMER-SPECIFIC INTEGRAND**, two-sector metric | **COMPILED, UNEXERCISED** | none |
| 11 | `MixedStar.cpp:254-264` | dark `B`, same shape, same `m_tot` | same | COMPILED, UNEXERCISED | none |
| 12 | `MixedStar.cpp:659-670` (constructor, `:400`) | second copy of #10 | same | COMPILED, UNEXERCISED | none |
| 13 | `MixedStar.cpp:675-686` | second copy of #11 | same | COMPILED, UNEXERCISED | none |
| 14 | `MixedStar.cpp:1158-1169` `BaryonNumIntegrand(double)` | `r·r·4π ρ_vis / sqrt(1 − 2 m_tot/r)`; `1e54` **commented out** (`:1163-1164`); **no domain guard at all** | CONSUMER-SPECIFIC INTEGRAND + same missing-`1e54` shape | COMPILED, UNEXERCISED | none |
| 15 | `MixedStar.cpp:1173-1184` `BaryonNumIntegrand_Dark(double)` | same | same | COMPILED, UNEXERCISED | none |

### 4.3 Sites found beyond the documented three — and one correction

| # | Site | Expression | Class | Reachability |
|---|---|---|---|---|
| 16 | `Extensions/MixedStar/DarkCore_Analysis.cpp:97-100` | `4π r².pow(2)` → `*= (1 − 2M_tot/r).pow(-0.5)` → `× exp(ν)` | **PROPER VOLUME** *and* **REDSHIFTED PROPER VOLUME** | COMPILED, UNEXERCISED |
| 17 | `Microphysics/BNV/Analysis/BNV_Analysis.cpp:64-67` | identical shape | PROPER VOLUME + REDSHIFTED | COMPILED, UNEXERCISED · **CANDIDATE** |
| 18 | `BNV/Analysis/Decay_Analysis.cpp:482-484`, `:577-580` | `4πr² (1−2M/r)^{-0.5} e^{ν}` | REDSHIFTED PROPER VOLUME | COMPILED, UNEXERCISED · CANDIDATE |
| 19 | `BNV/Internal/BNV_Chi.cpp:2035-2038`, `:2173-2176`, `:2331-2334` | `4πr²·rate/√(1−2M/r)·e^{ν}` (each preceded by its own `1e54`) | REDSHIFTED PROPER VOLUME | COMPILED, UNEXERCISED · CANDIDATE |
| 20 | `BNV/Channels/BNV_B_Chi_Photon.cpp:305-311`, `:348-350` | `4πr²/√(1−2M/r)·e^{ν}` | REDSHIFTED PROPER VOLUME | COMPILED, UNEXERCISED · CANDIDATE |
| 21 | `BNV_B_Chi_Photon.cpp:1216-1219`, `:1864-1867` | `4πr²/√(1−2M/r)·e^{2ν}` | REDSHIFTED PROPER VOLUME (`e^{2ν}` variant) | COMPILED, UNEXERCISED · CANDIDATE |
| 22 | `Physics/Evolution/RotochemicalCache.cpp:30`, `:65` | reads `geo.WV()` | **ALREADY-CONFORMANT CONSUMER** | **NOT COMPILED** · CANDIDATE |
| 23 | `TOVSolver.cpp:1976-1980` | `4π n_B r².pow(2)·(1−2M/r).pow(-0.5)`, then `× exp(ν)` | PROPER VOLUME | **PREPROCESSOR-DISABLED** — inside `#if DarkCore_Analysis` with `#define DarkCore_Analysis 0` at `TOVSolver.cpp:33` |

> **Correction to `PHASE3_CONSOLIDATION_PLAN.md` §5.** That table records
> `DarkCore_Analysis.cpp:97` as *"`4π r²` — coordinate volume, correctly not proper."* Line `:97`
> alone is `4π r²`, but line `:98` divides by `√(1 − 2M_tot/r)`. Site #16 **is** a proper-volume
> weight, and `:100` builds the redshifted variant. The Phase-3-entry classification of that row
> is superseded by this ADR.

### 4.4 Expressions that share the token `4πr²` but are NOT this measure

| Site | Quantity | Class |
|---|---|---|
| `RotationSolver.cpp:302,309,319,329,339,349` | `16π(p+ε)/(1−2m/r)`, `4π(p+ε)r/(1−2m/r)` | **HARTLE ODE COEFFICIENT — NOT PROPER VOLUME** (LIVE, O(Ω)) |
| `RotationSolver.cpp:1054-1097` | O(Ω²) source terms | HARTLE ODE COEFFICIENT — NOT PROPER VOLUME (UNREACHABLE SCAFFOLDING · CANDIDATE) |
| `TOVSolver::ODE` mass equation | `4πr²ρ` | **COORDINATE AREA — NOT PROPER VOLUME** |
| `NStar.cpp:825-831` (`EvaluateNu` surface BC), `MixedStar.cpp:781` | `1 − 2M/R`, clamp `1e-15`, `ν_R = ½ ln x` | **OUT OF SCOPE** — ν boundary condition, not a volume measure |
| `SurfaceGravity.cpp:205` | `c² M_cm/(R_cm² e^{ν})` | **OUT OF SCOPE** — uses `e^{ν}`, never `e^{Λ}` |

**Sharing the token `4πr²` does not make an expression the proper-volume measure.** Nothing in
§4.4 may be swept into a consolidation.

---

## 5. Mathematical definition proposed for single ownership

```
    f(r, m)   =  1 − 2 m / r                     [metric denominator, dimensionless]
    Λ(r, m)   =  −½ ln f                         [INV-03]
    e^{Λ}     =  f^{−1/2}                        [g_rr = e^{2Λ}]
    w_V(r)    =  4π r² e^{Λ}                     [INV-04, km² per km of coordinate radius]
```

with `r`, `m` both in geometric length units (km), per `CURRENT_ARCHITECTURE.md` §1
(`m [km] = GM/c²`). Domain and failure semantics are §11 and are **not** settled by this ADR.

---

## 6. Measured equivalence

All measurements were made with a **scratch program that links the built library and changes
nothing** (`NStar::SolveTOV_Profile` on `DS(CMF)-1_with_crust`, the same recipe as
`tests/thermal/heat_capacity_real_star.cpp:132-142`). No production or test file was modified and
no baseline was regenerated. Stars: 1.0, 1.4, 1.6, 2.0 M☉ (achieved 1.000044, 1.400022, 1.599976,
2.000095), `N` = 2629 / 2646 / 2635 / 2527 radial nodes.

### 6.1 Metric factor — algebraic forms compared node by node

`A = 1/√f` · `B = exp(Λ_stored)` · `C = exp(−½ ln f)` · `D = pow(f, −½)`

| Star | `A` vs `B` | `A` vs `C` | **`B` vs `C`** | `A` vs `D` |
|---|---|---|---|---|
| 1.0 M☉ | max_rel 2.220e-16, rms 8.21e-17, **1 ULP**, 2236/2629 bitwise | identical to the `A`–`B` column | **0.0 — 2629/2629 bitwise** | 2.220e-16, 1 ULP |
| 1.4 M☉ | 2.220e-16, rms 8.30e-17, 1 ULP, 2223/2646 | " | **0.0 — 2646/2646** | 2.220e-16, 1 ULP |
| 1.6 M☉ | 2.220e-16, rms 7.79e-17, 1 ULP, 2250/2635 | " | **0.0 — 2635/2635** | 2.220e-16, 1 ULP |
| 2.0 M☉ | 2.220e-16, rms 7.88e-17, 1 ULP, 2111/2527 | " | **0.0 — 2527/2527** | 2.220e-16, 1 ULP |

**Core vs surface.** Worst-case relative difference is *larger in the core*
(2.220e-16 for `r ≤ 0.5R`) than at the surface (1.63e-16 – 1.96e-16 for `r ≥ 0.9R`), and the
worst node sits at `r ≈ 1e-4 km` — the innermost decade, where `f → 1` and the two forms round
`1/√f` differently near unity. The surface is **not** the sensitive region: `min f` is
0.762 / 0.678 / 0.633 / 0.519 across the four stars, so the metric factor never exceeds 1.39.

**`B ≡ C` bitwise on every node of every star** is not a coincidence: `NStar.cpp:181` computes
`m_km` once and uses that same value both for the stored `m(km)` column (`:183`) and for the
`Λ` column (`:211-227`), with exactly the expression `−0.5·std::log(denom)`. On this construction
path, "stored `Λ`" *is* "derived `Λ`".

### 6.2 Controlled analytic star

Uniform-density (Schwarzschild interior) `m(r) = M (r/R)³`, 20 001 nodes, compactness
`2M/R ∈ [0.05, 0.98]`:

| `2M/R` | `A` vs `C` | `A` vs `D` |
|---|---|---|
| 0.05 – 0.60 | max_rel 2.220e-16, **max 1 ULP** | 2.220e-16, 1 ULP |
| 0.70 | 2.579e-16, 2 ULP | 2.220e-16, 1 ULP |
| 0.80 | 2.482e-16, 2 ULP | 2.220e-16, 1 ULP |
| 0.90 | 3.175e-16, 2 ULP | 2.220e-16, 1 ULP |
| 0.98 | 3.067e-16, 2 ULP | 2.220e-16, 1 ULP |

The `exp(−½ ln f)` route degrades from 1 to 2 ULP as `f → 0`; `pow(f, −½)` does not. Neither
exceeds 2 ULP even at `2M/R = 0.98`, four times more compact than any physical star in §6.1.

### 6.3 Weight-level comparison

`w_V` from `GeometryCache` versus the scalar inline product `4π·r·r·(1/√f)`:

| Star | max_rel | rms_rel | max ULP | bitwise |
|---|---|---|---|---|
| 1.0 M☉ | 5.922e-16 | 1.288e-16 | **3** | 1504/2629 |
| 1.4 M☉ | 5.515e-16 | 1.302e-16 | 3 | 1524/2646 |
| 1.6 M☉ | 6.378e-16 | 1.288e-16 | 3 | 1553/2635 |
| 2.0 M☉ | 6.378e-16 | 1.279e-16 | 3 | 1488/2527 |

**Pointwise, the canonical weight and the inline product are not bit-identical — they agree to
3 ULP.** This confirms the floating-point warning in `PHASE3_CONSOLIDATION_PLAN.md` §5 and
quantifies it for the first time.

---

## 7. Integrated baryon-number sensitivity, and the predeclared tolerance

### 7.1 Derivation *before* measurement

The two candidate compositions evaluate the same exact real number and differ only in the sequence
of correctly-rounded IEEE-754 binary64 operations:

```
 V1 (current, NStar.cpp:277-281)   r.pow(2) → ×(4π n_B) → ÷ f.sqrt() → ×1e54
 V2 (canonical)                    (4π)·r.pow(2) → ×exp(Λ) → ×n_B → ×1e54
```

The chains differ by at most **six** correctly-rounded operations (the `pow`/`sqrt`/`log`/`exp`
calls and the intervening products). Each contributes at most `½ ulp = 2^{−53} ≈ 1.11e-16`
relative. Worst-case **pointwise** discrepancy: `6 × 1.11e-16 ≈ 6.7e-16`.

The integral is a **positive-weighted sum of non-negative terms** — `n_B ≥ 0`, `w_V > 0`, and the
trapezoid weights of `DataSet::Integrate` are positive (INV-13: the integrator is the exact
integral of a *linear* interpolant). For such a sum the relative perturbation of the total is
bounded by the maximum relative perturbation of the terms; and since both variants use the same
quadrature on the same grid, the summation's own rounding is common-mode to first order:

```
 |ΔB| / B  ≤  6.7e-16 × (1 + O(N·u))  with N ≈ 2.6e3, N·u ≈ 2.9e-13   ⇒  ≈ 6.7e-16
```

**PREDECLARED TOLERANCE: `|ΔB|/B ≤ 1.0e-15` on `SeqPoint::b`.** The ~1.5× margin over the derived
`6.7e-16` covers `std::log` and `std::exp`, which IEEE-754 does not require to be correctly
rounded and whose platform implementations may cost an extra ULP each (§6.2 shows exactly that at
high compactness).

This bound is derived from the algebra and the operation counts alone. It does not depend on any
measurement, and it is recorded here *before* any implementation exists.

### 7.2 Measured — consistent with the bound, far inside it

Every variant was integrated through the **same** `DataSet::Interpolate`/`Integrate` calls as
production. The V1 replication reproduced `seq.b` **bitwise on all four stars**, which
authenticates the replication.

| Star | `seq.b` | V2 − V1 (canonical `w_V·n_B·1e54`) | V2b − V1 (`(n_B·1e54)·w_V`) | V3 − V1 (derived-Λ geometry) | V4 − V1 (`pow(f,−½)`) | V5 − V1 (`×exp(Λ_stored)`) |
|---|---|---|---|---|---|---|
| 1.0 M☉ | 1.27388873109354535e+57 | **1.368e-16** (1 ULP) | 1.368e-16 | 1.368e-16 | 1.368e-16 | 1.368e-16 |
| 1.4 M☉ | 1.83218336257875150e+57 | **0** (0 ULP) | 0 | 0 | 0 | 0 |
| 1.6 M☉ | 2.12457569547972117e+57 | **0** | 0 | 0 | 1.640e-16 (1 ULP) | 1.640e-16 |
| 2.0 M☉ | 2.74576306114479063e+57 | **0** | 0 | 0 | 0 | 0 |

**Maximum observed `|ΔB|/B` = 1.64e-16 — a single ULP of the double result, and 6× inside the
predeclared `1.0e-15`.** The integral is far less sensitive than the pointwise weight (§6.3, 3
ULP) because the trapezoidal sum averages the ±1 ULP jitter over ~2600 nodes.

**Consequence for the Phase-3-entry expectation.** `PHASE3_CONSOLIDATION_PLAN.md` §10 predeclared
*"NEW TOLERANCE PREDECLARED"* for the `NStar` baryon path. That expectation is **correct and
necessary**: the change is not bit-identical (the 1.0 M☉ star moves by 1 ULP). The tolerance is
now derived and stated: **1.0e-15 relative on `seq.b`**.

---

## 8. Layering analysis

### 8.1 Measured dependency directions at `0dd11a8`

Include edges between top-level modules, derived by scanning every `#include "CompactStar/…"`:

```
   EOS        ──→ Core (Prog.hpp)        ──→ Units.hpp
   Core       ──→ EOS                    ──→ Physics/State (Pulsar.hpp:49-51)
              ──→ Physics/Driver/Spin    (Pulsar.cpp:21-22)
              ──→ Physics/Spin.hpp       (Pulsar.cpp:22)
              ──→ Extensions/MixedStar   (TaskManager.cpp:12)
   Extensions ──→ Core                   (DarkCore_Analysis.hpp:44, LightDM_Scalar_Density.hpp:47)
   Physics/Evolution ──→ Core            *** .cpp ONLY ***  (StarContext.cpp:9,
                                          RotochemicalCache.cpp:14-15)
   Physics/Driver    ──→ Physics/Evolution, EOS, Units.hpp
```

Two facts matter and neither was previously recorded:

1. **`Core → Physics` and `Core ↔ Extensions` edges already exist.** So a `Core → Evolution`
   include would not be the first inversion in the tree.
2. **`Physics/Evolution`'s *headers* are deliberately Core-free.** `StarContext.hpp:53`
   forward-declares `class StarProfile`; `ProfileProvenance.hpp:47` does the same. Only `.cpp`
   files reach into Core. That is a property worth preserving, and Option D (§9) would destroy it.

Everything links into **one** static library (`CMakeLists.txt:122`), so no option here creates a
*link* cycle. The question is logical layering, not linkage.

### 8.2 Would `Core::NStar` / `Core::MixedStar` consuming `Evolution::GeometryCache` introduce a bad dependency?

**Compilation: no.** `GeometryCache.hpp` includes only `Zaki/Vector/DataColumn.hpp` and
`ProfileProvenance.hpp`; neither reaches Core. `NStar.cpp` could include it today.

**Architecture and semantics: yes — four independent obstructions, three of them fatal.**

1. **Construction requires a `StarContext`, hence a `StarProfile`.** `GeometryCache.hpp:104` —
   the only constructor is `explicit GeometryCache(const StarContext &)`. `MixedStar` has **no
   `StarProfile`**: it stores two raw `Zaki::Vector::DataSet`s, `ds_vis` and `ds_dar`
   (`MixedStar.hpp:269,272`). Sites #10–#15 therefore cannot consume `GeometryCache` at all
   without first migrating `MixedStar` to `StarProfile` — a change far outside proper-volume
   consolidation. **Fatal for `MixedStar`.**

2. **`NStar` builds its integrand *inside an open `EditScope`*.** `BuildFromTOV` opens
   `auto edit = prof_.Edit();` at `NStar.cpp:90` and holds it through the integrand construction
   at `:277-281`; `FinalizeSurface` does the same at `:533` through `:562-569`. A `StarContext`
   constructed there would bind columns (`StarContext.cpp:184-203`) and snapshot provenance
   (`GeometryCache.cpp:28`) at a **revision that is about to be bumped** when the `EditScope`
   destructs (`StarProfile.hpp:162-166`). Under **ADR-0003** that geometry is stale the instant it
   is created, and its provenance names a profile revision that never existed as a committed
   state. **Fatal, and it is an accepted-ADR conflict, not a style preference.**

3. **`GeometryCache` requires `ν`, which the proper-volume measure does not.**
   `GeometryCache.cpp:178-179` throws if the `nu` column is missing or mis-sized. `w_V = 4πr² e^{Λ}`
   contains no `ν`. Making `GeometryCache` the *mathematical* owner therefore couples the measure
   to a redshift input it does not use, and forbids any consumer that legitimately has `(r, m)`
   but no `ν`.

4. **`GeometryCache` has no scalar evaluate-at-`r` API.** Site #9 (`NStar::BaryonNumIntegrand`)
   and sites #14–#15 are *pointwise interpolated* queries at arbitrary `r`, not node lookups.
   `GeometryCache` exposes only whole `DataColumn`s (`GeometryCache.hpp:141-163`). Serving them
   would require a new interpolating API on a class whose documented contract is
   *"an immutable snapshot"* of node arrays.

**Answer to the audit's load-bearing question: yes.** Making Core consume
`Physics::Evolution::GeometryCache` introduces a reversed module dependency *and*, more decisively,
is blocked on its own terms by (1), (2) and (3). **The Phase-3-entry wording
`GeometryCache canonical; retire NStar/MixedStar inline forms` cannot be implemented as written.**

---

## 9. `GeometryCache` stored-Λ vs derived-Λ — measured

`GeometryCache.cpp:196-206` selects the profile's `Λ` column when present and otherwise derives it
from `(m, r)`. The audit asked: **can two `GeometryCache` instances over the same profile data
produce different `WV` depending on which route ran?**

Measured by building, for each star, a second `StarProfile` copy with the `MetricLambda` column
index cleared (`StarProfile.hpp:484-527`), so `StarContext::Lambda()` returns `nullptr`
(`StarContext.cpp:197`) and the `DeriveLambdaFromMR_` branch is forced:

| Quantity | 1.0 M☉ | 1.4 M☉ | 1.6 M☉ | 2.0 M☉ |
|---|---|---|---|---|
| `Λ` stored vs derived | **bitwise, 2629/2629** | 2646/2646 | 2635/2635 | 2527/2527 |
| `ExpLambda()` | bitwise | bitwise | bitwise | bitwise |
| `WV()` | **bitwise** | **bitwise** | **bitwise** | **bitwise** |
| `WVExpNu()` | bitwise | bitwise | bitwise | bitwise |
| `WVExp2Nu()` | bitwise | bitwise | bitwise | bitwise |

**Answer: no — on any profile produced by the `NStar` TOV path, the two routes are bit-identical.**
The reason is structural, not luck: `NStar.cpp:181,211-227` writes a `Λ` column computed by exactly
the expression `DeriveLambdaFromMR_` uses, from exactly the `m` and `r` values it stores.

**Scope of that answer.** It is a property of `NStar`'s construction, **not** a guarantee of the
contract. A `StarProfile` loaded from a file (`StarProfile.cpp:63-64` reads a `lambda(r)` column;
`StarBuilder` is LIVE on the file-reading path) could carry a `Λ` produced elsewhere, at a
different precision or from a different `m`. The proposed contract (§10) therefore keeps the
route-selection rule explicit rather than relying on the measured equality.

**ADR-0003 is not changed by this ADR.** The provenance contract already forbids mixing a
`GeometryCache` with a context it was not built from (`StarContext.cpp:763-771`); the equality
above adds evidence, not a new rule.

---

## 10. Proposed ownership

### 10.1 Proposed mathematical owner — a dependency-neutral geometry primitive

A header-only, module-neutral primitive owning exactly:

```
    MetricDenominator(r, m)         →  f = 1 − 2m/r  with the governed domain rule
    ExpLambda(r, m)  /  Lambda(r, m)
    ProperVolumeWeight(r, m)        →  4π r² e^{Λ}
```

**Proposed home: the `CompactStar/` top level**, alongside `CompactStar/Units.hpp`. The repository
layering dictates this, and the precedent is explicit and recent — `Units.hpp:35-39`, added by
**Phase 3A of this same program**, records the reasoning verbatim:

> *"Layering. This header includes nothing … It therefore adds no edge to the dependency graph and
> may be included from any layer. That is why it sits at the top level rather than inside `EOS` or
> `Physics`."*

A geometry primitive at the same level is includable by `Core`, `EOS`, `Physics`, `Extensions` and
`Microphysics` **without creating a single new inter-module edge**, and without disturbing the
Core-free `Physics/Evolution` headers (§8.1 fact 2). Final file name and namespace are left to the
owner; `CompactStar/Geometry.hpp` under `namespace CompactStar::Geometry` is the natural analogue
of `CompactStar::Units`. Unlike `Units.hpp` it will need `<cmath>`; that is a standard-library
include, not a dependency edge.

### 10.2 Proposed cached-representation owner — `GeometryCache`, unchanged

`GeometryCache` remains the **canonical cached owner** of `ExpLambda`, `WV`, `WVExpNu`, `WVExp2Nu`
and the mixed metric products, and remains the ADR-0003 provenance-carrying immutable snapshot.
It changes only in *where the formula comes from*: `DeriveLambdaFromMR_` delegates to the
primitive. Its `DataColumn` composition (`GeometryCache.cpp:211,223,241,244,247`) is **kept
verbatim**, which is what makes bit-identity for its consumers achievable by construction rather
than by measurement (§13).

### 10.3 Why this is the better reading of "single owner" (§12 of the audit brief)

INV-04's requirement is that the *measure* have one definition. Forcing every caller to instantiate
a `GeometryCache` is a different and stronger requirement, and §8.2 shows it is unattainable for
`MixedStar` and unsafe for `NStar`. Splitting A from B satisfies INV-04 **more** completely than
Option A does, because it can reach the six `MixedStar` sites and the scalar accessors that
Option A structurally cannot.

---

## 11. Domain and failure semantics — **NOT settled by this ADR**

INV-03 records *one* clamp: `denom ≤ 0 → 1e-15`. The audit found **six mutually inconsistent
behaviors** for the same degenerate input, all live or compiled today.

| Behavior | Site | `r ≤ 0` | `f = 0` | `f < 0` |
|---|---|---|---|---|
| **B1** clamp `f` to `1e-15` | `GeometryCache.cpp:129-144`; `NStar.cpp:211-227`, `:702-721` | `denom = 1` ⇒ `e^{Λ} = 1` | `e^{Λ} = 3.162278e+07` | `e^{Λ} = 3.162278e+07` |
| **B2** clamp with a logged error | `NStar.cpp:826-830` (`EvaluateNu` BC) | n/a | `x = 1e-15`, `Z_LOG_ERROR` | same |
| **B3** **silent divisor substitution** | every inline `/= (…).sqrt()` — `NStar.cpp:280`, `:567`; `MixedStar.cpp:243,245,259,261,664,666,680,682`; `BNV_Chi.cpp:2037,2175,2333`; `BNV_B_Chi_Photon.cpp:306,310,349,1217,1865` | `f` → `−inf`/NaN | **divides by 1** — i.e. `e^{Λ} = 1`, flat space — with only `Error(1) in [DataColumn.cpp:operator/=:1070] … is zero, dividing by 1 instead!` | **NaN** |
| **B4** fail-to-zero | `NStar.cpp:1055-1065` | returns `0` | returns `0` | returns `0` |
| **B5** `pow(f, −0.5)`, no guard | `DarkCore_Analysis.cpp:98`; `BNV_Analysis.cpp:65`; `Decay_Analysis.cpp:483,578`; `TOVSolver.cpp:1977` | ±inf | `+inf` | **NaN** |
| **B6** no guard at all | `MixedStar.cpp:1166,1181` | inf/NaN | `+inf` | **NaN** |

**B1 and B3 disagree by a factor of 3.16e7 at `f = 0`, and B4 disagrees with both by returning
zero.** B3 was measured directly: `Zaki::Vector::DataColumn::operator/=` substitutes `1` for a
zero divisor and logs, so the inline sites do not produce `inf` at `f = 0` — they **silently drop
the metric factor at that node**. This is a *library* behavior, not a CompactStar one, and it is
inherited by nineteen call sites that never asked for it.

**No current output moves under any choice.** Measured `min(1 − 2m/r)` across the four authenticated
stars is 0.762 / 0.678 / 0.633 / 0.519, so `max 2m/r = 0.481`; **zero nodes** in any star reach the
degenerate branch. The clamp is unreachable on every validated path.

The candidate contracts, none of which the evidence selects:

| Contract | Consequence |
|---|---|
| **clamp `f` to `ε`** (status quo B1) | finite, large, physically meaningless weight; preserves `GeometryCache` exactly; hides the input error |
| **return a finite regularized value** | as above with an explicit, documented `ε` rather than an inherited magic number |
| **return zero** (B4) | the shell contributes nothing; silently under-counts `B` |
| **throw / fail closed** | matches `GOVERNANCE.md` §3 in spirit; changes `GeometryCache` from total to partial and can abort a run mid-integration |

**This is an owner/ADR decision (Q3, §21), not an implementation detail.** The implementation must
not pick one by copying whichever site it happens to touch first.

---

## 12. Redshifted variants

Proposed contract, and the audit found no component with a competing claim:

- The **primitive owns `w_V` only**. It does not know about `ν`.
- **`GeometryCache` composes and caches** `w_V e^{ν}` (`:244`) and `w_V e^{2ν}` (`:247`), because
  the redshift columns are already its business (`:214-221`) and ADR-0003 already governs their
  provenance.
- **Consumers own their own physics factor**, never a second copy of the measure.

The nineteen candidate sites in §4.3 each re-derive `w_V e^{ν}` or `w_V e^{2ν}` inline. Under this
contract their future obligation is to obtain `w_V` from the single owner and apply their own
`e^{ν}`/`e^{2ν}` — **not** to invent a third canonical measure. `WVExpNu()` currently has **no
production consumer** and must not be deleted on that basis while it is the designated composition
point for those sites.

---

## 13. Preservation requirements

Per-target, as required by the audit brief. **No blanket tolerance.**

| Migration target | Standard | Basis |
|---|---|---|
| `GeometryCache::WV/WVExpNu/WVExp2Nu/ExpLambda` and their consumers — `StarContext.cpp:824`, `NeutrinoCooling_Details.cpp:68` | **BIT-IDENTICAL REQUIRED** | the `DataColumn` composition at `GeometryCache.cpp:211,223,241,244,247` is kept verbatim; only `Λ`'s *source* is delegated |
| `GeometryCache::DeriveLambdaFromMR_` fallback | **BIT-IDENTICAL REQUIRED** | primitive reproduces `GeometryCache.cpp:129-144` exactly, including `r ≤ 0 → denom = 1` and `ε = 1e-15` — unless §11/Q3 changes it, in which case the change is scientific-semantic and needs its own evidence |
| `NStar.cpp:277-281` — `BuildFromTOV` baryon integrand | **PREDECLARED TOLERANCE — `\|ΔB\|/B ≤ 1.0e-15` on `SeqPoint::b`** | §7.1 derivation; §7.2 measured max 1.64e-16 |
| `NStar.cpp:562-569` — `FinalizeSurface` baryon integrand | **PREDECLARED TOLERANCE — same 1.0e-15 — *plus* new coverage first** | textually and numerically identical algebra to the above, but on the **unvalidated** TOV Path 1; §4.2 #8 has zero coverage |
| `NStar.cpp:1053-1068` — scalar `BaryonNumIntegrand` | **NOT YET MIGRATABLE** | zero callers and carries the separate INV-14 `1e54` defect (§14) |
| `MixedStar.cpp` sites #10–#15 | **NOT YET MIGRATABLE** | no `StarProfile`, zero coverage, two-sector mass semantics (§15) |
| `DarkCore_Analysis`, `BNV_*`, `Decay_Analysis` (§4.3 #16–#21) | **NOT YET MIGRATABLE** | `GOVERNANCE.md` §5 candidates; contract only (§16) |
| `RotochemicalCache.cpp:30,65` | **no change** — already conformant, NOT COMPILED | reads `geo.WV()` today |
| §4.4 — Hartle ODE coefficients, TOV mass ODE, `EvaluateNu` BC, `SurfaceGravity` | **OUT OF SCOPE — must not change** | not this measure |

---

## 14. Relationship to INV-14 — and an explicit guard against scope creep

**INV-14 governs `B = ∫ 4πr² n_B (1−2m/r)^{−1/2} dr × 10^{54}` in full — the baryon-density
semantics *and* the unit conversion. This ADR governs `dV` only.**

Three statements are normative for any implementation of this ADR:

1. **Baryon-density semantics remain INV-14 and ADR-0001.** Any *per-species* integral must apply
   `n_i = Y_i n_B`; this ADR does not touch that.
2. **`FM3_TO_KM3` / `1e54` is not part of the proper-volume definition.** It appears at
   `NStar.cpp:23,281,569`, `MixedStar.cpp:248,264`, `RotochemicalCache.cpp:22,43`,
   `BNV_Chi.cpp:2031,2170,2328`, `BNV_B_Chi_Photon.cpp:302-303,1213,1861`. The proposed primitive
   must **not** own it. It is a consumer-side unit conversion (§3, meaning C), and
   `Units.hpp:33` already records `fm^-3 <-> km^-3` among the factors deliberately excluded from the
   3A consolidation.
3. **The zero-caller scalar accessor's missing `1e54` is a separately recorded defect and MUST NOT
   be repaired under this ADR.** `NStar::BaryonNumIntegrand` (`NStar.cpp:1053-1068`,
   `NStar.hpp:393`) computes the same geometric expression with no `1e54`, giving fm⁻³·km² rather
   than km⁻¹ — recorded at `SCIENTIFIC_INVARIANTS.md:407-409` as a ⚠ Defect under INV-14. It is
   harmless only because it has zero callers.

> **Instruction to any implementing agent.** Adding the missing `1e54` while "centralizing the
> proper-volume factor" would be a **scientific-semantic change smuggled inside a structural one**,
> in violation of `GOVERNANCE.md` §2 (*a change spanning classes takes the strictest applicable
> requirement*) and `AGENTS.md` §7. If the site is touched at all, its arithmetic must be left
> exactly as it is. Repairing it is a separate change under INV-14, requiring its own adjudication
> and its own evidence.

The same shape exists in `MixedStar::BaryonNumIntegrand` / `_Dark` (`MixedStar.cpp:1163-1164`,
`:1178-1179`, where the `1e54` is commented out). It is recorded here for the first time and
carries the same prohibition.

---

## 15. `MixedStar` disposition

**Classification: `MIXEDSTAR SHOULD USE THE SHARED PRIMITIVE BUT NEEDS NEW COVERAGE FIRST`.**

Findings:

- **The scalar primitive *can* serve it cleanly.** `w_V(r, m)` is a pure function of `(r, m)`;
  `MixedStar` simply supplies the **total** enclosed mass. There is no architectural obstruction —
  in contrast to `GeometryCache`, which cannot serve it at all (§8.2 obstruction 1).
- **`GeometryCache` cannot serve it without migrating `MixedStar` to `StarProfile`.**
  `MixedStar.hpp:269,272` — two raw `DataSet`s, no profile, no version, no provenance.
- **Its geometry is genuinely two-sector.** `mass_tot_dc` is assembled on a master grid chosen by
  which component extends further (`MixedStar.cpp:169-231`, `:596-653`), with the shorter sector's
  mass held at its surface value beyond its own radius (`:191-197`, `:221-227`). The metric factor
  is then evaluated on the **visible** grid for the visible integrand and the **dark** grid for the
  dark one, using `GetSubSet` when the grids differ in length (`:245`, `:261`, `:666`, `:682`).
  `GeometryCache` assumes one `(r, m, ν, Λ)` grid and has no representation for this.
- **Total mass is the physically correct choice and must be preserved.** `SCIENTIFIC_INVARIANTS.md:405`
  already records *"both sectors using `mass_tot_dc` for the metric factor."* The metric is sourced
  by all the mass inside `r`, not by one sector's. Any migration must keep this, and a detector
  must exist for the substitution (§18).
- **Coverage is zero and reachability is nil.** `CURRENT_ARCHITECTURE.md:91` — COMPILED,
  UNEXERCISED. Re-authenticated: the only `main/` reference is `tov_debug_main.cpp:51`, a
  `AlwaysExportMixed` predicate whose registration at `:83` is **commented out**. There is no
  live path and no test.
- **Six sites, in two duplicated build paths** — `SurfaceIsReached` (`:238-264`) and the
  constructor (`:659-686`) — plus two scalar accessors.

**Recommended 3D scope: govern `MixedStar` conceptually here; defer the source migration.**
See Q2 (§21). The ADR README already lists *"MixedStar modernization scope"* as an anticipated
Phase-3 ADR; that is the correct place for the migration, once coverage exists.

**`MixedStar` is legacy but live in the build.** Nothing in this ADR authorizes deleting or
modifying it.

---

## 16. Candidate scientific code — explicit non-scope

Sites #16–#21 (`DarkCore_Analysis`, `BNV_Analysis`, `Decay_Analysis`, `BNV_Chi`,
`BNV_B_Chi_Photon`) and #22 (`RotochemicalCache`) each reproduce their own volume factor. They are
**COMPILED, UNEXERCISED** or **NOT COMPILED**, unratified, and protected by `GOVERNANCE.md` §5 and
`CURRENT_ARCHITECTURE.md:110,132`.

**Do not migrate candidate scientific code merely to achieve textual deduplication.** No test
covers them; a refactor there would be unobserved, and touching unratified physics under a
structural ADR is exactly what §5 exists to prevent.

This ADR establishes only the **future contract** they must satisfy *if and when* they are
activated: obtain `w_V` from the single owner, apply their own `e^{ν}`/`e^{2ν}` and their own unit
conversion, and never re-derive the measure. Site #23 (`TOVSolver.cpp:1976-1980`) is
preprocessor-disabled (`TOVSolver.cpp:33`) and is out of scope for the same reason.

---

## 17. Alternatives considered

### Option A — `GeometryCache` is both the mathematical and the cached owner
Core consumers construct or receive a `GeometryCache`. **This is the Phase-3-entry wording.**
Blocked by §8.2 obstructions 1–3: unreachable for `MixedStar`, ADR-0003-conflicting for `NStar`,
and it couples the measure to a `ν` requirement it does not have. Not recommended — and *not*
rejected merely because it is inconvenient, but because it cannot be implemented for the sites
INV-04 names.

### Option B — dependency-neutral primitive owns the mathematics; `GeometryCache` owns the cached arrays  ← **RECOMMENDED**
§10. One formula; no new dependency edge; reaches `NStar`, `MixedStar` and the scalar accessors;
leaves `GeometryCache` canonical and bit-identical.
*Cost:* one new low-level abstraction, and the clamp/centre semantics must be defined explicitly
rather than inherited (§11 — which this ADR treats as a benefit, since they are currently
inconsistent six ways).

### Option C — `StarProfile` owns the metric / proper-volume calculation
`StarProfile` is a storage and versioning type (`CURRENT_ARCHITECTURE.md:88`); it already *stores* a
`MetricLambda` column, so this is superficially attractive. Rejected: `MixedStar` has no
`StarProfile` at all, so the largest cluster of duplicates remains unreachable; and it would move
derived-quantity computation into the type whose whole ADR-0003 role is to be the *source* whose
revisions derived caches key on.

### Option D — `Core` owns a geometry helper; `Evolution` depends on it
Workable for `NStar`, and `MixedStar` is also in Core. But it forces
`Physics/Evolution/GeometryCache.hpp` to include a `Core` header, destroying the deliberate
property that Evolution's *headers* are Core-free (§8.1 fact 2) — a strictly worse layering than
Option B for no compensating benefit. It also leaves `EOS` (which has no `Physics` edge at all)
unable to use the primitive without a new `EOS → Core` header edge.

### Option E — declare `GeometryCache` canonical by documentation only; leave the inline forms
Zero behavior risk; INV-04's fail-closed condition is *"no document says which is authoritative"*,
which a document can technically discharge. Rejected as insufficient: it leaves nineteen
independent copies of the formula, six inconsistent degenerate-input behaviors (§11), and a
`4πr²` in `DarkCore_Analysis` that a previous audit already misclassified (§4.3) precisely because
the definition is scattered. It converts a structural problem into a documentation obligation
that will decay.

### Decision table

| | **A** GeometryCache owns both | **B** primitive + GeometryCache | **C** StarProfile owns | **D** Core helper | **E** docs only |
|---|---|---|---|---|---|
| **Correctness** (one formula) | yes, where reachable | **yes, everywhere** | yes for `NStar` only | yes | no — 19 copies remain |
| **Dependency direction** | **Core → Evolution (reversed)** | **no new edge** | none new | Evolution → Core *header* (regression) | none |
| **Works for `NStar`** | conflicts with ADR-0003 (§8.2-2) | **yes** | yes | yes | n/a |
| **Works for `MixedStar`** | **no** (no `StarProfile`) | **yes** (scalar `(r,m)`) | **no** | yes | n/a |
| **Works for `GeometryCache`** | trivially | **yes, bit-identical** | yes | yes | yes |
| **Behavior-preservation risk** | high — forced ctor in an `EditScope` | **low** — composition kept verbatim | medium | low | zero |
| **API surface** | needs a new scalar interpolating API | 3 free functions | grows a storage type | 3 free functions in Core | none |
| **Future Hartle / rotochemical use** | needs `StarProfile` + `ν` | **usable from any layer** | limited | usable, with a Core edge | none |

---

## 18. Validation and detector requirements for a future implementation

An accepted implementation must satisfy all of the following. **These tests do not exist yet and
this ADR does not create them.**

### 18.1 Geometry-primitive identity (new, self-contained)
Analytic checks over controlled compactness `2M/R ∈ [0.05, 0.98]` on the uniform-density interior:
`e^{Λ}` vs `1/√f` vs `exp(−½ ln f)`; `w_V` vs `4πr² e^{Λ}`; exact Schwarzschild-exterior values.
Predeclared agreement: **≤ 2 ULP** (§6.2 measured exactly that, so the bound is evidenced, not
invented).

### 18.2 `GeometryCache` conformance (new, self-contained)
`WV`, `WVExpNu`, `WVExp2Nu`, `ExpLambda` must equal the primitive **bitwise** on a synthetic
profile, and must be **bitwise unchanged** from the pre-change build on the authenticated
1.6 M☉ star. Both the stored-Λ and derived-Λ routes must be exercised (§9 shows they agree today;
the test must keep it that way).

### 18.3 `NStar` baryon count
`SeqPoint::b` on the 1.0 / 1.4 / 1.6 / 2.0 M☉ stars must stay within the **predeclared
`1.0e-15`** relative bound of the values recorded in §7.2. Note that **no golden artifact
currently asserts `B`** — capturing those four values is part of the implementation, not optional.

### 18.4 Thermal and structural consumers — must not move at all
`heat_capacity_v1`, `heat_capacity_real_star`, `passive_cooling_regression`, the neutrino
diagnostics, `grid_convergence_cmf`: **byte-identical artifacts**, since §13 requires bit-identity
for every `w_V` consumer.

### 18.5 Hartle
`hartle_moment_inertia_analytic` and `hartle_moment_inertia_cmf` must be **unchanged**.
`hartle_I_dscmf1_debug.tsv` must stay byte-identical. Hartle does **not** consume this measure
(§4.4); if it moves, the change reached something it should not have.

### 18.6 Full suite
13/13 with the data root, 8/8 self-contained, plus the new 3D contracts. No new test may be added
that merely re-asserts an existing one.

### 18.7 Detector proof plan — the tests must *fail* under each of these mutations

| # | Injected fault | Expected detector |
|---|---|---|
| **D1** | wrong metric denominator (`1 − m/r`, or `1 + 2m/r`) | §18.1 analytic identity; §18.3 `B` moves ≫ 1e-15 |
| **D2** | `e^{Λ}` omitted from `w_V` (coordinate volume substituted) | §18.2 bitwise conformance; `heat_capacity_v1` (`C_⋆` shifts by ~10 %) |
| **D3** | `1e54` moved into the primitive, or added to the zero-caller accessor | §18.3 `B` off by 10⁵⁴; a static check that the primitive's translation unit contains no `1e54` |
| **D4** | component mass substituted for total mass in `MixedStar` | **requires new `MixedStar` coverage — does not exist today.** This is the concrete reason §15 says coverage must precede migration |
| **D5** | degenerate-input contract violated (clamp changed, or a `NaN`/`inf` escapes) | a domain test at `f → 0⁺`, `f = 0`, `f < 0`, `r = 0`, non-finite `r`/`m`, asserting the §11 contract the owner selects |
| **D6** | redshift variant composed from a second measure | `WVExpNu`/`WVExp2Nu` bitwise equality with `WV × ExpNu` / `WV × Exp2Nu` |

**No source mutation was executed by this audit.** D1–D6 are a plan.

---

## 19. Proposed migration sequence

Each step is independently revertible and independently evidenced.

1. **3D-1** — add the primitive; **no caller changes**. Prove it bitwise-reproduces
   `GeometryCache::DeriveLambdaFromMR_` and the analytic identities (§18.1). Expected impact:
   **none**, nothing calls it. `BIT-IDENTICAL` trivially.
2. **3D-2** — `GeometryCache::DeriveLambdaFromMR_` delegates to the primitive; composition kept
   verbatim. `BIT-IDENTICAL REQUIRED`; all five golden artifacts unchanged.
3. **3D-3** — `NStar.cpp:277-281` (`BuildFromTOV`) consumes the primitive. **Predeclared
   `1.0e-15`** on `seq.b`; the four §7.2 values become the new recorded reference.
4. **3D-4** — `NStar.cpp:562-569` (`FinalizeSurface`) — **only after** TOV Path-1 coverage exists.
   Same tolerance. May legitimately be deferred to 3E, which already owes Path-1 an equivalence
   harness (`PHASE3_CONSOLIDATION_PLAN.md` §11).
5. **Not scheduled by this ADR:** the scalar accessor (§14), `MixedStar` (§15), the candidates
   (§16).

---

## 20. Consequences

**If accepted:**
- INV-04's *ownership* question is decided; INV-04 moves to `GOVERNED (ACCEPTED)` only when the
  implementation lands and validates, exactly as INV-12 did under ADR-0003.
- `GOVERNANCE.md` §3 condition 5 is discharged for the proper-volume measure.
- One formula, one clamp, one degenerate-input rule, reachable from every layer.
- `GeometryCache` is confirmed canonical **for the cached representation** and its consumers do
  not move at all.
- The `DarkCore_Analysis` misclassification is corrected in the record.

**Costs:**
- One new low-level header at the `CompactStar/` root.
- `seq.b` may move by up to 1 ULP on some stars (§7.2). That is a real, predeclared numeric change
  on a quantity with **no golden baseline today** — which the implementation must create.
- The degenerate-input contract must be chosen (Q3), and whichever is chosen, at least four of the
  six current behaviors stop being reproduced. No validated output moves (§11), but the change is
  semantic and must be recorded as such.

**If rejected or deferred:** INV-04 stays fail-closed; the nineteen copies and six clamp behaviors
remain; and Phase 3E (canonical TOV path), which `PHASE3_CONSOLIDATION_PLAN.md` §11 lists as
depending on 3D, stays blocked.

---

## 21. Explicit non-scope

This ADR does **not**:

- authorize any implementation;
- change INV-03, INV-04, INV-12, INV-14 or INV-15, or any ADR-0001/0002/0003 provision;
- change `NStar::BaryonNumIntegrand`'s missing `1e54` (§14) or `MixedStar`'s equivalent;
- migrate, modify or delete `MixedStar` (§15);
- migrate, modify or delete candidate scientific code (§16);
- touch the Hartle ODE coefficients, the TOV mass ODE, `EvaluateNu`'s surface BC or
  `SurfaceGravity` (§4.4);
- reopen the solar-mass authority question (`PHASE3_CONSOLIDATION_PLAN.md` §4, U-4), which changes
  `m [km]` and hence every metric factor, and is deferred out of Phase 3 entirely;
- reopen ADR-0002 §6 (Pattern A vs B);
- select a numerical tolerance after observing migrated output — §7.1 is derived in advance and is
  binding.

---

## 22. Owner adjudication questions

> **RESOLVED 2026-09-01.** Q1 = Option B; Q2 = governed now, migration deferred; Q3 = hybrid
> physical-domain contract. See §0. The questions are retained below as posed, because the
> reasoning that framed them is part of the record.

### Q1 — What does "single owner" mean?

Should INV-04's single owner be

- **(a)** one **dependency-neutral mathematical primitive** owning `f`, `e^{Λ}` and `w_V`, with
  **`GeometryCache` retained as the canonical cached representation** — Option B; or
- **(b)** **`GeometryCache` itself**, requiring `Core::NStar` and `Core::MixedStar` to depend on
  `Physics::Evolution` — Option A, the Phase-3-entry wording?

**Recommendation: (a).** Not on taste: (b) cannot reach `MixedStar` at all
(`MixedStar.hpp:269,272` — no `StarProfile`), would construct a provenance-bearing
`GeometryCache` inside an open `EditScope` in conflict with **ADR-0003**
(`NStar.cpp:90,277-281`; `StarProfile.hpp:162-166`; `StarContext.cpp:763-771`), and would couple
the measure to a `ν` column it does not use (`GeometryCache.cpp:178-179`). (a) satisfies INV-04
**more** completely, adds no dependency edge — the layering is already settled by the `Units.hpp`
precedent from Phase 3A (`Units.hpp:35-39`) — and keeps every current `w_V` consumer bit-identical.

### Q2 — `MixedStar` scope for 3D

Should the 3D implementation **(a)** migrate `MixedStar`'s six sites now, after first adding
focused coverage; or **(b)** have this ADR govern `MixedStar` conceptually while the source
migration waits for the already-anticipated MixedStar modernization ADR?

**Recommendation: (b), on risk.** `MixedStar` is **COMPILED, UNEXERCISED** — its only `main/`
reference is a predicate whose registration is commented out (`tov_debug_main.cpp:51,83`) — has
**zero** test coverage, carries six duplicated sites across two build paths, and has two-sector
mass semantics (`MixedStar.cpp:169-231`) that no existing detector protects. Detector **D4**
(component mass substituted for total mass) has **nothing to fire in today**. Migrating unobserved
code inside a structural increment is precisely the risk `AGENTS.md` §5 forbids. (a) is acceptable
only if the coverage is built and merged *first*, as its own increment — at which point it is (b)
with extra steps.

### Q3 — Degenerate-input semantics

The `1e-15` clamp of INV-03 is **not** uniform. §11 measures **six** mutually inconsistent
behaviors for `f ≤ 0`, differing by a factor of 3.16e7 at `f = 0` and including a **silent
divisor substitution inherited from `Zaki::Vector::DataColumn::operator/=`** that nineteen call
sites never asked for. A single owner must pick one, and at least four of the six will stop being
reproduced.

Which contract governs `r ≤ 0`, `f = 0`, `f < 0`, and non-finite `r`/`m`?

- **(a)** clamp `f` to `ε = 1e-15` — status quo B1, preserves `GeometryCache` exactly;
- **(b)** return a finite regularized value with an explicit, documented `ε`;
- **(c)** return zero — B4, the scalar accessor's current behavior;
- **(d)** throw / fail closed, per `GOVERNANCE.md` §3.

**No recommendation is offered, and none should be inferred.** The evidence establishes that the
implementations disagree; it does not establish which is right, and the choice is a scientific
convention, not a numerical detail. **No currently validated output moves under any option** —
measured `max 2m/r = 0.481` across the four authenticated stars, so no node reaches the degenerate
branch. What the choice determines is what happens when a future EOS, a mixed-star configuration,
or a malformed input first does.

*(A supplementary note for whichever option is chosen: option (a) or (b) leaves `GeometryCache`
bit-identical and is the only pair compatible with §13's `BIT-IDENTICAL REQUIRED` row without
further evidence. (c) and (d) would make the primitive's adoption in `GeometryCache` a
scientific-semantic change in its own right.)*

---

## 23. Implementation record (Phase 3D-I)

Appended after implementation. **The accepted semantic decisions of §0 were not altered by
implementation**; nothing in the evidence below contradicts them.

**Commits.** Acceptance `f92df9a` (`docs: accept proper-volume ownership contract`);
implementation `refactor: establish canonical proper-volume measure`. Acceptance precedes
implementation in history, as required.

**Primitive.** `CompactStar/Geometry.hpp`, namespace `CompactStar::Geometry` — four free scalar
functions: `MetricDenominator`, `Lambda`, `ExpLambda`, `ProperVolumeWeight`. Standard library
only (`<cmath>`, `<stdexcept>`, `<string>`); no state, no registry, no hierarchy, no templates;
no dependency edge. `Core` was **not** made to depend on `Physics/Evolution`.

**Hybrid domain contract (§0-Q3) as built.** Regular centre returns `f = 1`, `Λ = 0` (as `-0.0`,
which is bit-identical to what the legacy `r ≤ 0 → denom = 1` branch produced and composes
correctly), `e^Λ = 1`, `w_V = 0`. `r < 0`, `r = 0` with `m ≠ 0`, `f ≤ 0`, and non-finite input
all throw `std::runtime_error` — the project's convention. No clamp, no epsilon, no tolerance
band around `r = 0`.

**Migrated — the validated path only.** `NStar::BuildFromTOV` Λ production (bit-identical) and
its baryon-number integrand; `GeometryCache::DeriveLambdaFromMR_` (delegates; the `eps`
parameter is removed). `GeometryCache`'s `DataColumn` composition is kept verbatim, so its
cached arrays are bit-identical by construction rather than by measurement — confirmed bitwise
96/96 on both Λ routes.

**Baryon-number result.** The §7.2 prediction, made from a scratch replication **before any
implementation existed**, was reproduced exactly: `|ΔB|/B = 1.368e-16` (one ULP) on the 1.0 M☉
star, **bitwise** on 1.4 / 1.6 / 2.0, with `M`, `R` and `ec` bitwise on all four. The pre-mutation
capture reproduced all four §7.2 values **bitwise**, authenticating that evidence on the
implementing machine. The `1.0e-15` bound of §7.1 was **not widened**.

**Behavior preservation.** All five protected artifacts **byte-identical**. The conditional
B-only exception the implementing brief permitted for `grid_convergence_cmf_1p6_debug.tsv` was
**not needed** — that artifact carries a `B` column and did not move at any resolution. No
golden was re-baselined.

**Tests and detectors.** New: `proper_volume_contract`, `geometry_cache_measure_contract`
(self-contained) and `baryon_number_cmf` (external-data) with the durable canonical reference
`tests/baselines/baryon_number_dscmf1_reference.tsv`. Detectors **D1, D2, D2b, D3, D5, D6 all
fire**; every mutation was reverted byte-identically by SHA-256. **D4 was correctly not
executed** — `MixedStar` is unmigrated and uncovered, so per §18.7 there was nothing for it to
fire in. Three initial non-firings were diagnosed rather than accepted: `heat_capacity_v1`'s
fixture is flat space (`m = 0` at every node, so `e^Λ ≡ 1`); `heat_capacity_real_star` contains
**no assertions at all** and is a diagnostic harness rather than a detector; and the first D6
mutation was bit-identical to the composition it replaced. The material detector for the cached
measure is `passive_cooling_regression`, under which dropping `e^Λ` moves `C_⋆` by **15–17 %**.

**Deferred, and still governed by this ADR (§13, §15, §16).** TOV Path 1 (`NStar::Append`,
`NStar::FinalizeSurface`); the zero-caller scalar `NStar::BaryonNumIntegrand(double)` with its
separate INV-14 defect, deliberately **not repaired**; all six `MixedStar` sites; and the §5
candidate code. None was touched.

**Invariants.** INV-03 keeps its headline; its clamp wording is superseded by this ADR and its
entry now separates governed contract, canonical conformance and deferred legacy nonconformance.
**INV-04 is NOT resolved** — it is `GOVERNED (ACCEPTED) — CANONICAL VALIDATED PATH CONFORMED;
LEGACY MIGRATIONS DEFERRED`. INV-14 untouched.

**Validation.** 16/16 authenticated, 10/10 self-contained. No test deleted, no tolerance widened.

Full record: [`docs/validation/PHASE3D_PROPER_VOLUME.md`](../validation/PHASE3D_PROPER_VOLUME.md).

---

## 24. Implementation record — Phase 3E-I2 (Path-1 conformance)

§13 deferred TOV Path 1 *"pending new coverage first."* Phase 3E-0 supplied that coverage
(`docs/validation/TOV_PATH_EQUIVALENCE.md`), so the deferral ended. **The accepted decisions of
§0 are unchanged.**

**Scope clarification.** ADR-0005's migration shorthand named only the `FinalizeSurface`
proper-volume factor, but the authenticated source held a **second** Path-1 nonconformance:
`NStar::Append` independently computed `f = 1 − 2m/r`, clamped `f ≤ 0 → 1e-15`, and formed
`Λ = −½ ln f`. That is governed by this ADR and was deferred for the same reason. Phase 3E-I2
migrated **both** sites; it did not broaden further.

| Site | Before | After |
|---|---|---|
| `NStar::Append` — Λ | local `denom`/clamp/`-0.5*log` block | `CompactStar::Geometry::Lambda(r_km, m_km)` |
| `NStar::FinalizeSurface` — baryon integrand | `r².pow(2) * 4π n_B /(1−2m/r).sqrt() * 1e54` | `w_V = Geometry::ProperVolumeWeight(r,m)` per node, then `w_V · n_B · 1e54` |

`FM3_TO_KM3` was **not** moved into `Geometry` — it belongs to INV-14, not to `dV`.

**Λ: bit-identical.** Phase 3E-0 measured Path-1 Λ ≡ Path-2 Λ bitwise *before* this migration;
after it, Λ is still bitwise between the paths (`max_profile_ulp = 0` across 14 central densities
and `radial_res` 5000/10000/20000). Path 2 is provably untouched — `BuildFromTOV` calls neither
migrated function — so **Path-1 Λ is bit-identical pre/post**.

**B: the governed movement, and a stronger-than-required outcome.** The composition now mirrors
`BuildFromTOV`, so both ordinary paths share one mathematical owner *and* one arithmetic
ordering. Result: **Path-1 `B` is now bitwise identical to Path-2 `B` at all 17 measured
comparisons.** The five comparisons that previously carried a one-ULP gap (`1.387e-16`,
`1.637e-16`, `1.269e-16`, `1.640e-16`, `1.640e-16`) are now exactly `0`. Path-1 `B` therefore
moved by at most **1.640e-16**, well inside the `1.0e-15` predeclared in §7.1 before any
implementation existed. The bound was not widened; the gap is closed rather than merely bounded.

At the four M☉ anchors Path-1 `B` was already bitwise-equal to Path 2 before I2, so those values
did not move at all: `1.27388873109076839e+57`, `1.83218336257875150e+57`,
`2.12457569547972117e+57`, `2.74576306114479063e+57`.

**`baryon_number_dscmf1_reference.tsv` is unchanged and was not re-baselined** — it is produced
through Path 2, which I2 does not touch. `baryon_number_cmf` reports worst `|ΔB|/B = 0.000e+00`.

**Detector D4.** Removing the proper-volume factor from the migrated `FinalizeSurface` fires
**8 assertions, every one of them B-related**, at `|ΔB|/B` = 8.4 %–19.5 %, while **all 7
radial-column assertions stay green**. That selectivity is what makes the
RADIAL-INTEGRATION vs BARYON-INTEGRATION classification trustworthy. Reverted byte-identically.

**Still nonconformant, and deliberately untouched:** `NStar::BaryonNumIntegrand(double)` (INV-14
defect, zero callers), all six `MixedStar` sites (no coverage; §0-Q2), the `GOVERNANCE.md` §5
candidate code, and `NStar::EvaluateNu`'s boundary-condition `x = 1e-15` (not this measure,
§4.4).

Full record: [`docs/validation/PHASE3E_I2_PATH1_GEOMETRY.md`](../validation/PHASE3E_I2_PATH1_GEOMETRY.md).

**Phase 3F boundary note (owner clarification, no change to §0 decisions).** The "candidate
scientific code" this ADR lists in §16 spans two different categories: the `GOVERNANCE.md` §5
core candidates from `675b4a9` (rotochemical), and the **project-specific extension modules**
`DarkCore_Analysis`, `BNV_*`, `Decay_Analysis`, which consume the core rather than form it. The
§16 *contract* (obtain `w_V` from the single owner, never a third measure) applies to both when
they are next touched; the extension modules are **not** a prerequisite for core closure. See
`docs/validation/PHASE3_CLOSEOUT.md` §7a.
