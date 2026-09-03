# Phase 4C-I0 — the EOS thermodynamic-derivative authority `dε/dp`

> **FORMAL STATUS: `PHASE-4C-I0 EOS DERIVATIVE AUTHORITY COMPLETE`**
>
> ADR-0007 was **ACCEPTED** by the project owner on 2026-09-02 (commit `13111de`), and this
> increment implements exactly one clause of it — **P5**, the `dε/dp` authority — and nothing
> else. **No O(Ω²) physics was added, activated, repaired or validated.** The four second-order
> candidate functions and `HartleResult` are **byte-identical** to `master` (`df859b5`);
> `RotationSolver.hpp` and `RotationSolver.cpp` were not touched at all. No monopole baseline
> exists and none may be created before 4D. `GOVERNANCE.md` §3.1 is **AUTHORIZED** by ADR-0007
> §9 and **NOT YET EXERCISED** — this increment changes no scientific output.

| Field | Value |
|---|---|
| **Starting HEAD** | `bcef5b57e3267fe26d1269efaa472d314818d367` — upstream equal, tree clean, **5 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Acceptance commit** | `13111de010375fe2ccbb371784ef2d27446a90ab` — `docs: accept Hartle second-order monopole contract` (documentation only; acceptance precedes implementation) |
| **Change class** | **scientific-semantic** (a new physical quantity acquires a single authority and a unit contract) **and structural** (EOS/TOV layer gains ownership; `StarProfile` gains a profile-attached carrier) — strictest applies, `GOVERNANCE.md` §2 |
| **Governing authority** | **ADR-0007 P5** (ACCEPTED 2026-09-02) — *"`dε/dp` is the derivative of the SAME `ε(p)` interpolant that constructed the star, owned by the EOS/TOV layer."* Also INV-02 (cgs inside TOV, geometric outside), INV-13 (interpolation method), ADR-0003 (`Version()` provenance), ADR-0001 (species are fractions), ADR-0005 (`I` and the goldens bitwise) |
| **Evidence** | `tests/eos/eos_derivative_contract.cpp` (self-contained, 15 checks), `tests/eos/eos_derivative_cmf.cpp` (external data, 17 checks), detectors D1–D3 (§11) |
| **Baseline** | pre: **23/23** (209.58 s) + **12/12** (14.10 s), seven artifact hashes as `PHASE3_CLOSEOUT.md` §1; post: **25/25** + **13/13**, same seven hashes (§12) |

---

## 1. Starting state

```
worktree  /Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotation
branch    physics/rotation-correctness      HEAD = bcef5b5   upstream equal   clean
master    df859b5                            5 ahead / 0 behind
```

Pre-task suites on that exact tree: **23/23 PASS** (209.58 s) and **12/12 PASS** (14.10 s),
confirmed from CTest's own `LastTest.log` (23 and 12 `Test Passed.`, zero `Test Failed.`, no
`LastTestsFailed.log`). Seven durable artifacts, `shasum -a 256`, unchanged from Phase 3:

| Artifact | SHA-256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |
| `tov_path_equivalence_dscmf1.tsv` | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |

## 2. ADR acceptance (precedes implementation)

Commit `13111de` moved ADR-0007 from PROPOSED to **ACCEPTED**, recording the owner's eight
adjudications against the ADR's thirteen decision questions, plus one modification:

| Owner decision | Effect here |
|---|---|
| canonical variable = Hartle `p₀*` | not exercised in I0 (4C-I1) |
| governed family = fixed `ε_c`, no surface condition | not exercised in I0 |
| products = seed-free coefficients per `Ω_geom²`, materialized only at an explicit `AngularVelocity` | not exercised in I0 |
| **`dε/dp` = derivative of the SAME `ε(p)` interpolant, EOS/TOV-owned** | **this increment** |
| surface `SURFACE ADEQUATE AS-IS`; `R_*` never labelled the exact `p = 0` surface | not exercised in I0 |
| validated `l = 0` is the Phase-4 deliverable unlocking Phase 5; `l = 2` a separate future extension, **not** thereby validated | recorded in the roadmap and invariants |
| candidate API replaced atomically in **4C-I1**, not deleted in I0 | candidate left byte-identical (§12) |
| homogeneous sequence-derivative response **not** public in 4C (validation use only, 4D) | nothing exposed here |
| `GOVERNANCE.md` §3.1 record accepted | **AUTHORIZED, NOT YET EXERCISED** (§13) |

Acceptance fixes a contract; it certifies no number. INV-08 is now
`GOVERNED (ADR-0007 ACCEPTED) — REPLACEMENT CONTRACT ESTABLISHED; CURRENT CANDIDATE
NONCONFORMANT / REPLACEMENT PENDING`.

## 3. What the derivative is

Hartle's `l = 0` mass equation (H67 eq. 97; HT68 eq. 15a) carries the factor `dE/dP`: the
derivative of the **one-parameter equation of state** `E = E(P)` (HT68 eq. 1) that defines the
star — the barotropic

```
dε/dp = 1/c_s²        [dimensionless]
```

with `c_s` the adiabatic sound speed of the cold EOS. It is **not** `dρ/dp` (rest-mass density;
related by `dN/N = dE/(E+P)`, HT68 eq. 27), and **not** a derivative of the radial profile.
Physically `dε/dp > 0` for matter with a finite sound speed; causal matter has `dε/dp ≥ 1`, and
incompressible matter has `dε/dp = 0` as a formal limit — so **no causal bound is asserted in
either direction**, because both `0` (the exact analytic value for a constant-density star) and
values below `1` (admissible for a synthetic test EOS) are legitimate inputs. What *is* asserted
on real matter is positivity and finiteness (§9).

## 4. Unit derivation — from the code's own constants, not from prose

The governed input tabulates `ε` as a **mass** density in g/cm³ against `p` as an **energy**
density in dyne/cm² = erg/cm³ (authenticated from the file header of
`DS(CMF)-1_with_crust.eos`: `e(g/cm^3)`, `p(dyne/cm^2)`, `rho(1/fm^3)`, then 17
species-fraction columns). The raw spline derivative therefore carries

```
[g cm⁻³] / [dyn cm⁻²] = [g cm⁻³][cm² s² g⁻¹ cm⁻¹] = s² cm⁻²
```

The profile conversions in `NStar::BuildFromTOV` (`NStar.cpp:180-186`) are

```
ε_geom = ε_cgs · (INV_FM4_2_INV_KM2 / INV_FM4_2_G_CM3)
p_geom = p_cgs · (INV_FM4_2_INV_KM2 / INV_FM4_2_Dyn_CM2)
```

so, dividing, the geometric (hence dimensionless) derivative is

```
dε_geom/dp_geom = (INV_FM4_2_Dyn_CM2 / INV_FM4_2_G_CM3) · dε_cgs/dp_cgs
```

**derived from the same two constants the profile itself uses**, not assumed. That ratio *is*
`c²` in cgs, because `INV_FM4_2_Dyn_CM2` takes an energy density to erg/cm³ while
`INV_FM4_2_G_CM3` takes the same energy density to its mass-density equivalent. The identity is
**asserted, not assumed**: check `Ua` compares it against an independent literal
`c = 2.99792458e10 cm/s` (exact by definition) and measures

```
INV_FM4_2_Dyn_CM2 / INV_FM4_2_G_CM3 = 8.987551787368e20        c² = 8.987551787368e20
relative difference = 1.458e-16                                bound 1e-13
```

The factor is applied in exactly one place, `EosDerivCgsToDimensionless()`
(`TOVSolver.cpp`, file-static), which is the sole owner of the conversion on this path.

## 5. API ownership — the narrowest architecture that satisfies A–E

The task's five requirements and how the chosen design meets each:

| | Requirement | How |
|---|---|---|
| A | EOS/TOV layer owns evaluating `dε/dp` | `TOVSolver::GetEDensDeriv` (`TOVSolver.hpp:754`) plus `HasEDensDeriv`, `EDensDerivPressMin/Max`. No other class evaluates it |
| B | `RotationSolver` never owns an EOS object or differentiates a profile | `RotationSolver.{hpp,cpp}` are **untouched by this increment** (`git diff` empty). The future solver reads the profile |
| C | a completed `StarProfile` can supply authoritative `dε/dp` at its nodes | `StarProfile::HasEosDEdP()` / `GetEosDEdP()` |
| D | point-constructed analytic stars can explicitly supply their own | `TOVPoint`'s trailing `dedp` parameter (§10) |
| E | absence fails closed when a consumer requiring it is invoked | `HasEosDEdP()` is false and `GetEosDEdP()` returns `nullptr`; the future solver checks. Nothing that does **not** need it is affected |

**The flow.** `Hidden_ImportEOS_Vis` builds `visi_eps_p_spline` → `SingleStarSolveToTOVPoints`
evaluates `GetEDensDeriv(y[1])` on the same iteration that evaluates `GetEDens(y[1])`, and
stores it on the `TOVPoint` → `NStar::BuildFromTOV` (bulk) or `NStar::Append` +
`FinalizeSurface` (incremental) installs it on the profile. Both ordinary-`NStar` construction
paths are supported, deliberately: a scientific authority that depended on which constructor was
used would be exactly the "uncertain authoritative code path" `GOVERNANCE.md` §3 condition 3
fails closed on.

**Why not a `radial` column** (`StarProfile`'s other option). The derivative is a property of the
*equation of state*, not an output of the radial integration, and keeping it out of `radial`
means the profile's column count, its species-column indices (still `8 … 8+n_species-1`), its
`Export` layout and every file derived from them are untouched — which is what lets this
increment guarantee byte-identical scientific output (§12) rather than merely hope for it. Two
registered tests read `Radial().Dim().size()` (`tov_path_equivalence_cmf.cpp:185`,
`grid_convergence_cmf.cpp:539`); a new column would have walked into both. There is still exactly
**one** authoritative copy per profile, it is cleared by `Reset()`, it travels with the profile
under the compiler-generated copy/move, and every mutator calls `Touch()` so any cache keyed on
`Version()` sees it (ADR-0003).

**Acceptance is all-or-nothing.** `SetEosDEdP` validates size against the radius column and
finiteness of every entry, and refuses anything else; `AppendEosDEdP` only stages, and
`FinalizeEosDEdP` applies the same test before publishing. A partial or corrupt set can never be
published as authoritative, so a consumer may trust `HasEosDEdP()` without inspecting nodes.

## 6. Interpolation authority — authenticated, not assumed

`TOV_gsl_interp_type = gsl_interp_steffen` (`TOVSolver.hpp:488`) is the **single** interpolation
type used for every EOS spline: `visi_eps_p_spline`, `visi_rho_p_spline`,
`visi_rho_i_p_spline[j]` and their dark counterparts, allocated at `TOVSolver.cpp:684-712`.
`GetEDensDeriv` differentiates **that same object** through **the same accelerator**
(`visi_p_accel`) that `GetEDens` uses. No second interpolant is constructed anywhere.

Steffen is a monotonicity-preserving cubic Hermite: `C¹`, exact on collinear data, `O(h²)` on
smooth data, with one-sided slopes at the endpoints and a limiter that engages where monotonicity
would otherwise be violated. Its second derivative jumps at knots, which is why the finite-difference
cross-check of §8 is taken away from knots. *(Aside, not repaired here: the comment block at
`TOVSolver.hpp:494-498` still describes a natural cubic spline, and the `PressureCutoff` comment
at `:549-554` still says `1e-5` where the code uses `1e-15` — both stale comments, no behaviour
difference, and outside this increment's scope.)*

Checks `Se` (§8) and detector `D3` (§11) exist precisely to hold this clause: they fail if the
returned value comes from any interpolant other than the star's own.

## 7. Derivative-domain semantics — stricter than a value lookup, deliberately

`GetEDens` goes through `SafeSplineEval`, which **clamps** an out-of-range pressure to the table
edge and logs a warning. That is defensible for a value. It is **not** defensible for a
derivative: a derivative is a *local* property of the state, so a clamped answer would describe
a different state than the one asked about, silently.

`SafeSplineEvalDeriv` (`TOVSolver.cpp`, file-static) therefore **fails closed** and returns NaN,
with a logged error, in every case it cannot answer exactly:

| Condition | Result |
|---|---|
| no interpolant (before `ImportEOS`, or a table too short/non-monotonic) | NaN |
| non-finite argument (NaN, ±inf) | NaN |
| `p < xmin` or `p > xmax` — **no clamping** | NaN |
| GSL error, or a non-finite derivative | NaN |
| `p == xmin` or `p == xmax` exactly | **evaluated** (endpoints are inside the domain) |

NaN rather than `0.0` is the deliberate choice: `0.0` is the physically meaningful value for
incompressible matter, so returning it on failure would make an error indistinguishable from a
legitimate answer. Callers test `std::isfinite`. This is asserted by checks `Da`, `Db`, `Dc`,
`Dd`, `De` (§8).

## 8. Analytic contract — `tests/eos/eos_derivative_contract.cpp` (self-contained)

The oracle is a closed-form `ε(p)` that the test both tabulates and differentiates itself, so
the comparison never asks CompactStar for the answer. Bounds were **predeclared from the
interpolation order before any comparison was run** and are recorded verbatim in the file header.

| Bound | Value | Basis |
|---|---|---|
| `L` linear reproduction | rel ≤ `1e-12` | a cubic Hermite with the common secant slope reproduces collinear data exactly; only floating point separates them |
| `S1` smooth interior, `N = 200` | rel ≤ `1e-3` | `O(h²)` with `h/p = ln(100)/200 = 2.3e-2` |
| `S2` refinement `N → 2N` | factor ≥ `3` | `O(h²)` predicts 4; 3 leaves room for the limiter |
| `S3` endpoints (one-sided) | rel ≤ `1e-2` | strictly weaker than the interior, reported separately rather than hidden |
| `F` vs central FD of the **same** spline, away from knots | rel ≤ `1e-6` | FD step `1e-5·p`, ~`1e-3` of the local knot spacing, so the stencil never straddles a knot |
| `U` unit identity | rel ≤ `1e-13` | exact-by-definition `c` vs the repository's constants |

**Measured — 15/15 checks pass.**

| Check | Result |
|---|---|
| `Ua` `INV_FM4_2_Dyn_CM2/INV_FM4_2_G_CM3 == c²` (independent literal) | `1.458e-16` |
| `Da` no EOS ⇒ `HasEDensDeriv()` false and NaN returned | pass |
| `La` interpolant exists and reports its own domain | `[1.0e33, 1.0e35]` |
| `Lb` 239 probes — every knot, midpoint and both endpoints — return the exact analytic value | worst rel `5.534e-14` (bound `1e-12`) |
| `Db` below domain ⇒ NaN, **not** the boundary derivative | pass |
| `Dc` above domain ⇒ NaN | pass |
| `Dd` exact endpoints evaluate | pass |
| `De` NaN / +inf / −inf ⇒ NaN | pass |
| `Sa` both resolutions evaluate everywhere, no non-finite value | pass |
| `Sb` `N=200` interior | max rel `8.568e-05` (bound `1e-3`) |
| `Sc` refinement `N=200 → 400` | `8.568e-05 → 2.131e-05`, factor **`4.02`** (required ≥ 3) — `O(h²)` confirmed |
| `Sd` endpoints | `6.948e-03` / `3.464e-03` (bound `1e-2`) |
| `Se` **the value IS the derivative of the star's own `ε(p)` spline** | `1.311e-09` / `2.615e-09` (bound `1e-6`) |
| `Pa`–`Pe` profile ownership | §10 |

`Se` is the strongest statement in the file: it pins the returned number to the *same*
interpolant `GetEDens` reads, not merely to the analytic EOS. A second interpolant or a profile
finite difference fails it — demonstrated by detectors D3 and D1 (§11).

The test EOSs are `ε = E₀ + S(p − p_min)` (linear, exact reproduction) and
`ε = E₀(p/P₀)^0.4` (power law, `dε/dp` decreasing with pressure, a genuinely nonzero third
derivative — a quadratic would have been reproduced almost exactly by Steffen's slope formula and
would not have tested discretization at all).

## 9. Real-EOS results — `tests/eos/eos_derivative_cmf.cpp` (external data)

DS(CMF)-1_with_crust, four stars built through `TOVSolver::SolveToProfile`. **17/17 checks
pass**: every star publishes the derivative, one value per radial node, **zero non-finite
values** (hence **zero unsupported-domain nodes** — every profile pressure lay inside the
interpolant's domain, which is what the fail-closed policy of §7 demands) and **zero
non-positive values**.

| M [M☉] | R [km] | nodes | min | median | max | centre | surface | max node step (at r) |
|---|---|---|---|---|---|---|---|---|
| 1.00 | 13.42632 | 2629 | 5.2825e0 | 5.9380e0 | 7.4243e3 | 5.2825e0 | 3.2423e3 | 1.525e0 (12.4603 km) |
| 1.40 | 13.54532 | 2646 | 3.7908e0 | 4.1953e0 | 5.4656e3 | 3.7908e0 | 3.2430e3 | 2.288e0 (12.9013 km) |
| 1.60 | 13.46832 | 2635 | 3.2786e0 | 3.5938e0 | 5.7302e3 | 3.2786e0 | 2.9703e3 | 1.538e0 (12.9503 km) |
| 2.00 | 12.71232 | 2527 | 2.5624e0 | 2.7455e0 | 4.6669e3 | 2.5624e0 | 2.8897e3 | 2.453e0 (12.4043 km) |

**Reading.** The centre value is the minimum on every star and falls monotonically with mass
(`5.28 → 2.56`): denser matter is stiffer, `c_s² = 1/(dε/dp)` rising from `0.19` to `0.39` — and
every value is causal (`dε/dp ≥ 1`) without that ever being asserted. The derivative then rises
by three orders of magnitude through the crust to `≈ 3e3` at the surface, which is the soft outer
crust, and the sharpest features sit at `12.4–13.0 km`, i.e. in the crust, not at the centre and
not at the EOS-floor node. The range reproduces the `3.8–4.9e3` recorded for this EOS in the
Phase-4A-0 audit (`PHASE4_ROTATION_ENTRY.md` §11C).

## 9.1 Diagnostic — why the profile finite difference is not the authority

The retired estimate `(dε/dr)/(dp/dr)`, computed by centred differences on the profile's own
geometric columns exactly as the unratified candidate does (`RotationSolver.cpp:1254-1277`), is
*mathematically the same quantity* by the chain rule on a barotrope. Every deviation below is
therefore **numerical noise in the retired method**, not physics. Reported, never asserted, never
used as an oracle.

| M [M☉] | median \|FD−EOS\|/EOS | max \|FD−EOS\|/EOS | nodes off by > 1 % | max FD node step |
|---|---|---|---|---|
| 1.00 | 5.3556e-07 | **4.8988e0** | 27 / 2629 | 4.506e0 |
| 1.40 | 4.5671e-07 | **3.1875e0** | 26 / 2646 | 3.548e0 |
| 1.60 | 4.2774e-07 | **2.6550e0** | 23 / 2635 | 3.076e0 |
| 2.00 | 3.4236e-07 | **1.5455e0** | 24 / 2527 | 2.699e0 |

**This is the case against it in one table.** In the smooth core the two agree to `~4e-7` — which
is why the defect is easy to miss — but at roughly two dozen crust nodes per star the finite
difference is wrong by **155 % to 490 %**, and its worst node-to-node step is `2.7–4.5` against
the authoritative `1.5–2.5`. A method that is excellent on average and catastrophic exactly where
the EOS is hardest is not an authority; it is a hazard. ADR-0007 P5 removes it, and detector D1
(§11) proves the removal is load-bearing.

## 10. Point-constructed and analytic stars

Existing analytic harnesses build `NStar` from explicit `TOVPoint` data with no live EOS object.
`dε/dp` is **never** reconstructed from those radial profiles. Instead `TOVPoint` gained one
trailing member,

```cpp
double dedp;   // DIMENSIONLESS; NaN (the default) means "not supplied"
```

with a defaulted trailing constructor parameter, so **every existing construction site compiles
unchanged** and no first-order test was forced to supply anything. The NaN default is what makes
the two states impossible to confuse:

- **explicitly supplied by the analytic EOS** — e.g. `0.0` at every node of a constant-density
  star, which is the correct formal derivative for incompressible matter;
- **not supplied** — the profile publishes nothing and any consumer requiring it fails closed.

`TOVPoint::Str()` was deliberately left unchanged, so no exported TOV-point stream gains a column.

**Measured** (`Pa`–`Pe`, on the exact Schwarzschild constant-density interior, `M = 2 km`,
`R = 13 km`, `N = 601`):

| Check | Result |
|---|---|
| `Pa` points supplying nothing ⇒ `HasEosDEdP()` false, `GetEosDEdP()` nullptr | pass |
| `Pb` an analytic star may supply its own | 601 values published |
| `Pc` the supplied constant-density value is **exactly `0.0`** at every node — incompressible, *not* the causal limit `1` | pass |
| `Pd` one interior NaN out of 601 withdraws the whole set | pass |
| `Pe` first-order results **bit-identical** with and without the derivative | `I = 1.571329385870e+02` both; `s(r)` and `s'(r)` bitwise equal at every node |

`Pe` is the guarantee that 4C-I0 added an authority and changed no physics.

## 11. Detectors

Each mutation: one unique production site, rebuilt, measured, then restored and verified
**byte-identically by SHA-256** (`TOVSolver.cpp` → `7a0aa20b…`, `NStar.cpp` → `b2b7948e…`, both
MATCH), with `os.utime` after every write so `make` could not skip a rebuild on a preserved
mtime. Script: scratchpad, untracked.

| | Mutation | Result — fired |
|---|---|---|
| **D1** | publish a radial-profile finite difference instead of the EOS-owned derivative (`NStar::BuildFromTOV`) | contract test **FAILS**: `Pa` (a star that supplied no derivative now publishes one) and `Pd` (the all-or-nothing rule is bypassed). The CMF test still passes — which is itself the finding: on *average* the FD looks fine, and only the fail-closed contract catches it |
| **D2** | omit the `c²` conversion (`EosDerivCgsToDimensionless()` returns `1.0`) | contract test **FAILS** on 5 checks — `Lb`, `Sb`, `Sc`, `Sd`, `Se` — each at relative error `1.000e0`, i.e. saturated: the value is wrong by the factor `c² ≈ 9e20`, exactly the predicted scale |
| **D3** | return the derivative of a **different** interpolation scheme (`gsl_interp_cspline`) over the same table | contract test **FAILS**: `Se` at `1.03e-3` against its `1e-6` bound (the check designed for exactly this), plus `Sb`, `Sc`. The CMF test **also fails on all four stars**: the non-monotone cspline produces **negative** `dε/dp` at 22–66 crust nodes — a physically impossible sound speed that the monotonicity-preserving Steffen authority never produces |

D3 is worth stating plainly: the choice of interpolant is not cosmetic. On this EOS a smoother
scheme yields `1/c_s² < 0` at dozens of crust nodes. `Se` catches the source-of-truth mismatch
analytically; the CMF positivity check catches its physical consequence.

## 12. Artifact protection

| Protected object | Result |
|---|---|
| `RotationSolver.hpp`, `RotationSolver.cpp` | **not touched** — `git diff HEAD` empty for both |
| `ODE_Hartle2_N_Fast`, `SolveHartle2_N`, `ODE_Hartle2_Mixed_Fast`, `SolveHartle2_Mixed` (254 lines) | **byte-identical to `master` (`df859b5`)** — extracted block SHA-256 `3e085cb4932cf6785673a614427c4d3c0f2168cdb81f7e96cd011e1bba6018c8` on both sides |
| `HartleResult` candidate fields | unchanged (same file, untouched) |
| the seven durable artifacts | **byte-identical** — all seven hashes as §1 |
| `I` | bit-identical (`Pe`; both Hartle CTests pass) |
| first-order response tests, TOV reference, passive cooling, path equivalence, grid convergence | pass, unchanged |
| `StarProfile` `radial` schema | unchanged — no column added, species indices still `8 … 8+n_species-1`, `Export` layout unchanged |

No existing baseline was re-established and no scientific output moved.

## 13. `dn_i/dp` — deferred

ADR-0007 P5 applies the same ownership principle to the species derivative `dn_i/dp`, and it is
**not implemented here**. Reasons, all substantive: species semantics are ADR-0001-governed
(`n_i = Y_i n_B`, and the stored columns are *fractions*, so `dn_i/dp` is not simply a spline
derivative of a stored column); `Y_i` itself depends on pressure; Phase 5 must decide the
derivative/composition interface; and the Hartle monopole ODE does not need it — only the
Phase-5 particle-number integral does (`PHASE4C_HARTLE2_DERIVATION.md` §15).

> **`dn_i/dp`: GOVERNED PRINCIPLE ESTABLISHED — IMPLEMENTATION DEFERRED TO PHASE 5.**

`GOVERNANCE.md` §3.1 status: **AUTHORIZED BY ADR-0007 §9 — CORRECTION NOT YET EXECUTED.** No
`hartle_monopole` baseline exists, and none may be created before 4D.

## 14. Readiness gate for 4C-I1

4C-I1 may begin when, and only when, all of the following hold — every one is true as of this
increment:

1. ADR-0007 is ACCEPTED with the Decision section recorded — ✅ `13111de`;
2. `dε/dp` has a single authority, is dimensionless, and its unit conversion is asserted against
   an independent constant — ✅ §4, check `Ua`;
3. that authority is the derivative of the star's own `ε(p)` interpolant, proven against the
   spline itself and not merely against an analytic EOS — ✅ check `Se`, detector `D3`;
4. its domain semantics fail closed, with NaN distinguishable from the physical value `0` —
   ✅ §7, checks `Da`–`De`;
5. a completed `StarProfile` can supply it at every node, all-or-nothing — ✅ §5, checks
   `Pa`–`Pd`, real stars in §9;
6. point-constructed analytic stars can supply their own explicitly, and cannot be confused with
   a reconstruction — ✅ §10;
7. first-order and TOV physics are unchanged and do not require it — ✅ `Pe`, §12;
8. the O(Ω²) candidate is still byte-identical and still marked nonconformant — ✅ §12, INV-08.

**What 4C-I1 must then do** (ADR-0007 §4, §9): implement the governed `l = 0` system (N1)–(N2)
in `p̂₀*` from the seed-free `s`, `s'` and this derivative; fixed-`ε_c` centre condition with the
regular-series start; `δM̂` with the exterior and surface-shell terms; the seed-free result type
with `At(AngularVelocity)` materialization; the `NStar` accessor; **fail closed when
`HasEosDEdP()` is false**; and delete the candidate atomically in the same commit. Then 4D.

## 15. Exact non-scope

No O(Ω²) equation was implemented, repaired, executed or baselined. No candidate was deleted or
modified. No second-order result type or `NStar` accessor was created. No `dn_i/dp`. No `l = 2`.
No MixedStar work (`FindMixedMomInertia`, `Solve_Mixed`, the dark-sector splines and
`fast_dEdP_v/_d` are untouched; the dark EOS has no derivative accessor and does not need one
yet). No change to the TOV surface convention, INV-06, `Geometry`, thermal drivers,
`RotochemicalCache`, or the BNV/Decay/DarkCore modules. No baseline re-established.

## 16. Status

> **`PHASE-4C-I0 EOS DERIVATIVE AUTHORITY COMPLETE`**

**Exact next task.** Phase 4C-I1 — atomically replace the invalid public Hartle2 candidate with
the ADR-0007 fixed-`ε_c`, `Ω²`-normalized monopole response using the governed EOS derivative,
while preserving all first-order results bitwise.
