# Phase 4D-RG — the EOS energy-density source of the Hartle monopole as a measure

> **FORMAL STATUS: `PHASE-4D EOS-MEASURE CORRECTION ADR PROPOSED — OWNER ADJUDICATION REQUIRED`**
>
> This record derives, from Hartle's own field equation, why the Phase-4D validation failed
> (`PHASE4D_MONOPOLE_VALIDATION.md`: `HARTLE MONOPOLE VALIDATION FAILED`, δM̂ ≈ −4.6 % on
> DS(CMF)-1), evaluates three numerical authorities for the correction in scratch, and grounds
> `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md` (**PROPOSED**). Nothing in
> production, tests or CMake is changed; no baseline is created; ADR-0007's accepted Decision is
> untouched.

| Field | Value |
|---|---|
| **Starting HEAD** | `246f3f24f84079e3fb928c2180e5fef7388a989d` (4D), branch `physics/rotation-correctness`, upstream equal, clean, **9 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Change class** | **documentation** (governance / primary derivation); scratch numerics outside the tree |
| **Pre-task suites** | full **30/30** (221.43 s), self-contained **16/16** (16.17 s); seven artifact hashes as `PHASE3_CLOSEOUT.md` §1 |
| **Primary sources re-read** | Hartle (1967) ApJ 150, 1005 — pp. 1009–1010 (eqs. 8–19), 1019 (83–91), 1020 (92–98), 1021 (99–107), 1022 (108); Chandrasekhar & Miller (1974) MNRAS 167, 63 — p. 69 (eqs. 50–59) |
| **Scratch** | `scratchpad/4dr/scratch_4dr.cpp` (DS(CMF)-1, HW EOS, analytic star; five source formulations), `scratchpad/4dr/twolayer.cpp` (exact two-layer star; jump operator) — not tracked |

---

## 1. The failure being corrected

Phase 4D established that `RotationSolver::ODE_HartleMonopole_` integrates ADR-0007 P2 exactly
(independent `(m₀,h₀)` and continuum solvers, published tables, detectors) but that ADR-0007 §7
item 11 is **not met** on DS(CMF)-1: the homogeneous `δM̂` misses the non-rotating sequence
derivative by `1.04e-3` at every resolution, the nodal `dε/dp` column integrated over the crust
misses ≈ 17 % of the crust's own `Δε`, and the omission on the *sourced* `δM̂` was estimated at
≈ 4.6 %. The dominant feature is the table's crust–core transition: `ε` rises 36 % over a 1.5 %
pressure interval — in the star a layer ≈ 0.44 m thick, narrower than any production node spacing
(§5). ADR-0007 P2 carries the EOS term as `4πr²(ε+p)(dε/dp)p̂₀*` with `dε/dp` a **column sampled
at the profile nodes and linearly interpolated** (INV-13). That representation is not
measure-complete: whatever `Δε` falls between two samples of a *derivative* is lost.

Two things are kept apart throughout, as the task requires: the **EOS derivative authority**
(`TOVSolver::GetEDensDeriv`, 4C-I0) is correct pointwise — scratch measured its nodal values
against the profile to `3e-14` (§5) — and is retained (§10 of ADR-0008, Q8); what fails is the use
of a **nodal derivative column as the representation of the source measure**.

## 2. The mathematical measure derivation (primary source)

Hartle's `l = 0` field equation (H67 eq. 93, p. 1020) carries the matter source as

```
(ΔG_t^t)_{l=0} = −(2/R²) dm₀/dR + [rotational terms] − 8π ξ dE/dR ,
```

i.e. the Eulerian density change of the displaced constant-density surfaces enters as
`−ξ₀ dE/dR`, and with (95) the `t`–`t` Einstein equation gives

```
dm₀/dR = −4πR² ξ₀ (dE/dR) + (1/12) j² R⁴ ω̄_R² − (1/3) R³ (j²)_R ω̄² .          (M1)
```

Hartle then substitutes (88)–(99), `p₀* = −ξ₀ (dP/dR)/(E+P)`, together with
`dE/dR = (dE/dP)(dP/dR)`, to write the first term as `4πR²(dE/dP)(E+P)p₀*` — his eq. (97). That
rewrite is an identity **only where `E(R)` is differentiable**. The primary statement is (M1),
whose matter term is a *measure*:

```
dm̂₀ |_EOS = −4π r² ξ̂₀ dε        (Lebesgue–Stieltjes in ε along the background)              (M2)
```

in CompactStar's normalization (`ξ̂₀ = ξ₀/Ω²`, `m̂₀ = m₀/Ω²`). **Sign conventions:** `r`
increases outward, `ε` decreases outward (`dε ≤ 0`), `ξ̂₀ > 0` is an outward displacement of the
isobar; the increment `−4πr²ξ̂₀ dε` is therefore **positive** for a layer moved outward across a
falling density — the extra mass enclosed by a fixed coordinate sphere. The same statement is
Hartle's Newtonian eq. (18), p. 1010: `δM = 4π ∫₀^a (−ξ₀ dρ/dR) R² dR = M^{[2]}(a)`, with his own
warning on the same page that "near the surface where the Eulerian expansion is no longer valid
the equations remain formally the same, but the identification of `p*` with the ratio of the
change in pressure to the unperturbed density can no longer be made".

Derivation of (M2) from ADR-0007's own variables, for the record: with the TOV relation
`dp/dr = −(ε+p)ν'` and `ξ̂₀ = p̂₀*/ν'`,

```
4πr²(ε+p)(dε/dp) p̂₀* dr = 4πr² p̂₀* (dε/dp) [−dp/(ν')]·… = −4πr² (p̂₀*/ν') (dε/dp) dp = −4πr² ξ̂₀ dε .
```

So the differential term of P2 and the measure (M2) coincide **wherever `dp/dr = −(ε+p)ν'` holds
pointwise and `ε` is absolutely continuous**. Neither holds inside a sub-node feature of a
piecewise-linear background — which is exactly where the two forms were found to differ (§6).

## 3. Smooth-equivalence proof

Let `ε(r)` be absolutely continuous on a segment `[a,b]` and let the background satisfy the TOV
relation there. Then `∫_a^b −4πr²ξ̂₀ dε = ∫_a^b 4πr²(ε+p)(dε/dp)p̂₀* dr` exactly (the substitution
of §2 is invertible). On the tabulated background the relevant statement is discrete: for a
segment on which `ε` is represented as **linear in `r`** between its end values `ε_a, ε_b`,

```
∫_a^b −4πr² ξ̂₀ dε = −(ε_b − ε_a)/(b − a) · ∫_a^b 4πr² ξ̂₀ dr ,                             (M3)
```

an ordinary ODE source `−4πr²(p̂₀*/ν') · Δε/Δr` that is (i) a bona fide term of the coupled
system (no operator splitting), (ii) **exact in the total measure** `Δε` of the segment whatever
happens inside it, and (iii) `O(h²)` accurate in the weight on smooth segments. Scratch
confirmation (§8): on the smooth HW EOS all five formulations agree with the differential column
form to `≤ 1.25e-5`; on the analytic constant-density star the secant form reproduces the column
form **exactly** (`0.0`), the interior source vanishing identically and the surface being the
terminal atom.

## 4. The jump condition

At a background radius `r_t` where `ε` is discontinuous with `p` continuous (a first-order
transition at pressure `p_t`; outward `ε⁻ → ε⁺`, `ε⁻ > ε⁺`), (M2) integrates across the atom to

```
Δm̂₀(r_t) = 4π r_t² (ε⁻ − ε⁺) ξ̂₀(r_t) ,                                                    (M4)
```

with **`p̂₀*` continuous** (its equation (100)/(98) contains `m̂₀`, which jumps finitely, but no
delta), **`ĥ₀` continuous** (98), hence `ξ̂₀ = p̂₀*/ν'` continuous because `ν' = (m + 4πr³p)/(r(r−2m))`
is continuous through a density jump (only `dp/dr` jumps). Multiple transitions compose
additively, each applied when the integration reaches its radius. The governed surface shell of
ADR-0007 P6, `4πR_*²ε_*ξ̂₀(R_*)`, is (M4) with `ε⁻ = ε_*`, `ε⁺ = 0` — the **terminal atom of the same
measure**, not a separate rule. Certification: §4.1 (two-layer exact star, `2e-10 … 4e-9`).

### 4.1 True first-order transition — certified on an exact two-layer star

Incompressible core `ρ₁ = 1e-3 km⁻²`, mantle `ρ₂ = 0.6ρ₁`, transition pressure `p_t = 0.05ρ₁`
(a genuine constant-pressure jump); core closed-form (Schwarzschild interior), mantle integrated
in pressure to `1e-13`; oracle `dM/dp_c` by central differences (`±1e-4 p_c`); homogeneous
`(m₀,h₀)` solve with `p̂₀*(0) = 1`, `dε/dp = 0` inside both layers, (M4) applied at `r_t`, terminal
atom at `R`:

| `p_c/ρ₁` | `M/R` | jump / `δM̂` | `δM̂` with (M4) vs sequence | without (M4) |
|---|---|---|---|---|
| 0.1 | 0.149 | 0.49 | **`4.1e-9`** | `5.2e-1` |
| 0.2 | 0.212 | 0.46 | **`2.3e-10`** | `4.9e-2` |
| 0.5 | 0.311 | 0.47 | **`3.9e-9`** | `2.5e-1` |

`p̂₀*` and `ĥ₀` continuous across `r_t`, `m̂₀` jumping by (M4), reproduce the exact sequence
derivative to integration accuracy. Smoothing the jump is neither needed nor admissible.

## 5. DS(CMF)-1 diagnosis, reproduced independently

Fresh scratch (`scratch_4dr`, 1.6 M☉, `ε_c = 7.312533e14 g/cm³`), radial resolution
1250/5000/10000/20000/40000 (res 2500 excluded, see below):

| res | nodes | spacing | column integrability deficit whole / crust | crust–core layer in `r` | knots in the densest profile interval |
|---|---|---|---|---|---|
| 1250 | 331 | 41 m | `1.6e-2` / `1.3e-1` | 0.44 m | 404 |
| 5000 | 1319 | 10 m | `1.9e-2` / `1.4e-1` | 0.44 m | 275 |
| 10000 | 2635 | 5.1 m | `2.4e-2` / `1.8e-1` | 0.48 m | 167 |
| 20000 | 5268 | 2.6 m | `2.3e-2` / `1.7e-1` | 0.41 m | 104 |
| 40000 | 10535 | 1.3 m | `2.3e-2` / `1.7e-1` | 0.44 m | 58 |

- The profile's `ε_i` equal the EOS interpolant's `ε(p_i)` to **`5e-16`**, and the profile's
  `dε/dp` column equals `GetEDensDeriv(p_i)` to **`3e-14`** — the star's own EOS representation,
  no second EOS.
- The table has **no vertical interval** (1191 rows, pressure strictly increasing); the
  crust–core near-jump (`ε` 4.93e13 → 6.71e13 g/cm³ over `Δp/p = 1.5 %`) is a *continuous*
  Steffen segment ≈ 0.44 m thick in the star at every resolution — case B of the task's
  taxonomy, not case C.
- Sequence derivative (oracle, the solver's own central `p_c` read back): `dM/dp_c = 9671.2 km³`
  at every valid resolution; column-form homogeneous identity `2.2e-5, 5.8e-4, 1.17e-3, 1.02e-3,
  1.04e-3` (erratic — node placement relative to the layer), reproducing 4D.
- **Res 2500 anomaly (recorded, out of scope):** `SetRadialRes(2500)` produced a *different*
  non-rotating star (`dM/dp_c = 9734`, no crust knot bracketed by any node pair, `m̂₀(R_*) = 698`
  vs ≈ 816–855) — a TOV-layer resampling artefact at that one resolution, not a Hartle matter;
  flagged for the TOV layer, excluded from every statement below.

## 6. Option A — pointwise EOS derivative at the actual RHS state

`term₁ = 4πr²(ε_lerp+p_lerp) · dε/dp_EOS(p_lerp(r)) · p̂₀*` evaluated inside the ODE (variant
POINT), at production tolerance (`1e-10`) and tighter (`1e-13`), with and without every mapped
EOS-knot radius as a mandatory step boundary (POINT+K), and with `ε` also from the EOS at the RHS
state (POINT-E):

| res | COLUMN `δM̂` | POINT(1e-10) | POINT(1e-13) | POINT+K | POINT-E | SECANT (§8) |
|---|---|---|---|---|---|---|
| 1250 | 855.475 | 914.908 | 914.908 | 914.908 | 932.496 | 866.084 |
| 5000 | 839.278 | 876.182 | 876.182 | 876.182 | 875.946 | 865.868 |
| 10000 | 825.970 | 870.734 | 870.734 | 870.734 | 869.930 | 865.866 |
| 20000 | 828.577 | 859.804 | 859.804 | 859.804 | 863.507 | 865.836 |
| 40000 | 827.999 | 865.590 | 865.590 | 865.590 | 866.235 | 865.846 |

Findings. (a) The GSL rk8pd driver **does resolve the 0.44-m layer**: tolerance and knot
boundaries change nothing (`≤ 1e-8`). (b) Yet Option A is **not convergent**: successive
differences `13, 5.5, 11, 5.8` (≈ 1–2 %), erratic in node placement — worse than the column
form. (c) POINT-E (both `ε` and `dε/dp` from the EOS) is equally erratic, so the defect is **not**
the lerped `ε`: it is the lerped `dp/dr`. Inside a sub-node layer the piecewise-linear background
has `dp_lerp/dr ≠ −(ε+p)ν'` by `O(1)` (the true slope is 36 % larger where `ε` has jumped), and
(97)'s `(ε+p)(dε/dp)dr` integrates to `Δε·[(ε+p)ν'/(−dp/dr)]` rather than `Δε`. (d) The
homogeneous identity improves with resolution (`3.7e-3 → 6.6e-5`) only because the layer's share
of the interval shrinks. **Verdict:** Option A recovers the transition numerically but is
formulation-inconsistent with INV-13's background wherever a feature is narrower than the node
spacing; it cannot represent a true discontinuity at all (a delta in `dε/dp` is not evaluable).

## 7. Option B — smooth derivative plus explicit discontinuity shells

Identification of "discontinuities" is the whole difficulty. DS(CMF)-1 has none in the strict
sense (§5): the importer requires strictly increasing pressure (`TOVSolver.cpp:666–675`) and the
Steffen interpolant is C¹, so a rule keyed on duplicate pressures or vertical table intervals
fires **nowhere** and leaves the 4D failure intact; a rule keyed on a steepness threshold is
arbitrary and forbidden by the task. Where genuine transition metadata exists (a future EOS-layer
authority: `p_t, ε⁻, ε⁺`), the jump operator (M4) is the correct treatment — and it is *also* what
the measure form does at an atom. Option B is therefore **subsumed by Option C**: its shells are
the singular part of the measure; its differential part is the absolutely-continuous part, which
on a tabulated background must itself be measure-complete (§8), not a sampled derivative.

## 8. Option C — the measure form, two realizations

**SECANT (profile partition):** on each profile interval the source is `−4πr²(p̂₀*/ν')·Δε_i/Δr_i`
with `Δε_i = ε_{i+1} − ε_i` from the profile's own `ε` nodes — (M3) on INV-13's own linear
representation of `ε(r)`; no EOS access. **SECANT+K (knot partition):** the same with every
EOS-table knot mapped into radius through the interval's linear `p(r)` and `ε` at the knot from
the interpolant (which passes through the table values), so the sub-node location of the layer
is honoured.

| res | SECANT `δM̂` | SECANT+K | SECANT J-identity | SECANT+K | COLUMN J | POINT J |
|---|---|---|---|---|---|---|
| 1250 | 866.084 | 866.447 | `4.6e-4` | `8.8e-4` | `2.2e-5` | `3.7e-3` |
| 5000 | 865.868 | 865.896 | `1.1e-4` | `1.4e-4` | `5.8e-4` | `6.5e-4` |
| 10000 | 865.866 | 865.856 | `1.0e-4` | `9.1e-5` | `1.17e-3` | `2.0e-4` |
| 20000 | 865.836 | 865.854 | `5.7e-5` | `7.8e-5` | `1.02e-3` | `1.1e-4` |
| 40000 | 865.846 | 865.848 | `7.2e-5` | `7.4e-5` | `1.04e-3` | `6.6e-5` |

- Measure completeness: the secant form captures the transition at every resolution, including
  41-m spacing; the corrected sourced `δM̂` at the production resolution is **865.87 vs 825.97,
  +4.83 %** (4D's post-hoc estimate was +4.6 %); the shell-excluded interior `m̂₀(R_*)` moves
  from 815.5 to 855.4.
- Convergence: relative spread `4e-5` over 5000 → 40000 and **no node-placement behaviour** (the
  column form's spread is `1.6 %`, Option A's `2 %`); the homogeneous identity sits at the
  oracle's floor (`≈ 7e-5`, set by the central-difference sequence derivative) from res 5000 on.
- Knot refinement: agrees with the profile-partition form to `3e-5` at ≥ 5000; at 41-m spacing
  both carry an `O(h)` weight-location error (`4.6e-4`, `8.8e-4`) — the refinement fixes the
  *location* but not the lerped background's own `O(h)` inconsistency inside a 41-m interval, so
  it helps only when the weight `4πr²ξ̂₀/ν'` varies materially across an interval. It needs an
  immutable EOS snapshot (§10) and is optional at production resolutions.

## 9. Numerical convergence comparison

| Formulation | Total measure captured | Order on smooth segments | Behaviour at a sub-node feature | Node-placement sensitivity | True discontinuity |
|---|---|---|---|---|---|
| COLUMN (ADR-0007 as implemented) | **no** (17 % of crust `Δε` lost) | `O(h²)` | loses the feature | erratic, 1.6 % | impossible |
| POINT (Option A) | yes, but mis-weighted | `O(h²)` | `O(1)` local inconsistency of `dp_lerp/dr` | erratic, 2 % | impossible |
| SECANT (Option C, profile) | **yes** | `O(h²)` | exact `Δε`, `O(h)` weight location | none measured (`4e-5`) | as a steep segment (exact `Δε`) or as an explicit atom |
| SECANT+K (Option C, knots) | **yes** | `O(h²)` | exact `Δε` and location | none measured | atom |
| jump operator (M4) | — | — | — | — | **exact** (`2e-10…4e-9`, §4.1) |

## 10. EOS ownership — what object supplies the measure

Measured facts: the profile's `ε` column **is** the EOS interpolant at the node pressures
(`5e-16`); the knots are the table rows the importer holds (`TOVSolver::eos_tab`, protected).
Consequences:

- **Candidate A (profile `ε` nodes)** supplies the absolutely-continuous part on the profile
  partition with no new object, no new ownership and no provenance change: it is the same
  representation that built the star. Sufficient at production resolutions (§8).
- **Candidate B (the `ε(p)` table)** is needed only for the knot refinement: an **immutable
  value-type snapshot** `(p_k, ε_k)` attached to `StarProfile` at construction (like the `dε/dp`
  column, all-or-nothing, `Touch()` on any later supply) — no GSL handles, no `TOVSolver` lifetime,
  no owning EOS object in `RotationSolver`. Optional; point-constructed stars omit it.
- **Candidate C (transition metadata `p_t, ε⁻, ε⁺`)** does not exist in the EOS layer today (the
  importer rejects non-increasing pressure; `gsl_interp_steffen` requires strictly increasing
  abscissae). A true first-order transition therefore **cannot yet be represented** by the star
  builder either; ADR-0008 fixes the semantics (M4) so that when the EOS layer gains such a
  representation the Hartle side is already governed. Stated plainly: *authoritative jump metadata
  is absent*.

## 11. Surface-shell unification

The explicit ADR-0007 P6 shell is the terminal atom `ε_* → 0` of (M2) at `R_*`. Proposed
semantics (ADR-0008 Q7): the interior measure runs over `[r₀, R_*)`; the terminal atom is added
**once**, computed by the same operator (M4) with `ε⁺ = 0`, and reported as the existing
`surface_shell_mass_over_Omega2` field — one formula, no double counting, byte-identical to today
on stars whose interior measure is smooth (analytic star: `0.0`; HW: `≤ 1.3e-5`).

## 12. Phase-5 consequence

H67 (109)–(114) build the baryon number the same way; ADR-0007 P13's
`(dn_i/dp)δp̂₀ = −(dn_i/dr)ξ̂₀` is the smooth rewrite of a measure `−4πr²e^{λ}ξ̂₀ dn_i`. A nodal
`dn_i/dp` column would repeat this failure at every composition discontinuity (and the crust is
made of them). ADR-0008 Q11 requires the scalar particle-number response to be measure-complete
from the outset; nothing of Phase 5 is implemented here.

## 13. Provenance

The corrected response consumes only profile-attached data (`ε`, `p`, `ν'`, optionally the knot
snapshot and, later, transition metadata) supplied at construction or through `Touch()`-ing
setters (`StarProfile.hpp:926–990`: `SetEosDEdP`, `AppendEosDEdP`, `ClearEosDEdP` all call
`Touch()`). The existing `HartleMonopoleResponse` provenance — profile identity + `Version()`
(`RotationSolver.hpp:236–237, 317`) — therefore remains sufficient **provided every new measure
datum follows the same rule**; no second cache authority is introduced.

## 14. Validation plan for the correction (predeclared)

See ADR-0008 §Validation: A–J of the task, with the homogeneous DS(CMF)-1 identity target
**`≤ 2e-4`** (twice the worst measure-form value at production resolution, `1.03e-4`; the oracle
floor is ≈ `7e-5`), a sourced Stieltjes/source identity, H re-run with no node-placement
behaviour, the new detector "omit the interior measure" that must reproduce the percent-level
failure, and the seven artifacts byte-identical. Only then the first monopole baseline.

## 15. Owner questions

ADR-0008 §11 Q1–Q12 (canonical form; source of truth; algorithm; sharp segments; true
discontinuities; transition location; shell unification; role of `dε/dp`; provenance;
point-constructed stars; Phase-5 measure completeness; correction/validation/baseline policy).

## 16. Exact non-scope

No production change (`RotationSolver`, `NStar`, `StarProfile`, `TOVSolver` untouched); no test
or CMake change; no baseline; no Phase-5 code; no `dn_i/dp`; no `l = 2`; no repair of the
res-2500 TOV anomaly; no edit of `GOVERNANCE.md` (its "currently only use" wording for §3.1 is
stale — ADR-0007 is a second use — recorded as documentation debt, not corrected here); ADR-0007's
Decision unchanged.
