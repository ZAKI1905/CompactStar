# ADR-0008 — Measure-complete EOS energy-density source for the Hartle monopole response

| Field | Value |
|---|---|
| **Status** | **PROPOSED** (2026-09-03) — owner adjudication required before any production correction |
| **Date** | 2026-09-03 |
| **Change class** | scientific-semantic (source measure of H67 eq. 97) and numerical-method (its integration); the production correction it authorizes proceeds under `GOVERNANCE.md` §3.1 as ADR-0007 §9 already does |
| **Governing authority** | Hartle (1967) ApJ 150, 1005, eqs. (18), (88), (93), (97)–(100), (107)–(108); Hartle & Thorne (1968); Chandrasekhar & Miller (1974) eqs. (50)–(56); ADR-0007 (ACCEPTED); ADR-0003 (provenance); ADR-0005 (TOV primitive); INV-06, INV-08, INV-09, INV-13 |
| **Affected invariants** | INV-08 (status), INV-09 (Phase-5 consequence), INV-13 (a tabulated background's representation of a measure) |
| **Blocks** | the Phase-4D-J correction increment, the first monopole baseline (§3.1 condition 7), roadmap 4E, Phase 5 |
| **Evidence record** | `docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md`; failure: `docs/validation/PHASE4D_MONOPOLE_VALIDATION.md` |

> **Relation to ADR-0007.** This ADR **amends / supersedes ADR-0007 P2, P5 and P6 only where EOS
> energy-density variation is represented by a measure rather than a pointwise radial
> derivative.** It leaves unaffected: `p₀*` as the canonical variable (P1); the fixed-`ε_c` family
> and the regular-centre start (P3, P4); the `Ω²` normalization and explicit-spin materialization
> (P9, P10); the first-order response and `I`; the surface `R_*` convention (P7, P8); the `l = 0`
> roadmap scope (P12); the exterior term `I²/R_*³` (P6). ADR-0007's accepted Decision is not
> edited; it receives only a post-validation note pointing here.

## Context

**What the code does.** `RotationSolver::ODE_HartleMonopole_` (`CompactStar/Core/src/RotationSolver.cpp:1251–1292`)
evaluates the EOS term of H67 (97) as `4π r² (ε+p) · dedp · p̂₀*` with `dedp` linearly interpolated
from the profile column that 4C-I0 fills from `TOVSolver::GetEDensDeriv` at the nodes
(`NStar.cpp:284, 608, 806`). Phase 4D verified that integration against two independent solvers,
the continuum limit, analytic limits and two published second-order tables, and found ADR-0007 §7
item 11 **unmet on DS(CMF)-1** (`1.04e-3` vs `1e-3`, resolution-independent) with the sourced `δM̂`
deficient by ≈ 4.6 %.

**What the physics says.** Hartle's field equation (93) carries the matter source as
`−8πξ dE/dR`: the increment of `m₀` from the displacement of constant-density surfaces is the
**measure** `−4πr²ξ₀ dε` (evidence record §2). Eq. (97)'s `(dE/dP)(E+P)p₀*` is that measure
rewritten through (88) where `E` is differentiable. ADR-0007 P2 adopted the rewrite, P5 gave the
derivative a single authority, and P6 handled the *surface* atom of the measure explicitly — but
**internal** variation narrower than the profile spacing has no representation: on DS(CMF)-1 the
crust–core transition (`Δε/ε = 36 %` over `Δp/p = 1.5 %`, ≈ 0.44 m in the star, a continuous
Steffen segment) lies between nodes at every resolution, and the nodal derivative column
integrated over the crust misses ≈ 17 % of the crust's own `Δε`.

**What the scratch experiments established** (evidence record §4.1, §5–§9; none of these numbers is
a contract value):

| Formulation of the EOS term | DS(CMF)-1 `δM̂`, res 5000 … 40000 | homogeneous identity | smooth HW EOS vs today | analytic star |
|---|---|---|---|---|
| nodal derivative column (today) | 839 → 826 → 829 → 828 (erratic, 1.6 %) | `5.8e-4 … 1.2e-3` | — | — |
| pointwise EOS derivative at the RHS state (Option A) | 876 → 871 → 860 → 866 (erratic, 2 %); tolerance- and knot-independent | `6.5e-4 → 6.6e-5` | `1.2e-5` | — |
| **measure form, profile partition (Option C)** | **865.87 → 865.87 → 865.84 → 865.85** (`4e-5`) | **`1.1e-4 → 7e-5`** (oracle floor) | `1.2e-5` | **identical (`0.0`)** |
| measure form, EOS-knot partition | 865.90 → 865.86 → 865.85 → 865.85 | `1.4e-4 → 7e-5` | `1.3e-5` | — |
| jump operator on an exact two-layer star | — | **`2e-10 … 4e-9`** (`5–52 %` without it) | — | — |

Option A resolves the layer numerically but is inconsistent with INV-13's piecewise-linear
background inside it (`dp_lerp/dr ≠ −(ε+p)ν'` by `O(1)` there, so `(ε+p)(dε/dp)dr` no longer
integrates to `Δε`); Option B has no threshold-free way to name the DS(CMF)-1 transition a
"discontinuity" (it is continuous) and is subsumed by the measure form at genuine atoms.

**Distinction kept.** The EOS derivative authority (`GetEDensDeriv`, 4C-I0) is correct — the
profile's `dε/dp` nodes match it to `3e-14` and the profile's `ε` nodes match the interpolant to
`5e-16`. What is not measure-complete is the *nodal column as the representation of the source*.

## Decision

*Empty while PROPOSED.* Filled in only on ratification (owner adjudication of Q1–Q12).

## Decision questions

- **Q1 — Canonical mathematical form.** (a) the differential derivative form (ADR-0007 P2 as is);
  (b) hybrid: derivative on continuous regions + explicit shells at declared discontinuities;
  (c) **the Lebesgue–Stieltjes measure** `dm̂₀|_EOS = −4πr²ξ̂₀ dε` for all `ε` variation, reducing
  to the derivative form on absolutely-continuous segments where the background is
  TOV-consistent, to a shell at an atom, and to a finite exact contribution across a very sharp
  continuous table segment. **Recommended: (c)** — from H67 (93) directly, threshold-free, and the
  only form that converged in scratch.
- **Q2 — Source of truth for the measure.** (a) the radial `StarProfile` `ε` nodes (= the star's
  own EOS interpolant at the node pressures, `5e-16`); (b) an immutable snapshot of the `ε(p)`
  table `(p_k, ε_k)` attached to the profile at construction, for knot refinement; (c) explicit
  transition metadata `(p_t, ε⁻, ε⁺)` — **absent from the EOS layer today** (importer rejects
  non-increasing pressure; Steffen needs strictly increasing abscissae). **Recommended: (a) as the
  mandatory representation, (b) optional and immutable, (c) governed here for the day the EOS
  layer can express it.** `RotationSolver` acquires no owning EOS object.
- **Q3 — Integration algorithm.** The measure's absolutely-continuous part as a genuine ODE term
  on each partition segment, `−4πr²(p̂₀*/ν')·Δε_k/Δr_k` (the derivative of the segment's linear
  `ε(r)`), integrated by the existing rk8pd driver with segment boundaries as mandatory step
  limits — **no operator splitting**; the singular part as the jump operator (M4) at declared
  atoms; the terminal atom at `R_*`. Measured: `O(h²)` on smooth segments, exact segment measure,
  no node-placement behaviour (evidence record §8–§9).
- **Q4 — Sharp continuous table segments.** Captured exactly in total measure by Q3 whatever their
  width; the weight's sub-segment location error is `O(h)` and vanishes with the optional knot
  partition (Q2b). No steepness threshold exists anywhere in the contract.
- **Q5 — True discontinuities.** `Δm̂₀ = 4πr_t²(ε⁻ − ε⁺)ξ̂₀(r_t)`; `p̂₀*`, `ĥ₀`, `ξ̂₀` continuous;
  multiple transitions compose additively; no smoothing. Certified `2e-10 … 4e-9` on an exact
  two-layer star.
- **Q6 — Transition-location ownership.** `r_t` is the radius at which the profile's (INV-13
  linear) `p(r)` equals the EOS-declared `p_t` — the same mapping that places table knots; the
  EOS layer owns `p_t, ε⁻, ε⁺`, the profile owns `r_t`.
- **Q7 — Surface shell / double-count rule.** **One unified measure:** interior over `[r₀, R_*)`,
  terminal atom `ε_* → 0` at `R_*` added exactly once by the same operator and still reported as
  `surface_shell_mass_over_Omega2`. Byte-identical to today on smooth interiors.
- **Q8 — Role retained by `dε/dp`.** The 4C-I0 authority stays the EOS derivative for: the
  regular-centre series (`b₅` needs `(dε/dp)_c` pointwise — correct), per-segment smooth
  cross-checks (`Δε_k/Δr_k` vs `(dε/dp)(dp/dr)` where the EOS is smooth), independent
  diagnostics, and future physics needing a pointwise derivative. **It is no longer the
  integrator's source input.** 4C-I0 is not wrong; the nodal column as a measure representation
  was.
- **Q9 — Provenance / cache.** All measure data are profile-attached and supplied at construction
  or through `Touch()`-ing setters (as `SetEosDEdP` already is), so profile identity + `Version()`
  (ADR-0003) remains sufficient; **any new measure datum must follow the same `Touch()` rule** —
  no second authority.
- **Q10 — Point-constructed / analytic stars.** The profile `ε` nodes define the measure (constant
  density ⇒ zero interior source, terminal atom = shell), so nothing new must be supplied; the
  explicit `dε/dp` column remains required for the centre series. Verified: secant ≡ column
  exactly on the analytic star.
- **Q11 — Phase-5 measure completeness.** The scalar particle-number response must be
  measure-complete from the outset: `δN̂_i` carries `−4πr²e^{λ}ξ̂₀ dn_i` (Stieltjes in `n_i`) with
  atoms at composition discontinuities; a nodal `dn_i/dp` column is **forbidden** as the
  integrator's source. Not implemented here.
- **Q12 — Correction / validation / baseline policy.** The production correction is a
  scientific-semantic change executed under §3.1 (this ADR, once ACCEPTED, is its authority; the
  current `246f3f2` response must not be baselined — it omits a ≈ 4.8 % physical term on the
  project's own EOS); re-validation per §Validation; **only then** the first
  `hartle_monopole_dscmf1_debug.tsv`.

## Alternatives

### Alternative A — pointwise EOS derivative at the actual ODE state

- **Statement** — keep H67 (97) as the ODE, but evaluate the same Steffen `dε/dp` (and `ε`) at the
  RHS pressure instead of interpolating a nodal column.
- **Required code changes** — `RotationSolver` needs a dependency-neutral immutable evaluator of
  `dε/dp(p)` (a frozen copy of the interpolant), `ODE_HartleMonopole_` calls it.
- **Migration risk** — none on smooth EOS (`1.2e-5` from today on HW); but on DS(CMF)-1 the answer
  wanders by 2 % with node placement and cannot represent an atom.
- **Validation needed** — would fail the convergence line (H) of §Validation by construction.
- **Implications** — rejected: numerically resolved, formulation-inconsistent (evidence §6).

### Alternative B — derivative column plus explicit shells at declared discontinuities

- **Statement** — keep the column; add `4πr_i²Δε_iξ̂₀(r_i)` at discontinuities identified by
  the EOS layer.
- **Required code changes** — EOS-layer transition metadata (does not exist), a jump operator in
  the solver.
- **Migration risk** — with no metadata the DS(CMF)-1 failure is untouched; any steepness
  threshold is arbitrary.
- **Validation needed** — the two-layer certification (passes) plus DS(CMF)-1 (fails).
- **Implications** — its shell is the atom of Alternative C; its smooth part must itself be
  measure-complete. Subsumed.

### Alternative C — measure-complete Stieltjes source (recommended)

- **Statement** — Q1(c), Q3, Q5, Q7 above.
- **Required code changes** (for the later correction increment, not this ADR) —
  `ODE_HartleMonopole_` term 1 replaced by the per-segment secant source; the driver loop made
  segment-aware (profile nodes; optional knot radii; declared atoms with the jump operator);
  `StarProfile` optionally carrying an immutable `(p_k, ε_k)` snapshot and, when the EOS layer can
  supply it, a transition list; the shell computed by the same operator. `TOVSolver` untouched
  unless the snapshot is exported from it.
- **Migration risk** — smooth-interior stars change by `≲ 1e-5` (HW) or not at all (analytic);
  DS(CMF)-1 `δM̂` changes by ≈ +4.8 %, `ξ̂₀(R_*)` by ≈ −1 %, which is the point.
- **Validation needed** — §Validation A–J.
- **Implications for existing outputs** — no monopole output has ever been baselined; the 4C-I1 /
  4D CMF diagnostic values are superseded and must not be used as expected values.

## Consequences (once accepted)

- The EOS energy-density source of the monopole is a measure; the nodal `dε/dp` column is
  demoted to diagnostics and the centre series; the surface shell becomes the terminal atom of
  one operator.
- The Phase-5 `δN̂_i` interface (ADR-0007 P13) must be re-expressed in measure form before
  implementation.
- The EOS layer gains a governed representation target for first-order transitions; until it
  exists, stars built from tables with duplicate-pressure rows remain unrepresentable (unchanged
  behaviour, now stated).
- ADR-0007 §7 item 11's bound is replaced by the tighter §Validation B target.
- `GOVERNANCE.md` §3.1 gains a third use (this ADR) after ADR-0002 and ADR-0007 — its "currently
  only use" sentence is stale documentation debt, not modified by any rank-below document.

## Validation (predeclared; measured before any baseline)

| # | Line | Bound / criterion |
|---|---|---|
| A | all Phase-4D analytic checks (`hartle_monopole_physics_analytic`) | green, values unchanged to `1e-12` (constant density: interior measure zero) |
| B | homogeneous DS(CMF)-1 sequence identity, res 10000 and 20000 | **`≤ 2e-4`** (scratch measure form `1.03e-4`, `5.7e-5`; oracle floor ≈ `7e-5`) |
| C | sourced DS(CMF)-1 source identity: integrated `m̂₀` vs the Stieltjes sum on the same partition | `≤ 1e-6` |
| D | radial convergence 5000/10000/20000/40000 of `δM̂` | successive differences monotone, relative spread `≤ 1e-4`; no node-placement behaviour |
| E | EOS sensitivity (Experiment I) re-run | the retired-FD substitution moves `δM̂` by `< 1e-3` (its 5 % was the missing measure) |
| F | HT68 published comparison | unchanged within `2e-5` |
| G | constant-density surface shell | byte-identical |
| H | detectors M1–M9 | still fire where applicable (M5 re-targeted to the measure) |
| I | new detector M10: omit the interior measure (use the nodal column) | must recreate the percent-level DS(CMF)-1 deficit (`≥ 3 %` on `δM̂`) |
| J | seven legacy artifacts, `I`, first-order tests | byte-identical |

Only after A–J: `tests/baselines/hartle_monopole_dscmf1_debug.tsv` with a same-build repeatability
tolerance, per `GOVERNANCE.md` §3.1 condition 7.

## Provenance

Drafted by the AI agent (Phase 4D-RG, 2026-09-03) from the primary source re-read
(H67 pp. 1009–1010, 1019–1022; C&M 1974 p. 69) and the scratch experiments recorded in
`docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md`. **The owner decides Q1–Q12.** No
production, test or CMake file changed; no baseline created.
