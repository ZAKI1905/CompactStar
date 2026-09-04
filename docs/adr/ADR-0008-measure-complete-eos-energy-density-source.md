# ADR-0008 — Measure-complete EOS energy-density source for the Hartle monopole response

| Field | Value |
|---|---|
| **Status** | **ACCEPTED** — 2026-09-03 (owner adjudication of Q1–Q12 recorded in the Decision below) |
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

**ACCEPTED 2026-09-03.** The owner adjudicated every decision question. The bindings below are
normative; where they differ from ADR-0007 they amend it, and only as scoped in the banner.

| # | Binding decision |
|---|---|
| **Q1** | **Option C.** The EOS energy-density contribution to `dm̂₀` is the measure `dm̂₀\|_EOS = −4πr²ξ̂₀ dε`, for **all** EOS energy-density variation. The `dε/dp` differential form is a smooth-region **rewrite**, not the canonical numerical source representation. |
| **Q2** | **Mandatory source of truth: the `StarProfile` energy-density values at the radial nodes** (which already equal the EOS interpolant evaluated on the star). An immutable EOS `(p_k, ε_k)` snapshot for sub-segment/knot refinement is **optional and future**, and is **not** required by this correction unless the accepted profile-partition contract cannot otherwise be met. True-transition metadata is governed here but does not exist in the current EOS layer and must not be invented. |
| **Q3** | **Per-segment ODE source.** On profile interval `[r_i, r_{i+1}]`: `dε/dr\|_seg = (ε_{i+1} − ε_i)/(r_{i+1} − r_i)` and `source_EOS(r) = −4πr² ξ̂₀(r) · dε/dr\|_seg` throughout that segment. **Profile-node boundaries are mandatory integration boundaries**; the integrator may never carry one segment's measure density into another. No operator splitting of the continuous part. Declared true discontinuities (when the EOS layer supports them) use the exact jump operator `Δm̂₀ = 4πr_t²(ε⁻ − ε⁺) ξ̂₀(r_t)`. **No steepness threshold.** |
| **Q4** | Sharp continuous segments are handled by the measure regardless of steepness. A rule of the form `if Δε/ε > X → transition` is **forbidden**. |
| **Q5** | True discontinuities: no smoothing; exact jump operator; `p̂₀*`, `ĥ₀`, `ξ̂₀` remain continuous through the material interface. The current EOS representation cannot express constant-pressure jumps; this ADR does not broaden into an importer redesign. |
| **Q6** | Transition location: the EOS layer owns `p_t, ε⁻, ε⁺`; the profile owns `r_t` through the governed `p(r)` mapping. **No transition detector exists or is created now.** |
| **Q7** | **One unified measure.** The interior measure covers `[r₀, R_*)`; the terminal `ε_* → 0` atom at `R_*` is applied **exactly once** with the same measure/jump semantics and continues to be exposed as `surface_shell_mass_over_Omega2`. It must not be double counted inside the final interior segment. |
| **Q8** | 4C-I0 stands. `dε/dp` remains authoritative for the regular-centre series, smooth-region cross-checks, EOS diagnostics and future pointwise-derivative physics. It is **no longer the monopole mass-source integrator input away from the centre**. The profile `dε/dp` authority is **not** removed and its fail-closed semantics are unchanged. |
| **Q9** | Provenance unchanged: `StarProfile` identity + `Version()` remains sufficient. Any future measure metadata must be profile-attached and `Touch()` the profile when changed. **No new cache authority.** |
| **Q10** | Point-constructed stars: the profile `ε` nodes define the measure. On a constant-density analytic star the interior `Δε = 0` and the terminal surface atom carries the whole density drop; an explicit `dε/dp` is still required for the centre series. |
| **Q11** | Phase-5 principle: the particle-number structural response **must** be measure-complete (`dn_i` with atoms at composition/species-density discontinuities). A nodal `dn_i/dp` column alone is **forbidden** as the Phase-5 integrator source. Nothing of Phase 5 is implemented under this ADR. |
| **Q12** | The `246f3f2` monopole output is **not** a reference. Required chronology: **ADR-0008 acceptance → production correction → independent revalidation → first monopole baseline.** No baseline may be created before that revalidation succeeds. |

The correction this Decision authorizes proceeds under `GOVERNANCE.md` §3.1 (third use, after
ADR-0002 and ADR-0007): the behaviour being replaced is adjudicated here as physically incomplete,
capturing it as the golden baseline would enshrine it, the minimum correction is named above, the
independent evidence that substitutes for regression is
`docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md` (measure derivation from H67 eq. (93),
smooth-region equivalence, exact-jump certification on a two-layer star, convergence comparison),
the scope is the ordinary-`NStar` `l = 0` EOS energy-density source only, all pre-correction
monopole numbers are recorded as non-references, and the baseline follows immediately after the
independent revalidation of §Validation.

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

## Implementation record — Phase 4D-RI (2026-09-03)

The correction is **implemented**; evidence: `docs/validation/PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md`.
Acceptance (`cc4bec4`) preceded every production edit.

| Clause | As implemented |
|---|---|
| Q1/Q3 source | `RotationSolver::ODE_HartleMonopole_` term 1 is `-4 pi r^2 xi0_hat * eps_slope` with `eps_slope = (eps_{i+1}-eps_i)/(r_{i+1}-r_i)` of the segment the driver is inside; `dedp` is no longer read by the right-hand side, so there is exactly one active EOS mass source |
| Q3 boundaries | `ComputeMonopoleResponse()` advances one governed segment per `gsl_odeiv2_driver_apply`, installing that segment's measure density first; the partition is validated (strictly increasing `r`, finite `Delta eps`) and fails closed otherwise |
| Q4 | no threshold of any kind exists in the implementation |
| Q5/Q6 | no internal atom runs and no detector exists — the current EOS layer declares none; the per-segment structure admits one at a boundary without touching the equation |
| Q7 | interior measure over `[r0, R_*)`; the terminal `eps_* -> 0` atom applied once, still published as `surface_shell_mass_over_Omega2`; verified to telescope correctly (`M4c`) |
| Q8 | the regular-centre series still consumes `(deps/dp)_c`; absence still fails closed (`M6`) |
| Q9 | provenance untouched; `(source_profile, source_version)` still guards the cache |
| Q10 | point-constructed constant-density stars: interior measure exactly `0`, output **bitwise unchanged** |
| Q2 | the optional EOS-knot snapshot was **not** added: the profile partition meets the accepted contract (same-partition accounting `1.4e-7 ... 3.8e-7`) |

Measured: constant-density star unchanged (`deltaM_hat = 1.4674047059e+03`, EOS channel exactly
`0`); smooth HT68 Harrison-Wheeler EOS agrees with the superseded differential form on `deltaM_hat`
to `5.5e-6 / 1.25e-5 / 1.15e-5` (bound `2e-5`); per-segment measure identity `~4e-13` of the total
EOS integral; same-partition accounting vs production `<= 3.8e-7` (Validation C, `1e-6`); DS(CMF)-1
`deltaM_hat` moves `+6.5 / +5.4 / +4.8 / +3.2 %` at 1.0/1.4/1.6/2.0 Msun, in the predicted
direction and size; EOS-derivative sensitivity now exactly `0.0` (Validation E); radial spread of
`deltaM_hat` `3.7e-5` (Validation D's spread half, `1e-4`) though **its monotonicity half is NOT
met** — the residual is the TOV background's own resolution dependence and is deferred to the
revalidation increment. Detectors D1-D4 all fire and were reverted byte-identically. Seven durable
artifacts and `I` unchanged.

> **Status: `CORRECTION IMPLEMENTED — INDEPENDENT REVALIDATION REQUIRED — NO MONOPOLE BASELINE`.**
> The independent `(m0,h0)` and continuum revalidation of §Validation was deliberately not
> repeated here. No baseline exists or may be created before it succeeds.

---


## Revalidation record — Phase 4D-RV (2026-09-03)

The corrected independent revalidation of §Validation was executed on the migrated ADR-0009
backgrounds; evidence: `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md`. The `(m₀,h₀)` oracle now carries the
EOS measure by its own Stieltjes route (midpoint atoms on an independently refined partition,
optionally at the EOS-table knots), never production's segment density.

| Line | Outcome |
|---|---|
| A | met — analytic values unchanged (`4.7e-15` between the Stieltjes and differential oracles on the zero-measure interior) |
| B | met — `1.03e-4` (10000), `5.7e-5` (20000) with the independent route; `1.17e-3` for the superseded form |
| C | met — `≤ 8.2e-8` |
| D | **spread met (`3.7e-5 ≤ 1e-4`); monotonicity NOT MET** — with `R_*` fixed to `3e-11`, `δM̂` over 5000/10000/20000/40000 moves `−2.1e-3, −3.0e-2, +9.5e-3 km³`; the EOS-knot oracle makes the values monotone (`865.896 → 865.856 → 865.854 → 865.848 → 865.845`), and the first-order `I` carries the same `−9.5e-3 km³` dip at 20000 — the residual is O(h) node placement of the profile-partition weight (Q4) and of the validated sampled first-order background, not monopole physics; the 4D-RI attribution to a moving `R_*` is superseded |
| E | met — `0.0` |
| F | met at its authenticated scope (`δM̂` measure vs differential at `3e14/1e15/3e15`: `1.25e-5`); on the denser HT68 models the superseded form's deficit grows to `1.7e-4` while the independent oracles agree with production to `≤ 1.9e-6` (recorded) |
| G | met — bitwise |
| H, I | M1–M9 fire; M10 recreates the ≥ 3 % deficit — see the record §13 |
| J | met — seven artifacts and `I` unchanged |

> **Status: `CORdocs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.mdTION IMPLEMENTED — INDEPENDENT REVALIDATION EXECUTED — VALIDATION D MONOTONICITY NOT MET — NO MONOPOLE BASELINE`.**
> Scientific status `HARTLE O(OMEGA^2) MONOPOLE RESPONSE CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE`. The accepted Decision is
> unchanged. The owner's adjudication of D's monotonicity clause against the measured node-placement
> floor (record §20) precedes any baseline.


## Post-validation owner clarification — Validation D monotonicity (Phase 4D-DA, 2026-09-04)

This clarification leaves the accepted Decision Q1–Q12 unchanged and preserves the original wording
of Validation line D above. Evidence: `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md` §21.

- **What D said.** "successive differences monotone, relative spread `≤ 1e-4`; no node-placement
  behaviour". The monotonicity phrase was never operationally defined, and the evidence this ADR cited
  for the corrected form (`PHASE4D_R_EOS_MEASURE_DERIVATION.md:224`, `865.868 / 865.866 / 865.836 /
  865.846`) already had non-monotone successive differences while being accepted as showing "no
  node-placement behaviour (`4e-5`)". Its evidenced purpose was to exclude the percent-level pathology
  of the nodal-derivative column.
- **Why it cannot be required of the accepted representation.** Q3/Q4 accept a profile-partition
  measure whose sub-segment weight location carries an `O(h)` error. For a sharp feature at phase `φ`
  within a cell that error is `Δε·W′·h·(1/2 − φ)`, with `φ → frac(2φ)` under grid doubling; a correct
  implementation therefore fails strict monotonicity in most phase orbits (toy model: values monotone
  25 %, difference magnitudes decreasing 55 % with a sampled background present), while the envelope
  `|e| ≤ Δε|W′|h/2` always holds. On DS(CMF)-1 the predicted sign matches production minus the
  EOS-knot oracle at 6/6 resolutions and the residual is 3–13 % of the computed envelope, which halves
  per doubling. The validated first-order sampled background (`I`) carries the same phase term.
  Literal "no node-placement behaviour" would contradict Q4; it is read as the absence of the material
  pathology, with the remaining representation term required to be bounded and vanishing.
- **Validation D as clarified (D′).** D′1 spread `≤ 1e-4` (unchanged); D′2 production vs the
  independent same-representation Stieltjes oracle `≤ 1e-6` at every ladder resolution (Validation C's
  bound); D′3 production vs the independent EOS-knot Stieltjes oracle `≤ 1e-4` (ADR-0007 §7-4); D′4
  `|production − knot oracle|_N ≤ E_crust(N) = Σ_crust |Δε_i||W′_i|h_i/2` from the production profile,
  with `E_crust ∝ 1/N`; D′5 `R_*` spread `≤ 1e-8`, M10 `≥ 3 %`, sequence identity `≤ 2e-4`. Measured:
  `3.71e-5`; `≤ 5.8e-7`; `≤ 3.2e-5`; ratios `≤ 0.131` with the envelope halving exactly; `3.4e-11`,
  `4.61 %`, `1.03e-4`. Continuum `δM̂(1.6 M☉) = 865.845 ± 0.001 km³` (diagnostic).
- **EOS-knot partition (Q2b)** remains optional: it removes a `≤ 3e-5` term below the scientific bound
  and cannot restore monotone differences while the first-order background is sampled.

> **Status: `CORdocs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.mdTION IMPLEMENTED — INDEPENDENTLY VERIFIED (Phase 4D-RV + 4D-DA) — FIRST MONOPOLE BASELINE NEXT`.**
> Scientific status `HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED`. No baseline was created by the adjudication; condition 7 of `GOVERNANCE.md`
> §3.1 is discharged by the separate first-baseline task.

## Provenance

Drafted by the AI agent (Phase 4D-RG, 2026-09-03) from the primary source re-read
(H67 pp. 1009–1010, 1019–1022; C&M 1974 p. 69) and the scratch experiments recorded in
`docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md`. **The owner decides Q1–Q12.** No
production, test or CMake file changed; no baseline created.
