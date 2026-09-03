# TOV-SURF-G — the ordinary-star surface event and integration-termination contract

> **FORMAL STATUS: `TOV SURFACE CONTRACT ADR PROPOSED — OWNER ADJUDICATION REQUIRED`**
>
> This record derives, from the TOV-RR-01 audit and from scratch measurements on the
> **byte-identical production physics**, what the ordinary-star surface `R_*` should mean, how it
> should be located, and what an integration that never reaches it may publish. It grounds
> `docs/adr/ADR-0009-tov-surface-event-and-termination.md` (**PROPOSED**, not accepted). No
> production, test, CMake or baseline file is changed.

| Field | Value |
|---|---|
| **Starting commit** | `b93be5f1df73eec0df1db1fc0ec39c524551e0b0` (TOV-RR-01 audit), on rotation source `e2fe0adf53d7975d3a5c57e84a8d3481ffd2ce41`; both in ancestry; master `df859b5a…` unchanged |
| **Branch / worktree** | `governance/tov-surface-contract` / `worktrees/CompactStar-tov-surface` |
| **Change class** | **documentation** (governance / numerical-physics derivation); scratch outside the tree |
| **Suites at entry** | full **31/31**, self-contained **17/17** (§16); seven artifact hashes as `PHASE3_CLOSEOUT.md` §1 |
| **Governing authority** | `GOVERNANCE.md` §2–§3; ADR-0005 (the primitive owns "the termination test"; §7.3 froze `SolveToProfile`'s fallback policy for Phase 3E only); ADR-0003 (provenance); ADR-0004 (`EvaluateNu` surface BC out of its scope); ADR-0007 P7 / ADR-0008 Q7 (surface = EOS-floor node, terminal atom); INV-06; INV-13; TOV-RR-01 |
| **Scratch** | `scratchpad/ts/ts_surface.cpp` — production `GetEDens`, production TOV right-hand side re-expressed **without** the fatal trial-state guard, scratch orchestration only; not tracked |

---

## 1. The inherited defect (binding evidence, not re-litigated)

TOV-RR-01 (`docs/validation/TOV_RADIAL_RES_2500_AUDIT.md`) established that
`TOVSolver::ODE` returns `GSL_EBADFUNC` whenever a **trial** Runge–Kutta state has
`p < PressureCutoff()`; that `gsl_odeiv2_driver_apply` treats a user-function error as fatal, not
as a rejected step; that consequently the accepted-state test `if (y[1] <= p_cut) break;` fired in
**0 of 87** audited solves; that every ordinary star therefore stops **one requested output target
short** of the cutoff, with `R_*` scattering 33 m (0.24 %) across resolutions; and that at the
DS(CMF)-1 crust–core near-discontinuity the same guard ends the integration 0.56 km early for a
scattered set of resolutions and, at the **production default 10000**, for ≈3 % of central densities
near 2 M☉ — including both `radial_res = 2500` rows of the durable `grid_convergence_cmf` artifacts.

Two earlier records had already seen the symptom without the diagnosis: `TOV_REFERENCE.md` §5
found production `R` always *smaller* than CompOSE and attributed the whole residual to the
table-floor convention; `GRID_CONVERGENCE.md` §7 attributed its anomalous order 0.70 to "the
driver aborts a step whose internal trial evaluation dips below the cutoff even though the recorded
grid point is still above it, so the truncation point depends on the requested step size" — and left
it unrepaired because the source was frozen. This record is the governance that both deferred.

## 2. What the surface physically is

The TOV solution `(m(r), p(r))` is a continuous trajectory; the EOS supplied to CompactStar ends
at a table floor `p_floor = eos_tab.pre[0]`, below which the star's matter is **not described** by
any governed input. `PressureCutoff() = max(1e-15 p_c, p_floor)` therefore marks the **edge of the
governed physical domain**: on DS(CMF)-1 it is the floor (`3.351885e25` dyn cm⁻², `n_B = 1e-7 fm⁻³`,
still inside the outer crust; `1e-15 p_c = 9.3e19` is far below it), on the HT68 Harrison–Wheeler
table it is the `1e-15 p_c` term (the table reaches 7.8 g cm⁻³). Neither is the `p = 0` vacuum
surface, and ADR-0007 P7 / ADR-0008 Q7 deliberately keep it that way (`R_*` is "the EOS-floor
surface", `ξ₀(R_*)` the displacement of the floor isobar, the shell the terminal `ε_* → 0` atom).

Nothing in this record changes that **value**. What it changes is the **locator**: today `R_*` is
the last radial output target the integrator happened to reach before a trial state strayed
below the domain edge — a number that depends on the output grid, on the inherited step size,
and on where a stiff feature sits inside a trial step. Physically, `R_*` should be **the unique
radius at which the accepted TOV solution reaches the domain edge, `p(R_*) = p_cut`** — a property
of the star, not of how it was sampled.

## 3. Trial states versus accepted states in adaptive Runge–Kutta

An embedded RK step evaluates the right-hand side at internal stages
`y + h Σ a_ij k_j`. Those stage values are **not** states of the solution; near a stiff feature they
can overshoot by orders of magnitude, and the error controller's whole purpose is to *reject* such
a step and shrink `h`. The audit's step-size trace showed exactly that sequence in the surviving
runs (`h` collapsing from ~1400 cm to 40–70 cm across the crust–core feature) and its interruption
in the dying ones (a trial stage at `p < p_cut`, `GSL_EBADFUNC`, state never advanced).

So three things must be kept apart:

| | meaning | correct response |
|---|---|---|
| **physical domain violation** | an *accepted* state outside the governed EOS domain | must not be published as an interior point |
| **integration event** | the accepted trajectory crosses `p = p_cut` | terminate here; this *is* `R_*` |
| **internal trial state** | a stage value of a step not yet accepted | evaluate a bounded continuation; let the controller decide |

The current guard conflates all three. A trial stage below `p_cut` is neither a domain violation
nor the event; making it fatal converts an ordinary step rejection into the end of the star.

## 4. Cutoff value versus cutoff locator

Two decisions that ADR-0009 keeps separate:

- **Threshold value** `p_cut = max(1e-15 p_c, p_floor)` — authenticated on both governed EOSs
  (§2); **unchanged** by this ADR. Whether `1e-15` is scientifically ideal, or whether the floor
  should be extrapolated toward vacuum, are separate questions (`TOV_REFERENCE.md` §5) and are
  **not** reopened.
- **Locator** — how the radius satisfying `p(R_*) = p_cut` is obtained (§6). This is what is
  defective and what the ADR governs.

## 5. Trial-state continuation below the cutoff

For the controller to *locate* the crossing it must be able to evaluate the right-hand side at
trial states slightly below `p_cut`. Three bounded continuations were evaluated (scratch,
production `GetEDens`):

| policy | EOS at trial `p < p_cut` | `p` in the TOV terms | result at a raised `p_cut = 10 p_floor` (so the policies differ) |
|---|---|---|---|
| **0** | `ε(max(p, p_cut))` | as is | `R_D = 13.438454857`, `M_D = 1.59997575293` |
| **1** | `ε(max(p, p_floor))` — the real table down to its floor, then clamp | as is | `R_D = 13.438454857`, `M_D = 1.59997575293` — **identical to 0 to all printed digits** |
| 2 | freeze `ε(p_cut)` and `p := p_cut` in every term | frozen | scratch run did not reach the event (loop exhausted); not diagnosed, **not recommended** — freezing `p` in the force term makes the right-hand side non-differentiable at the domain edge |

At the true DS(CMF)-1 cutoff (`p_cut = p_floor`) policies 0 and 1 coincide by construction, and
**policy 1 is exactly what `SafeSplineEval` already does** (`x_use = max(x, xmin)`, with a warning).
The accepted trajectory ends at the crossing, so the continuation never appears in output; the
measured sensitivity of `R_*` and `M` to the choice is **below `1e-10`**. Recommendation: the
narrowest continuation — evaluate EOS quantities at `max(p, p_floor)` (the existing clamp), use
`p` as is elsewhere, restrict it to trial states by construction (the event ends the accepted
trajectory), and **never publish** a sub-cutoff state as an interior point.

## 6. Locating the event — four algorithms measured

All on DS(CMF)-1, `ε_c = 7.312533e14` (1.6 M☉), production right-hand side without the fatal
guard, GSL rk8pd `1e-10/1e-10` as production, eight output partitions (the five production
`radial_res` values and three uniform grids of 1.75 / 7 / 28 m):

| N or grid | `r_k` (last node above) | `p_k/p_cut` | **A** first node ≤ cutoff | **B** linear in `p` | **B'** linear in `log p` | **C** bracketed re-integration | **D** pressure-coordinate landing | `M_D` |
|---|---|---|---|---|---|---|---|---|
| 2500 | 13.464228 | 2.10 | 13.492228 | 13.475283 | 13.464256 | **13.473358** | **13.473358** | 1.599975772 |
| 5000 | 13.463458 | 2.22 | 13.477458 | 13.474226 | 13.472369 | **13.473358** | **13.473358** | 1.599975772 |
| 10000 | 13.468323 | 1.53 | 13.475323 | 13.473583 | 13.473140 | **13.473358** | **13.473358** | 1.599975772 |
| 20000 | 13.471018 | 1.23 | 13.474518 | 13.473417 | 13.473295 | **13.473358** | **13.473358** | 1.599975772 |
| 40000 | 13.472681 | 1.06 | 13.474431 | 13.473368 | 13.473335 | **13.473358** | **13.473358** | 1.599975772 |
| uniform 1.75 m | 13.473260 | 1.01 | 13.475010 | 13.473359 | 13.473351 | **13.473358** | **13.473358** | 1.599975772 |
| uniform 7 m | 13.468010 | 1.57 | 13.475010 | 13.473569 | 13.473180 | **13.473358** | **13.473358** | 1.599975772 |
| uniform 28 m | 13.468010 | 1.57 | 13.496010 | 13.474170 | 13.468027 | **13.473358** | **13.473358** | 1.599975772 |
| **spread over 8 partitions** | | | **1.6e-3** | 1.4e-4 | 6.8e-4 | **5.6e-11** | **7.8e-11** | **1.1e-11** |

- **A** (today's semantics, made non-fatal) still quantizes `R_*` to the output grid: 21 m of
  scatter, and it always *overshoots* the crossing by up to one step (13.492 vs 13.473 at 2500).
- **B** (interpolation) is partition-dependent at `1e-4`; log-interpolation is worse, because
  `p(r)` near the table floor is not log-linear.
- **C** (bisection in `r` with re-integration from the last accepted state, fresh driver,
  tolerance 1 µm; 21–25 bisections) and **D** (integrate the terminal interval with `p` as the
  independent variable, `dr/dp = 1/f_p`, `dm/dp = f_m/f_p`, landing on `p_cut` exactly) agree with
  each other to **`< 1e-10`** at every partition and reproduce the same `R_*` and `M` on all eight.
  `p(R_C)/p_cut − 1 ≤ 7e-8`.

| method | accuracy | deterministic | EOS-domain safe | cost | depends on `radial_res` | changes artifacts | complexity |
|---|---|---|---|---|---|---|---|
| A | one output step | yes | yes | none | **yes (1.6e-3)** | no (today) | none |
| B | `O(h²)` interp. | yes | yes | none | yes (1.4e-4) | yes | trivial |
| **C** | integrator tolerance | yes (fixed tolerance) | yes (bracket stays ≥ floor except the final probe) | ~25 short integrations | **no (6e-11)** | yes | small |
| **D** | integrator tolerance | yes | yes (`p` monotone, never below `p_cut`) | one short integration | **no (8e-11)** | yes | small; needs `dp/dr < 0` (true inside the star) |

**Recommendation:** a deterministic bracketed event locator — **D** as canonical (single
integration, exact landing, no tolerance loop), **C** as the admissible equivalent/cross-check —
with `R_*` independent of the output partition to the solver's stated tolerance. A and B are
rejected: the partition would still define the star.

## 7. Is the output target partition scientific?

With the event contract in place the answer is measured, not argued: eight partitions from 1.75 m
to 28 m give `M` to `1.1e-11` and `R_*` to `7.8e-11` (§6), and the sequence derivative `dM/dp_c` to
`3e-6` (§12). **The partition is a sampling device.** Its only legitimate effects are on *sampled*
quantities — integrals over stored nodes (baryon number, `I`, thermal integrals) — which is the
proper subject of a grid-convergence study, not on `M`, `R_*`, `z_surf` or the central-density root.

## 8. Driver reset (audit option R2)

Under the corrected event contract, inherited versus per-segment-reset step size:

| N | inherited `h`: `R_D`, `M_D` | reset per segment: `R_D`, `M_D` | rel. difference |
|---|---|---|---|
| 2500 | 13.473358235, 1.59997577197 | 13.473358234, 1.59997577201 | `3.3e-11`, `2.3e-11` |
| 3000 | 13.473358235, 1.59997577203 | 13.473358235, 1.59997577198 | `4.0e-11`, `2.9e-11` |
| 10000 | 13.473358235, 1.59997577199 | 13.473358235, 1.59997577199 | `2.9e-11`, `8.2e-13` |

The reset that cured the audit's catastrophe is **unnecessary once the guard is gone** — it changes
nothing above `4e-11`. Recommendation: **C — not part of the physical contract**; the partition
must not be allowed to control the integrator's adaptive history *as a matter of contract*, and a
reset may be adopted only as an optional stabilization if a post-fix invariance test demonstrates a
need. Adopting it now would only have masked the guard.

## 9. Publication after an integration that did not reach the surface

Today `SingleStarSolveToTOVPoints` returns whatever points existed before a GSL failure, and
`SolveToProfile` (a) bisects on `tmp.back().m` of *any* returned profile, complete or truncated, and
(b) on failing to bracket the target mass, "falls back to the closest coarse sample" and publishes
it. That is how the `grid_convergence_cmf` `B_fixed_mass` row acquired an `ε_c` 0.22 % too high: the
root-finder silently paid for a missing crust with central density. Sequence derivatives inherit
the same silence (audit §14: a factor-4 error from one truncated neighbour).

Recommendation: **fail closed.** An integration that does not reach the governed surface event is
not a stellar solution and must publish **no** authoritative profile; the status must propagate
through `SingleStarSolveToTOVPoints`, `SolveToProfile` (both the bisection and the coarse
fallback), the sequence workflow (`Solve`), the sequence derivatives, and every test builder
(`BuildAtResolution`, `BuildAtDensity`, …). No root-finder may compensate for a missing crust by
changing `ε_c`.

## 10. Impact map — what moves when `R_*` sits at the crossing

Measured at `radial_res = 10000`, production last node → event (method D):

| star | production `R`, `M`, `p_last`, `ε_last` | corrected `R_*`, `M`, `ε(p_cut)` | `ΔR` | `ΔM/M` | `Δ(M/R)/(M/R)` | `Δz/z` |
|---|---|---|---|---|---|---|---|
| 1.0 M☉ | 13.426323, 1.000043795, 3.54e25, 1.73e8 | 13.427522, 1.000043795, 1.6588e8 | **+1.2 m** | `2.3e-10` | `8.9e-5` | `1.1e-4` |
| 1.4 M☉ | 13.545323, 1.400021843, 3.61e25, 1.75e8 | 13.546350, 1.400021843, 1.6588e8 | **+1.0 m** | `1.5e-10` | `7.6e-5` | `1.0e-4` |
| 1.6 M☉ | 13.468323, 1.599975771, 5.13e25, 2.28e8 | 13.473358, 1.599975772, 1.6588e8 | **+5.0 m** | `7.1e-10` | `3.7e-4` | `5.2e-4` |
| 2.0 M☉ | 12.712323, 2.000095227, 5.63e25, 2.44e8 | 12.715983, 2.000095227, 1.6588e8 | **+3.7 m** | `3.8e-10` | `2.9e-4` | `4.7e-4` |

and for the five truncated `radial_res = 10000` densities near 2 M☉ the correction is *not* small:
`R` +0.30…0.34 km, `M` +`6e-4` relative (e.g. `ε_c = 1.321825e15`: `12.348 → 12.684` km,
`2.004846 → 2.006117` M☉), i.e. the outer crust is restored.

Downstream consumers of the final node (no formula changed here; each will move as shown):

| consumer | site | moves by |
|---|---|---|
| `R`, `M` surface scalars | `NStar::BuildFromTOV` → `SetSurfaceScalars` | `ΔR` above; `M` `≤ 7e-10` |
| `ν(R) = ½ ln(1 − 2M/R)` outward boundary of the ν integration | `NStar::EvaluateNu` (ADR-0004 lists it out of its scope) | `Δz/z` = `1e-4…5e-4`, propagates to every `e^{ν}` |
| surface redshift `z_surf = e^{ν(R)}`, photon emitting area `4πR²e^{2ν(R)}` | `BuildFromTOV`, `PhotonCooling_Details.cpp:310` (INV-06) | `1e-4…1e-3` on luminosity normalization |
| compactness `M/R` | derived | `7.6e-5…3.7e-4` |
| baryon number `B = ∫₀^{R_*}` | `BuildFromTOV` `B_integrand.Integrate` | the omitted 1–5 m at `ρ ≈ 1.7e8` g cm⁻³: `≲ 1e-9` relative |
| moment of inertia `I` (first order) | `RotationSolver::FindNMomInertia` to the last node | same order as `B`'s tail, `≲ 1e-8`; the golden `hartle_I_dscmf1_debug.tsv` will **not** stay bitwise |
| Hartle `δM̂ = m̂₀(R_*) + shell + I²/R_*³` | `ComputeMonopoleResponse` | `I²/R_*³` by `3ΔR/R ≈ 1e-3`; the shell uses `ε_* = ε(p_cut)` instead of the last node's `ε_last` (`1.7–2.4e8 → 1.66e8`) |
| grid-convergence quantities | `grid_convergence_cmf` | the partition dependence of `R`, `e^{ν(R)}`, `L_ν` (order 0.70) largely **disappears**; what remains is the convergence of the sampled integrals |
| `TOV_REFERENCE` residual vs CompOSE | `TOV_REFERENCE.md` §4–§5 | production moves toward CompOSE by `ΔR` (1–5 m of the 27–47 m residual); the floor-convention part remains |

## 11. Baseline and artifact impact (nothing regenerated here)

| artifact | class | reason |
|---|---|---|
| `grid_convergence_cmf_1p6_debug.tsv` | **A — contaminated** | `radial_res = 2500` rows are a truncated star (audit §16) |
| `grid_convergence_cmf_1p6_trajectory.tsv` | **A — contaminated** | same rows |
| `tov_dscmf1_reference.tsv` | **B — valid star, surface values move** | carries `R` (+1–5 m), `M` (`≤ 7e-10`) |
| `passive_cooling_cmf_1p6_debug.tsv` | **B** | `T_inf` trajectory depends on `R`, `z_surf`, emitting area (`Δz/z ≈ 5e-4` at 1.6 M☉) |
| `hartle_I_dscmf1_debug.tsv` | **B (small)** | `R_km` column moves; `I` integral tail `≲ 1e-8`; **not** bitwise |
| `baryon_number_dscmf1_reference.tsv` | **B (small)** | `B` upper limit moves; `≲ 1e-9` |
| `tov_path_equivalence_dscmf1.tsv` | **D — unknown pending implementation** | records the equivalence of two paths that will move together; bitwise only if it stores differences |
| **C — byte-identical expected** | none | every artifact carries `R_*` or an integral to it |

**Migration procedure (predeclared):** old artifacts remain historical evidence; after the
correction passes §12 in full, each class-A/B/D artifact is regenerated **once**, in a dedicated
commit that records the old → new delta of every changed column against the impact map above; any
delta outside the map is a STOP. Class-A rows are replaced by the corrected star at the same
resolution (§13).

## 12. Post-fix validation (predeclared; plan only)

| # | line | predeclared bound |
|---|---|---|
| V1 | no ordinary star terminates via GSL failure on the sweep | 0 of all solves |
| V2 | every successful star reaches the governed event | status `SURFACE_REACHED` for all |
| V3 | `p(R_*)/p_cut − 1` | `≤ 1e-8` (locator; D lands exactly, C `≤ 7e-8` measured at 1 µm) |
| V4 | target-partition invariance of `M`, `R_*` over the five production `radial_res` values and three uniform grids | `M` `≤ 1e-9`, `R_*` `≤ 1e-8` (measured `1.1e-11`, `7.8e-11` — two decades of margin) |
| V5 | no scattered truncation around the crust–core transition | all 29 audit resolutions and the ≥150-density scan complete with `SURFACE_REACHED` |
| V6 | `dM/dp_c` stability | partition spread `≤ 1e-5` at fixed step (measured `3e-6`); step pair (5e-4, 1e-3) agree to `≤ 1e-5` (measured `5e-7`) |
| V7 | 1.0 / 1.4 / 1.6 / 2.0 M☉ canonical stars | `M` within `1e-9`, `R_*` within the impact map (+1.0/+1.2/+5.0/+3.7 m ± 0.5 m) |
| V8 | `z_surf`, `M/R` convergence with resolution | monotone, spread `≤ 1e-8` (they no longer depend on the partition) |
| V9 | artifact regeneration | only after V1–V8, per §11 |
| V10 | first-order / Hartle / thermal suites | reviewed; `I` and `δM̂` shifts within the impact map |
| V11 | detector: restore the fatal trial-state guard | must recreate the scattered truncation (`{2500, 2510, 2525, 3000}` at 1.6 M☉; the five 10000 densities) |
| V12 | detector: publish a GSL-truncated profile | must be caught by the fail-closed tests of every orchestrator |

**Sweep** (V1/V5): central densities `3.0e14 … 1.6e15 g cm⁻³`, ≥ 150 log-spaced, **plus** the five
failing densities exactly; resolutions `{2500, 5000, 10000, 20000, 40000}` and the three uniform
grids; both governed EOSs (DS(CMF)-1, HT68 HW).

**Tolerances, separated:** *scientific validation* bounds are the V-lines above (set from the
integrator tolerance `1e-10`, the measured partition spreads, and `TOV_REFERENCE` §9's existing
accuracy); *regression* tolerances for the regenerated artifacts are set **after** regeneration from
same-build repeatability (bitwise expected, `1e-12` fallback), never from the scientific bounds.

## 13. Grid-convergence artifact after the fix

The corrected solver makes `M`, `R_*`, `z_surf` partition-independent to `1e-11`, so the study's
subject becomes the convergence of *sampled* quantities (`B`, `I`, `C_*`, `L_ν`, `T_inf`). `2500`
becomes a valid coarse point of the same branch and should be **retained**; the sequence should be
**extended** to `{1250, 2500, 5000, 10000, 20000, 40000, 80000}` so that an order can be fitted
from at least four asymptotic pairs, and its order re-measured — it will no longer be 0.70. The
file shape is not a constraint.

## 14. Interaction with ADR-0008 / Phase 4D

Phase 4D, 4D-RG and 4D-RI used `radial_res ∈ {5000, 10000, 20000, 40000}` at the 1.6 M☉ `ε_c`,
all normal (audit §16), so their **measure-physics conclusions stand and ADR-0008 is not
reopened**. The corrected independent revalidation and the first monopole baseline must run on the
corrected background: `R_*` moves +5 m at 1.6 M☉ (`I²/R_*³` by `≈ 1e-3`), the shell's `ε_*` becomes
`ε(p_cut)`, and the residual non-monotone `δM̂` convergence recorded in
`PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md` §15 is this surface artifact and should vanish. **Recorded
dependency: ADR-0009 implementation + validation precede the Phase-4D corrected revalidation and
the first monopole baseline.**

## 15. Change class and governance path

The repair is **scientific-semantic** (the definition of `R_*`, hence INV-06 and every surface
scalar) **and numerical-method** (event location, right-hand-side domain contract) **and
structural** (completion status through every orchestrator, fail-closed publication). The strictest
class governs: scientific-semantic, ADR required — **this ADR**.

`GOVERNANCE.md` §3.1 does **not** apply as the licence: durable TOV baselines exist
(`tov_dscmf1_reference.tsv` and six more), so this is not pre-baseline work; the change proceeds
under the normal scientific-semantic path with the §11 migration procedure. The single §3.1-like
element — that the `grid_convergence_cmf` `2500` rows are adjudicated invalid and would enshrine a
defect if kept as reference — is handled by that same migration, not by an exception.

## 16. Owner questions and non-scope

The questions are ADR-0009 Q1–Q14; the recommendations are its §"Recommendations". Non-scope:
no production or test change; no repair; no baseline regeneration; no reopening of the `p_cut`
value, the floor-to-vacuum extrapolation, ADR-0008's physics, or Phase 5; no merge. Suites and
hashes are re-authenticated at exit (§17 of the final report).
