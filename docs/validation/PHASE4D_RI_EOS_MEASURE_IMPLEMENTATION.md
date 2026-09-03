# Phase 4D-RI — implementing the measure-complete EOS energy-density source (ADR-0008)

> **FORMAL STATUS: `CORRECTION IMPLEMENTED — INDEPENDENT REVALIDATION REQUIRED — NO MONOPOLE BASELINE`**
>
> ADR-0008 was **ACCEPTED** (2026-09-03, Q1–Q12 adjudicated) and its correction implemented: the
> EOS energy-density contribution to `dm₀/dr` is now the measure `−4πr²ξ̂₀ dε` of Hartle's eq. (93),
> integrated one governed profile segment at a time, with the surface shell the terminal
> `ε_* → 0` atom of that same measure. The nodal `dε/dp` column is retained for the regular-centre
> series and diagnostics (ADR-0008 Q8) and is no longer the radial mass source.
>
> **Nothing here is a validated physical number.** This increment proves conformance to an
> accepted contract; the corrected **independent** revalidation is the next increment, and the
> first monopole baseline stays blocked until it succeeds (ADR-0008 Q12).

| Field | Value |
|---|---|
| **Starting HEAD** | `8abbab48df792b33c87bcf87d52f639ed5ff7fc6` (Phase 4D-RG), branch `physics/rotation-correctness`, upstream equal, clean, **10 ahead / 0 behind** `master` = `df859b5a73c4cac0c115f240744d89ce9f830b8d` |
| **Change class** | **scientific-semantic** (the source measure of H67 eq. 97) **and numerical-method** (its integration), executed under `GOVERNANCE.md` §3.1 — third use, authorized by ADR-0008 |
| **Governing authority** | ADR-0008 (ACCEPTED 2026-09-03) Q1–Q12; ADR-0007 P1, P3, P4, P7–P14 unchanged; ADR-0003 (provenance); ADR-0004 (metric factor); ADR-0005; ADR-0006; INV-06, INV-08, INV-09, INV-13 |
| **Pre-task suites** | full **30/30** (224.43 s), self-contained **16/16** (16.38 s); seven artifact hashes as `PHASE3_CLOSEOUT.md` §1 |
| **Post-task suites** | full **31/31** (240.55 s), self-contained **17/17** (20.81 s) — §19 |
| **Acceptance commit** | `cc4bec4` — `docs: accept measure-complete Hartle EOS source` (precedes every production edit) |

---

## 1. Starting state

```
worktree  /Users/keeper/Documents/CompactStar/worktrees/CompactStar-rotation
branch    physics/rotation-correctness   HEAD = 8abbab4   upstream equal   clean
master    df859b5                        10 ahead / 0 behind
```

Seven durable artifacts, recorded before any edit and unchanged afterwards (§18):

| Artifact | SHA-256 |
|---|---|
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `tov_path_equivalence_dscmf1.tsv` | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |

## 2. ADR-0008 acceptance

`PROPOSED → ACCEPTED — 2026-09-03`, in its own commit `cc4bec4`, **before** any production edit.
ADR-0007's accepted Decision was not rewritten; it carries only the post-validation note added in
Phase 4D-RG pointing at ADR-0008.

## 3. Q1–Q12 as adjudicated

| # | Decision |
|---|---|
| Q1 | **Option C** — the EOS energy-density contribution is the measure `dm̂₀\|_EOS = −4πr²ξ̂₀ dε`, for all energy-density variation; the `dε/dp` form is a smooth-region rewrite, not the numerical source. |
| Q2 | **Profile `ε` nodes are the mandatory source of truth**; an immutable EOS `(p_k, ε_k)` snapshot is optional/future and **was not needed** (§14 shows the profile partition meets the accepted contract); no invented transition metadata. |
| Q3 | Per-segment source `−4πr²ξ̂₀(r)·(ε_{i+1}−ε_i)/(r_{i+1}−r_i)`; **profile-node boundaries are mandatory integration boundaries**; no operator splitting; exact jump operator for declared discontinuities; **no steepness threshold**. |
| Q4 | Sharp continuous segments handled by the measure; no `Δε/ε > X` rule exists anywhere in the implementation. |
| Q5 | True discontinuities: no smoothing; `p̂₀*`, `ĥ₀`, `ξ̂₀` continuous; current EOS metadata cannot express them, and none was invented (§10). |
| Q6 | EOS layer owns `p_t, ε⁻, ε⁺`; profile owns `r_t`; **no transition detector exists**. |
| Q7 | One unified measure: interior over `[r₀, R_*)`, terminal `ε_* → 0` atom applied **exactly once**, still exposed as `surface_shell_mass_over_Omega2` (§9). |
| Q8 | `dε/dp` retained for the regular-centre series, cross-checks and diagnostics; removed as the radial mass source; its fail-closed semantics unchanged (§8, M6). |
| Q9 | Provenance unchanged: profile identity + `Version()`; no new cache (§11). |
| Q10 | Point-constructed constant-density stars unchanged; explicit `dε/dp` still required for the centre series (M0, M1). |
| Q11 | Phase-5 particle-number response must be measure-complete; **nothing of Phase 5 implemented**. |
| Q12 | `246f3f2` is not a reference; acceptance → correction → independent revalidation → first baseline. **No baseline here.** |

## 4. `GOVERNANCE.md` §3.1 — restated before implementation

1. **Existing physically incomplete behaviour** — the monopole mass source used the
   nodal/interpolated `dε/dp` representation and lost EOS energy-density variation narrower than
   the profile spacing (`PHASE4D_MONOPOLE_VALIDATION.md` §14).
2. **Why the current output cannot be baselined** — DS(CMF)-1 `δM̂` omitted ≈ 4.8 % of the governed
   EOS-measure contribution; freezing it would enshrine the defect.
3. **Minimum authorized correction** — ADR-0008 only: the EOS energy-density source of the `l = 0`
   system, plus unifying the already-existing terminal shell under the same measure semantics.
4. **Independent verification substituting for regression** — the ADR-0008 derivation and scratch
   evidence (`PHASE4D_R_EOS_MEASURE_DERIVATION.md`: measure identity from H67 eq. (93), the
   smooth-region equivalence proof, the exact two-layer jump certification at `2e-10…4e-9`, the
   convergence comparison) — evidence that stands with no prior run of this code. The **corrected
   independent revalidation is still required** and is not claimed here.
5. **Narrow scope** — ordinary `NStar`, `l = 0`, the EOS energy-density source. Nothing else.
6. **Superseded outputs** — every pre-correction monopole number, including all 4C-I1 and 4D
   DS(CMF)-1 diagnostics, is recorded as **not a reference result**.
7. **Baseline** — forbidden until the corrected independent revalidation succeeds.

**Status: `CORRECTION IMPLEMENTED — INDEPENDENT REVALIDATION REQUIRED — NO MONOPOLE BASELINE`.**

## 5. The retired source representation

`RotationSolver::ODE_HartleMonopole_`, term 1 as ADR-0007 P2 stated it:

```
term1 = 4 pi r^2 (eps + p) * dedp * p0*_hat ,      dedp = the profile column, linearly interpolated
```

It was the only radial EOS mass source, and it is gone as a general source. `b.dedp` is no longer
read anywhere in the right-hand side, so there are not two active EOS mass sources.

## 6. The new measure source

```
term1 = 0                                             if the segment's Delta eps is exactly zero
term1 = -4 pi r^2 * xi0_hat(r) * eps_slope            otherwise,     xi0_hat = p0*_hat / nu'
```

with `eps_slope = (eps_{i+1} - eps_i)/(r_{i+1} - r_i)` the measure density of the segment the
driver is inside. The segment's **total** contribution is therefore exactly `Δε_i` whatever the
EOS does between the two nodes, which is the property the nodal derivative lacked. The slope is
never reconstructed from `dedp * dp/dr` (that product is what fails inside a sub-node feature,
where the piecewise-linear background violates `dp/dr = −(ε+p)ν'` at `O(1)`).

`ξ̂₀ = p̂₀*/ν'` is formed directly; where that division is ill-conditioned the right-hand side
applies **the same regular-centre limit the derived `ξ₀` column already uses**,
`ξ̂₀ → j_c² s_c² r/[4π(ε+3p)]` (`MonopoleBackground_::centre_xi_num`), and **fails closed**
(`GSL_EBADFUNC`) if that is unavailable too. No epsilon regularization, no denominator floor.

## 7. Segment integration architecture

`ComputeMonopoleResponse()` first validates the partition — strictly increasing radii, finite
energy-density increments — and fails closed otherwise (a new fail-closed path, ADR-0008 Q3).
The driver then advances **one governed segment per `gsl_odeiv2_driver_apply` call**, with that
segment's `eps_slope` installed before the call:

```
for i in 0 .. n-1:
    eps_slope = (E[hi] - E[lo]) / (R[hi] - R[lo])      # hi = i, lo = i-1  (i = 0: segment 0)
    driver_apply(driver, &r, R[i], y)                  # integrates exactly to R[i], never past it
    record m0_hat[i], p0*_hat[i]
```

`gsl_odeiv2_driver_apply` integrates to its target and never overshoots, so adaptive substeps stay
inside the segment and no step can carry one segment's `Δε` into another. The loop structure
(one call per node) is unchanged from 4C-I1; what is new is that the source is now a per-segment
quantity, which makes those boundaries load-bearing rather than incidental. Ordinary background
fields (`p`, `ε`, `m`, `ν`, `ν'`, `s`, `s'`) are still sampled at the actual RHS radius through
the existing shared-bracket linear interpolation — INV-13 is untouched; only the EOS
energy-density **measure** changed representation.

## 8. The centre keeps the pointwise EOS derivative

ADR-0007 P4's regular-centre series is unchanged and still consumes
`(dε/dp)_c` from the 4C-I0 authority:

```
p̂₀*(r₀) = (1/3) j² s² r₀² ,     m̂₀(r₀) = (4π/15)(ε+p)[(dε/dp)(p_c) + 2] j² s² r₀⁵
```

It is a local, well-resolved property of the central state, not an integrated measure (ADR-0008
Q8). A star without an authoritative `dε/dp` still fails closed (M6).

## 9. The terminal atom

The interior measure runs over `[r₀, R_*)`; the remaining jump `ε_* → 0` is applied **once**, by
the same operator a declared internal discontinuity would use:

```
surface_shell_hat = 4 pi R_*^2 (eps_* - 0) xi0_hat(R_*)
deltaM_hat        = m0_hat(R_*^-) + surface_shell_hat + I^2 / R_*^3
```

It is still published as `surface_shell_mass_over_Omega2`, and it is never also folded into the
last interior segment, whose measure is `ε[n−1] − ε[n−2]` (M4c verifies the interior measure
telescopes to `ε_c − ε_*`, leaving `ε_*` to the atom). `I²/R_*³` and the `δM̂` formula are
otherwise unchanged.

## 10. Internal atoms — contract only

No governed EOS today can declare a constant-pressure jump (the importer requires strictly
increasing pressure; `gsl_interp_steffen` requires strictly increasing abscissae), so **no
internal jump operator runs and no transition detector exists**. The implementation is
nonetheless shaped so a future declared transition is inserted at a segment boundary — the source
is already per-segment and `m̂₀` already accumulates across boundaries — without touching the
scientific equation. No EOS-import redesign was attempted.

## 11. Provenance

The response depends only on the profile's radial fields, the first-order response, the central
`dε/dp` and `Version()`. No new data and no new cache authority were introduced, so
`(source_profile, source_version)` still guards the cached `HartleMonopoleResponse`; the
contract test's provenance and staleness checks (`C6`, `C7a`) pass unchanged.

## 12. Analytic constant-density star — unchanged, by construction

`tests/rotation/hartle_monopole_measure_contract.cpp` (new, self-contained):

| Check | Result |
|---|---|
| M0 point-constructed star with explicit `dε/dp = 0` still computes | pass |
| M1a every interior segment has `Δε = 0` | `max \|Δε_i\| = 0.0` exactly |
| M1b the EOS channel of the source integrates to zero | `0.0` exactly |
| M1c production `m̂₀` equals the rotation-only integral at every node | worst rel `0.0` (bound `1e-14`) |
| M4a shell = terminal atom `4πR_*²ε_*ξ̂₀(R_*)` | `1.31420887e+03` both |
| M4b `δM̂ = m̂₀(R_*) + atom + I²/R_*³` | rel `0.0` |
| M4c interior measure telescopes to `ε_c − ε_*`, excluding the surface drop | interior sum `0.0`, `ε_* = 2.173e-4` |
| M6 no authoritative `dε/dp` ⇒ fails closed | pass |

`δM̂ = 1.4674047059e+03 km³`, `shell/δM̂ = 0.89560` — **identical to the pre-correction value**, and
the whole 4D analytic harness reproduces its 4D numbers digit for digit (`Cb 9.724e-09`,
`Cc 8.956e-01`, `Dc` intercepts `1.065e-06` / `4.915e-07`, `Ja 3.046e-09`). The published
Chandrasekhar–Miller comparison is likewise unmoved (`K1c 7.290e-04`, `K1d 2.699e-04`).

## 13. Smooth-EOS equivalence

Self-contained HT68 Harrison–Wheeler fixture, three central densities:

| `ε_c` [g/cm³] | nodes | `δM̂` rel vs the SUPERSEDED differential oracle | node-wise peak-relative (reported) |
|---|---|---|---|
| 3.00e14 | 3682 | `5.544e-06` | `1.140e-05` |
| 1.00e15 | 2730 | `1.254e-05` | `1.965e-05` |
| 3.00e15 | 2164 | `1.148e-05` | `2.166e-05` |

Bound `2e-5` on `δM̂` — met on all three. **Chronology, recorded:** a first run of this file also
asserted the bound on the node-wise metric and measured the third column; the two source
representations differ at `O(h²)` in the weight even where the EOS is smooth, and that difference
grows with the density gradient. The bound was **not widened** — it is asserted on `δM̂`, the
quantity ADR-0008's evidence table measured (`1.2e-5` on this EOS) and set it from, and the
node-wise number is reported. At the published level the same equivalence shows as
`K2h 4.937e-03 → 4.940e-03` and `K2i 1.118e-02` unchanged against HT68 Table 5 (bound `2e-2`).

## 14. Same-partition source accounting

A four-variable segment-aware accounting integrator (own driver, `rtol 1e-13`, `atol 1e-16`)
carries `m̂₀`, `p̂₀*`, the EOS channel of `m̂₀`, and the per-segment weight `∫4πr²ξ̂₀ dr`:

| `ε_c` | per-segment identity `Δm̂₀\|_EOS,i = −slope_i · W_i` | accounting vs production `m̂₀` |
|---|---|---|
| 3.00e14 | worst residual `4.245e-10 km³` = `3.769e-13` of the total EOS integral | worst node `3.757e-07`, `m̂₀(R_*)` `3.757e-07` |
| 1.00e15 | `1.647e-10 km³` = `3.761e-13` | `2.111e-07` |
| 3.00e15 | `5.732e-11 km³` = `3.760e-13` | `1.407e-07` |

Bounds `1e-10` (segment identity, measured against the star's total EOS integral because a
per-segment *relative* metric is meaningless on segments whose `Δε ≈ 0`) and `1e-6`
(ADR-0008 Validation C) — both met. This is an **implementation identity**, not a physical
validation. The profile-partition contract is therefore met without the optional EOS-knot
snapshot, so none was added (ADR-0008 Q2).

## 15. DS(CMF)-1 correction diagnostic — not a baseline, not frozen

Four standard stars, production resolution, finite responses, no solver failure, no stale
provenance:

| M☉ | `δM̂` before | `δM̂` after | change | `ξ̂₀(R_*)` before → after | `m̂₀(R_*)` after |
|---|---|---|---|---|---|
| 1.0 | 934.11 | **995.26** | **+6.54 %** | 3364.7 → 3303.5 | 992.13 |
| 1.4 | 892.93 | **940.92** | **+5.37 %** | 2111.8 → 2087.0 | 933.52 |
| 1.6 | 825.97 | **865.87** | **+4.83 %** | 1617.1 → 1601.7 | 855.44 |
| 2.0 | 556.15 | **573.74** | **+3.16 %** | 704.52 → 700.65 | 555.47 |

The direction and size are what ADR-0008 predicted (its scratch prototype gave `+4.8 %` and
`δM̂ = 865.866` at 1.6 M☉ against production's `865.866`); `R_*`, `I²/R_*³` and `(dε/dp)_c` are
unchanged, because the background is untouched. The 4D CMF harness, with the accepted source on
**both** sides, now agrees at `≤ 3.3e-7` (second-order-isolated chain) and `≤ 5.5e-6` (fully
independent chain) against the `1e-4` bound; the **superseded** differential oracle is reported
beside it and disagrees by `6.1e-2 / 5.1e-2 / 4.6e-2 / 3.1e-2` — that gap is the sub-node
energy-density variation the correction recovers. The near-vacuum identity tightened to
`2.1–2.7e-7` (bound `1e-6`).

**These numbers are diagnostics of this run.** None is frozen, none is a reference, and none is
claimed as independently validated physics.

The line that failed Phase 4D now meets its bound. The homogeneous solution's `δM̂` against the
non-rotating sequence derivative `(dM/dp_c) δp_c` — ADR-0007 §7 item 11, measured at
`1.17e-3 / 1.02e-3 / 1.04e-3` in Phase 4D — is now **`1.03e-4 / 5.7e-5 / 7.2e-5`** at radial
resolution 10000/20000/40000, inside both the original `1e-3` and ADR-0008's tightened Validation-B
target of `2e-4`; the Stieltjes cross-check against the profile's own `ε` steps agrees with the
ODE to `1e-6` (it was `1.1e-3` before, and that gap **was** the defect). This is measured with the
accepted source on the oracle side as well, and it lives inside the very harness the corrected
revalidation must re-run in full — so it is recorded as a strong implementation-level indicator,
**not** as the independent revalidation.

Two further measured consequences, recorded:

- **EOS-derivative sensitivity (Experiment I)** — substituting the retired profile finite
  difference for the Steffen authority now moves `δM̂` by **exactly `0.0`**: `dε/dp` no longer
  enters the radial source at all, and its only remaining influence is the `O(r₀⁵) ≈ 1e-25`
  centre-series initial value. ADR-0008 Validation line E (`< 1e-3`) is met with room to spare.
- **Radial convergence (Experiment H)** — `δM̂` at res 5000/10000/20000 is
  `865.868157 / 865.866034 / 865.836078`: relative **spread `3.7e-5`** (ADR-0008 Validation D asks
  `≤ 1e-4`) against `1.6 %` before the correction. Its successive differences are **not**
  monotone (`2.1e-3` then `3.0e-2 km³`) and `R_*` itself moves with the grid
  (`13.463458 → 13.471018 km`), so the residual is the TOV background's own resolution dependence.
  **The monotonicity half of line D is therefore NOT met at 4D-RI**; separating the background's
  contribution is the corrected revalidation increment's task, and the CMF harness now records
  this rather than asserting a 4D-era criterion.

## 16. Detectors

Each mutation applied at one production site of `RotationSolver.cpp` in a separate build tree,
the focused targets rebuilt and run, the source restored and verified byte-identical.

| Detector | Mutation | Fires on | Signature |
|---|---|---|---|
| **D1** | revert the EOS source to the nodal `dε/dp` differential form | `hartle_monopole_measure_contract` M3b (3 stars); `hartle_monopole_physics_cmf` G (8 lines) | same-partition accounting `5.9e-6 … 1.2e-5` against `1e-6`; CMF `δM̂` `6.1e-2 / 5.1e-2 / 4.6e-2 / 3.1e-2` against `1e-4` |
| **D2** | omit the EOS measure source entirely | measure contract M2 + M3b (6 lines); CMF G and the near-vacuum identity (12 lines) | accounting off by `1.7e1 … 3.5e1`; CMF `δM̂` off by `0.86 … 0.91` |
| **D3** | double-count the terminal surface atom | measure contract M4b; `hartle_monopole_contract` C13a; CMF matching arithmetic (4 stars) | `δM̂` identity `8.956e-1` — exactly the constant-density star's shell fraction |
| **D4** | use the previous segment's `Δε/Δr` for every interval | measure contract M2 + M3b (6 lines); CMF G (8 lines) | accounting `2.5e-3 … 3.2e-3` against `1e-6`; CMF `δM̂` `2.4e-3 … 2.6e-3` |

All four build, fire, and were reverted byte-identically (SHA-256 `c4f2c2f3898312d1…` before and
after the sweep). Two behaviours worth recording: `hartle_monopole_contract` is **blind** to D1, D2
and D4 — its star is constant-density, where every source form agrees, which is precisely why the
new measure-contract test exists; and D4's CMF signature (`2.5e-3`) is two orders of magnitude
smaller than D1's, so the segment-boundary discipline is a finer property than the source form and
needs the dedicated same-partition accounting to detect.

## 17. `radial_res = 2500` — excluded

The isolated TOV-layer anomaly recorded in `PHASE4D_R_EOS_MEASURE_DERIVATION.md` §5
(`SetRadialRes(2500)` producing a physically different DS(CMF)-1 star) is **out of scope here and
was not investigated or repaired**. No resolution statement in this record uses 2500; the
resolutions used are 5000, 10000, 20000 and 40000, and every convergence statement is
**CONDITIONAL ON THE TOV BACKGROUND; THE 2500 ANOMALY IS UNAUDITED**.

## 18. Legacy protection

First-order ODE, `HartleFirstOrderResponse` and `I` are untouched by construction — the change is
confined to the `l = 0` EOS source term and its integration loop. `hartle_moment_inertia_cmf`,
both first-order tests, the TOV and thermal suites pass unchanged, and the seven durable
artifacts are byte-identical to §1, in particular
`hartle_I_dscmf1_debug.tsv = ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`.

## 19. Suites and no-baseline status

| Configuration | Before | After | New tests |
|---|---|---|---|
| full (`COMPACTSTAR_EOS_DATA_ROOT` set) | 30/30, 224.43 s | **31/31, 240.55 s** | `hartle_monopole_measure_contract` (3.27 s) |
| self-contained | 16/16, 16.38 s | **17/17, 20.81 s** | same |

`git diff --check` clean; no scratch file is tracked. (One earlier suite run, deliberately made
concurrent with the detector sweep, aborted `heat_capacity_v1`; it passes in isolation in 10.6 s
and passes in both clean sequential runs above — a contention flake of that run, not a finding.)

**No baseline was created.** `tests/baselines/hartle_monopole_dscmf1_debug.tsv` does not exist, no
monopole regression is registered, and per ADR-0008 Q12 none may be created until the corrected
independent revalidation succeeds.

## 20. The revalidation gate

The next increment must satisfy ADR-0008 §Validation in full — A analytic unchanged; B homogeneous
DS(CMF)-1 sequence identity `≤ 2e-4`; C same-partition source identity `≤ 1e-6` (met here as an
implementation identity, §14); D radial convergence, **including the monotonicity half not met
here** (§15); E EOS sensitivity (met, §15); F HT68 published unchanged (met, §13); G
constant-density shell unchanged (met, §12); H detectors M1–M9 still firing; I the new detector
M10 "omit the interior measure" reproducing the percent-level deficit (D2 here is its
implementation-level form, §16); J seven artifacts byte-identical (met, §18) — **with the
independent `(m₀,h₀)` and continuum oracles re-run in full**, which this task deliberately did not
repeat. Only then the first monopole baseline.
