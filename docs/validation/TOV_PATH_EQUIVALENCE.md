# Phase 3E-0 — Equivalence measurement of the two live TOV paths

> **STATUS: `TOV-PATH EQUIVALENCE VERIFIED`**
>
> **No canonical TOV owner has been selected by this task.** This is measurement only. The
> `GOVERNANCE.md` fail-closed condition *"uncertain authoritative code path"* remains **OPEN**;
> resolving it requires the Phase-3E ownership ADR, which this task does not draft.

| Field | Value |
|---|---|
| **Starting HEAD** | `1ff63ee34b8c15aa78c55e43301cdaef69ecab7f` (10 ahead / 0 behind `master` `d2040d89…`, upstream equal, clean) |
| **Baseline** | **16/16** authenticated (193.50 s), **10/10** self-contained |
| **Change class** | test + documentation only; **no production source modified** |
| **Governing authority** | `GOVERNANCE.md`; ADR-0003, ADR-0004 (ACCEPTED); INV-03, INV-04, INV-13, INV-14 |

---

## 1. The two paths, as authenticated in source

**Path 1 — sequence / legacy-live** (`TOVSolver.cpp:1648`):

```
TOVSolver::Solve(Axis, dir, file)
  for idx = 0 .. Axis.res:                       # inclusive
    ec_req = Axis[idx]
    clamp ec to [central_eps_floor_factor*eps.front(), 0.999*eps.back()]
    init_press = p_of_e(ec)
    r = r_min ; y[1] = init_press ; y[0] = (4/3)*pi*r^3*GetEDens(y[1])
    RadiusLoop(r, y)                             # NStar::Append per step
    SurfaceIsReached()                           # -> NStar::FinalizeSurface()
    analysis ? analysis->Analyze(&n_star)
    n_exp_cond_f && n_exp_cond_f(n_star) ? ExportNStarProfile(idx, ...)
    Sequence::Add(n_star)
    n_star.Reset()
  ExportSequence(dir + file + "_Sequence.tsv")   # ALWAYS, unconditional
```

**Path 2 — validated profile** (`TOVSolver.cpp:2500`):

```
TOVSolver::SingleStarSolveToTOVPoints(ec, out_points)
  clamp ec  (same two expressions)
  init_press = p_of_e(ec)
  same y[] initialization
  radial loop -> out_points.emplace_back(TOVPoint)
NStar(points, species_labels) -> NStar::BuildFromTOV(...)
```

`SolveToProfile(target_M, …)` wraps Path 2 in a mass search. **It is not part of the equivalence
question** — Path 1 has no such orchestration, and comparing a sequence scan against a
root-finder would be a category error.

## 2. Source duplication map

| Item | Classification |
|---|---|
| `TOVSolver::ODE` | **SHARED FUNCTION** — one definition, both paths pass `{ODE, nullptr, 2, this}` |
| `p_of_e`, `GetEDens`, `GetNuDer`, `GetRho`, `GetRho_i`, `PressureCutoff` | **SHARED FUNCTION** |
| central-density clamp (`floor_e`, `ceil_e`) | **TEXTUALLY DUPLICATED BUT IDENTICAL** |
| initial radius / mass / pressure | **TEXTUALLY DUPLICATED BUT IDENTICAL** |
| GSL stepper `rk8pd`, initial step `1e-1`, abs `1e-10`, rel `1e-10` | **TEXTUALLY DUPLICATED BUT IDENTICAL** |
| `r_max`, `radial_res`, nominal step `(r_max−r_min)/radial_res` | **TEXTUALLY DUPLICATED BUT IDENTICAL** (Path 2 names the bounds `min_log_r`/`max_log_r`, but assigns them `r_min`/`r_max`) |
| step_scale ladder `100/1000/10000/100000 → .005/.025/.05/.25/1` | **TEXTUALLY DUPLICATED BUT IDENTICAL** |
| pressure cutoff evaluation | **TEXTUALLY DUPLICATED BUT IDENTICAL** — Path 1 calls `PressureCutoff()` per iteration, Path 2 hoists it to `p_cut`. `PressureCutoff()` reads only `init_press` and `eos_tab.pre[0]`, neither of which changes inside the loop, so hoisting is value-preserving. |
| append-then-test-cutoff **order** | **TEXTUALLY DUPLICATED BUT IDENTICAL** — both append the point, then break on `p ≤ cutoff`, so both include the first sub-cutoff node |
| `TOVPoint` construction | **TEXTUALLY DUPLICATED BUT IDENTICAL** — Path 1 builds it inline in a braced-init-list, Path 2 through named locals; same expressions, same inputs |
| **radial-loop body itself** | **DUPLICATED AND DIFFERENT (in code, not in value)** — `RadiusLoop` vs the inline loop in `SingleStarSolveToTOVPoints`; two texts, one algorithm |
| unit conversions (km, M☉, km⁻²) | **TEXTUALLY DUPLICATED BUT IDENTICAL** — `NStar::Append` and `BuildFromTOV` |
| `Lambda` construction | **TEXTUALLY DUPLICATED BUT IDENTICAL** expression; both retain/retained the legacy `1e-15` clamp shape, unreachable at `max 2m/r = 0.481` |
| `nu` reconstruction | **SHARED FUNCTION** — both call `EvaluateNu()` |
| species mapping | **SHARED SOURCE** — both take `eos_tab.extra_labels`; Path 1 via `NStar::Reset`, Path 2 via the `NStar(points, labels)` constructor |
| baryon-number integrand | **DUPLICATED AND DIFFERENT** — Path 2 uses the ADR-0004 primitive; Path 1's `FinalizeSurface` keeps the legacy inline `/(1−2m/r).sqrt()`. **Governed, intentional** (ADR-0004 §13). |
| `Find_MomInertia()` | **SHARED FUNCTION**, invoked at the same stage by both finalizers |
| `StarProfile` mirror surface scalars | **POSTPROCESSING DIFFERENCE** — `BuildFromTOV` calls `SetSurfaceScalars(M,R,z)`; `FinalizeSurface` **never does** (§7) |
| `seq.clear()` before fill | **POSTPROCESSING DIFFERENCE** — `BuildFromTOV` clears; `FinalizeSurface` does not (Path 1 relies on `NStar::Reset`) |
| `InitInterpolantsFromProfile_()` | **PATH-SPECIFIC ORCHESTRATION** — `FinalizeSurface` calls it; `BuildFromTOV` deliberately does not (`NStar.cpp:1199`). No measured numerical consequence (§4). |
| sequence accumulation, `Reset()`, TSV export | **PATH-SPECIFIC ORCHESTRATION** — Path 1 only |

## 3. Capture method

`tests/core/tov_path_equivalence_cmf.cpp` defines a **test-only** `CapturingTOVSolver : public
TOVSolver` that:

- sets the **existing** `n_exp_cond_f` callback to accept every star;
- overrides the **existing virtual** `ExportNStarProfile(idx, dir)` — which `Solve()` already
  calls after `SurfaceIsReached()` and before `Sequence::Add` / `n_star.Reset()`, the only
  window in which a finished Path-1 star exists;
- copies scalar and column **values** out immediately (never a pointer into `n_star`, which is
  reset; never a copy of `NStar`, which is non-copyable);
- does **not** call the base implementation, so no profile file is written;
- exposes `eos_tab.extra_labels` read-only, so both paths use one label ordering (ADR-0001).

It does **not** override `RadiusLoop`, `SurfaceIsReached` or `ODE`; does not touch tolerances,
the pressure cutoff or the step schedule; does not fill `NStar` directly; does not bypass
`Solve()`. **Nothing it does executes before or during the integration**, so it cannot perturb
the values it reads. No production API was added or changed. All filesystem side effects —
including `Solve()`'s unconditional `ExportSequence` — are directed to a scratch directory that
the test removes on exit.

**Axis semantics, authenticated not assumed.** `Zaki::Math::Axis::operator[](i)` returns
`min + i*(max−min)/res`, and `Solve()` iterates `idx = 0..res` **inclusive**. `res = 0` would
divide by zero, so a single central density is expressed as a degenerate range with `res = 1`;
the test asserts `Axis[0] == Axis[1] == ec` exactly before relying on it. This also yields two
consecutive identical stars, which doubles as a state-leakage probe.

## 4. Results — Experiment A, fixed central density

Anchors frozen from `tests/baselines/baryon_number_dscmf1_reference.tsv` (the achieved central
densities of the validated 1.0/1.4/1.6/2.0 M☉ reference stars). The **same frozen number feeds
both paths**; no mass search runs.

| Anchor | `ec` (g/cm³) | N₁ = N₂ | radial columns | ec/pc/M/R | I | \|ΔB\|/B |
|---|---|---|---|---|---|---|
| 1.0 M☉ | `4.545504e14` | 2629 | **25/25 bitwise** | bitwise | bitwise | **0 (bitwise)** |
| 1.4 M☉ | `6.164883e14` | 2646 | **25/25 bitwise** | bitwise | bitwise | **0 (bitwise)** |
| 1.6 M☉ | `7.312533e14` | 2635 | **25/25 bitwise** | bitwise | bitwise | **0 (bitwise)** |
| 2.0 M☉ | `1.298349e15` | 2527 | **25/25 bitwise** | bitwise | bitwise | **0 (bitwise)** |

The 25 columns are `r`, `m`, `nu'`, `p`, `eps`, `nB`, `nu`, `Lambda` and **all 17 species
columns**. `max_profile_ulp = 0` everywhere: not "agree to roundoff" — **identical bit patterns
at every node of every column**.

Classification of the (empty) difference set: no RADIAL INTEGRATION, UNIT CONVERSION, NU
RECONSTRUCTION, LAMBDA CONSTRUCTION or SPECIES MAPPING difference exists at these anchors.

## 5. Species semantics (ADR-0001)

17 species, and the **label list and its order are identical** on both paths at every anchor,
every resolution and every sequence node — necessarily so, since both read
`eos_tab.extra_labels`. Every species column is bitwise equal. Both paths therefore populate
species columns with the **same semantics**; no scientific-semantic discrepancy exists, and the
canonicalization claim is not blocked on this axis.

## 6. Surface termination

Both loops append the point and *then* break on `p ≤ PressureCutoff()`, so both include the
first sub-cutoff node. Measured at every anchor: identical node count, and bitwise-identical
last `r`, last `p` and last `m`. Example (1.6 M☉): `N = 2635`, last `r = 1.346832307596e+01` km,
last `p = 4.237e-14`. **No path includes or excludes a terminal node the other does not.**

## 7. Derived `SeqPoint` scalars, and one interface asymmetry

`SeqPoint`'s authenticated schema is exactly six fields: `ec`, `m`, `r`, `pc`, `b`, `I`. At
every anchor, `ec`, `m`, `r`, `pc` and `I` are **bit-identical** and `b` is within its
predeclared bound.

**`I` is bit-identical on both paths** (e.g. `8.699575833361e+01`, `1.356161305669e+02`,
`1.595871413252e+02`, `1.937231442198e+02` km³). Both finalizers invoke the same
`Find_MomInertia()` → `RotationSolver::FindNMomInertia()` at the same stage, on
bit-identical profile data, and obtain the same answer. Hartle was not modified or normalized;
**INV-07 remains unresolved** and no claim about physical Ω or J is made here.

**One genuine difference, predicted by the source audit before measurement and confirmed:**

| | Path 1 | Path 2 |
|---|---|---|
| `StarProfile::MassSurface()` | **0** | `1.476690e+00` km (1.0 M☉ anchor) |
| `StarProfile::RadiusSurface()` | **0** | `1.342632e+01` km |
| `StarProfile::ExpNuSurface()` | **0** | `8.831934e-01` |

`NStar::Reset` zeroes these via `SetSurfaceScalars(0,0,0)`; `BuildFromTOV` sets them; **
`FinalizeSurface` never calls `SetSurfaceScalars` at all**. The `SeqPoint` `m`/`r` are correct
on both paths, so this is **PROFILE CONSTRUCTION**, not radial integration — but any consumer
reading the profile's mirror scalars gets zero on Path 1. The test asserts this as the current
characterized state so a change is caught. **A canonical-path ADR must decide whether to
populate them on the legacy path or retire the mirror fields.**

## 8. Experiment B — the real multi-star sequence loop

Ten central densities from `4.5e14` to `1.30e15` g/cm³, solved in **one continuous `Solve()`
invocation**, each compared against an independent Path-2 solve at the same Axis value.

- **10/10 nodes equivalent.** All eight structural columns bitwise at every node.
- Worst `|ΔB|/B` across the sweep: **1.637e-16** (one ULP), inside the `1.0e-15` bound.
- `M` swept `9.867195e-01 → 2.000533e+00` M☉, strictly distinct — the loop really solved ten
  different stars.
- **No cross-star state leakage.** Two consecutive *identical* Axis nodes produce bit-identical
  stars at all four anchors, and every node in the ten-star sweep matches a solve performed in a
  freshly constructed solver. `n_star.Reset()`, `Sequence::Add` and EOS-accelerator reuse
  introduce no dependence of star *i+1* on star *i*.

## 9. Experiment C — radial-resolution cross-check

`Path1(h)` vs `Path2(h)` at the same `h`, 1.6 M☉ anchor. **Not** a convergence study.

| `radial_res` | N₁ = N₂ | radial columns | ec/pc/M/R/I | \|ΔB\|/B |
|---|---|---|---|---|
| 5000 | 1319 | **25/25 bitwise** | bitwise | 1.640e-16 |
| 10000 | 2635 | **25/25 bitwise** | bitwise | 0 (bitwise) |
| 20000 | 5268 | **25/25 bitwise** | bitwise | 1.640e-16 |

Equivalence is structural, not an accidental match at the default resolution.

## 10. Central-density clamp contract

`floor = 1.65880787e+09`, `ceil = 1.10576387e+16` g/cm³.

| Probe | Path 1 achieved `ec` | Path 2 achieved `ec` | `M` |
|---|---|---|---|
| below floor (`0.5×floor`) | `1.658807870000e+09` | `1.658807870000e+09` | `1.18357064e-03` both |
| above ceiling (`1.5×ceil`) | `1.105763869260e+16` | `1.105763869260e+16` | `1.68800975e+00` both |

Both paths clamp to the identical value and produce the identical star. These clamped
configurations are **excluded from the scientific equivalence baseline**; this is an
API-contract check only.

## 11. Baryon number — the known, governed 3D difference

Phase 3D migrated **only** Path 2's proper-volume integrand to the ADR-0004 primitive; Path 1's
`FinalizeSurface` still divides by `(1−2m/r).sqrt()` inline. Measured across all 17 comparisons:
**maximum `|ΔB|/B` = 1.640e-16**, one ULP, against the `1.0e-15` predeclared in ADR-0004 §7.1
*before* that implementation existed. Twelve of seventeen comparisons are **bitwise**.

This is an **intentional, governed legacy conformance gap** — not a TOV-integrator failure, and
not repaired here.

## 12. Caller and interface audit

The task brief named three Path-1 callers. **There are six live ones** (plus this test); the
count was an undercount, which matters for any migration plan.

| Caller | Numerical result consumed? | File side effect consumed? | Analysis callback? | Migration concern |
|---|---|---|---|---|
| `main/Test/spin_therm_evol_main.cpp:67` | no | **yes** — relies on `Solve()`'s built-in `ExportSequence` | no | sequence TSV must survive |
| `main/Test/tov_debug_main.cpp:165` | no — its `GetSequence()`/`Profile()` reads are **commented out** (`:153-160`) | **yes** | no | sequence TSV must survive |
| `main/Examples/sig_omega.cpp:64` | no | **yes** — calls `ExportSequence` again at `:70` | no | sequence TSV must survive |
| `main/Examples/Table_5-8_Glenn.cpp:50` | no | **yes** — `ExportSequence` at `:57`, then re-reads the TSV at `:63` and plots it | no | **strongest file dependency**: it consumes its own output as input |
| `main/Examples/coulatt.cpp:50` | no | `ExportSequence` call is **commented out** (`:53`); `Solve()`'s unconditional export still fires | no | white-dwarf EOS, not CMF |
| `main/Examples/polytrope.cpp:66` | no | `ExportSequence` commented out (`:69`); unconditional export still fires | no | polytrope EOS |

`main/Examples/rotating_ns.cpp:62` is `r_solver.Solve(...)` — a **`RotationSolver`**, not a TOV
path. `sig_omega_rho.cpp:73` and `sig_omega_rho_nstar.cpp:72` are commented out.

**No caller anywhere in `main/` attaches an `Analysis` or sets `n_exp_cond_f`** (verified by
grep across the whole directory). So today the per-star profile export and the analysis hook are
**dead on all production callers** — this test is the only user of the export hook.

**Every live caller consumes only the file side effect.** None reads in-memory sequence or
profile state after `Solve()`.

## 13. Sequence export interface

`Solve()` ends with an **unconditional** `ExportSequence(dir + file + "_Sequence.tsv")`.
`Sequence::Export` (`TOVSolver.cpp:110-127`) writes a `Zaki::File::VecSaver` 1-D table with the
fixed header

```
ec(g/cm^3)      M(Sun)          R(km)           pc(dyne/cm^2)   B               I(km^3)
```

— the six `SeqPoint` fields in declaration order, `%-14s` columns, tab-separated.

**Classification: PUBLIC / WORKFLOW INTERFACE.** This is not a debug artifact: it is the *only*
thing every live Path-1 caller consumes, and `Table_5-8_Glenn.cpp` re-reads its own emitted TSV
to produce a figure. A future canonical path must preserve this file, its name, its schema and
its column order, or explicitly replace every caller.

The exported text was **not** used as numerical truth anywhere in this measurement; all
comparisons used in-memory full-precision values.

## 14. Architecture hypothesis

**H2 — one TOV algorithm copied into two radial-loop implementations, with different
orchestration and output ownership.**

The evidence is decisive rather than suggestive:

- H1 (*two independent implementations*) is **refuted**. Independent implementations do not
  produce bit-identical values in 25 columns across 2500–5300 nodes, at three resolutions, at
  fourteen distinct central densities. They share `ODE`, `p_of_e`, `GetEDens`, `GetNuDer`,
  `GetRho`, `GetRho_i`, `PressureCutoff` and `EvaluateNu` outright, and the loops that are
  duplicated are duplicated *textually*, down to the step_scale ladder and the
  append-then-break order.
- H3 (*equivalent at the ODE layer but materially different in postprocessing*) is **too
  strong**. Postprocessing differs in *code* (`Append`+`FinalizeSurface` vs `BuildFromTOV`) but
  agrees in *value*: `ec`, `pc`, `M`, `R` and `I` are bit-identical, and `B` differs only by the
  one-ULP amount ADR-0004 predeclared. The only material postprocessing differences are the
  **unset mirror surface scalars** (§7) and the **sequence/export orchestration** (§13) — real,
  but confined to interface, not to the physics.

So the duplication is a maintenance and ownership problem, not a numerical one. **That is the
finding that should shape ADR-0005**: the canonical primitive can be designated without any
numerical migration risk on the radial solve, and the genuine work is in the *interface* —
the sequence TSV contract, the unset mirror scalars, the dead `Analysis`/export hooks, and
Path 1's unmigrated ADR-0004 conformance.

## 15. Detector proof

Both mutations were temporary, applied to a **single verified-unique** source site so they could
not move both paths, and reverted with byte-identity confirmed by SHA-256.

| ID | Mutation | Result |
|---|---|---|
| **D1** — radial integration | `RadiusLoop`'s GSL driver relative tolerance `1e-10 → 1e-8` (Path 1 only; `SingleStarSolveToTOVPoints` untouched) | **FIRES** — 37 assertions fail: 11 radial columns diverge, `max\|rel\| = 9.134e-07`, `max ULP = 7.5e9`; `M` diverges at `1.117e-12` |
| **D2** — Path-1 postprocessing | remove the proper-volume metric factor from `FinalizeSurface` (`BuildFromTOV` untouched) | **FIRES** — exactly **8** assertions fail, **all of them the `B` check**, at `\|ΔB\|/B` = 8.4 %–13.8 %; every radial column stays bitwise |

D2's selectivity is the point: a postprocessing-only fault produces a postprocessing-only
failure, which is what makes the §11 classification (RADIAL INTEGRATION vs BARYON INTEGRATION)
trustworthy rather than assumed.

Neither mutation is committed; `git status CompactStar/` is clean.

## 16. Protected artifacts

All six unchanged, byte-identical, and **not regenerated**:

| Artifact | SHA-256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |

New: `tests/baselines/tov_path_equivalence_dscmf1.tsv` (17 rows, all `EQUIVALENT`).

## 17. The exact boundary of what has and has not been validated

**Validated by this task.** On `DS(CMF)-1_with_crust`, at fourteen distinct central densities
spanning 0.99–2.00 M☉ and at `radial_res` 5000/10000/20000, the two live TOV paths produce
**bit-identical** radial structure (all 25 columns), **bit-identical** `ec`/`pc`/`M`/`R`/`I`,
identical species labels and ordering, identical surface termination, identical clamp behavior,
and `B` agreeing to at most one ULP — the governed ADR-0004 gap. Path 1's sequence loop shows no
cross-star state leakage.

**NOT validated, and explicitly out of scope.**

- **Only one EOS.** `coulatt` (white dwarf) and `polytrope` exercise Path 1 on EOS families this
  measurement never touched.
- **Only the nonrotating visible-sector path.** `Solve_Mixed` and the dark-sector paths are
  untouched.
- **Physical Ω and J.** `I` is the validated scale-free Hartle quantity only; INV-07 is
  unresolved and this task makes no claim about it.
- **Path 1's ADR-0004 conformance.** Its baryon integrand remains unmigrated by design.
- **The mirror surface scalars.** Path 1 leaves them zero; this is measured and asserted, not
  fixed.
- **The `Analysis` and per-star-profile-export hooks.** Dead on every production caller;
  exercised here only by the test's own capture. Their behavior under a canonical path is
  unexamined.
- **Ownership.** No canonical path is designated. The fail-closed condition stays open.
