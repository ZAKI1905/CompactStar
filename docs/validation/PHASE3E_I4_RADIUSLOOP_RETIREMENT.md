# Phase 3E-I4 — `GenTestSequence` coverage and `RadiusLoop` retirement

> **STATUS: `PHASE-3E-I4 DUPLICATE TOV RADIAL LOOP RETIRED`**
> and, on the evidence below, **`PHASE-3E CANONICAL TOV OWNERSHIP COMPLETE`**.
>
> **`GOVERNANCE.md` fail-closed condition #3 is CLOSED for ordinary visible-sector TOV radial
> numerical ownership.** Scope is stated precisely in §14 — it does **not** cover `MixedStar`,
> the two-fluid path, or dark-sector TOV.
>
> This document is written in two stages, matching the two commits. Commit A establishes
> coverage against the **legacy** implementation and proves the coverage can detect it. Commit B
> performs the migration and retirement. The ordering is the point: the test must exist against
> the old code before it is used as evidence for the new code.

| Field | Value |
|---|---|
| **Starting HEAD** | `d4f62da42750c752dea5421c3fea373003516459` (15 ahead / 0 behind `master`, upstream equal, clean) |
| **Baseline** | 18/18 authenticated (195.22 s), 10/10 self-contained |
| **Governing authority** | **ADR-0005 (ACCEPTED 2026-09-02)** §8, §15; ADR-0004; `GOVERNANCE.md` §3 condition 3 |

---

## 1. Why `RadiusLoop` survived I1

ADR-0005 §8 asserted that `Solve()` was `RadiusLoop`'s only live caller. **That was wrong**, and
Phase 3E-I1 found it: `TOVSolver::GenTestSequence` also calls it. Because that routine is
**public, compiled, and unexercised** — its sole repository reference is commented out at
`main/Test/tov_debug_main.cpp:196-203` — deleting `RadiusLoop` in I1 would have meant migrating
uncovered code inside a structural increment, which `AGENTS.md` forbids and which ADR-0004 §0-Q2
had already declined to do for `MixedStar`. I1 therefore retained it, marked it
non-authoritative, and deferred retirement to I4.

**I4 removes the blocker in the correct order: coverage first, migration second.**

## 2. `GenTestSequence` — authenticated behavior (pre-migration)

Call graph, as it exists at `d4f62da`:

```
TOVSolver::GenTestSequence(ec, dir, file)
  print TOV_gsl_interp_type
  radial_res_test test          # 16 entries
  p_of_e_prec = 1e-9            # SET, never restored
  for idx = 0 .. 15:
      test.Modify(this, idx)    # -> SetRadialRes(radial_res_set[idx]); never restored
      PrintStatus(idx, 16)
      init_press = p_of_e(ec)   # NOTE: no floor/ceiling clamp
      r = r_min ; y[1] = init_press ; y[0] = (4/3) pi r^3 GetEDens(y[1])
      RadiusLoop(r, y)          # <- the duplicate ordinary-star radial implementation
      SurfaceIsReached()        # -> NStar::FinalizeSurface()
      analysis ? analysis->Analyze(&n_star)
      n_exp_cond_f ? ( n_exp_cond_f(n_star) ? ExportNStarProfile(idx, dir+"/profiles"+file) )
      Sequence::Add(n_star)
      n_star.Reset()
  ExportSequence(dir + file + "_TestSequence.tsv")
  analysis ? analysis->Export(wrk_dir_ + dir)
```

**The authenticated 16-resolution list** (`radial_res_test::radial_res_set`, private to
`TOVSolver.cpp`), in order:

```
10000, 15000, 20000, 25000, 30000, 40000, 50000, 55000,
60000, 65000, 70000, 75000, 80000, 85000, 90000, 100000
```

**Persistent side effects, deliberately preserved and not "improved" in I4:**

- `p_of_e_prec` is set to `1e-9` and **never restored**;
- `radial_res` is left at `100000` after the call;
- `sequence` is **not** cleared on entry — repeated calls, or a call after `Solve()`, accumulate
  rows. `ClearSequence()` exists but `GenTestSequence` does not call it.

**Output.** `<file>_TestSequence.tsv`, written through the same `Sequence::Export` as the
ordinary workflow, so the schema is the identical six fields in the identical order:

```
ec(g/cm^3)    	 M(Sun)        	 R(km)         	 pc(dyne/cm^2) 	 B             	 I(km^3)       
```

**The file carries no radial-resolution column.** Row *i* is only mappable to a resolution
through the authenticated ordering above, which is why the coverage test duplicates that list
explicitly as an interface assertion.

## 3. Frozen central density

`ec = 731253342677476.12` g/cm³ — the achieved central density of the 1.6 M☉ DS(CMF)-1 reference
star, taken from the committed `tests/baselines/baryon_number_dscmf1_reference.tsv`. Measured EOS
bounds are `floor = 1.65880787e+09` and `ceil = 1.10576387e+16` g/cm³, so it is comfortably
in-domain and the clamp difference of §7 is a no-op for it.

No mass search runs inside the test; the value is frozen from existing authoritative evidence.

## 4. Commit A — the coverage test

`tests/core/tov_gen_test_sequence_cmf.cpp` (external-data). For the frozen `ec`, at all sixteen
production resolutions, it compares

> the **real** `GenTestSequence` (whatever radial engine it currently uses)
> **vs** `SingleStarSolveToTOVPoints` + the **same** `Append`/`FinalizeSurface` postprocessing.

Holding postprocessing fixed on both sides isolates the radial implementation — the only thing
I4 changes. The reference is built from real production APIs; **no TOV equations are
reimplemented and no independent solver exists in the test**. The reference also sets
`p_of_e_prec = 1e-9` so it matches the legacy routine's own precision.

In-memory values are captured through the **existing** virtual `ExportNStarProfile` hook enabled
via the **existing** `n_exp_cond_f` callback — the observational mechanism established in 3E-0.
No production API was added.

**Result against the legacy `RadiusLoop` implementation:**

| Assertion | Result |
|---|---|
| G1 one star per authenticated resolution | 16/16 |
| G2 `<file>_TestSequence.tsv` at the requested location | pass |
| G3a/b six fields, exact names and order; raw header byte-exact | pass |
| G4 exported row count = 16 | pass |
| G5 six numeric fields per row | pass |
| **G6 every resolution BIT-IDENTICAL to the canonical primitive** | **16/16 bitwise** on `(ec, M, R, pc, B, I)` |
| G7 the sweep varies structure | `R` 13.46832307596 → 13.47335607524 km |
| G8 exported rows match in-memory `SeqPoint` at export precision | worst rel **3.487e-09** |

> **`GEN_TEST_SEQUENCE EQUIVALENT TO CANONICAL PRIMITIVE`** — established **before** any
> production migration.

**Runtime: 4 s** (§7 pilot), well inside the ≤30 s guidance. The 16 production resolutions were
**not** reduced to make the test cheap.

**Pre-migration reference state**, captured for the post-migration comparison:

| | |
|---|---|
| `gts_TestSequence.tsv` SHA-256 | `e43bd05c1b22f7dbc26a39a10fe68d013bde5440eb3d5ce502eb83555245e511` |
| full-precision in-memory dump | 16 rows, `%.17g` per field (`--dump`) |
| res 10000 row | `M = 1.5999758341716293`, `R = 13.468323075955098`, `B = 2.1245756954797212e+57`, `I = 159.58714132523306` |

## 5. Detector D1 — before migration

Temporary mutation of a **`RadiusLoop`-only** numerical setting: GSL relative tolerance
`1e-10 → 1e-8`. `SingleStarSolveToTOVPoints` untouched; anchor verified unique.

| | Result |
|---|---|
| coverage test | **FAILS** (exit 1) |
| failing assertion | G6 — `0/16` bitwise |
| affected resolutions | **all 16** |
| worst relative discrepancy | **2.690e-10**, at `I(km^3)` @ res 15000 |
| example (res 10000) | `B` `2.12457569558348480e+57` vs canonical `2.12457569547972117e+57` |

Reverted; byte-identity confirmed by SHA-256.

**This is the last time a legitimate Path-1-only ordinary-star radial numerical detector can be
constructed.** After Commit B there is no second ordinary-star radial driver to perturb, and
that impossibility is the proof that ownership is resolved.

## 6. Commit A gates

| Gate | Result |
|---|---|
| coverage passes on old code | ✅ 9/9 assertions |
| D1 fires | ✅ |
| mutation reverted byte-identically | ✅ |
| production tree restored | ✅ no `CompactStar/` diff |
| no existing artifact changed | ✅ no `tests/baselines/` diff |
| `git diff --check` | ✅ clean |
| suite | **19/19** authenticated (198.97 s), **10/10** self-contained |

Commit A contains only the test, its CMake registration, and this evidence document.


---

# Commit B — migration and retirement

## 7. The out-of-domain audit (§18), measured before migrating

The legacy routine called `p_of_e(in_e_c)` **directly**, bypassing the central-density
floor/ceiling clamp that `Solve()` and `SolveToProfile()` both apply. `p_of_e` has its own,
*different* clamp — to `[eos_e_min, eos_e_max]`, returning `p_min`/`p_max` at the ends — whereas
the canonical primitive clamps `ec` to `[10·eos_e_min, 0.999·eos_e_max]` **before** inverting.
So the two agree inside that inner band and can differ outside it.

That difference was **measured, not reasoned about**, with a temporary probe (reverted
byte-identically):

> `GenTestSequence(ec = 1.0e8)` — below `eos_e_min = 1.65880787e+08` — **aborted the process,
> exit 134 (SIGABRT)**, before producing any star.

**Conclusion: the legacy out-of-domain behavior is not a well-defined public behavior; it is a
hard abort.** §18's first branch therefore applies, and no STOP condition arises. I4 guarantees
compatibility **on the ordinary valid domain only**, and the out-of-domain change — a clamped
star instead of a crash — is recorded as a deliberate improvement rather than claimed as
preservation. It affects no caller: `GenTestSequence` has none.

## 8. Migration

Inside `GenTestSequence`'s loop, this replaced the inline setup and the `RadiusLoop` call:

```cpp
std::vector<TOVPoint> tov_points;
SingleStarSolveToTOVPoints(in_e_c, tov_points);

for (const auto &tp : tov_points)
    n_star.Append(tp);
SurfaceIsReached();
```

**Preserved exactly**: `p_of_e_prec = 1e-9` (still set before the loop, still never restored);
`test.Modify(this, idx)` setting `radial_res` per iteration; the 16 resolutions and their order;
one integration per row; and the hook order `SurfaceIsReached → analysis->Analyze →
n_exp_cond_f → ExportNStarProfile → Sequence::Add → n_star.Reset`, followed by
`ExportSequence(..._TestSequence.tsv)` and `analysis->Export`. `Append`+`FinalizeSurface` was
**not** replaced by `BuildFromTOV` — I3 remains optional and untaken.

## 9. Post-migration result — bit-identical

| Check | Result |
|---|---|
| full-precision in-memory `SeqPoint` at all 16 resolutions, `(ec, M, R, pc, B, I)` at `%.17g` | **BIT-IDENTICAL** to the pre-migration dump |
| `gts_TestSequence.tsv` | **BYTE-IDENTICAL** — `e43bd05c1b22f7dbc26a39a10fe68d013bde5440eb3d5ce502eb83555245e511` before and after |
| coverage test | 9/9, `GEN_TEST_SEQUENCE EQUIVALENT TO CANONICAL PRIMITIVE` |

No tolerance was selected after the migration; none was needed.

## 10. `RadiusLoop` reference audit and deletion

Post-migration audit found **zero production callers**. Removed:

- the definition in `TOVSolver.cpp` (**155 lines**, including the I1 non-authoritative banner);
- the declaration in `TOVSolver.hpp`;
- the stale comment `"(by construction in RadiusLoop)"`, retargeted to the canonical primitive;
- the commented-out `RadiusLoop(r, y)` call left in an older disabled block;
- the `"identical to RadiusLoop"` / `"copy of RadiusLoop"` / `"exactly as in RadiusLoop"`
  comments inside `SingleStarSolveToTOVPoints`, which described a function that no longer exists.

The only surviving mentions are two comments that state, accurately, that the duplicate was
retired in 3E-I4.

## 11. One-owner structural proof (§21)

| Site | Function | Classification |
|---|---|---|
| `TOVSolver.cpp:2410` (driver), `:2457` (step ladder) | `SingleStarSolveToTOVPoints` | **the ONE ordinary visible-sector radial implementation** |
| `:2801`, `:2807` (drivers), `:2833` (ladder) | `RadiusLoopMixed` | two-fluid mixed-star core/mantle — **distinct physics, out of scope** |

> **ONE LIVE ORDINARY-STAR RADIAL NUMERICAL IMPLEMENTATION: `SingleStarSolveToTOVPoints`.**

All three orchestrators — `Solve(Axis)`, `SolveToProfile(target_M)`, `GenTestSequence(ec)` —
delegate to it.

## 12. Detector D1 — before and after

| | Result |
|---|---|
| **D1 BEFORE I4** | **FIRES** — RadiusLoop-only rel tol `1e-10 → 1e-8` breaks G6 to `0/16` bitwise, worst rel `2.690e-10` |
| **D1 AFTER I4** | **IMPOSSIBLE — ONE OWNER.** There is no second ordinary-star radial driver to perturb. Any mutation of `SingleStarSolveToTOVPoints` moves every orchestrator identically, so it is not a *Path-1-only* mutation and cannot express the duplicate-implementation fault. |

No duplicate was manufactured to force a mutation. **The impossibility is the proof that the
ownership problem is resolved** — the detector's own disappearance is the result.

## 13. Preservation results

| Check | Result |
|---|---|
| `tov_sequence_workflow_cmf` | **10/10** — `_Sequence.tsv` filename, header, order, units, row count, hooks unchanged |
| `tov_path_equivalence_cmf` | **67/67** |
| `baryon_number_cmf` | **14/14**, worst `\|ΔB\|/B = 0.000e+00` |
| six protected scientific artifacts | **byte-identical** |
| `baryon_number_dscmf1_reference.tsv` | **byte-identical** |
| `tov_path_equivalence_dscmf1.tsv` | **byte-identical** (numerically unchanged, as expected — I4 alters neither compared workflow's postprocessing) |
| Authenticated suite | **19/19** (199.68 s vs 198.97 s at Commit A — no regression) |
| Self-contained suite | **10/10** |
| Production diff | `TOVSolver.cpp`, `TOVSolver.hpp` **only** |

`NStar`, `Geometry.hpp`, `StarProfile`, `RotationSolver`, `MixedStar`, thermal and EOS code:
**untouched**. `SolveToProfile` untouched. The six `Solve()` callers untouched.

## 14. Governance closure — precise scope

`GOVERNANCE.md:70` condition 3: *"Two or more implementations are live and no document
establishes which is canonical."*

- ADR-0005 is **ACCEPTED** and names the authority — the "no document" half was discharged at I1.
- `Solve()`, `SolveToProfile()` and `GenTestSequence()` all delegate to that authority.
- `RadiusLoop` is **deleted**; the reference audit confirms no ordinary-star duplicate remains.
- All validation passes.

> **CONDITION #3: CLOSED for ordinary visible-sector TOV radial numerical ownership.**

**Explicitly NOT covered by this closure**, and not claimed:

- `MixedStar` / `Solve_Mixed` / `RadiusLoopMixed` two-fluid radial integration — a **distinct
  physics problem**, not a competing implementation of the ordinary-star path;
- dark-sector TOV (`ODE_Dark_Core`, `ODE_Dark_Mantle`);
- any future alternative solver;
- the two `NStar` construction styles (`BuildFromTOV` vs `Append`+`FinalizeSurface`), which are
  **value-equivalent and geometry-conformant but not textually consolidated** — that is optional
  increment I3, and it is a *postprocessing* question, not an authoritative-path question.

## 15. Explicit non-scope and deferred items

| Item | State |
|---|---|
| **3E-I3** — converge `Append`+`FinalizeSurface` onto `BuildFromTOV` | **OPTIONAL, not taken.** Reassess in the Phase-3 closeout; it does not block 3E. |
| mirror `M`/`R`/`z_surf` zeros | preserved (M1) — still `INTERNAL STATE ASYMMETRY — CURRENTLY UNOBSERVED` |
| `NStar::BaryonNumIntegrand(double)` | untouched — separate INV-14 defect |
| INV-04 | unchanged from its post-I2 status; I4 does not touch proper-volume ownership |
| INV-07 / Hartle | untouched; `I` bit-identical |
| `MixedStar`, candidates, ZakiLib, CONFIND | untouched |
