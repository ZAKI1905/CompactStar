# Phase 3E-I4 — `GenTestSequence` coverage and `RadiusLoop` retirement

> **STATUS (Commit A): `GEN_TEST_SEQUENCE COVERAGE ESTABLISHED — DUPLICATE RadiusLoop STILL
> LIVE`.**
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
