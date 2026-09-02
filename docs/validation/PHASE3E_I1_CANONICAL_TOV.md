# Phase 3E-I1 — Canonical TOV radial ownership: implementation and validation

> **STATUS: `PHASE-3E-I1 CANONICAL TOV RADIAL OWNERSHIP COMPLETE`**
>
> **`GOVERNANCE.md` fail-closed condition #3 is DISCHARGED for the ordinary visible-sector
> `Solve()` workflow, but NOT fully closed.** A second radial implementation remains reachable
> through the public, unexercised `GenTestSequence`. See §7 — this is the honest disposition, not
> a wording choice.

| Field | Value |
|---|---|
| **Starting HEAD** | `fbeb229fe9c54f67d3db4c4e0a5d43d1c9cfa19e` (12 ahead / 0 behind `master` `d2040d89…`, upstream equal, clean) |
| **Acceptance commit** | `f86618a` — `docs: accept canonical TOV ownership` |
| **Governing authority** | **ADR-0005 (ACCEPTED 2026-09-02)**; ADR-0004, ADR-0003, ADR-0001; INV-03, INV-04 |
| **Change class** | **STRUCTURAL / architecture** |
| **Baseline at start** | 17/17 authenticated, 10/10 self-contained (re-authenticated at this HEAD during 3E-G) |

---

## 1. Owner decisions, as implemented

| | Decision | Implemented as |
|---|---|---|
| **Q1** | `Solve()` preserved, subordinate | still public, still exported; no longer owns a radial integration |
| **Q2** | `_Sequence.tsv` contract preserved | unchanged filename, header, order, units; now guarded by a dedicated test |
| **Q3** | **P3 staged** — radial only | `Append`+`FinalizeSurface` untouched; no `BuildFromTOV` convergence |
| **Q4** | hooks preserved | `Analysis`, `n_exp_cond_f`, `virtual ExportNStarProfile` all unchanged |

## 2. What changed

**One function.** `TOVSolver::Solve` — the per-`Axis`-node body — now reads:

```cpp
// The ONE radial integration for this sequence member.
std::vector<TOVPoint> tov_points;
SingleStarSolveToTOVPoints(in_ax[idx], tov_points);

// Feed the points through the EXISTING Path-1 postprocessing, unchanged.
for (const auto &tp : tov_points)
    n_star.Append(tp);

SurfaceIsReached();
```

Removed from `Solve()`: its duplicate central-density clamp, its `p_of_e` conversion, its `y[]`
initialization, and its `RadiusLoop` call. All four now belong to the canonical primitive.

**Plus one documentation block** marking `TOVSolver::RadiusLoop` non-authoritative (§7).

Nothing else. `NStar` was not modified; `StarProfile`, `Geometry`, `RotationSolver`, `MixedStar`
and the thermal code were not touched.

## 3. `n_star` initialization is preserved, and why it needed checking

`Solve()` never initializes `n_star` inside its loop, so delegation could only be safe if the
initialization survives elsewhere. Authenticated:

- `TOVSolver::ImportEOS` calls `n_star.InitFromTOVSolver(this)` once
  (`TOVSolver.cpp:639`), which lays out the eight fixed columns, calls `ResetSpecies`, registers
  species labels from `eos_tab.extra_labels`, and calls `SetSurfaceScalars(0,0,0)`.
- `NStar::Reset()` (`NStar.cpp:783`) calls `prof_.Reset()`, clears `B_integrand` and sets
  `surface_ready = false`.
- `StarProfile::Reset()` (`StarProfile.hpp:875`) clears rows and `seq_point` and zeroes
  `M`/`R`/`z_surf`, but **deliberately does not clear `species_labels` or `species_idx`** — those
  two lines are commented out — so the column schema and species mapping survive every sequence
  member.

Delegation therefore changes only *when* `Append` is called (after the integration rather than
during it), not what state it is called against. **No `NStar` semantics were changed, and none
needed to be.**

## 4. The clamp is no longer duplicated

Before I1, `Solve()` and `SingleStarSolveToTOVPoints` each carried the same clamp. `Solve()`'s
copy is gone; the primitive's remains and is now the only one.

`central_eps_floor_factor`, the floor formula (`factor * eos_tab.eps.front()`) and the ceiling
formula (`0.999 * eos_tab.eps.back()`) are **unchanged**, and both were reading the same
`eos_tab`, so the effective central density is identical.

**The one observable difference is a log line.** The clamp warning now originates in
`SingleStarSolveToTOVPoints` and reads `"Requested eps(...) < floor(...) -> clamping."` instead of
`"Solve: requested eps(...) < floor(...) -> clamping."`. The warning still fires, at the same
threshold, for the same inputs. **No numerical behavior changes** — this is recorded here rather
than left for someone to discover in a diff.

## 5. Bit-identity evidence

Path 2 is provably untouched: no source change reaches `SingleStarSolveToTOVPoints`,
`NStar(points, labels)` or `BuildFromTOV`. The equivalence harness records, for every
experiment, the *relative differences between the two paths*. Re-emitting it after I1 gives:

> **`tov_path_equivalence_dscmf1.tsv` is BYTE-IDENTICAL to its pre-I1 baseline**
> (`a6b4cd799515b7dd261b7da6088fc06ea677137a993f44d46aafb1e508bf5786`).

Since Path 2 did not move and every `rel_M`, `rel_R`, `rel_B`, `rel_I`, `max_profile_rel` and
`max_profile_ulp` is unchanged, **every Path-1 value is bit-identical to pre-I1** — across all 25
radial columns (`r`, `m`, `nu'`, `p`, `eps`, `nB`, `nu`, `Lambda`, 17 species), all four
fixed-`ec` anchors, the ten-star sequence sweep, and `radial_res` 5000/10000/20000.

`67/67` assertions pass in the harness, unchanged from pre-I1.

| Quantity | I1 standard | Result |
|---|---|---|
| all 25 radial columns | BIT-IDENTICAL | **met** |
| `ec`, `pc`, `M`, `R`, `I` | BIT-IDENTICAL | **met** |
| Path-1 `B` | BIT-IDENTICAL (the ADR-0004 `1e-15` allowance belongs to I2, not I1) | **met** — unchanged; **none of the allowance was spent** |
| multi-star sequence | no state leakage | **met** |
| six protected artifacts | byte-identical | **met**, not regenerated |
| `tov_path_equivalence_dscmf1.tsv` | numerically unchanged | **byte-identical** |

## 6. `_Sequence.tsv` workflow contract

New test `tests/core/tov_sequence_workflow_cmf.cpp` (external-data). One test protects the
contract; deliberately **not** one test per caller.

| | Assertion | Result |
|---|---|---|
| W1 | file appears at exactly `<file>_Sequence.tsv` | pass |
| W2a/b | exactly six fields, exact names, **exact order** | pass |
| W3 | raw header line byte-exact | pass |
| W4 | row count `= Axis.res + 1` | pass |
| W5 | six numeric fields per row | pass |
| W6 | exported values match in-memory `SeqPoint` to export precision `%.8e` | pass, worst rel **3.137e-09** |
| W7 | `ec` column increasing and spanning the requested Axis | pass |

The header contract, verbatim:

```
ec(g/cm^3)    	 M(Sun)        	 R(km)         	 pc(dyne/cm^2) 	 B             	 I(km^3)       
```

**The rounded TSV is never used as a numerical oracle.** W6 runs the other way: the in-memory
full-precision `SeqPoint` is the reference and the file is the thing checked. In-memory values
are read through the existing `ExportNStarProfile` hook, the same observational mechanism 3E-0
established — no production API was added.

## 7. `RadiusLoop` disposition — and a correction to ADR-0005 §8

ADR-0005 §8 states: *"Its only live caller is `Solve()` itself."* **That is factually wrong, and
implementation found it.** `TOVSolver::GenTestSequence` (`TOVSolver.cpp:3244`) also calls
`RadiusLoop` (`:3308`).

`GenTestSequence` is a **public API that is compiled but unexercised**: its only repository
reference is commented out (`main/Test/tov_debug_main.cpp:196-203`), and it has zero test
coverage. It is the same category as `MixedStar`.

**Consequence — deletion is not trivially safe in I1.** Deleting `RadiusLoop` would require also
migrating `GenTestSequence`, which is uncovered code. Migrating unobserved code inside a
structural increment is precisely the risk `AGENTS.md` forbids, and the risk ADR-0004 §0-Q2
declined to take for `MixedStar`. §9 of the implementing brief provides for exactly this case.

**Taken: retain, mark non-authoritative, defer to I4.** `RadiusLoop` now carries a block comment
recording that it is non-authoritative as of ADR-0005, that `Solve()` no longer calls it, that
its sole remaining caller is the unexercised `GenTestSequence`, that no new callers may be added,
and that retirement is Phase 3E-I4.

**The accepted decision is unaffected.** Option A remains canonical and the staged retirement was
already ADR-0005's own preference (§8: *"no longer authoritative" is distinct from "deleted in
the first implementation commit"*). Only the **closure timing** changes.

## 8. Detector results

Every mutation was applied to a **verified-unique** source site and reverted with byte-identity
confirmed by SHA-256. None is committed.

| ID | Mutation | Result |
|---|---|---|
| **D2** | drop one `Sequence::Add` in the `Solve()` wrapper | **FIRES** — W4 (2 rows, expected 3) and W6a |
| **D3a** | rename `M(Sun)` → `M(Msun)` in `Sequence::Export` (not `MixedSequence`) | **FIRES** — W2b and W3 |
| **D3b** | `_Sequence.tsv` → `_Seq.tsv` at the `Solve()` site only | **FIRES** — W1 |

Two mutations were initially written against **non-unique anchors** (`"ec(g/cm^3)", "M(Sun)", …`
also matches `MixedSequence::Export`; `ExportSequence(in_dir + in_file + "_Sequence.tsv")` also
matches two mixed-star sites). Both were caught by an explicit uniqueness assertion and were
**never applied** — the resulting clean runs were not counted as detector passes. Re-targeted and
re-run.

### D1 — the impossibility proof, and its honest limit

ADR-0005 §16 predicted that after I1 the old "perturb a Path-1-only radial tolerance" mutation
would become **impossible to construct**. Structural audit of production:

| Site | Function | Scope |
|---|---|---|
| `TOVSolver.cpp:1760` | `RadiusLoop` | ordinary visible-sector — **still present** |
| `TOVSolver.cpp:2566` | `SingleStarSolveToTOVPoints` | ordinary visible-sector — **the canonical primitive** |
| `TOVSolver.cpp:2955`, `:2961` | `RadiusLoopMixed` | mixed-star core/mantle — out of scope |

**For the `Solve()` workflow: D1 IS impossible.** `Solve()` contains no GSL driver, no step
ladder and no `RadiusLoop` call; there is exactly **one** radial integration per sequence member,
and no Path-1-only tolerance exists to perturb.

**Globally: D1 is still constructible**, by perturbing `RadiusLoop` and reaching it through
`GenTestSequence`. The prediction holds for the workflow it was written about and does not hold
for the whole translation unit. Recorded as measured, not as predicted.

## 9. Fail-closed condition #3 — precise disposition

`GOVERNANCE.md:70`: *"Two or more implementations are live and no document establishes which is
canonical."*

- **The "no document" half is discharged**: ADR-0005 is ACCEPTED and names the authority.
- **The "two implementations" half is discharged for the ordinary visible-sector `Solve()`
  workflow**: it now has exactly one radial owner, and all six live callers reach it through that
  one owner.
- **It is NOT discharged globally**: `RadiusLoop` survives as an independent algorithm reachable
  through the public `GenTestSequence`.

> **Condition #3: DISCHARGED for the ordinary visible-sector `Solve()` workflow; PARTIALLY OPEN
> overall, pending 3E-I4.**

Explicitly **not** claimed: `MixedStar` canonicalization, dark-sector TOV canonicalization, or
any target-mass-search validation beyond existing evidence.

## 10. Deferred, unchanged, and untouched

| Item | State |
|---|---|
| **3E-I2** — Path-1 `FinalizeSurface` proper-volume → `CompactStar::Geometry` under `\|ΔB\|/B ≤ 1.0e-15` | **deferred**; the allowance is unspent |
| **3E-I3** — converge `Append`+`FinalizeSurface` onto `BuildFromTOV` | **deferred** (Q3 = P3) |
| **3E-I4** — migrate `GenTestSequence`, then delete `RadiusLoop` | **deferred**; now also the closure gate for condition #3 |
| mirror `M`/`R`/`z_surf` zeros | **preserved (M1)** — still `INTERNAL STATE ASYMMETRY — CURRENTLY UNOBSERVED`, not declared correct |
| `NStar::BaryonNumIntegrand(double)` (INV-14) | untouched |
| `Find_MomInertia` side effect, `RotationSolver` (INV-07) | untouched; `I` bit-identical |
| `SolveToProfile` — `N_coarse`, grid, stable-branch test, bisection, mass tolerance, fallback | untouched |
| `ImportEOS`, EOS splines, `p_of_e`, `GetEDens`, `GetRho`, `GetRho_i`, `PressureCutoff` | untouched |
| the six `Solve()` callers | **unmodified**; they compile and behave as before |
| **INV-04** | **NOT resolved** — Path 1 still uses the legacy inline proper-volume expression |

## 11. Validation summary

| Check | Result |
|---|---|
| Authenticated suite | **18/18 passed** |
| Self-contained suite | **10/10 passed** |
| Equivalence harness | 67/67, artifact **byte-identical** |
| Six protected artifacts | **byte-identical**, not regenerated |
| Radial integrations per sequence member | **exactly one** |
| Production files changed | `TOVSolver.cpp` only (plus tests/docs) |
| `NStar` production source | **unchanged** |
| Detector mutations | reverted **byte-identically** |
