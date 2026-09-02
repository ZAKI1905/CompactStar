# Phase 3E-I2 — Path-1 ADR-0004 geometry conformance

> **STATUS: `PHASE-3E-I2 PATH-1 GEOMETRY CONFORMANCE COMPLETE`**
>
> Both ordinary visible-sector `NStar` construction paths now use the canonical
> `CompactStar::Geometry` mathematics. **`GOVERNANCE.md` fail-closed condition #3 is unchanged
> by this increment** — still discharged for the `Solve()` workflow, still partially open
> overall, pending I4.

| Field | Value |
|---|---|
| **Starting HEAD** | `a31643f2a8fb083f4dec62a949e0e99838c564bb` (14 ahead / 0 behind `master`, upstream equal, clean) |
| **Governing authority** | **ADR-0004 (ACCEPTED 2026-09-01)** §13; **ADR-0005 (ACCEPTED 2026-09-02)** §15 |
| **Change class** | structural / architecture, with a governed numerical component |
| **Production diff** | `CompactStar/Core/src/NStar.cpp` **only** |

---

## 1. Scope, including one correction to ADR-0005's shorthand

ADR-0005 §15 described I2 as migrating *"`FinalizeSurface`'s proper-volume factor."* The
authenticated source held a **second** Path-1 ADR-0004 nonconformance that the shorthand omitted:

`NStar::Append` independently computed

```
f = 1 - 2m/r ;  if (f <= 0) f = 1e-15 ;  Lambda = -0.5 * log(f)
```

That is governed by ADR-0004 and was deferred out of Phase 3D for exactly the same reason as the
baryon integrand — Path 1 had no coverage. Phase 3E-0 supplied the coverage, so both deferrals
ended together. **I2 migrated both sites and nothing else.** ADR-0005's accepted ownership
decision is unaffected.

## 2. The two migrations

| Site | Before | After |
|---|---|---|
| `NStar::Append` — Λ | local `denom` / `1e-15` clamp / `-0.5*log(denom)` | `CompactStar::Geometry::Lambda(r_km, m_km)` |
| `NStar::FinalizeSurface` — baryon integrand | `r.pow(2)` × `4π n_B` ÷ `(1−2m/r).sqrt()` × `1e54` | per-node `w_V = Geometry::ProperVolumeWeight(r,m)`, then `w_V · n_B`, then `× 1e54` |

The `FinalizeSurface` composition deliberately **mirrors the already-validated `BuildFromTOV`
form**, so the two ordinary paths now share one mathematical owner *and* one arithmetic ordering.

`FM3_TO_KM3` was **not** moved into `Geometry`: it is INV-14 baryon-number semantics, not `dV`.

Unchanged at both sites: the M☉→km conversion, `TOVPoint` contents, column order, species
handling, `Edit`/`Version` semantics, append order, interpolation setup, surface-guess updates,
and the integration bounds.

## 3. Search proof — Path 1 no longer owns either formula

| Query, scoped to the function | Result |
|---|---|
| `NStar::Append` (lines 652–790) for `1e-15`, `denom`, `-0.5 * std::log`, `1.0 - 2.0` | **no code matches** — only explanatory comment text |
| `NStar::FinalizeSurface` (lines 540–660) for `.sqrt()`, `1.0 - 2.0`, `pow(2)` | **no code matches** — only explanatory comment text |
| `NStar.cpp` for `1e-15` in code | one match: `:861`, `EvaluateNu`'s boundary condition — **not this measure** (ADR-0004 §4.4), deliberately untouched |

`git diff` touches exactly two hunks in `FinalizeSurface` and one in `Append`. **`BuildFromTOV`
calls neither migrated function**, which is the pivot for §5.

## 4. Λ — bit-identical

Phase 3E-0 measured Path-1 Λ ≡ Path-2 Λ **bitwise** before this migration. After it, Λ is still
bitwise between the paths: `max_profile_ulp = 0` over all 25 radial columns, at 14 central
densities and `radial_res` 5000/10000/20000. Path 2 is untouched, so **Path-1 Λ is bit-identical
pre/post**. No tolerance was selected, and none was needed.

## 5. B — the governed movement, and a stronger outcome than required

Path 2 is provably untouched, so the recorded *path-to-path* differences pin the Path-1 movement
exactly.

| Comparison | `rel_B` pre-I2 | `rel_B` post-I2 |
|---|---|---|
| `A:1.0Msun`, `A:1.4Msun`, `A:1.6Msun`, `A:2.0Msun` | 0 (bitwise) | **0 (bitwise)** |
| `B:node0` | `1.387e-16` | **0** |
| `B:node3` | `1.637e-16` | **0** |
| `B:node9` | `1.269e-16` | **0** |
| `C:res5000` | `1.640e-16` | **0** |
| `C:res20000` | `1.640e-16` | **0** |
| all other rows | 0 | **0** |

**Maximum Path-1 `B` movement: 1.640e-16**, against the `1.0e-15` predeclared in ADR-0004 §7.1
*before* any implementation existed. The bound was not widened.

**Path-1 and Path-2 `B` are now bitwise identical at all 17 comparisons.** That is stronger than
the gate requires, and it is *structural* rather than coincidental: both paths now call the same
`ProperVolumeWeight` and compose `w_V · n_B · 1e54` in the same order. The former ADR-0004
Path-1 conformance gap is **closed**, not merely bounded.

Path-1 `B` at the four anchors — unchanged, because those were already bitwise pre-I2:

```
1.0 M☉  1.27388873109076839e+57
1.4 M☉  1.83218336257875150e+57
1.6 M☉  2.12457569547972117e+57
2.0 M☉  2.74576306114479063e+57
```

**No assertion was strengthened to bit identity.** ADR-0005 §16 permits that only if the source
structure makes it a stable contract; the governed `1.0e-15` bound is retained, and the bitwise
result is reported as evidence rather than encoded as a new requirement.

## 6. Artifacts

| Artifact | Required | Result |
|---|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | byte-identical | **unchanged** |
| `tov_dscmf1_reference.tsv` | byte-identical | **unchanged** |
| `grid_convergence_cmf_1p6_debug.tsv` | byte-identical | **unchanged** |
| `grid_convergence_cmf_1p6_trajectory.tsv` | byte-identical | **unchanged** |
| `hartle_I_dscmf1_debug.tsv` | byte-identical | **unchanged** |
| `baryon_number_dscmf1_reference.tsv` | byte-identical, **not** re-baselined | **unchanged**; `baryon_number_cmf` reports worst `\|ΔB\|/B = 0.000e+00` |

All six are produced through Path 2, which I2 does not touch — so their invariance is expected by
construction, not luck.

**`tov_path_equivalence_dscmf1.tsv` — re-emitted, and here is why that is correct.** Its role was
authenticated before deciding: the CTest runs `tov_path_equivalence_cmf <root>` with **no
`--compare`**, the test never loads a baseline, and the file is written only under `--emit`. Its
header states `rel_* are |path1-path2|/|path2|`. It is therefore a **measurement record of
path-to-path differences**, not a compared golden — and leaving it asserting a conformance gap
that no longer exists would make it false.

| | |
|---|---|
| old SHA-256 | `a6b4cd799515b7dd261b7da6088fc06ea677137a993f44d46aafb1e508bf5786` |
| new SHA-256 | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |
| fields changed | **`rel_B` only**, in 5 of 17 rows, each `≤1.640e-16 → 0` |
| fields unchanged | `rel_M`, `rel_R`, `rel_I`, `max_profile_rel`, `max_profile_ulp`, node counts, `status` |

## 7. Workflow interface

`tov_sequence_workflow_cmf`: **10/10**, unchanged. Filename, header, column order, units, row
count and hook order all as before; `W6b` reports the identical worst relative agreement
(`3.137e-09`) because it compares the exported file against the *current* in-memory `SeqPoint`
rather than against frozen text — exactly the property §23 requires when the governed `B`
intentionally migrates.

## 8. Detector D4

Temporary mutation, verified-unique anchor, reverted with byte-identity confirmed by SHA-256.

| Mutation | Result |
|---|---|
| replace `Geometry::ProperVolumeWeight(r,m)` with the coordinate weight `4π r²` in the **migrated** `FinalizeSurface` (`BuildFromTOV` untouched) | **FIRES — 8 assertions, every one B-related**: the four `A:` anchors, `C:res5000/10000/20000`, and the `B1` sequence aggregate, at `\|ΔB\|/B` = **8.388e-02, 1.181e-01, 1.379e-01, 1.951e-01** |
| radial-column assertions during that mutation | **all 7 still PASS** |

The selectivity is the point: a proper-volume-only fault produces a proper-volume-only failure,
which is what makes the RADIAL-INTEGRATION vs BARYON-INTEGRATION classification trustworthy.

**No Λ-domain detector was added.** ADR-0004's `proper_volume_contract` already pins the
canonical domain behavior (14 fail-closed cases, regular-centre limit, `e^Λ` vs two independent
forms), and `Append` now routes through that same primitive — proven by the §3 search. Adding a
second mutation would have inflated the detector count without testing anything new.

## 9. Deferred and untouched

| Item | State |
|---|---|
| **3E-I3** — converge `Append`+`FinalizeSurface` onto `BuildFromTOV` | **deferred, optional** (ADR-0005 Q3 = P3) |
| **3E-I4** — `GenTestSequence` coverage → migration → `RadiusLoop` deletion | **deferred, required** for full fail-closed closure |
| mirror `M`/`R`/`z_surf` zeros | **preserved (M1)** — still `INTERNAL STATE ASYMMETRY — CURRENTLY UNOBSERVED` |
| `NStar::BaryonNumIntegrand(double)` | untouched — separate INV-14 defect, zero callers |
| `GenTestSequence`, `RadiusLoop` | untouched; `RadiusLoop` remains marked non-authoritative |
| `MixedStar`, `Solve_Mixed` | untouched |
| `DarkCore_Analysis`, BNV, decay, `RotochemicalCache`, O(Ω²) Hartle | untouched |
| `TOVSolver.cpp/.hpp`, `Geometry.hpp`, `StarProfile`, `RotationSolver`, thermal, EOS | untouched |

## 10. Invariant disposition

- **INV-03** — metric convention unchanged. Canonical-path conformance now covers **both**
  ordinary visible-sector `NStar` construction paths. Legacy nonconformance narrowed to the
  scalar accessor, `MixedStar`, the candidates, and `EvaluateNu`'s out-of-scope BC.
- **INV-04** — headline updated to `GOVERNED (ADR-0004 ACCEPTED) — BOTH ORDINARY VISIBLE-SECTOR
  NStar PATHS CONFORMED; MixedStar / CANDIDATE / SCALAR-ACCESSOR MIGRATIONS DEFERRED`.
  **Still not globally resolved** — the deferred debt is retained verbatim, not erased.
- **INV-14** — unchanged. **INV-07** — unchanged; `I` bit-identical throughout.

## 11. Validation summary

| Check | Result |
|---|---|
| Authenticated suite | **18/18 passed** (197.47 s vs 198.28 s at I1 — no regression) |
| Self-contained suite | **10/10 passed** |
| Equivalence harness | **67/67**, runtime 2 s |
| `baryon_number_cmf` | **14/14**, worst `\|ΔB\|/B = 0.000e+00` |
| `tov_sequence_workflow_cmf` | **10/10** |
| Path-1 Λ | **bit-identical** pre/post |
| Path-1 `B` movement | **≤1.640e-16** (bound `1.0e-15`, unwidened) |
| Path-1 vs Path-2 `B` | **bitwise, 17/17** |
| Six protected artifacts | **byte-identical**, not regenerated |
| Production files changed | `NStar.cpp` only |
| Detector D4 | **fires, 8 B-only assertions**; revert byte-identical |
| Fail-closed #3 | **unchanged** — discharged for `Solve()`, partially open overall |
