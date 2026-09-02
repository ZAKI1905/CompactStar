# ADR-0005 — Canonical TOV numerical primitive and sequence workflow interface

| | |
|---|---|
| **Status** | **PROPOSED** |
| **Date drafted** | 2026-09-02 |
| **Change class** | **Structural / architecture** (`GOVERNANCE.md:51`) |
| **Drafted at** | `f3497ea82e5c19f92d89cacb76e179f085f4c53b` |
| **Roadmap increment** | Phase 3E (`docs/architecture/PHASE3_CONSOLIDATION_PLAN.md`) |
| **Fail-closed condition addressed** | **#3 — uncertain authoritative code path** (`GOVERNANCE.md:70`) |
| **Governing invariants** | INV-03, INV-04, INV-13, INV-14 (non-scope), INV-07 (non-scope) |
| **Prerequisite** | Phase 3E-0 — `TOV-PATH EQUIVALENCE VERIFIED` (`docs/validation/TOV_PATH_EQUIVALENCE.md`) |

> **This ADR proposes a contract. It writes no code and authorizes no implementation.**
> No production source, test, baseline, ADR status or invariant status was changed by the task
> that produced it.
>
> **The fail-closed condition remains OPEN.** Drafting an ADR does not close it; see §19.

---

## 1. Context

`GOVERNANCE.md:70` fail-closed condition #3 fires when *"two or more implementations are live and
no document establishes which is canonical."* CompactStar has had two live TOV integration paths
since before Phase 0, and no accepted document names either as the authority. Every scientific
harness built in Phase 2 protects one of them; until Phase 3E-0 the other had **no coverage at
all**.

Phase 3D deliberately declined to migrate the legacy path's proper-volume code for exactly this
reason (ADR-0004 §13): it is governed but was unobserved, and migrating unobserved code inside a
structural increment is the risk `AGENTS.md` §5 forbids. Phase 3E-0 built that coverage. This ADR
proposes the authority the coverage now makes safe to name.

## 2. Phase-3E-0 evidence (binding)

Measured on `DS(CMF)-1_with_crust` at identical EOS, identical central density and identical
numerical controls. **Do not reopen this numerical question unless the source changes.**

| Quantity | Result |
|---|---|
| all 25 radial columns (`r`, `m`, `nu'`, `p`, `eps`, `nB`, `nu`, `Lambda` + 17 species) | **bit-identical**, `max_profile_ulp = 0` |
| species labels and ordering | identical (both read `eos_tab.extra_labels`) |
| surface termination — node count, last `r`, `p`, `m` | identical |
| `SeqPoint` `ec`, `pc`, `M`, `R`, `I` | **bit-identical** |
| `SeqPoint` `B` | max `\|ΔB\|/B` = **1.640e-16**, inside ADR-0004's predeclared `1.0e-15` |
| `radial_res` 5000 / 10000 / 20000 | equivalent at each |
| ten-star Path-1 sequence sweep | 10/10 equivalent, **no state leakage** |
| central-density floor / ceiling clamp | identical achieved `ec` and `M` |

**Architecture classification: H2 VERIFIED — one TOV algorithm copied into two radial loops, with
different orchestration and output ownership.** H1 (two independent implementations) is refuted:
independent implementations do not agree bitwise across 2500–5300 nodes at three resolutions.

**The consequence that shapes this ADR:** canonicalizing the radial solve carries **no numerical
migration risk**. The genuine work is the *interface*.

## 3. Four architectural concepts, deliberately kept separate

Prior documents have called every public wrapper a "TOV implementation." That conflation is what
made the ownership question look harder than it is. This ADR proposes fixed terminology:

| Term | Component | Owns |
|---|---|---|
| **CANONICAL NUMERICAL PRIMITIVE** | `TOVSolver::SingleStarSolveToTOVPoints(ec, out_tov)` | the radial integration, and nothing else |
| **TARGET-MASS ORCHESTRATOR** | `TOVSolver::SolveToProfile(target_M, …)` | root-finding in central density over the primitive |
| **SEQUENCE / WORKFLOW ORCHESTRATOR** | `TOVSolver::Solve(Axis, dir, file)` | sweeping an Axis, accumulating a `Sequence`, and file output |
| **STAR-PROFILE CONSTRUCTOR** | `NStar::BuildFromTOV`, or `Append`+`FinalizeSurface` | turning `TOVPoint`s into a finished `StarProfile` |
| **WORKFLOW OUTPUT CONTRACT** | `<file>_Sequence.tsv` | the observable every live caller consumes |

**One numerical solver, several orchestration APIs** is not the same thing as "two TOV solvers."

## 4. Current call graph (authenticated at `f3497ea`)

```
Path 1  TOVSolver::Solve(Axis, dir, file)
          for idx = 0..Axis.res:                    # inclusive
            clamp ec -> p_of_e -> init y[]
            RadiusLoop(r, y)                        # NStar::Append per step
            SurfaceIsReached() -> NStar::FinalizeSurface()
            analysis? -> n_exp_cond_f? -> ExportNStarProfile
            Sequence::Add(n_star) -> n_star.Reset()
          ExportSequence(dir + file + "_Sequence.tsv")     # UNCONDITIONAL

Path 2  TOVSolver::SingleStarSolveToTOVPoints(ec, out_tov)
          clamp ec -> p_of_e -> init y[] -> radial loop -> out_tov
        NStar(points, labels) -> NStar::BuildFromTOV(...)

        TOVSolver::SolveToProfile(target_M, out_tov, out_labels)
          repeatedly calls SingleStarSolveToTOVPoints
```

**Shared outright** (one definition, both paths): `TOVSolver::ODE`, `p_of_e`, `GetEDens`,
`GetNuDer`, `GetRho`, `GetRho_i`, `PressureCutoff`, `EvaluateNu`, `Find_MomInertia`.

**Duplicated but measured identical**: central-density clamp, initial conditions, RK8PD setup,
initial GSL step, absolute/relative tolerance, radial grid, step_scale ladder,
append-before-break ordering, `TOVPoint` construction.

**Genuinely different**: the sequence/export orchestration; the star-profile construction style;
Path 1's unmigrated ADR-0004 baryon integrand; the mirror surface scalars.

## 5. The six live `Solve()` callers

Re-authenticated at `f3497ea`; unchanged since 3E-0.

| Caller | Numerical result consumed? | File side effect consumed? | Analysis / export hook? |
|---|---|---|---|
| `main/Test/spin_therm_evol_main.cpp:67` | no | **yes** (built-in `ExportSequence`) | no |
| `main/Test/tov_debug_main.cpp:165` | no — in-memory reads are commented out (`:153-160`) | **yes** | no |
| `main/Examples/sig_omega.cpp:64` | no | **yes** (`ExportSequence` at `:70`) | no |
| `main/Examples/Table_5-8_Glenn.cpp:50` | no | **yes** — `ExportSequence` at `:57`, then **re-reads its own TSV** at `:63` to plot | no |
| `main/Examples/coulatt.cpp:50` | no | unconditional export only | no |
| `main/Examples/polytrope.cpp:66` | no | unconditional export only | no |

`main/Examples/rotating_ns.cpp:62` is `r_solver.Solve(...)` on a **`RotationSolver`**
(`:41`) — not a TOV path. `sig_omega_rho.cpp` / `sig_omega_rho_nstar.cpp` are commented out.

**No caller anywhere in `main/` attaches an `Analysis` or sets `n_exp_cond_f`.** Every live caller
consumes only the file. None reads in-memory `Sequence` or `StarProfile` state after `Solve()`.

## 6. Alternatives for the canonical numerical primitive

### Option A — `SingleStarSolveToTOVPoints` is canonical ← **RECOMMENDED**

`Solve(Axis, …)` and `SolveToProfile(target_M, …)` both delegate to it. `RadiusLoop` ceases to be
an independent numerical authority.

- Already the most extensively validated path (Phase 2B-2 against the exact Schwarzschild
  interior to `3.5e-16` and the official CompOSE `eos.mr`; Phase 2B-4A grid convergence; Phase 3D).
- Already consumed by `SolveToProfile`, so one of the three orchestrators delegates to it today.
- Returns a neutral `std::vector<TOVPoint>` — it owns no file policy, no `NStar`, no `Sequence`.
- Proven exactly equivalent to `RadiusLoop` by 3E-0, so delegation is bit-identical by
  construction rather than by hope.
- Lets `Solve()` survive unchanged as a workflow facade.

**Source inspection found no material obstruction.** The one mechanical consequence: `Solve()`
currently appends to `n_star` *inside* the loop via `RadiusLoop`, so after delegation it must feed
the returned `TOVPoint`s to `NStar::Append` in a short loop. The clamp and `init_press` assignment
become the primitive's responsibility, which is safe — `PressureCutoff()` reads `init_press`, and
the primitive sets it before integrating exactly as `Solve()` does now. Memory cost of
materializing ~2600 `TOVPoint`s is ~0.35 MB per star.

### Option B — a new private `IntegrateSingleStar` primitive under both

Both `RadiusLoop` and `SingleStarSolveToTOVPoints` become wrappers around newly factored code.

**Rejected.** It creates a third implementation to review and validate in order to delete a
duplicate that Option A deletes directly. `SingleStarSolveToTOVPoints` is *already* the neutral,
side-effect-free primitive Option B would invent; extracting a fourth layer beneath it buys
symmetry, not safety. `AGENTS.md` forbids abstraction for aesthetic reasons.

### Option C — `SolveToProfile` is canonical

**Rejected as the primitive.** It combines radial integration with a target-mass root search, so
it cannot serve a fixed-`ec` sequence scan without running a search nobody asked for. It is an
*orchestrator*, and §3 keeps that concept separate.

### Option D — `NStar::SolveTOV_Profile` is canonical

**Rejected as the primitive.** It further combines EOS loading, the target-mass search and `NStar`
construction. It is the top of the convenience stack, not the bottom.

### Option E — keep both loops; name one authoritative by documentation only

**Rejected.** Zero migration cost, but it does not actually close fail-closed condition #3 in
substance: two live radial implementations remain, and the next divergence between them would be
silent. 3E-0's detector D1 exists precisely to catch that, and it should have something to catch.

## 7. Proposed contract

### 7.1 Canonical numerical primitive

**`TOVSolver::SingleStarSolveToTOVPoints(double ec, std::vector<TOVPoint>&)` is the single
authoritative TOV radial integration in CompactStar.** It owns the clamp, the initial conditions,
the GSL driver, the radial grid, the step schedule, the termination test and `TOVPoint`
construction. **It must not own file output, `NStar` construction, `Sequence` accumulation, or any
target-mass search.**

### 7.2 `Solve()` remains a supported workflow API, subordinate to the primitive

```
Solve(Axis, dir, file):
    for each ec in Axis:
        points = SingleStarSolveToTOVPoints(ec)      # the ONE numerical authority
        feed points through the star-finalization path
        analysis / export hooks as today
        Sequence::Add ; Reset
    ExportSequence(dir + file + "_Sequence.tsv")     # unchanged
```

`Solve()` stops being a numerical authority and becomes a **sequence/workflow orchestrator**.

### 7.3 `SolveToProfile` stays orchestration

Unchanged. It continues to call the primitive repeatedly. **`N_coarse`, the stable-branch
definition, the mass tolerance and the fallback policy must not change in Phase 3E** — those are
separate numerical-method questions with their own evidence requirements.

### 7.4 `_Sequence.tsv` is a compatibility contract

Preserve exactly, unless the owner later authorizes an interface change:

- **filename**: `<file>_Sequence.tsv`, written unconditionally at the end of `Solve()`
- **schema and order**: `ec(g/cm^3)`, `M(Sun)`, `R(km)`, `pc(dyne/cm^2)`, `B`, `I(km^3)`
- **format**: `Zaki::File::VecSaver` 1-D table, `%-14s` header fields, tab-separated
- **units and precision semantics** as emitted today

**Classification: PUBLIC / WORKFLOW INTERFACE**, not a debug artifact. It is the only observable
every live caller consumes, and `Table_5-8_Glenn.cpp` re-reads its own emitted TSV as program
input. The canonical primitive must **not** own this output; the workflow wrapper does.

### 7.5 Do not migrate the six callers

Replacing all six `Solve()` calls with bespoke loops around the primitive would duplicate sequence
orchestration in six places, force six caller changes, risk the TSV contract, and deliver **no
numerical benefit**. `Solve()` is retained as a compatibility facade. A modern API, if ever
wanted, can be added alongside. **Caller churn is not justified by the age of an API.**

## 8. `RadiusLoop` disposition

Once `Solve()` delegates, `TOVSolver::RadiusLoop` is no longer an independent numerical authority.
Its only live caller is `Solve()` itself.

The ADR distinguishes two things that are easy to conflate:

- **"no longer authoritative"** — required by this ADR, achieved the moment `Solve()` delegates;
- **"deleted"** — a separate, staged step (3E-I4) taken only after a reference audit proves no
  remaining use.

Preferred long-term state: **one radial integration implementation in the codebase.** The
implementation may combine deletion with I1 if the audit shows it is trivially safe, or stage it
for reviewability. Either is acceptable; leaving a permanently unreachable duplicate is not.

## 9. Star-profile postprocessing — P1 / P2 / P3

Numerical consolidation does **not** automatically resolve the two `NStar` construction styles.
3E-0 proved them value-equivalent for every physical sequence field except the governed `B`
difference and one interface asymmetry (§10).

| | **P1 — minimal** | **P2 — full convergence** | **P3 — staged** ← **RECOMMENDED** |
|---|---|---|---|
| What | `Solve()` delegates the radial solve; keeps `Append`+`FinalizeSurface` | `Solve()` also constructs each star via `BuildFromTOV` | I1 = radial only; a later increment decides postprocessing |
| Radial duplication removed | ✅ | ✅ | ✅ (at I1) |
| `NStar` construction duplication removed | ❌ | ✅ | deferred, explicitly |
| Mirror scalars | stay zero | become populated | unchanged at I1 |
| Hooks (`Analysis`, export) | trivially preserved | must be re-centred; they act on the solver's internal `n_star` | preserved at I1 |
| Migration risk | lowest | highest surgery in one step | lowest per step |

**P3 is recommended.** It resolves the ownership question — which is what fail-closed condition #3
actually asks — in the smallest step that genuinely resolves it, and leaves the *second*,
independent question (is converging `Append`+`FinalizeSurface` onto `BuildFromTOV` worth the API
surgery?) to be decided on its own evidence. **Total textual deduplication is not required to
establish authoritative ownership**, and treating it as required is how a structural increment
turns into a rewrite.

## 10. Mirror surface scalars

3E-0 measured: `NStar::Reset` zeroes `StarProfile`'s `M`, `R`, `z_surf` via
`SetSurfaceScalars(0,0,0)`; `BuildFromTOV` sets them; **`FinalizeSurface` never calls
`SetSurfaceScalars` at all.** So Path 1 leaves them at zero while its `SeqPoint` `m`/`r` are
correct. No current `Solve()` caller reads them.

**Classification: INTERNAL STATE ASYMMETRY — CURRENTLY UNOBSERVED.** It is *not* declared correct
merely because nothing observes it.

| | |
|---|---|
| **M1** preserve the zero state through 3E for strict behavioral minimality | keeps I1 bit-identical in every field |
| **M2** allow them to become populated when/if `Solve()` adopts `BuildFromTOV`, governed as an intentional internal-consistency correction | couples the fix to the increment that already touches this code and has coverage |
| **M3** add `SetSurfaceScalars` to `FinalizeSurface` independently | adjacent cleanup with no covering increment |

**Recommendation: M1 during 3E-I1, then M2 if and when postprocessing is unified.** M3 is
explicitly **not** recommended — not because the zeros are right, but because a correction should
land inside an increment that already has the coverage to detect what it changes. Known-incomplete
metadata should not be preserved forever for bit-identity's sake; it should be corrected
deliberately, once, under observation.

## 11. `Analysis` and profile-export hooks

3E-0 established that no `main/` caller attaches an `Analysis` or sets `n_exp_cond_f`, and that
the equivalence test is the only current user of the profile-export hook.

**Classification: LEGACY PUBLIC HOOK — CURRENTLY UNUSED.**

**Recommendation: preserve during 3E.** Canonical-TOV consolidation must not broaden into public
API deletion; retirement is a separate classification exercise for Phase 3F. Note that the
equivalence harness depends on `ExportNStarProfile` remaining virtual and `n_exp_cond_f`
remaining settable — removing them would remove the test that protects this migration.

## 12. Path-1 ADR-0004 conformance

ADR-0004 §13 deferred Path 1's proper-volume migration *"pending new coverage first."* **That
coverage now exists.** The deferral should therefore end.

The 3E implementation should migrate `NStar::FinalizeSurface`'s baryon-number metric factor —
today the legacy inline `/= (1.0 - 2.0*m/r).sqrt()` — to `CompactStar::Geometry`, under the
**already-predeclared and unwidened** `|ΔB|/B ≤ 1.0e-15` (ADR-0004 §7.1).

Must **not** change: `FM3_TO_KM3`, baryon-density semantics (ADR-0001), the interpolation method,
or the integration bounds.

**INV-14 explicit non-scope.** `NStar::BaryonNumIntegrand(double)` — the zero-caller scalar
accessor whose missing `1e54` is a separate INV-14 defect — **must not be repaired here.**
Repairing it would silently convert its output units inside a structural increment.

## 13. Hartle-`I` side effect

Both finalization paths currently compute `SeqPoint::I` via `Find_MomInertia()` →
`RotationSolver::FindNMomInertia()`. 3E-0 measured `I` **bit-identical** across both paths.

Do not remove or relocate this side effect silently. The implementation must keep `I`
bit-identical. **No claim about physical Ω or J is made; INV-07 remains unresolved.**

## 14. Preservation standards (predeclared)

| Target | Standard | Basis |
|---|---|---|
| all 25 radial columns | **BIT-IDENTICAL REQUIRED** | 3E-0 measured exact equality over 14 central densities, 3 resolutions, 17 species |
| `SeqPoint` `ec`, `pc`, `M`, `R` | **BIT-IDENTICAL REQUIRED** | same |
| `SeqPoint` `I` | **BIT-IDENTICAL REQUIRED** | §13 |
| `SeqPoint` `B`, while bringing Path 1 into ADR-0004 conformance | **`\|ΔB\|/B ≤ 1.0e-15`** | ADR-0004 §7.1, predeclared before that implementation existed; **must not be widened** |
| `_Sequence.tsv` filename, header, order, units | **UNCHANGED** | §7.4 |
| species labels and ordering | **UNCHANGED** | ADR-0001 |
| the six existing protected artifacts | **byte-identical**, *except* any artifact that genuinely contains a migrated Path-1 `B` value — **authenticate before specifying any exception** | none is currently known to carry a Path-1 `B`; every existing artifact is produced through Path 2 |
| `tov_path_equivalence_dscmf1.tsv` | remains authoritative; must still report `EQUIVALENT` | 3E-0 |
| `N_coarse`, stable-branch rule, mass tolerance, fallback policy | **UNCHANGED** | §7.3 |
| `ImportEOS`, EOS splines, `gsl_interp_accel`, `p_of_e`, `GetEDens`, `GetRho`, `GetRho_i`, `PressureCutoff` | **UNCHANGED** | already shared; measured equivalence gives no reason to touch them |

## 15. Proposed migration sequence

**3E-I1 — canonical radial solve.** `Solve(Axis)` delegates to `SingleStarSolveToTOVPoints`;
Path-1 `Append`+`FinalizeSurface` preserved; `Sequence` accumulation, hooks and `_Sequence.tsv`
preserved; `RadiusLoop` ceases to be a numerical authority. **All raw output bit-identical.**

**3E-I2 — Path-1 ADR-0004 conformance.** Migrate `FinalizeSurface`'s proper-volume factor to
`CompactStar::Geometry`. `B` within `1.0e-15`; everything else bit-identical.

**3E-I3 — optional profile-construction consolidation.** Only if the analysis at that point shows
it is worth doing. If it would create unnecessary API churn, defer it while clearly subordinating
the legacy postprocessing path.

**3E-I4 — `RadiusLoop` retirement.** Delete only after a reference audit proves no use remains.
May be combined with I1 if trivially safe.

## 16. Validation and detector plan (to be executed by the implementation, not here)

**Validation.** The existing `tov_path_equivalence_cmf` remains the primary authority and must
still report `EQUIVALENT`. Add a **file-interface test** that runs `Solve()` in a scratch
directory and verifies the exact output filename, header text, column ordering, row count, and
that the emitted values match the in-memory `Sequence` to export precision. That test is about
**interface compatibility, not TOV physics** — rounded TSV text must never be the scientific
numerical oracle.

**Caller compatibility.** One focused test protecting the API/file contract is sufficient;
**do not write one test per caller.** The six-caller audit is the evidence for *why* the contract
must survive, not a specification for six tests.

**Detectors** (plan; each must fire, and none may be a mutation that cannot legitimately fire):

| ID | Mutation | Must fail |
|---|---|---|
| **D1** | after consolidation, alter a Path-1-only radial tolerance | the equivalence test — **if any duplicate radial implementation survives, this is what catches it** |
| **D2** | break the `Solve()` wrapper's sequence row accumulation | the workflow test |
| **D3** | change `_Sequence.tsv` header, column order, or filename | the interface test |
| **D4** | remove the proper-volume factor from the migrated `FinalizeSurface` | the `B` equivalence assertion |
| **D5** | remove the `Find_MomInertia` side effect from one path | the `I` equivalence assertion |

D1 is the load-bearing one: after I1 it should become **impossible to construct**, because there
is only one radial implementation left to perturb. That impossibility is the proof that ownership
was actually resolved.

## 17. Consequences

**Positive.** One radial TOV implementation. Fail-closed condition #3 becomes closable. The
ownership boundary between numerics, star construction, orchestration and file output becomes
explicit and named. Path 1 gains ADR-0004 conformance. No caller changes; no file-format change.

**Negative / accepted.** `Solve()` remains a large legacy API surface. Two `NStar` construction
styles persist past I1 under P3. The mirror surface scalars stay wrong-but-unobserved until a
covered increment fixes them. Two currently-dead public hooks are retained.

**Neutral.** `SolveToProfile` and `NStar::SolveTOV_Profile` keep their present shape.

## 18. Explicit non-scope

`MixedStar` and `Solve_Mixed`; dark-sector paths; `RotationSolver` and the O(Ω²) Hartle work
(INV-07); the rotochemical and BNV candidates (`GOVERNANCE.md` §5); `NStar::BaryonNumIntegrand`'s
INV-14 defect; the solar-mass authority; the target-mass search parameters; EOS import and
interpolation; ZakiLib; CONFIND; ROOT. **No invariant status changes as a result of this ADR.**

## 19. When fail-closed condition #3 can close

Explicitly staged, because drafting is not deciding and deciding is not doing:

| Stage | Condition #3 |
|---|---|
| **ADR-0005 PROPOSED** (now) | **OPEN** — a proposal carries no authority (`docs/adr/README.md`) |
| **ADR-0005 ACCEPTED** | **STILL OPEN** — an authority is named, but two live radial implementations still exist in source |
| **3E-I1 lands and validates** — `Solve()` delegates, the duplicate radial authority is retired or subordinated, the equivalence and interface tests pass | **CAN CLOSE** |

## 20. Owner adjudication questions

Only genuine design decisions are surfaced. Implementation trivia is not put to the owner.

### Q1 — Preserve `Solve()`?

Should `TOVSolver::Solve(Axis, dir, file)` remain a supported public workflow API after
canonicalization, implemented as a wrapper around the one numerical primitive?

**Recommendation: YES.** All six live callers rely on its sequence-file workflow. Removing it
forces six caller migrations for zero numerical benefit, and one caller (`Table_5-8_Glenn`) reads
the emitted TSV back as program input.

### Q2 — Preserve the `_Sequence.tsv` contract?

Should the filename, schema, column order and units remain backward compatible?

**Recommendation: YES.** It is the only observable consumed by every live `Solve()` caller.

### Q3 — How far should postprocessing consolidation go in Phase 3E?

**P1** minimal / **P2** full / **P3** staged.

**Recommendation: P3.** Remove the duplicate radial implementation first while preserving the
wrapper; then decide separately whether converging `Append`+`FinalizeSurface` onto `BuildFromTOV`
buys enough architectural simplification to justify touching the hooks and internal state. This
also determines the mirror-surface-scalar disposition (§10): **M1** under I1, **M2** if P2 is
later taken.

### Q4 — The dead `Analysis` and profile-export hooks

Remove now, deprecate, or preserve?

**Recommendation: PRESERVE during 3E; classify in 3F.** They are unused by production callers but
the equivalence harness depends on the export hook remaining virtual and `n_exp_cond_f` remaining
settable.
