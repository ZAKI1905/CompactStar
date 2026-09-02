# Phase 3 closeout — audit, debt classification, and merge readiness

> **FORMAL STATUS: `PHASE 3 COMPLETE — MERGE READY`**
> **3E-I3: `DEFERRED — OPTIONAL`** (not justified now; explicit re-trigger condition in §8)
>
> This record is the durable authority for deciding the merge. It is an **audit**: no production
> source, test, baseline or CMake file was changed by the task that produced it, and no
> architectural cleanup was performed merely because it was available. The branch was **not**
> merged.

| Field | Value |
|---|---|
| **Audited HEAD** | `99d14884a371486ca9a66be14cf3496edb9a5265` |
| **Canonical `master`** | `d2040d899d10a3d0da54a5a2facd7d9cbec34850` |
| **Relationship** | 17 ahead / 0 behind; merge-base **is** `master` (fast-forwardable); upstream equal; clean tree |
| **Suite at this HEAD** | **19/19** authenticated, **197.89 s**; **10/10** self-contained, **13.44 s** |
| **Authority read** | `GOVERNANCE.md`; ADR-0001…0005 (all ACCEPTED); `SCIENTIFIC_INVARIANTS.md`; roadmap; `CURRENT_ARCHITECTURE.md`; `PHASE3_CONSOLIDATION_PLAN.md`; all eight Phase-3 validation records; `PHASE2B_CLOSURE.md`; `AGENTS.md` |

Per-test runtimes at this HEAD: `grid_convergence_cmf` 106.20 s, `heat_capacity_real_star`
51.34 s, `passive_cooling_regression` 11.03 s, `heat_capacity_v1` 9.83 s,
`hartle_moment_inertia_cmf` 4.84 s, `tov_gen_test_sequence_cmf` 3.24 s, `tov_reference_cmf`
3.13 s, `baryon_number_cmf` 2.92 s, `tov_path_equivalence_cmf` 1.92 s,
`cache_thermal_contract` 1.84 s, `photon_cooling_conformance` 1.02 s, all others ≤ 0.42 s.

---

## 1. Authenticated durable artifacts (unchanged by this audit)

| Artifact | SHA-256 |
|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `tov_dscmf1_reference.tsv` | `ba9f6ee51e501e5e5a2133f72d3d16f351e5c721eb3f7a7c04e4d922fbc13e28` |
| `grid_convergence_cmf_1p6_debug.tsv` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |
| `hartle_I_dscmf1_debug.tsv` | `ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15` |
| `baryon_number_dscmf1_reference.tsv` | `8da5799d21da2017dd7dc49dfec8571ade6efba22846a652796118f248d4a646` |
| `tov_path_equivalence_dscmf1.tsv` | `bbf61e5fddb4709500f22a1eb11b1e20554f7463376619e86e96ea0a2540d871` |

## 2. Phase-3 increment map — reconstructed from the validation records, not commit messages

| Increment | Governed goal | What actually landed | Class | Evidence | Status |
|---|---|---|---|---|---|
| **3A** | one owner for exactly-duplicated unit constants | `CompactStar/Units.hpp` (`KM3_TO_CM3`, `MEV_FM3_TO_ERG_CM3`), 4 sites | **engineering / bit-identical** | `PHASE3A_UNIT_DUPLICATES.md` | ✅ |
| **3B** | one cache-invalidation rule; provenance | ADR-0003 accepted; `ProfileProvenance`, `GeometryCache` provenance, `(identity, version)` keys, `C_⋆`/neutrino geometry keys, `StarContext` column re-binding; five hazards → enforced CTests | **structural** (fail-closed #4) | `CACHE_CORRECTNESS.md`, ADR-0003 record | ✅ |
| **3C** | one authority for `k_B` | four consumers → `Zaki::Physics::K_BOLTZ_EV`; three thermal goldens re-baselined with hashes | **numerical-method / constant-authority** | `PHASE3C_BOLTZMANN_AUTHORITY.md` | ✅ |
| **3D** | one owner for the proper-volume measure | ADR-0004 accepted; `CompactStar/Geometry.hpp`; `BuildFromTOV` + `GeometryCache` conformed; hybrid domain contract | **structural + scientific-semantic** (domain behavior governed) | `PHASE3D_PROPER_VOLUME.md` | ✅ (validated path) |
| **3E-0** | measure the two live TOV paths | 25 columns bit-identical at 14 densities × 3 resolutions; H2 verified; six callers audited | **measurement** (test/doc only) | `TOV_PATH_EQUIVALENCE.md` | ✅ |
| **3E-G** | propose TOV authority | ADR-0005 drafted (PROPOSED) | **governance** | ADR-0005 | ✅ |
| **3E-I1** | `Solve()` subordinate to the primitive | ADR-0005 accepted; `Solve()` delegates; `_Sequence.tsv` contract test | **structural** (fail-closed #3) | `PHASE3E_I1_CANONICAL_TOV.md` | ✅ |
| **3E-I2** | Path-1 geometry conformance | `Append` Λ and `FinalizeSurface` proper volume → `Geometry`; Path-1/Path-2 `B` bitwise | **structural** (governed 1e-15 numerical component) | `PHASE3E_I2_PATH1_GEOMETRY.md` | ✅ |
| **3E-I4** | retire the duplicate radial loop | `GenTestSequence` covered *before* migration, then migrated; `RadiusLoop` deleted; one ordinary-star radial owner | **structural** (fail-closed #3 closure) | `PHASE3E_I4_RADIUSLOOP_RETIREMENT.md` | ✅ |
| **3E-I3** | (optional) converge the two `NStar` constructors | not taken | — | §8 | deferred, optional |
| **3F** | classify dead/unreachable/legacy; closeout | this record | **documentation** | this document | ✅ |

Phase 3 was **not** uniformly "refactoring": 3C changed a physical constant's authority
(re-baselining three goldens under a predeclared, unwidened bound), and 3D replaced six mutually
inconsistent degenerate-input behaviors with one governed contract. Both were classified and
validated as such at the time.

## 3. Roadmap goal disposition (verified against the final source)

The ratified Phase-3 body lists five items; the ratified exit criterion is one sentence
(`MODERNIZATION_ROADMAP.md`: *"One authoritative owner per quantity; baselines still pass."*).

| Goal | Disposition | Basis |
|---|---|---|
| **A. Canonical TOV numerical ownership** | **SATISFIED** | one ordinary-star radial implementation (`SingleStarSolveToTOVPoints`); `Solve`, `SolveToProfile`, `GenTestSequence` delegate; `RadiusLoop` deleted; verified structurally in §14 |
| **B. Unit-conversion ownership** | **PARTIALLY SATISFIED BY DESIGN** | exact duplicates (3A) and `k_B` (3C) have one owner each. **Solar mass in km remains dual-authority** (`SUN_M_KM` in 9 files vs `GSL_CONST_CGSM_SOLAR_MASS` in 2) — **deferred out of Phase 3 by the ratified roadmap text itself** ("pending owner or ADR adjudication"). Not silently promoted; see §12 |
| **C. Proper-volume ownership** | **SATISFIED for ordinary visible-sector paths; DEFERRED BY GOVERNED DECISION elsewhere** | ADR-0004 §0-Q2 (accepted) defers `MixedStar`; `GOVERNANCE.md` §5 protects the two *core* candidates; the BNV/Decay/DarkCore modules are project-specific extensions outside core scope (owner clarification, §7a); INV-14 accessor separately governed |
| **D. Cache invalidation / provenance** | **SATISFIED** (scope: profile-derived caches) | ADR-0003 implemented and enforced; INV-12 resolved in that scope |
| **E. Dead / unreachable-code classification** | **SATISFIED** (classification), with **no deletions** beyond `RadiusLoop`, which had a reference audit and coverage first | §7 |

The roadmap also lists *"open for separate consideration in this phase"*: centralizing the
thermal-energy balance. It is explicitly **not** a prerequisite and requires its own ADR; it was
not taken up and is **not** an exit criterion.

## 4. ADR state

| ADR | Status | Implemented | Post-acceptance decision changes |
|---|---|---|---|
| 0001 species semantics | ACCEPTED 2026-08-31 | yes (Phase 2A) | none |
| 0002 heat-capacity ownership | ACCEPTED 2026-08-31 | yes (Phase 2A, Pattern A) | none |
| 0003 cache provenance | ACCEPTED 2026-09-01 | yes (3B) | none |
| 0004 proper-volume ownership | ACCEPTED 2026-09-01 | yes (3D + I2) | none — Q1/Q2/Q3 unchanged; two implementation records appended |
| 0005 canonical TOV primitive | ACCEPTED 2026-09-02 | yes (I1, I2, I4) | none — Q1–Q4 unchanged; one **factual** correction to §8 (the `GenTestSequence` caller ADR-0005 missed), marked in place; three implementation records appended |

No accepted decision was rewritten. Every implementation record was appended, not edited into
the decision text.

## 5. Invariant audit

| INV | Current headline | Changed in Phase 3? | Phase-3 blocker? | Later phase |
|---|---|---|---|---|
| 01 species semantics | GOVERNED (ACCEPTED) | no | no | — |
| 02 unit boundaries | VERIFIED CURRENT BEHAVIOR | **note only** (3C, at owner direction; headline deliberately unchanged) | no | solar-mass authority — deferred |
| 03 metric convention | VERIFIED CURRENT BEHAVIOR (domain governed by ADR-0004) | **yes** (3D, I2) | no | — |
| 04 proper volume | GOVERNED — both ordinary paths conformed; MixedStar/candidate/accessor deferred | **yes** (3D, I2) | no | MixedStar modernization |
| 05 radial centre | VERIFIED ⚠ (rotation half to re-audit) | no | no | Phase 4 |
| 06 surface | VERIFIED CURRENT BEHAVIOR | no | no | — (EOS-floor surface is a recorded 2B limitation) |
| 07 Hartle normalization | **UNRESOLVED** | no | no | **Phase 4** |
| 08 O(Ω²) order | INTENDED BUT UNVERIFIED | no | no | Phase 4 |
| 09 fixed-ε_c vs equilibrium derivative | INTENDED BUT UNVERIFIED | no | no | Phase 4/5 |
| 10 thermal redshift | VERIFIED CURRENT BEHAVIOR | no | no | — |
| 11 chemical-imbalance frame | **UNRESOLVED** | no | no | **Phase 5** |
| 12 cache invalidation | ✅ RESOLVED for profile-derived caches | **yes** (3B) | no | — |
| 13 interpolation | VERIFIED CURRENT BEHAVIOR | no | no | (observed order ≈0.70 recorded in 2B; not a Phase-3 item) |
| 14 baryon number | VERIFIED CURRENT BEHAVIOR (accessor defect recorded) | no | no | accessor: deferred debt (§9) |
| 15 heat capacity | GOVERNED (ACCEPTED) | no | no | — |
| 16 Direct Urca | VERIFIED CURRENT BEHAVIOR | no | no | Phase 5 (muon channel) |

Phase 3 changed exactly three headlines (INV-03, INV-04, INV-12) and one note (INV-02), each in
the increment that earned it. **No invariant was silently promoted**, and Phase 3 is not asked to
resolve Phase-4/5 physics.

## 6. Fail-closed condition audit (`GOVERNANCE.md` §3)

| # | Condition | Disposition | Scope / evidence |
|---|---|---|---|
| 1 | Ambiguous units | **DISCHARGED** for Phase-3 scope | INV-02 verified; the solar-mass split is *known and recorded*, not ambiguous |
| 2 | Ambiguous state meaning | **ACTIVE — OUT OF SCOPE** | INV-11 (chemical-imbalance frame) is a Phase-5 quantity; nothing in Phase 3 touched it |
| 3 | Uncertain authoritative code path | **DISCHARGED** for ordinary visible-sector TOV radial integration (I4) | `MixedStar`/`RadiusLoopMixed` is **NOT APPLICABLE** to #3: it is a single implementation of a *different* two-fluid problem, not a competing implementation of the same path |
| 4 | Uncertain cache validity | **DISCHARGED** for profile-derived caches (3B) | INV-12 scope is precise; no universal claim about every cache |
| 5 | Uncertain ownership | **DISCHARGED** for proper volume (ordinary paths), heat capacity, `k_B`; **ACTIVE — OUT OF PHASE-3 SCOPE** for the solar mass, which the ratified roadmap deferred | no document names a solar-mass authority; Phase 3 did not touch it; §12 |
| 6 | Absent validation for a scientific-semantic change | **DISCHARGED** | 3C and 3D each had baselines and validation; §5 candidates were neither activated nor modified |
| 7 | Source-of-truth disagreement | **DISCHARGED** | authenticated at every task; upstream equal at this HEAD |

**Does any active condition block the Phase-3 merge? No.** Both active conditions (#2 for
INV-11, #5 for the solar mass) attach to quantities Phase 3 neither touched nor claimed, and both
are explicitly assigned elsewhere by ratified documents.

## 7. Dead / unreachable / legacy / candidate classification (the 3F deliverable)

Controlled vocabulary: LIVE — CANONICAL · LIVE — SUPPORTED WRAPPER · LIVE — LEGACY ·
COMPILED — UNEXERCISED · PUBLIC — CURRENTLY UNUSED · **SUPPORTED EXTENSION HOOK — CURRENTLY
UNUSED BY CORE WORKFLOWS** · **PROJECT-SPECIFIC EXTENSION — UNEXERCISED** · UNREACHABLE
SCIENTIFIC CANDIDATE — PROTECTED · NOT COMPILED SCIENTIFIC CANDIDATE — PROTECTED · DEAD CODE ·
GENERATED ARTIFACT · DOCUMENTATION DEBT.

> **Owner architectural clarification, received during this audit and binding for
> classification (not a change to any accepted ADR):** (1) `Analysis` is an *intentional
> extension mechanism* — it lets a project-specific calculation run on each `NStar` **while**
> `TOVSolver` is looping, before the profile is reset, precisely to avoid the
> write-to-disk / re-read / re-interpolate cycle; the absence of a current `main/` caller does
> not make it obsolete. (2) The BNV / decay / dark-core directories are **project-specific
> research modules** that *consume* CompactStar's core machinery; they are not automatically
> future core-physics candidates, must not be required to validate or migrate for **core**
> closure, and must not be deleted for being unexercised. (3) No restructuring of those modules
> is authorized here. Evidence agrees with the intent: `DarkCore_Analysis`, `BNV_Analysis` and
> `Decay_Analysis` all **derive from `Core::Analysis`** — they are the hook's designed consumers —
> and `GOVERNANCE.md` §5's named candidate commit `675b4a9` touched only `RotationSolver` and
> `Rotochemical*`, never BNV/Decay/DarkCore. An earlier draft of this audit, and the INV-03/INV-04
> wording written in 3D/I2, conflated the two categories; both are corrected below and in §19.

| Item | Classification | Evidence | Delete? | Blocks Phase 3? |
|---|---|---|---|---|
| `SingleStarSolveToTOVPoints` | **LIVE — CANONICAL** | ADR-0005; sole ordinary-star radial driver | DO NOT DELETE | — |
| `Solve(Axis)`, `SolveToProfile` | **LIVE — SUPPORTED WRAPPER** | orchestrators over the primitive; six live `Solve` callers | DO NOT DELETE | — |
| `GenTestSequence` | **LIVE — SUPPORTED WRAPPER** (public diagnostic, now covered) | delegates since I4; `tov_gen_test_sequence_cmf`; sole repo reference still commented out | DO NOT DELETE | no |
| `TOVSolver::RadiusLoop` | **DELETED** (I4) | zero references remain except two historical comments | done | — |
| `Analysis` hook (`SetAnalysis`, `analysis->Analyze/Export`) | **SUPPORTED EXTENSION HOOK — CURRENTLY UNUSED BY CORE WORKFLOWS** | owner intent: in-memory per-star analysis during a sequence loop; three project modules derive from `Core::Analysis`; no current `main/` caller attaches one | **DO NOT DELETE** (ADR-0005 Q4; owner intent) | no |
| `n_exp_cond_f` / `virtual ExportNStarProfile` | **PUBLIC — USED BY TESTS ONLY** | only the three 3E harnesses set it; it is their observational capture mechanism | DO NOT DELETE (tests depend on it) | no |
| `NStar::Append` + `FinalizeSurface` | **LIVE — LEGACY** (value-equivalent, geometry-conformant, not consolidated) | Path-1 constructor; both wrappers use it | LATER — only via optional I3 | no |
| `NStar::BuildFromTOV` | **LIVE — CANONICAL** (Path-2 constructor) | all Phase-2 harnesses | DO NOT DELETE | — |
| `NStar` "LEGACY DS PATH" commented blocks (4) | **DEAD CODE** (commented) | `NStar.cpp` | LATER (cosmetic) | no |
| `NStar::BaryonNumIntegrand(double)` | **PUBLIC — CURRENTLY UNUSED, KNOWN DEFECT** | §9 | LATER, after INV-14 adjudication | no |
| `MixedStar`, `Solve_Mixed`, `RadiusLoopMixed` | **COMPILED — UNEXERCISED** (distinct two-fluid physics) | §10 | DO NOT DELETE | no |
| `MixedStar::BaryonNumIntegrand[_Dark]` | **PUBLIC — CURRENTLY UNUSED** (referenced only in commented-out code) | `MixedStar.cpp:31,34,421,424` | LATER | no |
| `DarkCore_Analysis` (`Extensions/MixedStar/`) | **PROJECT-SPECIFIC EXTENSION — UNEXERCISED** (an `Analysis` subclass) | compiled; no caller; inline proper-volume form governed *as a contract* by ADR-0004 §16 | DO NOT DELETE (owner intent) | **no — not core** |
| `BNV_Analysis`, `Decay_Analysis` (`Microphysics/BNV/Analysis/`) | **PROJECT-SPECIFIC EXTENSION — UNEXERCISED** (`Analysis` subclasses) | compiled; no caller | DO NOT DELETE (owner intent) | **no — not core** |
| `BNV_Chi`, `BNV_B_Chi_Photon` (`Microphysics/BNV/Internal`, `/Channels`) | **PROJECT-SPECIFIC EXTENSION — UNEXERCISED** | compiled; no evolution driver connects them | DO NOT DELETE (owner intent) | **no — not core** |
| `RotationSolver` O(Ω²) | **UNREACHABLE SCIENTIFIC CANDIDATE — PROTECTED** | `rot_solver` private, no accessor; `675b4a9` | DO NOT DELETE (§5) | no |
| `RotochemicalCache`, `Rotochemical`, `HeatingFromChem`, `Rates/Urca.hpp` | **NOT COMPILED SCIENTIFIC CANDIDATE — PROTECTED** | absent from CMake lists | DO NOT DELETE (§5) | no |
| `AccretionTorque`, `BNVSpinTorque`, `BNVSource`, `WeakRestoration`, `Coupling` | **DEAD CODE** (0–1 byte files, some in `install(FILES)`) | `CURRENT_ARCHITECTURE.md` component table | LATER | no |
| commented `RotationSolver` fragments | **DEAD CODE** (387 of 1279 lines commented) | `RotationSolver.cpp` | LATER — Phase 4 owns this file | no |
| commented `main/Examples` `Solve` calls (`sig_omega_rho*`, `coulatt`/`polytrope` exports) | **DEAD CODE** (commented) | 3E-0 caller audit | LATER | no |
| `build_xcode/` (123 files) | **GENERATED ARTIFACT** | tracked on `master` with identical count — **inherited, not introduced by Phase 3** | LATER (anticipated ADR "Generated artifacts in version control") | no |
| `docs/doxygen_output/` (2204 files) | **GENERATED ARTIFACT** | tracked on `master`, identical count — inherited | LATER (same ADR) | no |
| `main/Test/results/` (35 files incl. logs) | **GENERATED ARTIFACT** | tracked on `master`, identical count — inherited | LATER | no |
| historical statements in `PHASE2B_CLOSURE.md` (e.g. "five cache hazards remain") | **DOCUMENTATION DEBT** (historical, now scoped) | Phase-2B record; resolved in 3B | scope-marked in this task | no |

Nothing in this table blocks Phase 3. `GOVERNANCE.md` §5 applies to the two *core* candidate
rows (O(Ω²) Hartle; rotochemical) and neither was activated, modified, or ratified by Phase 3.
The project-extension rows are governed by owner intent and, for their inline geometry, by
ADR-0004's contract-only §16.

### 7a. Core library vs project-specific extensions — a boundary issue, recorded not resolved

CompactStar's **core** is the reusable neutron-star machinery: TOV, EOS, geometry, rotation,
thermal, evolution, caches. Around it sit **research-project modules** that inherit from or
consume that machinery to answer a specific paper's question. Today they are interleaved in the
tree:

| Module | Location today | Relationship to core |
|---|---|---|
| `DarkCore_Analysis` | `Extensions/MixedStar/` | `Analysis` subclass over `MixedStar` |
| `BNV_Analysis`, `Decay_Analysis` | `Microphysics/BNV/Analysis/` | `Analysis` subclasses over `NStar` |
| `BNV_Chi`, `BNV_B_Chi_Photon` | `Microphysics/BNV/Internal/`, `/Channels/` | rate/channel code consumed by the above |

Consequences for this closeout: these modules are **not** required to validate or migrate for
core Phase-3 completion; their inline proper-volume forms and architectural debt are **not** core
blockers; and they are **not deleted**. Phase 3 changed none of them.

**Recommendation for a future organisation task (not authorised here):** separate the two
tiers explicitly — e.g. an `extensions/` or `projects/` tree, or external project repositories
that depend on the core library — with the `Analysis` hook as the sanctioned extension seam.
That is a repository-organisation decision for the owner, later; the only obligation this audit
places on Phase 3 is to record the boundary accurately.

## 8. 3E-I3 decision — **NOT JUSTIFIED NOW; DEFERRED — OPTIONAL**

| Criterion | Finding |
|---|---|
| A. active correctness hazard from the duplication? | **No.** Both constructors are value-equivalent (3E-0) and geometry-conformant (I2); Path-1/Path-2 `B` is bitwise at all 17 comparisons |
| B. uncertain ownership under GOVERNANCE? | **No.** ADR-0005 §3 names the STAR-PROFILE CONSTRUCTOR role and Q3 = P3 explicitly leaves both in place; a document *does* say |
| C. current or imminent Phase-4 code depends on choosing one? | **No.** Phase 4 is rotation; `RotationSolver` consumes a finished profile identically from either constructor |
| D. would consolidation simplify the `RotationSolver` work? | **No.** That work touches `RotationSolver.cpp`, not profile construction |
| E. would it remove a known bug likely to contaminate Phase 4? | **No.** The only asymmetry — Path-1 mirror `M`/`R`/`z_surf` at zero — is confined to `n_star` inside `Solve`/`GenTestSequence`, which is `Reset()` immediately; no consumer reads those mirror fields on a Path-1 star |
| F. requires touching hooks / public workflow contracts? | **Yes.** P2 must re-centre the `Analysis`/`ExportNStarProfile` hooks (which act on the solver's `n_star`) and preserve `Sequence` — exactly the API surgery ADR-0005 §9 flagged |
| G. detector/validation need not otherwise satisfiable? | **No.** `tov_path_equivalence_cmf` and `tov_gen_test_sequence_cmf` already assert both constructors' outputs |

No criterion fires, F argues against, and *"textual duplication alone is not sufficient reason
for another structural migration."* → **§7 classification: I3 NOT JUSTIFIED — SKIP** for this
merge.

Why **DEFERRED — OPTIONAL** rather than **SHOULD NOT BE DONE**: the mirror-scalar asymmetry is a
real latent inconsistency (classified `INTERNAL STATE ASYMMETRY — CURRENTLY UNOBSERVED`, never
declared correct), and ADR-0005 §10 already fixes the governed place to correct it — **M2, inside
P2, under coverage**. The explicit **re-trigger condition** is: *a consumer that reads
`MassSurface()`, `RadiusSurface()` or `ExpNuSurface()` on a Path-1-constructed star appears.*
That would make criterion E true. Until then, I3 is not worth its API surgery. **Not implemented
here.**

## 9. INV-14 scalar accessor — **KNOWN DEFECT — UNEXERCISED — DEFERRED**

Re-authenticated at `NStar.cpp:1085-1100`: the public `NStar::BaryonNumIntegrand(double)`
computes `4π r² n_B / sqrt(1 − 2M/r)` — the same geometric shape as the integrand — **without
the `1e54` conversion**, and returns `0` for `f ≤ 0` (the ADR-0004 "B4" behavior, itself
nonconformant with the accepted domain contract). Declared in `NStar.hpp:393`; **zero callers**
in `CompactStar/`, `main/` or `tests/`. `MixedStar::BaryonNumIntegrand[_Dark]` are referenced
only from commented-out code.

- Phase-3 blocker? **No** — never executed on any path.
- Phase-4 blocker? **No** — Phase 4 does not integrate baryon number through this accessor.
- Repair now? **No.** Repairing it inside a closeout would silently convert a public function's
  output units with no covering test, which is the pattern every Phase-3 increment declined.
  It needs a small INV-14 adjudication (fix, deprecate, or delete) with a test.

## 10. `MixedStar` / two-fluid boundary

Re-authenticated: `MixedStar.cpp` **is compiled** (`Core/CMakeLists.txt:22`); the only `main/`
reference is a predicate whose registration is **commented out**
(`tov_debug_main.cpp:51,83`); **zero test coverage**; two-sector mass semantics
(`mass_tot_dc`) with six inline proper-volume sites (ADR-0004 §15); `RadiusLoopMixed` integrates
a **core + mantle two-fluid system with separate GSL drivers** — a genuinely different algorithm,
not a copy of the ordinary-star loop.

**Nothing about `MixedStar` prevents declaring the ordinary visible-sector modernization
complete.** It is a separate, unvalidated legacy physics scope. Its future work is already
governed: ADR-0004 §0-Q2 requires **focused coverage first**, then migration to the same
`Geometry` primitive with total enclosed mass; the ADR index carries an anticipated
"MixedStar modernization scope" ADR. Sharing the `TOVSolver` class does not make it Phase-3
scope.

## 11. Candidate scientific code

| Candidate | Compiled? | Reachable? | Validated? | Governed? | Target |
|---|---|---|---|---|---|
| O(Ω²) Hartle (`RotationSolver`) | yes | **no** (private, no accessor) | no | §5, INV-08 | Phase 4 |
| `RotochemicalCache` / `Rotochemical` | **no** | no | no | §5, INV-01 nonconformance recorded | Phase 5 |

**Project-specific extensions (owner clarification — not core candidates; listed for
completeness, not as closure prerequisites):**

| Module | Compiled? | Reachable? | Validated? | Governed? | Target |
|---|---|---|---|---|---|
| `BNV_Analysis`, `Decay_Analysis`, `BNV_Chi`, `BNV_B_Chi_Photon` | yes | no | no | owner intent (§7a); ADR-0004 §4.3/§16 contract only for their inline measures | project-scoped; future organisation task |
| `DarkCore_Analysis` | yes | no | no | owner intent (§7a); ADR-0004 §16 | project-scoped; depends on `MixedStar` |

**Phase 3 did not activate, modify, or ratify any of them — core candidates or project
extensions.** The only Phase-3 touch was *documentary*: ADR-0004 recorded a future proper-volume
contract for the inline sites. Scientific meaning is unchanged throughout.

## 12. Units / constants closeout

| Item | State | Disposition |
|---|---|---|
| exact-duplicate domain conversions (`KM3_TO_CM3`, MeV fm⁻³→erg cm⁻³) | one owner, `Units.hpp` (3A) | ✅ |
| `k_B` | one authority, `Zaki::Physics::K_BOLTZ_EV` (3C); zero literals, zero GSL Boltzmann on production paths | ✅ |
| **solar mass in km** | **two authorities**: `Zaki::Physics::SUN_M_KM` (9 files) vs `GSL_CONST_CGSM_SOLAR_MASS` (2 files), disagreeing at `6.2e-5` | **governed / deferred debt** — the ratified roadmap defers it "out of Phase 3 pending owner or ADR adjudication"; **not reopened here**; not a Phase-3 blocker |
| other local compact-object / geometric conversions (`INV_FM4_2_*`, pressure/energy geometric) | local, not bit-identical duplicates | deferred; owner intent: these belong in CompactStar, generic constants in ZakiLib |

Owner intent is respected as recorded: generic constants → ZakiLib; domain conversions →
CompactStar; solar mass → unadjudicated.

## 13. Cache closeout (ADR-0003)

Re-authenticated implementation: `ProfileProvenance{const StarProfile*, uint64_t}`;
`GeometryCache` carries and exposes it (immutable snapshot, no `Refresh()`);
`ProfileVersionedCache` keys on `(identity, version)`; `C_⋆` key includes geometry provenance and
**fails closed** on a foreign geometry; the `NeutrinoCooling` payload keys on identity + geometry;
`StarContext` re-binds column views before invalidating and advances its revision last.

**INV-12 scope is precise: RESOLVED for profile-derived caches.** It is not a universal claim
about every cache in CompactStar (e.g. EOS spline accelerators are not profile-derived and were
never in scope).

**No cache regression from 3D/3E:** `geometry_cache_measure_contract` G5 confirms provenance is
undisturbed by the Λ delegation; `cache_contract` and `cache_thermal_contract` pass at this HEAD;
3E touched no cache code.

## 14. Proper-volume closeout (ADR-0004)

Ordinary visible-sector conformance is complete: `BuildFromTOV` (3D), `Append` and
`FinalizeSurface` (I2), `GeometryCache::DeriveLambdaFromMR_` (3D). Path-1 and Path-2 `B` are
**bitwise identical**. Deferred, governed, and recorded in INV-04: `MixedStar` (§0-Q2), the §5
candidates (§16), the INV-14 scalar accessor. `EvaluateNu`'s boundary-condition `x = 1e-15` is
**not this measure** (§4.4).

**Correct INV-04 wording (already in place since I2, retained):**
`GOVERNED (ADR-0004 ACCEPTED) — BOTH ORDINARY VISIBLE-SECTOR NStar PATHS CONFORMED; MixedStar /
CANDIDATE / SCALAR-ACCESSOR MIGRATIONS DEFERRED`. Not "RESOLVED" — governed-but-nonconformant
sites remain — and not a blocker: unvalidated candidate code does not gate ordinary-star closure,
and governance nowhere says it should.

## 15. TOV closeout (ADR-0005)

Structural verification at this HEAD:

| Check | Result |
|---|---|
| ordinary visible-sector radial GSL drivers | **1** — `SingleStarSolveToTOVPoints` (`TOVSolver.cpp:2410`) |
| ordinary step_scale ladders | **1** (`:2457`) |
| `RadiusLoop` references (excl. `RadiusLoopMixed`) | **0 in code**; two historical comments |
| orchestrators delegating | `Solve` ✓ `SolveToProfile` ✓ `GenTestSequence` ✓ |
| workflow contracts | `_Sequence.tsv` guarded (`tov_sequence_workflow_cmf`); `_TestSequence.tsv` guarded (`tov_gen_test_sequence_cmf`) |
| six `Solve()` callers | unchanged since Phase-3 entry |
| `RadiusLoopMixed` drivers (`:2801/:2807`) | distinct two-fluid problem; excluded by design |

**Fail-closed #3: CLOSED for ordinary visible-sector TOV radial numerical ownership.** Not
claimed for `MixedStar`, dark-sector TOV, or future solvers.

## 16. Test-suite quality (19 authenticated / 10 self-contained)

| # | Test | Purpose | Data | Assertions | Detector-backed |
|---|---|---|---|---|---|
| 1 | `compactstar_library_smoke` | build/link smoke | self | exit code (intentional) | — |
| 2 | `heat_capacity_v1` | analytic/reference + contract (ADR-0002 V1) | self | yes; **independent long-double SI `k_B` oracle** | 3C (1-ULP agreement) |
| 3 | `tov_reference_analytic` | analytic/reference (Schwarzschild interior) | self | yes | — |
| 4 | `cache_contract` | contract (ADR-0003 P1–P4) | self | yes | 3B hazards stopped reproducing |
| 5 | `cache_thermal_contract` | contract | self | yes | 3B |
| 6 | `hartle_moment_inertia_analytic` | analytic/reference | self | yes | 2B-4B |
| 7 | `evolution_stepper_contract` | contract | self | yes | — |
| 8 | `photon_cooling_conformance` | conformance | self | yes; independent oracle | 3C |
| 9 | `proper_volume_contract` | contract (ADR-0004 domain + identities) | self | yes | 3D D1/D2b/D3/D5 |
| 10 | `geometry_cache_measure_contract` | contract (cache ≡ primitive; composition) | self | yes | 3D D1/D2/D6 |
| 11 | `heat_capacity_real_star` | **diagnostic harness** — ADR-0002 Tier-B evidence emitter | ext | **none** (fails only on load/solve) | — |
| 12 | `passive_cooling_regression` | regression vs golden | ext | yes (its own `REGRESSION` checks; `1e-5`/`1e-4`) | 3C, 3D D2 |
| 13 | `tov_reference_cmf` | reference (CompOSE `eos.mr`) | ext | yes | 2B-2 |
| 14 | `grid_convergence_cmf` | convergence | ext | yes | 2B-4A |
| 15 | `hartle_moment_inertia_cmf` | scale-free reference | ext | yes | 2B-4B |
| 16 | `baryon_number_cmf` | reference/regression (predeclared 1e-15) | ext | yes | 3D D1/D2b/D3 |
| 17 | `tov_path_equivalence_cmf` | equivalence (two orchestrators over one primitive) | ext | yes | 3E-0 D1/D2; I2 D4 |
| 18 | `tov_sequence_workflow_cmf` | interface/workflow | ext | yes | I1 D2/D3a/D3b |
| 19 | `tov_gen_test_sequence_cmf` | equivalence + interface | ext | yes | I4 D1 (pre-migration) |

Findings:

- **No test imports a production constant as its own oracle.** The three thermal tests derive
  `k_B` independently in `long double` from the SI-exact definitions; `geometry_cache_measure_contract`
  calls the `Geometry` primitive *by design* to pin cache ≡ primitive, while
  `proper_volume_contract` separately validates the primitive against two independent forms.
- **One diagnostic masquerades as a test: `heat_capacity_real_star`** has no assertions. This
  was found during 3D's detector work and is recorded there. It is intentionally an
  evidence-emitting harness (ADR-0002 Tier-B) and should be **classified, not deleted**; adding
  assertions is Phase-2A-adjacent coverage debt, not a Phase-3 blocker.
- **No redundancy.** `tov_path_equivalence_cmf` and `tov_gen_test_sequence_cmf` cover different
  orchestrators and different file contracts. After I4, the former's documented role as a
  "duplicate radial tolerance detector" is obsolete by construction (recorded in I1); its other
  roles — finalization, species, `B`, `I`, leakage, resolution — remain live.
- **External-data guard is correct**: all 9 external tests sit inside the
  `if(COMPACTSTAR_EOS_DATA_ROOT …)` block; 10 self-contained tests run in a clean checkout.
- **Sufficient to protect the merge?** Yes: every Phase-3 structural change has either a
  contract test with a demonstrated detector or a structural-impossibility proof (§17).

## 17. Detector / falsifiability audit

| Claimed guardrail | Falsifiable by | Demonstrated |
|---|---|---|
| exact-duplicate constants centralized bit-identically (3A) | artifact hash identity | five hashes unchanged |
| `k_B` authority (3C) | `1e-4` perturbation → `passive_cooling_regression` fails | fired; reverted |
| cache identity / provenance (3B) | five hazard reproductions stop; `cache_contract` P1–P4 | reproductions ceased |
| geometry factor & domain (3D) | D1 denominator, D2/D2b coordinate volume, D3 unit leak, D5 clamp, D6 composition | all fire; D4 correctly not executed |
| TOV path equivalence (3E-0) | D1 Path-1-only tolerance; D2 Path-1 metric factor | both fire, D2 B-only |
| workflow file contract (I1) | D2 row drop; D3a header; D3b filename | all fire |
| Path-1 geometry conformance (I2) | D4 proper-volume removal | fires, 8 B-only assertions, radial columns green |
| duplicate radial loop eliminated (I4) | D1 fires **before**; **structurally impossible after** | one-owner proof |
| mirror-scalar asymmetry is *asserted*, not assumed | `tov_path_equivalence_cmf` asserts current zeros | yes |
| no cross-star state leakage | identical-Axis-node bitwise assertion | yes |

**Phase-3 claims lacking falsifiability: none.** The one item without a *mutation* detector —
the Hartle-`I` side effect (ADR-0005 planned D5) — is protected by bitwise `I` equality
assertions in two tests, which any removal from one path would fail. The suite-level gap is
`heat_capacity_real_star` (§16), which is a Phase-2A evidence emitter rather than a Phase-3
claim.

## 18. Golden / baseline audit

| Artifact | Kind | Producer | Tolerance | Phase-3 change | Justification |
|---|---|---|---|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | **regression golden** (current output) | `passive_cooling_regression --emit` | `1e-5` / `1e-4` | **re-baselined in 3C** `edaa5e3b…` → `831744b0…` | `k_B` authority change; static response `1.7e-11`; trajectory `1.3e-6` shown to be step placement; owner adjudicated; tolerances **not** widened |
| `grid_convergence_cmf_1p6_debug.tsv` | static-diagnostic golden | `grid_convergence_cmf --emit-dir` | byte / predeclared | **re-baselined in 3C** `3be2005f…` → `61d84ddc…` | only `Cstar_1e8`, `dlnT_dt_1e8` moved, `1.7e-11`; 14 columns byte-identical |
| `grid_convergence_cmf_1p6_trajectory.tsv` | trajectory golden | same | `1e-5` | **re-baselined in 3C** `e04f7485…` → `ca32863d…` | same authority change |
| `tov_dscmf1_reference.tsv` | reference golden | `tov_reference_cmf --emit` | as 2B-2 | **unchanged** all phase | k_B-independent; a falsifiable prediction that held in 3C and 3D |
| `hartle_I_dscmf1_debug.tsv` | reference golden | `hartle_moment_inertia_cmf --emit` | as 2B-4B | **unchanged** all phase | same |
| `baryon_number_dscmf1_reference.tsv` | **canonical value table** (new, 3D) | `baryon_number_cmf --emit` | `1.0e-15` predeclared (ADR-0004 §7.1) | created 3D; unchanged since | B had no artifact before; pre-3D values in its header |
| `tov_path_equivalence_dscmf1.tsv` | **measurement record — NOT a compared golden** (new, 3E-0) | `tov_path_equivalence_cmf --emit` | n/a (CTest runs without `--compare`) | created 3E-0; **re-emitted in I2** `a6b4cd79…` → `bbf61e5f…` (only `rel_B` moved) | records path-to-path differences; leaving it asserting a closed gap would be false |

**No unjustified tolerance widening anywhere in Phase 3.** The 3C predeclared `1.7e-11` was left
unwidened and annotated; ADR-0004's `1.0e-15` was never widened and was not even spent in I1; the
regression tolerances are unchanged from Phase 2B.

## 19. Documentation consistency — findings and corrections made in this task

| Location | Finding | Action |
|---|---|---|
| `MODERNIZATION_ROADMAP.md` (Phase-3 status paragraph) | "**Phase 3 is not complete.**" | **corrected** (current-state) |
| `MODERNIZATION_ROADMAP.md` (entry-state blockquote) | "two fail-closed conditions are already active … two live TOV paths … five INV-12 hazards" | **scope-marked** as Phase-3-entry state; both since discharged |
| `PHASE3_CONSOLIDATION_PLAN.md` header note | "3E is governance only so far … ADR-0005 PROPOSED … no implementation" | **corrected** |
| `PHASE3_CONSOLIDATION_PLAN.md` §12 tail | "Phase 3 is not complete … no canonical TOV owner … That ADR is now drafted … PROPOSED" | **corrected** |
| `PHASE3_CONSOLIDATION_PLAN.md` 3F rows | predeclared "full 13/13" and "document only" | **marked executed** (documentation-only, no deletion) |
| `CURRENT_ARCHITECTURE.md` tests paragraph | "13 registered CTests … 13/13 / 8/8" | **corrected** to 19 / 10 / 9-external |
| `TARGET_ARCHITECTURE.md` gap table | "Single uniform cache-invalidation rule — ABSENT"; "`GeometryCache` version-gated — ABSENT" | **corrected** to IMPLEMENTED (3B) |
| `docs/adr/README.md` anticipated table | "MixedStar modernization scope — Phase 3"; solar-mass authority ADR not listed | **re-targeted** post-Phase-3; **solar-mass ADR added** as anticipated |
| `PHASE2B_CLOSURE.md` §12 | "Cache. … five hazards remain" | **scope-marked** historical (resolved in 3B); text otherwise preserved |
| `SCIENTIFIC_INVARIANTS.md` INV-03, INV-04 | "the §5-protected candidate code" used for BNV/Decay/DarkCore — a mislabel introduced in 3D/I2 | **corrected** to distinguish §5 core candidates from project-specific extensions (owner clarification) |
| ADR-0004 §24 / ADR-0005 §23 implementation records | provisional hook classification "LEGACY PUBLIC HOOK — CURRENTLY UNUSED" (ADR-0005 §11 itself deferred classification to 3F) | **3F note appended** recording the owner's clarification; **no accepted decision text touched** |
| `PHASE3D_PROPER_VOLUME.md`, `PHASE3E_I2_PATH1_GEOMETRY.md` | same "§5 candidate" phrasing | **left as historical**; superseded by this record |
| ADR-0003/4/5 records; Phase-2 validation docs; `PLAN` §1 "13/13 at entry"; roadmap Phase-2B "13/13" | historical statements about the state *at that time* | **left as historical** — accurate for their date |

## 20. Phase-2B parallel debt — not laundered into Phase 3

`PHASE2B_CLOSURE.md` is explicit and this audit re-authenticates it:

| Debt | Phase-2B status | Blocks Phase-3 merge? | Basis |
|---|---|---|---|
| **CI absent** | `CI IS A PHASE-2B COMPLETION BLOCKER` | **No** | roadmap keys Phase 3 off "baselines exist", which is satisfied; governance does not make CI a Phase-3 criterion; making it one now would be inventing an exit criterion |
| full coupled TOV→Hartle→cooling convergence | `OWNER ADJUDICATION REQUIRED` — a roadmap dependency inconsistency (needs Phase-4/5 physics) | **No** | not a Phase-3 item |
| continuum accuracy adequacy | unresolved; no governed accuracy budget; regression tolerance is not a continuum claim | **No** | recorded 2B limitation |

**Phase 2B remains formally NOT COMPLETE.** Phase 3 closing does not change that and does not
claim otherwise.

## 21. Phase-4 entry state (assuming Phase 3 closes)

| Prerequisite | State at entry |
|---|---|
| INV-07 physical Ω/J normalization | **UNRESOLVED** — Phase 4's central item |
| O(Ω²) Hartle candidate | **unvalidated, unratified** (§5), unreachable |
| scale-free `I = J/Ω` | **verified** (2B-4B), bit-identical across both constructors |
| TOV radial primitive | **canonical, sole owner** |
| geometry mathematics | **canonical** on all ordinary-star paths |
| caches | **governed** (ADR-0003) |
| thermal | **baseline-protected** (three goldens, tolerances unchanged) |
| INV-05 rotation-half re-audit | pending (Phase 4) |

**No remaining Phase-3 item blocks Phase 4.** Phase 4 is not begun here.

## 22. Merge mechanics

| Check | Result |
|---|---|
| ancestor relationship | `master` **is** the merge-base → fast-forward possible |
| behind count | 0 |
| conflicts with `master` | **none possible** (fast-forward) |
| working tree | clean before this audit; documentation-only diff after |
| upstream | equal |
| tracked scratch / test outputs introduced by Phase 3 | **none** — the tracked `build_xcode/`, `docs/doxygen_output/`, `main/Test/results/` all pre-exist on `master` with identical file counts |
| untracked scientific data committed | **none** (CompOSE datasets are external) |
| dependency binaries / external trees changed | **none** in `master..HEAD` |
| detector mutations left behind | **none** — every increment verified its revert by SHA-256 |
| temporary backups | none tracked |
| absolute local paths | **none in source or tests**; two evidence records (`PHASE3_CONSOLIDATION_PLAN.md`, `PHASE3C_BOLTZMANN_AUTHORITY.md`) cite the ZakiLib inspection path as provenance — evidence, not code |
| generated artifacts newly tracked | **none** |
| baseline changes lacking provenance | **none** — every change carries old/new SHA-256 and a reason |

**Classification: `PHASE-3 BRANCH MERGE READY`.** No blockers.

## 23. Phase-3 exit-criterion evaluation (read literally)

The ratified roadmap states exactly one exit criterion: **"One authoritative owner per quantity;
baselines still pass."**

| Clause | Result | Basis |
|---|---|---|
| one authoritative owner — TOV radial solve | **SATISFIED** | §15 |
| one authoritative owner — proper-volume measure | **SATISFIED** in ordinary visible-sector scope; remainder **deferred by accepted ADR-0004 §0-Q2 / §5** | §14 |
| one authoritative owner — cache invalidation | **SATISFIED** (profile-derived scope) | §13 |
| one authoritative owner — `k_B`, exact-duplicate conversions | **SATISFIED** | §12 |
| one authoritative owner — solar mass | **NOT satisfied — deferred by the ratified roadmap's own text**, which excludes it from Phase 3 pending adjudication | §12 |
| baselines still pass | **SATISFIED** — 19/19, 10/10; all seven artifacts byte-identical at this HEAD | §1 |

The criterion is **SATISFIED in its governed reading**: every quantity Phase 3 was ratified to
consolidate has one owner, and the one quantity without one (solar mass) was excluded from Phase
3 by the same ratified document that states the criterion. The roadmap contains no further
exit conditions; none was invented here. The "open for separate consideration" thermal-balance
item is explicitly not a prerequisite.

## 24. Formal Phase-3 status

> **`PHASE 3 COMPLETE — MERGE READY`**

Tied to §23 (ratified criterion satisfied), §6 (no active fail-closed condition attaches to
Phase-3 scope), §22 (no merge blocker), and §16–§18 (suite and baselines sufficient and
untampered). "Complete" is claimed on the ratified criterion, not on tests passing.

## 25. Known deferred debt carried across the merge (recorded, not smoothed)

1. **Solar-mass authority** — dual, unadjudicated; anticipated ADR added to the index.
2. **`MixedStar` / two-fluid** — compiled, unexercised, unvalidated, ADR-0004-governed; needs
   coverage before migration.
3. **§5 candidates** — O(Ω²) Hartle (Phase 4), rotochemical (Phase 5), BNV, DarkCore.
4. **INV-14 scalar accessor** — known defect, zero callers, needs adjudication + test.
5. **3E-I3** — optional; re-trigger condition in §8.
6. **Mirror surface scalars** on Path-1 stars — asserted zero; governed by ADR-0005 §10 (M2).
7. **`Analysis` extension hook** — supported, currently unused by core workflows; preserved by
   ADR-0005 Q4 and owner intent; its modern API may be reconsidered later, but the in-memory
   per-star sequence-analysis capability is to be kept. `n_exp_cond_f`/`ExportNStarProfile`
   remain the test harnesses' observational capture mechanism.
7a. **Core-vs-project boundary** — BNV / Decay / DarkCore modules are project-specific
   extensions interleaved in the core tree (§7a); a future organisation task should separate the
   tiers. No migration authorised in Phase 3.
8. **`heat_capacity_real_star` has no assertions** — evidence emitter; coverage debt.
9. **Generated artifacts tracked on `master`** (`build_xcode/`, `docs/doxygen_output/`,
   `main/Test/results/`) — inherited; anticipated ADR.
10. **Phase-2B**: CI absent; coupled-convergence adjudication; continuum accuracy — unchanged.

## 26. Exact next action

> Owner independently reviews this closeout record, then performs a controlled
> fast-forward/merge of `modernization/behavior-preserving-consolidation` into `master`,
> followed by post-merge verification (19/19, 10/10, seven artifact hashes byte-identical at the
> new `master`).

The merge was **not** performed by this task.
