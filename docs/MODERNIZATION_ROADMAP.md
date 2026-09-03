# CompactStar Modernization Roadmap

> **TOV-SURF-I-R (2026-09-03):** ADR-0009 source conformed and numerically validated.
> V7a/V7b clarification and the prior stop are retained. All seven artifacts remain
> unchanged; their dedicated migration is next, before corrected Phase-4D revalidation
> and the first monopole baseline. No Phase 5 work. Evidence: `docs/validation/TOV_SURFACE_IMPLEMENTATION.md:339`.

> **STATUS: RATIFIED — 2026-08-31.** The phase ordering and prerequisites are accepted
> (`GOVERNANCE.md` §7). Ordered by **dependency, not by date.** A phase begins only when its
> prerequisites hold.
>
> Governing principle: **do not propose scientific implementation before its prerequisites are
> governed and validated.** "Validated" means *evidence adequate to the change*, which is normally
> a regression baseline — and, in the single case `GOVERNANCE.md` §3.1 authorizes, independent
> physical verification standing in for a baseline that cannot legitimately exist yet (Phase 2A).

---

## Phase 0.5 — Source-of-truth reconciliation and governance foundation

**Status: ✅ COMPLETE — 2026-08-31.** All exit criteria satisfied; the governance package was
ratified by the owner at commit `617bb0e` (`GOVERNANCE.md` §7).

- Reconciled the source of truth. Owner commit `3639d71` ("updates", 2026-04-07) reworked
  `RotationSolver` and `MixedStar` and is **not** an ancestor of the previously audited
  `d91c31b`. Authenticated base: **`9f70f14`** (`master`), which contains both `3639d71` and the
  April candidate work `675b4a9`.
- Persisted the Phase-0 reconnaissance as a dated, non-normative snapshot with an explicit
  supersession appendix.
- Established `GOVERNANCE.md`, `docs/SCIENTIFIC_INVARIANTS.md`, the ADR system, the
  current/target architecture split, `AGENTS.md`, and this roadmap.
- Raised ADR-0001 (species profile semantics) as PROPOSED.
- **ADR-0001 ACCEPTED 2026-08-31** by owner adjudication: `n_B` in fm⁻³ in the
  `BaryonDensity` column; species columns are dimensionless fractions `Y_i = n_i/n_B`;
  `n_i = Y_i n_B` derived at the point of use; **no normalization on import**. INV-01 moves to
  GOVERNED (ACCEPTED).
- **ADR-0002 ACCEPTED 2026-08-31** by owner adjudication: one canonical physical heat capacity
  `C_⋆(T∞)` — the GR-integrated EOS/CompOSE-based stellar heat capacity, designated to
  `StarContext::HeatCapacityStar_Tinf` — for the evolved thermal degree of freedom; every energy
  channel divides by it; a driver-local constant is not acceptable as a production denominator.
  INV-15 moves to GOVERNED (ACCEPTED). **The software question of where the division occurs is
  explicitly deferred** (ADR-0002 §6). **No source change is authorized**; conformance is Phase 2A.
- **Repaired the Phase-2 / Phase-3 circular dependency** that ADR-0002 exposed. Phase 2 is split
  into **2A** (pre-baseline correctness prerequisites) and **2B** (the validation baseline);
  Phase 3 loses the heat-capacity ownership item.
- **Governance coherence audit before ratification.** Resolved two defects that ADR-0002 created
  or exposed. (a) **Authority inversion:** the Phase-2A pre-baseline exception existed only in the
  rank-8 skill plan while rank-1 `GOVERNANCE.md` §3 still made it fail-closed, so the accepted
  ADR-0002 workflow was procedurally invalid. The exception now lives once, at rank 1, as
  `GOVERNANCE.md` **§3.1**, with seven binding conditions; `AGENTS.md` and `AI_SKILL_PLAN.md`
  reference it rather than redefining it. (b) **Plumbing before mutation:** Phase 2A required
  durable scientific verification while the test mechanism was not introduced until Phase 2B.
  Minimal CTest plumbing moves to **Phase 1**. Also corrected the categorical convergence-order
  claim (INV-13) and stated the scope of invariant-register ratification.

**Exit criteria.**

| Criterion | Status |
|---|---|
| ADR-0001 adjudicated | ✅ **SATISFIED** — ACCEPTED 2026-08-31 |
| ADR-0002 adjudicated | ✅ **SATISFIED** — ACCEPTED 2026-08-31 |
| Governance package internally coherent | ✅ **SATISFIED** — coherence audit complete; no known contradiction between authority ranks |
| Owner reviews the remaining governance documents | ✅ **SATISFIED** — `GOVERNANCE.md`, `SCIENTIFIC_INVARIANTS.md`, `AI_SKILL_PLAN.md`, and this roadmap **RATIFIED 2026-08-31** at commit `617bb0e`, subject to the recorded status distinctions (`GOVERNANCE.md` §7) |

The two ADR acceptances and the package ratification are **separate decisions**. Accepting
ADR-0001 settled the species-semantics contract; accepting ADR-0002 settled heat-capacity
ownership; the 2026-08-31 ratification accepted the governance package as a set. **None of the
three** ratifies the Hartle O(Ω²) / rotochemical candidate code, validates any
`INTENDED BUT UNVERIFIED` item, resolves INV-07 or INV-11, excuses the `PhotonCooling` or
`RotochemicalCache` nonconformance, or certifies the numerical correctness of
`StarContext::HeatCapacityStar_Tinf`. `GOVERNANCE.md` §7 records the full non-scope.

---

## Phase 1 — Reproducible macOS build **and minimal validation plumbing**

**Status: ✅ COMPLETE.** Delivered in four increments — 1A (clean-checkout configure), 1B
(`RotationSolver` merge repair, library builds), 1C (CTest plumbing and a link smoke test), and
1D (default build type and warning policy). Every item below is satisfied and both exit-criteria
clauses are demonstrated. Evidence: `docs/build/MACOS_BUILD.md`.

**Prerequisite:** Phase 0.5 reviewed. ✅ **SATISFIED** — governance ratified 2026-08-31.

Nothing below Phase 1 could be verified, because at the start of this phase the project could not
be configured from a clean clone on any platform. It now can, on macOS.

**Build reproducibility**

- ✅ Guard the seven absent optional `main/` subdirectories so a clean clone configures. *(1A)*
- ✅ Stop writing generated configuration into the source tree — `configure_file` output moved to
  the binary directory. *(1A)*
- ✅ Define and document canonical macOS configure/build commands. *(1A)*
- ✅ Record exact dependency versions: GSL, OpenMP, Python3/NumPy, and the vendored Zaki and
  Confind archives. *(1A)*
- ✅ Adopt a warning policy and a default build type. *(Engineering class — must not change results.)* *(1D)* `Debug` is the default only when the user supplies none; `-Wall -Wextra` are `PRIVATE` to the `CompactStar` target; dependency headers are `SYSTEM`. The resulting 175 warnings are **inventoried, not repaired** — several touch scientific source and need their own change class. `-Werror` is deliberately not enabled.

**Minimal validation plumbing — generic infrastructure only**

At the start of Phase 1 the project had **no test mechanism of any kind**: `enable_testing`,
`include(CTest)`, and `add_test` appeared in no `CMakeLists.txt` anywhere in the tree. Phase 1
establishes only the generic ability to *run* a test — delivered by increment 1C:

- ✅ `enable_testing()` / CTest plumbing — `include(CTest)` in the top-level `CMakeLists.txt`, gated on standard `BUILD_TESTING`. *(1C)*
- ✅ A canonical location for test executables and fixtures — `tests/`, distinct from the manual demo programs in `main/Test/`. *(1C)*
- ✅ A trivial smoke test that runs through the standard test command — `compactstar_library_smoke`, which resolves a real out-of-line symbol from `libCompactStar.a`. *(1C)*
- ✅ The canonical test command, documented alongside the configure/build commands. *(1C)*

**No scientific baseline and no physics check is created here** — only the mechanism. This is
dependency/build class and must not change any number the code produces.

**Why this belongs in Phase 1 rather than Phase 2A.** Three reasons, all from repository evidence.
(a) `enable_testing()` is a top-level directive and must live in the same root `CMakeLists.txt`
that Phase 1 is already opening — alongside `add_subdirectory` at `:92` and `:153`, and the
`configure_file` relocation at `:84-87`. (b) It is dependency/build class, which is Phase 1's
class; Phase 2A is scientific-semantic. (c) Decisively, `GOVERNANCE.md` §3.1 condition 5 requires
Phase-2A work to be **narrowly scoped to the defect ADR-0002 governs**. Standing up a test
framework inside that phase would violate the scoping condition that authorizes the phase at all.

**Exit criteria.** A clean clone configures and builds on macOS with documented commands, **and**
a trivial test runs through the documented test command.

**Exit-criteria status. ✅ SATISFIED.** A clean checkout configures and the `CompactStar` library
builds and links on macOS with documented commands, and `ctest` runs one passing automated test.
The full phase item list is satisfied as well, including the warning policy and default build type.

**What Phase 1 does not establish.** It produced **no scientific validation whatsoever.** The
smoke test proves only that the library links and runs; the 175 warnings are recorded, not fixed;
and `Debug` versus `Release` is **not** claimed to be numerically equivalent — any future baseline
must state the configuration it was produced under.

**Decision gate — platform.** **Mac-first development is acceptable initially.** Cross-platform
dependency work moves earlier **only if** Linux or cloud compilation becomes required, CI must
run on Linux, or clean macOS dependency reproduction proves impossible. **A full Zaki or
Confind modernization is explicitly out of scope for now.**

---

## Phase 2A — Pre-baseline correctness prerequisites

**Status: ✅ COMPLETE.** Delivered in three increments — 2A-1 (Tier-A `C_⋆` verification),
2A-2 (authenticated real-star Tier-B verification, **V1 VERIFIED**), and 2A-3 (`PhotonCooling`
conformance plus the adjacent null-check repair). Evidence:
`docs/validation/HEAT_CAPACITY_V1.md`.

**Prerequisite:** Phase 1 complete — the project builds reproducibly **and a trivial test runs
through the documented test command.** ✅ **SATISFIED** — Phase 1 completed in increments 1A–1D.

**Authorized by `GOVERNANCE.md` §3.1** under **ADR-0002**. This is the first and currently only
invocation of the pre-baseline correctness exception. Every §3.1 condition must be satisfied and
reported, including its four report items.

**Durability requirement.** The independent verification below must be added as **tests running
under the Phase-1 plumbing**, not as one-off manual checks. §3.1 condition 7 requires a regression
baseline immediately afterward, and condition 4's independent evidence must remain re-runnable
after the correction lands — a check that cannot be repeated cannot support the exception.

**This phase is deliberately narrow.** It exists for one reason: a small number of known defects
would make a validation baseline *scientifically misleading* if frozen into it. Those, and only
those, are corrected first.

**This is not a general cleanup phase.** Behavior-preserving consolidation stays in Phase 3;
unrelated scientific corrections stay in their own phases. Admission requires showing that a
baseline captured *without* the correction would encode physics known in advance to be wrong.

- ✅ **Replace `PhotonCooling`'s production use of the constant `C_eff` with the governed
  `C_⋆(T∞)`** (ADR-0002; INV-15). *(2A-3)* The driver now divides by
  `StarContext::HeatCapacityStar_Tinf`; `Options::C_eff` was **removed**, not deprecated; the
  driver Doxygen states the governed equation; and both manual call-site assignments are gone.
  ADR-0002 **Pattern A** preserved — no central thermal-balance owner introduced.
- ✅ **Correct the immediately adjacent safety defect required to exercise that path** *(2A-3)* —
  the `NeutrinoCooling_Details.cpp` null-check ordering issue. The `if (!ctx.star)` guard now
  precedes the first dereference. Engineering class; no neutrino microphysics changed.
- ✅ **Add narrowly targeted verification of the stellar heat-capacity calculation itself**
  *(2A-1, 2A-2)* — ADR-0002 §V1, as durable CTest under the Phase-1 plumbing. All seven items
  addressed: exact unit chain, measured radial order **2.000**, `NT=160` interpolation error
  `6.6e-4`, endpoint clamping characterized, cache rebuilds verified, and on an authenticated
  `1.4 M☉` CMF star `C_⋆(10⁸ K) = 2.17e38 erg K⁻¹`. **V1 VERIFIED.**

**Evidence standard.** Phase-2A changes are governed by their own scientific/numerical class
under `GOVERNANCE.md` §2 and **require independent physical checks** — analytic limits, dimensional
analysis, convergence, and comparison against published values. They **must not** be validated by
comparison against the existing passive-cooling curve, which encodes the behavior ADR-0002
rejects. That is precisely the circularity this split removes.

**Exit criteria. ✅ SATISFIED.** The thermal energy equation now uses one heat capacity for both
channels; `C_⋆(T∞)` has passed independent verification (V1 VERIFIED); the adjacent null-check
defect is repaired. **`GOVERNANCE.md` §3.1 condition 7 is now outstanding**: the regression
baseline must follow immediately, with no unrelated work in between.

---

## Phase 2B — Validation baseline

**Status: ACTIVE — MERGE GATE SATISFIED.** `GOVERNANCE.md` §3.1 **condition 7 is SATISFIED** —
the passive-cooling regression baseline landed in increment 2B-1R, together with the integrator
repair it required. Two items remain outstanding — **CI** and the **full convergence scope** — so
the phase stays ACTIVE and is **not** COMPLETE.

**The exit criteria below are nevertheless SATISFIED, and so is the Phase-3 prerequisite**
("Phase 2B baselines exist"), which is what downstream phases key off. The Phase-2B closure audit
found the branch merge-ready on evidence: `docs/validation/PHASE2B_CLOSURE.md`. Measured
**13/13 CTests** with the authenticated CompOSE data root and **8/8** without it (the five
external-data tests are *excluded* by the CMake guard, not skipped).

**Prerequisite:** Phase 2A complete. ✅ **SATISFIED** — increments 2A-1 through 2A-3.

The codebase has zero assertions and zero CI. Until baselines exist, no numerical change can be
shown correct. Phase 1 supplied the *mechanism* and Phase 2A the first physics checks; Phase 2B
expands that plumbing into the actual scientific baseline.

- ◐ **Expand the Phase-1 test plumbing into a full harness; add CI.** The harness half is
  **done** — `tests/` went from the single Phase-1C smoke test to **13 registered CTests** (at Phase-2B closure; **19** after Phase 3) across
  Phase 2A–2B, each with a durable record under `docs/validation/`. **CI is the remaining
  Phase-2B item and is absent**: no `.github/`, no CI configuration of any kind. It blocks the
  COMPLETE label; it does **not** block merge or Phase-3 entry, neither of which references it.
  Minimum viable scope (self-contained CI over the 8 data-free tests, versus a full external-data
  scientific CI that has no reproducible provisioning path for the ~220 MB authenticated tables)
  is analysed in `docs/validation/PHASE2B_CLOSURE.md` §5.
- ✅ **TOV reference checks against known solutions — COMPLETE (2B-2).** The production TOV
  path is validated against references that are independent of CompactStar. Tier A: the exact
  Schwarzschild constant-density interior solution, evaluated directly against
  `TOVSolver::ODE` at two relativistic compactness values (`2GM/Rc² = 0.30, 0.50`) — max
  relative deviation **`3.5e-16`**, i.e. the production right-hand side is correct to roundoff.
  Tier B: the official CompOSE `eos.mr` for DS(CMF)-1_with_crust — `M_max` to **`2.8e-4`**,
  radii to **0.20–0.35 %**, plus published CMF anchors. The radius residual is systematic and
  **fully attributed** to the outer-boundary convention: production terminates at the EOS
  table floor (`n_B = 1e-7 fm⁻³`, still inside the outer crust), and the residual is a
  constant fraction (0.46–0.49, 8 % spread) of the hydrostatically-estimated omitted layer at
  every mass. Two pre-existing characteristics documented and deliberately **not repaired**
  (source frozen): that boundary convention, and a default `r_max = 70 km`/`radial_res = 10000`
  grid that costs ~2.7 m of radius accuracy while leaving the mass converged to nine figures.
  CTests `tov_reference_analytic` (self-contained) and `tov_reference_cmf` (external data,
  SKIPs at 77 without it); artifact `tests/baselines/tov_dscmf1_reference.tsv`. Evidence:
  `docs/validation/TOV_REFERENCE.md`.
- ✅ **First-order Hartle moment-of-inertia validation — COMPLETE (2B-4B).** The scale-free
  observable `I = J/Ω` is independently validated. The production first-order equation is
  derived term-by-term to be the published Hartle equation (`EQUATION MATCH`), and the
  arbitrary central normalization cancels from `I` **analytically** and numerically to
  **3.3e-15 over six decades** of central amplitude. An independent solver — conservative
  form, different state variables, built from the metric rather than production's coefficient
  helpers — agrees with production to **9.5e-9** on an exact constant-density background and
  to **1.3e-5 – 2.1e-5** across the CMF sequence at 1.0/1.4/1.6/2.0 M☉. Surface and volume
  extractions agree to `<8.1e-7`; both production and reference reproduce `I/(MR²) → 2/5`
  monotonically in the weak field; `I = 1.17–2.61e45 g cm²` sits within 7 % of the
  Breu–Rezzolla and Lattimer–Schutz relations (used as **sanity checks only** — the official
  DS(CMF)-1 `eos.mr` carries no `I`, and no EOS-specific published sequence exists). `I`
  inherits the 2B-4A surface jitter but the `r⁴` weighting does **not** amplify it (spread
  `6.6e-5` vs the radius's `6.9e-4`), and the production/reference agreement is
  resolution-independent. Three controlled detector mutations fire; a fourth (centre boundary)
  provably cannot, for a documented reason. CTests `hartle_moment_inertia_analytic`
  (self-contained) and `hartle_moment_inertia_cmf` (external data), labels
  `rotation;hartle;scientific`. **This does NOT resolve INV-07:** absolute `ω̄`, physical `Ω`
  and `J`, the `init_omega_bar` seed, and the `Ω [s⁻¹]` mislabel all remain unresolved, and
  **no claim is made about O(Ω²)**. Phase 4 remains blocked. Evidence:
  `docs/validation/HARTLE_MOMENT_INERTIA.md`.
- ◐ **Grid-convergence harness — PARTIAL (2B-4A).** The **nonrotating TOV → cooling slice is
  measured and CHARACTERIZED**. The first-order Hartle validation it was waiting on **has since
  landed** (2B-4B, `HARTLE-I VERIFIED` for the scale-free `I`), so that half of the dependency is
  discharged — but the full ratified TOV → Hartle → cooling convergence **still cannot be
  performed**, because no Hartle output reaches the thermal right-hand side at all today:
  `DriverContext` exposes only `star`, `geo`, `envelope`, `thermo`, `cfg`, and the thermal drivers
  reference no Hartle quantity. Creating that coupling needs the physical Ω normalization
  (INV-07, Phase 4) and either O(Ω²) response (Phase 4) or rotochemical heating (Phase 5).
  **The scope of this item is therefore left UNCHANGED pending owner adjudication** — relocating
  it would alter ratified scientific meaning. Evidence and exact proposed wording:
  `docs/validation/PHASE2B_CLOSURE.md` §6.
  Per INV-13, interpolation is **linear** and `DataSet::Integrate` is the trapezoid rule, whose
  nominal accuracy is O(Δr²) for sufficiently smooth integrands — so fourth-order behavior is
  not available from these components. The convergence order of a *complete* coupled observable
  is **measured by this harness, not assumed from the interpolation scheme** — and the
  measurement contradicts the nominal expectation: over `radial_res = 5000…40000` (measured
  Δr ratios 2.00), the observed order is **≈ 0.70** for `R`, `e^ν(R)` and `L_ν`, and **0.742**
  for the cooling-trajectory norm — **not 2** — because the stellar radius is fixed by a
  step-dependent *surface-termination event* at the EOS table floor rather than by a quadrature.
  `C_star` converges separately (core-dominated, `~1e-6` by the default grid). The
  variable-resolution builder is **bit-exact** against canonical `NStar::SolveTOV_Profile` at
  `radial_res = 10000` (all twelve checks, including the full 9-epoch trajectory, rel `0`).
  Thermal-integrator error is **subdominant by 2.7e3**; the target-mass search contributes
  **zero** differential error above 5000 (bit-identical `ε_c`). Default-grid continuum error:
  `R ~5.2e-4`, `T_inf ~4.5e-4`. **Accuracy adequacy is UNRESOLVED** — no governed downstream
  discretization budget exists, and the regression tolerances are deliberately not
  reinterpreted as one. Coarse-grid detector: `radial_res = 2500` is caught at 130x the
  default-grid radius error. **⚠ Reinterpreted 2026-09-03 (TOV-RR-01,
  `docs/validation/TOV_RADIAL_RES_2500_AUDIT.md`):** that 130x is **not** coarse-grid
  discretization error. At `radial_res = 2500` the TOV integration terminates at the crust–core
  transition — `R = 12.9042` km, `M = 1.59768` M☉ against `13.4635` / `1.59998` at 5000 — so both
  `radial_res = 2500` rows of the artifacts are a star missing its outer crust, and the
  `B_fixed_mass` row's `ε_c` is 0.22 % high because the mass root-finder compensated. The other
  four resolutions are clean and monotone, so the convergence conclusion stands on them; the 2500
  point is not a coarse-grid value of the same star. Not repaired — every candidate repair changes
  `R_*` and needs an ADR. **That ADR is accepted: `docs/adr/ADR-0009-tov-surface-event-and-termination.md` (ACCEPTED, TOV-SURF-I 2026-09-03; source conformed and validated in TOV-SURF-I-R; evidence `docs/validation/TOV_SURFACE_CONTRACT_DERIVATION.md`).** After the completed implementation/validation, this artifact is regenerated on the corrected solver with `2500` retained as a valid coarse point and the refinement set extended (evidence §13). CTest `grid_convergence_cmf` (labels
  `thermal;scientific;external-data;convergence`, ~106 s); artifacts
  `tests/baselines/grid_convergence_cmf_1p6_{debug,trajectory}.tsv`. Golden baselines
  unchanged. Evidence: `docs/validation/GRID_CONVERGENCE.md`.
- ✅ **Cache-correctness checks (INV-12) — COMPLETE (2B-3).** The supported same-star cache
  contracts are durably verified under CTest (`cache_contract`, `cache_thermal_contract`,
  label `cache`): mass-density, `Y_q` and Direct-Urca rebuild exactly on a sanctioned profile
  mutation and are stable on repeat; `ProfileVersionedCache` honors its same-star rebuild
  contract; the `NeutrinoCooling` payload rebuilds with `L_nu` scaling exactly with `rho` and
  drops the DU channel when composition closes it. **Five known INV-12 hazards are reproduced
  and quantified** — stale `GeometryCache` with no provenance (51.6 % geometry divergence,
  undetectable by any caller), the `C_star` key omitting the geometry (50 % error), the
  version-only generic key colliding across equal-version profiles (85.7 %), the concrete
  cross-star `NeutrinoCooling` collision (80 %), and `StarContext` column pointers never
  re-bound after a structural change. They live behind `--audit-known-hazards` and are
  deliberately **not** CTests: a known defect is never laundered into a green assertion.
  **The canonical passive-cooling baseline provably reaches none of them** — asserted over 602
  driver-context observations (one profile version, one `GeometryCache`, one `StarContext`,
  one thermo). Three controlled regressions were shown to be caught, then reverted exactly.
  **INV-12 architectural repair — DEFERRED TO PHASE 3**; the cache system is *not* corrected.
  Evidence: `docs/validation/CACHE_CORRECTNESS.md`.
- ✅ **Passive cooling regression — COMPLETE (2B-1, blocked → 2B-1R, established).** Captures the
  coherent `C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞` on the authenticated 1.6 M☉ CMF star over
  100 yr → 1 Myr: nine log-spaced epochs, one continuous integration, energy identity closing to
  2.1e-16, state tolerance 1e-5 derived from measurement, and a 1% photon perturbation caught
  1000× over. Golden values in `tests/baselines/passive_cooling_cmf_1p6_debug.tsv`; CTest
  `passive_cooling_regression` (labels `thermal;scientific;external-data;regression`), registered
  only when `COMPACTSTAR_EOS_DATA_ROOT` is supplied. **`GOVERNANCE.md` §3.1 condition 7 is
  SATISFIED.** Required first repairing two pre-existing defects — an unusable `MSBDF` default
  with no Jacobian, and an `EvolutionSystem` block that made a Spin state mandatory. Evidence:
  `docs/validation/PASSIVE_COOLING_BASELINE.md`.
  It remains a regression baseline, **not** a physics validation: the neutrino emissivity
  normalizations are still self-labeled placeholders (`NeutrinoCooling_Details.cpp:100-102`).

**Exit criteria. ✅ SATISFIED.** Baselines exist and run; a regression is detectable. Six
independent detector proofs were demonstrated and reverted — passive cooling (a 1 % photon
perturbation caught 1000× over), cache invalidation (×3), the coarse radial grid, and the Hartle
J/Ω/ODE extractions (×3). **Satisfying the exit criteria is not the same as completing the phase:**
CI and the convergence scope above remain outstanding items.

---

## Phase 3 — Behavior-preserving consolidation

**Prerequisite:** Phase 2B baselines exist. ✅ **SATISFIED** — `docs/validation/PHASE2B_CLOSURE.md`.

**Status: IN PROGRESS — increments 3A, 3B, 3C, 3D and 3E complete (3E-0, I1, I2, I4; I3
optional and not taken). `SingleStarSolveToTOVPoints` is the canonical TOV numerical primitive
and, after I4, the ONLY ordinary-star radial implementation — `RadiusLoop` is deleted and all
three orchestrators delegate to the one primitive. Fail-closed condition #3 is CLOSED for
ordinary visible-sector TOV radial ownership (not MixedStar, not dark-sector).**
**3A** centralized the exactly-duplicated unit constants (bit-identical). **3B** implemented the
ADR-0003 cache-provenance contract (goldens bit-identical; INV-12 resolved for profile-derived
caches). **3C** adopted the authoritative ZakiLib Boltzmann constant. **3D** established the
canonical proper-volume owner per ADR-0004 — on the validated path; **3E-I2 then conformed TOV
Path-1 as well**, so both ordinary visible-sector `NStar` paths share one geometry owner.
`MixedStar`, the §5 core candidates, the project-specific extension modules and the INV-14
scalar accessor remain deferred by governed decision. **Phase 3 is COMPLETE on its ratified exit
criterion** ("one authoritative owner per quantity; baselines still pass") — see
`docs/validation/PHASE3_CLOSEOUT.md`, the durable authority for the merge decision. The branch
was merged to `master` at **`df859b5`** and post-merge verified (19/19, 10/10, seven artifact
hashes unchanged; re-authenticated again at Phase-4 entry).
The durable plan —
ownership maps for all five targets, the cross-target dependency graph, the baseline-protection
matrix, the per-increment preservation standard and the recommended sequence — is
[`docs/architecture/PHASE3_CONSOLIDATION_PLAN.md`](architecture/PHASE3_CONSOLIDATION_PLAN.md).

> *(Scope note, Phase 3F: the blockquote below records the state **at Phase-3 entry**. Both
> fail-closed conditions it names have since been discharged — #3 by ADR-0005 / 3E-I4 for the
> ordinary visible sector, #4 by ADR-0003 / 3B for profile-derived caches.)*
>
> **Class correction from that audit.** The line below states that every item is engineering
> class. That is **not** what the higher-ranked documents require. Under `GOVERNANCE.md:51` a
> change that "moves ownership, changes boundaries, promotes or retires a code path" is
> **structural** and needs an ADR plus a same-change `CURRENT_ARCHITECTURE.md` update; two
> `GOVERNANCE.md` fail-closed conditions are already active (`:70` uncertain authoritative code
> path — the two live TOV paths; `:72` uncertain cache validity — the five INV-12 hazards); and
> `SCIENTIFIC_INVARIANTS.md:147` (INV-04) states outright that single proper-volume ownership is
> a "structural change; requires ADR". `GOVERNANCE.md` (rank 1) and `SCIENTIFIC_INVARIANTS.md`
> (rank 3) outrank this roadmap (rank 4), so **the canonical-TOV, proper-volume and cache items
> each require an ADR** and cannot be executed as pure refactoring. The recommended sequence also
> deliberately departs from the bullet order below — TOV is the highest-risk item and the only
> one with a live but wholly unvalidated path, so it is scheduled last rather than first. The
> bullet list is left unchanged; nothing here alters Phase-3 scientific scope.

Every item is **engineering class** and must produce bit-identical output, or a documented
tolerance.

- ✅ **Establish the canonical TOV path; retire or clearly subordinate the duplicate — COMPLETE
  (Phase 3E).** 3E-0 measured the two live paths **bit-identical** in all 25 radial columns
  (H2: one algorithm copied into two loops). **ADR-0005 ACCEPTED**: Q1 retain `Solve()` as a
  subordinate workflow orchestrator, Q2 preserve the `_Sequence.tsv` contract, Q3 = P3 staged,
  Q4 preserve the hooks. **I1** made `Solve()` delegate to `SingleStarSolveToTOVPoints`.
  **I2** conformed Path-1 `Append`/`FinalizeSurface` to `CompactStar::Geometry`, after which
  Path-1 and Path-2 `B` are **bitwise identical**. **I4** covered `GenTestSequence` *before*
  migrating it (16/16 bitwise, output byte-identical), then **deleted `RadiusLoop`**.
  `SingleStarSolveToTOVPoints` is now the **only** ordinary-star radial implementation, and
  `Solve`, `SolveToProfile` and `GenTestSequence` all delegate to it. Both file contracts are
  preserved and guarded; none of the callers changed.
  **`GOVERNANCE.md` fail-closed condition #3 is CLOSED for ordinary visible-sector TOV radial
  numerical ownership** — explicitly *not* `MixedStar`/`RadiusLoopMixed` two-fluid integration or
  the dark-sector paths, which are distinct physics scopes rather than competing implementations.
  **`3E-I3` remains OPTIONAL**: the two `NStar` construction styles are value-equivalent and
  geometry-conformant but not textually consolidated. Evidence:
  `docs/validation/TOV_PATH_EQUIVALENCE.md`, `docs/validation/PHASE3E_I1_CANONICAL_TOV.md`,
  `docs/validation/PHASE3E_I2_PATH1_GEOMETRY.md`,
  `docs/validation/PHASE3E_I4_RADIUSLOOP_RETIREMENT.md`; contract:
  `docs/adr/ADR-0005-canonical-tov-numerical-primitive.md`.
- ◐ **Single owner for unit conversions — PARTIAL: exact duplicate constants completed in 3A.**
  `KM3_TO_CM3` (2 production sites) and the MeV fm⁻³→erg cm⁻³ factor (2 sites) now have one
  owner, `CompactStar/Units.hpp` — a dependency-free header that adds no edge to the layer
  graph. **Bit-identical**: 352 lines of deterministic scientific output unchanged, all five
  golden artifact hashes unchanged, 13/13 and 8/8 green. Evidence:
  `docs/validation/PHASE3A_UNIT_DUPLICATES.md`. **k_B — ✅ RESOLVED in 3C.** The two local
  precisions (`8.617333262145e-11`, `8.617333262e-11`) are gone; all four production consumers
  now read the authoritative `Zaki::Physics::K_BOLTZ_EV * 1.0e-6`, and no GSL Boltzmann constant
  remains on a production path. ZakiLib was **not** modified and the archive was **not** rebuilt —
  it already owned the constant and CompactStar already linked the symbol. Fixed-state response
  is the analytic `1.68e-11`; the TOV and Hartle-I goldens stayed **byte-identical** as predicted,
  and only the three thermal artifacts were re-baselined. Evidence:
  `docs/validation/PHASE3C_BOLTZMANN_AUTHORITY.md`. **Still outstanding:** the **solar-mass authority** conflict (`SUN_M_KM` vs
  `GSL_CONST_CGSM_SOLAR_MASS`, differing at `6.2e-5`) remains a scientific/unit authority
  question **deferred out of Phase 3** pending owner or ADR adjudication.
- ◐ **Single owner for the proper-volume measure (INV-04) — ADR-0004 ACCEPTED; canonical
  validated path CONFORMED; legacy migrations DEFERRED.**
  `docs/adr/ADR-0004-proper-volume-and-metric-measure-ownership.md` was **ACCEPTED 2026-09-01**.
  Q1 = **Option B**: `CompactStar/Geometry.hpp` is the dependency-neutral *mathematical* owner of
  `f = 1 − 2m/r`, `Λ`, `e^{Λ}` and `w_V = 4πr² e^{Λ}` **including domain semantics**, while
  **`GeometryCache` is retained as the canonical *cached-representation* owner** and consumers
  keep their own physics factors and unit conversions. `Core` was **not** made to depend on
  `Physics/Evolution`. Q2 = `MixedStar` **governed now, source migration deferred** until focused
  coverage exists. Q3 = the **hybrid physical-domain contract**: exact regular-centre limit at
  `r = m = 0`, fail closed for `r < 0`, `r = 0` with `m ≠ 0`, `f ≤ 0` and non-finite input, and
  **no `1e-15` clamp** — replacing six mutually inconsistent legacy behaviors.
  This supersedes the Phase-3-entry wording *"`GeometryCache` canonical; retire the inline
  forms"*, which the 3D audit found not implementable as written.
  **Conformed:** `NStar::BuildFromTOV` (Λ bit-identical; `|ΔB|/B = 1.368e-16` against the
  `1.0e-15` predeclared before implementation) and `GeometryCache::DeriveLambdaFromMR_` (cached
  arrays bit-identical). **Not conformed, deferred:** TOV **Path 1**, the scalar
  `NStar::BaryonNumIntegrand(double)` (separate INV-14 defect, deliberately unrepaired),
  `MixedStar`'s six sites, and the §5 candidate code.
  **INV-04 is therefore `GOVERNED (ACCEPTED) — CANONICAL VALIDATED PATH CONFORMED; LEGACY
  MIGRATIONS DEFERRED`, not resolved.** Evidence: `docs/validation/PHASE3D_PROPER_VOLUME.md`.
- ✅ **One uniform cache-invalidation rule; version gate on `GeometryCache`; `StarContext`
  re-binds column pointers on invalidation (INV-12) — ADR-0003 ACCEPTED and IMPLEMENTED.**
  `docs/adr/ADR-0003-profile-cache-provenance-and-invalidation.md` was **ACCEPTED 2026-09-01**
  (Q1 = S1, Q2 = Option A) and implemented in increment 3B: a `(profile identity, version)`
  provenance token, `GeometryCache` carrying and exposing it, dependency-complete keys for `C_⋆`
  (+ geometry) and the `NeutrinoCooling` payload (+ identity, geometry), and `StarContext`
  re-binding its column views on a version change. **All five former hazards are now enforced
  CTests and no longer reproduce**; the goldens were byte-identical. **INV-12 is RESOLVED for
  profile-derived caches.** Evidence: `docs/validation/` and the ADR's implementation record.
- Classify dead and unreachable code; retire only after dependency review.

**Heat capacity is no longer a Phase-3 item.** Its physical ownership is governed by **ADR-0002**,
and its minimum source conformance is a **Phase-2A** prerequisite.

**Open for separate consideration in this phase:** whether the thermal-energy balance should be
**centralized** — drivers exposing power contributions (`−L_ν`, `−L_γ`, `+L_H`, …) with one
thermal-balance owner performing `dT∞/dt = L_net / C_⋆(T∞)` — instead of each driver dividing by
the shared `C_⋆` itself. ADR-0002 §6 states both patterns and deliberately decides neither.
Centralizing changes architectural ownership and the driver RHS contract, so it requires its own
ADR and must be evaluated **after** baselines exist. It is not a prerequisite for conformance.

**Exit criteria.** One authoritative owner per quantity; baselines still pass.

---

## Phase 4 — Rotation correctness

**Prerequisite:** Phase 3 ✅ (merged at `df859b5`); ADR on Hartle normalization accepted
✅ — **ADR-0006 ACCEPTED 2026-09-02** (Q1 = A + D, Q2 = A, Q3 = A, Q4 = A; binding clarification:
no implicit physical spin on `NStar` construction). **The Phase-4 prerequisite is satisfied.**

**Status: ENTRY AUDIT COMPLETE (increment 4A-0, 2026-09-02) —
`PHASE-4 ROTATION ENTRY AUDIT COMPLETE — ADR-0006 OWNER ADJUDICATION REQUIRED`.**
Evidence: `docs/validation/PHASE4_ROTATION_ENTRY.md`. Documentation only; no production change.
What the audit established: (a) `RotationSolver` provenance is fully resolved to four lineages
(owner first-order, AI candidate `675b4a9`, merge repair `57334d8`, later validation-only) with
no UNKNOWN block, and the source is byte-identical to the Phase 2B-4B-validated source;
(b) the INV-07 normalization contract is **derived** — `A = (Ω_phys/c)/Ω_raw`, `ω̄_phys = Aω̄_raw`,
`J_phys = AJ_raw`, `I` unchanged — and checked through the public API to `2e-16`; the hard-coded
seed corresponds to `Ω ≈ 2.2–2.4×10³ s⁻¹`; (c) the O(Ω²) candidate is **publicly callable with
zero repository callers** (the "structurally unreachable" wording was wrong and is corrected in
INV-08 and `CURRENT_ARCHITECTURE.md`), executes, and every output scales as the seed squared;
its equations are mapped term by term against Hartle's l = 0 system and found defective beyond the
`j²` factor (dimensionally inconsistent homogeneous `p0` equation, wrong source terms, a
non-Hartle boundary condition `δp(R) = 0` that shifts the central density, incomplete `δM`);
its `p0` is the Eulerian `δp`, not Hartle's `p₀*`.

**Proposed internal sequence (dependency order; the ratified item list below is unchanged):**

| Increment | Content | Gate |
|---|---|---|
| **4A** | ✅ **COMPLETE 2026-09-02.** ADR-0006 accepted, then implemented: typed physical spin input (`Ω [rad s⁻¹]`) with geometric internals and one conversion owner; seed internalized behind a non-public seam; `HartleResult` stripped of its five seed-normalized first-order fields (the `[s⁻¹]` mislabel is gone); seed-free `HartleFirstOrderResponse` and explicit `NStar::RotationAt(AngularVelocity)`; export headers corrected. **No implicit physical spin on construction.** V1–V9 pass at the predeclared bounds, four detectors fired and were reverted byte-identically, `I` bit-identical, seven goldens byte-identical, O(Ω²) byte-identical. Evidence: `docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md` | ADR-0006 ✅; `I` and the seven goldens bitwise ✅ |
| **4B** | ✅ **COMPLETE 2026-09-02.** Independent *physical* validation of the normalized first-order response: node-by-node comparison of `ω̄/Ω` and `ω̄'/Ω` against an independently derived and independently normalized profile (`2.9e-9` analytic vs a predeclared `1e-7`; `≤ 2.3e-5` on the four CMF stars vs `1e-4`), exterior-matching identities against an independent `I`, the volume-integral identity from production's own shape, and two **derived** weak-field coefficients. Detector D1 proves the coverage is new. **Zero production change.** Evidence: `docs/validation/PHASE4B_FIRST_ORDER_PHYSICS.md` | 4A ✅ |
| **4C** | **✅ COMPLETE — 4C-G + ADR-0007 ACCEPTED + 4C-I0 (2026-09-02) + 4C-I1 (2026-09-03). 4D executed 2026-09-03 — see the 4D row.** The O(Ω²) monopole contract was derived from the **primary source** (Hartle 1967 eqs. 87–88, 97, 99–100, 105–108, read from the journal scan; Hartle & Thorne 1968 §II secondary), not from the candidate: variable `p₀*`; the `l = 0` system per `Ω_geom²` from the verified `s`, `s'` (dimensionally audited); fixed-`ε_c` centre condition with a regular-series start and **no** surface condition; EOS-owned `dε/dp` (an API prerequisite — none exists); `δM = m₀(R_*) + 4πR_*²ε_*ξ₀(R_*) + J²/R_*³`; surface classified **`SURFACE ADEQUATE AS-IS`** with explicit `R_*` semantics; `l = 0` sufficient for all scalar counts (derived from the angular integration); the complete Phase-5 `δN_i` interface (five contributions); §3.1 conditions satisfiable; atomic replacement of the public candidate recommended. Every 4A-0 defect hypothesis classified against the primary text. Evidence: `docs/validation/PHASE4C_HARTLE2_DERIVATION.md`; `docs/adr/ADR-0007-hartle-second-order-monopole-response.md`. **Owner adjudication 2026-09-02:** `p₀*`; fixed `ε_c`; coefficients per `Ω_geom²` materialized only at an explicit `AngularVelocity`; EOS/TOV owns `dε/dp`; surface adequate as-is; **validated `l = 0` is the Phase-4 second-order deliverable that unlocks Phase 5, with `l = 2` shape/quadrupole a separate future rotation extension that blocks neither Phase-4 completion, Phase 5, nor the BNV thermal program (and is not thereby validated)**; atomic candidate replacement in 4C-I1; §3.1 authorized, not yet exercised; the homogeneous sequence-derivative response is not public in 4C. Implementation splits into **4C-I0** and **4C-I1**. **4C-I0 ✅ COMPLETE 2026-09-02:** `TOVSolver::GetEDensDeriv` is the sole `dε/dp` authority — the analytic derivative of the star's own `gsl_interp_steffen` `ε(p)` interpolant, dimensionless (conversion asserted against an independent `c` to `1.5e-16`), fail-closed off-domain, carried to consumers by `StarProfile::HasEosDEdP()`/`GetEosDEdP()` all-or-nothing, with an explicit supply path for point-constructed analytic stars. 15/15 self-contained + 17/17 CMF checks, detectors D1–D3 fired and reverted byte-identically; the retired profile FD measured at 155–490 % error on ~25 crust nodes per star. **`RotationSolver` untouched, candidate byte-identical, seven artifacts and `I` byte-identical, no monopole baseline.** `dn_i/dp` deferred to Phase 5. Evidence: `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md`. **4C-I1 ✅ COMPLETE 2026-09-03 — the `GOVERNANCE.md` §3.1 correction is EXECUTED:** the AI candidate (`SolveHartle2_N`, `ODE_Hartle2_N_Fast`, the MixedStar stubs, `HartleResult`, `GetHartleResult`, all candidate-only members) was **deleted** and the governed fixed-`ε_c`, `Ω²`-normalized `HartleMonopoleResponse` replaced it in the same commit — regular-centre series start, no shooting, background interpolated at the actual ODE radius, EOS-owned `dε/dp` with fail-closed absence, `δM̂` with shell and `I²/R_*³`, provenance-guarded caching, explicit `NStar::ComputeHartleMonopoleResponse()` (never automatic), materialization by scaling. Seed invariance `7.85e-15`; `I` bit-identical; seven artifacts unchanged; 27/27 + 14/14; detectors D1–D4 fired. **NOT validated physics — no monopole baseline.** Evidence: `docs/validation/PHASE4C_I1_MONOPOLE_IMPLEMENTATION.md` | ADR-0007 ✅ ACCEPTED; 4A/4B first-order layer ✅ |
| **4D** | **✅ EXECUTED 2026-09-03 — status `HARTLE MONOPOLE VALIDATION FAILED` (implementation verified; contract gap on stepped crusts).** The governed monopole response was recomputed by two test-only routes that never call production's ODE — Hartle's `(m₀, h₀)` system (97)+(98)+(90) on the tabulated background and a continuum solver on the closed-form constant-density interior — and compared node by node: analytic agreement `9.7e-9` (bound `1e-7`) converging as `h²` to the continuum; centre series `≤ 5e-10`; shell identity exact (shell = 90 % of `δM̂` on the homogeneous star); Newtonian intercepts `1e-6`; continuum first integral `6e-15`; homogeneous `δM̂` vs exact `dM/dp_c` `3e-9`; DS(CMF)-1 four stars `≤ 3.8e-7` (isolated) / `≤ 3.7e-5` (fully independent; bound `1e-4`); **published:** Chandrasekhar & Miller 1974 Table I 19/19 to `≤ 7.3e-4` (incl. the `ξ₀` sign change and the shell-excluded `δM/M` — C&M's own boundary condition drops the surface shell), Hartle & Thorne 1968 Tables 3/5 on the printed HW EOS 8/8 (`δR/R ≤ 4.9e-3`, `δM/M ≤ 1.1e-2`; bound `2e-2`); detectors M1–M9 all fire, reverted byte-identically; production diff **NONE**. **Not met:** ADR-0007 §7 item 11 on DS(CMF)-1 — `1.04e-3` vs `1e-3`, resolution-independent (10000/20000/40000), diagnosed as **crust density steps the nodal `dε/dp` column cannot represent** (17 % of the crust's `Δε` missing at every resolution). **The Stieltjes evaluation against the profile's own density steps reproduces the TOV-sequence derivative to `7e-5` and, applied to the sourced solution, puts the omitted internal delta-function shells at ≈ `4.6 %` of `δM̂` on DS(CMF)-1 — a substantive physical discrepancy of the accepted contract, reported and NOT repaired in 4D**; item 10 has no independent source (`c_s²` unavailable; FD moves `δM̂` by 5 %); item 9 order `2.35`, Richardson residual `7.7e-4`. Two labelled re-scopings (Newtonian intercept per the ADR's wording; tabulated first integral taken at the continuum level). **No monopole baseline created.** Evidence: `docs/validation/PHASE4D_MONOPOLE_VALIDATION.md` | 4C ✅ |
| **4E** | **◄ BLOCKED on the corrected revalidation.** (a) **ADR-0008 ACCEPTED and IMPLEMENTED (Phase 4D-RI, 2026-09-03)** — `docs/validation/PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md`: per-segment measure source with profile nodes as mandatory integration boundaries, terminal `ε_*→0` atom applied once, `dε/dp` retained for the regular-centre series; constant-density output bitwise unchanged, smooth-EOS equivalence `≤1.3e-5`, same-partition accounting `≤3.8e-7`, DS(CMF)-1 `δM̂` `+3.2…+6.5 %`, radial spread `1.6 %→3.7e-5`, detectors D1–D4 fire, seven artifacts and `I` unchanged, **no baseline**. **◄ NEXT — sequenced behind ADR-0009:** the TOV-RR-01 audit (`docs/validation/TOV_RADIAL_RES_2500_AUDIT.md`) found the ordinary-star surface located by a fatal trial-state guard (`R_*` one output step short, scattered truncations incl. ≈3 % of central densities at the default 10000), and `docs/adr/ADR-0009-tov-surface-event-and-termination.md` (ACCEPTED; source conformed and validated in TOV-SURF-I-R) governs the correction; **the corrected independent revalidation** (ADR-0008 §Validation A–J with the `(m₀,h₀)` and continuum oracles and M1–M10 re-run in full — line D's non-monotone residual is this surface artifact) **and the first monopole baseline wait for the artifact migration following ADR-0009's completed implementation and validation**, because `R_*` (+5 m at 1.6 M☉), `I²/R_*³` and the terminal atom's `ε_*` move. Superseded plan text: (a) **ADR-0008 PROPOSED (2026-09-03, Phase 4D-RG)** — `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md`: the EOS source of (97) as the measure `−4πr²ξ̂₀ dε` (H67 eq. 93), per-segment secant source on the profile partition (optional EOS-knot refinement), exact jump operator at declared discontinuities, surface shell as the terminal atom, `dε/dp` retained for the centre series/diagnostics, Phase-5 measure completeness, validation lines A–J, baseline only after re-validation. Scratch evidence in `docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md` (DS(CMF)-1 `δM̂` converged to `4e-5`, +4.8 %; Option A non-convergent; two-layer jump operator `2e-10…4e-9`). **Owner adjudicates Q1–Q12; then the correction increment (§3.1), re-validation, then the first monopole baseline** (§3.1 condition 7). (b) The Phase-5 structural fields themselves — `m₀/Ω²`, `δp₀/Ω²`, `ξ₀/Ω²`, `δM/Ω²`, `ω̄/Ω`, `I` — are **already carried by `HartleMonopoleResponse` / `PhysicalHartleMonopole` (4C-I1) and verified at the implementation level (4D)**; 4E reduces to the closeout plus a one-line ratification that those fields are the Phase-5 interface. **Phase 5 must not begin before (a).** | 4D ✅ (characterized) |

`MixedStar` two-fluid rotation is a separate track (compiled, unexercised, unvalidated; ADR-0004
§0-Q2 coverage first) and blocks no ordinary-`NStar` item. **Phase 5 must not begin before 4E.**

- **Re-audit `RotationSolver` and `MixedStar` against `9f70f14`.** The Phase-0 findings were made
  against the superseded `d91c31b` and must not be carried forward.
- ✅ **Physical Hartle normalization — COMPLETE (4A).** The seed is now private with no public
  setter (its value is deliberately unchanged at `5e-3`, since it carries no physical meaning
  and changing it would move `I` for no gain); a physical Ω is supplied through the typed
  `AngularVelocity` and the solution is scaled to it; the `[s^-1]` annotation is gone with the
  field that carried it (INV-07 CONFORMED).
- Verify true RHS-radius interpolation and cached bracket indices as introduced by `3639d71`.
- **Verify the second-order equations against a cited derivation** — restore the j² factor,
  complete δM with the exterior term, source `dε/dp` from the EOS rather than profile
  differences, impose proper central series expansions (INV-08). *(4C-G, 2026-09-02: derived
  from the primary source and governed in ADR-0007 PROPOSED — the derivation is done; the
  correction awaits adjudication and 4C-I.)*
- Make `HartleResult` reachable from `NStar`. *(4A-0 re-scoping: a raw `HartleResult` is already
  reachable through an external `RotationSolver`; what is needed is a **normalized** response
  accessor on `NStar` — ADR-0006 Q4, increment 4E.)*
- Convergence tests and validation against published moment-of-inertia and δM values.

**Exit criteria.** O(Ω) and O(Ω²) validated and reachable. The candidate status of `675b4a9` is
resolved — ratified or replaced. *(**Owner adjudication 2026-09-02, binding:** validated `l = 0` O(Ω²) structural response is the
Phase-4 second-order deliverable required to unlock Phase 5. `l = 2` shape/quadrupole physics is a
**separate future rotation extension** and blocks neither Phase-4 completion, Phase 5, nor the BNV
thermal program. No claim is made that `l = 2` physics is itself validated — it is out of scope.)*

---

## Phase 5 — Standard non-superfluid rotochemical heating

**Prerequisites:** Phase 4 · **ADR-0001 accepted ✅** · **ADR-0002 accepted ✅** · ADR on η conventions accepted ☐.

> **Species-semantics prerequisite: SATISFIED** (ADR-0001, 2026-08-31).
> **Phase 5 remains blocked** by every other prerequisite below.

- **Correct `RotochemicalCache` for ADR-0001 conformance** — construct `n_i = Y_i · n_B` before
  the `N_i`, `A_i`, and `B_i` species number-density integrations
  (`RotochemicalCache.cpp:147`, `:25-44`, `:47-104`). *Recorded in Phase 0.5; deliberately not
  implemented there.*
- Correct `A_i` (divide by Ω²) and `B_i` (geometry-consistent finite difference).
- Confirm the `Z_i` reduction under the ratified species semantics.
- Define chemical state: η_npe and η_npμ, ordering, redshift frame, units (INV-11).
- **Implement out-of-equilibrium weak rates** — ΔΓ(η,T) and the F/H(ξ = η/k_BT) functions. *This
  is the terminal blocker; nothing downstream exists without it.*
- `WeakRestoration` (currently a 0-byte file).
- Neutrino-rate modification for chemical disequilibrium.
- `HeatingFromChem`, with **single-source Γ** so heating and neutrino losses cannot double count.
  Its `+L_H,∞` term is subject to ADR-0002: it divides by the same governed `C_⋆(T∞)` as every
  other channel.
- Add both files to the build — they are still absent from every CMake source list.
- Fernández–Reisenegger regression.

**Exit criteria.** Standard rotochemical heating reproduces published results.

---

## Phase 6 — BNV integration

**Prerequisite:** Phase 5 validated.

BNV generalizes the rotochemical relation. It requires a governed, validated standard baseline to
generalize *from*. Beginning BNV work earlier would extend an unvalidated formalism and make both
unauditable.

---

## Dependency summary

```
0.5 governance ─► 1 build ─► 2A pre-baseline ─► 2B baseline ─► 3 consolidation ─► 4 rotation ─► 5 rotochemical ─► 6 BNV
  ✅ RATIFIED    ✅ COMPLETE      ✅ COMPLETE      ◐ ACTIVE (gate met)   ✅ COMPLETE (merged df859b5)   ◐ 4A–4D DONE — ADR-0008 ACCEPTED + IMPLEMENTED — CORRECTED REVALIDATION NEXT
                                    │                                 │                │              │
   ADR-0001 species semantics  ✅ ACCEPTED ──────────────────────────────────────────────────────────►│  (gate cleared)
   ADR-0002 heat capacity      ✅ ACCEPTED ─►│  (conformance is 2A work, not a gate)
   ADR thermal-balance arch.   ☐ open, deferred ────────────────────►│  (optional; not a gate)
   ADR-0006 Hartle normalization  ✅ ACCEPTED 2026-09-02 ──────────────────────────────►│
   ADR Hartle O(Ω²) equations  ☐ anticipated (after ADR-0006) ─────────────────────────►│
   ADR η conventions           ☐ open ────────────────────────────────────────────────────────────────►│
```

**The former Phase-2 / Phase-3 circularity is gone.** It ran:

```
Phase 2 baseline  ──blocked by──►  INV-15  ──fixed in──►  Phase 3  ──requires──►  Phase 2 baseline
```

ADR-0002 breaks it by deciding the physical ownership **now**, ahead of any baseline, and by
splitting the correction out of Phase 3 into **Phase 2A**, which precedes the baseline. Nothing in
Phase 2A depends on a passive-cooling baseline: it is validated by independent physical checks.

**One** unresolved invariant still gates the chain: **INV-11** (η conventions). **INV-07** is
**no longer a gate** — resolved as a contract by **ADR-0006 (ACCEPTED 2026-09-02)**; what remains
from it is Phase-4A implementation work, tracked as work rather than as an open decision.

**INV-01** (species semantics) is **no longer a gate** — resolved by ADR-0001. What remains from
it is a single Phase-5 implementation task: `RotochemicalCache` must construct `n_i = Y_i n_B`
before its species-density integrals. That is tracked as work, not as an open decision.

**INV-15** (heat-capacity ownership) is **no longer a gate** — resolved by ADR-0002. What remains
from it is **Phase-2A work** (`PhotonCooling` conformance plus `C_⋆(T∞)` verification) and an
**optional, deferred** Phase-3 architectural question (ADR-0002 §6). Neither is an open decision
about the physics.
