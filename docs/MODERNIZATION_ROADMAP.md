# CompactStar Modernization Roadmap

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

**Status: ACTIVE.** `GOVERNANCE.md` §3.1 **condition 7 is SATISFIED** — the passive-cooling
regression baseline landed in increment 2B-1R, together with the integrator repair it required.
Other Phase-2B items remain outstanding, so the phase stays active.

**Prerequisite:** Phase 2A complete. ✅ **SATISFIED** — increments 2A-1 through 2A-3.

The codebase has zero assertions and zero CI. Until baselines exist, no numerical change can be
shown correct. Phase 1 supplied the *mechanism* and Phase 2A the first physics checks; Phase 2B
expands that plumbing into the actual scientific baseline.

- Expand the Phase-1 test plumbing into a full harness; add CI.
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
- First-order Hartle moment-of-inertia checks against published values.
- ◐ **Grid-convergence harness — PARTIAL (2B-4A).** The **nonrotating TOV → cooling slice is
  measured and CHARACTERIZED**; the full ratified TOV → Hartle → cooling convergence remains
  pending the first-order Hartle validation and the unresolved INV-07 normalization.
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
  default-grid radius error. CTest `grid_convergence_cmf` (labels
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

**Exit criteria.** Baselines exist and run; a regression is detectable.

---

## Phase 3 — Behavior-preserving consolidation

**Prerequisite:** Phase 2B baselines exist.

Every item is **engineering class** and must produce bit-identical output, or a documented
tolerance.

- Establish the canonical TOV path; retire or clearly subordinate the duplicate.
- Single owner for unit conversions — remove the ~15 local re-derivations, including **k_B at
  two precisions**.
- Single owner for the proper-volume measure (INV-04).
- One uniform cache-invalidation rule; add a version gate to `GeometryCache`; re-bind
  `StarContext` column pointers on invalidation (INV-12).
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

**Prerequisite:** Phase 3; ADR on Hartle normalization accepted.

- **Re-audit `RotationSolver` and `MixedStar` against `9f70f14`.** The Phase-0 findings were made
  against the superseded `d91c31b` and must not be carried forward.
- Physical Hartle normalization — replace `init_omega_bar = 5e-3` and scale to a physical Ω
  (INV-07); correct the `[s^-1]` unit annotation.
- Verify true RHS-radius interpolation and cached bracket indices as introduced by `3639d71`.
- **Verify the second-order equations against a cited derivation** — restore the j² factor,
  complete δM with the exterior term, source `dε/dp` from the EOS rather than profile
  differences, impose proper central series expansions (INV-08).
- Make `HartleResult` reachable from `NStar`.
- Convergence tests and validation against published moment-of-inertia and δM values.

**Exit criteria.** O(Ω) and O(Ω²) validated and reachable. The candidate status of `675b4a9` is
resolved — ratified or replaced.

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
  ✅ RATIFIED    ✅ COMPLETE      ✅ COMPLETE       ◄ NEXT
                                    │                                 │                │              │
   ADR-0001 species semantics  ✅ ACCEPTED ──────────────────────────────────────────────────────────►│  (gate cleared)
   ADR-0002 heat capacity      ✅ ACCEPTED ─►│  (conformance is 2A work, not a gate)
   ADR thermal-balance arch.   ☐ open, deferred ────────────────────►│  (optional; not a gate)
   ADR Hartle normalization    ☐ open ─────────────────────────────────────────────────►│
   ADR η conventions           ☐ open ────────────────────────────────────────────────────────────────►│
```

**The former Phase-2 / Phase-3 circularity is gone.** It ran:

```
Phase 2 baseline  ──blocked by──►  INV-15  ──fixed in──►  Phase 3  ──requires──►  Phase 2 baseline
```

ADR-0002 breaks it by deciding the physical ownership **now**, ahead of any baseline, and by
splitting the correction out of Phase 3 into **Phase 2A**, which precedes the baseline. Nothing in
Phase 2A depends on a passive-cooling baseline: it is validated by independent physical checks.

**Two** unresolved invariants still gate the chain: **INV-07** (Hartle normalization) ·
**INV-11** (η conventions).

**INV-01** (species semantics) is **no longer a gate** — resolved by ADR-0001. What remains from
it is a single Phase-5 implementation task: `RotochemicalCache` must construct `n_i = Y_i n_B`
before its species-density integrals. That is tracked as work, not as an open decision.

**INV-15** (heat-capacity ownership) is **no longer a gate** — resolved by ADR-0002. What remains
from it is **Phase-2A work** (`PhotonCooling` conformance plus `C_⋆(T∞)` verification) and an
**optional, deferred** Phase-3 architectural question (ADR-0002 §6). Neither is an open decision
about the physics.
