# Phase 2B closure audit — obligations, CI, deferred convergence, merge readiness

> **THREE SEPARATE ANSWERS. They are deliberately not collapsed into one.**
>
> | Question | Answer |
> |---|---|
> | Is Phase 2B complete? | **`PHASE 2B ACTIVE — MERGE GATE SATISFIED`** |
> | Is the branch merge-ready? | **`PHASE-2 BRANCH MERGE READY`** |
> | Is Phase 3's prerequisite met? | **`PHASE 3 PREREQUISITE SATISFIED`** |
>
> Phase 2B is **not complete** — CI is absent and the full convergence item is `PARTIAL` — yet
> the ratified **exit criterion** ("Baselines exist and run; a regression is detectable") and the
> ratified **Phase-3 prerequisite** ("Phase 2B baselines exist") are both independently
> satisfied. Nothing about the branch's content blocks a fast-forward.
>
> **This is an audit. No production source was modified.** The only changes are this record, a
> stale-fact correction in `CURRENT_ARCHITECTURE.md`, a corrected figure in
> `HARTLE_MOMENT_INERTIA.md`, and roadmap status clarifications.

| Field | Value |
|---|---|
| **Starting commit** | `10c397d770f864ff20bd1c25a81249966de3fcbb` |
| **Canonical `master`** | `e2cc9270f5bdae677e659233398c360cfbb318a8` |
| **Relationship** | 10 ahead, 0 behind, 0 merges, upstream equal, clean tree |
| **Merge base** | `e2cc927…` (= `master`, so fast-forward is possible) |
| **Change class** | documentation / governance sync |
| **Governing authority** | `GOVERNANCE.md` → ACCEPTED ADRs → `SCIENTIFIC_INVARIANTS.md` → `MODERNIZATION_ROADMAP.md` |

---

## 1. Test results, both configurations

| Configuration | Registered | Result |
|---|---|---|
| **A** — full, with authenticated `COMPACTSTAR_EOS_DATA_ROOT` | 13 | **13/13 PASS** |
| **B** — clean checkout, no data root | 8 | **8/8 PASS** (13.9 s) |

The five external-data tests are **not registered at all** without the data root — the CMake
guard at `tests/CMakeLists.txt:71` excludes them; they do not appear as CTest *skips*:
`heat_capacity_real_star`, `passive_cooling_regression`, `tov_reference_cmf`,
`grid_convergence_cmf`, `hartle_moment_inertia_cmf`.

**A clean checkout therefore reproduces no external-data science.** Configuration B validates
units, analytic limits, cache contracts and the exact Schwarzschild TOV right-hand side — but
every check against the authenticated CompOSE authority (the passive-cooling golden curve, the
`eos.mr` comparison, the real-star heat capacity, grid convergence, the CMF Hartle sequence)
is absent. **An external-data exclusion is not a scientific PASS.**

## 2. Binding Phase-2B obligation matrix

Roadmap Phase 2B, exact item list at `docs/MODERNIZATION_ROADMAP.md:199–290`.

| # | Roadmap item (exact) | Status | Evidence | Blocks 2B completion? | Blocks merge? | Blocks Phase 3? | Remaining work |
|---|---|---|---|---|---|---|---|
| 1a | "Expand the Phase-1 test plumbing into a full harness…" | **DONE** (unticked in text) | 1 → 13 CTests across Phase 2A–2B | no | no | no | tick the item |
| 1b | "…**add CI**." | **NOT DONE** | no `.github/`, no CI config of any kind | **YES** | **no** | **no** | minimal macOS CI (§5) |
| 2 | "TOV reference checks against known solutions." | ✅ COMPLETE (2B-2) | `TOV_REFERENCE.md`; `3.5e-16` analytic, `2.8e-4` `M_max` | no | no | no | — |
| 3 | "First-order Hartle moment-of-inertia checks against published values." | ✅ COMPLETE (2B-4B) | `HARTLE_MOMENT_INERTIA.md`; `HARTLE-I VERIFIED` scale-free | no | no | no | — |
| 4 | "Grid-convergence harness … (TOV → Hartle → cooling)" | ◐ **PARTIAL** (2B-4A) | `GRID_CONVERGENCE.md`; nonrotating slice `CHARACTERIZED` | **YES** | **no** | **no** | not performable in Phase 2B — §6 |
| 5 | "Cache-correctness checks (INV-12)." | ✅ COMPLETE (2B-3) | `CACHE_CORRECTNESS.md` | no | no | no | Phase-3 repair |
| 6 | "Passive cooling regression." | ✅ COMPLETE (2B-1R) | `PASSIVE_COOLING_BASELINE.md`; §3.1 cond. 7 | no | no | no | — |
| 7 | **Exit criteria:** "Baselines exist and run; a regression is detectable." | ✅ **SATISFIED** | 13/13; four independent detector proofs | — | — | **this is the gate** | annotate |

A green check is used only where durable evidence exists — never because code is present.

## 3. Why "Phase 2B complete" and "exit criteria satisfied" are different things

The roadmap demonstrates its own convention twice. **Phase 1** (`:72`) and **Phase 2A** (`:135`)
are each marked `✅ COMPLETE`, and in both cases **every listed item carries a ✅ *and* the exit
criteria carry an explicit `✅ SATISFIED` annotation** (`:114`, `:180`). Phase 2B has neither:
two items are unticked and `:290` carries no annotation.

Read as a whole, the roadmap therefore treats a phase as COMPLETE when its items are discharged,
while using the **Exit criteria** line as the substantive gate that downstream phases key off.
That is exactly how Phase 3's prerequisite is written — "Phase 2B **baselines exist**"
(`:296`), not "Phase 2B complete".

## 4. Does a TOV → Hartle → cooling coupled observable exist today? **No.**

Established from source, not inferred:

- `DriverContext` (`CompactStar/Physics/Evolution/DriverContext.hpp:106,117,130,141,148`) exposes
  exactly five pointers to drivers: `star`, `geo`, `envelope`, `thermo`, `cfg`. **No Hartle
  quantity is reachable by any driver.**
- `grep` over `CompactStar/Physics/Driver/Thermal/` for `MomI`, `MomInertia`, `GetHartleResult`,
  `hartle` returns **nothing**.
- The canonical baseline runs `couple_spin = false`, `n_eta = 0`, drivers = {PhotonCooling,
  NeutrinoCooling} (`tests/thermal/passive_cooling_regression.cpp:158,203-204,259-260`).
- `SeqPoint::I` *is* computed on every star build (`NStar::BuildFromTOV → Find_MomInertia`), but
  nothing downstream of it reads it in the thermal path.

**Answer to the concrete question — "What Hartle output enters the present passive-cooling RHS?"
— none.** Phase 2B-4A measured TOV → cooling and Phase 2B-4B measured the scale-free `I`; those
are two disjoint validations, and this audit explicitly declines to present them as a coupled
observable.

## 5. CI adjudication — **`CI IS A PHASE-2B COMPLETION BLOCKER`**

CI is a **listed Phase-2B item** (`:199`), it is **absent** (no `.github/`, no `.gitlab-ci.yml`,
no equivalent), and the phase's own preamble names the problem it exists to fix: *"The codebase
has zero assertions and zero CI."* Under the convention of §3, an unticked item bars the COMPLETE
label.

It is **not** a merge blocker and **not** a Phase-3 blocker: neither the exit criterion nor the
Phase-3 prerequisite mentions CI, and both are satisfied without it.

### Minimum viable CI — scope only, deliberately not implemented

| Constraint | Finding |
|---|---|
| Vendored `Zaki`/`Confind` | static archives of unproven provenance for any non-macOS-arm64 target; source repos exist but are **not** demonstrated to reproduce them |
| Toolchain | CMake ≥ 4.x, AppleClang, GSL 2.7.1 (MacPorts `/opt/local`), OpenMP |
| Python | `find_package(Python3 … NumPy REQUIRED)`; needs `-DPython3_EXECUTABLE` pinned (Phase 1A finding) |
| Runner | a GitHub-hosted **macOS** runner is the only plausible host; GSL/OpenMP would need Homebrew or MacPorts provisioning |
| Data | the 5 external-data tests need ~220 MB of authenticated CompOSE tables that are **not** in the repository and have no reproducible provisioning path |

Two distinct things, which must not be conflated:

- **Self-contained CI** — configure, build, run the 8 data-free CTests. Achievable now. This is
  what "add CI" should mean for Phase 2B.
- **Full external-data scientific CI** — would additionally need the CompOSE tables provisioned
  reproducibly *with their recorded SHA-256s verified*. **No such provisioning path exists**, so
  this cannot be claimed today.

## 6. Disposition of the full convergence item — **`OWNER ADJUDICATION REQUIRED`**

The item's remaining scope is a coupled `TOV → Hartle → cooling` observable. §4 shows no such
coupling exists in production. Making it exist requires, at minimum, one of:

| Ingredient | Where it lives | Status |
|---|---|---|
| scale-free `I` | Phase 2B-4B | ✅ done — but inert for cooling |
| physical `Ω`/`J` normalization | **Phase 4** (`:325` — "ADR on Hartle normalization accepted") | INV-07 **UNRESOLVED** |
| O(Ω²) structural response | **Phase 4** | INV-08 unverified, candidate unratified |
| rotochemical spin-compression heating | **Phase 5** | blocked on Phase 4 + η ADR |

**§9 classification: `ROADMAP DEPENDENCY INCONSISTENCY`.** A Phase-2B item requires physics that
Phase 4 and Phase 5 must first define and validate. It is *not* a circular blockage, because
Phase 3 keys off "baselines exist" rather than Phase-2B completion — the ordering still runs.

**§10 disposition: `OWNER ADJUDICATION REQUIRED`.** Relocating a scientifically meaningful
requirement from Phase 2B to Phase 4/5 changes ratified meaning, which this audit has no
authority to do. **The roadmap item is therefore left unchanged.** Proposed wording, for the
owner to accept or reject:

> ◐ **Grid-convergence harness — PARTIAL (2B-4A).** The nonrotating TOV → cooling slice is
> measured and CHARACTERIZED. **The remaining `TOV → Hartle → cooling` scope is RELOCATED TO
> PHASE 4** and re-scoped as: *"convergence of the coupled observable, once a physically
> normalized Hartle quantity actually enters the thermal right-hand side."* Rationale: no Hartle
> output reaches the thermal RHS today (`DriverContext` exposes no Hartle field), and none can
> until INV-07 is resolved. Phase 2B retains no outstanding convergence obligation.

## 7. Ten-commit audit (`e2cc927..10c397d`)

| Commit | Increment | Class | Production changed? | Evidence | Merge concern |
|---|---|---|---|---|---|
| `f467b01` | 2A-1 | verification | no | `HEAT_CAPACITY_V1.md` | none |
| `a7beff4` | 2A-2 | verification | no | `HEAT_CAPACITY_V1.md` Tier B | none |
| `56aee7c` | 2A-3 | **scientific-semantic** | **YES** — PhotonCooling ×3, NeutrinoCooling_Details, 2 `main/` programs | ADR-0002 §3.1 chain; `HEAT_CAPACITY_V1.md` **V1 VERIFIED**; `photon_cooling_conformance` CTest | **cleared** — §8 |
| `800f245` | 2B-1 | verification (blocked) | no | honest BLOCKED record | none |
| `a2c1f62` | 2B-1R | **numerical-method + engineering** | **YES** — `EvolutionConfig.hpp`, `GSLIntegrator.cpp`, `EvolutionSystem.cpp` | `PASSIVE_COOLING_BASELINE.md`; GSL contract cited; `evolution_stepper_contract` CTest | **cleared** — §8 |
| `b8a26cb` | 2B-1R | verification | no | golden baseline + detector proof | none |
| `fc0f1ae` | 2B-2 | verification | no | `TOV_REFERENCE.md` | none |
| `e1d369f` | 2B-3 | verification | no | `CACHE_CORRECTNESS.md` + 3 detectors | none |
| `a24fe95` | 2B-4A | verification | no | `GRID_CONVERGENCE.md` | none |
| `10c397d` | 2B-4B | verification | no | `HARTLE_MOMENT_INERTIA.md` + 3 detectors | none |

**Eight of ten commits changed no production source at all.**

### 7.1 The two production-changing commits

**`56aee7c` — PhotonCooling conformance + NeutrinoCooling null-check.** Scientific-semantic.
Authority: ACCEPTED ADR-0002 via the `GOVERNANCE.md` §3.1 exception (§8 below). Independent
evidence that does not depend on the superseded output: `HEAT_CAPACITY_V1.md` **V1 VERIFIED**
(unit chain, analytic trapezoid, low-`T` scaling, radial order 2.000 at a cache node). Durable
conformance test: `tests/thermal/photon_cooling_conformance.cpp`. The `main/` edits are the
necessary consequence of deleting `PhotonCooling::Options::C_eff`.

**`a2c1f62` — stepper contract.** Numerical-method (default stepper `MSBDF` → `RKF45`) plus
engineering (fail-closed guard; removal of a dead block that made a Spin state mandatory).
Rationale recorded in the change and in `PASSIVE_COOLING_BASELINE.md`, citing the GSL 2.7.1
contract (`gsl_odeiv2.h:48`, `:69` — `GSL_ODEIV_JA_EVAL` dereferences a null jacobian with no
guard). Covered durably by `evolution_stepper_contract`. `NUMERICAL_POLICY.md` does not exist, and
`GOVERNANCE.md` §2 expressly defers that citation "once that document exists".

**Neither rests on "tests pass now."** Both carry independent physical or documented-contract
evidence predating any baseline.

## 8. `GOVERNANCE.md` §3.1 chain — all seven conditions SATISFIED

| # | Condition | Verdict |
|---|---|---|
| 1 | ACCEPTED ADR identifies a specific invalid behavior | ✅ ADR-0002 — `PhotonCooling` dividing by a driver-local `C_eff = 1e40` |
| 2 | Capturing it would enshrine the rejected behavior | ✅ stated in ADR-0002 and INV-15 |
| 3 | ADR explicitly identifies the minimum correction | ✅ divide by the canonical `C_⋆(T∞)`; **Pattern A preserved, Pattern B deferred** |
| 4 | Independent verification not depending on superseded output | ✅ `HEAT_CAPACITY_V1.md` — analytic limits, unit chain, closed-form quadrature |
| 5 | Narrowly scoped, no adjacent cleanup | ✅ only PhotonCooling + the null-check the roadmap explicitly authorized as "the immediately adjacent safety defect required to exercise that path", plus the two `main/` programs forced by the option's removal |
| 6 | Records which historical outputs are unsuitable as reference | ✅ ADR-0002, INV-15, `CURRENT_ARCHITECTURE.md:194` — retained, not deleted |
| 7 | Baseline created immediately after | ✅ `b8a26cb`, `passive_cooling_cmf_1p6_debug.tsv` — the two-commit gap is the honest BLOCKED report and its repair |

## 9. ADR-0002 post-acceptance diff audit — **no normative change**

Diff vs `master`: **+18 lines, 0 deletions**, all inserted at one point in the Provenance region.

The **Decision** section is **byte-identical** to `master`
(`sha256 48bafa2b838e204526dc7cd60d7792e45ff8dbc4799a9a800a7820a9dc74fcfe` on both sides), and the
status header remains `**ACCEPTED**` / `2026-08-31` on both sides.

Classification of every added line: **implementation-status / provenance update.** The addition
states explicitly that "The Context, Evidence, Rejected-current-behavior and Consequences sections
above are **left unchanged as the historical record**". **No post-acceptance normative semantic
change exists.** Not a merge blocker.

## 10. Invariant-status audit — no silent promotion

`diff` of every `## INV-` headline against `master`: **identical**. The file changed only inside
INV-15's sub-status rows (+17 / −12), each traceable:

| Change | Authority |
|---|---|
| INV-15 *Source conformance* ❌ → ✅ | commit `56aee7c` + `photon_cooling_conformance` |
| INV-15 *Numerical validation (V1)* ❌ → ✅ | `HEAT_CAPACITY_V1.md` **V1 VERIFIED** |
| INV-15 new row: *passive-cooling baseline* ✅ | `PASSIVE_COOLING_BASELINE.md`, §3.1 cond. 7 |
| Null-check hazard → ✅ CORRECTED | `56aee7c`, engineering class |

Verified still unresolved / unratified: **INV-07 UNRESOLVED**, **INV-11 unresolved**, **INV-12
hazards recorded honestly** (Phase 2B-3 deliberately did not alter its status — the four
identity/geometry/pointer defects remain), **INV-08 O(Ω²) unverified and unratified**,
**rotochemical candidate unratified**. The added INV-15 text itself preserves the limitation:
*"Regression only — the neutrino `Q0` normalizations remain placeholders."*

## 11. Branch validation matrix

| Domain | Reference / baseline | Independent? | CTest | External data | Detector shown | Status | Known limitations |
|---|---|---|---|---|---|---|---|
| Build smoke | real out-of-line symbol | n/a | `compactstar_library_smoke` | no | n/a | ✅ | asserts no science |
| Heat capacity (Tier A) | closed-form trapezoid, analytic limits | ✅ | `heat_capacity_v1` | no | ✅ | ✅ V1 VERIFIED | synthetic fixture |
| Heat capacity (real star) | authenticated CMF tables | ✅ | `heat_capacity_real_star` | **yes** | — | ✅ | magnitude-level |
| Photon conformance | ADR-0002 Pattern A | n/a | `photon_cooling_conformance` | no | — | ✅ conformant | not a physics validation |
| Evolution stepper | GSL 2.7.1 documented contract | ✅ | `evolution_stepper_contract` | no | ✅ | ✅ | — |
| Passive cooling | own golden curve | ✗ (regression) | `passive_cooling_regression` | **yes** | ✅ 1000× | ✅ baseline | `Q0_DU`/`Q0_MU` are **placeholders** |
| TOV analytic | exact Schwarzschild interior | ✅ | `tov_reference_analytic` | no | ✅ (GR vs Newtonian) | ✅ `3.5e-16` | — |
| TOV CMF | official CompOSE `eos.mr` | ✅ | `tov_reference_cmf` | **yes** | ✅ | ✅ | 0.20–0.35 % radius offset, attributed |
| Cache contracts | analytic rebuild factors | ✅ | `cache_contract`, `cache_thermal_contract` | no | ✅ ×3 | ✅ supported only | 5 hazards remain, Phase 3 |
| Grid convergence | self-refinement + Richardson | ✅ | `grid_convergence_cmf` | **yes** | ✅ coarse-grid | ◐ CHARACTERIZED | order ≈0.70, **not** O(Δr²) |
| Hartle-I analytic | exact interior + Newtonian limit | ✅ | `hartle_moment_inertia_analytic` | no | ✅ ×3 | ✅ VERIFIED | scale-free only |
| Hartle-I CMF | independent conservative solver | ✅ | `hartle_moment_inertia_cmf` | **yes** | — | ✅ VERIFIED | scale-free only |

## 12. Limitations that must survive the merge — recorded, not smoothed

**TOV.** The surface is the **EOS table floor**, not vacuum (`n_B = 1e-7 fm⁻³`, still inside the
outer crust). Radii run 0.20–0.35 % below the official reference; the residual is a constant
0.46–0.49 fraction of the hydrostatically estimated omitted layer at every mass. Default
`r_max = 70 km` / `radial_res = 10000` leaves ~80 % of the grid outside the star.

**Grid convergence.** Measured order **≈0.70** for `R`, `e^ν(R)`, `L_ν`; **0.742** for the
trajectory norm — **not** the nominal `O(Δr²)`. The order **drifts downward** with refinement.
Default-grid continuum error `T_inf ≈ 4.5e-4`, `R ≈ 5.2e-4`. **No governed downstream accuracy
budget exists**; adequacy is **unresolved**, and the regression tolerance `1e-5` is explicitly
*not* a continuum-accuracy claim.

**Cache.** Same-star contracts verified; **five hazards remain** — unversioned `GeometryCache`
(51.6 % divergence, unaskable), `C_star` key omitting geometry (50 %), version-only generic key
(85.7 %), cross-star `NeutrinoCooling` collision (80 %), unrebound column pointers. The canonical
baseline provably reaches none of them; **that is a claim about one procedure, not the API.**
Repair is Phase 3.

**Hartle.** `I` validated **scale-free only**. Physical `Ω`/`J` and the `init_omega_bar` seed
remain **INV-07 UNRESOLVED**; `HartleResult::Omega` is annotated `[s⁻¹]` but stores `km⁻¹`.
**No O(Ω²) claim.**

**Thermal.** `Q0_DU = 1e27`, `Q0_MU = 1e21` remain **self-described placeholders**. The golden
cooling curve is **regression evidence, not microphysical validation**.

## 13. Unit-constant inconsistency — **`PHASE-3 UNIT-CONSOLIDATION DEBT`**

Phase 2B-4B recorded this; the audit **corrects the figure**.

- `Zaki::Physics::SUN_M_KM = 1.476625038050 km` is **exactly the IAU nominal `GM_sun/c²`**
  (`1.32712440018e26 / c² = 1.476625038 km`) — the accurate, standard constant. It converts
  `M_sun → km` in `NStar::BuildFromTOV`, which sets the geometric `m [km]` column and hence `λ`.
- `GSL_CONST_CGSM_SOLAR_MASS = 1.98892e33 g` with GSL's own `G = 6.673e-8` gives
  `GM_sun/c² = 1.476716 km`. `TOVSolver` divides `m [g]` by that `M_sun` to report solar masses.
- **Internally consistent discrepancy: `6.2e-5`** — not the `2.8e-4` first reported, which came
  from comparing against CODATA `G = 6.67430e-8` rather than the `G` the build uses. It is the
  familiar consequence of `G` being far less precisely known than `GM_sun`.

**Affects:** the `g ↔ km` round trip, hence `m [km]` and `λ` at the `~6e-5` level — below every
tolerance currently in use and below the measured grid-convergence spread (`6.6e-5`).
**Does not affect:** any comparison performed on this branch. The TOV integration is entirely CGS
with GSL's `G`, `c`; the `eos.mr` comparison used CGS quantities; the Hartle production/reference
comparison used the same profile columns on both sides. **It is pre-existing on `master` and is
not introduced by this branch** (no commit here touches `TOVSolver` or `NStar`).

Phase 3 already carries "Single owner for unit conversions — remove the ~15 local re-derivations"
(`:298`). This belongs there. **Not repaired here.**

## 14. Merge mechanics

| Check | Result |
|---|---|
| `master` is an ancestor of HEAD | ✅ |
| fast-forward possible (0 behind) | ✅ |
| `git merge-tree` conflicts | **0** |
| local == remote | ✅ |
| working tree clean, no untracked files | ✅ |
| `git diff --check` | ✅ clean |
| external EOS files tracked | **0** |
| detector-mutation residue in production | **none** |
| generated output added by this branch | **none** |

**One pre-existing condition, not introduced here:** `build_xcode/` is tracked on `master`
(123 files). `.gitignore:4-7` documents it explicitly — `/build/` and `/build-*/` are anchored and
deliberately do not match `build_xcode` (hyphen vs underscore). **Zero of those files are added or
modified by this branch.** Under `GOVERNANCE.md` the "generated artifact" class requires an ADR on
what may be tracked; that remains outstanding **pre-existing debt**, not a Phase-2 merge blocker.

## 15. Final classifications

| | |
|---|---|
| **Phase 2B formal status** | **`PHASE 2B ACTIVE — MERGE GATE SATISFIED`** |
| **Branch readiness** | **`PHASE-2 BRANCH MERGE READY`** |
| **Phase 3 prerequisite** | **`PHASE 3 PREREQUISITE SATISFIED`** — ratified wording is "Phase 2B baselines exist" (`:296`); they exist and run |
| **CI** | **`CI REQUIRED BEFORE PHASE-2B COMPLETE BUT NOT BEFORE MERGE`** |
| **Full convergence item** | **`OWNER ADJUDICATION REQUIRED`** — roadmap left unchanged; proposed wording in §6 |
| **Unit-constant gap** | **`PHASE-3 UNIT-CONSOLIDATION DEBT`** |

No new ADR is required by this audit. The one item that would need owner authority — relocating
the convergence scope — is reported, not enacted.

## 16. Recommended next task

**The governed fast-forward merge of `modernization/thermal-prebaseline` into `master`.**

Every merge-safety and Phase-3 prerequisite is satisfied, all 10 commits are evidenced, the
accepted ADR is unchanged in substance, no invariant was silently promoted, and both test
configurations are green. The two outstanding Phase-2B items (CI, and the convergence scope
awaiting owner adjudication) do not bear on merge safety and are better resolved on `master`.
