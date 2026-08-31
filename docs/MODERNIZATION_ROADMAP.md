# CompactStar Modernization Roadmap

> **STATUS: DRAFT.** Ordered by **dependency, not by date.** A phase begins only when its
> prerequisites hold.
>
> Governing principle: **do not propose scientific implementation before its prerequisites are
> governed and validated.**

---

## Phase 0.5 — Source-of-truth reconciliation and governance foundation

**Status: this change.**

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

**Exit criteria.**

| Criterion | Status |
|---|---|
| ADR-0001 adjudicated | ✅ **SATISFIED** — ACCEPTED 2026-08-31 |
| ADR-0002 adjudicated | ✅ **SATISFIED** — ACCEPTED 2026-08-31 |
| Owner reviews the remaining governance documents | ☐ Outstanding — `GOVERNANCE.md`, `SCIENTIFIC_INVARIANTS.md`, `AI_SKILL_PLAN.md`, and this roadmap remain **DRAFT** |

Accepting ADR-0001 ratifies the species-semantics contract only; accepting ADR-0002 ratifies the
heat-capacity ownership convention only. Neither ratifies any other DRAFT document, nor the
Hartle O(Ω²) / rotochemical candidate code, nor the numerical correctness of
`StarContext::HeatCapacityStar_Tinf`.

---

## Phase 1 — Reproducible macOS build

**Prerequisite:** Phase 0.5 reviewed.

Nothing below Phase 1 can be verified, because the project currently cannot be configured from a
clean clone on any platform.

- Guard the seven absent optional `main/` subdirectories so a clean clone configures.
- Stop writing generated configuration into the source tree — move `configure_file` output to
  the binary directory.
- Define and document canonical macOS configure/build commands.
- Record exact dependency versions: GSL, OpenMP, Python3/NumPy, and the vendored Zaki and
  Confind archives.
- Adopt a warning policy and a default build type. *(Engineering class — must not change results.)*

**Exit criteria.** A clean clone configures and builds on macOS with documented commands.

**Decision gate — platform.** **Mac-first development is acceptable initially.** Cross-platform
dependency work moves earlier **only if** Linux or cloud compilation becomes required, CI must
run on Linux, or clean macOS dependency reproduction proves impossible. **A full Zaki or
Confind modernization is explicitly out of scope for now.**

---

## Phase 2A — Pre-baseline correctness prerequisites

**Prerequisite:** Phase 1 complete — the project builds reproducibly.

**This phase is deliberately narrow.** It exists for one reason: a small number of known defects
would make a validation baseline *scientifically misleading* if frozen into it. Those, and only
those, are corrected first.

**This is not a general cleanup phase.** Behavior-preserving consolidation stays in Phase 3;
unrelated scientific corrections stay in their own phases. Admission requires showing that a
baseline captured *without* the correction would encode physics known in advance to be wrong.

- **Replace `PhotonCooling`'s production use of the constant `C_eff` with the governed
  `C_⋆(T∞)`** (ADR-0002; INV-15). `PhotonCooling_Details.cpp:320` must divide by the canonical
  stellar heat capacity, not `drv.GetOptions().C_eff`
  (`PhotonCooling.hpp:229`). Update the driver's Doxygen (`PhotonCooling.hpp:55-62`, `:120-123`,
  `:214-229`; `PhotonCooling.cpp:27,36`) in the same change, and stop hand-setting `C` in
  `spin_therm_evol_2_main.cpp:245` and `spin_therm_evol_main.cpp:178`.
- **Correct the immediately adjacent safety defect required to exercise that path** — the
  confirmed null-check ordering issue in `NeutrinoCooling_Details.cpp`, which dereferences
  `ctx.star` at `:889` while its guard sits at `:901`. Engineering class, tracked separately from
  the INV-15 decision, admitted here only because routing a second driver through the same
  context path exercises the same unguarded pattern.
- **Add narrowly targeted verification of the stellar heat-capacity calculation itself** —
  ADR-0002 §V1: dimensional check through `KM3_TO_CM3`; the degenerate `c_V ∝ T` low-temperature
  slope; order-of-magnitude comparison against published total heat capacities; second-order
  convergence in Δr (INV-13); insensitivity to the `NT = 160` temperature grid; explicit
  statement of the endpoint-clamping behavior (`StarContext.cpp:725-728`); cache-rebuild
  correctness (INV-12).

**Evidence standard.** Phase-2A changes are governed by their own scientific/numerical class
under `GOVERNANCE.md` §2 and **require independent physical checks** — analytic limits, dimensional
analysis, convergence, and comparison against published values. They **must not** be validated by
comparison against the existing passive-cooling curve, which encodes the behavior ADR-0002
rejects. That is precisely the circularity this split removes.

**Exit criteria.** The thermal energy equation uses one heat capacity; `C_⋆(T∞)` has passed
independent verification; no known scientifically misleading defect remains on the thermal path.

---

## Phase 2B — Validation baseline

**Prerequisite:** Phase 2A complete.

The codebase has zero tests, zero assertions, and zero CI. Until baselines exist, no numerical
change can be shown correct.

- Introduce CTest and a test framework.
- TOV reference checks against known solutions.
- First-order Hartle moment-of-inertia checks against published values.
- Grid-convergence harness — noting INV-13: interpolation is **linear**, so expect
  second-order-in-Δr behavior, not fourth.
- Cache-correctness checks (INV-12).
- **Passive cooling regression.** Because Phase 2A has landed, this now captures a **physically
  coherent energy equation** — `C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞` — rather than deliberately
  preserving a known placeholder. It remains a regression baseline, not a physics validation: the
  neutrino emissivity normalizations are still self-labeled placeholders
  (`NeutrinoCooling_Details.cpp:100-102`).

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
