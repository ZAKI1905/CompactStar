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

**Exit criteria.**

| Criterion | Status |
|---|---|
| ADR-0001 adjudicated | ✅ **SATISFIED** — ACCEPTED 2026-08-31 |
| Owner reviews the remaining governance documents | ☐ Outstanding — `GOVERNANCE.md`, `SCIENTIFIC_INVARIANTS.md`, `AI_SKILL_PLAN.md`, and this roadmap remain **DRAFT** |

Accepting ADR-0001 ratifies the species-semantics contract only. It does not ratify any other
DRAFT document, nor the Hartle O(Ω²) / rotochemical candidate code.

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

## Phase 2 — Validation baseline

**Prerequisite:** Phase 1 complete.

The codebase has zero tests, zero assertions, and zero CI. Until baselines exist, no numerical
change can be shown correct.

- Introduce CTest and a test framework.
- TOV reference checks against known solutions.
- First-order Hartle moment-of-inertia checks against published values.
- Grid-convergence harness — noting INV-13: interpolation is **linear**, so expect
  second-order-in-Δr behavior, not fourth.
- Passive cooling regression capturing today's behavior as a baseline.

**Blocked by:** INV-15 (heat-capacity ownership) — a thermal baseline captured while two
different `C` values feed one energy balance would enshrine the inconsistency.

**Exit criteria.** Baselines exist and run; a regression is detectable.

---

## Phase 3 — Behavior-preserving consolidation

**Prerequisite:** Phase 2 baselines exist.

Every item is **engineering class** and must produce bit-identical output, or a documented
tolerance.

- Establish the canonical TOV path; retire or clearly subordinate the duplicate.
- Single owner for unit conversions — remove the ~15 local re-derivations, including **k_B at
  two precisions**.
- Single owner for the proper-volume measure (INV-04).
- One uniform cache-invalidation rule; add a version gate to `GeometryCache`; re-bind
  `StarContext` column pointers on invalidation (INV-12).
- Single owner for heat capacity (INV-15) — **requires an ADR**.
- Classify dead and unreachable code; retire only after dependency review.

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

**Prerequisites:** Phase 4 · **ADR-0001 accepted ✅** · ADR on η conventions accepted ☐.

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
0.5 governance ─► 1 build ─► 2 baseline ─► 3 consolidation ─► 4 rotation ─► 5 rotochemical ─► 6 BNV
                                               │                  │              │
   ADR-0001 species semantics  ✅ ACCEPTED ────────────────────────────────────►│  (gate cleared)
   ADR heat-capacity           ☐ open ───────►│                                  │
   ADR Hartle normalization    ☐ open ──────────────────────────►│               │
   ADR η conventions           ☐ open ──────────────────────────────────────────►│
```

**Three** unresolved invariants still gate the chain: **INV-07** (Hartle normalization) ·
**INV-11** (η conventions) · **INV-15** (heat-capacity ownership).

**INV-01** (species semantics) is **no longer a gate** — resolved by ADR-0001. What remains from
it is a single Phase-5 implementation task: `RotochemicalCache` must construct `n_i = Y_i n_B`
before its species-density integrals. That is tracked as work, not as an open decision.
