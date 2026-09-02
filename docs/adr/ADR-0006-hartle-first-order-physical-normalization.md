# ADR-0006 — Physical normalization and unit contract for Hartle first-order rotation

| Field | Value |
|---|---|
| **Status** | **PROPOSED** — carries no authority (`docs/adr/README.md` lifecycle) |
| **Date drafted** | 2026-09-02 |
| **Drafted at** | `df859b5a73c4cac0c115f240744d89ce9f830b8d` (branch `physics/rotation-correctness`) |
| **Change class** | **scientific-semantic** (units, physical meaning of `Ω`, `J`, `ω̄`) **and structural** (public rotation API, `HartleResult` ownership) — strictest applies (`GOVERNANCE.md` §2) |
| **Governing authority** | INV-07 (UNRESOLVED, fail-closed); INV-02 (unit boundary: geometric outside TOV); ADR-0005 §13 (`Find_MomInertia` side effect and `I` must stay bit-identical); Hartle (1967) ApJ 150, 1005 — validated first-order equation and exterior matching (`docs/validation/HARTLE_MOMENT_INERTIA.md`, `EQUATION MATCH`) |
| **Evidence record** | `docs/validation/PHASE4_ROTATION_ENTRY.md` (Phase 4A-0), §7–§9, §12–§13, §17–§18 |
| **Affected invariants** | **INV-07** (primary); INV-05 (rotation half); INV-02 (a new physical/geometric boundary); consumers INV-08, INV-09 (second order inherits the normalization) |
| **Blocks** | roadmap Phase 4A (implementation), 4B (physical validation); transitively 4C–4E and Phase 5 |
| **Explicitly does not decide** | anything at O(Ω²) — see §9; the cgs conversion of `J`, `I` (needs the unadjudicated `G`/solar-mass authority) |

---

## 1. Context

Every `NStar` construction runs `RotationSolver::FindNMomInertia()` (`NStar.cpp:335`, `:659`
→ `:1109`), which integrates the first-order Hartle frame-dragging equation from a hard-coded
seed `init_omega_bar = 5e-3` (`RotationSolver.cpp:701`) and extracts

```
J_raw = R⁴ ω̄'(R)/6 ,   Ω_raw = ω̄(R) + R ω̄'(R)/3 ,   I = J_raw / Ω_raw          (RotationSolver.cpp:737-739)
```

The equation is linear and homogeneous, so `ω̄_raw`, `J_raw`, `Ω_raw` are **arbitrary up to a
common factor**; only `I` is physical, and only `I` has been validated (Phase 2B-4B). The
arbitrary solution is nevertheless stored as if it were a result: `HartleResult::Omega` is
annotated `[s^-1]` (`RotationSolver.hpp:105`) but holds `Ω_raw` in km⁻¹ (`:744`); the export
header advertises `omega_bar_c (1/s)` and `Omega (1/s)` (`:637-638`) for seed-normalized values;
`stored_omega_bar_` feeds the O(Ω²) candidate, whose every output therefore scales as the seed
squared (entry record §12). On the audited stars the seed corresponds to `Ω_raw·c ≈ 2.2–2.4×10³
s⁻¹` — a fast millisecond pulsar — purely by accident of the number `5e-3`.

No document defines what physical spin the rotation solver represents, how a caller requests
one, or in which units first-order results are stored and exported. `GOVERNANCE.md` §3
condition 1 (ambiguous units) and condition 2 (ambiguous state meaning) are active for
`ω̄`, `Ω`, `J`. INV-07 explicitly requires an ADR. This ADR proposes the contract. **It writes no
code and authorizes no source change until accepted.**

### 1.1 Repository usage that constrains the choice (audited, not assumed)

| Where | Spin representation |
|---|---|
| `Physics/State/SpinState.hpp:65,117` | evolved spin DOF `Ω` in **rad/s** |
| `main/Test/spin_therm_evol_2_main.cpp:205` | `spin.Omega() = 100.0; // rad/s` |
| `Physics/Driver/Spin/MagneticDipole.hpp:17,50` | `Ω̇ = −K Ωⁿ` on that `Ω` |
| `Physics/Driver/Chem/src/Rotochemical.cpp:68,92,117` (candidate, not compiled) | reads `Ω`, `Ω̇` from `SpinState`; forcing `2ΩΩ̇` |
| `Core/Pulsar.hpp:272-284` | observed period `P` [s], `Ṗ` [s/s] |
| `RotationSolver::Solve_Mixed` (`:389`, `:540-542`) and the historical `Solve()` (`f6bbfea`, now undefined) | caller supplies **`ω̄_c` in km⁻¹** (a *seed*, not `Ω`) and receives `Ω·c` in s⁻¹ — `main/Examples/rotating_ns.cpp:55-65` passes log-spaced seeds `{5.66e-4, 5.66e-3}` km⁻¹ while quoting the target `Ω = 4.501e3 s⁻¹` only in a comment |
| `TOVSolver.cpp:119,179` | sequence export `I(km^3)` |

The living convention for physical spin in this repository is **`Ω` in rad s⁻¹**; the geometric
`km⁻¹` seed was only ever an internal integration variable that leaked into a legacy API.

## 2. Decision questions

1. What physical quantity does a caller supply to obtain a first-order rotating configuration,
   and in what units?
2. What is the canonical internal representation, and where exactly does the physical ↔
   geometric conversion happen?
3. Is the arbitrary numerical seed part of the scientific contract or an internal detail?
4. What do `HartleResult` (and any successor), `OmegaSeqPoint` and `ExportResults` store and
   label?
5. How is the first-order response exposed to consumers (Phase 5), and how does the contract
   propagate to O(Ω²)?

## 3. Definitions proposed (normative once accepted)

| Symbol | Definition | Unit (geometric / physical) |
|---|---|---|
| `ω(r)` | angular velocity of the local inertial frames (frame dragging) | km⁻¹ / s⁻¹ |
| `Ω` | uniform angular velocity of the fluid as seen at infinity | km⁻¹ / **rad s⁻¹ ≡ s⁻¹** |
| `ω̄(r) ≡ Ω − ω(r)` | Hartle's variable; the ODE unknown (`RotationSolver.hpp:72`) | km⁻¹ / s⁻¹ |
| `J` | total angular momentum, from exterior matching `ω̄ = Ω − 2J/r³` | km² / (g cm² s⁻¹ via `10¹⁰ c³/G`) |
| `I ≡ J/Ω` | moment of inertia (slow-rotation, `Ω → 0` limit) | km³ / (g cm² via `10¹⁵ c²/G`) |
| `c` | `Zaki::Physics::LIGHT_C_KM_S = 2.99792458e5 km s⁻¹` (exact by definition) | — |
| `Ω_geom = Ω_phys / c` | the **only** physical → geometric conversion for angular quantities | — |
| `Ω_K ≡ (M/R³)^{1/2}` | Keplerian-scale diagnostic in geometric units (no `G` needed) | km⁻¹ |
| `ω̄_raw`, `J_raw`, `Ω_raw` | the solution produced from the internal seed | geometric |
| `A ≡ Ω_geom / Ω_raw` | rescaling factor | dimensionless |

## 4. Proposed contract

- **P1 — Public input is physical.** The scientific API accepts a physical angular velocity
  `Ω_phys` in rad s⁻¹ (equivalently s⁻¹), the same quantity `SpinState::Omega()` evolves, through
  an explicitly named parameter or small typed wrapper; named helpers convert from spin frequency
  (`Ω = 2πf`) and period (`Ω = 2π/P`). **No public entry point accepts `km⁻¹`.**
- **P2 — Internal representation is geometric.** The ODE, the stored profile columns and the
  scalar results are geometric (km-based), consistent with `StarProfile` (INV-02). Exactly one
  function performs `Ω_geom = Ω_phys/c` and one performs the inverse; they are the sole owners of
  `c` on the rotation path.
- **P3 — Rescaling rule.** The physical first-order solution is `ω̄_phys = A ω̄_raw`,
  `ω̄'_phys = A ω̄'_raw`, `J_phys = A J_raw`, `Ω_geom = A Ω_raw` with `A = Ω_geom/Ω_raw`
  (derivation: entry record §8.1). `I = J_raw/Ω_raw` is computed from the raw solution and is
  **identical** before and after rescaling.
- **P4 — Seed is internal.** `init_omega_bar` is a numerical normalization: it may be retained
  (or replaced by any positive constant) but it must not be exposed as, or mistaken for, a
  physical quantity. `GetInitOmegaBar()` and the `omega_bar_c (1/s)` export column lose any
  physical reading (deprecate or relabel as an internal diagnostic). Its magnitude is
  constrained only by the integrator's tolerance behaviour (entry record §12): the seed
  invariance test (§7 item 1) is the guard.
- **P5 — Zero spin is well defined.** `Ω_phys = 0 ⇒ ω̄_phys ≡ 0`, `J_phys = 0`, `I` unchanged.
  Implementations must never form `I` as `J_phys/Ω_phys`.
- **P6 — Single canonical stored representation.** The result object stores geometric
  quantities only (`Ω [km⁻¹]`, `J [km²]`, `I [km³]`, `ω̄(r) [km⁻¹]`, `ω̄'(r) [km⁻²]`), with the
  `[s^-1]` annotation on `HartleResult::Omega` corrected and explicitly named accessors for
  physical `Ω` and `ω̄` (`× c`). **No duplicate `Ω_geom`/`Ω_phys` state.** cgs `J` and `I`
  accessors are deferred to the `G`/solar-mass authority ADR (they need `G`; `Ω` does not).
- **P7 — Export semantics.** Every exported column carries the unit of the value written;
  seed-normalized quantities are not exported under physical labels. A physical export carries
  the requested `Ω_phys`, `Ω_geom`, `J [km²]`, `I [km³]`, `M`, `R` with units in the header.
- **P8 — Response exposure.** Consumers receive **normalized response** — `I` and the
  dimensionless shape `ω̄(r)/Ω` (seed-free; `ω(r)/Ω = 1 − ω̄/Ω`) — plus the physical solution at a
  requested `Ω`. `NStar` exposes this through a public accessor of the normalized product, not
  through its private `RotationSolver`.
- **P9 — Propagation to second order.** Any O(Ω²) quantity must be either (a) computed from
  `ω̄_phys` at a requested `Ω`, or (b) exposed as a coefficient per unit `Ω_geom²` computed from
  the raw solution and divided by `Ω_raw²`. Either way it must be seed-free to the bound of §7
  item 1. This ADR does **not** decide the O(Ω²) equations (§9).
- **P10 — Regime diagnostic.** Results report `Ω/Ω_K`; the applicability threshold is a
  documented diagnostic, governed later with the O(Ω²) work.
- **P11 — `I` protection.** Nothing in this contract changes the ODE, the extraction formulas,
  or the value of `I`; the Hartle golden and both Hartle CTests remain bitwise (ADR-0005 §13).
- **P12 — MixedStar.** The same mathematics applies to `FindMixedMomInertia` and `Solve_Mixed`,
  but conformance of the two-fluid path is deferred to the MixedStar track (ADR-0004 §0-Q2
  precedent); this ADR governs it as a contract only.

## 5. Alternatives

### Q1 — public spin input

| | Option | Assessment |
|---|---|---|
| **A** | **`Ω` in rad s⁻¹ (= s⁻¹)** | matches `SpinState`, the spin drivers, the rotochemical candidate's forcing `2ΩΩ̇`, and the existing `(1/s)` export labels; the factor `2π` never enters the API. **Recommended** |
| B | spin frequency `f` in Hz | natural for observers, used nowhere in the repository; introduces the classic `2π` ambiguity at the boundary with `SpinState` |
| C | geometric `Ω` (or `ω̄_c`) in km⁻¹ | the historical `Solve()`/`Solve_Mixed` convention; silently mixes physical and geometric quantities and offers no way to request a target spin directly (the seed is not `Ω`). **Rejected** for the public API; retained only as the internal representation (P2) |
| D | a strongly typed quantity (e.g. an `AngularVelocity` value type with `FromRadPerS`, `FromHz`, `FromPeriod`, `GeomKmInv()`) | not an alternative to A but its safest **form**: it makes the unit part of the type and the conversion explicit. **Recommended together with A** |

**Recommendation: A, in form D.** Physical `Ω` [rad s⁻¹] at the API, geometric `km⁻¹` inside,
one named conversion each way.

### Q2 — the arbitrary seed

| | Option | Assessment |
|---|---|---|
| **A** | **strictly internal numerical normalization** (P4) | the equation is homogeneous; the seed carries no information. **Recommended** |
| B | expose it as an input "central `ω̄`" | perpetuates the legacy confusion between a seed and a spin; the physical `ω̄_c` follows from `Ω` and the star, not the other way round. Rejected |
| C | remove the seed by non-dimensionalizing the ODE (`ω̄/ω̄_c`) | equivalent to A with seed `1`; an implementation choice, not a contract difference |

### Q3 — `HartleResult` storage

| | Option | Assessment |
|---|---|---|
| **A** | **one canonical geometric representation + named physical accessors** (P6) | consistent with every other field and with `I = J/Ω`; one source of truth. **Recommended** |
| B | store both `Ω_geom` and `Ω_phys` | duplicate state that can drift; violates one-owner-per-quantity (`TARGET_ARCHITECTURE.md` §5 principle 2) |
| C | store physical only | breaks the geometric arithmetic of the struct and INV-02's boundary |

### Q4 — exposure to Phase 5

| | Option | Assessment |
|---|---|---|
| **A** | **explicit normalized response fields** (`ω̄/Ω`, `I`; later `Ω²` coefficients at fixed `ε_c`) through a public `NStar` accessor (P8, P9) | seed-free by construction; gives Phase 5 the derivatives it needs. **Recommended** |
| B | make `NStar::rot_solver` / raw `HartleResult` public | exposes seed-dependent profiles and the unvalidated candidate; already reachable via an external solver, which is the hazard, not the solution |
| C | require callers to run their own `RotationSolver` | today's de-facto state (entry record §13); no contract, no normalization |

## 6. Consequences (once accepted)

- INV-07 moves from UNRESOLVED to **GOVERNED (ACCEPTED)** as a contract; conformance is Phase 4A
  work and is tracked separately, exactly as INV-15 was.
- `GOVERNANCE.md` §3 conditions 1 and 2 are discharged for `ω̄`, `Ω`, `J` at first order.
- Forbidden: any public entry point taking a geometric spin; any physical label on a
  seed-normalized value; any O(Ω²) output that is not seed-free (P9); forming `I` from physical
  `J`/`Ω`.
- Required in the same change: `HartleResult` annotation and accessors, export headers,
  `CURRENT_ARCHITECTURE.md` rotation rows, INV-07 status, and the tests of §7.
- Not changed: the ODE, the extraction, `SeqPoint::I`, the seven goldens.
- Existing outputs: every historical `Omega`/`J`/`omega_bar` export from `Solve_Mixed` or the
  legacy `Solve()` is a seed-normalized quantity and is **not** a physical reference; `I`
  outputs are unaffected.

## 7. Validation requirements (predeclared; to be implemented as CTests under the Phase-1 plumbing, reusing `tests/rotation/`)

| # | Requirement | Bound |
|---|---|---|
| 1 | seed invariance: internal seed varied over `[10⁻³, 10³]`; `ω̄_phys(R)`, `J_phys`, `Ω` unchanged | rel `≤ 1e-10` |
| 2 | requested `Ω_phys` recovered from the physical surface values | rel `≤ 1e-13` |
| 3 | `J_phys = I Ω_phys / c` | rel `≤ 1e-13` |
| 4 | `I` and `hartle_I_dscmf1_debug.tsv` bit-identical; `hartle_moment_inertia_{analytic,cmf}` unchanged | bitwise |
| 5 | `ω̄_phys(r)` linear in requested `Ω`, node by node | rel `≤ 1e-14` |
| 6 | `Ω_geom · c = Ω_phys` against an independent literal `c = 299792.458 km/s` (the test must not import `LIGHT_C_KM_S`) | rel `≤ 1e-15` |
| 7 | unit annotations and exported header tokens match stored values (parsed and asserted) | exact |
| 8 | zero spin: `ω̄_phys ≡ 0`, `J_phys = 0`, `I` finite and equal to raw `I` | exact |
| 9 | regime diagnostic `Ω/Ω_K` reported for `1.0–2.0 M☉` at `Ω_phys ∈ {100, 2π·100, 2π·716} s⁻¹` | diagnostic |

Tolerances are fixed here, before implementation, from the entry record's measurements
(§8.3: `2e-16`; §12: `1e-5`–`1e-3` where the driver's absolute tolerance engages). A failure of
item 1 is to be resolved by governing the integrator tolerance, not by widening the bound.

## 8. Detectors (each applied, measured, reverted byte-identically)

| | Mutation | Must fail |
|---|---|---|
| D1 | omit `A` | items 2, 3 |
| D2 | multiply by `c` instead of dividing | item 6 (by `c² ≈ 9e10`), items 2, 3 |
| D3 | rescale `ω̄` but not `ω̄'`/`J` | items 2, 3 |
| D4 | let the seed leak into the normalized output | item 1 |

## 9. Governance boundary with second order — deliberately separate

This ADR certifies **nothing** at O(Ω²). The candidate (`RotationSolver.cpp:1024-1110`,
`:1126-1270`) has independent scientific questions — the perturbation variable (`p₀*` vs `δp`),
the `j²` factor, the fixed-central-density boundary condition versus the candidate's
`δp(R) = 0`, the exterior `J²/R³` term, the `dε/dp` authority — recorded in the entry record
§10–§12. They require **a separate ADR** (anticipated in `docs/adr/README.md`) which, because the
candidate's outputs cannot serve as a baseline, is expected to invoke `GOVERNANCE.md` §3.1 with
independent verification. Accepting ADR-0006 does not ratify, activate, or normalize the
candidate; P9 only fixes how any *future* second-order product relates to this contract.

## 10. Non-scope

Not decided here: the O(Ω²) equations and variables; the cgs conversion of `J`, `I` (`G`
authority); rotochemical semantics (INV-09, INV-11); `MixedStar` conformance; the regime
threshold; the integrator tolerance policy (`NUMERICAL_POLICY.md` does not exist —
`GOVERNANCE.md` §1); the Phase-2B coupled-convergence relocation.

## 11. Owner questions

| | Question | Recommendation |
|---|---|---|
| **Q1** | Public spin input: `Ω [rad s⁻¹]`, `f [Hz]`, `Ω_geom [km⁻¹]`, or a typed quantity? | **`Ω` in rad s⁻¹, in typed form** (5-Q1 A + D) |
| **Q2** | Does the arbitrary seed remain strictly internal? | **Yes** (5-Q2 A) |
| **Q3** | Store both `Ω_geom` and `Ω_phys`, or one canonical quantity with named accessors? | **One geometric canonical field + accessors** (5-Q3 A) |
| **Q4** | How is the first-order response exposed to Phase 5? | **Explicit normalized response fields** (5-Q4 A) |

## Decision

*Empty while PROPOSED.*

## Provenance

Drafted by an AI agent under `GOVERNANCE.md` §5 and `AGENTS.md` during roadmap increment 4A-0
at `df859b5`, from the audit in `docs/validation/PHASE4_ROTATION_ENTRY.md`: the first-order
equation and extraction were re-authenticated against the Phase 2B-4B validation (source
byte-identical), the rescaling contract was derived from linear homogeneity and exterior
matching and checked numerically on two stars through the public API (`J_phys/Ω_geom − I`
within `2e-16`), the repository's spin conventions were audited file by file, and every unit
annotation on the rotation path was tabulated. The agent recommends; **only the project owner
may move this ADR to ACCEPTED.**
