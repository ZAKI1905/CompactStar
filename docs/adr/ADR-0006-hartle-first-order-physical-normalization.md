# ADR-0006 — Physical normalization and unit contract for Hartle first-order rotation

| Field | Value |
|---|---|
| **Status** | **ACCEPTED** |
| **Date accepted** | 2026-09-02 — by project-owner adjudication |
| **Date drafted** | 2026-09-02 |
| **Authority** | **Project-owner adjudication.** The owner selected Q1 = A + D, Q2 = A, Q3 = A, Q4 = A and supplied the binding no-implicit-spin clarification recorded in the Decision. The public-API convention is owner-supplied; it is not an inference from repository evidence. |
| **Drafted at** | `df859b5a73c4cac0c115f240744d89ce9f830b8d` (branch `physics/rotation-correctness`) |
| **Accepted at** | `59714a3dab0037a77b100c27192d76ca2a14d030` (acceptance commit; implementation follows separately) |
| **Change class** | **scientific-semantic** (units, physical meaning of `Ω`, `J`, `ω̄`) **and structural** (public rotation API, `HartleResult` ownership) — strictest applies (`GOVERNANCE.md` §2) |
| **Governing authority** | INV-07 (UNRESOLVED and fail-closed at drafting; **this ADR is now its normative authority**); INV-02 (unit boundary: geometric outside TOV); ADR-0005 §13 (`Find_MomInertia` side effect and `I` must stay bit-identical); Hartle (1967) ApJ 150, 1005 — validated first-order equation and exterior matching (`docs/validation/HARTLE_MOMENT_INERTIA.md`, `EQUATION MATCH`) |
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

**All four were adjudicated by the owner on 2026-09-02 exactly as recommended — see the
Decision.**

## Decision

**Adjudicated by the project owner on 2026-09-02.** The proposed contract of §4 (P1–P12) is
**ACCEPTED**, with the four owner questions of §11 decided as follows and with one binding
clarification.

| Question | Decision |
|---|---|
| **Q1 — public spin input** | **A + D.** The public scientific spin input is a **physical angular velocity `Ω` in rad s⁻¹**, represented by a **small explicit typed quantity** rather than a bare `double`. Named factories convert from spin frequency and period; no implicit Hz ↔ rad s⁻¹ conversion exists. No public entry point accepts a geometric `km⁻¹` spin. |
| **Q2 — the arbitrary seed** | **A.** The central `ω̄` seed is **strictly internal numerical normalization**. It is not a scientific parameter, gets no public setter and no public constructor argument, and must never appear in a result or an export under a physical label. |
| **Q3 — result storage** | **A.** Store **one canonical geometric representation** (`Ω [km⁻¹]`, `J [km²]`, `I [km³]`, `ω̄ [km⁻¹]`, `ω̄' [km⁻²]`) and provide **named physical accessors**. **No duplicated `Ω_geom` / `Ω_phys` state.** |
| **Q4 — exposure** | **A.** Expose **explicit seed-free normalized first-order response through `NStar`** — never the raw seed-normalized profiles, and never the private `RotationSolver` or the O(Ω²) candidate. |

### Binding clarification — no implicit physical spin

> An `NStar` constructed without an explicit physical spin does **not** acquire an implicit
> physical `Ω`.
>
> Automatic star construction may compute
>
> ```
> I ,      ω̄(r)/Ω ,      [d ω̄/dr](r)/Ω
> ```
>
> because these are **seed-free**: the arbitrary normalization cancels analytically from each of
> them.
>
> A physical
>
> ```
> Ω ,      J ,      ω̄(r) ,      dω̄/dr(r)
> ```
>
> may be materialized **only after an explicit physical angular velocity is supplied.**

This is a consequence of Q1, Q2 and Q4 taken together, **not a new independent scientific
decision**: if construction supplied a default `Ω`, the seed would re-enter the public surface
(against Q2) and a star would silently carry a spin nobody requested (against Q1).

**Scope of acceptance.** This settles the first-order normalization and unit contract. Per §9 it
ratifies **nothing** at O(Ω²), does not activate or normalize the `SolveHartle2_N` candidate, and
does not decide the cgs conversion of `J` or `I` (which needs the unadjudicated `G` / solar-mass
authority). Acceptance fixes the convention; it does not certify that any implementation meets it
— conformance is roadmap increment 4A and is tracked separately, exactly as ADR-0002's was.

---

## 12. Implementation record (Phase 4A, 2026-09-02)

**The accepted decision was not modified during implementation.** Nothing in the evidence
contradicted it, so no return to owner adjudication was required. Q1 = A + D, Q2 = A, Q3 = A,
Q4 = A and the no-implicit-spin clarification were implemented as written.

**Components changed**

| Component | Change |
|---|---|
| `CompactStar/AngularVelocity.hpp` | **new** — typed physical angular velocity (rad s⁻¹) with `FromRadPerSecond` / `FromHz` / `FromPeriodSeconds`; **no factory accepts km⁻¹**. `GeomKmInverse()` and the free `AngularVelocityGeomToRadPerSecond()` are the two conversion owners. Top level beside `Units.hpp` / `Geometry.hpp`; depends on no CompactStar component. Fails closed on non-finite input and on a zero *period*; a zero angular velocity or frequency is valid. Negative values are accepted, because `MagneticDipole` writes `−K sign(Ω)|Ω|ⁿ` |
| `Core/RotationSolver.hpp` | `HartleFirstOrderResponse` (seed-free `I`, `ω̄/Ω`, `ω̄'/Ω`) and `PhysicalFirstOrderRotation` (one canonical geometric `Ω`, `J`, `I`, profiles, `OmegaRadPerSecond()`); `HartleResult` **stripped** of `Omega`, `J`, `I`, `omega_bar`, `domega_bar` and re-scoped to the second-order candidate; private `kDefaultInitOmegaBar = 5e-3`, `seed_omega_bar_`, `first_order_solve_count_`, `first_order_response_`; public `FirstOrderResponse()`; declared-never-defined `RotationSolverTestSeam` befriended |
| `Core/src/RotationSolver.cpp` | seed read from `seed_omega_bar_` (bit-identical default); solve-completion flag; seed-free response built **after** the untouched extraction; `HartleFirstOrderResponse::At()` — the single materialization owner; `ExportResults` header tokens corrected |
| `Core/NStar.hpp`, `Core/src/NStar.cpp` | `RotationResponse()` and `RotationAt(AngularVelocity)`; `RotationSolverTestSeam` befriended. `rot_solver` stays **private** |
| `tests/rotation/hartle_normalization_contract.cpp`, `..._cmf.cpp`, `tests/CMakeLists.txt` | **new** — V1–V9 as two CTests (self-contained + external data) |

**Raw / normalized separation.** The raw `ω̄_raw`, `J_raw`, `Ω_raw` never leave the solver. The
response is `I`, `ω̄/Ω`, `ω̄'/Ω`; the physical solution is that response times `Ω_geom`, with `J`
taken from the exterior matching applied to the scaled surface derivative and `I` carried
through rather than recomputed as `J/Ω` (so zero spin is not a `0/0`).

**Seed handling.** Value unchanged at `5e-3`, private, no public setter, no constructor
argument. Reached through `RotationSolverTestSeam`, declared in production and defined by the
validation harnesses. *(Phase-4B wording correction: the accurate classification is a*
***`PRIVILEGED TEST BACKDOOR — NOT SUPPORTED SCIENTIFIC API`***; *calling it reachable only by
the harnesses is too strong in C++, since any translation unit could define the befriended type.
What Q2 requires does hold: no supported public seed setter, no supported public seed constructor
argument, and no production consumer of the seam.)*

**Unit corrections.** The `HartleResult::Omega` `[s^-1]`/km⁻¹ mislabel is gone with the field;
the legacy export header now reads `omega_bar_c_seed (1/s)`, `M (M_sun)`, `R (km)`,
`J_seednorm (km^2)`, `Omega_seednorm (1/s)`. No cgs `J`/`I` accessor was added (§10).

**Validation.** All of §7 passes, at the bounds fixed above before implementation: seed
invariance `≤ 4.3e-15` (bound `1e-10`), requested-`Ω` recovery `2.1e-16` (`1e-13`),
`J = IΩ_phys/c` `3.1e-16` (`1e-13`), conversion exact against an independent SI literal `c`
(`1e-15`), linearity `0` at every node (`1e-14`), zero spin exact, annotations and header tokens
exact. **`I` bit-identical; the seven Phase-3 artifacts byte-identical**, with
`hartle_I_dscmf1_debug.tsv` additionally **re-emitted** and compared byte-for-byte.

**Detectors.** All four of §8 fire and were reverted byte-identically (SHA-256 verified), with
the contract test green at the reverted state. D2 fires at `8.99e10` relative error — the `c²`
signature — as predicted.

**Second-order non-scope.** `ODE_Hartle2_N_Fast`, `SolveHartle2_N`, `ODE_Hartle2_Mixed_Fast` and
`SolveHartle2_Mixed` are **byte-identical** to the pre-change source. Nothing at O(Ω²) was
repaired, normalized, executed for a result or baselined; INV-08 is unchanged. P9 governs how a
*future* second-order product must relate to this contract, and that is Phase-4C work.

**Known residual, recorded not repaired.** `Solve_Mixed` retains four inline `× c` conversions
and its own seed literal on the legacy two-fluid path. P12 governs it as a contract and defers
conformance to the `MixedStar` track; only its export *labels* were corrected here.

Full record: `docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md`.

### Evidence update — Phase 4B (2026-09-02). **No decision text is changed.**

The contract this ADR fixes has since been validated as *physics*, not only as an API contract.
Phase 4B compared the seed-free response `ω̄/Ω` and `ω̄'/Ω` **node by node** against an
independently derived and independently normalized profile, and added three closed-form physical
checks. Against bounds predeclared from the Phase-2B record: profile agreement `2.9e-9` on the
exact analytic star (bound `1e-7`) and `≤ 2.3e-5` across the four authenticated CMF stars (bound
`1e-4`); the exterior identities `s(R) = 1 − 2I/R³` and `s'(R) = 6I/R⁴` against the
**independent** `I`; the volume-integral identity from production's own shape (`1.1e-7` analytic,
`≤ 3.1e-5` CMF); and the derived weak-field coefficients `ω(0)/Ω → 2(M/R)`, `ω(R)/Ω → 0.8(M/R)`.
The independent oracle's own numerical floor is four to ten decades below the difference it
measures.

**§7 item 9 (the `Ω/Ω_K` diagnostic) is explicitly NOT promoted by this.** Phase 4B records the
distinction in the terms this ADR should carry: *normalization correctness is not slow-rotation
truncation accuracy*. `Ω/Ω_K ≈ 0.6` at 716 Hz says nothing about the size of the neglected O(Ω²)
terms; §10's deferral of a governed regime threshold stands.

Phase 4B changed **no production source** and left the O(Ω²) candidate byte-identical.
Record: `docs/validation/PHASE4B_FIRST_ORDER_PHYSICS.md`.

---

## Provenance

Drafted by an AI agent under `GOVERNANCE.md` §5 and `AGENTS.md` during roadmap increment 4A-0
at `df859b5`, from the audit in `docs/validation/PHASE4_ROTATION_ENTRY.md`: the first-order
equation and extraction were re-authenticated against the Phase 2B-4B validation (source
byte-identical), the rescaling contract was derived from linear homogeneity and exterior
matching and checked numerically on two stars through the public API (`J_phys/Ω_geom − I`
within `2e-16`), the repository's spin conventions were audited file by file, and every unit
annotation on the rotation path was tabulated. The agent recommended; the agent did not select
among the alternatives.

**Adjudicated by the project owner on 2026-09-02**, who selected Q1 = A + D, Q2 = A, Q3 = A,
Q4 = A and supplied the binding no-implicit-spin clarification recorded in the Decision. Per
`GOVERNANCE.md` §5, acceptance of this ADR ratifies **only** the first-order normalization and
unit contract. It does **not** ratify the Hartle O(Ω²) or rotochemical candidate code, does not
certify that any implementation of this contract is numerically correct, and confers accepted
status on no other document.
