# Phase 4A — first-order physical normalization: implementation and validation

> **FORMAL STATUS: `PHASE-4A FIRST-ORDER PHYSICAL NORMALIZATION COMPLETE`**
>
> `ADR-0006` is **ACCEPTED** and the source now conforms. A star no longer carries an arbitrary
> seed-normalized angular velocity behind a `[s^-1]` label; it carries a **seed-free response**,
> and a **physical** first-order solution exists only where a caller asked for one.
>
> **`I` is bit-identical and all seven Phase-3 artifacts are byte-identical.** The O(Ω²)
> candidate is **untouched**: all four of its functions are byte-identical to the pre-change
> source, and nothing here validates, activates or normalizes it (INV-08 unchanged).

| Field | Value |
|---|---|
| **Starting HEAD** | `59714a3dab0037a77b100c27192d76ca2a14d030` (Phase 4A-0 entry audit) |
| **Canonical `master`** | `df859b5a73c4cac0c115f240744d89ce9f830b8d` — branch 1 ahead / 0 behind at entry, upstream equal, clean tree |
| **Acceptance commit** | `a9c423c40b44686cf90623fd69af4fad39a8b46e` — *docs: accept Hartle first-order normalization contract* |
| **Change class** | **scientific-semantic** (units, the meaning of `Ω`, `J`, `ω̄`) **and structural** (public rotation API) — strictest applies (`GOVERNANCE.md` §2) |
| **Governing authority** | **ADR-0006 (ACCEPTED 2026-09-02)**, Q1 = A + D, Q2 = A, Q3 = A, Q4 = A, plus the binding no-implicit-spin clarification; INV-07; INV-02; ADR-0005 §13 (`I` stays bit-identical) |
| **Evidence it builds on** | `docs/validation/PHASE4_ROTATION_ENTRY.md` (4A-0); `docs/validation/HARTLE_MOMENT_INERTIA.md` (2B-4B, `EQUATION MATCH`, `HARTLE-I VERIFIED`) |
| **Build / toolchain** | `Debug`, AppleClang 17, CMake 4.2, GSL 2.7.1, macOS arm64 |

---

## 1. What changed, in one paragraph

The frame-dragging equation is linear and homogeneous, so its solution is fixed only up to a
constant. Production seeded it with `init_omega_bar = 5e-3` and then published the resulting
`Ω_raw` and `J_raw` as if they were physical — `HartleResult::Omega` was annotated `[s^-1]` while
storing km⁻¹, and on the audited stars that seed happened to correspond to a 350–385 Hz pulsar
that nobody had requested. Phase 4A separates the three things that were conflated: the **raw
internal solution** (an implementation detail), the **seed-free response** of the star
(`I`, `ω̄/Ω`, `ω̄'/Ω` — what construction may legitimately produce), and a **physical solution**
at an explicitly requested angular velocity (what a caller must ask for). The arbitrary seed no
longer reaches any public surface, and the physical spin input is now a typed quantity in
rad s⁻¹ rather than a bare `double` whose unit lived in a comment.

## 2. Owner decisions implemented

| Question | Decision | Where it lands |
|---|---|---|
| **Q1** public spin input | **A + D** — physical `Ω` in rad s⁻¹, carried by an explicit typed quantity | `CompactStar/AngularVelocity.hpp`; `NStar::RotationAt(AngularVelocity)` |
| **Q2** the arbitrary seed | **A** — strictly internal numerical normalization | `RotationSolver::seed_omega_bar_` (private, no public setter); reachable only through `RotationSolverTestSeam` |
| **Q3** result storage | **A** — one canonical geometric representation, named physical accessors, no duplicated state | `PhysicalFirstOrderRotation::Omega_geom` + `OmegaRadPerSecond()` |
| **Q4** exposure | **A** — seed-free normalized response through `NStar` | `NStar::RotationResponse()` returning `HartleFirstOrderResponse` |
| **Binding clarification** | no implicit physical spin on construction | construction publishes only `I`, `ω̄/Ω`, `ω̄'/Ω`; `Ω`, `J`, `ω̄(r)` require an explicit `AngularVelocity` |

## 3. The three first-order concepts, now distinct in code

| | Concept | Where it lives | Seed-dependent? | Reachable? |
|---|---|---|---|---|
| **A** | raw numerical solution `ω̄_raw`, `dω̄_raw/dr`, `J_raw`, `Ω_raw` | `RotationSolver` privates (`stored_omega_bar_`, locals) | yes, by construction | **no** — no public API returns any of it |
| **B** | seed-free response: `I [km³]`, `ω̄/Ω` (dimensionless), `ω̄'/Ω [km⁻¹]` | `HartleFirstOrderResponse`, via `NStar::RotationResponse()` | **no** — analytic cancellation, measured to `4.2e-15` | yes, for every built star |
| **C** | physical solution at a requested `Ω`: `Ω_geom [km⁻¹]`, `J [km²]`, `I`, `ω̄(r) [km⁻¹]`, `ω̄'(r) [km⁻²]` | `PhysicalFirstOrderRotation`, via `NStar::RotationAt(AngularVelocity)` | **no** | only after an explicit `AngularVelocity` |

The materialization is `ω̄_phys(r) = [ω̄/Ω](r) · Ω_geom`, with `J` then taken from the exterior
matching applied to the **scaled** surface derivative, `J = R⁴ ω̄'_phys(R)/6`. Deriving `J` from
the published array rather than from `I·Ω_geom` keeps `J = I Ω` a genuine consistency check on
the arrays a consumer actually receives — which is why detectors D1 and D3 fire on it (§7).

`I` is carried through unchanged rather than recomputed as `J/Ω`, so zero spin is an ordinary
case and not a `0/0`.

## 4. The conversion, and its single owner

```
Omega_geom [km^-1]  =  Omega_phys [rad s^-1] / c ,     c = Zaki::Physics::LIGHT_C_KM_S
```

Forward conversion: `AngularVelocity::GeomKmInverse()`. Inverse rendering:
`CompactStar::AngularVelocityGeomToRadPerSecond(double)`. **These two functions are the only
places on the governed first-order path where an angular velocity meets the speed of light.**
The inverse deliberately returns a plain `double` rather than an `AngularVelocity`, so it cannot
be used to smuggle a geometric number back into the scientific input path (ADR-0006 P1).

`c` is taken from ZakiLib rather than re-derived locally, per the owner's recorded
constant-authority intent and the Phase-3C `k_B` precedent: a second local `c` would recreate
exactly the dual-authority defect that increment eliminated. The **test** oracle is independent
(§5, V6).

**Residual, recorded not repaired.** `RotationSolver::Solve_Mixed` still multiplies by
`LIGHT_C_KM_S` inline at four sites when it fills the legacy two-fluid sequence. That is
`MixedStar` first-order code, explicitly out of Phase-4A scope; ADR-0006 P12 governs it as a
contract and defers conformance to the `MixedStar` track. Its export labels were corrected (§6).

## 5. Validation — ADR-0006 §7, every tolerance predeclared before implementation

Two new CTests. Neither imports `Zaki::Physics::LIGHT_C_KM_S` for a conversion check; both use
the SI-exact literal `c = 299792.458 km/s`, so a wrong *conversion* is detectable and not masked
by agreement with the constant production uses.

| # | Requirement | Bound | Measured | Where |
|---|---|---|---|---|
| **V1** | internal seed varied over `[1e-3, 1e3]`; `Ω`, `J_phys`, `ω̄_phys(R)`, `I` unchanged | `≤ 1e-10` | `Ω` **0**, `J` `2.07e-15`, `ω̄(R)` `4.05e-16`, `I` `1.99e-15`; real stars `4.25e-15` | V1a–V1f |
| **V2** | requested `Ω` recovered from the surface of the physical solution | `≤ 1e-13` | `1.97e-16` analytic, `2.07e-16` CMF | V2a–V2c |
| **V3** | `J = I · Ω_phys / c` | `≤ 1e-13` | `1.88e-16` analytic, `3.07e-16` CMF | V3a, V3b |
| **V4** | `I` bit-identical; golden artifact byte-identical; 2B-4B harnesses unchanged and green | bitwise | `I == SeqPoint::I` **bitwise** on the analytic star and all four CMF stars; artifact hash unchanged; both 2B-4B tests untouched and passing | V4a, V4b, §8 |
| **V5** | `ω̄_phys(r)` linear in requested `Ω`, node by node | `≤ 1e-14` | **0** over 4001 nodes, for `2Ω` and for `−Ω/2`; derivative likewise **0** | V5a, V5b |
| **V6** | `Ω_geom · c = Ω_phys` against the independent oracle | `≤ 1e-15` | **0** analytic and CMF | V6a–V6c |
| **V7** | annotations and exported header tokens match the stored values | exact | all seven checks pass | V7a–V7f |
| **V8** | zero spin: `Ω = 0`, `J = 0`, `ω̄ ≡ 0`, `I` finite and unchanged | exact | **exactly** zero; `I` bit-identical to the scale-free value | V8a–V8c |
| **V9** | slow-rotation diagnostic `Ω/Ω_K` reported for 1.0–2.0 M☉ | diagnostic | §9 | V9a |

Three checks beyond the ADR's list, because they make the contract falsifiable rather than
merely satisfied:

- **C1c** — the response's own surface combination `[ω̄/Ω](R) + R[ω̄'/Ω](R)/3` is **exactly 1**.
  A normalized quantity that is not normalized fails here, and it is simultaneously a unit
  check: a dimensionless quantity must come out a pure number. D4 fires on it.
- **V7e1** — the legacy sequence export is **inert** on the ordinary `NStar` path: no rows, no
  file. That is the separation ADR-0006 P7 asks to be made explicit, asserted rather than
  claimed. (`VecSaver::Export1D` returns immediately on an empty vector.)
- **V7e3** — the corrected header tokens are parsed against a written row, so a token cannot
  drift onto the wrong column.

## 6. Unit corrections

| Site | Before | After |
|---|---|---|
| `HartleResult::Omega` | annotated `[s^-1]`, stored km⁻¹ — an outright mislabel published through `GetHartleResult()` | **field removed.** `HartleResult` is now the second-order candidate's result only; first-order values live in the two new seed-free types |
| `HartleResult::J`, `::I`, `::omega_bar`, `::domega_bar` | seed-normalized values behind a public accessor | **removed** for the same reason. No consumer existed: `RotochemicalCache` reads only `hartle.p0`, and it is not compiled |
| `ExportResults` header | `omega_bar_c (1/s)`, `M`, `R`, `J`, `Omega (1/s)` — a seed and a seed-normalized `Ω` advertised as physical; three columns with no unit at all | `omega_bar_c_seed (1/s)`, `M (M_sun)`, `R (km)`, `J_seednorm (km^2)`, `Omega_seednorm (1/s)` — every token names the value written, and the seed-normalized quantities are labelled as such |
| `init_omega_bar` | no documented unit; a bare `5e-3` literal in the solve | documented as an arbitrary internal normalization; sourced from the private `seed_omega_bar_` (default unchanged at `5e-3`) |

**No cgs `J` or `I` accessor was added.** Those need `G`, whose authority is still unadjudicated
(`Zaki::Physics::SUN_M_KM` vs `GSL_CONST_CGSM_SOLAR_MASS`, ~6.2e-5 apart). ADR-0006 §10 defers
them deliberately; `c` alone suffices for `Ω`.

## 7. Detectors — applied, measured, reverted byte-identically

Each mutation was verified to have **exactly one** matching site before being applied, and the
file's SHA-256 was confirmed identical to its pre-mutation value after revert. The contract test
is **green** at the reverted state.

| | Mutation | ADR-0006 §8 requires | Fired |
|---|---|---|---|
| **D1** | omit the rescaling factor entirely (publish the normalized shape as if physical) | items 2, 3 | **FIRES** — V2a, V3a, V5a, V5b, V7b, V7b2, V7c, V7d, V8a |
| **D2** | multiply by `c` where the conversion must divide | item 6, catastrophically | **FIRES** — V6a at `8.99e10` relative error (`c²`), plus V2b, V3a, V7a, V7b2 |
| **D3** | rescale `ω̄` but leave `ω̄'` (and hence `J`) unscaled | items 2, 3 | **FIRES** — V2a, V3a, V5b, V7b, V7b2, V7c, V8a |
| **D4** | let the internal seed leak into the normalized response | item 1 | **FIRES** — V1b, V1c at `1.0e6` spread (six decades of seed), plus C1c, V2a, V3a, V7b, V7b2, V7c |

An initial detector run produced contaminated D3/D4 results: the revert restored the original
file *mtime* (`shutil.copy2` preserves it), so `make` skipped the rebuild and the D2 mutation
persisted in the built objects. The harness was corrected to stamp the mtime on restore and the
whole suite was re-run; the table above is the clean run. Recorded because a stale-build
false positive is exactly the kind of thing a detector audit exists to catch.

**No O(Ω²) detector was defined.** Its equations are not governed yet (ADR-0006 §9).

## 8. `I`, the artifacts, and the O(Ω²) candidate

- **`I` bit-identical.** `RotationResponse().I == GetSequence().I` bitwise on the analytic star
  and on all four CMF stars. The ODE, its coefficients, the seed value, the extraction formulas
  and `nstar_ptr->MomI` are untouched; the only change inside `FindNMomInertia` is that the seed
  is read from a private member whose default is the same `5e-3`, plus the construction of the
  normalized response *after* the extraction.
- **Seven artifacts byte-identical** (§10), and `hartle_I_dscmf1_debug.tsv` **re-emitted** from
  the rebuilt binary compares byte-identical to the committed golden — a stronger check than the
  hash of an untouched file.
- **O(Ω²) untouched.** `ODE_Hartle2_N_Fast` (87 lines), `SolveHartle2_N` (145),
  `ODE_Hartle2_Mixed_Fast` (7) and `SolveHartle2_Mixed` (4) are **byte-identical** to the
  pre-implementation source. The `j²` factor, the `p0` equation, `dε/dp`, the centre conditions,
  the surface shooting, `ξ₀` and `δM` are all exactly as INV-08 records them. The candidate
  still consumes the private `stored_omega_bar_` / `stored_domega_bar_` arrays, which still
  carry the arbitrary seed — **its outputs remain quadratic in that seed and are not
  normalized.** No candidate output was baselined.

## 9. V9 — slow-rotation diagnostic (reported, never enforced)

`Ω_K = (M/R³)^{1/2}` in geometric units, so no `G` is needed.

| `M` [M☉] | `R` [km] | `Ω_K` [km⁻¹] | `Ω/Ω_K` at 100 rad/s | at 2π·100 Hz | at 2π·716 Hz |
|---|---|---|---|---|---|
| 1.0000 | 13.4263 | `2.470068e-2` | `1.350e-2` | `8.485e-2` | `6.075e-1` |
| 1.4000 | 13.5453 | `2.884154e-2` | `1.157e-2` | `7.267e-2` | `5.203e-1` |
| 1.6000 | 13.4683 | `3.109722e-2` | `1.073e-2` | `6.740e-2` | `4.826e-1` |
| 2.0001 | 12.7123 | `3.791605e-2` | `8.797e-3` | `5.528e-2` | `3.958e-1` |

ADR-0006 §10 defers governing a slow-rotation threshold to the O(Ω²) work, so nothing fails on
these values. They are recorded because the 716 Hz column sits at `0.4–0.6` of the Keplerian
scale, which is where a *second-order* treatment would start to matter — a fact Phase 4C will
need and this increment must not pretend to have settled.

## 10. Final validation

| Configuration | Registered | Result | Wall time | Before 4A |
|---|---|---|---|---|
| Full, authenticated `COMPACTSTAR_EOS_DATA_ROOT` | **21** | **21/21 PASS** | 209.95 s | 19/19, 200.25 s |
| Self-contained, no data root | **11** | **11/11 PASS** | 15.20 s | 10/10, 13.97 s |

Both counts are up by one: `hartle_normalization_contract` (self-contained, 0.42 s) and
`hartle_normalization_cmf` (external data, 3.30 s). The `+9.7 s` and `+1.2 s` are those two
tests; no pre-existing test slowed measurably, consistent with the per-star work being one
`O(n)` division loop added to a solve that already integrates ~2600 nodes.

**Artifacts.** All seven Phase-3 goldens are byte-identical, and
`hartle_I_dscmf1_debug.tsv` **re-emitted from the rebuilt binary** compares byte-for-byte with
the committed file (`ddf018579364e9b3bf24f1e3a3e2577e70fb09b8b84eb718237d9da683aa9d15`).
No baseline was regenerated or re-baselined.

**Performance.** One first-order ODE solve per built star, unchanged: asserted directly through
the solve counter (P1a, P1b, P1c). Twenty-five materializations at different angular velocities
perform **zero** additional integrations — the first-order problem is linear, so a physical
solution is a scaling, not a solve.

## 11. Scope — what this increment did not do

No O(Ω²) equation was altered, repaired, normalized, executed for a result, or baselined. No
rotochemical work. `MixedStar` scientific equations untouched (only its legacy export *labels*
were corrected). `TOVSolver`, `Geometry.hpp`, the thermal drivers, `RotochemicalCache` and the
BNV / Decay / DarkCore modules untouched. `NStar::rot_solver` remains private. No cgs `J`/`I`
accessor. No slow-rotation threshold imposed. No golden re-baselined. Phase 4B is not begun.

## 12. Status

> **`PHASE-4A FIRST-ORDER PHYSICAL NORMALIZATION COMPLETE`**

INV-07 moves to **GOVERNED (ADR-0006 ACCEPTED) — FIRST-ORDER PHYSICAL NORMALIZATION CONFORMED**.
This is **not** a claim that rotation is validated: O(Ω²) remains an unverified candidate
(INV-08), and the physical first-order response has not yet been checked against independent
physical references across analytic and real-star cases. That is Phase 4B.

**Exact next task.** *Phase 4B — independently validate the physically normalized first-order
response across analytic and real-star cases before governing or replacing the O(Ω²) candidate.*
