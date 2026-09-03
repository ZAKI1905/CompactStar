# CompactStar Scientific Invariants

> **STATUS: RATIFIED — 2026-08-31.** Normative at authority rank 3 (`GOVERNANCE.md` §1, §7).
> Ratified at commit `617bb0e`, subject to the status distinctions below: ratification
> accepts **the register**, not the correctness of the code or physics it records. See
> *What ratification of this register means*.
>
> **Evidence base:** `docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md`, audited at
> commit `d91c31b`. Entries touching `RotationSolver` or `MixedStar` are marked ⚠ where
> commit `3639d71` (reachable from `9f70f14`) may supersede the cited evidence.
> **Rotation entries re-audited at `df859b5` (2026-09-02, Phase 4A-0):** INV-05 (rotation half),
> INV-07 and INV-08 now cite current source lines; evidence in
> `docs/validation/PHASE4_ROTATION_ENTRY.md`.

## Status vocabulary

| Status | Meaning |
|---|---|
| **GOVERNED (ACCEPTED)** | Ratified by an ACCEPTED ADR. **Normative** — this is what the convention *is*, independent of whether all code conforms. Conformance is tracked separately. |
| **VERIFIED CURRENT BEHAVIOR** | The code demonstrably does this. Says nothing about correctness. |
| **INTENDED BUT UNVERIFIED** | Code or documentation asserts this; not confirmed by execution or test. |
| **PROPOSED** | Suggested convention. Not yet adopted. |
| **UNRESOLVED** | Sources conflict, or the convention is absent. **Fail-closed.** |
| **LEGACY / COMPATIBILITY** | True of a retained older path, not of the canonical one. |

**GOVERNED (ACCEPTED) is a statement about the contract, not about the code.** An entry may be
governed while some component still violates it; such violations are recorded explicitly as
implementation nonconformance.

**A VERIFIED CURRENT BEHAVIOR entry is a statement about the code, not an endorsement of the
physics.** Promotion to intended physics requires an ADR.

## What ratification of this register means

The owner ratified this register on 2026-08-31 (`GOVERNANCE.md` §7). That decision means the owner accepts **the register**: its status vocabulary; the
ADR-backed contracts it records as GOVERNED (ACCEPTED); the fail-closed designation of its
UNRESOLVED entries; and its descriptions of verified current behavior as **the best available
evidence** about what the code does today.

It does **not** mean any of the following, and no later document may read it that way:

| Ratification did **not** | Because |
|---|---|
| Validate `INTENDED BUT UNVERIFIED` entries | They stay unverified until executed, tested, or derived — INV-08, INV-09 |
| Resolve `UNRESOLVED` entries | They remain fail-closed — INV-07, INV-11, and sub-items of INV-06 and INV-16 |
| Accept the Hartle O(Ω²) or rotochemical candidate code | `GOVERNANCE.md` §5 — merging is not ratification, and neither is the passage of time |
| Permit implementation nonconformance | `RotochemicalCache` (INV-01) and `PhotonCooling` (INV-15) remain nonconformant and must be corrected, not excused |
| Certify numerical correctness | No entry here is a validation result. A designated implementation is the accepted *owner* of a quantity, never a claim that it computes it correctly |
| Endorse the physics of a `VERIFIED CURRENT BEHAVIOR` entry | That label describes the code. Placeholder emissivities and arbitrary normalizations are recorded, not blessed |

Ratifying the register is what makes it *citable as authority* at rank 3. It settles what the
conventions **are**, not whether the code meets them or whether the numbers are right.

---

## INV-01 — Species profile semantics — **GOVERNED (ACCEPTED)**

**Statement.**

> The `StarProfile::BaryonDensity` column stores total baryon number density `n_B` in fm⁻³.
> Per-species profile columns store dimensionless composition fractions `Y_i = n_i / n_B`.
> Physical species number densities are derived as `n_i = Y_i · n_B`.

This applies identically whether the species label is human-readable (`"n"`, `"p"`, `"e"`,
`"mu"`) or a CompOSE numeric particle code.

**Normative decision.** `docs/adr/ADR-0001-species-profile-semantics.md` — **ACCEPTED
2026-08-31 by project-owner adjudication.** The owner supplied the authoritative EOS schema
contract; this is not an inference from code.

**Ingestion contract.** The input EOS table already carries `n_B` in its baryon-density column
and `Y_i` in its extra composition columns. `TOVSolver` and `NStar` **preserve** these values —
**no normalization is applied or wanted** on import. Verified: `TOVSolver.cpp:552` copies extras
verbatim; `NStar.cpp:199` and `:236` copy `n_B` and `Y_i` directly into the profile.

**Storage contract: RESOLVED. Implementation conformance: INCOMPLETE.**

| Conformance | Component | Evidence |
|---|---|---|
| ✅ Conformant | Direct-Urca number-density reconstruction | `StarContext.cpp:544-546` — `nn = Yn*nB` etc., with the in-code comment *"Convert fractions to number densities in fm^-3."* |
| ✅ Conformant | Charge-fraction construction | `StarContext.cpp:691-696` — `Y_q = Σ_i q_i Y_i`, dimensionless |
| ✅ Conformant | Legacy per-species density reconstruction | `TOVSolver.cpp:1971` — species column × baryon density |
| ❌ **Nonconformant** | `RotochemicalCache` | `RotochemicalCache.cpp:147` passes the raw `Y_i` column into `ComputeEnclosedNumber` (`:25-44`) and `ComputeStructuralDerivative` (`:47-104`), which document and treat it as `n_i` in fm⁻³ (`RotochemicalCache.hpp:116,133`). No `× n_B` is applied. |

The nonconformant component is **not compiled** (absent from every CMake source list) and
`Build()` is never called, so it has never produced output. This is pre-existing unvalidated
candidate code from `675b4a9`, not a regression.

**Confidence.** High. The contract is owner-supplied; conformance was verified file-by-file
against `9f70f14`.

**Impact.** The storage contract no longer blocks Phase 5. The `RotochemicalCache` correction —
construct `n_i = Y_i n_B` before the `N_i`, `A_i`, `B_i` integrals — is recorded as a Phase-5
implementation task. **No currently compiled behavior changes**, and no previously generated
output is invalidated.

**Documentation debt.** `StarProfile.hpp:45,296` still describe species columns as densities and
must be corrected; the name `rho_i` is semantically misleading under this invariant. See
ADR-0001 Consequences.

---

## INV-02 — Unit boundaries: CGS inside TOV, geometric outside — **VERIFIED CURRENT BEHAVIOR**

**Statement.** The TOV integrator works in CGS with explicit `G` and `c`. Conversion to
geometric units happens once, at `NStar::Append`. All stored profile data is geometric:
r, m in km; p, ε in km⁻²; ν′ in km⁻¹; baryon density in fm⁻³.

**Evidence.** `TOVSolver.cpp:1365-1406` (CGS RHS with
`GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT`); conversion at `NStar.cpp:179-202`, `:687-696`;
column labels `"r(km)"`, `"m(km)"`, `"p(km^-2)"`, `"eps(km^-2)"`.

**Confidence.** High. **Documented?** Storage units yes (ARCHITECTURE.md); **the CGS/geometric
split is undocumented.**

**Note.** Local conversion constants are re-derived rather than drawn from
`Zaki::Physics::Constants`. **The k_B part of that debt is closed (Phase 3C):** production
formerly carried **k_B at two different precisions** — `8.617333262145e-11` in the thermal
drivers and `8.617333262e-11` in `CompOSE_Thermo` — and all four consumers now delegate to the
authoritative generic constant `Zaki::Physics::K_BOLTZ_EV * 1.0e-6`. Adopting it was an
owner-authorized **numerical-method / constant-authority change**, not behavior-preserving
work: the direct effect on `CvDensity_cgs` and `C_⋆` is `+1.68e-11` (matching the analytic
expectation), while `L_ν`, `L_γ`, `M` and `R` are bit-identical at fixed state. Evidence:
`docs/validation/PHASE3C_BOLTZMANN_AUTHORITY.md`.

**Still outstanding local conversion debt:** the exact-duplicate constants centralized in
Phase 3A (`CompactStar/Units.hpp`) covered only `km³→cm³` and `MeV fm⁻³→erg cm⁻³`; the
**solar-mass authority** (`Zaki::Physics::SUN_M_KM` vs `GSL_CONST_CGSM_SOLAR_MASS`, differing at
`6.2e-5`) remains **unadjudicated and deferred**, and other geometric/compact-object conversions
are still re-derived locally. See INV-13.

**This invariant's own status is unchanged:** INV-02 governs the CGS-inside-TOV /
geometric-outside unit boundary, which is VERIFIED CURRENT BEHAVIOR and is not affected by the
k_B authority change.

**Angular quantities — a second governed boundary (Phase 4A, ADR-0006).** Rotation adds a
*physical/geometric* boundary alongside the CGS/geometric one: the public first-order rotation
API takes and returns a physical `Ω` in rad s⁻¹ while everything stored stays geometric
(`Ω [km⁻¹]`, `J [km²]`, `I [km³]`, `ω̄ [km⁻¹]`). It is crossed in exactly **two** functions,
`AngularVelocity::GeomKmInverse()` (`Ω_phys/c`) and `AngularVelocityGeomToRadPerSecond()`
(`× c`), both in `CompactStar/AngularVelocity.hpp`, which uses the authoritative
`Zaki::Physics::LIGHT_C_KM_S` rather than a local literal — the Phase-3C constant-authority
precedent. The legacy `MixedStar` sequence path (`RotationSolver::Solve_Mixed`) retains its own
inline `× c` and is governed but not yet conformant (ADR-0006 P12).

---

## INV-03 — Metric convention — **VERIFIED CURRENT BEHAVIOR** (domain governed by ADR-0004)

**Statement.** `ds² = −e^{2ν} dt² + e^{2Λ} dr² + r² dΩ²`, so `g_tt = −e^{2ν}`, `g_rr = e^{2Λ}`,
with `Λ = −½ ln(1 − 2m/r)`. **The metric convention itself is unchanged.**

**Domain and failure semantics are governed by ADR-0004** (ACCEPTED 2026-09-01, rank-2
authority), which supersedes this entry's former wording *"degenerate argument clamped:
`denom ≤ 0 → 1e-15`"*. That wording recorded one of **six** mutually inconsistent behaviors the
Phase-3D audit measured for the same degenerate input, disagreeing by a factor of `3.16e7` at
`f = 0` (ADR-0004 §11). The accepted contract is:

| Case | Behavior |
|---|---|
| `r = 0`, `m = 0` | exact regular-centre limit: `f = 1`, `Λ = 0`, `e^Λ = 1`, `w_V = 0` |
| `r < 0`; `r = 0` with `m ≠ 0`; `f = 1 − 2m/r ≤ 0`; non-finite `r` or `m` | **fail closed** |
| `r > 0`, `f > 0` | evaluate, **no clamp, no epsilon** |

**Three distinct statuses, deliberately not merged:**

- **Governed contract** — the above, for all code under ADR-0004.
- **Canonical-path conformance** — `CompactStar/Geometry.hpp` implements it and is the sole
  definition used by **both ordinary visible-sector `NStar` construction paths**:
  `NStar::BuildFromTOV` (Λ production and the baryon integrand, Phase 3D) and, since
  **Phase 3E-I2**, `NStar::Append` (Λ production) and `NStar::FinalizeSurface` (baryon
  proper-volume factor). `GeometryCache::DeriveLambdaFromMR_` also delegates.
  **No canonical `1e-15` clamp remains on any of those paths.**
- **Legacy nonconformance, still deferred** — the scalar `NStar::BaryonNumIntegrand(double)`
  (separate INV-14 defect, zero callers), all six `MixedStar` sites, and the inline forms in the
  **project-specific extension modules** (`DarkCore_Analysis`, `BNV_*`, `Decay_Analysis` — owner
  clarification, Phase 3F: research-project code consuming the core, not core candidates) still
  carry their own inline forms and their own degenerate behavior. **These have not migrated.**
  `MixedStar` waits on focused coverage (ADR-0004 §0-Q2); the project modules are outside core
  closure scope. (`TOVSolver::RadiusLoop` was deleted in Phase 3E-I4.) `NStar::EvaluateNu`'s boundary-condition `x = 1e-15` is **not this measure**
  (ADR-0004 §4.4) and is deliberately untouched.

**Evidence.** Canonical: `CompactStar/Geometry.hpp`; `NStar.cpp` (BuildFromTOV);
`GeometryCache.cpp:DeriveLambdaFromMR_`. Contract tests: `proper_volume_contract`,
`geometry_cache_measure_contract`. Legacy, unmigrated: `StarProfile.hpp:246-247`;
`SurfaceGravity.cpp:38-42`; `NStar.cpp` Path-1 blocks; `MixedStar.cpp`.

**No validated output moved.** Measured `max 2m/r = 0.481` across the four authenticated stars,
so no node on any validated path reaches a fail-closed branch (ADR-0004 §11).

**Confidence.** High. **Documented?** Yes.

---

## INV-04 — Proper-volume measure — **GOVERNED (ADR-0004 ACCEPTED) — BOTH ORDINARY VISIBLE-SECTOR NStar PATHS CONFORMED; MixedStar / CANDIDATE / SCALAR-ACCESSOR MIGRATIONS DEFERRED**

**Statement.** The canonical measure is `w_V(r) = 4πr² e^{Λ}`, with redshifted variants
`w_V e^{ν}` and `w_V e^{2ν}`.

**Ownership, resolved by ADR-0004 §0-Q1 (Option B) into three distinct roles:**

| Role | Owner | Scope |
|---|---|---|
| **A. Mathematical** | `CompactStar/Geometry.hpp` (`CompactStar::Geometry`) | `f`, `Λ`, `e^Λ`, `w_V`, **and the domain/failure semantics**. Dependency-neutral: standard library only, no layer edge. |
| **B. Cached representation** | `Physics::Evolution::GeometryCache` | `ExpLambda`, `WV`, `WVExpNu`, `WVExp2Nu`; ADR-0003 provenance. Obtains the *formula* from A; keeps its own `DataColumn` composition. |
| **C. Consumer integrands** | `NStar`, thermal drivers, … | their own physics factor (`n_B`, `c_V`, `Q_ν`) and their own unit conversions (`1e54` is INV-14, **not** part of `dV`). |

`Core` was **not** made to depend on `Physics/Evolution` to obtain the measure — that was the
Phase-3-entry wording (Option A) and ADR-0004 §8.2 shows it is unattainable for `MixedStar` and
unsafe for `NStar` (it would construct a provenance-bearing cache inside an open `EditScope`,
against ADR-0003).

**Conformed — BOTH ordinary visible-sector `NStar` construction paths.**

- *Path 2* (`NStar::BuildFromTOV`, Phase 3D): Λ production bit-identical; baryon integrand
  `|ΔB|/B = 1.368e-16` on 1.0 M☉, bitwise on 1.4/1.6/2.0, against the `1.0e-15` predeclared in
  ADR-0004 §7.1 **before** implementation.
- *Path 1* (`NStar::Append` + `NStar::FinalizeSurface`, **Phase 3E-I2**): Λ production
  **bit-identical**; baryon proper-volume factor migrated under the same governed bound, with
  measured movement ≤ **1.640e-16**. Because the composition now mirrors `BuildFromTOV`, Path-1
  and Path-2 `B` are **bitwise identical at all 17 measured comparisons** — the former
  conformance gap is closed, not merely bounded.
- `GeometryCache::DeriveLambdaFromMR_` delegates; all cached arrays bit-identical.

**NOT conformed — deferred, and deliberately still recorded here:**

| Site | Blocking |
|---|---|
| ~~TOV Path 1 — `NStar::Append` Λ block, `NStar::FinalizeSurface` baryon integrand~~ | **CONFORMED in Phase 3E-I2** once Phase 3E-0 supplied the coverage ADR-0004 §13 was waiting for |
| `NStar::BaryonNumIntegrand(double)` scalar accessor | separate INV-14 defect (missing `1e54`), zero callers — **must not be repaired opportunistically** |
| `MixedStar.cpp` — six sites, two build paths | zero coverage, two-sector mass semantics; ADR-0004 §0-Q2, §15 |
| `DarkCore_Analysis`, `BNV_*`, `Decay_Analysis` | **project-specific extension modules** (owner clarification, Phase 3F) — not core candidates; ADR-0004 §16 contract only; not a core closure prerequisite |

**Evidence.** `CompactStar/Geometry.hpp`; `GeometryCache.cpp`; `NStar.cpp` (BuildFromTOV).
Tests: `proper_volume_contract`, `geometry_cache_measure_contract`, `baryon_number_cmf`.

**Redshifted variants.** The primitive owns `w_V` only. `GeometryCache` composes `w_V e^{ν}` and
`w_V e^{2ν}` from the **single** `WV` array; `geometry_cache_measure_contract` G3 pins that
bitwise so a second measure cannot appear (ADR-0004 §12).

**`ARCHITECTURE.md` claim resolved.** `WV()` is the canonical *cached* measure; the canonical
*mathematical* measure is the primitive, and `Core`'s validated path now uses it. The claim that
`WV()` is "the universal radial integration measure" remains inaccurate for the deferred legacy
sites above, and `CURRENT_ARCHITECTURE.md` says so explicitly.

**Remaining action.** Migrate the deferred sites as their blocking coverage lands. **INV-04 is
not fully resolved until then.**

---

## INV-05 — Radial center convention — **VERIFIED CURRENT BEHAVIOR**

**Statement (as of `d91c31b`).** Integration starts at `r_min = 1 cm = 1e-5 km`, never r = 0.
Initial conditions `y[0] = (4/3)πr³ε(p_c)`, `y[1] = p_c`. Central ε clamped to
`[10·ε_min, 0.999·ε_max]`. **No series expansion is used** in either the TOV or the rotation
solver.

**Evidence.** `TOVSolver.hpp:559,629`; `TOVSolver.cpp:1708-1712`, `:2545-2550` (duplicated).

**Rotation half — re-audited at `df859b5` (Phase 4A-0, supersedes the `9f70f14` ⚠).**
`RotationSolver` carries `kR_EPS_KM = 1e-6 km` (`RotationSolver.cpp:31`) and `SafeR0()`
(`:33-38`, used by `Solve_Mixed` `:375` and `FindMixedMomInertia` `:947`), with `r_safe` guards in
the four first-order right-hand sides (`:244`, `:260`, `:276`, `:292`). `FindNMomInertia` starts
at the first strictly-positive profile radius, falling back to `kR_EPS_KM` (`:684-690`) — on
production profiles `r₀ = 1e-5 km` — with `ω̄(r₀) = init_omega_bar`, `ω̄'(r₀) = 0` (`:704-706`).
**No series expansion**; the regular series `ω̄ = ω̄_c[1 + (8π/5)(ε_c+p_c)r² + …]` bounds the
truncation at `≲ 1e-12` relative, and `I` is immune to the centre condition by construction
(`HARTLE_MOMENT_INERTIA.md` §15.1). The O(Ω²) candidate starts at `r[0]` with literal zeros
(`:1170-1171`, `:1204-1205`) — assessed in INV-08 (item D). Evidence:
`docs/validation/PHASE4_ROTATION_ENTRY.md` §7.3.

**Confidence.** High for TOV and for rotation (at `df859b5`). **Documented?** No — implicit only.

---

## INV-06 — Stellar surface convention — **VERIFIED CURRENT BEHAVIOR**

**Statement.** Surface is where `p(r) ≤ max(1e-15·p_c, eos_tab.pre[0])`. Surface quantities are
the last grid point: `R = r[-1]`, `M = m[-1]`, `z_surf = exp(ν[-1])`. The ν integration uses
`ν(R) = ½ ln(1 − 2M/R)` as its outward boundary condition. Photon emitting area is
`A_∞ = 4πR² e^{2ν(R)}`.

**Evidence.** `TOVSolver.cpp:1204-1208`, `:1372-1385`, `:1873`, `:2657`; `NStar.cpp:297,301,336`,
`:823-842`; `PhotonCooling_Details.cpp:310`.

**Phase 4C-G note (2026-09-02) — O(Ω²) surface semantics, no convention change.** The last node
`R_*` is the EOS-table-floor surface (`p_* = 3.351885e25 dyn cm⁻²`, `n_B = 1e-7 fm⁻³` for
DS(CMF)-1; `TOV_REFERENCE.md` §5). For the Hartle monopole sector this was audited and
classified **`SURFACE ADEQUATE AS-IS`**: exterior matching identities are exact in vacuum and
hold at `R_*` to `≈ 1e-6`; the surface mass-shell and number-boundary terms are evaluated
explicitly from `ε_*`, `n_*`; `ξ₀(R_*)` is the displacement of the floor isobar and differs from
the `p = 0` surface displacement by `O(ΔR/R_*) ≈ 4–7e-3` relative
(`docs/validation/PHASE4C_HARTLE2_DERIVATION.md` §11; ADR-0007 P7). The doc-comment at
`TOVSolver.hpp:549-554` still says "`1e-5` times smaller than the central pressure" where the
code uses `1e-15` (`TOVSolver.cpp:1206`) — a stale comment, not a behaviour difference.

**⚠ UNRESOLVED sub-item — heat-blanket base.** Two competing thresholds coexist:
`TbDefinition.hpp:60` uses ρ_b = 1e10 g/cm³ (located by inward scan,
`TbDefinition.cpp:91-95`), while `StarBuilder.hpp` carries
`Options::blanket_energy_density_km2 = 7.4237e-9`. **No document says which is authoritative.**

---

## INV-07 — Hartle first-order normalization — **GOVERNED (ADR-0006 ACCEPTED) — FIRST-ORDER PHYSICAL NORMALIZATION CONFORMED AND PHYSICAL RESPONSE INDEPENDENTLY VERIFIED**

**Statement.** The frame-dragging equation for ω̄ is linear and homogeneous, so its solution is
determined only up to scale. The code fixes that scale with a hard-coded
`init_omega_bar = 5e-3`. Consequently `HartleResult::Omega` and `.J` carry an **arbitrary,
unphysical normalization**; only `I = J/Ω` is scale-invariant and physically meaningful.

**Evidence (current, `df859b5`).** Seed `init_omega_bar = 5e-3` at `RotationSolver.cpp:701`
(`FindNMomInertia`) and `:954` (`FindMixedMomInertia`); extraction `J = R⁴y[1]/6`,
`Ω = y[0] + Ry[1]/3`, `I = J/Ω` at `:737-739` (and `:992-994`, `:488-490`). Historical
(`d91c31b`) locations were `:1240` and `:1276-1288`.

**Re-audit, 2026-09-02 (Phase 4A-0 — `docs/validation/PHASE4_ROTATION_ENTRY.md` §8–§9).**
- The seed's physical meaning was measured: on the audited stars `Ω_raw·c ≈ 2.2–2.4×10³ s⁻¹`
  (a fast millisecond pulsar), purely by accident of the number `5e-3`.
- The normalization contract was **derived** from linear homogeneity and exterior matching:
  `A = Ω_target_geom/Ω_raw`, `Ω_target_geom = Ω_phys/c`; `ω̄_phys = Aω̄_raw`, `J_phys = AJ_raw`,
  `I` unchanged; checked numerically through the public API to `2e-16`.
- Every O(Ω²) candidate output was measured to scale as the seed **squared** (entry record §12),
  confirming the impact statement below.
- Unit audit: `HartleResult::Omega` annotated `[s^-1]` (`RotationSolver.hpp:105`) still stores
  km⁻¹ (`:744`); `ExportResults` advertises `omega_bar_c (1/s)`, `Omega (1/s)` (`:637-638`) for
  seed-normalized values populated only by `Solve_Mixed` (`:540-542`, `× LIGHT_C_KM_S`); `J`,
  `M`, `R` are exported without units; `init_omega_bar` itself has no documented unit
  (`hpp:212`); `p0_c` has none (`hpp:117`).
- The legacy public entry points `Solve(Axis,…)`, `Solve(double,…)`, `ODE`, `GetMass`,
  `GetPress` are **declared but undefined** (`hpp:370-375,403,315,318`; no symbol in the library).

**Impact.** Any O(Ω²) quantity is quadratic in ω̄ and therefore inherits the square of this
arbitrary factor (measured). **This blocks Phase-4 and, transitively, Phase-5.**

**Normative decision.** **`docs/adr/ADR-0006-hartle-first-order-physical-normalization.md` —
ACCEPTED 2026-09-02 by project-owner adjudication** (Q1 = A + D, Q2 = A, Q3 = A, Q4 = A). The
governed contract: physical `Ω` in rad s⁻¹ at the public API carried by an explicit typed
quantity; geometric km⁻¹ internally with **one** named conversion `Ω_geom = Ω_phys/c`; rescaling
`A = Ω_geom/Ω_raw`; the seed **strictly internal**; one canonical geometric result with named
physical accessors and no duplicated state; unit-true exports; and a seed-free normalized
response (`I`, `ω̄/Ω`, `ω̄'/Ω`) exposed through `NStar`.

**Binding clarification (ADR-0006 Decision).** An `NStar` constructed without an explicit
physical spin does **not** acquire an implicit physical `Ω`. Construction may compute the
seed-free `I`, `ω̄/Ω` and `ω̄'/Ω`; physical `Ω`, `J`, `ω̄(r)` and `dω̄/dr` may be materialized only
after an explicit physical angular velocity is supplied.

**Source conformance: ✅ COMPLETE (Phase 4A, 2026-09-02).** Evidence:
`docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md`.

| Aspect | Status |
|---|---|
| **Governing convention** | ✅ **RESOLVED** — ADR-0006 (Q1 = A + D, Q2 = A, Q3 = A, Q4 = A) |
| **Source conformance** | ✅ **COMPLETE** — typed physical spin input, one conversion owner, seed internalized, seed-free response exposed, unit annotations corrected |
| **Contract validation** | ✅ **COMPLETE** — ADR-0006 §7 items V1–V9, two CTests, four detectors fired and reverted byte-identically |
| **Independent physical validation of the normalized response** | ✅ **COMPLETE (Phase 4B, 2026-09-02)** — `PHYSICAL FIRST-ORDER HARTLE RESPONSE VERIFIED`; evidence `docs/validation/PHASE4B_FIRST_ORDER_PHYSICS.md` |

What the code now does. The public scientific spin input is a physical `Ω` in rad s⁻¹ carried by
the typed `CompactStar::AngularVelocity` (`CompactStar/AngularVelocity.hpp`), whose
`GeomKmInverse()` and the companion `AngularVelocityGeomToRadPerSecond()` are the **only** two
places on the governed first-order path where an angular velocity meets `c`. The arbitrary seed
is the private `RotationSolver::seed_omega_bar_` (`RotationSolver.hpp:326`, default unchanged at
`5e-3`, `:321`), with **no public setter**; it is reached through the declared-never-defined
`RotationSolverTestSeam` (`:225`, `:551`), so seed invariance is *proved* rather than asserted.
That seam is precisely a **`PRIVILEGED TEST BACKDOOR — NOT SUPPORTED SCIENTIFIC API`** (Phase-4B
wording correction): calling it reachable *only* by the harnesses is too strong in C++, since any
translation unit could define the befriended type. What ADR-0006 Q2 requires does hold — no
supported public seed setter, no supported public seed constructor argument, and no production
consumer of the seam.
Construction publishes only the **seed-free** `HartleFirstOrderResponse` — `I`, `ω̄/Ω`, `ω̄'/Ω`
(`:159`) — through `NStar::RotationResponse()` (`NStar.hpp:387`); a physical `Ω`, `J` and `ω̄(r)`
exist only via `NStar::RotationAt(AngularVelocity)` (`:409`). **No implicit physical spin is
conferred by construction**, per ADR-0006's binding clarification.

Measured: seed invariance `≤ 4.3e-15` over six decades (bound `1e-10`); requested `Ω` recovered
to `2.1e-16` (bound `1e-13`); `J = I Ω_phys/c` to `3.1e-16` (bound `1e-13`); conversion exact
against an independent SI literal `c` (bound `1e-15`); `ω̄_phys` linear in `Ω` to `0` at every
node. **`I` is bit-identical and all seven Phase-3 artifacts are byte-identical.**

**Phase 4B — the response is physically right, not merely contract-conformant (2026-09-02).**
The seed-free shape `s(r) = ω̄(r)/Ω` and `s'(r) = ω̄'(r)/Ω` was compared **node by node** against
an independently derived and independently normalized profile (the conservative-form solver in
`tests/rotation/hartle_reference.hpp`, which normalizes by its own surface extraction and never
calls production's ODE, its coefficient helpers, or the materialization). Agreement is `2.9e-9`
(`s`) and `9.5e-9` (`s'`) on the exact analytic star against a predeclared bound of `1e-7`, and
`≤ 1.85e-5` and `≤ 2.28e-5` across the four authenticated CMF stars against `1e-4`. The response
further satisfies the exterior-matching identities `s(R) = 1 − 2I/R³` and `s'(R) = 6I/R⁴` against
the **independent** `I`, reproduces `I` through a volume integral that reads only the interior
(`1.1e-7` analytic, `≤ 3.1e-5` CMF), and reproduces two **derived** weak-field coefficients,
`ω(0)/Ω → 2(M/R)` and `ω(R)/Ω → 0.8(M/R)`, to `3.0e-4` and `1.7e-3` at `M/R = 0.002`. A detector
that corrupts the interior shape while leaving `I` and the surface untouched fails the Phase-4B
tests and passes every Phase-2B and Phase-4A test — the proof that this validation is not
redundant. **Zero production source changed.**

**The `[s^-1]` mislabel is gone**, together with the four other seed-normalized first-order
fields: `HartleResult` no longer carries `Omega`, `J`, `I`, `omega_bar` or `domega_bar` at all
(`RotationSolver.hpp:203`) and is now the second-order candidate's result only. The legacy
`ExportResults` header now labels its seed-normalized values as such
(`RotationSolver.cpp:646`).

**Not resolved by this.** The `MixedStar` two-fluid rotation path still carries its own inline
`× c` conversions and its own seed literal; ADR-0006 P12 governs it as a contract and defers
conformance to the `MixedStar` track. And **nothing at O(Ω²) is normalized** — see INV-08.

---

## INV-08 — Hartle perturbative order — **GOVERNED (ADR-0007 + ADR-0008 ACCEPTED) — MEASURE-COMPLETE O(Ω²) MONOPOLE SOURCE CONFORMED; INDEPENDENT REVALIDATION REQUIRED**

**Statement.** First order O(Ω) supplies frame dragging and I. Second order O(Ω²) supplies the
monopole structural perturbations `(m₀, p₀)` with isobar displacement `ξ₀ = −p₀/(dp/dr)`,
surface condition `p₀(R) = 0` enforced by exact superposition of a particular and a homogeneous
solution.

**Evidence (current, `df859b5`).** `ODE_Hartle2_N_Fast` `RotationSolver.cpp:1024-1110`;
`SolveHartle2_N` `:1126-1270`; superposition `:1233-1250`; `ξ₀` `:1252-1264`; `delta_M` `:1268`.
The candidate code is whitespace-identical to `675b4a9` (author "Claude", 2026-04-05); the merge
and the Phase-1B repair changed no formula. Full equation map:
`docs/validation/PHASE4_ROTATION_ENTRY.md` §10–§12.

**Current behaviour re-authenticated 2026-09-02 (Phase 4A-0). The historical `d91c31b`
findings below are retained, corrected where marked.**

- **Reachability — corrected.** `NStar::rot_solver` is private with no accessor (`NStar.hpp:100`),
  but `RotationSolver` is a public class: `AttachNStar`, `FindNMomInertia`, `SolveHartle2_N` and
  `GetHartleResult` are public and defined (`hpp:306,382,390,397`), and `NStar` grants
  `friend class RotationSolver` (`NStar.hpp:82`). An external solver attached to any `NStar`
  executes the candidate from user code; this was compiled, linked and run against the
  unmodified library (entry record §13). **Correct classification: `PUBLIC SCIENTIFIC CANDIDATE
  — ZERO REPOSITORY CALLERS — EXECUTABLE FROM USER CODE — UNVALIDATED`**, not "structurally
  unreachable". Zero repository callers remains true.
- **Variable.** The code's `p0` is, by annotation (`hpp:114`, `[km⁻²]`), by the one correct
  homogeneous term `dm0/dr = 4πr²(dε/dp)p0` (`:1054`) and by `ξ₀ = −p0/(dp/dr)` (`:1258-1261`),
  the **Eulerian monopole pressure perturbation `δp₀ = (ε+p)p₀*`** — not Hartle's dimensionless
  `p₀*`. Its evolution equation (`:1068`) is **dimensionally inconsistent under either reading**
  (the `p0` and `m0` coefficients are one power of `1/r` too many — the comment at `:1057` writes
  `r²(r−2m)` for the TOV factor `r(r−2m)`), and the GR factors `(1+dε/dp)ν'` and
  `(1+8πr²p)/(1−2m/r)` are absent.
- **`j²` factor — re-confirmed.** `S_m` (`:1088`) uses `1/(1−2m/r)`, the reciprocal of the
  `(1−2m/r)` part of `j² = e^{−2ν}(1−2m/r)`, and omits `e^{−2ν}`; the second `S_m` term is also
  wrong in sign, by a factor 2, and by one power of `r` (dimension km⁻¹ where dimensionless is
  required); `S_p` (`:1100`) has a dimensionally inconsistent second term. The comment
  `:1078-1080` claiming `ν` is unavailable is false — the profile carries `MetricNu`
  (`StarProfile.hpp:246`).
- **Boundary condition — corrected.** The superposition arithmetic is correct *for the condition
  it imposes*, `p0(R) = δp(R) = 0` (achieved to `1e-21`). That is **not** Hartle's condition: it
  forces `ξ₀(R) = 0` (no monopole surface displacement) and `p0_c = δp(r₀) ≠ 0` (measured
  `−0.47 p_c`, `−0.52 p_c`), i.e. the central density changes. Hartle's fixed-central-density
  solution has `p₀*(0) = 0` with no homogeneous admixture, `p₀*(R) ≠ 0`, `ξ₀(R) = p₀*(R)/ν'(R)`.
  The candidate confuses the Eulerian `δp` with the Lagrangian `Δp = δp + ξ₀ dp/dr`, whose
  vanishing already *defines* `ξ₀`. **Consequently it does not produce the fixed-`ε_c`
  derivatives INV-09 requires.**
- **`delta_M = m0[-1]` only** (`:1268`); the exterior term `J²/R³` is absent. Measured values at
  the raw seed are **negative** (`δM/M = −1.6` on the 1.4 M☉ CMF star, `−3.4e-6` on the
  constant-density star), where the physical monopole mass change is positive and small.
- **`dε/dp`** by centered profile finite differences (`:1140-1161`) — mathematically the EOS
  derivative by the chain rule on a barotrope, numerically inferior (range `3.8–4.9e3` on the
  CMF star); the correct authority is the EOS (`1/c_s²`). The `1.0` fallback (`:1148`) is the
  **causal limit `c_s² = 1`** — *correction:* incompressible matter has `dε/dp → 0`, not `→ ∞`
  (the constant-density star gives exactly `0`).
- **Centre — downgraded to a documented approximation.** Literal `{0,0}` / `{0,1}` starts at
  `r₀ = 1e-5 km` (`:1170-1171`, `:1204-1205`) versus the regular series
  `p₀* = (1/3)j_c²ω̄_c²r² + …`, `m₀ = (4π/15)(ε_c+p_c)[(dε/dp)_c+2]j_c²ω̄_c²r⁵ + …`: relative
  truncation `O((r₀/R)²) ≈ 6e-13` — acceptable for the correct equations; no series expansion.
- **Seed dependence — measured.** `SolveHartle2_N` consumes `stored_omega_bar_` /
  `stored_domega_bar_` (`:1184-1185`) written at seed `5e-3`; every output (`m0`, `p0`, `ξ₀`,
  `p0_c`, `δM`) scales as the seed **squared** to `≤ 1.5e-3` over six decades, exactly where the
  driver's fixed absolute tolerance `1e-10` is not engaged (entry record §12). Quadratic scaling
  is a normalization diagnostic, **not** correctness.
- Shipped source still carries `[FIX: confirm exact from textbook]` (`:1019`) and `???`
  (`:1091`). MixedStar second order is a no-op stub (`:1114-1120`, `:1274-1277`).
- *(Historical, `d91c31b` numbering: `:1554-1700`, `:1591-1678`, `:1502`, `:1504-1507`, `:1514`,
  `:1696`, `:1444`, `:1518`, `:1542-1548`, `:1702-1705`.)*

**Phase 4A note (2026-09-02) — the candidate was deliberately NOT touched.** ADR-0006
conformance changed the *first-order* surface only. All four second-order functions are
**byte-identical** to the pre-change source: `ODE_Hartle2_N_Fast` (87 lines,
`RotationSolver.cpp:1140`), `SolveHartle2_N` (145 lines, `:1242`), `ODE_Hartle2_Mixed_Fast`
(7 lines) and `SolveHartle2_Mixed` (4 lines). The `j²` factor, the `p0` equation, the `dε/dp`
reconstruction, the literal centre conditions, the surface shooting, `ξ₀` and `δM` are exactly
as recorded above. The candidate still reads the private `stored_omega_bar_` /
`stored_domega_bar_` arrays (`:1300`), which still carry the arbitrary seed, so **its outputs
remain quadratic in that seed and are not normalized** — ADR-0006 P9 requires a future
second-order product to be seed-free, and that is Phase-4C work. Two changes are visible from
here and neither is a repair: `HartleResult` lost its five seed-normalized *first-order* fields
(`RotationSolver.hpp:203`), so the struct is now the candidate's result only; and line numbers
moved. **No candidate output was executed for a result or baselined.** Evidence:
`docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md` §8.

**Phase 4C-G note (2026-09-02) — governed replacement derived and proposed; candidate untouched.**
Primary-source access was authenticated (the journal scan of Hartle 1967, ApJ 150, 1005, read
page by page; Hartle & Thorne 1968 §II as secondary) and the `l = 0` system CompactStar should
implement was derived from it: variable `p₀* = −ξ₀(dp/dr)/(ε+p) = ν'ξ₀` (H67 87–88, 99),
equations H67 (97) and (100) written per `Ω_geom²` from the verified `s`, `s'`
(`docs/validation/PHASE4C_HARTLE2_DERIVATION.md` §6–§7, dimensionally audited term by term),
fixed-`ε_c` centre condition `m₀(0) = p₀*(0) = 0` with **no** surface condition (H67 108, p. 1009;
HT68 §II f), regular series `p₀* → (1/3)j_c²ω̄_c²r²`, `m₀ → (4π/15)(ε_c+p_c)[(dε/dp)_c+2]j_c²ω̄_c²r⁵`
(H67 108, re-derived), EOS-owned `dε/dp` (no API exists — an implementation prerequisite),
`δM = m₀(R_*) + 4πR_*²ε_*ξ₀(R_*) + J²/R_*³` (H67 105–107 plus the explicit surface-shell term),
`ξ₀(R_*)` generally nonzero. Every 4A-0 defect hypothesis above was classified against the
primary text: all **CONFIRMED** (the centre-start item **PARTLY CORRECT** — bound confirmed,
series start now recommended); the candidate's own premises — its equation citations
"(30)–(33)" and "no `e^{−2ν}` available" — **REFUTED** (record §17). The surface convention was
audited for O(Ω²): **`SURFACE ADEQUATE AS-IS`** — exterior identities hold at the EOS-floor node
to `≈ 1e-6`, the shell and number-boundary terms are `≲ 1e-5`, and `ξ₀(R_*)` is the
floor-isobar displacement with an `O(ΔR/R) ≈ 4–7e-3` relative systematic (record §11). The
angular integration of H67 (109) shows `l = 0` suffices for every scalar count (record §14).
**`GOVERNANCE.md` §3.1 is applicable once ADR-0007 is accepted; all seven conditions are
recorded in ADR-0007 §9.** The candidate is byte-identical; no candidate output was executed.

**Normative decision (2026-09-02).**
**`docs/adr/ADR-0007-hartle-second-order-monopole-response.md` — ACCEPTED 2026-09-02 by
project-owner adjudication.** The governed contract: integrate Hartle's `p₀*`; the fixed-`ε_c`
family (`m₀(0) = p₀*(0) = 0`) with **no surface condition**; the `l = 0` equations H67 (97),
(100) per `Ω_geom²` from the verified `s`, `s'`; a regular-series start at `r₀`; `dε/dp` owned by
the EOS/TOV layer as the derivative of the star's own `ε(p)` interpolant; `δM = m₀(R_*) +
4πR_*²ε_*ξ₀(R_*) + J²/R_*³`; surface `SURFACE ADEQUATE AS-IS` with explicit `R_*` semantics
(`R_*` is never labelled the exact `p = 0` surface); seed-free coefficients materialized only at
an explicit `AngularVelocity`. `GOVERNANCE.md` §3.1 is **AUTHORIZED** by that ADR (§9) —
**the correction is not yet executed.** Modified at acceptance: the homogeneous
sequence-derivative response is **not** a public API in Phase 4C (validation use only).

**Phase 4C-I0 note (2026-09-02) — the EOS derivative authority of ADR-0007 P5 is in place; the
candidate is untouched.** `TOVSolver::GetEDensDeriv` now owns `dε/dp` as the derivative of the
same `ε(p)` Steffen interpolant that builds the star, delivered dimensionless (the cgs→geometric
factor `INV_FM4_2_Dyn_CM2/INV_FM4_2_G_CM3` matches `c²` to `1.5e-16` against an independent
literal) and **fail-closed** off-domain, and `StarProfile::HasEosDEdP()`/`GetEosDEdP()` carry it
to consumers all-or-nothing. `RotationSolver.{hpp,cpp}` were **not touched**; the four
second-order functions and `HartleResult` are byte-identical to `master`. The retired profile
finite difference was measured against the new authority on DS(CMF)-1: median agreement `~4e-7`
in the core but **155–490 % error at ~25 crust nodes per star**, which is why P5 removes it from
authority. Evidence: `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md`. **No O(Ω²) physics was
implemented, executed or baselined.**

**Phase 4C-I1 note (2026-09-03) — the candidate is GONE and the governed source is in place.**
`GOVERNANCE.md` §3.1 was **executed**: `SolveHartle2_N`, `ODE_Hartle2_N_Fast`, the two MixedStar
stubs, `HartleResult`, `GetHartleResult` and every proven candidate-only member were **deleted**,
and the governed fixed-`ε_c`, `Ω²`-normalized response replaced them in the same commit
(zero live definitions, zero public declarations, zero compiled callers remain). The new
`HartleMonopoleResponse` integrates `p₀*` from the regular-centre series with no shooting and no
surface condition, interpolates every background input at the **actual** ODE radius through one
shared bracket, takes `dε/dp` only from the Phase-4C-I0 authority (absence fails closed), and
composes `δM̂ = m̂₀(R_*) + 4πR_*²ε_*ξ̂₀(R_*) + I²/R_*³`. Measured: seed invariance **7.85e-15**
(bound 1e-10), exact quadratic materialization, `I` bit-identical, seven artifacts unchanged,
zero monopole solves on ordinary construction. Evidence:
`docs/validation/PHASE4C_I1_MONOPOLE_IMPLEMENTATION.md`.

**Phase 4D evidence (2026-09-03, `docs/validation/PHASE4D_MONOPOLE_VALIDATION.md`).** The governed response was recomputed by
two test-only routes that never call production's ODE — Hartle's `(m₀, h₀)` system (97)+(98)+(90) on
the tabulated background, and a continuum solver on the closed-form constant-density interior — and
compared node by node. **Implementation level:** analytic profile agreement `9.7e-9` (bound `1e-7`),
converging to the continuum solution as `h²`; centre series `≤ 5.0e-10` over the first ten nodes
(bound `1e-8`); shell identity exact, shell = 90 % of `δM̂` on the homogeneous star; Newtonian
intercepts `|δM̂/R³ − 1| = 1.1e-6`, `|3Mξ̂/R⁴ − 1| = 4.9e-7` (bound `5e-3`, monotone); continuum
first integral `6.1e-15` (bound `1e-9`; tabulated form is an exact-`h²` floor, `1.9e-9` at
N = 16001, `4.8e-10` at 32001); homogeneous `δM̂` vs the exact `dM/dp_c` `3.0e-9`; DS(CMF)-1 four
stars, second-order-isolated `≤ 3.8e-7` and fully independent `≤ 3.7e-5` (bound `1e-4`); near-vacuum
identity `≤ 1.7e-8`; **published:** Chandrasekhar & Miller (1974) Table I, 19/19 homogeneous
configurations, `ξ₀` and the shell-excluded `δM/M` to `≤ 7.3e-4`; Hartle & Thorne (1968) Tables 3/5
on the printed HW EOS, 8/8, `δR/R ≤ 4.9e-3`, `δM/M ≤ 1.1e-2` (bound `2e-2`); detectors M1–M9 all
fire and were reverted byte-identically. **Not met:** §7 item 11 on the tabulated DS(CMF)-1
background — the homogeneous `δM̂` vs `(dM/dp_c)δp_c` is `1.04e-3` (bound `1e-3`), resolution-
independent (10000/20000/40000), diagnosed as **density steps in the crust EOS that the nodal
`dε/dp` column cannot represent** (the column integrated over the crust misses 17 % of the crust's own
`Δε` at every resolution; the dominant feature is the table's crust–core transition, `Δε/ε = 36 %` over
`Δp/p = 1.5 %`, a layer ≈ 0.6 m thick that falls between profile nodes at every tested resolution). **The same Stieltjes evaluation against the profile's own density steps
reproduces the independent TOV-sequence derivative to `7e-5`, and applied to the SOURCED solution it
quantifies the contribution production omits at ≈ `4.6 %` of `δM̂`** — Hartle's `dE/dP` delta functions at
internal density discontinuities, which the smooth Steffen authority cannot carry (the 4C-I0 "FD crust
noise" was those steps). §7 item 10 has no independent derivative source (`c_s²`: CONDITIONAL CHECK
UNAVAILABLE; the retired FD moves `δM̂` by 5 %); §7 item 9 measured order `2.35` on `δM̂`, Richardson
residual `7.7e-4` at res 20000. The two 4D-entry re-scopings (Newtonian *intercept* per the ADR's own
wording; the tabulated first integral taken at the continuum level) are recorded verbatim with their
original measurements.

**Status.** **GOVERNED (ADR-0007 ACCEPTED) — O(Ω²) MONOPOLE SOURCE CONFORMED; IMPLEMENTATION INDEPENDENTLY VERIFIED; PHYSICAL VALIDATION FAILED ON TABULATED CRUSTS WITH DENSITY STEPS (`δM̂` ≈ −4.6 % ON DS(CMF)-1) — ADR-0007 AMENDMENT REQUIRED**. The implementation is
verified against independent formulations, analytic limits, continuum convergence and two published
second-order calculations; **on a tabulated crust with density discontinuities the accepted contract
omits Hartle's internal delta-function shells, and on DS(CMF)-1 that omission is ≈ `4.6 %` of `δM̂`**
(≈ `1e-3` of the homogeneous `δM̂`) — a substantive physical discrepancy that 4D was required to
report and not repair. Consequently **no monopole baseline was created** (`GOVERNANCE.md` §3.1
condition 7 remains deferred: **CORRECTION EXECUTED — IMPLEMENTATION INDEPENDENTLY VERIFIED —
PHYSICAL VALIDATION FAILED ON STEPPED CRUSTS — NO BASELINE YET**), **no O(Ω²) `δM̂` from a stepped
tabulated EOS may be cited as a result**, and the smallest next action is an **ADR-0007 amendment
adjudicating internal density discontinuities** (internal shells `4πr_i²Δε_i ξ₀(r_i)`, or term 1
integrated against `dε` on the EOS table), its implementation, a re-run of 4D Experiment J, then the
baseline. The `l = 2` sector remains out of scope, not verified.

**Phase 4D-RG (2026-09-03).** The amendment is drafted: `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md` (**PROPOSED**) —
the EOS source of Hartle's (97) is the measure `−4πr²ξ̂₀ dε` of his eq. (93); a nodal `dε/dp` column
cannot represent it (the 4C-I0 derivative authority itself is correct to `3e-14` at the nodes). In
scratch the measure form on the profile's own `ε` nodes converges DS(CMF)-1's `δM̂` to `4e-5` across
res 5000–40000 (today: 1.6 % erratic), sits at the sequence-derivative oracle's floor (`≈ 7e-5`), and
reduces exactly to today's result on smooth or constant-density stars; the internal jump operator is
certified to `2e-10…4e-9` on an exact two-layer star. Evidence: `docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md`.
**Phase 4D-RI (2026-09-03).** ADR-0008 was **ACCEPTED** (Q1–Q12) and its correction **implemented**:
the EOS energy-density contribution to `dm₀/dr` is now the measure `−4πr²ξ̂₀ dε`, integrated one
governed profile segment at a time with the profile nodes as mandatory integration boundaries, and
the surface shell is the terminal `ε_* → 0` atom of that same measure. `dε/dp` is retained for the
regular-centre series and diagnostics and is no longer the radial mass source. Measured: the
constant-density analytic star is **bitwise unchanged** (its interior measure is identically zero);
the smooth HT68 EOS agrees with the superseded differential form on `δM̂` to `≤ 1.3e-5`; the
same-partition source accounting reproduces production to `≤ 3.8e-7`; DS(CMF)-1 `δM̂` moves
`+6.5 / +5.4 / +4.8 / +3.2 %` at 1.0/1.4/1.6/2.0 M☉ — the direction and size ADR-0008 predicted —
and its radial spread falls from `1.6 %` to `3.7e-5`; the retired-FD substitution now moves `δM̂` by
exactly `0.0`. Detectors D1–D4 fire and revert byte-identically. Evidence: `docs/validation/PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md`.

**Status.** **GOVERNED (ADR-0007 + ADR-0008 ACCEPTED) — MEASURE-COMPLETE O(Ω²) MONOPOLE SOURCE CONFORMED; INDEPENDENT REVALIDATION REQUIRED**. The
implementation conforms to the amended contract; **the corrected independent revalidation — the
`(m₀,h₀)` and continuum oracles, the homogeneous sequence campaign, the published comparisons and
M1–M10 in full — has NOT been repeated**, so no O(Ω²) number here is validated physics, none may be
cited as a result, and **no monopole baseline exists or may be created** (ADR-0008 Q12).
`GOVERNANCE.md` §3.1: **CORRECTION IMPLEMENTED — INDEPENDENT REVALIDATION REQUIRED — NO BASELINE
YET.** ADR-0008 Validation D's monotonicity half is recorded as **not met** (the residual is the TOV
background's own resolution dependence). The `l = 2` sector remains out of scope, not verified.

---

## INV-09 — Fixed-ε_c versus equilibrium-sequence derivatives — **INTENDED BUT UNVERIFIED**

**Statement.** Hartle perturbation theory gives the rotational response at **fixed central
energy density**; physical spin-down follows an **equilibrium sequence** along which ε_c changes.
The two are reconciled by
`A_i = (∂N_i/∂Ω²)|_{ε_c}`, `B_i = (∂N_i/∂ε_c)|_Ω`, and the baryon-conserving reduction
`Z_i = A_i − B_i (A_B / B_B)`.

**Evidence.** `RotochemicalCache.hpp:59-76`; reduction implemented at
`RotochemicalCache.cpp:167-171`.

**Verified sub-claim.** The Z_i reduction is **correct in form** — the hardest conceptual piece
of the formalism is right.

**Defects.** `A_i` is computed as `∫ (dn_i/dp)·p₀·w_V dr` and **never divided by Ω²**
(`RotochemicalCache.cpp:93-103`) — so it is δN_i at one spin rate, not a derivative with respect
to Ω². `B_i` integrates perturbed-star radial grids against the **unperturbed** star's `WV()`
weight (`:135-136`, `:162-163`) — a metric/grid mismatch with a possible out-of-range read.

**Compilation status.** `RotochemicalCache.cpp` is **not in any CMake source list** in either
`d91c31b` or `9f70f14` (re-checked at `df859b5`). Never compiled; `Build()` has zero callers.

**Phase 4A-0 note (2026-09-02).** `A_i` requires the O(Ω²) solution of the **fixed-central-density**
family (`p₀*(0) = 0`) per unit `Ω²`. The current O(Ω²) candidate imposes `δp(R) = 0` instead,
which shifts the central density (INV-08), and its `p0` is seed-normalized; so even after the
missing `÷Ω²`, `A_i` computed from it would not be `(∂N_i/∂Ω²)|_{ε_c}`. The Phase-4 output for
Phase 5 must be `Ω²`-normalized response coefficients at fixed `ε_c`
(`docs/validation/PHASE4_ROTATION_ENTRY.md` §16; ADR-0006 P8–P9).

**Phase 4C-I1 note (2026-09-03).** The fixed-`ε_c` structural response is now **implemented**
(`NStar::ComputeHartleMonopoleResponse()` → `HartleMonopoleResponse`, coefficients per
`Ω_geom²`), so the `A_i` side of this invariant has a governed source to consume once Phase 4D
validates it. **This invariant is not resolved:** the baryon-conserving equilibrium-sequence
reduction — `B_i`, `A_B/B_B`, `Z_i` — remains Phase 5, the public homogeneous/sequence-derivative
response is deliberately not exposed, `dn_i/dp` is not implemented, and the defects recorded below
stand.

**Phase 4D note (2026-09-03).** The fixed-`ε_c` `A_i`-side source is now **independently verified at
the implementation level** (`docs/validation/PHASE4D_MONOPOLE_VALIDATION.md`), and the sequence-derivative identity
was exercised test-side: exact (`3e-9`) on the analytic star, `1.04e-3` on DS(CMF)-1 against the
predeclared `1e-3` because the tabulated crust's density steps are not represented by the nodal
`dε/dp` column — an omission worth ≈ `4.6 %` of the sourced `δM̂`. That is precisely the `B_i`-side
machinery this invariant will need (`δN_i` integrals over displaced layers carry the same internal
steps). *Phase 4D-RI (2026-09-03):* ADR-0008 is ACCEPTED and the `A_i`-side source is corrected, so
the fixed-`ε_c` fields Phase 5 will consume are measure-complete; the baryon-conserving reduction
(`B_i`, `A_B/B_B`, `Z_i`), `dn_i/dp` and its own measure completeness (ADR-0008 Q11) remain
unimplemented. **This invariant is not resolved.**
*Phase 4D-RG (2026-09-03):* `docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md` (PROPOSED) Q11 requires the scalar particle-number response to be
**measure-complete** — `δN̂_i` carries `−4πr²e^{λ}ξ̂₀ dn_i` with atoms at composition discontinuities; a nodal
`dn_i/dp` column is forbidden as the integrator's source. The `B_i` side inherits the same rule.

**Phase 4C-G / ADR-0007 note (2026-09-02).** ADR-0007 (ACCEPTED) makes the **fixed-`ε_c` Hartle
response the governed Phase-4 structural family**: Phase 4 delivers `(∂/∂Ω²)|_{ε_c}` coefficients
and nothing else. The baryon-conserving equilibrium-sequence reduction — `B_i`, `A_B/B_B`, and
`Z_i = A_i − B_i(A_B/B_B)` — **remains Phase 5**, and constant `B` is never imposed inside the
Hartle solve. At acceptance the owner further deferred *public* ownership of `B_i` /
`∂N_i/∂ε_c`: the regular homogeneous solution may be computed internally or test-side in 4D to
validate the sequence-derivative identity, but it is not a public Phase-4C API. **This invariant
is not resolved** — the defects below stand until Phase 5 addresses them.

**Phase 4C-G note (2026-09-02).** The complete O(Ω²) expression for a scalar count was derived
from H67 (109) (`docs/validation/PHASE4C_HARTLE2_DERIVATION.md` §14–§15): per unit `Ω_geom²`,
`A_i = ∫ w_V {(dn_i/dp)δp̂₀ + n_i[m̂₀/(r−2m) + (1/3)r²e^{−2ν}s²]} dr + w_V(R_*) n_i(R_*) ξ̂₀(R_*)`
with `n_i = Y_i n_B`. The `∫(dn_i/dp)p₀ w_V dr` integrand above is therefore **incomplete** as
well as un-normalized: the proper-volume (`m₀`) and time-dilation (`ω̄²`) terms and the boundary
term are of the same order; `h₀` cancels; the `l = 2` sector integrates to zero. The fixed-`ε_c`
family is the particular solution with `p₀*(0) = 0`; the regular homogeneous solution is the
TOV sequence derivative (`B_i`'s source) and ADR-0007 P11 proposes exposing it. `Z_i` stays
Phase 5.

---

## INV-10 — Thermal redshift convention — **VERIFIED CURRENT BEHAVIOR**

**Statement.** The thermal ODE variable is `x = ln(T∞ / T_ref)` with `T_ref = 1e8 K`, so
`T∞ = T_ref · exp(x)`. Local temperature is `T(r) = T∞ · e^{−ν(r)}`. Drivers contribute
`dx/dt = (1/T∞)(dT∞/dt)`.

**Evidence.** `ThermalState.hpp:52-66,166`; `ThermalState.cpp:76`; `StarContext.cpp:804-805`;
`TbDefinition.cpp:138`; `PhotonCooling_Details.cpp:324`; `NeutrinoCooling_Details.cpp:969`.

**Confidence.** High. **Documented?** Yes — `ThermalState.hpp:52-66`.

**Undocumented sub-item.** The heat-capacity cache is tabulated on a log-spaced grid over
`[1e-5, 1] MeV` with 160 points and is **clamped at both ends** (`StarContext.cpp:725-728`,
`:775-777`). Silent clamping outside that range is a numerical convention that should be governed.

---

## INV-11 — Chemical-imbalance redshift convention — **UNRESOLVED**

**Statement.** *No convention exists.* `ChemState` states that the meaning of η_i "is defined by
the drivers/microphysics layer, not by this class," and that units "are typically energy
(e.g. erg), but this is not enforced." No driver defines it.

**Evidence.** `ChemState.hpp:46-57`. `η_npe` and `η_npμ` appear **only in comments**
(`Rotochemical.cpp:104-111`, `ChemState.hpp:18`) — no symbol, enum, index constant, or accessor.

**What is undefined.** (a) Local versus redshifted frame — the F&R convention `η^∞ = e^ν η_local`
is unrepresented; (b) DOF count and ordering — `n_eta = 0` in `RunBuilder.cpp:39` and all four
main programs; (c) units — the driver writes `Z·2ΩΩ̇` where Z is a particle count, with no
susceptibility `∂μ/∂N` to convert number to energy; (d) whether components are per-species or
per-reaction-channel — `Rotochemical.cpp:110-111` documents the channel combination
`(Z_n − Z_p − Z_e)` in prose, then writes raw per-species `Z_i` into slot `i`.

**Confidence.** High that nothing is defined. **Fail-closed. Blocks Phase-5.**

---

## INV-12 — Cache invalidation — ✅ **RESOLVED for profile-derived caches** (ADR-0003)

**Resolution, 2026-09-01 (Phase 3B).** ADR-0003 (**ACCEPTED**) established a
`(StarProfile identity, Version())` provenance contract, and it is implemented and enforced by
CTest. All five measured hazards are closed:

| Hazard | Now |
|---|---|
| `GeometryCache` had no provenance | carries `(profile, version)`; `Matches()` answers staleness; still an immutable snapshot with no `Refresh()` |
| `C_⋆` key omitted the geometry | key is `(profile version, thermo identity, geometry provenance)`; a caller-supplied foreign geometry **fails closed** |
| `ProfileVersionedCache` keyed on version alone | keys on `(identity, version)` |
| `NeutrinoCooling` collided across equal-version stars | payload depends on profile **and** geometry provenance; mismatch fails closed via the driver's `ok=false` |
| `StarContext` never re-bound column views | contract **S1**: re-binds before use on a revision change, advances the cached revision only on success, and throws on an unusable schema |

Enforced by `cache_contract` (P1–P4), `cache_thermal_contract` (T3–T4) and `heat_capacity_v1`
(U7.d). The former `--audit-known-hazards` mode is **removed** — no known-bug output is retained
as expected behavior. Evidence, including the superseded pre-repair measurements:
`docs/validation/CACHE_CORRECTNESS.md`; contract: `docs/adr/ADR-0003-*.md`.

**Scope of this resolution — precise.** It covers **profile-derived** caches only. Deliberately
**outside** the ADR-0003 contract and **unchanged**: `RotationSolver`'s bracket/`omega_bar`
acceleration caches (algorithm-local, cannot outlive their own solve), `TOVSolver`'s EOS splines
and `gsl_interp_accel` (keyed to the imported EOS table, not to a profile), and the
`TimeSeriesObserver` name→pointer maps (bookkeeping). No provenance hazard was demonstrated for
any of them; this invariant makes no claim about future caches of other kinds.

---

### Historical record — the pre-repair measurement (Phase 2B-3, superseded)

**Statement (as it stood before Phase 3B).** Five structurally different invalidation rules
coexist.

| Rule | Mechanism | Users |
|---|---|---|
| 0 | **No gate at all** — deep copy at construction, never rebuilt | `GeometryCache` (`GeometryCache.cpp:174-237`) |
| 1 | `ProfileVersionedCache`: rebuild iff `version != last` | **NeutrinoCooling only** (`NeutrinoCooling_Cache.hpp:116`) |
| 2 | Shared snapshot drops *all* derived caches, then eagerly rebuilds two | `StarContext::RefreshDerivedCachesIfNeeded_` (`:241-281`) |
| 3 | Rule 2 **plus** own `prof_version` **plus** raw-pointer `thermo_tag` identity | Heat capacity (`:712-717`) |
| 4 | Bare null-check outside the version gate | `ChargeFractionYq` (`:636-642`) |

**Two structural hazards.**
- **`GeometryCache` carries no version stamp and exposes no `Invalidate()`.** If a profile is
  re-solved after construction, every downstream integral silently uses stale geometry while
  `StarContext`'s own caches correctly rebuild.
- **`StarContext` binds raw column pointers once, in its constructor** (`:177-196`), and
  `RefreshDerivedCachesIfNeeded_` never re-binds them. A profile mutation that reallocates
  column storage bumps the version — so payloads rebuild — while the seven cached pointers
  (`StarContext.hpp:219-225`) **dangle**. The version gate gives false confidence here.

**Correction to ARCHITECTURE.md:247**, which claims `ProfileVersionedCache` is "used by DUrca
mask, heat capacity, rotochemical coefficients." **False on all three counts.**

---

## INV-13 — Interpolation method — **VERIFIED CURRENT BEHAVIOR**

**Statement.** Interpolation is **linear**, not spline. `DataSet::Integrate` is therefore the
exact integral of a linear interpolant — i.e. the trapezoid rule. Heat-capacity interpolation is
linear in log T. CompOSE derivative nodes are linearly interpolated.

**Evidence.** `dependencies/include/Zaki/Vector/DataSet.hpp:572` —
`kInterpType = gsl_interp_linear;` with `gsl_interp_steffen` commented out on the next line.
Also `StarContext.cpp:736`; `CompOSE_Thermo.cpp:549`.

**Confidence.** High. **Documented?** **No — and documented misleadingly.** ARCHITECTURE.md
refers to "GSL spline" and "interpolation"; a reader would reasonably assume cubic.

**Why this matters.** This sets the accuracy ceiling of every radial integral in the code —
enclosed baryon number, moment of inertia, neutrino luminosity, heat capacity. The
linear-interpolation and trapezoidal components have **nominal O(Δr²) accuracy for sufficiently
smooth data**, so fourth-order behavior is not available from them.

**Nominal order is not observed order.** The statement above is a property of the quadrature
scheme in isolation. It does **not** license the stronger claim that a complete coupled
observable — a TOV solve feeding a Hartle solve feeding a cooling integration, over a
non-uniform grid, through clamped and table-interpolated EOS quantities — will exhibit
second-order convergence. Smoothness is assumed rather than established; the grid is set by the
integrator, not chosen for a refinement study; and clamping (INV-10) is not smooth at all.
**The convergence order of the complete calculation must be measured, not inferred from the
interpolation scheme.** Designing that measurement is Phase-2B work.

**Phase 4C-I0 note (2026-09-02) — the EOS tables are a separate interpolation domain.** This
invariant governs `Zaki::Vector::DataSet` (linear). The **EOS** splines inside `TOVSolver` are a
different mechanism and a different choice: `TOV_gsl_interp_type = gsl_interp_steffen`
(`TOVSolver.hpp:488`), a monotonicity-preserving C¹ cubic, used for `ε(p)`, `n_B(p)` and every
species column. ADR-0007 P5 makes the `dε/dp` authority the **analytic derivative of that same
interpolant** (`TOVSolver::GetEDensDeriv`), evaluated through the same accelerator — not a
re-interpolation under a smoother scheme, which would make the derivative inconsistent with the
`ε` the star was actually built from. Detector D3
(`docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md` §11) measured the consequence of violating that:
a `gsl_interp_cspline` derivative over the same table yields **negative** `dε/dp` at 22–66 crust
nodes on DS(CMF)-1, i.e. `c_s² < 0`. Steffen's monotonicity preservation is load-bearing here,
not cosmetic. *Phase 4D-RG note (2026-09-03):* a linearly interpolated background represents a
**measure** faithfully only through the *values* it carries at the nodes — the profile's `ε` nodes are the
EOS interpolant's own values (`5e-16`), so `Δε` between nodes is exact — while a *sampled derivative*
column loses whatever variation falls between samples (17 % of the crust's `Δε` on DS(CMF)-1). Inside a
sub-node feature the linear background also violates `dp/dr = −(ε+p)ν'` by `O(1)`, which is why
evaluating the EOS derivative at the actual ODE state does not converge either (`docs/adr/ADR-0008-measure-complete-eos-energy-density-source.md`). *(Unrepaired and unchanged: the comment block at `TOVSolver.hpp:494-498` still
describes a natural cubic spline, and `TOVSolver.hpp:549-554` still says `1e-5` where
`PressureCutoff` uses `1e-15` — stale comments, no behaviour difference.)*

---

## INV-14 — Baryon-number definition — **VERIFIED CURRENT BEHAVIOR**

**Statement.** `B = ∫ 4πr² n_B(r) (1 − 2m/r)^{−1/2} dr × 10⁵⁴`, evaluated as a linear-interpolant
integral over `[r₀, r_N]`. The `1e54` converts fm⁻³·km³ to a dimensionless count.

**Evidence.** `NStar.cpp:23` (`FM3_TO_KM3`), `:274-283`, `:314-316`. MixedStar equivalents at
`MixedStar.cpp:158-186`, both sectors using `mass_tot_dc` for the metric factor.

**⚠ Defect.** The public scalar accessor `NStar::BaryonNumIntegrand` (`NStar.cpp:1053-1068`,
declared `NStar.hpp:393`) computes the **same formula without the `1e54` factor** — same name,
different units (fm⁻³·km² vs km⁻¹). Harmless today only because it has zero callers.

**Relation to INV-01 (resolved).** This invariant integrates the `BaryonDensity` column `n_B`
directly, not a species column, so it is **unaffected in form** by ADR-0001. Any *per-species*
enclosed-number integral, however, must apply `n_i = Y_i n_B` — which the `RotochemicalCache`
integrand does not currently do (see INV-01 conformance table).

**Phase 4C-G note (2026-09-02).** This integral is Hartle & Thorne (1968) eq. (3e). Its O(Ω²)
monopole extension (density, proper-volume, time-dilation and boundary terms; `l = 2` integrates
out) is derived in `docs/validation/PHASE4C_HARTLE2_DERIVATION.md` §14 and reduces to this
expression at `Ω = 0`.

---

## INV-15 — Heat-capacity ownership — **GOVERNED (ACCEPTED)**

**Statement.**

> CompactStar has **one** physical heat capacity for the evolved isothermal-interior thermal
> degree of freedom: `C_⋆(T∞)`, the star-integrated, GR-weighted, temperature- and
> structure-dependent heat capacity
>
> ```
> C_⋆(T∞) = ∫₀^R c_V( T_local(r), n_B(r), Y_q(r) ) · 4π r² e^{Λ(r)} dr ,
>            T_local(r) = T∞ · e^{−ν(r)}
> ```
>
> with `c_V` supplied by the EOS/CompOSE thermodynamics, in the redshift convention of INV-10 and
> the proper-volume measure of INV-04. **Every** energy channel acting on the thermal degree of
> freedom divides by this same quantity:
>
> ```
> C_⋆(T∞) · dT∞/dt = − L_ν,∞ − L_γ,∞ + L_H,∞ + ⋯
> ```
>
> equivalently, in the evolved logarithmic variable,
>
> ```
> d/dt ln(T∞/T_ref) = ( − L_ν,∞ − L_γ,∞ + L_H,∞ + ⋯ ) / ( T∞ · C_⋆(T∞) )
> ```
>
> A driver-local constant such as `C_eff = 1e40 erg K⁻¹` is **not** the physical heat capacity of
> production thermal evolution and is not acceptable as a production physical denominator. It may
> be retained only as an **explicitly selected** test/debug approximation, and only if it applies
> to the whole thermal balance rather than to a single channel.

**Normative decision.** `docs/adr/ADR-0002-thermal-heat-capacity-ownership.md` — **ACCEPTED
2026-08-31 by project-owner adjudication.** The repository contained two mutually inconsistent
behaviors and no document ranked them; the physical convention is owner-supplied, not inferred.

**Designated implementation candidate.** `StarContext::HeatCapacityStar_Tinf(...)`
(`StarContext.hpp:151`; implemented `StarContext.cpp:704-818`) is the currently designated
implementation of `C_⋆(T∞)`. It integrates the CompOSE `c_V` (`StarContext.cpp:807-808`;
`CompOSE_Thermo.cpp:711-732`) at the Tolman local temperature `T∞ e^{−ν}`
(`StarContext.cpp:804-805`) against the canonical `GeometryCache::WV()` measure
(`StarContext.cpp:758`, `:810-811`), yielding erg K⁻¹ via `KM3_TO_CM3 = 1e15` (`:761`).

**Designation is not certification.** `C_⋆(T∞)` is the accepted physical *owner*; **this does not
assert that its present implementation is numerically correct.** It has never been tested — the
repository has no test suite. See *Numerical validation* below.

**Three statuses, tracked separately.**

| Aspect | Status |
|---|---|
| **Governing physics** | ✅ **RESOLVED** — ADR-0002, one canonical `C_⋆(T∞)` |
| **Source conformance** | ✅ **COMPLETE** — `PhotonCooling` divides by `C_⋆(T∞)`; the driver-local `C_eff` option is removed (roadmap Phase 2A-3) |
| **Numerical validation (V1)** | ✅ **COMPLETE** — `docs/validation/HEAT_CAPACITY_V1.md`, **V1 VERIFIED** |
| **Passive-cooling regression baseline** | ✅ **COMPLETE** — `GOVERNANCE.md` §3.1 condition 7 **SATISFIED**. Golden values at `tests/baselines/passive_cooling_cmf_1p6_debug.tsv`, nine epochs from 100 yr to 1 Myr, energy identity closing to 2.1e-16, state tolerance 1e-5. See `docs/validation/PASSIVE_COOLING_BASELINE.md`. **Regression only — the neutrino `Q0` normalizations remain placeholders and this validates no neutrino microphysics.** |

**Source conformance detail.**

| Conformance | Component | Evidence |
|---|---|---|
| ✅ Conformant | `NeutrinoCooling` | Obtains `ctx.star->HeatCapacityStar_Tinf(Tinf_MeV, *ctx.thermo, ctx.geo)` and divides `L_ν,∞` by it, converting per INV-10 |
| ✅ Conformant | `PhotonCooling` | Obtains the same canonical `C_⋆(T∞)` from `StarContext` and divides `L_γ,∞` by it (Phase 2A-3). `PhotonCooling::Options::C_eff` no longer exists |

The run now integrates the governed equation
`C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞`, both channels sharing one denominator — ADR-0002 **Pattern A**
(additive drivers, shared `C_⋆`). Pattern B, a centralized thermal-energy balance, remains the
deferred Phase-3 architectural question and was **not** introduced.

**Every passive-cooling output produced before this correction is superseded as a validation
reference** and must not be used as a regression baseline (ADR-0002). Those artifacts are retained
as historical products, not deleted.

**Numerical validation detail.** Required before the Phase-2B thermal baseline, and **not**
measurable against the existing passive-cooling curve, which encodes the rejected behavior:
dimensional check through `KM3_TO_CM3`; the degenerate `c_V ∝ T` low-temperature slope;
order-of-magnitude comparison against published total heat capacities (~10³⁷–10³⁸ erg K⁻¹ at
`T ≈ 10⁸ K` for a canonical star, which alone separates `C_⋆` from the `1e40` placeholder);
grid refinement in Δr — the quadrature is trapezoidal, of nominal order O(Δr²) for smooth data
(INV-13), with the observed order to be *measured* rather than assumed; insensitivity to the
`NT = 160` temperature tabulation (`StarContext.cpp:777`); explicit statement of the
endpoint-clamping behavior
(`StarContext.cpp:725-728`, see INV-10); and cache-rebuild correctness (INV-12). ADR-0002 §V1.

**Impact.** INV-15 no longer blocks a thermal validation baseline as a *decision*. The roadmap
circularity is repaired: photon-cooling conformance and targeted `C_⋆` verification are
**Phase 2A**, and the baseline is **Phase 2B**. No golden regression may be captured from the
existing passive-cooling outputs.

### Known implementation hazards — recorded separately, not part of INV-15

These are engineering/numerical defects on the heat-capacity path. They are **not** heat-capacity
semantics questions and are deliberately not folded into the invariant above.

- **Null-check ordering — ✅ CORRECTED (Phase 2A-3).** `NeutrinoCooling_Details.cpp` used to
  dereference `ctx.star` in the `HeatCapacityStar_Tinf` call twelve lines *before* its
  `if (!ctx.star)` guard, so the guard could never fire. The guard now precedes the first use.
  Engineering-class fix; no emissivity, rate, cache, option, or numerical constant changed.
- **Cache key omits the geometry.** `HeatCapacityStar_Tinf` accepts an optional `GeometryCache`
  and falls back to constructing one locally (`StarContext.cpp:754-755`), but the cache key is
  only `(profile version, thermo pointer identity)` (`:712-714`). A later call at the same profile
  version supplying a *different* `GeometryCache` reuses the earlier table. See INV-12 rule 3.
- **Silent endpoint clamping.** Outside the tabulated `[1e-5, 1] MeV` range (≈ `[1.16e5, 1.16e10] K`)
  the endpoint value is returned with no signal (`StarContext.cpp:725-728`). Already flagged as an
  ungoverned numerical convention in INV-10.

---

## INV-16 — Direct-Urca threshold — **VERIFIED CURRENT BEHAVIOR**

**Statement.** Direct Urca is permitted where the momentum triangle inequality
`k_Fn ≤ k_Fp + k_Fe` holds, with `k_F = (3π² n)^{1/3}`. Guards `n_B ≥ 1e-6 fm⁻³` and
`n ≥ 1e-12 fm⁻³` suppress a crust false positive. The boundary is the end of the **contiguous**
allowed region scanning outward from index 0.

**Species input (per INV-01 / ADR-0001).** The Fermi momenta are built from **derived number
densities**, not from the stored fraction columns:

> `n_n = Y_n n_B`,  `n_p = Y_p n_B`,  `n_e = Y_e n_B`

**The implementation already conforms** — `StarContext.cpp:544-546` performs exactly this
conversion, under the in-code comment *"Convert fractions to number densities in fm^-3."*

**Evidence.** `StarContext.cpp:412-591`, with 34 lines of in-code rationale at `:412-445`.

**Confidence.** High. **Documented?** Yes — among the best-documented physics in the codebase.

**⚠ Two items requiring adjudication.**
- **The muon channel is omitted.** Only the electron triangle is tested
  (`StarContext.cpp:564-569`). Whether `n → p + μ + ν̄_μ` is intentionally excluded is **not
  documented**. This is physically consequential for the DU boundary radius.
- The per-zone mask is built but never read: the only consumer uses just the boundary index and
  integrates `[0 .. durca_last]` (`NeutrinoCooling_Details.cpp:82,193`), never calling
  `DirectUrcaMask()`. Behavior is identical only if the allowed region is contiguous — which the
  construction enforces but the mask API implies need not hold.

**Note.** The threshold logic is real physics; the **rates gated by it are placeholders**
(`Q0_DU = 1.0e27`, `Q0_MU = 1.0e21`, linear in ρ, self-labeled
"Placeholder normalizations", `NeutrinoCooling_Details.cpp:100-102`). `K_PBF = 0.0` hardcoded.

---

## Summary

| Status | Entries |
|---|---|
| **GOVERNED (ACCEPTED)** | **INV-01** — ADR-0001, accepted 2026-08-31 · **INV-15** — ADR-0002, accepted 2026-08-31 · **INV-07** — ADR-0006, accepted 2026-09-02, **first-order source conformed and physical response independently verified 2026-09-02** · **INV-08** — ADR-0007, accepted 2026-09-02, **O(Ω²) monopole source conformed 2026-09-03; Phase 4D verified the implementation and FAILED the physics on tabulated crusts with density steps; **ADR-0008 ACCEPTED and the measure-complete source implemented 2026-09-03 (Phase 4D-RI)** — independent revalidation required, no baseline** |
| VERIFIED CURRENT BEHAVIOR | INV-02, 03, 04, 05, 06, 10, 12, 13, 14, 16 |
| INTENDED BUT UNVERIFIED | INV-09 |
| **UNRESOLVED (fail-closed)** | **INV-11** — and sub-items of INV-06, INV-16 |

**One unresolved invariant still blocks a downstream phase:**

- **INV-11** (η convention) blocks Phase 5.

**INV-07 is fully resolved for first order.** ADR-0006 (ACCEPTED 2026-09-02) settled the
contract, Phase 4A made the source conform, and **Phase 4B verified the normalized response
against independent evidence** — an independently normalized profile, the exterior-matching and
volume identities, and derived weak-field coefficients. What remains from it is `MixedStar`
conformance on its own track. **Rotation as a whole is still not validated:** O(Ω²) remains an
unverified candidate (INV-08), and nothing in the first-order work bears on the adequacy of the
slow-rotation truncation itself.

**INV-08 wording corrected (2026-09-02):** the O(Ω²) candidate is *publicly callable with zero
repository callers*, not "structurally unreachable".
**Phase 4C-G / ADR-0007 (2026-09-02) and 4C-I1 (2026-09-03):** the governed replacement was
derived from the primary source, **ADR-0007 was ACCEPTED**, and the `GOVERNANCE.md` §3.1
correction has now been **executed** — the AI-authored candidate is deleted and the governed
monopole response is in its place. INV-08 is GOVERNED and its **source conforms**, but
**no O(Ω²) number is validated physics**: that requires Phase-4D independent validation, and
**no monopole baseline exists or may be created** before it.
**Phase 4D (2026-09-03):** the independent validation was executed (`docs/validation/PHASE4D_MONOPOLE_VALIDATION.md`). Every
implementation-level line passed — independent `(m₀, h₀)` and continuum solvers, centre series,
Newtonian intercepts, shell, published Chandrasekhar–Miller and Hartle–Thorne second-order tables,
detectors M1–M9 — but ADR-0007 §7 item 11 was **not met** on DS(CMF)-1 (`1.04e-3` vs `1e-3`,
resolution-independent), diagnosed as crust density steps invisible to the nodal `dε/dp` column.
The Stieltjes evaluation against the profile's own density steps then showed the omission is
≈ `4.6 %` of the sourced `δM̂` on DS(CMF)-1 — **HARTLE MONOPOLE VALIDATION FAILED** on stepped crusts;
implementation verified, contract amendment required, still no baseline.

**INV-01 is resolved** as a storage contract (ADR-0001). One implementation nonconformance
remains — `RotochemicalCache` must construct `n_i = Y_i n_B` — tracked as a Phase-5 task, not as
an open invariant.

**INV-15 is resolved** as a physical ownership convention (ADR-0002): one canonical `C_⋆(T∞)`.
It **no longer blocks the thermal baseline as a decision.** What remains is work, tracked in two
places and not as an open invariant: `PhotonCooling` source conformance and targeted `C_⋆`
verification are roadmap **Phase 2A**; the baseline itself is **Phase 2B**. The `NeutrinoCooling`
null-check ordering defect is recorded as a separate implementation hazard, not as part of INV-15.

**GOVERNED (ACCEPTED) is a claim about the contract, not about the code.** Both governed
invariants currently have live implementation nonconformance — `RotochemicalCache` for INV-01,
`PhotonCooling` for INV-15 — and neither designated implementation has been numerically validated.

Accepting ADR-0001 resolves only INV-01; accepting ADR-0002 resolves only INV-15. **Neither
conferred accepted status on this document as a whole.** The register was ratified separately by
the owner on 2026-08-31 (`GOVERNANCE.md` §7) — a decision about the register, not about the
physics or the code it describes.
