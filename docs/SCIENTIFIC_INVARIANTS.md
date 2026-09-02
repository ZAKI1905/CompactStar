# CompactStar Scientific Invariants

> **STATUS: RATIFIED — 2026-08-31.** Normative at authority rank 3 (`GOVERNANCE.md` §1, §7).
> Ratified at commit `617bb0e`, subject to the status distinctions below: ratification
> accepts **the register**, not the correctness of the code or physics it records. See
> *What ratification of this register means*.
>
> **Evidence base:** `docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md`, audited at
> commit `d91c31b`. Entries touching `RotationSolver` or `MixedStar` are marked ⚠ where
> commit `3639d71` (reachable from `9f70f14`) may supersede the cited evidence.

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

## INV-05 — Radial center convention — **VERIFIED CURRENT BEHAVIOR** ⚠

**Statement (as of `d91c31b`).** Integration starts at `r_min = 1 cm = 1e-5 km`, never r = 0.
Initial conditions `y[0] = (4/3)πr³ε(p_c)`, `y[1] = p_c`. Central ε clamped to
`[10·ε_min, 0.999·ε_max]`. **No series expansion is used** in either the TOV or the rotation
solver.

**Evidence.** `TOVSolver.hpp:559,629`; `TOVSolver.cpp:1708-1712`, `:2545-2550` (duplicated).

**⚠ Supersession.** `9f70f14` adds explicit regularization in `RotationSolver`:
`kR_EPS_KM = 1e-6` km and `SafeR0()` (`RotationSolver.cpp:31-37`), with `r_safe` guards at
`:244,260,276,292`. **The rotation half of this entry must be re-audited against `9f70f14`.**

**Confidence.** High for TOV, superseded for rotation. **Documented?** No — implicit only.

---

## INV-06 — Stellar surface convention — **VERIFIED CURRENT BEHAVIOR**

**Statement.** Surface is where `p(r) ≤ max(1e-15·p_c, eos_tab.pre[0])`. Surface quantities are
the last grid point: `R = r[-1]`, `M = m[-1]`, `z_surf = exp(ν[-1])`. The ν integration uses
`ν(R) = ½ ln(1 − 2M/R)` as its outward boundary condition. Photon emitting area is
`A_∞ = 4πR² e^{2ν(R)}`.

**Evidence.** `TOVSolver.cpp:1204-1208`, `:1372-1385`, `:1873`, `:2657`; `NStar.cpp:297,301,336`,
`:823-842`; `PhotonCooling_Details.cpp:310`.

**⚠ UNRESOLVED sub-item — heat-blanket base.** Two competing thresholds coexist:
`TbDefinition.hpp:60` uses ρ_b = 1e10 g/cm³ (located by inward scan,
`TbDefinition.cpp:91-95`), while `StarBuilder.hpp` carries
`Options::blanket_energy_density_km2 = 7.4237e-9`. **No document says which is authoritative.**

---

## INV-07 — Hartle first-order normalization ⚠ **UNRESOLVED**

**Statement.** The frame-dragging equation for ω̄ is linear and homogeneous, so its solution is
determined only up to scale. The code fixes that scale with a hard-coded
`init_omega_bar = 5e-3`. Consequently `HartleResult::Omega` and `.J` carry an **arbitrary,
unphysical normalization**; only `I = J/Ω` is scale-invariant and physically meaningful.

**Evidence.** `RotationSolver.cpp:1240` (`d91c31b`). Extraction:
`J = R⁴y[1]/6`, `Ω = y[0] + Ry[1]/3`, `I = J/Ω` (`:1276-1288`).

**⚠ Persists.** `init_omega_bar` still appears 13 times in `9f70f14`. **Unresolved in both
states.**

**Secondary defect.** `HartleResult::Omega` is documented `[s^-1]` (`RotationSolver.hpp:105`)
but stores a geometric km⁻¹ value. Elsewhere the file multiplies by
`Zaki::Physics::LIGHT_C_KM_S` where s⁻¹ is wanted.

**Impact.** Any O(Ω²) quantity is quadratic in ω̄ and therefore inherits the square of this
arbitrary factor. **This blocks Phase-4 and, transitively, Phase-5.**

**Proposed resolution.** Rescale to a physical Ω before any second-order quantity is exported,
and correct the unit annotation. Requires an ADR.

---

## INV-08 — Hartle perturbative order — **INTENDED BUT UNVERIFIED** ⚠

**Statement.** First order O(Ω) supplies frame dragging and I. Second order O(Ω²) supplies the
monopole structural perturbations `(m₀, p₀)` with isobar displacement `ξ₀ = −p₀/(dp/dr)`,
surface condition `p₀(R) = 0` enforced by exact superposition of a particular and a homogeneous
solution.

**Evidence.** `RotationSolver.cpp:1554-1700`; superposition `:1591-1678`.

**Verified sub-claim.** The superposition construction is correctly implemented — the one
second-order component that matches its documentation.

**Unverified and defective sub-claims (as of `d91c31b`).**
- The factor `j² = e^{−2ν}(1 − 2m/r)` is **dropped** from the source terms. The code names it
  (`:1502`), states that ν is unavailable in the fast cache (`:1504-1507`), and omits
  `e^{−2ν}`. `S_m` (`:1514`) then uses `1/(1 − 2m/r)` — the **reciprocal** of the required factor.
- `delta_M = m0[-1]` only; the exterior matching term `J²/R³` is absent (`:1696`).
- `dε/dp` is reconstructed by centered finite differences of the **profile**, not obtained from
  the EOS, with a `1.0` fallback labeled "incompressible limit" (incompressible is dε/dp → ∞).
- Central boundary conditions are literal `{0,0}` and `{0,1}` — no series expansion.
- Shipped source carries `[FIX: confirm exact from textbook]` (`:1444`) and `???` (`:1518`).
- MixedStar second order is a no-op stub (`:1542-1548`, `:1702-1705`).
- Zero call sites; `rot_solver` is private with no accessor — **structurally unreachable**.

**Status.** **UNVERIFIED SCIENTIFIC CANDIDATE** under `GOVERNANCE.md` §5. Must not be cited as
implemented physics.

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
`d91c31b` or `9f70f14`. Never compiled; `Build()` has zero callers.

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
| **GOVERNED (ACCEPTED)** | **INV-01** — ADR-0001, accepted 2026-08-31 · **INV-15** — ADR-0002, accepted 2026-08-31 |
| VERIFIED CURRENT BEHAVIOR | INV-02, 03, 04, 05⚠, 06, 10, 12, 13, 14, 16 |
| INTENDED BUT UNVERIFIED | INV-08⚠, 09 |
| **UNRESOLVED (fail-closed)** | **INV-07⚠, 11** — and sub-items of INV-06, INV-16 |

**Two unresolved invariants block downstream phases:**

- **INV-07** (Hartle normalization) blocks Phase 4, transitively Phase 5.
- **INV-11** (η convention) blocks Phase 5.

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
