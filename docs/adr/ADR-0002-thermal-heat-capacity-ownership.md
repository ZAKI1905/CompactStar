# ADR-0002 — Heat-capacity ownership for the evolved thermal degree of freedom

| Field | Value |
|---|---|
| **Status** | **ACCEPTED** |
| **Date proposed** | 2026-08-31 |
| **Date accepted** | 2026-08-31 |
| **Authority** | **Project-owner adjudication.** The owner fixed the physical convention; this is not an inference from repository evidence, since the repository currently contains two mutually inconsistent behaviors and no document ranked them. |
| **Change class** | scientific-semantic |
| **Governing authority** | This ADR is now the normative authority for the heat capacity entering the thermal energy balance |
| **Affected invariants** | INV-15 (primary), INV-10 (redshift convention, unchanged), INV-12 (cache invalidation), INV-13 (interpolation order) |
| **Unblocks** | The heat-capacity prerequisite of the thermal validation baseline. Source conformance becomes roadmap **Phase 2A**; the baseline itself becomes **Phase 2B**. |
| **Explicitly does not decide** | Where in software the division by `C_⋆` is performed — see §6, *Deferred architectural question* |

## Context

CompactStar evolves a single scalar thermal degree of freedom. Per **INV-10**, the ODE variable is

```
x = ln(T∞ / T_ref),   T_ref = 1e8 K,   T(r) = T∞ · e^{−ν(r)}
```

and each driver contributes `dx/dt = (1/T∞)(dT∞/dt)` additively through `RHSAccumulator`
(`ThermalState.hpp:52-66`; `PhotonCooling_Details.cpp:324`; `NeutrinoCooling_Details.cpp:969`).

Two live drivers write into that one slot, and **each divides its luminosity by a different
heat capacity.** Because `RHSAccumulator` is additive and both contributions land in the same
element, the integrated equation is not a physical energy balance at all:

```
dT∞/dt  =  − L_ν,∞ / C_⋆(T∞)   −   L_γ,∞ / C_eff        with C_eff ≠ C_⋆
```

`C` is a property of the star — of its structure, composition, and temperature. It is not a
property of an emission channel. A configuration in which the same degree of freedom has two
different heat capacities depending on which term is being added is not a modeling choice with
a defensible physical reading; it is an inconsistency.

### Why this had to be adjudicated before, not after, a baseline

The roadmap as previously written was circular. Phase 2 declared the thermal validation baseline
**blocked by INV-15**, on the correct grounds that a passive-cooling regression captured while
two different `C` feed one energy equation would enshrine the inconsistency. Phase 3 then placed
the heat-capacity ownership correction **after** Phase-2 baselines existed. Each waited on the
other.

Resolving that by capturing the baseline anyway would have frozen a known-wrong energy equation
into a golden regression file, and every later correction would then have been measured as a
*regression against a placeholder*. `AGENTS.md` §5 requires a baseline before scientific
refactoring; it does not authorize a baseline that is known in advance to encode incorrect
physics. The circularity is repaired by ADR-0002 together with the roadmap change it mandates
(see Consequences).

## Evidence

Authenticated against `modernization/governance-foundation` at `ba49e10`, which carries the same
source tree as `9f70f14` for every file cited here.

### Path A — PhotonCooling: driver-local constant

| Claim | Evidence |
|---|---|
| `PhotonCooling::Options` carries a driver-local heat capacity as a plain configuration scalar | `PhotonCooling.hpp:229` — `double C_eff = 1.0e40;` inside `struct Options` |
| Documented default is 1.0e40 erg/K | `PhotonCooling.hpp:214-229` — *"Effective heat capacity … Units: erg/K … Here it is treated as a user-provided effective constant."* |
| The RHS divides by it | `PhotonCooling_Details.cpp:320` — `d.dTinf_dt_K_s = -d.L_gamma_inf_erg_s / drv.GetOptions().C_eff;` |
| …then converts to the log variable per INV-10 | `PhotonCooling_Details.cpp:324` — `d.dLnTinf_dt_1_s = d.dTinf_dt_K_s / d.Tinf_K;` |
| Only validation is positivity | `PhotonCooling_Details.cpp:229-234` — `if (!(drv.GetOptions().C_eff > 0.0)) { d.ok = false; … }` |
| The live end-to-end program sets it by hand and flags it as wrong | `spin_therm_evol_2_main.cpp:245` — `photonOpts.C_eff = 1.0e40; // Change this!` |
| A second program does the same | `spin_therm_evol_main.cpp:178` — `thermOpts.C_eff = 1.0e40;` with the comment *"big C_eff → slow cooling for demonstration"* |
| **No coupling to the stellar heat capacity exists in this driver** | `HeatCapacityStar_Tinf` appears nowhere in `PhotonCooling.hpp`, `PhotonCooling_Details.hpp`, `PhotonCooling.cpp`, or `PhotonCooling_Details.cpp`. The driver touches `ctx.star` only for the envelope `Tb` mapping and surface gravity (`PhotonCooling_Details.cpp:153-177`), never for `C`. |

The header itself states the gap plainly: *"It does not compute C_eff from EOS/microphysics
(yet)"* (`PhotonCooling.hpp:122`). `C_eff = 1e40 erg/K` is a demonstration knob chosen to make a
cooling curve visibly evolve, not a computed physical quantity.

### Path B — NeutrinoCooling: the GR-integrated stellar heat capacity

| Claim | Evidence |
|---|---|
| The driver obtains `C` from the star context, keyed on the evolved `T∞` | `NeutrinoCooling_Details.cpp:889` — `d.C_eff_erg_K = ctx.star->HeatCapacityStar_Tinf(Tinf_MeV, *ctx.thermo, ctx.geo);` |
| `T∞` is converted K → MeV first | `NeutrinoCooling_Details.cpp:880` with `MEV_PER_K = 8.617333262145e-11` (`:214`) |
| CompOSE thermodynamics is a hard requirement | `NeutrinoCooling_Details.cpp:871-877` — fails closed with `ctx.thermo == nullptr` |
| The result is validated for positivity and finiteness | `NeutrinoCooling_Details.cpp:891-896` |
| It becomes the RHS denominator | `NeutrinoCooling_Details.cpp:968` — `d.dTinf_dt_K_s = -d.L_nu_inf_erg_s / d.C_eff_erg_K;` |
| …then converts to the log variable per INV-10 | `NeutrinoCooling_Details.cpp:969` |

The earlier constant-`C` policy is present but commented out
(`NeutrinoCooling_Details.cpp:401-416`, `:786`, `:866`), including a comment that proposed
*"reuse the same C_eff as photon cooling in a shared thermal config object"* (`:403`). That is
the historical record of the very ambiguity this ADR closes — and it closes it the other way.

### The stellar heat-capacity calculation, documented rather than inferred

Traced through `StarContext::HeatCapacityStar_Tinf` (`StarContext.cpp:704-739`) and
`StarContext::BuildHeatCapacityCache_` (`StarContext.cpp:744-818`).

| Property | Present implementation | Evidence |
|---|---|---|
| **Governing integral** | `C_⋆(T∞) = ∫ c_V(T_local(r), n_B(r), Y_q(r)) dV` | `StarContext.hpp:135-143`; loop at `StarContext.cpp:788-817` |
| **Local temperature relation** | Tolman isothermal core: `T_local(r) = T∞ · e^{−ν(r)}`, via `GeometryCache::ExpMinusNu()` | `StarContext.cpp:804-805` |
| **Redshift convention** | Consistent with **INV-10** — the same `T∞ = T_local e^{ν}` convention the ODE variable uses | `StarContext.cpp:804-805` vs `ThermalState.hpp:52-66` |
| **Proper-volume weighting** | `dV = 4πr² e^{Λ(r)} dr`, taken from the canonical `GeometryCache::WV()` (INV-04) | `StarContext.cpp:758`, `:810-811`; `GeometryCache.hpp:132` |
| **Thermodynamic input** | `CompOSE_Thermo::CvDensity_cgs_ForCooling(T_MeV, n_B, Y_q)`, i.e. `c_V = T · n_B · dQ2/dT` with a constrained low-T slope fit blended to the table derivative | `StarContext.cpp:807-808`; `CompOSE_Thermo.cpp:711-732`, `:694-706` |
| **Composition input** | `Y_q(r)` from `ChargeFractionYq()` — under **ADR-0001** this is `Σ_i q_i Y_i`, correctly dimensionless | `StarContext.cpp:770-772`, `:691-696` |
| **Units** | `c_V` in erg cm⁻³ K⁻¹; `WV·dr` in km³ scaled by `KM3_TO_CM3 = 1e15` → cm³; product **erg K⁻¹**. Dimensionally consistent with `L` in erg s⁻¹ and `T` in K. | `CompOSE_Thermo.cpp:717-731`; `StarContext.cpp:761`, `:810-811` |
| **Temperature dependence** | Real. `C_⋆` is tabulated on a log-spaced `T∞` grid of 160 points over `[1e-5, 1] MeV` (≈ `[1.16e5, 1.16e10] K`) and interpolated linearly in `C` against `log T` | `StarContext.cpp:775-777`, `:785`, `:730-738` |
| **Radial quadrature** | Trapezoid over the profile grid — consistent with **INV-13**; nominal accuracy O(Δr²) for sufficiently smooth data, so fourth-order behavior is unavailable from this component | `StarContext.cpp:800-815` |
| **Caching / invalidation** | Rebuilt when the profile version changes **or** when the `CompOSE_Thermo` pointer identity changes; otherwise reused. This is INV-12's "rule 3". | `StarContext.cpp:709-717` |

Two properties of this implementation are recorded here so that acceptance cannot later be
misread as endorsement of them:

- **Silent clamping.** Outside the tabulated range, `C_⋆` is returned as the endpoint value
  (`StarContext.cpp:725-728`). This is already flagged as an ungoverned numerical convention in
  **INV-10**.
- **The cache key omits the geometry.** `geo` is an optional argument; when it is null a local
  `GeometryCache` is constructed from the profile (`StarContext.cpp:754-755`), but the cache key
  is only `(profile version, thermo pointer)` (`:712-714`). A subsequent call at the same profile
  version supplying a *different* `GeometryCache` reuses the earlier table. Related to **INV-12**.

### Separately: a confirmed implementation hazard, not part of this decision

`NeutrinoCooling_Details.cpp` **dereferences `ctx.star` at line 889** — the
`HeatCapacityStar_Tinf` call above — while its `if (!ctx.star)` guard sits **twelve lines later
at `:901`**. The guard cannot fire; a null `ctx.star` is a null dereference before it is reached.

This is an **engineering-class defect**, not a heat-capacity semantics question, and it is
recorded as such (`CURRENT_ARCHITECTURE.md` §4; INV-15 *Known implementation hazards*). It
appears here only because any change routing `PhotonCooling` through the same context path will
exercise the same ordering, which is why its correction is scoped into Phase 2A.

## Decision

**CompactStar has exactly one physical heat capacity for the evolved isothermal-interior thermal
degree of freedom:**

> ### `C_⋆(T∞)` — the star-integrated, GR-weighted, temperature-dependent heat capacity.

It is a function of temperature and of stellar structure. The canonical physical definition is

```
C_⋆(T∞)  =  ∫₀^R  c_V( T_local(r), n_B(r), Y_q(r) ) · 4π r² e^{Λ(r)} dr ,
              T_local(r) = T∞ · e^{−ν(r)}
```

with `c_V` supplied by the EOS/CompOSE thermodynamics, in the redshift convention already
governing `ThermalState` (**INV-10**) and the proper-volume measure of **INV-04**.

**The designated implementation of this quantity is the path represented by
`StarContext::HeatCapacityStar_Tinf(...)`** (`StarContext.hpp:151`, implemented at
`StarContext.cpp:704-818`).

### Governing thermal-energy invariant

Every energy channel acting on the thermal degree of freedom divides by **this same** `C_⋆(T∞)`:

```
C_⋆(T∞) · dT∞/dt  =  − L_ν,∞  −  L_γ,∞  +  L_H,∞  +  ⋯
```

subject to the exact redshift convention already governing `ThermalState`. Equivalently, since
the evolved variable is logarithmic (**INV-10**):

```
d/dt ln(T∞ / T_ref)  =  ( − L_ν,∞ − L_γ,∞ + L_H,∞ + ⋯ ) / ( T∞ · C_⋆(T∞) )
```

Terms may be added to the right-hand side — chemical heating `L_H,∞`, BNV heating, and others —
without altering this ADR. What the ADR fixes is that **there is one denominator, and it is
`C_⋆(T∞)`.**

### Status of the designated implementation

**`C_⋆(T∞)` is the accepted physical owner. Its current implementation is not thereby validated.**

This distinction is load-bearing. Acceptance settles *which quantity is physically authoritative*.
It makes no claim that `HeatCapacityStar_Tinf` is numerically correct, that its 160-point grid is
adequate, that its clamping behavior is right, that `CvDensity_cgs_ForCooling` reproduces a
published `c_V`, or that its cache invalidates soundly. None of that has been tested — the
repository has no test suite. Numerical validation is required and is specified in
*Validation requirements* below.

## Rejected current behavior

**A driver-local constant is not the physical heat capacity of production thermal evolution.**

Specifically rejected:

1. **`PhotonCooling::Options::C_eff = 1.0e40 erg/K` as a production physical denominator.**
   It is temperature-independent, structure-independent, EOS-independent, and hand-set at the
   call site with the comment `// Change this!` (`spin_therm_evol_2_main.cpp:245`). It is a
   demonstration knob. It has never been a physical claim, and this ADR forbids it from
   functioning as one.

2. **Any configuration in which one energy channel silently uses a heat capacity different from
   the others.** This is the actual defect: not that `C_eff` exists, but that it serves as the
   *default* normalization for the photon channel while the neutrino channel uses `C_⋆(T∞)`,
   with both summed into the same state element and nothing in the code or the output signalling
   the mismatch.

3. **The `NeutrinoCooling_Details.cpp:403` proposal** to *"reuse the same C_eff as photon cooling
   in a shared thermal config object."* Sharing a constant would make the two channels
   consistent with each other and consistently wrong. Consistency is necessary but not
   sufficient; the shared quantity must be the physical one.

### What is *not* rejected

A constant heat capacity may be **retained as an explicitly selected test/debug approximation** —
for infrastructure tests, for analytic-solution checks where a closed-form cooling curve is
wanted, or for isolating integrator behavior from EOS behavior. Two conditions bind such a
retention:

- it must be **explicitly selected**, never a silent default; and
- when selected it must apply to the **whole thermal balance**, not to one channel, so that the
  evolved equation remains internally consistent even when approximate.

The mechanism that enforces this (an option, a policy object, a context-supplied strategy) is an
architectural question and is deferred — see §6.

## Implementation implications

**No source change is authorized by this ADR.** Everything below is recorded as future work,
scoped to roadmap **Phase 2A**.

### Required for conformance

1. **`PhotonCooling` must stop using a driver-local constant as its production denominator.**
   `PhotonCooling_Details.cpp:320` must divide by the governed `C_⋆(T∞)` obtained through the
   shared context path, not by `drv.GetOptions().C_eff`.
2. **`PhotonCooling::Options::C_eff` must lose its status as the default physical normalization.**
   Whether the field is removed, renamed to mark it as a test override, or subordinated to an
   explicit approximation mode is an API question, settled with §6.
3. **Doxygen must follow.** `PhotonCooling.hpp:55-62`, `:120-123`, `:214-229` and
   `PhotonCooling.cpp:27,36` all document the constant-`C_eff` equation as the driver's physics.
   They become wrong on the day the source is corrected and must be updated in the same change.
4. **`spin_therm_evol_2_main.cpp:245` and `spin_therm_evol_main.cpp:178`** must stop hand-setting
   a physical heat capacity. The `// Change this!` marker is discharged by the correction.
5. **The null-check ordering defect at `NeutrinoCooling_Details.cpp:889` vs `:901`** must be
   corrected as part of the same phase, because routing a second driver through
   `ctx.star->HeatCapacityStar_Tinf(...)` exercises the same unguarded pattern. It is tracked as
   an engineering defect on its own merits, **not** as part of this decision.

### Consequences of correcting the source (recorded, not performed)

The correction is a **scientific-semantic change**: it changes numbers the code produces.
The photon cooling rate is currently `L_γ/1e40`; after conformance it becomes `L_γ/C_⋆(T∞)`.
Whether cooling accelerates or slows depends on `C_⋆` for the star and EOS in use, and the
magnitude is not predictable from this document — it must be measured. Under `GOVERNANCE.md` §2
the change therefore owes: a cited governing authority (this ADR), validation before and after,
and a provenance record.

Because no trustworthy passive-cooling baseline exists to measure against — which is precisely
the circularity this ADR breaks — that validation must be **independent physical verification of
`C_⋆(T∞)` itself**, not agreement with the previous curve. See below.

## Deferred architectural question

**This ADR governs the single heat-capacity denominator. It does not decide where in software
the division by `C_⋆(T∞)` occurs.** Both patterns below satisfy the governing invariant. Choosing
between them is a **structural/architecture** change requiring its own ADR, and it must be
evaluated **after** validation infrastructure exists — not now.

### Pattern A — shared denominator in additive drivers

Each thermal driver computes its own power contribution and obtains the identical canonical
`C_⋆(T∞)` through a shared context/helper:

```
NeutrinoCooling:   −L_ν     / C_⋆
PhotonCooling:     −L_γ     / C_⋆
HeatingFromChem:   +L_H     / C_⋆
```

Because `RHSAccumulator` is additive, these sum to the correct energy equation **provided every
consumer uses the same canonical heat capacity.** This is the smaller change: it preserves the
existing driver contract (`TARGET_ARCHITECTURE.md` §5, principle 3 — pure, additive drivers) and
matches what `NeutrinoCooling` already does. Its weakness is that the invariant is maintained by
convention across N call sites rather than structurally, and it queries `C_⋆` once per driver
per evaluation.

### Pattern B — centralized thermal-energy balance

Physics components expose power contributions only —

```
−L_ν ,  −L_γ ,  +L_H ,  ⋯
```

— and a single thermal-balance owner performs

```
dT∞/dt = L_net / C_⋆(T∞)
```

This makes the single-denominator invariant structural rather than conventional, removes
duplicate heat-capacity queries, and makes future energy-accounting audits substantially easier —
which will matter once chemical heating and BNV heating are present and terms must be shown not
to double-count (`MODERNIZATION_ROADMAP.md` Phase 5, *single-source Γ*). Its cost is real:
it changes architectural ownership, requires a new RHS contract distinguishing "power" from
"state derivative", and touches every thermal driver.

### Why the choice is deferred

Pattern B is a structural change to the driver contract, on the live thermal path, in a
repository with **zero tests and zero baselines** (`CURRENT_ARCHITECTURE.md` §5). Making it now
would violate `AGENTS.md` §5. Pattern A can be adopted immediately as the minimum conforming
change and does not foreclose Pattern B later — a centralized owner can be introduced afterwards
as a behavior-preserving consolidation, measured against baselines that by then exist.

**Recorded for a future ADR. Not decided here.** The roadmap carries it as an open Phase-3
consideration.

## Validation requirements

Acceptance of the *convention* required no numerical validation: it is a physical-ownership
adjudication by the owner. Validation is owed by the **implementation**, and it splits in two.

### V1 — Verification of `C_⋆(T∞)` itself (roadmap Phase 2A)

Narrowly targeted; **must not** be measured against the existing passive-cooling curve, which is
known to encode the rejected behavior.

1. **Dimensional and unit check.** Confirm end-to-end that `c_V` [erg cm⁻³ K⁻¹] × proper volume
   [cm³] yields `C_⋆` [erg K⁻¹], including the `KM3_TO_CM3 = 1e15` factor
   (`StarContext.cpp:761`).
2. **Analytic limit.** For a degenerate free Fermi gas, `c_V ∝ T` at fixed composition, so
   `C_⋆(T∞)` should be very nearly linear in `T∞` in the low-temperature regime. Verify the
   tabulated `C_⋆` reproduces that slope over the log-spaced grid.
3. **Order-of-magnitude check against literature.** For a canonical ~1.4 M⊙ npeμ star, published
   total heat capacities are of order `10³⁷–10³⁸ erg K⁻¹` at `T ≈ 10⁸ K`. Compare. This check
   alone distinguishes `C_⋆` from the `1e40` placeholder and is the single most informative test
   available before any baseline exists.
4. **Grid convergence in `r`.** The integral is trapezoidal over the profile grid (INV-13), whose
   nominal order is O(Δr²) for sufficiently smooth data. **Measure** the observed order under
   refinement; do not assume it. A departure from the nominal order is itself a finding — the
   profile grid is set by the TOV integrator rather than chosen for refinement, and `c_V` reaches
   the integrand through table interpolation and domain clamping.
5. **Grid convergence in `T`.** Vary `NT` from 160 (`StarContext.cpp:777`) and confirm the
   interpolated `C_⋆` is insensitive to the tabulation density over the range of interest.
6. **Clamping behavior stated, not assumed.** Record explicitly what happens for
   `T∞ < 1e-5 MeV` (≈ `1.16e5 K`) and `T∞ > 1 MeV` (≈ `1.16e10 K`), where the endpoint value is
   returned silently (`StarContext.cpp:725-728`). Governing that clamp remains INV-10 work.
7. **Cache correctness.** Confirm a rebuild on profile-version change; record the behavior when
   the same profile version is queried with a different `GeometryCache`, which the current key
   does not distinguish (`StarContext.cpp:712-714`). INV-12.

### V2 — Thermal baseline (roadmap Phase 2B, after V1)

Once `PhotonCooling` conforms and `C_⋆` has passed V1, capture the passive-cooling regression.
It then records a **physically coherent energy equation** rather than deliberately preserving a
placeholder — which is the whole point of ordering the two phases this way.

## Consequences for existing passive-cooling outputs

Stated plainly, because it affects committed results.

- **Every passive-cooling result produced to date used the rejected behavior.** Both live
  programs set `C_eff = 1.0e40` by hand (`spin_therm_evol_2_main.cpp:245`,
  `spin_therm_evol_main.cpp:178`), so every run in which photon cooling was active integrated
  `dT∞/dt = −L_γ/1e40 − L_ν/C_⋆(T∞)`.
- **Those outputs are not a validation reference and must not be treated as one.** They are
  retained as historical artifacts of infrastructure development. Under `GOVERNANCE.md` §1 rule 7,
  that the code produced them is evidence of behavior, not of correctness.
- **No golden regression may be captured from them.** This is the operative prohibition: it is
  what prevents the circularity from being resolved by freezing the placeholder.
- **The magnitude of the change is not asserted here.** How much the corrected curve differs is
  an empirical question for Phase 2A, not a claim this ADR is entitled to make.
- **Nothing outside the thermal path is invalidated.** TOV structure, the O(Ω) moment of inertia,
  and the Direct-Urca threshold do not depend on `C`. Neutrino *luminosity* is likewise
  unaffected — only its conversion to `dT∞/dt` was ever in question, and `NeutrinoCooling`
  already conforms.

## Consequences

Now in force:

- **INV-15 moves from UNRESOLVED to GOVERNED (ACCEPTED)** as to *governing physics*, citing this
  ADR. Its **source conformance** and **numerical validation** remain explicitly incomplete and
  are tracked separately in `docs/SCIENTIFIC_INVARIANTS.md`.
- **`CURRENT_ARCHITECTURE.md` §3 drops heat capacity from its unresolved-authority list** — the
  *ownership* question is decided. Its driver table and pipeline diagram continue to report the
  **live nonconformance** honestly, because no source correction has landed.
- **The roadmap circularity is repaired.** Phase 2 splits into **Phase 2A** (a deliberately narrow
  pre-baseline correctness phase: photon-cooling conformance, the adjacent null-check ordering
  defect, and targeted `C_⋆` verification) and **Phase 2B** (the validation baseline). Phase 3
  loses the heat-capacity ownership item, since ownership is now governed and minimum conformance
  is a Phase-2A prerequisite; Phase 3 retains only the *architectural* question of §6.
- **Phase 2A is not a general cleanup phase.** It admits only work without which a baseline would
  be scientifically misleading. Unrelated scientific corrections — INV-07 Hartle normalization,
  INV-16 muon channel, the placeholder emissivity normalizations, INV-04 proper-volume
  consolidation — stay in their own phases. Phase-2A changes are governed by their own
  scientific/numerical class and require **independent physical checks**, not comparison against
  a known-wrong passive-cooling curve.
- **A constant heat capacity is forbidden as a silent default** anywhere on the production
  thermal path, and forbidden in all cases from applying to one channel while others use `C_⋆`.
- **No source change is authorized by this ADR.**

## Alternatives

*Retained for the record. Alternative 1 was selected by owner adjudication.*

### Alternative 1 — `C_⋆(T∞)`, the GR-integrated EOS-based stellar heat capacity ✅ **SELECTED**

- **Statement.** One physical heat capacity, temperature- and structure-dependent, computed as
  the proper-volume integral of the EOS `c_V` at the Tolman local temperature. Every channel
  divides by it.
- **Required code changes.** `PhotonCooling_Details.cpp:320` routes to the canonical `C_⋆`;
  `PhotonCooling::Options::C_eff` loses default-physical status; Doxygen at
  `PhotonCooling.hpp:55-62,120-123,214-229` and `PhotonCooling.cpp:27,36` corrected; both test
  mains stop hand-setting `C`; `NeutrinoCooling_Details.cpp:889/:901` ordering corrected.
- **Migration risk.** Moderate and **visible**. It changes the photon cooling rate by a factor
  that must be measured, on the live thermal path, with no baseline. Mitigated by validating
  `C_⋆` independently (V1) rather than against the old curve.
- **Validation needed.** V1 in full, before V2.
- **Implications for existing outputs.** Every passive-cooling run to date is superseded as a
  reference. Nothing outside the thermal path is affected.

### Alternative 2 — a shared constant `C_eff` for all thermal channels

- **Statement.** Keep a constant heat capacity but give it one owner in `DriverContext` or
  `EvolutionConfig`, so all channels at least share it. Essentially
  `NeutrinoCooling_Details.cpp:403`.
- **Required code changes.** Smallest: hoist `C_eff` out of `PhotonCooling::Options`, have
  `NeutrinoCooling` read it instead of calling `HeatCapacityStar_Tinf`.
- **Migration risk.** Low mechanically — and that is the trap. It produces a *self-consistent*
  energy equation that is still not the star's physics, and it would look correct in every
  internal consistency check.
- **Validation needed.** None possible in principle: there is no physical reference a constant
  `C` can be validated against for a real star.
- **Implications for existing outputs.** Would make the neutrino channel *worse*, discarding a
  real GR integral in favor of a knob. **Rejected.**

### Alternative 3 — per-channel heat capacities, made explicit rather than removed

- **Statement.** Accept that each channel may carry its own `C`, and require only that each be
  documented.
- **Required code changes.** Documentation only; the code already behaves this way.
- **Migration risk.** None — because nothing changes.
- **Validation needed.** Not applicable.
- **Implications for existing outputs.** Preserves them. But `dT∞/dt` is a single scalar equation
  for a single thermal degree of freedom; there is no physical reading in which its denominator
  depends on which loss term is being summed. This alternative documents an inconsistency instead
  of resolving it. **Rejected.**

### Alternative 4 — defer until a baseline exists

- **Statement.** Capture the Phase-2 baseline as-is, then correct heat capacity in Phase 3 and
  measure the delta.
- **Required code changes.** None now.
- **Migration risk.** **High and structural.** It is the circular dependency: a known
  scientifically inconsistent placeholder would be frozen into a golden regression purely to
  satisfy process, and every subsequent correction would register as a regression against it.
- **Validation needed.** N/A.
- **Implications for existing outputs.** Would elevate known-wrong outputs to reference status.
  **Rejected** — this is the specific failure mode ADR-0002 exists to prevent.

## Provenance

The dual heat-capacity ownership was identified during Phase-0 reconnaissance (2026-08-31) and
recorded as INV-15 and as item 8 of the reconnaissance findings
(`docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md:240`). It was re-authenticated
independently for this ADR against `modernization/governance-foundation` at `ba49e10`: both
driver paths were re-read end to end, the `StarContext` heat-capacity implementation was traced
through `BuildHeatCapacityCache_` and its CompOSE `c_V` source, units were verified through the
`KM3_TO_CM3` scaling, and the previously reported `NeutrinoCooling` null-check ordering issue was
confirmed present at `:889` versus `:901`.

**Implementation-conformance update, 2026-08-31 (status only).** The nonconformance this ADR
recorded has been **corrected**. Roadmap Phase 2A-1/2A-2 produced the §V1 evidence
(`docs/validation/HEAT_CAPACITY_V1.md`, **V1 VERIFIED**), and Phase 2A-3 changed `PhotonCooling`
to divide by the canonical `C_⋆(T∞)` from `StarContext::HeatCapacityStar_Tinf`, removing
`PhotonCooling::Options::C_eff` outright. **Pattern A was preserved and the §6 architectural
question remains deferred.** The adjacent `NeutrinoCooling` null-check ordering defect was
repaired in the same increment.

The Context, Evidence, Rejected-current-behavior and Consequences sections above are **left
unchanged as the historical record of why the correction was required** — they describe the state
of the code at acceptance, not today's. **`GOVERNANCE.md` §3.1 condition 7 — SATISFIED.** Phase 2B-1 was blocked by two pre-existing
defects unrelated to this ADR (an unusable `MSBDF` default with no Jacobian, and an
`EvolutionSystem` block that made a Spin state mandatory). Both were repaired in Phase 2B-1R, and
the passive-cooling regression baseline was established immediately afterward:
`tests/baselines/passive_cooling_cmf_1p6_debug.tsv`, CTest `passive_cooling_regression`, evidence
in `docs/validation/PASSIVE_COOLING_BASELINE.md`. Every §3.1 condition for the Phase-2A-3
correction is now discharged.

**Post-acceptance clarification, 2026-08-31 (editorial).** During the Phase-0.5 governance
coherence audit, two statements of *numerical* wording in this ADR — the radial-quadrature row of
the evidence table and item 4 of §V1 — asserted second-order convergence in Δr as a requirement
rather than as the quadrature's nominal order to be measured. They were corrected to match the
amended INV-13. **No part of the Decision, the rejected behavior, the deferred architectural
question, or the consequences was altered**, and the correction was necessary because an ACCEPTED
ADR outranks `SCIENTIFIC_INVARIANTS.md` (`GOVERNANCE.md` §1) and would otherwise have re-imposed
the assumption the invariant register now rejects.

**Adjudicated by the project owner on 2026-08-31.** The owner fixed the physical convention —
one canonical `C_⋆(T∞)`, the GR-integrated EOS-based quantity — and simultaneously directed that
the *software* question of where the division occurs be held open pending validation
infrastructure. The agent drafted the alternatives, gathered and verified the evidence, and
recorded the decision; it did not select among the alternatives.

Per `GOVERNANCE.md` §5, acceptance of this ADR ratifies **only** the heat-capacity ownership
convention and the governing thermal-energy invariant. It does **not**:

- certify that `StarContext::HeatCapacityStar_Tinf` is numerically correct;
- ratify the placeholder neutrino emissivity normalizations that feed `L_ν`
  (`NeutrinoCooling_Details.cpp:100-102`);
- ratify the Hartle O(Ω²) or `RotochemicalCache` candidate code from `675b4a9`;
- decide the thermal-balance architecture (§6);
- confer accepted status on any other DRAFT governance document.
