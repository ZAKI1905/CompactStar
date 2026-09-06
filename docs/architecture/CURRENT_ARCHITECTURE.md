# CompactStar — Current Architecture

> **Unit-1D candidate (2026-09-06): ADR-0012 A1 production correction IMPLEMENTED + CANDIDATE VALIDATED.**
> Focused 42/42 and data-free 30/30 pass; full raw 50 PASS / 3 FAIL (rc=8),
> classified A=50, B=3, C=0. Every Category-B output equals its validated successor exactly.
> C1–C6 are repeatable successor candidates, not governed baselines; H7/H8 remain byte-identical.
> Independent review: PENDING. Human ratification: PENDING. Governed baseline supersession: PENDING.
> ADR-0011 Q4: NUMERIC UNIT PREREQUISITE SATISFIED IN CANDIDATE; NOT YET CLOSED; PB7 incomplete.
> INV-09 remains INTENDED BUT UNVERIFIED (unresolved); INV-11 remains UNRESOLVED.
> No Phase-5B production, Btilde, paper Z/W, evolution, or BNV work begins.
> Evidence: `docs/validation/RELATIVISTIC_UNIT_BOUNDARY_IMPLEMENTATION.md:1`.
> Earlier entries below retain their historical scope.


> **Track-R Structure-1 (2026-09-05): HUMAN-RATIFIED WITH QUALIFIED CLAIM — STATIC FREE-GAS STRUCTURE AND COMMON-STATE FR2005 TABLE-1 NUMBERS VERIFIED.** At the source's printed `rho_c = 1.10e15 g cm^-3`, one common central state reproduces `M`, `R`, and `R_infinity` with no source-bin widening, mass fitting, or EOS tuning and with independent TOV validation. The source-qualified `M_max` selection semantics remain unresolved; the continuous pre-Sigma supremum is not claimed as reproduced. INV-09/INV-11 and global response/evolution scope are unchanged. Evidence: `docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_RATIFICATION.md:1`.

> **Historical Phase 5A-5 state (2026-09-05): TRACK-R FREE-GAS LOCAL THERMODYNAMIC COVERAGE HUMAN-RATIFIED AND COMPLETE — WHOLE-STAR STRUCTURE REPRODUCTION READY.** Typed active results cover
> vacuum, smooth 1D p-e, value-only neutron onset, smooth 2D npe, value-only muon onset, and
> smooth 3D npe-mu through the pre-Sigma-minus ceiling, with explicit numerical-resolution
> failures. The code is ratified as scientific authority for this Track-R local source domain.
> No whole-star structure or global response has been computed; P-7 remains a required planning
> input for the separate structure task. Evidence: `docs/validation/PHASE5A5_TRACKR_PE_BRANCH.md:63-280`;
> `docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:771-828`;
> `docs/validation/PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:3-46`, `:94-135`.
>
> **Historical Phase 5A-3 state (2026-09-04): Track-R cold free-gas LOCAL model implemented and
> analytically validated.** `TrackRFreeGasThermodynamicProvider` supplies the source-model cold
> ideal `n,p,e,mu` energy, neutral conjugates, constrained Hessian, active-branch equilibrium
> recovery, muon onset, and the pre-`Sigma-minus` domain endpoint. It is not connected to a star,
> rotation, global susceptibility, coefficient, rate, or evolution path. Evidence:
> `docs/validation/PHASE5A3_TRACKR_FREEGAS_LOCAL_PROVIDER.md`.
>
> **Current Phase 5A-2 state (2026-09-04): local thermodynamic contract IMPLEMENTED AND
> VERIFIED.** `CompactStar/EOS/LocalThermodynamics.*` implements the accepted ADR-0010 generic
> cold charge-neutral boundary and analytic free electrons/muons. A test-only analytic toy and
> independent full-intrinsic projection fixture pass V1-V10. The provider is not connected to a
> star or evolution. No APR, DS(CMF) off-equilibrium provider, stellar susceptibility, paper
> `Z/W`, evolution, superfluidity, or BNV is implemented. Evidence:
> `docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_IMPLEMENTATION.md`.
>
> **Historical Phase 4E closeout (2026-09-04): PHASE 4 ROTATION CORRECTNESS COMPLETE —
> PHASE-5 STRUCTURAL INTERFACE RATIFIED.** The existing normalized
> `HartleFirstOrderResponse` + `HartleMonopoleResponse` are the supported structural inputs;
> no wrapper, new physics or baseline change. INV-08 is CLOSED/VERIFIED only for ordinary
> NStar, fixed-ε_c, l=0 O(Ω²) on ADR-0009 backgrounds. INV-09 remains unresolved.
> **PHASE 5 NOT YET BEGUN.** Normal reviewed branch integration precedes a subsequent
> Phase-5 task; no merge here. Contract, eight unchanged hashes and final serial suites:
> `docs/validation/PHASE4_CLOSEOUT.md:64`, `docs/validation/PHASE4_CLOSEOUT.md:158`, `docs/validation/PHASE4_CLOSEOUT.md:183`, `docs/validation/PHASE4_CLOSEOUT.md:219`.
> Earlier dated entries below preserve their historical scope, including past blockers.


> **Phase 4D-BL (2026-09-04): first verified Hartle monopole regression baseline established.**
> Scientific status remains `HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED`.
> Complete serial suites: **45/45 full, 679.30 s; 22/22 self-contained, 91.38 s**.
> Eight current artifact hashes: `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:250`; producer/repeatability: `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:1`.
> Scientific blocker for Phase 4E **CLEARED**; closeout/field ratification remains the next task.
> Phase 5 has not begun. Earlier dated entries below retain their historical scope.

> **Phase 4D-DA (2026-09-04):** owner adjudication of ADR-0008 Validation D — the monotonicity clause
> was a heuristic the accepted profile-partition measure cannot satisfy (`O(h)` phase-dependent location
> term, confirmed 6/6 on the real star); clarified into D′1–D′5, all met by the 4D-RV evidence. Status
> `HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED` (`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md` §21). **No baseline yet — the first monopole baseline is the next, separate task.**
>
> **Phase 4D-RV (2026-09-03):** the corrected independent Hartle-monopole revalidation ran on the
> migrated ADR-0009 backgrounds — status `HARTLE O(OMEGA^2) MONOPOLE RESPONSE CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE`
> (`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md`). Every physics line passes (independent Stieltjes-measure oracle `≤ 5.3e-5`, published
> tables reproduced, M1–M10 load-bearing, production diff NONE); ADR-0008 Validation D's
> monotonicity clause is NOT MET at the node-placement floor (`δM̂` spread `3.7e-5`, differences
> non-monotone). **No monopole baseline; owner adjudication of D required; Phase 5 remains blocked.**
>
> **TOV-SURF-MR (2026-09-03):** ADR-0009 is accepted, source conformed and
> numerically validated; the seven surface artifacts are now migrated. The ordinary
> primitive retains its accepted-solution `p=p_cut` event and fail-closed completion;
> cutoff, target-mass tolerance and inherited driver history are unchanged.
> Current hashes and producer reproduction: `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:889`.
> Main 41/41 and self-contained 20/20 pass; corrected Phase-4D independent revalidation
> is ready but has not run. The first monopole baseline and Phase 5 remain blocked.
> Evidence: `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:960`, `docs/validation/TOV_SURFACE_ARTIFACT_MIGRATION.md:1003`;
> source validation: `docs/validation/TOV_SURFACE_IMPLEMENTATION.md:339`.

> **STATUS: DESCRIPTIVE.** Authoritative for component boundaries and ownership
> (`GOVERNANCE.md` authority rank 6). Describes **only** behavior that is compiled and reachable.
>
> **Current snapshot:** `11ffe45` — the head of roadmap Phase 1 (increments 1A–1D).
>
> **Evidence scope is deliberately mixed, and the layers must not be conflated:**
>
> | Claim area | Authenticated at | How |
> |---|---|---|
> | Build, configuration, warning policy, test plumbing (§5) | **`11ffe45`** | Re-authenticated during this Phase-1 sync against the actual CMake files and a clean configure/build/test |
> | `RotationSolver` merge integrity and buildability | **`57334d8`** | Re-authenticated during Phase 1B against the three-way merge `9f70f14` = `3639d71` + `e60e656` |
> | `RotationSolver` API reachability, provenance, unit annotations | **`df859b5`** | Phase 4A-0 entry audit — `docs/validation/PHASE4_ROTATION_ENTRY.md` |
> | O(Ω²) monopole system — primary-source derivation, surface audit, governed replacement (ADR-0007) | **`bb073c8`** | Phase 4C-G — `docs/validation/PHASE4C_HARTLE2_DERIVATION.md` |
> | EOS thermodynamic derivative `dε/dp` — authority, units, domain semantics, profile delivery | **`13111de`** | Phase 4C-I0 — `docs/validation/PHASE4C_I0_EOS_DERIVATIVE.md` |
> | O(Ω²) monopole — candidate retirement and governed replacement, provenance, materialization | **`a97f9c5`** | Phase 4C-I1 — `docs/validation/PHASE4C_I1_MONOPOLE_IMPLEMENTATION.md` |
> | O(Ω²) monopole — independent physical validation (analytic, continuum, DS(CMF)-1, Chandrasekhar–Miller 1974, Hartle–Thorne 1968) | **`377bc4a`** | Phase 4D — `docs/validation/PHASE4D_MONOPOLE_VALIDATION.md` — `HARTLE MONOPOLE VALIDATION FAILED` on stepped crusts (implementation verified; ADR-0007 amendment required) |
> | O(Ω²) monopole — EOS energy-density source as a measure (derivation, scratch evaluation of three numerical authorities) | **`246f3f2`** | Phase 4D-RG — `docs/validation/PHASE4D_R_EOS_MEASURE_DERIVATION.md` — ADR-0008 **PROPOSED** |
> | O(Ω²) monopole — measure-complete EOS source implemented (per-segment measure, terminal atom, `dε/dp` retained for the centre series) | **`8abbab4`** | Phase 4D-RI — `docs/validation/PHASE4D_RI_EOS_MEASURE_IMPLEMENTATION.md` — ADR-0008 **ACCEPTED** |
> | O(Ω²) monopole — corrected independent revalidation on the migrated ADR-0009 backgrounds (independent Stieltjes-measure oracle, ADR-0008 A–J, M1–M10) | **`218c9aa`** | Phase 4D-RV — `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md` — **CHARACTERIZED; ADR-0008 D monotonicity NOT MET; no baseline** |
> | O(Ω²) monopole — adjudication of ADR-0008 Validation D monotonicity (numerical analysis, phase prediction, computed envelope) | **`42b34ac`** | Phase 4D-DA — `docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md` §21 — **VERIFIED; first baseline next** |
> | First-order rotation API, units and normalization | **Phase 4A** | ADR-0006 implemented and validated — `docs/validation/PHASE4A_FIRST_ORDER_NORMALIZATION.md` |
> | First-order rotation **physics** (profile shape) | **Phase 4B** | Independently verified — `docs/validation/PHASE4B_FIRST_ORDER_PHYSICS.md` |
> | Everything else — EOS, TOV, Evolution, Drivers, Microphysics, caches | `9f70f14`, from the `d91c31b` audit | Phase-0 reconnaissance; see `docs/reconnaissance/2026-08-31-phase-0-reconnaissance.md` |
>
> **The repository was not re-audited in full at `11ffe45`.** Phase 1 was a build/test phase and
> touched no scientific source. Unchanged-subsystem claims still carry their original Phase-0
> evidence scope and its caveats.
>
> **Historical base commit:** `9f70f14` (`master` merge of owner work `3639d71` + April candidate
> `675b4a9`) — retained because most component-level evidence below still derives from it.
>
> For intended design, see [`TARGET_ARCHITECTURE.md`](TARGET_ARCHITECTURE.md). **Do not read
> intent from this document.**

## Scientific source library

- Shared source root is `/Users/keeper/Documents/CompactStar/literature`; it is outside Git and shared
  by all worktrees.
- Before re-deriving established scientific material, consult `/Users/keeper/Documents/CompactStar/literature/catalog.tsv` and
  `/Users/keeper/Documents/CompactStar/literature/README.md`.
- PDFs under this root are read-only scientific references; do not rename, move, or edit them in
  normal research tasks.
- `literature/_incoming` is staging/provenance only and is not scientific authority.
- `catalog.tsv` role labels define source authority; later entries may explicitly supersede earlier
  treatments.
- User-authored BNV papers are in `literature/bnv/zakeri`.
- Future scientific work should name the exact source PDF used when equations or decisions depend
  on it.

## Labels used here

| Label | Meaning |
|---|---|
| **LIVE** | Compiled, reachable, and exercised by a program that runs |
| **COMPILED, UNEXERCISED** | Built into the library; no program calls it |
| **UNREACHABLE SCAFFOLDING** | Compiled but cannot be called — no public accessor |
| **PUBLIC, ZERO CALLERS** | Compiled and callable through the public API; no program or test in the repository calls it |
| **DECLARED, UNDEFINED** | Declared in a public header with no definition in any translation unit; any caller fails to link |
| **NOT COMPILED** | Source exists; absent from every CMake source list |
| **EMPTY** | Zero-byte or body-less file |
| **CANDIDATE** | Present but scientifically unverified (`GOVERNANCE.md` §5) |

---

## 1. What actually executes

One program exercises the full pipeline end to end: `main/Test/spin_therm_evol_2_main.cpp`.

```
CompOSE EOS table
   ↓
TOVSolver::SolveToProfile → SingleStarSolveToTOVPoints     LIVE
   │   integrates in CGS with explicit G, c  (INV-02)
   ↓   unit conversion happens once, at NStar::Append
NStar::prof_ : StarProfile        LIVE — canonical star representation
   │   geometric units: km, km⁻²  ·  versioned (m_version + EditScope)
   ↓
RotationSolver::FindNMomInertia → SeqPoint::I              LIVE, O(Ω) only
   │   raw internal solve from an arbitrary seed (private; never published)
   └─► HartleFirstOrderResponse   { I, ω̄/Ω, ω̄'/Ω }          LIVE — SEED-FREE (ADR-0006)
          via NStar::RotationResponse()
          NO implicit physical spin: construction confers no Ω
             │
             └─ + explicit AngularVelocity [rad/s]  ──►  PhysicalFirstOrderRotation
                   { Ω [km⁻¹], J [km²], I [km³], ω̄(r), ω̄'(r) }   via NStar::RotationAt()
                   a scaling of the response — NOT a second ODE solve
   ↓
StarContext(prof_)                LIVE
   ├─→ GeometryCache               LIVE — deep copy, NO version gate (INV-12)
   ├─→ DirectUrcaMask              LIVE — hand-rolled gate
   ├─→ HeatCapacity C(T∞)          LIVE — third, different gate
   └─→ MassDensity g/cm³           LIVE
   ↓
DriverContext {star, geo, envelope, thermo, cfg}    5 raw non-owning pointers
   ↓
StateLayout{Spin:1, Thermal:1}    ← n_eta = 0; the Chem block is never configured
   ↓
EvolutionSystem                   LIVE — plain registration-order loop
   ├─ MagneticDipole      LIVE      dΩ/dt = −K·sign(Ω)·|Ω|ⁿ
   ├─ NeutrinoCooling     LIVE      L_ν ÷ C_⋆(T∞)  ← governed denominator (ADR-0002)
   │                                placeholder Q₀ normalizations
   └─ PhotonCooling       LIVE      L_γ ÷ C_⋆(T∞)  ← same denominator
                                    CONFORMING with ADR-0002 (Pattern A)
   ↓
GSLIntegrator (RKF45, rtol 1e-6 / atol 1e-10)       LIVE
   │   default was MSBDF and unusable — implicit, no Jacobian (Phase 2B-1R)
   │   MSBDF now REJECTED with a diagnostic until Jacobian support exists
   ↓
TimeSeriesObserver + DiagnosticsObserver            LIVE
```

**Everything not on this diagram is scaffolding.**

---

## 2. Component status

### EOS — Phase 5A local thermodynamics

| Component | Status | Note |
|---|---|---|
| `ILocalThermodynamicProvider` and value types | **LIVE — ADR-0010 CONFORMED LOCALLY** | Independent cold charge-neutral `x=(n_B,n_e,n_mu)` interface returning validated species reconstruction, energy density, `g=(mu_n,-eta_npe,-eta_npmu)`, and `ChargeNeutralChemicalHessian H=partial g/partial x`. Metadata declares identity, composition, conventions, and smooth domain. No star dependency, intrinsic charged-potential requirement, projector, inverse, or paper-`B` API |
| `ActiveLocalThermodynamicEvaluation`, active-chart and boundary value types | **LIVE — HUMAN-RATIFIED FOR TRACK-R LOCAL DOMAIN** | Variant of full 3x3, explicit npe 2x2, explicit p-e 1x1, value-only neutron/muon thresholds, and value-only vacuum. Each alternative declares active particles, response dimension, and smooth/threshold/vacuum status; boundary types have no Hessian. `Evaluate` retains its full-chart meaning; `EvaluateActive`/`EquilibriumAt` expose typed alternatives. `CompactStar/EOS/LocalThermodynamics.hpp`; `docs/validation/PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:27-46`, `:54-74` |
| `ColdRelativisticIdealFermion` | **LIVE — ANALYTICALLY VALIDATED; HUMAN-RATIFIED WITHIN TRACK-R SOURCE MODEL** | Source-independent n,p,e,mu spin-1/2 primitive; zero-density chemical potential is rest mass and derivative is unavailable. `CompactStar/EOS/src/LocalThermodynamics.cpp:104`; `docs/validation/PHASE5A5_TRACKR_PE_INDEPENDENT_REVIEW.md:134-220`; `docs/validation/PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:27-46` |
| `TrackRFreeGasThermodynamicProvider` | **LIVE — HUMAN-RATIFIED WHOLE-STAR LOCAL COVERAGE AND QUALIFIED STATIC STRUCTURE** | Cold ideal vacuum/p-e/npe/npe-mu source model; separate 1D/2D/3D active charts, value-only thresholds, explicit near-onset numerical failures, and the pre-Sigma-minus endpoint. The additive `BarotropeAt` capability supplies equilibrium values independently of chemical H. The Structure-1 generated-table/canonical-TOV path and common-state FR2005 Table-1 numbers are human-ratified at Levels 1 and 2, without global response or evolution; the existing response guards remain unchanged. Level 3 is not ratified: source-qualified `M_max` selection semantics remain unresolved. `docs/validation/TRACKR_FREEGAS_WHOLESTAR_STRUCTURE1_RATIFICATION.md:1`; `CompactStar/EOS/TrackRFreeGasThermodynamics.hpp`; `docs/validation/PHASE5A5_TRACKR_LOCAL_RATIFICATION.md:27-46`, `:94-135` |
| `ColdRelativisticFreeLepton` | **LIVE — ANALYTICALLY VALIDATED** | T=0 relativistic free electrons/muons using repository constant authority; returns total chemical potential and energy density, with analytic `dmu/dn` only for positive density. At zero density the value limit is explicit and the derivative is unavailable/fail-closed |
| Phase 5A-2 analytic toy | **TEST-ONLY FIXTURE** | Positive stable charge-neutral potential with nonzero cross derivatives; validates V1-V10 and creates no scientific baseline. It is not APR, DS(CMF), BPAL, or a neutron-star EOS |

### Core

| Component | Status | Note |
|---|---|---|
| `StarProfile` | **LIVE** | Canonical for `NStar`. Column enum + `m_version` + `EditScope` RAII |
| `StarProfileView` | **COMPILED, UNEXERCISED** | `NStar::View()` has zero callers |
| `NStar` | **LIVE** | Migration to `prof_` complete; legacy `ds` is commented out, not dual-live |
| `MixedStar` | **COMPILED, UNEXERCISED** | No surviving `main/` uses it. Master-grid totals added by `3639d71` |
| `TOVSolver` | **LIVE** | Two live integration paths — see §3 |
| `TOVSolver_Thread` | **COMPILED, UNEXERCISED** | Bookkeeping subclass, 124 lines |
| `RotationSolver` — O(Ω) | **LIVE** | Runs on every star build; feeds `SeqPoint::I`. Its profile-backed interpolation path was restored in Phase 1B — see below. **`I = J/Ω` VALIDATED as a scale-free observable (Phase 2B-4B)**. **Phase 4A: the first-order normalization is GOVERNED and CONFORMED (ADR-0006).** The arbitrary seed is private (`RotationSolver.hpp:326`, default `5e-3`) with no public setter; construction publishes the **seed-free** `HartleFirstOrderResponse` (`I`, `ω̄/Ω`, `ω̄'/Ω`) and **no implicit physical spin**; a physical solution requires an explicit `AngularVelocity` in rad s⁻¹. `I` is bit-identical across the change. **Phase 4B: the normalized response is INDEPENDENTLY VERIFIED as physics** — node-by-node against an independently normalized profile (`2.9e-9` analytic, `≤ 2.3e-5` on the CMF sequence), plus exterior-matching and volume identities and derived weak-field coefficients |
| `CompactStar::AngularVelocity` | **LIVE** | Top-level, dependency-neutral typed physical angular velocity (`CompactStar/AngularVelocity.hpp`). Factories `FromRadPerSecond` / `FromHz` / `FromPeriodSeconds`; **no factory accepts km⁻¹**. Together with `AngularVelocityGeomToRadPerSecond()` it is the **sole owner of the physical ↔ geometric angular conversion** on the governed path, using `Zaki::Physics::LIGHT_C_KM_S` |
| `TOVSolver` — EOS derivative `dε/dp` | **LIVE — SOLE AUTHORITY (ADR-0007 P5)** | `GetEDensDeriv` (`TOVSolver.hpp`) returns the **dimensionless** barotropic `dε/dp = 1/c_s²` as the analytic derivative of the **same** `visi_eps_p_spline` (`gsl_interp_steffen`) that `GetEDens` reads, through the same accelerator; the cgs→geometric factor `INV_FM4_2_Dyn_CM2/INV_FM4_2_G_CM3` (= `c²`, asserted to `1.5e-16` against an independent literal) is applied in exactly one place. **Fail-closed**: out-of-domain, non-finite input, or no interpolant returns NaN — never a clamped boundary value, and never `0.0`, which is the physical value for incompressible matter. Profile finite differencing is **not** an authority (measured 155–490 % error at ~25 crust nodes per star). Phase 4C-I0 |
| `StarProfile` — `HasEosDEdP()` / `GetEosDEdP()` | **LIVE** | Profile-attached carrier for the above, filled through `TOVPoint::dedp` by **both** ordinary-`NStar` construction paths (`BuildFromTOV`, and `Append` + `FinalizeSurface`). Deliberately **not** a `radial` column: the schema, the species-column indices and the export layout are unchanged. All-or-nothing — a set is published only if every node carries a finite value, so a consumer may trust the flag without inspecting nodes. Absent by default on point-constructed stars, which may supply their own explicitly (`0.0` for constant density); first-order and TOV physics neither read it nor require it |
| `HartleFirstOrderResponse`, `PhysicalFirstOrderRotation` | **LIVE — NORMALIZED RESPONSE RATIFIED AS PHASE-5 INPUT (4E)** | `docs/validation/PHASE4_CLOSEOUT.md:64`, `docs/validation/PHASE4_CLOSEOUT.md:115`. The two first-order result types (`RotationSolver.hpp:159`, `:110`). The response is seed-free by construction; the physical solution stores **one canonical geometric `Ω`** with a named `OmegaRadPerSecond()` accessor and no duplicated state (ADR-0006 Q3) |
| `RotationSolver` — O(Ω²) monopole | **LIVE — MEASURE-COMPLETE GOVERNED SOURCE (ADR-0007 + ADR-0008); INDEPENDENTLY PHYSICALLY VERIFIED ON THE ADR-0009 SURFACE (4D-RV + 4D-DA); FIRST TRUSTED REGRESSION BASELINE ESTABLISHED (4D-BL); PHASE-5 STRUCTURAL INTERFACE RATIFIED (4E)** | **Phase 4C-I1 (2026-09-03) executed the `GOVERNANCE.md` §3.1 correction ADR-0007 §9 authorized:** the AI-authored candidate of `675b4a9` — `SolveHartle2_N`, `ODE_Hartle2_N_Fast`, the two MixedStar stubs, `HartleResult`, `GetHartleResult` and every proven candidate-only member — was **deleted**, and the governed response replaced it in the same commit (zero live definitions, zero public declarations, zero compiled callers remain). `ODE_HartleMonopole_` integrates the ADR-0007 P2 system in `m̂₀ = m₀/Ω²` and `p̂₀* = p₀*/Ω²`, driven by the **verified seed-free** `ω̄/Ω` and `ω̄'/Ω`; every background input is interpolated at the **actual** ODE radius through one shared bracket (INV-13); `1 − 2m/r` comes from `Geometry::MetricDenominator` (ADR-0004, no clamp); `dε/dp` comes only from the Phase-4C-I0 profile authority and its absence **fails closed**; the start is the regular-centre series (no shooting, no surface condition, no homogeneous admixture); `δM̂ = m̂₀(R_*) + 4πR_*²ε_*ξ̂₀(R_*) + I²/R_*³` with `R_*` the EOS-floor node, never labelled `p = 0`. Published atomically with source-profile identity and `Version()` provenance; a stale response is never returned. **Measured:** seed invariance `7.85e-15` (bound `1e-10`), exact quadratic materialization, `I` bit-identical, seven artifacts unchanged. **Phase 4D (2026-09-03):** recomputed by an independent `(m₀, h₀)` solver and a continuum solver — analytic agreement `9.7e-9` converging as `h²`, centre series `5e-10`, Newtonian intercepts `1e-6`, DS(CMF)-1 `≤ 3.8e-7`, Chandrasekhar–Miller 1974 19/19 to `7.3e-4`, Hartle–Thorne 1968 8/8 to `1.1e-2`, detectors M1–M9 fire. **Failed:** ADR-0007 §7 item 11 on the tabulated crust (`1.04e-3` vs `1e-3`), diagnosed and quantified — the smooth `dε/dp` authority omits Hartle's internal delta-function shells at the crust's density steps, worth ≈ `4.6 %` of `δM̂` on DS(CMF)-1 (Stieltjes evaluation on the profile's own `ε` steps, corroborated by the TOV-sequence derivative to `7e-5`); **not repaired in 4D**; **Phase 4D-RI (2026-09-03) implemented the ADR-0008 correction**: term 1 of `ODE_HartleMonopole_` is now the measure `−4πr²ξ̂₀·(Δε_i/Δr_i)` of the segment the driver is inside, the driver advances one profile segment per call with the partition validated and fail-closed, the surface shell is the terminal `ε_* → 0` atom of that same measure, and `dε/dp` is kept only for the regular-centre series (ADR-0008 Q8). Constant-density output bitwise unchanged; DS(CMF)-1 `δM̂` `+3.2…+6.5 %`; radial spread `1.6 % → 3.7e-5`; detectors D1–D4 fire. **Historical 4D-RI disposition:** corrected revalidation was outstanding and no baseline existed. **Current 4E:** INV-08 CLOSED/VERIFIED within ordinary-NStar fixed-ε_c l=0 scope; existing response interface ratified; `docs/validation/PHASE4_CLOSEOUT.md:158`. Independently VERIFIED and first baseline established; `docs/validation/PHASE4D_MONOPOLE_BASELINE.md:1` (INV-08) |
| `HartleMonopoleResponse`, `PhysicalHartleMonopole` | **LIVE — MEASURE-COMPLETE SOURCE (ADR-0008); INDEPENDENTLY PHYSICALLY VERIFIED (4D-RV + 4D-DA); FIRST TRUSTED REGRESSION BASELINE ESTABLISHED (4D-BL); NORMALIZED RESPONSE RATIFIED AS PHASE-5 INPUT (4E); PHYSICAL TYPE FOR EXPLICIT SPIN ONLY** | `docs/validation/PHASE4_CLOSEOUT.md:85`, `docs/validation/PHASE4_CLOSEOUT.md:115`. The two second-order result types (`RotationSolver.hpp`), mirroring the accepted first-order architecture. Every canonical field carries the `_over_Omega2` suffix because that is what it is — a coefficient, not a perturbation. `NStar::ComputeHartleMonopoleResponse()` **performs work and says so**; it is deliberately **not** run by star construction, so no existing workflow pays for an integration it never uses (measured: zero solves on construction, one per explicit request, none per materialization). `NStar::MonopoleAt(AngularVelocity)` materializes by scaling — zero spin gives exact zeros, `+Ω` and `−Ω` give bit-identical results. No implicit physical spin (ADR-0006 clarification, extended by ADR-0007 P9). The homogeneous sequence-derivative response is deliberately **not** public in Phase 4C |
| `RotationSolver::Solve(Axis,…)`, `Solve(double,…)`, `ODE`, `GetMass`, `GetPress` | **DECLARED, UNDEFINED** | Declared at `RotationSolver.hpp:370-375,403,315,318`; no definition, no symbol in `libCompactStar.a` (removed by the owner's rework `3639d71`). Sole reference `main/Examples/rotating_ns.cpp:62` is not in `main/Examples/CMakeLists.txt` and would not link |
| `StarBuilder` | **LIVE** | On the file-reading path only |
| `SeqPoint`, `Prog` | **LIVE** | |

### Evolution

**LIVE:** `State`, `StateVector`, `StateLayout`, `StatePacking`, `RHSAccumulator`, `StarContext`,
`GeometryCache`, `DriverContext`, `EvolutionSystem`, `EvolutionConfig`, `GSLIntegrator`,
`IObserver`, `TimeSeriesObserver`, `DiagnosticsObserver`, `DiagnosticCatalog`,
`DiagnosticPacket`, `RunPaths`, `RunBuilder`, `RunObservers`.

| Component | Status | Note |
|---|---|---|
| `ProfileVersionedCache` | **LIVE — one client** | Only `NeutrinoCooling`. Not used by DU mask, heat capacity, or rotochemical |
| `UnitContract` / `UnitVocabulary` | **COMPILED, UNEXERCISED** | Plumbing real; every producer returns an empty contract |
| `RotochemicalCache` | **NOT COMPILED · CANDIDATE** | Absent from `Physics/Evolution/CMakeLists.txt`; `Build()` has zero callers. **Known nonconformant with ADR-0001** — see below |
| `CheckpointObserver` | **EMPTY** | 0 bytes; not in any CMakeLists |

### Drivers

| Driver | Status | Note |
|---|---|---|
| `MagneticDipole` | **LIVE** | `use_moment_of_inertia` is a logged no-op |
| `NeutrinoCooling` | **LIVE** | Divides `L_ν,∞` by the governed `C_⋆(T∞)` — **conformant** with ADR-0002 (`NeutrinoCooling_Details.cpp:889`, `:968`). Emissivity normalizations self-labeled placeholders; `K_PBF = 0.0`. Carries a null-deref ordering defect — see §4 |
| `PhotonCooling` | **LIVE** | Divides `L_γ,∞` by the canonical `C_⋆(T∞)` from `StarContext::HeatCapacityStar_Tinf` — **conformant** with ADR-0002 (Phase 2A-3). `Options::C_eff` removed. Requires `ctx.star` and `ctx.thermo` when actively cooling; a disabled driver still needs neither |
| `Rotochemical` | **NOT COMPILED · CANDIDATE** | `Driver/Chem/CMakeLists.txt` sources list is empty |
| `HeatingFromChem` | **NOT COMPILED** | Header only; no `.cpp`; commented out of CMakeLists |
| `AccretionTorque`, `BNVSpinTorque`, `BNVSource`, `WeakRestoration`, `Coupling` | **EMPTY** | 0–1 byte files; five appear in `install(FILES …)` rules |
| Envelope: `IEnvelope`, `SurfaceGravity` | **LIVE** | `SurfaceGravity` is the most complete file in the layer |
| `EnvelopePotekhin1997` | **LIVE** | Iron fit real; Accreted variant is a self-declared placeholder |
| `EnvelopePotekhin2003` | **LIVE, MISLEADING** | Forwards verbatim to the 1997 fit |

### Microphysics

| Component | Status |
|---|---|
| `Rates/Urca.hpp` | **NOT COMPILED** — mislabeled copy of `RateChannels.hpp`; no includer; `# add_subdirectory(Rates)` commented out |
| BNV channels / analysis | **COMPILED, UNEXERCISED** — no evolution driver connects them |

---

## 3. Ambiguous authority — unresolved

These are live conflicts. Under `GOVERNANCE.md` §3 they are fail-closed until adjudicated.

1. **TOV integration — RESOLVED for the ordinary visible sector (ADR-0005, Phase 3E).**
   **`TOVSolver::SingleStarSolveToTOVPoints` is the CANONICAL NUMERICAL PRIMITIVE** and, as of
   Phase 3E-I4, the **only** ordinary-star radial implementation in the codebase. It owns the
   central-density clamp, the `p_of_e` conversion, the initial conditions, the GSL RK8PD driver,
   the radial grid, the step ladder, the pressure-cutoff termination and `TOVPoint` construction
   — and nothing else.

   | Role | Component |
   |---|---|
   | **CANONICAL NUMERICAL PRIMITIVE** | `TOVSolver::SingleStarSolveToTOVPoints` |
   | **TARGET-MASS ORCHESTRATOR** | `TOVSolver::SolveToProfile` |
   | **SEQUENCE / WORKFLOW ORCHESTRATOR** | `TOVSolver::Solve(Axis, dir, file)` |
   | **RESOLUTION-SWEEP DIAGNOSTIC ORCHESTRATOR** | `TOVSolver::GenTestSequence(ec, …)` |
   | **WORKFLOW OUTPUT CONTRACTS** | `<file>_Sequence.tsv`, `<file>_TestSequence.tsv` |

   > **ADR-0009 source conformed (TOV-SURF-I-R):** ordinary trial pressure below
   > the cutoff is no longer a fatal surface signal. The final accepted crossing
   > is located by `dr/dp=1/f_p`, `dm/dp=f_m/f_p`, reusing the radial RHS;
   > status `SURFACE_REACHED` is required for publication. `SolveToProfile` skips
   > failed coarse samples, rejects failed bisection members and has no nearest
   > fallback. `Solve`, `GenTestSequence` and the NStar wrapper guard completion.
   > No MixedStar event migration. Source: `CompactStar/Core/src/TOVSolver.cpp:1490`,
   > `:2493`, `:2516`, `:2613`; validation: `docs/validation/TOV_SURFACE_IMPLEMENTATION.md:339`.

   **All three orchestrators delegate to the one primitive.** `TOVSolver::RadiusLoop` — the
   duplicate ordinary-star radial loop — was **REMOVED** in Phase 3E-I4 after `GenTestSequence`,
   its last caller, was covered and migrated. Both output contracts are preserved and guarded by
   tests. ADR-0009 subsequently added completion guards without changing these output contracts.

   The primitive path is **VALIDATED** as of Phase 2B-2 against the exact Schwarzschild interior
   solution (to `3.5e-16`) and the official CompOSE `eos.mr` (`M_max` to `2.8e-4`, radii to
   0.20–0.35 %); see `docs/validation/TOV_REFERENCE.md`. That document also records two
   deliberately unrepaired characteristics: the surface is the EOS table floor rather than
   vacuum, and the default `r_max = 70 km` with `radial_res = 10000` leaves ~80 % of the radial
   grid outside the star.

   **`GOVERNANCE.md` fail-closed condition #3 is CLOSED for ordinary visible-sector TOV radial
   numerical ownership.** It does **not** cover `MixedStar`/`RadiusLoopMixed` two-fluid
   integration or the dark-sector TOV paths — those are distinct physics scopes, still
   uncanonicalized, and remain listed below.

2. **Two `NStar` profile-construction blocks** — `BuildFromTOV` and
   `InitFromTOVSolver`+`Append`+`FinalizeSurface`, with duplicated hardcoded column layouts.
   **Still two after 3E-I1/I2/I4**, deliberately: ADR-0005 Q3 = **P3 staged**, so I1 unified the
   radial numerics and I2 unified the *geometry mathematics* — but the two construction styles
   themselves remain. Converging them is **optional increment I3**, a postprocessing question
   that does not reopen the authoritative-path question. **As of Phase 3E-I2 both ordinary visible-sector paths use the canonical
   `CompactStar::Geometry` owner**: `BuildFromTOV` (3D) and now `Append` (Λ) plus
   `FinalizeSurface` (proper volume). Path-1 and Path-2 `B` are consequently **bitwise
   identical**, closing the last ADR-0004 conformance gap on these paths. Path 1 still leaves the
   profile's mirror `M`/`R`/`z_surf` at **zero** (preserved under M1; classified
   `INTERNAL STATE ASYMMETRY — CURRENTLY UNOBSERVED`, not declared correct). The `Analysis` and
   profile-export hooks are preserved, pending classification in Phase 3F.
   Phase 3E-0 measured the two blocks as **value-equivalent** but found two real interface asymmetries:
   `FinalizeSurface` never calls `SetSurfaceScalars`, so Path 1 leaves the profile's mirror
   `M`/`R`/`z_surf` at **zero** while Path 2 populates them; and Path 1 alone owns the
   sequence-accumulation and unconditional `_Sequence.tsv` export, which is the **only** thing
   its six live `main/` callers consume. No caller anywhere attaches an `Analysis` or an export
   condition, so those two hooks are currently dead on all production callers.
3. **Proper volume defined in several places** (INV-04) — **partially resolved in Phase 3D**;
   the validated path now has one owner, the legacy sites do not. See the Phase-3D entry below.

**Resolved since the Phase-0 audit:**

- **Species semantics** are no longer ambiguous. Per **ADR-0001 (ACCEPTED 2026-08-31)**,
  `StarProfile::BaryonDensity` stores `n_B` in fm⁻³ and species columns store dimensionless
  fractions `Y_i = n_i/n_B`, with `n_i = Y_i n_B` derived at the point of use. `TOVSolver` and
  `NStar` preserve the EOS-supplied values without normalization. See INV-01.
- **Heat-capacity ownership** is no longer ambiguous. Per **ADR-0002 (ACCEPTED 2026-08-31)**,
  the thermal degree of freedom has exactly one physical heat capacity — `C_⋆(T∞)`, the
  GR-integrated EOS/CompOSE-based stellar heat capacity, designated to
  `StarContext::HeatCapacityStar_Tinf` — and every energy channel divides by it. See INV-15.
  **The ownership question is decided; the code does not yet conform** — see the next section.

---

### `PhotonCooling` — ADR-0002 conformance (Phase 2A-3)

**Corrected.** Previously the driver divided `L_γ,∞` by a driver-local constant
`Options::C_eff = 1e40`, hand-set at the call site, while `NeutrinoCooling` used the canonical
`C_⋆(T∞)` — two different heat capacities summed into one state element.

Now:

- `PhotonCooling` obtains `C_⋆(T∞)` from `StarContext::HeatCapacityStar_Tinf`, the same call
  `NeutrinoCooling` makes, so the live equation is
  `C_⋆(T∞) dT∞/dt = −L_ν,∞ − L_γ,∞`.
- `PhotonCooling::Options::C_eff` **no longer exists**; it was removed rather than deprecated,
  because ADR-0002 permits a constant heat capacity only as an explicit *whole-balance*
  approximation, which a per-channel option can never be.
- The driver requires `ctx.star` and `ctx.thermo` when actively cooling and **fails closed** with a
  diagnostic if either is absent. A deliberately disabled driver
  (`radiating_fraction <= 0` or `global_scale <= 0`) still returns zero and needs neither.
- A new `C_star_erg_K` diagnostic exposes the denominator so later regressions can audit it.
- **ADR-0002 Pattern A is preserved.** No central thermal-energy manager, power-only driver
  interface, or new RHS ownership was introduced; Pattern B remains the deferred Phase-3 question.

Verified by `tests/thermal/photon_cooling_conformance.cpp`, which asserts
`C_star_erg_K == HeatCapacityStar_Tinf(...)`, `dT∞/dt == −L_γ/C_⋆`, and
`d ln T∞/dt == (dT∞/dt)/T∞`, and which **fails against the old `1e40` denominator** (confirmed by
a controlled temporary regression).

**Consequence for historical outputs:** every passive-cooling product generated before this change
is superseded as a validation reference. They are retained, not deleted, and must not be used as a
regression baseline.

**`main/Test/spin_therm_evol_main.cpp` is no longer compatible** with the governed contract: it
builds an empty `StarContext` and sets no `ctx.thermo`, so its photon cooling now fails closed and
contributes zero. That program is an infrastructure demo that relied on the constant precisely
because it has no thermodynamics; wiring a real EOS into it is outside this correction and was not
attempted. `spin_therm_evol_2_main.cpp` supplies both and is unaffected.

### `RotochemicalCache` — ADR-0001 nonconformance

Recorded for Phase 5. **Not repaired here.**

- Status is unchanged: **NOT COMPILED · CANDIDATE**. It is absent from every CMake source list,
  `Build()` has zero callers, and it has never produced output.
- Under the accepted ADR-0001 semantics, its per-species integration path uses the **wrong
  primitive**: `RotochemicalCache.cpp:147` passes the raw stored `Y_i` column into
  `ComputeEnclosedNumber` (`:25-44`) and `ComputeStructuralDerivative` (`:47-104`), both of which
  document and treat the argument as `n_i` in fm⁻³ (`RotochemicalCache.hpp:116,133`). No `× n_B`
  is applied, so `N_i`, `A_i`, and `B_i` are all computed from fractions where densities are
  required.
- This is **pre-existing unvalidated candidate code** from `675b4a9`, not a regression introduced
  by ADR-0001. The ADR made an existing latent defect explicit.

## 4. Known structural hazards

> **Cache provenance is governed by ADR-0003 (ACCEPTED) and was implemented in Phase 3B.**
> The five hazards that Phase 2B-3 measured are **repaired and now enforced by CTest**; the
> `--audit-known-hazards` mode no longer exists. Historical measurements are preserved, labelled
> superseded, in `docs/validation/CACHE_CORRECTNESS.md`.

**Current cache contract.**

- **Runtime provenance is `(StarProfile identity, StarProfile::Version())`**, carried by the
  typed `ProfileProvenance` in `CompactStar/Physics/Evolution/ProfileProvenance.hpp`. There is no
  generated UUID, no global registry and no persistent/serialized profile ID.
- **Lifetime rule:** a profile-derived context, cache or snapshot **must not outlive the
  `StarProfile` it was built from**. Pointer identity is meaningful for exactly that long. In
  practice `StarContext` is constructed only by callers — `main/Test/spin_therm_evol_2_main.cpp`
  and the test harnesses; no library code constructs one — so the caller already owns and
  outlives both.
- **`GeometryCache` is an immutable provenance-carrying snapshot.** It records
  `(profile, version)` at construction and exposes `Provenance()`, `SourceProfile()`,
  `SourceVersion()` and `Matches(ctx)`. There is deliberately **no `Refresh()`**: a changed
  profile means the caller constructs a new one. Its geometry arrays (`R`, `Mass`, `Area`,
  `ExpNu`, `ExpLambda`, `WV`, `WVExpNu`, `WVExp2Nu`, …) are **numerically unchanged** by that
  work, and remain bit-identical after Phase 3D (below) — only the *source of the Λ formula*
  changed, not the composition.
- **The proper-volume measure has one mathematical owner** — `CompactStar/Geometry.hpp`
  (`CompactStar::Geometry`), per **ADR-0004 (ACCEPTED 2026-09-01)**. It owns `f = 1 − 2m/r`,
  `Λ = −½ ln f`, `e^Λ` and `w_V = 4π r² e^Λ`, **including the domain and failure semantics**. It
  is header-only, depends on the standard library alone, holds no state and adds **no edge to
  the dependency graph** — which is why it sits at the top level beside `Units.hpp` and can be
  included from `Core`, `EOS`, `Physics`, `Extensions` and `Microphysics` alike.

  Ownership is deliberately **split three ways** (ADR-0004 §0-Q1, Option B):
  the primitive owns the *mathematics*; **`GeometryCache` remains the canonical *cached*
  representation** (`ExpLambda`, `WV`, `WVExpNu`, `WVExp2Nu`, with its ADR-0003 provenance and
  its `DataColumn` composition kept verbatim); and **consumers keep their own physics factor and
  unit conversions** — `n_B`, `c_V`, `Q_ν`, and the `1e54` baryon conversion, which belongs to
  INV-14 and is deliberately *not* in the measure. `Core` was **not** made to depend on
  `Physics/Evolution`; `NStar.cpp` includes no `GeometryCache`.

  **Degenerate inputs no longer disagree.** Before 3D the repository held six mutually
  inconsistent behaviors for the same invalid input, differing by `3.16e7` at `f = 0` and
  including a silent divisor substitution inherited from `Zaki::Vector::DataColumn::operator/=`.
  The accepted contract is: exact regular-centre limit at `r = m = 0`; **fail closed** for
  `r < 0`, for `r = 0` with `m ≠ 0`, for `f ≤ 0`, and for non-finite input; no clamp and no
  epsilon anywhere. `SurfaceGravity` and the Hartle/TOV coefficients are **not** this measure
  and are untouched (ADR-0004 §4.4).

  **What actually conformed, and what did not.** Conformed: `NStar::BuildFromTOV` — Λ production
  (bit-identical) and the baryon-number integrand (`|ΔB|/B = 1.368e-16`, one ULP, against
  `1.0e-15` predeclared before implementation) — and `GeometryCache::DeriveLambdaFromMR_`.
  **Still carrying their own inline forms, and not migrated:** TOV **Path 1**
  (`NStar::Append`, `NStar::FinalizeSurface`), pending Phase 3E equivalence coverage; the scalar
  `NStar::BaryonNumIntegrand(double)`, which has a separate INV-14 defect and zero callers and
  was deliberately **not** repaired here; all six `MixedStar` sites, pending focused coverage
  (`MixedStar` is COMPILED but UNEXERCISED, with zero tests); and the `GOVERNANCE.md` §5
  candidate scientific code. **These are governed by ADR-0004 but are not yet conformant.**

- **`StarContext` follows contract S1** — it stays valid across a sanctioned in-place profile
  mutation. On a revision change `RefreshDerivedCachesIfNeeded_` **re-binds the column views
  first**, then invalidates derived payloads, then advances the cached revision **last**. A
  profile that no longer satisfies the required schema makes `BindColumnsOrThrow_` throw; because
  the revision is advanced only on success, the context is never left falsely marked current, and
  a later repair restores it. The revision check deliberately precedes the `IsValid()` guard, and
  the four derived accessors refresh before testing `IsValid()`, so a failed re-bind fails closed
  rather than degrading to silent nulls.
- **`ProfileVersionedCache<T>` keys on `(identity, version)`**, not the version alone. It knows
  profile provenance **only**; consumer-specific dependencies stay typed at the consumer, and
  there is no generic `void*` dependency list or dependency graph.
- **Heat capacity `C_⋆(T∞)` is keyed on `(profile version, thermo identity, geometry
  provenance)`.** A caller-supplied `GeometryCache` whose provenance does not match the context's
  current profile **throws** — silently substituting a locally built geometry would discard what
  the caller asked for. An *equivalent* geometry rebuilt for the same `(profile, version)` is
  interchangeable and does not force a rebuild: the key is provenance, not object address.
  `CompOSE_Thermo` needs no version because it exposes no setter or reload — pointer identity
  suffices.
- **The `NeutrinoCooling` payload depends on `(profile identity, version, geometry
  provenance)`.** Driver `Options` are **not** in the key: `SetOptions` exists, but the payload
  builder reads no option — `include_*` and `global_scale` are applied at evaluation time. A
  `DriverContext` whose `geo` does not belong to its `star` fails closed through the driver's
  ordinary `ok=false` diagnostic, without terminating.

**Out of the ADR-0003 contract, unchanged:** `RotationSolver`'s bracket/`omega_bar` acceleration
caches (algorithm-local, cannot outlive their own solve), `TOVSolver`'s EOS splines and
`gsl_interp_accel` (keyed to the imported EOS table, not to a profile), and the
`TimeSeriesObserver` name→pointer maps (bookkeeping, not scientific state). No provenance hazard
was demonstrated for any of them.

- **No driver dependency graph.** `IDriver::DependsOn()` / `Updates()` are declared, overridden
  by every driver, and **never called**. `EvolutionSystem.cpp:125-134` iterates in registration
  order. The comment at `IDriver.hpp:102-103` referring to an "Evolution graph" is stale.
- **Declared-but-undefined symbols reachable from compiled code** —
  `Physics::Spin::DipoleFieldEstimate` and `CharacteristicAge` are called at `Pulsar.cpp:204,212`
  with no definition anywhere.
- **Null-deref ordering — ✅ RESOLVED (Phase 2A-3).** `NeutrinoCooling_Details.cpp` used to
  dereference `ctx.star` in its `HeatCapacityStar_Tinf` call twelve lines before the
  `if (!ctx.star)` guard, which therefore could never fire. The guard now precedes first use.
  No emissivity, rate, cache, option, or numerical constant changed.
- **Stiff early transient, explicit stepper.** `T∞` falls from `1e9 K` to `1.2e7 K` between 100
  and ~300 yr under the placeholder neutrino normalizations, so the passive system is effectively
  stiff there. The default `RKF45` is adequate on measured evidence, but one output cadence
  (`samples_per_decade = 300`) produces a non-finite state in that window while 50–250 and 400 all
  complete and agree to 4.3e-6. Cause undiagnosed; recorded in
  `docs/validation/PASSIVE_COOLING_BASELINE.md`. A stiff method with a real Jacobian is the likely
  long-term answer, and `MSBDF` remains unavailable until one exists.
- **Heat-capacity cache key — ✅ REPAIRED (Phase 3B, ADR-0003).** The key was
  `(profile version, thermo pointer)` while the supplied `GeometryCache` was a genuine input, so a
  later call at the same profile version with a different geometry silently reused the earlier
  table — measured at **50 %** error in Phase 2B-3. The key now includes the geometry's
  provenance, and a supplied geometry belonging to a different profile or revision **fails
  closed**. See the cache contract above.
- **Canonical-baseline exposure — NONE, and now additionally protected by the contract above.** `passive_cooling_regression` now asserts, on every one
  of 602 driver-context observations in the canonical run, that exactly one profile version, one
  `GeometryCache`, one `StarContext` and one `CompOSE_Thermo` are in play from start to finish.
  Every hazard above needs either a structural mutation during evolution or a cache/driver reused
  across contexts; the canonical run does neither.

---

## 5. Build reality

Re-authenticated at **`11ffe45`** after roadmap Phase 1. Full evidence and commands:
[`docs/build/MACOS_BUILD.md`](../build/MACOS_BUILD.md).

### Configuration and build

- **A clean macOS checkout configures out of source, and the library builds and links.**
  The seven historical/private program directories (`BNV_2022`, `EOS`, `SpinDown_2022`,
  `SpinDown_2023`, `BNV_2023`, `BNV_2024`, `LightDM_2024`) are gitignored, absent from a clean
  checkout, and are now **optional**: each is added only if it carries a `CMakeLists.txt`,
  otherwise a `STATUS` message is emitted (`main/CMakeLists.txt:7-23`). Tracked
  `main/Examples` and `main/Test` remain **required** and unconditional (`:28-29`), so a
  misspelled required directory still fails configure.
- **Python discovery still needs help on the reference Mac.**
  `find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)`
  (`CMakeLists.txt:82`) fails when the first Python CMake finds lacks NumPy. The canonical
  configure command therefore passes `-DPython3_EXECUTABLE`, or a project `.venv` supplies one.
- **Builds are genuinely out of source.** `configure_file` writes the generated
  `CompactStarConfig.h` into the **binary tree** at
  `${CMAKE_BINARY_DIR}/generated/include/CompactStar/Core/CompactStarConfig.h`
  (`CMakeLists.txt:110-116`), and that generated root precedes the source root on the include
  path (`:137`). A fresh configure leaves tracked source state unchanged.
  **Artifact debt remains:** a stale generated copy is still tracked at
  `CompactStar/Core/CompactStarConfig.h`, alongside an editor duplicate
  `CompactStar/Core/CompactStarConfig 2.h`. Both are inert for the build — the binary-tree copy
  wins — and neither was removed; retiring them is a generated-artifact decision.

### Platform support — unchanged

- **macOS is still the only authenticated development platform.** This repository ships Zaki and
  Confind as prebuilt static archives for `Darwin/{arm64,x86_64}` only
  (`dependencies/lib/{Zaki,Confind}/Darwin/…`), and `CMakeLists.txt:95-101` hard-fails when the
  archive for the host platform is absent. Configuration therefore still fails on other platforms.
- **Dependency source exists externally but is not integrated here.** Per project-owner authority,
  ZakiLib (private) and CONFIND (public) exist as external version-controlled repositories. That
  is **not** the same as authenticated binary equivalence or cross-platform support: their current
  source revisions have **not** been shown to reproduce the archives CompactStar consumes, and
  their build contracts have not been authenticated. This repository contains no dependency
  source, no submodule, no `FetchContent`, and no package-manager integration.
  **Cross-platform status is unchanged.**

### Build configuration, warnings, diagnostics

- **A default development build type now exists.** `Debug` is selected when the user supplies
  none, and only for single-configuration generators (`CMakeLists.txt:35-43`); a user-supplied
  `-DCMAKE_BUILD_TYPE` is never overridden and multi-config generators are untouched.
  `Debug` is a **development** default, not a production or scientific configuration, and
  **`Debug` and `Release` are not claimed to be numerically equivalent** — no comparison exists.
- **A warning policy exists for CompactStar-owned code.** `-Wall -Wextra`, applied **`PRIVATE`**
  to the `CompactStar` target for AppleClang/Clang/GNU (`CMakeLists.txt:161-166`), so the flags
  reach neither consumers of the installed library, nor `main/`, nor `tests/`.
- **Dependency headers are `SYSTEM`-classified** (`CMakeLists.txt:144-149`), so GSL and the
  vendored Zaki/Confind headers compile with `-isystem` and their diagnostics are excluded from
  CompactStar's inventory. This suppresses reporting; it does not fix those dependencies.
- **`-Werror` is not enabled.** A clean `Debug` build succeeds **with warnings**; they were
  inventoried during Phase 1D and **none was repaired**. Several classes touch scientific/EOS
  code — hidden overloaded virtuals in the EOS particle hierarchy, VLAs in `SigmaOmegaRho*`, and
  set-but-unused variables in `TOVSolver` — and require classified review before any change.
  See [`docs/build/MACOS_BUILD.md`](../build/MACOS_BUILD.md) for the full table.
- **No sanitizer policy and no assertion policy exist.** Neither has been established. Production
  source still contains no assertions; the checks introduced in Phase 2A–2B live in `tests/`, not
  in the library.

### Automated tests

**Current Phase 5A-2 inventory:** 46 full / 23 self-contained, all green without skips or
exclusions in final serial runs (686.34 s / 91.14 s). The added
`rotochemical_local_thermodynamics` test is a self-contained analytic validation fixture, not a
scientific baseline. All eight governed artifacts are unchanged. Evidence:
`docs/validation/PHASE5A2_LOCAL_THERMODYNAMIC_IMPLEMENTATION.md`. The dated Phase-3 inventory
below is historical.

- **CTest infrastructure exists.** `include(CTest)` at top level provides standard `BUILD_TESTING`
  and `enable_testing()`; `tests/` is added only when testing is on (`CMakeLists.txt:222-226`).
- **`tests/` is the canonical automated-test root**, holding **19 registered CTests** as of the
  Phase-3 closeout (`tests/CMakeLists.txt`). Ten are self-contained; nine are registered only
  when `COMPACTSTAR_EOS_DATA_ROOT` supplies the authenticated CompOSE tables (the guard
  *excludes* them, it does not skip them): `heat_capacity_real_star`,
  `passive_cooling_regression`, `tov_reference_cmf`, `grid_convergence_cmf`,
  `hartle_moment_inertia_cmf`, `baryon_number_cmf`, `tov_path_equivalence_cmf`,
  `tov_sequence_workflow_cmf`, `tov_gen_test_sequence_cmf`. Measured at the closeout HEAD:
  **19/19 pass** with the data root (197.89 s), **10/10 pass** without it (13.44 s). Purpose
  classification and the one assertion-free diagnostic (`heat_capacity_real_star`) are recorded
  in `docs/validation/PHASE3_CLOSEOUT.md` §16.
- **`main/Test/` remains manual demo/debug programs**, not a suite. None of its eight executables
  is registered with CTest.
- **The smoke test is infrastructure validation, not a scientific baseline.** It establishes that
  CompactStar headers compile, that a test executable links a genuine out-of-line symbol from
  `libCompactStar.a` (`CompactStar::Core::Prog::GetName()`, verified present as a defined `T`
  symbol rather than an unresolved import), that transitive dependencies link, and that the binary
  runs and exits zero under CTest. **It asserts no scientific value of any kind.**
- **A scientific test suite now exists**, built across Phase 2A–2B: heat capacity
  (`tests/thermal/heat_capacity_v1.cpp`, `heat_capacity_real_star.cpp`), photon conformance
  (`photon_cooling_conformance.cpp`), the stepper contract
  (`evolution_stepper_contract.cpp`), the passive-cooling regression baseline
  (`passive_cooling_regression.cpp`), TOV references (`tests/core/tov_reference_*.cpp`), cache
  contracts (`tests/evolution/cache_contract.cpp`, `tests/thermal/cache_thermal_contract.cpp`),
  grid convergence (`tests/thermal/grid_convergence_cmf.cpp`) and the scale-free Hartle moment
  of inertia (`tests/rotation/hartle_moment_inertia_*.cpp`). Each has a durable evidence record
  under `docs/validation/`. **Scope limits are recorded in those records and are not superseded
  by the existence of a green suite.**
- **No CI exists.** No `.github/`, `.gitlab-ci.yml`, or equivalent is present.

### Tracked-artifact debt — unchanged

- **Generated artifacts still dominate the tracked tree** — Doxygen output, a committed Xcode
  build directory, and ~400 MiB of committed run results including three ~90 MB profiler traces.
  The **~84%** figure comes from the Phase-0 reconnaissance and is retained with that evidence
  scope; it was not recomputed here. Phase 1 added **3** tracked files to a tree of ~2,650, so it
  cannot have moved the figure materially.

---

## 6. What this document does **not** claim

- It does **not** claim the rotochemical pipeline is operational. Only the independent local
  ADR-0010 thermodynamic provider boundary and analytic free leptons are compiled; no star,
  coefficient, reaction, heating, or evolution connection exists.
- It does **not** claim second-order Hartle is validated. It is publicly callable, has zero
  repository callers, is unverified, and its equations are recorded as defective (INV-08;
  `docs/validation/PHASE4_ROTATION_ENTRY.md` §10–§12). Phase 4A left it byte-identical.
  Phase 4C-G derived the governed replacement from the primary source, **ADR-0007 was ACCEPTED
  on 2026-09-02**, and **Phase 4C-I1 (2026-09-03) replaced the candidate with a conforming
  implementation**. **Phase 4D (2026-09-03) then validated that implementation independently**
  (`docs/validation/PHASE4D_MONOPOLE_VALIDATION.md`): two test-only solvers in different variables,
  analytic and continuum limits, and the published Chandrasekhar–Miller (1974) and Hartle–Thorne
  (1968) second-order tables all agree with production. What 4D did **not** establish is ADR-0007
  §7 item 11 on a tabulated crust with density discontinuities, and the diagnosis quantified the
  omitted internal delta-function shells at ≈ `4.6 %` of `δM̂` on DS(CMF)-1 — so the status is
  **HARTLE MONOPOLE VALIDATION FAILED** (implementation verified; accepted contract incomplete for
  stepped crusts). **ADR-0008 was then ACCEPTED and its correction implemented (Phase 4D-RI,
  2026-09-03):** the EOS source is the measure `−4πr²ξ̂₀ dε`, integrated one profile segment at a
  time. The implementation conforms, but **the corrected independent revalidation has not been
  repeated**, so still: no O(Ω²) number is validated physics, none may be cited, **no monopole
  baseline exists**, and the `l = 2` sector is out of scope rather than verified. **Phase 4D-RV
  (2026-09-03) then executed that corrected independent revalidation on the migrated ADR-0009
  backgrounds** (`docs/validation/PHASE4D_CORRECTED_MONOPOLE_REVALIDATION.md`): every physics line passes with an independent Stieltjes-measure oracle
  (`≤ 5.3e-5` on DS(CMF)-1, sequence identity `≤ 1.03e-4`, published tables reproduced, M1–M10
  fire), but ADR-0008 Validation D's monotonicity clause is NOT MET at the node-placement floor —
  status **CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE** at 4D-RV; **Phase 4D-DA (2026-09-04)
  adjudicated that clause as a heuristic the accepted representation cannot satisfy and clarified it into
  D′1–D′5, all met — status `HARTLE O(OMEGA^2) MONOPOLE RESPONSE VERIFIED`; the first monopole baseline is the next separate task.**
- It **does** now claim the physically normalized first-order response is verified as physics
  (Phase 4B): the shape agrees node-by-node with an independently derived and independently
  normalized profile, satisfies the exterior and volume identities, and reproduces derived
  weak-field coefficients. It does **not** claim anything about the adequacy of truncating the
  structure at Hartle order — `Ω/Ω_K` reaches ~0.6 at 716 Hz, and the size of the neglected
  O(Ω²) terms is an open O(Ω²) question, not a first-order one.
- It claims of the O(Ω) solver **only** that its scale-free observable `I = J/Ω` is validated
  (Phase 2B-4B: equation match against published Hartle, analytic and numerical cancellation
  of the arbitrary normalization, agreement with an independent solver to 9.5e-9 analytic /
  2.1e-5 on the CMF sequence, and the correct Newtonian limit). It does **not** claim the
  absolute first-order normalization is correct: `ω̄`, `Ω`, `J`, the `init_omega_bar` seed and
  the `Ω [s⁻¹]` annotation remain **unresolved** (INV-07).
- It does **not** claim placeholder emissivities represent real microphysics.
- It does **not** claim the passive-cooling regression validates the physics. A regression
  baseline now exists (`tests/baselines/passive_cooling_cmf_1p6_debug.tsv`, CTest
  `passive_cooling_regression`, requiring the authenticated external CMF data), and
  `GOVERNANCE.md` §3.1 condition 7 is satisfied — but the neutrino emissivity normalizations it
  freezes are self-labelled placeholders, so it detects change, not correctness.
- It does **not** claim CI exists. The regression must be run by hand.
- It does **not** claim the corrected cooling *trajectory* has been validated. Only the denominator
  identity and the verified `C_⋆(T∞)` are established; the neutrino emissivity normalizations
  remain self-labelled placeholders.
- It does **not** claim the Phase-1B `RotationSolver` repair validated any physics. The library
  compiles; INV-07 and INV-08 are exactly as unresolved as before.
- It does **not** claim a scientific test suite exists. One infrastructure smoke test exists, and
  it asserts no scientific value. Scientific validation is roadmap Phase 2A/2B.
- It does **not** claim `Debug` is a production or scientifically validated configuration, nor
  that `Debug` and `Release` agree numerically.
- It does **not** claim the external ZakiLib/CONFIND repositories constitute cross-platform
  support. Their revisions have not been shown to reproduce the archives consumed here.
