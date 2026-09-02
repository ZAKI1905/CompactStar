# Phase 3C — Boltzmann-constant authority and ownership audit

> **STATUS: `PHASE-3C IMPLEMENTED — AUTHORITY ADOPTED`.**
>
> The **ownership** question was settled by owner intent plus objective authority: **ZakiLib
> already defines Boltzmann's constant, correctly, and CompactStar already links the symbol.**
> No ZakiLib change and no dependency-archive rebuild were needed.
>
> The owner adjudicated the open tolerance question in the Phase 3C brief and **directed
> adoption with re-baselining** of the `k_B`-dependent artifacts, with the predeclared
> `1.7e-11` bound left **unwidened** and the shift **not** characterized as a regression.
> §13–§15 record that adjudication and the implementation evidence. The `1.3e-6` trajectory
> figure is now *explained*, not merely observed: an identical-magnitude shift is produced by
> changing the integrator tolerance alone with `k_B` held fixed (§14).

| Field | Value |
|---|---|
| **CompactStar HEAD** | `9a0da377f5bb378ff729d0ddba4b4dae7a53f49e` (5 ahead / 0 behind `master` `d2040d89…`) |
| **CompactStar baseline** | 13/13 authenticated, 8/8 self-contained (Phase 3B, unchanged by this audit) |
| **ZakiLib inspected** | `/Users/keeper/Documents/CompactStar/external/Zaki`, branch `master`, HEAD `b9ddebaded24962468954846f47238aec2726fd4`, **clean**, remote `git@github.com:ZAKI1905/ZakiLib.git`, `b9ddeba Initial import of Zaki library` |
| **Change class of this audit** | documentation only |

---

## 1. CompactStar `k_B` occurrence map

**Production — four declarations, two values, two distinct roles.**

| Location | Value | Role | Live? | Validated by | Numerical consequence |
|---|---|---|---|---|---|
| `Physics/Driver/Thermal/src/PhotonCooling_Details.cpp:347` | `8.617333262145e-11` (**LONG**) | `MEV_PER_K`; `T_inf[MeV] = T_inf[K] × MEV_PER_K` (`:348`) | yes | `passive_cooling_regression`, `photon_cooling_conformance` | per-step evaluation |
| `Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp:215` | `8.617333262145e-11` (**LONG**) | `MEV_PER_K`; same conversion (`:881`) | yes | `passive_cooling_regression` | per-step evaluation |
| `EOS/src/CompOSE_Thermo.cpp:567` | `8.617333262e-11` (**SHORT**) | `kB_MeV_per_K` in `CvDensity_cgs`; `c_V,nat × k_B × MeV fm⁻³→erg cm⁻³` (`:571`) | yes | `heat_capacity_v1`, `heat_capacity_real_star` | enters the **cached** `C_⋆` table |
| `EOS/src/CompOSE_Thermo.cpp:721` | `8.617333262e-11` (**SHORT**) | same, in `CvDensity_cgs_ForCooling` (`:725`) | yes | as above | enters the cached `C_⋆` table |

**Tests — four literals, all oracles, all out of scope for change:**
`photon_cooling_conformance.cpp:58` (LONG, "matches the drivers"), `heat_capacity_v1.cpp:61`
(SHORT — and `:57-58` explicitly documents the split as INV-02's recorded issue),
`heat_capacity_real_star.cpp:56` (SHORT). **`main/` contains none.**

The two roles are physically the same constant: one converts a temperature, the other converts a
heat-capacity density. INV-02 already records this split by file and line.

## 2. ZakiLib already owns Boltzmann's constant

| Quantity | Zaki symbol | Value | Unit | Location | Used by CompactStar today? |
|---|---|---|---|---|---|
| Boltzmann | `Zaki::Physics::K_BOLTZ_SI` | `1.380649e-23` | J/K | decl `Constants.hpp:150`, def `Physics/src/Constants.cpp:114` | no |
| Boltzmann | `Zaki::Physics::K_BOLTZ_EV` | `8.61733326214518e-5` | **eV/K** | decl `Constants.hpp:151`, def `Constants.cpp:115` | no |
| Kelvin→GeV | `Zaki::Physics::KEL_2_GEV` | `K_BOLTZ_EV*1e-9` | GeV/K | `Constants.cpp:56` | no |

`K_BOLTZ_SI = 1.380649e-23` **is the SI-exact defining value**. `K_BOLTZ_EV` is the exact
quotient rounded to 15 significant figures — **3.0e-16** from the correctly-rounded double.

**CompactStar already consumes this exact header and archive**: 90 uses of
`Zaki::Physics::MEV_2_INV_FM`, 36 of `SUN_M_KM`, 34 of `INV_FM4_2_INV_KM2`, 23 of `YR_2_SEC`,
14 of `Q_E`, 12 of `LIGHT_C_KM_S`, and more. Adopting `K_BOLTZ_EV` therefore adds **no new
dependency edge and no new component boundary**.

## 3. Independent SI derivation

Both inputs are exact by definition under the post-2019 SI:

```
k_B = 1.380649e-23 J/K          (exact)
1 eV = 1.602176634e-19 J        (exact)

k_B[eV/K]  = 8.61733326214517743366365933408063920…e-5
k_B[MeV/K] = 8.61733326214517743366365933408063920…e-11
```

Correctly rounded to `double`: `8.61733326214517686e-11`, bits `0x3dd7afe8c1362e29`.
No external table is needed — the value follows from the definitions.

## 4. Bit-level comparison — nothing is bit-identical to what exists

| Candidate | double | bits | ULP vs exact | rel vs exact |
|---|---|---|---|---|
| **exact SI reference** | `8.61733326214517686e-11` | `0x3dd7afe8c1362e29` | 0 | 0 |
| `Zaki::K_BOLTZ_EV * 1e-6` | `8.61733326214517944e-11` | `0x3dd7afe8c1362e2b` | **+2** | `3.0e-16` |
| `Zaki::K_BOLTZ_EV / 1e6` | `8.61733326214518074e-11` | `0x3dd7afe8c1362e2c` | +3 | `4.5e-16` |
| CompactStar **LONG** | `8.61733326214499979e-11` | `0x3dd7afe8c1362da0` | **−137** | `2.06e-14` |
| CompactStar **SHORT** | `8.61733326200000020e-11` | `0x3dd7afe8c1347764` | **−112 325** | `1.69e-11` |
| GSL-derived | `8.61734200025788339e-11` | `0x3dd7afea542f9b53` | **+6.76e9** | `1.01e-06` |

`LONG` is `K_BOLTZ_EV` truncated to 13 significant figures. `SHORT` is a 10-figure rounding.
**LONG vs SHORT: `1.68e-11` relative.**

**No candidate reproduces either current literal bit-for-bit.** Phase 3C therefore cannot be
bit-identical engineering under any option.

## 5. GSL is the weakest authority (Model D rejected)

`GSL_CONST_MKSA_BOLTZMANN = 1.3806504e-23` (`gsl_const_mksa.h:38`) is a **pre-2019 CODATA-era**
value, **`1.01e-6` away from the SI-exact definition** — five orders of magnitude worse than
either CompactStar literal and ten orders worse than ZakiLib's. GSL is also inconsistent with
its own `Q_E` for this purpose. The owner's instruction not to assume GSL supersedes ZakiLib is
independently confirmed by measurement: here GSL is simply obsolete.

## 6. Measured numerical impact — the analytic estimate was badly wrong

Analytic sensitivity (`k_B` enters multiplicatively) predicts:

| quantity | analytic |
|---|---|
| `CvDensity_cgs` ∝ `k_B` (SHORT) | `1.68e-11` |
| `C_⋆` ∝ `CvDensity` | `1.68e-11` |
| `T_inf[MeV]` ∝ `k_B` (LONG) | `2.08e-14` |
| `L_ν` — `K_DU`/`K_MU` contain no `k_B` | **0, independent** |

A controlled mutation pointing all four production sites at `Zaki::K_BOLTZ_EV * 1e-6`, built and
re-emitted, then **fully reverted** (production is byte-identical to `HEAD`), measures the
trajectory instead:

| column | max &#124;rel&#124; change |
|---|---|
| `T_inf` | `1.68e-7` |
| `C_⋆` | `1.72e-7` |
| `L_ν` | `1.01e-6` |
| `L_ν,DU` | `1.01e-6` |
| `L_ν,MU` | **`1.34e-6`** |
| `L_γ` | `4.06e-7` |
| `T_surf` | `1.02e-7` |
| `T_b` | `1.68e-7` |
| `d ln T_inf/dt` | `6.65e-7` |

**The observed effect is ~10⁵ times the input perturbation.** The cause is not error growth in
the physics: the RKF45 driver runs at `rtol = 1e-6`, so a `1.7e-11` change in `C_⋆` perturbs
adaptive step selection and the output lands at the integrator's own tolerance scale. `L_ν` is
analytically independent of `k_B` yet moves `1.01e-6` — which is only explicable as step-choice
sensitivity, and confirms the mechanism.

**Consequences measured, not assumed:**

- `passive_cooling_regression --compare` still **PASSES** (`1.3e-6` ≪ state `1e-5`, luminosity `1e-4`).
- The golden TSV is written at `%.12e`, so it is **NOT byte-identical** — the goldens would have to be re-emitted.
- `heat_capacity_v1` **FAILS** (U1.a, U1.b): its oracle hard-codes the SHORT literal, so it must be updated in lockstep.
- `photon_cooling_conformance` **passes** — its `kMEV_PER_K` only feeds a call, and its assertion is a conformance equality, not a value check.

## 7. Dependency and build path

- Headers resolve from the **vendored** `dependencies/include/` (`CMakeLists.txt:90,149`), **not**
  from `external/Zaki`. The vendored `Zaki/Physics/Constants.hpp` is **byte-identical** to the
  source-repo header.
- The archive is `dependencies/lib/Zaki/Darwin/arm64/libZaki.a` (`CMakeLists.txt:92`).
- `K_BOLTZ_SI` and `K_BOLTZ_EV` are `extern const double` **defined in `Constants.cpp`** — not
  header-only — so they live in the archive. `nm` confirms both:
  `__ZN4Zaki7Physics10K_BOLTZ_EVE`, `__ZN4Zaki7Physics10K_BOLTZ_SIE`.
- The **compiled values were authenticated** by locating their IEEE-754 byte patterns inside the
  extracted `Constants.cpp.o`: `1.380649e-23` → `0x3b30b0e6d55e647c` **found**;
  `8.61733326214518e-5` → `0x3f1696fe94e2cfa0` **found**. The archive matches the source exactly;
  there is no stale-archive risk for these symbols.

**Therefore: no ZakiLib source change, no archive rebuild, and no dependency-version update are
required.** The symbol is already present, correct, and linkable. (Had ZakiLib needed a change,
it would have been a two-repository action — a correct ZakiLib commit is *not* automatically an
authorized CompactStar dependency update.)

## 8. Ownership alternatives

| Model | Assessment |
|---|---|
| **A** — ZakiLib owns SI `k_B`; CompactStar derives MeV/K from a generic eV/J conversion | Viable, but redundant: ZakiLib **already** publishes the eV/K form, so re-deriving it in CompactStar adds arithmetic and ULP noise for nothing. |
| **B** — CompactStar consumes `Zaki::Physics::K_BOLTZ_EV` and scales to MeV/K | **RECOMMENDED.** eV/K is a general particle/nuclear representation, squarely inside the owner's ZakiLib boundary; the symbol exists, is authoritative to `3e-16`, and is already linked. |
| **C** — CompactStar owns a named literal documented as SI-derived | Defeats the owner's stated reusable-library boundary and re-creates the duplication INV-02 records. Rejected. |
| **D** — use GSL | Rejected on evidence: GSL's value is pre-2019 and `1.01e-6` from exact (§5). |

## 9. Classification of `k_B [MeV/K]`

**`ZAKILIB-GENERIC`.** Boltzmann's constant is a fundamental constant of thermal physics, not a
compact-object convention. The owner's boundary explicitly names it. The MeV/K *scaling* is a
trivial decimal shift of ZakiLib's eV/K form, not a domain conversion — unlike, say, MeV fm⁻³ →
erg cm⁻³, which is genuinely nuclear-astrophysics-specific and correctly lives in
`CompactStar/Units.hpp` (Phase 3A).

## 10. Change class and ADR requirement

**Change class: NUMERICAL-METHOD / CONSTANT-AUTHORITY CHANGE.** Not engineering — no option is
bit-identical (§4), and the validated goldens move (§6).

**Not structural**: CompactStar already depends on this header and archive and on ten-plus
sibling symbols; no ownership boundary or component boundary moves.

**ADR: not required.** `GOVERNANCE.md:51` requires an ADR for *structural* changes; the
numerical-method row requires "convergence evidence · regression against baseline ·
`NUMERICAL_POLICY.md` citation once that document exists; until then the numerical rationale is
recorded in the change itself" — which this record provides. No ZakiLib governance record is
needed either, since ZakiLib is not modified.

## 11. Proposed implementation plan (NOT executed)

1. Replace all four production declarations with `Zaki::Physics::K_BOLTZ_EV * 1e-6`, keeping the
   local names (`MEV_PER_K`, `kB_MeV_per_K`) so expression structure is untouched. Use `* 1e-6`
   rather than `/ 1e6`: measured at **+2 ULP** vs **+3 ULP** from the exact reference (§4).
2. Update the three test oracles to the same authority, with a comment stating they now assert
   against ZakiLib rather than an independent literal — and record the loss of independence
   explicitly, since `heat_capacity_v1` exists partly to catch a wrong `k_B`.
3. **Re-emit all affected goldens** and record old→new hashes with this document as rationale.
4. Re-run 13/13 and 8/8; confirm every change is bounded by the measured `1.3e-6`.
5. Update INV-02 to record that the dual-precision `k_B` split is resolved.

Steps 1–2 must land together; step 3 is the governance-sensitive one.

## 12. The question for the owner

Everything above is settled by evidence except one thing, and it is not implementation trivia.

**The Phase-3 plan predeclared `≤1.7e-11` as 3C's preservation tolerance**
(`PHASE3_CONSOLIDATION_PLAN.md`, derived from the constants before any run). The measured
trajectory impact is **`1.3e-6`**, about **10⁵ times larger**, because of adaptive-stepper
amplification the predeclaration did not anticipate. The rule that a tolerance may not be chosen
after seeing the output means this cannot be waved through by widening it now.

> **Q1.** Adopt the authoritative `k_B` and **re-baseline** the five golden artifacts — accepting
> a one-time, physically negligible shift of up to `1.3e-6` (well inside the existing `1e-5` /
> `1e-4` regression tolerances, and *toward* the exact SI value) — recording the old and new
> hashes and this audit as justification?
>
> **Or** defer 3C, leaving the two-precision split in place, on the grounds that the Phase-2B
> goldens should not move for a `1.7e-11` correctness improvement?
>
> A third option exists and should be named: adopt the authority **only** in `CompOSE_Thermo`
> (the SHORT, less accurate site) and leave the LONG driver literal, which is already within
> `2.1e-14` of exact. That halves the diff but leaves INV-02's split formally unresolved.

**Recommendation: the first option.** The change moves every consumer *toward* the SI-exact
value, eliminates a documented invariant defect, costs nothing in dependency terms, and the
resulting shift is four orders of magnitude inside the tolerances that Phase 2B derived from
measurement. But re-baselining validated scientific evidence is the owner's decision to take
knowingly, not a side effect for me to absorb.

---

## 13. Owner adjudication (Phase 3C brief)

The owner selected the first option of §12: **adopt the authoritative ZakiLib constant and
re-baseline only the artifacts that actually depend on it.** The brief attached four binding
constraints, all of which are honoured here:

1. The predeclared `1.7e-11` bound is **not** widened retroactively. It stands in
   `PHASE3_CONSOLIDATION_PLAN.md` as written, annotated with the reason it was not met.
2. The re-baselined goldens are **not** to be characterized as a regression. They are a
   one-time authority correction moving every consumer *toward* SI-exact.
3. No test may be deleted to accommodate re-baselining. None was; all 13 remain.
4. **INV-02 keeps its `VERIFIED CURRENT BEHAVIOR` headline.** The owner explicitly corrected
   the audit's own recommendation here: INV-02 records verified behaviour, not an unresolved
   defect, so only its *note* was updated. This document's §11 recommendation to "update INV-02
   to record that the split is resolved" was therefore executed as a note-only edit.

## 14. Why `1.3e-6` is step selection, not physics

The audit (§6) measured `1.3e-6` on the passive-cooling trajectory and attributed it to
adaptive-stepper amplification. That attribution was an inference. It has now been **tested
directly**, by holding the authority fixed and changing only the integrator tolerance:

| Comparison | `k_B` authority | `rtol` | max &#124;rel&#124; |
|---|---|---|---|
| old vs new | **changed** | `1e-6` (default) | `1.343e-6` |
| old vs new | **changed** | `1e-8` | `3.517e-7` (**3.8×** smaller) |
| old vs old | **identical** | `1e-6` → `1e-8` | **`1.146e-6`** |
| new vs new | **identical** | `1e-6` → `1e-8` | `1.543e-7` |

The third row is decisive. Changing *only* the step-size tolerance, with the constant held
fixed, moves the trajectory by essentially the same amount as changing the constant does. The
`1.3e-6` is therefore a property of where the adaptive driver places its steps, not a physical
displacement — and it shrinks when the integrator is tightened, which a physical shift would
not do.

The **static** diagnostics confirm the complementary half. With no time integration at all, the
fixed-state response is the analytic ratio and nothing else:

| Quantity (fixed state) | measured | analytic |
|---|---|---|
| `CvDensity_cgs` | `+1.684741e-11` | `+1.684735e-11` |
| `CvDensity_cgs_ForCooling` | `+1.684736e-11` | `+1.684735e-11` |
| `C_⋆(10^8 K)` | `+1.686970e-11` | `≈1.685e-11` |
| `L_ν`, `L_ν,DU`, `L_ν,MU`, `L_γ` | **exactly 0** | `k_B`-independent at fixed `T` |
| `M`, `R` | **exactly 0** | `k_B`-independent |

`grid_convergence_cmf_1p6_debug.tsv` is a static (non-evolved) artifact and shows the same
thing at file level: of its sixteen columns, **fourteen are byte-identical** — including
`ec_gcm3`, `achieved_M`, `R_km`, `B`, `z_surf`, `compactness`, `Lnu_1e8` and `Lgamma_1e8` — and
only the two heat-capacity-derived columns move:

| Column | max &#124;rel&#124; |
|---|---|
| `Cstar_1e8` | `1.6983e-11` |
| `dlnT_dt_1e8` | `1.7136e-11` |

**Both are inside the predeclared `1.7e-11` bound.** The bound was correct for the physics; it
was applied to an evolved trajectory whose step placement is not a physical observable.

## 15. Implementation and evidence

**Production sites changed — four, all to the same authority:**

| File | Before | After |
|---|---|---|
| `EOS/src/CompOSE_Thermo.cpp` (`CvDensity_cgs`) | local `8.617333262e-11` | `Zaki::Physics::K_BOLTZ_EV * 1.0e-6` |
| `EOS/src/CompOSE_Thermo.cpp` (`CvDensity_cgs_ForCooling`) | local `8.617333262e-11` | `Zaki::Physics::K_BOLTZ_EV * 1.0e-6` |
| `Physics/Driver/Thermal/src/PhotonCooling_Details.cpp:347` | local `8.617333262145e-11` | `Zaki::Physics::K_BOLTZ_EV * 1.0e-6` |
| `Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp:215` | local `8.617333262145e-11` | `Zaki::Physics::K_BOLTZ_EV * 1.0e-6` |

No `k_B` literal and no GSL Boltzmann constant remain on any production path.
`CompactStar/Units.hpp` deliberately does **not** own `k_B` — ZakiLib does — and its note
records that.

**Test oracles are independent of the production authority.** The three thermal conformance
tests derive `k_B [MeV/K]` themselves in `long double` from the two SI-exact defining
constants, and never import `Zaki::Physics::K_BOLTZ_EV`:

```cpp
static constexpr long double kKbJ_per_K_L = 1.380649e-23L;   // SI exact
static constexpr long double kMeV_in_J_L  = 1.602176634e-13L; // SI exact
constexpr double kB_MeV_per_K_SI = static_cast<double>(kKbJ_per_K_L / kMeV_in_J_L);
static constexpr double kKbTolUlp = 8.0;  // fixed from representation, not from output
```

Measured production-vs-oracle agreement: **1 ULP (`1.5e-16`)**, inside the 8-ULP allowance,
which was fixed from the floating-point representation before the comparison was run.

**Artifact dependency classification (all five goldens re-emitted and compared):**

| Golden | `k_B`? | Result |
|---|---|---|
| `tov_dscmf1_reference.tsv` | no | **byte-identical** |
| `hartle_I_dscmf1_debug.tsv` | no | **byte-identical** |
| `passive_cooling_cmf_1p6_debug.tsv` | yes | re-baselined |
| `grid_convergence_cmf_1p6_debug.tsv` | yes | re-baselined |
| `grid_convergence_cmf_1p6_trajectory.tsv` | yes | re-baselined |

The two `k_B`-independent artifacts remaining byte-identical is a **falsifiable prediction that
held**: TOV structure and the Hartle moment of inertia contain no thermal constant, so a
genuine `k_B`-only change must leave them untouched. Had either moved, the change would have
been touching something it should not.

**Re-baseline provenance:**

| Golden | old SHA-256 | new SHA-256 |
|---|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `edaa5e3bc5984e2d8cf5acee6664816083ecf83415a9355582614cc3baf42d6e` | `831744b0a206541fd0e24adc67876cc1ee4d02d89a580942a9fb0c6749999453` |
| `grid_convergence_cmf_1p6_debug.tsv` | `3be2005f8cfdae2798637e4d51674461c9f56dc36ea48d79ad9459109dcc3c88` | `61d84ddcb87645197c5406c880b648fdf3bb9b0ed8c58350800ca2f2d296ff40` |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `e04f748536d27331bbd383ce4aa11547d5a4f12f927b5df9192d36522a986194` | `ca32863dabaa28fad63d5c36b287a3b94e9b6b85f11980bf2be4e65499d9a0c6` |

Changed columns, `passive_cooling_cmf_1p6_debug.tsv`: `Tinf_K` `1.679e-7`, `C_star_erg_K`
`1.720e-7`, `L_nu_inf_erg_s` `1.007e-6`, `L_nu_DU_inf_erg_s` `1.007e-6`, `L_nu_MU_inf_erg_s`
`1.343e-6`, `L_gamma_inf_erg_s` `4.063e-7`, `Tsurf_K` `1.016e-7`, `Tb_K` `1.679e-7`,
`dLnTinf_dt_1_s` `6.651e-7`, `Lnu_over_Lgamma` `6.010e-7`. That `L_ν` moves at all — it is
analytically independent of `k_B` at fixed temperature, and moved **exactly zero** in the
fixed-state test — is itself the signature of step-placement, not of physics.

Trajectory artifact: `Tinf_K` `1.679e-7`, `L_nu_inf` `1.007e-6`, `L_gamma_inf` `4.063e-7`,
`C_star` `1.720e-7`.

All changes are one to two orders of magnitude inside the `1e-5`/`1e-4` regression tolerances
that Phase 2B derived from measurement, and those tolerances were **not** altered.

## 16. Detector sensitivity is retained (the re-baseline did not blind the test)

Re-baselining a golden always raises the question of whether the test still detects anything.
It does. With the new baselines in place, an artificial `+1.0e-4` relative perturbation was
injected into `k_B` at **all four** production sites, the library rebuilt, and the regression
re-run:

```
The following tests FAILED:
    10 - passive_cooling_regression (Failed)
```

The perturbation was then reverted and the revert verified **byte-identical by SHA-256** against
the pre-perturbation copies of all three source files, the library rebuilt, and the test
returned to `Passed`. So `passive_cooling_regression` still discriminates a `1e-4` thermal
perturbation while accepting the `1.3e-6` step-placement shift — the detector band sits between
the two, which is where it belongs.

## 17. Validation summary

| Check | Result |
|---|---|
| Full authenticated suite | **13/13 passed** (196.68 s) |
| Self-contained suite | **8/8 passed** (13.39 s) |
| Production `k_B` literals remaining | **0** |
| GSL Boltzmann on a production path | **0** |
| ZakiLib modified / archive rebuilt | **no / no** |
| Production-vs-oracle agreement | **1 ULP** (allowance 8 ULP, predeclared) |
| `k_B`-independent goldens | **byte-identical** (TOV, Hartle-I) |
| `k_B`-dependent goldens | 3 re-baselined, hashes recorded (§15) |
| Tests deleted or tolerances widened | **none** |
| Predeclared `1.7e-11` bound | **left unwidened**; satisfied by every static observable |
| Detector sensitivity after re-baseline | **retained** (§16) |
| INV-02 headline | **unchanged** (`VERIFIED CURRENT BEHAVIOR`); note only |
