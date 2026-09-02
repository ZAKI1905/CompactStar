# Phase 3C — Boltzmann-constant authority and ownership audit

> **STATUS: `PHASE-3C OWNER ADJUDICATION REQUIRED`.**
>
> The **ownership** question is settled by owner intent plus objective authority: **ZakiLib
> already defines Boltzmann's constant, correctly, and CompactStar already links the symbol.**
> No ZakiLib change and no dependency-archive rebuild are needed.
>
> What requires adjudication is a **numerical consequence that the Phase-3 plan's own
> predeclared tolerance did not survive.** The plan predeclared `≤1.7e-11` for increment 3C,
> derived from the two constants. The measured effect on the validated passive-cooling
> trajectory is **up to 1.3e-6** — roughly **10⁵ times larger** — because the adaptive ODE
> stepper amplifies the perturbation to its own `rtol = 1e-6` scale. The governing rule is that
> a tolerance may not be widened after observing the output, so this is the owner's call, not
> mine. See §9.

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
