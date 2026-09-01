# ADR-0002 §V1 — Verification of the canonical stellar heat capacity `C_⋆(T∞)`

> **STATUS: V1 BLOCKED.**
> Tier-A analytic verification of the production algorithm is complete and passes.
> **Tier B — the real canonical-star magnitude check — could not be run: no authenticated
> CompOSE EOS assets exist on this machine.** ADR-0002 §V1 item 3 is therefore NOT satisfied, and
> a second, independent blocker was found: the thermodynamic *content* of the CompOSE entropy
> cannot be established from the code alone. **Do not proceed to the `PhotonCooling` conformance
> change on this evidence.**

| Field | Value |
|---|---|
| **Governing authority** | `docs/adr/ADR-0002-thermal-heat-capacity-ownership.md` §V1; INV-10, INV-12, INV-13, INV-15 |
| **Code commit** | `e2cc9270f5bdae677e659233398c360cfbb318a8` |
| **Build configuration** | `Debug` (Phase-1D default) |
| **Toolchain** | AppleClang 17.0.0.17000604, CMake 4.2.1, macOS 15.6.1 (24G90), arm64 |
| **Verification program** | `tests/thermal/heat_capacity_v1.cpp` (CTest: `heat_capacity_v1`) |
| **Production code changed** | **None.** |

## Accuracy criterion — fixed before results were seen

**Relative error ≤ 1 %** for the `T∞`-cache interpolation and for cache-vs-direct agreement.

Justification, stated in advance: `C_⋆(T∞)` enters the thermal RHS as a denominator alongside
neutrino emissivities whose normalisations are self-labelled placeholders
(`NeutrinoCooling_Details.cpp:100-102`, `Q0_DU = 1.0e27`, `Q0_MU = 1.0e21`, linear in ρ) and
alongside composition from an EOS whose own uncertainty is far larger. A sub-percent numerical
error in the denominator is not the limiting term in any Phase-2A use. This is a *numerical
adequacy* criterion for the interpolation scheme only; it is **not** a claim about physical
accuracy.

## 1. Production path, re-authenticated at `e2cc927`

```
StarContext::HeatCapacityStar_Tinf      StarContext.cpp:704-739
  -> StarContext::BuildHeatCapacityCache_   StarContext.cpp:744-820
       -> CompOSE_Thermo::CvDensity_cgs_ForCooling   CompOSE_Thermo.cpp:718-732
            -> CvDensity_Natural_ForCooling  (= T · n_B · dQ2/dT)      :711-715
                 -> dQ2dT_ForCooling         (low-T fit blended to table)  :673-707
                      -> LowT_Slope_dQ2dT0_  (least squares through origin) :625-669
                      -> dQ2dT               (table finite differences)
```

| Property | Authenticated behavior | Evidence |
|---|---|---|
| `T∞` cache range | `[1e-5, 1] MeV` ≈ `[1.16e5, 1.16e10] K` | `StarContext.cpp:775-776` |
| `NT` | `160`, log-spaced | `:777`, `:785` |
| Interpolation | linear in `C` against **log T** | `:730-738` |
| Radial quadrature | trapezoid over the profile grid | `:794-814` |
| Proper volume | `WV = 4π r² e^{Λ}` from `GeometryCache` | `:758`; `GeometryCache.cpp:230` |
| Tolman factor | `T_local = T∞ · e^{−ν}` via `ExpMinusNu` | `:804-805`; `GeometryCache.cpp:206` |
| `KM3_TO_CM3` | `1.0e15` | `:761`, applied `:810-811` |
| Cache key | `(profile version, thermo pointer identity)` — **geometry excluded** | `:712-714` |
| Endpoint clamping | returns endpoint value, silently | `:725-728` |
| Low-T policy | `dQ2/dT` = constant `a0` below the blend window | `CompOSE_Thermo.cpp:700-706` |
| Units | `erg cm⁻³ K⁻¹` × km³ × `1e15` → **erg K⁻¹** | verified numerically, §3 below |

## 2. CompOSE thermodynamic content — **INTERPRETATION BLOCKER**

`Q2` is read as CompOSE quantity 2, entropy per baryon `s/n_B`
(`CompOSE_Thermo.cpp:326`, `:354`), and `c_V = T n_B (∂Q2/∂T)` follows from `s = n_B Q2`.
The algebra is correct. **What the entropy physically contains is not established.**

Two findings:

- **The lepton flag is parsed and then never used.** The `eos.thermo` header is read as
  `mn mp Il` into `m_Il`, documented "Lepton flag (1 if leptons included)"
  (`CompOSE_Thermo.hpp:208`; `CompOSE_Thermo.cpp:291`). **`m_Il` appears nowhere else in the
  class.** Nothing verifies, branches on, or reports whether the loaded table's entropy includes
  leptons. `CvDensity_cgs_ForCooling` returns the heat capacity *of whatever the table contains*.
  ADR-0002 requires the **total** heat capacity of the evolved degree of freedom; if the table is
  baryon-only, the electron and muon contributions are silently missing.
- **Crust/lattice contributions are not represented at all** in this path — the integral is
  `c_V(T, n_B, Y_q)` from a single bulk table over the whole radial range.

Neither question can be settled without an authenticated table. **This is an interpretation
blocker independent of the missing magnitude check**, and it is exactly the failure mode ADR-0002
warns against: a designated implementation returning a plausible number is not evidence that the
number is the governed quantity.

### The low-temperature model is an extrapolation of five orders of magnitude

`dQ2dT_ForCooling` blends a fitted low-`T` slope `a0` into the tabulated derivative
(`CompOSE_Thermo.cpp:700-706`). `a0` is a least-squares slope through the origin over the first
`lowT_fit_points` temperature grid points — **default 3**, with the in-code comment "uses
T = 0, 2, 4 MeV by default" (`CompOSE_Thermo.hpp:71`). The switch is `lowT_switch_MeV = 2.0` with
blend width `1.0`, and `BlendLowT_` returns **exactly 0** for `T ≤ 1 MeV`.

The live thermal program does not override any of these (`spin_therm_evol_2_main.cpp:153-156`
sets only `use_central_difference`, `clamp_to_domain`, `Tmin_for_derivative_MeV`).

Two consequences, both material:

1. **Across the entire neutron-star cooling range** (`10⁶–10⁹ K` ≈ `10⁻⁵–10⁻⁴ MeV`) the blend
   weight is 0, so `dQ2/dT ≡ a0`, a constant, and therefore **`c_V ∝ T` and `C_⋆ ∝ T∞` hold by
   construction, not by physics.** Any test of the low-`T` slope is close to tautological. It
   verifies plumbing; it cannot corroborate the physics. ADR-0002 §V1 item 2 must be read with
   this caveat.
2. **`a0` is sampled at MeV-scale temperatures and applied at `10⁻⁵` MeV** — a 5-decade
   extrapolation. For ideal degenerate matter `s ∝ T`, so the *form* is defensible; whether the
   *coefficient* measured near a few MeV still describes matter at `10⁸ K` is a physics question
   this verification cannot answer without a real table.

## 3. Tier A — self-contained analytic verification

Runs in a clean checkout with **no external EOS assets**. It exercises the real production path
`StarProfile → StarContext → GeometryCache → HeatCapacityStar_Tinf`, and loads its fixture
through the **real `CompOSE_Thermo` parser** — the thermo class is not mocked.

### Fixture (synthetic, deliberately non-physical)

Written in the genuine on-disk CompOSE format derived from the production parser:
axis files `"<dim> <N> v…"`, and `eos.thermo` rows `iT inb iYq Q1 Q2 … Q7 Nadd` with 1-based
indices.

```
Q2(T, n_B, Y_q) = 0.25 · T      (entropy per baryon; slope 1/MeV, independent of n_B, Y_q)
T grid  = {0, 0.25, 0.5, 0.75, 1.0} MeV     (entirely below the 1 MeV blend threshold)
n_B grid = {0.05, 0.30, 0.60, 1.00} fm⁻³
Y_q grid = {0, 0.10, 0.30, 0.50}
```

Because `Q2` is exactly linear in `T`, the least-squares slope reproduces `a0 = 0.25` exactly and
the tabulated derivative agrees, so `c_V = T · n_B · 0.25` holds analytically everywhere tested.

Star: `N` uniform points on `r ∈ [0, 10] km`, `ν = 0`, `Λ = 0`, `n_B = 0.30 fm⁻³`, one proton
species (CompOSE code `11`, `q = +1`) so `Y_q = 0.10`. Then `e^{−ν} = 1`, `WV = 4πr²`, and

```
C_⋆(T∞) = T∞ · n_B · a · k_B[MeV/K] · (MeV fm⁻³ → erg cm⁻³) · 1e15 · (4/3)π R³
```

### Unit chain — ADR-0002 §V1 item 1 · **PASS**

At `T∞ = 1e-5 MeV` (the first cache node, so no interpolation is involved):

| Quantity | Value (erg K⁻¹) |
|---|---|
| Independent direct integral | `4.337437e+35` |
| Production `HeatCapacityStar_Tinf` | `4.337437e+35` |
| Closed-form trapezoid | `4.337437e+35` |
| Closed-form continuum | `4.337437e+35` |

Direct-vs-closed-form and production-vs-closed-form agree to **< 1e-12 relative**. The test also
asserts the result is *not* `1e15`× smaller than expected, so **removing or mis-powering
`KM3_TO_CM3` fails the test.** The `k_B` value used for the reference is `8.617333262e-11`, matching
`CompOSE_Thermo.cpp` — deliberately, since `NeutrinoCooling_Details.cpp:214` carries a different
precision (`8.617333262145e-11`); that split is INV-02's recorded two-precision `k_B` issue and is
not adjudicated here.

### Low-`T` scaling — item 2 · **PASS, with the §2 caveat**

`d ln C_⋆ / d ln T∞ = 1.000180` over `1e-5 → 1e-4 MeV`. The `1.8e-4` deviation is the `T∞`-cache
interpolation error, consistent with §V1 item 5 below. **This confirms correct plumbing only** —
linearity is enforced by the production low-`T` model (§2).

### Direct integral vs production cache · **PASS**

The reference re-implements the ADR integral from public inputs (`R`, `WV`, `ExpMinusNu`, `n_B`,
and `Y_q` reconstructed independently as `Σ qᵢ Yᵢ`) and **never calls
`HeatCapacityStar_Tinf`**.

| `T∞` (MeV) | cached | direct | rel. err |
|---|---|---|---|
| 1.0e-5 | 4.337437e+35 | 4.337437e+35 | 0 (cache node) |
| 3.0e-5 | 1.301726e+36 | 1.301231e+36 | 3.80e-4 |
| 1.0e-4 | 4.339231e+36 | 4.337437e+36 | 4.13e-4 |
| 3.0e-4 | 1.301321e+37 | 1.301231e+37 | 6.86e-5 |
| 1.0e-3 | 4.340153e+37 | 4.337437e+37 | 6.26e-4 |
| 1.0e-2 | 4.340180e+38 | 4.337437e+38 | 6.32e-4 |
| 1.0e-1 | 4.339283e+39 | 4.337437e+39 | 4.26e-4 |

Max **6.3e-4**, within the 1 % criterion.

### Radial convergence — item 4 · **PASS, order measured**

Evaluated at the first cache node so the measurement is the radial quadrature alone.

| N | rel. err vs continuum | observed order |
|---|---|---|
| 101 | 5.000e-05 | — |
| 201 | 1.250e-05 | **2.000** |
| 401 | 3.125e-06 | **2.000** |
| 801 | 7.813e-07 | **2.000** |

Observed order **exactly 2.000**, matching the nominal `O(Δr²)` of the trapezoid on smooth data
(INV-13). The order was measured, not assumed.

> **Methodological finding worth carrying forward.** At *non-node* temperatures the observed order
> collapsed to ≈ 0 — errors plateaued at ~4.1e-4 for every `N`. That is the `T∞`-interpolation
> floor (~6.6e-4) masking the radial error, not a defect in the quadrature. Any future convergence
> study of this integral must evaluate at cache nodes or bypass the cache, or it will measure the
> wrong thing.

### `T∞`-grid sensitivity — item 5 · **PASS, production `NT` unmodified**

Surrogate caches were built in *test* code — log-spaced nodes, values from the direct integral,
production's interpolation rule — so the production `NT = 160` constant was never touched.

| `NT` | max interpolation rel. err |
|---|---|
| 40 | 1.09e-2 |
| 80 | 2.54e-3 |
| 160 | 6.55e-4 |
| 320 | 1.63e-4 |
| **production (160)** | **6.55e-4** |

Second-order in node spacing, as expected. `NT = 160` meets the 1 % criterion with roughly 15×
margin; `NT = 40` would **not**. The error arises because `C_⋆` is linear in `T∞` while the scheme
interpolates linearly against `log T`.

### Endpoint clamping — item 6 · **CHARACTERIZED, NOT ENDORSED**

> **CURRENT NUMERICAL BEHAVIOR — NOT YET GOVERNED AS CORRECT.**

- `T∞ < 1e-5 MeV` → returns the low endpoint value exactly.
- `T∞ > 1 MeV` → returns the high endpoint value exactly.
- At `T∞ = 1e-9 MeV` the clamp **overstates `C_⋆` by a factor of 1.0e4** relative to the direct
  integral — precisely as expected for a quantity linear in `T∞` that is frozen at `1e-5 MeV`.

`1e-5 MeV ≈ 1.16e5 K` is below the neutron-star cooling range, so ordinary Phase-2A use should not
reach the clamp. It is recorded, not fixed; INV-10 retains adjudication.

### Cache correctness — item 7 · **PASS on version/thermo; HAZARD CONFIRMED on geometry**

| Case | Result |
|---|---|
| Repeated identical query | Bit-identical. **PASS** |
| Profile-version mutation (`n_B` doubled via `EditScope`) | Version `1 → 2`; cache rebuilt; `C_⋆` ratio `2.000000` — correct direction and magnitude. **PASS** |
| Distinct `CompOSE_Thermo` object, slope doubled | Cache rebuilt; ratio `2.000000`. **PASS** |
| **Same profile version + same thermo, different `GeometryCache`** | **HAZARD CONFIRMED** |

For the last case a second `GeometryCache` built with `e^{Λ} = 2` — which should double the proper
volume and hence `C_⋆` — was passed on a later call. Production returned
`4.339231e+36` both times, **bit-identical**: the stale table was reused because the cache key is
only `(profile version, thermo pointer)` (`StarContext.cpp:712-714`). This is the INV-12 hazard,
now demonstrated rather than inferred. **Not fixed here.**

Practical exposure is limited while a single `GeometryCache` is built once per star and reused —
which is what `spin_therm_evol_2_main.cpp` does — but nothing in the API enforces that.

## 4. Tier B — real canonical star · **BLOCKED**

**No authenticated CompOSE EOS assets exist on this machine.** Searched: the worktree, the main
checkout, and `$HOME` to depth 8. No `eos.thermo`, `eos.t`, `eos.nb`, or `eos.yq`; no
`DS(CMF)-1_with_crust`, no `CMF-1_general`, and no `EOS/CompOSE/` directory at the location the
live program expects (`<repo>/EOS/CompOSE/`, from
`spin_therm_evol_2_main.cpp:115-119`). The directory is not gitignored — it is simply absent.

No replacement EOS was downloaded and no substitute table was used, per the task constraints.
**No Tier-B executable was written**, because a scientific integration test that cannot be run
even once is unverified candidate code (`GOVERNANCE.md` §5); the next task should add it together
with the data.

Consequently **not performed**: achieved mass/radius/central density for a canonical `1.4 M☉`
star; `C_⋆(10⁸ K)` on a real profile; the literature magnitude comparison; the real-star
`C_⋆/T∞` trend; and the CompOSE table-domain coverage audit (`n_B(r)`, `Y_q(r)`, `T_local(r)`
versus the table grids, and how much of the integral comes from clamped thermodynamics).

**No SHA256 hashes are recorded because no data was used.**

## 5. ADR-0002 §V1 checklist

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | Dimensional/unit check end to end, incl. `KM3_TO_CM3` | ✅ **PASS** | Exact to < 1e-12 against closed form; a missing `1e15` fails the test |
| 2 | Analytic limit — `c_V ∝ T`, so `C_⋆` linear in `T∞` | ⚠️ **PASS (weak evidence)** | Slope `1.000180`. But production enforces linearity by construction below 1 MeV (§2) — this cannot corroborate the physics |
| 3 | Order-of-magnitude vs published `10³⁷–10³⁸ erg K⁻¹` at `10⁸ K` | ⛔ **BLOCKED** | No EOS assets; Tier B not run |
| 4 | Grid convergence in `r`, second order expected | ✅ **PASS** | Observed order **exactly 2.000** across three refinements |
| 5 | Grid convergence in `T` (`NT` sensitivity) | ✅ **PASS** | `NT=160` → 6.55e-4, within the 1 % criterion; `NT=40` would fail |
| 6 | Clamping behavior stated explicitly | ✅ **PASS (recorded, not endorsed)** | Both endpoints characterized; 1e4× overstatement at 1e-9 MeV |
| 7 | Cache correctness | ⚠️ **PARTIAL** | Version and thermo rebuilds correct; **GeometryCache-key hazard CONFIRMED** |
| — | CompOSE thermodynamic content (lepton/crust) | ⛔ **BLOCKED** | `m_Il` parsed, never used; content unverifiable without a real table |

## 6. Conclusion

**V1 BLOCKED.**

The production *algorithm* is sound where it can be checked without real data: the unit chain is
exact, the radial quadrature converges at its nominal second order, the `T∞` cache is adequate at
`NT = 160` under a criterion fixed in advance, endpoint clamping behaves as documented, and the
cache rebuilds correctly on profile-version and thermo-object changes.

That is **not** sufficient to certify `C_⋆(T∞)` as the governed physical quantity. Two independent
blockers remain:

1. **The magnitude check cannot be performed** — the single most informative discriminator
   available (ADR-0002's own words), and the one that would distinguish a real `C_⋆` from the
   `1e40` placeholder it is meant to replace.
2. **The thermodynamic content of the CompOSE entropy is unestablished** — the lepton flag is
   read and discarded, so nothing confirms the integral is the *total* heat capacity ADR-0002
   requires rather than a baryon-only subset.

Per ADR-0002 and the Phase-2A gate, **the `PhotonCooling` conformance change must not proceed on
this evidence.**
