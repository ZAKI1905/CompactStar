# ADR-0002 §V1 — Verification of the canonical stellar heat capacity `C_⋆(T∞)`

> **STATUS: V1 VERIFIED.**
> Tier-A analytic verification (Phase 2A-1) and Tier-B verification on an authenticated
> canonical `1.4 M☉` star (Phase 2A-2) are both complete. All seven ADR-0002 §V1 items are
> adequately addressed. `C_⋆(10⁸ K) = 2.17e38 erg K⁻¹`, within the same order of magnitude as the
> conventional neutron-star estimate (`~1e38`), differing by a factor of ~2.2.
>
> Two documented hazards remain and do **not** invalidate the Phase-2A use: the endpoint clamp
> (below the cooling range) and the `GeometryCache` cache-key omission (INV-12), which cannot fire
> when one `GeometryCache` is built per star, as production does.
>
> **Supersedes the Phase-2A-1 status of `V1 BLOCKED`.** Both blockers are resolved: the historical
> data were located and authenticated, and the finite-`T` table's lepton content is established.

| Field | Value |
|---|---|
| **Governing authority** | `docs/adr/ADR-0002-thermal-heat-capacity-ownership.md` §V1; INV-10, INV-12, INV-13, INV-15 |
| **Code commit** | Tier A `e2cc927`; Tier B `f467b01` |
| **Build configuration** | `Debug` (Phase-1D default) |
| **Toolchain** | AppleClang 17.0.0.17000604, CMake 4.2.1, macOS 15.6.1 (24G90), arm64 |
| **Verification programs** | `tests/thermal/heat_capacity_v1.cpp` (CTest `heat_capacity_v1`, Tier A, no external data) · `tests/thermal/heat_capacity_real_star.cpp` (CTest `heat_capacity_real_star`, Tier B, labels `thermal;external-data;scientific`) |
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

## 2. CompOSE thermodynamic content — **RESOLVED**

### Correction to the Phase-2A-1 claim

Phase 2A-1 stated that "nothing establishes whether entropy includes leptons". **As a
format-level claim that was too strong and is withdrawn.** The official CompOSE specification
defines the first `eos.thermo` row as `m_n  m_p  I_l`, where `I_l = 1` means the EoS includes
electrons and/or muons, and quantity `Q2 = s/n_B` is the **total** entropy per baryon of the
tabulated EoS. The format does answer the question; Phase 2A-1 simply had no table to read.

What remains correct from Phase 2A-1 is the narrower engineering observation: **`m_Il` is parsed
into a member (`CompOSE_Thermo.cpp:291`) and never read anywhere else in the class.** The code
does not assert, branch on, or report the lepton content of whatever table it is handed. That is
a robustness gap, not a physics defect — the table must be authenticated externally, which this
document now does.

### The authenticated finite-`T` table

`eos.thermo` first line, read directly:

```
938.91871259  938.91871250     1
```

so `m_n = 938.91871259 MeV`, `m_p = 938.91871250 MeV`, **`I_l = 1` → leptons are included.**

Its bundled `README` (official CompOSE distribution, dated 2017):

> "Eos table, CMF hadronic EoS (with electrons), downloaded from http://compose.obspm.fr"

with references Dexheimer & Schramm (ApJ 683, 943, 2008); Schurhoff, Schramm & Dexheimer
(ApJ 724, L74, 2010); Dexheimer, Negreiros & Schramm (PRC 92, 012801, 2015).

`eos.init` records the grids and particle list:

| Property | Value |
|---|---|
| `T` | 81 points, `0 – 160 MeV`, uniform step **2 MeV** |
| `n_B` | 301 points, `1.0e-2 – 3.01 fm⁻³` |
| `Y_q` | 54 points, `0 – 0.53` |
| Particle codes | `0, 1, 11, 10, 100, 112, 111, 110, 121, 120, 500, 501, 502` |

Decoding the CompOSE codes: **`0` = electron, `1` = muon**, `10` = neutron, `11` = proton,
`100` = Λ, `110/111/112` = Σ⁻/Σ⁰/Σ⁺, `120/121` = Ξ⁻/Ξ⁰, `500/501/502` = mesons.

**Both electrons and muons appear in the particle list.** Note the mild tension worth recording:
the table *title* says "with electrons" while the particle list also carries the muon code. `I_l = 1`
covers "electrons and/or muons" either way, so the entropy is lepton-inclusive; the exact muon
treatment is not further resolved here and is not material at the 0.2 % level that separates the
components audited below.

### Total vs partial — classification per ADR-0002

**PARTIAL `C_V`, by a quantified and small margin.**

The tabulated entropy covers baryons, hyperons and leptons over `n_B ≥ 0.01 fm⁻³`. It does **not**
cover the crust: no ion lattice, no nuclei. The cold structural star extends to
`n_B ≈ 1.06e-7 fm⁻³`, five decades below the finite-`T` floor, and those zones are evaluated at the
table edge under `clamp_to_domain = true`.

Measured on the canonical star (§5): the sub-floor region is **92 of 2646 radial zones (3.5 %)**,
confined to the outer **0.64 km** of a `13.55 km` star, and contributes **0.196 %** of `C_⋆(10⁸ K)`.
The missing crust physics therefore cannot move the result at the order-of-magnitude level being
tested. It is recorded as a caveat on calling the quantity `C_⋆,total`, not as a blocker.

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

## 4. Data provenance and hashes

### How the data were located

Phase 2A-1 found no CompOSE assets and declared Tier B blocked. The identity of the local
directory name `CMF-1_general` was then recovered from **repository evidence**: the committed run
log `main/Test/results/spin_therm_evol_2/spin_therm_evol_2_main.log:27` records the absolute path
the historical run actually read —

```
.../GoogleDrive-m.zakeri@eku.edu/My Drive/Research/Tools/Coding/CompactStar/main/EOS/CompOSE/DS(CMF)-1_with_crust/...
```

That directory still exists and contains both `CMF-1_general/` and `DS(CMF)-1_with_crust/`, each
with its **original official CompOSE `README` and `eos.init`**, file-dated 2017–2021. **Nothing was
downloaded**; the authenticated originals were used.

### Identity — what is and is not established

| Claim | Basis | Status |
|---|---|---|
| The finite-`T` table is the official CompOSE **CMF hadronic EoS (with electrons)** | its own bundled `README`: "downloaded from http://compose.obspm.fr" | **AUTHENTICATED** (official distribution file) |
| Its grids are `81 × 301 × 54`, `T ≤ 160 MeV`, `n_B 0.01–3.01`, `Y_q 0–0.53` | `eos.init` + parser output | **AUTHENTICATED** (read directly) |
| The cold table is the official **CMF #1 EoS with crust** | its `README`: "CMF #1 EoS with crust, downloaded from https://compose.obspm.fr", Dexheimer et al. refs | **AUTHENTICATED** |
| `CMF-1_general` (the local name in source) denotes that finite-`T` table | the historical run log path + directory contents | **AUTHENTICATED** as owner-supplied local naming |
| These correspond to CompOSE catalogue **EoS ID 115** and **ID 188** | dimensions and titles match the IDs named in the task | **PLAUSIBLE, NOT INDEPENDENTLY CONFIRMED** — the live CompOSE catalogue was not queried and no official checksum file was compared |

The numeric EoS IDs are therefore **inference**, recorded as such. The physical identity of the
tables does not depend on them: it rests on the official README and `eos.init` shipped with the
data.

### Governed data root and SHA256

Copied (not symlinked, since the source is a cloud mount) to a local root outside tracked source.
`eos.compo` (408 MB) was not copied — it is not consumed by this path.

```
/Users/keeper/Documents/CompactStar/data/compose/
├── DNS-CMF-Hadronic-with-electrons/   (220 MB)
└── DS-CMF-1-with-crust/               (404 KB)
```

| SHA256 | File |
|---|---|
| `a456fb8595208ddf3119350a856fbf2b906c0a0e19bb7c716571748d0aa0724b` | `DNS-CMF-Hadronic-with-electrons/eos.thermo` |
| `2e4c6ec1feb85b16d0ee7036dce183782a9f681577e79c72315171069aa8513d` | `DNS-CMF-Hadronic-with-electrons/eos.t` |
| `3f79dbcc6f8b519696377f89ebc86464bc55cd61d9e2459f6e21e2d9e00f380d` | `DNS-CMF-Hadronic-with-electrons/eos.nb` |
| `d98fcd2f7752039c552c2ef2d04ab485b75db47a61f8ae1740875b54bf9824fd` | `DNS-CMF-Hadronic-with-electrons/eos.yq` |
| `48af68b1b4f6727252ae0051fc35c6445240241d654566973a068d20e1f35222` | `DNS-CMF-Hadronic-with-electrons/README` |
| `412f6739c769df650b3238a6e4b6f0d0f2d7d4a5df43e7d16c80e913bcaddbfb` | `DNS-CMF-Hadronic-with-electrons/eos.init` |
| `5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd` | `DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos` |
| `4e69b9193e0f075584239d818e1b459791da4d12427531914c86cdd209c898a8` | `DS-CMF-1-with-crust/README` |
| `8b4472405295655cf530572af7edd7448efd0393d0bf9ad86ead3e4c87228c90` | `DS-CMF-1-with-crust/eos.init` |

No data file is committed. The path appears in **no** committed source: the test receives it from
the CMake cache variable `COMPACTSTAR_EOS_DATA_ROOT`.

## 5. Tier B — authenticated canonical star · **COMPLETE**

```bash
cmake -S . -B build -DPython3_EXECUTABLE="$(command -v python3)"       -DCOMPACTSTAR_EOS_DATA_ROOT=/path/to/data/compose
cmake --build build -j8
ctest --test-dir build -L scientific --output-on-failure
```

Without the cache variable the Tier-B test is **not registered** and Tier A still runs; with it,
the test **exits non-zero if the data are missing** — it never silently skips.

### Canonical star, built through the production TOV path

`NStar::SolveTOV_Profile` → `StarProfile` → `StarContext` → `GeometryCache`.

| Property | Value |
|---|---|
| Requested mass | `1.400000 M☉` |
| **Achieved mass** | **`1.400022 M☉`** |
| **Radius** | **`13.545 km`** |
| Central energy density | `6.1649e14` (profile units) |
| Radial points | `2646` |
| Profile version | `5` |
| `n_B` centre → surface | `3.447e-1 → 1.056e-7 fm⁻³` |
| Strong-sector species found | 12 |
| `Y_q` range over the star | `[0, 0.1511]` vs table `[0, 0.53]` — fully inside |

### Structural / thermodynamic domain overlap

| Quantity | Value |
|---|---|
| Finite-`T` table `n_B` floor | `1.00e-2 fm⁻³` |
| Zones below the floor | **92 / 2646 (3.48 %)** |
| First such radius | `12.908 km` of `13.545 km` — the outer `0.64 km` |
| `Y_q` coverage | fully native, no clamping |
| `T_local` coverage | fully native (`T_local ≤ 0.09 MeV` ≪ `160 MeV`) |
| **Clamped contribution to `C_⋆(10⁸ K)`** | **`4.249e35` erg/K = 0.196 %** |
| Native-domain contribution | `2.167e38` erg/K = 99.80 % |

**Domain clamping neither dominates nor invalidates the result.** The mismatch is real — the cold
EoS carries a crust the finite-`T` table does not describe — but it is confined to the outermost
5 % of the radius and a fifth of a percent of the integral.

### `C_⋆(T∞)` on the real star

| `T∞` [K] | `T∞` [MeV] | `C_⋆` [erg K⁻¹] | `C_⋆/T∞` [erg K⁻²] | `d ln C/d ln T` |
|---|---|---|---|---|
| 1e6 | 8.617e-5 | 2.1725e+36 | 2.1725e+30 | 0.99973 |
| 3e6 | 2.585e-4 | 6.5156e+36 | 2.1719e+30 | 1.00038 |
| 1e7 | 8.617e-4 | 2.1729e+37 | 2.1729e+30 | 0.99989 |
| 3e7 | 2.585e-3 | 6.5178e+37 | 2.1726e+30 | 1.00006 |
| **1e8** | 8.617e-3 | **2.1727e+38** | 2.1727e+30 | 1.00005 |
| 3e8 | 2.585e-2 | 6.5186e+38 | 2.1729e+30 | 0.99973 |
| 1e9 | 8.617e-2 | 2.1722e+39 | 2.1722e+30 | — |

`C_⋆/T∞` is constant to 5 significant figures, as the low-`T` model enforces (§2 caveat).

### Literature magnitude — order-of-magnitude diagnostic

```
C_⋆(10⁸ K) = 2.17e38 erg K⁻¹        C_⋆/T∞ = 2.17e30 erg K⁻²
broad neutron-star expectation      ~1e37 – 1e38 erg K⁻¹
conventional estimate C ~ 1e39·T9   ~1e38 erg K⁻¹ at 10⁸ K
ratio to 1e38                        2.17
```

`2.17e38 erg K⁻¹` is **within the same order of magnitude** as the conventional `~1e38 erg K⁻¹`
estimate, differing by a factor of **~2.2**. Note it lies numerically *above* `1e38`, i.e. above
the quoted `10³⁷–10³⁸` band rather than at its upper end.

A factor of ~2 is the expected quality of agreement for a diagnostic of this kind applied to a
different EoS: the reference scaling is a generic non-superfluid estimate, while the value here
is specific to this CMF model's composition, mass and radius. **No attempt was made to decompose
the difference into individual causes** — attributing it to particular species or to the stellar
radius would require a separate controlled study that has not been performed. It is recorded
generically as EoS/composition/model dependence.

**No tolerance was tuned and no old cooling curve was used** — the comparison is against published
scaling only.

This is the check ADR-0002 called "the single most informative test available", and it is the one
that separates a genuine `C_⋆` from the `1e40 erg K⁻¹` `PhotonCooling` placeholder: the placeholder
is **46× larger** than the computed value at `10⁸ K`.

### Low-`T` fit sensitivity — the Phase-2A-1 concern, now measured

Phase 2A-1 flagged that `dQ2/dT` is fitted over the first few table temperatures and applied five
decades lower. The real table's grid is `0, 2, 4, …, 160 MeV`, so the default 3-point stencil is
exactly `T = 0, 2, 4 MeV` — as `CompOSE_Thermo.hpp:71` claims. Varying `lowT_fit_points` **in test
code only** (production defaults untouched):

| `lowT_fit_points` | fit spans `T` up to | `C_⋆(10⁸ K)` [erg K⁻¹] |
|---|---|---|
| 2 | 2 MeV | 2.1545e+38 |
| **3 (production default)** | 4 MeV | **2.1727e+38** |
| 4 | 6 MeV | 2.1486e+38 |
| 5 | 8 MeV | 2.1201e+38 |

**Spread: 1.025× (2.5 %).** The extrapolation is *methodologically* aggressive but *numerically*
benign for this table: `s/n_B` is close to linear in `T` from 0 to 8 MeV at these densities, so the
inferred coefficient barely depends on how many points anchor it. This resolves the Phase-2A-1
concern in the direction of confidence, and the criterion was that instability would block V1 —
2.5 % does not.

The residual caveat stands and is not removed by this measurement: the table has **no data
between 0 and 2 MeV**, so linearity below 2 MeV is assumed, not sampled. The stability above only
shows the fit is insensitive to which MeV-scale points are used.

## 6. ADR-0002 §V1 checklist

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | Dimensional/unit check end to end, incl. `KM3_TO_CM3` | ✅ **PASS** | Tier A: exact to < 1e-12 against closed form; a missing `1e15` fails the test |
| 2 | Analytic limit — `c_V ∝ T`, `C_⋆` linear in `T∞` | ✅ **PASS** | Tier A slope `1.000180`; real star `d ln C/d ln T = 1.0000` to 5 s.f. **Caveat: linearity is enforced by the low-`T` model below 1 MeV, so the meaningful check is the coefficient — see item 2b** |
| 2b | Low-`T` coefficient robustness | ✅ **PASS** | Spread across 2/3/4/5-point fit stencils: **1.025×** |
| 3 | Order-of-magnitude vs published `10³⁷–10³⁸ erg K⁻¹` at `10⁸ K` | ✅ **PASS** | **`2.17e38 erg K⁻¹`**, factor **2.2** from the conventional `1e39·T9` estimate |
| 4 | Grid convergence in `r`, second order expected | ✅ **PASS** | Tier A observed order **exactly 2.000** across three refinements |
| 5 | Grid convergence in `T` (`NT` sensitivity) | ✅ **PASS** | `NT=160` → `6.55e-4`, within the 1 % criterion; `NT=40` would fail |
| 6 | Clamping behavior stated explicitly | ✅ **PASS (recorded, not endorsed)** | Both endpoints characterized; `1e4×` overstatement at `1e-9 MeV`, which is below the cooling range and never reached in normal use |
| 7 | Cache correctness | ✅ **PASS with a documented hazard** | Version and thermo rebuilds correct; **`GeometryCache`-key omission CONFIRMED** but unreachable when one `GeometryCache` is built per star, as production does |
| — | CompOSE thermodynamic content | ✅ **RESOLVED** | `I_l = 1`; official README "with electrons"; particle list includes `e` and `μ` plus hyperons. **PARTIAL** only w.r.t. crust lattice/ions, quantified at **0.196 %** of the integral |
| — | Structural / thermo domain overlap | ✅ **QUANTIFIED** | 3.48 % of zones clamped, outer 0.64 km, 0.196 % of `C_⋆` |

## 7. Conclusion

**V1 VERIFIED.**

The production implementation `StarContext::HeatCapacityStar_Tinf` is verified as an adequate
realisation of the ADR-0002 canonical `C_⋆(T∞)` for the purposes of Phase 2A:

- the unit chain is exact and the `km³ → cm³` factor is detectably correct;
- the radial quadrature converges at its nominal second order;
- the `T∞` cache at `NT = 160` is accurate to `6.6e-4`, well inside the pre-declared 1 % criterion;
- on an authenticated canonical `1.4 M☉` CMF star, `C_⋆(10⁸ K) = 2.17e38 erg K⁻¹` — the same
  order of magnitude as the conventional estimate, differing by a factor of ~2.2, with the
  difference recorded generically as EoS/composition dependence rather than attributed to a
  specific cause — and **46× smaller than the `1e40 erg K⁻¹` placeholder it replaced**, which is
  the practically decisive result;
- the low-`T` coefficient is stable to 2.5 % against the fit stencil, resolving the Phase-2A-1
  concern;
- domain clamping is confined to 0.196 % of the integral.

### Caveats carried forward — none blocking Phase 2A

1. **Crust thermodynamics are not represented.** The finite-`T` table floors at `0.01 fm⁻³`; the
   outer `0.64 km` is evaluated at the table edge. Strictly the quantity is not `C_⋆,total`. At
   0.196 % this cannot affect Phase-2A conclusions, but it must not be forgotten if crust physics
   later matters (e.g. superfluid or magnetar-scale modelling).
2. **No table data between 0 and 2 MeV.** Linearity of `s/n_B` below 2 MeV is assumed, not
   sampled. The fit-stencil stability measured here does not test that assumption.
3. **`m_Il` is parsed and never used.** The code cannot detect being handed a lepton-free table.
   A cheap robustness improvement — read and report it — is recommended, but is a separate
   engineering change.
4. **`GeometryCache` cache-key omission (INV-12)** remains, demonstrated in Tier A. Safe only by
   convention; nothing in the API enforces one geometry per profile version.
5. **Endpoint clamping (INV-10)** remains ungoverned, though outside the cooling range.

### Consequence — conformance has since landed (Phase 2A-3)

This document supplied the `GOVERNANCE.md` §3.1 condition-4 evidence for the `PhotonCooling`
correction: independent verification that depends on no prior passive-cooling output.

**That correction is now in the tree.** `PhotonCooling` divides by the canonical `C_⋆(T∞)`
obtained from `StarContext::HeatCapacityStar_Tinf`, `PhotonCooling::Options::C_eff` has been
removed, and the ADR-0002 Pattern A additive-driver architecture is preserved. Measured on this
same authenticated star at `T∞ = 10⁸ K`, the instantaneous denominator ratio is **46.0**
(`1e40 / 2.1727e38`) — a *local* rate ratio only, not a claim about the cooling trajectory, which
is nonlinear in `T∞` and also driven by neutrino losses.

Remaining §3.1 obligation: **condition 7**, the regression baseline, which must follow immediately
(roadmap Phase 2B).
