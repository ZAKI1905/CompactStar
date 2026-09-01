# Phase 3A — centralizing exactly-duplicated unit constants

> **RESULT: `BIT-IDENTICAL ENGINEERING CONSOLIDATION`.**
>
> All 352 lines of captured deterministic scientific output are **bit-identical** before and
> after. All five golden artifacts are **byte-identical**. 13/13 and 8/8 CTests pass. No test
> source, baseline, ADR or invariant was touched.

| Field | Value |
|---|---|
| **Starting commit** | `3a2be24cec26ffb9426e286765380a58032ea72a` |
| **Change class** | engineering / behavior-preserving (roadmap Phase 3A, plan increment U-1) |
| **ADR** | **none required** — the entry audit classified exact-duplicate consolidation as engineering, and source inspection confirmed it |
| **Baseline at entry** | 13/13 authenticated, 8/8 self-contained |

---

## 1. The two conversions

| Conversion | Literal | Exact? |
|---|---|---|
| km³ → cm³ | `1.0e15` | exact by definition, `(10^5 cm)^3` |
| MeV fm⁻³ → erg cm⁻³ | `1.602176634e33` | fixed by the SI-exact elementary charge [CompOSE manual] |

Both were already **numerically identical** at every production site, which is what makes this
consolidation behavior-preserving by construction rather than by measurement.

## 2. Pre-change production occurrence map

| Value | Path | Function | Kind | Same value? | In scope |
|---|---|---|---|---|---|
| `1.0e15` | `Physics/Evolution/src/StarContext.cpp:761` | `BuildHeatCapacityCache_` | production | yes | **YES** |
| `1.0e15` | `Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp:94` | `BuildNeutrinoCoolingCache` | production | yes | **YES** |
| `1.602176634e33` | `EOS/src/CompOSE_Thermo.cpp:569` | `CvDensity_cgs` | production | yes | **YES** |
| `1.602176634e33` | `EOS/src/CompOSE_Thermo.cpp:726` | `CvDensity_cgs_ForCooling` | production | yes | **YES** |

**Excluded after inspection — same token, different physical meaning:**

| Site | Expression | Why excluded |
|---|---|---|
| `Microphysics/BNV/Analysis/src/Decay_Analysis.cpp:111,140` | `out *= LIGHT_C_M_S * 1e15` | a **length** conversion (m→fm), not km³→cm³ |
| `Microphysics/BNV/Internal/src/BNV_Chi.cpp:93` | `out *= LIGHT_C_M_S * 1e15` | same |
| `NeutrinoCooling_Details.cpp:140` (comment), `:264,283` (commented code) | `rho15 = rho/1e15` | a **density normalization** to 10¹⁵ g cm⁻³ |
| `NeutrinoCooling_Details.cpp:223` | commented-out declaration | not a declaration |
| `Core/src/TaskManager.cpp:463` | commented plot range | not a conversion |

Consolidating any of these under a km³→cm³ name would have been a semantic error, not a
deduplication.

## 3. The owner and why its layer is correct

**`CompactStar/Units.hpp`** — new, header-only, **includes nothing at all** (not the standard
library, not Zaki, not Confind, not another CompactStar component), `inline constexpr double`,
namespace `CompactStar::Units`.

Layer rationale, from the measured dependency graph:

| Candidate | Verdict |
|---|---|
| `CompactStar/Physics/...` | **rejected** — `EOS → Physics` is currently **0 files**; putting it there would create a new *reversed* dependency edge for `CompOSE_Thermo` |
| `CompactStar/EOS/...` | **rejected** — `KM3_TO_CM3` is a pure geometric factor with no EOS meaning; it would make `Physics → EOS` carry a non-EOS concern |
| `CompactStar/Core/...` | **rejected** — `Core` already includes both `EOS` and `Physics`; it is the *highest* layer, not a base |
| **`CompactStar/Units.hpp` (top level)** | **selected** — depends on nothing, so it adds **no edge** to the graph and any layer may include it |

**No CMake change was needed.** The repo root is already a `PUBLIC` include root
(`CMakeLists.txt:132-139`) and `CompactStar_SRC_Files` (`CompactStar/CMakeLists.txt:7`) lists
only source files, so a header-only addition needs no registration and there is no header
install rule to update.

**Deliberately not a units library**: no strong types, no dimensional analysis, no user-defined
literals, no conversion functions, no runtime state. Two constants, nothing more.

## 4. Explicitly out of scope

| Constant | Status |
|---|---|
| `k_B` — `8.617333262145e-11` (thermal) vs `8.617333262e-11` (`CompOSE_Thermo`) | **untouched.** Two precisions differing at `1.7e-11`; unifying changes numbers → roadmap **3C** |
| `Zaki::Physics::SUN_M_KM` vs `GSL_CONST_CGSM_SOLAR_MASS` | **untouched.** Disagree at `6.2e-5`; a scientific/unit authority decision → **deferred out of Phase 3** |
| `G`, `c`, other GSL constants; pressure/energy-density geometric conversions; fm⁻³↔km⁻³; year↔second | **untouched** — none was a bit-identical duplicate |
| `INV_1E9_POW6/8`, `Q0_DU`, `Q0_MU` | **untouched** — single-site, not duplicates |
| proper volume, cache architecture, TOV ownership | **untouched** — increments 3B/3D/3E |

Expression structure was preserved exactly: only the *spelling* of one operand changed at each
use site. No multiplication was reordered and no combined coefficient was precomputed — in
particular `kB_MeV_per_K * MEV_FM3_TO_ERG_CM3` remains two separate factors in the original
order.

## 5. Why test literals were retained

`tests/thermal/heat_capacity_v1.cpp:62-63`, `heat_capacity_real_star.cpp:57` and
`hartle_moment_inertia_cmf.cpp:66` keep their own literals, **unchanged**.

These tests are **oracles for exactly these factors**. `heat_capacity_v1` contains check `U1.c`,
*"magnitude excludes a missing 1e15"*, which exists specifically to catch a missing or
mis-powered `10^15` in production. If its reference calculation imported the production constant
it is meant to verify, the check would become circular and stop detecting the defect it was
written for. `hartle_moment_inertia_cmf.cpp:66` is a different composite
(`1.0e15 * c²/G` → g cm²), not the same quantity.

**No test source file was modified.** `git diff -- tests` is empty.

## 6. Behavior-preservation evidence

Deterministic scientific output was captured from six existing harnesses **before** the edit
(from the unmodified build) and again after, filtered to remove absolute paths only
(`s|/[A-Za-z0-9_./-]*||g`); no numeric field was filtered and no tolerance was applied.

| Harness | Lines compared | Result |
|---|---|---|
| `heat_capacity_v1` | 42 | **IDENTICAL** |
| `heat_capacity_real_star` | 247 | **IDENTICAL** |
| `photon_cooling_conformance` | 18 | **IDENTICAL** |
| `cache_thermal_contract` | 10 | **IDENTICAL** |
| `passive_cooling_regression` (9 epochs vs golden) | 17 | **IDENTICAL** |
| `hartle_moment_inertia_cmf` | 18 | **IDENTICAL** |
| **total** | **352** | **bit-identical** |

Compile-time confirmation that the literals themselves are unchanged:

```cpp
static_assert(CompactStar::Units::KM3_TO_CM3 == 1.0e15, ...);
static_assert(CompactStar::Units::MEV_FM3_TO_ERG_CM3 == 1.602176634e33, ...);
```

Both pass. No dedicated test executable was added — the existing scientific tests already
exercise both consumers.

## 7. Regression suites

| Configuration | Before | After |
|---|---|---|
| Full, authenticated data root | 13/13 PASS | **13/13 PASS** |
| Self-contained | 8/8 PASS | **8/8 PASS** |

## 8. Golden artifacts — byte-identical

| Artifact | SHA-256 before | after |
|---|---|---|
| `passive_cooling_cmf_1p6_debug.tsv` | `edaa5e3b…f42d6e` | **same** |
| `tov_dscmf1_reference.tsv` | `ba9f6ee5…c13e28` | **same** |
| `grid_convergence_cmf_1p6_debug.tsv` | `3be2005f…cc3c88` | **same** |
| `grid_convergence_cmf_1p6_trajectory.tsv` | `e04f7485…a986194` | **same** |
| `hartle_I_dscmf1_debug.tsv` | `ddf01857…aa9d15` | **same** |

`git diff -- tests/baselines` is empty.

## 9. Post-change duplicate search

```
km^3 -> cm^3 production declarations of the value : CompactStar/Units.hpp:54   (one)
MeV fm^-3 -> erg cm^-3 production declarations    : CompactStar/Units.hpp:58   (one)
```

Production consumers now reference the owner at `StarContext.cpp:810,811`,
`NeutrinoCooling_Details.cpp:185,197`, `CompOSE_Thermo.cpp:571,725`.

The claim is **one production owner per conversion** — not "one occurrence in the repository".
Tests legitimately keep independent literals (§5) and documentation quotes the values.

## 10. Conclusion

**`BIT-IDENTICAL ENGINEERING CONSOLIDATION`.** The 3A classification from the Phase-3 entry
audit held under implementation: nothing in the source contradicted it, no ADR was required, and
no numerical result moved.
