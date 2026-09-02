// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2026
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file Units.hpp
 * @brief Single owner for unit-conversion factors that were previously re-declared,
 *        with numerically identical values, in more than one production translation unit.
 *
 * Scope — deliberately narrow (roadmap Phase 3A, `docs/architecture/PHASE3_CONSOLIDATION_PLAN.md`
 * increment 3A). This header owns ONLY conversions whose duplicated declarations were already
 * bit-for-bit identical, so that centralizing them is behavior-preserving by construction.
 * It is **not** a units library: no strong types, no dimensional analysis, no user-defined
 * literals, no conversion functions, no runtime state.
 *
 * Explicitly NOT owned here, each for a stated reason:
 *   - k_B (MeV/K). RESOLVED in roadmap increment 3C, and deliberately NOT owned here: it is a
 *     fundamental constant, so its authority is the reusable library. All four production
 *     consumers now use `Zaki::Physics::K_BOLTZ_EV * 1.0e-6`. (Before 3C they carried two
 *     divergent literals, `8.617333262145e-11` in the thermal drivers and `8.617333262e-11`
 *     in `CompOSE_Thermo`.) See docs/validation/PHASE3C_BOLTZMANN_AUTHORITY.md.
 *   - The solar-mass conversions (`Zaki::Physics::SUN_M_KM` vs `GSL_CONST_CGSM_SOLAR_MASS`).
 *     These disagree at ~6.2e-5 and choosing one is a scientific/unit authority decision
 *     requiring owner or ADR adjudication — explicitly deferred out of Phase 3.
 *   - G, c and the other GSL physical constants; pressure/energy-density geometric conversions;
 *     fm^-3 <-> km^-3; year <-> second. None was a bit-identical duplicate.
 *
 * Layering. This header includes nothing — not the standard library, not Zaki, not Confind, and
 * no other CompactStar component. It therefore adds no edge to the dependency graph and may be
 * included from any layer. That is why it sits at the top level rather than inside `EOS` or
 * `Physics`: `EOS` currently includes nothing from `Physics`, so placing it under `Physics`
 * would create a new reversed dependency, and `Core` already depends on both.
 *
 * Test-side literals are intentionally NOT routed through this header. Several tests
 * (notably `tests/thermal/heat_capacity_v1.cpp`) act as oracles for exactly these factors —
 * `heat_capacity_v1` exists in part to detect a missing or mis-powered 10^15 — and would lose
 * their independence if they imported the production constant they are meant to verify.
 */

#ifndef CompactStar_Units_H
#define CompactStar_Units_H

namespace CompactStar::Units
{

/// Cubic kilometres to cubic centimetres: (1e5 cm)^3. Exact by definition.
/// Previously declared identically in `Physics/Evolution/src/StarContext.cpp` and
/// `Physics/Driver/Thermal/src/NeutrinoCooling_Details.cpp`.
inline constexpr double KM3_TO_CM3 = 1.0e15;

/// MeV fm^-3 to erg cm^-3 [CompOSE manual].
/// Previously declared identically twice in `EOS/src/CompOSE_Thermo.cpp`.
inline constexpr double MEV_FM3_TO_ERG_CM3 = 1.602176634e33;

} // namespace CompactStar::Units

#endif /* CompactStar_Units_H */
