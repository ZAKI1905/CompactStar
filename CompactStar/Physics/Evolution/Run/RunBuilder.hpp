#pragma once
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file RunBuilder.hpp
 * @brief Small orchestration helpers to reduce boilerplate in Evolution debug mains.
 *
 * This module provides narrowly-scoped “builder” utilities for:
 *  - configuring Evolution::Config with common defaults,
 *  - wiring DriverContext from already-built StarContext/GeometryCache,
 *  - building StateVector, StateLayout, RHSAccumulator consistently,
 *  - collecting diagnostics-capable drivers.
 *
 * Important: this is NOT intended to be a “simulation framework” that builds
 * stars, chooses EOS, or owns application concerns. Star construction belongs
 * to Core (e.g. NStar, SolveTOV_Profile) or higher-level application code.
 *
 * @ingroup Evolution
 */

#include <cstddef>
#include <memory>
#include <vector>

#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
#include "CompactStar/Physics/Evolution/RHSAccumulator.hpp"
#include "CompactStar/Physics/Evolution/StateLayout.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/Tags.hpp"

// Diagnostics interface (used only for collecting pointers).
#include "CompactStar/Physics/Driver/Diagnostics/DriverDiagnostics.hpp"

namespace CompactStar::Physics::Evolution
{
class StarContext;
class GeometryCache;
} // namespace CompactStar::Physics::Evolution

namespace CompactStar::Physics::Evolution::Run
{

/**
 * @brief Bundle of Evolution state wiring objects commonly needed by a main().
 *
 * A debug main typically needs all of these:
 * - StateVector: block registry and view into concrete state blocks
 * - StateLayout: packing order and total dimension
 * - RHSAccumulator: per-tag RHS accumulation buffers
 *
 * This struct keeps them together for ergonomic use.
 *
 * Ownership:
 * - StateVector stores references/copies depending on your implementation; this struct
 *   only stores the StateVector instance itself.
 * - The concrete state blocks (ThermalState, SpinState, etc.) are owned by your main.
 */
struct StateWiring
{
	Evolution::StateVector state_vec;
	Evolution::StateLayout layout;
	Evolution::RHSAccumulator rhs;

	/// Total packed dimension (layout.TotalSize()).
	std::size_t dim = 0;
};

/**
 * @brief Create a Config with safe/typical defaults for debug runs.
 *
 * @return Config with deterministic defaults:
 *   - RKF45
 *   - rtol=1e-6, atol=1e-10
 *   - dt_save=1e5
 *   - max_steps=1e6
 *
 * Caller can override fields after return.
 */
Evolution::Config MakeDefaultConfig();

/**
 * @brief Wire a DriverContext from pre-built StarContext/GeometryCache and Config.
 *
 * @param star StarContext (must outlive ctx).
 * @param geo  GeometryCache (must outlive ctx).
 * @param cfg  Evolution::Config (must outlive ctx).
 *
 * @return DriverContext with star/geo/cfg set and envelope set to nullptr.
 */
Evolution::DriverContext MakeDriverContext(Evolution::StarContext &star,
										   Evolution::GeometryCache &geo,
										   Evolution::Config &cfg);

/**
 * @brief Wire a DriverContext from pre-built StarContext/GeometryCache, Config, and Envelope.
 *
 * @param star StarContext (must outlive ctx).
 * @param geo  GeometryCache (must outlive ctx).
 * @param cfg  Evolution::Config (must outlive ctx).
 * @param env  Envelope model (must outlive ctx).
 * @return * Evolution::DriverContext
 */
Evolution::DriverContext MakeDriverContext(Evolution::StarContext &star,
										   Evolution::GeometryCache &geo,
										   Evolution::Config &cfg,
										   const Physics::Driver::Thermal::Boundary::IEnvelope *env);

/**
 * @brief
 *
 * @param star StarContext (must outlive ctx).
 * @param geo GeometryCache (must outlive ctx).
 * @param cfg Evolution::Config (must outlive ctx).
 * @param thermo EOS::CompOSE_Thermo (must outlive ctx).
 * @return Evolution::DriverContext
 */
Evolution::DriverContext MakeDriverContext(Evolution::StarContext &star,
										   Evolution::GeometryCache &geo,
										   Evolution::Config &cfg,
										   const EOS::CompOSE_Thermo *thermo);
/**
 * @brief
 *
 * @param star StarContext (must outlive ctx).
 * @param geo GeometryCache (must outlive ctx).
 * @param cfg Evolution::Config (must outlive ctx).
 * @param env Envelope model (must outlive ctx).
 * @param thermo EOS::CompOSE_Thermo (must outlive ctx).
 * @return Evolution::DriverContext
 */
Evolution::DriverContext MakeDriverContext(Evolution::StarContext &star,
										   Evolution::GeometryCache &geo,
										   Evolution::Config &cfg,
										   const Physics::Driver::Thermal::Boundary::IEnvelope *env,
										   const EOS::CompOSE_Thermo *thermo);

/**
 * @brief Configure StateLayout packing order for the given tags.
 *
 * @param wiring StateWiring with state_vec already registered.
 * @param order  Tag order for packing (authoritative).
 */
void ConfigureLayout(StateWiring &wiring,
					 const std::vector<Physics::State::StateTag> &order);

/**
 * @brief Configure RHSAccumulator buffers for a set of registered tags.
 *
 * This reads the sizes from the StateVector blocks and configures rhs accordingly.
 *
 * @param wiring StateWiring with state_vec already registered.
 * @param tags   Tags to configure for RHS accumulation.
 */
void ConfigureRHS(StateWiring &wiring,
				  const std::vector<Physics::State::StateTag> &tags);

/**
 * @brief Collect diagnostics-capable drivers into a non-owning pointer list.
 *
 * This is a tiny helper so mains don’t repeat the dynamic_cast boilerplate.
 *
 * @param drivers Vector of shared_ptr<IDriver> (or DriverPtr).
 * @return Vector of pointers to drivers that implement IDriverDiagnostics.
 *
 * @note The returned pointers are non-owning and remain valid only as long as
 *       the underlying shared_ptr instances remain alive.
 */
std::vector<const Physics::Driver::Diagnostics::IDriverDiagnostics *>
CollectDiagnosticsDrivers(const std::vector<Evolution::DriverPtr> &drivers);

} // namespace CompactStar::Physics::Evolution::Run