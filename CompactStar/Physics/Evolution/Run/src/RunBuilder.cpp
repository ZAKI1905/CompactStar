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
 * @file RunBuilder.cpp
 * @brief Implementation of Evolution run builder helpers.
 *
 * @ingroup Evolution
 */

#include "CompactStar/Physics/Evolution/Run/RunBuilder.hpp"

// #include "CompactStar/Physics/Driver/Diagnostics/DriverDiagnostics.hpp"
#include "CompactStar/Physics/Driver/IDriver.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/State/SpinState.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

#include <Zaki/Util/Logger.hpp>

namespace CompactStar::Physics::Evolution::Run
{
//--------------------------------------------------------------
Evolution::Config MakeDefaultConfig()
{
	Evolution::Config cfg;

	// Coupling flags: leave conservative defaults; caller can override.
	cfg.couple_spin = true;
	cfg.n_eta = 0;

	// Integrator defaults
	cfg.stepper = Evolution::StepperType::RKF45;
	cfg.rtol = 1e-6;
	cfg.atol = 1e-10;
	cfg.max_internal_steps = 1000000;
	cfg.max_samples = 1000000;
	cfg.dt_save = 1.0e5;

	return cfg;
}
//--------------------------------------------------------------
Evolution::DriverContext MakeDriverContext(Evolution::StarContext &star,
										   Evolution::GeometryCache &geo,
										   Evolution::Config &cfg)
{
	// Evolution::DriverContext ctx;
	// ctx.star = &star;
	// ctx.geo = &geo;
	// // ctx.envelope = nullptr; // reserved for future
	// ctx.cfg = &cfg;
	// return ctx;

	return MakeDriverContext(star, geo, cfg, nullptr);
}

//--------------------------------------------------------------
Evolution::DriverContext MakeDriverContext(Evolution::StarContext &star,
										   Evolution::GeometryCache &geo,
										   Evolution::Config &cfg,
										   const Physics::Driver::Thermal::Boundary::IEnvelope *env)
{
	Evolution::DriverContext ctx;
	ctx.star = &star;
	ctx.geo = &geo;
	ctx.envelope = env; // <-- enable envelope wiring
	ctx.cfg = &cfg;
	return ctx;
}

//--------------------------------------------------------------
void ConfigureLayout(StateWiring &wiring,
					 const std::vector<Physics::State::StateTag> &order)
{
	wiring.layout.Configure(wiring.state_vec, order);
	wiring.dim = wiring.layout.TotalSize();
}
//--------------------------------------------------------------
void ConfigureRHS(StateWiring &wiring,
				  const std::vector<Physics::State::StateTag> &tags)
{
	// Policy: configure rhs buffer sizes for each tag from the state vector.
	// ADAPT HERE if your StateVector API differs.
	for (const auto tag : tags)
	{
		switch (tag)
		{
		case Physics::State::StateTag::Thermal:
			wiring.rhs.Configure(tag, wiring.state_vec.GetThermal().Size());
			break;

		case Physics::State::StateTag::Spin:
			wiring.rhs.Configure(tag, wiring.state_vec.GetSpin().Size());
			break;

		default:
			// Conservative: if a tag isn’t handled here, warn once and skip.
			// If you add more tags (Chemical, BNV, etc.), extend this switch.
			Z_LOG_WARNING("RunBuilder::ConfigureRHS: unhandled tag; extend switch for tag.");
			break;
		}
	}
}
//--------------------------------------------------------------
std::vector<const Physics::Driver::Diagnostics::IDriverDiagnostics *>
CollectDiagnosticsDrivers(const std::vector<Evolution::DriverPtr> &drivers)
{
	std::vector<const Physics::Driver::Diagnostics::IDriverDiagnostics *> out;
	out.reserve(drivers.size());

	for (const auto &d : drivers)
	{
		if (!d)
			continue;

		auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(d.get());
		if (p)
			out.push_back(p);
	}

	return out;
}
//--------------------------------------------------------------
} // namespace CompactStar::Physics::Evolution::Run