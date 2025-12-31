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
 * @file EvolutionSystem.cpp
 * @brief Implementation of EvolutionSystem (RHS functor for dY/dt).
 */

#include <stdexcept>

#include "CompactStar/Physics/Driver/IDriver.hpp"
#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/Observers/IObserver.hpp"
#include "CompactStar/Physics/Evolution/RHSAccumulator.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateLayout.hpp"
#include "CompactStar/Physics/Evolution/StatePacking.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"

#include "Zaki/Util/Logger.hpp"

#include "CompactStar/Physics/State/SpinState.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"
namespace CompactStar::Physics::Evolution
{

//--------------------------------------------------------------
// EvolutionSystem::EvolutionSystem
//--------------------------------------------------------------
EvolutionSystem::EvolutionSystem(const DriverContext &ctx,
								 StateVector &state,
								 RHSAccumulator &rhs,
								 const StateLayout &layout,
								 std::vector<DriverPtr> drivers)
	: m_ctx(ctx),
	  m_state(state),
	  m_rhs(rhs),
	  m_layout(layout),
	  m_drivers(std::move(drivers))
{
	// Sanity-check the static context up front so misconfigurations
	// fail early and loudly.
	ValidateContext();

	// Ensure we have at least one physics driver.
	if (m_drivers.empty())
	{
		throw std::runtime_error(
			"EvolutionSystem: constructed with no physics drivers. "
			"At least one IDriver must be provided.");
	}
}

//--------------------------------------------------------------
// EvolutionSystem::ValidateContext
//--------------------------------------------------------------
void EvolutionSystem::ValidateContext() const
{
	// Check that all required context pointers are non-null.
	if (!m_ctx.star)
	{
		throw std::runtime_error("EvolutionSystem: Context.star is null.");
	}
	if (!m_ctx.geo)
	{
		throw std::runtime_error("EvolutionSystem: Context.geo is null.");
	}

	if (!m_ctx.cfg)
	{
		throw std::runtime_error("EvolutionSystem: Context.cfg is null.");
	}

	// NOTE: m_ctx.envelope is optional and must be validated by drivers
	// that require it (based on their Options/surface model).

	// Additional invariants could be checked here, e.g.:
	//  - that StateLayout::TotalSize() matches what the integrator expects,
	//  - consistency between Config and State/StateLayout sizes.
}

//--------------------------------------------------------------
// EvolutionSystem::operator()
//--------------------------------------------------------------
int EvolutionSystem::operator()(double t, const double y[], double dydt[]) const
{
	// 1) Unpack flat y[] into logical State blocks.
	//
	// This uses the free helper defined in StatePacking.hpp:
	//   UnpackStateVector(StateVector&, const StateLayout&, const double*)
	UnpackStateVector(m_state, m_layout, y);

	{
		const auto &th = m_state.GetThermal();
		const auto &sp = m_state.GetSpin();

		// Z_LOG_INFO("DEBUG t=" + std::to_string(t) +
		// 		   "  y[0]=" + std::to_string(y[0]) +
		// 		   "  y[1]=" + std::to_string(y[1]) +
		// 		   "  unpacked: Tinf=" + std::to_string(th.Tinf()) +
		// 		   "  Omega=" + std::to_string(sp.Omega()));
	}

	// 2) Clear RHS accumulator before driver contributions.
	m_rhs.Clear();

	// 3) Let each physics driver accumulate its contribution to dY/dt.
	//
	// Drivers are responsible for:
	//   - reading whatever parts of m_state they depend on,
	//   - using the StarContext (geometry, EOS, etc.) as needed,
	//   - adding their contributions into m_rhs via AddTo(...).
	// const StarContext &starCtx = *m_ctx.star;

	for (const auto &drv : m_drivers)
	{
		if (!drv)
		{
			throw std::runtime_error(
				"EvolutionSystem::operator(): encountered null driver pointer.");
		}

		drv->AccumulateRHS(t, m_state, m_rhs, m_ctx);
	}

	{
		// Temporarily add a getter if you don’t have one:
		// double RHSAccumulator::Peek(StateTag tag, size_t i) const;

		// Z_LOG_INFO("DEBUG rhs: Thermal[0]=" + std::to_string(m_rhs.Peek(State::StateTag::Thermal, 0)) +
		// 		   "  Spin[0]=" + std::to_string(m_rhs.Peek(State::StateTag::Spin, 0)));
	}
	// 4) Scatter accumulated RHS blocks back into the flat dydt[] array.
	//
	// This uses the helper from StatePacking.hpp:
	//   ScatterRHSFromAccumulator(const RHSAccumulator&, const StateLayout&, double*)
	ScatterRHSFromAccumulator(m_rhs, m_layout, dydt);

	// Z_LOG_INFO("DEBUG dydt[0]=" + std::to_string(dydt[0]) +
	// 		   "  dydt[1]=" + std::to_string(dydt[1]));
	// GSL convention: return 0 for success.
	return 0;
}

//--------------------------------------------------------------
// EvolutionSystem::AddObserver
//--------------------------------------------------------------
void EvolutionSystem::AddObserver(ObserverPtr obs)
{
	if (!obs)
	{
		throw std::runtime_error("EvolutionSystem::AddObserver: null observer pointer.");
	}
	m_observers.push_back(std::move(obs));
}

//--------------------------------------------------------------
// EvolutionSystem::NotifyStart
//--------------------------------------------------------------
void EvolutionSystem::NotifyStart(double t0, double t1, const double *y0) const
{
	if (m_observers.empty())
		return;

	// Unpack y0 into logical StateVector so observers see consistent Y.
	UnpackStateVector(m_state, m_layout, y0);

	Observers::RunInfo run;
	run.t0 = t0;
	run.tf = t1;

	for (const auto &obs : m_observers)
	{
		if (obs)
			obs->OnStart(run, m_state, m_ctx);
	}
}

//--------------------------------------------------------------
// EvolutionSystem::NotifySample
//--------------------------------------------------------------
void EvolutionSystem::NotifySample(double t, const double *y, std::size_t sample_index) const
{
	if (m_observers.empty())
		return;

	UnpackStateVector(m_state, m_layout, y);

	Observers::SampleInfo s;
	s.t = t;
	s.sample_index = sample_index;
	s.step_index = sample_index; // current integrator emits one sample per chunk

	for (const auto &obs : m_observers)
	{
		if (obs)
			obs->OnSample(s, m_state, m_ctx);
	}
}

//--------------------------------------------------------------
// EvolutionSystem::NotifyFinish
//--------------------------------------------------------------
void EvolutionSystem::NotifyFinish(double t, const double *y, bool ok) const
{
	if (m_observers.empty())
		return;

	UnpackStateVector(m_state, m_layout, y);

	Observers::FinishInfo fin;
	fin.t_final = t;
	fin.ok = ok;

	for (const auto &obs : m_observers)
	{
		if (obs)
			obs->OnFinish(fin, m_state, m_ctx);
	}
}
//--------------------------------------------------------------
} // namespace CompactStar::Physics::Evolution