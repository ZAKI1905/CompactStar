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
 * @file GSLIntegrator.cpp
 * @brief Implementation of GSLIntegrator (GSL driver wrapper).
 */

#include "CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp"

#include <algorithm> // std::min
#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include <Zaki/Util/Instrumentor.hpp> // PROFILE_FUNCTION
#include <Zaki/Util/Logger.hpp>		  // Z_LOG_INFO, Z_LOG_ERROR

#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"

namespace CompactStar::Physics::Evolution
{

//--------------------------------------------------------------
//  Local helpers
//--------------------------------------------------------------

/**
 * @brief Map StepperType → GSL stepper type.
 *
 * Throws std::runtime_error if the enum value is not recognized. This is
 * intentionally strict so invalid configurations fail fast and loudly.
 */
static const gsl_odeiv2_step_type *SelectStepper(StepperType type)
{
	using ST = StepperType;

	switch (type)
	{
	case ST::RKF45:
		return gsl_odeiv2_step_rkf45;
	case ST::RKCK:
		return gsl_odeiv2_step_rkck;
	case ST::RK8PD:
		return gsl_odeiv2_step_rk8pd;
	case ST::RK2:
		return gsl_odeiv2_step_rk2;
	case ST::MSBDF:
		return gsl_odeiv2_step_msbdf;
	default:
		throw std::runtime_error("GSLIntegrator::SelectStepper: unsupported StepperType value.");
	}
}

/**
 * @brief Human-readable name for logging.
 */
static const char *StepperTypeName(StepperType type)
{
	using ST = StepperType;

	switch (type)
	{
	case ST::RKF45:
		return "RKF45";
	case ST::RKCK:
		return "RKCK";
	case ST::RK8PD:
		return "RK8PD";
	case ST::RK2:
		return "RK2";
	case ST::MSBDF:
		return "MSBDF";
	default:
		return "UnknownStepper";
	}
}

//--------------------------------------------------------------
/**
 * @brief GSL-compatible RHS callback forwarding to EvolutionSystem.
 *
 * This function matches the `gsl_odeiv2_system::function` signature and
 * simply forwards to `EvolutionSystem::operator()`.
 */
static int GslRHS(double t, const double y[], double dydt[], void *params)
{
	auto *sys = static_cast<EvolutionSystem *>(params);
	return sys->operator()(t, y, dydt);
}

//--------------------------------------------------------------
//  GSLIntegrator::GSLIntegrator
//--------------------------------------------------------------
GSLIntegrator::GSLIntegrator(const EvolutionSystem &sys,
							 const Config &cfg,
							 std::size_t dim)
	: m_sys(&sys),
	  m_cfg(&cfg),
	  m_dim(dim)
{
	if (!m_sys)
	{
		throw std::runtime_error("GSLIntegrator: EvolutionSystem pointer must not be null.");
	}
	if (!m_cfg)
	{
		throw std::runtime_error("GSLIntegrator: Config pointer must not be null.");
	}
	if (m_dim == 0)
	{
		throw std::runtime_error("GSLIntegrator: dimension must be > 0.");
	}
}

//--------------------------------------------------------------
//  GSLIntegrator::Integrate
//--------------------------------------------------------------
bool GSLIntegrator::Integrate(double t0, double t1, double *y) const
{
	PROFILE_FUNCTION();

	if (!y)
	{
		throw std::runtime_error("GSLIntegrator::Integrate: y pointer must not be null.");
	}

	// ---------------------------------------------------------------------
	//  Build GSL system description
	// ---------------------------------------------------------------------
	gsl_odeiv2_system sys;
	sys.function = &GslRHS;
	sys.jacobian = nullptr; // analytic Jacobian not implemented (let GSL approximate / not use)
	sys.dimension = m_dim;
	sys.params = const_cast<EvolutionSystem *>(m_sys);

	// Select stepper; will throw if configuration is invalid.
	const gsl_odeiv2_step_type *step_type = SelectStepper(m_cfg->stepper);

	// ---------------------------------------------------------------------
	//  Initial step size heuristic
	// ---------------------------------------------------------------------
	//
	// We pick an initial step smaller than dt_save so the driver has room
	// to adapt before the first requested output.
	double h_init = 0.0;
	if (t1 > t0)
	{
		const double span = t1 - t0;
		const double dt_save = (m_cfg->dt_save > 0.0) ? m_cfg->dt_save : span;

		// h_init = 0.1 * dt_save; // start with 10% of dt_save
		// More conservative: also limit to 1% of total span
		h_init = std::min(0.1 * dt_save, 0.01 * span);

		if (h_init <= 0.0)
		{
			h_init = span * 0.01;
		}
	}
	else
	{
		// Degenerate interval; we will exit immediately below, but avoid
		// constructing a driver with zero/negative h.
		h_init = 1.0;
	}

	// ---------------------------------------------------------------------
	//  Allocate the GSL driver
	// ---------------------------------------------------------------------
	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
		&sys,
		step_type,
		h_init,
		m_cfg->atol,
		m_cfg->rtol);

	if (!driver)
	{
		throw std::runtime_error("GSLIntegrator::Integrate: failed to allocate GSL driver.");
	}

	// ---------------------------------------------------------------------
	//  Log configuration once at the start
	// ---------------------------------------------------------------------
	std::ostringstream oss;
	oss << "Using GSL stepper '"
		<< StepperTypeName(m_cfg->stepper)
		<< "' (rtol=" << m_cfg->rtol
		<< ", atol=" << m_cfg->atol
		<< ", max_samples=" << m_cfg->max_samples
		<< ", dt_save=" << m_cfg->dt_save
		<< ")";

	Z_LOG_INFO(oss.str());

	// If there is nothing to do (t0 >= t1), we can return early.
	if (t0 >= t1)
	{
		gsl_odeiv2_driver_free(driver);
		return true;
	}

	// Notify observers at start (t0 snapshot)
	m_sys->NotifyStart(t0, t1, y);
	// ---------------------------------------------------------------------
	//  Main integration loop
	// ---------------------------------------------------------------------
	double t = t0;
	// std::size_t steps_used = 0;

	// Treat max_steps as max output samples (since we sample once per dt_save).
	std::size_t samples_written = 0;

	// dt_save constant for this integration call
	const double dt_save = (m_cfg->dt_save > 0.0) ? m_cfg->dt_save : (t1 - t0);

	// Cap internal GSL substeps per driver_apply
	if (m_cfg->max_internal_steps > 0)
	{
		gsl_odeiv2_driver_set_nmax(driver, m_cfg->max_internal_steps);
	}

	// We advance in chunks of dt_save, letting the GSL driver adapt internally.
	while (t < t1)
	{
		const double t_target = std::min(t + dt_save, t1);

		int status = gsl_odeiv2_driver_apply(driver, &t, t_target, y);

		if (status != GSL_SUCCESS)
		{
			std::ostringstream oss;

			oss << "GSLIntegrator: GSL step failed at t=" << t
				<< " with status=" << status
				<< " (" << gsl_strerror(status) << ")";
			Z_LOG_ERROR(oss.str());

			gsl_odeiv2_driver_free(driver);

			m_sys->NotifyFinish(t, y, false);

			return false;
		}

		// Notify observers once per dt_save chunk (sample cadence)
		m_sys->NotifySample(t, y, samples_written);

		++samples_written;
		if (samples_written > m_cfg->max_samples)
		{
			std::ostringstream oss;

			oss << "GSLIntegrator: exceeded max_samples=" << m_cfg->max_samples
				<< " before reaching t1=" << t1
				<< " (t=" << t << ")";
			Z_LOG_ERROR(oss.str());

			gsl_odeiv2_driver_free(driver);

			m_sys->NotifyFinish(t, y, false);

			return false;
		}

		// At this point, in a future phase, we can notify Observers with (t,y).
	}

	gsl_odeiv2_driver_free(driver);

	m_sys->NotifyFinish(t, y, true);

	return true;
}

//--------------------------------------------------------------

} // namespace CompactStar::Physics::Evolution