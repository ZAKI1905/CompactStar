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
// static const char *StepperTypeName(StepperType type)
// {
// 	using ST = StepperType;

// 	switch (type)
// 	{
// 	case ST::RKF45:
// 		return "RKF45";
// 	case ST::RKCK:
// 		return "RKCK";
// 	case ST::RK8PD:
// 		return "RK8PD";
// 	case ST::RK2:
// 		return "RK2";
// 	case ST::MSBDF:
// 		return "MSBDF";
// 	default:
// 		return "UnknownStepper";
// 	}
// }

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
// bool GSLIntegrator::Integrate(double t0, double t1, double *y) const
// {
// 	PROFILE_FUNCTION();

// 	if (!y)
// 	{
// 		throw std::runtime_error("GSLIntegrator::Integrate: y pointer must not be null.");
// 	}

// 	// ---------------------------------------------------------------------
// 	//  Build GSL system description
// 	// ---------------------------------------------------------------------
// 	gsl_odeiv2_system sys;
// 	sys.function = &GslRHS;
// 	sys.jacobian = nullptr; // analytic Jacobian not implemented (let GSL approximate / not use)
// 	sys.dimension = m_dim;
// 	sys.params = const_cast<EvolutionSystem *>(m_sys);

// 	// Select stepper; will throw if configuration is invalid.
// 	const gsl_odeiv2_step_type *step_type = SelectStepper(m_cfg->stepper);

// 	auto compute_log_q = [&]() -> double
// 	{
// 		if (m_cfg->log_q > 1.0 && std::isfinite(m_cfg->log_q))
// 			return m_cfg->log_q;

// 		const double spd = m_cfg->samples_per_decade;
// 		if (!(spd > 0.0) || !std::isfinite(spd))
// 			return 10.0; // fallback: 1 sample per decade

// 		return std::pow(10.0, 1.0 / spd);
// 	};

// 	const double q = (m_cfg->save_cadence == SaveCadence::LogTime) ? compute_log_q() : 0.0;

// 	auto next_target_time = [&](double t_current) -> double
// 	{
// 		if (m_cfg->save_cadence == SaveCadence::LinearDt)
// 		{
// 			const double dt = (m_cfg->dt_save > 0.0) ? m_cfg->dt_save : (t1 - t0);
// 			return std::min(t_current + dt, t1);
// 		}

// 		// LogTime
// 		const double floor_t = (m_cfg->log_t_floor > 0.0) ? m_cfg->log_t_floor : 1.0;

// 		// ensure we're on a positive time for geometric stepping
// 		const double tpos = std::max(t_current, floor_t);

// 		// geometric step
// 		double tnext = tpos * q;

// 		// guard: if q is too close to 1, you can stall due to floating precision
// 		if (!(tnext > t_current))
// 			tnext = t_current + 1.0; // 1 second minimal progress (or something policy-based)

// 		return std::min(tnext, t1);
// 	};

// 	// ---------------------------------------------------------------------
// 	//  Initial step size heuristic
// 	// ---------------------------------------------------------------------
// 	//
// 	// We pick an initial step smaller than dt_save so the driver has room
// 	// to adapt before the first requested output.
// 	double h_init = 0.0;
// 	if (t1 > t0)
// 	{
// 		const double span = t1 - t0;

// 		if (m_cfg->save_cadence == SaveCadence::LogTime)
// 		{
// 			const double floor_t = (m_cfg->log_t_floor > 0.0) ? m_cfg->log_t_floor : 1.0;
// 			const double tpos = std::max(t0, floor_t);
// 			const double dt_first = tpos * (q - 1.0);
// 			h_init = std::min(0.1 * dt_first, 0.01 * span);
// 			if (!(h_init > 0.0) || !std::isfinite(h_init))
// 				h_init = 0.01 * span;
// 		}
// 		else
// 		{
// 			const double dt_save = (m_cfg->dt_save > 0.0) ? m_cfg->dt_save : span;
// 			h_init = std::min(0.1 * dt_save, 0.01 * span);
// 			if (h_init <= 0.0)
// 				h_init = 0.01 * span;
// 		}
// 	}
// 	else
// 		h_init = 1.0;

// 	// ---------------------------------------------------------------------
// 	//  Allocate the GSL driver
// 	// ---------------------------------------------------------------------
// 	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
// 		&sys,
// 		step_type,
// 		h_init,
// 		m_cfg->atol,
// 		m_cfg->rtol);

// 	if (!driver)
// 	{
// 		throw std::runtime_error("GSLIntegrator::Integrate: failed to allocate GSL driver.");
// 	}

// 	// ---------------------------------------------------------------------
// 	//  Log configuration once at the start
// 	// ---------------------------------------------------------------------
// 	Z_LOG_INFO(std::string("Using GSL integrator: ") + m_cfg->ToLogString());

// 	// If there is nothing to do (t0 >= t1), we can return early.
// 	if (t0 >= t1)
// 	{
// 		gsl_odeiv2_driver_free(driver);
// 		return true;
// 	}

// 	// Notify observers at start (t0 snapshot)
// 	m_sys->NotifyStart(t0, t1, y);
// 	// ---------------------------------------------------------------------
// 	//  Main integration loop
// 	// ---------------------------------------------------------------------
// 	double t = t0;
// 	// std::size_t steps_used = 0;

// 	// Treat max_steps as max output samples (since we sample once per dt_save).
// 	std::size_t samples_written = 0;

// 	// dt_save constant for this integration call
// 	// const double dt_save = (m_cfg->dt_save > 0.0) ? m_cfg->dt_save : (t1 - t0);

// 	// Cap internal GSL substeps per driver_apply
// 	if (m_cfg->max_internal_steps > 0)
// 	{
// 		gsl_odeiv2_driver_set_nmax(driver, m_cfg->max_internal_steps);
// 	}

// 	// We advance in chunks of dt_save, letting the GSL driver adapt internally.
// 	while (t < t1)
// 	{
// 		// const double t_target = std::min(t + dt_save, t1);
// 		const double t_target = next_target_time(t);

// 		int status = gsl_odeiv2_driver_apply(driver, &t, t_target, y);

// 		if (status != GSL_SUCCESS)
// 		{
// 			std::ostringstream oss;

// 			oss << "GSLIntegrator: GSL step failed at t=" << t
// 				<< " with status=" << status
// 				<< " (" << gsl_strerror(status) << ")";
// 			Z_LOG_ERROR(oss.str());

// 			gsl_odeiv2_driver_free(driver);

// 			m_sys->NotifyFinish(t, y, false);

// 			return false;
// 		}

// 		// Notify observers once per dt_save chunk (sample cadence)
// 		m_sys->NotifySample(t, y, samples_written);

// 		++samples_written;
// 		if (samples_written > m_cfg->max_samples)
// 		{
// 			std::ostringstream oss;

// 			oss << "GSLIntegrator: exceeded max_samples=" << m_cfg->max_samples
// 				<< " before reaching t1=" << t1
// 				<< " (t=" << t << ")";
// 			Z_LOG_ERROR(oss.str());

// 			gsl_odeiv2_driver_free(driver);

// 			m_sys->NotifyFinish(t, y, false);

// 			return false;
// 		}

// 		// At this point, in a future phase, we can notify Observers with (t,y).
// 	}

// 	gsl_odeiv2_driver_free(driver);

// 	m_sys->NotifyFinish(t, y, true);

// 	return true;
// }
//--------------------------------------------------------------
//  GSLIntegrator::Integrate
//--------------------------------------------------------------
bool GSLIntegrator::Integrate(double t0, double t1, double *y) const
{
	PROFILE_FUNCTION();

	if (!y)
		throw std::runtime_error("GSLIntegrator::Integrate: y pointer must not be null.");

	// ---------------------------------------------------------------------
	//  Build GSL system description
	// ---------------------------------------------------------------------
	gsl_odeiv2_system sys;
	sys.function = &GslRHS;
	sys.jacobian = nullptr; // analytic Jacobian not implemented
	sys.dimension = m_dim;
	sys.params = const_cast<EvolutionSystem *>(m_sys);

	// Select stepper; will throw if configuration is invalid.
	const gsl_odeiv2_step_type *step_type = SelectStepper(m_cfg->stepper);

	// ---------------------------------------------------------------------
	//  Early-out: degenerate interval
	// ---------------------------------------------------------------------
	// Z_LOG_INFO(std::string("Using GSL integrator: ") + m_cfg->ToLogString());

	if (t0 >= t1)
	{
		// Still notify start/finish? Up to you; current behavior returns early.
		return true;
	}

	const double span = t1 - t0;

	// For LogTime, validate q; if invalid, downgrade to LinearDt to avoid stalls/explosions.
	double q = 0.0;
	SaveCadence cadence = m_cfg->save_cadence;

	if (cadence == SaveCadence::LogTime)
	{
		q = m_cfg->EffectiveLogQ();

		if (!(q > 1.0))
		{
			Z_LOG_WARNING("Invalid LogTime cadence; falling back to LinearDt.");
			cadence = SaveCadence::LinearDt;
			q = 0.0;
		}
	}

	const auto effective_dt_save = [&]() -> double
	{
		// Only used for LinearDt mode.
		if (m_cfg->dt_save > 0.0 && std::isfinite(m_cfg->dt_save))
			return m_cfg->dt_save;
		return span;
	};

	const auto log_floor_t = [&]() -> double
	{
		const double f = m_cfg->log_t_floor;
		return (f > 0.0 && std::isfinite(f)) ? f : 1.0;
	};

	// Scale-aware minimal progress (prevents infinite loops when t_target == t due to precision).
	const auto min_progress = [&](double t_current) -> double
	{
		// A small step relative to |t|, but never zero.
		const double s = std::max(1.0, std::fabs(t_current));
		return 64.0 * std::numeric_limits<double>::epsilon() * s;
	};

	auto next_target_time = [&](double t_current) -> double
	{
		if (cadence == SaveCadence::LinearDt)
		{
			const double dt = effective_dt_save();
			return std::min(t_current + dt, t1);
		}

		// LogTime:
		// We want geometric spacing in time *once we're past floor*.
		// If we're below floor, we first jump to floor (or slightly above, if already there).
		const double floor_t = log_floor_t();

		double tnext = 0.0;
		if (t_current < floor_t)
		{
			tnext = floor_t;
		}
		else
		{
			tnext = t_current * q;
		}

		// Ensure strict progress in floating arithmetic.
		const double eps = min_progress(t_current);
		if (!(tnext > t_current + eps))
			tnext = t_current + std::max(eps, 1.0 * std::numeric_limits<double>::epsilon());

		return std::min(tnext, t1);
	};

	// ---------------------------------------------------------------------
	//  Initial step size heuristic
	// ---------------------------------------------------------------------
	double h_init = 0.0;
	{
		// We pick an initial step comfortably smaller than the first requested output interval.
		double dt_first = 0.0;

		if (cadence == SaveCadence::LinearDt)
		{
			dt_first = std::min(effective_dt_save(), span);
		}
		else
		{
			// LogTime: the first target is either floor_t (if t0<floor) or t0*q (if t0>=floor).
			const double floor_t = log_floor_t();
			const double t_first_target = (t0 < floor_t) ? floor_t : (t0 * q);
			dt_first = std::min(std::max(0.0, t_first_target - t0), span);

			// If dt_first collapses (can happen with extreme params), fall back to a small fraction of span.
			if (!(dt_first > 0.0) || !std::isfinite(dt_first))
				dt_first = 0.01 * span;
		}

		// Conservative: 10% of first output interval, but also cap at 1% of total span.
		h_init = std::min(0.1 * dt_first, 0.01 * span);

		if (!(h_init > 0.0) || !std::isfinite(h_init))
			h_init = 0.01 * span;

		// Final guard: never zero.
		if (!(h_init > 0.0))
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
		throw std::runtime_error("GSLIntegrator::Integrate: failed to allocate GSL driver.");

	// Cap internal GSL substeps per driver_apply
	if (m_cfg->max_internal_steps > 0)
		gsl_odeiv2_driver_set_nmax(driver, m_cfg->max_internal_steps);

	// ---------------------------------------------------------------------
	//  Notify observers at start (t0 snapshot)
	// ---------------------------------------------------------------------
	m_sys->NotifyStart(t0, t1, y);

	// ---------------------------------------------------------------------
	//  Main integration loop: advance to each requested output time
	// ---------------------------------------------------------------------
	double t = t0;
	std::size_t samples_written = 0;

	while (t < t1)
	{
		const double t_target = next_target_time(t);

		// Defensive: if we somehow fail to advance, break with error to avoid infinite loops.
		if (!(t_target > t))
		{
			Z_LOG_ERROR("GSLIntegrator: non-advancing target time; aborting to avoid infinite loop.");
			gsl_odeiv2_driver_free(driver);
			m_sys->NotifyFinish(t, y, false);
			return false;
		}

		const int status = gsl_odeiv2_driver_apply(driver, &t, t_target, y);

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

		// Sample (one per requested output time)
		m_sys->NotifySample(t, y, samples_written);

		++samples_written;
		if (samples_written >= m_cfg->max_samples)
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
	}

	gsl_odeiv2_driver_free(driver);
	m_sys->NotifyFinish(t, y, true);
	return true;
}

//--------------------------------------------------------------

} // namespace CompactStar::Physics::Evolution