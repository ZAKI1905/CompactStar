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

#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"

#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

namespace CompactStar::Physics::Evolution
{

/**
 * @brief Convert StepperType to a stable, human-readable name.
 */
const char *StepperTypeName(StepperType s) noexcept
{
	switch (s)
	{
	case StepperType::RK2:
		return "RK2";
	case StepperType::RKF45:
		return "RKF45";
	case StepperType::RKCK:
		return "RKCK";
	case StepperType::RK8PD:
		return "RK8PD";
	case StepperType::MSBDF:
		return "MSBDF";
	default:
		return "Unknown";
	}
}

/**
 * @brief Convert SaveCadence to a stable, human-readable name.
 */
const char *SaveCadenceName(SaveCadence c) noexcept
{
	switch (c)
	{
	case SaveCadence::LinearDt:
		return "LinearDt";
	case SaveCadence::LogTime:
		return "LogTime";
	default:
		return "Unknown";
	}
}

/**
 * @brief Compute the effective geometric factor q for LogTime sampling.
 */
double Config::EffectiveLogQ() const
{

	if (save_cadence != SaveCadence::LogTime)
		return 0.0;

	if (log_q > 1.0 && std::isfinite(log_q))
		return log_q;

	if (samples_per_decade > 0.0 && std::isfinite(samples_per_decade))
		return std::pow(10.0, 1.0 / samples_per_decade);

	return 0.0;
}

/**
 * @brief Build a compact, unambiguous configuration string for logs.
 */
std::string Config::ToLogString() const
{
	std::ostringstream oss;
	oss.setf(std::ios::scientific);
	oss << std::setprecision(6);

	oss << "stepper=" << StepperTypeName(stepper)
		<< ", rtol=" << rtol
		<< ", atol=" << atol
		<< ", max_samples=" << max_samples
		<< ", max_internal_steps=" << max_internal_steps
		<< ", save_cadence=" << SaveCadenceName(save_cadence);

	if (save_cadence == SaveCadence::LinearDt)
	{
		oss << ", dt_save_s=" << dt_save;
	}
	else
	{
		const double q = EffectiveLogQ();

		oss << ", log_t_floor_s=" << log_t_floor
			<< ", samples_per_decade=" << samples_per_decade
			<< ", q=" << q;

		if (log_q > 1.0 && std::isfinite(log_q))
			oss << " (override)";
	}

	if (!run_label.empty())
		oss << ", run_label=\"" << run_label << "\"";

	return oss.str();
}

} // namespace CompactStar::Physics::Evolution