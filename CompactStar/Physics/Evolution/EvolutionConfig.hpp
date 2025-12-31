// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file EvolutionConfig.hpp
 * @brief User-configurable options for chemical/thermal/spin evolution runs.
 *
 * Includes integrator tolerances, enabled physics channels, envelope and gap
 * choices, output cadence, and initial-condition sizes (e.g. n_eta).
 *
 * @ingroup PhysicsEvolution
 */

#ifndef CompactStar_Physics_Evolution_Config_H
#define CompactStar_Physics_Evolution_Config_H

#include <cstddef>
#include <string>

namespace CompactStar
{
namespace Physics
{
namespace Evolution
{

//==============================================================
/**
 * @enum StepperType
 * @brief Available ODE steppers (GSL backends).
 *
 * These map directly onto `gsl_odeiv2_step_type` implementations.
 *
 * Typical guidance:
 *  - Use an explicit RK method (RKF45 / RKCK / RK8PD) for non-stiff
 *    or exploratory runs.
 *  - Use MSBDF for stiff late-time thermal/chemical evolution.
 */
enum class StepperType
{
	RKF45, /*!< Runge–Kutta–Fehlberg 4(5) — robust general-purpose non-stiff stepper. */
	RKCK,  /*!< Cash–Karp RK45 — similar to RKF45, sometimes slightly more stable.    */
	RK8PD, /*!< Dormand–Prince 8(5,3) — high-accuracy explicit RK (more expensive).   */
	RK2,   /*!< Simple RK2 (midpoint) — mainly for debugging and sanity checks.       */

	MSBDF /*!< Multistep BDF — stiff solver; good default for late-time evolution.   */
};

//==============================================================
// /**
//  * @enum EnvelopeModel
//  * @brief Surface boundary models mapping \f$T_b \to T_s\f$.
//  */
// enum class EnvelopeModel
// {
// 	Iron,
// 	Accreted,
// 	Custom /*!< Provide a custom callable in Phase 2+. */
// };
//==============================================================

/**
 * @struct Config
 * @brief Configuration for a single evolution run.
 *
 * This struct is intentionally small and POD-like; it is passed by const
 * reference into `EvolutionSystem` and `GSLIntegrator`.
 *
 * ### Integrator-related fields
 *  - `stepper`    : choice of GSL backend (`StepperType`).
 *  - `rtol, atol` : relative/absolute tolerances for adaptive stepping.
 *  - `max_steps`  : hard cap on total number of internal steps.
 *  - `dt_save`    : cadence at which we *request* output samples.
 *
 * The actual integrator is implemented in `GSLIntegrator.cpp` and uses
 * these settings to construct a `gsl_odeiv2_driver`.
 */
struct Config
{
	// ---- Integrator ------------------------------------------------------
	StepperType stepper = StepperType::MSBDF; /*!< Time stepper choice (GSL backend). */
	double rtol = 1e-6;						  /*!< Relative tolerance for adaptive stepping.  */
	double atol = 1e-10;					  /*!< Absolute tolerance for adaptive stepping.  */
	std::size_t max_samples = 1000000;		  /*!< Cap on total number of samples.            */
	std::size_t max_internal_steps = 10000;	  /*!< Cap on internal steps per driver_apply.   */

	// ---- Output ----------------------------------------------------------
	double dt_save = 1.0e5; /*!< Spacing of requested saved samples (s). */
	bool save_intermediate = true;

	// ---- Physics toggles -------------------------------------------------
	bool use_isothermal_core = true;	 /*!< Assume isothermal interior (standard). */
	bool enable_MU = true;				 /*!< Modified Urca emissivity/reactions.    */
	bool enable_DU = true;				 /*!< Direct Urca emissivity/reactions.      */
	bool enable_PBF = false;			 /*!< Pair-breaking/formation emission.      */
	bool enable_BNV = false;			 /*!< Baryon-number violating processes.     */
	bool enable_rotochem_driver = false; /*!< Spin-down driver for \f$\eta\f$.       */
	bool couple_spin = false;			 /*!< Include \f$\Omega\f$ in the state.     */

	// // ---- Envelope and gaps -----------------------------------------------
	// EnvelopeModel envelope = EnvelopeModel::Iron;
	// double envelope_xi = 0.0; /*!< Light-element column parameter (if used). */
	// bool superfluid_n = false;
	// bool superfluid_p = false;

	// ---- Chemical imbalances ---------------------------------------------
	std::size_t n_eta = 1; /*!< Number of \f$\eta_i\f$ components. */

	// ---- Units policy (documentation only; conversions handled elsewhere) -
	std::string unit_policy = "cgs_with_Gc1";

	// ---- Misc ------------------------------------------------------------
	std::string run_label; /*!< Free-form label for outputs. */
};
//==============================================================

} // namespace Evolution
} // namespace Physics
} // namespace CompactStar

#endif /* CompactStar_Physics_Evolution_Config_H */