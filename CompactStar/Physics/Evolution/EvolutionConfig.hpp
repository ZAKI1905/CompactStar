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
 * @file EvolutionConfig.hpp
 * @brief User-configurable options for chemical/thermal/spin evolution runs.
 *
 * This header defines:
 *  - ODE stepper selection (GSL backends),
 *  - output sampling schedule (linear vs logarithmic),
 *  - integrator tolerances and limits,
 *  - physics-channel toggles,
 *  - small run metadata.
 *
 * The configuration is passed by const reference into the evolution system
 * and integrator. It is intended to remain simple and POD-like.
 *
 * @ingroup PhysicsEvolution
 */

#ifndef CompactStar_Physics_Evolution_Config_H
#define CompactStar_Physics_Evolution_Config_H

#include <cstddef>
#include <cstdint>
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
 * These map onto `gsl_odeiv2_step_type` implementations. Choice typically
 * depends on stiffness:
 *  - Explicit RK methods are suitable for non-stiff problems.
 *  - BDF methods are suitable for stiff problems (common in late-time cooling).
 */
enum class StepperType
{
	/**
	 * @brief Runge–Kutta–Fehlberg 4(5) explicit method (general-purpose non-stiff).
	 */
	RKF45,

	/**
	 * @brief Cash–Karp RK45 explicit method (similar class to RKF45).
	 */
	RKCK,

	/**
	 * @brief Dormand–Prince 8(5,3) explicit method (high-order, higher cost).
	 */
	RK8PD,

	/**
	 * @brief Simple RK2 (midpoint), mainly for debugging / sanity checks.
	 */
	RK2,

	/**
	 * @brief Multistep BDF (stiff solver), typically a good default for stiff evolution.
	 */
	MSBDF
};

/**
 * @brief Convert StepperType to a stable, human-readable name.
 *
 * This function is intended for logging/metadata only.
 *
 * Requirements:
 *  - Must not allocate.
 *  - Must not throw.
 *
 * @param s Stepper type.
 * @return Stable string name (e.g., "MSBDF").
 */
const char *StepperTypeName(StepperType s) noexcept;

//==============================================================
/**
 * @enum SaveCadence
 * @brief Output sampling schedule requested by the integrator.
 *
 * This affects how the integrator chooses successive *target times* `t_target`
 * at which it requests output samples (i.e., calls `NotifySample(...)`).
 *
 * @note This does NOT control the internal GSL adaptive step sizes.
 *       GSL still chooses internal substeps to satisfy error tolerances;
 *       this only controls the externally visible sampling times.
 */
enum class SaveCadence : std::uint8_t
{
	/**
	 * @brief Linear time cadence: targets follow an arithmetic progression.
	 *
	 * Target times are requested as:
	 *   t_{k+1} = t_k + dt_save
	 */
	LinearDt = 0,

	/**
	 * @brief Logarithmic time cadence: targets follow a geometric progression.
	 *
	 * Target times are requested as:
	 *   t_{k+1} = t_k * q
	 *
	 * The geometric factor q is derived from:
	 *  - log_q if provided and valid (>1), otherwise
	 *  - samples_per_decade via q = 10^(1/samples_per_decade).
	 */
	LogTime = 1
};

/**
 * @brief Convert SaveCadence to a stable, human-readable name.
 *
 * This function is intended for logging/metadata only.
 *
 * Requirements:
 *  - Must not allocate.
 *  - Must not throw.
 *
 * @param c Save cadence.
 * @return Stable string name (e.g., "LogTime").
 */
const char *SaveCadenceName(SaveCadence c) noexcept;

//==============================================================
/**
 * @struct Config
 * @brief Configuration for a single evolution run.
 *
 * This struct is intentionally small and passed by const reference into
 * EvolutionSystem and GSLIntegrator. The integrator uses these settings
 * to construct and control a `gsl_odeiv2_driver`.
 */
struct Config
{
	// ---------------------------------------------------------------------
	// Integrator
	// ---------------------------------------------------------------------

	/**
	 * @brief Time stepper choice (GSL backend).
	 */
	// RKF45 — selected from the Phase-2B-1R convergence study on the authenticated
	// 1.6 Msun CMF passive-cooling trajectory (100 yr -> 1 Myr). RKF45, RKCK and RK8PD all
	// complete, are stable under tolerance tightening, and agree to ~3e-7 in T_inf, far
	// inside the 1e-3 gate; RKF45 additionally showed the smallest nominal-to-tightest
	// difference (6.4e-8). Being indistinguishable in accuracy, the lower-complexity
	// general-purpose method is preferred. Evidence: docs/validation/PASSIVE_COOLING_BASELINE.md.
	//
	// MSBDF was the previous default and was UNUSABLE: it is implicit and GSL dereferences
	// sys.jacobian unconditionally, which CompactStar does not supply, so every default run
	// died with SIGSEGV. MSBDF remains in the enum for future Jacobian support and is now
	// rejected with a clear error by GSLIntegrator rather than crashing.
	//
	// This is a non-stiff choice. A future rotochemical/chemical system may need a stiff
	// method and a real Jacobian; that is NOT decided here.
	StepperType stepper = StepperType::RKF45;

	/**
	 * @brief Relative tolerance for adaptive stepping (dimensionless).
	 *
	 * Used by the GSL driver to control local truncation error.
	 */
	double rtol = 1e-6;

	/**
	 * @brief Absolute tolerance for adaptive stepping (units of each state component).
	 *
	 * Used by the GSL driver to control local truncation error.
	 */
	double atol = 1e-10;

	/**
	 * @brief Cap on total number of externally emitted samples.
	 *
	 * This is a safety cap on the number of times the integrator will call
	 * `NotifySample(...)` during a run.
	 */
	std::size_t max_samples = 1000000;

	/**
	 * @brief Cap on internal GSL substeps per driver_apply.
	 *
	 * This limits work per requested output interval to prevent runaway
	 * internal stepping in difficult regions.
	 */
	std::size_t max_internal_steps = 10000;

	// ---------------------------------------------------------------------
	// Output sampling requested by integrator
	// ---------------------------------------------------------------------

	/**
	 * @brief Sampling cadence policy for the integrator output schedule.
	 *
	 * Determines how the integrator selects successive `t_target` values.
	 */
	SaveCadence save_cadence = SaveCadence::LinearDt;

	/**
	 * @brief LinearDt: spacing of requested samples in seconds.
	 *
	 * Used only when save_cadence == SaveCadence::LinearDt.
	 */
	double dt_save = 1.0e6;

	/**
	 * @brief LogTime: minimum positive time used to avoid t=0 issues (seconds).
	 *
	 * When save_cadence == SaveCadence::LogTime, sampling requires a positive
	 * reference time. This value is used as a floor for early scheduling.
	 */
	double log_t_floor = 1.0;

	/**
	 * @brief LogTime: desired number of samples per decade in time.
	 *
	 * Used only when save_cadence == SaveCadence::LogTime and log_q is not set.
	 * The geometric factor is derived as:
	 *   q = 10^(1/samples_per_decade)
	 */
	double samples_per_decade = 50.0;

	/**
	 * @brief LogTime: explicit geometric ratio q (override).
	 *
	 * If > 1 and finite, this value overrides samples_per_decade.
	 * Used only when save_cadence == SaveCadence::LogTime.
	 */
	double log_q = 0.0;

	/**
	 * @brief Whether to save intermediate samples.
	 *
	 * This flag is reserved for higher-level output policies (e.g., whether
	 * observers write intermediate samples or only final snapshots). The
	 * integrator may still call NotifySample based on cadence.
	 */
	bool save_intermediate = true;

	// ---------------------------------------------------------------------
	// Physics toggles
	// ---------------------------------------------------------------------

	/**
	 * @brief Assume an isothermal interior (standard cooling approximation).
	 */
	bool use_isothermal_core = true;

	/**
	 * @brief Enable Modified Urca emissivity / reactions.
	 */
	bool enable_MU = true;

	/**
	 * @brief Enable Direct Urca emissivity / reactions.
	 */
	bool enable_DU = true;

	/**
	 * @brief Enable Pair-Breaking and Formation (PBF) neutrino emission.
	 */
	bool enable_PBF = false;

	/**
	 * @brief Enable baryon-number violating processes (BNV).
	 */
	bool enable_BNV = false;

	/**
	 * @brief Enable a rotochemical / spin-down driver term for chemical imbalance evolution.
	 */
	bool enable_rotochem_driver = false;

	/**
	 * @brief Couple spin (Omega) into the evolved state vector.
	 *
	 * If true, the integrator state includes Omega and evolves it.
	 * If false, Omega may be treated as an external parameter or fixed.
	 */
	bool couple_spin = false;

	// ---------------------------------------------------------------------
	// Chemical imbalances
	// ---------------------------------------------------------------------

	/**
	 * @brief Number of chemical imbalance components (eta_i) included in the state.
	 */
	std::size_t n_eta = 1;

	// ---------------------------------------------------------------------
	// Units policy (documentation-only)
	// ---------------------------------------------------------------------

	/**
	 * @brief Units policy label (documentation only).
	 *
	 * Conversions are handled elsewhere. This string is meant to be printed in
	 * metadata/logs so output files remain interpretable.
	 */
	std::string unit_policy = "cgs_with_Gc1";

	// ---------------------------------------------------------------------
	// Misc
	// ---------------------------------------------------------------------

	/**
	 * @brief Optional free-form label to tag outputs (e.g., parameter scan ID).
	 */
	std::string run_label;

	// ---------------------------------------------------------------------
	// Helper methods (implemented in EvolutionConfig.cpp)
	// ---------------------------------------------------------------------

	/**
	 * @brief Compute the effective geometric factor q for LogTime sampling.
	 *
	 * Behavior:
	 *  - If save_cadence != LogTime, returns 0.0.
	 *  - If log_q is valid (>1 and finite), returns log_q.
	 *  - Else if samples_per_decade is valid (>0 and finite),
	 *    returns 10^(1/samples_per_decade).
	 *  - Otherwise returns 0.0 to signal invalid configuration.
	 *
	 * @return Effective q for LogTime, or 0.0 if not applicable/invalid.
	 */
	double EffectiveLogQ() const;

	/**
	 * @brief Build a compact, unambiguous configuration string for logs.
	 *
	 * The returned string includes cadence-specific parameters so that the
	 * run configuration is fully reconstructible from logs.
	 *
	 * @return Human-readable configuration string.
	 */
	std::string ToLogString() const;
};

//==============================================================

} // namespace Evolution
} // namespace Physics
} // namespace CompactStar

#endif /* CompactStar_Physics_Evolution_Config_H */