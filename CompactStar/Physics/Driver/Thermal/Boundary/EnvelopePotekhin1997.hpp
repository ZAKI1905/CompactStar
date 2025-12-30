#ifndef COMPACTSTAR_PHYSICS_DRIVER_THERMAL_BOUNDARY_ENVELOPEPOTEKHIN1997_HPP
#define COMPACTSTAR_PHYSICS_DRIVER_THERMAL_BOUNDARY_ENVELOPEPOTEKHIN1997_HPP

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
 * @file EnvelopePotekhin1997.hpp
 * @brief Potekhin-style Tb→Ts envelope fits (1997-era / classic blanketing envelope models).
 *
 * @ingroup PhysicsDriverThermalBoundary
 *
 * This header declares concrete implementations of the Boundary::IEnvelope interface
 * using Potekhin-style analytic fits. These provide the **thermal boundary condition**
 * relating the local temperature at the base of the blanketing envelope (Tb) to the
 * local effective surface temperature (Ts).
 *
 * ## Definitions and conventions
 * - **Tb**: local temperature at the base of the envelope [K].
 *   In CompactStar, Tb is typically evaluated at a density-defined boundary
 *   (e.g., \f$\rho_b \sim 10^{10}\,\mathrm{g\,cm^{-3}}\f$), chosen via
 *   Boundary::TbDefinition / Boundary::FindTbIndex().
 * - **Ts**: local effective surface temperature [K].
 * - **g14**: surface gravity in units of \f$10^{14}\,\mathrm{cm\,s^{-2}}\f$,
 *   typically computed via Boundary::SurfaceGravity_g14().
 * - **xi**: light-element column / accreted-envelope parameter (dimensionless
 *   “knob” used by accreted-envelope fits). For pure iron envelopes, xi is ignored.
 *
 * ## Scope
 * - These classes are intentionally small and stable; physics details live in the
 *   corresponding .cpp implementation (EnvelopePotekhin1997.cpp).
 * - No star-structure queries occur here; callers supply Tb and g14 explicitly.
 *
 * ## Notes
 * - The Potekhin family of fits appears in several “generations” (1997/1999-era,
 *   later updates, and review-parameterizations). This header is for the classic
 *   1997-style form; keep newer variants in separate headers (e.g., 2003) to avoid
 *   silent physics changes.
 */

#include "CompactStar/Physics/Driver/Thermal/Boundary/IEnvelope.hpp"

namespace CompactStar::Physics::Driver::Thermal::Boundary
{

/**
 * @class EnvelopePotekhin1997_Iron
 * @brief Potekhin-style heavy-element (iron) envelope fit.
 *
 * Returns \f$T_s(T_b, g_{14})\f$ for an iron envelope. The parameter @p xi is ignored.
 *
 * @note Input validation (e.g., Tb>0, g14>0) is performed in the .cpp implementation.
 */
class EnvelopePotekhin1997_Iron final : public IEnvelope
{
  public:
	~EnvelopePotekhin1997_Iron() override = default;

	/**
	 * @brief Map base-of-envelope temperature to surface temperature.
	 * @param Tb   Base-of-envelope local temperature [K].
	 * @param g14  Surface gravity in units of 1e14 cm s^-2.
	 * @param xi   Ignored for iron envelopes (kept for interface compatibility).
	 * @return Local effective surface temperature Ts [K].
	 */
	double Ts_from_Tb(double Tb, double g14, double xi) const override;

	/// Optional human-readable identifier for logs/diagnostics.
	const char *ModelName() override { return "Potekhin1997_Iron"; }
};

/**
 * @class EnvelopePotekhin1997_Accreted
 * @brief Potekhin-style light-element (accreted) envelope fit.
 *
 * Returns \f$T_s(T_b, g_{14}, \xi)\f$ for accreted/light-element envelopes.
 * The parameter @p xi controls the amount of light elements (model-dependent).
 *
 * @note Input validation (e.g., Tb>0, g14>0, finite xi) is performed in the .cpp implementation.
 */
class EnvelopePotekhin1997_Accreted final : public IEnvelope
{
  public:
	~EnvelopePotekhin1997_Accreted() override = default;

	/**
	 * @brief Map base-of-envelope temperature to surface temperature.
	 * @param Tb   Base-of-envelope local temperature [K].
	 * @param g14  Surface gravity in units of 1e14 cm s^-2.
	 * @param xi   Light-element column / accreted-envelope parameter.
	 * @return Local effective surface temperature Ts [K].
	 */
	double Ts_from_Tb(double Tb, double g14, double xi) const override;

	const char *ModelName() override { return "Potekhin1997_Accreted"; }
};

} // namespace CompactStar::Physics::Driver::Thermal::Boundary

#endif /* COMPACTSTAR_PHYSICS_DRIVER_THERMAL_BOUNDARY_ENVELOPEPOTEKHIN1997_HPP */