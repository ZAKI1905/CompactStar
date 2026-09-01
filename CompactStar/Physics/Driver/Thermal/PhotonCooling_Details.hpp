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
 * @file PhotonCooling_Details.hpp
 * @brief Shared computations for PhotonCooling (physics + diagnostics).
 *
 * This file centralizes the computation of:
 *  - surface temperature mapping (currently trivial / placeholder),
 *  - A_inf and A_eff_inf,
 *  - luminosity at infinity,
 *  - physical cooling rate dTinf/dt [K/s],
 *  - ODE RHS contribution d/dt ln(Tinf/Tref) [1/s].
 *
 * The goal is to avoid the “diagnostics drift” problem: if you change the physics
 * implementation, diagnostics automatically match, because they call the same helpers.
 */

#include <cmath>
#include <string>

// #include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
// #include "CompactStar/Physics/State/ThermalState.hpp"

namespace CompactStar::Physics::Evolution
{
class StateVector;
class DriverContext;
} // namespace CompactStar::Physics::Evolution
namespace CompactStar::Physics::Evolution::Diagnostics
{
class DiagnosticPacket;
}

namespace CompactStar::Physics::Driver::Thermal
{
class PhotonCooling;
} // namespace CompactStar::Physics::Driver::Thermal

namespace CompactStar::Physics::Driver::Thermal::Detail
{

/**
 * @brief Bundle of derived photon-cooling quantities for diagnostics/physics.
 *
 * All fields are *computed* outputs. If a field couldn't be computed (missing cache),
 * the helper sets `ok=false` and provides a message; callers decide whether to proceed.
 */
struct PhotonCooling_Details
{
	bool ok = true;
	std::string message;

	// Inputs (resolved)
	double Tinf_K = 0.0;
	double Tsurf_K = 0.0;

	// NEW: Envelope inputs (resolved)
	double Tb_K = 0.0; // local base-of-envelope temperature [K] at rho_b
	double g14 = 0.0;  // surface gravity in units of 1e14 cm s^-2

	// Geometry/redshift
	double R_surf_km = 0.0;
	double R_surf_cm = 0.0;
	double exp2nu_surf = 0.0;

	// Areas
	double A_inf_cm2 = 0.0;
	double A_eff_inf_cm2 = 0.0;

	// Luminosity and RHS
	double L_gamma_inf_erg_s = 0.0;

	/// Canonical GR-integrated stellar heat capacity C_*(T_inf) [erg/K] used as the
	/// denominator of the photon-cooling RHS, per ADR-0002. Obtained from
	/// StarContext::HeatCapacityStar_Tinf — never a driver-local constant.
	double C_star_erg_K = 0.0;

	double dTinf_dt_K_s = 0.0;	 // physical cooling rate [K/s]
	double dLnTinf_dt_1_s = 0.0; // ODE RHS contribution [1/s]
};

/**
 * @brief Compute derived quantities for PhotonCooling.
 *
 * @param drv PhotonCooling instance (options/constants live in drv.opts_).
 * @param Y State vector (ThermalState read here).
 * @param ctx Driver context (GeometryCache read here).
 */
PhotonCooling_Details ComputeDerived(const PhotonCooling &drv,
									 const Evolution::StateVector &Y,
									 const Evolution::DriverContext &ctx);

/**
 * @brief Write a diagnostics packet for PhotonCooling at a given snapshot.
 *
 * This uses ComputeDerived() so diagnostics and physics stay consistent.
 */
void Diagnose(const PhotonCooling &self,
			  double t,
			  const Evolution::StateVector &Y,
			  const Evolution::DriverContext &ctx,
			  Evolution::Diagnostics::DiagnosticPacket &out);

} // namespace CompactStar::Physics::Driver::Thermal::Detail