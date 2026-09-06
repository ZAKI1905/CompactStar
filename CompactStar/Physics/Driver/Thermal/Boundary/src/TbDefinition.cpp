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
 * @file TbDefinition.cpp
 * @brief Helpers for defining the base-of-envelope temperature Tb and locating the envelope boundary by density.
 *
 * This translation unit implements:
 *  - FindTbIndex(): selects the radial index i_b corresponding to a chosen envelope-base density threshold rho_b.
 *  - ComputeTb(): computes the local base-of-envelope temperature Tb from a redshifted isothermal-core temperature T_inf.
 *
 * ## Physical meaning
 * In neutron-star cooling models, the outer heat-blanketing envelope is often treated via a boundary condition
 * at a fiducial mass/energy density rho_b (commonly ~1e10 g/cm^3). The mapping Tb -> Ts is then supplied by an
 * envelope model (e.g., Potekhin fits), while the interior is assumed isothermal in the redshifted sense:
 *
 *   T_inf = exp(nu(r)) * T_local(r)  (constant in the isothermal interior),
 *
 * so that at the envelope base (index i_b):
 *
 *   Tb = T_local(r_b) = T_inf / exp(nu(r_b)).
 *
 * ## Units and conversions
 * CompactStar stores the stellar energy density in geometric units as eps [km^-2].
 * The envelope threshold rho_b is typically specified in cgs mass density units [g/cm^3].
 *
 * Accepted ADR-0012 owns both directions in CompactStar/RelativityUnits.hpp:
 * rho_b and profile epsilon use the same G,c as the canonical TOV solve.
 *
 * ## Algorithmic notes
 * FindTbIndex() assumes the stored energy density decreases monotonically with radius.
 * It returns the first index i where eps(r_i) <= eps_b. Depending on convention, users may prefer
 * the last index above threshold (i-1) or an interpolated boundary radius.
 *
 * ## Error handling
 * These functions throw std::runtime_error if required profile arrays are missing/inconsistent
 * or if the requested threshold is not reached within the profile.
 */

#include "CompactStar/Physics/Driver/Thermal/Boundary/TbDefinition.hpp"

#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"

#include "CompactStar/RelativityUnits.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace CompactStar::Physics::Driver::Thermal::Boundary
{
// -------------------------------------------------------
std::size_t FindTbIndex(const Evolution::StarContext &star, double rho_b)
{

	// Density in km^-2, radius in km
	const auto &rho = star.EnergyDensity();
	const auto &r = star.Radius();

	// Convert rho_b from g/cm^3 to km^-2
	const double rho_b_km2 = CompactStar::RelativityUnits::MassDensityGcm3ToEnergyKmMinus2(rho_b); // g/cm^3 to km^-2

	if (rho->Size() == 0 || r->Size() == 0 || rho->Size() != r->Size())
		throw std::runtime_error("FindTbIndex: missing or inconsistent R/Rho arrays.");

	// // We assume rho decreases outward; find the first index where rho <= rho_b_km2.
	// for (std::size_t i = 0; i < rho->Size(); ++i)
	// {
	// 	if (rho->operator[](i) <= rho_b_km2)
	// 		return i;
	// }

	// Find the outermost index i such that rho[i] >= rho_b_km2.
	// Scan from the surface inward for early exit (blanket is near the surface).
	for (std::size_t i = rho->Size(); i-- > 0;)
	{
		if ((*rho)[i] >= rho_b_km2)
			return i; // inside (deeper) side of the boundary
	}

	throw std::runtime_error("FindTbIndex: rho_b threshold not reached in profile.");

	// Alternatively, since DataColumn has a GetClosestIdx method:
	// return rho->GetClosestIdx(rho_b_km2);

	throw std::runtime_error("FindTbIndex: rho_b threshold not reached in profile.");
}

// -------------------------------------------------------
double ComputeTb(const Evolution::StarContext &star,
				 const Evolution::GeometryCache *geo,
				 double T_inf,
				 const TbDefinition &def)
{
	if (!(T_inf > 0.0))
		throw std::runtime_error("ComputeTb: T_inf must be > 0.");

	const std::size_t i_b = FindTbIndex(star, def.rho_b);

	if (!def.assume_isothermal_redshifted)
	{
		throw std::runtime_error(
			"ComputeTb: assume_isothermal_redshifted=false not implemented; "
			"provide Tb policy or local temperature profile.");
	}

	double expnu_b = 0.0;

	if (def.prefer_geometry_cache && geo && geo->ExpNu().Size() > i_b)
	{
		expnu_b = geo->ExpNu()[i_b];
	}
	else
	{
		expnu_b = exp(star.Nu()->operator[](i_b));
	}

	if (!(expnu_b > 0.0))
		throw std::runtime_error("ComputeTb: exp(nu) at base-of-envelope is invalid.");

	// T_local = T_inf / exp(nu)
	return T_inf / expnu_b;
}
// -------------------------------------------------------
} // namespace CompactStar::Physics::Driver::Thermal::Boundary