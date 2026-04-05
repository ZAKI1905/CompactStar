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
 * @file RotochemicalCache.hpp
 * @brief Precomputed rotochemical heating coefficients dN_i/dOmega^2.
 *
 * Implements the Fernandez & Reisenegger (2005) formalism:
 *
 *   dN_i/dOmega^2 = A_i - B_i * (A_B / B_B)
 *
 * where:
 *   A_i = (partial N_i / partial Omega^2)|_{eps_c}  (structural, from Hartle)
 *   B_i = (partial N_i / partial eps_c)|_Omega       (sequence, finite diff)
 *   A_B, B_B = same for total baryon number
 *
 * The baryon conservation constraint dN_B/dOmega^2 = 0 determines
 * d(eps_c)/d(Omega^2), eliminating one degree of freedom.
 *
 * @ingroup PhysicsEvolution
 */

#ifndef CompactStar_Physics_Evolution_RotochemicalCache_H
#define CompactStar_Physics_Evolution_RotochemicalCache_H

#include <string>
#include <vector>

#include <Zaki/Vector/DataColumn.hpp>

namespace CompactStar
{

// Forward declarations
namespace Core
{
struct HartleResult;
class NStar;
class MixedStar;
} // namespace Core

namespace Physics::Evolution
{

class StarContext;
class GeometryCache;

/**
 * @struct RotochemicalSpeciesData
 * @brief Per-species data for the rotochemical calculation.
 */
struct RotochemicalSpeciesData
{
	std::string label; ///< Species name ("n", "p", "e", "mu")

	double N_total = 0.0; ///< Total enclosed number (background TOV)

	/// A_i = (partial N_i / partial Omega^2)|_{eps_c}
	/// Computed from volume integral with Hartle O(Omega^2) perturbations
	double A = 0.0;

	/// B_i = (partial N_i / partial eps_c)|_Omega
	/// Computed from finite difference of equilibrium stars
	double B = 0.0;

	/// Z_i = dN_i/dOmega^2 = A_i - B_i * (A_B / B_B)
	/// The final reduced coefficient used by the rotochemical driver
	double Z = 0.0;
};

/**
 * @class RotochemicalCache
 * @brief Stores precomputed rotochemical heating coefficients.
 *
 * This cache is populated once for a given stellar configuration
 * and Hartle solution, then used by the Rotochemical IDriver
 * during time evolution.
 */
class RotochemicalCache
{
  public:
	RotochemicalCache() = default;

	// ------------------------------------------------------------------
	// Species data
	// ------------------------------------------------------------------
	std::vector<RotochemicalSpeciesData> species;

	// ------------------------------------------------------------------
	// Total baryon number derivatives
	// ------------------------------------------------------------------
	double A_B = 0.0; ///< (partial N_B / partial Omega^2)|_{eps_c}
	double B_B = 0.0; ///< (partial N_B / partial eps_c)|_Omega

	// ------------------------------------------------------------------
	// Validity
	// ------------------------------------------------------------------
	bool valid = false;

	// ------------------------------------------------------------------
	// Computation methods
	// ------------------------------------------------------------------

	/**
	 * @brief Compute enclosed particle number N_i for a single species.
	 *
	 * N_i = integral_0^R  n_i(r) * 4*pi*r^2 * e^{Lambda(r)} dr
	 *
	 * @param n_i Species number density column (fm^-3), on the radial grid
	 * @param geo GeometryCache providing WV() = 4*pi*r^2*e^Lambda
	 * @return Total enclosed number (in km^-3 * km^3 = dimensionless after FM3_TO_KM3)
	 */
	static double ComputeEnclosedNumber(
		const Zaki::Vector::DataColumn &n_i,
		const Zaki::Vector::DataColumn &r_col,
		const GeometryCache &geo);

	/**
	 * @brief Compute A_i = (partial N_i / partial Omega^2)|_{eps_c}
	 * for a single species.
	 *
	 * This involves a volume integral of the perturbation kernel:
	 *   A_i ~ integral [ (dn_i/dp) * p0 + n_i * (metric corrections) ]
	 *          * proper_volume * dr / Omega^2
	 *
	 * @param n_i Species density column (fm^-3)
	 * @param p_col Pressure column
	 * @param hartle HartleResult with xi0, m0, p0 profiles
	 * @param geo GeometryCache
	 * @return A_i (the structural derivative)
	 */
	static double ComputeStructuralDerivative(
		const Zaki::Vector::DataColumn &n_i,
		const Zaki::Vector::DataColumn &r_col,
		const Zaki::Vector::DataColumn &p_col,
		const Core::HartleResult &hartle,
		const GeometryCache &geo);

	/**
	 * @brief Compute all coefficients for an NStar.
	 *
	 * Steps:
	 *  1. For each species: compute N_i, A_i
	 *  2. Compute A_B (total baryons)
	 *  3. Compute B_i, B_B from finite difference (requires perturbed stars)
	 *  4. Combine: Z_i = A_i - B_i * (A_B / B_B)
	 *
	 * @param star The NStar (provides profile, species, EOS)
	 * @param hartle Second-order Hartle results
	 * @param geo GeometryCache
	 * @param star_plus NStar at eps_c + delta (for B_i computation)
	 * @param star_minus NStar at eps_c - delta (for B_i computation)
	 * @param delta_eps_c Central energy density step for finite difference
	 */
	void Build(const Core::NStar &star,
			   const Core::HartleResult &hartle,
			   const GeometryCache &geo,
			   const Core::NStar &star_plus,
			   const Core::NStar &star_minus,
			   double delta_eps_c);
};

} // namespace Physics::Evolution
} // namespace CompactStar

#endif /* CompactStar_Physics_Evolution_RotochemicalCache_H */
