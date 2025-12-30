// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file StarContext.hpp
 * @brief Read-only bridge to Core structures (StarProfile) for the evolution module.
 *
 * StarContext centralizes access to geometry, composition, and metric quantities
 * derived from a precomputed neutron-star profile. It does not run TOV/rotation;
 * it only exposes cached, interpolation-friendly views needed during RHS evaluations.
 *
 * @ingroup PhysicsEvolution
 */
#ifndef CompactStar_Physics_Evolution_StarContext_H
#define CompactStar_Physics_Evolution_StarContext_H

// #include <Zaki/Vector/DataSet.hpp>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace Zaki::Vector
{
class DataColumn; // dependencies/include/Zaki/Vector/DataColumn.hpp
}
namespace CompactStar::Core
{
class StarProfile; // Core/StarProfile.hpp
}

namespace CompactStar::Physics::Evolution
{

//==============================================================
//                      StarContext Class
//==============================================================
/**
 * @class StarContext
 * @brief Immutable, per-star adapter exposing cached geometry and composition.
 *
 * Constructed from an existing @c CompactStar::Core::StarProfile.
 * Caches non-owning pointers to frequently accessed columns (r, m, nu, lambda, nb).
 * Validates consistency (row counts) up front to fail fast before evolution.
 */
class StarContext
{
  public:
	/// Default constructor (invalid context).
	StarContext() = default;

	/**
	 * @brief Construct from a precomputed neutron-star profile.
	 * @param prof Reference to an initialized @c CompactStar::Core::StarProfile.
	 */
	explicit StarContext(const CompactStar::Core::StarProfile &prof);

	/// True iff bound to a profile and required columns were found.
	bool IsValid() const { return m_prof != nullptr && m_r != nullptr && m_m != nullptr; }

	// --------------------
	// Basic geometry / grid
	// --------------------
	std::size_t Size() const;										 ///< Number of radial samples.
	const Zaki::Vector::DataColumn *Radius() const { return m_r; }	 ///< r(km)
	const Zaki::Vector::DataColumn *Mass() const { return m_m; }	 ///< m(km)
	const Zaki::Vector::DataColumn *Nu() const { return m_nu; }		 ///< nu
	const Zaki::Vector::DataColumn *Lambda() const { return m_lam; } ///< lambda

	// --------------------
	// Thermodynamic background quantities
	// --------------------
	const Zaki::Vector::DataColumn *BaryonDensity() const { return m_nb; }	///< nB(fm^-3) or nullptr
	const Zaki::Vector::DataColumn *Pressure() const { return m_pre; }		///< p (km^-2)
	const Zaki::Vector::DataColumn *EnergyDensity() const { return m_eps; } ///< eps (km^-2)

	// --------------------
	// Derived cached columns
	// --------------------
	/**
	 * @brief Mass density in g/cm^3 (cached derived column).
	 *
	 * This is computed from the profile's energy density (km^-2) using:
	 *   rho = eps / c^2  (in geometrized units with G=c=1, eps has length^-2)
	 * and converted to cgs: g/cm^3.
	 *
	 * Cache invalidation is based on StarProfile versioning.
	 *
	 * @return Pointer to cached column, or nullptr if energy density is unavailable.
	 */
	const Zaki::Vector::DataColumn *MassDensity_gcm3() const;

	// --------------------
	// Global scalars (derived from columns)
	// --------------------
	double RadiusSurface() const; ///< r[-1] (km)
	double MassSurface() const;	  ///< m[-1] (km)
	double ExpNuSurface() const;  ///< exp(nu[-1]) if nu exists, else 0

  private:
	void BindColumnsOrThrow_(); // sets m_r/m_m/m_nu/m_lam/m_nb/m_pre/m_eps
	void ValidateOrThrow_();	// checks consistent row counts for required cols

	// Derived cache helpers
	void RefreshDerivedCachesIfNeeded_() const;
	void BuildMassDensityCache_() const;

	// Convenience: current version from profile
	std::uint64_t ProfileVersion_() const;

  private:
	const CompactStar::Core::StarProfile *m_prof = nullptr;

	// Non-owning cached pointers to frequently accessed columns (set in ctor).

	const Zaki::Vector::DataColumn *m_r = nullptr;	 ///< r(km)
	const Zaki::Vector::DataColumn *m_m = nullptr;	 ///< m(km)
	const Zaki::Vector::DataColumn *m_nu = nullptr;	 ///< nu
	const Zaki::Vector::DataColumn *m_lam = nullptr; ///< lambda
	const Zaki::Vector::DataColumn *m_nb = nullptr;	 ///< nB(fm^-3)
	const Zaki::Vector::DataColumn *m_pre = nullptr; ///< p (km^-2)
	const Zaki::Vector::DataColumn *m_eps = nullptr; ///< eps (km^-2)

	// -------- derived caches (owned by StarContext) --------
	mutable std::uint64_t m_cached_version = 0;

	mutable std::unique_ptr<Zaki::Vector::DataColumn> m_rho_gcm3; // cached derived column
};

} // namespace CompactStar::Physics::Evolution

#endif /* CompactStar_Physics_Evolution_StarContext_H */