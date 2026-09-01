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
#include "CompactStar/Physics/Evolution/ProfileProvenance.hpp"
#include <memory>
#include <vector>

namespace Zaki::Vector
{
class DataColumn; // dependencies/include/Zaki/Vector/DataColumn.hpp
}
namespace CompactStar::Core
{
class StarProfile; // Core/StarProfile.hpp
}

namespace CompactStar
{
namespace EOS
{
class CompOSE_Thermo;
}
} // namespace CompactStar

namespace CompactStar::Physics::Evolution
{

class GeometryCache;

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
	// Heat capacity cache (GR, for cooling in terms of T_infty)
	// --------------------

	// --------------------
	// Heat capacity cache (GR, for cooling in terms of T_infty)
	// --------------------

	/**
	 * @brief Star-integrated heat capacity C(T_infty) for an isothermal (Tolman) core.
	 *
	 * Implements:
	 *   C(Tinf) = ∫ cV(Tlocal, nB, Yq) dV
	 * with:
	 *   Tlocal(r) = Tinf * exp(-nu(r))
	 *   dV = 4*pi*r^2*exp(Lambda(r)) dr   (via GeometryCache::WV()).
	 *
	 * CV is in erg K^{-1}.
	 *
	 * The cache is built on-demand and invalidated when the profile version changes.
	 *
	 * @param Tinf_MeV Redshifted temperature at infinity (MeV).
	 * @param thermo   CompOSE thermodynamics table interface (used for cV).
	 * @param geo	   Pointer to GeometryCache
	 */
	double HeatCapacityStar_Tinf(double Tinf_MeV,
								 const CompactStar::EOS::CompOSE_Thermo &thermo,
								 const GeometryCache *geo = nullptr) const;

	// --------------------
	// Derived cached masks
	// --------------------
	/**
	 * @brief Direct Urca kinematic allowance mask (cached).
	 *
	 * mask[i] = 1 if direct Urca is kinematically allowed at radius index i,
	 * else 0. Computed from Fermi-momentum triangle condition:
	 *   kFn <= kFp + kFe  with kF = (3*pi^2*n)^(1/3).
	 *
	 * Number densities must be provided in fm^-3 for n_n, n_p, n_e.
	 *
	 * Cache invalidation is based on StarProfile versioning.
	 *
	 * @return Pointer to cached mask, or nullptr if required density columns are unavailable.
	 */
	const std::vector<std::uint8_t> *DirectUrcaMask() const;

	/**
	 * @brief Last index (largest r index) where direct Urca is allowed.
	 *
	 * Returns -1 if mask unavailable or no region allows DU.
	 */
	long DirectUrcaLastAllowedIndex() const;

	/**
	 * @brief Radius (km) at the last allowed index, or 0 if not available.
	 */
	double DirectUrcaBoundaryRadius_km() const;

	// --------------------
	// Global scalars (derived from columns)
	// --------------------
	double RadiusSurface() const; ///< r[-1] (km)
	double MassSurface() const;	  ///< m[-1] (km)
	double ExpNuSurface() const;  ///< exp(nu[-1]) if nu exists, else 0

	// Convenience: current version from profile
	std::uint64_t ProfileVersion() const;

	/**
	 * @brief Current provenance of the bound profile: (identity, version).
	 *
	 * ADR-0003. A derived snapshot records this at construction; a consumer compares its
	 * recorded token against this to decide whether the snapshot still applies.
	 */
	[[nodiscard]] ProfileProvenance Provenance() const;

  private:
	void BindColumnsOrThrow_() const; // sets m_r/m_m/m_nu/m_lam/m_nb/m_pre/m_eps
	void ValidateOrThrow_();	// checks consistent row counts for required cols

	// Derived cache helpers
	void RefreshDerivedCachesIfNeeded_() const;
	void BuildMassDensityCache_() const;

	// Direct Urca cache builder
	void BuildDirectUrcaMaskCache_() const;

	// Derived cache helpers (continued)
	void BuildYqCache_() const;
	const Zaki::Vector::DataColumn *ChargeFractionYq() const;

	// Builds the cache in cgs units: erg cm^-3 K^-1
	void BuildHeatCapacityCache_(const CompactStar::EOS::CompOSE_Thermo &thermo,
								 const CompactStar::Physics::Evolution::GeometryCache *geo) const;

  private:
	const CompactStar::Core::StarProfile *m_prof = nullptr;

	// Non-owning cached pointers to frequently accessed columns (set in ctor).

	// ADR-0003 (S1): these views are RE-BOUND when the bound profile's revision changes,
	// so they are mutable — the re-bind happens inside the const refresh path, before any
	// cached view or payload may be used.
	mutable const Zaki::Vector::DataColumn *m_r = nullptr;	 ///< r(km)
	mutable const Zaki::Vector::DataColumn *m_m = nullptr;	 ///< m(km)
	mutable const Zaki::Vector::DataColumn *m_nu = nullptr;	 ///< nu
	mutable const Zaki::Vector::DataColumn *m_lam = nullptr; ///< lambda
	mutable const Zaki::Vector::DataColumn *m_nb = nullptr;	 ///< nB(fm^-3)
	mutable const Zaki::Vector::DataColumn *m_pre = nullptr; ///< p (km^-2)
	mutable const Zaki::Vector::DataColumn *m_eps = nullptr; ///< eps (km^-2)

	// -------- derived caches (owned by StarContext) --------
	mutable std::uint64_t m_cached_version = 0;

	mutable std::unique_ptr<Zaki::Vector::DataColumn> m_rho_gcm3; // cached derived column

	mutable std::unique_ptr<std::vector<std::uint8_t>> m_durca_mask; // 0/1 mask along r
	mutable long m_durca_last_allowed = -1;							 // cached boundary index (or -1)
	mutable double m_durca_boundary_r_km = 0.0;						 // cached boundary radius (km)

	// -------- Heat capacity cache (owned by StarContext) --------
	mutable std::unique_ptr<Zaki::Vector::DataColumn> m_Yq_cache; // cached strong-sector charge fraction Yq(r)

	// Heat capacity cache (star-integrated C(Tinf) for isothermal core)
	struct HeatCapacityCache
	{
		bool loaded = false;

		std::uint64_t prof_version = 0;	  // profile version snapshot
		const void *thermo_tag = nullptr; // pointer identity of thermo used to build this cache

		// ADR-0003: the GeometryCache is a genuine input to BuildHeatCapacityCache_, so it
		// belongs in the validity condition. Keyed on the geometry's PROVENANCE, not on the
		// address of the GeometryCache object: an equivalent geometry rebuilt for the same
		// (profile, version) is interchangeable and must not force a rebuild (ADR-0003 §11).
		ProfileProvenance geo_prov{};

		std::vector<double> Tinf_MeV; // grid (MeV)
		std::vector<double> C_star;	  // integrated heat capacity on the grid

		mutable std::size_t last_i = 0; // accel for interpolation
	};

	mutable HeatCapacityCache m_cv_cache;
};

} // namespace CompactStar::Physics::Evolution

#endif /* CompactStar_Physics_Evolution_StarContext_H */