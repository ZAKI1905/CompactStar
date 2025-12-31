// -*- lsst-c++ -*-
/**
 * @file StarContext.cpp
 * @brief See StarContext.hpp. Read-only adapter over StarProfile.
 */

#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Core/StarProfile.hpp"

#include <Zaki/Physics/Constants.hpp> // unit conversions
#include <Zaki/Vector/DataColumn.hpp>

#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::Evolution
{
// --------------------
// Helpers
// --------------------
std::uint64_t StarContext::ProfileVersion_() const
{
	if (!m_prof)
		return 0;

	// Examples: m_prof->Version(), m_prof->GetVersion(), m_prof->version().
	return static_cast<std::uint64_t>(m_prof->Version());
}
//==============================================================
//                   StarContext Class
//==============================================================

//--------------------------------------------------------------
StarContext::StarContext(const CompactStar::Core::StarProfile &prof)
	: m_prof(&prof)
{
	BindColumnsOrThrow_();
	ValidateOrThrow_();

	// Initialize cache version snapshot to current (cache is empty until requested).
	m_cached_version = ProfileVersion_();
}

//--------------------------------------------------------------
std::size_t StarContext::Size() const
{
	return (m_r ? m_r->Size() : 0);
}

//--------------------------------------------------------------
double StarContext::RadiusSurface() const
{
	if (!m_r || m_r->Size() == 0)
		return 0.0;

	return (*m_r)[-1]; // km
}

//--------------------------------------------------------------
double StarContext::MassSurface() const
{
	if (!m_m || m_m->Size() == 0)
		return 0.0;

	return (*m_m)[-1]; // km
}

//--------------------------------------------------------------
double StarContext::ExpNuSurface() const
{
	if (!m_nu || m_nu->Size() == 0)
		return 0.0;

	return std::exp((*m_nu)[-1]);
}

//==============================================================
//                   Private helpers
//==============================================================

//--------------------------------------------------------------
void StarContext::BindColumnsOrThrow_()
{
	// Mandatory geometry
	m_r = m_prof->GetRadius();
	m_m = m_prof->GetMass();

	if (!m_r)
		throw std::runtime_error("StarContext: missing radius column r(km)");
	if (!m_m)
		throw std::runtime_error("StarContext: missing mass column m(km)");

	// Optional but strongly expected
	m_nu = m_prof->GetMetricNu();
	m_lam = m_prof->GetMetricLambda();

	// Optional thermodynamics
	m_nb = m_prof->GetBaryonDensity();
	m_pre = m_prof->GetPressure();
	m_eps = m_prof->GetEnergyDensity();
}

//--------------------------------------------------------------
void StarContext::ValidateOrThrow_()
{
	const std::size_t n = m_r->Size();
	if (n == 0)
		throw std::runtime_error("StarContext: profile has zero rows");

	auto check = [&](const Zaki::Vector::DataColumn *c,
					 const char *name)
	{
		if (c && c->Size() != n)
		{
			throw std::runtime_error(
				std::string("StarContext: column '") + name +
				"' has size " + std::to_string(c->Size()) +
				", expected " + std::to_string(n));
		}
	};

	check(m_m, "m");
	check(m_nu, "nu");
	check(m_lam, "lambda");
	check(m_nb, "nB");
	check(m_pre, "p");
	check(m_eps, "eps");
}

//--------------------------------------------------------------
const Zaki::Vector::DataColumn *StarContext::MassDensity_gcm3() const
{
	if (!IsValid())
		return nullptr;

	// Need energy density to build rho
	if (!m_eps)
		return nullptr;

	RefreshDerivedCachesIfNeeded_();

	return m_rho_gcm3.get();
}

//--------------------------------------------------------------
void StarContext::RefreshDerivedCachesIfNeeded_() const
{
	// If the profile is not valid, do not attempt to build derived caches.
	if (!IsValid())
		return;

	const auto v = ProfileVersion_();

	// If profile changed since last snapshot, invalidate derived caches.
	if (v != m_cached_version)
	{
		// Invalidate mass density cache
		m_rho_gcm3.reset();

		// Invalidate direct Urca caches
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;

		// Update cached version
		m_cached_version = v;
	}

	// Build on demand
	if (!m_rho_gcm3 && m_eps)
		BuildMassDensityCache_();

	// Build direct-Urca mask only when needed, not every time.
	// Guard on required inputs
	if (!m_durca_mask)
	{
		if (m_nb && m_prof)
			BuildDirectUrcaMaskCache_();
	}
}

//--------------------------------------------------------------
// Build mass density cache in g/cm^3 from energy density eps (km^-2)
void StarContext::BuildMassDensityCache_() const
{
	// Defensive: ensure profile still valid
	if (!m_prof || !m_eps)
		return;

	const std::size_t n = static_cast<std::size_t>(m_eps->Size());
	if (n == 0)
	{
		m_rho_gcm3.reset();
		return;
	}

	// ---- Unit conversion factor ----
	// eps is stored as (km^-2). Convert to rho [g/cm^3].
	const double kKmMinus2_to_gcm3 = Zaki::Physics::MEV_FM3_2_G_CM3 /
									 Zaki::Physics::MEV_FM3_2_INV_KM2;

	// Build values
	std::vector<double> rho(n);
	for (std::size_t i = 0; i < n; ++i)
	{
		const double eps_km2 = (*m_eps)[i]; // adjust if DataColumn uses operator[]
		rho[i] = eps_km2 * kKmMinus2_to_gcm3;
	}

	// Construct column
	m_rho_gcm3 = std::make_unique<Zaki::Vector::DataColumn>("rho(g/cm^3)", rho);
}

//--------------------------------------------------------------
// Build direct Urca mask cache
void StarContext::BuildDirectUrcaMaskCache_() const
{
	// Defensive: ensure profile still valid
	if (!m_prof || !m_nb || !m_r)
		return;

	const std::size_t n = static_cast<std::size_t>(m_nb->Size());
	if (n == 0)
	{
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;
		return;
	}

	// Species fractions Yi = n_i / nB are stored as "species" columns.
	// Per CompOSE convention:
	//   10: neutron, 11: proton, 0: electron
	const auto *Yn_col = m_prof->GetSpeciesPtr("10");
	const auto *Yp_col = m_prof->GetSpeciesPtr("11");
	const auto *Ye_col = m_prof->GetSpeciesPtr("0");

	if (!Yn_col || !Yp_col || !Ye_col)
	{
		// If any required species fraction is missing, we cannot build the mask.
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;
		return;
	}

	// Sanity check on sizes
	if (static_cast<std::size_t>(Yn_col->Size()) != n ||
		static_cast<std::size_t>(Yp_col->Size()) != n ||
		static_cast<std::size_t>(Ye_col->Size()) != n ||
		static_cast<std::size_t>(m_r->Size()) != n)
	{
		m_durca_mask.reset();
		m_durca_last_allowed = -1;
		m_durca_boundary_r_km = 0.0;
		return;
	}

	// Store as 0/1 bytes.
	std::vector<uint8_t> mask(n, 0);

	// kF = (3*pi^2*n)^(1/3).
	const double three_pi2 = 3.0 * M_PI * M_PI;

	// Compute mask and boundary in ONE reverse scan to avoid a second pass.
	long last_allowed = -1;
	double boundary_r_km = 0.0;

	for (long il = static_cast<long>(n) - 1; il >= 0; --il)
	{
		const std::size_t i = static_cast<std::size_t>(il);

		const double nB = (*m_nb)[i];
		const double Yn = (*Yn_col)[i];
		const double Yp = (*Yp_col)[i];
		const double Ye = (*Ye_col)[i];

		// If nB>0 and Y>=0 then nn/np/ne are automatically >=0.
		if (!(nB > 0.0) || Yn < 0.0 || Yp < 0.0 || Ye < 0.0)
		{
			mask[i] = 0;
			continue;
		}

		const double nn = Yn * nB;
		const double np = Yp * nB;
		const double ne = Ye * nB;

		// If any is exactly 0, cbrt handles it fine and DU is not allowed anyway.
		const double kFn = std::cbrt(three_pi2 * nn);
		const double kFp = std::cbrt(three_pi2 * np);
		const double kFe = std::cbrt(three_pi2 * ne);

		const bool allowed = (kFn <= (kFp + kFe));
		mask[i] = allowed ? 1 : 0;

		// First allowed we encounter in reverse order is the outermost allowed point.
		if (allowed && last_allowed < 0)
		{
			last_allowed = il;
			boundary_r_km = (*m_r)[i];
			// We do NOT break because we still need to fill mask for smaller i.
		}
	}

	m_durca_mask = std::make_unique<std::vector<std::uint8_t>>(std::move(mask));
	m_durca_last_allowed = last_allowed;
	m_durca_boundary_r_km = boundary_r_km;
}

//--------------------------------------------------------------
// Direct Urca kinematic allowance mask (cached).
// mask[i] = 1 if direct Urca is kinematically allowed at radius index i,
// else 0. Computed from Fermi-momentum triangle condition:
//   kFn <= kFp + kFe  with kF = (3*pi^2*n)^(1/3).
// Number densities must be provided in fm^-3 for n_n, n_p, n_e.
// Cache invalidation is based on StarProfile versioning.
const std::vector<std::uint8_t> *StarContext::DirectUrcaMask() const
{
	if (!IsValid())
		return nullptr;

	RefreshDerivedCachesIfNeeded_();

	return m_durca_mask.get();
}

//--------------------------------------------------------------
// Last index (largest r index) where direct Urca is allowed.
// Returns -1 if mask unavailable or no region allows DU.
long StarContext::DirectUrcaLastAllowedIndex() const
{
	if (!IsValid())
		return -1;

	RefreshDerivedCachesIfNeeded_();

	return m_durca_last_allowed;
}

//--------------------------------------------------------------
// Radius (km) at the last allowed index, or 0 if not available.
double StarContext::DirectUrcaBoundaryRadius_km() const
{
	if (!IsValid())
		return 0.0;

	RefreshDerivedCachesIfNeeded_();

	return m_durca_boundary_r_km;
}

//==============================================================
} // namespace CompactStar::Physics::Evolution