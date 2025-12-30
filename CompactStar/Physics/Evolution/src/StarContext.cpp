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

	// Adjust name to your actual API.
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
	const auto v = ProfileVersion_();

	// If profile changed since last snapshot, invalidate derived caches.
	if (v != m_cached_version)
	{
		m_rho_gcm3.reset();
		m_cached_version = v;
	}

	// Build on demand
	if (!m_rho_gcm3 && m_eps)
		BuildMassDensityCache_();
}
//--------------------------------------------------------------

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
	// Adjust to your DataColumn API: name/units/metadata.
	m_rho_gcm3 = std::make_unique<Zaki::Vector::DataColumn>("rho(g/cm^3)", rho);
}

//==============================================================
} // namespace CompactStar::Physics::Evolution