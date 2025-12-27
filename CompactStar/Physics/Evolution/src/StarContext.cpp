// -*- lsst-c++ -*-
/**
 * @file StarContext.cpp
 * @brief See StarContext.hpp. Read-only adapter over StarProfile.
 */

#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Core/StarProfile.hpp"

#include <Zaki/Vector/DataColumn.hpp>

#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::Evolution
{
//==============================================================
//                   StarContext Class
//==============================================================

//--------------------------------------------------------------
StarContext::StarContext(const CompactStar::Core::StarProfile &prof)
	: m_prof(&prof)
{
	BindColumnsOrThrow_();
	ValidateOrThrow_();
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

//==============================================================
} // namespace CompactStar::Physics::Evolution