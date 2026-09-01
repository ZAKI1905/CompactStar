// -*- lsst-c++ -*-
/**
 * @file GeometryCache.cpp
 * @brief See GeometryCache.hpp.
 *
 * Implementation policy:
 *  - Prefer DataColumn algebra (exp/log/pow, operator*, etc.)
 *  - Only do explicit loops when we must apply a per-row safety clamp
 *    (e.g., deriving Lambda from 1-2m/r).
 */

#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"

#include <cmath>
#include <stdexcept>

namespace CompactStar::Physics::Evolution
{

// -------------------------------------------------------
// GeometryCache::GeometryCache
// -------------------------------------------------------
GeometryCache::GeometryCache(const StarContext &ctx)
{
	// ADR-0003: record provenance before building. Purely additive — no geometry
	// array below is affected.
	m_prov = ctx.Provenance();

	Build_(ctx);
}

// -------------------------------------------------------
// GeometryCache::Size
// -------------------------------------------------------
std::size_t GeometryCache::Size() const
{
	return m_r.Size();
}

//--------------------------------------------------------------
bool GeometryCache::Matches(const StarContext &ctx) const
{
	return m_prov.IsSet() && m_prov == ctx.Provenance();
}


// -------------------------------------------------------
// GeometryCache::R
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::R() const { return m_r; }

// -------------------------------------------------------
// GeometryCache::Mass in km
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::Mass() const { return m_mass; }

// -------------------------------------------------------
// GeometryCache::Area
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::Area() const { return m_area; }

// -------------------------------------------------------
// GeometryCache::ExpNu
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::ExpNu() const { return m_expNu; }

// -------------------------------------------------------
// GeometryCache::ExpMinusNu
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::ExpMinusNu() const { return m_expMinusNu; }

// -------------------------------------------------------
// GeometryCache::Exp2Nu
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::Exp2Nu() const { return m_exp2Nu; }

// -------------------------------------------------------
// GeometryCache::ExpLambda
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::ExpLambda() const { return m_expLam; }

// -------------------------------------------------------
// GeometryCache::ExpMinusLambda
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::ExpMinusLambda() const { return m_expMinusLam; }

// -------------------------------------------------------
// GeometryCache::ExpNuMinusLambda
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::ExpNuMinusLambda() const { return m_expNuMinusLam; }

// -------------------------------------------------------
// GeometryCache::ExpMinusNuMinusLambda
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::ExpMinusNuMinusLambda() const { return m_expMinusNuMinusLam; }

// -------------------------------------------------------
// GeometryCache::WV
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::WV() const { return m_wV; }

// -------------------------------------------------------
// GeometryCache::WVExpNu
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::WVExpNu() const { return m_wVExpNu; }

// -------------------------------------------------------
// GeometryCache::WVExp2Nu
// -------------------------------------------------------
const Zaki::Vector::DataColumn &GeometryCache::WVExp2Nu() const { return m_wVExp2Nu; }

// -------------------------------------------------------
// GeometryCache::DeriveLambdaFromMR_
//
// Only place we explicitly loop: we need denom clamp.
// Everything else is DataColumn algebra.
// -------------------------------------------------------
Zaki::Vector::DataColumn
GeometryCache::DeriveLambdaFromMR_(const Zaki::Vector::DataColumn &r,
								   const Zaki::Vector::DataColumn &m,
								   double eps)
{
	const std::size_t N = std::min(r.Size(), m.Size());
	Zaki::Vector::DataColumn Lambda;
	Lambda.SetLabel("Lambda(derived)");
	Lambda.Reserve(N);

	for (std::size_t i = 0; i < N; ++i)
	{
		const double r_km = r[i];
		const double m_km = m[i];

		double denom = 1.0;
		if (r_km > 0.0)
		{
			denom = 1.0 - 2.0 * m_km / r_km;
			if (denom <= 0.0)
				denom = eps;
		}

		// Lambda = -0.5 * ln(denom)
		Lambda.PushBack(-0.5 * std::log(denom));
	}

	return Lambda;
}

// -------------------------------------------------------
// GeometryCache::Build_
//
// Build with DataColumn expressions:
//   Mass      = m(km)
//   area      = 4*pi * r^2
//   expNu     = exp(nu)
//   exp2Nu    = expNu * expNu
//   expLam    = exp(Lambda)
//   exp(-Lam) = 1 / expLam  (via operator/(double, DataColumn))
//   wV        = area * expLam
//   wV*eNu    = wV * expNu
//   wV*e2Nu   = wV * exp2Nu
// -------------------------------------------------------
void GeometryCache::Build_(const StarContext &ctx)
{
	const auto *r_col = ctx.Radius();
	const auto *nu_col = ctx.Nu();
	const auto *lam_col = ctx.Lambda();
	const auto *m_col = ctx.Mass();

	// -----------------------------
	// Validate inputs
	// -----------------------------
	if (!r_col || r_col->Size() == 0)
		throw std::runtime_error("GeometryCache: missing/empty radius column");

	const std::size_t N = r_col->Size();

	if (!nu_col || nu_col->Size() != N)
		throw std::runtime_error("GeometryCache: missing/invalid nu column (size mismatch)");

	// -----------------------------
	//    Copy primitive columns
	//	 (we cache values locally)
	// -----------------------------
	m_r = *r_col;
	m_r.SetLabel("r(km)");

	m_mass = *m_col;
	m_mass.SetLabel("m(km)");
	// -----------------------------

	// -----------------------------
	// Lambda: use profile if present, else derive from m,r
	// -----------------------------
	Zaki::Vector::DataColumn Lambda;
	if (lam_col && lam_col->Size() == N)
	{
		Lambda = *lam_col;
		Lambda.SetLabel("Lambda");
	}
	else
	{
		if (!m_col || m_col->Size() != N)
			throw std::runtime_error("GeometryCache: Lambda missing and cannot derive it (need m(km) column)");
		Lambda = DeriveLambdaFromMR_(m_r, *m_col, 1e-15);
	}

	// -----------------------------
	// Area and metric exponentials (DataColumn algebra)
	// -----------------------------
	m_area = (4.0 * M_PI) * m_r.pow(2);
	m_area.SetLabel("4*pi*r^2");

	m_expNu = exp(*nu_col);
	m_expNu.SetLabel("exp(nu)");

	m_expMinusNu = 1.0 / m_expNu;
	m_expMinusNu.SetLabel("exp(-nu)");

	m_exp2Nu = m_expNu * m_expNu;
	m_exp2Nu.SetLabel("exp(2*nu)");

	m_expLam = exp(Lambda);
	m_expLam.SetLabel("exp(Lambda)");

	m_expMinusLam = 1.0 / m_expLam;
	m_expMinusLam.SetLabel("exp(-Lambda)");

	// -----------------------------
	// Mixed metric products (common in transport)
	// -----------------------------
	m_expNuMinusLam = m_expNu * m_expMinusLam;
	m_expNuMinusLam.SetLabel("exp(nu - Lambda)");

	m_expMinusNuMinusLam = m_expMinusNu * m_expMinusLam;
	m_expMinusNuMinusLam.SetLabel("exp(-(nu + Lambda))");

	// -----------------------------
	// Canonical weights
	// -----------------------------
	m_wV = m_area * m_expLam;
	m_wV.SetLabel("wV = 4*pi*r^2*exp(Lambda)");

	m_wVExpNu = m_wV * m_expNu;
	m_wVExpNu.SetLabel("wV*exp(nu)");

	m_wVExp2Nu = m_wV * m_exp2Nu;
	m_wVExp2Nu.SetLabel("wV*exp(2nu)");
}
// -------------------------------------------------------
} // namespace CompactStar::Physics::Evolution