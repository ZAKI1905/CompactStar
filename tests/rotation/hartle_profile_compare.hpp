// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_profile_compare.hpp
 * @brief Phase 4B — shared fixtures and profile-comparison utilities for the independent
 *        validation of the physically normalized first-order Hartle response.
 *
 * This header holds ONLY test-side plumbing: the exact constant-density fixture, the
 * background adapter, and the node-by-node comparison metrics. **The independent solver itself
 * lives in `hartle_reference.hpp` and is not duplicated here.**
 *
 * THE OBJECT UNDER TEST. Phase 4A exposed the seed-free first-order response
 *
 *     s(r)  = omega_bar(r) / Omega            [dimensionless]
 *     s'(r) = [d omega_bar/dr](r) / Omega     [km^-1]
 *
 * from which the frame-dragging fraction follows as `omega(r)/Omega = 1 - s(r)`. Phase 4A
 * proved the *contract* around these (the seed does not leak, a requested spin comes back,
 * the units are right). Phase 4B asks the different and harder question: **is the shape
 * itself right?** — which only an independently derived profile can answer.
 *
 * HOW s' IS COMPARED. `s'` vanishes at the centre by regularity, so a pointwise *relative*
 * error is meaningless there: both sides are zero and the ratio is 0/0. Two honest measures
 * are reported instead, and both are asserted:
 *
 *   - `rel_to_scale_sp` — the largest absolute difference divided by the profile's own peak
 *     |s'|, i.e. an error relative to the scale of the quantity;
 *   - `max_rel_sp` — the largest pointwise relative difference restricted to nodes where
 *     |s'_ref| is at least 1 % of that peak, so the ratio is taken only where it means
 *     something.
 *
 * `s` is O(1) everywhere on a physical star, so it is compared by ordinary relative error.
 */

#ifndef CompactStar_Tests_HartleProfileCompare_H
#define CompactStar_Tests_HartleProfileCompare_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/TOVSolver.hpp" // TOVPoint is only forward-declared in NStar.hpp
#include "hartle_reference.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace hartle_4b
{

using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;

// ---------------------------------------------------------------------------
//  Exact Schwarzschild constant-density interior (Schwarzschild 1916; MTW §23.7).
//  The same fixture Phase 2B-4B used, restated so the validated harnesses stay untouched.
// ---------------------------------------------------------------------------
struct UniformStar
{
	double M_km = 0, R_km = 0, rho0 = 0; // rho0 in km^-2
	double yR() const { return std::sqrt(1.0 - 2.0 * M_km / R_km); }
	double y(double r) const { return std::sqrt(1.0 - 2.0 * M_km * r * r / (R_km * R_km * R_km)); }
	double m(double r) const { return M_km * r * r * r / (R_km * R_km * R_km); }
	double lambda(double r) const { return -std::log(y(r)); }
	double nu(double r) const { return std::log(1.5 * yR() - 0.5 * y(r)); }
	double p(double r) const { return rho0 * (y(r) - yR()) / (3.0 * yR() - y(r)); }
	double nuprime(double r) const
	{
		return 2.0 * M_km * r / (R_km * R_km * R_km * y(r) * (3.0 * yR() - y(r)));
	}
};

inline UniformStar MakeUniform(double M_km, double R_km)
{
	UniformStar u;
	u.M_km = M_km;
	u.R_km = R_km;
	u.rho0 = 3.0 * M_km / (4.0 * M_PI * R_km * R_km * R_km);
	return u;
}

/// Production TOV profiles begin at r_min = 1e-5 km, never at r = 0; the analytic tabulation
/// mirrors that so the star is an admissible production profile.
inline constexpr double kRStart_km = 1.0e-5;

inline double UniformGridR(const UniformStar &u, std::size_t i, std::size_t N)
{
	return kRStart_km +
		   static_cast<double>(i) * (u.R_km - kRStart_km) / static_cast<double>(N - 1);
}

/// The analytic star tabulated as a reference-solver background.
inline hartle_ref::Background TabulateUniform(const UniformStar &u, std::size_t N)
{
	hartle_ref::Background b;
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = UniformGridR(u, i, N);
		b.r.push_back(r);
		b.m.push_back(u.m(r));
		b.p.push_back(u.p(r));
		b.eps.push_back(u.rho0);
		b.nu.push_back(u.nu(r));
		b.lambda.push_back(u.lambda(r));
	}
	return b;
}

/// The same analytic star as a production `NStar`, through the public TOVPoint constructor.
/// TOVPoint units: r [km], m [Msun], nu_der [1/cm], p [dyne/cm^2], e [g/cm^3].
inline std::unique_ptr<NStar> UniformProductionStar(const UniformStar &u, std::size_t N)
{
	const double km2_to_gcm3 = Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;
	const double km2_to_dyn = Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_INV_KM2;

	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = UniformGridR(u, i, N);
		pts.emplace_back(r, u.m(r) / Zaki::Physics::SUN_M_KM, u.nuprime(r) / 1.0e5, 0.0,
						 u.p(r) * km2_to_dyn, u.rho0 * km2_to_gcm3, 0.1, std::vector<double>{});
	}
	return std::make_unique<NStar>(pts);
}

/// Adapt a production `StarProfile` into the reference solver's background. The reference then
/// builds its right-hand side from `nu`, `lambda` and `(eps + p)` — never from production's
/// `1/(1 - 2m/r)` coefficient helpers — which is what keeps it an independent oracle.
inline hartle_ref::Background BackgroundFromProfile(const NStar &ns)
{
	hartle_ref::Background b;
	const auto &P = ns.Profile();
	const auto *r = P.GetRadius();
	const auto *m = P.GetMass();
	const auto *p = P.GetPressure();
	const auto *e = P.GetEnergyDensity();
	const auto *nu = P.GetMetricNu();
	const auto *lam = P.GetMetricLambda();
	if (!r || !m || !p || !e || !nu || !lam)
		return b;
	const int n = static_cast<int>(r->Size());
	b.r.reserve(static_cast<std::size_t>(n));
	for (int i = 0; i < n; ++i)
	{
		b.r.push_back((*r)[i]);
		b.m.push_back((*m)[i]);
		b.p.push_back((*p)[i]);
		b.eps.push_back((*e)[i]);
		b.nu.push_back((*nu)[i]);
		b.lambda.push_back((*lam)[i]);
	}
	return b;
}

/// Production's seed-free normalized response, as plain vectors on the profile grid.
inline void ProductionShape(const NStar &ns, std::vector<double> &s, std::vector<double> &sp)
{
	const auto &R = ns.RotationResponse();
	const int n = static_cast<int>(R.omega_bar_over_Omega.Size());
	s.assign(static_cast<std::size_t>(n), 0.0);
	sp.assign(static_cast<std::size_t>(n), 0.0);
	for (int i = 0; i < n; ++i)
	{
		s[static_cast<std::size_t>(i)] = R.omega_bar_over_Omega[i];
		sp[static_cast<std::size_t>(i)] = R.domega_bar_over_Omega_dr[i];
	}
}

/// Node-by-node comparison of two independently normalized shapes.
struct ShapeCmp
{
	std::size_t n = 0;
	double max_abs_s = 0, max_rel_s = 0, rms_s = 0, r_worst_s = 0;
	double max_abs_sp = 0, scale_sp = 0, rel_to_scale_sp = 0, max_rel_sp = 0, rms_sp = 0,
		   r_worst_sp = 0;
};

inline ShapeCmp CompareShapes(const std::vector<double> &sA, const std::vector<double> &spA,
							  const std::vector<double> &sB, const std::vector<double> &spB,
							  const std::vector<double> &rr)
{
	ShapeCmp c;
	c.n = std::min({sA.size(), spA.size(), sB.size(), spB.size(), rr.size()});
	auto rel = [](double a, double b) {
		return std::fabs(b) > 0.0 ? std::fabs(a - b) / std::fabs(b) : std::fabs(a - b);
	};

	for (std::size_t i = 0; i < c.n; ++i)
		c.scale_sp = std::max(c.scale_sp, std::fabs(spB[i]));

	double acc_s = 0.0, acc_sp = 0.0;
	for (std::size_t i = 0; i < c.n; ++i)
	{
		const double ds = std::fabs(sA[i] - sB[i]);
		if (ds > c.max_abs_s)
		{
			c.max_abs_s = ds;
			c.r_worst_s = rr[i];
		}
		c.max_rel_s = std::max(c.max_rel_s, rel(sA[i], sB[i]));
		acc_s += ds * ds;

		const double dsp = std::fabs(spA[i] - spB[i]);
		if (dsp > c.max_abs_sp)
		{
			c.max_abs_sp = dsp;
			c.r_worst_sp = rr[i];
		}
		// A pointwise ratio is only meaningful where s' is not itself ~0 (it vanishes at the
		// centre by regularity).
		if (std::fabs(spB[i]) >= 1e-2 * c.scale_sp)
			c.max_rel_sp = std::max(c.max_rel_sp, rel(spA[i], spB[i]));
		acc_sp += dsp * dsp;
	}
	if (c.n > 0)
	{
		c.rms_s = std::sqrt(acc_s / static_cast<double>(c.n));
		c.rms_sp = std::sqrt(acc_sp / static_cast<double>(c.n));
	}
	c.rel_to_scale_sp = c.scale_sp > 0.0 ? c.max_abs_sp / c.scale_sp : 0.0;
	return c;
}

/**
 * @brief The moment of inertia recovered from a normalized shape by volume quadrature:
 *
 *     I = (8 pi / 3) int_0^R r^4 (eps + p) e^{lambda - nu} [omega_bar(r)/Omega] dr .
 *
 * Derived by integrating the conservative form of the Hartle equation
 * (`docs/validation/HARTLE_MOMENT_INERTIA.md` §8). It reads the INTERIOR of the shape, so it
 * is sensitive to exactly what a surface-only check cannot see — and because Phase 4A made
 * `omega_bar/Omega` the actual production response object, this now tests production's own
 * published profile rather than a reference copy of it.
 *
 * Trapezoid on the background grid, matching the interpolation order the profile itself
 * carries (INV-13).
 */
inline double VolumeIntegralI(const hartle_ref::Background &bg, const std::vector<double> &s)
{
	const std::size_t n = std::min(bg.N(), s.size());
	if (n < 2)
		return 0.0;
	auto integrand = [&](std::size_t i) {
		const double ri = bg.r[i];
		return ri * ri * ri * ri * (bg.eps[i] + bg.p[i]) *
			   std::exp(bg.lambda[i] - bg.nu[i]) * s[i];
	};
	double acc = 0.0;
	for (std::size_t i = 0; i + 1 < n; ++i)
		acc += 0.5 * (integrand(i) + integrand(i + 1)) * (bg.r[i + 1] - bg.r[i]);
	return (8.0 * M_PI / 3.0) * acc;
}

} // namespace hartle_4b

#endif /* CompactStar_Tests_HartleProfileCompare_H */
