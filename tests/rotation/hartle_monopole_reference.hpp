// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_reference.hpp
 * @brief Phase 4D — INDEPENDENT test-only reference solvers for the Hartle O(Omega^2)
 *        monopole (l = 0) response. Never compiled into production.
 *
 * WHY THIS EXISTS. Production (`RotationSolver::ODE_HartleMonopole_`, ADR-0007 P2) integrates
 * Hartle's (m0, p0*) system — H67 eqs. (97) and (100). To validate it, something must compute
 * the same physics by a DIFFERENT route. This header does so with Hartle's *other* pair of
 * variables, (m0, h0):
 *
 *     H67 (97)   dm0/dr  = 4 pi r^2 (E+P)(dE/dP) p0*  + (1/12) j^2 r^4 (omegabar')^2
 *                          - (1/3) r^3 (dj^2/dr) omegabar^2
 *     H67 (98)   dh0/dr  = m0 (1 + 8 pi r^2 P)/(r-2M)^2 + 4 pi (E+P) r^2 p0* /(r-2M)
 *                          - (1/12) r^4 j^2 (omegabar')^2/(r-2M)
 *     H67 (90)   p0* = gamma - h0 + (1/3) e^{-nu_H} r^2 omegabar^2        (first integral)
 *
 * with p0* an ALGEBRAIC by-product of the first integral, never an integrated variable. In
 * CompactStar's convention (nu_H = 2 nu, j^2 = e^{-2nu}(1-2m/r), j^2/(r-2m) = e^{-2nu}/r,
 * dj^2/dr = -8 pi r (eps+p) e^{-2nu}) and per unit Omega_geom^2 (s = omegabar/Omega,
 * s' = omegabar'/Omega, mhat = m0/Omega^2, hhat = h0/Omega^2, phat = p0* /Omega^2):
 *
 *     dmhat/dr = 4 pi r^2 (eps+p)(deps/dp) phat + (1/12) r^4 e^{-2nu}(1-2m/r) s'^2
 *                + (8 pi/3) r^4 (eps+p) e^{-2nu} s^2
 *     dhhat/dr = mhat (1 + 8 pi r^2 p)/(r-2m)^2 + 4 pi (eps+p) r^2 phat/(r-2m)
 *                - (1/12) r^3 e^{-2nu} s'^2
 *     phat     = gammahat - hhat + (1/3) r^2 e^{-2nu} s^2
 *
 * Production's (100) is what one gets by differentiating (90) and substituting (98) — Hartle
 * says so on p. 1021 — so the two formulations are equivalent for the exact background but
 * have different state vectors and a different pressure-side equation, which is exactly what an
 * independent check needs. The only assumption the equivalence rests on is that the tabulated
 * nu' is the derivative of the tabulated nu; on a linearly interpolated background that holds
 * to O(h^2), and Phase 4D measures that floor explicitly rather than assuming it away.
 *
 * WHAT THIS HEADER DOES NOT DO. It never calls `RotationSolver::ODE_HartleMonopole_`,
 * `NStar::ComputeHartleMonopoleResponse`, `HartleMonopoleResponse::At`,
 * `Geometry::MetricDenominator` or any production helper. Its interpolation is its own
 * (binary search, not a cached bracket walk). Its centre initialisation is its own: the
 * leading power of each right-hand side is used numerically rather than the closed-form series
 * being transcribed. Its exterior extraction is its own arithmetic on its own solution.
 *
 * The background (profile columns) and the EOS derivative are INPUTS to both solvers, not the
 * thing under test; Phase 4B verified the first-order shape and Phase 4C-I0 the derivative.
 *
 * Also here: a CONTINUUM solver for the exact Schwarzschild constant-density interior with no
 * tabulation at all — first order and second order integrated together on the closed-form
 * background, with p0* obtained BOTH from the first integral (90) and by integrating (100)
 * directly, so that Hartle's first integral can be tested at the continuum level.
 */

#ifndef CompactStar_Tests_HartleMonopoleReference_H
#define CompactStar_Tests_HartleMonopoleReference_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "hartle_profile_compare.hpp" // hartle_4b::UniformStar — the analytic background

namespace hartle_mono_ref
{

// ===========================================================================
//  Tabulated background for the (m0, h0) solver. Geometric units (km).
// ===========================================================================
struct Background2
{
	std::vector<double> r;	  // km, increasing
	std::vector<double> p;	  // km^-2
	std::vector<double> eps;  // km^-2
	std::vector<double> m;	  // km
	std::vector<double> nu;	  // dimensionless  (g_tt = -e^{2 nu})
	std::vector<double> nup;  // km^-1
	std::vector<double> dedp; // dimensionless  (EOS authority, Phase 4C-I0)
	std::vector<double> s;	  // omegabar/Omega
	std::vector<double> sp;	  // omegabar'/Omega, km^-1

	std::size_t N() const { return r.size(); }
	double R() const { return r.back(); }

	bool Consistent() const
	{
		const std::size_t n = r.size();
		return n >= 4 && p.size() == n && eps.size() == n && m.size() == n && nu.size() == n &&
			   nup.size() == n && dedp.size() == n && s.size() == n && sp.size() == n;
	}

	/// Linear interpolation by binary search (INV-13: the background is tabulated; nothing
	/// smoother is justified). Independent of production's cached-bracket walk.
	static double Lerp(const std::vector<double> &x, const std::vector<double> &y, double xq)
	{
		if (xq <= x.front())
			return y.front();
		if (xq >= x.back())
			return y.back();
		std::size_t lo = 0, hi = x.size() - 1;
		while (hi - lo > 1)
		{
			const std::size_t mid = (lo + hi) / 2;
			if (x[mid] <= xq)
				lo = mid;
			else
				hi = mid;
		}
		const double t = (xq - x[lo]) / (x[hi] - x[lo]);
		return y[lo] + t * (y[hi] - y[lo]);
	}

	struct Point
	{
		double p, eps, m, nu, nup, dedp, s, sp;
	};
	Point At(double rq) const
	{
		return Point{Lerp(r, p, rq),   Lerp(r, eps, rq),  Lerp(r, m, rq),
					 Lerp(r, nu, rq),  Lerp(r, nup, rq),  Lerp(r, dedp, rq),
					 Lerp(r, s, rq),   Lerp(r, sp, rq)};
	}
};

// ===========================================================================
//  The (m0, h0) solver on a tabulated background.
// ===========================================================================
struct MHOptions
{
	bool sources_on = true;	   ///< false = the regular HOMOGENEOUS solution (no omegabar terms)
	double gammahat = 0.0;	   ///< the first-integral constant; see Solve for how it is set
	double phat_at_r0 = std::numeric_limits<double>::quiet_NaN(); ///< NaN = regular series
	double r0 = std::numeric_limits<double>::quiet_NaN();		  ///< NaN = first grid radius
	double rtol = 1e-13;	   ///< tighter than production's 1e-10
	double atol = 1e-16;
	double h_init = 1e-6;
	double I_exterior = 0.0;   ///< the I to use in I^2/R^3 (independent or production, caller's choice)
	/// ADR-0008 (ACCEPTED 2026-09-03): evaluate the EOS energy-density contribution to dm0/dr as
	/// the MEASURE -4 pi r^2 xi0_hat d(eps), one profile segment at a time, instead of the
	/// pointwise 4 pi r^2 (eps+p)(deps/dp) p0* rewrite. false keeps the ADR-0007 differential
	/// form, which is what a smooth-EOS equivalence check compares against.
	bool eos_measure = false;
	const gsl_odeiv2_step_type *step = gsl_odeiv2_step_rk8pd; ///< admissibility sweeps may vary this
};

struct MHResult
{
	bool ok = false;
	std::string message;
	std::vector<double> mhat, hhat, phat, dphat, xihat; // on the background grid
	double gammahat = 0.0;
	double R_star = 0.0, eps_star = 0.0;
	double mhat_R = 0.0, phat_R = 0.0, xihat_R = 0.0;
	double shell_hat = 0.0, exterior_hat = 0.0, deltaM_hat = 0.0;
	double r0_used = 0.0;
};

namespace detail
{
struct Params
{
	const Background2 *bg;
	bool sources;
	double gammahat;
	bool measure = false;	 ///< ADR-0008 source form
	double eps_slope = 0.0;	 ///< d(eps)/dr of the segment the driver is inside [km^-3]
};

/// y[0] = mhat, y[1] = hhat.
inline int RHS(double r, const double y[], double f[], void *pv)
{
	const Params *P = static_cast<const Params *>(pv);
	const auto b = P->bg->At(r);

	const double D = 1.0 - 2.0 * b.m / r; // metric denominator, plain expression
	const double r_2m = r * D;
	const double e2 = std::exp(-2.0 * b.nu);
	const double ep = b.eps + b.p;
	const double r2 = r * r, r3 = r2 * r, r4 = r3 * r;

	const double s = P->sources ? b.s : 0.0;
	const double sp = P->sources ? b.sp : 0.0;

	// H67 (90): p0* from the first integral
	const double phat = P->gammahat - y[1] + (1.0 / 3.0) * r2 * e2 * s * s;

	// H67 (97). The EOS term is either Hartle's smooth rewrite (ADR-0007 P2) or the measure
	// -4 pi r^2 xi0_hat d(eps) of his eq. (93) on the current segment (ADR-0008 Q1/Q3).
	double term1;
	if (P->measure)
	{
		term1 = 0.0;
		if (P->eps_slope != 0.0)
		{
			const double xi = phat / b.nup;
			if (!std::isfinite(xi))
				return GSL_EBADFUNC;
			term1 = -4.0 * M_PI * r2 * xi * P->eps_slope;
		}
	}
	else
		term1 = 4.0 * M_PI * r2 * ep * b.dedp * phat;

	f[0] = term1 + (1.0 / 12.0) * r4 * e2 * D * sp * sp +
		   (8.0 * M_PI / 3.0) * r4 * ep * e2 * s * s;

	// H67 (98)
	f[1] = y[0] * (1.0 + 8.0 * M_PI * r2 * b.p) / (r_2m * r_2m) +
		   4.0 * M_PI * ep * r2 * phat / r_2m - (1.0 / 12.0) * r3 * e2 * sp * sp;

	if (!std::isfinite(f[0]) || !std::isfinite(f[1]) || !(D > 0.0))
		return GSL_EBADFUNC;
	return GSL_SUCCESS;
}
} // namespace detail

/**
 * @brief Integrate the (m0, h0) system from the centre to the surface.
 *
 * Centre: fixed-eps_c is p0*(0) = 0 and m0(0) = 0 (H67 108, p. 1009; HT68 §II f). With
 * hhat(r0) = 0 the first-integral constant is chosen so that phat(r0) takes the regular-series
 * value; mhat(r0) is set from the leading power of its own right-hand side,
 * mhat(r0) = r0 f_m(r0)/(k+1) with k = 4 for the rotational source and k = 2 for the
 * homogeneous family — a numerical use of the series rather than a transcription of the
 * closed-form coefficient, so the initialisation is genuinely this solver's own.
 *
 * Surface: R_* is the last background node (INV-06). delta M is assembled EXACTLY as ADR-0007
 * P6 defines it, from this solver's own mhat and xihat and the caller's I.
 */
inline MHResult Solve(const Background2 &bg, const MHOptions &opt)
{
	MHResult out;
	if (!bg.Consistent())
	{
		out.message = "inconsistent background";
		return out;
	}
	const std::size_t n = bg.N();

	double r0 = std::isfinite(opt.r0) ? opt.r0 : bg.r.front();
	if (!(r0 > 0.0))
		r0 = 1e-6;
	out.r0_used = r0;

	const auto b0 = bg.At(r0);
	const double e2_0 = std::exp(-2.0 * b0.nu);
	const double j2_0 = e2_0 * (1.0 - 2.0 * b0.m / r0);

	// phat(r0): regular series (1/3) j^2 s^2 r0^2 for the rotational family, or the caller's
	// value (e.g. 1 for the homogeneous family).
	double phat0;
	if (std::isfinite(opt.phat_at_r0))
		phat0 = opt.phat_at_r0;
	else
		phat0 = opt.sources_on ? (1.0 / 3.0) * j2_0 * b0.s * b0.s * r0 * r0 : 0.0;

	// hhat(r0) = 0 (its expansion begins at O(r^4)); gammahat follows from (90).
	const double s0 = opt.sources_on ? b0.s : 0.0;
	const double gammahat = phat0 - (1.0 / 3.0) * r0 * r0 * e2_0 * s0 * s0;

	detail::Params P{&bg, opt.sources_on, gammahat, opt.eos_measure, 0.0};
	if (opt.eos_measure && n > 1)
		P.eps_slope = (bg.eps[1] - bg.eps[0]) / (bg.r[1] - bg.r[0]);

	// mhat(r0) from the leading power of its own RHS.
	double y[2] = {0.0, 0.0};
	{
		double f[2];
		if (detail::RHS(r0, y, f, &P) != GSL_SUCCESS)
		{
			out.message = "RHS failed at r0";
			return out;
		}
		const int k = opt.sources_on ? 4 : 2;
		y[0] = r0 * f[0] / static_cast<double>(k + 1);
	}

	gsl_odeiv2_system sys = {detail::RHS, nullptr, 2, &P};
	gsl_odeiv2_driver *d =
		gsl_odeiv2_driver_alloc_y_new(&sys, opt.step, opt.h_init, opt.atol, opt.rtol);

	out.mhat.assign(n, 0.0);
	out.hhat.assign(n, 0.0);
	out.phat.assign(n, 0.0);
	out.dphat.assign(n, 0.0);
	out.xihat.assign(n, 0.0);

	double r = r0;
	for (std::size_t i = 0; i < n; ++i)
	{
		const double rt = bg.r[i];
		// ADR-0008 Q3: one governed segment per driver call, with that segment's own measure
		// density installed first. Node boundaries are hard integration boundaries.
		if (opt.eos_measure && i > 0)
			P.eps_slope = (bg.eps[i] - bg.eps[i - 1]) / (bg.r[i] - bg.r[i - 1]);
		if (rt > r0)
		{
			if (gsl_odeiv2_driver_apply(d, &r, rt, y) != GSL_SUCCESS)
			{
				out.message = "GSL failed at r = " + std::to_string(rt);
				gsl_odeiv2_driver_free(d);
				return out;
			}
		}
		const double e2 = std::exp(-2.0 * bg.nu[i]);
		const double s = opt.sources_on ? bg.s[i] : 0.0;
		const double phat = gammahat - y[1] + (1.0 / 3.0) * rt * rt * e2 * s * s;
		out.mhat[i] = y[0];
		out.hhat[i] = y[1];
		out.phat[i] = phat;
		out.dphat[i] = (bg.eps[i] + bg.p[i]) * phat;
		double xi = phat / bg.nup[i];
		if (!std::isfinite(xi))
		{
			// regular-centre limit of p0* /nu' (PHASE4C_HARTLE2_DERIVATION.md §9.2)
			const double den = 4.0 * M_PI * (bg.eps[i] + 3.0 * bg.p[i]);
			xi = (den != 0.0 && opt.sources_on) ? j2_0 * b0.s * b0.s * rt / den : 0.0;
		}
		out.xihat[i] = xi;
	}
	gsl_odeiv2_driver_free(d);

	const std::size_t last = n - 1;
	out.gammahat = gammahat;
	out.R_star = bg.r[last];
	out.eps_star = bg.eps[last];
	out.mhat_R = out.mhat[last];
	out.phat_R = out.phat[last];
	out.xihat_R = out.xihat[last];
	out.shell_hat = 4.0 * M_PI * out.R_star * out.R_star * out.eps_star * out.xihat_R;
	out.exterior_hat = opt.I_exterior * opt.I_exterior / (out.R_star * out.R_star * out.R_star);
	out.deltaM_hat = out.mhat_R + out.shell_hat + out.exterior_hat;

	for (std::size_t i = 0; i < n; ++i)
		if (!std::isfinite(out.mhat[i]) || !std::isfinite(out.phat[i]) ||
			!std::isfinite(out.xihat[i]))
		{
			out.message = "non-finite output";
			return out;
		}
	out.ok = true;
	return out;
}

// ===========================================================================
//  The CONTINUUM solver on the exact constant-density interior — no tabulation.
//
//  State: y[0] = omegabar (raw), y[1] = q = r^4 j omegabar', y[2] = m0 (raw),
//         y[3] = h0 (raw), y[4] = p0* (raw) integrated DIRECTLY by H67 (100).
//  Everything raw is divided by Omega_raw^2 (or Omega_raw) at the end, so no seed survives.
//  The first integral (90) then gives a second p0* from h0; the difference between the two is
//  a pure continuum test of the transcription of (98), (100) and (90) against each other.
// ===========================================================================
struct ContinuumNode
{
	double r = 0, s = 0, sp = 0, mhat = 0, hhat = 0, phat_from_90 = 0, phat_from_100 = 0,
		   xihat = 0, dphat = 0;
};

struct ContinuumResult
{
	bool ok = false;
	std::string message;
	std::vector<ContinuumNode> nodes;
	double Omega_raw = 0, J_raw = 0, I = 0;
	double R = 0, eps_R = 0;
	double mhat_R = 0, phat_R = 0, xihat_R = 0, shell_hat = 0, exterior_hat = 0, deltaM_hat = 0;
	double first_integral_max_residual = 0; ///< max |phat_100 - phat_90| / scale over nodes
	double first_integral_scale = 0;
};

namespace detail
{
struct UParams
{
	const hartle_4b::UniformStar *u;
	double gamma_raw;
};

inline int URHS(double r, const double y[], double f[], void *pv)
{
	const UParams *P = static_cast<const UParams *>(pv);
	const auto &u = *P->u;
	const double eps = u.rho0, p = u.p(r), m = u.m(r), nu = u.nu(r), nup = u.nuprime(r);
	const double lam = u.lambda(r);
	const double D = 1.0 - 2.0 * m / r;
	const double r_2m = r * D;
	const double e2 = std::exp(-2.0 * nu);
	const double j = std::exp(-(nu + lam));
	const double ep = eps + p;
	const double r2 = r * r, r3 = r2 * r, r4 = r3 * r;

	const double ob = y[0];
	const double dob = y[1] / (r4 * j); // omegabar' from the conserved flux
	// first order (conservative form, as hartle_reference.hpp)
	f[0] = dob;
	f[1] = 16.0 * M_PI * r4 * ep * std::exp(lam - nu) * ob;

	// second order: constant density => deps/dp = 0, so the p0* term of (97) vanishes
	const double m0 = y[2], h0 = y[3], p0_100 = y[4];
	const double p0_90 = P->gamma_raw - h0 + (1.0 / 3.0) * r2 * e2 * ob * ob;

	f[2] = (1.0 / 12.0) * r4 * e2 * D * dob * dob + (8.0 * M_PI / 3.0) * r4 * ep * e2 * ob * ob;
	f[3] = m0 * (1.0 + 8.0 * M_PI * r2 * p) / (r_2m * r_2m) + 4.0 * M_PI * ep * r2 * p0_90 / r_2m -
		   (1.0 / 12.0) * r3 * e2 * dob * dob;
	// H67 (100) for p0* directly — the production-side equation, integrated here at tight
	// tolerance so the first integral can be tested at the continuum level
	f[4] = -m0 * (1.0 + 8.0 * M_PI * r2 * p) / (r_2m * r_2m) - 4.0 * M_PI * ep * r2 * p0_100 / r_2m +
		   (1.0 / 12.0) * r3 * e2 * dob * dob + (2.0 / 3.0) * r * e2 * ob * (ob + r * dob - r * nup * ob);

	for (int k = 0; k < 5; ++k)
		if (!std::isfinite(f[k]))
			return GSL_EBADFUNC;
	return GSL_SUCCESS;
}
} // namespace detail

/// Integrate on the closed-form background from r0 to the surface, recording at @p r_out.
inline ContinuumResult SolveUniformContinuum(const hartle_4b::UniformStar &u,
											 const std::vector<double> &r_out, double r0,
											 double rtol = 1e-13, double atol = 1e-18)
{
	ContinuumResult out;
	if (r_out.empty() || !(r0 > 0.0))
	{
		out.message = "bad inputs";
		return out;
	}
	const double R = u.R_km;

	// Centre: omegabar(r0) = 1 (raw seed), q(r0) = 0 (regular); m0, h0 from the series with
	// gamma_raw chosen so that p0*(r0) is the regular-series value.
	const double nu0 = u.nu(r0), lam0 = u.lambda(r0);
	const double e2_0 = std::exp(-2.0 * nu0);
	const double j2_0 = std::exp(-2.0 * (nu0 + lam0));
	const double ob0 = 1.0;
	const double p0_series = (1.0 / 3.0) * j2_0 * ob0 * ob0 * r0 * r0;
	const double gamma_raw = p0_series - (1.0 / 3.0) * r0 * r0 * e2_0 * ob0 * ob0; // h0(r0)=0

	detail::UParams P{&u, gamma_raw};
	double y[5] = {ob0, 0.0, 0.0, 0.0, p0_series};
	{
		double f[5];
		if (detail::URHS(r0, y, f, &P) != GSL_SUCCESS)
		{
			out.message = "RHS failed at r0";
			return out;
		}
		y[2] = r0 * f[2] / 5.0; // m0 ~ r^5
	}

	gsl_odeiv2_system sys = {detail::URHS, nullptr, 5, &P};
	gsl_odeiv2_driver *d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-7, atol, rtol);

	std::vector<double> targets = r_out;
	std::sort(targets.begin(), targets.end());
	if (targets.back() < R)
		targets.push_back(R);

	std::vector<double> raw_r, raw_ob, raw_dob, raw_m0, raw_h0, raw_p100;
	double r = r0;
	for (double rt : targets)
	{
		if (rt > r0)
		{
			if (gsl_odeiv2_driver_apply(d, &r, rt, y) != GSL_SUCCESS)
			{
				out.message = "GSL failed at r = " + std::to_string(rt);
				gsl_odeiv2_driver_free(d);
				return out;
			}
		}
		const double rr = std::max(rt, r0);
		const double jj = std::exp(-(u.nu(rr) + u.lambda(rr)));
		raw_r.push_back(rt);
		raw_ob.push_back(y[0]);
		raw_dob.push_back(rt > r0 ? y[1] / (rr * rr * rr * rr * jj) : 0.0);
		raw_m0.push_back(y[2]);
		raw_h0.push_back(y[3]);
		raw_p100.push_back(y[4]);
	}
	gsl_odeiv2_driver_free(d);

	// Own exterior extraction at the surface
	const double obR = raw_ob.back(), dobR = raw_dob.back();
	out.Omega_raw = obR + R * dobR / 3.0;
	out.J_raw = R * R * R * R * dobR / 6.0;
	out.I = out.J_raw / out.Omega_raw;
	const double W2 = out.Omega_raw * out.Omega_raw;

	double scale = 0.0, maxres = 0.0;
	for (std::size_t i = 0; i < raw_r.size(); ++i)
	{
		ContinuumNode nd;
		nd.r = raw_r[i];
		const double rr = std::max(nd.r, r0);
		const double e2 = std::exp(-2.0 * u.nu(rr));
		nd.s = raw_ob[i] / out.Omega_raw;
		nd.sp = raw_dob[i] / out.Omega_raw;
		nd.mhat = raw_m0[i] / W2;
		nd.hhat = raw_h0[i] / W2;
		const double p90 = gamma_raw - raw_h0[i] + (1.0 / 3.0) * rr * rr * e2 * raw_ob[i] * raw_ob[i];
		nd.phat_from_90 = p90 / W2;
		nd.phat_from_100 = raw_p100[i] / W2;
		nd.dphat = (u.rho0 + u.p(rr)) * nd.phat_from_90;
		const double nup = u.nuprime(rr);
		double xi = nd.phat_from_90 / nup;
		if (!std::isfinite(xi))
			xi = j2_0 * nd.s * nd.s * rr / (4.0 * M_PI * (u.rho0 + 3.0 * u.p(rr)));
		nd.xihat = xi;
		scale = std::max(scale, std::fabs(nd.phat_from_90));
		out.nodes.push_back(nd);
	}
	for (const auto &nd : out.nodes)
		maxres = std::max(maxres, std::fabs(nd.phat_from_100 - nd.phat_from_90));
	out.first_integral_scale = scale;
	out.first_integral_max_residual = (scale > 0.0) ? maxres / scale : maxres;

	const auto &S = out.nodes.back();
	out.R = R;
	out.eps_R = u.rho0;
	out.mhat_R = S.mhat;
	out.phat_R = S.phat_from_90;
	out.xihat_R = S.xihat;
	out.shell_hat = 4.0 * M_PI * R * R * u.rho0 * out.xihat_R;
	out.exterior_hat = out.I * out.I / (R * R * R);
	out.deltaM_hat = out.mhat_R + out.shell_hat + out.exterior_hat;
	out.ok = true;
	return out;
}

// ===========================================================================
//  Comparison metric shared by the 4D harnesses (the same rule Phase 4B used).
// ===========================================================================
struct Cmp
{
	std::size_t n = 0;
	double max_abs = 0, max_rel = 0, rms = 0, r_worst = 0, scale = 0, rel_to_scale = 0;
	double surface_rel = 0;
};

/// max_rel is taken only where |ref| >= 1 % of the profile's peak; rel_to_scale is the largest
/// absolute difference over that peak — the two honest measures for a field that vanishes at
/// the centre.
inline Cmp Compare(const std::vector<double> &a, const std::vector<double> &ref,
				   const std::vector<double> &rr)
{
	Cmp c;
	c.n = std::min({a.size(), ref.size(), rr.size()});
	for (std::size_t i = 0; i < c.n; ++i)
		c.scale = std::max(c.scale, std::fabs(ref[i]));
	double acc = 0.0;
	for (std::size_t i = 0; i < c.n; ++i)
	{
		const double d = std::fabs(a[i] - ref[i]);
		if (d > c.max_abs)
		{
			c.max_abs = d;
			c.r_worst = rr[i];
		}
		if (std::fabs(ref[i]) >= 1e-2 * c.scale && c.scale > 0.0)
			c.max_rel = std::max(c.max_rel, d / std::fabs(ref[i]));
		acc += d * d;
	}
	if (c.n)
	{
		c.rms = std::sqrt(acc / static_cast<double>(c.n));
		const double rl = ref[c.n - 1];
		c.surface_rel = std::fabs(rl) > 0 ? std::fabs(a[c.n - 1] - rl) / std::fabs(rl)
										 : std::fabs(a[c.n - 1] - rl);
	}
	c.rel_to_scale = c.scale > 0.0 ? c.max_abs / c.scale : 0.0;
	return c;
}

} // namespace hartle_mono_ref

#endif /* CompactStar_Tests_HartleMonopoleReference_H */
