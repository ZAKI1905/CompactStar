// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_reference.hpp
 * @brief INDEPENDENT test-side first-order Hartle frame-dragging solver.
 *
 * Phase 2B-4B. This is a reference implementation written to validate production; it is
 * NOT a copy of it and deliberately does not call
 *   RotationSolver::ODE_N_Fast,
 *   RotationSolver::GetHartleOmegaCoeff_N_Fast,
 *   RotationSolver::GetHartleDOmegaCoeff_N_Fast
 * or reuse their expanded second-order form.
 *
 * INDEPENDENT REPRESENTATION — the CONSERVATIVE Hartle system.
 * Hartle (1967), "Slowly Rotating Relativistic Stars. I. Equations of Structure",
 * ApJ 150, 1005, writes the frame-dragging equation in the self-adjoint form
 *
 *      (1/r^4) d/dr [ r^4 j d(omegabar)/dr ] + (4/r) (dj/dr) omegabar = 0 ,
 *      j(r) = exp[-(nu + lambda)] .
 *
 * Introducing the conserved-flux variable
 *
 *      q(r) = r^4 j d(omegabar)/dr                     [km^2]
 *
 * turns it into a FIRST-ORDER system in (omegabar, q):
 *
 *      d(omegabar)/dr = q / (r^4 j)
 *      dq/dr          = -4 r^3 j' omegabar
 *                     = +16 pi r^4 (eps + p) exp(lambda - nu) omegabar ,
 *
 * where the second line uses the exact identity derived in
 * docs/validation/HARTLE_MOMENT_INERTIA.md §6:
 *
 *      j'/j = -(nu' + lambda') = -4 pi r (eps + p) / (1 - 2m/r) ,
 *      so    j' = -4 pi r (eps + p) exp(lambda - nu) .
 *
 * The RHS is therefore built from the METRIC columns (nu, lambda) and (eps + p) directly,
 * never from production's 1/(1 - 2m/r) coefficient helpers. Different state variables,
 * different right-hand side, same physics.
 *
 * UNITS. Geometric, exactly as StarProfile stores them:
 *   r, m [km];  p, eps [km^-2];  nu, lambda dimensionless;
 *   omegabar [km^-1];  q [km^2];  J [km^2];  Omega [km^-1];  I [km^3].
 *
 * SCALE FREEDOM. The system is linear and homogeneous, so omegabar -> A*omegabar maps
 * q -> A*q, J -> A*J, Omega -> A*Omega and leaves I = J/Omega invariant. The central
 * amplitude is a free parameter of this solver precisely so that invariance can be
 * measured rather than assumed.
 */

#ifndef CompactStar_Tests_HartleReference_H
#define CompactStar_Tests_HartleReference_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

namespace hartle_ref
{

/// Tabulated spherical background in geometric units, as a StarProfile provides it.
struct Background
{
	std::vector<double> r;	   // km, increasing, r[0] >= 0
	std::vector<double> m;	   // km
	std::vector<double> p;	   // km^-2
	std::vector<double> eps;   // km^-2
	std::vector<double> nu;	   // dimensionless
	std::vector<double> lambda; // dimensionless

	std::size_t N() const { return r.size(); }
	double R() const { return r.back(); }

	/// Linear interpolation with clamping — the background is tabulated, so nothing
	/// smoother is justified (INV-13).
	static double Interp(const std::vector<double> &x, const std::vector<double> &y, double xq)
	{
		const std::size_t n = x.size();
		if (xq <= x.front()) return y.front();
		if (xq >= x.back()) return y.back();
		std::size_t lo = 0, hi = n - 1;
		while (hi - lo > 1)
		{
			const std::size_t mid = (lo + hi) / 2;
			if (x[mid] <= xq) lo = mid; else hi = mid;
		}
		const double t = (xq - x[lo]) / (x[hi] - x[lo]);
		return y[lo] + t * (y[hi] - y[lo]);
	}

	/// j(r) = exp[-(nu + lambda)]
	double j_at(double rq) const
	{
		return std::exp(-(Interp(r, nu, rq) + Interp(r, lambda, rq)));
	}
	/// The source weight  w(r) = (eps + p) * exp(lambda - nu)   [km^-2]
	double w_at(double rq) const
	{
		const double ep = Interp(r, eps, rq) + Interp(r, p, rq);
		return ep * std::exp(Interp(r, lambda, rq) - Interp(r, nu, rq));
	}
};

struct Result
{
	bool ok = false;
	std::string message;
	double seed = 0;			 // omegabar(r0)
	double r0 = 0;				 // start radius [km]
	double omega_bar_R = 0;		 // omegabar(R)      [km^-1]
	double domega_bar_R = 0;	 // omegabar'(R)     [km^-2]
	double J = 0;				 // R^4 omegabar'(R)/6   [km^2]
	double Omega = 0;			 // omegabar(R) + R omegabar'(R)/3  [km^-1]
	double I_surface = 0;		 // J/Omega          [km^3]
	double I_volume = 0;		 // (8pi/3) int r^4 (eps+p) e^{lam-nu} (omegabar/Omega) dr
	std::vector<double> ombar;	 // solution on the background grid
};

/// GSL right-hand side for the conservative system. y[0] = omegabar, y[1] = q.
inline int RHS(double r, const double y[], double f[], void *params)
{
	const Background *bg = static_cast<const Background *>(params);
	const double rs = (r > 1e-12) ? r : 1e-12;
	const double r4 = rs * rs * rs * rs;
	f[0] = y[1] / (r4 * bg->j_at(rs));
	f[1] = 16.0 * M_PI * r4 * bg->w_at(rs) * y[0];
	return GSL_SUCCESS;
}

/**
 * @brief Solve the conservative Hartle system on @p bg from @p r0 to the surface.
 *
 * Regular centre: omegabar'(0) = 0 implies q = r^4 j omegabar' = O(r^5), so q(r0) = 0 is
 * the regular start and its truncation error is O(r0^5) — far smaller than production's
 * equivalent condition at the same radius. @p r0 is a free parameter so that centre
 * sensitivity can be measured (§8 of the task).
 */
inline Result Solve(const Background &bg, double seed, double r0,
					double rtol = 1e-12, double atol = 1e-14)
{
	Result out;
	out.seed = seed;
	out.r0 = r0;
	if (bg.N() < 4 || !(bg.R() > 0.0))
	{
		out.message = "degenerate background";
		return out;
	}

	gsl_odeiv2_system sys = {RHS, nullptr, 2,
							 const_cast<void *>(static_cast<const void *>(&bg))};
	gsl_odeiv2_driver *d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, atol, rtol);

	double y[2] = {seed, 0.0};
	double r = r0;
	out.ombar.assign(bg.N(), 0.0);

	for (std::size_t i = 0; i < bg.N(); ++i)
	{
		const double rt = bg.r[i];
		if (rt <= r0)
		{
			out.ombar[i] = y[0];
			continue;
		}
		if (gsl_odeiv2_driver_apply(d, &r, rt, y) != GSL_SUCCESS)
		{
			out.message = "GSL failed at r = " + std::to_string(rt);
			gsl_odeiv2_driver_free(d);
			return out;
		}
		out.ombar[i] = y[0];
	}
	gsl_odeiv2_driver_free(d);

	const double R = bg.R();
	const double jR = bg.j_at(R);
	out.omega_bar_R = y[0];
	out.domega_bar_R = y[1] / (R * R * R * R * jR);

	// Vacuum exterior: (r^4 omegabar')' = 0 gives omegabar' = 6J/r^4 and
	// omegabar = Omega - 2J/r^3, hence the two surface relations below.
	out.J = R * R * R * R * out.domega_bar_R / 6.0;
	out.Omega = out.omega_bar_R + R * out.domega_bar_R / 3.0;
	out.I_surface = out.J / out.Omega;

	// Independent volume quadrature of the SAME solution:
	//   I = (8pi/3) int_0^R r^4 (eps+p) e^{lambda-nu} (omegabar/Omega) dr
	// Trapezoid on the background grid. This reads the interior solution rather than the
	// ODE state at R, so a corrupted J extraction shows up as surface/volume disagreement.
	double acc = 0.0;
	auto integrand = [&](std::size_t i) {
		const double ri = bg.r[i];
		const double ep = bg.eps[i] + bg.p[i];
		return ri * ri * ri * ri * ep *
			   std::exp(bg.lambda[i] - bg.nu[i]) * out.ombar[i];
	};
	for (std::size_t i = 0; i + 1 < bg.N(); ++i)
		acc += 0.5 * (integrand(i) + integrand(i + 1)) * (bg.r[i + 1] - bg.r[i]);
	out.I_volume = (8.0 * M_PI / 3.0) * acc / out.Omega;

	out.ok = true;
	return out;
}

} // namespace hartle_ref

#endif
