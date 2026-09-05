// TEST SIDE ONLY. Enthalpy TOV in (r^2,m/r^3), RKF45, direct phase-space EOS.
#pragma once
#include "local_oracle.hpp"
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <stdexcept>
namespace structure1
{
struct OracleStar
{
	double M, R, Hc, Hs;
};
struct EnthalpyOracle
{
	LocalOracle eos;
	double h0 = eos.masses[1] + eos.masses[2];
	double c = GSL_CONST_CGSM_SPEED_OF_LIGHT, G = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT;
	// Numerical unit boundary independent of the production table adapter.
	double to_geo = 1.602176634e33 * G / std::pow(c, 4) * 1e10;
	double solar_length = G * GSL_CONST_CGSM_SOLAR_MASS / (c * c) / 1e5;
	double Hc = 0;
	OraclePoint at_H(double H) const
	{
		if (H <= 0)
			return {};
		double dh = h0 * std::expm1(H), h = h0 + dh;
		const double mp = eos.masses[1], me = eos.masses[2];
		if (h < eos.masses[0])
		{
			double k2 = dh * (h + h0) * (h - (mp - me)) * (h + (mp - me)) / (4 * h * h);
			double n = std::pow(k2, 1.5) / (3 * M_PI * M_PI * std::pow(eos.b, 3));
			return eos.from_ns({0, n, n, 0});
		}
		double t = ((h - mp) * (h + mp) + me * me) / (2 * h);
		if (t > eos.masses[3])
		{
			double a = eos.masses[3], b = h - mp;
			for (int j = 0; j < 64; ++j)
			{
				t = (a + b) / 2;
				if (t == a || t == b)
					break;
				double p = eos.density(t, me) + eos.density(t, eos.masses[3]);
				if (eos.potential(p, mp) + t < h)
					a = t;
				else
					b = t;
			}
			t = (a + b) / 2;
		}
		double ne = eos.density(t, me), nm = eos.density(t, eos.masses[3]);
		return eos.from_ns({eos.density(h, eos.masses[0]), ne + nm, ne, nm});
	}
	double H_for_rho(double rho) const
	{
		double a = .05, b = .3;
		for (int j = 0; j < 75; ++j)
		{
			double h = (a + b) / 2;
			if (h == a || h == b)
				break;
			double r = at_H(h).e * 1.602176634e33 / (c * c);
			if (r < rho)
				a = h;
			else
				b = h;
		}
		return (a + b) / 2;
	}
	double H_for_pressure(double p_cgs) const
	{
		if (p_cgs == 0)
			return 0;
		double a = 0, b = .3;
		for (int j = 0; j < 100; ++j)
		{
			double h = (a + b) / 2;
			if (h == a || h == b)
				break;
			if (at_H(h).p * 1.602176634e33 < p_cgs)
				a = h;
			else
				b = h;
		}
		return (a + b) / 2;
	}
	static int rhs(double x, const double y[], double f[], void *ptr)
	{
		auto &o = *static_cast<EnthalpyOracle *>(ptr);
		auto q = o.at_H(o.Hc - x);
		double e = q.e * o.to_geo, p = q.p * o.to_geo;
		double A = (1 - 2 * y[0] * y[1]) / (y[1] + 4 * M_PI * p);
		f[0] = 2 * A;
		f[1] = (4 * M_PI * e - 3 * y[1]) * A / y[0];
		return std::isfinite(f[0]) && std::isfinite(f[1]) ? GSL_SUCCESS : GSL_EBADFUNC;
	}
	OracleStar solve(double requested_rho, double cut_cgs = 0, double tol = 1e-11, double start = 1e-9)
	{
		Hc = H_for_rho(requested_rho);
		double Hs = H_for_pressure(cut_cgs);
		auto q = at_H(Hc);
		double e = q.e * to_geo, p = q.p * to_geo;
		double delta = start * Hc, x = delta;
		double y[2] = {3 * delta / (2 * M_PI * (e + 3 * p)), 4 * M_PI * e / 3 - 4 * M_PI * (e + p) * q.slope * delta / 5};
		gsl_odeiv2_system sys = {rhs, nullptr, 2, this};
		auto *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, delta, tol * 1e-3, tol);
		int status = gsl_odeiv2_driver_apply(d, &x, Hc - Hs, y);
		gsl_odeiv2_driver_free(d);
		if (status != GSL_SUCCESS)
			throw std::runtime_error("independent enthalpy TOV failure");
		double r = std::sqrt(y[0]);
		return {y[1] * r * r * r / solar_length, r, Hc, Hs};
	}
};
} // namespace structure1
