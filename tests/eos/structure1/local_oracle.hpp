// Test-only phase-space quadratures and common-lepton equilibrium inversion.
// No production thermodynamic evaluation, chemical Hessian, or EOS spline.
#pragma once
#include <Zaki/Physics/Constants.hpp>
#include <array>
#include <cmath>
#include <gsl/gsl_integration.h>
namespace structure1
{
struct OraclePoint
{
	double n, e, p, h, slope;
	std::array<double, 4> ns;
};
struct LocalOracle
{
	std::array<double, 4> masses = {Zaki::Physics::NEUTRON_M_MEV, Zaki::Physics::PROTON_M_MEV, Zaki::Physics::ELECTRON_M_MEV, Zaki::Physics::MUON_M_MEV};
	double b = 1 / Zaki::Physics::MEV_2_INV_FM;
	gsl_integration_glfixed_table *rule = gsl_integration_glfixed_table_alloc(96);
	~LocalOracle() { gsl_integration_glfixed_table_free(rule); }
	double density(double u, double m) const { return u <= m ? 0 : std::pow((u - m) * (u + m), 1.5) / (3 * M_PI * M_PI * b * b * b); }
	double potential(double n, double m) const { return std::hypot(m, b * std::cbrt(3 * M_PI * M_PI * n)); }
	std::array<double, 4> common(double t) const
	{
		double e = density(t, masses[2]), mu = density(t, masses[3]), p = e + mu;
		return {density(potential(p, masses[1]) + t, masses[0]), p, e, mu};
	}
	OraclePoint from_ns(std::array<double, 4> ns) const
	{
		OraclePoint v{};
		v.ns = ns;
		v.n = ns[0] + ns[1];
		std::array<double, 4> S{}, u{};
		for (int i = 0; i < 4; ++i)
		{
			double k = b * std::cbrt(3 * M_PI * M_PI * ns[i]);
			u[i] = std::hypot(masses[i], k);
			S[i] = u[i] * k / (M_PI * M_PI * b * b * b);
			double e = 0, p = 0;
			for (size_t j = 0; j < 96; ++j)
			{
				double q, w;
				gsl_integration_glfixed_point(0, 1, j, &q, &w, rule);
				double z = std::hypot(masses[i], q * k);
				e += w * q * q * z;
				p += w * q * q * q * q / z;
			}
			v.e += 3 * ns[i] * e;
			v.p += ns[i] * k * k * p;
		}
		double C;
		if (ns[0] == 0)
		{
			v.h = u[1] + u[2];
			C = S[1] * S[2] / (S[1] + S[2]);
		}
		else
		{
			v.h = u[0];
			double L = S[2] + S[3];
			C = S[0] + S[1] * L / (S[1] + L);
		}
		v.slope = v.h * C / v.n;
		return v;
	}
	OraclePoint at_n(double n) const
	{
		const double t0 = ((masses[0] - masses[1]) * (masses[0] + masses[1]) + masses[2] * masses[2]) / (2 * masses[0]);
		if (n <= density(t0, masses[2]))
			return from_ns({0, n, n, 0});
		double a = t0, z = 150;
		for (int j = 0; j < 90; ++j)
		{
			double t = (a + z) / 2;
			if (t == a || t == z)
				break;
			auto ns = common(t);
			if (ns[0] + ns[1] < n)
				a = t;
			else
				z = t;
		}
		return from_ns(common((a + z) / 2));
	}
};
} // namespace structure1
