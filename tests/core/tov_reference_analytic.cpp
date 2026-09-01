// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file tov_reference_analytic.cpp
 * @brief Equation-level validation of the production TOV right-hand side against the
 *        EXACT Schwarzschild constant-density interior solution.
 *
 * Self-contained: builds a synthetic constant-density EOS in a temp directory and loads it
 * through the production `TOVSolver::ImportEOS` path. No external CompOSE assets.
 *
 * WHAT IS TESTED. `TOVSolver::ODE` — the actual production right-hand side — is evaluated at
 * exact (r, m(r), p(r)) taken from the analytic solution, and its returned dm/dr and dp/dr are
 * compared against the analytic derivatives. This exercises the GR factors, the CGS unit
 * constants G and c, the signs, the (rho + p/c^2) and (m + 4 pi r^3 p/c^2) terms, and the
 * (1 - 2Gm/(c^2 r)) denominator.
 *
 * ANALYTIC REFERENCE (Schwarzschild 1916 interior solution; see e.g. Misner, Thorne & Wheeler,
 * "Gravitation", §23.7, or Shapiro & Teukolsky, "Black Holes, White Dwarfs and Neutron Stars",
 * §5.5). For constant mass density rho0, with
 *
 *     m(r) = (4 pi / 3) rho0 r^3 ,   M = (4 pi / 3) rho0 R^3 ,   y = 2 G M / (c^2 R) ,
 *     A(r) = sqrt(1 - y r^2 / R^2) ,  B = sqrt(1 - y) ,
 *
 * the pressure and its derivative are
 *
 *     p(r)     = rho0 c^2 (A - B) / (3B - A) ,
 *     dp/dr    = -2 rho0 c^2 y r B / ( R^2 A (3B - A)^2 ) ,
 *
 * the latter obtained by differentiating p(r) in closed form. The Buchdahl bound requires
 * y < 8/9; both test cases stay well below it and are genuinely relativistic.
 *
 * The fixture is SYNTHETIC and NON-PHYSICAL: a strictly incompressible star. It asserts no
 * neutron-star property.
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_const_cgsm.h>

#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::TOVSolver;

static const double G_cgs = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT;
static const double c_cgs = GSL_CONST_CGSM_SPEED_OF_LIGHT;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}

/// Synthetic incompressible EOS: eps = rho0 (constant, g/cm^3), p strictly increasing.
static void WriteConstantDensityEOS(const fs::path &f, double rho0)
{
	std::ofstream o(f);
	o << "e(g/cm^3)\tp(dyne/cm^2)\trho(1/fm^3)\n";
	o << std::setprecision(17) << std::scientific;
	const int N = 400;
	const double lp0 = std::log(1.0e22), lp1 = std::log(1.0e37);
	for (int i = 0; i < N; ++i)
	{
		const double p = std::exp(lp0 + (lp1 - lp0) * double(i) / double(N - 1));
		o << rho0 << "\t" << p << "\t" << 0.3 << "\n"; // n_B unused by the TOV RHS
	}
}

struct Exact
{
	double rho0, R, M, y;
	double A(double r) const { return std::sqrt(1.0 - y * r * r / (R * R)); }
	double B() const { return std::sqrt(1.0 - y); }
	double m(double r) const { return (4.0 / 3.0) * M_PI * rho0 * r * r * r; }
	double p(double r) const
	{
		const double a = A(r), b = B();
		return rho0 * c_cgs * c_cgs * (a - b) / (3.0 * b - a);
	}
	double dpdr(double r) const
	{
		const double a = A(r), b = B();
		return -2.0 * rho0 * c_cgs * c_cgs * y * r * b / (R * R * a * (3.0 * b - a) * (3.0 * b - a));
	}
	double dmdr(double r) const { return 4.0 * M_PI * r * r * rho0; }
};

static Exact MakeExact(double rho0, double y_target)
{
	// y = (8 pi G / 3 c^2) rho0 R^2  ->  solve for R
	const double k = 8.0 * M_PI * G_cgs / (3.0 * c_cgs * c_cgs);
	const double R = std::sqrt(y_target / (k * rho0));
	Exact e;
	e.rho0 = rho0;
	e.R = R;
	e.M = (4.0 / 3.0) * M_PI * rho0 * R * R * R;
	e.y = 2.0 * G_cgs * e.M / (c_cgs * c_cgs * R);
	return e;
}

static double Rel(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout << "TOV equation-level reference: exact Schwarzschild constant-density interior\n"
			  << "(synthetic incompressible fixture; asserts no neutron-star property)\n\n";

	const fs::path root = fs::temp_directory_path() / "compactstar_tov_analytic";
	fs::remove_all(root);
	fs::create_directories(root);
	const fs::path eosf = root / "const_density.eos";

	const double rho0 = 5.0e14; // g/cm^3
	WriteConstantDensityEOS(eosf, rho0);

	TOVSolver tov;
	tov.SetWrkDir(root.string());
	tov.ImportEOS(eosf.string(), true);

	// Sanity: the production EOS lookup must return the constant density.
	{
		const double e1 = tov.GetEDens(1.0e30);
		const double e2 = tov.GetEDens(1.0e34);
		Report("E0 production GetEDens returns the constant fixture density",
			   Rel(e1, rho0) < 1e-10 && Rel(e2, rho0) < 1e-10,
			   "eps(1e30)=" + std::to_string(e1) + ", eps(1e34)=" + std::to_string(e2));
	}

	// Two genuinely relativistic compactness values, both safely under Buchdahl 8/9.
	for (double y_target : {0.30, 0.50})
	{
		const Exact ex = MakeExact(rho0, y_target);
		std::cout << "\ncompactness 2GM/(Rc^2) = " << ex.y
				  << "   R = " << ex.R / 1.0e5 << " km"
				  << "   M = " << ex.M / GSL_CONST_CGSM_SOLAR_MASS << " Msun"
				  << "   (Buchdahl bound 8/9 = " << (8.0 / 9.0) << ")\n";

		double worst_m = 0.0, worst_p = 0.0;
		// center-near, mid-radius, near-surface. r=0 is avoided: production uses a
		// finite-radius center convention.
		for (double frac : {0.02, 0.10, 0.25, 0.50, 0.75, 0.90, 0.98})
		{
			const double r = frac * ex.R;
			double yv[2] = {ex.m(r), ex.p(r)};
			double f[2] = {0.0, 0.0};
			const int rc = TOVSolver::ODE(r, yv, f, &tov);
			if (rc != 0)
			{
				Report("ODE returned success", false,
					   "GSL status " + std::to_string(rc) + " at r/R=" + std::to_string(frac));
				continue;
			}
			const double em = Rel(f[0], ex.dmdr(r));
			const double ep = Rel(f[1], ex.dpdr(r));
			worst_m = std::max(worst_m, em);
			worst_p = std::max(worst_p, ep);
			std::cout << "    r/R=" << frac
					  << "  dm/dr rel=" << em << "  dp/dr rel=" << ep << "\n";
		}
		Report("A1 production dm/dr == 4 pi r^2 rho (y=" + std::to_string(y_target) + ")",
			   worst_m < 1e-10, "max rel " + std::to_string(worst_m));
		Report("A2 production dp/dr == exact Schwarzschild derivative (y=" +
				   std::to_string(y_target) + ")",
			   worst_p < 1e-8, "max rel " + std::to_string(worst_p));
	}

	// Sign/GR sanity: dp/dr must be negative, and the Newtonian limit must NOT reproduce
	// the relativistic value (i.e. the GR correction factors are actually active).
	{
		const Exact ex = MakeExact(rho0, 0.50);
		const double r = 0.5 * ex.R;
		double yv[2] = {ex.m(r), ex.p(r)};
		double f[2] = {0.0, 0.0};
		TOVSolver::ODE(r, yv, f, &tov);
		const double newtonian = -G_cgs * ex.m(r) * rho0 / (r * r);
		Report("A3 dp/dr is negative", f[1] < 0.0, "dp/dr = " + std::to_string(f[1]));
		Report("A4 GR corrections are active (differs from Newtonian by >10%)",
			   Rel(f[1], newtonian) > 0.10,
			   "GR/Newtonian = " + std::to_string(f[1] / newtonian));
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "production TOV RHS matches the exact interior solution"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	fs::remove_all(root);
	return g_fail == 0 ? 0 : 1;
}
