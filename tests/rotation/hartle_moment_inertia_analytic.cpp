// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_moment_inertia_analytic.cpp
 * @brief Phase 2B-4B — self-contained validation of the SCALE-FREE first-order Hartle
 *        moment of inertia I = J/Omega. No external EOS assets.
 *
 * WHAT IS AND IS NOT CLAIMED. This validates I as a scale-free observable. It makes NO
 * claim about the absolute normalization of omegabar, about physical Omega or J, about the
 * meaning of `init_omega_bar`, or about O(Omega^2) Hartle physics. INV-07 remains
 * UNRESOLVED for all of those. See docs/validation/HARTLE_MOMENT_INERTIA.md.
 *
 * BACKGROUND. The exact Schwarzschild constant-density interior, in geometric units (km):
 *
 *   rho0 = 3M/(4 pi R^3),      m(r) = M r^3/R^3
 *   y(r) = sqrt(1 - 2 M r^2/R^3),        y_R = sqrt(1 - 2M/R)
 *   e^lambda = 1/y
 *   e^nu     = (3/2) y_R - (1/2) y
 *   p(r)     = rho0 (y - y_R)/(3 y_R - y)
 *   nu'(r)   = 2 M r / ( R^3 y (3 y_R - y) )
 *
 * (Schwarzschild 1916; see Misner, Thorne & Wheeler §23.7 or Shapiro & Teukolsky §5.5.)
 * The fixture is synthetic and incompressible; it asserts no neutron-star property. Its
 * purpose is that the NEWTONIAN limit of a uniform sphere is known exactly:
 *
 *   I / (M R^2) -> 2/5   as   M/R -> 0.
 *
 * The reference solver (tests/rotation/hartle_reference.hpp) is independent of production:
 * different state variables (conservative flux q = r^4 j omegabar'), different right-hand
 * side, built from the metric rather than from production's 1/(1-2m/r) coefficient helpers.
 */

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "hartle_reference.hpp"

#include <Zaki/Physics/Constants.hpp>

using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}
static double Rel(double a, double b)
{
	return std::fabs(b) > 0.0 ? std::fabs(a - b) / std::fabs(b) : std::fabs(a - b);
}
static std::string Sci(double v, int p = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", p, v);
	return b;
}

// ---------------------------------------------------------------------------
//  Exact constant-density interior, tabulated in the geometric units StarProfile uses.
// ---------------------------------------------------------------------------
struct Uniform
{
	double M_km = 0, R_km = 0, rho0 = 0; // rho0 in km^-2
	double yR() const { return std::sqrt(1.0 - 2.0 * M_km / R_km); }
	double y(double r) const { return std::sqrt(1.0 - 2.0 * M_km * r * r / (R_km * R_km * R_km)); }
	double m(double r) const { return M_km * r * r * r / (R_km * R_km * R_km); }
	double lambda(double r) const { return -std::log(y(r)); }
	double expnu(double r) const { return 1.5 * yR() - 0.5 * y(r); }
	double nu(double r) const { return std::log(expnu(r)); }
	double p(double r) const { return rho0 * (y(r) - yR()) / (3.0 * yR() - y(r)); }
	double nuprime(double r) const
	{
		return 2.0 * M_km * r / (R_km * R_km * R_km * y(r) * (3.0 * yR() - y(r)));
	}
};

static Uniform MakeUniform(double M_km, double R_km)
{
	Uniform u;
	u.M_km = M_km;
	u.R_km = R_km;
	u.rho0 = 3.0 * M_km / (4.0 * M_PI * R_km * R_km * R_km);
	return u;
}

/// Production TOV profiles begin at r_min = 1e-5 km, never at r = 0 (dividing by r(0)
/// is a hard error inside the production profile machinery). The analytic tabulation
/// mirrors that so both solvers see the same admissible grid.
static constexpr double kRStart_km = 1.0e-5;

static double GridR(const Uniform &u, std::size_t i, std::size_t N)
{
	const double h = (u.R_km - kRStart_km) / static_cast<double>(N - 1);
	return kRStart_km + static_cast<double>(i) * h;
}

static hartle_ref::Background TabulateAnalytic(const Uniform &u, std::size_t N)
{
	hartle_ref::Background b;
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = GridR(u, i, N);
		b.r.push_back(r);
		b.m.push_back(u.m(r));
		b.p.push_back(u.p(r));
		b.eps.push_back(u.rho0);
		b.nu.push_back(u.nu(r));
		b.lambda.push_back(u.lambda(r));
	}
	return b;
}

/// Build the SAME analytic star as a production NStar, through the public
/// NStar(TOVPoint vector) constructor, so production Hartle can be run on it.
/// TOVPoint units: r [km], m [Msun], nu_der [1/cm], p [dyne/cm^2], e [g/cm^3].
static std::unique_ptr<NStar> ProductionStar(const Uniform &u, std::size_t N)
{
	const double km2_to_gcm3 =
		Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;
	const double km2_to_dyn =
		Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_INV_KM2;

	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = GridR(u, i, N);
		pts.emplace_back(r,
						 u.m(r) / Zaki::Physics::SUN_M_KM, // Msun
						 u.nuprime(r) / 1.0e5,			   // 1/km -> 1/cm
						 0.0,							   // nu: rebuilt by BuildFromTOV
						 u.p(r) * km2_to_dyn,			   // dyne/cm^2
						 u.rho0 * km2_to_gcm3,			   // g/cm^3
						 0.1,							   // n_B: inert here
						 std::vector<double>{});
	}
	return std::make_unique<NStar>(pts);
}

int main()
{
	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 2B-4B — scale-free first-order Hartle I = J/Omega\n"
				 "Self-contained: exact Schwarzschild constant-density interior.\n"
				 "Validates I ONLY. Absolute omegabar / Omega / J normalization (INV-07)\n"
				 "is NOT addressed here. No O(Omega^2) claim.\n\n";

	// A canonical moderately relativistic case for the structural checks.
	const Uniform ref_u = MakeUniform(/*M_km=*/2.0, /*R_km=*/13.0); // M/R ~ 0.154
	const auto bg = TabulateAnalytic(ref_u, 4001);

	// -----------------------------------------------------------------------
	// A0 — the background itself
	// -----------------------------------------------------------------------
	std::cout << "A0 analytic background sanity\n";
	{
		Report("A0a j(R) = 1, i.e. nu(R) + lambda(R) = 0 (Schwarzschild matching)",
			   Rel(bg.j_at(bg.R()), 1.0) < 1e-12,
			   "j(R) = " + Sci(bg.j_at(bg.R()), 12));
		Report("A0b surface pressure vanishes", std::fabs(bg.p.back()) < 1e-14 * ref_u.rho0,
			   "p(R)/rho0 = " + Sci(bg.p.back() / ref_u.rho0));
		Report("A0c sub-Buchdahl", 2.0 * ref_u.M_km / ref_u.R_km < 8.0 / 9.0,
			   "2M/R = " + Sci(2.0 * ref_u.M_km / ref_u.R_km));
	}

	// -----------------------------------------------------------------------
	// A1 — multi-seed invariance of I  (the central scale-freedom claim)
	// -----------------------------------------------------------------------
	std::cout << "\nA1 normalization invariance of I over six decades of central amplitude\n";
	{
		const std::vector<double> seeds = {1.0e-3, 5.0e-3, 1.0, 1.0e3};
		std::vector<hartle_ref::Result> rs;
		const double r0 = bg.r[1];
		for (double s : seeds)
			rs.push_back(hartle_ref::Solve(bg, s, r0));
		bool all_ok = true;
		std::cout << "      seed        omegabar(R)     J             Omega         I=J/Omega\n";
		for (std::size_t i = 0; i < seeds.size(); ++i)
		{
			all_ok = all_ok && rs[i].ok;
			if (!rs[i].ok) continue;
			std::cout << "   " << Sci(seeds[i]) << "  " << Sci(rs[i].omega_bar_R)
					  << "  " << Sci(rs[i].J) << "  " << Sci(rs[i].Omega) << "  "
					  << Sci(rs[i].I_surface, 10) << "\n";
		}
		if (all_ok)
		{
			double worst_I = 0, worst_J = 0, worst_O = 0;
			for (std::size_t i = 1; i < seeds.size(); ++i)
			{
				const double k = seeds[i] / seeds[0];
				worst_I = std::max(worst_I, Rel(rs[i].I_surface, rs[0].I_surface));
				worst_J = std::max(worst_J, Rel(rs[i].J, k * rs[0].J));
				worst_O = std::max(worst_O, Rel(rs[i].Omega, k * rs[0].Omega));
			}
			Report("A1a J scales linearly with the central amplitude", worst_J < 1e-12,
				   "max deviation from exact proportionality " + Sci(worst_J));
			Report("A1b Omega scales linearly with the central amplitude", worst_O < 1e-12,
				   "max deviation from exact proportionality " + Sci(worst_O));
			Report("A1c I = J/Omega is INVARIANT across 6 decades of amplitude",
				   worst_I < 1e-12, "max relative spread " + Sci(worst_I));
		}
		else
			Report("A1 multi-seed solve", false, "a solve failed");
	}

	// -----------------------------------------------------------------------
	// A2 — centre-radius sensitivity of the reference solver
	// -----------------------------------------------------------------------
	std::cout << "\nA2 reference-solver centre sensitivity (regular start q(r0)=0, error O(r0^5))\n";
	double centre_spread = 0.0;
	{
		const std::vector<double> r0s = {bg.r[1], 1e-3, 1e-4, kRStart_km};
		std::vector<double> Is;
		for (double r0 : r0s)
		{
			const auto res = hartle_ref::Solve(bg, 1.0, r0);
			if (!res.ok) continue;
			Is.push_back(res.I_surface);
			std::cout << "      r0 = " << Sci(r0) << " km   I = " << Sci(res.I_surface, 10)
					  << " km^3\n";
		}
		double lo = Is.empty() ? 0 : Is[0], hi = lo;
		for (double v : Is) { lo = std::min(lo, v); hi = std::max(hi, v); }
		centre_spread = lo > 0 ? (hi - lo) / lo : 1.0;
		Report("A2a I is insensitive to the reference start radius over 4 decades",
			   Is.size() >= 3 && centre_spread < 1e-9,
			   "relative spread " + Sci(centre_spread));
	}

	// -----------------------------------------------------------------------
	// A3 — surface extraction vs independent volume quadrature
	// -----------------------------------------------------------------------
	std::cout << "\nA3 surface J = R^4 omegabar'(R)/6  vs  independent volume integral\n";
	{
		const auto res = hartle_ref::Solve(bg, 1.0, bg.r[1]);
		const double d = Rel(res.I_volume, res.I_surface);
		std::cout << "      I_surface = " << Sci(res.I_surface, 10)
				  << "   I_volume = " << Sci(res.I_volume, 10) << " km^3\n";
		Report("A3a the two independent extractions of I agree", res.ok && d < 1e-5,
			   "relative difference " + Sci(d));
	}

	// -----------------------------------------------------------------------
	// A4 — the Newtonian limit:  I/(M R^2) -> 2/5
	// -----------------------------------------------------------------------
	std::cout << "\nA4 weak-field limit of a uniform sphere: I/(M R^2) -> 2/5\n";
	std::cout << "      M/R        I_ref/(M R^2)   dev from 0.4    I_prod/(M R^2)   "
				 "dev from 0.4\n";
	{
		const std::vector<double> CofR = {0.15, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002};
		std::vector<double> dev_ref, dev_prod;
		for (double C : CofR)
		{
			const double R_km = 13.0;
			const Uniform u = MakeUniform(C * R_km, R_km);
			const auto b = TabulateAnalytic(u, 4001);
			const auto res = hartle_ref::Solve(b, 1.0, b.r[1]);
			const double norm = u.M_km * u.R_km * u.R_km;
			const double q_ref = res.ok ? res.I_surface / norm : 0.0;

			double q_prod = 0.0;
			auto ns = ProductionStar(u, 4001);
			if (ns)
				q_prod = ns->GetSequence().I / norm;

			dev_ref.push_back(std::fabs(q_ref - 0.4) / 0.4);
			dev_prod.push_back(std::fabs(q_prod - 0.4) / 0.4);
			std::cout << "   " << Sci(C) << "  " << Sci(q_ref, 8) << "   "
					  << Sci(dev_ref.back()) << "   " << Sci(q_prod, 8) << "   "
					  << Sci(dev_prod.back()) << "\n";
		}
		// The approach must be monotone in C and must actually reach 2/5.
		bool mono_ref = true, mono_prod = true;
		for (std::size_t i = 1; i < dev_ref.size(); ++i)
		{
			mono_ref = mono_ref && dev_ref[i] < dev_ref[i - 1];
			mono_prod = mono_prod && dev_prod[i] < dev_prod[i - 1];
		}
		Report("A4a reference I/(M R^2) approaches 2/5 monotonically as M/R -> 0",
			   mono_ref, "deviation " + Sci(dev_ref.front()) + " -> " + Sci(dev_ref.back()));
		Report("A4b reference reaches 2/5 at the weakest field", dev_ref.back() < 5e-3,
			   "deviation at M/R = 0.002 is " + Sci(dev_ref.back()));
		Report("A4c PRODUCTION I/(M R^2) approaches 2/5 monotonically", mono_prod,
			   "deviation " + Sci(dev_prod.front()) + " -> " + Sci(dev_prod.back()));
		Report("A4d PRODUCTION reaches 2/5 at the weakest field", dev_prod.back() < 5e-3,
			   "deviation at M/R = 0.002 is " + Sci(dev_prod.back()));
	}

	// -----------------------------------------------------------------------
	// A5 — production vs reference on the same analytic background
	// -----------------------------------------------------------------------
	std::cout << "\nA5 production Hartle vs independent reference on the analytic star\n";
	{
		auto ns = ProductionStar(ref_u, 4001);
		const auto res = hartle_ref::Solve(bg, 1.0, bg.r[1]);
		if (!ns || !res.ok)
			Report("A5 solve", false, "production or reference solve failed");
		else
		{
			// First confirm production rebuilt the SAME metric from nu'.
			const auto *nu_col = ns->Profile().GetMetricNu();
			double worst_nu = 0.0;
			for (std::size_t i = 0; i < static_cast<std::size_t>(nu_col->Size()); ++i)
				worst_nu = std::max(worst_nu,
									std::fabs((*nu_col)[i] - ref_u.nu(bg.r[i])));
			Report("A5a production reconstructs the analytic metric nu(r) from nu'",
				   worst_nu < 1e-9, "max |nu_prod - nu_exact| = " + Sci(worst_nu));

			const double Ip = ns->GetSequence().I;
			const double d = Rel(Ip, res.I_surface);
			std::cout << "      I_production = " << Sci(Ip, 10)
					  << "   I_reference = " << Sci(res.I_surface, 10) << " km^3\n";
			Report("A5b production I agrees with the independent reference", d < 1e-4,
				   "relative difference " + Sci(d));
		}
	}

	// -----------------------------------------------------------------------
	// A6 — reference-solver numerical floor
	// -----------------------------------------------------------------------
	std::cout << "\nA6 reference numerical floor (must be below the production/reference gap)\n";
	{
		const auto loose = hartle_ref::Solve(bg, 1.0, bg.r[1], 1e-10, 1e-12);
		const auto tight = hartle_ref::Solve(bg, 1.0, kRStart_km, 1e-14, 1e-16);
		const double floor = Rel(loose.I_surface, tight.I_surface);
		auto ns = ProductionStar(ref_u, 4001);
		const double gap = ns ? Rel(ns->GetSequence().I, tight.I_surface) : 1.0;
		std::cout << "      reference floor = " << Sci(floor)
				  << "   production/reference gap = " << Sci(gap) << "\n";
		Report("A6a reference numerical error is subdominant to the production gap",
			   floor < 0.1 * gap,
			   "floor " + Sci(floor) + " vs gap " + Sci(gap) + " (ratio " +
				   Sci(floor > 0 ? gap / floor : 0.0, 1) + ")");
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "analytic Hartle-I checks passed"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
