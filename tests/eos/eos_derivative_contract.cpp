#include "tests/relativity/fixture_units.hpp"
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file eos_derivative_contract.cpp
 * @brief Phase 4C-I0 — contract for the EOS thermodynamic derivative `d(eps)/dp`
 *        (ADR-0007 P5, ACCEPTED 2026-09-02). Self-contained: builds its own synthetic EOS
 *        tables in a temp directory and needs NO external assets.
 *
 * WHAT IS UNDER TEST. ADR-0007 P5 gives the EOS/TOV layer sole authority over `d(eps)/dp`,
 * and fixes it as the derivative of the *same* `eps(p)` interpolant that constructed the star.
 * `TOVSolver::GetEDensDeriv` is that authority. This file asserts the four things a caller
 * must be able to rely on:
 *
 *   1. it really is the derivative of that interpolant (not a re-interpolation, not a
 *      finite difference of anything);
 *   2. it is delivered DIMENSIONLESS, with the cgs -> geometric factor applied exactly once;
 *   3. its domain semantics FAIL CLOSED — a derivative is a local property of the state, so
 *      an out-of-range query must not be answered with a boundary value;
 *   4. a profile publishes the derivative only when every node has an authoritative one, and
 *      a point-constructed analytic star can supply its own explicitly.
 *
 * WHY A SYNTHETIC EOS. The oracle has to be independent of production. Here it is a closed-form
 * `eps(p)` that this file both tabulates and differentiates analytically, so the comparison
 * never consults CompactStar for the answer.
 *
 * ============================== PREDECLARED BOUNDS ==============================
 * Fixed BEFORE any comparison was run, from the interpolation order alone. GSL's Steffen
 * interpolant (`TOV_gsl_interp_type`, TOVSolver.hpp:488) is a monotonicity-preserving cubic
 * Hermite: C^1, with slopes accurate to O(h^2) on smooth data and one-sided at the endpoints.
 *
 *   L  linear reproduction        rel <= 1e-12   a cubic Hermite with the common secant slope
 *                                                reproduces collinear data exactly; only
 *                                                floating point separates the two.
 *   S1 smooth interior, N = 200   rel <= 1e-3    O(h^2) with h/p = ln(100)/200 = 2.3e-2.
 *   S2 refinement N -> 2N         factor >= 3    O(h^2) predicts 4; 3 leaves room for the
 *                                                limiter, which is not second-order accurate
 *                                                where it engages.
 *   S3 endpoints (one-sided)      rel <= 1e-2    strictly weaker than the interior, and
 *                                                reported separately rather than hidden.
 *   F  vs central FD of the SAME  rel <= 1e-6    the FD is a TEST ORACLE ONLY, never an
 *      spline, away from knots                   authority; step 1e-5*p, ~1e-3 of the local
 *                                                knot spacing, so it never straddles a knot.
 *   U  unit identity              rel <= 1e-13   exact-by-definition constant vs the repo's.
 * ===============================================================================
 *
 * NOT TESTED HERE: any O(Omega^2) physics. 4C-I0 adds the derivative authority and nothing
 * else; the monopole solver is 4C-I1 and its validation is 4D (INV-08).
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;

using CompactStar::Core::NStar;
using CompactStar::Core::StarProfile;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

// ---------------------------------------------------------------------------
//  Predeclared bounds — see the header comment. Nothing below revises them.
// ---------------------------------------------------------------------------
static constexpr double kBoundLinearExact = 1.0e-12;
static constexpr double kBoundSmoothInterior = 1.0e-3;
static constexpr double kRefinementFactorMin = 3.0;
static constexpr double kBoundSmoothEndpoint = 1.0e-2;
static constexpr double kBoundVsFiniteDiff = 1.0e-6;
static constexpr double kBoundUnitIdentity = 1.0e-13;

/// The exact defining value of the speed of light, in cm/s. An INDEPENDENT literal: this test
/// must not take the repository's own constant as its oracle for the unit conversion.
static constexpr double kLightC_cm_per_s = 2.99792458e10;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &detail)
{
	std::cout << (ok ? "  [ OK ] " : "  [FAIL] ") << id;
	if (!detail.empty())
		std::cout << "   " << detail;
	std::cout << "\n";
	if (!ok)
		++g_fail;
}

static std::string Sci(double v, int p = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", p, v);
	return b;
}

static double Rel(double a, double b)
{
	const double d = std::fabs(a - b);
	return std::fabs(b) > 0.0 ? d / std::fabs(b) : d;
}

// ---------------------------------------------------------------------------
//  The synthetic equations of state. Both are written in the table's own units:
//  eps in g/cm^3 against p in dyne/cm^2, exactly like the authenticated CompOSE input.
// ---------------------------------------------------------------------------

/// Linear: eps = E0 + S*(p - P_MIN). Collinear data, so the interpolant must reproduce it.
struct LinearEos
{
	double E0 = 1.0e14;	  // g/cm^3 at p = P_MIN
	double S = 5.0e-21;	  // (g/cm^3)/(dyne/cm^2)
	double p_min = 1.0e33; // dyne/cm^2
	double Eps(double p) const { return E0 + S * (p - p_min); }
	double DEpsDp(double) const { return S; }
};

/// Power law: eps = E0 * (p/P0)^gamma, gamma < 1, so d(eps)/dp decreases with pressure — the
/// qualitative behaviour of real cold matter, and with a genuinely nonzero third derivative,
/// which is what makes the O(h^2) statement testable at all.
struct PowerEos
{
	double E0 = 1.0e15;	 // g/cm^3
	double P0 = 1.0e35;	 // dyne/cm^2
	double gamma = 0.40; // dimensionless
	double Eps(double p) const { return E0 * std::pow(p / P0, gamma); }
	double DEpsDp(double p) const { return E0 * gamma * std::pow(p / P0, gamma - 1.0) / P0; }
};

/// Baryon density is irrelevant to this contract but the importer requires a third column;
/// a smooth monotone stand-in keeps the table well-formed.
static double RhoOf(double p, double p_min, double p_max)
{
	const double t = std::log(p / p_min) / std::log(p_max / p_min);
	return 1.0e-7 * std::pow(1.0e7, t); // 1e-7 .. 1.0 fm^-3
}

template <class Eos>
static std::vector<double> LogGrid(double p_min, double p_max, std::size_t n)
{
	std::vector<double> p(n);
	for (std::size_t i = 0; i < n; ++i)
		p[i] = p_min * std::pow(p_max / p_min,
								static_cast<double>(i) / static_cast<double>(n - 1));
	return p;
}

template <class Eos>
static bool WriteEosFile(const fs::path &path, const Eos &eos,
						 const std::vector<double> &p_nodes)
{
	std::ofstream out(path);
	if (!out)
		return false;
	out << "e(g/cm^3)\tp(dyne/cm^2)\trho(1/fm^3)\n";
	out << std::scientific << std::setprecision(16);
	const double p_min = p_nodes.front(), p_max = p_nodes.back();
	for (double p : p_nodes)
		out << eos.Eps(p) << "\t" << p << "\t" << RhoOf(p, p_min, p_max) << "\n";
	return true;
}

// ---------------------------------------------------------------------------
//  An exact constant-density interior, used only for the profile-side contract.
//  Schwarzschild 1916; the same fixture the Phase-2B/4B rotation harnesses use.
// ---------------------------------------------------------------------------
struct UniformStar
{
	double M_km = 2.0, R_km = 13.0, rho0 = 0.0; // rho0 in km^-2
	double yR() const { return std::sqrt(1.0 - 2.0 * M_km / R_km); }
	double y(double r) const { return std::sqrt(1.0 - 2.0 * M_km * r * r / (R_km * R_km * R_km)); }
	double m(double r) const { return M_km * r * r * r / (R_km * R_km * R_km); }
	double p(double r) const { return rho0 * (y(r) - yR()) / (3.0 * yR() - y(r)); }
	double nuprime(double r) const
	{
		return 2.0 * M_km * r / (R_km * R_km * R_km * y(r) * (3.0 * yR() - y(r)));
	}
};

/// Build the analytic star as TOVPoints. `dedp` semantics:
///   mode 0 — do not supply one (the default NaN): the profile must publish nothing;
///   mode 1 — supply the analytic value 0.0 at every node (incompressible matter);
///   mode 2 — supply it everywhere except one interior node, to prove all-or-nothing.
static std::vector<TOVPoint> UniformPoints(const UniformStar &u, std::size_t N, int mode)
{
	const double km2_to_gcm3 =
		relativity_fixture::eps_to_rho;
	const double km2_to_dyn =
		relativity_fixture::pressure_to_cgs;
	const double r0 = 1.0e-5;

	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = r0 + static_cast<double>(i) * (u.R_km - r0) /
									static_cast<double>(N - 1);
		double dedp = std::numeric_limits<double>::quiet_NaN();
		if (mode == 1)
			dedp = 0.0; // exact: constant density has d(eps)/dp = 0
		else if (mode == 2)
			dedp = (i == N / 2) ? std::numeric_limits<double>::quiet_NaN() : 0.0;

		pts.emplace_back(r, u.m(r) / relativity_fixture::solar_km, u.nuprime(r) / 1.0e5, 0.0,
						 u.p(r) * km2_to_dyn, u.rho0 * km2_to_gcm3, 0.1,
						 std::vector<double>{}, dedp);
	}
	return pts;
}

//==============================================================
int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4C-I0 — EOS thermodynamic-derivative contract (ADR-0007 P5).\n"
				 "Self-contained. Bounds predeclared from the interpolation order.\n"
				 "No O(Omega^2) physics is exercised or claimed (INV-08).\n\n";

	const fs::path wrk = fs::temp_directory_path() / "compactstar_eos_deriv_contract";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	// ------------------------------------------------------------------
	//  U — the unit identity, checked against an independent literal.
	// ------------------------------------------------------------------
	std::cout << "U. Unit contract: the cgs -> dimensionless factor\n";
	{
		const double repo_ratio =
			Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_G_CM3;
		const double c2 = kLightC_cm_per_s * kLightC_cm_per_s;
		const double rel = Rel(repo_ratio, c2);
		Report("Ua INV_FM4_2_Dyn_CM2 / INV_FM4_2_G_CM3 == c^2 (independent literal)",
			   rel <= kBoundUnitIdentity,
			   "ratio=" + Sci(repo_ratio, 12) + "  c^2=" + Sci(c2, 12) +
				   "  rel=" + Sci(rel) + "  bound=" + Sci(kBoundUnitIdentity));
		std::cout << "     (eps is tabulated as a MASS density g/cm^3 while p is an ENERGY\n"
					 "      density dyne/cm^2 = erg/cm^3, so the raw spline derivative carries\n"
					 "      s^2/cm^2 and the dimensionless value is that times c^2.)\n";
	}

	// ------------------------------------------------------------------
	//  Pre-import fail-closed behaviour.
	// ------------------------------------------------------------------
	std::cout << "\nD. Domain semantics — fail closed\n";
	{
		TOVSolver bare;
		const bool has = bare.HasEDensDeriv();
		const double v = bare.GetEDensDeriv(1.0e34);
		Report("Da a solver with no EOS reports no derivative and returns NaN",
			   !has && std::isnan(v), "HasEDensDeriv=" + std::string(has ? "true" : "false"));
	}

	// ------------------------------------------------------------------
	//  L — linear EOS: exact reproduction.
	// ------------------------------------------------------------------
	std::cout << "\nL. Linear EOS — the interpolant must reproduce collinear data exactly\n";
	{
		LinearEos eos;
		const double p_min = 1.0e33, p_max = 1.0e35;
		eos.p_min = p_min;
		const auto nodes = LogGrid<LinearEos>(p_min, p_max, 120);
		const fs::path f = wrk / "linear.eos";
		if (!WriteEosFile(f, eos, nodes))
			Report("L fixture written", false, "cannot write " + f.string());
		else
		{
			TOVSolver tov;
			tov.SetWrkDir(wrk.string());
			tov.ImportEOS(f.string(), true);

			Report("La the interpolant exists and reports its own domain",
				   tov.HasEDensDeriv() &&
					   Rel(tov.EDensDerivPressMin(), p_min) <= 1e-15 &&
					   Rel(tov.EDensDerivPressMax(), p_max) <= 1e-15,
				   "[" + Sci(tov.EDensDerivPressMin()) + ", " +
					   Sci(tov.EDensDerivPressMax()) + "]");

			const double expect = eos.DEpsDp(0.0) * kLightC_cm_per_s * kLightC_cm_per_s;
			double worst = 0.0, worst_p = 0.0;
			// knots, inter-knot midpoints, and both exact endpoints
			std::vector<double> probes;
			for (std::size_t i = 0; i < nodes.size(); ++i)
			{
				probes.push_back(nodes[i]);
				if (i + 1 < nodes.size())
					probes.push_back(std::sqrt(nodes[i] * nodes[i + 1]));
			}
			bool all_finite = true;
			for (double p : probes)
			{
				const double got = tov.GetEDensDeriv(p);
				if (!std::isfinite(got))
				{
					all_finite = false;
					break;
				}
				const double rel = Rel(got, expect);
				if (rel > worst)
				{
					worst = rel;
					worst_p = p;
				}
			}
			Report("Lb every knot, midpoint and endpoint returns the exact analytic value",
				   all_finite && worst <= kBoundLinearExact,
				   "expect=" + Sci(expect, 12) + "  worst rel=" + Sci(worst) +
					   " at p=" + Sci(worst_p) + "  bound=" + Sci(kBoundLinearExact) +
					   "  n_probes=" + std::to_string(probes.size()));

			// --- domain rejection, on a live interpolant ---
			const double below = p_min * (1.0 - 1.0e-9);
			const double above = p_max * (1.0 + 1.0e-9);
			Report("Db below the domain -> NaN (NOT the boundary derivative)",
				   std::isnan(tov.GetEDensDeriv(below)), "p=" + Sci(below));
			Report("Dc above the domain -> NaN", std::isnan(tov.GetEDensDeriv(above)),
				   "p=" + Sci(above));
			Report("Dd the exact endpoints are IN the domain and evaluate",
				   std::isfinite(tov.GetEDensDeriv(p_min)) &&
					   std::isfinite(tov.GetEDensDeriv(p_max)),
				   "xmin=" + Sci(tov.GetEDensDeriv(p_min)) +
					   "  xmax=" + Sci(tov.GetEDensDeriv(p_max)));
			const double nan_in = std::numeric_limits<double>::quiet_NaN();
			const double inf_in = std::numeric_limits<double>::infinity();
			Report("De non-finite input -> NaN (NaN, +inf, -inf)",
				   std::isnan(tov.GetEDensDeriv(nan_in)) &&
					   std::isnan(tov.GetEDensDeriv(inf_in)) &&
					   std::isnan(tov.GetEDensDeriv(-inf_in)),
				   "");
			std::cout << "     (GetEDens CLAMPS and warns; the derivative REFUSES. A derivative\n"
						 "      is a local property of the state, so a clamped answer would\n"
						 "      describe a different state than the one asked about.)\n";
		}
	}

	// ------------------------------------------------------------------
	//  S — smooth power-law EOS: discretization behaviour and FD cross-check.
	// ------------------------------------------------------------------
	std::cout << "\nS. Smooth power-law EOS — discretization behaviour of the SAME interpolant\n";
	{
		PowerEos eos;
		const double p_min = 1.0e33, p_max = 1.0e35;
		const double c2 = kLightC_cm_per_s * kLightC_cm_per_s;

		auto measure = [&](std::size_t n, double &interior, double &endpoint,
						   double &fd_worst) -> bool {
			const auto nodes = LogGrid<PowerEos>(p_min, p_max, n);
			const fs::path f = wrk / ("power_" + std::to_string(n) + ".eos");
			if (!WriteEosFile(f, eos, nodes))
				return false;
			TOVSolver tov;
			tov.SetWrkDir(wrk.string());
			tov.ImportEOS(f.string(), true);
			if (!tov.HasEDensDeriv())
				return false;

			interior = 0.0;
			endpoint = 0.0;
			fd_worst = 0.0;
			for (std::size_t i = 0; i < nodes.size(); ++i)
			{
				const double p = nodes[i];
				const double got = tov.GetEDensDeriv(p);
				if (!std::isfinite(got))
					return false;
				const double rel = Rel(got, eos.DEpsDp(p) * c2);
				if (i == 0 || i + 1 == nodes.size())
					endpoint = std::max(endpoint, rel);
				else
					interior = std::max(interior, rel);
			}

			// Central FD of the SAME spline, at knot midpoints so the stencil never straddles
			// a knot (Steffen is only C^1, so a straddling stencil would measure the kink, not
			// the derivative). TEST ORACLE ONLY — never a production authority.
			for (std::size_t i = 0; i + 1 < nodes.size(); ++i)
			{
				const double pm = std::sqrt(nodes[i] * nodes[i + 1]);
				const double d = 1.0e-5 * pm;
				const double e_hi = tov.GetEDens(pm + d);
				const double e_lo = tov.GetEDens(pm - d);
				const double fd = ((e_hi - e_lo) / (2.0 * d)) * c2;
				const double got = tov.GetEDensDeriv(pm);
				fd_worst = std::max(fd_worst, Rel(got, fd));
			}
			return true;
		};

		double i200 = 0, e200 = 0, fd200 = 0, i400 = 0, e400 = 0, fd400 = 0;
		const bool ok200 = measure(200, i200, e200, fd200);
		const bool ok400 = measure(400, i400, e400, fd400);

		Report("Sa both resolutions evaluate everywhere without a non-finite value",
			   ok200 && ok400, "");
		Report("Sb N=200 interior nodes are within the predeclared O(h^2) bound",
			   ok200 && i200 <= kBoundSmoothInterior,
			   "max rel=" + Sci(i200) + "  bound=" + Sci(kBoundSmoothInterior));
		const double factor = (i400 > 0.0) ? i200 / i400
										   : std::numeric_limits<double>::infinity();
		Report("Sc doubling the node count reduces the interior error like O(h^2)",
			   ok200 && ok400 && factor >= kRefinementFactorMin,
			   "N=200: " + Sci(i200) + "  N=400: " + Sci(i400) +
				   "  factor=" + Sci(factor, 2) + "  required>=" +
				   Sci(kRefinementFactorMin, 1));
		Report("Sd endpoints (one-sided slopes) stay within their own weaker bound",
			   ok200 && ok400 && e200 <= kBoundSmoothEndpoint &&
				   e400 <= kBoundSmoothEndpoint,
			   "N=200: " + Sci(e200) + "  N=400: " + Sci(e400) +
				   "  bound=" + Sci(kBoundSmoothEndpoint));
		Report("Se the returned value IS the derivative of the star's own eps(p) spline "
			   "(central FD of that spline, away from knots)",
			   ok200 && ok400 && fd200 <= kBoundVsFiniteDiff &&
				   fd400 <= kBoundVsFiniteDiff,
			   "N=200: " + Sci(fd200) + "  N=400: " + Sci(fd400) +
				   "  bound=" + Sci(kBoundVsFiniteDiff));
		std::cout << "     (Se is the strongest statement here: it pins the value to the SAME\n"
					 "      interpolant GetEDens reads, not merely to the analytic EOS. A second\n"
					 "      interpolant, or a profile finite difference, would fail it.)\n";
	}

	// ------------------------------------------------------------------
	//  P — the profile side: all-or-nothing publication, explicit analytic supply.
	// ------------------------------------------------------------------
	std::cout << "\nP. Profile ownership — point-constructed stars and fail-closed absence\n";
	{
		UniformStar u;
		u.rho0 = 3.0 * u.M_km / (4.0 * M_PI * u.R_km * u.R_km * u.R_km);
		const std::size_t N = 601;

		auto pts_none = UniformPoints(u, N, 0);
		auto pts_full = UniformPoints(u, N, 1);
		auto pts_hole = UniformPoints(u, N, 2);

		NStar ns_none(pts_none);
		NStar ns_full(pts_full);
		NStar ns_hole(pts_hole);

		Report("Pa a star whose points supply no derivative publishes none (fail closed)",
			   !ns_none.Profile().HasEosDEdP() &&
				   ns_none.Profile().GetEosDEdP() == nullptr,
			   "HasEosDEdP=false, GetEosDEdP=nullptr");

		const auto *col = ns_full.Profile().GetEosDEdP();
		bool full_ok = ns_full.Profile().HasEosDEdP() && col != nullptr;
		bool all_zero = full_ok;
		if (full_ok)
		{
			full_ok = (col->Size() == N);
			for (std::size_t i = 0; full_ok && i < col->Size(); ++i)
				if (!((*col)[static_cast<int>(i)] == 0.0))
					all_zero = false;
		}
		Report("Pb an analytic star may supply its own EOS derivative explicitly",
			   full_ok, "size=" + std::to_string(col ? col->Size() : 0) +
							"  expected=" + std::to_string(N));
		Report("Pc the supplied constant-density value is exactly 0 at every node "
			   "(incompressible: d(eps)/dp = 0, NOT the causal limit 1)",
			   full_ok && all_zero, "");

		Report("Pd one missing node withdraws the whole set — a partial set is no authority",
			   !ns_hole.Profile().HasEosDEdP() &&
				   ns_hole.Profile().GetEosDEdP() == nullptr,
			   "one interior NaN out of " + std::to_string(N));

		// First-order behaviour must be untouched by the presence or absence of the derivative.
		const double I_none = ns_none.GetSequence().I;
		const double I_full = ns_full.GetSequence().I;
		const auto &R_none = ns_none.RotationResponse();
		const auto &R_full = ns_full.RotationResponse();
		bool shape_bitwise = R_none.valid && R_full.valid &&
							 R_none.omega_bar_over_Omega.Size() ==
								 R_full.omega_bar_over_Omega.Size();
		for (int i = 0; shape_bitwise && i < static_cast<int>(R_none.omega_bar_over_Omega.Size());
			 ++i)
		{
			if (R_none.omega_bar_over_Omega[i] != R_full.omega_bar_over_Omega[i] ||
				R_none.domega_bar_over_Omega_dr[i] != R_full.domega_bar_over_Omega_dr[i])
				shape_bitwise = false;
		}
		Report("Pe first-order results are BIT-IDENTICAL with and without the derivative",
			   I_none == I_full && shape_bitwise,
			   "I=" + Sci(I_none, 12) + " both; s and s' bitwise equal");
		std::cout << "     (4C-I0 adds an authority; it changes no physics. First-order and TOV\n"
					 "      construction must never begin to require a derivative they do not\n"
					 "      use — only the future O(Omega^2) solver may demand it.)\n";
	}

	// ------------------------------------------------------------------
	std::cout << "\nc_s^2 cross-check: NOT AVAILABLE IN CURRENT GOVERNED INPUT.\n"
				 "  The governed import path reads e(g/cm^3), p(dyne/cm^2), rho(1/fm^3) and\n"
				 "  species-fraction columns only (TOVSolver.cpp:509-528). No sound-speed or\n"
				 "  thermodynamic-derivative column is imported, and 4C-I0 does not broaden the\n"
				 "  import path to add one (ADR-0007 P5: one source of truth).\n";

	fs::remove_all(wrk);

	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail
			  << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
