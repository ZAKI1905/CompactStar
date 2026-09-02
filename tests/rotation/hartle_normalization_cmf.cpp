// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_normalization_cmf.cpp
 * @brief Phase 4A — the ADR-0006 first-order normalization contract on the authenticated
 *        DS(CMF)-1_with_crust neutron-star sequence, plus the V9 slow-rotation diagnostic.
 *
 *   usage: hartle_normalization_cmf <EOS_DATA_ROOT>
 *
 * The self-contained companion `hartle_normalization_contract` proves the contract on an exact
 * analytic background. This harness re-proves the parts that could conceivably depend on the
 * background — seed invariance, requested-Omega recovery, `J = I Omega`, and the bit-identity
 * of `I` — on the four real stars the Phase-2B/3 baselines are built from, and records the
 * regime diagnostic ADR-0006 §7 item 9 asks for.
 *
 * Tolerances are ADR-0006 §7, fixed before implementation. `c = 299792.458 km/s` is the
 * independent SI-exact oracle; `Zaki::Physics::LIGHT_C_KM_S` — the constant production converts
 * with — is deliberately not imported for the conversion checks.
 *
 * V9 IS A DIAGNOSTIC, NOT A GATE. `Omega/Omega_K` is reported so the slow-rotation regime is
 * visible; ADR-0006 §10 explicitly defers governing a threshold to the O(Omega^2) work, so
 * nothing here fails on its value.
 *
 * SCOPE. The normalization CONTRACT only. The first-order physics is Phase 2B-4B
 * (`HARTLE_MOMENT_INERTIA.md`); O(Omega^2) remains an unverified candidate (INV-08) and is not
 * touched, executed or claimed here.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;

// Non-public test seam — see hartle_normalization_contract.cpp for the rationale (ADR-0006 Q2).
namespace CompactStar::Core
{
struct RotationSolverTestSeam
{
	static void SetSeed(RotationSolver &rs, double seed) { rs.seed_omega_bar_ = seed; }
	static std::size_t Solves(const NStar &ns) { return ns.rot_solver.first_order_solve_count_; }
};
} // namespace CompactStar::Core

using CompactStar::AngularVelocity;
using CompactStar::Core::NStar;
using CompactStar::Core::PhysicalFirstOrderRotation;
using CompactStar::Core::RotationSolver;
using CompactStar::Core::RotationSolverTestSeam;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

static constexpr double kIndependentLightSpeed_km_s = 299792.458; // SI-exact

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

static double SurfaceOmegaFrom(const PhysicalFirstOrderRotation &ph)
{
	const int n = static_cast<int>(ph.omega_bar.Size());
	const double R = (*ph.r_grid)[n - 1];
	return ph.omega_bar[n - 1] + R * ph.domega_bar[n - 1] / 3.0;
}

/// Build a star at a target mass through the production construction, exactly as the Phase-2B
/// Hartle harness does. No production API is added.
static std::unique_ptr<NStar> Build(const fs::path &cold, const fs::path &wrk, double M)
{
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	std::vector<TOVPoint> pts;
	std::vector<std::string> labels;
	const int n = tov.SolveToProfile(M, pts, &labels);
	if (n <= 0 || pts.size() < 4)
		return nullptr;
	return labels.empty() ? std::make_unique<NStar>(pts)
						  : std::make_unique<NStar>(pts, labels);
}

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: hartle_normalization_cmf <EOS_DATA_ROOT>\n";
		return 2;
	}
	const fs::path root = argv[1];
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "authenticated EOS missing: " << cold.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / "compactstar_hartle_norm_cmf";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 4A — ADR-0006 normalization contract on DS(CMF)-1_with_crust\n"
				 "Contract only: normalization, units, seed isolation. No O(Omega^2) claim.\n\n";

	// The three angular velocities ADR-0006 §7 item 9 names, plus a slow case.
	const std::vector<std::pair<std::string, double>> spins = {
		{"100 rad/s", 1.0e2},
		{"2 pi x 100 Hz", 2.0 * M_PI * 1.0e2},
		{"2 pi x 716 Hz", 2.0 * M_PI * 716.0}};

	struct Row
	{
		double M = 0, R = 0, I = 0, Omega_K = 0;
		std::vector<double> ratio;
	};
	std::vector<Row> rows;

	double worst_recov = 0.0, worst_J = 0.0, worst_conv = 0.0, worst_seed = 0.0;
	bool all_I_bitwise = true, all_valid = true, all_one_solve = true;

	std::cout << "B1 sequence: seed-free response, bit-identical I, contract at each spin\n";
	for (double Mt : {1.0, 1.4, 1.6, 2.0})
	{
		auto ns = Build(cold, wrk / ("M" + std::to_string(int(Mt * 10))), Mt);
		if (!ns)
		{
			Report("B1 build M=" + Sci(Mt, 1), false, "solve failed");
			continue;
		}
		const auto &resp = ns->RotationResponse();
		all_valid = all_valid && resp.valid;
		all_I_bitwise = all_I_bitwise && (resp.I == ns->GetSequence().I);
		all_one_solve = all_one_solve && (RotationSolverTestSeam::Solves(*ns) == 1);

		Row w;
		w.M = ns->GetSequence().m;
		w.R = ns->GetSequence().r;
		w.I = resp.I;
		// Keplerian scale in geometric units: Omega_K = sqrt(M/R^3), no G needed.
		const double M_km = w.M * Zaki::Physics::SUN_M_KM;
		w.Omega_K = std::sqrt(M_km / (w.R * w.R * w.R));

		for (const auto &sp : spins)
		{
			const auto ph = ns->RotationAt(AngularVelocity::FromRadPerSecond(sp.second));
			worst_recov = std::max(worst_recov, Rel(SurfaceOmegaFrom(ph), ph.Omega_geom));
			worst_J = std::max(
				worst_J, Rel(ph.J, resp.I * sp.second / kIndependentLightSpeed_km_s));
			worst_conv = std::max(
				worst_conv, Rel(ph.Omega_geom * kIndependentLightSpeed_km_s, sp.second));
			w.ratio.push_back(ph.Omega_geom / w.Omega_K);
		}

		// V1 on a real star: the seed must not reach the physical answer.
		{
			const auto omega = AngularVelocity::FromRadPerSecond(2.0 * M_PI * 716.0);
			std::vector<double> Js, Ws;
			for (double s : {1.0e-3, 5.0e-3, 1.0, 1.0e3})
			{
				RotationSolver rs;
				rs.AttachNStar(ns.get());
				RotationSolverTestSeam::SetSeed(rs, s);
				rs.FindNMomInertia();
				const auto ph = rs.FirstOrderResponse().At(omega);
				Js.push_back(ph.J);
				Ws.push_back(rs.FirstOrderResponse().I);
			}
			for (std::size_t i = 1; i < Js.size(); ++i)
			{
				worst_seed = std::max(worst_seed, Rel(Js[i], Js[0]));
				worst_seed = std::max(worst_seed, Rel(Ws[i], Ws[0]));
			}
		}

		std::cout << "   M = " << std::fixed << std::setprecision(6) << w.M << " Msun  R = "
				  << w.R << " km  I = " << std::scientific << std::setprecision(8) << w.I
				  << " km^3\n"
				  << std::scientific << std::setprecision(8);
		rows.push_back(std::move(w));
	}

	Report("B1a every star has a valid seed-free first-order response", all_valid && !rows.empty(),
		   std::to_string(rows.size()) + " stars");
	Report("V4b response.I is BIT-IDENTICAL to SeqPoint::I on every star", all_I_bitwise,
		   "the governed observable is untouched by ADR-0006 conformance");
	Report("P1c each star build performs exactly one first-order ODE solve", all_one_solve,
		   "no extra integration was introduced");
	Report("V1f the internal seed does not reach J or I on real stars", worst_seed <= 1e-10,
		   "worst relative spread over seeds [1e-3, 1e3] " + Sci(worst_seed) + " (bound 1e-10)");
	Report("V2c the requested Omega is recovered at the surface on every star and spin",
		   worst_recov <= 1e-13, "worst relative difference " + Sci(worst_recov) + " (bound 1e-13)");
	Report("V3b J = I * Omega_phys / c on every star and spin, with the independent c",
		   worst_J <= 1e-13, "worst relative difference " + Sci(worst_J) + " (bound 1e-13)");
	Report("V6c Omega_geom * c_independent reproduces the requested rad/s", worst_conv <= 1e-15,
		   "worst relative difference " + Sci(worst_conv) + " (bound 1e-15)");

	// =======================================================================
	//  V9 — slow-rotation regime diagnostic (REPORTED, never enforced)
	// =======================================================================
	std::cout << "\nV9 slow-rotation diagnostic: Omega / Omega_K with Omega_K = sqrt(M/R^3) "
				 "[geometric]\n";
	std::cout << "   M [Msun]   R [km]     Omega_K [1/km]";
	for (const auto &sp : spins)
		std::cout << "   " << std::setw(14) << sp.first;
	std::cout << "\n";
	for (const auto &w : rows)
	{
		std::cout << "   " << std::fixed << std::setprecision(4) << w.M << "     " << w.R
				  << "   " << std::scientific << std::setprecision(6) << w.Omega_K;
		for (double q : w.ratio)
			std::cout << "   " << std::setw(14) << Sci(q, 6);
		std::cout << "\n" << std::scientific << std::setprecision(8);
	}
	std::cout << "   (ADR-0006 §10 defers governing a slow-rotation threshold to the O(Omega^2)\n"
				 "    work; this table is reported, and nothing above fails on its value.)\n";
	Report("V9a the regime diagnostic is finite and ordered on every star",
		   !rows.empty() && std::all_of(rows.begin(), rows.end(),
										[](const Row &w) {
											return w.Omega_K > 0.0 && std::isfinite(w.Omega_K) &&
												   w.ratio.size() == 3 &&
												   w.ratio[0] < w.ratio[1] &&
												   w.ratio[1] < w.ratio[2];
										}),
		   "Omega/Omega_K increases with the requested spin, as it must");

	fs::remove_all(wrk);
	std::cout << "\n"
			  << (g_fail == 0 ? "CMF normalization contract checks passed"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
