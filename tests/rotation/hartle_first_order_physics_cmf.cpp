// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_first_order_physics_cmf.cpp
 * @brief Phase 4B — INDEPENDENT physical validation of the normalized first-order Hartle
 *        response on the authenticated DS(CMF)-1_with_crust sequence.
 *
 *   usage: hartle_first_order_physics_cmf <EOS_DATA_ROOT>
 *
 * The self-contained companion proves the profile on an exact analytic background, where the
 * background itself carries no discretization error. This harness repeats the comparison on the
 * four real stars the Phase-2B/3 baselines are built from, where the background is a tabulated
 * TOV solution — the case that actually matters for published physics.
 *
 * The independent oracle is `hartle_reference.hpp`, integrating the CONSERVATIVE system in
 * different state variables from the metric columns, normalized by its OWN surface extraction.
 * It never calls production's ODE, its coefficient helpers, or `HartleFirstOrderResponse::At`.
 *
 * ============================ PREDECLARED BOUNDS ============================
 * Fixed BEFORE the comparison, from the published Phase-2B record only
 * (`docs/validation/HARTLE_MOMENT_INERTIA.md` §10, §8, §11), by the same rule the analytic
 * harness states:
 *
 *     profile bound = 5 x (Phase-2B production/reference `I` agreement for this star class),
 *                     rounded up to the next power of ten.
 *
 *   CMF sequence : 2B measured 1.32e-5 … 2.15e-5 -> 5 x 2.15e-5 = 1.08e-4 -> BOUND 1e-4
 *   volume id    : 2B surf/vol 6.62e-7 … 8.05e-7, plus the same I agreement
 *                                                -> 5 x (8.05e-7 + 2.15e-5) -> BOUND 1e-4
 *   exterior     : algebraic identity, roundoff only               -> BOUND 1e-12
 *
 * None of these depends on a Phase-4B measurement, and none is widened.
 *
 * SCOPE. First order only. No O(Omega^2) claim (INV-08). The `Omega/Omega_K` diagnostic is
 * context, not a validated slow-rotation error budget.
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
#include "hartle_profile_compare.hpp"
#include "hartle_reference.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;
using CompactStar::AngularVelocity;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;
using namespace hartle_4b;

// ---- predeclared bounds (see file header) ---------------------------------------------
static constexpr double kBoundProfileCMF = 1.0e-4;
static constexpr double kBoundVolumeCMF = 1.0e-4;
static constexpr double kBoundExteriorIdentity = 1.0e-12;
static constexpr double kRefFloorMargin = 1.0e-3;

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

struct Row
{
	double M = 0, R = 0, I_prod = 0, I_ref = 0, I_vol = 0;
	double max_rel_s = 0, rms_s = 0, r_worst_s = 0;
	double sp_abs_peak = 0, sp_max_rel = 0, rms_sp = 0, r_worst_sp = 0;
	double drag_c = 0, drag_R = 0, Omega_K = 0;
	std::size_t n = 0;
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: hartle_first_order_physics_cmf <EOS_DATA_ROOT>\n";
		return 2;
	}
	const fs::path root = argv[1];
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "authenticated EOS missing: " << cold.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / "compactstar_hartle_4b_cmf";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 4B — independent physical validation of the normalized first-order\n"
				 "Hartle response on DS(CMF)-1_with_crust. Bounds predeclared from Phase 2B.\n"
				 "First order only — no O(Omega^2) claim (INV-08).\n\n";

	std::vector<Row> rows;
	double worst_profile_s = 0.0, worst_profile_sp = 0.0, worst_vol = 0.0, worst_ext = 0.0,
		   worst_ext_cross = 0.0;
	bool ref_admissible = true;
	double ref_floor_reported = 0.0, ref_ratio_reported = 0.0;

	std::cout << "B full normalized-profile comparison, star by star\n";
	for (double Mt : {1.0, 1.4, 1.6, 2.0})
	{
		auto ns = Build(cold, wrk / ("M" + std::to_string(int(Mt * 10))), Mt);
		if (!ns || !ns->RotationResponse().valid)
		{
			Report("B build M = " + Sci(Mt, 1), false, "solve or response unavailable");
			continue;
		}
		const auto bg = BackgroundFromProfile(*ns);
		if (bg.N() < 4)
		{
			Report("B background M = " + Sci(Mt, 1), false, "profile lacks metric columns");
			continue;
		}
		std::vector<double> sP, spP;
		ProductionShape(*ns, sP, spP);

		// Same centre convention as production: both begin at the first profile node.
		const auto ref0 = hartle_ref::Solve(bg, 5.0e-3, bg.r[0], 1e-12, 1e-14);
		if (!ref0.ok)
		{
			Report("B reference solve M = " + Sci(Mt, 1), false, ref0.message);
			continue;
		}
		const ShapeCmp c = CompareShapes(sP, spP, ref0.s, ref0.s_prime, bg.r);

		Row w;
		w.n = c.n;
		w.M = ns->GetSequence().m;
		w.R = ns->GetSequence().r;
		w.I_prod = ns->RotationResponse().I;
		w.I_ref = ref0.I_surface;
		w.I_vol = VolumeIntegralI(bg, sP);
		w.max_rel_s = c.max_rel_s;
		w.rms_s = c.rms_s;
		w.r_worst_s = c.r_worst_s;
		w.sp_abs_peak = c.rel_to_scale_sp;
		w.sp_max_rel = c.max_rel_sp;
		w.rms_sp = c.rms_sp;
		w.r_worst_sp = c.r_worst_sp;
		w.drag_c = 1.0 - sP.front();
		w.drag_R = 1.0 - sP.back();
		w.Omega_K = std::sqrt(w.M * Zaki::Physics::SUN_M_KM / (w.R * w.R * w.R));

		worst_profile_s = std::max(worst_profile_s, c.max_rel_s);
		worst_profile_sp = std::max({worst_profile_sp, c.rel_to_scale_sp, c.max_rel_sp});
		worst_vol = std::max(worst_vol, Rel(w.I_vol, w.I_prod));

		// Exterior identities on this star.
		const double s_R = sP.back(), sp_R = spP.back();
		worst_ext = std::max({worst_ext, Rel(s_R, 1.0 - 2.0 * w.I_prod / std::pow(w.R, 3)),
							  Rel(sp_R, 6.0 * w.I_prod / std::pow(w.R, 4))});
		worst_ext_cross =
			std::max({worst_ext_cross, Rel(s_R, 1.0 - 2.0 * w.I_ref / std::pow(w.R, 3)),
					  Rel(sp_R, 6.0 * w.I_ref / std::pow(w.R, 4))});

		// Reference admissibility, demonstrated on the 1.6 Msun star as the task requires.
		if (std::fabs(Mt - 1.6) < 1e-9)
		{
			double floor_s = 0.0, floor_sp = 0.0;
			struct V
			{
				const char *tag;
				double seed, rtol, atol;
			};
			const std::vector<V> variants = {{"tol 1e-13 / 1e-15", 5.0e-3, 1e-13, 1e-15},
											 {"tol 1e-14 / 1e-16", 5.0e-3, 1e-14, 1e-16},
											 {"seed x 1e3", 5.0, 1e-12, 1e-14},
											 {"seed x 1e-3", 5.0e-6, 1e-12, 1e-14}};
			std::cout << "   reference numerical floor on the 1.6 Msun star:\n";
			for (const auto &v : variants)
			{
				const auto rr = hartle_ref::Solve(bg, v.seed, bg.r[0], v.rtol, v.atol);
				if (!rr.ok)
					continue;
				const ShapeCmp cc = CompareShapes(rr.s, rr.s_prime, ref0.s, ref0.s_prime, bg.r);
				floor_s = std::max(floor_s, cc.max_rel_s);
				floor_sp = std::max(floor_sp, cc.rel_to_scale_sp);
				std::cout << "      " << std::setw(20) << std::left << v.tag << std::right
						  << "  s " << Sci(cc.max_rel_s) << "   s' " << Sci(cc.rel_to_scale_sp)
						  << "\n";
			}
			ref_floor_reported = floor_s;
			ref_ratio_reported = c.max_rel_s > 0 ? floor_s / c.max_rel_s : 0.0;
			ref_admissible = floor_s <= kRefFloorMargin * c.max_rel_s &&
							 floor_sp <= kRefFloorMargin * std::max(c.rel_to_scale_sp, 1e-300);
		}

		rows.push_back(w);
	}

	std::cout << "\n   M[Msun]  R[km]     N     max rel s   rms s      worst r    max rel s'  "
				 "rms s'\n";
	for (const auto &w : rows)
		std::cout << "   " << std::fixed << std::setprecision(4) << w.M << "  " << w.R << "  "
				  << std::setw(5) << w.n << "  " << std::scientific << std::setprecision(3)
				  << w.max_rel_s << "  " << w.rms_s << "  " << w.r_worst_s << "  "
				  << w.sp_max_rel << "  " << w.rms_sp << "\n"
				  << std::scientific << std::setprecision(8);

	Report("B0a the independent reference is numerically admissible as the oracle",
		   ref_admissible,
		   "1.6 Msun floor " + Sci(ref_floor_reported) + " vs measured difference (ratio " +
			   Sci(ref_ratio_reported, 1) + ", bound 1e-3)");
	Report("Ba s(r) agrees with the independent profile at every node of every star",
		   rows.size() == 4 && worst_profile_s <= kBoundProfileCMF,
		   "worst relative " + Sci(worst_profile_s) + " over " + std::to_string(rows.size()) +
			   " stars (predeclared bound " + Sci(kBoundProfileCMF, 0) + ")");
	Report("Bb s'(r) agrees with the independent profile at every node of every star",
		   worst_profile_sp <= kBoundProfileCMF,
		   "worst " + Sci(worst_profile_sp) + " (bound " + Sci(kBoundProfileCMF, 0) + ")");

	std::cout << "\n   I comparison (Phase 2B measured 1.32e-5 … 2.15e-5 for these stars)\n";
	std::cout << "   M[Msun]   I_prod[km^3]    I_ref[km^3]     rel        I_vol[km^3]     "
				 "vol rel\n";
	for (const auto &w : rows)
		std::cout << "   " << std::fixed << std::setprecision(4) << w.M << "   "
				  << std::scientific << std::setprecision(8) << w.I_prod << "  " << w.I_ref
				  << "  " << Sci(Rel(w.I_prod, w.I_ref)) << "  " << w.I_vol << "  "
				  << Sci(Rel(w.I_vol, w.I_prod)) << "\n";

	Report("Ca exterior identities s(R) = 1 - 2I/R^3 and s'(R) = 6I/R^4 hold to roundoff "
		   "(internal consistency between the published profile and the published I)",
		   worst_ext <= kBoundExteriorIdentity,
		   "worst relative " + Sci(worst_ext) + " (bound " + Sci(kBoundExteriorIdentity, 0) +
			   ")");
	Report("Cb the same identities hold against the INDEPENDENT solver's I",
		   worst_ext_cross <= kBoundProfileCMF,
		   "worst relative " + Sci(worst_ext_cross) + " (bound " + Sci(kBoundProfileCMF, 0) +
			   ")");
	Report("Da I recovered by volume quadrature over production's own s(r) matches its I",
		   worst_vol <= kBoundVolumeCMF,
		   "worst relative " + Sci(worst_vol) + " (bound " + Sci(kBoundVolumeCMF, 0) + ")");

	// =======================================================================
	//  Spin reversal on a real star
	// =======================================================================
	std::cout << "\nF spin reversal on the 1.4 Msun star\n";
	{
		auto ns = Build(cold, wrk / "rev", 1.4);
		if (ns && ns->RotationResponse().valid)
		{
			const double w = 1.0e2;
			const auto pp = ns->RotationAt(AngularVelocity::FromRadPerSecond(w));
			const auto pm = ns->RotationAt(AngularVelocity::FromRadPerSecond(-w));
			const int n = static_cast<int>(pp.omega_bar.Size());
			bool exact = (pm.Omega_geom == -pp.Omega_geom) && (pm.J == -pp.J) && (pm.I == pp.I);
			for (int i = 0; i < n && exact; ++i)
				exact = (pm.omega_bar[i] == -pp.omega_bar[i]) &&
						(pm.domega_bar[i] == -pp.domega_bar[i]);
			Report("Fa reversing the spin negates Omega, J and both profiles EXACTLY", exact,
				   "bit-exact antisymmetry at all " + std::to_string(n) + " nodes; I unchanged");
		}
		else
			Report("Fa spin reversal", false, "1.4 Msun star unavailable");
	}

	// =======================================================================
	//  Diagnostics
	// =======================================================================
	std::cout << "\nG diagnostics — frame dragging and slow-rotation context\n";
	std::cout << "   M[Msun]  omega(0)/Omega  omega(R)/Omega   Omega_K[1/km]   "
				 "Omega/Omega_K @ 2pi*716Hz\n";
	const double w716_geom = 2.0 * M_PI * 716.0 / 299792.458;
	for (const auto &w : rows)
		std::cout << "   " << std::fixed << std::setprecision(4) << w.M << "   "
				  << std::scientific << std::setprecision(6) << w.drag_c << "    " << w.drag_R
				  << "    " << w.Omega_K << "    " << (w716_geom / w.Omega_K) << "\n";
	{
		bool ordered = !rows.empty();
		for (const auto &w : rows)
			ordered = ordered && w.drag_c > 0.0 && w.drag_c < 1.0 && w.drag_R > 0.0 &&
					  w.drag_R < w.drag_c;
		Report("Ga the frame-dragging fraction lies in 0 < omega/Omega < 1 and decreases "
			   "outward (DIAGNOSTIC ordering, not a derived inequality)",
			   ordered, "measured on all four validated stars");
	}
	std::cout << "\n   NOTE, recorded deliberately: normalization correctness is NOT "
				 "slow-rotation\n"
				 "   truncation accuracy. Omega/Omega_K reaching ~0.6 at 716 Hz says nothing "
				 "here about\n"
				 "   the size of the neglected O(Omega^2) terms; that question belongs with "
				 "the\n"
				 "   O(Omega^2) validation and regime work (ADR-0006 SS10, INV-08).\n";

	fs::remove_all(wrk);
	std::cout << "\n"
			  << (g_fail == 0 ? "CMF first-order physics checks passed"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
