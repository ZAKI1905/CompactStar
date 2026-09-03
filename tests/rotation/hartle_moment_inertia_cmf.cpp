// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_moment_inertia_cmf.cpp
 * @brief Phase 2B-4B — the scale-free first-order Hartle moment of inertia on the
 *        authenticated DS(CMF)-1_with_crust neutron-star sequence.
 *
 *   usage: hartle_moment_inertia_cmf <EOS_DATA_ROOT> [--emit <file.tsv>]
 *
 * PRIMARY REFERENCE is the independent test-side solver in hartle_reference.hpp, not any
 * published fit. Breu & Rezzolla (2016) and Lattimer & Schutz (2005) appear ONLY as
 * secondary plausibility checks with their own EOS scatter, never as numerical truth.
 *
 * The official CompOSE DS(CMF)-1 distribution ships `eos.mr` with exactly two columns,
 * R [km] and M [Msun] — it contains NO moment of inertia. No EOS-specific published I
 * sequence for this model was located; see docs/validation/HARTLE_MOMENT_INERTIA.md §18.
 *
 * SCOPE: I = J/Omega only. No claim about absolute omegabar, Omega, J, physical spin, or
 * O(Omega^2). INV-07 remains UNRESOLVED.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "hartle_reference.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

// ---------------------------------------------------------------------------
//  Geometric <-> physical conversion.
//
//  A geometric moment of inertia has dimensions of length^3; I_phys = I_geom * c^2/G.
//    c^2/G = (2.99792458e10 cm/s)^2 / 6.67430e-8 cgs = 1.346590922e28 g/cm
//    I[g cm^2] = I[km^3] * 1e15 (cm^3/km^3) * 1.346590922e28 (g/cm)
//
//  NOTE, recorded rather than silently absorbed: the repository's own
//  Zaki::Physics::SUN_M_KM = 1.476625038050 km corresponds to M_sun = 1.98835e33 g under
//  this c^2/G, while GSL_CONST_CGSM_SOLAR_MASS = 1.98892e33 g. That is a pre-existing
//  ~2.8e-4 inconsistency between two constants in the build. It affects ONLY the cgs
//  number printed for reference; every validation comparison below is performed in
//  consistent geometric units (km) and is untouched by it.
// ---------------------------------------------------------------------------
static const double kC2_over_G = 1.346590922e28;	   // g/cm
static const double kKm3_to_gcm2 = 1.0e15 * kC2_over_G; // 1.346590922e43

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

static hartle_ref::Background FromProfile(const NStar &ns)
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
	const std::size_t n = static_cast<std::size_t>(r->Size());
	b.r.reserve(n);
	for (std::size_t i = 0; i < n; ++i)
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

static std::vector<std::string> g_labels;

/// Build a star at a chosen radial resolution through the production construction.
/// Identical pattern to Phase 2B-4A, which proved it bit-exact against
/// NStar::SolveTOV_Profile at radial_res = 10000. No production API is added.
static std::unique_ptr<NStar> Build(const fs::path &cold, const fs::path &wrk,
									std::size_t res, bool fixed_ec, double ec, double M)
{
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	if (res)
		tov.SetRadialRes(res);
	std::vector<TOVPoint> pts;
	std::vector<std::string> labels;
	int n = 0;
	if (fixed_ec)
	{
		n = tov.SingleStarSolveToTOVPoints(ec, pts);
		labels = g_labels;
	}
	else
		n = tov.SolveToProfile(M, pts, &labels);
	if (n <= 0 || pts.size() < 4)
		return nullptr;
	return labels.empty() ? std::make_unique<NStar>(pts)
						  : std::make_unique<NStar>(pts, labels);
}

struct Row
{
	double M_target = 0, M = 0, R = 0, C = 0, ec = 0;
	int N = 0;
	double I_prod = 0, I_ref = 0, I_vol = 0;
	double d_prod_ref = 0, d_surf_vol = 0;
	double I_over_MR2 = 0, Ibar = 0, I_cgs = 0;
	double Ibar_BR = 0, I_LS_km3 = 0;
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: hartle_moment_inertia_cmf <EOS_DATA_ROOT> [--emit <file>]\n";
		return 2;
	}
	const fs::path root = argv[1];
	std::string emit;
	for (int i = 2; i < argc; ++i)
		if (std::strcmp(argv[i], "--emit") == 0 && i + 1 < argc)
			emit = argv[++i];

	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "authenticated EOS missing: " << cold.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / "compactstar_hartle_cmf";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 2B-4B — scale-free Hartle I on DS(CMF)-1_with_crust\n"
				 "Primary reference: independent conservative-form solver.\n"
				 "Universal relations are SECONDARY sanity checks only.\n"
				 "Validates I = J/Omega only; INV-07 (absolute normalization) untouched.\n\n";

	// Capture the EOS species labels once, from the canonical production call.
	{
		TOVSolver probe;
		probe.SetWrkDir((wrk / "labels").string());
		probe.ImportEOS(cold.string(), true);
		std::vector<TOVPoint> tmp;
		if (probe.SolveToProfile(1.4, tmp, &g_labels) <= 0 ||
			probe.LastSolveStatus() != CompactStar::Core::TOVSolveStatus::SURFACE_REACHED)
			return 4;
	}

	// =======================================================================
	//  B1 — the real neutron-star sequence
	// =======================================================================
	std::cout << "B1 stable-branch sequence: production vs independent reference\n";
	std::vector<Row> rows;
	for (double Mt : {1.0, 1.4, 1.6, 2.0})
	{
		Row w;
		w.M_target = Mt;
		auto ns = Build(cold, wrk / ("M" + std::to_string(int(Mt * 10))), 0, false, 0.0, Mt);
		if (!ns)
		{
			Report("B1 build M=" + std::to_string(Mt), false, "solve failed");
			continue;
		}
		const auto &s = ns->GetSequence();
		w.M = s.m;
		w.R = s.r;
		w.ec = s.ec;
		w.N = static_cast<int>(ns->Profile().GetRadius()->Size());
		w.I_prod = s.I;

		const auto bg = FromProfile(*ns);
		const auto res = hartle_ref::Solve(bg, 1.0, bg.r[1]);
		if (!res.ok)
		{
			Report("B1 reference solve M=" + std::to_string(Mt), false, res.message);
			continue;
		}
		w.I_ref = res.I_surface;
		w.I_vol = res.I_volume;
		w.d_prod_ref = Rel(w.I_prod, w.I_ref);
		w.d_surf_vol = Rel(w.I_vol, w.I_ref);

		const double M_km = w.M * Zaki::Physics::SUN_M_KM;
		w.C = M_km / w.R;
		w.I_over_MR2 = w.I_prod / (M_km * w.R * w.R);
		w.Ibar = w.I_prod / (M_km * M_km * M_km);
		w.I_cgs = w.I_prod * kKm3_to_gcm2;

		// Breu & Rezzolla 2016 (MNRAS 459, 646), eq. (7): Ibar = I/M^3 fitted as a
		// polynomial in 1/C with C = M/R, coefficients a1..a4 below. Quoted as a fit
		// with a few-percent EOS scatter, NOT as an exact relation.
		const double a1 = 8.134e-1, a2 = 2.101e-1, a3 = 3.175e-3, a4 = -2.717e-4;
		w.Ibar_BR = a1 / w.C + a2 / (w.C * w.C) + a3 / std::pow(w.C, 3) +
					a4 / std::pow(w.C, 4);

		// Lattimer & Schutz 2005 (ApJ 629, 979): I ~ 0.237 M R^2 (1 + 4.2 x + 90 x^4),
		// x = (M/Msun)(km/R). Quoted to ~10%.
		const double x = w.M / w.R;
		w.I_LS_km3 = 0.237 * M_km * w.R * w.R *
					 (1.0 + 4.2 * x + 90.0 * std::pow(x, 4));
		rows.push_back(w);
	}

	std::cout << "   M[Msun]   R[km]      C        I_prod[km^3]  I_ref[km^3]   "
				 "prod/ref     surf/vol\n";
	for (const auto &w : rows)
		std::cout << "   " << std::setprecision(6) << std::fixed << w.M << "  " << w.R
				  << "  " << w.C << "   " << std::scientific << std::setprecision(6)
				  << w.I_prod << "  " << w.I_ref << "  " << Sci(w.d_prod_ref) << "  "
				  << Sci(w.d_surf_vol) << "\n"
				  << std::scientific << std::setprecision(8);

	{
		double worst_pr = 0, worst_sv = 0;
		for (const auto &w : rows)
		{
			worst_pr = std::max(worst_pr, w.d_prod_ref);
			worst_sv = std::max(worst_sv, w.d_surf_vol);
		}
		Report("B1a production I agrees with the independent reference across the sequence",
			   !rows.empty() && worst_pr < 1e-3, "worst relative difference " + Sci(worst_pr));
		Report("B1b surface and volume extractions of I agree",
			   !rows.empty() && worst_sv < 1e-4, "worst relative difference " + Sci(worst_sv));
		Report("B1c every model is sub-Buchdahl and physically ordered",
			   !rows.empty() && std::all_of(rows.begin(), rows.end(), [](const Row &w) {
				   return 2.0 * w.C < 8.0 / 9.0 && w.I_prod > 0.0 && w.R > 0.0;
			   }),
			   "compactness range " + Sci(rows.front().C) + " … " + Sci(rows.back().C));
	}

	// =======================================================================
	//  B2 — physical magnitude and published universal relations (SANITY ONLY)
	// =======================================================================
	std::cout << "\nB2 magnitude and universal relations (secondary sanity checks)\n";
	std::cout << "   M[Msun]  I[g cm^2]      I/(M R^2)   Ibar=I/M^3   Breu-Rezzolla  "
				 "ratio    Lattimer-Schutz ratio\n";
	for (const auto &w : rows)
		std::cout << "   " << std::setprecision(3) << std::fixed << w.M << "   "
				  << std::scientific << std::setprecision(4) << w.I_cgs << "   "
				  << w.I_over_MR2 << "   " << w.Ibar << "   " << w.Ibar_BR << "   "
				  << std::setprecision(4) << std::fixed << (w.Ibar / w.Ibar_BR) << "        "
				  << (w.I_prod / w.I_LS_km3) << "\n"
				  << std::scientific << std::setprecision(8);
	{
		bool mag_ok = true, br_ok = true, ls_ok = true;
		for (const auto &w : rows)
		{
			mag_ok = mag_ok && w.I_cgs > 5e44 && w.I_cgs < 5e45;
			br_ok = br_ok && Rel(w.Ibar, w.Ibar_BR) < 0.15;
			ls_ok = ls_ok && Rel(w.I_prod, w.I_LS_km3) < 0.20;
		}
		Report("B2a I is of the expected neutron-star magnitude (0.5–5 x 1e45 g cm^2)",
			   mag_ok, "range " + Sci(rows.front().I_cgs) + " … " + Sci(rows.back().I_cgs));
		Report("B2b Ibar sits within the Breu-Rezzolla fit's EOS scatter (<15%)", br_ok,
			   "a fit, not truth — deviation is EOS-dependent by construction");
		Report("B2c I sits within the Lattimer-Schutz relation's ~10-20% band", ls_ok,
			   "an empirical relation, not truth");
	}

	// =======================================================================
	//  B3 — radial convergence of I
	// =======================================================================
	std::cout << "\nB3 radial convergence of I (fixed ec*, the 2B-4A resolution matrix)\n";
	const double ec_star = 7.312533426775e14; // canonical 1.6 Msun central density (2B-4A)
	std::cout << "      res     N      dr_eff[km]   I_prod[km^3]    I_ref[km^3]     "
				 "prod/ref\n";
	std::vector<double> Ip_seq, Ir_seq, dr_seq, R_seq;
	for (std::size_t res : {2500, 5000, 10000, 20000, 40000})
	{
		auto ns = Build(cold, wrk / ("C" + std::to_string(res)), res, true, ec_star, 0.0);
		if (!ns) { std::cout << "      res=" << res << " FAILED\n"; continue; }
		const auto bg = FromProfile(*ns);
		const auto rr = hartle_ref::Solve(bg, 1.0, bg.r[1]);
		const double Ip = ns->GetSequence().I;
		const double dr = bg.R() / double(bg.N() - 1);
		Ip_seq.push_back(Ip);
		Ir_seq.push_back(rr.ok ? rr.I_surface : 0.0);
		dr_seq.push_back(dr);
		R_seq.push_back(bg.R());
		std::cout << "   " << std::setw(7) << res << " " << std::setw(6) << bg.N() << "  "
				  << Sci(dr) << "  " << Sci(Ip, 8) << "  " << Sci(rr.I_surface, 8) << "  "
				  << Sci(rr.ok ? Rel(Ip, rr.I_surface) : -1.0) << "\n";
	}
	{
		// Differences and observed order on the three finest, with the ACTUAL dr ratio.
		std::string diffs;
		for (std::size_t i = 1; i + 1 < Ip_seq.size(); ++i)
			diffs += " |d| " + Sci(std::fabs(Ip_seq[i] - Ip_seq[i + 1]));
		double p = std::nan("");
		const std::size_t n = Ip_seq.size();
		if (n >= 3)
		{
			const double d12 = Ip_seq[n - 3] - Ip_seq[n - 2];
			const double d23 = Ip_seq[n - 2] - Ip_seq[n - 1];
			const double scale = std::fabs(Ip_seq[n - 1]);
			if (d12 * d23 > 0.0 && std::fabs(d23) / scale > 1e-9)
			{
				const double r = 0.5 * (dr_seq[n - 3] / dr_seq[n - 2] +
										dr_seq[n - 2] / dr_seq[n - 1]);
				p = std::log(std::fabs(d12 / d23)) / std::log(r);
			}
		}
		std::cout << "      successive production |dI|:" << diffs << "\n";
		std::cout << "      observed order (3 finest): "
				  << (std::isfinite(p) ? Sci(p) : std::string("not reliably measurable"))
				  << "\n";
		// Phase 2B-4A established that the stellar RADIUS itself does not converge
		// smoothly: it is fixed by a step-dependent surface-termination event at the EOS
		// table floor, with a sub-linear, drifting order. I is an r^4-weighted integral
		// over the same star, so §15's question is whether that weighting AMPLIFIES the
		// inherited surface behaviour. The criterion is therefore derived from the
		// measured radius behaviour on this very matrix, not invented: I's relative
		// spread must be no larger than R's over the same resolutions.
		auto spread = [](const std::vector<double> &v, std::size_t from) {
			if (v.size() <= from + 1) return 1.0;
			double lo = v[from], hi = v[from];
			for (std::size_t i = from; i < v.size(); ++i)
			{
				lo = std::min(lo, v[i]);
				hi = std::max(hi, v[i]);
			}
			return lo > 0.0 ? (hi - lo) / lo : 1.0;
		};
		// index 1 onwards = 5000..40000; radial_res = 2500 is outside the asymptotic
		// regime (2B-4A) and is reported but excluded, never silently dropped.
		const double sI = spread(Ip_seq, 1);
		const double sR = spread(R_seq, 1);
		std::cout << "      relative spread over res 5000-40000:  I = " << Sci(sI)
				  << "   R = " << Sci(sR) << "\n";
		std::cout << "      (radial_res = 2500 is outside the asymptotic regime per 2B-4A "
					 "and is excluded from the spread, not dropped: I = "
				  << Sci(Ip_seq.empty() ? 0.0 : Ip_seq[0], 8) << " km^3)\n";
		Report("B3a the r^4 weighting does NOT amplify the inherited surface behaviour: "
			   "I is no more grid-sensitive than R itself",
			   Ip_seq.size() == 5 && sI <= sR,
			   "I spread " + Sci(sI) + " vs R spread " + Sci(sR));

		double worst = 0, best = 1.0;
		for (std::size_t i = 0; i < Ip_seq.size(); ++i)
			if (Ir_seq[i] > 0)
			{
				const double d = Rel(Ip_seq[i], Ir_seq[i]);
				worst = std::max(worst, d);
				best = std::min(best, d);
			}
		Report("B3b production tracks the reference at every resolution", worst < 1e-3,
			   "worst production/reference difference " + Sci(worst));
		Report("B3c the production/reference agreement is resolution-INDEPENDENT, so the "
			   "residual grid sensitivity belongs to the TOV background, not to Hartle",
			   worst > 0.0 && worst / best < 3.0,
			   "prod/ref difference stays in [" + Sci(best) + ", " + Sci(worst) + "]");
	}

	// =======================================================================
	//  B4 — reference numerical floor on the real star
	// =======================================================================
	std::cout << "\nB4 reference numerical floor vs the production/reference gap\n";
	{
		auto ns = Build(cold, wrk / "floor", 0, false, 0.0, 1.4);
		if (ns)
		{
			const auto bg = FromProfile(*ns);
			const auto loose = hartle_ref::Solve(bg, 1.0, bg.r[1], 1e-10, 1e-12);
			const auto tight = hartle_ref::Solve(bg, 1.0, 1e-6, 1e-14, 1e-16);
			const double floor = Rel(loose.I_surface, tight.I_surface);
			const double gap = Rel(ns->GetSequence().I, tight.I_surface);
			std::cout << "      reference floor = " << Sci(floor)
					  << "   production/reference gap = " << Sci(gap) << "\n";
			Report("B4a reference numerical error is subdominant to the production gap",
				   floor < 0.1 * gap,
				   "floor " + Sci(floor) + " vs gap " + Sci(gap) + " (ratio " +
					   Sci(floor > 0 ? gap / floor : 0.0, 1) + ")");
		}
	}

	if (!emit.empty())
	{
		std::ofstream o(emit);
		o << "# Phase 2B-4B — scale-free first-order Hartle I on DS(CMF)-1_with_crust\n"
		  << "# I in geometric km^3; I_cgs = I_km3 * 1.346590922e43 g cm^2\n"
		  << "# Ibar = I/M^3 and C = M/R in geometric units (SUN_M_KM = 1.476625038050 km)\n"
		  << "# Breu-Rezzolla and Lattimer-Schutz are SANITY FITS, not reference truth.\n"
		  << "M_target\tM_achieved\tR_km\tcompactness\tepsilon_c\tradial_res\tN_profile"
			 "\tI_production_km3\tI_reference_km3\tI_volume_km3\tprod_ref_rel\tsurf_vol_rel"
			 "\tI_cgs\tI_over_MR2\tIbar\tBreu_Rezzolla\tLattimer_Schutz_km3\n";
		o << std::setprecision(12) << std::scientific;
		for (const auto &w : rows)
			o << w.M_target << "\t" << w.M << "\t" << w.R << "\t" << w.C << "\t" << w.ec
			  << "\t" << 10000 << "\t" << w.N << "\t" << w.I_prod << "\t" << w.I_ref << "\t"
			  << w.I_vol << "\t" << w.d_prod_ref << "\t" << w.d_surf_vol << "\t" << w.I_cgs
			  << "\t" << w.I_over_MR2 << "\t" << w.Ibar << "\t" << w.Ibar_BR << "\t"
			  << w.I_LS_km3 << "\n";
		std::cout << "\n  artifact written: " << emit << "\n";
	}

	fs::remove_all(wrk);
	std::cout << "\n"
			  << (g_fail == 0 ? "CMF Hartle-I checks passed"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
