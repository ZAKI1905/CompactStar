// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file tov_reference_cmf.cpp
 * @brief Validate the PRODUCTION TOV path against the OFFICIAL CompOSE mass-radius
 *        table `eos.mr` shipped with the DS(CMF)-1_with_crust distribution.
 *
 *   usage: tov_reference_cmf <EOS_DATA_ROOT> [--emit <baseline.tsv>]
 *
 * The reference is an EXTERNAL, third-party result: `eos.mr` was produced by the CompOSE
 * project's own TOV solver from `eos.nb` / `eos.thermo`. No CompactStar-generated output is
 * used as a reference anywhere in this file.
 *
 * KNOWN REPRESENTATION GAP (measured, documented, NOT repaired here — the production
 * scientific source is frozen for this phase). The reference `eos.mr` was computed from the
 * CompOSE tables `eos.nb`/`eos.thermo`, whereas production TOV integrates the tabulated
 * `DS(CMF)-1_with_crust.eos`. Those two representations of the same EOS agree to ~5e-9 in
 * n_B and p at every grid point, and in energy density above the crust/core splice, but
 * BELOW n_B = 4.0e-2 fm^-3 the `.eos` energy density is uniformly larger by
 *
 *      delta = 6.8866e-04  =  m_n(free, 939.5654 MeV) / m_n(CMF header, 938.9187125 MeV) - 1
 *
 * The origin is a deliberate branch in CompactStar's own converter,
 * CompactStar/EOS/src/CompOSE_EOS.cpp:294-311, which normalizes the crust to a hardcoded
 * 939.5653 MeV and the core to the header mass m_n. Measured and found NOT to explain the
 * radius residual (see B7 below and docs/validation/TOV_REFERENCE.md §3.2, §5).
 *
 * TOLERANCE BUDGET (fixed a priori from the numerics below, NOT tuned to observed results):
 *   - production radial grid quantization of R: step = r_max/radial_res = 70 km/10000 = 7 m,
 *     so dR/R <~ 5.2e-4 for a ~13.5 km star;
 *   - crust energy-density convention offset above: ~6.9e-4 fractional;
 *   - linear interpolation of a 179-point reference curve in M;
 *   - the target-mass bisection tolerance inside SolveToProfile.
 * Summed with margin: 0.5 % on R, 1 % on M_max. Tier A already certifies the differential
 * equation itself to machine precision; Tier B's job is to catch gross defects in EOS
 * ingestion, unit handling, the center/surface treatment and the central-density search.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;
using namespace CompactStar::Core;

static const double kTolR = 5.0e-3;	 // 0.5 % on radius
static const double kTolMmax = 1.0e-2; // 1 % on maximum mass

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}
static double Rel(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

// --------------------------------------------------------------------------
// Official CompOSE mass-radius reference
// --------------------------------------------------------------------------
struct MR
{
	std::vector<double> R, M; // km, Msun — as tabulated (ordered by decreasing R)
	size_t i_max = 0;

	bool Load(const fs::path &f)
	{
		std::ifstream in(f);
		if (!in)
			return false;
		double r, m;
		while (in >> r >> m)
		{
			R.push_back(r);
			M.push_back(m);
		}
		if (M.size() < 10)
			return false;
		i_max = size_t(std::max_element(M.begin(), M.end()) - M.begin());
		return true;
	}
	double Mmax() const { return M[i_max]; }
	double R_at_Mmax() const { return R[i_max]; }

	/// R(M) on the STABLE branch (everything up to and including the maximum-mass point),
	/// by linear interpolation in M. Returns <0 if M is outside the stable branch.
	double R_of_M(double m) const
	{
		for (size_t i = 0; i + 1 <= i_max; ++i)
		{
			const double a = M[i], b = M[i + 1];
			if ((m >= a && m <= b) || (m <= a && m >= b))
			{
				if (a == b)
					return R[i];
				const double t = (m - a) / (b - a);
				return R[i] + t * (R[i + 1] - R[i]);
			}
		}
		return -1.0;
	}
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: tov_reference_cmf <EOS_DATA_ROOT> [--emit <baseline.tsv>]\n";
		return 2;
	}
	const fs::path root = argv[1];
	std::string emit;
	for (int i = 2; i < argc; ++i)
		if (std::string(argv[i]) == "--emit" && i + 1 < argc)
			emit = argv[++i];

	const fs::path dist = root / "DS-CMF-1-with-crust";
	const fs::path cold = dist / "DS(CMF)-1_with_crust.eos";
	const fs::path mrf = dist / "eos.mr";

	std::cout << std::fixed << std::setprecision(6);
	std::cout << "Production TOV path vs official CompOSE eos.mr (DS(CMF)-1_with_crust)\n\n";

	if (!fs::exists(cold) || !fs::exists(mrf))
	{
		std::cout << "  [SKIP] authenticated EOS distribution incomplete under "
				  << dist.string() << " (need the .eos table and the official eos.mr)\n";
		return 77; // CTest SKIP_RETURN_CODE
	}

	MR ref;
	if (!ref.Load(mrf))
	{
		Report("R0 official eos.mr parsed", false, "could not read " + mrf.string());
		return 1;
	}
	std::cout << "  reference eos.mr: " << ref.M.size() << " rows"
			  << "   M_max = " << ref.Mmax() << " Msun @ R = " << ref.R_at_Mmax() << " km"
			  << "\n                    R(1.4) = " << ref.R_of_M(1.4)
			  << "   R(1.6) = " << ref.R_of_M(1.6) << " km\n\n";
	Report("R0 official eos.mr is a usable stable-branch reference",
		   ref.M.size() > 100 && ref.Mmax() > 1.5 && ref.R_of_M(1.4) > 0.0,
		   std::to_string(ref.M.size()) + " rows, M_max=" + std::to_string(ref.Mmax()));

	std::vector<std::string> emit_rows;

	// ----------------------------------------------------------------------
	// B1 — target-mass path (the one production actually uses)
	// ----------------------------------------------------------------------
	std::cout << "\nB1 target-mass path: NStar::SolveTOV_Profile\n";
	struct B1Row { double M, R_prod, R_ref; };
	std::vector<B1Row> b1;
	double R16_production = 0.0;
	for (double Mt : {1.0, 1.4, 1.6})
	{
		NStar ns;
		ns.SetWrkDir(fs::temp_directory_path().string());
		const int n = ns.SolveTOV_Profile(cold.string(), Mt, "tov_ref");
		if (n <= 0)
		{
			Report("B1 solve M=" + std::to_string(Mt), false, "SolveTOV_Profile failed");
			continue;
		}
		const auto &s = ns.GetSequence();
		const double Rref = ref.R_of_M(s.m); // compare at the ACHIEVED mass
		const double dR = Rel(s.r, Rref);
		if (std::fabs(Mt - 1.6) < 1e-9)
			R16_production = s.r;
		if (Rref > 0.0)
			b1.push_back({s.m, s.r, Rref});
		std::cout << "    target " << Mt << " -> M=" << s.m << " R=" << s.r
				  << " km   ref R(M)=" << Rref << " km   dR=" << std::scientific << dR
				  << std::fixed << "   ec=" << std::scientific << s.ec << std::fixed << "\n";
		emit_rows.push_back(std::to_string(Mt) + "\t" + std::to_string(s.m) + "\t" +
							std::to_string(s.r) + "\t" + std::to_string(Rref));
		Report("B1 R matches official eos.mr at M=" + std::to_string(Mt),
			   Rref > 0.0 && dR < kTolR, "rel " + std::to_string(dR) +
										 " (budget " + std::to_string(kTolR) + ")");
		Report("B1 target mass actually achieved at M=" + std::to_string(Mt),
			   Rel(s.m, Mt) < 1.0e-3, "achieved " + std::to_string(s.m));
	}

	// ----------------------------------------------------------------------
	// B2 — fixed-central-density path and the maximum mass
	// ----------------------------------------------------------------------
	std::cout << "\nB2 fixed-central-density path: TOVSolver::SingleStarSolveToTOVPoints\n";
	TOVSolver tov;
	tov.SetWrkDir(fs::temp_directory_path().string());
	tov.ImportEOS(cold.string(), true);

	double Mmax = 0.0, R_at_Mmax = 0.0, ec_at_Mmax = 0.0;
	std::vector<std::pair<double, double>> curve; // (M, R)
	for (int i = 0; i <= 60; ++i)
	{
		const double ec = std::pow(10.0, std::log10(3.0e14) +
										 (std::log10(5.0e15) - std::log10(3.0e14)) * i / 60.0);
		std::vector<TOVPoint> pts;
		if (tov.SingleStarSolveToTOVPoints(ec, pts) <= 0 || pts.empty())
			continue;
		curve.emplace_back(pts.back().m, pts.back().r);
		if (pts.back().m > Mmax)
		{
			Mmax = pts.back().m;
			R_at_Mmax = pts.back().r;
			ec_at_Mmax = ec;
		}
	}
	std::cout << "    swept " << curve.size() << " central densities;  M_max = " << Mmax
			  << " Msun @ R = " << R_at_Mmax << " km  (ec = " << std::scientific
			  << ec_at_Mmax << std::fixed << " g/cm^3)\n";
	Report("B2 maximum mass matches the official sequence", Rel(Mmax, ref.Mmax()) < kTolMmax,
		   "production " + std::to_string(Mmax) + " vs official " +
			   std::to_string(ref.Mmax()) + ", rel " + std::to_string(Rel(Mmax, ref.Mmax())));
	Report("B2 radius at maximum mass matches", Rel(R_at_Mmax, ref.R_at_Mmax()) < 2.0 * kTolR,
		   "production " + std::to_string(R_at_Mmax) + " vs official " +
			   std::to_string(ref.R_at_Mmax()) + ", rel " +
			   std::to_string(Rel(R_at_Mmax, ref.R_at_Mmax())));

	// The two production entry points must agree with each other.
	{
		double Rc = -1.0;
		for (size_t i = 0; i + 1 < curve.size(); ++i)
		{
			const double a = curve[i].first, b = curve[i + 1].first;
			if (a <= 1.4 && 1.4 <= b)
			{
				const double t = (1.4 - a) / (b - a);
				Rc = curve[i].second + t * (curve[i + 1].second - curve[i].second);
				break;
			}
		}
		Report("B2 fixed-ec and target-mass entry points agree at 1.4 Msun",
			   Rc > 0.0 && Rel(Rc, ref.R_of_M(1.4)) < kTolR,
			   "ec-path R(1.4)=" + std::to_string(Rc) + " vs official " +
				   std::to_string(ref.R_of_M(1.4)));
	}

	// ----------------------------------------------------------------------
	// B3 — radial-resolution sensitivity of the production discretization
	// ----------------------------------------------------------------------
	std::cout << "\nB3 radial-resolution sensitivity at fixed ec = 9e14 g/cm^3\n";
	std::cout << "    (production default is r_max = 70 km, radial_res = 10000 -> 7 m step)\n";
	double R_coarse = 0.0, R_fine = 0.0;
	for (auto cfgp : std::vector<std::pair<double, size_t>>{
			 {70.0, 10000}, {70.0, 40000}, {20.0, 10000}, {20.0, 80000}})
	{
		TOVSolver t2;
		t2.SetWrkDir(fs::temp_directory_path().string());
		t2.ImportEOS(cold.string(), true);
		t2.SetMaxRadius(cfgp.first);
		t2.SetRadialRes(cfgp.second);
		std::vector<TOVPoint> pts;
		if (t2.SingleStarSolveToTOVPoints(9.0e14, pts) <= 0)
			continue;
		const double step_m = cfgp.first * 1000.0 / double(cfgp.second);
		std::cout << "    r_max=" << std::setw(5) << cfgp.first << " km  res=" << std::setw(6)
				  << cfgp.second << "  (step " << std::setprecision(2) << step_m
				  << " m)  ->  M=" << std::setprecision(8) << pts.back().m
				  << "  R=" << pts.back().r << " km\n"
				  << std::setprecision(6);
		if (cfgp.first == 70.0 && cfgp.second == 10000)
			R_coarse = pts.back().r;
		if (cfgp.first == 20.0 && cfgp.second == 80000)
			R_fine = pts.back().r;
	}
	Report("B3 production default resolution is within its own quantization budget",
		   R_coarse > 0.0 && R_fine > 0.0 && Rel(R_coarse, R_fine) < kTolR,
		   "default R=" + std::to_string(R_coarse) + " vs 80x-finer R=" +
			   std::to_string(R_fine) + ", rel " + std::to_string(Rel(R_coarse, R_fine)));

	// ----------------------------------------------------------------------
	// B4/B5 — center and surface condition audits
	// ----------------------------------------------------------------------
	std::cout << "\nB4/B5 center and surface conditions\n";
	{
		std::vector<TOVPoint> pts;
		tov.SingleStarSolveToTOVPoints(9.0e14, pts);
		const auto &c = pts.front();
		const auto &s = pts.back();
		const double r_cm = c.r * 1.0e5;
		const double m_g = c.m * 1.98892e33;
		// Regular series at the center: m(r) = (4/3) pi r^3 eps_c.
		const double m_expect = (4.0 / 3.0) * M_PI * r_cm * r_cm * r_cm * c.e;
		std::cout << "    first point:  r=" << std::scientific << r_cm << " cm  m="
				  << m_g << " g  eps=" << c.e << " g/cm^3  p=" << c.p << "\n";
		std::cout << "    last  point:  r=" << std::fixed << s.r << " km  p="
				  << std::scientific << s.p << " dyne/cm^2  eps=" << s.e
				  << " g/cm^3\n" << std::fixed;
		Report("B4 center uses the regular series m(r)=(4/3)pi r^3 eps_c",
			   Rel(m_g, m_expect) < 1.0e-6,
			   "m(r_min)=" + std::to_string(m_g) + " vs (4/3)pi r^3 eps_c=" +
				   std::to_string(m_expect));
		Report("B4 central offset r_min is negligible against the stellar radius",
			   r_cm / (s.r * 1.0e5) < 1.0e-5,
			   "r_min/R = " + std::to_string(r_cm / (s.r * 1.0e5)));
		// The DS(CMF)-1 table stops at n_B = 1e-7 fm^-3; there is no vacuum boundary in
		// it. Production therefore terminates the integration at the table's own lowest
		// pressure, which is the same outer boundary the reference eos.mr was built with.
		// Read that floor straight from the .eos file — independently of the solver.
		double p_floor = 0.0, e_floor = 0.0;
		{
			std::ifstream in(cold);
			std::string line;
			std::getline(in, line); // header
			if (std::getline(in, line))
			{
				std::istringstream ls(line);
				ls >> e_floor >> p_floor;
			}
		}
		std::cout << "    EOS table floor: p=" << std::scientific << p_floor
				  << " dyne/cm^2  eps=" << e_floor << " g/cm^3\n" << std::fixed;
		Report("B5 the surface is set by the EOS table floor, not by anything else",
			   p_floor > 0.0 && s.p >= p_floor && s.p / p_floor < 2.0,
			   "p_surf=" + std::to_string(s.p) + " vs table floor " +
				   std::to_string(p_floor) + " (ratio " +
				   std::to_string(p_floor > 0.0 ? s.p / p_floor : -1.0) + ")");
		Report("B5 surface pressure is negligible against the center",
			   s.p > 0.0 && s.p / c.p < 1.0e-8,
			   "p_surf/p_c = " + std::to_string(s.p / c.p));
		Report("B5 mass is essentially converged at the surface",
			   Rel(s.m, pts[pts.size() - 2].m) < 1.0e-6,
			   "last-step mass change rel " +
				   std::to_string(Rel(s.m, pts[pts.size() - 2].m)));
	}

	// ----------------------------------------------------------------------
	// B6 — consistency with the star used by the passive-cooling baseline
	// ----------------------------------------------------------------------
	std::cout << "\nB6 cooling-fingerprint consistency\n";
	{
		// The passive-cooling regression baseline was established on a 1.6 Msun star.
		// Its radius must be reproduced here and must sit inside the same budget
		// against the official reference. This checks that the thermal baseline rests
		// on a TOV solution consistent with the external reference — it does NOT use
		// CompactStar output as the reference.
		const double Rref16 = ref.R_of_M(1.6);
		Report("B6 the 1.6 Msun cooling star agrees with the official reference",
			   R16_production > 0.0 && Rel(R16_production, Rref16) < kTolR,
			   "cooling-path R=" + std::to_string(R16_production) + " km vs official " +
				   std::to_string(Rref16) + " km, rel " +
				   std::to_string(Rel(R16_production, Rref16)));
	}

	// ----------------------------------------------------------------------
	// B-anchor — published sanity anchor, independent of eos.mr
	// ----------------------------------------------------------------------
	// The CMF (Dexheimer & Schramm) family is published with a maximum mass near
	// 2.1 Msun and a 1.4 Msun radius near 13-14 km; see the distribution README:
	//   V. Dexheimer, Publ. Astron. Soc. Australia 34 (2017);
	//   V. Dexheimer, R. O. Gomes, T. Klaehn, S. Han, M. Salinas,
	//   Phys. Rev. C 103 (2021) 025808.
	// These are deliberately WIDE brackets: their job is to catch a gross unit or
	// scaling error that a self-consistent comparison against eos.mr could not,
	// not to certify agreement at the percent level.
	std::cout << "\nB-anchor published sanity brackets (CMF family)\n";
	{
		double R14 = 0.0;
		for (const auto &row : b1)
			if (std::fabs(row.M - 1.4) < 0.01)
				R14 = row.R_prod;
		Report("BA maximum mass is inside the published CMF bracket [1.9, 2.3] Msun",
			   Mmax > 1.9 && Mmax < 2.3, "M_max = " + std::to_string(Mmax));
		Report("BA R(1.4) is inside the published CMF bracket [12, 15] km",
			   R14 > 12.0 && R14 < 15.0, "R(1.4) = " + std::to_string(R14));
		Report("BA the star is causal and sub-Buchdahl at maximum mass",
			   2.0 * 6.67430e-8 * Mmax * 1.98892e33 /
					   (R_at_Mmax * 1.0e5 * 2.99792458e10 * 2.99792458e10) < 8.0 / 9.0,
			   "2GM/Rc^2 = " +
				   std::to_string(2.0 * 6.67430e-8 * Mmax * 1.98892e33 /
								  (R_at_Mmax * 1.0e5 * 2.99792458e10 * 2.99792458e10)));
	}

	// ----------------------------------------------------------------------
	// B7 — attribution of the systematic radius residual
	// ----------------------------------------------------------------------
	// The B1 residuals are systematic, single-signed (production R is always smaller than
	// the reference) and grow toward low mass. That is the signature of a surface
	// boundary, not of a solver error: production terminates at the lowest pressure the
	// DS(CMF)-1 table contains (n_B = 1e-7 fm^-3, eps ~ 1.7e8 g/cm^3), which is still
	// inside the outer crust. The layer below it carries negligible MASS but tens of
	// metres of RADIUS.
	//
	// Estimate that omitted height hydrostatically. Near the surface p << rho c^2, so
	//     dp/dr = -rho g_eff ,   g_eff = G M / ( r^2 sqrt(1 - 2GM/(r c^2)) ) ,
	// and with a local polytrope rho = rho0 (p/p0)^(1/Gamma) fitted to the two lowest
	// tabulated rows,
	//     delta_r = (1/g_eff) * (p0/rho0) / (1 - 1/Gamma) .
	// If the residual is a boundary-convention effect it must be BOUNDED by this height
	// and must be the SAME FRACTION of it at every mass. A solver or unit error would
	// show no such scaling.
	std::cout << "\nB7 attribution of the systematic radius residual\n";
	{
		double e0 = 0, p0 = 0, e1 = 0, p1 = 0;
		{
			std::ifstream in(cold);
			std::string line;
			std::getline(in, line);
			std::getline(in, line);
			{ std::istringstream ls(line); ls >> e0 >> p0; }
			std::getline(in, line);
			{ std::istringstream ls(line); ls >> e1 >> p1; }
		}
		const double Gamma = std::log(p1 / p0) / std::log(e1 / e0);
		const double I_erg_g = (p0 / e0) / (1.0 - 1.0 / Gamma);
		std::cout << "    EOS floor: p0=" << std::scientific << p0 << " dyne/cm^2, rho0="
				  << e0 << " g/cm^3, local Gamma=" << std::fixed << Gamma << "\n";

		const double G = 6.67430e-8, c2 = 2.99792458e10 * 2.99792458e10,
					 Msun = 1.98892e33;
		std::vector<double> ratios;
		std::cout << "      M      R_prod      R_ref    residual   omitted-layer   fraction\n";
		for (const auto &row : b1)
		{
			const double M_g = row.M * Msun, R_cm = row.R_prod * 1.0e5;
			const double g_eff = G * M_g / (R_cm * R_cm *
											std::sqrt(1.0 - 2.0 * G * M_g / (R_cm * c2)));
			const double dr_km = (I_erg_g / g_eff) / 1.0e5;
			const double resid_km = row.R_ref - row.R_prod;
			const double frac = resid_km / dr_km;
			ratios.push_back(frac);
			std::cout << "   " << std::setprecision(3) << row.M << "   " << std::setprecision(6)
					  << row.R_prod << "  " << row.R_ref << "   " << std::setprecision(5)
					  << resid_km << " km    " << dr_km << " km     "
					  << std::setprecision(4) << frac << "\n" << std::setprecision(6);
		}
		bool bounded = !ratios.empty();
		for (double f : ratios)
			bounded = bounded && f > 0.0 && f < 1.0;
		Report("B7 the residual is bounded by the omitted outer-crust layer",
			   bounded, "fractions all in (0,1)");

		double lo = 1e9, hi = -1e9;
		for (double f : ratios) { lo = std::min(lo, f); hi = std::max(hi, f); }
		const double spread = (ratios.empty() || lo <= 0.0) ? 1.0 : (hi - lo) / lo;
		Report("B7 the residual is the same fraction of that layer at every mass "
			   "(a boundary convention, not a solver error)",
			   spread < 0.15,
			   "fraction spans [" + std::to_string(lo) + ", " + std::to_string(hi) +
				   "], relative spread " + std::to_string(spread));
	}

	if (!emit.empty())
	{
		std::ofstream o(emit);
		o << "# CompactStar TOV reference comparison — DS(CMF)-1_with_crust\n"
		  << "# reference: official CompOSE eos.mr (external, third-party)\n"
		  << "# target_M\tachieved_M\tR_production_km\tR_official_km\n";
		for (const auto &r : emit_rows)
			o << r << "\n";
		o << "# M_max_production\t" << Mmax << "\tM_max_official\t" << ref.Mmax() << "\n";
		std::cout << "\n  wrote " << emit << "\n";
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "production TOV path agrees with the official sequence"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
