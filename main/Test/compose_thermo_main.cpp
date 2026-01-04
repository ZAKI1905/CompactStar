//------------------------------------------------------------------------------
// File: tests/compose_thermo_main.cpp
//
// Minimal sanity-test driver for CompactStar::EOS::CompOSE_Thermo.
//
// This program:
//  1) Resolves the CompOSE EOS directory (EOS/CompOSE/CMF-1_general)
//  2) Loads eos.t / eos.nb / eos.yq / eos.thermo
//  3) Prints grid metadata and a handful of Q2, dQ2dT, and CvDensity samples
//
// Expected usage:
//   - Add this file to a small test executable target in CMake, or temporarily
//     compile it in your existing binary target for quick verification.
//
// NOTE:
//   This intentionally does not depend on star structure; it validates parsing,
//   interpolation, and derivative stencils in isolation.
//------------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CompactStar/EOS/CompOSE_Thermo.hpp"

// Your path utilities
#include "Zaki/String/Directory.hpp"

int main()
{
	try
	{
		//--------------------------------------------------------------------------
		// Resolve path roots relative to this translation unit.
		const Zaki::String::Directory this_file_dir(__FILE__);

		// Canonical run paths under "<...>/results/<out_dir>/"
		const Zaki::String::Directory base_results_dir = this_file_dir.ParentDir() + "/results";
		const Zaki::String::Directory out_dir = "compose_therm";

		// EOS root location
		Zaki::String::Directory eos_root =
			this_file_dir.ParentDir().ParentDir() + "/EOS/CompOSE/";
		const std::string eos_name = "CMF-1_general";

		const Zaki::String::Directory eos_dir = eos_root + "/" + eos_name;

		std::cout << "============================================================\n";
		std::cout << "CompOSE_Thermo sanity test\n";
		std::cout << "------------------------------------------------------------\n";
		std::cout << "this_file_dir     : " << this_file_dir.Str() << "\n";
		std::cout << "base_results_dir  : " << base_results_dir.Str() << "\n";
		std::cout << "out_dir           : " << out_dir.Str() << "\n";
		std::cout << "eos_dir           : " << eos_dir.Str() << "\n";
		std::cout << "============================================================\n\n";

		//--------------------------------------------------------------------------
		// Load thermo
		CompactStar::EOS::CompOSE_Thermo::Options opt;
		opt.use_central_difference = true;
		opt.clamp_to_domain = true;

		// For very low-T calls (e.g. T << 2 MeV), it can be helpful to anchor
		// the derivative to the first interval (0 -> T1). You can toggle this:
		// opt.Tmin_for_derivative_MeV = 2.0;
		opt.Tmin_for_derivative_MeV = 0.0;

		CompactStar::EOS::CompOSE_Thermo thermo(eos_dir.Str(), opt);

		if (!thermo.IsLoaded())
			throw std::runtime_error("Thermo table did not load (IsLoaded() == false).");

		//--------------------------------------------------------------------------
		// Print grid summary
		const auto &Tg = thermo.TGrid_MeV();
		const auto &Nbg = thermo.NbGrid_fm3();
		const auto &Yqg = thermo.YqGrid();

		std::cout << "Grid sizes:\n";
		std::cout << "  NT  = " << thermo.NT() << "\n";
		std::cout << "  NNb = " << thermo.NNb() << "\n";
		std::cout << "  NYq = " << thermo.NYq() << "\n\n";

		std::cout << "Grid ranges:\n";
		std::cout << "  T   : [" << Tg.front() << ", " << Tg.back() << "] MeV\n";
		std::cout << "  nB  : [" << Nbg.front() << ", " << Nbg.back() << "] fm^-3\n";
		std::cout << "  Yq  : [" << Yqg.front() << ", " << Yqg.back() << "]\n\n";

		//--------------------------------------------------------------------------
		// Choose a few sample points for sanity checks.
		//
		// Strategy:
		//  - midpoints in nB and Yq to avoid boundary artifacts
		//  - a few temperatures: 0, 2, 8, 50 MeV (or clamped if outside)
		//
		const auto pick_mid = [](const std::vector<double> &g) -> double
		{
			return g[g.size() / 2];
		};

		const double nb_mid = pick_mid(Nbg);
		const double yq_mid = pick_mid(Yqg);

		const int precision = 5;
		// Some temperature probes; will be clamped if outside domain
		const std::vector<double> T_probe = {0.0, 0.2, 0.5, 1.0, 2.0, 4.0, 10.0, 50.0, 80.0, 160.0};

		std::cout << "Sample evaluation point (mid-grid):\n";
		std::cout << "  nB = " << nb_mid << " fm^-3\n";
		std::cout << "  Yq = " << yq_mid << "\n\n";

		std::cout << std::scientific << std::setprecision(precision);

		std::cout << "Thermo samples at mid-grid:\n";
		std::cout << "  Columns:\n T[MeV], Q2=s/nB, dQ2/dT[1/MeV], CvDensity[erg K^-1 cm^-3], CvPerBaryon, dQ2_Cool/dT, Cv_Cool, Q2_Cool\n";
		for (double T : T_probe)
		{
			const double q2 = thermo.Q2(T, nb_mid, yq_mid);
			const double dq2 = thermo.dQ2dT(T, nb_mid, yq_mid);
			const double cv = thermo.CvDensity_cgs(T, nb_mid, yq_mid);
			const double cvb = thermo.CvPerBaryon(T, nb_mid, yq_mid);
			const double dq2_cool = thermo.dQ2dT_ForCooling(T, nb_mid, yq_mid);
			const double cv_cool = thermo.CvDensity_cgs_ForCooling(T, nb_mid, yq_mid);
			const double q2_cool = thermo.Q2_ForCooling(T, nb_mid, yq_mid);

			std::cout << "  "
					  << std::setw(precision) << T
					  << "  " << std::setw(precision) << q2
					  << "  " << std::setw(precision) << dq2
					  << "  " << std::setw(precision) << cv
					  << "  " << std::setw(precision) << cvb
					  << "  " << std::setw(precision) << dq2_cool
					  << "  " << std::setw(precision) << cv_cool
					  << "  " << std::setw(precision) << q2_cool

					  << "\n";
		}
		std::cout << "\n";

		//--------------------------------------------------------------------------
		// Additional sanity checks:
		//  1) Check that Q2(0,...) is near ~0 (not guaranteed for all models, but typical).
		//  2) Check that CvDensity >= 0 for a handful of points (should generally be true).
		//
		// We'll evaluate a small set across the grid (low, mid, high) for nB and Yq.
		//
		const std::vector<double> nb_probe = {Nbg.front(), nb_mid, Nbg.back()};
		const std::vector<double> yq_probe = {Yqg.front(), yq_mid, Yqg.back()};
		const double T_check = (Tg.size() > 1 ? Tg[1] : Tg.front()); // first nonzero if available

		std::cout << "Positivity spot-check at T = " << T_check << " MeV:\n";
		for (double nb : nb_probe)
		{
			for (double yq : yq_probe)
			{
				const double cv = thermo.CvDensity_cgs(T_check, nb, yq);
				std::cout << "  nB=" << std::setw(precision) << nb
						  << "  Yq=" << std::setw(precision) << yq
						  << "  Cv=" << std::setw(precision) << cv
						  << (cv < 0.0 ? "  [WARNING: negative]" : "")
						  << "\n";
			}
		}

		std::cout << "\nDone.\n";
		return 0;
	}
	catch (const std::exception &e)
	{
		std::cerr << "\n[ERROR] " << e.what() << "\n";
		return 1;
	}
}