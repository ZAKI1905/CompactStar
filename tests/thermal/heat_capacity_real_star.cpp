// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file heat_capacity_real_star.cpp
 * @brief ADR-0002 §V1 Tier-B verification on an authenticated canonical star.
 *
 * VERIFICATION ONLY — changes no production behavior.
 *
 * Requires authenticated official CompOSE data, supplied via the data root passed
 * as argv[1] (wired from the CMake cache variable COMPACTSTAR_EOS_DATA_ROOT). No
 * absolute path is baked into this file. If the data are absent the program exits
 * non-zero: it must not silently skip.
 *
 * Expected layout under the data root:
 *     DS-CMF-1-with-crust/DS(CMF)-1_with_crust.eos     (cold structure, T=0)
 *     DNS-CMF-Hadronic-with-electrons/eos.{t,nb,yq,thermo}   (finite-T thermodynamics)
 *
 * This program reports numbers. It deliberately does NOT assert a literature value
 * as a pass/fail threshold — the magnitude comparison is an order-of-magnitude
 * diagnostic and is adjudicated in docs/validation/HEAT_CAPACITY_V1.md.
 */

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::StarProfile;
using CompactStar::EOS::CompOSE_Thermo;
using CompactStar::Physics::Evolution::GeometryCache;
using CompactStar::Physics::Evolution::StarContext;

// k_B in MeV/K — the value CompOSE_Thermo itself uses (INV-02 records the split).
constexpr double kB_MeV_per_K = 8.617333262e-11;
constexpr double KM3_TO_CM3 = 1.0e15;

static double KtoMeV(double T_K) { return T_K * kB_MeV_per_K; }

int main(int argc, char **argv)
{
	std::cout << std::scientific << std::setprecision(6);
	if (argc < 2)
	{
		std::cerr << "usage: heat_capacity_real_star <COMPACTSTAR_EOS_DATA_ROOT>\n";
		return 2;
	}
	const fs::path root(argv[1]);
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	const fs::path finiteT = root / "DNS-CMF-Hadronic-with-electrons";

	if (!fs::exists(cold) || !fs::exists(finiteT / "eos.thermo"))
	{
		std::cerr << "REQUIRED AUTHENTICATED DATA MISSING under " << root << "\n"
				  << "  expected: " << cold << "\n"
				  << "            " << (finiteT / "eos.thermo") << "\n";
		return 3; // hard failure, never a silent skip
	}

	std::cout << "ADR-0002 V1 Tier-B — authenticated canonical star\n"
			  << "  cold structure : " << cold << "\n"
			  << "  finite-T thermo: " << finiteT << "\n\n";

	// ---------------------------------------------------------------
	// 1) Finite-T table header and grids
	// ---------------------------------------------------------------
	CompOSE_Thermo::Options opt; // production defaults, as the live program uses
	opt.use_central_difference = true;
	opt.clamp_to_domain = true;
	opt.Tmin_for_derivative_MeV = 0.0;
	CompOSE_Thermo thermo(finiteT.string(), opt);
	if (!thermo.IsLoaded())
	{
		std::cerr << "finite-T table failed to load\n";
		return 4;
	}
	const auto &Tg = thermo.TGrid_MeV();
	const auto &nbg = thermo.NbGrid_fm3();
	const auto &yqg = thermo.YqGrid();
	std::cout << "T1  finite-T grids\n"
			  << "      NT=" << Tg.size() << "  T in [" << Tg.front() << ", " << Tg.back() << "] MeV\n"
			  << "      Nnb=" << nbg.size() << "  nB in [" << nbg.front() << ", " << nbg.back() << "] fm^-3\n"
			  << "      NYq=" << yqg.size() << "  Yq in [" << yqg.front() << ", " << yqg.back() << "]\n"
			  << "      first 3 T points (the default low-T fit stencil): "
			  << Tg[0] << ", " << Tg[1] << ", " << Tg[2] << " MeV\n\n";

	// ---------------------------------------------------------------
	// 2) Canonical star through the production TOV path
	// ---------------------------------------------------------------
	const double target_M = 1.4;
	NStar ns;
	ns.SetWrkDir(fs::temp_directory_path().string());
	const int n_rows = ns.SolveTOV_Profile(cold.string(), target_M, "hc_tierB");
	if (n_rows <= 0)
	{
		std::cerr << "SolveTOV_Profile failed\n";
		return 5;
	}
	const auto &seq = ns.GetSequence();
	StarContext ctx(ns.Profile());
	GeometryCache geo(ctx); // ONE stable GeometryCache for all queries (INV-12 hazard avoided)

	const auto *rcol = ns.Profile().GetPtr(StarProfile::Column::Radius);
	const auto *nbcol = ns.Profile().GetPtr(StarProfile::Column::BaryonDensity);
	std::cout << "T2  canonical star (production TOV path)\n"
			  << "      requested M = " << target_M << " Msun\n"
			  << "      achieved  M = " << seq.m << " Msun,  R = " << seq.r << " km\n"
			  << "      central eps = " << seq.ec << "\n"
			  << "      radial points = " << n_rows << ",  profile version = "
			  << ctx.ProfileVersion() << "\n"
			  << "      nB range = [" << (*nbcol)[0] << ", " << (*nbcol)[-1] << "] fm^-3"
			  << "  (central -> surface)\n\n";

	// ---------------------------------------------------------------
	// 3) Domain overlap: cold star vs finite-T table
	// ---------------------------------------------------------------
	const double nb_min_tab = nbg.front();
	const std::size_t N = geo.Size();
	const auto &wv = geo.WV();
	const auto &r = geo.R();
	const auto &emnu = geo.ExpMinusNu();

	std::size_t n_below = 0;
	double r_first_below = -1.0;
	for (std::size_t i = 0; i < N; ++i)
		if ((*nbcol)[i] < nb_min_tab)
		{
			if (r_first_below < 0.0)
				r_first_below = r[i];
			++n_below;
		}
	std::cout << "T3  thermo-domain coverage\n"
			  << "      finite-T table starts at nB = " << nb_min_tab << " fm^-3\n"
			  << "      radial zones with nB below that: " << n_below << " / " << N
			  << " (" << (100.0 * double(n_below) / double(N)) << " %)\n"
			  << "      first such radius = " << r_first_below << " km of R = " << seq.r << " km\n"
			  << "      clamp_to_domain = true, so those zones are evaluated at the table edge\n\n";

	// Reconstruct Y_q = sum_i q_i Y_i over the strong sector, independently of
	// StarContext's private cache, using the same CompOSE codes/charges it uses.
	const std::pair<const char *, double> kStrong[] = {
		{"10", 0.0}, {"11", +1.0}, {"20", -1.0}, {"21", 0.0}, {"22", +1.0}, {"23", +2.0},
		{"100", 0.0}, {"110", -1.0}, {"111", 0.0}, {"112", +1.0}, {"120", -1.0}, {"121", 0.0}};
	std::vector<double> Yq(N, 0.0);
	int n_species_found = 0;
	for (const auto &kv : kStrong)
	{
		const auto *col = ns.Profile().GetSpeciesPtr(kv.first);
		if (!col)
			continue;
		++n_species_found;
		for (std::size_t i = 0; i < N; ++i)
			Yq[i] += kv.second * (*col)[i];
	}
	std::cout << "      strong-sector species found in the cold profile: "
			  << n_species_found << ",  Yq(center)=" << Yq[0]
			  << "  Yq(surface)=" << Yq[N - 1] << "\n";
	{
		double ylo = 1e300, yhi = -1e300;
		for (double v : Yq) { ylo = std::min(ylo, v); yhi = std::max(yhi, v); }
		std::cout << "      Yq range over the star = [" << ylo << ", " << yhi
				  << "] vs table [" << yqg.front() << ", " << yqg.back() << "]\n";
	}

	// Contribution of clamped vs native zones, using the reconstructed Y_q.
	auto integrate = [&](double Tinf_MeV, bool only_clamped, bool only_native) {
		double sum = 0.0;
		for (std::size_t i = 0; i + 1 < N; ++i)
		{
			const bool clamped = ((*nbcol)[i] < nb_min_tab) || ((*nbcol)[i + 1] < nb_min_tab);
			if (only_clamped && !clamped)
				continue;
			if (only_native && clamped)
				continue;
			const double dr = r[i + 1] - r[i];
			const double cv0 = thermo.CvDensity_cgs_ForCooling(Tinf_MeV * emnu[i], (*nbcol)[i], Yq[i]);
			const double cv1 = thermo.CvDensity_cgs_ForCooling(Tinf_MeV * emnu[i + 1], (*nbcol)[i + 1], Yq[i + 1]);
			sum += 0.5 * (cv0 * wv[i] + cv1 * wv[i + 1]) * KM3_TO_CM3 * dr;
		}
		return sum;
	};

	// ---------------------------------------------------------------
	// 4) C_star(T_inf) through production
	// ---------------------------------------------------------------
	std::cout << "T4  production C_star(T_inf)\n"
			  << "      T_inf [K]      T_inf [MeV]    C_star [erg/K]   C_star/T_inf   dlnC/dlnT\n";
	const std::vector<double> TK = {1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9};
	std::vector<double> Cs;
	for (double T : TK)
		Cs.push_back(ctx.HeatCapacityStar_Tinf(KtoMeV(T), thermo, &geo));
	for (std::size_t i = 0; i < TK.size(); ++i)
	{
		double slope = std::nan("");
		if (i + 1 < TK.size())
			slope = std::log(Cs[i + 1] / Cs[i]) / std::log(TK[i + 1] / TK[i]);
		std::cout << "      " << TK[i] << "   " << KtoMeV(TK[i]) << "   " << Cs[i]
				  << "   " << (Cs[i] / TK[i]) << "   " << slope << "\n";
	}

	const double C_1e8 = ctx.HeatCapacityStar_Tinf(KtoMeV(1e8), thermo, &geo);
	const double clamped_1e8 = integrate(KtoMeV(1e8), true, false);
	const double native_1e8 = integrate(KtoMeV(1e8), false, true);
	const double tot_1e8 = clamped_1e8 + native_1e8;
	std::cout << "\nT5  clamped-zone contribution at T_inf = 1e8 K"
			  << "\n"
			  << "      native-domain part = " << native_1e8 << " erg/K ("
			  << (100.0 * native_1e8 / tot_1e8) << " %)\n"
			  << "      clamped part       = " << clamped_1e8 << " erg/K ("
			  << (100.0 * clamped_1e8 / tot_1e8) << " %)\n";

	std::cout << "\nT6  literature magnitude diagnostic (order of magnitude only)\n"
			  << "      C_star(1e8 K) = " << C_1e8 << " erg/K\n"
			  << "      broad NS expectation ~1e37 - 1e38 erg/K\n"
			  << "      C ~ 1e39 * T9 erg/K  =>  ~1e38 at 1e8 K\n"
			  << "      ratio to 1e38 = " << (C_1e8 / 1.0e38) << "\n";

	// ---------------------------------------------------------------
	// 5) Low-T fit sensitivity — production defaults NOT changed
	// ---------------------------------------------------------------
	std::cout << "\nT7  low-T fit sensitivity at T_inf = 1e8 K"
			  << " (lowT_fit_points varied in TEST code only)\n";
	double lo = 1e300, hi = -1e300;
	for (int nfit : {2, 3, 4, 5})
	{
		CompOSE_Thermo::Options o = opt;
		o.lowT_fit_points = nfit;
		CompOSE_Thermo th(finiteT.string(), o);
		StarContext c2(ns.Profile());
		GeometryCache g2(c2);
		const double v = c2.HeatCapacityStar_Tinf(KtoMeV(1e8), th, &g2);
		lo = std::min(lo, v);
		hi = std::max(hi, v);
		std::cout << "      lowT_fit_points=" << nfit << " (T up to "
				  << Tg[std::min<std::size_t>(nfit - 1, Tg.size() - 1)] << " MeV)"
				  << "   C_star(1e8 K) = " << v << " erg/K\n";
	}
	std::cout << "      spread across fit choices = " << (hi / lo) << "x\n";

	// ---------------------------------------------------------------
	// 6) Post-conformance sanity check: the PhotonCooling denominator
	// ---------------------------------------------------------------
	std::cout << "\nT8  PhotonCooling denominator after ADR-0002 conformance\n";
	{
		CompactStar::Physics::State::ThermalState th;
		th.Resize(1);
		th.SetTinf(1e8);
		CompactStar::Physics::Evolution::StateVector Y;
		Y.Register(CompactStar::Physics::State::StateTag::Thermal, th);

		CompactStar::Physics::Driver::Thermal::PhotonCooling::Options po;
		po.surface_model =
			CompactStar::Physics::Driver::Thermal::PhotonCooling::Options::SurfaceModel::ApproxFromTinf;
		CompactStar::Physics::Driver::Thermal::PhotonCooling pc(po);

		CompactStar::Physics::Evolution::DriverContext dctx;
		dctx.star = &ctx;
		dctx.geo = &geo;
		dctx.thermo = &thermo;

		const auto pd =
			CompactStar::Physics::Driver::Thermal::Detail::ComputeDerived(pc, Y, dctx);
		if (!pd.ok)
		{
			std::cout << "      ComputeDerived not ok: " << pd.message << "\n";
		}
		else
		{
			const double old_rate = -pd.L_gamma_inf_erg_s / 1.0e40;
			std::cout << "      C_star used by PhotonCooling = " << pd.C_star_erg_K << " erg/K\n"
					  << "      L_gamma_inf                  = " << pd.L_gamma_inf_erg_s << " erg/s\n"
					  << "      dTinf/dt (governed C_star)   = " << pd.dTinf_dt_K_s << " K/s\n"
					  << "      dTinf/dt (old 1e40)          = " << old_rate << " K/s\n"
					  << "      instantaneous ratio          = " << (pd.dTinf_dt_K_s / old_rate)
					  << "   (expected 1e40 / C_star)\n"
					  << "      NOTE: this is a LOCAL denominator ratio only. The full corrected\n"
					  << "            cooling trajectory is NOT 46x faster -- C_star and L_gamma are\n"
					  << "            both temperature dependent and neutrino cooling acts too.\n";
		}
	}

	std::cout << "\nTier-B numbers produced. Adjudication is in "
			  << "docs/validation/HEAT_CAPACITY_V1.md.\n";
	return 0;
}
