// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file photon_cooling_conformance.cpp
 * @brief ADR-0002 source-conformance test for PhotonCooling (roadmap Phase 2A).
 *
 * Asserts that the photon channel divides by the SAME canonical C_*(T_inf) that
 * NeutrinoCooling uses — ADR-0002 Pattern A — rather than a driver-local constant.
 *
 * Self-contained: builds a synthetic CompOSE fixture in a temp directory and needs no
 * external EOS data. The fixture is deliberately non-physical; this asserts no
 * scientific value, only the governed denominator identity.
 *
 * This test FAILS against the pre-change implementation, which divided by
 * PhotonCooling::Options::C_eff = 1e40.
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::StarProfile;
using CompactStar::EOS::CompOSE_Thermo;
using CompactStar::Physics::Driver::Thermal::PhotonCooling;
using CompactStar::Physics::Evolution::DriverContext;
using CompactStar::Physics::Evolution::GeometryCache;
using CompactStar::Physics::Evolution::StarContext;
using CompactStar::Physics::Evolution::StateVector;
using CompactStar::Physics::State::ThermalState;
namespace Detail = CompactStar::Physics::Driver::Thermal::Detail;

constexpr double kSlope = 0.25;
constexpr double kR_km = 10.0;
constexpr double kNb = 0.30;
constexpr double kYq = 0.10;
constexpr double kTinf_K = 1.0e8;
constexpr double kMEV_PER_K = 8.617333262145e-11; // matches the drivers
constexpr double kOldPlaceholder = 1.0e40;		  // the rejected PhotonCooling C_eff

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}
static void Note(const std::string &id, const std::string &d)
{
	std::cout << "  [RECORD] " << id << " — " << d << "\n";
}

static void WriteAxis(const fs::path &p, const std::vector<double> &v)
{
	std::ofstream o(p);
	o << "1 " << v.size() << "\n" << std::setprecision(17);
	for (double x : v)
		o << x << "\n";
}

static void WriteFixture(const fs::path &dir, double slope)
{
	fs::create_directories(dir);
	const std::vector<double> T = {0.0, 0.25, 0.5, 0.75, 1.0};
	const std::vector<double> nb = {0.05, 0.30, 0.60, 1.00};
	const std::vector<double> yq = {0.0, 0.10, 0.30, 0.50};
	WriteAxis(dir / "eos.t", T);
	WriteAxis(dir / "eos.nb", nb);
	WriteAxis(dir / "eos.yq", yq);
	std::ofstream th(dir / "eos.thermo");
	th << "939.565 938.272 1\n" << std::setprecision(17);
	for (std::size_t it = 0; it < T.size(); ++it)
		for (std::size_t ib = 0; ib < nb.size(); ++ib)
			for (std::size_t iy = 0; iy < yq.size(); ++iy)
				th << (it + 1) << " " << (ib + 1) << " " << (iy + 1)
				   << " 0 " << (slope * T[it]) << " 0 0 0 0 0 0\n";
}

static void FillProfile(StarProfile &prof, std::size_t N)
{
	auto edit = prof.Edit();
	auto &rad = prof.RadialMutable();
	rad.ClearRows();
	rad.Reserve(9, N);
	const char *labels[] = {"r(km)", "m(km)", "nu_prime(km^-1)", "p(km^-2)",
							"eps(km^-2)", "nB(fm^-3)", "nu", "lambda"};
	using C = StarProfile::Column;
	const C cols[] = {C::Radius, C::Mass, C::MetricNuPrime, C::Pressure,
					  C::EnergyDensity, C::BaryonDensity, C::MetricNu, C::MetricLambda};
	for (int i = 0; i < 8; ++i)
	{
		rad[i].SetLabel(labels[i]);
		prof.SetColumnIndex(cols[i], i);
	}
	prof.ResetSpecies(1);
	rad[8].SetLabel("11");
	prof.SetSpeciesColumn("11", 8);

	const double h = kR_km / double(N - 1);
	for (std::size_t i = 0; i < N; ++i)
	{
		rad[0].PushBack(double(i) * h);
		for (int c = 1; c <= 4; ++c)
			rad[c].PushBack(0.0);
		rad[5].PushBack(kNb);
		rad[6].PushBack(0.0); // nu = 0
		rad[7].PushBack(0.0); // Lambda = 0
		rad[8].PushBack(kYq);
	}
}

static double RelErr(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout << "ADR-0002 PhotonCooling source-conformance (synthetic fixture)\n\n";

	const fs::path root = fs::temp_directory_path() / "compactstar_pc_conformance";
	fs::remove_all(root);
	WriteFixture(root / "linear", kSlope);
	CompOSE_Thermo::Options topt;
	CompOSE_Thermo thermo((root / "linear").string(), topt);

	StarProfile prof;
	FillProfile(prof, 2001);
	StarContext ctx_star(prof);
	GeometryCache geo(ctx_star);

	ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(kTinf_K);
	StateVector Y;
	Y.Register(CompactStar::Physics::State::StateTag::Thermal, thermal);

	PhotonCooling::Options opts;
	opts.surface_model = PhotonCooling::Options::SurfaceModel::ApproxFromTinf;
	opts.radiating_fraction = 1.0;
	opts.global_scale = 1.0;
	PhotonCooling drv(opts);

	DriverContext ctx;
	ctx.star = &ctx_star;
	ctx.geo = &geo;
	ctx.thermo = &thermo;

	// Independent expectation: the canonical heat capacity, obtained directly.
	const double C_expected =
		ctx_star.HeatCapacityStar_Tinf(kTinf_K * kMEV_PER_K, thermo, &geo);

	// ---------------------------------------------------------------
	std::cout << "P1  denominator identity\n";
	const auto d = Detail::ComputeDerived(drv, Y, ctx);
	Report("P1.a ComputeDerived ok", d.ok, d.message.empty() ? "(no message)" : d.message);
	std::cout << "      T_inf              = " << d.Tinf_K << " K\n"
			  << "      L_gamma_inf        = " << d.L_gamma_inf_erg_s << " erg/s\n"
			  << "      C_star_erg_K       = " << d.C_star_erg_K << " erg/K\n"
			  << "      C_expected         = " << C_expected << " erg/K\n"
			  << "      dTinf/dt           = " << d.dTinf_dt_K_s << " K/s\n"
			  << "      dlnTinf/dt         = " << d.dLnTinf_dt_1_s << " 1/s\n";
	Report("P1.b C_star_erg_K == HeatCapacityStar_Tinf", RelErr(d.C_star_erg_K, C_expected) < 1e-12,
		   "rel err " + std::to_string(RelErr(d.C_star_erg_K, C_expected)));
	Report("P1.c dTinf/dt == -L_gamma / C_star",
		   RelErr(d.dTinf_dt_K_s, -d.L_gamma_inf_erg_s / C_expected) < 1e-12, "exact");
	Report("P1.d dlnTinf/dt == dTinf/dt / T_inf",
		   RelErr(d.dLnTinf_dt_1_s, d.dTinf_dt_K_s / d.Tinf_K) < 1e-12, "exact");

	// ---------------------------------------------------------------
	std::cout << "\nP2  sensitivity: the old 1e40 denominator must NOT reproduce this\n";
	{
		const double old_rate = -d.L_gamma_inf_erg_s / kOldPlaceholder;
		const double ratio = d.dTinf_dt_K_s / old_rate;
		std::cout << "      dTinf/dt with C_star   = " << d.dTinf_dt_K_s << " K/s\n"
				  << "      dTinf/dt with 1e40     = " << old_rate << " K/s\n"
				  << "      ratio                  = " << ratio << "\n";
		Report("P2 corrected rate differs decisively from the 1e40 rate",
			   std::fabs(ratio - 1.0) > 1.0, "ratio " + std::to_string(ratio));
	}

	// ---------------------------------------------------------------
	std::cout << "\nP3  fail-closed on missing context (was a latent dereference)\n";
	{
		DriverContext c_nostar = ctx;
		c_nostar.star = nullptr;
		const auto r1 = Detail::ComputeDerived(drv, Y, c_nostar);
		Report("P3.a ctx.star == nullptr -> ok=false, no dereference", !r1.ok, r1.message);

		DriverContext c_nothermo = ctx;
		c_nothermo.thermo = nullptr;
		const auto r2 = Detail::ComputeDerived(drv, Y, c_nothermo);
		Report("P3.b ctx.thermo == nullptr -> ok=false", !r2.ok, r2.message);
	}

	// ---------------------------------------------------------------
	std::cout << "\nP4  disabled driver keeps zero semantics and needs no thermodynamics\n";
	{
		for (auto mode : {0, 1})
		{
			PhotonCooling::Options o = opts;
			if (mode == 0)
				o.radiating_fraction = 0.0;
			else
				o.global_scale = 0.0;
			PhotonCooling dis(o);
			DriverContext c = ctx;
			c.star = nullptr; // prove no star/thermo is demanded when contributing zero
			c.thermo = nullptr;
			const auto r = Detail::ComputeDerived(dis, Y, c);
			const std::string tag = (mode == 0) ? "radiating_fraction<=0" : "global_scale<=0";
			Report("P4 " + tag + " -> ok, all-zero, no star/thermo needed",
				   r.ok && r.L_gamma_inf_erg_s == 0.0 && r.dTinf_dt_K_s == 0.0 &&
					   r.dLnTinf_dt_1_s == 0.0,
				   r.message);
		}
		Note("P4 ordering",
			 "pre-existing: SurfaceModel::EnvelopeTbTs still requires ctx.star earlier, for the "
			 "Tb mapping — unchanged by this correction");
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "PhotonCooling conforms to ADR-0002"
							  : "CONFORMANCE FAILURES: " + std::to_string(g_fail))
			  << "\n";
	fs::remove_all(root);
	return g_fail == 0 ? 0 : 1;
}
