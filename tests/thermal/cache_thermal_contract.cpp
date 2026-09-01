// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file cache_thermal_contract.cpp
 * @brief Phase 2B-3 — INV-12 cache contracts on the THERMAL path: the StarContext charge
 *        fraction Y_q (observed through C_star) and the NeutrinoCooling profile-versioned
 *        cache. Self-contained; builds its own synthetic CompOSE fixture.
 *
 *   usage: cache_thermal_contract
 *
 * SCOPE. Deliberately does NOT duplicate what `heat_capacity_v1` (U7) already covers
 * durably: repeated-query stability of C_star, its rebuild on profile-version change, and
 * its rebuild on a different CompOSE_Thermo object. Those remain that test's contracts and
 * are re-authenticated by running it. This file adds the two thermal cache dependencies
 * that no existing test exercises — Y_q, and the NeutrinoCooling payload — plus a
 * quantified reproduction of the geometry-key hazard.
 *
 * All checks are ordinary pass/fail contracts. The two cases in `RunProvenanceContracts`
 * were previously defect reproductions behind `--audit-known-hazards`; ADR-0003 (ACCEPTED)
 * makes them requirements, so that mode is gone.
 *
 * FIXTURE. Q2 = s/n_B = slope * T * (1 + Yq), so
 *      c_V = T * n_B * dQ2/dT = T * n_B * slope * (1 + Yq)
 * exactly, making C_star strictly proportional to (1 + Yq) and therefore a clean probe of
 * the Y_q cache, which StarContext keeps private. With nu = Lambda = 0 the NeutrinoCooling
 * coefficients reduce to K ∝ ∫ rho * 4 pi r^2 dr, i.e. exactly linear in the energy density,
 * so a scaling mutation has an exact predicted effect. Synthetic throughout; asserts no
 * neutron-star property.
 *
 * Y_q is Σ q_i Y_i over STRONG-SECTOR species only (StarContext::BuildYqCache_); electrons
 * (code "0") are not summed. With species {10, 11, 0} this gives Y_q = Y_p exactly.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp"
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Details.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateVector.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::StarProfile;
using CompactStar::EOS::CompOSE_Thermo;
using CompactStar::Physics::Evolution::DriverContext;
using CompactStar::Physics::Evolution::GeometryCache;
using CompactStar::Physics::Evolution::StarContext;
using CompactStar::Physics::Evolution::StateVector;
namespace Th = CompactStar::Physics::Driver::Thermal;

static const double kSlope = 1.0e-2; // 1/MeV
static const double kR_km = 10.0;
static const double kNb = 0.30; // fm^-3, a grid node
static const std::size_t kN = 401;

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
/// Scientific formatting: std::to_string would print 1e-15 as 0.000000.
static std::string Sci(double v, int prec = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
}

// ---------------------------------------------------------------------------
//  Synthetic CompOSE fixture: Q2 = slope * T * (1 + yq)
// ---------------------------------------------------------------------------
static void WriteAxis(const fs::path &p, const std::vector<double> &v)
{
	std::ofstream out(p);
	out << "1 " << v.size() << "\n";
	out << std::setprecision(17);
	for (double x : v)
		out << x << "\n";
}

static void WriteThermo(const fs::path &dir, double slope)
{
	fs::create_directories(dir);
	const std::vector<double> T = {0.0, 0.25, 0.5, 0.75, 1.0};
	const std::vector<double> nb = {0.05, 0.30, 0.60, 1.00};
	const std::vector<double> yq = {0.0, 0.10, 0.30, 0.50};
	WriteAxis(dir / "eos.t", T);
	WriteAxis(dir / "eos.nb", nb);
	WriteAxis(dir / "eos.yq", yq);

	std::ofstream th(dir / "eos.thermo");
	th << "939.565 938.272 0\n";
	th << std::setprecision(17);
	for (std::size_t it = 0; it < T.size(); ++it)
		for (std::size_t ib = 0; ib < nb.size(); ++ib)
			for (std::size_t iy = 0; iy < yq.size(); ++iy)
			{
				const double Q2 = slope * T[it] * (1.0 + yq[iy]);
				th << (it + 1) << " " << (ib + 1) << " " << (iy + 1)
				   << " 0 " << Q2 << " 0 0 0 0 0 0\n";
			}
}

// ---------------------------------------------------------------------------
//  Synthetic star. nu = Lambda = 0; species 10 (n), 11 (p), 0 (e).
// ---------------------------------------------------------------------------
static void FillProfile(StarProfile &prof, double Yp, double Yn, double eps_km2)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	radial.ClearRows();
	radial.Reserve(11, kN);

	radial[0].SetLabel("r(km)");
	prof.SetColumnIndex(StarProfile::Column::Radius, 0);
	radial[1].SetLabel("m(km)");
	prof.SetColumnIndex(StarProfile::Column::Mass, 1);
	radial[2].SetLabel("nu_prime(km^-1)");
	prof.SetColumnIndex(StarProfile::Column::MetricNuPrime, 2);
	radial[3].SetLabel("p(km^-2)");
	prof.SetColumnIndex(StarProfile::Column::Pressure, 3);
	radial[4].SetLabel("eps(km^-2)");
	prof.SetColumnIndex(StarProfile::Column::EnergyDensity, 4);
	radial[5].SetLabel("nB(fm^-3)");
	prof.SetColumnIndex(StarProfile::Column::BaryonDensity, 5);
	radial[6].SetLabel("nu");
	prof.SetColumnIndex(StarProfile::Column::MetricNu, 6);
	radial[7].SetLabel("lambda");
	prof.SetColumnIndex(StarProfile::Column::MetricLambda, 7);

	prof.ResetSpecies(3);
	radial[8].SetLabel("10");
	prof.SetSpeciesColumn("10", 8);
	radial[9].SetLabel("11");
	prof.SetSpeciesColumn("11", 9);
	radial[10].SetLabel("0");
	prof.SetSpeciesColumn("0", 10);

	const double h = kR_km / static_cast<double>(kN - 1);
	for (std::size_t i = 0; i < kN; ++i)
	{
		radial[0].PushBack(static_cast<double>(i) * h);
		radial[1].PushBack(0.0);
		radial[2].PushBack(0.0);
		radial[3].PushBack(0.0);
		radial[4].PushBack(eps_km2);
		radial[5].PushBack(kNb);
		radial[6].PushBack(0.0);
		radial[7].PushBack(0.0);
		radial[8].PushBack(Yn);
		radial[9].PushBack(Yp);
		radial[10].PushBack(Yp); // charge neutrality: Y_e = Y_p
	}
}

static void SetSpeciesInPlace(StarProfile &prof, int col, double v)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	for (std::size_t i = 0; i < kN; ++i)
		radial[col][i] = v;
}

static void ScaleEpsInPlace(StarProfile &prof, double f)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	for (std::size_t i = 0; i < kN; ++i)
		radial[4][i] *= f;
}

/// Build a StateVector carrying one thermal DOF at T_inf.
struct ThermalHolder
{
	CompactStar::Physics::State::ThermalState thermal;
	StateVector Y;
	explicit ThermalHolder(double Tinf_K)
	{
		thermal.Resize(1);
		thermal.SetTinf(Tinf_K);
		Y.Register(CompactStar::Physics::State::StateTag::Thermal, thermal);
	}
};

static Th::NeutrinoCooling::Options CoolingOpts()
{
	Th::NeutrinoCooling::Options no;
	no.include_direct_urca = true;
	no.include_modified_urca = true;
	no.include_pair_breaking = false;
	no.global_scale = 1.0;
	return no;
}

// ===========================================================================
//  SUPPORTED CONTRACTS
// ===========================================================================
static void RunSupportedContracts(const fs::path &root)
{
	CompOSE_Thermo::Options opt;
	CompOSE_Thermo thermo((root / "yqdep").string(), opt);
	if (!thermo.IsLoaded())
	{
		Report("T0 synthetic thermo fixture loads", false, "CompOSE_Thermo not loaded");
		return;
	}
	Report("T0 synthetic thermo fixture loads", true, "Q2 = slope*T*(1+Yq)");

	const double T = 1.0e-4; // MeV

	// -----------------------------------------------------------------------
	// T1 — the Y_q cache rebuilds, observed through C_star
	// -----------------------------------------------------------------------
	std::cout << "\nT1 StarContext charge-fraction (Y_q) cache\n";
	{
		StarProfile prof;
		FillProfile(prof, /*Yp=*/0.10, /*Yn=*/0.85, /*eps=*/1.0e-4);
		StarContext sc(prof);

		const std::uint64_t v0 = sc.ProfileVersion();
		const double C0 = sc.HeatCapacityStar_Tinf(T, thermo, nullptr);
		Report("T1a C_star builds with Y_q = Y_p = 0.10", C0 > 0.0,
			   "C_star = " + std::to_string(C0) + " erg/K");

		// Stability at a fixed version.
		bool stable = true;
		for (int k = 0; k < 4; ++k)
			stable = stable && sc.HeatCapacityStar_Tinf(T, thermo, nullptr) == C0;
		Report("T1b repeated queries at a fixed version are bit-identical", stable,
			   "4 repeats");

		// Change ONE charged strong-sector fraction: Y_p 0.10 -> 0.30 (both yq grid nodes).
		// Keep Y_e = Y_p for charge neutrality; Y_q is the strong-sector sum, so Y_q = Y_p.
		SetSpeciesInPlace(prof, 9, 0.30);  // proton
		SetSpeciesInPlace(prof, 10, 0.30); // electron
		const std::uint64_t v1 = sc.ProfileVersion();
		Report("T1c the species mutation increments Version()", v1 > v0,
			   "version " + std::to_string(v0) + " -> " + std::to_string(v1));

		const double C1 = sc.HeatCapacityStar_Tinf(T, thermo, nullptr);
		const double expected = C0 * (1.0 + 0.30) / (1.0 + 0.10);
		Report("T1d C_star rebuilds and tracks (1 + Y_q) exactly",
			   Rel(C1, expected) < 1e-9,
			   "C_star " + std::to_string(C0) + " -> " + std::to_string(C1) +
				   " (expected " + std::to_string(expected) + ", rel dev " +
				   std::to_string(Rel(C1, expected)) + ")");

		// Back down, to rule out a one-way latch.
		SetSpeciesInPlace(prof, 9, 0.10);
		SetSpeciesInPlace(prof, 10, 0.10);
		Report("T1e C_star rebuilds back",
			   Rel(sc.HeatCapacityStar_Tinf(T, thermo, nullptr), C0) < 1e-9,
			   "returned to " + std::to_string(sc.HeatCapacityStar_Tinf(T, thermo, nullptr)));
	}

	// -----------------------------------------------------------------------
	// T2 — NeutrinoCooling profile-versioned cache, same-star semantics
	// -----------------------------------------------------------------------
	std::cout << "\nT2 NeutrinoCooling cache — same-star rebuild (the contract it was "
				 "designed for)\n";
	{
		StarProfile prof;
		FillProfile(prof, /*Yp=*/0.15, /*Yn=*/0.85, /*eps=*/1.0e-4); // DU allowed
		StarContext sc(prof);
		GeometryCache geo0(sc);

		DriverContext ctx;
		ctx.star = &sc;
		ctx.geo = &geo0;
		ctx.thermo = &thermo;

		Th::NeutrinoCooling nu(CoolingOpts());
		ThermalHolder st(1.0e8);

		const auto d1 = Th::Detail::NeutrinoCooling_Details::ComputeDerived(nu, st.Y, ctx);
		Report("T2a diagnostics compute cleanly", d1.ok && d1.L_nu_inf_erg_s > 0.0,
			   d1.ok ? "L_nu = " + std::to_string(d1.L_nu_inf_erg_s) + " erg/s, n_zones = " +
						   std::to_string(d1.n_zones)
					 : d1.message);

		// (2,3) deterministic repeat at the same state and context
		bool det = true;
		for (int k = 0; k < 4; ++k)
		{
			const auto dr = Th::Detail::NeutrinoCooling_Details::ComputeDerived(nu, st.Y, ctx);
			det = det && dr.L_nu_inf_erg_s == d1.L_nu_inf_erg_s;
		}
		Report("T2b repeated queries at the same state/context are bit-identical", det,
			   "4 repeats");

		// (4) sanctioned in-place mutation affecting the coefficient
		const double kFactor = 3.0;
		ScaleEpsInPlace(prof, kFactor);

		// (5,6) fresh GeometryCache consistent with the new profile, wired into the context
		GeometryCache geo1(sc);
		ctx.geo = &geo1;

		// (7,8) same driver instance -> version-driven rebuild
		const auto d2 = Th::Detail::NeutrinoCooling_Details::ComputeDerived(nu, st.Y, ctx);
		// With nu = Lambda = 0 the coefficients are K ∝ ∫ rho * 4 pi r^2 dr, exactly linear
		// in eps, so L_nu must scale by the same factor at fixed T_inf.
		Report("T2c the profile-versioned cache rebuilds and L_nu scales exactly with rho",
			   d2.ok && Rel(d2.L_nu_inf_erg_s, kFactor * d1.L_nu_inf_erg_s) < 1e-9,
			   "L_nu " + std::to_string(d1.L_nu_inf_erg_s) + " -> " +
				   std::to_string(d2.L_nu_inf_erg_s) + " (expected x" +
				   std::to_string(kFactor) + ", rel dev " +
				   std::to_string(Rel(d2.L_nu_inf_erg_s, kFactor * d1.L_nu_inf_erg_s)) + ")");

		// A composition mutation that closes the DU channel must also rebuild.
		SetSpeciesInPlace(prof, 9, 0.05);
		SetSpeciesInPlace(prof, 10, 0.05); // Y_n = 0.85 > 8*0.05 -> DU forbidden
		GeometryCache geo2(sc);
		ctx.geo = &geo2;
		const auto d3 = Th::Detail::NeutrinoCooling_Details::ComputeDerived(nu, st.Y, ctx);
		Report("T2d closing the DU channel rebuilds the cache and removes the DU luminosity",
			   d3.ok && d3.L_nu_DU_inf_erg_s == 0.0 && d2.L_nu_DU_inf_erg_s > 0.0,
			   "L_nu_DU " + std::to_string(d2.L_nu_DU_inf_erg_s) + " -> " +
				   std::to_string(d3.L_nu_DU_inf_erg_s) + ", MU unchanged at " +
				   std::to_string(d3.L_nu_MU_inf_erg_s));
	}
}

// ===========================================================================
//  ADR-0003 ENFORCED CONTRACTS — formerly the "known hazard" reproductions
//
//  Before Phase 3B these two cases reproduced measured DEFECTS (50 % and 80 % relative
//  error) and were kept out of CTest because they could not pass. ADR-0003 (ACCEPTED) makes
//  each one a requirement. The historical measurements survive in
//  docs/validation/CACHE_CORRECTNESS.md, labelled superseded.
// ===========================================================================
static void RunProvenanceContracts(const fs::path &root)
{
	CompOSE_Thermo::Options opt;
	CompOSE_Thermo thermo((root / "yqdep").string(), opt);
	if (!thermo.IsLoaded())
	{
		Report("T0b fixture loads for the provenance contracts", false, "not loaded");
		return;
	}
	const double T = 1.0e-4;

	std::cout << "\nADR-0003 ENFORCED PROVENANCE CONTRACTS\n";

	// -----------------------------------------------------------------------
	// T3 (was HAZARD B) — geometry participates in the C_star validity condition
	// -----------------------------------------------------------------------
	std::cout << "\nT3 C_star honours the supplied GeometryCache\n";
	{
		StarProfile probe;
		FillProfile(probe, 0.10, 0.85, 1.0e-4);
		StarContext sc(probe);
		GeometryCache geo_ok(sc);

		// A geometry built from a DIFFERENT profile (thicker: e^Lambda = 2).
		StarProfile fat;
		FillProfile(fat, 0.10, 0.85, 1.0e-4);
		{
			auto edit = fat.Edit();
			auto &radial = fat.RadialMutable();
			for (std::size_t i = 0; i < kN; ++i)
				radial[7][i] = std::log(2.0);
		}
		StarContext sc_fat(fat);
		GeometryCache geo_foreign(sc_fat);

		const double c_ok = sc.HeatCapacityStar_Tinf(T, thermo, &geo_ok);
		Report("T3a a matching geometry is accepted", c_ok > 0.0,
			   "C_star = " + Sci(c_ok));

		bool threw = false;
		try { (void)sc.HeatCapacityStar_Tinf(T, thermo, &geo_foreign); }
		catch (const std::exception &) { threw = true; }
		Report("T3b a geometry from a different profile FAILS CLOSED", threw,
			   threw ? "threw as required (previously returned a 50 %-wrong C_star)"
					 : "NO THROW — stale geometry accepted");

		// An equivalent geometry rebuilt for the SAME (profile, version) is interchangeable:
		// provenance, not object address, is the key (ADR-0003 §11).
		GeometryCache geo_equiv(sc);
		const double c_equiv = sc.HeatCapacityStar_Tinf(T, thermo, &geo_equiv);
		Report("T3c an equivalent geometry with identical provenance is accepted and agrees",
			   c_equiv == c_ok,
			   "same value; a different C++ object does not force a spurious rebuild");
	}

	// -----------------------------------------------------------------------
	// T4 (was HAZARD D) — a reused driver must not serve one star's payload for another
	// -----------------------------------------------------------------------
	std::cout << "\nT4 NeutrinoCooling across two equal-version stars\n";
	{
		StarProfile A, B;
		FillProfile(A, 0.15, 0.85, 1.0e-4);
		FillProfile(B, 0.15, 0.85, 5.0e-4); // 5x the energy density
		StarContext sa(A), sb(B);
		GeometryCache ga(sa), gb(sb);

		DriverContext ca;
		ca.star = &sa; ca.geo = &ga; ca.thermo = &thermo;
		DriverContext cb;
		cb.star = &sb; cb.geo = &gb; cb.thermo = &thermo;
		ThermalHolder st(1.0e8);

		Report("T4a precondition: the two stars share a numeric Version()",
			   A.Version() == B.Version(),
			   "both " + std::to_string(A.Version()));

		Th::NeutrinoCooling reused(CoolingOpts());
		const auto rA = Th::Detail::NeutrinoCooling_Details::ComputeDerived(reused, st.Y, ca);
		const auto rB = Th::Detail::NeutrinoCooling_Details::ComputeDerived(reused, st.Y, cb);

		Th::NeutrinoCooling fresh(CoolingOpts());
		const auto fB = Th::Detail::NeutrinoCooling_Details::ComputeDerived(fresh, st.Y, cb);

		Report("T4b the reused driver returns star B's OWN luminosity, not star A's",
			   rA.ok && rB.ok && fB.ok && rB.L_nu_inf_erg_s == fB.L_nu_inf_erg_s &&
				   rB.L_nu_inf_erg_s != rA.L_nu_inf_erg_s,
			   "reused B = " + Sci(rB.L_nu_inf_erg_s) + ", fresh B = " +
				   Sci(fB.L_nu_inf_erg_s) + ", A = " + Sci(rA.L_nu_inf_erg_s));

		// Returning to star A must rebuild again, not serve B's payload.
		const auto rA2 = Th::Detail::NeutrinoCooling_Details::ComputeDerived(reused, st.Y, ca);
		Report("T4c switching back rebuilds again", rA2.ok &&
			   rA2.L_nu_inf_erg_s == rA.L_nu_inf_erg_s,
			   "A revisited = " + Sci(rA2.L_nu_inf_erg_s));

		// A geometry that does not belong to ctx.star must fail closed, through the
		// driver's own diagnostic mechanism (no termination).
		DriverContext bad;
		bad.star = &sa; bad.geo = &gb; bad.thermo = &thermo;
		const auto rbad = Th::Detail::NeutrinoCooling_Details::ComputeDerived(reused, st.Y, bad);
		Report("T4d a mismatched profile/geometry pair fails closed rather than computing",
			   !rbad.ok, rbad.ok ? "COMPUTED ANYWAY" : "ok=false: " + rbad.message);
	}
}

int main(int argc, char **argv)
{
	if (argc < 1)
		return 2;
	(void)argv;

	std::cout << std::scientific << std::setprecision(6);
	const fs::path root = fs::temp_directory_path() / "compactstar_cache_thermal";
	fs::remove_all(root);
	WriteThermo(root / "yqdep", kSlope);

	std::cout << "THERMAL CACHE CONTRACTS (Y_q via C_star; NeutrinoCooling payload)\n"
				 "Synthetic fixture; asserts no neutron-star property.\n"
				 "C_star repeat/version/thermo-identity contracts live in heat_capacity_v1 "
				 "(U7).\n\n";
	RunSupportedContracts(root);
	RunProvenanceContracts(root);

	std::cout << "\n"
			  << (g_fail == 0
					  ? "supported thermal cache contracts and ADR-0003 provenance contracts hold"
					  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	fs::remove_all(root);
	return g_fail == 0 ? 0 : 1;
}
