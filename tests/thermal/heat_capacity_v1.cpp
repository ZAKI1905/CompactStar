// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2026
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file heat_capacity_v1.cpp
 * @brief ADR-0002 §V1 Tier-A verification of StarContext::HeatCapacityStar_Tinf.
 *
 * VERIFICATION ONLY. This program changes no production behavior. It exercises the
 * production path
 *
 *     StarProfile -> StarContext -> GeometryCache -> HeatCapacityStar_Tinf
 *
 * against a *synthetic, deliberately non-physical* fixture whose exact answer is known
 * in closed form, and against an independent re-implementation of the ADR-0002 integral
 *
 *     C_star(T_inf) = ∫ c_V(T_local, n_B, Y_q) · 4π r² e^{Λ} dr ,  T_local = T_inf e^{-ν}
 *
 * The fixture is loaded through the real CompOSE_Thermo parser — it is NOT mocked.
 *
 * NOTHING HERE IS A SCIENTIFIC BASELINE. No neutron-star mass, radius, luminosity, EOS
 * value, or literature number is asserted. Magnitude comparison against real neutron-star
 * heat capacities requires the authenticated CompOSE tables and is Tier B, which is not
 * part of this executable. See docs/validation/HEAT_CAPACITY_V1.md.
 */

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"

namespace fs = std::filesystem;

using CompactStar::Core::StarProfile;
using CompactStar::EOS::CompOSE_Thermo;
using CompactStar::Physics::Evolution::GeometryCache;
using CompactStar::Physics::Evolution::StarContext;

// ---------------------------------------------------------------------------
//  Unit constants — must match CompOSE_Thermo::CvDensity_cgs_ForCooling exactly.
//  (CompOSE_Thermo.cpp uses 8.617333262e-11; NeutrinoCooling_Details.cpp uses
//   8.617333262145e-11. That two-precision split is INV-02's recorded k_B issue.
//   We deliberately use the value the code under test uses.)
// ---------------------------------------------------------------------------
constexpr double kB_MeV_per_K = 8.617333262e-11;
constexpr double MeVfm3_to_ergcm3 = 1.602176634e33;
constexpr double KM3_TO_CM3 = 1.0e15;

// Synthetic thermodynamic slope: Q2(T,nB,Yq) = kFixtureSlope · T  (units 1/MeV)
constexpr double kFixtureSlope = 0.25;

// Synthetic star
constexpr double kR_km = 10.0;
constexpr double kNb_fm3 = 0.30;
constexpr double kYq = 0.10;

// ---------------------------------------------------------------------------
//  Pass/fail bookkeeping
// ---------------------------------------------------------------------------
static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &detail)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << detail << "\n";
	if (!ok)
		++g_fail;
}
static void Note(const std::string &id, const std::string &detail)
{
	std::cout << "  [RECORD] " << id << " — " << detail << "\n";
}

// ---------------------------------------------------------------------------
//  Synthetic CompOSE fixture, written in the real on-disk format.
//
//  Axis files:   "<dim> <N> v1 v2 ... vN"     (CompOSE header form)
//  eos.thermo:   header "mn mp Il", then rows
//                "iT inb iYq Q1 Q2 Q3 Q4 Q5 Q6 Q7 Nadd"  (1-based indices)
//  Only Q2 = s/n_B is consumed by CompOSE_Thermo; the rest are filler zeros.
//
//  Q2 = slope · T makes dQ2/dT = slope exactly, on BOTH the low-T least-squares
//  path (a0 = Σ TᵢQᵢ / Σ Tᵢ² = slope) and the tabulated-derivative path, so
//  c_V = T · n_B · slope holds analytically over the whole fixture domain.
// ---------------------------------------------------------------------------
static void WriteAxis(const fs::path &p, const std::vector<double> &v)
{
	std::ofstream out(p);
	out << "1 " << v.size() << "\n";
	out << std::setprecision(17);
	for (double x : v)
		out << x << "\n";
}

static void WriteSyntheticThermo(const fs::path &dir, double slope)
{
	fs::create_directories(dir);

	// T grid must cover the production cache range [1e-5, 1] MeV with room to spare,
	// and must stay entirely below CompOSE_Thermo's default low-T blend window
	// (lowT_switch 2.0 MeV, width 1.0 MeV) so the blend weight is exactly 0.
	const std::vector<double> T = {0.0, 0.25, 0.5, 0.75, 1.0};
	const std::vector<double> nb = {0.05, 0.30, 0.60, 1.00};
	const std::vector<double> yq = {0.0, 0.10, 0.30, 0.50};

	WriteAxis(dir / "eos.t", T);
	WriteAxis(dir / "eos.nb", nb);
	WriteAxis(dir / "eos.yq", yq);

	std::ofstream th(dir / "eos.thermo");
	th << "939.565 938.272 0\n"; // mn mp Il  (Il=0: see report — flag is parsed, never used)
	th << std::setprecision(17);
	for (std::size_t it = 0; it < T.size(); ++it)
		for (std::size_t ib = 0; ib < nb.size(); ++ib)
			for (std::size_t iy = 0; iy < yq.size(); ++iy)
			{
				const double Q2 = slope * T[it]; // entropy per baryon, s/n_B
				th << (it + 1) << " " << (ib + 1) << " " << (iy + 1)
				   << " 0 " << Q2 << " 0 0 0 0 0 0\n";
			}
}

// ---------------------------------------------------------------------------
//  Synthetic star: ν = 0, Λ = 0, n_B = const, one proton species (CompOSE code 11)
//  so ChargeFractionYq() yields Y_q = Y_p.
//
//  With ν = Λ = 0:  e^{-ν} = 1,  WV = 4π r² ,  so
//      C_star(T_inf) = T_inf · n_B · slope · K · 1e15 · ∫₀^R 4π r² dr
//  with K = k_B[MeV/K] · MeVfm3→ergcm3.
// ---------------------------------------------------------------------------
static void FillProfile(StarProfile &prof, std::size_t N, double R_km,
						double nb_fm3, double yq, double lambda = 0.0)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	radial.ClearRows();
	radial.Reserve(9, N);

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

	prof.ResetSpecies(1);
	radial[8].SetLabel("11"); // proton, q = +1
	prof.SetSpeciesColumn("11", 8);

	const double h = R_km / static_cast<double>(N - 1);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = static_cast<double>(i) * h;
		radial[0].PushBack(r);
		radial[1].PushBack(0.0);
		radial[2].PushBack(0.0);
		radial[3].PushBack(0.0);
		radial[4].PushBack(0.0);
		radial[5].PushBack(nb_fm3);
		radial[6].PushBack(0.0);	 // nu = 0
		radial[7].PushBack(lambda);	 // Lambda
		radial[8].PushBack(yq);
	}
}

// Closed-form continuum answer for the fixture (ν = Λ = 0).
static double AnalyticContinuum(double Tinf_MeV, double slope, double nb, double R_km)
{
	const double K = kB_MeV_per_K * MeVfm3_to_ergcm3;
	const double volume_km3 = (4.0 / 3.0) * M_PI * R_km * R_km * R_km;
	return Tinf_MeV * nb * slope * K * KM3_TO_CM3 * volume_km3;
}

// Closed-form value of the *trapezoid* of 4πr² on the uniform grid, so the
// discretization gap is predicted rather than tolerated.
static double AnalyticTrapezoid(double Tinf_MeV, double slope, double nb,
								double R_km, std::size_t N)
{
	const double K = kB_MeV_per_K * MeVfm3_to_ergcm3;
	const double h = R_km / static_cast<double>(N - 1);
	double s = 0.0;
	for (std::size_t i = 0; i + 1 < N; ++i)
	{
		const double r0 = static_cast<double>(i) * h;
		const double r1 = static_cast<double>(i + 1) * h;
		s += 0.5 * (r0 * r0 + r1 * r1) * h;
	}
	return Tinf_MeV * nb * slope * K * KM3_TO_CM3 * 4.0 * M_PI * s;
}

// ---------------------------------------------------------------------------
//  Independent direct evaluation of the ADR-0002 integral.
//  Uses only public production inputs; DOES NOT call HeatCapacityStar_Tinf and
//  does not reuse its 160-point T_inf cache.
// ---------------------------------------------------------------------------
static double DirectIntegral(const StarProfile &prof, const GeometryCache &geo,
							 const CompOSE_Thermo &thermo, double Tinf_MeV)
{
	const auto &r = geo.R();
	const auto &wv = geo.WV();
	const auto &emnu = geo.ExpMinusNu();
	const auto *nb = prof.GetPtr(StarProfile::Column::BaryonDensity);
	// Y_q reconstructed independently as sum_i q_i Y_i (ADR-0001 semantics).
	// The fixture carries a single charged strong-sector species: CompOSE code 11
	// (proton, q = +1). We do NOT call StarContext's private Yq cache.
	const auto *yq = prof.GetSpeciesPtr("11");

	const std::size_t N = geo.Size();
	double sum = 0.0;
	for (std::size_t i = 0; i + 1 < N; ++i)
	{
		const double dr = r[i + 1] - r[i];
		const double cv0 = thermo.CvDensity_cgs_ForCooling(Tinf_MeV * emnu[i], (*nb)[i], (*yq)[i]);
		const double cv1 = thermo.CvDensity_cgs_ForCooling(Tinf_MeV * emnu[i + 1], (*nb)[i + 1], (*yq)[i + 1]);
		sum += 0.5 * (cv0 * wv[i] + cv1 * wv[i + 1]) * KM3_TO_CM3 * dr;
	}
	return sum;
}

static double RelErr(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout << "ADR-0002 V1 Tier-A verification (synthetic fixture; NOT a scientific baseline)\n\n";

	const fs::path root = fs::temp_directory_path() / "compactstar_hcv1_fixture";
	fs::remove_all(root);
	WriteSyntheticThermo(root / "linear", kFixtureSlope);
	WriteSyntheticThermo(root / "linear2x", 2.0 * kFixtureSlope);

	CompOSE_Thermo::Options opt; // defaults: clamp_to_domain, low-T fit on
	CompOSE_Thermo thermo((root / "linear").string(), opt);
	if (!thermo.IsLoaded())
	{
		std::cerr << "fixture failed to load\n";
		return 2;
	}

	const std::size_t N0 = 2001;
	StarProfile prof;
	FillProfile(prof, N0, kR_km, kNb_fm3, kYq);
	StarContext ctx(prof);
	GeometryCache geo(ctx);

	// -------------------------------------------------------------------
	// U1 — unit chain, including the km³→cm³ factor
	// -------------------------------------------------------------------
	std::cout << "U1  unit chain (erg/K) and the 1e15 km^3->cm^3 factor\n";
	{
		// Evaluate the quadrature itself, free of the T_inf-cache interpolation:
		// the direct integral uses no cache, and T = 1e-5 MeV is exactly the first
		// cache node, where production returns the node value without interpolating.
		const double T = 1.0e-5; // MeV — first production cache node
		const double dir = DirectIntegral(prof, geo, thermo, T);
		const double prod = ctx.HeatCapacityStar_Tinf(T, thermo, &geo);
		const double trap = AnalyticTrapezoid(T, kFixtureSlope, kNb_fm3, kR_km, N0);
		const double cont = AnalyticContinuum(T, kFixtureSlope, kNb_fm3, kR_km);
		std::cout << "      direct integral    = " << dir << " erg/K\n"
				  << "      production (node)  = " << prod << " erg/K\n"
				  << "      analytic trapezoid = " << trap << " erg/K\n"
				  << "      analytic continuum = " << cont << " erg/K\n";
		Report("U1.a direct integral == closed-form trapezoid", RelErr(dir, trap) < 1e-12,
			   "rel err " + std::to_string(RelErr(dir, trap)));
		Report("U1.b production at a cache node == closed-form trapezoid",
			   RelErr(prod, trap) < 1e-12, "rel err " + std::to_string(RelErr(prod, trap)));
		// A missing (or wrong-power) 1e15 shows up as a factor 1e15 discrepancy.
		Report("U1.c magnitude excludes a missing 1e15", RelErr(prod, trap * 1e-15) > 1.0,
			   "production is not 1e15x smaller than expected");
		Note("U1.d", "trapezoid vs continuum discretization gap = " +
						 std::to_string(RelErr(trap, cont)) + " at N=" + std::to_string(N0));
	}

	// -------------------------------------------------------------------
	// U2 — low-T scaling  d ln C / d ln T
	// -------------------------------------------------------------------
	std::cout << "\nU2  low-temperature scaling d ln C_star / d ln T_inf\n";
	{
		const double Ta = 1.0e-5, Tb = 1.0e-4; // one decade, inside the cache domain
		const double Ca = ctx.HeatCapacityStar_Tinf(Ta, thermo, &geo);
		const double Cb = ctx.HeatCapacityStar_Tinf(Tb, thermo, &geo);
		const double slope = std::log(Cb / Ca) / std::log(Tb / Ta);
		std::cout << "      C(" << Ta << ") = " << Ca << ",  C(" << Tb << ") = " << Cb << "\n";
		Report("U2 slope ~= 1 for the degenerate fixture", std::fabs(slope - 1.0) < 1e-3,
			   "measured slope " + std::to_string(slope));
		Note("U2 caveat",
			 "the production low-T model fixes dQ2/dT to a constant below the blend window, "
			 "so c_V is linear in T BY CONSTRUCTION; this checks plumbing, not physics");
	}

	// -------------------------------------------------------------------
	// U3 — production cache vs independent direct integral
	// -------------------------------------------------------------------
	std::cout << "\nU3  production cache vs independent direct integral\n";
	{
		double worst = 0.0;
		for (double T : {1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 1.0e-2, 1.0e-1})
		{
			const double prod = ctx.HeatCapacityStar_Tinf(T, thermo, &geo);
			const double dir = DirectIntegral(prof, geo, thermo, T);
			const double e = RelErr(prod, dir);
			worst = std::max(worst, e);
			std::cout << "      T_inf=" << T << "  cached=" << prod
					  << "  direct=" << dir << "  rel=" << e << "\n";
		}
		// Criterion fixed in advance (see report §accuracy criterion): 1%.
		Report("U3 cache vs direct within 1%", worst < 1.0e-2,
			   "max rel err " + std::to_string(worst));
	}

	// -------------------------------------------------------------------
	// U4 — radial convergence: measured, not assumed
	// -------------------------------------------------------------------
	std::cout << "\nU4  radial convergence (order measured)\n";
	{
		// Use the first cache node so the measured error is the RADIAL quadrature
		// error only; at non-node temperatures the T_inf-interpolation floor
		// (~7e-4, see U5) dominates and flattens the observed order to ~0.
		const double T = 1.0e-5;
		const double truth = AnalyticContinuum(T, kFixtureSlope, kNb_fm3, kR_km);
		std::vector<std::size_t> Ns = {101, 201, 401, 801};
		std::vector<double> errs;
		for (std::size_t n : Ns)
		{
			StarProfile p;
			FillProfile(p, n, kR_km, kNb_fm3, kYq);
			StarContext c(p);
			GeometryCache g(c);
			const double v = c.HeatCapacityStar_Tinf(T, thermo, &g);
			const double e = RelErr(v, truth);
			errs.push_back(e);
			std::cout << "      N=" << n << "  rel err=" << e << "\n";
		}
		double worst_order = 1e9;
		for (std::size_t i = 1; i < errs.size(); ++i)
		{
			const double order = std::log(errs[i - 1] / errs[i]) / std::log(2.0);
			std::cout << "      observed order (" << Ns[i - 1] << "->" << Ns[i] << ") = "
					  << order << "\n";
			worst_order = std::min(worst_order, order);
		}
		Report("U4 observed order >= 1.8 (nominal O(dr^2))", worst_order > 1.8,
			   "min observed order " + std::to_string(worst_order));
	}

	// -------------------------------------------------------------------
	// U5 — T-grid (NT) sensitivity, WITHOUT changing production NT
	// -------------------------------------------------------------------
	std::cout << "\nU5  T_inf-grid sensitivity (surrogate caches; production NT untouched)\n";
	{
		auto surrogate_err = [&](std::size_t NT) {
			// Reproduce the production cache concept: log-spaced nodes over
			// [1e-5, 1] MeV, value from the DIRECT integral, then production's
			// interpolation rule (linear in C against log T).
			std::vector<double> Tn(NT), Cn(NT);
			const double lo = std::log(1.0e-5), hi = std::log(1.0);
			for (std::size_t k = 0; k < NT; ++k)
			{
				Tn[k] = std::exp(lo + (hi - lo) * double(k) / double(NT - 1));
				Cn[k] = DirectIntegral(prof, geo, thermo, Tn[k]);
			}
			double worst = 0.0;
			for (double T : {1.7e-5, 5.5e-5, 1.3e-4, 4.4e-4, 2.2e-3, 3.3e-2, 1.7e-1})
			{
				std::size_t i = 0;
				while (i + 2 < NT && Tn[i + 1] < T)
					++i;
				const double w = (std::log(T) - std::log(Tn[i])) / (std::log(Tn[i + 1]) - std::log(Tn[i]));
				const double interp = (1.0 - w) * Cn[i] + w * Cn[i + 1];
				worst = std::max(worst, RelErr(interp, DirectIntegral(prof, geo, thermo, T)));
			}
			return worst;
		};
		for (std::size_t NT : {40u, 80u, 160u, 320u})
			std::cout << "      NT=" << NT << "  max interp rel err=" << surrogate_err(NT) << "\n";

		double worst_prod = 0.0;
		for (double T : {1.7e-5, 5.5e-5, 1.3e-4, 4.4e-4, 2.2e-3, 3.3e-2, 1.7e-1})
			worst_prod = std::max(worst_prod,
								  RelErr(ctx.HeatCapacityStar_Tinf(T, thermo, &geo),
										 DirectIntegral(prof, geo, thermo, T)));
		std::cout << "      production NT=160 max rel err=" << worst_prod << "\n";
		Report("U5 production NT=160 within the 1% criterion", worst_prod < 1.0e-2,
			   "max rel err " + std::to_string(worst_prod));
	}

	// -------------------------------------------------------------------
	// U6 — endpoint clamping: CHARACTERIZED, NOT ENDORSED
	// -------------------------------------------------------------------
	std::cout << "\nU6  endpoint clamping — CURRENT NUMERICAL BEHAVIOR, NOT GOVERNED AS CORRECT\n";
	{
		const double lo = ctx.HeatCapacityStar_Tinf(1.0e-5, thermo, &geo);
		const double below = ctx.HeatCapacityStar_Tinf(1.0e-9, thermo, &geo);
		const double hi = ctx.HeatCapacityStar_Tinf(1.0, thermo, &geo);
		const double above = ctx.HeatCapacityStar_Tinf(50.0, thermo, &geo);
		Report("U6.a below-range returns the low endpoint", RelErr(below, lo) < 1e-12,
			   "C(1e-9 MeV) == C(1e-5 MeV)");
		Report("U6.b above-range returns the high endpoint", RelErr(above, hi) < 1e-12,
			   "C(50 MeV) == C(1 MeV)");
		const double true_below = DirectIntegral(prof, geo, thermo, 1.0e-9);
		Note("U6.c", "at T_inf=1e-9 MeV the clamp overstates C_star by factor " +
						 std::to_string(below / true_below) +
						 " vs the direct integral — recorded for INV-10, not fixed here");
	}

	// -------------------------------------------------------------------
	// U7 — cache correctness
	// -------------------------------------------------------------------
	std::cout << "\nU7  cache behavior\n";
	{
		const double T = 1.0e-4;
		// A: repeated identical query
		const double a1 = ctx.HeatCapacityStar_Tinf(T, thermo, &geo);
		const double a2 = ctx.HeatCapacityStar_Tinf(T, thermo, &geo);
		Report("U7.a repeated query stable", a1 == a2, "bit-identical");

		// B: profile-version mutation (geo=nullptr so geometry is rebuilt from the profile)
		StarProfile p2;
		FillProfile(p2, N0, kR_km, kNb_fm3, kYq);
		StarContext c2(p2);
		const std::uint64_t v_before = c2.ProfileVersion();
		const double b1 = c2.HeatCapacityStar_Tinf(T, thermo, nullptr);
		{
			auto edit = p2.Edit();
			auto &rad = p2.RadialMutable();
			for (std::size_t i = 0; i < N0; ++i)
				rad[5].ValuesMutable()[i] = 2.0 * kNb_fm3; // double n_B
		}
		const std::uint64_t v_after = c2.ProfileVersion();
		const double b2 = c2.HeatCapacityStar_Tinf(T, thermo, nullptr);
		Report("U7.b profile version advanced", v_after != v_before,
			   "version " + std::to_string(v_before) + " -> " + std::to_string(v_after));
		Report("U7.b cache rebuilt, C_star doubled with n_B", RelErr(b2, 2.0 * b1) < 1e-9,
			   "ratio " + std::to_string(b2 / b1));

		// C: thermo-object identity
		CompOSE_Thermo thermo2((root / "linear2x").string(), opt);
		const double c_1 = ctx.HeatCapacityStar_Tinf(T, thermo, &geo);
		const double c_2 = ctx.HeatCapacityStar_Tinf(T, thermo2, &geo);
		Report("U7.c cache rebuilds on a different thermo object",
			   RelErr(c_2, 2.0 * c_1) < 1e-9, "ratio " + std::to_string(c_2 / c_1));

		// D: KNOWN HAZARD — cache key omits the GeometryCache (INV-12)
		StarProfile p3;
		FillProfile(p3, N0, kR_km, kNb_fm3, kYq, /*lambda=*/std::log(2.0)); // e^Λ = 2
		StarContext c3(p3);
		GeometryCache geo_fat(c3); // semantically different geometry

		StarProfile p4;
		FillProfile(p4, N0, kR_km, kNb_fm3, kYq);
		StarContext c4(p4);
		GeometryCache geo_flat(c4);

		const double d1 = c4.HeatCapacityStar_Tinf(T, thermo, &geo_flat);
		const double d2 = c4.HeatCapacityStar_Tinf(T, thermo, &geo_fat); // same version+thermo
		const bool stale = (d1 == d2);
		Note("U7.d GeometryCache-key hazard (INV-12)",
			 stale ? "CONFIRMED: second call with a DIFFERENT GeometryCache reused the cached "
					 "table (identical result); the key is only (profile version, thermo ptr)"
				   : "not reproduced: results differed");
		std::cout << "      flat geometry = " << d1 << ",  fat geometry (e^Lambda=2) = " << d2
				  << "  (ratio 2 expected if the GeometryCache were honored)\n";
	}

	std::cout << "\n" << (g_fail == 0 ? "Tier-A checks passed" : "Tier-A FAILURES: " + std::to_string(g_fail))
			  << "\n";
	std::cout << "Tier B (real canonical star, literature magnitude) NOT part of this executable.\n";
	fs::remove_all(root);
	return g_fail == 0 ? 0 : 1;
}
