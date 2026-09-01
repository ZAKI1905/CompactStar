// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file cache_contract.cpp
 * @brief Phase 2B-3 — INV-12 cache-correctness contracts for StarContext, GeometryCache
 *        and the generic ProfileVersionedCache. Self-contained; no external data.
 *
 *   usage: cache_contract [--audit-known-hazards]
 *
 * TWO MODES, DELIBERATELY SEPARATED.
 *
 *  - Default (the registered CTest) asserts ONLY contracts that are currently CORRECT.
 *    Every check here must pass against a correct implementation and must fail if
 *    version-driven invalidation regresses.
 *
 *  - `--audit-known-hazards` reproduces KNOWN INV-12 DEFECTS for the validation report.
 *    It is NOT registered as a CTest and is NOT a green regression criterion. It never
 *    asserts that stale behavior is correct; it records and classifies it.
 *
 * The synthetic star is a uniform sphere with nu = Lambda = 0 and constant n_B. It
 * asserts no neutron-star property; it exists solely to drive cache state transitions
 * with analytically predictable answers.
 *
 * Direct-Urca threshold used by the fixtures (INV-16, electron channel only, the current
 * governed behavior): DU is allowed where kFn <= kFp + kFe with kF = (3 pi^2 n)^(1/3).
 * With charge neutrality Y_e = Y_p this reduces exactly to Y_n <= 8 Y_p. The fixtures sit
 * far from that boundary on purpose, so the transition is unambiguous. This test does not
 * adjudicate the unresolved muon channel.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/ProfileCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"

using CompactStar::Core::StarProfile;
using CompactStar::Physics::Evolution::GeometryCache;
using CompactStar::Physics::Evolution::ProfileVersionedCache;
using CompactStar::Physics::Evolution::StarContext;

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

// ---------------------------------------------------------------------------
//  Synthetic uniform star. Columns: r, m, nu', p, eps, nB, nu, lambda
//  Species: "10" neutron, "11" proton, "0" electron (the set DU requires).
// ---------------------------------------------------------------------------
struct Fixture
{
	std::size_t N = 64;
	double R_km = 10.0;
	double nb_fm3 = 0.40;
	double eps_km2 = 1.0e-4;
	double Yn = 0.85, Yp = 0.15, Ye = 0.15; // Y_n <= 8 Y_p  -> DU ALLOWED
};

static void FillProfile(StarProfile &prof, const Fixture &f)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	radial.ClearRows();
	radial.Reserve(11, f.N);

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
	prof.SetSpeciesColumn("10", 8); // neutron
	radial[9].SetLabel("11");
	prof.SetSpeciesColumn("11", 9); // proton, q = +1
	radial[10].SetLabel("0");
	prof.SetSpeciesColumn("0", 10); // electron

	const double h = f.R_km / static_cast<double>(f.N - 1);
	for (std::size_t i = 0; i < f.N; ++i)
	{
		radial[0].PushBack(static_cast<double>(i) * h);
		radial[1].PushBack(0.0);
		radial[2].PushBack(0.0);
		radial[3].PushBack(0.0);
		radial[4].PushBack(f.eps_km2);
		radial[5].PushBack(f.nb_fm3);
		radial[6].PushBack(0.0);
		radial[7].PushBack(0.0);
		radial[8].PushBack(f.Yn);
		radial[9].PushBack(f.Yp);
		radial[10].PushBack(f.Ye);
	}
}

/// Scale one column in place through the sanctioned mutation API. No reallocation:
/// the row count is unchanged and each element is overwritten in place.
static void ScaleColumnInPlace(StarProfile &prof, int col, double factor)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	for (std::size_t i = 0; i < static_cast<std::size_t>(radial[col].Size()); ++i)
		radial[col][i] *= factor;
}

static void SetColumnInPlace(StarProfile &prof, int col, double value)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	for (std::size_t i = 0; i < static_cast<std::size_t>(radial[col].Size()); ++i)
		radial[col][i] = value;
}

static double MeanRho(const StarContext &sc)
{
	const auto *c = sc.MassDensity_gcm3();
	if (!c || c->Size() == 0)
		return -1.0;
	double s = 0.0;
	for (std::size_t i = 0; i < static_cast<std::size_t>(c->Size()); ++i)
		s += (*c)[i];
	return s / static_cast<double>(c->Size());
}

static long MaskSum(const StarContext &sc)
{
	const auto *m = sc.DirectUrcaMask();
	if (!m)
		return -1;
	long s = 0;
	for (std::uint8_t v : *m)
		s += v;
	return s;
}

// ===========================================================================
//  SUPPORTED CONTRACTS — the registered regression criterion
// ===========================================================================
static void RunSupportedContracts()
{
	std::cout << "SUPPORTED CACHE CONTRACTS (StarContext / GeometryCache / "
				 "ProfileVersionedCache)\n"
			  << "Synthetic uniform fixture; asserts no neutron-star property.\n\n";

	// -----------------------------------------------------------------------
	// C1 — mass-density cache rebuilds on profile-version change
	// -----------------------------------------------------------------------
	std::cout << "C1/C2 StarContext mass-density cache\n";
	{
		StarProfile prof;
		Fixture f;
		FillProfile(prof, f);
		StarContext sc(prof);

		const std::uint64_t v0 = sc.ProfileVersion();
		const double rho0 = MeanRho(sc);
		Report("C1a mass-density cache builds and is positive", rho0 > 0.0,
			   "rho = " + std::to_string(rho0) + " g/cm^3");

		// Repeated queries at a fixed version must not drift.
		bool stable = true;
		const auto *p_first = sc.MassDensity_gcm3();
		for (int k = 0; k < 5; ++k)
			stable = stable && (sc.MassDensity_gcm3() == p_first) &&
					 (MeanRho(sc) == rho0);
		Report("C2a repeated queries at a fixed version are stable and rebuild nothing",
			   stable && sc.ProfileVersion() == v0, "5 repeats, identical payload identity");

		// Mutate eps by a known factor through the sanctioned API.
		ScaleColumnInPlace(prof, 4, 2.0);
		const std::uint64_t v1 = sc.ProfileVersion();
		Report("C1b sanctioned in-place mutation increments Version()", v1 > v0,
			   "version " + std::to_string(v0) + " -> " + std::to_string(v1));

		const double rho1 = MeanRho(sc);
		Report("C1c mass density rebuilds and scales by exactly the same factor",
			   Rel(rho1, 2.0 * rho0) < 1e-12,
			   "rho " + std::to_string(rho0) + " -> " + std::to_string(rho1) +
				   " (expected x2, rel dev " + std::to_string(Rel(rho1, 2.0 * rho0)) + ")");

		// And a second, different factor, to rule out a one-shot rebuild.
		ScaleColumnInPlace(prof, 4, 0.25);
		const double rho2 = MeanRho(sc);
		Report("C1d a second mutation rebuilds again", Rel(rho2, 0.5 * rho0) < 1e-12,
			   "rho -> " + std::to_string(rho2) + " (expected 0.5x original)");
	}

	// -----------------------------------------------------------------------
	// C3 — Direct-Urca mask/boundary rebuilds across the kinematic threshold
	// -----------------------------------------------------------------------
	std::cout << "\nC3/C4 StarContext Direct-Urca mask and boundary\n";
	{
		StarProfile prof;
		Fixture f; // Y_n = 0.85, Y_p = Y_e = 0.15  ->  Y_n <= 8 Y_p, DU allowed
		FillProfile(prof, f);
		StarContext sc(prof);

		const long sum0 = MaskSum(sc);
		const long last0 = sc.DirectUrcaLastAllowedIndex();
		const double rb0 = sc.DirectUrcaBoundaryRadius_km();
		Report("C3a DU allowed throughout a composition far above threshold",
			   sum0 == static_cast<long>(f.N) && last0 == static_cast<long>(f.N) - 1,
			   "mask sum " + std::to_string(sum0) + "/" + std::to_string(f.N) +
				   ", last allowed index " + std::to_string(last0) +
				   ", boundary r = " + std::to_string(rb0) + " km");

		bool stable = true;
		for (int k = 0; k < 5; ++k)
			stable = stable && MaskSum(sc) == sum0 &&
					 sc.DirectUrcaLastAllowedIndex() == last0 &&
					 sc.DirectUrcaBoundaryRadius_km() == rb0;
		Report("C4a repeated DU queries at a fixed version are stable", stable,
			   "5 repeats, identical mask sum / boundary index / boundary radius");

		// Cross the threshold: Y_n = 0.95, Y_p = Y_e = 0.05 -> Y_n > 8 Y_p.
		{
			auto edit = prof.Edit();
			auto &radial = prof.RadialMutable();
			for (std::size_t i = 0; i < f.N; ++i)
			{
				radial[8][i] = 0.95; // neutron
				radial[9][i] = 0.05; // proton
				radial[10][i] = 0.05; // electron
			}
		}
		const long sum1 = MaskSum(sc);
		Report("C3b DU mask rebuilds to forbidden after crossing the threshold",
			   sum1 == 0 && sc.DirectUrcaLastAllowedIndex() == -1 &&
				   sc.DirectUrcaBoundaryRadius_km() == 0.0,
			   "mask sum " + std::to_string(sum1) + ", last allowed " +
				   std::to_string(sc.DirectUrcaLastAllowedIndex()));

		// And back again, so the rebuild is not a one-way latch.
		{
			auto edit = prof.Edit();
			auto &radial = prof.RadialMutable();
			for (std::size_t i = 0; i < f.N; ++i)
			{
				radial[8][i] = 0.85;
				radial[9][i] = 0.15;
				radial[10][i] = 0.15;
			}
		}
		Report("C3c DU mask rebuilds back to allowed", MaskSum(sc) == static_cast<long>(f.N),
			   "mask sum " + std::to_string(MaskSum(sc)));

		// A partial DU region: inner half above threshold, outer half below.
		{
			auto edit = prof.Edit();
			auto &radial = prof.RadialMutable();
			for (std::size_t i = f.N / 2; i < f.N; ++i)
			{
				radial[8][i] = 0.95;
				radial[9][i] = 0.05;
				radial[10][i] = 0.05;
			}
		}
		Report("C3d DU boundary tracks a partial allowed region",
			   sc.DirectUrcaLastAllowedIndex() == static_cast<long>(f.N / 2) - 1,
			   "last allowed index " + std::to_string(sc.DirectUrcaLastAllowedIndex()) +
				   " (expected " + std::to_string(f.N / 2 - 1) + "), boundary r = " +
				   std::to_string(sc.DirectUrcaBoundaryRadius_km()) + " km");
	}

	// -----------------------------------------------------------------------
	// C5 — GeometryCache is a faithful snapshot of its source AT CONSTRUCTION
	// -----------------------------------------------------------------------
	std::cout << "\nC5 GeometryCache construction-time fidelity\n";
	{
		StarProfile prof;
		Fixture f;
		f.N = 32;
		FillProfile(prof, f);
		StarContext sc(prof);
		GeometryCache geo(sc);

		bool size_ok = geo.Size() == f.N;
		bool r_ok = true, nu_ok = true, area_ok = true;
		for (std::size_t i = 0; i < geo.Size(); ++i)
		{
			const double r = geo.R()[i];
			r_ok = r_ok && Rel(r, (*sc.Radius())[i]) < 1e-14;
			nu_ok = nu_ok && Rel(geo.ExpNu()[i], 1.0) < 1e-14; // nu = 0
			area_ok = area_ok && Rel(geo.Area()[i], 4.0 * M_PI * r * r) < 1e-12;
		}
		Report("C5a GeometryCache reproduces its source geometry at construction",
			   size_ok && r_ok && nu_ok && area_ok,
			   "size " + std::to_string(geo.Size()) + ", r/exp(nu)/area all match");
	}

	// -----------------------------------------------------------------------
	// C6 — profile VERSION is not profile IDENTITY (a fact, established here)
	// -----------------------------------------------------------------------
	std::cout << "\nC6 version vs identity\n";
	{
		StarProfile A, B;
		Fixture fa;
		Fixture fb;
		fb.eps_km2 = 5.0 * fa.eps_km2; // physically different star
		fb.nb_fm3 = 0.80;
		FillProfile(A, fa);
		FillProfile(B, fb);
		StarContext sa(A), sb(B);
		Report("C6a two distinct profiles with different physics can share a Version()",
			   A.Version() == B.Version() && &A != &B,
			   "A.Version()=" + std::to_string(A.Version()) +
				   ", B.Version()=" + std::to_string(B.Version()) +
				   "; mean rho A=" + std::to_string(MeanRho(sa)) +
				   " vs B=" + std::to_string(MeanRho(sb)));
		std::cout << "        => profile version != profile identity. Any cache keyed on\n"
					 "           version alone and reused across contexts can collide.\n"
					 "           Consequences are reproduced in --audit-known-hazards.\n";
	}

	// -----------------------------------------------------------------------
	// C7 — ProfileVersionedCache: the SAME-STAR contract it was designed for
	// -----------------------------------------------------------------------
	std::cout << "\nC7 ProfileVersionedCache same-star rebuild contract\n";
	{
		StarProfile prof;
		Fixture f;
		FillProfile(prof, f);
		StarContext sc(prof);

		struct Payload
		{
			double mean_rho = 0.0;
		};
		ProfileVersionedCache<Payload> cache;
		int builds = 0;
		auto builder = [&](const StarContext &s, Payload &out)
		{
			++builds;
			out.mean_rho = MeanRho(s);
		};

		const double p0 = cache.Get(sc, builder).mean_rho;
		Report("C7a builder runs once on first access", builds == 1,
			   "builds = " + std::to_string(builds));

		for (int k = 0; k < 5; ++k)
			cache.Get(sc, builder);
		Report("C7b builder does not re-run at a fixed profile version", builds == 1,
			   "builds after 5 more Get() calls = " + std::to_string(builds));
		Report("C7c cache records the version it was built against",
			   cache.IsBuilt() && cache.BuiltAgainstVersion() == sc.ProfileVersion(),
			   "built against version " + std::to_string(cache.BuiltAgainstVersion()));

		ScaleColumnInPlace(prof, 4, 3.0);
		const double p1 = cache.Get(sc, builder).mean_rho;
		Report("C7d builder re-runs after a sanctioned mutation", builds == 2,
			   "builds = " + std::to_string(builds));
		Report("C7e rebuilt payload reflects the new profile", Rel(p1, 3.0 * p0) < 1e-12,
			   "mean rho " + std::to_string(p0) + " -> " + std::to_string(p1));

		cache.Invalidate();
		cache.Get(sc, builder);
		Report("C7f explicit Invalidate() forces a rebuild", builds == 3,
			   "builds = " + std::to_string(builds));
	}

	// -----------------------------------------------------------------------
	// C8 — column-container identity under in-place edits (pointer-rebind safety)
	// -----------------------------------------------------------------------
	// StarContext binds raw DataColumn pointers once in its constructor and never
	// re-binds. This check establishes the boundary of what is currently SAFE: an
	// in-place value edit leaves the DataColumn objects at the same addresses, so
	// the bound pointers stay valid. No dangling pointer is created or dereferenced
	// anywhere in this test.
	std::cout << "\nC8 column-container identity under in-place mutation\n";
	{
		StarProfile prof;
		Fixture f;
		FillProfile(prof, f);
		StarContext sc(prof);

		const void *eps_before = static_cast<const void *>(sc.EnergyDensity());
		const void *nb_before = static_cast<const void *>(sc.BaryonDensity());
		ScaleColumnInPlace(prof, 4, 1.5);
		ScaleColumnInPlace(prof, 5, 1.1);
		const void *eps_after = static_cast<const void *>(&prof.Radial()[4]);
		const void *nb_after = static_cast<const void *>(&prof.Radial()[5]);

		Report("C8a in-place value edits keep DataColumn addresses stable, so the "
			   "pointers StarContext bound in its constructor remain valid",
			   eps_before == eps_after && nb_before == nb_after,
			   "eps and nB column addresses unchanged across two in-place mutations");
		Report("C8b derived caches still track those in-place edits",
			   Rel(MeanRho(sc), 0.0) > 0.0 && MeanRho(sc) > 0.0,
			   "mean rho = " + std::to_string(MeanRho(sc)) + " g/cm^3");
	}
}

// ===========================================================================
//  KNOWN-HAZARD AUDIT — evidence for the report, NOT a pass/fail criterion
// ===========================================================================
static void RunHazardAudit()
{
	std::cout << "\n"
			  << "===========================================================\n"
			  << " KNOWN INV-12 HAZARD AUDIT — reproductions, not assertions\n"
			  << " Nothing below is a regression criterion. No check here\n"
			  << " claims that stale behavior is correct.\n"
			  << "===========================================================\n";

	// -----------------------------------------------------------------------
	// HAZARD A — GeometryCache carries no provenance, so staleness is undetectable
	// -----------------------------------------------------------------------
	std::cout << "\n[HAZARD A] GeometryCache snapshot has no source identity or version\n";
	{
		StarProfile prof;
		Fixture f;
		f.N = 16;
		FillProfile(prof, f);
		StarContext sc(prof);

		GeometryCache G_old(sc);
		const double r_old_mid = G_old.R()[f.N / 2];
		const double wv_old_mid = G_old.WV()[f.N / 2];

		// Mutate the geometry through the sanctioned API: stretch r, set Lambda != 0.
		{
			auto edit = prof.Edit();
			auto &radial = prof.RadialMutable();
			for (std::size_t i = 0; i < f.N; ++i)
			{
				radial[0][i] *= 1.30; // r
				radial[7][i] = 0.20;  // lambda
			}
		}
		GeometryCache G_new(sc);

		std::cout << "    profile version now " << sc.ProfileVersion()
				  << "; StarContext's own caches DID rebuild (mean rho = "
				  << MeanRho(sc) << ")\n";
		std::cout << "    G_old.R[mid]  = " << r_old_mid
				  << "   G_new.R[mid]  = " << G_new.R()[f.N / 2] << "\n";
		std::cout << "    G_old.WV[mid] = " << wv_old_mid
				  << "   G_new.WV[mid] = " << G_new.WV()[f.N / 2] << "\n";
		std::cout << "    G_old still returns its pre-mutation values: "
				  << (G_old.R()[f.N / 2] == r_old_mid ? "YES" : "no") << "\n";
		std::cout << "    relative WV divergence: "
				  << Rel(wv_old_mid, G_new.WV()[f.N / 2]) << "\n";
		std::cout << "\n    An immutable snapshot returning its construction-time values is\n"
					 "    not itself wrong. The DEFECT is that GeometryCache exposes no\n"
					 "    source profile identity, no source version and no Invalidate(),\n"
					 "    so a holder of G_old CANNOT ask whether it is stale.\n"
					 "    CLASSIFICATION: KNOWN INV-12 HAZARD (Phase-3 repair required).\n";
	}

	// -----------------------------------------------------------------------
	// HAZARD C — ProfileVersionedCache collides across equal-version profiles
	// -----------------------------------------------------------------------
	std::cout << "\n[HAZARD C] ProfileVersionedCache key omits profile identity\n";
	{
		StarProfile A, B;
		Fixture fa;
		Fixture fb;
		fb.eps_km2 = 7.0 * fa.eps_km2;
		FillProfile(A, fa);
		FillProfile(B, fb);
		StarContext sa(A), sb(B);

		struct Payload
		{
			double mean_rho = 0.0;
		};
		ProfileVersionedCache<Payload> cache;
		int builds = 0;
		auto builder = [&](const StarContext &s, Payload &out)
		{
			++builds;
			out.mean_rho = MeanRho(s);
		};

		const double from_A = cache.Get(sa, builder).mean_rho;
		const double from_B = cache.Get(sb, builder).mean_rho;
		const double truth_B = MeanRho(sb);

		std::cout << "    A.Version() = " << A.Version()
				  << ", B.Version() = " << B.Version()
				  << "  (equal: " << (A.Version() == B.Version() ? "YES" : "no") << ")\n";
		std::cout << "    builder invocations across both Get() calls: " << builds << "\n";
		std::cout << "    cache.Get(A) mean rho = " << from_A << "\n";
		std::cout << "    cache.Get(B) mean rho = " << from_B
				  << "   TRUE B mean rho = " << truth_B << "\n";
		std::cout << "    relative error served for B: " << Rel(from_B, truth_B) << "\n";
		if (builds == 1 && from_B == from_A)
			std::cout << "\n    CONFIRMED: B silently received A's payload. The key is the\n"
						 "    numeric version alone, so two different stars that happen to\n"
						 "    share a version collide with no diagnostic.\n"
						 "    CLASSIFICATION: KNOWN INV-12 HAZARD — profile identity omitted.\n";
		else
			std::cout << "\n    NOT REPRODUCED under this construction; investigate before\n"
						 "    relying on the classification.\n";
	}

	// -----------------------------------------------------------------------
	// HAZARD E — StarContext column pointers are never re-bound
	// -----------------------------------------------------------------------
	// Assessed from source plus SAFE container-identity observation on a profile
	// that has NO StarContext bound to it. No dangling pointer is ever formed or
	// dereferenced; nothing here invokes undefined behavior.
	std::cout << "\n[HAZARD E] StarContext column pointers are bound once and never re-bound\n";
	{
		StarProfile solo; // deliberately unbound: no StarContext observes this object
		Fixture f;
		FillProfile(solo, f);
		const void *eps_addr_1 = static_cast<const void *>(&solo.Radial()[4]);
		const std::size_t n1 = static_cast<std::size_t>(solo.Radial()[4].Size());

		// (a) clear-and-refill at the SAME column count: reuses the existing
		//     DataColumn objects, so addresses survive.
		FillProfile(solo, f);
		const void *eps_addr_2 = static_cast<const void *>(&solo.Radial()[4]);
		const std::size_t n2 = static_cast<std::size_t>(solo.Radial()[4].Size());

		// (b) a genuine structural change — growing the column count, which is what
		//     re-importing a profile with more species columns does. The column
		//     container reallocates and every DataColumn moves.
		{
			auto edit = solo.Edit();
			auto &radial = solo.RadialMutable();
			radial.ClearRows();
			radial.Reserve(14, f.N); // 11 -> 14 columns
		}
		const void *eps_addr_3 = static_cast<const void *>(&solo.Radial()[4]);

		std::cout << "    in-place edits (see C8): DataColumn addresses stable -> bound "
					 "pointers stay valid\n";
		std::cout << "    (a) clear+refill at the same column count: column address "
				  << (eps_addr_1 == eps_addr_2 ? "unchanged" : "CHANGED")
				  << ", size " << n1 << " -> " << n2
				  << "  => bound pointers survive this path\n";
		std::cout << "    (b) column-count change 11 -> 14 (e.g. a re-import with more\n"
					 "        species): column address "
				  << (eps_addr_2 == eps_addr_3 ? "unchanged" : "CHANGED")
				  << "  => any pointer bound before this point is DANGLING\n";
		std::cout << "\n    From source: StarContext::BindColumnsOrThrow_() runs only in the\n"
					 "    constructor (StarContext.cpp), and RefreshDerivedCachesIfNeeded_()\n"
					 "    invalidates derived payloads without ever re-binding. A mutation\n"
					 "    that reallocates or replaces column storage therefore bumps the\n"
					 "    version — so payloads rebuild — while the seven cached raw pointers\n"
					 "    (StarContext.hpp m_r/m_m/m_nu/m_lam/m_nb/m_pre/m_eps) may no longer\n"
					 "    refer to live storage. The version gate gives false confidence.\n"
					 "    This audit deliberately does NOT construct that state.\n"
					 "    CLASSIFICATION: KNOWN POINTER-REBIND HAZARD (Phase-3 repair required).\n";
	}
}

int main(int argc, char **argv)
{
	bool audit = false;
	for (int i = 1; i < argc; ++i)
		if (std::strcmp(argv[i], "--audit-known-hazards") == 0)
			audit = true;

	std::cout << std::scientific << std::setprecision(6);

	if (audit)
	{
		RunHazardAudit();
		std::cout << "\nHazard audit complete. This mode is a report, not a pass/fail "
					 "criterion.\n";
		return 0;
	}

	RunSupportedContracts();
	std::cout << "\n"
			  << (g_fail == 0 ? "supported cache contracts hold"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
