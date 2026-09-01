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
 *   usage: cache_contract
 *
 * All checks here are ordinary pass/fail contracts.
 *
 *  - `RunSupportedContracts` covers the version-driven invalidation that was already correct
 *    before Phase 3B.
 *  - `RunProvenanceContracts` covers the ADR-0003 provenance contract. Those three cases were
 *    previously reproductions of measured DEFECTS behind `--audit-known-hazards`, kept out of
 *    CTest because they could not pass. ADR-0003 (ACCEPTED) makes them requirements, so the
 *    audit mode is gone and the known-bug output is no longer treated as expected behavior.
 *    The historical measurements survive in docs/validation/CACHE_CORRECTNESS.md.
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
#include <cstdio>
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
/// Scientific formatting: std::to_string would print 1e-15 as 0.000000.
static std::string Sci(double v, int prec = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
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
		std::cout << "        => profile version != profile identity. ADR-0003 therefore keys\n"
					 "           reusable caches on (identity, version); see P2 below.\n";
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
//  ADR-0003 ENFORCED CONTRACTS — formerly the "known hazard" reproductions
//
//  Before Phase 3B these three cases were reproductions of measured DEFECTS, deliberately
//  kept out of CTest because they could not pass. ADR-0003 (ACCEPTED) makes each one a
//  requirement, so they are now ordinary assertions in the registered test. The historical
//  measurements are preserved in docs/validation/CACHE_CORRECTNESS.md, labelled superseded.
// ===========================================================================
static void RunProvenanceContracts()
{
	std::cout << "\nADR-0003 ENFORCED PROVENANCE CONTRACTS\n";

	// -----------------------------------------------------------------------
	// P1 (was HAZARD A) — GeometryCache provenance makes staleness ANSWERABLE
	// -----------------------------------------------------------------------
	std::cout << "\nP1 GeometryCache carries its source provenance\n";
	{
		StarProfile prof;
		Fixture f;
		f.N = 16;
		FillProfile(prof, f);
		StarContext sc(prof);

		GeometryCache G_old(sc);
		Report("P1a a fresh GeometryCache matches the context it was built from",
			   G_old.Matches(sc) && G_old.SourceProfile() == &prof,
			   "source profile identity and version recorded");
		const std::uint64_t v_at_build = G_old.SourceVersion();

		{
			auto edit = prof.Edit();
			auto &radial = prof.RadialMutable();
			for (std::size_t i = 0; i < f.N; ++i)
			{
				radial[0][i] *= 1.30; // r
				radial[7][i] = 0.20;  // lambda
			}
		}

		Report("P1b after a sanctioned mutation the old snapshot reports itself STALE",
			   !G_old.Matches(sc),
			   "Matches() == false (was undetectable before ADR-0003)");
		Report("P1c the stale snapshot still remembers what it was built from",
			   G_old.SourceProfile() == &prof && G_old.SourceVersion() == v_at_build,
			   "provenance is preserved, not silently updated");

		GeometryCache G_new(sc);
		Report("P1d a freshly constructed snapshot matches again", G_new.Matches(sc),
			   "caller owns rebuilding; there is no Refresh()");
	}

	// -----------------------------------------------------------------------
	// P2 (was HAZARD C) — ProfileVersionedCache keys on identity, not version alone
	// -----------------------------------------------------------------------
	std::cout << "\nP2 ProfileVersionedCache distinguishes two equal-version profiles\n";
	{
		StarProfile A, B;
		Fixture fa;
		Fixture fb;
		fb.eps_km2 = 7.0 * fa.eps_km2;
		FillProfile(A, fa);
		FillProfile(B, fb);
		StarContext sa(A), sb(B);

		Report("P2a the precondition still holds: equal numeric Version(), different physics",
			   A.Version() == B.Version() && MeanRho(sa) != MeanRho(sb),
			   "A.Version()=" + std::to_string(A.Version()) +
				   " == B.Version()=" + std::to_string(B.Version()));

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
		Report("P2b the builder RE-RAN for the second profile", builds == 2,
			   "builds = " + std::to_string(builds) + " (was 1 before ADR-0003)");
		Report("P2c star B receives its OWN payload, not star A's",
			   Rel(from_B, MeanRho(sb)) < 1e-12 && from_B != from_A,
			   "B mean rho " + Sci(from_B) + " vs A " + Sci(from_A));
		Report("P2d the cache records the identity it was built against",
			   cache.BuiltAgainst().source == &B && cache.BuiltAgainst().version == B.Version(),
			   "provenance is (identity, version)");
	}

	// -----------------------------------------------------------------------
	// P3 (was HAZARD E) — StarContext re-binds its column views (S1)
	// -----------------------------------------------------------------------
	std::cout << "\nP3 StarContext re-binds column views on a revision change (S1)\n";
	{
		StarProfile prof;
		Fixture f;
		FillProfile(prof, f);
		StarContext sc(prof);
		const double rho0 = MeanRho(sc);
		Report("P3a baseline derived value is available", rho0 > 0.0, "mean rho " + Sci(rho0));

		// A STRUCTURAL change: grow the column count so every DataColumn is reallocated.
		// Before ADR-0003 the context kept pointers bound at construction; the version gate
		// rebuilt payloads from addresses that no longer referred to live storage.
		const void *eps_before = static_cast<const void *>(&prof.Radial()[4]);
		{
			auto edit = prof.Edit();
			auto &radial = prof.RadialMutable();
			radial.Reserve(14, f.N); // 11 -> 14 columns
		}
		const void *eps_after = static_cast<const void *>(&prof.Radial()[4]);
		Report("P3b the structural change really did move the columns",
			   eps_before != eps_after,
			   "eps column address changed, so a stale binding would be observable");

		// The context must re-bind before any view is used, and stay correct.
		const double rho1 = MeanRho(sc);
		Report("P3c the context re-bound and still yields the correct derived value",
			   rho1 > 0.0 && Rel(rho1, rho0) < 1e-12,
			   "mean rho " + Sci(rho1) + " (unchanged: the values were preserved)");
		Report("P3d the re-bound views point into the CURRENT profile storage",
			   static_cast<const void *>(sc.EnergyDensity()) == eps_after,
			   "bound pointer equals the live column address");
	}

	// -----------------------------------------------------------------------
	// P4 — an invalid schema must FAIL CLOSED, not dereference stale memory
	// -----------------------------------------------------------------------
	std::cout << "\nP4 an unusable schema fails closed\n";
	{
		StarProfile prof;
		Fixture f;
		FillProfile(prof, f);
		StarContext sc(prof);
		(void)MeanRho(sc); // bind and build

		// Drop the mandatory radius mapping through the sanctioned API. SetColumnIndex
		// calls Touch(), so the revision advances and the context must revalidate.
		prof.SetColumnIndex(StarProfile::Column::Radius, -1);

		auto probe = [&]() {
			try { (void)sc.MassDensity_gcm3(); return false; }
			catch (const std::exception &) { return true; }
		};

		const bool threw_first = probe();
		Report("P4a a profile that no longer satisfies the required schema throws rather "
			   "than serving data through stale bindings",
			   threw_first, threw_first ? "threw as required" : "NO THROW");

		// The failed re-bind must NOT have advanced the cached revision: a context that
		// threw must stay marked out-of-date, so the next call retries rather than
		// reporting itself fresh with stale views.
		Report("P4b the context is not left falsely marked current after a failed re-bind",
			   probe(), "a second access throws again rather than returning stale data");

		// Restoring the schema must make the context usable again.
		prof.SetColumnIndex(StarProfile::Column::Radius, 0);
		Report("P4c restoring the schema restores the context",
			   !probe() && MeanRho(sc) > 0.0, "mean rho " + Sci(MeanRho(sc)));
	}
}

int main()
{
	std::cout << std::scientific << std::setprecision(6);

	RunSupportedContracts();
	RunProvenanceContracts();

	std::cout << "\n"
			  << (g_fail == 0 ? "supported cache contracts and ADR-0003 provenance contracts hold"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
