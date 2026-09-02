// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file geometry_cache_measure_contract.cpp
 * @brief Phase 3D — `GeometryCache` conformance to the ADR-0004 proper-volume contract.
 *
 * Self-contained: a synthetic uniform-density star built directly as a `StarProfile`. No EOS
 * assets and no TOV solve.
 *
 * ADR-0004 keeps `GeometryCache` as the canonical **cached** owner (role B) while the
 * mathematics moves to `CompactStar/Geometry.hpp` (role A). That split is only safe if three
 * things are true, and this test asserts each of them:
 *
 *   G1. The two Λ routes agree. `GeometryCache` uses the profile's `Λ` column when present and
 *       otherwise derives Λ from `(m, r)`. ADR-0004 §9 measured these bit-identical on all four
 *       authenticated stars — but that was a property of how `NStar` happens to build profiles,
 *       **not** a guarantee. A profile read from a file could carry a Λ produced elsewhere. This
 *       pins the derived route to the canonical formula so the two cannot silently diverge.
 *
 *   G2/G4. The cached `WV` is the same measure the primitive defines — `4π r² e^Λ`, node by
 *       node. If `GeometryCache` ever grew a second definition, this fires.
 *
 *   G3. The redshifted measures are **composed from the one `WV` array**, not independently
 *       re-derived (ADR-0004 §12). This is the D6 detector's target: composing `WVExpNu` from a
 *       second, separately re-derived measure must fail here even when it agrees to many digits.
 *
 * Bitwise equality is required throughout, not a tolerance. Every quantity compared is produced
 * from the same inputs by the same operations; anything less than bitwise means a second
 * definition has appeared somewhere.
 *
 * ADR-0003 provenance is NOT re-tested here — `cache_contract` owns that — beyond confirming
 * this change did not disturb it.
 */

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Geometry.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"

using CompactStar::Core::StarProfile;
using CompactStar::Physics::Evolution::GeometryCache;
using CompactStar::Physics::Evolution::StarContext;
namespace Geo = CompactStar::Geometry;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}
static std::string Sci(double v, int prec = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
}

// ---------------------------------------------------------------------------
//  Synthetic uniform-density star: m(r) = M (r/R)^3, so 2m/r rises smoothly to
//  2M/R = 0.40 at the surface — comparable to the 0.481 maximum measured on the
//  authenticated CMF stars (ADR-0004 §11), and safely inside the domain.
//  nu(r) varies so that exp(nu) != 1 and the composition checks have real content.
// ---------------------------------------------------------------------------
namespace
{
constexpr std::size_t kN = 96;
constexpr double kR = 12.0; // km
constexpr double kM = 2.4;  // km  => 2M/R = 0.40

/// Independent oracle for the stored Lambda column. Deliberately does NOT call the primitive:
/// it is the expression the profile-writing path has always used, written out here.
double StoredLambdaOracle(double r, double m)
{
	if (r == 0.0)
		return 0.0;
	return -0.5 * std::log(1.0 - 2.0 * m / r);
}

void FillProfile(StarProfile &prof, bool write_lambda_column)
{
	auto edit = prof.Edit();
	auto &radial = prof.RadialMutable();
	radial.ClearRows();
	radial.Reserve(8, kN);

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

	// The derived-Lambda route is forced by leaving the MetricLambda mapping absent, which is
	// how a profile without that column presents itself to StarContext::Lambda().
	prof.SetColumnIndex(StarProfile::Column::MetricLambda, write_lambda_column ? 7 : -1);

	const double h = kR / static_cast<double>(kN - 1);
	for (std::size_t i = 0; i < kN; ++i)
	{
		const double r = static_cast<double>(i) * h;
		const double x = r / kR;
		const double m = kM * x * x * x;

		radial[0].PushBack(r);
		radial[1].PushBack(m);
		radial[2].PushBack(0.0);
		radial[3].PushBack(0.0);
		radial[4].PushBack(1.0e-4);
		radial[5].PushBack(0.40);
		radial[6].PushBack(-0.25 + 0.10 * x * x); // non-trivial nu
		radial[7].PushBack(StoredLambdaOracle(r, m));
	}
}
} // namespace

int main()
{
	std::cout << "\n=== ADR-0004 GeometryCache measure conformance ===\n";

	StarProfile prof_stored, prof_derived;
	FillProfile(prof_stored, true);
	FillProfile(prof_derived, false);

	StarContext sc_stored(prof_stored), sc_derived(prof_derived);
	GeometryCache g_stored(sc_stored), g_derived(sc_derived);

	// -----------------------------------------------------------------------
	// G1 — the two Lambda routes must produce identical cached arrays.
	// -----------------------------------------------------------------------
	std::cout << "\nG1 stored-Lambda route vs derived-Lambda route\n";
	{
		Report("G1a both routes produced a cache of the same size",
			   g_stored.Size() == kN && g_derived.Size() == kN,
			   "stored " + std::to_string(g_stored.Size()) + ", derived " +
				   std::to_string(g_derived.Size()));

		struct Arr
		{
			const char *name;
			const Zaki::Vector::DataColumn *a;
			const Zaki::Vector::DataColumn *b;
		};
		const Arr arrays[] = {
			{"ExpLambda", &g_stored.ExpLambda(), &g_derived.ExpLambda()},
			{"WV", &g_stored.WV(), &g_derived.WV()},
			{"WVExpNu", &g_stored.WVExpNu(), &g_derived.WVExpNu()},
			{"WVExp2Nu", &g_stored.WVExp2Nu(), &g_derived.WVExp2Nu()},
		};

		for (const auto &A : arrays)
		{
			std::size_t bitwise = 0;
			double worst = 0.0;
			for (std::size_t i = 0; i < kN; ++i)
			{
				const double x = (*A.a)[i], y = (*A.b)[i];
				if (x == y)
					++bitwise;
				else
					worst = std::max(worst, std::fabs(x - y) / std::fabs(x));
			}
			Report(std::string("G1 ") + A.name + " bitwise on both routes", bitwise == kN,
				   std::to_string(bitwise) + "/" + std::to_string(kN) +
					   (bitwise == kN ? " bitwise" : " bitwise, worst rel " + Sci(worst)));
		}
	}

	// -----------------------------------------------------------------------
	// G2 / G4 — the cached measure IS the primitive's measure.
	// -----------------------------------------------------------------------
	std::cout << "\nG2 cached WV vs the canonical primitive\n";
	{
		std::size_t bitwise_prim = 0, bitwise_comp = 0;
		double worst_prim = 0.0;
		for (std::size_t i = 0; i < kN; ++i)
		{
			const double r = g_stored.R()[i];
			const double m = (*sc_stored.Mass())[i];
			const double want = Geo::ProperVolumeWeight(r, m);
			const double got = g_stored.WV()[i];
			if (got == want)
				++bitwise_prim;
			else if (std::fabs(want) > 0.0)
				worst_prim = std::max(worst_prim, std::fabs(got - want) / std::fabs(want));

			// The composition the cache itself performs: (4*pi*r^2) * exp(Lambda).
			if (got == g_stored.Area()[i] * g_stored.ExpLambda()[i])
				++bitwise_comp;
		}
		Report("G2a cached WV equals Geometry::ProperVolumeWeight node by node",
			   bitwise_prim == kN,
			   std::to_string(bitwise_prim) + "/" + std::to_string(kN) + " bitwise" +
				   (bitwise_prim == kN ? "" : ", worst rel " + Sci(worst_prim)));
		Report("G2b cached WV equals Area * ExpLambda", bitwise_comp == kN,
			   std::to_string(bitwise_comp) + "/" + std::to_string(kN) + " bitwise");
	}

	// -----------------------------------------------------------------------
	// G3 — redshifted measures are COMPOSED from the one WV array (ADR-0004 §12).
	//      This is detector D6's target.
	// -----------------------------------------------------------------------
	std::cout << "\nG3 redshifted measures are composed from the single WV array\n";
	{
		std::size_t n1 = 0, n2 = 0;
		double worst1 = 0.0, worst2 = 0.0;
		for (std::size_t i = 0; i < kN; ++i)
		{
			const double wv = g_stored.WV()[i];
			const double e1 = g_stored.ExpNu()[i];
			const double e2 = g_stored.Exp2Nu()[i];

			const double p1 = g_stored.WVExpNu()[i], w1 = wv * e1;
			const double p2 = g_stored.WVExp2Nu()[i], w2 = wv * e2;
			if (p1 == w1)
				++n1;
			else
				worst1 = std::max(worst1, std::fabs(p1 - w1) / std::fabs(w1));
			if (p2 == w2)
				++n2;
			else
				worst2 = std::max(worst2, std::fabs(p2 - w2) / std::fabs(w2));
		}
		Report("G3a WVExpNu == WV * ExpNu bitwise", n1 == kN,
			   std::to_string(n1) + "/" + std::to_string(kN) + " bitwise" +
				   (n1 == kN ? "" : ", worst rel " + Sci(worst1)));
		Report("G3b WVExp2Nu == WV * Exp2Nu bitwise", n2 == kN,
			   std::to_string(n2) + "/" + std::to_string(kN) + " bitwise" +
				   (n2 == kN ? "" : ", worst rel " + Sci(worst2)));

		// exp(nu) must actually vary, or G3 would be vacuous.
		Report("G3c the fixture's exp(nu) is non-trivial",
			   std::fabs(g_stored.ExpNu()[0] - g_stored.ExpNu()[kN - 1]) > 1e-3,
			   "exp(nu): " + Sci(g_stored.ExpNu()[0]) + " -> " + Sci(g_stored.ExpNu()[kN - 1]));
	}

	// -----------------------------------------------------------------------
	// G5 — ADR-0003 provenance is undisturbed by the ADR-0004 delegation.
	// -----------------------------------------------------------------------
	std::cout << "\nG5 ADR-0003 provenance unaffected\n";
	{
		Report("G5a cache carries the provenance of the context it was built from",
			   g_stored.Provenance() == sc_stored.Provenance(), "provenance matches");
		Report("G5b caches over distinct profiles carry distinct provenance",
			   !(g_stored.Provenance() == g_derived.Provenance()),
			   "distinct source profiles are distinguishable");
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "GeometryCache conforms to the ADR-0004 measure contract"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
