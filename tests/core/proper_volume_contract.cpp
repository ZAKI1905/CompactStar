// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file proper_volume_contract.cpp
 * @brief Contract test for the canonical metric/proper-volume primitive, `CompactStar/Geometry.hpp`.
 *
 * Self-contained: no EOS assets, no TOV solve, no profile. It exercises the primitive directly,
 * which is the point — the mathematical owner must be verifiable without any of the machinery
 * that consumes it.
 *
 * WHAT IS TESTED — the accepted contract of **ADR-0004** (2026-09-01):
 *
 *   A. the exact regular-centre limit at `(r, m) = (0, 0)`;
 *   B. fail-closed behavior on every input the contract declares invalid;
 *   C. the metric factor against two independent algebraic forms on the physical domain;
 *   D. the proper-volume composition `w_V = 4π r² e^Λ`.
 *
 * ORACLE INDEPENDENCE. Group C first pins `f` itself: the primitive's `MetricDenominator` must
 * be bitwise the `1 - 2m/r` subtraction, written out independently here. It then compares the
 * production `ExpLambda` against `1/sqrt(f)` and `exp(-0.5*log(f))` — the three algebraic forms
 * of the SAME `f`, which is the comparison ADR-0004 §6.2 measured. Comparing instead against a
 * separately constructed `1 - c` would measure the fixture's own rounding (`2*(0.5*c*R)/R` is not
 * bitwise `c`), not the primitive. The `2` ULP allowance is the bound §6.2 derived from the
 * operation counts *before* this test existed; it is not tuned to observed output.
 *
 * WHAT IS NOT TESTED HERE. Baryon number (see `baryon_number_cmf`), the cached `GeometryCache`
 * arrays (see `geometry_cache_measure_contract`), and anything requiring a real EOS.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "CompactStar/Geometry.hpp"

namespace Geo = CompactStar::Geometry;

static int g_fail = 0;

/// Scientific formatting: std::to_string prints -0.0 and 1e-15 alike as "0.000000".
static std::string Sci(double v, int prec = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
}

static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}

/// Distance in representable doubles between `a` and `b`. Both must be finite and same-signed.
static long UlpDistance(double a, double b)
{
	if (a == b)
		return 0;
	std::int64_t ia = 0, ib = 0;
	std::memcpy(&ia, &a, sizeof ia);
	std::memcpy(&ib, &b, sizeof ib);
	if ((ia < 0) != (ib < 0))
		return std::numeric_limits<long>::max();
	return static_cast<long>(ia > ib ? ia - ib : ib - ia);
}

/// Runs `fn`, returning true iff it threw `std::runtime_error` — the project's fail-closed
/// convention. A silently-returned value, or any other exception type, is a contract violation.
template <typename F>
static bool FailsClosed(F &&fn, std::string *what = nullptr)
{
	try
	{
		const double v = fn();
		if (what)
			*what = "returned " + std::to_string(v) + " instead of throwing";
		return false;
	}
	catch (const std::runtime_error &e)
	{
		if (what)
			*what = std::string("threw: ") + e.what();
		return true;
	}
	catch (...)
	{
		if (what)
			*what = "threw a non-runtime_error exception";
		return false;
	}
}

int main()
{
	std::cout << "\n=== ADR-0004 proper-volume / metric primitive contract ===\n\n";
	std::cout << std::setprecision(17);

	// -----------------------------------------------------------------------
	// A. Regular centre — the EXACT (0,0) point, no tolerance band.
	// -----------------------------------------------------------------------
	std::cout << "A. regular-centre limit (r = 0, m = 0)\n";
	{
		const double f = Geo::MetricDenominator(0.0, 0.0);
		const double lam = Geo::Lambda(0.0, 0.0);
		const double el = Geo::ExpLambda(0.0, 0.0);
		const double wv = Geo::ProperVolumeWeight(0.0, 0.0);

		Report("A1 f == 1 exactly", f == 1.0, "f = " + std::to_string(f));
		Report("A2 Lambda == 0 exactly", lam == 0.0, "Lambda = " + std::to_string(lam));
		Report("A3 expLambda == 1 exactly", el == 1.0, "exp(Lambda) = " + std::to_string(el));
		Report("A4 wV == 0 exactly", wv == 0.0, "wV = " + std::to_string(wv));
		// Lambda is -0.0 here: -0.5 * log(1.0) == -0.0. That is exactly what the legacy
		// `r <= 0 -> denom = 1` branch produced too, so it is bit-identical to pre-3D
		// behavior, and it composes correctly. Signed zero is asserted to be harmless
		// rather than asserted away.
		Report("A5 the centre's signed zero composes correctly",
			   lam == 0.0 && std::exp(lam) == 1.0 && 1.0 / el == 1.0,
			   "Lambda = " + Sci(lam) + " (signbit " + std::to_string(std::signbit(lam)) +
				   "), exp(Lambda) = 1, 1/exp(Lambda) = 1");
	}

	// -----------------------------------------------------------------------
	// B. Fail-closed domain. Every one of these silently produced a number before ADR-0004.
	// -----------------------------------------------------------------------
	std::cout << "\nB. invalid input fails closed\n";
	{
		const double inf = std::numeric_limits<double>::infinity();
		const double nan = std::numeric_limits<double>::quiet_NaN();
		std::string w;

		Report("B1 negative radius", FailsClosed([] { return Geo::ExpLambda(-1.0, 0.5); }, &w), w);
		Report("B2 negative radius, zero mass",
			   FailsClosed([] { return Geo::ExpLambda(-1.0, 0.0); }, &w), w);
		Report("B3 zero radius with non-zero mass",
			   FailsClosed([] { return Geo::ExpLambda(0.0, 1.0); }, &w), w);
		Report("B4 zero radius with tiny non-zero mass (no epsilon band)",
			   FailsClosed([] { return Geo::ExpLambda(0.0, 1e-300); }, &w), w);
		Report("B5 non-finite radius (inf)",
			   FailsClosed([inf] { return Geo::ExpLambda(inf, 1.0); }, &w), w);
		Report("B6 non-finite radius (NaN)",
			   FailsClosed([nan] { return Geo::ExpLambda(nan, 1.0); }, &w), w);
		Report("B7 non-finite mass (inf)",
			   FailsClosed([inf] { return Geo::ExpLambda(10.0, inf); }, &w), w);
		Report("B8 non-finite mass (NaN)",
			   FailsClosed([nan] { return Geo::ExpLambda(10.0, nan); }, &w), w);
		// f == 0 exactly: the Schwarzschild radius, r = 2m.
		Report("B9 f == 0 exactly (horizon)",
			   FailsClosed([] { return Geo::ExpLambda(10.0, 5.0); }, &w), w);
		Report("B10 f < 0 (inside the horizon)",
			   FailsClosed([] { return Geo::ExpLambda(10.0, 8.0); }, &w), w);

		// The contract must reach every entry point, not only the one the caller happened to use.
		Report("B11 MetricDenominator fails closed too",
			   FailsClosed([] { return Geo::MetricDenominator(10.0, 8.0); }, &w), w);
		Report("B12 Lambda fails closed too",
			   FailsClosed([] { return Geo::Lambda(10.0, 8.0); }, &w), w);
		Report("B13 ProperVolumeWeight fails closed too",
			   FailsClosed([] { return Geo::ProperVolumeWeight(10.0, 8.0); }, &w), w);

		// The old behavior this replaces: a 1e-15 clamp returned e^Lambda = 3.162278e+07.
		// If any clamp were reintroduced, B9/B10 would return that instead of throwing.
		Report("B14 no 1e-15 clamp survives (would yield expLambda = 3.162278e+07)",
			   FailsClosed([] { return Geo::ExpLambda(10.0, 5.0); }), "horizon throws, not clamped");
	}

	// -----------------------------------------------------------------------
	// C. Physical domain — against two independent algebraic forms.
	// -----------------------------------------------------------------------
	std::cout << "\nC. metric factor on the physical domain, vs independent forms (<= 2 ULP)\n";
	{
		// 2M/R from 0.05 to 0.98. The upper end is past the Buchdahl bound (8/9) and past any
		// neutron star, deliberately: the primitive must stay correct right up to the horizon.
		const std::vector<double> compactness = {0.05, 0.10, 0.20, 0.30, 0.40, 0.481, 0.50,
												 0.60, 0.70, 0.80, 0.888888888888889, 0.95, 0.98};
		const double R = 12.0; // km; any positive radius exercises the same algebra
		long worst_inv_sqrt = 0, worst_exp_log = 0, worst_f = 0;
		double worst_c = 0.0;

		for (const double c : compactness)
		{
			const double m = 0.5 * c * R;

			// Independence check on f itself: the same subtraction, written out here rather
			// than obtained from the primitive. (Comparing against `1 - c` instead would be
			// comparing against a DIFFERENT quantity -- 2*(0.5*c*R)/R is not bitwise c -- and
			// would measure the fixture's rounding, not the primitive's.)
			const double f_indep = 1.0 - 2.0 * m / R;
			const double f_prim = Geo::MetricDenominator(R, m);
			worst_f = std::max(worst_f, UlpDistance(f_prim, f_indep));

			const double el = Geo::ExpLambda(R, m);

			// The three algebraic forms of the same f, per ADR-0004 SS6.2.
			const long d1 = UlpDistance(el, 1.0 / std::sqrt(f_prim));
			const long d2 = UlpDistance(el, std::exp(-0.5 * std::log(f_prim)));
			if (d1 > worst_inv_sqrt)
			{
				worst_inv_sqrt = d1;
				worst_c = c;
			}
			worst_exp_log = std::max(worst_exp_log, d2);
		}

		Report("C0 f is bitwise the canonical 1 - 2m/r subtraction", worst_f == 0,
			   "worst = " + std::to_string(worst_f) + " ULP");
		Report("C1 expLambda vs 1/sqrt(f) within 2 ULP", worst_inv_sqrt <= 2,
			   "worst = " + std::to_string(worst_inv_sqrt) + " ULP at 2m/R = " +
				   std::to_string(worst_c));
		Report("C2 expLambda vs exp(-0.5 log f) is bitwise", worst_exp_log == 0,
			   "worst = " + std::to_string(worst_exp_log) + " ULP");

		// Monotonicity: e^Lambda must grow without bound as the horizon is approached.
		double prev = 0.0;
		bool monotone = true;
		for (const double c : compactness)
		{
			const double el = Geo::ExpLambda(R, 0.5 * c * R);
			if (el <= prev)
				monotone = false;
			prev = el;
		}
		Report("C3 expLambda increases monotonically with compactness", monotone,
			   "e^Lambda(0.98) = " + std::to_string(prev));
	}

	// -----------------------------------------------------------------------
	// D. Proper-volume composition.
	// -----------------------------------------------------------------------
	std::cout << "\nD. proper-volume weight composition\n";
	{
		long worst = 0;
		for (const double c : {0.05, 0.20, 0.481, 0.70, 0.95})
			for (const double r : {1e-5, 0.5, 5.0, 12.0, 70.0})
			{
				const double m = 0.5 * c * r;
				const double got = Geo::ProperVolumeWeight(r, m);
				const double want = (4.0 * M_PI) * (r * r) * Geo::ExpLambda(r, m);
				worst = std::max(worst, UlpDistance(got, want));
			}
		Report("D1 wV == 4*pi*r^2*expLambda bitwise", worst == 0,
			   "worst = " + std::to_string(worst) + " ULP");

		// Flat-space limit: with m = 0 the weight is exactly the Euclidean shell area.
		const double r = 7.5;
		Report("D2 m = 0 gives exactly 4*pi*r^2 (flat space)",
			   Geo::ProperVolumeWeight(r, 0.0) == (4.0 * M_PI) * (r * r),
			   "wV(7.5, 0) = " + std::to_string(Geo::ProperVolumeWeight(r, 0.0)));

		// The measure must be strictly larger than the coordinate one whenever m > 0.
		Report("D3 wV > coordinate 4*pi*r^2 for m > 0",
			   Geo::ProperVolumeWeight(12.0, 2.0) > (4.0 * M_PI) * 144.0,
			   "proper > coordinate");

		// The primitive owns no unit conversion. If 1e54 ever migrated in, wV at r = 1 km,
		// m = 0 would be 1e54 times too large; this pins the magnitude to O(4*pi).
		const double w_unit = Geo::ProperVolumeWeight(1.0, 0.0);
		Report("D4 primitive carries no baryon-density unit conversion",
			   w_unit > 12.0 && w_unit < 13.0,
			   "wV(1 km, 0) = " + std::to_string(w_unit) + " (must be 4*pi, not 4*pi*1e54)");
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "ADR-0004 primitive contract holds"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
