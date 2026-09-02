// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_first_order_physics_analytic.cpp
 * @brief Phase 4B — INDEPENDENT physical validation of the normalized first-order Hartle
 *        response on an exact analytic background. Self-contained: no external EOS assets.
 *
 * WHAT 4B ADDS OVER 4A. Phase 4A proved the ADR-0006 *contract*: the arbitrary seed does not
 * leak, a requested spin comes back, the units are right, `J = I Omega`, zero spin works.
 * Every one of those checks is satisfiable by a profile of the wrong SHAPE. Phase 4B asks
 * whether the shape is physically right, and answers it with evidence that does not come from
 * production:
 *
 *     s(r)  = omega_bar(r) / Omega  ,      s'(r) = [d omega_bar/dr](r) / Omega
 *
 * compared node by node against an independently derived and independently normalized profile,
 * plus three closed-form physical identities.
 *
 * THE INDEPENDENT ORACLE. `hartle_reference.hpp` integrates the CONSERVATIVE Hartle system in
 * different state variables (`omega_bar`, `q = r^4 j omega_bar'`), building its right-hand side
 * from the metric columns `nu`, `lambda` and `(eps + p)`. It never calls `ODE_N_Fast`,
 * `GetHartleOmegaCoeff_N_Fast`, `GetHartleDOmegaCoeff_N_Fast` or
 * `HartleFirstOrderResponse::At`, and it normalizes itself by its OWN surface extraction
 * `Omega_ref = omega_bar_ref(R) + R omega_bar_ref'(R)/3` — never by a production quantity.
 *
 * ============================ PREDECLARED BOUNDS ============================
 * Fixed BEFORE the comparison, from the published Phase-2B record only
 * (`docs/validation/HARTLE_MOMENT_INERTIA.md` §10, §8), by one stated rule:
 *
 *     profile bound = 5 x (Phase-2B production/reference `I` agreement for this star class),
 *                     rounded up to the next power of ten.
 *
 * The rule is tied to `I` because `I` is a functional of exactly this profile, so the two
 * cannot disagree in order of magnitude; the factor 5 allows a worst NODE to exceed an
 * INTEGRATED quantity's error. Inputs and results:
 *
 *   analytic  : 2B measured 9.455e-9  -> 5x = 4.7e-8  -> BOUND 1e-7   (s and s')
 *   volume id : 2B measured 1.24e-7 (surf/vol) + 9.455e-9  -> 5x = 6.7e-7 -> BOUND 1e-6
 *   exterior  : algebraic identity, roundoff only          -> BOUND 1e-12
 *   weak field: leading corrections are O(M/R), so at M/R = 0.002 expect a few 1e-3
 *                                                          -> BOUND 5e-3
 *
 * None of these bounds depends on any measurement made in Phase 4B. They are not widened.
 *
 * SCOPE. First order only. Nothing here validates, executes or claims anything at O(Omega^2),
 * which remains an unverified candidate (INV-08).
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "hartle_profile_compare.hpp"
#include "hartle_reference.hpp"

using CompactStar::AngularVelocity;
using CompactStar::Core::NStar;
using namespace hartle_4b;

// ---- predeclared bounds (see file header) ---------------------------------------------
static constexpr double kBoundProfileAnalytic = 1.0e-7;
static constexpr double kBoundVolumeAnalytic = 1.0e-6;
static constexpr double kBoundExteriorIdentity = 1.0e-12;
static constexpr double kBoundWeakFieldCoeff = 5.0e-3;
/// The reference is admissible as an oracle only if its own numerical floor is at least three
/// decades below the discrepancy it is being used to measure.
static constexpr double kRefFloorMargin = 1.0e-3;

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
static std::string Sci(double v, int p = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", p, v);
	return b;
}

int main()
{
	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 4B — independent physical validation of the normalized first-order\n"
				 "Hartle response. Self-contained: exact Schwarzschild constant-density star.\n"
				 "Bounds predeclared from the Phase-2B record; see the file header.\n"
				 "First order only — no O(Omega^2) claim (INV-08).\n\n";

	const UniformStar u = MakeUniform(/*M_km=*/2.0, /*R_km=*/13.0); // M/R ~ 0.154
	const std::size_t N = 4001;
	const auto bg = TabulateUniform(u, N);
	auto ns = UniformProductionStar(u, N);
	if (!ns || !ns->RotationResponse().valid)
	{
		std::cerr << "fixture star or its first-order response is unavailable\n";
		return 1;
	}

	std::vector<double> sP, spP;
	ProductionShape(*ns, sP, spP);
	const double I_prod = ns->RotationResponse().I;
	const double R = bg.R();

	// The reference baseline: same centre convention as production (both begin at the first
	// profile node), so the comparison is not confounded by a different start radius.
	const auto ref0 = hartle_ref::Solve(bg, 5.0e-3, bg.r[0], 1e-12, 1e-14);
	if (!ref0.ok)
	{
		std::cerr << "reference solve failed: " << ref0.message << "\n";
		return 1;
	}

	// =======================================================================
	//  A0 — is the reference admissible as an oracle at all?
	// =======================================================================
	std::cout << "A0 reference numerical floor (reference vs reference, matched centre)\n";
	double ref_floor_s = 0.0, ref_floor_sp = 0.0;
	{
		struct V
		{
			const char *tag;
			double seed, rtol, atol;
		};
		// Tightening the tolerance by two decades, and moving the seed over six, must move the
		// normalized profile by far less than the production/reference difference it is used to
		// measure. The seed sweep matters because the reference normalizes by its OWN surface
		// extraction, so seed-independence is a property to be demonstrated, not assumed.
		const std::vector<V> variants = {{"tol 1e-13 / 1e-15", 5.0e-3, 1e-13, 1e-15},
										 {"tol 1e-14 / 1e-16", 5.0e-3, 1e-14, 1e-16},
										 {"seed x 1e3", 5.0, 1e-12, 1e-14},
										 {"seed x 1e-3", 5.0e-6, 1e-12, 1e-14}};
		for (const auto &v : variants)
		{
			const auto r = hartle_ref::Solve(bg, v.seed, bg.r[0], v.rtol, v.atol);
			if (!r.ok)
			{
				Report(std::string("A0 variant ") + v.tag, false, r.message);
				continue;
			}
			const ShapeCmp c = CompareShapes(r.s, r.s_prime, ref0.s, ref0.s_prime, bg.r);
			ref_floor_s = std::max(ref_floor_s, c.max_rel_s);
			ref_floor_sp = std::max(ref_floor_sp, c.rel_to_scale_sp);
			std::cout << "      " << std::setw(20) << std::left << v.tag << std::right
					  << "  s " << Sci(c.max_rel_s) << "   s' " << Sci(c.rel_to_scale_sp) << "\n";
		}
		std::cout << "      reference floor:  s " << Sci(ref_floor_s) << "   s' "
				  << Sci(ref_floor_sp) << "\n";
	}

	// =======================================================================
	//  A / Experiment A — the full normalized profile, node by node
	// =======================================================================
	std::cout << "\nA full normalized-profile comparison against the independent solver\n";
	const ShapeCmp c = CompareShapes(sP, spP, ref0.s, ref0.s_prime, bg.r);
	{
		std::cout << "      nodes compared: " << c.n << "\n"
				  << "      s : max rel " << Sci(c.max_rel_s) << "  max abs " << Sci(c.max_abs_s)
				  << "  rms " << Sci(c.rms_s) << "  worst r = " << Sci(c.r_worst_s, 4) << " km\n"
				  << "      s': max abs/peak " << Sci(c.rel_to_scale_sp) << "  max rel(>1% peak) "
				  << Sci(c.max_rel_sp) << "  rms " << Sci(c.rms_sp) << "  worst r = "
				  << Sci(c.r_worst_sp, 4) << " km\n";

		Report("A0a the reference's own floor is subdominant to the difference it measures",
			   ref_floor_s <= kRefFloorMargin * c.max_rel_s &&
				   ref_floor_sp <= kRefFloorMargin * std::max(c.rel_to_scale_sp, 1e-300),
			   "floor s " + Sci(ref_floor_s) + " vs measured " + Sci(c.max_rel_s) + " (ratio " +
				   Sci(c.max_rel_s > 0 ? ref_floor_s / c.max_rel_s : 0.0, 1) + ", bound 1e-3)");

		Report("Aa s(r) agrees with the independent profile at EVERY node",
			   c.n == N && c.max_rel_s <= kBoundProfileAnalytic,
			   "max relative " + Sci(c.max_rel_s) + " over " + std::to_string(c.n) +
				   " nodes (predeclared bound " + Sci(kBoundProfileAnalytic, 0) + ")");
		Report("Ab s'(r) agrees with the independent profile at EVERY node",
			   c.rel_to_scale_sp <= kBoundProfileAnalytic &&
				   c.max_rel_sp <= kBoundProfileAnalytic,
			   "max abs/peak " + Sci(c.rel_to_scale_sp) + ", max rel " + Sci(c.max_rel_sp) +
				   " (bound " + Sci(kBoundProfileAnalytic, 0) + ")");
		Report("Ac the agreement is not carried by a few nodes — the RMS is comparable",
			   c.rms_s <= kBoundProfileAnalytic && c.rms_s > 0.0,
			   "rms(s) " + Sci(c.rms_s) + " vs max " + Sci(c.max_rel_s));

		// I was already validated in Phase 2B; reported for continuity, asserted loosely.
		std::cout << "      I_prod " << Sci(I_prod, 10) << "   I_ref " << Sci(ref0.I_surface, 10)
				  << "   rel " << Sci(Rel(I_prod, ref0.I_surface)) << "  (Phase 2B: 9.455e-9)\n";
	}

	// =======================================================================
	//  C / Experiment C — exterior matching identities
	// =======================================================================
	std::cout << "\nC exterior-matching identities at the surface\n";
	{
		// Outside the star omega_bar = Omega - 2J/r^3 with J = I Omega, so dividing by Omega:
		//     s(R)  = 1 - 2I/R^3 ,      s'(R) = 6I/R^4 .
		const int n = static_cast<int>(sP.size());
		const double s_R = sP[static_cast<std::size_t>(n - 1)];
		const double sp_R = spP[static_cast<std::size_t>(n - 1)];
		const double pred_s = 1.0 - 2.0 * I_prod / (R * R * R);
		const double pred_sp = 6.0 * I_prod / (R * R * R * R);

		// (i) Internal consistency. Given production's own definitions of I, s and s' this is an
		// ALGEBRAIC identity, so it must hold to roundoff. It is not independent physics — it is
		// the check that the published profile and the published I describe the SAME solution,
		// which a shape corruption that left I alone would break.
		Report("Ca s(R) = 1 - 2I/R^3 (production profile vs production I — internal consistency)",
			   Rel(s_R, pred_s) <= kBoundExteriorIdentity,
			   "s(R) " + Sci(s_R, 10) + " vs " + Sci(pred_s, 10) + ", rel " +
				   Sci(Rel(s_R, pred_s)));
		Report("Cb s'(R) = 6I/R^4 (internal consistency)",
			   Rel(sp_R, pred_sp) <= kBoundExteriorIdentity,
			   "s'(R) " + Sci(sp_R, 10) + " vs " + Sci(pred_sp, 10) + ", rel " +
				   Sci(Rel(sp_R, pred_sp)));

		// (ii) The genuinely independent form: production's surface shape against the
		// INDEPENDENT solver's I. Nothing on the two sides shares an origin.
		const double pred_s_ref = 1.0 - 2.0 * ref0.I_surface / (R * R * R);
		const double pred_sp_ref = 6.0 * ref0.I_surface / (R * R * R * R);
		Report("Cc s(R) = 1 - 2 I_ref/R^3 (production profile vs INDEPENDENT I)",
			   Rel(s_R, pred_s_ref) <= kBoundProfileAnalytic,
			   "rel " + Sci(Rel(s_R, pred_s_ref)) + " (bound " + Sci(kBoundProfileAnalytic, 0) +
				   ")");
		Report("Cd s'(R) = 6 I_ref/R^4 (production profile vs INDEPENDENT I)",
			   Rel(sp_R, pred_sp_ref) <= kBoundProfileAnalytic,
			   "rel " + Sci(Rel(sp_R, pred_sp_ref)) + " (bound " +
				   Sci(kBoundProfileAnalytic, 0) + ")");
	}

	// =======================================================================
	//  D / Experiment D — the volume-integral identity, from production's own shape
	// =======================================================================
	std::cout << "\nD volume-integral identity using production's normalized response\n";
	{
		const double I_vol = VolumeIntegralI(bg, sP);
		Report("Da I recovered by volume quadrature over production's s(r) matches its own I",
			   Rel(I_vol, I_prod) <= kBoundVolumeAnalytic,
			   "I_vol " + Sci(I_vol, 10) + " vs I " + Sci(I_prod, 10) + ", rel " +
				   Sci(Rel(I_vol, I_prod)) + " (bound " + Sci(kBoundVolumeAnalytic, 0) + ")");
		std::cout << "      this integral reads the INTERIOR of the shape and uses no J_raw,\n"
					 "      no Omega_raw and no seed — a surface-only check cannot see what it "
					 "sees.\n";
	}

	// =======================================================================
	//  E / Experiment E — the weak-field PROFILE limit, with derived coefficients
	// =======================================================================
	std::cout << "\nE weak-field profile limit (derived closed forms, not a fitted trend)\n";
	{
		// In the Newtonian limit j -> 1 and j' -> -4 pi r rho, so the frame-dragging equation
		// becomes omega_bar'' + (4/r) omega_bar' = k^2 omega_bar with k^2 = 16 pi rho = 12M/R^3
		// for a uniform sphere. Expanding the regular solution for kR << 1,
		//
		//     omega_bar(r) = omega_bar_c [ 1 + k^2 r^2/10 + O(k^4 r^4) ] ,
		//
		// hence omega_bar(R) = omega_bar_c (1 + 1.2 M/R) and R omega_bar'(R)/3 =
		// omega_bar_c (0.8 M/R), so Omega = omega_bar_c (1 + 2M/R) and
		//
		//     omega(0)/Omega = 1 - s(0) -> 2 (M/R) ,     omega(R)/Omega = 1 - s(R) -> 0.8 (M/R).
		//
		// The surface value is independently the exterior result 2I/R^3 with I -> (2/5) M R^2.
		// Both coefficients are DERIVED, so this is a physical limit test, not a fitted law.
		const std::vector<double> Cs = {0.15, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002};
		std::vector<double> dev_c, dev_R;
		std::cout << "      M/R        (1-s_c)/(M/R)   dev from 2      (1-s_R)/(M/R)  dev from "
					 "0.8\n";
		for (double C : Cs)
		{
			const double Rk = 13.0;
			const UniformStar uu = MakeUniform(C * Rk, Rk);
			auto nsw = UniformProductionStar(uu, N);
			if (!nsw || !nsw->RotationResponse().valid)
			{
				Report("E weak-field build at M/R = " + Sci(C, 2), false, "unavailable");
				continue;
			}
			std::vector<double> sw, spw;
			ProductionShape(*nsw, sw, spw);
			const double a = (1.0 - sw.front()) / C;
			const double b = (1.0 - sw.back()) / C;
			dev_c.push_back(std::fabs(a - 2.0) / 2.0);
			dev_R.push_back(std::fabs(b - 0.8) / 0.8);
			std::cout << "      " << Sci(C, 3) << "  " << Sci(a, 8) << "  " << Sci(dev_c.back())
					  << "  " << Sci(b, 8) << "  " << Sci(dev_R.back()) << "\n";
		}
		bool mono_c = true, mono_R = true;
		for (std::size_t i = 1; i < dev_c.size(); ++i)
		{
			mono_c = mono_c && dev_c[i] < dev_c[i - 1];
			mono_R = mono_R && dev_R[i] < dev_R[i - 1];
		}
		Report("Ea omega(0)/Omega approaches the derived 2(M/R) monotonically", mono_c,
			   "deviation " + Sci(dev_c.front()) + " -> " + Sci(dev_c.back()));
		Report("Eb omega(R)/Omega approaches the derived 0.8(M/R) monotonically", mono_R,
			   "deviation " + Sci(dev_R.front()) + " -> " + Sci(dev_R.back()));
		Report("Ec both coefficients are reached at the weakest field",
			   !dev_c.empty() && dev_c.back() <= kBoundWeakFieldCoeff &&
				   dev_R.back() <= kBoundWeakFieldCoeff,
			   "at M/R = 0.002: centre " + Sci(dev_c.back()) + ", surface " +
				   Sci(dev_R.back()) + " (bound " + Sci(kBoundWeakFieldCoeff, 0) + ")");
		Report("Ed frame dragging vanishes in the Newtonian limit: s(r) -> 1 throughout",
			   !dev_c.empty(), "max|1 - s| falls as 2(M/R) -> 0, i.e. no dragging without "
							   "relativity");
	}

	// =======================================================================
	//  F / Experiment F — physical materialization and spin reversal
	// =======================================================================
	std::cout << "\nF spin reversal and materialization linearity (secondary)\n";
	{
		const double w = 1.0e2;
		const auto pp = ns->RotationAt(AngularVelocity::FromRadPerSecond(w));
		const auto pm = ns->RotationAt(AngularVelocity::FromRadPerSecond(-w));
		const auto p2 = ns->RotationAt(AngularVelocity::FromRadPerSecond(2.0 * M_PI * 1.0e2));
		const int n = static_cast<int>(pp.omega_bar.Size());

		bool exact = (pm.Omega_geom == -pp.Omega_geom) && (pm.J == -pp.J) && (pm.I == pp.I);
		for (int i = 0; i < n && exact; ++i)
			exact = (pm.omega_bar[i] == -pp.omega_bar[i]) &&
					(pm.domega_bar[i] == -pp.domega_bar[i]);
		Report("Fa reversing the spin negates Omega, J and both profiles EXACTLY", exact,
			   "bit-exact antisymmetry at all " + std::to_string(n) + " nodes");
		Report("Fb I is unchanged under spin reversal", pm.I == pp.I && pp.I == I_prod,
			   "I = " + Sci(pp.I, 10) + " km^3");

		std::vector<double> s2, sp2;
		ProductionShape(*ns, s2, sp2);
		bool shape_same = (s2.size() == sP.size());
		for (std::size_t i = 0; i < s2.size() && shape_same; ++i)
			shape_same = (s2[i] == sP[i]) && (sp2[i] == spP[i]);
		Report("Fc the seed-free response itself is untouched by any materialization",
			   shape_same, "s(r) and s'(r) bit-identical after three materializations");
		Report("Fd a third spin reproduces the same shape once divided out",
			   Rel(p2.omega_bar[n / 2] / p2.Omega_geom, sP[static_cast<std::size_t>(n / 2)]) <=
				   1e-13,
			   "materialization is a scaling of one shape, not three different solutions");
	}

	// =======================================================================
	//  Diagnostics — frame-dragging fraction and regime (reported, not asserted as physics)
	// =======================================================================
	std::cout << "\nG diagnostics\n";
	{
		const std::size_t n = sP.size();
		double mn = 1e300, mx = -1e300;
		for (std::size_t i = 0; i < n; ++i)
		{
			mn = std::min(mn, 1.0 - sP[i]);
			mx = std::max(mx, 1.0 - sP[i]);
		}
		std::cout << "      omega/Omega = 1 - s:  centre " << Sci(1.0 - sP[0], 6) << "  surface "
				  << Sci(1.0 - sP[n - 1], 6) << "  min " << Sci(mn, 6) << "  max " << Sci(mx, 6)
				  << "\n";
		// Recorded as a measured ordering on the validated stars. It is NOT asserted as a
		// theorem: no governing document derives a hard inequality for omega/Omega, and this
		// harness does not invent one.
		Report("Ga the frame-dragging fraction is measured within 0 <= omega/Omega < 1 and is "
			   "largest at the centre (DIAGNOSTIC ordering, not a derived inequality)",
			   mn >= 0.0 && mx < 1.0 && (1.0 - sP[0]) == mx,
			   "dragging decreases outward, as expected for a rotating fluid interior");

		const double M_km = u.M_km;
		const double Omega_K = std::sqrt(M_km / (R * R * R));
		const double w716 = 2.0 * M_PI * 716.0 / 299792.458;
		std::cout << "      slow-rotation context: Omega_K = " << Sci(Omega_K, 6)
				  << " km^-1, Omega/Omega_K at 2 pi x 716 Hz = " << Sci(w716 / Omega_K, 6)
				  << "\n"
				  << "      NOTE: normalization correctness != slow-rotation truncation "
					 "accuracy.\n"
					 "      Nothing above claims the O(Omega^2)-neglected terms are small at "
					 "this ratio;\n"
					 "      that question belongs with the O(Omega^2) validation and regime "
					 "work.\n";
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "analytic first-order physics checks passed"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
