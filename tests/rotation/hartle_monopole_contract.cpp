#include "tests/relativity/fixture_units.hpp"
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_contract.cpp
 * @brief Phase 4C-I1 — implementation contract for the governed O(Omega^2) monopole response
 *        (ADR-0007, ACCEPTED 2026-09-02). Self-contained: no external EOS assets.
 *
 * WHAT THIS IS, AND IS NOT.
 *
 * This file protects the *contract* of the replacement that retired the `675b4a9` candidate:
 * fail-closed inputs, the regular-centre start, provenance and staleness, seed invariance,
 * the materialization algebra, and the one-solve-per-request performance rule.
 *
 * **It is NOT independent physical validation of the ODE.** No oracle here computes the
 * monopole solution by another route, so nothing here can tell you the numbers are right.
 * That is Phase 4D — regular-centre series against an independent solver in different
 * variables, Newtonian homogeneous-star limits, exterior identities, published Hartle-Thorne
 * results, convergence, EOS-derivative sensitivity, and the M1-M9 detectors. Until 4D,
 * INV-08 stands at "source conformed, independent physical validation pending", and no
 * baseline of these numbers may be created.
 *
 * THE FIXTURE. The exact Schwarzschild constant-density interior (Schwarzschild 1916; MTW
 * §23.7) — the same star the validated Phase-2B/4A/4B rotation harnesses use. Its EOS
 * derivative is known in closed form: incompressible matter has `d(eps)/dp = 0` exactly, and
 * it is supplied explicitly through the governed `TOVPoint` mechanism from Phase 4C-I0. That
 * exercises the derivative path with a value no reconstruction could have invented.
 *
 * ======================== PREDECLARED BOUNDS ========================
 * Fixed BEFORE any measurement, and inherited from the accepted contract rather than chosen
 * here:
 *   seed invariance of every coefficient   rel <= 1e-10   ADR-0007 §7 item 1 (= ADR-0006 §7-1;
 *                                                         the production driver's absolute
 *                                                         tolerance is 1e-10)
 *   quadratic materialization              rel <= 1e-14   ADR-0007 §7 item 2 (pure arithmetic)
 *   delta_M internal identity              rel <= 1e-14   pure arithmetic on published fields
 * ====================================================================
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#include <Zaki/Physics/Constants.hpp>

// ---------------------------------------------------------------------------
//  PRIVILEGED TEST BACKDOOR — NOT SUPPORTED SCIENTIFIC API.
//
//  `RotationSolverTestSeam` is declared, never defined, in production; a harness defines it
//  to reach internal numerical state that must not become public API (ADR-0006 Q2). Here it
//  reaches the arbitrary first-order seed — so seed invariance of the SECOND-order
//  coefficients can be *proved* rather than asserted — and the two solve counters, so
//  "one integration per explicit request, none per materialization" is measured.
// ---------------------------------------------------------------------------
namespace CompactStar::Core
{
struct RotationSolverTestSeam
{
	static void SetSeed(RotationSolver &rs, double seed) { rs.seed_omega_bar_ = seed; }
	static double DefaultSeed() { return RotationSolver::kDefaultInitOmegaBar; }
	static void SetSeed(NStar &ns, double seed) { ns.rot_solver.seed_omega_bar_ = seed; }

	static std::size_t FirstOrderSolves(const NStar &ns)
	{
		return ns.rot_solver.first_order_solve_count_;
	}
	static std::size_t MonopoleSolves(const NStar &ns)
	{
		return ns.rot_solver.monopole_solve_count_;
	}

	/// Rebuild a star's profile in place. `BuildFromTOV` is private and the public constructors
	/// build immediately, so this is the only way to install a NON-DEFAULT seed *before* the
	/// first-order solve runs — which is what proving seed invariance requires.
	static void Rebuild(NStar &ns, const std::vector<TOVPoint> &pts)
	{
		ns.BuildFromTOV(pts, nullptr);
	}

	/// Bump the profile version without touching a single datum, so staleness can be tested
	/// as the pure provenance question it is.
	static void TouchProfile(NStar &ns) { ns.prof_.Touch(); }
};
} // namespace CompactStar::Core

using CompactStar::AngularVelocity;
using CompactStar::Core::HartleMonopoleResponse;
using CompactStar::Core::NStar;
using CompactStar::Core::PhysicalHartleMonopole;
using CompactStar::Core::RotationSolverTestSeam;
using CompactStar::Core::StarProfile;
using CompactStar::Core::TOVPoint;

static constexpr double kBoundSeedInvariance = 1.0e-10;
static constexpr double kBoundQuadratic = 1.0e-14;
static constexpr double kBoundDeltaMIdentity = 1.0e-14;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id;
	if (!d.empty())
		std::cout << " — " << d;
	std::cout << "\n";
	if (!ok)
		++g_fail;
}
static double Rel(double a, double b)
{
	return std::fabs(b) > 0.0 ? std::fabs(a - b) / std::fabs(b) : std::fabs(a - b);
}
static std::string Sci(double v, int p = 4)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", p, v);
	return b;
}

// ---------------------------------------------------------------------------
//  Exact Schwarzschild constant-density interior.
// ---------------------------------------------------------------------------
struct UniformStar
{
	double M_km = 2.0, R_km = 13.0, rho0 = 0.0; // rho0 in km^-2
	double yR() const { return std::sqrt(1.0 - 2.0 * M_km / R_km); }
	double y(double r) const
	{
		return std::sqrt(1.0 - 2.0 * M_km * r * r / (R_km * R_km * R_km));
	}
	double m(double r) const { return M_km * r * r * r / (R_km * R_km * R_km); }
	double p(double r) const { return rho0 * (y(r) - yR()) / (3.0 * yR() - y(r)); }
	double nuprime(double r) const
	{
		return 2.0 * M_km * r / (R_km * R_km * R_km * y(r) * (3.0 * yR() - y(r)));
	}
};

static UniformStar MakeUniform(double M_km = 2.0, double R_km = 13.0)
{
	UniformStar u;
	u.M_km = M_km;
	u.R_km = R_km;
	u.rho0 = 3.0 * M_km / (4.0 * M_PI * R_km * R_km * R_km);
	return u;
}

/// `supply_dedp == false` reproduces a star that carries no EOS-derivative authority.
static std::vector<TOVPoint> UniformPoints(const UniformStar &u, std::size_t N,
										   bool supply_dedp)
{
	const double km2_to_gcm3 =
		relativity_fixture::eps_to_rho;
	const double km2_to_dyn =
		relativity_fixture::pressure_to_cgs;
	const double r0 = 1.0e-5;

	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r =
			r0 + static_cast<double>(i) * (u.R_km - r0) / static_cast<double>(N - 1);
		// Incompressible matter: d(eps)/dp = 0 EXACTLY. Supplied explicitly through the
		// Phase-4C-I0 mechanism, never reconstructed from the radial profile.
		const double dedp = supply_dedp ? 0.0 : std::numeric_limits<double>::quiet_NaN();
		pts.emplace_back(r, u.m(r) / relativity_fixture::solar_km, u.nuprime(r) / 1.0e5, 0.0,
						 u.p(r) * km2_to_dyn, u.rho0 * km2_to_gcm3, 0.1,
						 std::vector<double>{}, dedp);
	}
	return pts;
}

static std::unique_ptr<NStar> BuildStar(const UniformStar &u, std::size_t N, bool supply_dedp)
{
	return std::make_unique<NStar>(UniformPoints(u, N, supply_dedp));
}

/// A star whose solver carries a NON-DEFAULT internal first-order seed. The seed must be in
/// place before the first-order solve runs, and the public constructors solve immediately, so
/// the profile is installed afterwards through the privileged seam.
static std::unique_ptr<NStar> BuildSeeded(const UniformStar &u, std::size_t N, double seed)
{
	const auto pts = UniformPoints(u, N, /*supply_dedp=*/true);
	auto ns = std::make_unique<NStar>();
	RotationSolverTestSeam::SetSeed(*ns, seed);
	RotationSolverTestSeam::Rebuild(*ns, pts);
	return ns;
}

int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout
		<< "Phase 4C-I1 — governed O(Omega^2) monopole IMPLEMENTATION CONTRACT (ADR-0007).\n"
		   "Self-contained. This protects contracts, NOT the physics: independent scientific\n"
		   "validation of the ODE is Phase 4D (INV-08). No baseline is produced.\n\n";

	const UniformStar u = MakeUniform();
	const std::size_t N = 601;

	// ------------------------------------------------------------------
	//  C1 — fail closed without an authoritative d(eps)/dp
	// ------------------------------------------------------------------
	std::cout << "C1/C2. EOS-derivative gate (ADR-0007 P5)\n";
	{
		auto ns = BuildStar(u, N, /*supply_dedp=*/false);
		const bool has = ns->Profile().HasEosDEdP();
		const bool ok_compute = ns->ComputeHartleMonopoleResponse();
		const auto *resp = ns->MonopoleResponse();
		bool threw = false;
		try
		{
			(void)ns->MonopoleAt(AngularVelocity::FromRadPerSecond(100.0));
		}
		catch (const std::exception &)
		{
			threw = true;
		}
		Report("C1 a star with no authoritative d(eps)/dp fails closed",
			   !has && !ok_compute && resp == nullptr && threw,
			   "HasEosDEdP=false, Compute=false, MonopoleResponse=nullptr, MonopoleAt threw");
	}

	auto ns = BuildStar(u, N, /*supply_dedp=*/true);
	Report("C2 an explicit d(eps)/dp = 0 (incompressible) is accepted",
		   ns->Profile().HasEosDEdP(), "supplied through the governed TOVPoint mechanism");

	// ------------------------------------------------------------------
	//  C14 — computing confers no spin; C11/§31 — solve accounting
	// ------------------------------------------------------------------
	std::cout << "\nC11/C14. Explicit computation, no implicit spin, solve accounting\n";
	{
		const std::size_t mono_before = RotationSolverTestSeam::MonopoleSolves(*ns);
		Report("C11a ordinary construction runs ZERO monopole solves", mono_before == 0,
			   "monopole solves after construction = " + std::to_string(mono_before));
		Report("C11b ordinary construction did run the first-order solve",
			   RotationSolverTestSeam::FirstOrderSolves(*ns) == 1, "first-order solves = 1");
		Report("C14a before computing, no response is published",
			   ns->MonopoleResponse() == nullptr, "MonopoleResponse() == nullptr");
	}

	const bool computed = ns->ComputeHartleMonopoleResponse();
	Report("C2b the governed response computes on the analytic star", computed, "");
	if (!computed)
	{
		std::cout << "\nFAIL — cannot continue without a response\n";
		return 1;
	}
	const HartleMonopoleResponse *R0 = ns->MonopoleResponse();
	Report("C11c exactly ONE monopole integration was performed",
		   RotationSolverTestSeam::MonopoleSolves(*ns) == 1,
		   "monopole solves = " + std::to_string(RotationSolverTestSeam::MonopoleSolves(*ns)));
	{
		const bool again = ns->ComputeHartleMonopoleResponse();
		Report("C11d recomputing on an UNCHANGED profile reuses the cache (still one solve)",
			   again && RotationSolverTestSeam::MonopoleSolves(*ns) == 1,
			   "monopole solves = " +
				   std::to_string(RotationSolverTestSeam::MonopoleSolves(*ns)));
	}
	Report("C14b the normalized response carries no angular velocity — it is per Omega^2",
		   R0 != nullptr && R0->valid,
		   "fields are *_over_Omega2; a physical spin exists only via MonopoleAt()");

	// ------------------------------------------------------------------
	//  C3/C4 — shape and finiteness
	// ------------------------------------------------------------------
	std::cout << "\nC3/C4. Response shape and finiteness\n";
	{
		const std::size_t n = ns->Profile().GetRadius()->Size();
		const bool sizes = R0->m0_over_Omega2.Size() == n &&
						   R0->p0star_over_Omega2.Size() == n &&
						   R0->delta_p0_over_Omega2.Size() == n &&
						   R0->xi0_over_Omega2.Size() == n;
		Report("C3 every coefficient column has one value per radial node", sizes,
			   "n = " + std::to_string(n));

		bool finite = std::isfinite(R0->deltaM_over_Omega2) &&
					  std::isfinite(R0->surface_shell_mass_over_Omega2) &&
					  std::isfinite(R0->surface_xi0_over_Omega2) &&
					  std::isfinite(R0->R_surface) && std::isfinite(R0->I);
		for (int i = 0; finite && i < static_cast<int>(n); ++i)
			finite = std::isfinite(R0->m0_over_Omega2[i]) &&
					 std::isfinite(R0->p0star_over_Omega2[i]) &&
					 std::isfinite(R0->delta_p0_over_Omega2[i]) &&
					 std::isfinite(R0->xi0_over_Omega2[i]);
		Report("C4 every published value is finite", finite, "");
	}

	// ------------------------------------------------------------------
	//  C5 — the regular-centre start (ADR-0007 P4)
	// ------------------------------------------------------------------
	std::cout << "\nC5. Regular-centre initialization (ADR-0007 P4)\n";
	{
		const auto &P = ns->Profile();
		const auto *rc = P.GetRadius();
		const auto *ec = P.GetEnergyDensity();
		const auto *pc = P.GetPressure();
		const auto *mc = P.GetMass();
		const auto *nuc = P.GetMetricNu();
		const auto *dc = P.GetEosDEdP();
		const auto &fo = ns->RotationResponse();

		const double r0 = (*rc)[0];
		const double j2 = std::exp(-2.0 * (*nuc)[0]) * (1.0 - 2.0 * (*mc)[0] / r0);
		const double s0 = fo.omega_bar_over_Omega[0];

		const double phat_expect = (1.0 / 3.0) * j2 * s0 * s0 * r0 * r0;
		const double mhat_expect = (4.0 * M_PI / 15.0) * ((*ec)[0] + (*pc)[0]) *
								   ((*dc)[0] + 2.0) * j2 * s0 * s0 * std::pow(r0, 5);

		const double rel_p = Rel(R0->p0star_over_Omega2[0], phat_expect);
		const double rel_m = Rel(R0->m0_over_Omega2[0], mhat_expect);
		Report("C5a p0*_hat(r0) equals the governed regular series (1/3) j_c^2 s_c^2 r0^2",
			   rel_p <= 1e-12,
			   "got " + Sci(R0->p0star_over_Omega2[0]) + " expect " + Sci(phat_expect) +
				   "  rel=" + Sci(rel_p));
		Report("C5b m0_hat(r0) equals (4pi/15)(eps+p)[(deps/dp)+2] j_c^2 s_c^2 r0^5",
			   rel_m <= 1e-12,
			   "got " + Sci(R0->m0_over_Omega2[0]) + " expect " + Sci(mhat_expect) +
				   "  rel=" + Sci(rel_m));
		std::cout << "     (NOT the retired candidate's literal {0,0} start: the series makes\n"
					 "      the fixed-eps_c family exact to rounding, ADR-0007 P4.)\n";
	}

	// ------------------------------------------------------------------
	//  C13 — delta_M internal identity (ADR-0007 P6)
	// ------------------------------------------------------------------
	std::cout << "\nC13. delta_M composition (ADR-0007 P6)\n";
	{
		const int last = static_cast<int>(R0->m0_over_Omega2.Size()) - 1;
		const double Rs = R0->R_surface;
		const double expect = R0->m0_over_Omega2[last] +
							  R0->surface_shell_mass_over_Omega2 +
							  (R0->I * R0->I) / (Rs * Rs * Rs);
		const double rel = Rel(R0->deltaM_over_Omega2, expect);
		Report("C13a deltaM_hat = m0_hat(R_*) + shell_hat + I^2/R_*^3", rel <= kBoundDeltaMIdentity,
			   "deltaM_hat=" + Sci(R0->deltaM_over_Omega2) + "  rel=" + Sci(rel) +
				   "  bound=" + Sci(kBoundDeltaMIdentity));
		Report("C13b the surface shell term is present and non-zero on a constant-density star",
			   R0->surface_shell_mass_over_Omega2 != 0.0 &&
				   std::isfinite(R0->surface_shell_mass_over_Omega2),
			   "shell_hat=" + Sci(R0->surface_shell_mass_over_Omega2) + "  m0_hat(R_*)=" +
				   Sci(R0->m0_over_Omega2[last]) + "  I^2/R_*^3=" +
				   Sci((R0->I * R0->I) / (Rs * Rs * Rs)));
		std::cout << "     (The shell is negligible on an EOS-floor star and DOMINANT here, "
					 "which is why\n      ADR-0007 P6 computes it rather than assuming it "
					 "small.)\n";
	}

	// ------------------------------------------------------------------
	//  C6/C7 — provenance and staleness
	// ------------------------------------------------------------------
	std::cout << "\nC6/C7. Provenance and staleness (ADR-0003)\n";
	{
		const auto &P = ns->Profile();
		const bool match = R0->MatchesSource(static_cast<const void *>(&P), P.Version());
		Report("C6 the response records its source profile identity and Version()", match,
			   "version = " + std::to_string(P.Version()));

		const std::uint64_t v_before = P.Version();
		RotationSolverTestSeam::TouchProfile(*ns); // a bare version bump; data untouched
		const std::uint64_t v_after = ns->Profile().Version();
		const auto *stale = ns->MonopoleResponse();
		Report("C7a a profile version bump makes the cached response stale, and it is NOT "
			   "returned",
			   v_after != v_before && stale == nullptr,
			   "version " + std::to_string(v_before) + " -> " + std::to_string(v_after) +
				   ", MonopoleResponse() == nullptr");

		bool threw = false;
		try
		{
			(void)ns->MonopoleAt(AngularVelocity::FromRadPerSecond(100.0));
		}
		catch (const std::exception &)
		{
			threw = true;
		}
		Report("C7b materializing from a stale response is refused, not silently served", threw,
			   "MonopoleAt threw");

		const std::size_t before = RotationSolverTestSeam::MonopoleSolves(*ns);
		const bool re = ns->ComputeHartleMonopoleResponse();
		const std::size_t after = RotationSolverTestSeam::MonopoleSolves(*ns);
		Report("C7c explicit recomputation restores a current response, with ONE new solve",
			   re && ns->MonopoleResponse() != nullptr && after == before + 1,
			   "solves " + std::to_string(before) + " -> " + std::to_string(after));
	}

	// ------------------------------------------------------------------
	//  C8/C9/C10/C12 — materialization algebra
	// ------------------------------------------------------------------
	std::cout << "\nC8/C9/C10/C12. Materialization at an explicit AngularVelocity\n";
	{
		// Read through a CONST reference deliberately: `NStar::GetSequence()` has a non-const
		// overload that returns `prof_.SeqMutable()`, which calls `Touch()` — so reading the
		// sequence through a mutable star bumps the profile version and, quite correctly,
		// makes every version-keyed cache stale. That is pre-existing ADR-0003 behaviour, not
		// something this increment introduces, and it is exactly what the provenance check is
		// for. It is recorded in the validation record; no first-order code is changed here.
		const NStar &cns = *ns;
		const double I_before = cns.GetSequence().I;
		const std::size_t solves_before = RotationSolverTestSeam::MonopoleSolves(*ns);

		const auto zero = ns->MonopoleAt(AngularVelocity::FromRadPerSecond(0.0));
		bool all_zero = zero.delta_M == 0.0 && zero.surface_shell_mass == 0.0 &&
						zero.surface_xi0 == 0.0 && zero.Omega_geom == 0.0;
		for (int i = 0; all_zero && i < static_cast<int>(zero.m0.Size()); ++i)
			all_zero = zero.m0[i] == 0.0 && zero.p0star[i] == 0.0 &&
					   zero.delta_p0[i] == 0.0 && zero.xi0[i] == 0.0;
		Report("C8 zero spin materializes EXACT zeros in every field", all_zero, "");

		const double w = 2.0 * M_PI * 300.0;
		const auto pp = ns->MonopoleAt(AngularVelocity::FromRadPerSecond(w));
		const auto mm = ns->MonopoleAt(AngularVelocity::FromRadPerSecond(-w));
		bool identical = pp.delta_M == mm.delta_M &&
						 pp.surface_shell_mass == mm.surface_shell_mass &&
						 pp.surface_xi0 == mm.surface_xi0;
		for (int i = 0; identical && i < static_cast<int>(pp.m0.Size()); ++i)
			identical = pp.m0[i] == mm.m0[i] && pp.p0star[i] == mm.p0star[i] &&
						pp.delta_p0[i] == mm.delta_p0[i] && pp.xi0[i] == mm.xi0[i];
		Report("C9 +Omega and -Omega materialize BIT-IDENTICAL perturbations", identical,
			   "Omega = +/-" + Sci(w, 3) + " s^-1; the response is quadratic in Omega");

		const auto p2 = ns->MonopoleAt(AngularVelocity::FromRadPerSecond(2.0 * w));
		double worst = 0.0;
		for (int i = 0; i < static_cast<int>(pp.m0.Size()); ++i)
		{
			worst = std::max(worst, Rel(p2.m0[i], 4.0 * pp.m0[i]));
			worst = std::max(worst, Rel(p2.p0star[i], 4.0 * pp.p0star[i]));
			worst = std::max(worst, Rel(p2.delta_p0[i], 4.0 * pp.delta_p0[i]));
			worst = std::max(worst, Rel(p2.xi0[i], 4.0 * pp.xi0[i]));
		}
		worst = std::max(worst, Rel(p2.delta_M, 4.0 * pp.delta_M));
		Report("C10 doubling Omega quadruples every field: Q(2W) = 4 Q(W)",
			   worst <= kBoundQuadratic,
			   "worst rel=" + Sci(worst) + "  bound=" + Sci(kBoundQuadratic));

		for (int k = 0; k < 25; ++k)
			(void)ns->MonopoleAt(AngularVelocity::FromHz(50.0 + 10.0 * k));
		Report("C11e 25 further materializations perform ZERO additional integrations",
			   RotationSolverTestSeam::MonopoleSolves(*ns) == solves_before,
			   "solves still " + std::to_string(solves_before));

		Report("C12 the first-order moment of inertia is unchanged, bitwise",
			   cns.GetSequence().I == I_before, "I = " + Sci(I_before, 12));
	}

	// ------------------------------------------------------------------
	//  Seed invariance (§25) — the hard requirement of ADR-0006 P9
	// ------------------------------------------------------------------
	std::cout << "\nS. Seed invariance of the O(Omega^2) coefficients (ADR-0006 P9, "
				 "ADR-0007 §7-1)\n";
	{
		auto ref = BuildSeeded(u, N, RotationSolverTestSeam::DefaultSeed());
		const bool ref_ok = ref->ComputeHartleMonopoleResponse();
		const auto *A = ref->MonopoleResponse();

		bool all_ok = ref_ok && A != nullptr;
		double worst_all = 0.0;
		std::string worst_tag;

		for (double seed : {1.0e-3, 1.0e-1, 1.0e1, 1.0e3})
		{
			auto other = BuildSeeded(u, N, seed);
			const bool ok = other->ComputeHartleMonopoleResponse();
			const auto *B = other->MonopoleResponse();
			if (!ok || B == nullptr || A == nullptr)
			{
				all_ok = false;
				continue;
			}
			double worst = Rel(B->deltaM_over_Omega2, A->deltaM_over_Omega2);
			worst = std::max(worst, Rel(B->surface_xi0_over_Omega2, A->surface_xi0_over_Omega2));
			worst = std::max(worst, Rel(B->surface_shell_mass_over_Omega2,
										A->surface_shell_mass_over_Omega2));
			for (int i = 0; i < static_cast<int>(A->m0_over_Omega2.Size()); ++i)
			{
				worst = std::max(worst, Rel(B->m0_over_Omega2[i], A->m0_over_Omega2[i]));
				worst = std::max(worst, Rel(B->p0star_over_Omega2[i], A->p0star_over_Omega2[i]));
				worst = std::max(worst, Rel(B->xi0_over_Omega2[i], A->xi0_over_Omega2[i]));
			}
			if (worst > worst_all)
			{
				worst_all = worst;
				worst_tag = Sci(seed, 1);
			}
			std::cout << "     seed " << Sci(seed, 1) << " : worst rel = " << Sci(worst) << "\n";
		}
		Report("Sa every monopole coefficient is invariant under the internal first-order seed "
			   "over [1e-3, 1e3]",
			   all_ok && worst_all <= kBoundSeedInvariance,
			   "worst rel=" + Sci(worst_all) + " (seed " + worst_tag + ")  bound=" +
				   Sci(kBoundSeedInvariance));
		std::cout << "     (The solver consumes only omega_bar/Omega and its derivative. The raw\n"
					 "      stored_omega_bar_ arrays are unreachable from it by construction.)\n";
	}

	// ------------------------------------------------------------------
	//  Diagnostics — reported, never asserted as physics
	// ------------------------------------------------------------------
	std::cout << "\nDiagnostics (analytic star, M = " << Sci(u.M_km, 2)
			  << " km, R = " << Sci(u.R_km, 2) << " km) — NOT scientific claims:\n";
	{
		const auto *Rr = ns->MonopoleResponse();
		if (Rr != nullptr)
		{
			const int last = static_cast<int>(Rr->m0_over_Omega2.Size()) - 1;
			std::cout << "     R_*                 = " << Sci(Rr->R_surface) << " km"
					  << "   (production surface node, NOT the exact p = 0 radius)\n"
					  << "     I                   = " << Sci(Rr->I) << " km^3\n"
					  << "     m0_hat(R_*)         = " << Sci(Rr->m0_over_Omega2[last])
					  << " km^3\n"
					  << "     p0*_hat(R_*)        = " << Sci(Rr->p0star_over_Omega2[last])
					  << " km^2\n"
					  << "     xi0_hat(R_*)        = " << Sci(Rr->surface_xi0_over_Omega2)
					  << " km^3\n"
					  << "     shell_hat           = "
					  << Sci(Rr->surface_shell_mass_over_Omega2) << " km^3\n"
					  << "     deltaM_hat          = " << Sci(Rr->deltaM_over_Omega2)
					  << " km^3\n";
		}
	}

	std::cout << "\nNON-CLAIM: nothing above verifies the O(Omega^2) PHYSICS. These are\n"
				 "implementation contracts. Independent scientific validation is Phase 4D;\n"
				 "until then INV-08 remains 'source conformed, validation pending' and no\n"
				 "monopole baseline may be created.\n";

	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail
			  << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
