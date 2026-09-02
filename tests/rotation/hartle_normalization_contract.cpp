// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_normalization_contract.cpp
 * @brief Phase 4A — conformance tests for the ADR-0006 first-order physical-normalization and
 *        unit contract. Self-contained: no external EOS assets.
 *
 * ADR-0006 (ACCEPTED 2026-09-02) requires these checks, and fixes every tolerance below
 * **before** the implementation existed (ADR-0006 §7). They are reproduced here verbatim:
 *
 *   V1 internal-seed invariance over [1e-3, 1e3]        rel <= 1e-10
 *   V2 requested Omega recovered at the surface         rel <= 1e-13
 *   V3 J = I * Omega_phys / c                           rel <= 1e-13
 *   V4 I bit-identical (also: the golden artifact)      bitwise
 *   V5 omega_bar_phys linear in requested Omega         rel <= 1e-14
 *   V6 Omega_geom * c = Omega_phys, independent oracle  rel <= 1e-15
 *   V7 annotations and exported header tokens agree     exact
 *   V8 zero spin well defined, no 0/0                   exact
 *
 * WHAT IS CLAIMED. That the *contract* holds: the arbitrary seed does not reach any public
 * result, a requested physical spin is the spin that comes back, the units are what they say,
 * and `I` is untouched. It does **not** revalidate the first-order physics — that is Phase
 * 2B-4B (`docs/validation/HARTLE_MOMENT_INERTIA.md`, `EQUATION MATCH`, `HARTLE-I VERIFIED`) —
 * and it makes **no** claim about O(Omega^2), which remains an unverified candidate (INV-08).
 *
 * INDEPENDENT UNIT ORACLE. V6 uses a literal `c = 299792.458 km/s`, exact by the SI definition
 * of the metre and the second. It deliberately does **not** import
 * `Zaki::Physics::LIGHT_C_KM_S`, which is the constant production converts with: a test that
 * imported it could not detect a wrong conversion, only a wrong constant.
 *
 * BACKGROUND. The exact Schwarzschild constant-density interior, in geometric units (km) — the
 * same fixture Phase 2B-4B used, restated here so that the validated 2B-4B harnesses are left
 * byte-untouched (ADR-0006 §7 item 4):
 *
 *   rho0 = 3M/(4 pi R^3),   m(r) = M r^3/R^3,   y(r) = sqrt(1 - 2 M r^2/R^3)
 *   p(r) = rho0 (y - y_R)/(3 y_R - y),   nu'(r) = 2 M r / ( R^3 y (3 y_R - y) )
 *
 * (Schwarzschild 1916; MTW §23.7 / Shapiro & Teukolsky §5.5.) The fixture is synthetic and
 * incompressible; it asserts no neutron-star property. It is used because the contract is a
 * statement about scaling and units, which any admissible background exercises equally.
 */

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
//  Non-public test seam (ADR-0006 Q2, task §14).
//
//  The arbitrary omega_bar seed must NOT become public scientific API, yet seed invariance has
//  to be *proved* rather than asserted. `RotationSolverTestSeam` is declared — never defined —
//  in `RotationSolver.hpp` and befriended by `RotationSolver` and `NStar`; defining it here
//  gives this harness the access it needs.
//
//  CLASSIFICATION (Phase 4B, precise): this is a
//
//      PRIVILEGED TEST BACKDOOR — NOT SUPPORTED SCIENTIFIC API.
//
//  It is *not* a mechanism only a test can use: a friend declaration names a type, and any
//  translation unit could define that type and obtain the same access. What ADR-0006 Q2
//  actually requires, and what holds, is narrower and sufficient: there is no supported public
//  seed setter, no supported public seed constructor argument, and no production consumer of
//  the seam.
// ---------------------------------------------------------------------------
namespace CompactStar::Core
{
struct RotationSolverTestSeam
{
	static void SetSeed(RotationSolver &rs, double seed) { rs.seed_omega_bar_ = seed; }
	static double Seed(const RotationSolver &rs) { return rs.seed_omega_bar_; }
	static double DefaultSeed() { return RotationSolver::kDefaultInitOmegaBar; }
	static std::size_t Solves(const RotationSolver &rs) { return rs.first_order_solve_count_; }
	static std::size_t Solves(const NStar &ns) { return ns.rot_solver.first_order_solve_count_; }

	/// Test-only: the legacy sequence container is populated in production ONLY by
	/// `Solve_Mixed`, and `VecSaver::Export1D` returns immediately on an empty vector — so an
	/// ordinary NStar exports no file at all. Pushing one row here lets the corrected header be
	/// parsed against the values it labels, without giving production a new writer.
	static void PushSeqPoint(RotationSolver &rs, const OmegaSeqPoint &pt)
	{
		rs.omega_seq_pts.push_back(pt);
	}
	static std::size_t SeqPointCount(const RotationSolver &rs)
	{
		return rs.omega_seq_pts.size();
	}
};
} // namespace CompactStar::Core

using CompactStar::AngularVelocity;
using CompactStar::Core::HartleFirstOrderResponse;
using CompactStar::Core::NStar;
using CompactStar::Core::PhysicalFirstOrderRotation;
using CompactStar::Core::RotationSolver;
using CompactStar::Core::RotationSolverTestSeam;
using CompactStar::Core::TOVPoint;

// The independent oracle: exact by the SI definitions of the metre and the second.
static constexpr double kIndependentLightSpeed_km_s = 299792.458;

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

// ---------------------------------------------------------------------------
//  Exact constant-density interior, tabulated in the geometric units StarProfile uses.
// ---------------------------------------------------------------------------
struct Uniform
{
	double M_km = 0, R_km = 0, rho0 = 0; // rho0 in km^-2
	double yR() const { return std::sqrt(1.0 - 2.0 * M_km / R_km); }
	double y(double r) const { return std::sqrt(1.0 - 2.0 * M_km * r * r / (R_km * R_km * R_km)); }
	double m(double r) const { return M_km * r * r * r / (R_km * R_km * R_km); }
	double p(double r) const { return rho0 * (y(r) - yR()) / (3.0 * yR() - y(r)); }
	double nuprime(double r) const
	{
		return 2.0 * M_km * r / (R_km * R_km * R_km * y(r) * (3.0 * yR() - y(r)));
	}
};

static Uniform MakeUniform(double M_km, double R_km)
{
	Uniform u;
	u.M_km = M_km;
	u.R_km = R_km;
	u.rho0 = 3.0 * M_km / (4.0 * M_PI * R_km * R_km * R_km);
	return u;
}

/// Production TOV profiles begin at r_min = 1e-5 km, never at r = 0; the analytic tabulation
/// mirrors that so the star is an admissible production profile.
static constexpr double kRStart_km = 1.0e-5;

static double GridR(const Uniform &u, std::size_t i, std::size_t N)
{
	return kRStart_km + static_cast<double>(i) * (u.R_km - kRStart_km) / static_cast<double>(N - 1);
}

/// Build the analytic star as a production `NStar`, through the public TOVPoint constructor.
/// TOVPoint units: r [km], m [Msun], nu_der [1/cm], p [dyne/cm^2], e [g/cm^3].
static std::unique_ptr<NStar> ProductionStar(const Uniform &u, std::size_t N)
{
	const double km2_to_gcm3 =
		Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;
	const double km2_to_dyn =
		Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_INV_KM2;

	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = GridR(u, i, N);
		pts.emplace_back(r, u.m(r) / Zaki::Physics::SUN_M_KM, u.nuprime(r) / 1.0e5, 0.0,
						 u.p(r) * km2_to_dyn, u.rho0 * km2_to_gcm3, 0.1, std::vector<double>{});
	}
	return std::make_unique<NStar>(pts);
}

/// Surface angular velocity recomputed from a physical solution by exterior matching:
/// `Omega = omega_bar(R) + R omega_bar'(R)/3`. This is the *independent* recovery of the
/// requested spin — it reads the published arrays rather than the stored scalar.
static double SurfaceOmegaFrom(const PhysicalFirstOrderRotation &ph)
{
	const int n = static_cast<int>(ph.omega_bar.Size());
	const double R = (*ph.r_grid)[n - 1];
	return ph.omega_bar[n - 1] + R * ph.domega_bar[n - 1] / 3.0;
}

int main()
{
	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 4A — ADR-0006 first-order physical-normalization contract\n"
				 "Self-contained: exact Schwarzschild constant-density interior.\n"
				 "Validates the CONTRACT (normalization, units, seed isolation), not the\n"
				 "first-order physics (Phase 2B-4B) and nothing at O(Omega^2) (INV-08).\n\n";

	const Uniform u = MakeUniform(/*M_km=*/2.0, /*R_km=*/13.0); // M/R ~ 0.154
	auto ns = ProductionStar(u, 4001);
	if (!ns)
	{
		std::cerr << "fixture star construction failed\n";
		return 1;
	}

	// =======================================================================
	//  C0 — the AngularVelocity value type (ADR-0006 Q1 = A + D)
	// =======================================================================
	std::cout << "C0 typed physical spin input\n";
	{
		Report("C0a default construction is exactly zero spin",
			   AngularVelocity{}.RadPerSecond() == 0.0 && AngularVelocity{}.IsZero(),
			   "Omega = " + Sci(AngularVelocity{}.RadPerSecond()));

		const double w = 4.5e3;
		Report("C0a2 FromRadPerSecond round-trips exactly",
			   AngularVelocity::FromRadPerSecond(w).RadPerSecond() == w, "Omega = " + Sci(w));

		// Omega = 2 pi f, applied in exactly one place.
		const double f = 716.0;
		Report("C0b FromHz applies Omega = 2 pi f exactly",
			   AngularVelocity::FromHz(f).RadPerSecond() == 2.0 * M_PI * f,
			   "f = " + Sci(f) + " Hz -> " + Sci(AngularVelocity::FromHz(f).RadPerSecond()) +
				   " rad/s");

		// Omega = 2 pi / P, and consistent with FromHz(1/P).
		const double P = 1.0 / f;
		Report("C0c FromPeriodSeconds agrees with FromHz(1/P)",
			   Rel(AngularVelocity::FromPeriodSeconds(P).RadPerSecond(),
				   AngularVelocity::FromHz(f).RadPerSecond()) < 1e-15,
			   "P = " + Sci(P) + " s");

		// Negative spin is valid: MagneticDipole writes -K sign(Omega) |Omega|^n.
		Report("C0d negative angular velocity is accepted (sign convention not restricted)",
			   AngularVelocity::FromRadPerSecond(-w).RadPerSecond() == -w, "Omega = " + Sci(-w));

		auto throws = [](auto &&fn) {
			try
			{
				fn();
				return false;
			}
			catch (const std::exception &)
			{
				return true;
			}
		};
		const double nan_v = std::numeric_limits<double>::quiet_NaN();
		const double inf_v = std::numeric_limits<double>::infinity();
		Report("C0e non-finite inputs fail closed",
			   throws([&] { AngularVelocity::FromRadPerSecond(nan_v); }) &&
				   throws([&] { AngularVelocity::FromRadPerSecond(inf_v); }) &&
				   throws([&] { AngularVelocity::FromHz(nan_v); }) &&
				   throws([&] { AngularVelocity::FromPeriodSeconds(inf_v); }),
			   "NaN and infinity rejected by every factory");
		Report("C0f a zero PERIOD fails closed, while a zero frequency does not",
			   throws([&] { AngularVelocity::FromPeriodSeconds(0.0); }) &&
				   !throws([&] { AngularVelocity::FromHz(0.0); }),
			   "a zero period is undefined; a zero frequency is the non-rotating limit");
	}

	// =======================================================================
	//  V6 — the conversion, against an INDEPENDENT c (must precede everything that uses it)
	// =======================================================================
	std::cout << "\nV6 physical <-> geometric conversion vs an independent c = 299792.458 km/s\n";
	{
		double worst = 0.0;
		for (double w : {1.0e2, 2.0 * M_PI * 1.0e2, 2.0 * M_PI * 716.0, -3.3e2})
		{
			const auto omega = AngularVelocity::FromRadPerSecond(w);
			worst = std::max(worst,
							 Rel(omega.GeomKmInverse() * kIndependentLightSpeed_km_s, w));
		}
		Report("V6a Omega_geom * c_independent reproduces the requested rad/s", worst <= 1e-15,
			   "worst relative difference " + Sci(worst) + " (bound 1e-15)");

		// The inverse rendering must be the exact inverse of the forward conversion.
		const double wg = AngularVelocity::FromRadPerSecond(4.5e3).GeomKmInverse();
		Report("V6b the geometric -> physical rendering is the inverse of the forward conversion",
			   Rel(CompactStar::AngularVelocityGeomToRadPerSecond(wg) /
					   kIndependentLightSpeed_km_s,
				   wg) <= 1e-15,
			   "round trip through both conversion owners");
	}

	// =======================================================================
	//  C1 — construction yields a seed-free response and NO implicit spin
	// =======================================================================
	std::cout << "\nC1 ordinary construction: seed-free response, no implicit physical spin\n";
	const auto &resp = ns->RotationResponse();
	{
		Report("C1a a built star has a valid first-order response", resp.valid,
			   "response.valid = " + std::string(resp.valid ? "true" : "false"));

		Report("C1b construction performs exactly ONE first-order ODE solve",
			   RotationSolverTestSeam::Solves(*ns) == 1,
			   "solves = " + std::to_string(RotationSolverTestSeam::Solves(*ns)));

		// V4, locally: the governed observable is untouched by this increment.
		Report("V4a response.I is BIT-IDENTICAL to SeqPoint::I", resp.I == ns->GetSequence().I,
			   "I = " + Sci(resp.I, 14) + " km^3");

		// The normalization identity: the response is a ratio, so its own surface combination
		// is exactly the dimensionless 1. This is simultaneously a unit check (V7) — a
		// dimensionless quantity must come out as a pure number.
		const int n = static_cast<int>(resp.omega_bar_over_Omega.Size());
		const double R = (*resp.r_grid)[n - 1];
		const double unit_omega =
			resp.omega_bar_over_Omega[n - 1] + R * resp.domega_bar_over_Omega_dr[n - 1] / 3.0;
		Report("C1c the response is normalized: [omega_bar/Omega](R) + R [omega_bar'/Omega](R)/3 = 1",
			   Rel(unit_omega, 1.0) <= 1e-13,
			   "value = " + Sci(unit_omega, 14) + ", i.e. dimensionless as documented");

		// The frame-dragging ratio is a property of the star, available with no spin at all.
		const double drag_c = 1.0 - resp.omega_bar_over_Omega[0];
		Report("C1d the frame-dragging ratio omega(0)/Omega is available without any spin",
			   drag_c > 0.0 && drag_c < 1.0, "omega_c/Omega = " + Sci(drag_c));
	}

	// =======================================================================
	//  V1 — internal-seed invariance (the central ADR-0006 Q2 claim)
	// =======================================================================
	std::cout << "\nV1 internal-seed invariance over six decades, at a FIXED requested Omega\n";
	{
		const double w_req = 2.0 * M_PI * 716.0; // the fastest known pulsar, as a stress case
		const auto omega = AngularVelocity::FromRadPerSecond(w_req);
		const std::vector<double> seeds = {1.0e-3, 5.0e-3, 1.0e-1, 1.0, 1.0e1, 1.0e3};

		std::vector<double> Is, Js, Ws, OBs;
		std::cout << "      seed        I [km^3]        Omega_geom      J [km^2]        "
					 "omega_bar(R)\n";
		for (double s : seeds)
		{
			RotationSolver rs;
			rs.AttachNStar(ns.get());
			RotationSolverTestSeam::SetSeed(rs, s);
			rs.FindNMomInertia();
			const auto ph = rs.FirstOrderResponse().At(omega);
			const int n = static_cast<int>(ph.omega_bar.Size());
			Is.push_back(rs.FirstOrderResponse().I);
			Js.push_back(ph.J);
			Ws.push_back(ph.Omega_geom);
			OBs.push_back(ph.omega_bar[n - 1]);
			std::cout << "   " << Sci(s, 1) << "  " << Sci(Is.back(), 8) << "  "
					  << Sci(Ws.back(), 8) << "  " << Sci(Js.back(), 8) << "  "
					  << Sci(OBs.back(), 8) << "\n";
		}

		auto spread = [](const std::vector<double> &v) {
			double w = 0.0;
			for (std::size_t i = 1; i < v.size(); ++i)
				w = std::max(w, Rel(v[i], v[0]));
			return w;
		};
		const double sI = spread(Is), sJ = spread(Js), sW = spread(Ws), sOB = spread(OBs);
		Report("V1a Omega is invariant to the internal seed", sW <= 1e-10,
			   "relative spread " + Sci(sW) + " (bound 1e-10)");
		Report("V1b J_phys is invariant to the internal seed", sJ <= 1e-10,
			   "relative spread " + Sci(sJ) + " (bound 1e-10)");
		Report("V1c omega_bar_phys(R) is invariant to the internal seed", sOB <= 1e-10,
			   "relative spread " + Sci(sOB) + " (bound 1e-10)");
		Report("V1d I is invariant to the internal seed", sI <= 1e-10,
			   "relative spread " + Sci(sI) + " (bound 1e-10)");

		// The production seed is unchanged by this increment, and is not public API.
		Report("V1e the production default seed is unchanged",
			   RotationSolverTestSeam::DefaultSeed() == 5e-3,
			   "default internal seed = " + Sci(RotationSolverTestSeam::DefaultSeed()));
	}

	// =======================================================================
	//  V2 / V3 / V5 / V8 — the materialized physical solution
	// =======================================================================
	std::cout << "\nV2/V3/V5 requested-Omega recovery, J = I Omega, and linearity\n";
	{
		const std::vector<double> ws = {1.0e2, 2.0 * M_PI * 1.0e2, 2.0 * M_PI * 716.0, -3.3e2};
		double worst_recov = 0.0, worst_J = 0.0, worst_rad = 0.0;

		std::cout << "      Omega_phys [1/s]  Omega_geom [1/km]  J [km^2]        "
					 "surface recovery  J vs I*Omega\n";
		for (double w : ws)
		{
			const auto ph = ns->RotationAt(AngularVelocity::FromRadPerSecond(w));
			const double recov = Rel(SurfaceOmegaFrom(ph), ph.Omega_geom);
			const double dJ = Rel(ph.J, resp.I * w / kIndependentLightSpeed_km_s);
			worst_recov = std::max(worst_recov, recov);
			worst_J = std::max(worst_J, dJ);
			worst_rad = std::max(worst_rad, Rel(ph.OmegaRadPerSecond(), w));
			std::cout << "   " << Sci(w, 8) << "  " << Sci(ph.Omega_geom, 8) << "  "
					  << Sci(ph.J, 8) << "  " << Sci(recov) << "  " << Sci(dJ) << "\n";
		}
		Report("V2a the requested Omega is recovered from the surface of the physical solution",
			   worst_recov <= 1e-13,
			   "worst relative difference " + Sci(worst_recov) + " (bound 1e-13)");
		Report("V2b OmegaRadPerSecond() returns the requested rad/s", worst_rad <= 1e-13,
			   "worst relative difference " + Sci(worst_rad) + " (bound 1e-13)");
		Report("V3a J = I * Omega_phys / c, with c from the INDEPENDENT oracle", worst_J <= 1e-13,
			   "worst relative difference " + Sci(worst_J) + " (bound 1e-13)");

		// V5 — linearity node by node.
		const double w0 = 2.0 * M_PI * 300.0;
		const auto p1 = ns->RotationAt(AngularVelocity::FromRadPerSecond(w0));
		const auto p2 = ns->RotationAt(AngularVelocity::FromRadPerSecond(2.0 * w0));
		const auto p3 = ns->RotationAt(AngularVelocity::FromRadPerSecond(-0.5 * w0));
		const int n = static_cast<int>(p1.omega_bar.Size());
		double worst_lin = 0.0, worst_lin_d = 0.0;
		for (int i = 0; i < n; ++i)
		{
			worst_lin = std::max(worst_lin, Rel(p2.omega_bar[i], 2.0 * p1.omega_bar[i]));
			worst_lin = std::max(worst_lin, Rel(p3.omega_bar[i], -0.5 * p1.omega_bar[i]));
			worst_lin_d = std::max(worst_lin_d, Rel(p2.domega_bar[i], 2.0 * p1.domega_bar[i]));
		}
		Report("V5a omega_bar_phys(r) scales linearly with the requested Omega at EVERY node",
			   worst_lin <= 1e-14 && n > 100,
			   "worst node relative difference " + Sci(worst_lin) + " over " +
				   std::to_string(n) + " nodes (bound 1e-14)");
		Report("V5b d(omega_bar)/dr scales with it — the derivative is not left behind",
			   worst_lin_d <= 1e-14,
			   "worst node relative difference " + Sci(worst_lin_d) + " (bound 1e-14)");
	}

	std::cout << "\nV8 zero spin\n";
	{
		const auto z = ns->RotationAt(AngularVelocity{});
		const int n = static_cast<int>(z.omega_bar.Size());
		double max_ob = 0.0, max_dob = 0.0;
		for (int i = 0; i < n; ++i)
		{
			max_ob = std::max(max_ob, std::fabs(z.omega_bar[i]));
			max_dob = std::max(max_dob, std::fabs(z.domega_bar[i]));
		}
		Report("V8a Omega, J and omega_bar are EXACTLY zero at zero spin",
			   z.Omega_geom == 0.0 && z.J == 0.0 && max_ob == 0.0 && max_dob == 0.0,
			   "Omega = " + Sci(z.Omega_geom) + ", J = " + Sci(z.J) +
				   ", max|omega_bar| = " + Sci(max_ob));
		Report("V8b I survives zero spin unchanged and finite — it is never formed as J/Omega",
			   z.I == resp.I && std::isfinite(z.I) && z.I > 0.0,
			   "I = " + Sci(z.I, 14) + " km^3, bit-identical to the scale-free value");
		Report("V8c a zero-spin solution is still valid", z.valid, "valid = true");
	}

	// =======================================================================
	//  §24 — a physical solution is a scaling, never a new integration
	// =======================================================================
	std::cout << "\nP1 performance contract: materialization does not re-integrate\n";
	{
		RotationSolver rs;
		rs.AttachNStar(ns.get());
		rs.FindNMomInertia();
		const std::size_t after_solve = RotationSolverTestSeam::Solves(rs);
		for (int k = 1; k <= 25; ++k)
			(void)rs.FirstOrderResponse().At(AngularVelocity::FromHz(10.0 * k));
		const std::size_t after_scalings = RotationSolverTestSeam::Solves(rs);
		Report("P1a 25 materializations at different Omega perform ZERO extra ODE solves",
			   after_solve == 1 && after_scalings == 1,
			   "solve count " + std::to_string(after_solve) + " -> " +
				   std::to_string(after_scalings));

		// The star's own response object is stable: NStar caches, it does not recompute.
		const auto *a = &ns->RotationResponse();
		(void)ns->RotationAt(AngularVelocity::FromHz(500.0));
		const auto *b = &ns->RotationResponse();
		Report("P1b NStar::RotationResponse() is a stable cached object, not a recomputation",
			   a == b && RotationSolverTestSeam::Solves(*ns) == 1,
			   "same address, star solve count still " +
				   std::to_string(RotationSolverTestSeam::Solves(*ns)));
	}

	// =======================================================================
	//  V7 — annotations and exported header tokens describe the stored values
	// =======================================================================
	std::cout << "\nV7 unit annotations and export header\n";
	{
		const double w = 2.0 * M_PI * 400.0;
		const auto ph = ns->RotationAt(AngularVelocity::FromRadPerSecond(w));
		const int n = static_cast<int>(ph.omega_bar.Size());
		const double R = (*ph.r_grid)[n - 1];

		// Omega_geom is documented [km^-1]: it must equal Omega_phys/c with the independent c.
		Report("V7a Omega_geom really is km^-1",
			   Rel(ph.Omega_geom, w / kIndependentLightSpeed_km_s) <= 1e-15,
			   "stored " + Sci(ph.Omega_geom) + " vs Omega_phys/c " +
				   Sci(w / kIndependentLightSpeed_km_s));

		// omega_bar is documented [km^-1] in the SAME system as Omega_geom. Two statements are
		// needed: the surface combination reproduces the geometric Omega, and multiplying it by
		// c — and only by c — turns it into the requested rad/s. The ratio, not a relative
		// difference, is what separates the two unit systems.
		const double surf = SurfaceOmegaFrom(ph);
		Report("V7b omega_bar is km^-1, in the same system as Omega_geom",
			   Rel(surf, ph.Omega_geom) <= 1e-13,
			   "surface combination " + Sci(surf) + " vs stored Omega_geom " +
				   Sci(ph.Omega_geom));
		Report("V7b2 omega_bar is NOT in s^-1: it is exactly c times smaller than the request",
			   Rel(surf * kIndependentLightSpeed_km_s, w) <= 1e-13 &&
				   Rel(w / surf, kIndependentLightSpeed_km_s) <= 1e-13,
			   "requested/surface = " + Sci(w / surf) + " = c, so the stored profile is "
			   "geometric as annotated");

		// J is documented [km^2]: I [km^3] * Omega [km^-1] is km^2.
		Report("V7c J is km^2, consistent with I [km^3] x Omega [km^-1]",
			   Rel(ph.J, ph.I * ph.Omega_geom) <= 1e-13,
			   "J = " + Sci(ph.J) + " km^2");

		// The response's shape column is documented dimensionless.
		Report("V7d omega_bar_over_Omega is dimensionless (unchanged by the requested Omega)",
			   Rel(ph.omega_bar[n / 2] / ph.Omega_geom, resp.omega_bar_over_Omega[n / 2]) <= 1e-13,
			   "the physical profile divided by Omega returns the stored shape");
		(void)R;

		// ---- the legacy sequence export header ----------------------------------------
		// First, the separation ADR-0006 P7 requires: the ordinary NStar first-order path does
		// not export through here at all. `omega_seq_pts` is populated only by Solve_Mixed, and
		// VecSaver::Export1D returns immediately on an empty vector, so no file is produced.
		const fs::path wrk = fs::temp_directory_path() / "compactstar_phase4a_export";
		fs::remove_all(wrk);
		fs::create_directories(wrk);
		{
			RotationSolver rs_empty;
			rs_empty.SetWrkDir(wrk.string());
			rs_empty.ExportResults("nstar_path_probe.tsv");
			Report("V7e1 the legacy sequence export is inert on the ordinary NStar path",
				   RotationSolverTestSeam::SeqPointCount(rs_empty) == 0 &&
					   !fs::exists(wrk / "nstar_path_probe.tsv"),
				   "no sequence rows and no file: physical NStar results come from "
				   "NStar::RotationAt(), not from this legacy export");
		}

		// Now parse the corrected header against a row whose values are known, so that each
		// token is checked against the quantity it actually labels.
		RotationSolver rs;
		rs.SetWrkDir(wrk.string());
		const double v_seed = 1.5e3, v_M = 1.4, v_R = 13.5, v_J = 2.5, v_Omega = 4.5e3;
		RotationSolverTestSeam::PushSeqPoint(
			rs, CompactStar::Core::OmegaSeqPoint(v_seed, v_M, v_R, v_J, v_Omega));
		rs.ExportResults("rotation_seq_header_probe.tsv");

		std::ifstream in((wrk / "rotation_seq_header_probe.tsv").string());
		std::string header, row;
		std::getline(in, header);
		std::getline(in, row);
		in.close();

		auto has = [&header](const std::string &tok) {
			return header.find(tok) != std::string::npos;
		};
		const bool tokens_ok = has("omega_bar_c_seed (1/s)") && has("M (M_sun)") &&
							   has("R (km)") && has("J_seednorm (km^2)") &&
							   has("Omega_seednorm (1/s)");
		Report("V7e2 every exported header token carries the unit of the value written",
			   tokens_ok, "header: \"" + header + "\"");

		// The tokens must be in the order the values are written, or the units label the wrong
		// columns. Parse the row and match it against what was pushed.
		std::vector<double> vals;
		{
			std::string r = row;
			for (char &ch : r)
				if (ch == '\t')
					ch = ' ';
			std::istringstream ss(r);
			double v;
			while (ss >> v)
				vals.push_back(v);
		}
		const bool order_ok = vals.size() == 5 && Rel(vals[0], v_seed) < 1e-7 &&
							  Rel(vals[1], v_M) < 1e-7 && Rel(vals[2], v_R) < 1e-7 &&
							  Rel(vals[3], v_J) < 1e-7 && Rel(vals[4], v_Omega) < 1e-7;
		Report("V7e3 the header tokens are in the order of the values they label", order_ok,
			   "parsed " + std::to_string(vals.size()) + " values matching the written row");

		Report("V7f the export no longer advertises a seed-normalized value as a physical spin",
			   header.find("omega_bar_c (1/s)") == std::string::npos &&
				   header.find("Omega (1/s)") == std::string::npos,
			   "the pre-ADR-0006 tokens 'omega_bar_c (1/s)' and 'Omega (1/s)' are gone; the "
			   "values are labelled as the seed-normalized legacy quantities they are");
		fs::remove_all(wrk);
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "ADR-0006 first-order normalization contract checks passed"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
