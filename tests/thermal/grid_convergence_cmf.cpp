// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file grid_convergence_cmf.cpp
 * @brief Phase 2B-4A — measure how the NONROTATING production pipeline
 *
 *        TOV structure -> StarProfile/GeometryCache -> C_star, L_nu, L_gamma
 *                      -> passive thermal evolution -> T_inf(t)
 *
 *        converges under refinement of the TOV radial resolution.
 *
 *   usage: grid_convergence_cmf <EOS_DATA_ROOT> [--emit-dir <dir>]
 *
 * SCOPE BOUNDARY (binding). This is the NONROTATING SLICE only. It does not exercise the
 * Hartle observable, does not claim Hartle convergence, does not resolve INV-07, does not
 * touch `init_omega_bar`, and does not complete the ratified TOV -> Hartle -> cooling
 * roadmap item. The moment of inertia is never reported as a convergence observable here.
 *
 * NOTE ON RotationSolver. `NStar::BuildFromTOV` calls `Find_MomInertia()` ->
 * `RotationSolver::FindNMomInertia()` as a side effect of constructing ANY StarProfile.
 * That is pre-existing production behavior on the canonical path, and the established
 * passive-cooling baseline already incurs it. It cannot be avoided without deviating from
 * the production construction this study exists to measure. It is inert for the observable
 * studied here: `couple_spin = false`, no spin driver is registered, and `I` never enters
 * `T_inf(t)`. No statement about its accuracy or convergence is made or implied.
 *
 * INV-13 governs the expectation: interpolation is linear and `DataSet::Integrate` is the
 * trapezoid rule, nominally O(dr^2) for smooth integrands — but "the convergence order of
 * the complete calculation must be measured, not inferred from the interpolation scheme."
 * Nothing here assumes second order.
 *
 * TWO EXPERIMENTS, deliberately not conflated:
 *   A  fixed central energy density  -> isolates radial discretization of ONE physical star
 *   B  fixed target mass 1.6 Msun    -> the actual production workflow, mass search included
 *
 * The scientific setup (initial state, interval, epochs, drivers, envelope, stepper,
 * tolerances, cadence) is copied from `passive_cooling_regression.cpp` so this is the same
 * cooling problem, not a subtly different one.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Driver/Thermal/Boundary/EnvelopePotekhin1997.hpp"
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp"
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling_Details.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling_Details.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp"
#include "CompactStar/Physics/Evolution/Observers/IObserver.hpp"
#include "CompactStar/Physics/Evolution/Run/RunBuilder.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StatePacking.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;
namespace P = CompactStar::Physics;
namespace Th = CompactStar::Physics::Driver::Thermal;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;
using CompactStar::EOS::CompOSE_Thermo;

// ---------------------------------------------------------------------------
//  Frozen scientific configuration — identical to passive_cooling_regression.cpp
// ---------------------------------------------------------------------------
constexpr double kTargetM = 1.6;   // Msun
constexpr double kTinf0_K = 1.0e9; // K
constexpr double kT0_yr = 1.0e2;
constexpr double kT1_yr = 1.0e6;
constexpr double kRtol = 1e-6;
constexpr double kAtol = 1e-10;
constexpr double kSamplesPerDecade = 150.0;
constexpr P::Evolution::StepperType kStepper = P::Evolution::StepperType::RKF45;
static const std::vector<double> kEpochsYr = {
	1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6};

/// Static-diagnostic temperature for Experiment A.
constexpr double kTstatic_K = 1.0e8;

/// Refinement sequence. Contains the production default (10000).
static const std::vector<std::size_t> kResolutions = {2500, 5000, 10000, 20000, 40000};

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
/// Scientific formatting for report strings — std::to_string would print 1e-7 as 0.000000.
static std::string Sci(double v, int prec = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
}

// ---------------------------------------------------------------------------
//  Per-resolution record
// ---------------------------------------------------------------------------
struct Case
{
	std::size_t res = 0;
	int n_profile = 0;
	double dr_min = 0, dr_max = 0, dr_med = 0, dr_eff = 0; // km
	double ec = 0, pc = 0, M = 0, R = 0, B = 0;
	double z_surf = 0;	   // e^{nu(R)}
	double compactness = 0; // 2GM/(Rc^2)
	// static thermal diagnostics at kTstatic_K
	double C_star = 0, L_nu = 0, L_nu_DU = 0, L_nu_MU = 0, L_gamma = 0, dlnT_dt = 0;
	// trajectory (Experiment B only)
	std::vector<double> Tinf, Lnu_t, Lgam_t, Cstar_t;
	bool ok = false;
	double seconds = 0;
};

// ---------------------------------------------------------------------------
//  Build a star at a given radial resolution.
//
//  This is the ONLY structural deviation from NStar::SolveTOV_Profile: the internal
//  TOVSolver is constructed here so SetRadialRes() can be called on it. Everything
//  downstream — SolveToProfile and BuildFromTOV, via the public NStar(points, labels)
//  constructor — is the identical production code. Section EQ below proves the
//  res = 10000 case reproduces the canonical path.
// ---------------------------------------------------------------------------
static std::vector<std::string> g_species_labels; // captured once from the EOS table

static std::unique_ptr<NStar> BuildAtResolution(const fs::path &cold, const fs::path &wrk,
												std::size_t res, bool fixed_ec,
												double ec_fixed, Case &c)
{
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	tov.SetRadialRes(res); // the only knob varied in this study

	std::vector<TOVPoint> pts;
	std::vector<std::string> labels;
	int n = 0;
	if (fixed_ec)
	{
		n = tov.SingleStarSolveToTOVPoints(ec_fixed, pts);
		// SingleStarSolveToTOVPoints does not emit species labels, but they are a
		// property of the EOS TABLE, not of the resolution. Reuse the list captured
		// from the canonical SolveToProfile call so BuildFromTOV registers species
		// exactly as production does.
		labels = g_species_labels;
	}
	else
	{
		n = tov.SolveToProfile(kTargetM, pts, &labels);
	}
	if (n <= 0 || pts.size() < 4)
		return nullptr;
	// Reject a degenerate profile rather than feeding it to BuildFromTOV.
	if (!(pts.back().r > 0.0) || !std::isfinite(pts.back().m) || !(pts.back().p > 0.0))
		return nullptr;

	// Measured radial spacing of the ACTUAL grid — not 1/res, which is not the grid.
	std::vector<double> dr;
	dr.reserve(pts.size());
	for (std::size_t i = 0; i + 1 < pts.size(); ++i)
		dr.push_back(pts[i + 1].r - pts[i].r);
	std::vector<double> sorted = dr;
	std::sort(sorted.begin(), sorted.end());
	c.res = res;
	c.n_profile = n;
	c.dr_min = sorted.empty() ? 0.0 : sorted.front();
	c.dr_max = sorted.empty() ? 0.0 : sorted.back();
	c.dr_med = sorted.empty() ? 0.0 : sorted[sorted.size() / 2];
	c.dr_eff = pts.back().r / static_cast<double>(pts.size() - 1);

	auto ns = labels.empty() ? std::make_unique<NStar>(pts)
							 : std::make_unique<NStar>(pts, labels);
	const auto &seq = ns->GetSequence();
	c.ec = seq.ec;
	c.pc = seq.pc;
	c.M = seq.m;
	c.R = seq.r;
	c.B = seq.b;
	c.z_surf = ns->Profile().ExpNuSurface();
	const double G = 6.67430e-8, c2 = 2.99792458e10 * 2.99792458e10, Msun = 1.98892e33;
	c.compactness = 2.0 * G * c.M * Msun / (c.R * 1.0e5 * c2);
	return ns;
}

// ---------------------------------------------------------------------------
//  Static thermal diagnostics at a fixed T_inf, through the production drivers.
// ---------------------------------------------------------------------------
static bool StaticDiagnostics(NStar &ns, CompOSE_Thermo &thermo, double Tinf_K, Case &c)
{
	P::Evolution::StarContext sc(ns.Profile());
	P::Evolution::GeometryCache geo(sc);
	P::Evolution::Config cfg;
	cfg.couple_spin = false;
	cfg.n_eta = 0;
	Th::Boundary::EnvelopePotekhin1997_Iron env;
	P::Evolution::DriverContext ctx =
		P::Evolution::Run::MakeDriverContext(sc, geo, cfg, &env, &thermo);

	P::State::ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(Tinf_K);
	P::Evolution::StateVector Y;
	Y.Register(P::State::StateTag::Thermal, thermal);

	Th::PhotonCooling::Options po;
	po.surface_model = Th::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;
	po.radiating_fraction = 1.0;
	po.global_scale = 1.0;
	Th::PhotonCooling photon(po);

	Th::NeutrinoCooling::Options no;
	no.include_direct_urca = true;
	no.include_modified_urca = true;
	no.include_pair_breaking = false;
	no.global_scale = 1.0;
	Th::NeutrinoCooling neutrino(no);

	const auto pd = Th::Detail::ComputeDerived(photon, Y, ctx);
	const auto nd = Th::Detail::NeutrinoCooling_Details::ComputeDerived(neutrino, Y, ctx);
	if (!pd.ok || !nd.ok)
		return false;
	c.C_star = pd.C_star_erg_K;
	c.L_gamma = pd.L_gamma_inf_erg_s;
	c.L_nu = nd.L_nu_inf_erg_s;
	c.L_nu_DU = nd.L_nu_DU_inf_erg_s;
	c.L_nu_MU = nd.L_nu_MU_inf_erg_s;
	c.dlnT_dt = pd.dLnTinf_dt_1_s + nd.dLnTinf_dt_1_s;
	return true;
}

// ---------------------------------------------------------------------------
//  Passive-cooling trajectory — the same problem as passive_cooling_regression.
// ---------------------------------------------------------------------------
static bool CoolingTrajectory(NStar &ns, CompOSE_Thermo &thermo, double rtol, double atol,
							  std::vector<double> &Tinf, std::vector<double> &Lnu,
							  std::vector<double> &Lgam, std::vector<double> &Cst,
							  std::string &err)
{
	P::Evolution::StarContext sc(ns.Profile());
	P::Evolution::GeometryCache geo(sc);

	P::Evolution::Config cfg;
	cfg.couple_spin = false;
	cfg.n_eta = 0;
	cfg.rtol = rtol;
	cfg.atol = atol;
	cfg.max_internal_steps = 1'000'000;
	cfg.max_samples = 1'000'000;
	cfg.stepper = kStepper;
	cfg.save_cadence = P::Evolution::SaveCadence::LogTime;
	cfg.samples_per_decade = kSamplesPerDecade;

	Th::Boundary::EnvelopePotekhin1997_Iron env;
	P::Evolution::DriverContext ctx =
		P::Evolution::Run::MakeDriverContext(sc, geo, cfg, &env, &thermo);

	P::State::ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(kTinf0_K);

	P::Evolution::Run::StateWiring w;
	std::vector<P::State::StateTag> tags{P::State::StateTag::Thermal};
	w.state_vec.Register(P::State::StateTag::Thermal, thermal);
	P::Evolution::Run::ConfigureLayout(w, tags);
	P::Evolution::Run::ConfigureRHS(w, tags);

	Th::PhotonCooling::Options po;
	po.surface_model = Th::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;
	po.radiating_fraction = 1.0;
	po.global_scale = 1.0;
	auto photon = std::make_shared<Th::PhotonCooling>(po);

	Th::NeutrinoCooling::Options no;
	no.include_direct_urca = true;
	no.include_modified_urca = true;
	no.include_pair_breaking = false;
	no.global_scale = 1.0;
	auto neutrino = std::make_shared<Th::NeutrinoCooling>(no);

	std::vector<P::Evolution::DriverPtr> drivers{photon, neutrino};
	P::Evolution::EvolutionSystem system(ctx, w.state_vec, w.rhs, w.layout,
										 std::move(drivers));

	struct Capture final : P::Evolution::Observers::IObserver
	{
		const Th::PhotonCooling *ph;
		const Th::NeutrinoCooling *nu;
		std::vector<double> *T, *Ln, *Lg, *C;
		std::string *err;
		std::size_t next = 0;

		void Record(double t_s, const P::Evolution::StateVector &Y,
					const P::Evolution::DriverContext &c)
		{
			if (next >= kEpochsYr.size())
				return;
			const double target = kEpochsYr[next] * Zaki::Physics::YR_2_SEC;
			if (t_s < target * (1.0 - 1e-9))
				return;
			const auto pd = Th::Detail::ComputeDerived(*ph, Y, c);
			const auto nd = Th::Detail::NeutrinoCooling_Details::ComputeDerived(*nu, Y, c);
			if (!pd.ok || !nd.ok)
			{
				*err = "driver diagnostics not ok at t=" + std::to_string(t_s);
				next = kEpochsYr.size();
				return;
			}
			T->push_back(pd.Tinf_K);
			Ln->push_back(nd.L_nu_inf_erg_s);
			Lg->push_back(pd.L_gamma_inf_erg_s);
			C->push_back(pd.C_star_erg_K);
			++next;
		}
		void OnStart(const P::Evolution::Observers::RunInfo &r,
					 const P::Evolution::StateVector &Y,
					 const P::Evolution::DriverContext &c) override { Record(r.t0, Y, c); }
		void OnSample(const P::Evolution::Observers::SampleInfo &s_,
					  const P::Evolution::StateVector &Y,
					  const P::Evolution::DriverContext &c) override { Record(s_.t, Y, c); }
		std::string Name() const override { return "GridConvergenceCapture"; }
	};

	Tinf.clear(); Lnu.clear(); Lgam.clear(); Cst.clear();
	auto cap = std::make_shared<Capture>();
	cap->ph = photon.get();
	cap->nu = neutrino.get();
	cap->T = &Tinf; cap->Ln = &Lnu; cap->Lg = &Lgam; cap->C = &Cst;
	cap->err = &err;
	system.AddObserver(cap);

	P::Evolution::GSLIntegrator integrator(system, cfg, w.dim);
	std::vector<double> y(w.dim);
	P::Evolution::PackStateVector(w.state_vec, w.layout, y.data());

	const double YR = Zaki::Physics::YR_2_SEC;
	if (!integrator.Integrate(kT0_yr * YR, kT1_yr * YR, y.data()))
	{
		err = "integration failed";
		return false;
	}
	if (Tinf.size() != kEpochsYr.size())
	{
		err = "captured " + std::to_string(Tinf.size()) + " of " +
			  std::to_string(kEpochsYr.size()) + " epochs";
		return false;
	}
	return true;
}

// ---------------------------------------------------------------------------
//  Trajectory norms
// ---------------------------------------------------------------------------
struct Norm
{
	double dmax = 0;   // max_k |ln(T1/T2)|
	double drms = 0;   // RMS_k |ln(T1/T2)|
	double dfinal = 0; // |ln| at the last epoch
	std::size_t k_max = 0;
};
static Norm TrajectoryNorm(const std::vector<double> &a, const std::vector<double> &b)
{
	Norm n;
	double s2 = 0.0;
	for (std::size_t k = 0; k < a.size() && k < b.size(); ++k)
	{
		const double d = std::fabs(std::log(a[k] / b[k]));
		if (d > n.dmax) { n.dmax = d; n.k_max = k; }
		s2 += d * d;
		n.dfinal = d;
	}
	n.drms = std::sqrt(s2 / static_cast<double>(a.size()));
	return n;
}

/// Observed order from three successive values, using the ACTUAL spacing ratio.
///
/// Returns NaN — reported as "not reliably measurable" — when the differences are not
/// interpretable: zero, sign-flipped (no monotone approach), or already at the
/// floating-point roundoff floor of the quantity itself. Producing a number in those
/// cases would be meaningless, so none is produced.
static constexpr double kRoundoffFloor = 1.0e-9; // relative to |q|
static double ObservedOrder(double q1, double q2, double q3, double h1, double h2, double h3)
{
	const double d12 = q1 - q2, d23 = q2 - q3;
	if (!(std::fabs(d23) > 0.0) || !(std::fabs(d12) > 0.0))
		return std::nan("");
	if (d12 * d23 <= 0.0) // sign flip: not a monotone approach, order meaningless
		return std::nan("");
	const double scale = std::fabs(q3) > 0.0 ? std::fabs(q3) : 1.0;
	if (std::fabs(d23) / scale < kRoundoffFloor) // differences are roundoff, not truncation
		return std::nan("");
	const double r1 = h1 / h2, r2 = h2 / h3;
	const double r = 0.5 * (r1 + r2);
	if (!(r > 1.0))
		return std::nan("");
	return std::log(std::fabs(d12 / d23)) / std::log(r);
}
static std::string OrderStr(double p)
{
	if (!std::isfinite(p))
		return "not reliably measurable";
	char b[64];
	snprintf(b, sizeof(b), "%.3f", p);
	return b;
}

// ===========================================================================
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: grid_convergence_cmf <EOS_DATA_ROOT> [--emit-dir <dir>]\n";
		return 2;
	}
	const fs::path root = argv[1];
	fs::path emit_dir;
	for (int i = 2; i < argc; ++i)
		if (std::strcmp(argv[i], "--emit-dir") == 0 && i + 1 < argc)
			emit_dir = argv[++i];

	const fs::path dist = root / "DS-CMF-1-with-crust";
	const fs::path cold = dist / "DS(CMF)-1_with_crust.eos";
	const fs::path fineT = root / "DNS-CMF-Hadronic-with-electrons";
	if (!fs::exists(cold) || !fs::exists(fineT / "eos.thermo"))
	{
		std::cerr << "authenticated EOS data missing under " << root.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / "compactstar_grid_convergence";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	CompOSE_Thermo::Options thOpt;
	thOpt.use_central_difference = true;
	thOpt.clamp_to_domain = true;
	thOpt.Tmin_for_derivative_MeV = 0.0;
	CompOSE_Thermo thermo(fineT.string(), thOpt);
	if (!thermo.IsLoaded())
	{
		std::cerr << "CompOSE_Thermo failed to load\n";
		return 3;
	}

	std::cout << std::scientific << std::setprecision(8);
	std::cout << "Phase 2B-4A — nonrotating TOV -> cooling grid convergence\n"
			  << "DS(CMF)-1_with_crust, target " << kTargetM << " Msun, "
			  << kT0_yr << " -> " << kT1_yr << " yr, RKF45, rtol " << kRtol
			  << ", atol " << kAtol << "\n"
			  << "Scope: NONROTATING slice. No Hartle claim; I is not an observable here.\n\n";

	// =======================================================================
	//  EQ — production equivalence of the variable-resolution builder at 10000
	// =======================================================================
	std::cout << "EQ  production equivalence at radial_res = 10000\n";
	Case eq_canon, eq_test;
	std::unique_ptr<NStar> ns_canon, ns_test;
	{
		// Canonical: exactly what production calls.
		ns_canon = std::make_unique<NStar>();
		ns_canon->SetWrkDir((wrk / "canon").string());
		const int n = ns_canon->SolveTOV_Profile(cold.string(), kTargetM, "eq");
		if (n <= 0)
		{
			std::cerr << "canonical SolveTOV_Profile failed\n";
			return 3;
		}
		const auto &s = ns_canon->GetSequence();
		eq_canon.n_profile = n;
		eq_canon.ec = s.ec; eq_canon.pc = s.pc; eq_canon.M = s.m; eq_canon.R = s.r;
		eq_canon.B = s.b;
		eq_canon.z_surf = ns_canon->Profile().ExpNuSurface();

		// Capture the EOS species-label list once, from the canonical production call.
		{
			TOVSolver probe;
			probe.SetWrkDir((wrk / "labels").string());
			probe.ImportEOS(cold.string(), true);
			std::vector<TOVPoint> tmp;
			probe.SolveToProfile(kTargetM, tmp, &g_species_labels);
		}

		// Test-side builder at the production default resolution.
		ns_test = BuildAtResolution(cold, wrk / "eqtest", 10000, false, 0.0, eq_test);
		if (!ns_test)
		{
			std::cerr << "test-side builder failed at res=10000\n";
			return 3;
		}
		Report("EQ1 radial point count identical",
			   eq_test.n_profile == eq_canon.n_profile,
			   std::to_string(eq_test.n_profile) + " vs " + std::to_string(eq_canon.n_profile));
		Report("EQ2 achieved mass identical", Rel(eq_test.M, eq_canon.M) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.M, eq_canon.M)));
		Report("EQ3 radius identical", Rel(eq_test.R, eq_canon.R) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.R, eq_canon.R)));
		Report("EQ4 central density identical", Rel(eq_test.ec, eq_canon.ec) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.ec, eq_canon.ec)));
		Report("EQ5 central pressure identical", Rel(eq_test.pc, eq_canon.pc) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.pc, eq_canon.pc)));
		Report("EQ6 baryon number identical", Rel(eq_test.B, eq_canon.B) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.B, eq_canon.B)));
		Report("EQ7 surface redshift factor identical",
			   Rel(eq_test.z_surf, eq_canon.z_surf) < 1e-12,
			   "e^nu(R) " + std::to_string(eq_test.z_surf) + " vs " +
				   std::to_string(eq_canon.z_surf));
		Report("EQ8 species labels present and equal in count",
			   ns_test->Profile().Radial().Dim().size() ==
				   ns_canon->Profile().Radial().Dim().size(),
			   "columns " + std::to_string(ns_test->Profile().Radial().Dim().size()));

		// Thermal diagnostics and the full trajectory must agree too.
		StaticDiagnostics(*ns_canon, thermo, kTstatic_K, eq_canon);
		StaticDiagnostics(*ns_test, thermo, kTstatic_K, eq_test);
		Report("EQ9 C_star identical", Rel(eq_test.C_star, eq_canon.C_star) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.C_star, eq_canon.C_star)));
		Report("EQ10 L_nu identical", Rel(eq_test.L_nu, eq_canon.L_nu) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.L_nu, eq_canon.L_nu)));
		Report("EQ11 L_gamma identical", Rel(eq_test.L_gamma, eq_canon.L_gamma) < 1e-12,
			   "rel " + std::to_string(Rel(eq_test.L_gamma, eq_canon.L_gamma)));

		std::string err;
		std::vector<double> Ta, La, Ga, Ca, Tb, Lb, Gb, Cb;
		const bool oka = CoolingTrajectory(*ns_canon, thermo, kRtol, kAtol, Ta, La, Ga, Ca, err);
		const bool okb = CoolingTrajectory(*ns_test, thermo, kRtol, kAtol, Tb, Lb, Gb, Cb, err);
		Norm nq;
		if (oka && okb)
			nq = TrajectoryNorm(Ta, Tb);
		Report("EQ12 full 9-epoch cooling trajectory identical",
			   oka && okb && nq.dmax < 1e-12,
			   oka && okb ? "max |dln T_inf| = " + std::to_string(nq.dmax) : err);
		eq_canon.Tinf = Ta; eq_test.Tinf = Tb;
	}
	if (g_fail)
	{
		std::cout << "\nSTOP: the variable-resolution builder does not reproduce the "
					 "canonical construction.\nA convergence study on a non-equivalent "
					 "construction would be meaningless.\n";
		return 4;
	}
	const double ec_star = eq_canon.ec; // Experiment A authority
	std::cout << "  => equivalent. Experiment A will hold ec* = " << ec_star
			  << " g/cm^3 fixed (the canonical 1.6 Msun central density).\n";

	// =======================================================================
	//  EXPERIMENT A — fixed central energy density
	// =======================================================================
	std::cout << "\nEXPERIMENT A  fixed ec* — isolates radial discretization\n";
	std::vector<Case> A;
	for (std::size_t res : kResolutions)
	{
		Case c;
		const auto t0 = std::chrono::steady_clock::now();
		auto ns = BuildAtResolution(cold, wrk / ("A" + std::to_string(res)), res, true,
									ec_star, c);
		if (!ns) { std::cout << "    res=" << res << " FAILED\n"; A.push_back(c); continue; }
		c.ok = StaticDiagnostics(*ns, thermo, kTstatic_K, c);
		c.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		A.push_back(c);
		std::cout << "    res=" << std::setw(6) << res << "  N=" << std::setw(6) << c.n_profile
				  << "  dr_eff=" << c.dr_eff << " km  M=" << c.M << "  R=" << c.R
				  << "  C*=" << c.C_star << "  Lnu=" << c.L_nu << "  (" << std::fixed
				  << std::setprecision(1) << c.seconds << "s)\n" << std::scientific
				  << std::setprecision(8);
	}

	// =======================================================================
	//  EXPERIMENT B — fixed target mass, the production workflow
	// =======================================================================
	std::cout << "\nEXPERIMENT B  fixed target mass " << kTargetM
			  << " Msun — the production workflow\n";
	std::vector<Case> B;
	for (std::size_t res : kResolutions)
	{
		Case c;
		const auto t0 = std::chrono::steady_clock::now();
		auto ns = BuildAtResolution(cold, wrk / ("B" + std::to_string(res)), res, false, 0.0, c);
		if (!ns) { std::cout << "    res=" << res << " FAILED\n"; B.push_back(c); continue; }
		StaticDiagnostics(*ns, thermo, kTstatic_K, c);
		std::string err;
		c.ok = CoolingTrajectory(*ns, thermo, kRtol, kAtol, c.Tinf, c.Lnu_t, c.Lgam_t,
								 c.Cstar_t, err);
		c.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		if (!c.ok)
			std::cout << "    res=" << res << " cooling FAILED: " << err << "\n";
		B.push_back(c);
		std::cout << "    res=" << std::setw(6) << res << "  N=" << std::setw(6) << c.n_profile
				  << "  dr_eff=" << c.dr_eff << " km  M=" << c.M << "  ec=" << c.ec
				  << "  R=" << c.R
				  << (c.ok ? "  Tinf(1Myr)=" + std::to_string(c.Tinf.back()) : "")
				  << "  (" << std::fixed << std::setprecision(1) << c.seconds << "s)\n"
				  << std::scientific << std::setprecision(8);
	}

	// =======================================================================
	//  Convergence analysis
	// =======================================================================
	// Order fitting uses the four finest resolutions; 2500 is retained and reported
	// as a coarse diagnostic (and as the detector proof) but is excluded from order
	// fitting if it is not in the asymptotic regime. That exclusion is stated, never silent.
	auto fit = [&](const std::vector<Case> &S, auto get, const char *name)
	{
		std::cout << "    " << std::setw(22) << std::left << name << std::right;
		// successive differences over the whole sequence
		for (std::size_t i = 0; i + 1 < S.size(); ++i)
			std::cout << "  |d" << S[i].res << "-" << S[i + 1].res
					  << "|=" << std::setprecision(3) << std::fabs(get(S[i]) - get(S[i + 1]));
		std::cout << std::setprecision(8) << "\n";
		// order from the three finest (indices n-3, n-2, n-1)
		const std::size_t n = S.size();
		if (n < 3)
			return std::nan("");
		const double p = ObservedOrder(get(S[n - 3]), get(S[n - 2]), get(S[n - 1]),
									   S[n - 3].dr_eff, S[n - 2].dr_eff, S[n - 1].dr_eff);
		const double q3 = get(S[n - 1]);
		const double rel23 = std::fabs(q3) > 0.0
								 ? std::fabs(get(S[n - 2]) - q3) / std::fabs(q3)
								 : 0.0;
		std::cout << "        observed order (3 finest): " << OrderStr(p)
				  << "   [finest relative difference " << Sci(rel23) << "]";
		if (!std::isfinite(p) && rel23 < kRoundoffFloor)
			std::cout << "  <- at the roundoff floor of this quantity";
		std::cout << "\n";
		return p;
	};

	std::cout << "\nA-CONV  structural and thermal convergence at fixed ec*\n";
	const double pA_M = fit(A, [](const Case &c) { return c.M; }, "M [Msun]");
	const double pA_R = fit(A, [](const Case &c) { return c.R; }, "R [km]");
	fit(A, [](const Case &c) { return c.B; }, "B [baryon number]");
	fit(A, [](const Case &c) { return c.z_surf; }, "e^nu(R)");
	const double pA_C = fit(A, [](const Case &c) { return c.C_star; }, "C_star(1e8 K)");
	const double pA_Lnu = fit(A, [](const Case &c) { return c.L_nu; }, "L_nu(1e8 K)");
	fit(A, [](const Case &c) { return c.L_gamma; }, "L_gamma(1e8 K)");
	fit(A, [](const Case &c) { return c.dlnT_dt; }, "dlnT_inf/dt(1e8 K)");

	std::cout << "\nB-CONV  structural convergence at fixed target mass\n";
	fit(B, [](const Case &c) { return c.M; }, "achieved M [Msun]");
	fit(B, [](const Case &c) { return c.ec; }, "selected ec [g/cm^3]");
	fit(B, [](const Case &c) { return c.R; }, "R [km]");
	const double pB_C = fit(B, [](const Case &c) { return c.C_star; }, "C_star(1e8 K)");

	std::cout << "\nB-TRAJ  passive-cooling trajectory norms  D = max_k |ln(T1/T2)|\n";
	std::vector<Norm> Bn;
	for (std::size_t i = 0; i + 1 < B.size(); ++i)
	{
		if (!B[i].ok || !B[i + 1].ok) { Bn.push_back(Norm{}); continue; }
		const Norm nn = TrajectoryNorm(B[i].Tinf, B[i + 1].Tinf);
		Bn.push_back(nn);
		std::cout << "    " << std::setw(6) << B[i].res << " -> " << std::setw(6)
				  << B[i + 1].res << "   Dmax=" << nn.dmax << "  RMS=" << nn.drms
				  << "  final=" << nn.dfinal << "  worst epoch=" << std::fixed
				  << std::setprecision(0) << kEpochsYr[nn.k_max] << " yr\n"
				  << std::scientific << std::setprecision(8);
	}
	double pB_traj = std::nan("");
	if (Bn.size() >= 3 && Bn[Bn.size() - 1].dmax > 0.0 && Bn[Bn.size() - 2].dmax > 0.0)
	{
		const std::size_t n = B.size();
		const double r = 0.5 * (B[n - 3].dr_eff / B[n - 2].dr_eff +
								B[n - 2].dr_eff / B[n - 1].dr_eff);
		pB_traj = std::log(Bn[Bn.size() - 2].dmax / Bn[Bn.size() - 1].dmax) / std::log(r);
		std::cout << "    trajectory-norm contraction order (2 finest pairs): "
				  << OrderStr(pB_traj) << "\n";
	}

	// Final-epoch T_inf convergence
	const double pB_Tfin = fit(B, [](const Case &c) { return c.ok ? c.Tinf.back() : 0.0; },
							   "T_inf(1 Myr) [K]");

	// =======================================================================
	//  Richardson estimates — only where the order is stable and positive
	// =======================================================================
	std::cout << "\nRICHARDSON  continuum estimates (only where the measured order supports it)\n";
	auto richardson = [&](const std::vector<Case> &S, auto get, double p, const char *name)
	{
		const std::size_t n = S.size();
		if (!std::isfinite(p) || p <= 0.1 || n < 2)
		{
			std::cout << "    " << std::setw(22) << std::left << name << std::right
					  << "  not justified (order " << OrderStr(p) << ")\n";
			return std::nan("");
		}
		const double r = S[n - 2].dr_eff / S[n - 1].dr_eff;
		const double qf = get(S[n - 1]), qc = get(S[n - 2]);
		const double est = qf + (qf - qc) / (std::pow(r, p) - 1.0);
		std::cout << "    " << std::setw(22) << std::left << name << std::right
				  << "  finest = " << qf << "   extrapolated = " << est
				  << "   |diff|/|finest| = " << Rel(qf, est) << "\n";
		return est;
	};
	richardson(A, [](const Case &c) { return c.R; }, pA_R, "A: R [km]");
	richardson(A, [](const Case &c) { return c.C_star; }, pA_C, "A: C_star(1e8 K)");
	richardson(A, [](const Case &c) { return c.M; }, pA_M, "A: M [Msun]");
	richardson(B, [](const Case &c) { return c.ok ? c.Tinf.back() : 0.0; }, pB_Tfin,
			   "B: T_inf(1 Myr) [K]");

	// =======================================================================
	//  Thermal-integrator floor
	// =======================================================================
	std::cout << "\nFLOOR  thermal-integrator error vs radial-refinement error\n";
	std::cout << "    (test configuration only; production defaults untouched)\n";
	double floor_10k = -1.0, floor_fine = -1.0;
	{
		const std::vector<std::pair<const char *, std::size_t>> probes = {
			{"res=10000", 10000}, {"res=40000", 40000}};
		for (const auto &pr : probes)
		{
			Case c;
			auto ns = BuildAtResolution(cold, wrk / ("F" + std::to_string(pr.second)),
										pr.second, false, 0.0, c);
			if (!ns) continue;
			std::vector<double> T1, L1, G1, C1, T2, L2, G2, C2;
			std::string err;
			const bool o1 = CoolingTrajectory(*ns, thermo, kRtol, kAtol, T1, L1, G1, C1, err);
			// tighten BOTH tolerances by 100x; the physical problem is untouched
			const bool o2 = CoolingTrajectory(*ns, thermo, kRtol * 1e-2, kAtol * 1e-2,
											  T2, L2, G2, C2, err);
			if (!o1 || !o2) { std::cout << "    " << pr.first << " floor probe failed\n"; continue; }
			const Norm nf = TrajectoryNorm(T1, T2);
			if (pr.second == 10000) floor_10k = nf.dmax; else floor_fine = nf.dmax;
			std::cout << "    " << pr.first << "  rtol " << kRtol << "->" << kRtol * 1e-2
					  << ", atol " << kAtol << "->" << kAtol * 1e-2
					  << "   Dmax=" << nf.dmax << "  (worst epoch " << std::fixed
					  << std::setprecision(0) << kEpochsYr[nf.k_max] << " yr)\n"
					  << std::scientific << std::setprecision(8);
		}
	}

	// =======================================================================
	//  Default-grid discretization error and the coarse detector
	// =======================================================================
	std::cout << "\nDEFAULT-GRID ERROR  res=10000 against the finest computed grid\n";
	double def_err_T = -1.0, def_err_R = -1.0, def_err_C = -1.0;
	std::size_t i10 = 0, ifin = B.size() - 1;
	for (std::size_t i = 0; i < B.size(); ++i)
		if (B[i].res == 10000) i10 = i;
	if (B[i10].ok && B[ifin].ok)
	{
		def_err_T = TrajectoryNorm(B[i10].Tinf, B[ifin].Tinf).dmax;
		def_err_R = Rel(B[i10].R, B[ifin].R);
		def_err_C = Rel(B[i10].C_star, B[ifin].C_star);
		std::cout << "    R                  rel = " << def_err_R << "\n"
				  << "    C_star(1e8 K)      rel = " << def_err_C << "\n"
				  << "    T_inf trajectory   Dmax = " << def_err_T << "\n";
	}

	std::cout << "\nDETECTOR  can the analysis distinguish a materially coarse grid?\n";
	double coarse_T = -1.0, coarse_R = -1.0;
	if (B[0].ok && B[ifin].ok)
	{
		coarse_T = TrajectoryNorm(B[0].Tinf, B[ifin].Tinf).dmax;
		coarse_R = Rel(B[0].R, B[ifin].R);
		std::cout << "    res=" << B[0].res << " vs finest:  R rel = " << coarse_R
				  << ",  T_inf Dmax = " << coarse_T << "\n";
		std::cout << "    res=10000  vs finest:  R rel = " << def_err_R
				  << ",  T_inf Dmax = " << def_err_T << "\n";
	}

	// =======================================================================
	//  Assertions — convergence EVIDENCE only. No invented accuracy threshold.
	// =======================================================================
	std::cout << "\nCRITERIA\n";
	Report("G1 the coarse grid is detected as materially different from the finest",
		   coarse_R > 0.0 && def_err_R > 0.0 && coarse_R > 10.0 * def_err_R,
		   "coarse R error " + Sci(coarse_R) + " vs default " + Sci(def_err_R) +
			   "  (ratio " + Sci(def_err_R > 0 ? coarse_R / def_err_R : 0.0, 1) + ")");
	Report("G2 refinement contracts the trajectory norm",
		   Bn.size() >= 3 && Bn[Bn.size() - 1].dmax < Bn[Bn.size() - 2].dmax &&
			   Bn[Bn.size() - 2].dmax < Bn[Bn.size() - 3].dmax,
		   Bn.size() >= 3 ? "D: " + Sci(Bn[Bn.size() - 3].dmax) + " -> " +
								Sci(Bn[Bn.size() - 2].dmax) + " -> " +
								Sci(Bn[Bn.size() - 1].dmax)
						  : "insufficient pairs");
	Report("G3 Experiment A radius contracts under refinement",
		   A.size() >= 4 && std::fabs(A[A.size() - 1].R - A[A.size() - 2].R) <
								std::fabs(A[A.size() - 2].R - A[A.size() - 3].R),
		   "successive |dR| decreasing on the finest pairs");
	Report("G4 thermal-integrator error is subdominant to radial-refinement error",
		   floor_10k >= 0.0 && def_err_T > 0.0 && floor_10k < 0.1 * def_err_T,
		   "integrator floor " + Sci(floor_10k) + " vs default-grid radial error " +
			   Sci(def_err_T) + "  (floor is " +
			   Sci(floor_10k > 0 ? def_err_T / floor_10k : 0.0, 1) + "x smaller)");
	Report("G5 the target-mass search stays inside its own tolerance at every resolution",
		   std::all_of(B.begin(), B.end(),
					   [](const Case &c) { return c.n_profile > 0 &&
											std::fabs(c.M - kTargetM) < 1e-4; }),
		   "max |M - 1.6| = " +
			   Sci(std::max_element(B.begin(), B.end(),
									[](const Case &a, const Case &b) {
										return std::fabs(a.M - kTargetM) <
											   std::fabs(b.M - kTargetM);
									})
					   ->M - kTargetM) +
			   " (tolerance 1e-4 Msun)");

	// =======================================================================
	//  Artifacts
	// =======================================================================
	if (!emit_dir.empty())
	{
		fs::create_directories(emit_dir);
		std::ofstream o(emit_dir / "grid_convergence_cmf_1p6_debug.tsv");
		o << "# Phase 2B-4A convergence matrix — DS(CMF)-1_with_crust, nonrotating slice\n"
		  << "# static thermal diagnostics evaluated at T_inf = " << kTstatic_K << " K\n"
		  << "# NOT a golden cooling curve; the golden baseline remains "
			 "passive_cooling_cmf_1p6_debug.tsv\n"
		  << "experiment\tradial_res\tN_profile\tdr_eff_km\tdr_min_km\tdr_max_km"
			 "\tec_gcm3\tachieved_M\tR_km\tB\tz_surf\tcompactness"
			 "\tCstar_1e8\tLnu_1e8\tLgamma_1e8\tdlnT_dt_1e8\n";
		o << std::setprecision(12) << std::scientific;
		auto row = [&](const char *tag, const Case &c) {
			o << tag << "\t" << c.res << "\t" << c.n_profile << "\t" << c.dr_eff << "\t"
			  << c.dr_min << "\t" << c.dr_max << "\t" << c.ec << "\t" << c.M << "\t"
			  << c.R << "\t" << c.B << "\t" << c.z_surf << "\t" << c.compactness << "\t"
			  << c.C_star << "\t" << c.L_nu << "\t" << c.L_gamma << "\t" << c.dlnT_dt << "\n";
		};
		for (const auto &c : A) row("A_fixed_ec", c);
		for (const auto &c : B) row("B_fixed_mass", c);

		std::ofstream t(emit_dir / "grid_convergence_cmf_1p6_trajectory.tsv");
		t << "# Phase 2B-4A passive-cooling trajectories per radial resolution "
			 "(Experiment B)\n"
		  << "# NOT a golden baseline. Reference values live in "
			 "passive_cooling_cmf_1p6_debug.tsv\n"
		  << "radial_res\tepoch_yr\tTinf_K\tL_nu_inf\tL_gamma_inf\tC_star\n";
		t << std::setprecision(12) << std::scientific;
		for (const auto &c : B)
			if (c.ok)
				for (std::size_t k = 0; k < kEpochsYr.size(); ++k)
					t << c.res << "\t" << kEpochsYr[k] << "\t" << c.Tinf[k] << "\t"
					  << c.Lnu_t[k] << "\t" << c.Lgam_t[k] << "\t" << c.Cstar_t[k] << "\n";
		std::cout << "\n  artifacts written to " << emit_dir.string() << "\n";
	}

	fs::remove_all(wrk);
	std::cout << "\n"
			  << (g_fail == 0 ? "convergence criteria satisfied"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
