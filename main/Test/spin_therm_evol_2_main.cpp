//==============================================================
// -*- lsst-c++ -*-
/*
 * CompactStar test main: coupled Spin + Thermal evolution (refactored)
 *
 * NEW PATH:
 *   NStar::SolveTOV_Profile(...) -> StarProfile -> StarContext(StarProfile) -> GeometryCache
 *
 * Then:
 *   ThermalState + PhotonCooling + NeutrinoCooling
 *   SpinState    + MagneticDipole
 *   EvolutionSystem + GSLIntegrator
 *
 * Refactor goal:
 *  - Keep star-building in main (Core responsibility)
 *  - Use Evolution/Run helpers to eliminate boilerplate:
 *      RunPaths     : canonical output paths
 *      RunBuilder   : cfg defaults, ctx wiring, layout/rhs wiring, diag-driver collection
 *      RunObservers : standard observer presets/factories
 */

#include <cstdio>
#include <gsl/gsl_errno.h>

#include <iostream>
#include <memory>
#include <vector>

#include <Zaki/String/Directory.hpp>
#include <Zaki/Util/Logger.hpp>

// -------------------------------------------------------
// Evolution infrastructure
// -------------------------------------------------------
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StatePacking.hpp"

// Run helpers (NEW)
#include "CompactStar/Physics/Evolution/Run/RunBuilder.hpp"
#include "CompactStar/Physics/Evolution/Run/RunObservers.hpp"
#include "CompactStar/Physics/Evolution/Run/RunPaths.hpp"

// -------------------------------------------------------
// State blocks
// -------------------------------------------------------
#include "CompactStar/Physics/State/SpinState.hpp"
#include "CompactStar/Physics/State/Tags.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

// -------------------------------------------------------
// Drivers
// -------------------------------------------------------
#include "CompactStar/Physics/Driver/Spin/MagneticDipole.hpp"
#include "CompactStar/Physics/Driver/Thermal/Boundary/EnvelopePotekhin1997.hpp"
#include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"

// -------------------------------------------------------
// Core: build the star profile the NEW way
// -------------------------------------------------------
#include "CompactStar/Core/NStar.hpp"

// -------------------------------------------------------
// GSL error handler (so GSL does not abort the program)
// -------------------------------------------------------
static void
my_gsl_error_handler(const char *reason,
					 const char *file,
					 int line,
					 int gsl_errno)
{
	std::fprintf(stderr,
				 "GSL ERROR: %s:%d: %s (%s)\n",
				 file, line, reason, gsl_strerror(gsl_errno));
	// keep going — this is a debug main
}

//==============================================================
int main()
{
	// Resolve path roots relative to this translation unit.
	const Zaki::String::Directory this_file_dir(__FILE__);

	// Canonical run paths under "<...>/results/<out_dir>/"
	const Zaki::String::Directory base_results_dir = this_file_dir.ParentDir() + "/results";
	const Zaki::String::Directory out_dir = "spin_therm_evol_2";

	// Zaki::Util::Instrumentor::BeginSession("spin-therm-evol-2-main", base_results_dir + out_dir + "/timing_profile.json");

	PROFILE_FUNCTION();

	using namespace CompactStar;

	// -------------------------------------------------------
	// 0) GSL + logging
	// -------------------------------------------------------
	gsl_set_error_handler(&my_gsl_error_handler);

	Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Info);
	Zaki::Util::LogManager::SetBlackWhite(false);

	const Physics::Evolution::Run::RunPaths paths =
		Physics::Evolution::Run::MakeRunPaths(base_results_dir, out_dir, "spin_therm_evol_2_main.log");

	Zaki::Util::LogManager::SetLogFile(paths.log_file);

	// -------------------------------------------------------
	// 1) Build a star the NEW way: SolveTOV_Profile -> StarProfile
	// -------------------------------------------------------
	Zaki::String::Directory eos_root =
		this_file_dir.ParentDir().ParentDir() + "/EOS/CompOSE/";
	const std::string eos_name = "DS(CMF)-1_with_crust";
	const Zaki::String::Directory eos_file =
		eos_root + eos_name + "/" + eos_name + ".eos";

	Core::NStar ns;
	ns.SetWrkDir(base_results_dir);

	const double target_M = 1.6; // Mass in Msun
	const int n_rows = ns.SolveTOV_Profile(eos_file, target_M, out_dir);

	if (n_rows <= 0)
	{
		Z_LOG_ERROR("SolveTOV_Profile failed (n_rows <= 0).");
		return 1;
	}

	// Optional: export profile for sanity checks
	ns.Export(out_dir + "/NStar_Profile.tsv");

	// -------------------------------------------------------
	// 2) StarContext + GeometryCache from *StarProfile*
	// -------------------------------------------------------
	Physics::Evolution::StarContext starCtx(ns.Profile());
	Physics::Evolution::GeometryCache geo(starCtx);

	if (geo.R().Size() == 0 || geo.Exp2Nu().Size() == 0)
	{
		Z_LOG_ERROR("GeometryCache is empty (R or Exp2Nu size is 0).");
		return 1;
	}

	// -------------------------------------------------------
	// 3) Load CompOSE thermo tables for Cv and neutrino cooling
	// -------------------------------------------------------
	const std::string general_eos_name = "CMF-1_general";
	const Zaki::String::Directory eos_dir = eos_root + general_eos_name;

	CompactStar::EOS::CompOSE_Thermo::Options thOpt;
	thOpt.use_central_difference = true;
	thOpt.clamp_to_domain = true;
	thOpt.Tmin_for_derivative_MeV = 0.0;

	CompactStar::EOS::CompOSE_Thermo thermo(eos_dir.Str(), thOpt);
	if (!thermo.IsLoaded())
	{
		Z_LOG_ERROR("CompOSE_Thermo did not load (IsLoaded() == false). Check eos_dir and required tables.");
		return 1;
	}

	// -------------------------------------------------------
	// 4) Evolution configuration (defaults + overrides)
	// -------------------------------------------------------
	Physics::Evolution::Config cfg = Physics::Evolution::Run::MakeDefaultConfig();
	cfg.couple_spin = true;
	cfg.n_eta = 0;

	// (optional overrides)
	// Explicit: the canonical scientific run must not depend on an implicit default.
	cfg.stepper = Physics::Evolution::StepperType::RKF45;
	cfg.rtol = 1e-6;
	cfg.atol = 1e-10;
	cfg.max_internal_steps = 1'000'000;
	cfg.max_samples = 1'000'000;
	// cfg.dt_save = 1.0e2 * Zaki::Physics::YR_2_SEC; // save every 100 years

	cfg.save_cadence = Physics::Evolution::SaveCadence::LogTime;
	// cfg.log_t_floor = 1.0 * Zaki::Physics::YR_2_SEC; // start log spacing at 1 year (or 1 day)
	cfg.log_t_floor = 1.0;			//
	cfg.samples_per_decade = 150.0; //  per decade
	// -------------------------------------------------------
	// 5) DriverContext wiring
	// -------------------------------------------------------
	// Physics::Evolution::DriverContext ctx =
	// 	Physics::Evolution::Run::MakeDriverContext(starCtx, geo, cfg);
	Physics::Driver::Thermal::Boundary::EnvelopePotekhin1997_Iron env97_fe;

	Physics::Evolution::DriverContext ctx =
		Physics::Evolution::Run::MakeDriverContext(starCtx, geo, cfg,
												   &env97_fe, &thermo);

	// -------------------------------------------------------
	// 6) Dynamic states: Thermal + Spin
	// -------------------------------------------------------
	Physics::State::ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(1e9); // K

	Physics::State::SpinState spin;
	spin.Resize(1);
	spin.Omega() = 100.0; // rad/s

	// -------------------------------------------------------
	// 7) State wiring (StateVector + Layout + RHS) via RunBuilder helpers
	// -------------------------------------------------------
	Physics::Evolution::Run::StateWiring wiring;

	// Register blocks into StateVector
	wiring.state_vec.Register(Physics::State::StateTag::Thermal, thermal);
	wiring.state_vec.Register(Physics::State::StateTag::Spin, spin);

	// Packing order
	Physics::Evolution::Run::ConfigureLayout(
		wiring,
		std::vector<Physics::State::StateTag>{
			Physics::State::StateTag::Thermal,
			Physics::State::StateTag::Spin});

	// RHS buffers
	Physics::Evolution::Run::ConfigureRHS(
		wiring,
		std::vector<Physics::State::StateTag>{
			Physics::State::StateTag::Thermal,
			Physics::State::StateTag::Spin});

	// -------------------------------------------------------
	// 8) Drivers
	// -------------------------------------------------------
	Physics::Driver::Spin::MagneticDipole::Options spinOpts;
	spinOpts.braking_index = 3.0;
	spinOpts.K_prefactor = 1e-15;
	spinOpts.use_moment_of_inertia = false;

	auto spinDriver =
		std::make_shared<Physics::Driver::Spin::MagneticDipole>(spinOpts);

	Physics::Driver::Thermal::PhotonCooling::Options photonOpts;
	// photonOpts.surface_model =
	// 	Physics::Driver::Thermal::PhotonCooling::Options::SurfaceModel::ApproxFromTinf;
	photonOpts.surface_model =
		Physics::Driver::Thermal::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;
	photonOpts.radiating_fraction = 1.0;
	// Heat capacity is no longer a PhotonCooling option: per ADR-0002 the denominator is
	// the canonical C_*(T_inf) from StarContext + CompOSE_Thermo, shared with NeutrinoCooling.
	photonOpts.global_scale = 1.0;

	auto photonDriver =
		std::make_shared<Physics::Driver::Thermal::PhotonCooling>(photonOpts);

	Physics::Driver::Thermal::NeutrinoCooling::Options nuOpts;
	nuOpts.include_direct_urca = true;
	nuOpts.include_modified_urca = true;
	nuOpts.include_pair_breaking = false;
	nuOpts.global_scale = 1.0;

	auto neutrinoDriver =
		std::make_shared<Physics::Driver::Thermal::NeutrinoCooling>(nuOpts);

	std::vector<Physics::Evolution::DriverPtr> drivers;
	drivers.push_back(spinDriver);
	drivers.push_back(photonDriver);
	drivers.push_back(neutrinoDriver);

	// -------------------------------------------------------
	// 9) Build EvolutionSystem
	// -------------------------------------------------------
	Physics::Evolution::EvolutionSystem system(
		ctx,
		wiring.state_vec,
		wiring.rhs,
		wiring.layout,
		std::move(drivers));

	// -------------------------------------------------------
	// 10) Observers (default presets + factories)
	// -------------------------------------------------------
	// IMPORTANT: CollectDiagnosticsDrivers must be called on the *drivers owned by system*.
	// Since we moved drivers into system above, we cannot reuse the local vector.
	// Policy: query system for its drivers if you have an accessor; if not, keep a copy.
	//
	// Minimal, robust approach: collect diagnostics pointers BEFORE std::move(drivers).
	//
	// We already moved drivers. Therefore we must rebuild the diag list directly from the
	// shared_ptr objects we still hold (spinDriver/photonDriver/neutrinoDriver).
	// That keeps this main compile-safe regardless of EvolutionSystem API.

	std::vector<Physics::Evolution::DriverPtr> driver_refs;
	driver_refs.push_back(spinDriver);
	driver_refs.push_back(photonDriver);
	driver_refs.push_back(neutrinoDriver);

	const auto diag_drivers =
		Physics::Evolution::Run::CollectDiagnosticsDrivers(driver_refs);

	// DiagnosticsObserver
	{
		auto diag = Physics::Evolution::Run::MakeDiagnosticsObserver(paths, diag_drivers);
		system.AddObserver(diag);
	}

	// TimeSeriesObserver (catalog-driven)
	{
		auto ts = Physics::Evolution::Run::MakeTimeSeriesObserver(paths, diag_drivers);
		system.AddObserver(ts);
	}

	// -------------------------------------------------------
	// 11) Pack initial state -> y[]
	// -------------------------------------------------------
	std::vector<double> y(wiring.dim);
	Physics::Evolution::PackStateVector(wiring.state_vec, wiring.layout, y.data());

	std::cout << "Initial conditions:\n";
	std::cout << "  Tinf  = " << thermal.Tinf() << " K\n";
	std::cout << "  Omega = " << spin.Omega() << " rad/s\n";

	// -------------------------------------------------------
	// 12) Integrate
	// -------------------------------------------------------
	Physics::Evolution::GSLIntegrator integrator(system, cfg, wiring.dim);

	const double t0 = 1.0e2 * Zaki::Physics::YR_2_SEC; // 100 yr
	const double t1 = 1.0e6 * Zaki::Physics::YR_2_SEC; // 1 Myr

	const bool ok = integrator.Integrate(t0, t1, y.data());
	if (!ok)
	{
		std::cerr << "GSLIntegrator: integration failed or max_steps exceeded.\n";
		return 1;
	}

	// -------------------------------------------------------
	// 13) Unpack final state and print
	// -------------------------------------------------------
	Physics::Evolution::UnpackStateVector(wiring.state_vec, wiring.layout, y.data());

	std::cout << "Final conditions (t = " << t1 / Zaki::Physics::YR_2_SEC << " yr):\n";
	std::cout << "  Tinf  = " << thermal.Tinf() << " K\n";
	std::cout << "  Omega = " << spin.Omega() << " rad/s\n";

	std::cout << "[debug] done.\n";

	// Zaki::Util::Instrumentor::EndSession();
	return 0;
}
//==============================================================