// // -*- lsst-c++ -*-
// /*
//  * CompactStar test main: coupled Spin + Thermal evolution
//  *
//  * NEW PATH:
//  *   NStar::SolveTOV_Profile(...) -> StarProfile -> StarContext(StarProfile) -> GeometryCache
//  *
//  * Then:
//  *   ThermalState + PhotonCooling
//  *   SpinState    + MagneticDipole
//  *   EvolutionSystem + GSLIntegrator
//  */

// #include <cstdio>
// #include <gsl/gsl_errno.h>
// // #include <gsl/gsl_strerror.h>
// #include <iostream>
// #include <memory>
// #include <vector>

// #include <Zaki/String/Directory.hpp>
// #include <Zaki/Util/Logger.hpp>

// // -------------------------------------------------------
// // Evolution infrastructure
// // -------------------------------------------------------
// #include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
// #include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
// #include "CompactStar/Physics/Evolution/GeometryCache.hpp"
// #include "CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp"
// #include "CompactStar/Physics/Evolution/Observers/DiagnosticsObserver.hpp"
// #include "CompactStar/Physics/Evolution/Observers/TimeSeriesObserver.hpp"
// #include "CompactStar/Physics/Evolution/RHSAccumulator.hpp"
// #include "CompactStar/Physics/Evolution/StarContext.hpp"
// #include "CompactStar/Physics/Evolution/StateLayout.hpp"
// #include "CompactStar/Physics/Evolution/StatePacking.hpp"
// #include "CompactStar/Physics/Evolution/StateVector.hpp"
// // -------------------------------------------------------
// // State blocks
// // -------------------------------------------------------
// #include "CompactStar/Physics/State/SpinState.hpp"
// #include "CompactStar/Physics/State/Tags.hpp"
// #include "CompactStar/Physics/State/ThermalState.hpp"

// // -------------------------------------------------------
// // Drivers
// // -------------------------------------------------------
// #include "CompactStar/Physics/Driver/Spin/MagneticDipole.hpp"
// #include "CompactStar/Physics/Driver/Thermal/NeutrinoCooling.hpp"
// #include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"
// // -------------------------------------------------------
// // Core: build the star profile the NEW way
// // -------------------------------------------------------
// #include "CompactStar/Core/NStar.hpp"

// // -------------------------------------------------------
// // GSL error handler (so GSL does not abort the program)
// // -------------------------------------------------------
// static void
// my_gsl_error_handler(const char *reason,
// 					 const char *file,
// 					 int line,
// 					 int gsl_errno)
// {
// 	std::fprintf(stderr,
// 				 "GSL ERROR: %s:%d: %s (%s)\n",
// 				 file, line, reason, gsl_strerror(gsl_errno));
// 	// keep going — this is a debug main
// }

// //==============================================================
// int main()
// {
// 	using namespace CompactStar;

// 	Zaki::String::Directory dir(__FILE__);
// 	std::cout << "[debug] this file dir        = " << dir << "\n";
// 	std::cout << "[debug] this file parent dir = " << dir.ParentDir() << "\n";

// 	// -------------------------------------------------------
// 	// 0) GSL + logging
// 	// -------------------------------------------------------
// 	gsl_set_error_handler(&my_gsl_error_handler);
// 	Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Info);
// 	Zaki::Util::LogManager::SetBlackWhite(false);
// 	Zaki::Util::LogManager::SetLogFile(dir + "spin_therm_evol_2_main.log");
// 	// -------------------------------------------------------
// 	// 1) Build a star the NEW way: SolveTOV_Profile -> StarProfile
// 	// -------------------------------------------------------
// 	Zaki::String::Directory eos_root =
// 		dir.ParentDir().ParentDir() + "/EOS/CompOSE/";
// 	std::string eos_name = "DS(CMF)-1_with_crust";

// 	Zaki::String::Directory eos_file =
// 		eos_root + eos_name + "/" + eos_name + ".eos";

// 	Zaki::String::Directory base_results_dir = dir.ParentDir() + "/results";
// 	Zaki::String::Directory out_dir = "spin_therm_evol_2";

// 	Core::NStar ns;
// 	ns.SetWrkDir(base_results_dir);

// 	const double target_M = 1.8; // Msun
// 	const int n_rows = ns.SolveTOV_Profile(eos_file, target_M, out_dir);

// 	if (n_rows <= 0)
// 	{
// 		Z_LOG_ERROR("SolveTOV_Profile failed (n_rows <= 0).");
// 		return 1;
// 	}

// 	// Optional: export the profile (handy for sanity checks)
// 	ns.Export(out_dir + "/NStar_Profile.tsv");

// 	// -------------------------------------------------------
// 	// 2) StarContext + GeometryCache from *StarProfile*
// 	// -------------------------------------------------------
// 	Physics::Evolution::StarContext starCtx(ns.Profile());
// 	Physics::Evolution::GeometryCache geo(starCtx);

// 	// Validate GeometryCache contents
// 	if (geo.R().Size() == 0 || geo.Exp2Nu().Size() == 0)
// 	{
// 		Z_LOG_ERROR("GeometryCache is empty (R or Exp2Nu size is 0).");
// 		return 1;
// 	}

// 	// -------------------------------------------------------
// 	// 3) Evolution configuration
// 	// -------------------------------------------------------
// 	Physics::Evolution::Config cfg;

// 	cfg.couple_spin = true;
// 	cfg.n_eta = 0;

// 	cfg.stepper = Physics::Evolution::StepperType::RKF45;
// 	cfg.rtol = 1e-6;
// 	cfg.atol = 1e-10;
// 	cfg.max_steps = 1000000;
// 	cfg.dt_save = 1.0e5;

// 	// -------------------------------------------------------
// 	// 4) EvolutionSystem::Context wiring
// 	// -------------------------------------------------------
// 	Physics::Evolution::DriverContext ctx;
// 	ctx.star = &starCtx;
// 	ctx.geo = &geo;
// 	ctx.envelope = nullptr; // later
// 	ctx.cfg = &cfg;

// 	// -------------------------------------------------------
// 	// 5) Dynamic states: Thermal + Spin
// 	// -------------------------------------------------------
// 	Physics::State::ThermalState thermal;
// 	thermal.Resize(1);
// 	thermal.SetTinf(1.0e8); // K

// 	Physics::State::SpinState spin;
// 	spin.Resize(1);
// 	spin.Omega() = 100.0; // rad/s

// 	// -------------------------------------------------------
// 	// 6) Register blocks into StateVector
// 	// -------------------------------------------------------
// 	Physics::Evolution::StateVector stateVec;
// 	stateVec.Register(Physics::State::StateTag::Thermal, thermal);
// 	stateVec.Register(Physics::State::StateTag::Spin, spin);

// 	// -------------------------------------------------------
// 	// 7) Define packing order
// 	// -------------------------------------------------------
// 	Physics::Evolution::StateLayout layout;
// 	layout.Configure(
// 		stateVec,
// 		{Physics::State::StateTag::Thermal,
// 		 Physics::State::StateTag::Spin});

// 	const std::size_t dim = layout.TotalSize();

// 	// -------------------------------------------------------
// 	// 8) RHS accumulator
// 	// -------------------------------------------------------
// 	Physics::Evolution::RHSAccumulator rhs;
// 	rhs.Configure(Physics::State::StateTag::Thermal, thermal.Size());
// 	rhs.Configure(Physics::State::StateTag::Spin, spin.Size());

// 	// -------------------------------------------------------
// 	// 9) Drivers
// 	// -------------------------------------------------------
// 	Physics::Driver::Spin::MagneticDipole::Options spinOpts;
// 	spinOpts.braking_index = 3.0;
// 	spinOpts.K_prefactor = 1e-15;
// 	spinOpts.use_moment_of_inertia = false;

// 	auto spinDriver =
// 		std::make_shared<Physics::Driver::Spin::MagneticDipole>(spinOpts);

// 	Physics::Driver::Thermal::PhotonCooling::Options thermOpts;
// 	thermOpts.surface_model =
// 		Physics::Driver::Thermal::PhotonCooling::Options::SurfaceModel::ApproxFromTinf;

// 	// PhotonCooling uses GeometryCache (4πR^2 e^{2ν}) internally when ctx.geo is set.
// 	// radiating_fraction is a *dimensionless* prefactor (e.g. hot-spot fraction).
// 	thermOpts.radiating_fraction = 1.0;
// 	thermOpts.C_eff = 1.0e40;
// 	thermOpts.global_scale = 1.0;

// 	auto thermalDriver =
// 		std::make_shared<Physics::Driver::Thermal::PhotonCooling>(thermOpts);

// 	// ----------------------------
// 	// Neutrino cooling (core)
// 	// ----------------------------
// 	Physics::Driver::Thermal::NeutrinoCooling::Options nuOpts;
// 	nuOpts.include_direct_urca = true;
// 	nuOpts.include_modified_urca = true;
// 	nuOpts.include_pair_breaking = false;
// 	nuOpts.global_scale = 1.0;

// 	auto neutrinoDriver =
// 		std::make_shared<Physics::Driver::Thermal::NeutrinoCooling>(nuOpts);
// 	// ----------------------------
// 	std::vector<Physics::Evolution::DriverPtr> drivers;
// 	drivers.push_back(spinDriver);
// 	drivers.push_back(thermalDriver);
// 	drivers.push_back(neutrinoDriver);

// 	// -------------------------------------------------------
// 	// 10) Build EvolutionSystem
// 	// -------------------------------------------------------
// 	Physics::Evolution::EvolutionSystem system(
// 		ctx,
// 		stateVec,
// 		rhs,
// 		layout,
// 		std::move(drivers));

// 	// -------------------------------------------------------
// 	// 11) Observer: Diagnostics
// 	// -------------------------------------------------------
// 	namespace PhyEvolObs = Physics::Evolution::Observers;
// 	PhyEvolObs::DiagnosticsObserver::Options dopts;
// 	dopts.output_path = base_results_dir + "/" + out_dir + "/diagnostics.jsonl";
// 	dopts.record_every_n_steps = 1000; // Record every 1000 steps
// 	dopts.record_at_start = true;	   // common choice
// 	dopts.write_catalog = true;
// 	dopts.catalog_output_path = base_results_dir + "/" + out_dir + "/diagnostics.catalog.json";

// 	std::vector<const Physics::Driver::Diagnostics::IDriverDiagnostics *> diag_drivers;

// 	if (auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(spinDriver.get()))
// 		diag_drivers.push_back(p);

// 	if (auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(thermalDriver.get()))
// 		diag_drivers.push_back(p);

// 	if (auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(neutrinoDriver.get()))
// 		diag_drivers.push_back(p);

// 	auto diag = std::make_shared<PhyEvolObs::DiagnosticsObserver>(dopts, diag_drivers);

// 	system.AddObserver(diag);
// 	// -------------------------------------------------------
// 	// 11b) Observer: TimeSeries (compact table for plotting)
// 	// -------------------------------------------------------
// 	{
// 		namespace PhyEvolObs = Physics::Evolution::Observers;
// 		using TS = PhyEvolObs::TimeSeriesObserver;

// 		TS::Options topts;
// 		topts.output_path = base_results_dir + "/" + out_dir + "/timeseries.csv";
// 		topts.format = TS::OutputFormat::CSV;
// 		topts.append = false;
// 		topts.record_at_start = true;
// 		topts.record_every_n_samples = 1; // record every OnSample (dt_save cadence)
// 		topts.record_every_dt = 0.0;	  // disable time-trigger (optional)
// 		topts.write_header = true;
// 		topts.write_sidecar_metadata = true;
// 		topts.float_precision = 17;

// 		// NEW: schema-driven columns from diagnostics catalog
// 		topts.use_catalog = true;
// 		topts.catalog_path = base_results_dir + "/" + out_dir + "/diagnostics.catalog.json";

// 		// Pick the profile(s) you want (must exist in ProducerCatalog::profiles)
// 		topts.catalog_profiles = {"timeseries_default"};

// 		// Optional convenience
// 		topts.include_builtin_time = true;
// 		topts.include_builtin_sample_index = true;

// 		// If (later) you add DriverScalar columns, pass diagnostics-capable drivers here.
// 		std::vector<const Physics::Driver::Diagnostics::IDriverDiagnostics *> ts_drivers;
// 		if (auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(spinDriver.get()))
// 			ts_drivers.push_back(p);
// 		if (auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(thermalDriver.get()))
// 			ts_drivers.push_back(p);

// 		if (auto *p = dynamic_cast<const Physics::Driver::Diagnostics::IDriverDiagnostics *>(neutrinoDriver.get()))
// 			ts_drivers.push_back(p);

// 		auto ts = std::make_shared<TS>(topts, ts_drivers);
// 		system.AddObserver(ts);
// 	}

// 	// -------------------------------------------------------
// 	// 12) Pack initial state -> y[]
// 	// -------------------------------------------------------
// 	std::vector<double> y(dim);
// 	Physics::Evolution::PackStateVector(stateVec, layout, y.data());

// 	std::cout << "Initial conditions:\n";
// 	std::cout << "  Tinf  = " << thermal.Tinf() << "\n";
// 	std::cout << "  Omega = " << spin.Omega() << "\n";

// 	// -------------------------------------------------------
// 	// 13) Integrate
// 	// -------------------------------------------------------
// 	Physics::Evolution::GSLIntegrator integrator(system, cfg, dim);

// 	const double t0 = 0.0;
// 	const double t1 = 1.0e10; // s

// 	const bool ok = integrator.Integrate(t0, t1, y.data());
// 	if (!ok)
// 	{
// 		std::cerr << "GSLIntegrator: integration failed or max_steps exceeded.\n";
// 		return 1;
// 	}

// 	// -------------------------------------------------------
// 	// 14) Unpack final state and print
// 	// -------------------------------------------------------
// 	Physics::Evolution::UnpackStateVector(stateVec, layout, y.data());

// 	std::cout << "Final conditions (t = " << t1 << " s):\n";
// 	std::cout << "  Tinf  = " << thermal.Tinf() << "\n";
// 	std::cout << "  Omega = " << spin.Omega() << "\n";

// 	std::cout << "[debug] done.\n";
// 	return 0;
// }

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

	Zaki::Util::Instrumentor::BeginSession("spin-therm-evol-2-main", base_results_dir + out_dir + "/timing_profile.json");

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

	const double target_M = 1.8; // Msun
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
	// 3) Evolution configuration (defaults + overrides)
	// -------------------------------------------------------
	Physics::Evolution::Config cfg = Physics::Evolution::Run::MakeDefaultConfig();
	cfg.couple_spin = true;
	cfg.n_eta = 0;

	// (optional overrides)
	cfg.rtol = 1e-6;
	cfg.atol = 1e-10;
	cfg.max_steps = 1'000'000;
	cfg.dt_save = 1.0e6;

	// -------------------------------------------------------
	// 4) DriverContext wiring
	// -------------------------------------------------------
	// Physics::Evolution::DriverContext ctx =
	// 	Physics::Evolution::Run::MakeDriverContext(starCtx, geo, cfg);
	Physics::Driver::Thermal::Boundary::EnvelopePotekhin1997_Iron env97_fe;

	Physics::Evolution::DriverContext ctx =
		Physics::Evolution::Run::MakeDriverContext(starCtx, geo, cfg, &env97_fe);

	// -------------------------------------------------------
	// 5) Dynamic states: Thermal + Spin
	// -------------------------------------------------------
	Physics::State::ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(1.0e8); // K

	Physics::State::SpinState spin;
	spin.Resize(1);
	spin.Omega() = 100.0; // rad/s

	// -------------------------------------------------------
	// 6) State wiring (StateVector + Layout + RHS) via RunBuilder helpers
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
	// 7) Drivers
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
	photonOpts.C_eff = 1.0e40;
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
	// 8) Build EvolutionSystem
	// -------------------------------------------------------
	Physics::Evolution::EvolutionSystem system(
		ctx,
		wiring.state_vec,
		wiring.rhs,
		wiring.layout,
		std::move(drivers));

	// -------------------------------------------------------
	// 9) Observers (default presets + factories)
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
	// 10) Pack initial state -> y[]
	// -------------------------------------------------------
	std::vector<double> y(wiring.dim);
	Physics::Evolution::PackStateVector(wiring.state_vec, wiring.layout, y.data());

	std::cout << "Initial conditions:\n";
	std::cout << "  Tinf  = " << thermal.Tinf() << " K\n";
	std::cout << "  Omega = " << spin.Omega() << " rad/s\n";

	// -------------------------------------------------------
	// 11) Integrate
	// -------------------------------------------------------
	Physics::Evolution::GSLIntegrator integrator(system, cfg, wiring.dim);

	const double t0 = 0.0;
	const double t1 = 1.0e9;

	const bool ok = integrator.Integrate(t0, t1, y.data());
	if (!ok)
	{
		std::cerr << "GSLIntegrator: integration failed or max_steps exceeded.\n";
		return 1;
	}

	// -------------------------------------------------------
	// 12) Unpack final state and print
	// -------------------------------------------------------
	Physics::Evolution::UnpackStateVector(wiring.state_vec, wiring.layout, y.data());

	std::cout << "Final conditions (t = " << t1 << " s):\n";
	std::cout << "  Tinf  = " << thermal.Tinf() << " K\n";
	std::cout << "  Omega = " << spin.Omega() << " rad/s\n";

	std::cout << "[debug] done.\n";

	Zaki::Util::Instrumentor::EndSession();
	return 0;
}
//==============================================================