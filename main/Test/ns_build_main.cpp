// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 *
 * Debug / smoke-test main for the TOV solver.
 *
 * Goal:
 *  - make sure TOVSolver.hpp / TOVSolver.cpp actually link
 *  - exercise: ImportEOS -> Solve(...) -> ExportSequence(...)
 *  - touch both NS and (optionally) mixed-star paths
 *
 * This file is intentionally verbose and defensive: it prints
 * out what it’s doing so we can see where it fails (path,
 * EOS import, GSL spline, surface reached, export, ...).
 */

#include <gsl/gsl_errno.h>
#include <iostream>

#include <Zaki/String/Directory.hpp>
#include <Zaki/Util/Logger.hpp>

#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/Physics/Evolution/EvolutionConfig.hpp>
#include <CompactStar/Physics/Evolution/EvolutionSystem.hpp>
#include <CompactStar/Physics/Evolution/GeometryCache.hpp>
#include <CompactStar/Physics/Evolution/StarContext.hpp>
//--------------------------------------------------------------
// GSL error handler (so GSL does not abort the program)
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
//--------------------------------------------------------------

int main()
{
	// 1) GSL error handler
	gsl_set_error_handler(&my_gsl_error_handler);

	// 2) Logging: make it a bit chatty for debugging
	Zaki::Util::LogManager::SetLogLevels(Zaki::Util::LogLevel::Info);

	// 3) Working directory based on this file
	//    dir      -> .../CompactStar/Core/main/Test/
	//    dir.ParentDir() -> .../CompactStar/Core/main/
	Zaki::String::Directory dir(__FILE__);
	std::cout << "[debug] this file dir = " << dir << "\n";
	std::cout << "[debug] this file parent dir = " << dir.ParentDir() << "\n";

	Zaki::String::Directory eos_root =
		dir.ParentDir().ParentDir() + "/EOS/CompOSE/";
	std::string eos_name = "DS(CMF)-1_with_crust";

	std::cout << "[debug] assuming EOS root = " << eos_root << "\n";
	std::cout << "[debug] assuming EOS name = " << eos_name << "\n";

	// Build the EOS file path:
	Zaki::String::Directory eos_file =
		eos_root + eos_name + "/" + eos_name + ".eos";

	Zaki::String::Directory base_results_dir = dir.ParentDir() + "/results";

	// 8) Solve the sequence
	//    Output directory: results/ns_build/
	Zaki::String::Directory out_dir = "ns_build";

	CompactStar::Core::NStar ns;
	ns.SetWrkDir(base_results_dir);
	const double target_M = 1.8; // Msun
	const int n = ns.SolveTOV_Profile(eos_file, target_M, out_dir);

	ns.Export(out_dir + "/NStar_Profile.tsv");

	CompactStar::Physics::Evolution::StarContext starCtx(ns.Profile());
	CompactStar::Physics::Evolution::GeometryCache geo(starCtx);
	std::cout << "\n geo.Exp2Nu().Size() = " << geo.Exp2Nu().Size() << "\n";
	std::cout << "[debug] done.\n";

	// --------------------------------------------------------------
	// 1) Evolution configuration
	// --------------------------------------------------------------
	CompactStar::Physics::Evolution::Config cfg;

	// Enable spin in the state vector (so MagneticDipole has something to act on).
	cfg.couple_spin = true;

	// No chemical imbalances for this test.
	cfg.n_eta = 0;

	// Simple integrator settings for a short test run.
	cfg.stepper = CompactStar::Physics::Evolution::StepperType::RKF45;
	cfg.rtol = 1e-6;
	cfg.atol = 1e-10;
	cfg.max_internal_steps = 1000000;
	cfg.max_samples = 1000000;
	cfg.dt_save = 1.0e5; // not used directly by GSLIntegrator yet, but kept for consistency

	CompactStar::Physics::Evolution::DriverContext ctx;
	ctx.star = &starCtx; // e.g. pointer to StarContext built from NStar + EOS
	ctx.geo = &geo;		 // e.g. GeometryCache
	// ctx.envelope = nullptr; // e.g. Thermal::IEnvelope implementation
	ctx.cfg = &cfg;

	return 0;
}