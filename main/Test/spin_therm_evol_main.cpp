// -*- lsst-c++ -*-
/*
 * CompactStar test main: coupled Spin + Thermal evolution
 *
 * Uses:
 *  - SpinState + MagneticDipole driver  (dΩ/dt)
 *  - ThermalState + PhotonCooling driver (dT∞/dt)
 *
 * This is primarily an infrastructure test:
 *  - multi-block StateLayout ({Thermal, Spin}),
 *  - multiple drivers writing into RHSAccumulator,
 *  - packing/unpacking via StatePacking,
 *  - GSLIntegrator for time-stepping.
 */

#include <iostream>
#include <memory>
#include <vector>

#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp"
#include "CompactStar/Physics/Evolution/RHSAccumulator.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StateLayout.hpp"
#include "CompactStar/Physics/Evolution/StatePacking.hpp"

#include "CompactStar/Physics/State/SpinState.hpp"
#include "CompactStar/Physics/State/Tags.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

#include "CompactStar/Physics/Driver/Spin/MagneticDipole.hpp"
#include "CompactStar/Physics/Driver/Thermal/PhotonCooling.hpp"

#include "CompactStar/Core/StarBuilder.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

int main()
{
	using namespace CompactStar;

	// --------------------------------------------------------------
	// 1) Solve TOV to get a neutron-star profile
	// --------------------------------------------------------------
	Zaki::String::Directory dir(__FILE__);
	std::cout << "[debug] this file dir = " << dir << "\n";

	// Create a TOV solver
	CompactStar::Core::TOVSolver tov;

	Zaki::String::Directory eos_root =
		dir.ParentDir().ParentDir() + "/EOS/CompOSE/";
	// Zaki::String::Directory eos_root = "EOS/CompOSE/";
	std::string eos_name = "DS(CMF)-1_with_crust";

	tov.ImportEOS(eos_root + eos_name + "/" + eos_name + ".eos");
	tov.SetWrkDir(dir.ParentDir() + "/results");
	Zaki::Math::Axis ec_axis({{5.0e+14, 1.913e15}, 20, "Log"});

	// Solve the sequence
	//    Output directory: results/out_dir/
	//    Output file name: out_file_Sequence.tsv
	Zaki::String::Directory out_dir = "spin_therm_evol/";
	Zaki::String::Directory out_file = "spin_therm_evol";

	tov.Solve(ec_axis, out_dir, out_file);

	Core::StarBuilder::Options sbOpts; // defaults are fine for now
	Core::StarBuilder::Output sbOut;
	Core::StarBuilder::BuildFromSequence(tov.GetWrkDir(), out_dir, out_file.Str(), 1.8, sbOut, sbOpts);

	Core::NStar ns;

	Physics::Evolution::StarContext starCtx;
	Physics::Evolution::GeometryCache geo(starCtx);
	// --------------------------------------------------------------
	// 1) Evolution configuration
	// --------------------------------------------------------------
	Physics::Evolution::Config cfg;

	// Enable spin in the state vector (so MagneticDipole has something to act on).
	cfg.couple_spin = true;

	// No chemical imbalances for this test.
	cfg.n_eta = 0;

	// Simple integrator settings for a short test run.
	cfg.stepper = Physics::Evolution::StepperType::RKF45;
	cfg.rtol = 1e-6;
	cfg.atol = 1e-10;
	cfg.max_internal_steps = 1000000;
	cfg.max_samples = 1000000;
	cfg.dt_save = 1.0e5; // not used directly by GSLIntegrator yet, but kept for consistency

	// --------------------------------------------------------------
	// 2) Static context (star, geometry, envelope)
	// --------------------------------------------------------------
	//
	// For now we don't need any of these in the drivers we’re using, so we
	// can leave them as nullptr. Later, when we wire in real microphysics
	// and envelope models, we’ll fill these.
	Physics::Evolution::DriverContext ctx;
	ctx.star = &starCtx;	// e.g. pointer to StarContext built from NStar + EOS
	ctx.geo = &geo;			// e.g. GeometryCache
	ctx.envelope = nullptr; // e.g. Thermal::IEnvelope implementation
	ctx.cfg = &cfg;

	// --------------------------------------------------------------
	// 3) Dynamic state blocks: Thermal + Spin
	// --------------------------------------------------------------

	// --- Thermal block: 1 DOF = T∞ ---------------------------------
	Physics::State::ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(1.0e8); // K, rough initial redshifted internal temperature

	// Optional: set T_surf if we use DirectTSurf model later.
	// For now we use ApproxFromTinf in PhotonCooling::Options, so this is not needed.
	// thermal.T_surf = 1.0e6;

	// --- Spin block: 1 DOF = Ω -------------------------------------
	Physics::State::SpinState spin;
	spin.Resize(1);
	spin.Omega() = 100.0; // initial spin frequency [rad/s], arbitrary

	// --------------------------------------------------------------
	// 4) StateVector: register both blocks with tags
	// --------------------------------------------------------------
	Physics::Evolution::StateVector stateVec;
	stateVec.Register(Physics::State::StateTag::Thermal, thermal);
	stateVec.Register(Physics::State::StateTag::Spin, spin);

	// --------------------------------------------------------------
	// 5) Layout: decide packing order in y[]
	// --------------------------------------------------------------
	//
	// Here we choose y = [ Thermal block, then Spin block ].
	// We could swap the order if we prefer.
	Physics::Evolution::StateLayout layout;
	layout.Configure(
		stateVec,
		{Physics::State::StateTag::Thermal,
		 Physics::State::StateTag::Spin});

	const std::size_t dim = layout.TotalSize();

	// --------------------------------------------------------------
	// 6) RHSAccumulator: configure blocks for both tags
	// --------------------------------------------------------------
	Physics::Evolution::RHSAccumulator rhs;
	rhs.Configure(Physics::State::StateTag::Thermal, thermal.Size());
	rhs.Configure(Physics::State::StateTag::Spin, spin.Size());

	// --------------------------------------------------------------
	// 7) Drivers: MagneticDipole (Spin) + PhotonCooling (Thermal)
	// --------------------------------------------------------------

	// --- Spin driver: MagneticDipole --------------------------------
	Physics::Driver::Spin::MagneticDipole::Options spinOpts;
	spinOpts.braking_index = 3.0;
	spinOpts.K_prefactor = 1e-15; // arbitrary; tweak to see visible evolution
	spinOpts.use_moment_of_inertia = false;

	auto spinDriver =
		std::make_shared<Physics::Driver::Spin::MagneticDipole>(spinOpts);

	// --- Thermal driver: PhotonCooling ------------------------------
	Physics::Driver::Thermal::PhotonCooling::Options thermOpts;

	// Use simple T_surf ≈ T∞ mapping for now.
	thermOpts.surface_model =
		Physics::Driver::Thermal::PhotonCooling::Options::SurfaceModel::ApproxFromTinf;

	// Effective area and heat capacity are just toy values here; adjust to keep
	// dT/dt in a reasonable range.
	thermOpts.radiating_fraction = 1.0; // in code units; real code would use ~4πR^2 z^2
	thermOpts.C_eff = 1.0e40;			// big C_eff → slow cooling for demonstration
	thermOpts.global_scale = 1.0;

	auto thermalDriver =
		std::make_shared<Physics::Driver::Thermal::PhotonCooling>(thermOpts);

	// --- Driver list -------------------------------------------------
	std::vector<Physics::Evolution::DriverPtr> drivers;
	drivers.push_back(spinDriver);
	drivers.push_back(thermalDriver);

	// --------------------------------------------------------------
	// 8) EvolutionSystem: glue together context, state, RHS, layout, drivers
	// --------------------------------------------------------------
	Physics::Evolution::EvolutionSystem system(
		ctx,
		stateVec,
		rhs,
		layout,
		std::move(drivers));

	// --------------------------------------------------------------
	// 9) Pack initial state into flat y[]
	// --------------------------------------------------------------
	std::vector<double> y(dim);
	Physics::Evolution::PackStateVector(stateVec, layout, y.data());

	std::cout << "Initial conditions:\n";
	std::cout << "  Tinf  = " << thermal.Tinf() << "\n";
	std::cout << "  Omega = " << spin.Omega() << "\n";

	// --------------------------------------------------------------
	// 10) GSL integrator and evolve
	// --------------------------------------------------------------
	Physics::Evolution::GSLIntegrator integrator(system, cfg, dim);

	const double t0 = 0.0;
	const double t1 = 1.0e10; // s; adjust as desired

	const bool ok = integrator.Integrate(t0, t1, y.data());
	if (!ok)
	{
		std::cerr << "GSLIntegrator: integration failed or max_steps exceeded.\n";
		return 1;
	}

	// --------------------------------------------------------------
	// 11) Unpack final state back into blocks and print
	// --------------------------------------------------------------
	Physics::Evolution::UnpackStateVector(stateVec, layout, y.data());

	std::cout << "Final conditions (t = " << t1 << " s):\n";
	std::cout << "  Tinf  = " << thermal.Tinf() << "\n";
	std::cout << "  Omega = " << spin.Omega() << "\n";

	return 0;
}