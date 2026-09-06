#include "tests/relativity/candidate_capture.hpp"
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file passive_cooling_regression.cpp
 * @brief Phase-2B passive-cooling regression baseline (GOVERNANCE.md §3.1 condition 7).
 *
 * Freezes the end-to-end behavior of the now-coherent thermal equation
 *
 *     C_*(T_inf) dT_inf/dt = -L_nu,inf - L_gamma,inf
 *
 * with BOTH channels dividing by the same canonical C_*(T_inf) (ADR-0002 Pattern A).
 *
 * ============================ WHAT THIS IS NOT ============================
 * This is a REGRESSION baseline, not a validation of the physics. The neutrino
 * emissivity normalizations it freezes are self-labelled placeholders in source
 * (NeutrinoCooling_Details.cpp: Q0_DU = 1.0e27, Q0_MU = 1.0e21, "Placeholder
 * normalizations (must match your emissivity model)"). A future scientifically
 * justified change to them is EXPECTED to move these values and must be
 * separately reviewed and re-baselined.
 * =========================================================================
 *
 * Uses only production APIs on authenticated external CMF data; no mocked star,
 * no mocked heat capacity, no simplified cooling law, no historical result file.
 *
 * Usage:
 *   passive_cooling_regression <EOS_DATA_ROOT> --compare <baseline.tsv>
 *   passive_cooling_regression <EOS_DATA_ROOT> --emit    <baseline.tsv>
 *   passive_cooling_regression <EOS_DATA_ROOT> --study            (numerics study)
 */

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <chrono>
#include <map>
#include <sstream>
#include <set>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/EOS/CompOSE_Thermo.hpp"
#include "CompactStar/Physics/Driver/Spin/MagneticDipole.hpp"
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
#include "CompactStar/Physics/State/SpinState.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;
namespace P = CompactStar::Physics;
namespace Th = CompactStar::Physics::Driver::Thermal;
using CompactStar::EOS::CompOSE_Thermo;

// ---------------------------------------------------------------------------
//  Frozen configuration — authenticated from main/Test/spin_therm_evol_2_main.cpp
//  at commit 56aee7c. Any change here changes the baseline by definition.
// ---------------------------------------------------------------------------
constexpr double kTargetM = 1.6;   // Msun
constexpr double kTinf0_K = 1.0e9; // K
constexpr double kT0_yr = 1.0e2;
constexpr double kT1_yr = 1.0e6;
constexpr double kRtol = 1e-6;
constexpr double kAtol = 1e-10;
constexpr double kSamplesPerDecade = 150.0;
static const char *kSourceCommit = "94165f359ccfadf34bf64ede62e7bef9c581a067 (ADR-0009 validated science)";
// Selected from the Phase-2B-1R convergence study; see PASSIVE_COOLING_BASELINE.md.
constexpr P::Evolution::StepperType kBaselineStepper = P::Evolution::StepperType::RKF45;

// Logarithmically separated diagnostic checkpoints [yr].
// ---------------------------------------------------------------------------
//  Tolerance policy — DERIVED FROM MEASUREMENT, not chosen to make the first run pass.
//
//  Inputs measured in the Phase-2B-1R study (--study), RKF45 / thermal-only:
//     repeat-run variation                 0            (bit-identical over 3 runs)
//     nominal 1e-6 vs tighter 3e-7 / 1e-7  1.19e-7      (T_inf)
//     cadence 150 -> 75 samples/decade     2.04e-6      (T_inf),  1.22e-5 (luminosity)
//     floating-point floor                 1e-14
//  driver = max(...) = 2.04e-6 for state, 1.22e-5 for luminosities.
//  Fixed safety factor x5, chosen in advance and applied uniformly:
//     state       5 x 2.04e-6 = 1.02e-5  -> 1e-5
//     luminosity  5 x 1.22e-5 = 6.1e-5   -> 1e-4
//  Both are far tighter than percent level, so the baseline is numerically mature.
//  See docs/validation/PASSIVE_COOLING_BASELINE.md.
// ---------------------------------------------------------------------------
constexpr double kTolState = 1e-5; // Tinf, C_star, Tsurf, Tb, dLnTinf/dt
constexpr double kTolLumin = 1e-4; // luminosities (steep powers of T amplify dT)
constexpr double kTolEnergyIdentity = 1e-10;
constexpr double kTolPatternA = 1e-12; // C_star must be IDENTICAL across both channels

static const std::vector<double> kCheckpointsYr = {
	1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6};

struct Checkpoint
{
	double t_s = 0, t_yr = 0, Tinf_K = 0, C_star_erg_K = 0;
	double L_nu_inf = 0, L_nu_DU = 0, L_nu_MU = 0, L_gamma_inf = 0;
	double dLnTinf_dt = 0, Tsurf_K = 0, Tb_K = 0;
	double energy_identity_rel = 0; // |RHS_sum - (-(Lnu+Lg)/C)/T| / |expected|
	double C_star_match_rel = 0;    // |C_nu - C_gamma| / C_gamma  (ADR-0002 Pattern A)
};

struct StarInfo
{
	double requested_M = 0, achieved_M = 0, R_km = 0, ec = 0;
	int n_rows = 0;
	std::uint64_t prof_version = 0;
};

// ---------------------------------------------------------------------------
//  Phase-2B-3 (INV-12): observations proving the canonical run never enters a
//  known cache-hazard state. Everything here is recorded test-side, from the
//  DriverContext the observer already receives; no production API was added.
// ---------------------------------------------------------------------------
struct CacheSafety
{
	std::uint64_t version_at_start = 0;
	std::uint64_t version_at_finish = 0;
	std::size_t distinct_versions_during_run = 0;
	std::size_t distinct_geo_ptrs = 0;
	std::size_t distinct_star_ptrs = 0;
	std::size_t distinct_thermo_ptrs = 0;
	std::size_t observations = 0;
	bool geo_is_the_run_geo = false;   // ctx.geo == &geo constructed from this run's StarContext
	bool star_is_the_run_star = false; // ctx.star == &starCtx
};

struct RunOpts
{
	double rtol = kRtol;
	double atol = kAtol;
	double spd = kSamplesPerDecade;
	P::Evolution::StepperType stepper = kBaselineStepper; // always explicit
	bool segmented = false; // false = one continuous Integrate(t0,t1) — the canonical procedure
	bool include_spin = false; // Configuration B (thermal-only) is the canonical baseline
	// Detector proof only (--detector). NEVER used by the canonical run.
	double photon_global_scale = 1.0;
};

// ---------------------------------------------------------------------------
static bool RunTrajectory(const fs::path &root, const RunOpts &ro,
						  std::vector<Checkpoint> &out, StarInfo &si,
						  std::string &err, CacheSafety *cs = nullptr)
{
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	const fs::path fineT = root / "DNS-CMF-Hadronic-with-electrons";
	if (!fs::exists(cold) || !fs::exists(fineT / "eos.thermo"))
	{
		err = "authenticated EOS data missing under " + root.string();
		return false;
	}

	CompactStar::Core::NStar ns;
	ns.SetWrkDir(fs::temp_directory_path().string());
	const int n_rows = ns.SolveTOV_Profile(cold.string(), kTargetM, "pc_regression");
	if (n_rows <= 0)
	{
		err = "SolveTOV_Profile failed";
		return false;
	}
	const auto &seq = ns.GetSequence();
	si = {kTargetM, seq.m, seq.r, seq.ec, n_rows, 0};

	P::Evolution::StarContext starCtx(ns.Profile());
	P::Evolution::GeometryCache geo(starCtx);
	si.prof_version = starCtx.ProfileVersion();

	CompOSE_Thermo::Options thOpt;
	thOpt.use_central_difference = true;
	thOpt.clamp_to_domain = true;
	thOpt.Tmin_for_derivative_MeV = 0.0;
	CompOSE_Thermo thermo(fineT.string(), thOpt);
	if (!thermo.IsLoaded())
	{
		err = "CompOSE_Thermo failed to load";
		return false;
	}

	P::Evolution::Config cfg;
	cfg.couple_spin = ro.include_spin;
	cfg.n_eta = 0;
	cfg.rtol = ro.rtol;
	cfg.atol = ro.atol;
	cfg.max_internal_steps = 1'000'000;
	cfg.max_samples = 1'000'000;
	cfg.stepper = ro.stepper; // explicit: the canonical run must not rely on the default
	cfg.save_cadence = P::Evolution::SaveCadence::LogTime;
	cfg.samples_per_decade = ro.spd;

	// Exactly the live wiring from spin_therm_evol_2_main.cpp:188-192 — including the
	// Potekhin-1997 iron envelope, which the Tb->Ts mapping needs.
	Th::Boundary::EnvelopePotekhin1997_Iron env97_fe;
	P::Evolution::DriverContext ctx =
		P::Evolution::Run::MakeDriverContext(starCtx, geo, cfg, &env97_fe, &thermo);

	P::State::ThermalState thermal;
	thermal.Resize(1);
	thermal.SetTinf(kTinf0_K);
	P::State::SpinState spin;
	spin.Resize(1);
	spin.Omega() = 100.0;

	P::Evolution::Run::StateWiring w;
	std::vector<P::State::StateTag> tags{P::State::StateTag::Thermal};
	w.state_vec.Register(P::State::StateTag::Thermal, thermal);
	if (ro.include_spin)
	{
		w.state_vec.Register(P::State::StateTag::Spin, spin);
		tags.push_back(P::State::StateTag::Spin);
	}
	P::Evolution::Run::ConfigureLayout(w, tags);
	P::Evolution::Run::ConfigureRHS(w, tags);

	Th::PhotonCooling::Options po;
	po.surface_model = Th::PhotonCooling::Options::SurfaceModel::EnvelopeTbTs;
	po.radiating_fraction = 1.0;
	po.global_scale = ro.photon_global_scale; // 1.0 canonical; perturbed only by --detector
	auto photon = std::make_shared<Th::PhotonCooling>(po);

	Th::NeutrinoCooling::Options no;
	no.include_direct_urca = true;
	no.include_modified_urca = true;
	no.include_pair_breaking = false;
	no.global_scale = 1.0;
	auto neutrino = std::make_shared<Th::NeutrinoCooling>(no);

	std::vector<P::Evolution::DriverPtr> drivers;
	if (ro.include_spin)
	{
		P::Driver::Spin::MagneticDipole::Options so;
		so.braking_index = 3.0;
		so.K_prefactor = 1e-15;
		so.use_moment_of_inertia = false;
		drivers.push_back(std::make_shared<P::Driver::Spin::MagneticDipole>(so));
	}
	drivers.push_back(photon);
	drivers.push_back(neutrino);

	P::Evolution::EvolutionSystem system(ctx, w.state_vec, w.rhs, w.layout,
										 std::move(drivers));

	// -----------------------------------------------------------------------
	//  Test-side observer: capture checkpoints DURING one continuous integration.
	//  Production observer architecture is untouched.
	// -----------------------------------------------------------------------
	struct Capture final : P::Evolution::Observers::IObserver
	{
		const Th::PhotonCooling *ph;
		const Th::NeutrinoCooling *nu;
		std::vector<Checkpoint> *out;
		std::string *err;
		std::size_t next = 0;

		// INV-12 canonical-path observations. Recorded on EVERY callback, not just
		// at the checkpoints we keep, so a mid-run swap could not slip through.
		std::set<const void *> geo_seen, star_seen, thermo_seen;
		std::set<std::uint64_t> version_seen;
		std::size_t n_obs = 0;

		void Observe(const P::Evolution::DriverContext &c)
		{
			++n_obs;
			geo_seen.insert(static_cast<const void *>(c.geo));
			star_seen.insert(static_cast<const void *>(c.star));
			thermo_seen.insert(static_cast<const void *>(c.thermo));
			if (c.star)
				version_seen.insert(c.star->ProfileVersion());
		}

		void Record(double t_s, const P::Evolution::StateVector &Y,
					const P::Evolution::DriverContext &c)
		{
			if (next >= kCheckpointsYr.size())
				return;
			const double target = kCheckpointsYr[next] * Zaki::Physics::YR_2_SEC;
			// Geometric LogTime sampling lands on decade boundaries exactly when
			// t0 is a decade and samples_per_decade divides evenly; accept a small
			// relative window so we never restart the integrator to hit a time.
			if (t_s < target * (1.0 - 1e-9))
				return;

			const auto pd = Th::Detail::ComputeDerived(*ph, Y, c);
			const auto nd = Th::Detail::NeutrinoCooling_Details::ComputeDerived(*nu, Y, c);
			if (!pd.ok || !nd.ok)
			{
				*err = "driver diagnostics not ok at t=" + std::to_string(t_s) +
					   " s: [photon] " + pd.message + " [neutrino] " + nd.message;
				next = kCheckpointsYr.size();
				return;
			}

			Checkpoint ck;
			ck.t_s = t_s;
			ck.t_yr = t_s / Zaki::Physics::YR_2_SEC;
			ck.Tinf_K = pd.Tinf_K;
			ck.C_star_erg_K = pd.C_star_erg_K;
			ck.L_nu_inf = nd.L_nu_inf_erg_s;
			ck.L_nu_DU = nd.L_nu_DU_inf_erg_s;
			ck.L_nu_MU = nd.L_nu_MU_inf_erg_s;
			ck.L_gamma_inf = pd.L_gamma_inf_erg_s;
			ck.dLnTinf_dt = pd.dLnTinf_dt_1_s + nd.dLnTinf_dt_1_s;
			ck.Tsurf_K = pd.Tsurf_K;
			ck.Tb_K = pd.Tb_K;

			const double dT_exp = -(ck.L_nu_inf + ck.L_gamma_inf) / ck.C_star_erg_K;
			const double dLn_exp = dT_exp / ck.Tinf_K;
			ck.energy_identity_rel = std::fabs(ck.dLnTinf_dt - dLn_exp) / std::fabs(dLn_exp);

			// ADR-0002 Pattern A: both channels MUST share one denominator.
			ck.C_star_match_rel =
				std::fabs(nd.C_eff_erg_K - pd.C_star_erg_K) / pd.C_star_erg_K;

			out->push_back(ck);
			++next;
		}

		void OnStart(const P::Evolution::Observers::RunInfo &r,
					 const P::Evolution::StateVector &Y,
					 const P::Evolution::DriverContext &c) override
		{
			Observe(c);
			Record(r.t0, Y, c);
		}
		void OnSample(const P::Evolution::Observers::SampleInfo &s_,
					  const P::Evolution::StateVector &Y,
					  const P::Evolution::DriverContext &c) override
		{
			Observe(c);
			Record(s_.t, Y, c);
		}
		std::string Name() const override { return "PassiveCoolingCapture"; }
	};

	out.clear();
	auto cap = std::make_shared<Capture>();
	cap->ph = photon.get();
	cap->nu = neutrino.get();
	cap->out = &out;
	cap->err = &err;
	system.AddObserver(cap);

	P::Evolution::GSLIntegrator integrator(system, cfg, w.dim);
	std::vector<double> y(w.dim);
	P::Evolution::PackStateVector(w.state_vec, w.layout, y.data());

	// INV-12: the profile must not change across the integration.
	const std::uint64_t version_at_start = starCtx.ProfileVersion();

	const double YR = Zaki::Physics::YR_2_SEC;
	const double t0 = kT0_yr * YR;
	const double t1 = kT1_yr * YR;

	if (!ro.segmented)
	{
		// CANONICAL: one continuous integration. The integrator keeps its adaptive
		// step history for the whole trajectory.
		if (!integrator.Integrate(t0, t1, y.data()))
		{
			err = "integration failed (continuous)";
			return false;
		}
	}
	else
	{
		// Diagnostic only: restart at every checkpoint, which resets adaptive step
		// history and the initial-step heuristic. Used once to quantify the
		// difference against the continuous procedure.
		double t_cur = t0;
		for (double ck_yr : kCheckpointsYr)
		{
			const double t_ck = ck_yr * YR;
			if (t_ck > t_cur)
			{
				if (!integrator.Integrate(t_cur, t_ck, y.data()))
				{
					err = "integration failed (segmented) before " +
						  std::to_string(ck_yr) + " yr";
					return false;
				}
				t_cur = t_ck;
			}
		}
	}

	if (!err.empty())
		return false;
	if (out.size() != kCheckpointsYr.size())
	{
		err = "captured " + std::to_string(out.size()) + " checkpoints, expected " +
			  std::to_string(kCheckpointsYr.size());
		return false;
	}
	if (cs)
	{
		cs->version_at_start = version_at_start;
		cs->version_at_finish = starCtx.ProfileVersion();
		cs->distinct_versions_during_run = cap->version_seen.size();
		cs->distinct_geo_ptrs = cap->geo_seen.size();
		cs->distinct_star_ptrs = cap->star_seen.size();
		cs->distinct_thermo_ptrs = cap->thermo_seen.size();
		cs->observations = cap->n_obs;
		cs->geo_is_the_run_geo =
			cap->geo_seen.size() == 1 &&
			*cap->geo_seen.begin() == static_cast<const void *>(&geo);
		cs->star_is_the_run_star =
			cap->star_seen.size() == 1 &&
			*cap->star_seen.begin() == static_cast<const void *>(&starCtx);
	}

	return true;
}

// ---------------------------------------------------------------------------
static const char *kCols =
	"t_yr\tTinf_K\tC_star_erg_K\tL_nu_inf_erg_s\tL_nu_DU_inf_erg_s\t"
	"L_nu_MU_inf_erg_s\tL_gamma_inf_erg_s\tTsurf_K\tTb_K\tdLnTinf_dt_1_s\t"
	"Lnu_over_Lgamma";

static void WriteBaseline(const fs::path &p, const std::vector<Checkpoint> &cks,
						  const StarInfo &si)
{
	std::ofstream o(p);
	o << std::setprecision(12) << std::scientific;
	o << "# CompactStar passive-cooling regression baseline\n"
	  << "# schema_version\t2\n"
	  << "#\n"
	  << "# ============================ NOT A PHYSICS VALIDATION ======================\n"
	  << "# REGRESSION baseline only. The neutrino emissivity normalizations frozen here\n"
	  << "# are SELF-LABELLED PLACEHOLDERS in source: Q0_DU = 1.0e27, Q0_MU = 1.0e21\n"
	  << "# (NeutrinoCooling_Details.cpp, \"Placeholder normalizations\"). A scientifically\n"
	  << "# justified change to them is EXPECTED to move every value below and must be\n"
	  << "# separately reviewed and re-baselined. Nothing here validates neutrino physics.\n"
	  << "# ===========================================================================\n"
	  << "#\n"
	  << "# --- provenance ---\n"
	  << "# source_commit\t" << kSourceCommit << "\n"
	  << "# thermal_configuration_origin\t800f24522bf7fc56387c2288e38c7906c1b54fcc (+ Phase-2B-1R repair)\n"
	  << "# build_configuration\tDebug\n"
	  << "# compiler\tAppleClang 17.0.0.17000604\n"
	  << "# gsl_version\t2.7.1\n"
	  << "#\n"
	  << "# --- numerical method (the metadata hole the Jan-2026 run exposed) ---\n"
	  << "# stepper\tRKF45\n"
	  << "# integration\tsingle continuous Integrate(t0,t1)\n"
	  << "# rtol\t" << kRtol << "\n"
	  << "# atol\t" << kAtol << "\n"
	  << "# save_cadence\tLogTime\n"
	  << "# samples_per_decade\t" << kSamplesPerDecade << "\n"
	  << "#\n"
	  << "# --- EOS provenance (dual representation) ---\n"
	  << "# eos_structural\tCMF #1 EoS with crust (official CompOSE), PROCESSED .eos\n"
	  << "# sha256_processed_eos\t5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd\n"
	  << "# sha256_cold_raw_eos_thermo\t416444999ccac569e2c9b34808888949c36d759f30cce25dab0d42c13e900ce3\n"
	  << "# sha256_cold_raw_eos_nb\td9c8e78c2fcf37fe770fecfc2d3a211d840a28299821a56c77e66f9ff74edef8\n"
	  << "# sha256_cold_raw_eos_t\t1a37b9563c40962b203e7bca1aa3b41e8c8b1427953df68095a51dd2cc17ff96\n"
	  << "# sha256_cold_raw_eos_yq\t1a37b9563c40962b203e7bca1aa3b41e8c8b1427953df68095a51dd2cc17ff96\n"
	  << "# eos_thermo\tCMF hadronic EoS with electrons (official CompOSE), RAW tables\n"
	  << "# sha256_thermo_eos_thermo\ta456fb8595208ddf3119350a856fbf2b906c0a0e19bb7c716571748d0aa0724b\n"
	  << "# sha256_thermo_eos_nb\t3f79dbcc6f8b519696377f89ebc86464bc55cd61d9e2459f6e21e2d9e00f380d\n"
	  << "# sha256_thermo_eos_t\t2e4c6ec1feb85b16d0ee7036dce183782a9f681577e79c72315171069aa8513d\n"
	  << "# sha256_thermo_eos_yq\td98fcd2f7752039c552c2ef2d04ab485b75db47a61f8ae1740875b54bf9824fd\n"
	  << "#\n"
	  << "# --- star (structural fingerprint) ---\n"
	  << "# requested_mass_Msun\t" << kTargetM << "\n"
	  << "# achieved_mass_Msun\t" << si.achieved_M << "\n"
	  << "# radius_km\t" << si.R_km << "\n"
	  << "# central_eps\t" << si.ec << "\n"
	  << "# radial_points\t" << si.n_rows << "\n"
	  << "# profile_version\t" << si.prof_version << "\n"
	  << "#\n"
	  << "# --- thermal configuration ---\n"
	  << "# Tinf_initial_K\t" << kTinf0_K << "\n"
	  << "# t0_yr\t" << kT0_yr << "\n"
	  << "# t1_yr\t" << kT1_yr << "\n"
	  << "# envelope\tEnvelopeTbTs / Potekhin1997 Iron / xi=0 / rho_b=1e10\n"
	  << "# photon\tradiating_fraction=1, global_scale=1, C_star from StarContext\n"
	  << "# neutrino\tDU=on, MU=on, PBF=off, global_scale=1\n"
	  << "# state\tTHERMAL ONLY (no SpinState). Spin verified bit-identically decoupled.\n"
	  << "#\n"
	  << "# --- tolerances applied by the regression ---\n"
	  << "# tol_state\t" << kTolState << "\n"
	  << "# tol_luminosity\t" << kTolLumin << "\n"
	  << kCols << "\n";
	for (const auto &c : cks)
		o << c.t_yr << "\t" << c.Tinf_K << "\t" << c.C_star_erg_K << "\t" << c.L_nu_inf
		  << "\t" << c.L_nu_DU << "\t" << c.L_nu_MU << "\t" << c.L_gamma_inf << "\t"
		  << c.Tsurf_K << "\t" << c.Tb_K << "\t" << c.dLnTinf_dt << "\t"
		  << (c.L_gamma_inf > 0.0 ? c.L_nu_inf / c.L_gamma_inf : 0.0) << "\n";
}

static bool ReadBaseline(const fs::path &p, std::vector<Checkpoint> &out)
{
	std::ifstream in(p);
	if (!in)
		return false;
	std::string line;
	out.clear();
	while (std::getline(in, line))
	{
		if (line.empty() || line[0] == '#')
			continue;
		if (line.rfind("t_yr", 0) == 0)
			continue;
		std::istringstream is(line);
		Checkpoint c;
		double ratio = 0.0;
		is >> c.t_yr >> c.Tinf_K >> c.C_star_erg_K >> c.L_nu_inf >> c.L_nu_DU >>
			c.L_nu_MU >> c.L_gamma_inf >> c.Tsurf_K >> c.Tb_K >> c.dLnTinf_dt >> ratio;
		out.push_back(c);
	}
	return !out.empty();
}

static double Rel(double a, double b)
{
	if (b == 0.0)
		return (a == 0.0) ? 0.0 : 1.0;
	return std::fabs(a - b) / std::fabs(b);
}


int main(int argc, char **argv)
{
	std::cout << std::scientific << std::setprecision(8);
	if (argc < 2)
	{
		std::cerr << "usage: passive_cooling_regression <EOS_DATA_ROOT> [--compare|--emit <f>|--study]\n";
		return 2;
	}
	const fs::path root(argv[1]);
	std::string mode = (argc > 2) ? argv[2] : "--compare";
	fs::path bfile = (argc > 3) ? fs::path(argv[3]) : fs::path();

	auto StepName = [](P::Evolution::StepperType t) {
		switch (t)
		{
		case P::Evolution::StepperType::RKF45: return "RKF45";
		case P::Evolution::StepperType::RKCK: return "RKCK";
		case P::Evolution::StepperType::RK8PD: return "RK8PD";
		case P::Evolution::StepperType::RK2: return "RK2";
		case P::Evolution::StepperType::MSBDF: return "MSBDF";
		default: return "?";
		}
	};

	if (mode == "--spdscan")
	{
		StarInfo si;
		std::string err;
		for (double spd : {50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0})
		{
			RunOpts ro;
			ro.stepper = kBaselineStepper;
			ro.spd = spd;
			std::vector<Checkpoint> r;
			const bool ok = RunTrajectory(root, ro, r, si, err);
			std::cout << "  spd=" << spd << "  " << (ok ? "OK   Tinf(1e6yr)=" : "FAIL ") ;
			if (ok) std::cout << r.back().Tinf_K;
			else std::cout << err;
			std::cout << "\n" << std::flush;
			err.clear();
		}
		return 0;
	}

	if (mode == "--steppers")
	{
		// Candidate comparison + convergence, Configuration A (spin+thermal).
		const std::vector<P::Evolution::StepperType> cands = {
			P::Evolution::StepperType::RKF45,
			P::Evolution::StepperType::RKCK,
			P::Evolution::StepperType::RK8PD};
		const std::vector<double> rtols = {1e-6, 3e-7, 1e-7};
		std::map<std::string, std::vector<Checkpoint>> nominal;
		StarInfo si;
		std::string err;

		for (auto st : cands)
		{
			std::cout << "\n=== " << StepName(st) << " ===\n";
			std::vector<Checkpoint> ref;
			for (double rt : rtols)
			{
				RunOpts ro;
				ro.stepper = st;
				ro.rtol = rt;
				ro.atol = (rt <= 1e-7) ? 1e-12 : ((rt <= 3e-7) ? 1e-11 : 1e-10);
				ro.include_spin = true; // Configuration A: thermal-only not yet repaired
				std::vector<Checkpoint> r;
				const auto t_start = std::chrono::steady_clock::now();
				const bool ok = RunTrajectory(root, ro, r, si, err);
				const double secs = std::chrono::duration<double>(
										std::chrono::steady_clock::now() - t_start).count();
				if (!ok)
				{
					std::cout << "  rtol=" << rt << "  FAILED: " << err << "\n";
					err.clear();
					continue;
				}
				std::cout << "  rtol=" << rt << "  runtime=" << secs << " s"
						  << "  Tinf(1e6yr)=" << r.back().Tinf_K << "\n";
				if (rt == 1e-6)
					nominal[StepName(st)] = r;
				ref = r; // last (tightest) becomes the local reference
			}
			if (!ref.empty() && nominal.count(StepName(st)))
			{
				double worst = 0;
				for (std::size_t k = 0; k < ref.size(); ++k)
					worst = std::max(worst, Rel(nominal[StepName(st)][k].Tinf_K, ref[k].Tinf_K));
				std::cout << "  nominal(1e-6) vs tightest(1e-7): max rel dTinf = " << worst
						  << (worst < 1e-3 ? "   [within 1e-3 gate]" : "   [EXCEEDS 1e-3 GATE]")
						  << "\n";
			}
		}

		std::cout << "\n=== cross-stepper agreement at rtol=1e-6 ===\n";
		if (nominal.size() >= 2)
		{
			auto it = nominal.begin();
			const auto &base = it->second;
			const std::string bname = it->first;
			for (++it; it != nominal.end(); ++it)
			{
				double worst = 0;
				for (std::size_t k = 0; k < base.size(); ++k)
					worst = std::max(worst, Rel(it->second[k].Tinf_K, base[k].Tinf_K));
				std::cout << "  " << bname << " vs " << it->first
						  << ": max rel dTinf = " << worst << "\n";
			}
		}
		return 0;
	}

	if (mode == "--study")
	{
		RunOpts base;
		base.stepper = kBaselineStepper;
		StarInfo si;
		std::string err;
		std::vector<std::vector<Checkpoint>> reps;
		std::cout << "Numerics study for " << StepName(kBaselineStepper) << "\n";

		for (int i = 0; i < 3; ++i)
		{
			std::vector<Checkpoint> r;
			if (!RunTrajectory(root, base, r, si, err)) { std::cerr << err << "\n"; return 3; }
			reps.push_back(r);
		}
		double rep = 0;
		bool bitwise = true;
		for (std::size_t k = 0; k < reps[0].size(); ++k)
			for (int i = 1; i < 3; ++i)
			{
				rep = std::max(rep, Rel(reps[i][k].Tinf_K, reps[0][k].Tinf_K));
				if (reps[i][k].Tinf_K != reps[0][k].Tinf_K) bitwise = false;
			}
		std::cout << "\nA. repeatability (3 runs): max rel dTinf = " << rep
				  << (bitwise ? "  [BIT-IDENTICAL]" : "  [not bit-identical]") << "\n";

		auto compare = [&](const char *label, RunOpts ro) {
			std::vector<Checkpoint> r;
			if (!RunTrajectory(root, ro, r, si, err)) { std::cout << "  " << label << " FAILED: " << err << "\n"; err.clear(); return 0.0; }
			double wT = 0, wL = 0;
			for (std::size_t k = 0; k < r.size(); ++k)
			{
				wT = std::max(wT, Rel(reps[0][k].Tinf_K, r[k].Tinf_K));
				wL = std::max(wL, std::max(Rel(reps[0][k].L_nu_inf, r[k].L_nu_inf),
										   Rel(reps[0][k].L_gamma_inf, r[k].L_gamma_inf)));
			}
			std::cout << "  " << label << ": max rel dTinf = " << wT << ", dL = " << wL << "\n";
			return wT;
		};

		std::cout << "\nB. tolerance tightening\n";
		RunOpts t1o = base; t1o.rtol = 3e-7; t1o.atol = 1e-11;
		const double d_3e7 = compare("rtol 3e-7", t1o);
		RunOpts t2o = base; t2o.rtol = 1e-7; t2o.atol = 1e-12;
		const double d_1e7 = compare("rtol 1e-7", t2o);

		std::cout << "\nC. cadence sensitivity\n";
		RunOpts c1 = base; c1.spd = 75.0;  const double d_75  = compare("75 samples/decade", c1);
		RunOpts c2 = base; c2.spd = 300.0; const double d_300 = compare("300 samples/decade", c2);

		std::cout << "\nD. continuous vs segmented (restart at each checkpoint)\n";
		RunOpts seg = base; seg.segmented = true;
		const double d_seg = compare("segmented", seg);

		std::cout << "\nE. spin decoupling\n";
		RunOpts sp = base; sp.include_spin = true;
		const double d_spin = compare("spin+thermal vs thermal-only", sp);

		const double floor_fp = 1e-14;
		const double drive = std::max({rep, d_3e7, d_1e7, d_75, d_300, floor_fp});
		std::cout << "\nF. tolerance inputs\n"
				  << "   repeat=" << rep << "  tighten=" << std::max(d_3e7, d_1e7)
				  << "  cadence=" << std::max(d_75, d_300) << "  fp floor=" << floor_fp << "\n"
				  << "   max = " << drive << "   x5 safety = " << (5.0 * drive) << "\n"
				  << "   (segmented=" << d_seg << " and spin=" << d_spin
				  << " are reported findings, NOT tolerance inputs)\n";
		return 0;
	}

	if (mode == "--detector")
	{
		// Regression-detector proof: a MODEST 1% perturbation of the photon channel must
		// be caught by the committed tolerances. Not committed as a configuration.
		std::vector<Checkpoint> gold, pert;
		if (!ReadBaseline(bfile, gold)) { std::cerr << "cannot read baseline\n"; return 5; }
		RunOpts ro;
		ro.photon_global_scale = 1.01;
		StarInfo si;
		std::string err;
		if (!RunTrajectory(root, ro, pert, si, err)) { std::cerr << err << "\n"; return 3; }
		int caught = 0;
		for (std::size_t k = 0; k < gold.size(); ++k)
		{
			if (Rel(pert[k].Tinf_K, gold[k].Tinf_K) > kTolState) ++caught;
			if (Rel(pert[k].L_gamma_inf, gold[k].L_gamma_inf) > kTolLumin) ++caught;
		}
		std::cout << "detector: PhotonCooling global_scale 1.00 -> 1.01\n"
				  << "  checkpoints exceeding tolerance: " << caught << "\n"
				  << "  max rel dTinf = ";
		double w = 0;
		for (std::size_t k = 0; k < gold.size(); ++k)
			w = std::max(w, Rel(pert[k].Tinf_K, gold[k].Tinf_K));
		std::cout << w << "\n"
				  << (caught > 0 ? "  DETECTOR WORKS: a 1% photon perturbation fails the regression\n"
								 : "  DETECTOR FAILED: perturbation went unnoticed\n");
		return caught > 0 ? 0 : 1;
	}

	std::vector<Checkpoint> cks;
	StarInfo si;
	CacheSafety cs;
	std::string err;
	if (!RunTrajectory(root, {}, cks, si, err, &cs))
	{
		std::cerr << "FAILED: " << err << "\n";
		return 3;
	}

	std::cout << "star: requested " << si.requested_M << " Msun, achieved " << si.achieved_M
			  << " Msun, R = " << si.R_km << " km, " << si.n_rows << " points\n\n";
	std::cout << "t_yr\tTinf_K\tC_star\tL_nu\tL_gamma\tL_nu/L_gamma\tenergy_id_rel\n";
	for (const auto &c : cks)
		std::cout << c.t_yr << "\t" << c.Tinf_K << "\t" << c.C_star_erg_K << "\t"
				  << c.L_nu_inf << "\t" << c.L_gamma_inf << "\t"
				  << (c.L_gamma_inf > 0 ? c.L_nu_inf / c.L_gamma_inf : 0.0) << "\t"
				  << c.energy_identity_rel << "\n";

	int fails = 0;

	// -----------------------------------------------------------------------
	//  Phase-2B-3 (INV-12) — canonical-path cache-safety contract
	// -----------------------------------------------------------------------
	//  The known INV-12 hazards (stale GeometryCache, version-only cache identity,
	//  cross-star driver reuse, unrebound column pointers) all require either a
	//  structural mutation during evolution or a driver/context reused across
	//  stars. These assertions prove the canonical baseline does neither. They are
	//  observations of the DriverContext the observer already receives; they add
	//  no production API and cannot change any golden value.
	std::cout << "\nINV-12 canonical-path cache safety (" << cs.observations
			  << " driver-context observations during the run):\n";
	auto cache_chk = [&](const char *what, bool ok, const std::string &detail) {
		std::cout << (ok ? "  [ok]   " : "  [FAIL] ") << what << " — " << detail << "\n";
		if (!ok)
			++fails;
	};
	cache_chk("profile version is unchanged across the integration",
			  cs.version_at_start == cs.version_at_finish,
			  "start " + std::to_string(cs.version_at_start) + ", finish " +
				  std::to_string(cs.version_at_finish));
	cache_chk("profile version never changed mid-run",
			  cs.distinct_versions_during_run == 1,
			  std::to_string(cs.distinct_versions_during_run) +
				  " distinct version(s) observed");
	cache_chk("exactly one GeometryCache is used, and it is this run's",
			  cs.distinct_geo_ptrs == 1 && cs.geo_is_the_run_geo,
			  std::to_string(cs.distinct_geo_ptrs) +
				  " distinct GeometryCache object(s), identity matches the run's");
	cache_chk("exactly one StarContext is used, so no driver is reused across stars",
			  cs.distinct_star_ptrs == 1 && cs.star_is_the_run_star,
			  std::to_string(cs.distinct_star_ptrs) +
				  " distinct StarContext object(s), identity matches the run's");
	cache_chk("exactly one CompOSE_Thermo is used, so the heat-capacity key is stable",
			  cs.distinct_thermo_ptrs == 1,
			  std::to_string(cs.distinct_thermo_ptrs) + " distinct thermo object(s)");
	cache_chk("the run actually produced observations", cs.observations > 0,
			  std::to_string(cs.observations) + " observations");
	if (!fails)
		std::cout << "  => KNOWN CACHE HAZARDS ARE NOT REACHED BY THIS BASELINE.\n"
					 "     (Scope: this canonical procedure only. Not a claim about the\n"
					 "     general API — see docs/validation/CACHE_CORRECTNESS.md.)\n";

	for (const auto &c : cks)
		if (!(c.C_star_match_rel < kTolPatternA))
		{
			std::cerr << "\nADR-0002 PATTERN-A VIOLATION at t = " << c.t_yr
					  << " yr: photon and neutrino used different C_star (rel "
					  << c.C_star_match_rel << ")\n";
			++fails;
		}
	for (const auto &c : cks)
		if (!(c.energy_identity_rel < kTolEnergyIdentity))
		{
			std::cerr << "\nENERGY-ACCOUNTING DEFECT at t = " << c.t_yr
					  << " yr: identity residual " << c.energy_identity_rel << "\n";
			++fails;
		}
	if (fails)
		return 4; // not baseline drift — a physics/accounting defect

	unit_candidate_evidence::Capture("passive_cooling_cmf_1p6_debug.tsv", [&](const fs::path &p) { WriteBaseline(p, cks, si); });

	if (mode == "--emit")
	{
		WriteBaseline(bfile, cks, si);
		std::cout << "\nbaseline written: " << bfile << "\n";
		return 0;
	}

	std::vector<Checkpoint> gold;
	if (!ReadBaseline(bfile, gold))
	{
		std::cerr << "cannot read baseline: " << bfile << "\n";
		return 5;
	}
	if (gold.size() != cks.size())
	{
		std::cerr << "checkpoint count mismatch: baseline " << gold.size() << " vs run "
				  << cks.size() << "\n";
		return 6;
	}

	std::cout << "\ncomparison against golden baseline (state tol " << kTolState
			  << ", luminosity tol " << kTolLumin << "):\n";
	auto chk = [&](const char *name, double a, double b, double tol, double t) {
		const double r = Rel(a, b);
		if (!(r <= tol))
		{
			std::cerr << "  REGRESSION t=" << t << " yr  " << name << ": run=" << a
					  << " baseline=" << b << " rel=" << r << " > " << tol << "\n";
			++fails;
		}
	};
	for (std::size_t k = 0; k < cks.size(); ++k)
	{
		const auto &c = cks[k];
		const auto &g = gold[k];
		chk("t_yr", c.t_yr, g.t_yr, 1e-12, g.t_yr);
		chk("Tinf_K", c.Tinf_K, g.Tinf_K, kTolState, g.t_yr);
		chk("C_star_erg_K", c.C_star_erg_K, g.C_star_erg_K, kTolState, g.t_yr);
		chk("dLnTinf_dt_1_s", c.dLnTinf_dt, g.dLnTinf_dt, kTolLumin, g.t_yr);
		chk("Tsurf_K", c.Tsurf_K, g.Tsurf_K, kTolState, g.t_yr);
		chk("Tb_K", c.Tb_K, g.Tb_K, kTolState, g.t_yr);
		chk("L_nu_inf_erg_s", c.L_nu_inf, g.L_nu_inf, kTolLumin, g.t_yr);
		chk("L_nu_DU_inf_erg_s", c.L_nu_DU, g.L_nu_DU, kTolLumin, g.t_yr);
		chk("L_nu_MU_inf_erg_s", c.L_nu_MU, g.L_nu_MU, kTolLumin, g.t_yr);
		chk("L_gamma_inf_erg_s", c.L_gamma_inf, g.L_gamma_inf, kTolLumin, g.t_yr);
	}
	std::cout << (fails == 0 ? "  all checkpoints within tolerance\n"
							 : "  " + std::to_string(fails) + " regression(s)\n");
	return fails == 0 ? 0 : 1;
}
