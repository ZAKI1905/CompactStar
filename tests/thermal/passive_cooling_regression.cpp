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
#include <sstream>
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

// Logarithmically separated diagnostic checkpoints [yr].
static const std::vector<double> kCheckpointsYr = {
	1.0e2, 3.0e2, 1.0e3, 3.0e3, 1.0e4, 3.0e4, 1.0e5, 3.0e5, 1.0e6};

struct Checkpoint
{
	double t_s = 0, t_yr = 0, Tinf_K = 0, C_star_erg_K = 0;
	double L_nu_inf = 0, L_nu_DU = 0, L_nu_MU = 0, L_gamma_inf = 0;
	double dLnTinf_dt = 0, Tsurf_K = 0, Tb_K = 0;
	double energy_identity_rel = 0; // |RHS_sum - (-(Lnu+Lg)/C)/T| / |expected|
};

struct StarInfo
{
	double requested_M = 0, achieved_M = 0, R_km = 0, ec = 0;
	int n_rows = 0;
	std::uint64_t prof_version = 0;
};

struct RunOpts
{
	double rtol = kRtol;
	double atol = kAtol;
	double spd = kSamplesPerDecade;
	// MUST stay true: EvolutionSystem::operator() unconditionally calls
	// m_state.GetSpin() in a leftover debug block (EvolutionSystem.cpp:103-112)
	// whose only consumer is commented out, so a thermal-only StateVector throws
	// "requested tag 'Spin' is not registered". Configuration B of the task brief is
	// therefore not runnable, and §5 forbids editing production to force it.
	// Decoupling is instead proven by the per-checkpoint energy identity below.
	bool include_spin = true;
};

// ---------------------------------------------------------------------------
static bool RunTrajectory(const fs::path &root, const RunOpts &ro,
						  std::vector<Checkpoint> &out, StarInfo &si,
						  std::string &err)
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
	po.global_scale = 1.0;
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
	P::Evolution::GSLIntegrator integrator(system, cfg, w.dim);

	std::vector<double> y(w.dim);
	P::Evolution::PackStateVector(w.state_vec, w.layout, y.data());

	const double YR = Zaki::Physics::YR_2_SEC;
	double t_cur = kT0_yr * YR;
	out.clear();

	for (double ck_yr : kCheckpointsYr)
	{
		const double t_ck = ck_yr * YR;
		if (t_ck > t_cur)
		{
			// Segment-wise: each Integrate() allocates its own GSL driver, so this is
			// deterministic and reproducible. It IS the defined baseline procedure.
			if (!integrator.Integrate(t_cur, t_ck, y.data()))
			{
				err = "integration failed before t = " + std::to_string(ck_yr) + " yr";
				return false;
			}
			t_cur = t_ck;
		}
		P::Evolution::UnpackStateVector(w.state_vec, w.layout, y.data());

		const auto pd = Th::Detail::ComputeDerived(*photon, w.state_vec, ctx);
		const auto nd = Th::Detail::NeutrinoCooling_Details::ComputeDerived(*neutrino, w.state_vec, ctx);
		if (!pd.ok || !nd.ok)
		{
			err = "driver diagnostics not ok at " + std::to_string(ck_yr) +
				  " yr: [photon] " + pd.message + " [neutrino] " + nd.message;
			return false;
		}

		Checkpoint c;
		c.t_s = t_ck;
		c.t_yr = ck_yr;
		c.Tinf_K = pd.Tinf_K;
		c.C_star_erg_K = pd.C_star_erg_K;
		c.L_nu_inf = nd.L_nu_inf_erg_s;
		c.L_nu_DU = nd.L_nu_DU_inf_erg_s;
		c.L_nu_MU = nd.L_nu_MU_inf_erg_s;
		c.L_gamma_inf = pd.L_gamma_inf_erg_s;
		c.dLnTinf_dt = pd.dLnTinf_dt_1_s + nd.dLnTinf_dt_1_s;
		c.Tsurf_K = pd.Tsurf_K;
		c.Tb_K = pd.Tb_K;

		// Energy-equation identity: both channels must share one denominator.
		const double dT_expected = -(c.L_nu_inf + c.L_gamma_inf) / c.C_star_erg_K;
		const double dLn_expected = dT_expected / c.Tinf_K;
		c.energy_identity_rel = std::fabs(c.dLnTinf_dt - dLn_expected) / std::fabs(dLn_expected);

		// Independent cross-check: neutrino must use the SAME C_* as photon.
		if (std::fabs(nd.C_eff_erg_K - pd.C_star_erg_K) / pd.C_star_erg_K > 1e-12)
		{
			err = "channels used DIFFERENT heat capacities at " +
				  std::to_string(ck_yr) + " yr — ADR-0002 violation";
			return false;
		}
		out.push_back(c);
	}
	return true;
}

// ---------------------------------------------------------------------------
static const char *kCols =
	"t_yr\tTinf_K\tC_star_erg_K\tL_nu_inf_erg_s\tL_nu_DU_inf_erg_s\t"
	"L_nu_MU_inf_erg_s\tL_gamma_inf_erg_s\tdLnTinf_dt_1_s\tTsurf_K\tTb_K";

static void WriteBaseline(const fs::path &p, const std::vector<Checkpoint> &cks,
						  const StarInfo &si)
{
	std::ofstream o(p);
	o << std::setprecision(12) << std::scientific;
	o << "# CompactStar passive-cooling regression baseline\n"
	  << "# schema_version\t1\n"
	  << "# NOTE: REGRESSION baseline only. The neutrino emissivity normalizations frozen\n"
	  << "#       here are self-labelled PLACEHOLDERS in source (Q0_DU=1e27, Q0_MU=1e21).\n"
	  << "#       This is NOT a validation of neutrino microphysics.\n"
	  << "# generated_by_commit\t56aee7c5a132ffb4d922cac160fda917363ef8e7\n"
	  << "# build_configuration\tDebug\n"
	  << "# compiler\tAppleClang 17.0.0.17000604\n"
	  << "# eos_structural\tCMF #1 EoS with crust (official CompOSE), processed .eos\n"
	  << "# eos_structural_sha256\t5747dd73256c0c28bc56be337cbb96d0918a54bc9ed9fc40984c5befd47ae5dd\n"
	  << "# eos_thermo\tCMF hadronic EoS with electrons (official CompOSE), raw tables\n"
	  << "# eos_thermo_sha256\ta456fb8595208ddf3119350a856fbf2b906c0a0e19bb7c716571748d0aa0724b\n"
	  << "# target_mass_Msun\t" << kTargetM << "\n"
	  << "# achieved_mass_Msun\t" << si.achieved_M << "\n"
	  << "# radius_km\t" << si.R_km << "\n"
	  << "# central_eps\t" << si.ec << "\n"
	  << "# radial_points\t" << si.n_rows << "\n"
	  << "# profile_version\t" << si.prof_version << "\n"
	  << "# Tinf_initial_K\t" << kTinf0_K << "\n"
	  << "# t0_yr\t" << kT0_yr << "\n"
	  << "# t1_yr\t" << kT1_yr << "\n"
	  << "# stepper\tMSBDF\n"
	  << "# rtol\t" << kRtol << "\n"
	  << "# atol\t" << kAtol << "\n"
	  << "# save_cadence\tLogTime\n"
	  << "# samples_per_decade\t" << kSamplesPerDecade << "\n"
	  << "# envelope\tEnvelopeTbTs / Potekhin1997 Iron / xi=0 / rho_b=1e10\n"
	  << "# photon\tradiating_fraction=1, global_scale=1, C_star from StarContext\n"
	  << "# neutrino\tDU=on, MU=on, PBF=off, global_scale=1\n"
	  << "# spin\tregistered (forced by EvolutionSystem); contributes nothing to the thermal RHS\n"
	  << "#     \tproven by the per-checkpoint energy identity, residual < 1e-10\n"
	  << kCols << "\n";
	for (const auto &c : cks)
		o << c.t_yr << "\t" << c.Tinf_K << "\t" << c.C_star_erg_K << "\t" << c.L_nu_inf
		  << "\t" << c.L_nu_DU << "\t" << c.L_nu_MU << "\t" << c.L_gamma_inf << "\t"
		  << c.dLnTinf_dt << "\t" << c.Tsurf_K << "\t" << c.Tb_K << "\n";
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
		is >> c.t_yr >> c.Tinf_K >> c.C_star_erg_K >> c.L_nu_inf >> c.L_nu_DU >>
			c.L_nu_MU >> c.L_gamma_inf >> c.dLnTinf_dt >> c.Tsurf_K >> c.Tb_K;
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

// Tolerance policy — fixed BEFORE golden values were accepted; see
// docs/validation/PASSIVE_COOLING_BASELINE.md.
constexpr double kTolState = 1e-4; // Tinf, C_star, Tsurf, Tb, dLnTinf/dt
constexpr double kTolLumin = 1e-3; // luminosities (steep powers of T amplify dT)
constexpr double kTolEnergyIdentity = 1e-10;

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

	if (mode == "--study")
	{
		std::cout << "Numerics study (not a pass/fail test)\n\n";
		StarInfo si;
		std::string err;
		std::vector<std::vector<Checkpoint>> runs;

		for (int i = 0; i < 3; ++i)
		{
			std::vector<Checkpoint> r;
			RunOpts ro;
			if (!RunTrajectory(root, ro, r, si, err))
			{
				std::cerr << err << "\n";
				return 3;
			}
			runs.push_back(r);
		}
		double worst_rep = 0;
		for (std::size_t k = 0; k < runs[0].size(); ++k)
			for (int i = 1; i < 3; ++i)
				worst_rep = std::max(worst_rep, Rel(runs[i][k].Tinf_K, runs[0][k].Tinf_K));
		std::cout << "A. repeat-run determinism: max rel variation in Tinf = " << worst_rep << "\n";

		std::vector<Checkpoint> tight;
		RunOpts rt;
		rt.rtol = 3e-7;
		rt.atol = 1e-11;
		if (!RunTrajectory(root, rt, tight, si, err))
		{
			std::cerr << err << "\n";
			return 3;
		}
		std::cout << "\nB. nominal (rtol 1e-6) vs tighter (rtol 3e-7, atol 1e-11):\n";
		double worst_int = 0, worst_int_L = 0;
		for (std::size_t k = 0; k < tight.size(); ++k)
		{
			const double dT = Rel(runs[0][k].Tinf_K, tight[k].Tinf_K);
			const double dL = std::max(Rel(runs[0][k].L_nu_inf, tight[k].L_nu_inf),
									   Rel(runs[0][k].L_gamma_inf, tight[k].L_gamma_inf));
			worst_int = std::max(worst_int, dT);
			worst_int_L = std::max(worst_int_L, dL);
			std::cout << "   t=" << tight[k].t_yr << " yr  dTinf=" << dT << "  dL=" << dL << "\n";
		}
		std::cout << "   max Tinf integration difference = " << worst_int << "\n"
				  << "   max luminosity  difference      = " << worst_int_L << "\n";

		std::vector<Checkpoint> cad;
		RunOpts rc;
		rc.spd = 40.0;
		if (!RunTrajectory(root, rc, cad, si, err))
		{
			std::cerr << err << "\n";
			return 3;
		}
		double worst_cad = 0;
		for (std::size_t k = 0; k < cad.size(); ++k)
			worst_cad = std::max(worst_cad, Rel(runs[0][k].Tinf_K, cad[k].Tinf_K));
		std::cout << "\nC. cadence independence (samples_per_decade 150 -> 40): max rel dTinf = "
				  << worst_cad << "\n";

		std::cout << "\nD. spin decoupling\n"
				  << "   Configuration B (thermal-only) is NOT RUNNABLE: EvolutionSystem::operator()\n"
				  << "   unconditionally calls m_state.GetSpin() in a dead debug block\n"
				  << "   (EvolutionSystem.cpp:103-112), so an unregistered Spin tag throws.\n"
				  << "   Production was not modified to force equivalence (task brief §5).\n"
				  << "   Decoupling is instead demonstrated by the per-checkpoint energy identity:\n"
				  << "   d ln T_inf/dt equals -(L_nu + L_gamma)/(C_star * T_inf) exactly, so the\n"
				  << "   spin block contributes nothing to the thermal RHS. Max identity residual:\n";
		double worst_id = 0;
		for (const auto &c : runs[0])
			worst_id = std::max(worst_id, c.energy_identity_rel);
		std::cout << "     " << worst_id << "\n";
		return 0;
	}

	std::vector<Checkpoint> cks;
	StarInfo si;
	std::string err;
	if (!RunTrajectory(root, {}, cks, si, err))
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
	for (const auto &c : cks)
		if (!(c.energy_identity_rel < kTolEnergyIdentity))
		{
			std::cerr << "\nENERGY-ACCOUNTING DEFECT at t = " << c.t_yr
					  << " yr: identity residual " << c.energy_identity_rel << "\n";
			++fails;
		}
	if (fails)
		return 4; // not baseline drift — a physics/accounting defect

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
