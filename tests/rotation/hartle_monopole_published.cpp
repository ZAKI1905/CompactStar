// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_published.cpp
 * @brief Phase 4D Experiment K — the governed O(Omega^2) monopole response against PUBLISHED
 *        second-order results. Self-contained (historical fixture embedded).
 *
 * K1. Chandrasekhar & Miller (1974), MNRAS 167, 63, Table I (p. 73): slowly rotating
 *     HOMOGENEOUS configurations, R/R_S = 1.15 ... 100, tabulated to 4-5 significant figures.
 *     Compared on CompactStar's exact constant-density background (Phase 4D analytic fixture):
 *       I/MR^2                        first order
 *       varpi_1 = (Omega - omega)_1   first order:   s(R) = varpi_1 * I_tab       (unit GJ/R_S^3 c^2)
 *       xi_0                          SECOND order:  xi_hat(R)/R_S^3 = xi0_tab * I_tab^2  (unit G^2J^2/R_S^3 c^6)
 *       deltaM/M                      SECOND order, in C&M's OWN convention — see below.
 *     C&M impose "all of these functions join continuously with their exterior solutions"
 *     (their eq. 32), i.e. m0 CONTINUOUS at the boundary. For a density discontinuity that
 *     omits the surface-shell contribution 4 pi R^2 eps xi_0(R) that Hartle's (97) carries as a
 *     delta function (Phase 4C-G §9.3; ADR-0007 P4). Their tabulated deltaM/M is therefore
 *     [m0(R^-) + J^2/R^3]/M, whose Newtonian limit is I Omega^2/M — NOT the physical
 *     fixed-central-density mass change Omega^2 R^3/M that Experiment D verifies. The
 *     comparison is made like-for-like (shell excluded) and the shell-inclusive value is
 *     printed beside it so the reader sees the size of the omitted term.
 *
 * K2. Hartle & Thorne (1968), ApJ 153, 807, Tables 3 and 5 (Harrison-Wheeler EOS, Table 1):
 *     eight configurations, E_c = 1e14 ... 1e17 g/cm^3, at Omega^2 = M/R^3:
 *       Table 3: R, M                          (non-rotating, checks the EOS reconstruction)
 *       Table 5: R_g/R = sqrt(I/M)/R, omega_s/Omega = 1 - s(R), omega_c/Omega = 1 - s(0),
 *                deltaR/R = xi_hat(R) M/R^4 (SECOND order), deltaM/M = deltaM_hat/R^3 (SECOND order)
 *       and Omega [s^-1] = c sqrt(M/R^3) as a unit sanity check.
 *     HT68 quote "about 1 per cent or better" for their models; the EOS is printed to three
 *     significant figures. ADR-0007 §7 item 8 predeclares 2e-2.
 *
 * Neither comparison uses any CompactStar number as an expected answer.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "eos_table_knots.hpp"
#include "hartle_monopole_reference.hpp"
#include "hartle_profile_compare.hpp"
#include "hartle_thorne_1968_hw_eos.hpp"

#include <Zaki/Physics/Constants.hpp>

#include <unistd.h>

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolveStatus;
using CompactStar::Core::TOVSolver;
using hartle_4b::UniformStar;

/// Exposes the protected cutoff accessor for the ADR-0009 completeness check.
struct Probe : TOVSolver
{
	using TOVSolver::PressureCutoff;
};

static constexpr double kK_Published = 2.0e-2; // ADR-0007 §7 item 8
static constexpr double kF_Smooth = 2.0e-5;	// ADR-0008 Validation F (smooth-EOS measure-vs-differential equivalence on deltaM_hat)
// AUTHENTICATED SCOPE of the 2e-5 bound (Phase 4D-RV, recorded chronology): ADR-0008's Context
// table set "smooth HW EOS vs today <= 1.25e-5" from scratch stars at eps_c = 3e14, 1e15, 3e15
// g/cm^3 (PHASE4D_R_EOS_MEASURE_DERIVATION.md §8), and 4D-RI §13 asserted the bound on deltaM_hat at
// exactly those three densities (5.5e-6 / 1.25e-5 / 1.15e-5). A first 4D-RV run of this file
// asserted the same bound over all eight HT68 Table-5 configurations and measured 4.1e-5 (6e15),
// 1.7e-4 (1e16), 2.5e-5 (3e16), 1.4e-4 (1e17): the SUPERSEDED nodal-derivative form loses more of
// the densified table's sub-node slope structure on the denser, steeper stars — the ADR-0008
// mechanism itself — while the INDEPENDENT Stieltjes oracles agree with production to <= 2e-6 on
// all eight (F2). The bound is therefore asserted at its authenticated scope (F1) and every
// configuration is reported (F1r); nothing was widened.
static const double kF_AuthenticatedEc[] = {3.0e14, 1.0e15, 3.0e15};
static constexpr std::size_t kN_Grid = 4001;

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
static std::string Sci(double v, int p = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", p, v);
	return b;
}
static double Rel(double a, double b)
{
	return std::fabs(b) > 0.0 ? std::fabs(a - b) / std::fabs(b) : std::fabs(a - b);
}

// ---------------------------------------------------------------------------
//  K1 — Chandrasekhar & Miller 1974, Table I (transcribed from the 400-dpi scan, p. 73)
// ---------------------------------------------------------------------------
struct CMRow
{
	double x, I, I_MR2, dM_M, Q, w1, xi0, xi2, eps;
};
static const std::vector<CMRow> kCM = {
	// R/R_S   I          I/MR^2     dM/M       Q          varpi_1     xi_0        xi_2        eps
	{1.125, 5.057e-1, 7.991e-1, 3.442e0, 2.002e0, 5.7270e-1, -4.7390e-1, -9.2667e-1, 5.2307e0},
	{1.15, 5.105e-1, 7.720e-1, 3.454e0, 2.076e0, 6.4391e-1, -4.2403e-1, -1.3942e0, 6.0898e0},
	{1.2, 5.248e-1, 7.289e-1, 3.412e0, 2.169e0, 7.4818e-1, -2.7505e-1, -2.0869e0, 6.7278e0},
	{1.3, 5.657e-1, 6.695e-1, 3.225e0, 2.441e0, 8.5741e-1, +1.2203e-1, -3.5565e0, 7.9708e0},
	{1.4, 6.171e-1, 6.297e-1, 2.993e0, 2.813e0, 8.9174e-1, +5.8490e-1, -5.0402e0, 9.0327e0},
	{1.5, 6.758e-1, 6.007e-1, 2.757e0, 3.272e0, 8.8714e-1, +1.0740e0, -6.4841e0, 9.8930e0},
	{1.6, 7.406e-1, 5.786e-1, 2.533e0, 3.806e0, 8.6205e-1, +1.5681e0, -7.8651e0, 1.0569e1},
	{1.7, 8.106e-1, 5.610e-1, 2.327e0, 4.406e0, 8.2652e-1, +2.0553e0, -9.1734e0, 1.1088e1},
	{1.8, 8.856e-1, 5.467e-1, 2.140e0, 5.064e0, 7.8622e-1, +2.5287e0, -1.0407e1, 1.1478e1},
	{1.9, 9.653e-1, 5.348e-1, 1.971e0, 5.773e0, 7.4440e-1, +2.9847e0, -1.1566e1, 1.1761e1},
	{2.0, 1.049e0, 5.245e-1, 1.819e0, 6.525e0, 7.0294e-1, +3.4213e0, -1.2655e1, 1.1960e1},
	{2.5, 1.534e0, 4.909e-1, 1.259e0, 1.081e1, 5.2376e-1, +5.3097e0, -1.7183e1, 1.2131e1},
	{3.0, 2.123e0, 4.718e-1, 9.163e-1, 1.567e1, 3.9700e-1, +6.7734e0, -2.0549e1, 1.1658e1},
	{4.0, 3.604e0, 4.505e-1, 5.440e-1, 2.636e1, 2.4623e-1, +8.8387e0, -2.5151e1, 1.0300e1},
	{5.0, 5.487e0, 4.390e-1, 3.588e-1, 3.773e1, 1.6625e-1, +1.0206e1, -2.8126e1, 9.0297e0},
	{10.0, 2.091e1, 4.182e-1, 9.493e-2, 9.785e1, 4.5821e-2, +1.3235e1, -3.4564e1, 5.3517e0},
	{20.0, 8.177e1, 4.089e-1, 2.437e-2, 2.217e2, 1.1980e-2, +1.4900e1, -3.8031e1, 2.8966e0},
	{35.0, 2.481e2, 4.051e-1, 8.046e-3, 4.086e2, 3.9848e-3, +1.5644e1, -3.9569e1, 1.7106e0},
	{50.0, 5.043e2, 4.034e-1, 3.960e-3, 5.960e2, 1.9668e-3, +1.5948e1, -4.0192e1, 1.2131e0},
	{100.0, 2.009e3, 4.018e-1, 9.950e-4, 1.221e3, 4.9585e-4, +1.6305e1, -4.0926e1, 6.1574e-1},
};

/// Exact constant-density interior on a uniform grid with the exact EOS derivative (= 0).
static std::unique_ptr<NStar> HomogeneousStar(const UniformStar &u, std::size_t N)
{
	const double km2_to_gcm3 = Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;
	const double km2_to_dyn = Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_INV_KM2;
	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = hartle_4b::UniformGridR(u, i, N);
		pts.emplace_back(r, u.m(r) / Zaki::Physics::SUN_M_KM, u.nuprime(r) / 1.0e5, 0.0, u.p(r) * km2_to_dyn,
						 u.rho0 * km2_to_gcm3, 0.1, std::vector<double>{}, 0.0);
	}
	auto ns = std::make_unique<NStar>(pts);
	if (!ns->ComputeHartleMonopoleResponse())
		return nullptr;
	return ns;
}

static void RunChandrasekharMiller()
{
	std::cout << "K1. Chandrasekhar & Miller 1974, Table I — homogeneous configurations (R = 12 km, M = R/(2 R/R_S))\n"
				 "    Bound " << Sci(kK_Published, 0) << " on every compared quantity. R/R_S = 1.125 is the Buchdahl limit (p_c infinite) and is skipped.\n\n"
				 "    R/R_S   I/MR^2 prod  tab      rel     | s(R) prod    tab       rel     | xi0 prod     tab        rel     | dM/M(C&M) prod  tab      rel     | shell-incl.  ratio\n";
	double worst_I = 0, worst_s = 0, worst_xi = 0, worst_dM = 0;
	int n_ok = 0, n_tot = 0;
	for (const auto &row : kCM)
	{
		if (row.x <= 1.13)
			continue;
		const double R = 12.0, M = R / (2.0 * row.x), RS = 2.0 * M;
		const UniformStar u = hartle_4b::MakeUniform(M, R);
		auto ns = HomogeneousStar(u, kN_Grid);
		if (!ns)
		{
			Report("K1 build R/R_S=" + Sci(row.x, 3), false, "no response");
			continue;
		}
		const NStar &c = *ns;
		const auto *mono = c.MonopoleResponse();
		const auto &rot = c.RotationResponse();
		const int last = static_cast<int>(rot.omega_bar_over_Omega.Size()) - 1;
		const double I = mono->I, Rstar = mono->R_surface;
		const double I_MR2 = I / (M * R * R);
		const double sR = rot.omega_bar_over_Omega[last];
		const double sR_tab = row.w1 * row.I; // varpi_1 [J/R_S^3] * I [R_S^3] = omegabar(R)/Omega
		const double xi = mono->surface_xi0_over_Omega2 / (RS * RS * RS);
		const double xi_tab = row.xi0 * row.I * row.I; // xi0 [J^2/R_S^3] * I^2 -> xi_hat / R_S^3
		const double conv = RS * RS * RS * RS / (M * I * I); // deltaM_hat -> deltaM/M in units J^2/R_S^4
		const double dM_noshell = (mono->m0_over_Omega2[last] + I * I / (Rstar * Rstar * Rstar)) * conv;
		const double dM_full = mono->deltaM_over_Omega2 * conv;
		const double rI = Rel(I_MR2, row.I_MR2), rs = Rel(sR, sR_tab), rx = Rel(xi, xi_tab), rd = Rel(dM_noshell, row.dM_M);
		worst_I = std::max(worst_I, rI);
		worst_s = std::max(worst_s, rs);
		worst_xi = std::max(worst_xi, rx);
		worst_dM = std::max(worst_dM, rd);
		const bool ok = rI <= kK_Published && rs <= kK_Published && rx <= kK_Published && rd <= kK_Published;
		n_ok += ok;
		++n_tot;
		char b[512];
		snprintf(b, sizeof(b), "    %6.3f  %.5f  %.4f  %.1e | %.6f  %.6f  %.1e | %+.5e %+.5e %.1e | %.5e  %.4e  %.1e | %.4e  %.3e %s\n",
				 row.x, I_MR2, row.I_MR2, rI, sR, sR_tab, rs, xi, xi_tab, rx, dM_noshell, row.dM_M, rd, dM_full,
				 dM_full / dM_noshell, ok ? "" : "  <-- outside bound");
		std::cout << b;
	}
	Report("K1a I/MR^2 agrees with C&M Table I at every tabulated compactness", worst_I <= kK_Published, "worst rel=" + Sci(worst_I));
	Report("K1b s(R) = varpi_1 I agrees with C&M Table I (first order, 5 s.f.)", worst_s <= kK_Published, "worst rel=" + Sci(worst_s));
	Report("K1c xi_0(R) agrees with C&M Table I (SECOND order, incl. the sign change near R/R_S = 1.25)", worst_xi <= kK_Published,
		   "worst rel=" + Sci(worst_xi));
	Report("K1d deltaM/M agrees with C&M Table I in C&M's convention (m0 continuous at the boundary, shell excluded)",
		   worst_dM <= kK_Published, "worst rel=" + Sci(worst_dM));
	std::cout << "    " << n_ok << "/" << n_tot << " configurations inside the bound on all four quantities.\n"
			  << "    The shell-inclusive column is the PHYSICAL fixed-central-density deltaM (Experiment D: -> Omega^2 R^3/M);\n"
			  << "    C&M's boundary condition (their eq. 32) drops the shell for a density discontinuity, which is why the\n"
			  << "    ratio grows like ~ (5/2)(R/R_S) at large R/R_S. The agreement of the shell-EXCLUDED value with their\n"
			  << "    table validates the interior m0 integration; the shell itself was validated analytically in Experiments C, D and J.\n";
}

// ---------------------------------------------------------------------------
//  K2 — Hartle & Thorne 1968, Tables 3 and 5 on the Harrison-Wheeler EOS
// ---------------------------------------------------------------------------
static void RunHartleThorne(const fs::path &wrk)
{
	std::cout << "\nK2. Hartle & Thorne 1968, Tables 3 and 5 — Harrison-Wheeler EOS (Table 1), Omega^2 = M/R^3\n";
	const fs::path eos = wrk / "hw1968_densified.eos";
	if (!ht68::WriteDensifiedEos(eos.string(), 40))
	{
		Report("K2 fixture written", false, eos.string());
		return;
	}
	std::cout << "    fixture: " << ht68::HWTable1().size() << " printed rows, densified log-log to 40/decade -> " << eos.filename().string()
			  << "; G c^-2 = " << Sci(ht68::kGoverC2_cm_per_g, 3) << " cm/g (HT68 footnote)\n\n";

	struct Out
	{
		double ec, R, M, Rg_R, Omega, ws, wc, dR_R, dM_M, shell_over_dM;
		std::size_t n;
		bool complete = false;
		double f_diff = 0, f_stj = 0, f_knot = 0, f_sec = 0; // deltaM_hat rel: production vs the oracle forms
	};
	const eos_knots::Knots knots = eos_knots::Read(eos.string());
	std::vector<Out> outs;
	const auto &T3 = ht68::Table3();
	const auto &T5 = ht68::Table5();
	const double c_km_s = ht68::kC_cm_per_s * 1.0e-5;

	for (std::size_t k = 0; k < T3.size(); ++k)
	{
		Probe tov;
		const fs::path w = wrk / ("hw" + std::to_string(k));
		fs::create_directories(w);
		tov.SetWrkDir(w.string());
		tov.ImportEOS(eos.string(), true);
		std::vector<TOVPoint> pts;
		if (tov.SingleStarSolveToTOVPoints(T3[k].ec, pts) <= 0 || pts.size() < 4)
		{
			Report("K2 TOV at E_c=" + Sci(T3[k].ec, 2), false, "no profile");
			continue;
		}
		const bool complete = tov.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED && pts.back().p == tov.PressureCutoff() &&
							  std::isfinite(pts.back().e) && std::isfinite(pts.back().dedp);
		NStar ns(pts);
		if (!ns.ComputeHartleMonopoleResponse())
		{
			Report("K2 monopole at E_c=" + Sci(T3[k].ec, 2), false, "no response");
			continue;
		}
		const NStar &c = ns;
		const auto *mono = c.MonopoleResponse();
		const auto &rot = c.RotationResponse();
		const int last = static_cast<int>(rot.omega_bar_over_Omega.Size()) - 1;
		Out o;
		o.ec = T3[k].ec;
		o.n = pts.size();
		o.R = mono->R_surface;
		o.M = pts.back().m; // Msun
		const double Mkm = o.M * Zaki::Physics::SUN_M_KM;
		o.Rg_R = std::sqrt(mono->I / Mkm) / o.R;
		o.Omega = std::sqrt(Mkm / (o.R * o.R * o.R)) * c_km_s;
		o.ws = 1.0 - rot.omega_bar_over_Omega[last];
		o.wc = 1.0 - rot.omega_bar_over_Omega[0];
		o.dR_R = mono->surface_xi0_over_Omega2 * Mkm / (o.R * o.R * o.R * o.R);
		o.dM_M = mono->deltaM_over_Omega2 / (o.R * o.R * o.R);
		o.shell_over_dM = mono->surface_shell_mass_over_Omega2 / mono->deltaM_over_Omega2;
		o.complete = complete;
		// F — ADR-0008 Validation F on the smooth HW EOS: production (measure) vs the SUPERSEDED
		// differential oracle on deltaM_hat, plus the INDEPENDENT Stieltjes oracles (reported).
		{
			const auto &P = c.Profile();
			const auto *r = P.GetRadius();
			const int n = static_cast<int>(r->Size());
			hartle_mono_ref::Background2 bg;
			for (int i = 0; i < n; ++i)
			{
				bg.r.push_back((*r)[i]);
				bg.p.push_back((*P.GetPressure())[i]);
				bg.eps.push_back((*P.GetEnergyDensity())[i]);
				bg.m.push_back((*P.GetMass())[i]);
				bg.nu.push_back((*P.GetMetricNu())[i]);
				bg.nup.push_back((*P.GetMetricNuPrime())[i]);
				bg.dedp.push_back((*P.GetEosDEdP())[i]);
				bg.s.push_back(rot.omega_bar_over_Omega[i]);
				bg.sp.push_back(rot.domega_bar_over_Omega_dr[i]);
			}
			hartle_mono_ref::MHOptions od;
			od.I_exterior = mono->I;
			od.eos_measure = false;
			const auto diff = hartle_mono_ref::Solve(bg, od);
			hartle_mono_ref::MHOptions os = od;
			os.eos_measure = true;
			const auto sec = hartle_mono_ref::Solve(bg, os);
			hartle_mono_ref::StieltjesOptions sp;
			sp.refine = 4;
			const auto stj = hartle_mono_ref::SolveStieltjes(bg, od, sp);
			hartle_mono_ref::StieltjesOptions sk;
			sk.refine = 2;
			sk.knot_p = &knots.p;
			sk.knot_eps = &knots.eps;
			const auto knt = hartle_mono_ref::SolveStieltjes(bg, od, sk);
			o.f_diff = diff.ok ? Rel(mono->deltaM_over_Omega2, diff.deltaM_hat) : 1.0;
			o.f_sec = sec.ok ? Rel(mono->deltaM_over_Omega2, sec.deltaM_hat) : 1.0;
			o.f_stj = stj.ok ? Rel(mono->deltaM_over_Omega2, stj.deltaM_hat) : 1.0;
			o.f_knot = (knots.ok && knt.ok) ? Rel(mono->deltaM_over_Omega2, knt.deltaM_hat) : 1.0;
		}
		outs.push_back(o);
	}

	std::cout << "    E_c[g/cm^3]  nodes | R[km] prod  tab    rel    | M prod   tab    rel    | Rg/R prod tab   rel    | Omega prod  tab     rel\n";
	double wR = 0, wM = 0, wRg = 0, wOm = 0, wws = 0, wwc = 0, wdR = 0, wdM = 0;
	for (const auto &o : outs)
	{
		const auto t3 = *std::find_if(T3.begin(), T3.end(), [&](const ht68::NonRot &r) { return r.ec == o.ec; });
		const auto t5 = *std::find_if(T5.begin(), T5.end(), [&](const ht68::Rot &r) { return r.ec == o.ec; });
		const double rR = Rel(o.R, t3.R_km), rM = Rel(o.M, t3.M_sun), rRg = Rel(o.Rg_R, t5.Rg_over_R), rOm = Rel(o.Omega, t5.Omega_s);
		wR = std::max(wR, rR);
		wM = std::max(wM, rM);
		wRg = std::max(wRg, rRg);
		wOm = std::max(wOm, rOm);
		char b[400];
		snprintf(b, sizeof(b), "    %.2e    %5zu | %7.3f  %6.2f  %.1e | %.4f  %.3f  %.1e | %.4f  %.3f  %.1e | %.3e  %.2e  %.1e\n", o.ec, o.n, o.R,
				 t3.R_km, rR, o.M, t3.M_sun, rM, o.Rg_R, t5.Rg_over_R, rRg, o.Omega, t5.Omega_s, rOm);
		std::cout << b;
	}
	std::cout << "\n    E_c[g/cm^3] | w_s/W prod   tab      rel    | w_c/W prod  tab     rel    | dR/R prod  tab     rel    | dM/M prod  tab     rel    | shell/dM\n";
	for (const auto &o : outs)
	{
		const auto t5 = *std::find_if(T5.begin(), T5.end(), [&](const ht68::Rot &r) { return r.ec == o.ec; });
		const double rws = Rel(o.ws, t5.ws_over_W), rwc = Rel(o.wc, t5.wc_over_W), rdR = Rel(o.dR_R, t5.dR_over_R), rdM = Rel(o.dM_M, t5.dM_over_M);
		wws = std::max(wws, rws);
		wwc = std::max(wwc, rwc);
		wdR = std::max(wdR, rdR);
		wdM = std::max(wdM, rdM);
		char b[400];
		snprintf(b, sizeof(b), "    %.2e    | %.4e  %.2e  %.1e | %.4f  %.3f  %.1e | %.4f  %.3f  %.1e | %.4f  %.4f  %.1e | %.1e\n", o.ec, o.ws,
				 t5.ws_over_W, rws, o.wc, t5.wc_over_W, rwc, o.dR_R, t5.dR_over_R, rdR, o.dM_M, t5.dM_over_M, rdM, o.shell_over_dM);
		std::cout << b;
	}
	const bool all = outs.size() == T3.size();
	Report("K2a all eight HT68 configurations reconstructed from the printed EOS", all, std::to_string(outs.size()) + "/8");
	Report("K2z every HT68 star is a complete ADR-0009 star (SURFACE_REACHED, last node p == p_cut, finite surface data)",
		   all && std::all_of(outs.begin(), outs.end(), [](const Out &o) { return o.complete; }), "");
	{
		std::cout << "\n    F. ADR-0008 Validation F on the smooth HW EOS — production deltaM_hat vs the oracle forms (relative):\n"
					 "    E_c[g/cm^3] | superseded differential | Stieltjes profile K=4 | Stieltjes EOS-knot K=2 | secant\n";
		double wf = 0.0, ws_ = 0.0, wk = 0.0;
		for (const auto &o : outs)
		{
			char b[300];
			snprintf(b, sizeof(b), "    %.2e    |  %.3e              |  %.3e            |  %.3e             |  %.3e\n", o.ec, o.f_diff, o.f_stj, o.f_knot, o.f_sec);
			std::cout << b;
			wf = std::max(wf, o.f_diff);
			ws_ = std::max(ws_, o.f_stj);
			wk = std::max(wk, o.f_knot);
		}
		double wf_auth = 0.0;
		std::size_t n_auth = 0;
		for (const auto &o : outs)
			for (double ec : kF_AuthenticatedEc)
				if (o.ec == ec)
				{
					wf_auth = std::max(wf_auth, o.f_diff);
					++n_auth;
				}
		Report("F1 smooth-EOS measure-vs-differential equivalence on deltaM_hat at the AUTHENTICATED scope of ADR-0008 Validation F (eps_c = 3e14, 1e15, 3e15; bound 2e-5)",
			   all && n_auth == 3 && wf_auth <= kF_Smooth, "worst rel=" + Sci(wf_auth) + "  bound=" + Sci(kF_Smooth));
		std::cout << "     RECORD — F1r the same equivalence over all eight HT68 configurations: worst rel " << Sci(wf)
				  << (wf <= kF_Smooth ? " <= " : " > ") << Sci(kF_Smooth, 0)
				  << " (not asserted: outside the bound's authenticated scope; the superseded form's deficit grows with the density gradient,\n"
				  << "      production is the measure-complete side and the independent oracles confirm it at F2).\n";
		Report("F2 the INDEPENDENT Stieltjes oracles (profile and EOS-knot partitions) agree with production on the smooth EOS within the same bound",
			   all && ws_ <= kF_Smooth && wk <= kF_Smooth, "worst rel profile=" + Sci(ws_) + "  knots=" + Sci(wk) + "  (reported beside F1; same bound)");
	}
	Report("K2b Table 3 radii agree (EOS reconstruction, non-rotating)", all && wR <= kK_Published, "worst rel=" + Sci(wR));
	Report("K2c Table 3 masses agree (EOS reconstruction, non-rotating)", all && wM <= kK_Published, "worst rel=" + Sci(wM));
	Report("K2d Table 5 R_g/R = sqrt(I/M)/R agrees (first order)", all && wRg <= kK_Published, "worst rel=" + Sci(wRg));
	Report("K2e Table 5 Omega = c sqrt(M/R^3) agrees (unit sanity)", all && wOm <= kK_Published, "worst rel=" + Sci(wOm));
	Report("K2f Table 5 omega_s/Omega agrees (first order, surface frame dragging)", all && wws <= kK_Published, "worst rel=" + Sci(wws));
	Report("K2g Table 5 omega_c/Omega agrees (first order, central frame dragging)", all && wwc <= kK_Published, "worst rel=" + Sci(wwc));
	Report("K2h Table 5 deltaR/R = xi_0(R)/R agrees (SECOND order)", all && wdR <= kK_Published, "worst rel=" + Sci(wdR));
	Report("K2i Table 5 deltaM/M = [m0(R) + J^2/R^3]/M agrees (SECOND order)", all && wdM <= kK_Published, "worst rel=" + Sci(wdM));
	std::cout << "    (HT68 surface: P = 0 on the printed table; CompactStar's cutoff max(1e-15 p_c, table floor) omits an outer\n"
				 "     layer of metres on these models — the shell/dM column shows the surface term is negligible here, unlike\n"
				 "     the homogeneous case. Transcription uncertainty: see hartle_thorne_1968_hw_eos.hpp.)\n";
}

int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4D Experiment K — governed O(Omega^2) monopole response vs PUBLISHED second-order results\n\n";
	const fs::path wrk = fs::temp_directory_path() / ("compactstar_hartle_monopole_4d_published_" + std::to_string(::getpid()));
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	RunChandrasekharMiller();
	RunHartleThorne(wrk);

	fs::remove_all(wrk);
	std::cout << "\nSLOW-ROTATION DISCLAIMER: these comparisons establish the O(Omega^2) COEFFICIENTS; HT68's\n"
				 "Omega^2 = M/R^3 configurations are themselves slow-rotation truncations.\n";
	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
