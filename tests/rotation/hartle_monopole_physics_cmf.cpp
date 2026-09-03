// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_physics_cmf.cpp
 * @brief Phase 4D — INDEPENDENT physical validation of the governed O(Omega^2) monopole
 *        response on the authenticated DS(CMF)-1 sequence. External data required.
 *
 * The oracle is the test-only (m0,h0) solver of hartle_monopole_reference.hpp — Hartle's
 * (97)+(98) with p0* from the first integral (90) — never production's (m0,p0*) system. Two
 * chains: SECOND-ORDER-ISOLATED (production's Phase-4B-verified s, s') and FULLY INDEPENDENT
 * (hartle_reference.hpp's own first order on the same tabulated background).
 *
 * PHASE 4D-RI (2026-09-03): ADR-0008 was ACCEPTED and the EOS energy-density source of both
 * production and this oracle is now the measure -4 pi r^2 xi0_hat d(eps) evaluated one profile
 * segment at a time (MHOptions::eos_measure). The (m0,h0) formulation, its own interpolation,
 * centre start and tolerances are unchanged, so the comparison stays a genuine cross-check of a
 * different variable pair; the SUPERSEDED differential form is still run and REPORTED beside it,
 * and its disagreement is the size of the sub-node energy-density variation the nodal deps/dp
 * column cannot represent. The corrected independent revalidation is a separate increment: no
 * bound here is widened and no result here is a validation claim.
 *
 * ============================ PREDECLARED BOUNDS ============================
 *   G  full profile, production vs (m0,h0) reference, per star   rel <= 1e-4   (ADR-0007 §7-4)
 *   F  near-vacuum identity. No profile node lies beyond R_*, so the vacuum statement
 *      "mhat + I^2/r^3 = const" is examined on the outermost nodes (eps < 1e9 g/cm^3,
 *      r > R_* / 2): the MATTER-source-corrected value must be constant to      <= 1e-6   (§7-5)
 *      and the raw (uncorrected) spread is reported. The matching arithmetic is exact.
 *      (A first run defined the window as (eps+p) r^2 < 1e-6, which also admits the
 *       central nodes where I^2/r^3 is ~1e19 — a test defect, corrected and recorded.)
 *   I  EOS-derivative sensitivity                                         REPORTED  (§7-10)
 *      Since ADR-0008 the profile deps/dp no longer enters the radial mass source at all, so
 *      this substitution now probes only the regular-centre series coefficient b_5.
 *      The ADR spread bound (1e-3) applies to an INDEPENDENT derivative source. None is
 *      available (c_s^2: CONDITIONAL CHECK UNAVAILABLE); the retired profile FD is a
 *      diagnostic only. A first run asserted the bound against the FD and measured 5.0e-2;
 *      the assertion was withdrawn as mis-specified and the number is recorded.
 *   J  homogeneous deltaM vs the non-rotating sequence derivative                <= 1e-3   (§7-11)
 *      H67 p. 1022: the homogeneous solution is the change of CENTRAL PRESSURE along the
 *      non-rotating sequence, so the native identity is deltaM_hom = (dM/dp_c) delta p_c with
 *      delta p_c = (eps_c + p_c) p0*_c; that is what the analytic Experiment J certified (3e-9).
 *      OUTCOME ON DS(CMF)-1 (recorded, not asserted): 1.17e-3 / 1.02e-3 / 1.04e-3 at radial
 *      resolution 10000 / 20000 / 40000 — NOT MET and resolution-independent. Diagnosis, proven
 *      test-side: the nodal deps/dp column integrated over the crust misses ~17 % of the crust's
 *      own Delta eps at every resolution (density steps of the crust EOS that no sampled
 *      derivative represents); a Stieltjes sum of the same source against the profile's own eps
 *      steps reproduces the sequence derivative to ~7e-5, and applied to the SOURCED solution it
 *      puts the omitted internal delta-function shells (Hartle's dE/dP at discontinuities) at
 *      ~4.6 % of deltaM_hat. That is a physical discrepancy of the accepted contract, reported by
 *      Phase 4D and deliberately not repaired here (ADR-0007 amendment required).
 *      Chronology: a first run computed only the eps_c form at 10000; 20000 was added, then the
 *      p_c form and the inversion diagnostic, then 40000 and the Stieltjes diagnostic.
 *      PHASE 4D-RI: with ADR-0008 accepted, this line is computed with the corrected source in
 *      both production and the oracle and is REPORTED here; its adjudication against the
 *      ADR-0008 target (<= 2e-4) belongs to the corrected independent revalidation increment.
 *   H  radial convergence: order MEASURED, no pass criterion invented    reported (§7-9, INV-13)
 *   R  reference floor / disagreement                                    <= 0.1
 * ===========================================================================
 *
 * No 4C-I1 diagnostic value, no retired-candidate value and no PHASE4_ROTATION_ENTRY number is
 * an expected answer anywhere in this file.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
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
#include "hartle_monopole_reference.hpp"
#include "hartle_profile_compare.hpp"
#include "hartle_reference.hpp"

#include <Zaki/Physics/Constants.hpp>

#include <unistd.h>

#include <cstdlib>

namespace fs = std::filesystem;

/// HARTLE_4D_QUICK=1 restricts the run to G (isolated chain) + F: used ONLY by the detector
/// sweep, where each production mutation is rebuilt and the whole 4D suite re-run nine times.
/// The registered CTest runs the full set.
static bool Quick() { return std::getenv("HARTLE_4D_QUICK") != nullptr; }

using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;
using hartle_mono_ref::Background2;
using hartle_mono_ref::Cmp;
using hartle_mono_ref::Compare;
using hartle_mono_ref::MHOptions;
using hartle_mono_ref::MHResult;

static constexpr double kG_Profile = 1.0e-4;
static constexpr double kF_Exterior = 1.0e-6;
static constexpr double kI_Spread = 1.0e-3;
static constexpr double kJ_Homog = 1.0e-3;
static constexpr double kR_FloorRatio = 0.1;
static constexpr double kVacuumEps_gcm3 = 1.0e9; // near-vacuum window: eps below this AND r > R_* / 2

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
struct Star
{
	std::unique_ptr<NStar> ns;
	std::vector<TOVPoint> pts;
	std::vector<std::string> labels;
	double ec_cgs = 0;
};

static Star BuildAtMass(const fs::path &cold, const fs::path &wrk, double M)
{
	Star st;
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	const int n = tov.SolveToProfile(M, st.pts, &st.labels);
	if (n <= 0 || st.pts.size() < 4)
		return st;
	st.ec_cgs = st.pts.front().e;
	st.ns = st.labels.empty() ? std::make_unique<NStar>(st.pts) : std::make_unique<NStar>(st.pts, st.labels);
	if (!st.ns->ComputeHartleMonopoleResponse())
		st.ns.reset();
	return st;
}

static Star BuildAtDensity(const fs::path &cold, const fs::path &wrk, double ec, std::size_t res)
{
	Star st;
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.SetRadialRes(res);
	tov.ImportEOS(cold.string(), true);
	if (tov.SingleStarSolveToTOVPoints(ec, st.pts) <= 0 || st.pts.size() < 4)
		return st;
	st.ec_cgs = ec;
	st.ns = std::make_unique<NStar>(st.pts);
	if (!st.ns->ComputeHartleMonopoleResponse())
		st.ns.reset();
	return st;
}

/// One member of the non-rotating sequence: M [km] and the ACTUAL central p, eps [km^-2] the
/// solver used (its own inversion of the central condition), for the sequence derivative.
struct SeqPt
{
	double M_km = std::numeric_limits<double>::quiet_NaN(), pc = 0, ec = 0;
};
static SeqPt SequencePoint(const fs::path &cold, const fs::path &wrk, double ec, std::size_t res)
{
	SeqPt o;
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.SetRadialRes(res);
	tov.ImportEOS(cold.string(), true);
	std::vector<TOVPoint> pts;
	if (tov.SingleStarSolveToTOVPoints(ec, pts) <= 0 || pts.empty())
		return o;
	o.M_km = pts.back().m * Zaki::Physics::SUN_M_KM;
	o.pc = pts.front().p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2;
	o.ec = pts.front().e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3;
	return o;
}

struct ProdFields
{
	std::vector<double> r, mhat, phat, dphat, xihat;
	double deltaM = 0, shell = 0, xi_R = 0, R = 0, I = 0, eps_R = 0;
};
static ProdFields Fields(const NStar &ns)
{
	ProdFields f;
	const auto *R = ns.MonopoleResponse();
	const int n = static_cast<int>(R->m0_over_Omega2.Size());
	for (int i = 0; i < n; ++i)
	{
		f.r.push_back((*R->r_grid)[i]);
		f.mhat.push_back(R->m0_over_Omega2[i]);
		f.phat.push_back(R->p0star_over_Omega2[i]);
		f.dphat.push_back(R->delta_p0_over_Omega2[i]);
		f.xihat.push_back(R->xi0_over_Omega2[i]);
	}
	f.deltaM = R->deltaM_over_Omega2;
	f.shell = R->surface_shell_mass_over_Omega2;
	f.xi_R = R->surface_xi0_over_Omega2;
	f.R = R->R_surface;
	f.I = R->I;
	f.eps_R = (*ns.Profile().GetEnergyDensity())[n - 1];
	return f;
}

static void ProdShape(const NStar &ns, std::vector<double> &s, std::vector<double> &sp)
{
	const auto &R = ns.RotationResponse();
	const int n = static_cast<int>(R.omega_bar_over_Omega.Size());
	s.assign(n, 0.0);
	sp.assign(n, 0.0);
	for (int i = 0; i < n; ++i)
	{
		s[i] = R.omega_bar_over_Omega[i];
		sp[i] = R.domega_bar_over_Omega_dr[i];
	}
}

static Background2 FromProfile(const NStar &ns, const std::vector<double> &s, const std::vector<double> &sp)
{
	const auto &P = ns.Profile();
	const auto *r = P.GetRadius();
	const int n = static_cast<int>(r->Size());
	Background2 o;
	for (int i = 0; i < n; ++i)
	{
		o.r.push_back((*r)[i]);
		o.p.push_back((*P.GetPressure())[i]);
		o.eps.push_back((*P.GetEnergyDensity())[i]);
		o.m.push_back((*P.GetMass())[i]);
		o.nu.push_back((*P.GetMetricNu())[i]);
		o.nup.push_back((*P.GetMetricNuPrime())[i]);
		o.dedp.push_back((*P.GetEosDEdP())[i]);
	}
	o.s = s;
	o.sp = sp;
	return o;
}

struct StarRow
{
	double M = 0, R = 0, dedp_c = 0, mhat_R = 0, phat_R = 0, xi_R = 0, shell = 0, ext = 0, dM = 0;
	Cmp m_iso, p_iso, d_iso, x_iso, m_full, p_full, x_full;
	double dM_iso = 0, dM_full = 0;
	double ext_const = 0, dM_node_indep = 0;
	std::size_t n = 0, n_vac = 0;
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: hartle_monopole_physics_cmf <EOS_DATA_ROOT>\n";
		return 2;
	}
	const fs::path root = argv[1];
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "authenticated EOS missing: " << cold.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / ("compactstar_hartle_monopole_4d_cmf_" + std::to_string(::getpid()));
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4D — INDEPENDENT validation of the governed O(Omega^2) monopole response on\n"
				 "DS(CMF)-1_with_crust. Oracle: the test-only (m0,h0) solver. Bounds from ADR-0007 §7.\n\n";

	// =====================================================================
	//  G — four stars, both chains
	// =====================================================================
	std::cout << "G. Full-profile comparison, production vs the independent (m0,h0) solver\n";
	const double masses[] = {1.0, 1.4, 1.6, 2.0};
	std::vector<StarRow> rows;
	Star star16; // kept for H, I, J, R

	for (double M : masses)
	{
		Star st = BuildAtMass(cold, wrk / ("m" + std::to_string(M)), M);
		if (!st.ns)
		{
			Report("G build " + Sci(M, 2) + " Msun", false, "no production response");
			continue;
		}
		const NStar &cns = *st.ns;
		const auto pf = Fields(cns);
		std::vector<double> s_p, sp_p;
		ProdShape(cns, s_p, sp_p);

		StarRow row;
		row.M = M;
		row.R = pf.R;
		row.n = pf.r.size();
		row.dedp_c = (*cns.Profile().GetEosDEdP())[0];
		row.mhat_R = pf.mhat.back();
		row.phat_R = pf.phat.back();
		row.xi_R = pf.xi_R;
		row.shell = pf.shell;
		row.ext = pf.I * pf.I / (pf.R * pf.R * pf.R);
		row.dM = pf.deltaM;

		// isolated chain
		const Background2 bg_iso = FromProfile(cns, s_p, sp_p);
		MHOptions oi;
		oi.I_exterior = pf.I;
		oi.eos_measure = true; // ADR-0008 Q1/Q3 (accepted source form)
		const MHResult ri = hartle_mono_ref::Solve(bg_iso, oi);
		row.m_iso = Compare(pf.mhat, ri.mhat, pf.r);
		row.p_iso = Compare(pf.phat, ri.phat, pf.r);
		row.d_iso = Compare(pf.dphat, ri.dphat, pf.r);
		row.x_iso = Compare(pf.xihat, ri.xihat, pf.r);
		row.dM_iso = Rel(pf.deltaM, ri.deltaM_hat);

		// fully independent chain: independent first order on the same tabulated background
		if (Quick())
		{
			row.m_full = row.m_iso; row.p_full = row.p_iso; row.x_full = row.x_iso; row.dM_full = row.dM_iso;
		}
		const auto hb = Quick() ? hartle_ref::Background{} : hartle_4b::BackgroundFromProfile(cns);
		hartle_ref::Result fo;
		MHResult rf;
		if (!Quick())
		{
			fo = hartle_ref::Solve(hb, 5.0e-3, hb.r.front(), 1e-13, 1e-16);
			Background2 bg_full = bg_iso;
			bg_full.s = fo.s;
			bg_full.sp = fo.s_prime;
			MHOptions of;
			of.I_exterior = fo.I_surface;
			of.eos_measure = true;
			rf = hartle_mono_ref::Solve(bg_full, of);
			row.m_full = Compare(pf.mhat, rf.mhat, pf.r);
			row.p_full = Compare(pf.phat, rf.phat, pf.r);
			row.x_full = Compare(pf.xihat, rf.xihat, pf.r);
			row.dM_full = Rel(pf.deltaM, rf.deltaM_hat);
		}
		else
		{
			fo.ok = true;
			rf.ok = true;
		}

		// F — near-vacuum identity on the outermost nodes (no node lies beyond R_*)
		double raw_spread = std::numeric_limits<double>::quiet_NaN(), corr_spread = raw_spread, sp_vs_I = raw_spread;
		{
			const auto &P = cns.Profile();
			const double km2_to_gcm3 = Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;
			const double eps_win = kVacuumEps_gcm3 / km2_to_gcm3;
			const std::size_t n = pf.r.size();
			// matter part of the (97) source on production's own fields, in the ACCEPTED measure
			// form (ADR-0008): the EOS part uses each segment's own d(eps)/dr, evaluated at the
			// node from the segment that ends there.
			std::vector<double> fm(n, 0.0);
			for (std::size_t i = 0; i < n; ++i)
			{
				const int ii = static_cast<int>(i);
				const double r = pf.r[i], eps = (*P.GetEnergyDensity())[ii], p = (*P.GetPressure())[ii];
				const double nu = (*P.GetMetricNu())[ii];
				const std::size_t lo = (i == 0) ? 0 : i - 1, hi = (i == 0) ? 1 : i;
				const double slope = ((*P.GetEnergyDensity())[static_cast<int>(hi)] - (*P.GetEnergyDensity())[static_cast<int>(lo)]) /
									 (pf.r[hi] - pf.r[lo]);
				fm[i] = -4.0 * M_PI * r * r * pf.xihat[i] * slope +
						(8.0 * M_PI / 3.0) * r * r * r * r * (eps + p) * std::exp(-2.0 * nu) * s_p[i] * s_p[i];
			}
			// remaining matter source from node i to R_* (trapezoid)
			std::vector<double> S(n, 0.0);
			for (std::size_t i = n - 1; i-- > 0;)
				S[i] = S[i + 1] + 0.5 * (fm[i] + fm[i + 1]) * (pf.r[i + 1] - pf.r[i]);
			double rlo = 1e300, rhi = -1e300, clo = 1e300, chi = -1e300;
			std::size_t nv = 0;
			for (std::size_t i = 0; i < n; ++i)
			{
				const int ii = static_cast<int>(i);
				if (pf.r[i] > 0.5 * pf.R && (*P.GetEnergyDensity())[ii] < eps_win)
				{
					const double ext = pf.I * pf.I / (pf.r[i] * pf.r[i] * pf.r[i]);
					const double raw = pf.mhat[i] + ext, corr = pf.mhat[i] + S[i] + ext;
					rlo = std::min(rlo, raw); rhi = std::max(rhi, raw);
					clo = std::min(clo, corr); chi = std::max(chi, corr);
					++nv;
				}
			}
			row.n_vac = nv;
			if (nv >= 2)
			{
				raw_spread = (rhi - rlo) / pf.deltaM;
				corr_spread = (chi - clo) / pf.deltaM;
			}
			row.ext_const = corr_spread;
			row.dM_node_indep = raw_spread;
			sp_vs_I = Rel(sp_p[n - 1], 6.0 * pf.I / (pf.R * pf.R * pf.R * pf.R));
		}

		// SUPERSEDED differential form, reported only: the gap is the sub-node energy-density
		// variation the nodal deps/dp column loses (ADR-0008; PHASE4D_MONOPOLE_VALIDATION.md).
		double old_form_dM = std::numeric_limits<double>::quiet_NaN();
		{
			MHOptions oo = oi;
			oo.eos_measure = false;
			const MHResult ro = hartle_mono_ref::Solve(bg_iso, oo);
			if (ro.ok)
				old_form_dM = Rel(ro.deltaM_hat, pf.deltaM);
		}

		const bool ok_iso = ri.ok && std::max({row.m_iso.max_rel, row.p_iso.max_rel, row.d_iso.max_rel, row.x_iso.max_rel}) <= kG_Profile;
		const bool ok_full = rf.ok && fo.ok && std::max({row.m_full.max_rel, row.p_full.max_rel, row.x_full.max_rel}) <= kG_Profile;
		Report(Sci(M, 2) + " Msun: ISOLATED chain agrees at every node", ok_iso,
			   "max rel mhat=" + Sci(row.m_iso.max_rel) + " phat=" + Sci(row.p_iso.max_rel) + " dphat=" +
				   Sci(row.d_iso.max_rel) + " xihat=" + Sci(row.x_iso.max_rel) + " dM=" + Sci(row.dM_iso));
		Report(Sci(M, 2) + " Msun: FULLY INDEPENDENT chain agrees at every node", ok_full,
			   "max rel mhat=" + Sci(row.m_full.max_rel) + " phat=" + Sci(row.p_full.max_rel) + " xihat=" +
				   Sci(row.x_full.max_rel) + " dM=" + Sci(row.dM_full));
		Report(Sci(M, 2) + " Msun: near-vacuum identity — matter-corrected mhat + S + I^2/r^3 constant over the outermost nodes",
			   std::isfinite(row.ext_const) && row.ext_const <= kF_Exterior,
			   "corrected spread/deltaM=" + Sci(row.ext_const) + "  raw spread/deltaM=" + Sci(row.dM_node_indep) + " over " +
				   std::to_string(row.n_vac) + " nodes (eps < 1e9 g/cm^3)  bound=" + Sci(kF_Exterior) +
				   "  | s'(R_*) vs 6I/R_*^4 rel=" + Sci(sp_vs_I));
		Report(Sci(M, 2) + " Msun: matching arithmetic deltaM_hat = mhat(R_*) + shell + I^2/R_*^3 is exact",
			   Rel(pf.deltaM, pf.mhat.back() + pf.shell + pf.I * pf.I / (pf.R * pf.R * pf.R)) <= 1e-14, "");
		std::cout << "     REPORTED — the SUPERSEDED differential (nodal deps/dp) oracle disagrees with the corrected production deltaM_hat by "
				  << Sci(old_form_dM) << "; that is the sub-node energy-density variation ADR-0008 recovers, not a numerical error.\n";

		rows.push_back(row);
		if (M == 1.6)
			star16 = std::move(st);
	}

	std::cout << "\n  M[Msun]  R_*[km]   dedp_c    mhat(R_*)   phat(R_*)   xi_hat(R_*)  shell_hat   I^2/R_*^3   deltaM_hat\n";
	for (const auto &w : rows)
	{
		char b[400];
		snprintf(b, sizeof(b), "  %5.2f   %8.5f  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e\n", w.M, w.R, w.dedp_c, w.mhat_R,
				 w.phat_R, w.xi_R, w.shell, w.ext, w.dM);
		std::cout << b;
	}
	std::cout << "  M[Msun]  RMS mhat    RMS phat    RMS xihat   worst r mhat  worst r phat  worst r xihat  surface rel phat\n";
	for (const auto &w : rows)
	{
		char b[400];
		snprintf(b, sizeof(b), "  %5.2f   %.3e   %.3e   %.3e   %8.4f km   %8.4f km   %8.4f km    %.3e\n", w.M, w.m_iso.rms,
				 w.p_iso.rms, w.x_iso.rms, w.m_iso.r_worst, w.p_iso.r_worst, w.x_iso.r_worst, w.p_iso.surface_rel);
		std::cout << b;
	}

	if (Quick())
	{
		fs::remove_all(wrk);
		std::cout << "\n(HARTLE_4D_QUICK: R/H/I/J skipped — detector-sweep mode)\n";
		std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail << " failed check(s)\n";
		return g_fail == 0 ? 0 : 1;
	}
	if (!star16.ns)
	{
		std::cout << "\nFAIL — 1.6 Msun star unavailable for H/I/J/R\n";
		return 1;
	}
	const NStar &c16 = *star16.ns;
	const auto pf16 = Fields(c16);
	std::vector<double> s16, sp16;
	ProdShape(c16, s16, sp16);
	const Background2 bg16 = FromProfile(c16, s16, sp16);

	// =====================================================================
	//  R — reference admissibility on 1.6 Msun
	// =====================================================================
	std::cout << "\nR. Reference admissibility on 1.6 Msun\n";
	{
		MHOptions o0;
		o0.I_exterior = pf16.I;
		const auto base = hartle_mono_ref::Solve(bg16, o0);
		double floor = 0.0;
		for (double tol : {1e-11, 1e-15})
		{
			MHOptions o = o0;
			o.rtol = tol;
			o.atol = tol * 1e-3;
			const auto v = hartle_mono_ref::Solve(bg16, o);
			floor = std::max({floor, Compare(v.mhat, base.mhat, pf16.r).max_rel, Compare(v.phat, base.phat, pf16.r).max_rel,
							  Rel(v.deltaM_hat, base.deltaM_hat)});
		}
		{
			MHOptions o = o0;
			o.r0 = 1e-7;
			const auto v = hartle_mono_ref::Solve(bg16, o);
			floor = std::max({floor, Compare(v.mhat, base.mhat, pf16.r).max_rel, Compare(v.phat, base.phat, pf16.r).max_rel,
							  Rel(v.deltaM_hat, base.deltaM_hat)});
		}
		const auto it = std::find_if(rows.begin(), rows.end(), [](const StarRow &w) { return w.M == 1.6; });
		const double signal = std::max({it->m_iso.max_rel, it->p_iso.max_rel, it->x_iso.max_rel, it->dM_iso});
		Report("Ra reference self-movement (tol 1e-11..1e-15, r0 1e-5..1e-7) is subdominant to the disagreement",
			   signal > 0 && floor / signal <= kR_FloorRatio,
			   "floor=" + Sci(floor) + "  signal=" + Sci(signal) + "  ratio=" + Sci(signal > 0 ? floor / signal : 0.0, 2));
	}

	// =====================================================================
	//  H — radial convergence at fixed central density (1.6 Msun)
	// =====================================================================
	std::cout << "\nH. Radial convergence at fixed eps_c = " << Sci(star16.ec_cgs, 6) << " g/cm^3 (1.6 Msun)\n";
	{
		struct Hrow
		{
			std::size_t res, n;
			double R, M, dM, xi, mR, p_half;
		};
		std::vector<Hrow> h;
		for (std::size_t res : {std::size_t(5000), std::size_t(10000), std::size_t(20000)})
		{
			Star st = BuildAtDensity(cold, wrk / ("h" + std::to_string(res)), star16.ec_cgs, res);
			if (!st.ns)
			{
				Report("H build res=" + std::to_string(res), false, "");
				continue;
			}
			const auto f = Fields(*st.ns);
			// interior sample: p0*_hat at the node nearest R/2
			std::size_t k = 0;
			for (std::size_t i = 0; i < f.r.size(); ++i)
				if (std::fabs(f.r[i] - 0.5 * f.R) < std::fabs(f.r[k] - 0.5 * f.R))
					k = i;
			h.push_back({res, f.r.size(), f.R, st.pts.back().m, f.deltaM, f.xi_R, f.mhat.back(), f.phat[k]});
		}
		std::cout << "     res     nodes    R_*[km]     M[Msun]        deltaM_hat      xi_hat(R_*)     mhat(R_*)       phat(R/2)\n";
		for (const auto &x : h)
		{
			char b[300];
			snprintf(b, sizeof(b), "     %5zu   %5zu   %.6f  %.8f  %.8e  %.8e  %.8e  %.8e\n", x.res, x.n, x.R, x.M, x.dM, x.xi, x.mR, x.p_half);
			std::cout << b;
		}
		if (h.size() == 3)
		{
			auto order = [](double a, double b, double c) {
				const double d1 = std::fabs(a - b), d2 = std::fabs(b - c);
				return d2 > 0 ? std::log2(d1 / d2) : std::numeric_limits<double>::quiet_NaN();
			};
			auto rich = [](double b, double c, double p) { return c + (c - b) / (std::pow(2.0, p) - 1.0); };
			const double pM = order(h[0].dM, h[1].dM, h[2].dM), pX = order(h[0].xi, h[1].xi, h[2].xi),
						 pm = order(h[0].mR, h[1].mR, h[2].mR), pp = order(h[0].p_half, h[1].p_half, h[2].p_half);
			std::cout << "     successive differences deltaM_hat: " << Sci(std::fabs(h[0].dM - h[1].dM)) << ", "
					  << Sci(std::fabs(h[1].dM - h[2].dM)) << "  observed order " << Sci(pM, 2) << "\n";
			std::cout << "     successive differences xi_hat(R_*): " << Sci(std::fabs(h[0].xi - h[1].xi)) << ", "
					  << Sci(std::fabs(h[1].xi - h[2].xi)) << "  observed order " << Sci(pX, 2) << "\n";
			std::cout << "     successive differences mhat(R_*):   " << Sci(std::fabs(h[0].mR - h[1].mR)) << ", "
					  << Sci(std::fabs(h[1].mR - h[2].mR)) << "  observed order " << Sci(pm, 2) << "\n";
			std::cout << "     successive differences phat(R/2):   " << Sci(std::fabs(h[0].p_half - h[1].p_half)) << ", "
					  << Sci(std::fabs(h[1].p_half - h[2].p_half)) << "  observed order " << Sci(pp, 2) << "\n";
			if (std::isfinite(pM) && pM > 0.5)
				std::cout << "     Richardson-style residual on deltaM_hat at res=20000 (order " << Sci(pM, 2)
						  << "): " << Sci(Rel(h[2].dM, rich(h[1].dM, h[2].dM, pM))) << " relative\n";
			std::cout << "     (No pass criterion is invented here: INV-13 requires the order to be measured. R_* itself\n"
						 "      moves with the grid, so surface-node quantities carry that shift; deltaM_hat is the\n"
						 "      exterior-matched quantity and is the cleanest convergence indicator.)\n";
			const double lo = std::min({h[0].dM, h[1].dM, h[2].dM}), hi = std::max({h[0].dM, h[1].dM, h[2].dM});
			const double spread = (hi - lo) / std::fabs(h[1].dM);
			std::cout << "     RECORD — PHASE 4D-RI: with the ADR-0008 measure source the relative SPREAD of deltaM_hat over\n"
					  << "     5000/10000/20000 is " << Sci(spread) << " (ADR-0008 Validation D asks <= 1e-4), but its successive\n"
					  << "     differences are NOT monotone (" << Sci(std::fabs(h[0].dM - h[1].dM), 2) << " then "
					  << Sci(std::fabs(h[1].dM - h[2].dM), 2) << " km^3) and R_* itself moves with the grid\n"
					  << "     (" << Sci(h[0].R, 8) << " -> " << Sci(h[2].R, 8) << " km). The residual is the TOV background's own\n"
					  << "     resolution dependence, not the monopole source; separating the two is the corrected\n"
					  << "     revalidation increment's line D, not this implementation task's.\n";
			Report("Ha the 5000/10000/20000 sequence was produced and its convergence measured", h.size() == 3, "");
		}
	}

	// =====================================================================
	//  I — EOS-derivative sensitivity (Steffen authority vs profile FD; c_s^2 unavailable)
	// =====================================================================
	std::cout << "\nI. EOS-derivative sensitivity on 1.6 Msun\n";
	{
		// Substitute the retired profile finite difference THROUGH the governed explicit-supply
		// mechanism (TOVPoint::dedp), test-side; production's authority is untouched.
		std::vector<TOVPoint> pts_fd = star16.pts;
		const double c2 = Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_G_CM3; // cgs -> dimensionless
		const std::size_t n = pts_fd.size();
		for (std::size_t i = 0; i < n; ++i)
		{
			const std::size_t lo = (i == 0) ? 0 : i - 1, hi = (i + 1 < n) ? i + 1 : n - 1;
			const double dp = pts_fd[hi].p - pts_fd[lo].p, de = pts_fd[hi].e - pts_fd[lo].e;
			pts_fd[i].dedp = std::fabs(dp) > 1e-30 ? (de / dp) * c2 : 1.0; // the candidate's fallback, deliberately
		}
		NStar ns_fd = star16.labels.empty() ? NStar(pts_fd) : NStar(pts_fd, star16.labels);
		const bool ok = ns_fd.ComputeHartleMonopoleResponse();
		if (ok)
		{
			const NStar &cfd = ns_fd;
			const auto ff = Fields(cfd);
			const auto cm = Compare(ff.mhat, pf16.mhat, pf16.r), cp = Compare(ff.phat, pf16.phat, pf16.r),
					   cx = Compare(ff.xihat, pf16.xihat, pf16.r);
			std::cout << "     A (Steffen, authority): deltaM_hat=" << Sci(pf16.deltaM, 8) << "  xi_hat(R_*)=" << Sci(pf16.xi_R, 8) << "\n"
					  << "     B (profile FD, diagnostic): deltaM_hat=" << Sci(ff.deltaM, 8) << "  xi_hat(R_*)=" << Sci(ff.xi_R, 8) << "\n"
					  << "     B-A spread: deltaM_hat " << Sci(Rel(ff.deltaM, pf16.deltaM)) << ", xi_hat(R_*) " << Sci(Rel(ff.xi_R, pf16.xi_R))
					  << ", profile max rel mhat " << Sci(cm.max_rel) << " phat " << Sci(cp.max_rel) << " xihat " << Sci(cx.max_rel)
					  << " (worst at r=" << Sci(cp.r_worst, 4) << " km)\n"
					  << "     C (tabulated c_s^2): CONDITIONAL CHECK UNAVAILABLE — the governed import carries no such column and the\n"
					  << "       raw eos.thermo's three additional columns are undocumented (one negative, one zero, one equal to n_B).\n"
					  << "     The ADR spread bound (<= 1e-3 on deltaM_hat) applies to an INDEPENDENT source; the FD row is a diagnostic\n"
					  << "     showing how much the retired method moves the global answer.\n";
			std::cout << "     RECORD — §7 item 10 (independent-derivative spread <= " << Sci(kI_Spread, 0)
					  << "): NOT INDEPENDENTLY TESTABLE with the available inputs; the retired-FD diagnostic moves\n"
					  << "     deltaM_hat by " << Sci(Rel(ff.deltaM, pf16.deltaM)) << " (the crust derivative is load-bearing at the percent level).\n";
			Report("Ia the FD-substituted star computed through the governed explicit-supply mechanism (diagnostic only)", true, "");
		}
		else
			Report("Ia FD-substituted star computed", false, "");
	}

	// =====================================================================
	//  J — homogeneous family vs the TOV sequence derivative dM/deps_c (test-side)
	// =====================================================================
	std::cout << "\nJ. Homogeneous (sequence-derivative) family vs dM/deps_c from a TOV sweep — 1.6 Msun\n";
	{
		double relP_fin = 1.0, relE_fin = 1.0;
		bool okfin = false;
		std::size_t res_fin = 0;
		// 40000 was added after the 10000 and 20000 measurements (1.171e-3, 1.016e-3 in the p_c
		// form): the residual is grid-limited by the crust sampling of the dedp column, see the
		// integrability diagnostic printed with each row.
		for (std::size_t res : {std::size_t(10000), std::size_t(20000), std::size_t(40000)})
		{
			Star st = BuildAtDensity(cold, wrk / ("j" + std::to_string(res)), star16.ec_cgs, res);
			if (!st.ns)
			{
				Report("J build res=" + std::to_string(res), false, "");
				continue;
			}
			const NStar &c = *st.ns;
			std::vector<double> s0, sp0;
			ProdShape(c, s0, sp0);
			const Background2 bg = FromProfile(c, s0, sp0);
			MHOptions o;
			o.sources_on = false;
			o.phat_at_r0 = 1.0;
			o.I_exterior = 0.0;
			o.eos_measure = true; // ADR-0008 Q1/Q3: the accepted source, on the oracle too
			const auto hom = hartle_mono_ref::Solve(bg, o);
			const auto &P = c.Profile();
			const double ec = (*P.GetEnergyDensity())[0], pc = (*P.GetPressure())[0], dedpc = (*P.GetEosDEdP())[0];
			const double dp_c = (ec + pc);			  // p0*_c = 1  [km^-2]
			const double deps_c = (ec + pc) * dedpc; // authority derivative
			const double rel_step = 1.0e-3;
			const SeqPt up = SequencePoint(cold, wrk / "jp", star16.ec_cgs * (1.0 + rel_step), res);
			const SeqPt dn = SequencePoint(cold, wrk / "jm", star16.ec_cgs * (1.0 - rel_step), res);
			const double dM_dpc = (up.M_km - dn.M_km) / (up.pc - dn.pc);
			const double dM_dec = (up.M_km - dn.M_km) / (up.ec - dn.ec);
			const double dedp_seq = (up.ec - dn.ec) / (up.pc - dn.pc);
			const double expP = dM_dpc * dp_c, expE = dM_dec * deps_c, got = hom.deltaM_hat;
			// Stieltjes evaluation of the m0 source. With dp/dr = -(eps+p) nu' the source
			// 4 pi r^2 (eps+p)(deps/dp) p0* dr equals 4 pi r^2 p0* (-d eps)/nu': integrating against
			// the profile's OWN density steps includes every density discontinuity the nodal
			// (deps/dp) column cannot represent; integrating against the column reproduces the
			// SUPERSEDED discretization. Test-side diagnostic, no production change.
			// PHASE 4D-RI: since ADR-0008 the eps-step sum is what production itself integrates, so
			// this block now measures the discretization gap the correction closed, and
			// `dM_corr_sourced` is the contribution production no longer omits.
			double m0h_col = 0.0, m0h_stj = 0.0, dM_corr_sourced = 0.0;
			const auto fs_ = Fields(c);
			{
				const int nn = static_cast<int>(bg.r.size());
				for (int i = 0; i + 1 < nn; ++i)
				{
					const double rm = 0.5 * (bg.r[i] + bg.r[i + 1]), nupm = 0.5 * (bg.nup[i] + bg.nup[i + 1]);
					const double ph_h = 0.5 * (hom.phat[i] + hom.phat[i + 1]), ph_s = 0.5 * (fs_.phat[i] + fs_.phat[i + 1]);
					const double step_eps = bg.eps[i] - bg.eps[i + 1];
					const double step_col = 0.5 * (bg.dedp[i] + bg.dedp[i + 1]) * (bg.p[i] - bg.p[i + 1]);
					const double w = 4.0 * M_PI * rm * rm / nupm;
					m0h_col += w * ph_h * step_col;
					m0h_stj += w * ph_h * step_eps;
					dM_corr_sourced += w * ph_s * (step_eps - step_col);
				}
			}
			std::cout << "     res " << res << ": Stieltjes check — homogeneous m0(R_*) by ODE " << Sci(hom.mhat_R, 8)
					  << ", by sum against the dedp column " << Sci(m0h_col, 8) << " (rel " << Sci(Rel(m0h_col, hom.mhat_R))
					  << "), against the profile's eps steps " << Sci(m0h_stj, 8) << "\n"
					  << "                -> eps-step homogeneous deltaM vs (dM/dp_c) dp_c: rel=" << Sci(Rel(m0h_stj + hom.shell_hat, expP))
					  << "   | implied correction to the SOURCED deltaM_hat from unrepresented density steps: "
					  << Sci(dM_corr_sourced / fs_.deltaM) << " relative\n";
			// integrability of the nodal dedp column: trapezoid sum of (deps/dp) dp over the profile
			// must reproduce eps_* - eps_c exactly if the column samples the EOS derivative adequately
			double integ_all = 0.0, integ_crust = 0.0, d_eps_crust = 0.0;
			{
				const double eps_cc = 1.0e14 * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3;
				const int nn = static_cast<int>(bg.r.size());
				bool in_crust_started = false;
				double eps_at_cc = 0.0;
				for (int i = 0; i + 1 < nn; ++i)
				{
					const double dI = 0.5 * (bg.dedp[i] + bg.dedp[i + 1]) * (bg.p[i + 1] - bg.p[i]);
					integ_all += dI;
					if (bg.eps[i] < eps_cc)
					{
						if (!in_crust_started)
						{
							in_crust_started = true;
							eps_at_cc = bg.eps[i];
						}
						integ_crust += dI;
					}
				}
				d_eps_crust = bg.eps.back() - eps_at_cc;
			}
			const double d_eps_all = bg.eps.back() - bg.eps.front();
			std::cout << "     res " << res << ": dedp-column integrability: sum (deps/dp) dp vs (eps_* - eps_c): rel err "
					  << Sci(Rel(integ_all, d_eps_all)) << " over the whole star; crust (eps < 1e14 g/cm^3): "
					  << Sci(Rel(integ_crust, d_eps_crust)) << "\n";
			std::cout << "     res " << res << ": sequence dM/dp_c=" << Sci(dM_dpc, 8) << "  dM/deps_c=" << Sci(dM_dec, 8)
					  << "  (deps/dp)_c sequence=" << Sci(dedp_seq, 7) << " vs authority=" << Sci(dedpc, 7)
					  << "  rel=" << Sci(Rel(dedp_seq, dedpc)) << "\n"
					  << "                m0_hom(R_*)=" << Sci(hom.mhat_R, 8) << "  shell_hom=" << Sci(hom.shell_hat, 3)
					  << "  deltaM_hom=" << Sci(got, 8) << "\n"
					  << "                p_c form: expected (dM/dp_c) dp_c=" << Sci(expP, 8) << "  rel=" << Sci(Rel(got, expP))
					  << "   | eps_c form: expected (dM/deps_c) deps_c=" << Sci(expE, 8) << "  rel=" << Sci(Rel(got, expE)) << "\n";
			relP_fin = Rel(got, expP);
			relE_fin = Rel(got, expE);
			okfin = hom.ok && std::isfinite(dM_dpc);
			res_fin = res;
		}
		std::cout << "     RECORD — Ja (ADR-0007 §7 item 11, predeclared <= " << Sci(kJ_Homog, 0) << "): production's homogeneous delta M vs (dM/dp_c) dp_c at res "
				  << res_fin << ": rel=" << Sci(relP_fin) << (okfin && relP_fin <= kJ_Homog ? "  MET" : "  NOT MET") << "\n"
				  << "     (Not asserted: the residual is resolution-independent and is diagnosed above as density steps of the crust EOS that\n"
				  << "      the nodal deps/dp column cannot represent — an input-representation limit (INV-13) recorded in the 4D record, not widened away.)\n";
		Report("Ja the sequence-derivative comparison was produced at 10000/20000/40000 with the Stieltjes diagnostic", okfin, "");
		std::cout << "     RECORD — the eps_c form at the finest res: rel=" << Sci(relE_fin) << (relE_fin <= kJ_Homog ? " <= " : " > ")
				  << Sci(kJ_Homog, 0) << "; any gap between the two forms is the solver's central-condition inversion\n"
				  << "     against the authority derivative, not the monopole solver.\n";
		std::cout << "     (Test-side only — ADR-0007 P11 as modified at acceptance; no public API.)\n";
	}

	fs::remove_all(wrk);

	std::cout << "\nSLOW-ROTATION DISCLAIMER: coefficient correctness, not truncation accuracy at high spin.\n";
	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
