// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_physics_cmf.cpp
 * @brief Phase 4D-RV — corrected INDEPENDENT physical revalidation of the governed O(Omega^2)
 *        monopole response on the migrated (ADR-0009) DS(CMF)-1 backgrounds. External data.
 *
 * THE CLAIM UNDER TEST (ADR-0007 as amended by ADR-0008, on the ADR-0009 surface event). For
 * ordinary NStar, with R_* the accepted p = p_cut event stored as the final node, the fixed-eps_c
 * l = 0 response m0_hat, p0*_hat, delta_p0_hat, xi0_hat, deltaM_hat per Omega_geom^2 is
 * physically correct when the EOS energy-density source is the MEASURE -4 pi r^2 xi0_hat d(eps)
 * (smooth variation, sharp continuous tabulated variation, the terminal eps_* -> 0 atom). NOT
 * claimed: l = 2, high-spin truncation accuracy, baryon-conserving reduction, Phase-5 A_i/B_i/Z_i,
 * rotochemical heating, MixedStar rotation.
 *
 * THE INDEPENDENT ORACLE (hartle_monopole_reference.hpp). Hartle's (m0,h0) pair, (97)+(98) with
 * p0* an algebraic by-product of the first integral (90) — never production's (m0,p0*) system —
 * with its own interpolation, centre start, exterior arithmetic and tighter tolerances. Since
 * Phase 4D-RV the EOS measure is represented INDEPENDENTLY of production's per-segment
 * eps_slope density: `SolveStieltjes` accumulates the measure as midpoint ATOMS on an
 * independently constructed partition (K-fold refinement of every profile interval, optionally
 * with every EOS-table knot mapped into the interval first), while its driver integrates the
 * smooth part with no EOS term. The production-like secant realisation (`eos_measure`) and the
 * SUPERSEDED differential form are still run and REPORTED, never asserted.
 *
 * Two chains everywhere: SECOND-ORDER-ISOLATED (production's Phase-4B-verified s, s') and
 * FULLY INDEPENDENT (hartle_reference.hpp's own first order on the same tabulated background).
 *
 * ============================ PREDECLARED BOUNDS ============================
 *   Z  ADR-0009 background: every star SURFACE_REACHED, last node p == p_cut exactly, finite
 *      surface EOS data, no partial profile                                     exact
 *   G  full profile, production vs the INDEPENDENT Stieltjes oracle (profile partition K=4,
 *      and EOS-knot partition K=2), per star, both chains                       rel <= 1e-4   (ADR-0007 §7-4)
 *   F  near-vacuum identity on the outermost nodes (eps < 1e9 g/cm^3, r > R_* / 2), matter-
 *      corrected mhat + S + I^2/r^3 constant                                  <= 1e-6   (§7-5)
 *   C  ADR-0008 Validation C: sourced same-partition accounting vs production m0_hat  <= 1e-6
 *   R  reference floor: the EOS-knot oracle's own floor (tolerance, centre start, K-refinement)
 *      must be subdominant to its production disagreement                        ratio <= 0.1
 *      The profile-partition oracle converges in K to the SAME continuous-density measure
 *      production integrates, so its production disagreement sits AT its own floor; that
 *      agreement is reported as floor-limited and both must lie two decades under G  <= 1e-6
 *   H  ADR-0008 Validation D: deltaM_hat over 5000/10000/20000/40000 at fixed eps_c —
 *      relative spread <= 1e-4 (asserted, Hb) AND successive differences of decreasing
 *      magnitude |d1| >= |d2| >= |d3| (Hc). Hc is RECORDED, never waived: if it is not met the
 *      scientific status of the record is CHARACTERIZED, not VERIFIED (task §27), and the
 *      per-resolution decomposition (first-order I, EOS channel, rotational channel, both
 *      independent oracles) is printed so the cause is attributable. 80000 is a diagnostic
 *      extension; R_* spread <= 1e-8 (ADR-0009 floor) is asserted (Ha)
 *   I  ADR-0008 Validation E: the retired-FD substitution moves deltaM_hat by      < 1e-3
 *   J  ADR-0008 Validation B: homogeneous deltaM vs (dM/dp_c) dp_c, res 10000 and 20000,
 *      on complete ADR-0009 stars, INDEPENDENT Stieltjes oracle                <= 2e-4;  40000 reported
 * ===========================================================================
 * Chronology kept from the historical harness (recorded, not hidden): Phase 4D asserted §7-11 at
 * 1e-3 and measured 1.04e-3 (the nodal-derivative defect); Phase 4D-RI reported line J and H with
 * the corrected source on both sides and left their adjudication to this revalidation. The
 * near-vacuum window was corrected once in 4D ((eps+p)r^2 < 1e-6 admitted the central nodes).
 * No 4C-I1, 4D or 4D-RI diagnostic value is an expected answer anywhere in this file.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "eos_table_knots.hpp"
#include "hartle_monopole_reference.hpp"
#include "hartle_profile_compare.hpp"
#include "hartle_reference.hpp"

#include <Zaki/Physics/Constants.hpp>

#include <unistd.h>

namespace fs = std::filesystem;

/// HARTLE_4D_QUICK=1 restricts the run to Z + G (isolated chain) + F + C: used ONLY by the
/// detector sweep. The registered CTest runs the full set.
static bool Quick() { return std::getenv("HARTLE_4D_QUICK") != nullptr; }

using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolveStatus;
using CompactStar::Core::TOVSolver;
using hartle_mono_ref::Background2;
using hartle_mono_ref::Cmp;
using hartle_mono_ref::Compare;
using hartle_mono_ref::MHOptions;
using hartle_mono_ref::MHResult;
using hartle_mono_ref::StieltjesOptions;

static constexpr double kG_Profile = 1.0e-4;
static constexpr double kF_Exterior = 1.0e-6;
static constexpr double kC_Accounting = 1.0e-6;
static constexpr double kR_FloorRatio = 0.1;
static constexpr double kH_Spread = 1.0e-4;
static constexpr double kH_Rstar = 1.0e-8;
static constexpr double kI_FD = 1.0e-3;
static constexpr double kJ_Homog = 2.0e-4;
static constexpr double kVacuumEps_gcm3 = 1.0e9;
static constexpr int kK_Profile = 4; // Stieltjes refinement, profile partition
static constexpr int kK_Knots = 2;   // Stieltjes refinement, EOS-knot partition

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
//  ADR-0009-aware star construction. `Probe` only exposes the protected cutoff accessor.
// ---------------------------------------------------------------------------
struct Probe : TOVSolver
{
	using TOVSolver::PressureCutoff;
};

struct Star
{
	std::unique_ptr<NStar> ns;
	std::vector<TOVPoint> pts;
	std::vector<std::string> labels;
	double ec_cgs = 0, p_cut = 0;
	bool complete = false; ///< SURFACE_REACHED, last node exactly at the cutoff, finite surface data
	std::string why;
};

static void CheckComplete(Star &st, const Probe &tov)
{
	st.p_cut = tov.PressureCutoff();
	const bool status = tov.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED;
	const bool n_ok = st.pts.size() >= 4;
	const bool at_cut = n_ok && st.pts.back().p == st.p_cut;
	const bool finite = n_ok && std::isfinite(st.pts.back().e) && st.pts.back().e > 0.0 &&
						std::isfinite(st.pts.back().dedp) && std::isfinite(st.pts.back().m);
	bool increasing = n_ok;
	for (std::size_t i = 1; n_ok && i < st.pts.size(); ++i)
		increasing = increasing && st.pts[i].r > st.pts[i - 1].r && st.pts[i].p >= st.p_cut;
	st.complete = status && n_ok && at_cut && finite && increasing;
	if (!st.complete)
		st.why = std::string("status=") + (status ? "SURFACE_REACHED" : "other") + " nodes=" +
				 std::to_string(st.pts.size()) + " at_cut=" + (at_cut ? "1" : "0") + " finite=" +
				 (finite ? "1" : "0") + " increasing=" + (increasing ? "1" : "0");
}

static Star BuildAtMass(const fs::path &cold, const fs::path &wrk, double M)
{
	Star st;
	Probe tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	const int n = tov.SolveToProfile(M, st.pts, &st.labels);
	if (n <= 0 || st.pts.size() < 4)
	{
		st.why = "SolveToProfile returned no complete profile";
		return st;
	}
	CheckComplete(st, tov);
	st.ec_cgs = st.pts.front().e;
	st.ns = st.labels.empty() ? std::make_unique<NStar>(st.pts) : std::make_unique<NStar>(st.pts, st.labels);
	if (!st.ns->ComputeHartleMonopoleResponse())
		st.ns.reset();
	return st;
}

static Star BuildAtDensity(const fs::path &cold, const fs::path &wrk, double ec, std::size_t res)
{
	Star st;
	Probe tov;
	tov.SetWrkDir(wrk.string());
	tov.SetRadialRes(res);
	tov.ImportEOS(cold.string(), true);
	if (tov.SingleStarSolveToTOVPoints(ec, st.pts) <= 0 || st.pts.size() < 4)
	{
		st.why = "primitive returned no complete profile";
		return st;
	}
	CheckComplete(st, tov);
	st.ec_cgs = ec;
	st.ns = std::make_unique<NStar>(st.pts);
	if (!st.ns->ComputeHartleMonopoleResponse())
		st.ns.reset();
	return st;
}

/// One member of the non-rotating sequence: M [km] and the ACTUAL central p, eps [km^-2] the
/// solver used, for the sequence derivative. Only a COMPLETE (SURFACE_REACHED, p = p_cut) star
/// may enter the derivative (ADR-0009 Q8; task §21).
struct SeqPt
{
	double M_km = std::numeric_limits<double>::quiet_NaN(), pc = 0, ec = 0;
	bool complete = false;
};
static SeqPt SequencePoint(const fs::path &cold, const fs::path &wrk, double ec, std::size_t res)
{
	SeqPt o;
	Probe tov;
	tov.SetWrkDir(wrk.string());
	tov.SetRadialRes(res);
	tov.ImportEOS(cold.string(), true);
	std::vector<TOVPoint> pts;
	if (tov.SingleStarSolveToTOVPoints(ec, pts) <= 0 || pts.empty())
		return o;
	o.complete = tov.LastSolveStatus() == TOVSolveStatus::SURFACE_REACHED && pts.back().p == tov.PressureCutoff();
	if (!o.complete)
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

/// Linear interpolation of a field on the profile grid at a requested radius.
static double AtRadius(const std::vector<double> &r, const std::vector<double> &f, double rq)
{
	return Background2::Lerp(r, f, rq);
}

// ---------------------------------------------------------------------------
//  ADR-0008 Validation C — same-partition source accounting (own driver, rtol 1e-13). Carries
//  m0_hat, p0*_hat, the EOS channel of m0_hat and the per-segment weight, one governed segment
//  per driver call. An implementation identity, distinct from the independent oracle above.
// ---------------------------------------------------------------------------
struct AccParams
{
	const Background2 *bg;
	double slope;
};
static int AccRHS(double r, const double y[], double f[], void *pv)
{
	const AccParams *P = static_cast<const AccParams *>(pv);
	const auto b = P->bg->At(r);
	const double D = 1.0 - 2.0 * b.m / r, r_2m = r * D, e2 = std::exp(-2.0 * b.nu), ep = b.eps + b.p;
	const double r2 = r * r, r3 = r2 * r, r4 = r3 * r;
	const double mhat = y[0], phat = y[1];
	const double xi = phat / b.nup;
	if (!std::isfinite(xi) || !(D > 0.0))
		return GSL_EBADFUNC;
	const double term1 = (P->slope != 0.0) ? -4.0 * M_PI * r2 * xi * P->slope : 0.0;
	f[0] = term1 + (1.0 / 12.0) * r4 * e2 * D * b.sp * b.sp + (8.0 * M_PI / 3.0) * r4 * ep * e2 * b.s * b.s;
	f[1] = -mhat * (1.0 + 8.0 * M_PI * r2 * b.p) / (r_2m * r_2m) - 4.0 * M_PI * ep * r2 * phat / r_2m +
		   (1.0 / 12.0) * r3 * e2 * b.sp * b.sp + (2.0 / 3.0) * r * e2 * b.s * (b.s + r * b.sp - r * b.nup * b.s);
	f[2] = term1;
	f[3] = 4.0 * M_PI * r2 * xi;
	return (std::isfinite(f[0]) && std::isfinite(f[1])) ? GSL_SUCCESS : GSL_EBADFUNC;
}
struct AccResult
{
	bool ok = false;
	std::vector<double> mhat;
	double eos_total = 0.0, worst_segment_abs = 0.0;
};
static AccResult Account(const Background2 &b, double mhat0, double phat0)
{
	AccResult out;
	const std::size_t n = b.N();
	AccParams P{&b, 0.0};
	gsl_odeiv2_system sys = {AccRHS, nullptr, 4, &P};
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-16, 1e-13);
	double y[4] = {mhat0, phat0, 0.0, 0.0};
	double r = b.r.front();
	out.mhat.assign(n, 0.0);
	out.mhat[0] = y[0];
	for (std::size_t i = 1; i < n; ++i)
	{
		P.slope = (b.eps[i] - b.eps[i - 1]) / (b.r[i] - b.r[i - 1]);
		const double before = y[2];
		y[3] = 0.0;
		if (gsl_odeiv2_driver_apply(d, &r, b.r[i], y) != GSL_SUCCESS)
		{
			gsl_odeiv2_driver_free(d);
			return out;
		}
		out.worst_segment_abs = std::max(out.worst_segment_abs, std::fabs((y[2] - before) + P.slope * y[3]));
		out.mhat[i] = y[0];
	}
	gsl_odeiv2_driver_free(d);
	out.eos_total = y[2];
	out.ok = true;
	return out;
}

// ---------------------------------------------------------------------------
struct OracleSet
{
	MHResult secant, differential, stj_profile, stj_knots;
};
static OracleSet RunOracles(const Background2 &bg, double I_ext, const eos_knots::Knots *knots,
							bool sources_on = true, double phat_at_r0 = std::numeric_limits<double>::quiet_NaN())
{
	OracleSet o;
	MHOptions base;
	base.I_exterior = I_ext;
	base.sources_on = sources_on;
	base.phat_at_r0 = phat_at_r0;
	MHOptions sec = base;
	sec.eos_measure = true;
	o.secant = hartle_mono_ref::Solve(bg, sec);
	MHOptions dif = base;
	dif.eos_measure = false;
	o.differential = hartle_mono_ref::Solve(bg, dif);
	StieltjesOptions sp;
	sp.refine = kK_Profile;
	o.stj_profile = hartle_mono_ref::SolveStieltjes(bg, base, sp);
	if (knots != nullptr && knots->ok)
	{
		StieltjesOptions sk;
		sk.refine = kK_Knots;
		sk.knot_p = &knots->p;
		sk.knot_eps = &knots->eps;
		o.stj_knots = hartle_mono_ref::SolveStieltjes(bg, base, sk);
	}
	return o;
}

struct ChainCmp
{
	Cmp m, p, d, x;
	double dM = 0;
	bool ok = false;
	double Worst() const { return std::max({m.max_rel, p.max_rel, x.max_rel, dM}); }
};
static ChainCmp CompareFields(const ProdFields &pf, const MHResult &ref)
{
	ChainCmp c;
	c.ok = ref.ok;
	c.m = Compare(pf.mhat, ref.mhat, pf.r);
	c.p = Compare(pf.phat, ref.phat, pf.r);
	c.d = Compare(pf.dphat, ref.dphat, pf.r);
	c.x = Compare(pf.xihat, ref.xihat, pf.r);
	c.dM = Rel(pf.deltaM, ref.deltaM_hat);
	return c;
}

struct StarRow
{
	double target = 0, ec = 0, M = 0, R = 0, p_cut = 0, eps_R = 0, I = 0, dedp_c = 0;
	double mhat_R = 0, phat_R = 0, xi_R = 0, shell = 0, ext = 0, dM = 0;
	ChainCmp iso_prof, iso_knot, iso_sec, full_prof, full_knot;
	double old_form = 0, prof_vs_knot_dM = 0, prof_vs_knot_node = 0;
	double ext_const = 0, raw_spread = 0;
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
	const fs::path wrk = fs::temp_directory_path() / ("compactstar_hartle_monopole_4drv_cmf_" + std::to_string(::getpid()));
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4D-RV — corrected INDEPENDENT validation of the governed O(Omega^2) monopole response on\n"
				 "DS(CMF)-1_with_crust (ADR-0007 + ADR-0008 on the ADR-0009 surface event). Oracle: the test-only\n"
				 "(m0,h0) solver with the measure accumulated as Stieltjes atoms on an independent partition.\n\n";

	const eos_knots::Knots knots = eos_knots::Read(cold.string());
	Report("K0 EOS-table knots read for the knot-refined measure partition", knots.ok,
		   std::to_string(knots.rows) + " rows");
	const double km2_to_gcm3 = Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;

	// =====================================================================
	//  Z + G + F + C — four stars
	// =====================================================================
	std::cout << "\nG. Four canonical stars (governed target-mass workflow), production vs the independent oracle\n";
	const double masses[] = {1.0, 1.4, 1.6, 2.0};
	std::vector<StarRow> rows;
	Star star16;

	for (double M : masses)
	{
		Star st = BuildAtMass(cold, wrk / ("m" + std::to_string(M)), M);
		Report("Z " + Sci(M, 2) + " Msun: complete ADR-0009 star (SURFACE_REACHED, last node p == p_cut, finite surface data)",
			   st.complete, st.why);
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
		row.target = M;
		row.ec = st.ec_cgs;
		row.M = st.pts.back().m;
		row.R = pf.R;
		row.p_cut = st.p_cut;
		row.eps_R = pf.eps_R * km2_to_gcm3;
		row.I = pf.I;
		row.n = pf.r.size();
		row.dedp_c = (*cns.Profile().GetEosDEdP())[0];
		row.mhat_R = pf.mhat.back();
		row.phat_R = pf.phat.back();
		row.xi_R = pf.xi_R;
		row.shell = pf.shell;
		row.ext = pf.I * pf.I / (pf.R * pf.R * pf.R);
		row.dM = pf.deltaM;

		// ---- chain A: production s, s' ----
		const Background2 bg_iso = FromProfile(cns, s_p, sp_p);
		const OracleSet oa = RunOracles(bg_iso, pf.I, &knots);
		row.iso_prof = CompareFields(pf, oa.stj_profile);
		row.iso_knot = CompareFields(pf, oa.stj_knots);
		row.iso_sec = CompareFields(pf, oa.secant);
		row.old_form = oa.differential.ok ? Rel(oa.differential.deltaM_hat, pf.deltaM) : std::numeric_limits<double>::quiet_NaN();
		row.prof_vs_knot_dM = Rel(oa.stj_profile.deltaM_hat, oa.stj_knots.deltaM_hat);
		row.prof_vs_knot_node = Compare(oa.stj_profile.mhat, oa.stj_knots.mhat, pf.r).max_rel;

		// ---- chain B: independent first order ----
		hartle_ref::Result fo;
		if (!Quick())
		{
			const auto hb = hartle_4b::BackgroundFromProfile(cns);
			fo = hartle_ref::Solve(hb, 5.0e-3, hb.r.front(), 1e-13, 1e-16);
			Background2 bg_full = bg_iso;
			bg_full.s = fo.s;
			bg_full.sp = fo.s_prime;
			const OracleSet ob = RunOracles(bg_full, fo.I_surface, &knots);
			row.full_prof = CompareFields(pf, ob.stj_profile);
			row.full_knot = CompareFields(pf, ob.stj_knots);
			row.full_prof.ok = row.full_prof.ok && fo.ok;
			row.full_knot.ok = row.full_knot.ok && fo.ok;
		}

		// ---- F: near-vacuum identity on the outermost nodes ----
		double sp_vs_I = std::numeric_limits<double>::quiet_NaN();
		{
			const auto &P = cns.Profile();
			const double eps_win = kVacuumEps_gcm3 / km2_to_gcm3;
			const std::size_t n = pf.r.size();
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
			row.raw_spread = nv >= 2 ? (rhi - rlo) / pf.deltaM : std::numeric_limits<double>::quiet_NaN();
			row.ext_const = nv >= 2 ? (chi - clo) / pf.deltaM : std::numeric_limits<double>::quiet_NaN();
			sp_vs_I = Rel(sp_p[n - 1], 6.0 * pf.I / (pf.R * pf.R * pf.R * pf.R));
		}

		// ---- C: same-partition accounting (ADR-0008 Validation C) ----
		const auto acc = Account(bg_iso, pf.mhat.front(), pf.phat.front());
		const double acc_node = acc.ok ? Compare(acc.mhat, pf.mhat, pf.r).max_rel : 1.0;
		const double acc_R = acc.ok ? Rel(acc.mhat.back(), pf.mhat.back()) : 1.0;

		const bool ok_iso = row.iso_prof.ok && row.iso_knot.ok && row.iso_prof.Worst() <= kG_Profile &&
							row.iso_knot.Worst() <= kG_Profile && row.iso_prof.d.max_rel <= kG_Profile;
		Report(Sci(M, 2) + " Msun: SECOND-ORDER-ISOLATED chain — production agrees with the INDEPENDENT Stieltjes oracle (profile K=4 and EOS-knot K=2) at every node",
			   ok_iso,
			   "profile: mhat=" + Sci(row.iso_prof.m.max_rel) + " phat=" + Sci(row.iso_prof.p.max_rel) + " dphat=" +
				   Sci(row.iso_prof.d.max_rel) + " xihat=" + Sci(row.iso_prof.x.max_rel) + " dM=" + Sci(row.iso_prof.dM) +
				   " | knots: mhat=" + Sci(row.iso_knot.m.max_rel) + " phat=" + Sci(row.iso_knot.p.max_rel) + " xihat=" +
				   Sci(row.iso_knot.x.max_rel) + " dM=" + Sci(row.iso_knot.dM) + "  bound=" + Sci(kG_Profile));
		if (!Quick())
			Report(Sci(M, 2) + " Msun: FULLY INDEPENDENT chain agrees at every node (profile and knot oracles)",
				   row.full_prof.ok && row.full_knot.ok && row.full_prof.Worst() <= kG_Profile && row.full_knot.Worst() <= kG_Profile,
				   "profile: mhat=" + Sci(row.full_prof.m.max_rel) + " phat=" + Sci(row.full_prof.p.max_rel) + " xihat=" +
					   Sci(row.full_prof.x.max_rel) + " dM=" + Sci(row.full_prof.dM) + " | knots: mhat=" + Sci(row.full_knot.m.max_rel) +
					   " phat=" + Sci(row.full_knot.p.max_rel) + " xihat=" + Sci(row.full_knot.x.max_rel) + " dM=" + Sci(row.full_knot.dM));
		Report(Sci(M, 2) + " Msun: near-vacuum identity — matter-corrected mhat + S + I^2/r^3 constant over the outermost nodes",
			   std::isfinite(row.ext_const) && row.ext_const <= kF_Exterior,
			   "corrected spread/deltaM=" + Sci(row.ext_const) + "  raw spread/deltaM=" + Sci(row.raw_spread) + " over " +
				   std::to_string(row.n_vac) + " nodes (eps < 1e9 g/cm^3)  bound=" + Sci(kF_Exterior) +
				   "  | s'(R_*) vs 6I/R_*^4 rel=" + Sci(sp_vs_I));
		Report(Sci(M, 2) + " Msun: matching arithmetic deltaM_hat = mhat(R_*) + shell + I^2/R_*^3 is exact",
			   Rel(pf.deltaM, pf.mhat.back() + pf.shell + pf.I * pf.I / (pf.R * pf.R * pf.R)) <= 1e-14, "");
		Report(Sci(M, 2) + " Msun: ADR-0008 Validation C — same-partition source accounting reproduces production m0_hat",
			   acc.ok && acc_node <= kC_Accounting && acc_R <= kC_Accounting,
			   "worst node rel=" + Sci(acc_node) + "  m0_hat(R_*) rel=" + Sci(acc_R) + "  worst segment residual=" +
				   Sci(acc.worst_segment_abs) + " km^3 of an EOS total " + Sci(acc.eos_total) + "  bound=" + Sci(kC_Accounting));
		std::cout << "     REPORTED — production-like secant oracle (variable-pair cross-check only): mhat=" << Sci(row.iso_sec.m.max_rel)
				  << " phat=" << Sci(row.iso_sec.p.max_rel) << " dM=" << Sci(row.iso_sec.dM) << "\n"
				  << "     REPORTED — measure-discretization diagnostic, profile-partition vs EOS-knot Stieltjes oracle: deltaM_hat rel="
				  << Sci(row.prof_vs_knot_dM) << ", node-wise mhat rel=" << Sci(row.prof_vs_knot_node) << "\n"
				  << "     REPORTED — the SUPERSEDED differential (nodal deps/dp) oracle disagrees with production deltaM_hat by "
				  << Sci(row.old_form) << " (the sub-node energy-density variation ADR-0008 recovers; M10's target)\n";

		rows.push_back(row);
		if (M == 1.6)
			star16 = std::move(st);
	}

	std::cout << "\n  target  eps_c[g/cm^3]    M[Msun]         R_*[km]      p_cut[dyn/cm^2]  eps_*[g/cm^3]  I[km^3]       dedp_c\n";
	for (const auto &w : rows)
	{
		char b[400];
		snprintf(b, sizeof(b), "  %5.2f   %.9e  %.10f  %.9f  %.7e  %.7e  %.9e  %.7e\n", w.target, w.ec, w.M, w.R, w.p_cut, w.eps_R, w.I, w.dedp_c);
		std::cout << b;
	}
	std::cout << "  target  mhat(R_*)        phat(R_*)        xi_hat(R_*)      shell_hat        I^2/R_*^3        deltaM_hat        nodes\n";
	for (const auto &w : rows)
	{
		char b[400];
		snprintf(b, sizeof(b), "  %5.2f   %.9e  %.9e  %.9e  %.9e  %.9e  %.10e  %zu\n", w.target, w.mhat_R, w.phat_R, w.xi_R, w.shell, w.ext, w.dM, w.n);
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
		std::cout << "\nFAIL — 1.6 Msun star unavailable for R/H/I/J\n";
		return 1;
	}
	const NStar &c16 = *star16.ns;
	const auto pf16 = Fields(c16);
	std::vector<double> s16, sp16;
	ProdShape(c16, s16, sp16);
	const Background2 bg16 = FromProfile(c16, s16, sp16);

	// =====================================================================
	//  R — the reference's own floor on 1.6 Msun
	// =====================================================================
	std::cout << "\nR. Reference floor on 1.6 Msun — integrator tolerance, centre start, measure-partition refinement, EOS-knot refinement\n";
	{
		MHOptions o0;
		o0.I_exterior = pf16.I;
		StieltjesOptions p4;
		p4.refine = kK_Profile;
		const auto base = hartle_mono_ref::SolveStieltjes(bg16, o0, p4);
		StieltjesOptions k2;
		k2.refine = kK_Knots;
		k2.knot_p = &knots.p;
		k2.knot_eps = &knots.eps;
		const auto kbase = hartle_mono_ref::SolveStieltjes(bg16, o0, k2);
		auto move = [&](const MHResult &v, const MHResult &ref) {
			return std::max({Compare(v.mhat, ref.mhat, pf16.r).max_rel, Compare(v.phat, ref.phat, pf16.r).max_rel,
							 Compare(v.xihat, ref.xihat, pf16.r).max_rel, Rel(v.deltaM_hat, ref.deltaM_hat)});
		};
		double floor_tol = 0.0, floor_r0 = 0.0, floor_K = 0.0, floor_knotK = 0.0;
		for (double tol : {1e-11, 1e-15})
		{
			MHOptions o = o0;
			o.rtol = tol;
			o.atol = tol * 1e-3;
			floor_tol = std::max(floor_tol, move(hartle_mono_ref::SolveStieltjes(bg16, o, p4), base));
		}
		{
			MHOptions o = o0;
			o.r0 = 1e-7;
			floor_r0 = std::max(floor_r0, move(hartle_mono_ref::SolveStieltjes(bg16, o, p4), base));
		}
		// The reference's own discretization floor is its movement under FURTHER refinement
		// (K=8 vs K=4; knot K=4 vs K=2); the coarser points form the convergence ladder and are
		// reported, not counted as the floor of the K=4 / K=2 references.
		std::cout << "     profile-partition refinement ladder K (deltaM_hat, and movement vs K=4):\n";
		for (int K : {1, 2, 8})
		{
			StieltjesOptions s;
			s.refine = K;
			const auto v = hartle_mono_ref::SolveStieltjes(bg16, o0, s);
			const double mv = move(v, base);
			if (K > kK_Profile)
				floor_K = std::max(floor_K, mv);
			std::cout << "       K=" << K << ": deltaM_hat=" << Sci(v.deltaM_hat, 10) << "  movement=" << Sci(mv) << (K > kK_Profile ? "  (floor)" : "  (ladder)") << "\n";
		}
		std::cout << "       K=4: deltaM_hat=" << Sci(base.deltaM_hat, 10) << " (reference)\n";
		std::cout << "     EOS-knot partition refinement ladder K (deltaM_hat, and movement vs knot K=2):\n";
		for (int K : {1, 4})
		{
			StieltjesOptions s = k2;
			s.refine = K;
			const auto v = hartle_mono_ref::SolveStieltjes(bg16, o0, s);
			const double mv = move(v, kbase);
			if (K > kK_Knots)
				floor_knotK = std::max(floor_knotK, mv);
			std::cout << "       K=" << K << ": deltaM_hat=" << Sci(v.deltaM_hat, 10) << "  movement=" << Sci(mv) << (K > kK_Knots ? "  (floor)" : "  (ladder)") << "\n";
		}
		std::cout << "       K=2: deltaM_hat=" << Sci(kbase.deltaM_hat, 10) << " (reference)\n";
		const auto it = std::find_if(rows.begin(), rows.end(), [](const StarRow &w) { return w.target == 1.6; });
		const double signal = it->iso_prof.Worst();
		const double signal_k = it->iso_knot.Worst();
		const double floor = std::max({floor_tol, floor_r0, floor_K});
		std::cout << "     floor components: tolerance " << Sci(floor_tol) << ", centre start " << Sci(floor_r0)
				  << ", partition refinement " << Sci(floor_K) << ", knot refinement " << Sci(floor_knotK) << "\n"
				  << "     absolute floor on deltaM_hat: " << Sci(floor * base.deltaM_hat) << " km^3; relative floor " << Sci(floor)
				  << "; production-oracle disagreement (profile oracle) " << Sci(signal) << ", (knot oracle) " << Sci(signal_k) << "\n";
		const double floor_k = std::max(floor_knotK, std::max(floor_tol, floor_r0));
		Report("Ra the knot-oracle floor (tolerance, centre start, K-refinement) is subdominant to its production disagreement",
			   signal_k > 0 && floor_k / signal_k <= kR_FloorRatio,
			   "floor=" + Sci(floor_k) + "  signal=" + Sci(signal_k) + "  floor/signal=" + Sci(signal_k > 0 ? floor_k / signal_k : 0.0, 2));
		Report("Rb the profile-partition oracle and production agree down to that oracle's own floor, both two decades under the G bound",
			   std::max(floor, signal) <= 1.0e-6,
			   "floor=" + Sci(floor) + "  disagreement=" + Sci(signal) + "  (floor-limited agreement at " + Sci(std::max(floor, signal)) + ")");
	}

	// =====================================================================
	//  H — corrected radial convergence at fixed eps_c (ADR-0008 Validation D)
	// =====================================================================
	std::cout << "\nH. Radial convergence at fixed eps_c = " << Sci(star16.ec_cgs, 9) << " g/cm^3 (1.6 Msun), migrated ADR-0009 background\n";
	{
		struct Hrow
		{
			std::size_t res, n;
			double R, M, I, mR, xi, shell, ext, dM, p25, p50, p75, dM_prof, dM_knot;
			double eos_part, rot_part; // m0_hat(R_*) channels from the same-partition accounting integrator
			bool complete;
		};
		std::vector<Hrow> h;
		for (std::size_t res : {std::size_t(5000), std::size_t(10000), std::size_t(20000), std::size_t(40000), std::size_t(80000)})
		{
			Star st = BuildAtDensity(cold, wrk / ("h" + std::to_string(res)), star16.ec_cgs, res);
			if (!st.ns)
			{
				Report("H build res=" + std::to_string(res), false, st.why);
				continue;
			}
			const auto f = Fields(*st.ns);
			std::vector<double> s, sp;
			ProdShape(*st.ns, s, sp);
			const Background2 bg = FromProfile(*st.ns, s, sp);
			const OracleSet o = RunOracles(bg, f.I, &knots);
			const auto acc = Account(bg, f.mhat.front(), f.phat.front());
			h.push_back({res, f.r.size(), f.R, st.pts.back().m, f.I, f.mhat.back(), f.xi_R, f.shell, f.I * f.I / (f.R * f.R * f.R), f.deltaM,
						 AtRadius(f.r, f.phat, 0.25 * f.R), AtRadius(f.r, f.phat, 0.5 * f.R), AtRadius(f.r, f.phat, 0.75 * f.R),
						 o.stj_profile.deltaM_hat, o.stj_knots.deltaM_hat,
						 acc.ok ? acc.eos_total : std::numeric_limits<double>::quiet_NaN(),
						 acc.ok ? acc.mhat.back() - acc.mhat.front() - acc.eos_total : std::numeric_limits<double>::quiet_NaN(), st.complete});
		}
		std::cout << "     res    nodes   complete  R_*[km]         M[Msun]          I[km^3]          mhat(R_*)        xi_hat(R_*)      shell_hat     I^2/R_*^3      deltaM_hat\n";
		for (const auto &x : h)
		{
			char b[400];
			snprintf(b, sizeof(b), "     %5zu  %5zu   %d   %.9f  %.10f  %.9e  %.9e  %.9e  %.4e  %.9e  %.10e\n", x.res, x.n, x.complete ? 1 : 0, x.R, x.M, x.I, x.mR, x.xi, x.shell, x.ext, x.dM);
			std::cout << b;
		}
		std::cout << "     res    phat(R/4)        phat(R/2)        phat(3R/4)       | independent oracle deltaM_hat: profile K=4     EOS-knot K=2\n";
		for (const auto &x : h)
		{
			char b[400];
			snprintf(b, sizeof(b), "     %5zu  %.9e  %.9e  %.9e  |  %.10e  %.10e\n", x.res, x.p25, x.p50, x.p75, x.dM_prof, x.dM_knot);
			std::cout << b;
		}
		std::cout << "     res    m0_hat(R_*) EOS channel   rotational channel   | first-order I[km^3]   I vs 80000\n";
		for (const auto &x : h)
		{
			char b[400];
			snprintf(b, sizeof(b), "     %5zu  %.10e       %.10e     |  %.10e   %+.3e\n", x.res, x.eos_part, x.rot_part, x.I, h.back().I != 0 ? x.I / h.back().I - 1.0 : 0.0);
			std::cout << b;
		}
		const bool have4 = h.size() >= 4 && h[0].res == 5000 && h[3].res == 40000;
		if (have4)
		{
			auto diffs = [&](auto get) {
				std::vector<double> d;
				for (std::size_t i = 1; i < 4; ++i)
					d.push_back(get(h[i]) - get(h[i - 1]));
				return d;
			};
			const auto dP = diffs([](const Hrow &x) { return x.dM; });
			const auto dO = diffs([](const Hrow &x) { return x.dM_prof; });
			const auto dK = diffs([](const Hrow &x) { return x.dM_knot; });
			auto spread4 = [&](auto get) {
				double lo = 1e300, hi = -1e300;
				for (std::size_t i = 0; i < 4; ++i)
				{
					lo = std::min(lo, get(h[i]));
					hi = std::max(hi, get(h[i]));
				}
				return (hi - lo) / std::fabs(get(h[1]));
			};
			auto decreasing = [](const std::vector<double> &d) {
				return std::fabs(d[0]) >= std::fabs(d[1]) && std::fabs(d[1]) >= std::fabs(d[2]);
			};
			auto same_sign = [](const std::vector<double> &d) {
				return (d[0] >= 0) == (d[1] >= 0) && (d[1] >= 0) == (d[2] >= 0);
			};
			const double sprP = spread4([](const Hrow &x) { return x.dM; });
			const double sprO = spread4([](const Hrow &x) { return x.dM_prof; });
			const double sprK = spread4([](const Hrow &x) { return x.dM_knot; });
			double rspread = 0.0;
			for (std::size_t i = 0; i < 4; ++i)
				rspread = std::max(rspread, Rel(h[i].R, h[1].R));
			std::cout << "     successive differences 5000->10000->20000->40000 [km^3]:\n"
					  << "       production            : " << Sci(dP[0]) << ", " << Sci(dP[1]) << ", " << Sci(dP[2])
					  << "   magnitudes decreasing: " << (decreasing(dP) ? "YES" : "NO") << "   monotone values: " << (same_sign(dP) ? "YES" : "NO") << "\n"
					  << "       oracle, profile K=4   : " << Sci(dO[0]) << ", " << Sci(dO[1]) << ", " << Sci(dO[2])
					  << "   magnitudes decreasing: " << (decreasing(dO) ? "YES" : "NO") << "   monotone values: " << (same_sign(dO) ? "YES" : "NO") << "\n"
					  << "       oracle, EOS-knot K=2  : " << Sci(dK[0]) << ", " << Sci(dK[1]) << ", " << Sci(dK[2])
					  << "   magnitudes decreasing: " << (decreasing(dK) ? "YES" : "NO") << "   monotone values: " << (same_sign(dK) ? "YES" : "NO") << "\n";
			if (h.size() >= 5)
				std::cout << "       extension 40000->80000: production " << Sci(h[4].dM - h[3].dM) << ", profile oracle " << Sci(h[4].dM_prof - h[3].dM_prof)
						  << ", knot oracle " << Sci(h[4].dM_knot - h[3].dM_knot) << " km^3 (diagnostic)\n";
			std::cout << "     relative spread over 5000..40000: production " << Sci(sprP) << ", profile oracle " << Sci(sprO) << ", knot oracle " << Sci(sprK)
					  << "; R_* spread " << Sci(rspread) << "\n";
			Report("Ha every convergence star is a complete ADR-0009 star and R_* is partition-invariant at the validated floor",
				   std::all_of(h.begin(), h.end(), [](const Hrow &x) { return x.complete; }) && rspread <= kH_Rstar,
				   "R_* spread=" + Sci(rspread) + "  bound=" + Sci(kH_Rstar));
			Report("Hb ADR-0008 Validation D (spread): production deltaM_hat relative spread over 5000/10000/20000/40000 <= 1e-4",
				   sprP <= kH_Spread, "spread=" + Sci(sprP) + "  bound=" + Sci(kH_Spread));
			auto diffsI = diffs([](const Hrow &x) { return x.I; });
			auto diffsE = diffs([](const Hrow &x) { return x.eos_part; });
			auto diffsR = diffs([](const Hrow &x) { return x.rot_part; });
			std::cout << "       first-order I         : " << Sci(diffsI[0]) << ", " << Sci(diffsI[1]) << ", " << Sci(diffsI[2])
					  << "   magnitudes decreasing: " << (decreasing(diffsI) ? "YES" : "NO") << "   monotone values: " << (same_sign(diffsI) ? "YES" : "NO") << "\n"
					  << "       m0 EOS channel        : " << Sci(diffsE[0]) << ", " << Sci(diffsE[1]) << ", " << Sci(diffsE[2])
					  << "   magnitudes decreasing: " << (decreasing(diffsE) ? "YES" : "NO") << "   monotone values: " << (same_sign(diffsE) ? "YES" : "NO") << "\n"
					  << "       m0 rotational channel : " << Sci(diffsR[0]) << ", " << Sci(diffsR[1]) << ", " << Sci(diffsR[2])
					  << "   magnitudes decreasing: " << (decreasing(diffsR) ? "YES" : "NO") << "   monotone values: " << (same_sign(diffsR) ? "YES" : "NO") << "\n";
			const bool hc = decreasing(dP);
			std::cout << "     RECORD — Hc ADR-0008 Validation D (monotonicity of production deltaM_hat, 5000->10000->20000->40000): "
					  << (hc ? "MET" : "NOT MET") << " — |d| = " << Sci(std::fabs(dP[0])) << " -> " << Sci(std::fabs(dP[1])) << " -> " << Sci(std::fabs(dP[2]))
					  << " km^3; reading (ii) 'values monotone': " << (same_sign(dP) ? "MET" : "NOT MET") << ".\n"
					  << "     This line is RECORDED, NOT WAIVED: when it is not met the scientific status of the revalidation record is\n"
					  << "     CHARACTERIZED — INDEPENDENT VALIDATION INCOMPLETE (task §27), no monopole baseline is created, and the\n"
					  << "     decomposition above attributes the residual (first-order sampled background vs measure weight location).\n";
		}
		else
			Report("Ha the 5000/10000/20000/40000 sequence was produced", false, "");
	}

	// =====================================================================
	//  I — EOS derivative / measure sensitivity on 1.6 Msun
	// =====================================================================
	std::cout << "\nI. EOS-derivative / measure sensitivity on 1.6 Msun\n";
	{
		std::vector<TOVPoint> pts_fd = star16.pts;
		const double c2 = Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_G_CM3;
		const std::size_t n = pts_fd.size();
		for (std::size_t i = 0; i < n; ++i)
		{
			const std::size_t lo = (i == 0) ? 0 : i - 1, hi = (i + 1 < n) ? i + 1 : n - 1;
			const double dp = pts_fd[hi].p - pts_fd[lo].p, de = pts_fd[hi].e - pts_fd[lo].e;
			pts_fd[i].dedp = std::fabs(dp) > 1e-30 ? (de / dp) * c2 : 1.0; // the retired candidate's fallback, deliberately
		}
		NStar ns_fd = star16.labels.empty() ? NStar(pts_fd) : NStar(pts_fd, star16.labels);
		const bool ok = ns_fd.ComputeHartleMonopoleResponse();
		const OracleSet o = RunOracles(bg16, pf16.I, &knots);
		std::cout << "     A governed measure (production):            deltaM_hat=" << Sci(pf16.deltaM, 10) << "  xi_hat(R_*)=" << Sci(pf16.xi_R, 10) << "\n";
		if (ok)
		{
			const auto ff = Fields(ns_fd);
			std::cout << "     B retired profile-FD derivative substituted:  deltaM_hat=" << Sci(ff.deltaM, 10) << "  spread vs A=" << Sci(Rel(ff.deltaM, pf16.deltaM)) << "\n";
			Report("Ia ADR-0008 Validation E — the retired-FD substitution moves deltaM_hat by < 1e-3 (deps/dp enters only the centre series now)",
				   Rel(ff.deltaM, pf16.deltaM) < kI_FD, "spread=" + Sci(Rel(ff.deltaM, pf16.deltaM)) + "  bound=" + Sci(kI_FD));
		}
		else
			Report("Ia FD-substituted star computed", false, "");
		std::cout << "     C nodal deps/dp differential form (oracle):    deltaM_hat=" << Sci(o.differential.deltaM_hat, 10) << "  vs A=" << Sci(Rel(o.differential.deltaM_hat, pf16.deltaM))
				  << "  (the superseded representation; the historical ~5 % deficit as a diagnostic)\n"
				  << "     D EOS-knot-refined measure (oracle):           deltaM_hat=" << Sci(o.stj_knots.deltaM_hat, 10) << "  vs A=" << Sci(Rel(o.stj_knots.deltaM_hat, pf16.deltaM)) << "\n"
				  << "       profile-partition Stieltjes (oracle):        deltaM_hat=" << Sci(o.stj_profile.deltaM_hat, 10) << "  vs A=" << Sci(Rel(o.stj_profile.deltaM_hat, pf16.deltaM)) << "\n"
				  << "       production-like secant (oracle):             deltaM_hat=" << Sci(o.secant.deltaM_hat, 10) << "  vs A=" << Sci(Rel(o.secant.deltaM_hat, pf16.deltaM)) << "\n";
	}

	// =====================================================================
	//  J — homogeneous family vs the TOV sequence derivative (ADR-0008 Validation B)
	// =====================================================================
	std::cout << "\nJ. Homogeneous (sequence-derivative) family vs dM/dp_c from a TOV sweep of complete ADR-0009 stars — 1.6 Msun\n";
	{
		bool ok_bind = true;
		for (std::size_t res : {std::size_t(10000), std::size_t(20000), std::size_t(40000)})
		{
			Star st = BuildAtDensity(cold, wrk / ("j" + std::to_string(res)), star16.ec_cgs, res);
			if (!st.ns || !st.complete)
			{
				Report("J build res=" + std::to_string(res), false, st.why);
				ok_bind = false;
				continue;
			}
			const NStar &c = *st.ns;
			std::vector<double> s0, sp0;
			ProdShape(c, s0, sp0);
			const Background2 bg = FromProfile(c, s0, sp0);
			const OracleSet hom = RunOracles(bg, 0.0, &knots, false, 1.0);
			const auto &P = c.Profile();
			const double ec = (*P.GetEnergyDensity())[0], pc = (*P.GetPressure())[0], dedpc = (*P.GetEosDEdP())[0];
			const double dp_c = (ec + pc); // p0*_c = 1
			const double rel_step = 1.0e-3;
			const SeqPt up = SequencePoint(cold, wrk / "jp", star16.ec_cgs * (1.0 + rel_step), res);
			const SeqPt dn = SequencePoint(cold, wrk / "jm", star16.ec_cgs * (1.0 - rel_step), res);
			Report("J res=" + std::to_string(res) + ": both sequence neighbours are complete ADR-0009 stars", up.complete && dn.complete, "");
			const double dM_dpc = (up.M_km - dn.M_km) / (up.pc - dn.pc);
			const double dedp_seq = (up.ec - dn.ec) / (up.pc - dn.pc);
			const double expP = dM_dpc * dp_c;
			const double rP = Rel(hom.stj_profile.deltaM_hat, expP), rK = Rel(hom.stj_knots.deltaM_hat, expP),
						 rS = Rel(hom.secant.deltaM_hat, expP), rD = Rel(hom.differential.deltaM_hat, expP);
			std::cout << "     res " << res << ": dM/dp_c=" << Sci(dM_dpc, 9) << " km^3  (deps/dp)_c sequence=" << Sci(dedp_seq, 7) << " vs authority=" << Sci(dedpc, 7)
					  << " rel=" << Sci(Rel(dedp_seq, dedpc)) << "\n"
					  << "                expected (dM/dp_c) dp_c=" << Sci(expP, 9) << "\n"
					  << "                homogeneous deltaM: Stieltjes profile K=4 " << Sci(hom.stj_profile.deltaM_hat, 9) << " rel=" << Sci(rP)
					  << " | EOS-knot K=2 " << Sci(hom.stj_knots.deltaM_hat, 9) << " rel=" << Sci(rK)
					  << " | secant (reported) rel=" << Sci(rS) << " | superseded differential (reported) rel=" << Sci(rD) << "\n";
			const bool bind = res <= 20000;
			const bool ok = hom.stj_profile.ok && hom.stj_knots.ok && std::isfinite(dM_dpc) && up.complete && dn.complete && rP <= kJ_Homog && rK <= kJ_Homog;
			if (bind)
			{
				ok_bind = ok_bind && ok;
				Report("Ja ADR-0008 Validation B at res " + std::to_string(res) + ": independent homogeneous deltaM vs (dM/dp_c) dp_c <= 2e-4 (profile and knot oracles)",
					   ok, "rel profile=" + Sci(rP) + "  knot=" + Sci(rK) + "  bound=" + Sci(kJ_Homog));
			}
			else
				std::cout << "     RECORD — res 40000 (diagnostic, not asserted): rel profile=" << Sci(rP) << "  knot=" << Sci(rK) << "\n";
		}
		Report("Jb the sequence-derivative identity holds at both binding resolutions", ok_bind, "");
		std::cout << "     (Test-side only — ADR-0007 P11 as modified at acceptance; no public homogeneous API.)\n";
	}

	fs::remove_all(wrk);
	std::cout << "\nSLOW-ROTATION DISCLAIMER: coefficient correctness, not truncation accuracy at high spin.\n";
	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
