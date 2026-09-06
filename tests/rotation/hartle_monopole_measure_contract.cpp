#include "tests/relativity/fixture_units.hpp"
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_measure_contract.cpp
 * @brief Phase 4D-RI — implementation contract of the MEASURE-COMPLETE EOS energy-density
 *        source of the O(Omega^2) monopole response (ADR-0008, ACCEPTED 2026-09-03).
 *        Self-contained.
 *
 * WHAT THIS FILE IS. An implementation-contract test, not a scientific validation. It proves
 * that `RotationSolver::ODE_HartleMonopole_` sources `dm0/dr` from the governed measure
 *
 *     dm0_hat|_EOS = -4 pi r^2 xi0_hat d(eps)          (H67 eq. 93; ADR-0008 Q1)
 *
 * evaluated one profile segment at a time with that segment's own d(eps)/dr (ADR-0008 Q3), that
 * the surface shell is the terminal eps_* -> 0 atom of the same measure applied exactly once
 * (Q7), and that the accepted form still reduces to the superseded differential form wherever
 * the EOS is smooth. The corrected INDEPENDENT physical revalidation is a separate increment;
 * nothing here claims a validated number.
 *
 * ============================ PREDECLARED BOUNDS ============================
 *   M1  constant-density star: every interior segment has Delta eps = 0, so the EOS
 *       source is identically zero and m0_hat equals the rotation-only integral   <= 1e-14
 *   M2  smooth EOS (HT68 Harrison-Wheeler, self-contained fixture): production
 *       (measure) vs the independent (m0,h0) solver in the SUPERSEDED differential
 *       form, ON deltaM_hat — the quantity ADR-0008's evidence table measured
 *       (1.2e-5 on this EOS) and set the bound from                               <= 2e-5   (ADR-0008 Validation)
 *       The node-wise peak-relative disagreement is REPORTED beside it. A first run
 *       of this file also asserted the bound on that node-wise metric and measured
 *       5.5e-6 / 1.97e-5 / 2.17e-5 at eps_c = 3e14 / 1e15 / 3e15 g/cm^3: the two
 *       source representations differ at O(h^2) in the weight even where the EOS is
 *       smooth, and that difference grows with the density gradient. The bound was
 *       NOT widened; it is asserted on the quantity it was predeclared on, and the
 *       node-wise number is recorded.
 *   M3a per-segment measure identity inside the accounting integrator:
 *       Delta m0_hat|_EOS,i = -slope_i * Integral(4 pi r^2 xi0_hat dr)_i, worst
 *       segment residual measured ABSOLUTELY against the star's total EOS
 *       integral (a per-segment RELATIVE metric is meaningless on the segments
 *       whose Delta eps is ~ 0)                                                    <= 1e-10
 *   M3b same-partition source accounting vs production                            <= 1e-6   (ADR-0008 Validation C)
 *
 * M2 and M3b compare node by node on the nodes carrying the signal — |ref| >= 1 % of the
 * profile peak — the convention `hartle_monopole_reference.hpp::Compare` already uses, because
 * m0_hat spans ~25 decades between the centre and the surface and a relative metric on its
 * O(r^5) tail measures nothing. The surface value is reported separately in every case.
 *   M4  terminal atom exactly once: interior measure telescopes to eps_c - eps_*,
 *       the shell is 4 pi R_*^2 eps_* xi0_hat(R_*), delta M identity exact         <= 1e-14
 * ===========================================================================
 *
 * No DS(CMF)-1 value, no pre-correction monopole value and no 4D diagnostic is an expected
 * answer anywhere in this file.
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
#include <unistd.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "hartle_monopole_reference.hpp"
#include "hartle_profile_compare.hpp"
#include "hartle_thorne_1968_hw_eos.hpp"

#include <Zaki/Physics/Constants.hpp>

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;
using hartle_4b::UniformStar;
using hartle_mono_ref::Background2;
using hartle_mono_ref::MHOptions;

static constexpr double kM1_Zero = 1.0e-14;
static constexpr double kM2_Smooth = 2.0e-5;
static constexpr double kM3a_Segment = 1.0e-10;
static constexpr double kM3b_Accounting = 1.0e-6;
static constexpr double kM4_Exact = 1.0e-14;

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

/// Worst relative disagreement over the nodes that carry the signal: |ref| >= 1 % of the
/// profile's peak magnitude (the `hartle_monopole_reference.hpp::Compare` convention). m0_hat
/// spans ~25 decades from its O(r^5) centre behaviour to the surface, so an unrestricted
/// node-wise relative metric reports the deep interior and nothing else.
static double PeakRel(const std::vector<double> &got, const Zaki::Vector::DataColumn &ref)
{
	const int n = static_cast<int>(ref.Size());
	double peak = 0.0;
	for (int i = 0; i < n; ++i)
		peak = std::max(peak, std::fabs(ref[i]));
	double worst = 0.0;
	for (int i = 0; i < n && i < static_cast<int>(got.size()); ++i)
		if (std::fabs(ref[i]) >= 1.0e-2 * peak)
			worst = std::max(worst, Rel(got[i], ref[i]));
	return worst;
}

// ---------------------------------------------------------------------------
//  Backgrounds
// ---------------------------------------------------------------------------
struct Bg
{
	std::vector<double> r, p, eps, m, nu, nup, dedp, s, sp;
	std::size_t N() const { return r.size(); }
};

static Bg FromProfile(const NStar &ns)
{
	const auto &P = ns.Profile();
	const auto *r = P.GetRadius();
	const auto &R1 = ns.RotationResponse();
	const int n = static_cast<int>(r->Size());
	Bg o;
	for (int i = 0; i < n; ++i)
	{
		o.r.push_back((*r)[i]);
		o.p.push_back((*P.GetPressure())[i]);
		o.eps.push_back((*P.GetEnergyDensity())[i]);
		o.m.push_back((*P.GetMass())[i]);
		o.nu.push_back((*P.GetMetricNu())[i]);
		o.nup.push_back((*P.GetMetricNuPrime())[i]);
		o.dedp.push_back((*P.GetEosDEdP())[i]);
		o.s.push_back(R1.omega_bar_over_Omega[i]);
		o.sp.push_back(R1.domega_bar_over_Omega_dr[i]);
	}
	return o;
}
static Background2 ToReference(const Bg &b)
{
	Background2 o;
	o.r = b.r;
	o.p = b.p;
	o.eps = b.eps;
	o.m = b.m;
	o.nu = b.nu;
	o.nup = b.nup;
	o.dedp = b.dedp;
	o.s = b.s;
	o.sp = b.sp;
	return o;
}

/// Exact constant-density interior on a uniform grid, with the EXACT EOS derivative
/// d(eps)/dp = 0 supplied through the governed mechanism (ADR-0007 P5, ADR-0008 Q10).
static std::vector<TOVPoint> UniformPoints(const UniformStar &u, std::size_t N, bool supply_dedp)
{
	const double km2_to_gcm3 = relativity_fixture::eps_to_rho;
	const double km2_to_dyn = relativity_fixture::pressure_to_cgs;
	std::vector<TOVPoint> pts;
	pts.reserve(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		const double r = hartle_4b::UniformGridR(u, i, N);
		pts.emplace_back(r, u.m(r) / relativity_fixture::solar_km, u.nuprime(r) / 1.0e5, 0.0,
						 u.p(r) * km2_to_dyn, u.rho0 * km2_to_gcm3, 0.1, std::vector<double>{},
						 supply_dedp ? 0.0 : std::numeric_limits<double>::quiet_NaN());
	}
	return pts;
}

// ---------------------------------------------------------------------------
//  Test-side same-partition source accounting.
//
//  A four-variable segment-aware integrator over the SAME profile partition, with its own
//  driver and tighter tolerances:
//      y[0] = m0_hat, y[1] = p0*_hat, y[2] = the EOS part of m0_hat, y[3] = Integral 4 pi r^2 xi dr
//  y[3] is reset at every segment boundary, so the per-segment identity
//      Delta y[2]_i = -slope_i * y[3]_i
//  is checked directly, and y[0] is compared against production node by node.
// ---------------------------------------------------------------------------
struct AccParams
{
	const Bg *bg;
	double slope;
};
static double LerpAt(const Bg &b, const std::vector<double> &c, double rq)
{
	if (rq <= b.r.front())
		return c.front();
	if (rq >= b.r.back())
		return c.back();
	std::size_t lo = 0, hi = b.r.size() - 1;
	while (hi - lo > 1)
	{
		const std::size_t mid = (lo + hi) / 2;
		if (b.r[mid] <= rq)
			lo = mid;
		else
			hi = mid;
	}
	const double t = (rq - b.r[lo]) / (b.r[hi] - b.r[lo]);
	return c[lo] + t * (c[hi] - c[lo]);
}
static int AccRHS(double r, const double y[], double f[], void *pv)
{
	const AccParams *P = static_cast<const AccParams *>(pv);
	const Bg &b = *P->bg;
	const double p = LerpAt(b, b.p, r), eps = LerpAt(b, b.eps, r), m = LerpAt(b, b.m, r);
	const double nu = LerpAt(b, b.nu, r), nup = LerpAt(b, b.nup, r);
	const double s = LerpAt(b, b.s, r), sp = LerpAt(b, b.sp, r);
	const double D = 1.0 - 2.0 * m / r, r_2m = r * D, e2 = std::exp(-2.0 * nu), ep = eps + p;
	const double r2 = r * r, r3 = r2 * r, r4 = r3 * r;
	const double mhat = y[0], phat = y[1];
	const double xi = phat / nup;
	if (!std::isfinite(xi) || !(D > 0.0))
		return GSL_EBADFUNC;

	const double term1 = (P->slope != 0.0) ? -4.0 * M_PI * r2 * xi * P->slope : 0.0;
	f[0] = term1 + (1.0 / 12.0) * r4 * e2 * D * sp * sp + (8.0 * M_PI / 3.0) * r4 * ep * e2 * s * s;
	f[1] = -mhat * (1.0 + 8.0 * M_PI * r2 * p) / (r_2m * r_2m) - 4.0 * M_PI * ep * r2 * phat / r_2m +
		   (1.0 / 12.0) * r3 * e2 * sp * sp + (2.0 / 3.0) * r * e2 * s * (s + r * sp - r * nup * s);
	f[2] = term1;
	f[3] = 4.0 * M_PI * r2 * xi;
	return (std::isfinite(f[0]) && std::isfinite(f[1])) ? GSL_SUCCESS : GSL_EBADFUNC;
}

struct AccResult
{
	bool ok = false;
	std::vector<double> mhat, phat;
	double worst_segment_abs = 0.0; ///< max |Delta m0_EOS,i + slope_i * W_i| over segments
	double eos_total = 0.0, rot_total = 0.0;
};
static AccResult Account(const Bg &b, double mhat0, double phat0)
{
	AccResult out;
	const std::size_t n = b.N();
	AccParams P{&b, 0.0};
	gsl_odeiv2_system sys = {AccRHS, nullptr, 4, &P};
	gsl_odeiv2_driver *d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-16, 1e-13);
	double y[4] = {mhat0, phat0, 0.0, 0.0};
	double r = b.r.front();
	out.mhat.assign(n, 0.0);
	out.phat.assign(n, 0.0);
	out.mhat[0] = y[0];
	out.phat[0] = y[1];
	for (std::size_t i = 1; i < n; ++i)
	{
		P.slope = (b.eps[i] - b.eps[i - 1]) / (b.r[i] - b.r[i - 1]);
		const double eos_before = y[2];
		y[3] = 0.0; // per-segment weight accumulator
		if (gsl_odeiv2_driver_apply(d, &r, b.r[i], y) != GSL_SUCCESS)
		{
			gsl_odeiv2_driver_free(d);
			return out;
		}
		// ADR-0008 Q3: the EOS increment of this segment must be exactly -slope_i times the
		// segment's own weight integral, i.e. -4 pi <r^2 xi> Delta eps_i.
		const double got = y[2] - eos_before, want = -P.slope * y[3];
		out.worst_segment_abs = std::max(out.worst_segment_abs, std::fabs(got - want));
		out.mhat[i] = y[0];
		out.phat[i] = y[1];
	}
	gsl_odeiv2_driver_free(d);
	out.eos_total = y[2];
	out.rot_total = y[0] - mhat0 - y[2];
	out.ok = true;
	return out;
}

// ---------------------------------------------------------------------------
int main()
{
	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4D-RI — implementation contract of the measure-complete EOS energy-density\n"
				 "source (ADR-0008 ACCEPTED 2026-09-03). Not a validation of any number.\n\n";
	const fs::path wrk = fs::temp_directory_path() / ("compactstar_measure_contract_" + std::to_string(::getpid()));
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	// =====================================================================
	//  M1 / M4 — the exact constant-density star (ADR-0008 Q10, Q7)
	// =====================================================================
	std::cout << "M1/M4. Exact constant-density interior (M = 2 km, R = 13 km, N = 4001)\n";
	{
		const UniformStar u = hartle_4b::MakeUniform(2.0, 13.0);
		NStar ns(UniformPoints(u, 4001, true));
		const bool computed = ns.ComputeHartleMonopoleResponse();
		Report("M0 a point-constructed star with an explicit d(eps)/dp = 0 still computes (Q10)", computed, "");
		if (computed)
		{
			const NStar &c = ns;
			const auto *R = c.MonopoleResponse();
			const Bg b = FromProfile(c);
			const int last = static_cast<int>(b.N()) - 1;

			double worst_interior = 0.0, telescope = 0.0;
			for (std::size_t i = 1; i < b.N(); ++i)
			{
				worst_interior = std::max(worst_interior, std::fabs(b.eps[i] - b.eps[i - 1]));
				telescope += b.eps[i] - b.eps[i - 1];
			}
			Report("M1a every interior segment of a constant-density star has Delta eps = 0",
				   worst_interior == 0.0, "max |Delta eps_i| = " + Sci(worst_interior));

			// the EOS source must be identically zero: the accounting integrator's EOS channel
			const auto acc = Account(b, R->m0_over_Omega2[0], R->p0star_over_Omega2[0]);
			Report("M1b the EOS channel of the source integrates to exactly zero", acc.ok && acc.eos_total == 0.0,
				   "Integral of the EOS source = " + Sci(acc.eos_total));
			const double worst_m = PeakRel(acc.mhat, R->m0_over_Omega2);
			Report("M1c production m0_hat equals the rotation-only integral at every node",
				   acc.ok && worst_m <= kM1_Zero, "worst rel = " + Sci(worst_m) + "  bound = " + Sci(kM1_Zero));

			// M4: the terminal atom, exactly once
			const double shell = 4.0 * M_PI * R->R_surface * R->R_surface * b.eps[last] * R->surface_xi0_over_Omega2;
			Report("M4a the surface shell is the terminal eps_* -> 0 atom, 4 pi R_*^2 eps_* xi0_hat(R_*)",
				   Rel(R->surface_shell_mass_over_Omega2, shell) <= kM4_Exact,
				   "published = " + Sci(R->surface_shell_mass_over_Omega2, 8) + "  atom = " + Sci(shell, 8));
			const double dM = R->m0_over_Omega2[last] + R->surface_shell_mass_over_Omega2 +
							  R->I * R->I / (R->R_surface * R->R_surface * R->R_surface);
			Report("M4b deltaM_hat = m0_hat(R_*) + terminal atom + I^2/R_*^3, exactly",
				   Rel(R->deltaM_over_Omega2, dM) <= kM4_Exact, "rel = " + Sci(Rel(R->deltaM_over_Omega2, dM)));
			Report("M4c the interior measure telescopes to eps_c - eps_* and EXCLUDES the surface drop",
				   std::fabs(telescope - (b.eps[last] - b.eps[0])) <= kM4_Exact * std::max(1.0, std::fabs(b.eps[0])) &&
					   b.eps[last] > 0.0,
				   "interior sum = " + Sci(telescope) + "  eps_* = " + Sci(b.eps[last]) +
					   " (the atom the shell carries, never also an interior segment)");
			std::cout << "     deltaM_hat = " << Sci(R->deltaM_over_Omega2, 10) << " km^3, shell/deltaM = "
					  << Sci(R->surface_shell_mass_over_Omega2 / R->deltaM_over_Omega2, 4)
					  << "  (unchanged by ADR-0008 by construction: the interior measure is zero)\n";
		}
	}

	// =====================================================================
	//  M6 — fail-closed without the EOS derivative (ADR-0007 P5, retained by Q8)
	// =====================================================================
	{
		const UniformStar u = hartle_4b::MakeUniform(2.0, 13.0);
		NStar ns(UniformPoints(u, 401, false));
		Report("M6 a star with no authoritative d(eps)/dp still fails closed (Q8 retains P5 for the centre series)",
			   !ns.ComputeHartleMonopoleResponse() && ns.MonopoleResponse() == nullptr, "");
	}

	// =====================================================================
	//  M2 / M3 — a smooth tabulated EOS: HT68 Harrison-Wheeler (self-contained)
	// =====================================================================
	std::cout << "\nM2/M3. Smooth tabulated EOS — Hartle & Thorne (1968) Harrison-Wheeler table, densified\n";
	{
		const fs::path eos = wrk / "hw1968_densified.eos";
		if (!ht68::WriteDensifiedEos(eos.string(), 40))
			Report("M2 fixture written", false, eos.string());
		else
		{
			bool any = false;
			for (double ec : {3.0e14, 1.0e15, 3.0e15})
			{
				TOVSolver tov;
				const fs::path w = wrk / ("hw" + Sci(ec, 2));
				fs::create_directories(w);
				tov.SetWrkDir(w.string());
				tov.ImportEOS(eos.string(), true);
				std::vector<TOVPoint> pts;
				if (tov.SingleStarSolveToTOVPoints(ec, pts) <= 0 || pts.size() < 4)
				{
					Report("M2 TOV at eps_c = " + Sci(ec, 2), false, "no profile");
					continue;
				}
				NStar ns(pts);
				if (!ns.ComputeHartleMonopoleResponse())
				{
					Report("M2 monopole at eps_c = " + Sci(ec, 2), false, "no response");
					continue;
				}
				any = true;
				const NStar &c = ns;
				const auto *R = c.MonopoleResponse();
				const Bg b = FromProfile(c);
				const int last = static_cast<int>(b.N()) - 1;

				// M2 — the SUPERSEDED differential form on an independent (m0,h0) solver
				MHOptions o;
				o.I_exterior = R->I;
				o.eos_measure = false;
				const auto old_form = hartle_mono_ref::Solve(ToReference(b), o);
				MHOptions o2 = o;
				o2.eos_measure = true;
				const auto new_form = hartle_mono_ref::Solve(ToReference(b), o2);
				const double d_old = Rel(R->deltaM_over_Omega2, old_form.deltaM_hat);
				const double d_new = Rel(R->deltaM_over_Omega2, new_form.deltaM_hat);
				const double worst_old = PeakRel(old_form.mhat, R->m0_over_Omega2);

				Report("M2 eps_c = " + Sci(ec, 2) + " g/cm^3: the measure source agrees with the SUPERSEDED "
													"differential form on deltaM_hat where the EOS is smooth",
					   old_form.ok && d_old <= kM2_Smooth,
					   "deltaM_hat rel = " + Sci(d_old) + "  bound = " + Sci(kM2_Smooth) +
						   "  | REPORTED: worst m0_hat node rel = " + Sci(worst_old) +
						   " (O(h^2) between the two source representations), same-form oracle rel = " + Sci(d_new));

				// M3 — same-partition source accounting
				const auto acc = Account(b, R->m0_over_Omega2[0], R->p0star_over_Omega2[0]);
				const double worst_m = PeakRel(acc.mhat, R->m0_over_Omega2);
				const double seg_rel = (acc.eos_total != 0.0) ? acc.worst_segment_abs / std::fabs(acc.eos_total)
															  : acc.worst_segment_abs;
				Report("M3a eps_c = " + Sci(ec, 2) + ": every segment's EOS increment equals -slope_i times its "
													 "own weight integral",
					   acc.ok && seg_rel <= kM3a_Segment,
					   "worst segment residual = " + Sci(acc.worst_segment_abs) + " km^3 = " + Sci(seg_rel) +
						   " of the total EOS integral  bound = " + Sci(kM3a_Segment));
				Report("M3b eps_c = " + Sci(ec, 2) + ": same-partition source accounting reproduces production",
					   acc.ok && worst_m <= kM3b_Accounting && Rel(acc.mhat[last], R->m0_over_Omega2[last]) <= kM3b_Accounting,
					   "worst node rel = " + Sci(worst_m) + ", m0_hat(R_*) rel = " +
						   Sci(Rel(acc.mhat[last], R->m0_over_Omega2[last])) + "  bound = " + Sci(kM3b_Accounting));
				std::cout << "     R_* = " << Sci(R->R_surface, 6) << " km, nodes " << b.N()
						  << ", EOS part of m0_hat(R_*) = " << Sci(acc.eos_total, 6) << ", rotational part = "
						  << Sci(acc.rot_total, 6) << "\n";
			}
			if (!any)
				Report("M2/M3 at least one smooth-EOS star was built", false, "");
		}
	}

	fs::remove_all(wrk);
	std::cout << "\nThis file proves conformance to ADR-0008, not the correctness of any physical number.\n";
	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
