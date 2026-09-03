// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file eos_derivative_cmf.cpp
 * @brief Phase 4C-I0 — the EOS thermodynamic derivative `d(eps)/dp` on the authenticated
 *        DS(CMF)-1 sequence (ADR-0007 P5, ACCEPTED 2026-09-02). External data required.
 *
 * The self-contained companion (`eos_derivative_contract.cpp`) proves the derivative is the
 * derivative of the star's own `eps(p)` interpolant, dimensionless, and fail-closed off-domain.
 * This file asks the remaining question: does that authority actually deliver a complete,
 * finite, physically sane value at every node of a *real* star built from the governed CompOSE
 * table — including across the crust, where the table is at its most awkward?
 *
 * ASSERTED (physical statements, safe by construction of a barotropic EOS):
 *   - every star publishes the derivative (`HasEosDEdP()`), one value per radial node;
 *   - zero non-finite values, hence zero unsupported-domain nodes: every profile pressure
 *     lay inside the interpolant's domain, which is what the fail-closed policy demands;
 *   - every value is strictly positive and finite — `d(eps)/dp = 1/c_s^2 > 0` for matter with
 *     a finite sound speed.
 *
 * DIAGNOSTIC, NOT ASSERTED: the retired profile finite-difference estimate
 * `(d eps/dr)/(dp/dr)`, computed here on the same nodes purely to show *why* ADR-0007 P5
 * removed it from authority. It is mathematically the same quantity by the chain rule on a
 * barotrope, so any disagreement is numerical noise — and the point of printing it is to show
 * how much noise there is. It is never compared against a bound and never used as an oracle.
 *
 * NOT TESTED HERE: any O(Omega^2) physics. The monopole solver is 4C-I1 and its validation is
 * 4D; INV-08 remains "candidate nonconformant / replacement pending".
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
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;

using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &detail)
{
	std::cout << (ok ? "  [ OK ] " : "  [FAIL] ") << id;
	if (!detail.empty())
		std::cout << "   " << detail;
	std::cout << "\n";
	if (!ok)
		++g_fail;
}

static std::string Sci(double v, int p = 4)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", p, v);
	return b;
}

static std::unique_ptr<NStar> Build(const fs::path &cold, const fs::path &wrk, double M)
{
	TOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	std::vector<TOVPoint> pts;
	std::vector<std::string> labels;
	const int n = tov.SolveToProfile(M, pts, &labels);
	if (n <= 0 || pts.size() < 4)
		return nullptr;
	return labels.empty() ? std::make_unique<NStar>(pts)
						  : std::make_unique<NStar>(pts, labels);
}

struct Row
{
	double M = 0, R = 0;
	std::size_t n = 0, n_nonfinite = 0, n_nonpositive = 0;
	double dmin = 0, dmax = 0, dmed = 0, d_centre = 0, d_surface = 0;
	double max_jump = 0, r_max_jump = 0;	 // largest node-to-node relative step, EOS-owned
	double fd_max_rel = 0, fd_med_rel = 0;	 // profile-FD deviation from the EOS-owned value
	double fd_max_jump = 0;					 // largest node-to-node relative step, profile-FD
	std::size_t fd_bad = 0;					 // FD nodes deviating by more than 1 %
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: eos_derivative_cmf <EOS_DATA_ROOT>\n";
		return 2;
	}
	const fs::path root = argv[1];
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "authenticated EOS missing: " << cold.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / "compactstar_eos_deriv_cmf";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4C-I0 — EOS d(eps)/dp on DS(CMF)-1_with_crust (ADR-0007 P5).\n"
				 "Authority = derivative of the star's own eps(p) Steffen interpolant.\n"
				 "The profile finite difference is printed as a DIAGNOSTIC only.\n\n";

	const double masses[] = {1.0, 1.4, 1.6, 2.0};
	std::vector<Row> rows;

	for (double M : masses)
	{
		auto ns = Build(cold, wrk, M);
		if (!ns)
		{
			Report("build " + Sci(M, 2) + " Msun", false, "solver returned no profile");
			continue;
		}

		const auto &P = ns->Profile();
		const auto *r = P.GetRadius();
		const auto *p = P.GetPressure();
		const auto *e = P.GetEnergyDensity();
		const auto *d = P.GetEosDEdP();

		Row row;
		row.M = M;

		const bool published = P.HasEosDEdP() && d != nullptr && r != nullptr &&
							   p != nullptr && e != nullptr;
		Report(Sci(M, 2) + " Msun: the profile publishes an EOS derivative", published,
			   published ? ("nodes=" + std::to_string(d->Size())) : "absent");
		if (!published)
			continue;

		const std::size_t n = static_cast<std::size_t>(r->Size());
		row.n = n;
		row.R = (*r)[static_cast<int>(n) - 1];

		Report(Sci(M, 2) + " Msun: one derivative value per radial node",
			   d->Size() == n,
			   "dEdP=" + std::to_string(d->Size()) + "  nodes=" + std::to_string(n));

		std::vector<double> vals;
		vals.reserve(n);
		for (std::size_t i = 0; i < n; ++i)
		{
			const double v = (*d)[static_cast<int>(i)];
			if (!std::isfinite(v))
				++row.n_nonfinite;
			else if (!(v > 0.0))
				++row.n_nonpositive;
			else
				vals.push_back(v);
		}
		row.d_centre = (*d)[0];
		row.d_surface = (*d)[static_cast<int>(n) - 1];

		if (!vals.empty())
		{
			std::vector<double> s = vals;
			std::sort(s.begin(), s.end());
			row.dmin = s.front();
			row.dmax = s.back();
			row.dmed = s[s.size() / 2];
		}

		// Largest node-to-node relative step in the authoritative value: this is the crust
		// showing up in the EOS itself, not a numerical artefact.
		for (std::size_t i = 1; i < n; ++i)
		{
			const double a = (*d)[static_cast<int>(i) - 1], b = (*d)[static_cast<int>(i)];
			if (!std::isfinite(a) || !std::isfinite(b) || a == 0.0)
				continue;
			const double jump = std::fabs(b - a) / std::fabs(a);
			if (jump > row.max_jump)
			{
				row.max_jump = jump;
				row.r_max_jump = (*r)[static_cast<int>(i)];
			}
		}

		// ---- DIAGNOSTIC: the retired profile finite difference -------------------------
		// Centred differences of the profile's own geometric columns, exactly as the
		// unratified O(Omega^2) candidate computes them (RotationSolver.cpp:1254-1277).
		// Both columns are km^-2, so the ratio is directly comparable to the EOS value.
		std::vector<double> fd_rel;
		fd_rel.reserve(n);
		double fd_prev = 0.0;
		bool have_prev = false;
		for (std::size_t i = 1; i + 1 < n; ++i)
		{
			const double dp = (*p)[static_cast<int>(i) + 1] - (*p)[static_cast<int>(i) - 1];
			const double de = (*e)[static_cast<int>(i) + 1] - (*e)[static_cast<int>(i) - 1];
			if (std::fabs(dp) <= 1.0e-30)
				continue; // the candidate substitutes 1.0 here; we simply skip the node
			const double fd = de / dp;
			const double auth = (*d)[static_cast<int>(i)];
			if (!std::isfinite(fd) || !std::isfinite(auth) || auth == 0.0)
				continue;
			const double rel = std::fabs(fd - auth) / std::fabs(auth);
			fd_rel.push_back(rel);
			if (rel > 1.0e-2)
				++row.fd_bad;
			if (have_prev && fd_prev != 0.0)
				row.fd_max_jump = std::max(row.fd_max_jump,
										   std::fabs(fd - fd_prev) / std::fabs(fd_prev));
			fd_prev = fd;
			have_prev = true;
		}
		if (!fd_rel.empty())
		{
			std::vector<double> s = fd_rel;
			std::sort(s.begin(), s.end());
			row.fd_med_rel = s[s.size() / 2];
			row.fd_max_rel = s.back();
		}
		// --------------------------------------------------------------------------------

		Report(Sci(M, 2) + " Msun: zero non-finite values (hence zero unsupported-domain "
						   "nodes)",
			   row.n_nonfinite == 0,
			   "non-finite=" + std::to_string(row.n_nonfinite) + " of " +
				   std::to_string(n));
		Report(Sci(M, 2) + " Msun: every value is strictly positive (d(eps)/dp = 1/c_s^2 > 0)",
			   row.n_nonpositive == 0,
			   "non-positive=" + std::to_string(row.n_nonpositive) + " of " +
				   std::to_string(n));

		rows.push_back(row);
	}

	// ------------------------------------------------------------------
	std::cout << "\nAuthoritative d(eps)/dp (dimensionless) along each profile\n";
	std::cout << "  M[Msun]    R[km]     nodes      min         median        max        "
				 " centre      surface    max node step (at r)\n";
	for (const auto &w : rows)
	{
		char b[512];
		snprintf(b, sizeof(b),
				 "  %5.2f   %8.5f   %6zu   %.4e  %.4e  %.4e  %.4e  %.4e   %.3e (%.4f km)\n",
				 w.M, w.R, w.n, w.dmin, w.dmed, w.dmax, w.d_centre, w.d_surface,
				 w.max_jump, w.r_max_jump);
		std::cout << b;
	}

	std::cout << "\nDIAGNOSTIC — retired profile finite difference (d eps/dr)/(dp/dr)\n"
				 "  Same quantity by the chain rule on a barotrope, so every deviation below is\n"
				 "  numerical noise in the retired method, not physics. ADR-0007 P5 removed it\n"
				 "  from authority for exactly this reason; it is not compared to any bound.\n";
	std::cout << "  M[Msun]   median |FD-EOS|/EOS   max |FD-EOS|/EOS   nodes off by >1%   "
				 "max FD node step\n";
	for (const auto &w : rows)
	{
		char b[512];
		snprintf(b, sizeof(b), "  %5.2f       %.4e            %.4e          %6zu / %zu      "
							   "  %.3e\n",
				 w.M, w.fd_med_rel, w.fd_max_rel, w.fd_bad, w.n, w.fd_max_jump);
		std::cout << b;
	}
	std::cout << "  (Compare the last column with the authoritative 'max node step' above: the\n"
				 "   EOS-owned value follows the crust; the finite difference amplifies it.)\n";

	std::cout << "\nc_s^2 cross-check: NOT AVAILABLE IN CURRENT GOVERNED INPUT — the imported\n"
				 "  table carries e(g/cm^3), p(dyne/cm^2), rho(1/fm^3) and species fractions\n"
				 "  only. 4C-I0 does not broaden the import path (ADR-0007 P5).\n";

	Report("all four masses produced a profile with a published derivative",
		   rows.size() == 4, "stars=" + std::to_string(rows.size()) + " of 4");

	fs::remove_all(wrk);

	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail
			  << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
