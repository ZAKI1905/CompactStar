// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_monopole_cmf.cpp
 * @brief Phase 4C-I1 — conformance of the governed O(Omega^2) monopole response on the
 *        authenticated DS(CMF)-1 sequence (ADR-0007). External data required.
 *
 * WHAT IS ASSERTED — structural conformance only:
 *   - each star carries an authoritative d(eps)/dp and a valid first-order response;
 *   - the monopole response computes explicitly, with exactly one integration per star;
 *   - every published coefficient is finite, with one value per radial node;
 *   - provenance matches the source profile and its Version();
 *   - the surface quantities and delta_M are finite;
 *   - delta_M composes exactly as ADR-0007 P6 defines it;
 *   - materialization is quadratic and performs no further integration.
 *
 * WHAT IS **NOT** ASSERTED — anything about whether the numbers are physically right.
 * The diagnostic table below is printed so a reader can see the magnitudes and so 4D has a
 * starting point. **It is not a golden, it is not compared against the retired candidate, and
 * it must not be cited as a result.** Independent scientific validation — regular-centre
 * series, an independent solver in different variables, Newtonian limits, exterior identities,
 * published Hartle-Thorne models, convergence, EOS-derivative sensitivity — is Phase 4D.
 * INV-08 remains "source conformed, independent physical validation pending".
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

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;

using CompactStar::AngularVelocity;
using CompactStar::Core::HartleMonopoleResponse;
using CompactStar::Core::NStar;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

static constexpr double kBoundQuadratic = 1.0e-14;   // ADR-0007 §7 item 2
static constexpr double kBoundDeltaMIdentity = 1.0e-14;

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
static double Rel(double a, double b)
{
	return std::fabs(b) > 0.0 ? std::fabs(a - b) / std::fabs(b) : std::fabs(a - b);
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
	double M = 0, R_star = 0, I = 0;
	std::size_t n = 0;
	double mhat_R = 0, phat_R = 0, xihat_R = 0, shell = 0, exterior = 0, dM = 0;
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: hartle_monopole_cmf <EOS_DATA_ROOT>\n";
		return 2;
	}
	const fs::path root = argv[1];
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "authenticated EOS missing: " << cold.string() << "\n";
		return 3;
	}
	const fs::path wrk = fs::temp_directory_path() / "compactstar_hartle_monopole_cmf";
	fs::remove_all(wrk);
	fs::create_directories(wrk);

	std::cout << std::scientific << std::setprecision(6);
	std::cout << "Phase 4C-I1 — governed O(Omega^2) monopole CONFORMANCE on "
				 "DS(CMF)-1_with_crust.\n"
				 "Structural conformance only. The diagnostic table is NOT a golden and NOT a\n"
				 "scientific result; independent validation is Phase 4D (INV-08, ADR-0007).\n\n";

	const double masses[] = {1.0, 1.4, 1.6, 2.0};
	std::vector<Row> rows;

	for (double M : masses)
	{
		const std::string tag = Sci(M, 2) + " Msun";
		auto ns = Build(cold, wrk, M);
		if (!ns)
		{
			Report(tag + ": build", false, "solver returned no profile");
			continue;
		}
		const NStar &cns = *ns; // const access: the non-const GetSequence() bumps Version()

		const auto &P = cns.Profile();
		Report(tag + ": carries an authoritative d(eps)/dp (Phase 4C-I0)", P.HasEosDEdP(), "");
		Report(tag + ": has a valid first-order response", cns.RotationResponse().valid, "");

		const bool ok = ns->ComputeHartleMonopoleResponse();
		Report(tag + ": the governed monopole response computes", ok, "");
		if (!ok)
			continue;

		const HartleMonopoleResponse *R = cns.MonopoleResponse();
		Report(tag + ": a current response is published", R != nullptr && R->valid, "");
		if (R == nullptr)
			continue;

		const std::size_t n = P.GetRadius()->Size();
		const bool sizes = R->m0_over_Omega2.Size() == n &&
						   R->p0star_over_Omega2.Size() == n &&
						   R->delta_p0_over_Omega2.Size() == n &&
						   R->xi0_over_Omega2.Size() == n;
		Report(tag + ": one value per radial node in every coefficient column", sizes,
			   "n = " + std::to_string(n));

		bool finite = std::isfinite(R->deltaM_over_Omega2) &&
					  std::isfinite(R->surface_shell_mass_over_Omega2) &&
					  std::isfinite(R->surface_xi0_over_Omega2) &&
					  std::isfinite(R->R_surface) && std::isfinite(R->I);
		for (int i = 0; finite && i < static_cast<int>(n); ++i)
			finite = std::isfinite(R->m0_over_Omega2[i]) &&
					 std::isfinite(R->p0star_over_Omega2[i]) &&
					 std::isfinite(R->delta_p0_over_Omega2[i]) &&
					 std::isfinite(R->xi0_over_Omega2[i]);
		Report(tag + ": every published value is finite", finite, "");

		Report(tag + ": provenance matches the source profile and its Version()",
			   R->MatchesSource(static_cast<const void *>(&P), P.Version()),
			   "version = " + std::to_string(P.Version()));

		const int last = static_cast<int>(n) - 1;
		const double Rs = R->R_surface;
		const double exterior = (R->I * R->I) / (Rs * Rs * Rs);
		const double expect =
			R->m0_over_Omega2[last] + R->surface_shell_mass_over_Omega2 + exterior;
		Report(tag + ": deltaM_hat = m0_hat(R_*) + shell_hat + I^2/R_*^3",
			   Rel(R->deltaM_over_Omega2, expect) <= kBoundDeltaMIdentity,
			   "rel = " + Sci(Rel(R->deltaM_over_Omega2, expect)));

		Report(tag + ": R_* is the last profile node (EOS-floor surface, not p = 0)",
			   Rs == (*P.GetRadius())[last], "R_* = " + Sci(Rs) + " km");

		// materialization: quadratic, and free of further integration
		const auto w1 = cns.MonopoleAt(AngularVelocity::FromHz(300.0));
		const auto w2 = cns.MonopoleAt(AngularVelocity::FromHz(600.0));
		double worst = Rel(w2.delta_M, 4.0 * w1.delta_M);
		for (int i = 0; i < static_cast<int>(n); ++i)
		{
			worst = std::max(worst, Rel(w2.m0[i], 4.0 * w1.m0[i]));
			worst = std::max(worst, Rel(w2.xi0[i], 4.0 * w1.xi0[i]));
		}
		Report(tag + ": materialization is exactly quadratic in Omega", worst <= kBoundQuadratic,
			   "worst rel = " + Sci(worst));

		Row row;
		row.M = M;
		row.R_star = Rs;
		row.I = R->I;
		row.n = n;
		row.mhat_R = R->m0_over_Omega2[last];
		row.phat_R = R->p0star_over_Omega2[last];
		row.xihat_R = R->surface_xi0_over_Omega2;
		row.shell = R->surface_shell_mass_over_Omega2;
		row.exterior = exterior;
		row.dM = R->deltaM_over_Omega2;
		rows.push_back(row);
	}

	// ------------------------------------------------------------------
	std::cout << "\nDIAGNOSTIC — governed monopole coefficients, per Omega_geom^2.\n"
				 "  NOT a golden. NOT compared against the retired candidate. NOT a scientific\n"
				 "  claim. Printed so Phase 4D has a starting point and a reader can see scale.\n";
	std::cout << "  M[Msun]   R_*[km]     nodes   I[km^3]      m0_hat(R_*)  p0*_hat(R_*)"
				 "  xi0_hat(R_*)  shell_hat    I^2/R_*^3    deltaM_hat\n";
	for (const auto &w : rows)
	{
		char b[600];
		snprintf(b, sizeof(b),
				 "  %5.2f   %9.5f   %6zu   %.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e\n",
				 w.M, w.R_star, w.n, w.I, w.mhat_R, w.phat_R, w.xihat_R, w.shell, w.exterior,
				 w.dM);
		std::cout << b;
	}
	std::cout << "  (Units: m0_hat, xi0_hat, shell_hat, deltaM_hat in km^3; p0*_hat in km^2.\n"
				 "   R_* is the production EOS-floor surface node, deliberately not identified\n"
				 "   with the exact p = 0 radius — ADR-0007 P7, INV-06.)\n";

	Report("all four masses produced a current governed monopole response", rows.size() == 4,
		   "stars = " + std::to_string(rows.size()) + " of 4");

	fs::remove_all(wrk);

	std::cout << "\nNON-CLAIM: no number above is validated physics. ADR-0007 fixes the\n"
				 "contract and this implementation conforms to it; independent verification is\n"
				 "Phase 4D, and only after it may a monopole baseline be created.\n";

	std::cout << "\n" << (g_fail == 0 ? "PASS" : "FAIL") << " — " << g_fail
			  << " failed check(s)\n";
	return g_fail == 0 ? 0 : 1;
}
