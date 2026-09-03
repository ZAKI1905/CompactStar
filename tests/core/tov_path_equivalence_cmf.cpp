// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file tov_path_equivalence_cmf.cpp
 * @brief Phase 3E-0 — measure the observable numerical relationship between the two LIVE TOV
 *        paths, as they exist today. This designates no canonical owner.
 *
 * THE TWO PATHS.
 *
 *   Path 1 (sequence / legacy-live):
 *       TOVSolver::Solve(Axis, dir, file)
 *         -> per Axis node: clamp ec, p_of_e(ec), RadiusLoop(r, y)
 *              -> NStar::Append(TOVPoint) per radial step
 *         -> SurfaceIsReached() -> NStar::FinalizeSurface()
 *         -> analysis? -> n_exp_cond_f? -> ExportNStarProfile(idx, dir)
 *         -> Sequence::Add(n_star) -> n_star.Reset()
 *
 *   Path 2 (validated profile):
 *       TOVSolver::SingleStarSolveToTOVPoints(ec, out_points)
 *         -> vector<TOVPoint>
 *       NStar(points, species_labels) -> NStar::BuildFromTOV(...)
 *
 * Every Phase-2 scientific harness protects Path 2. Path 1 had no coverage before this test.
 *
 * WHAT THIS TEST DOES NOT DO. It does not reimplement or mock the TOV equations, does not
 * modify production source, and does not compare Path 1 against Path 2's *mass-search*
 * orchestration (`SolveToProfile`), which Path 1 does not possess. The equivalence authority is:
 * same EOS, same central energy density, same numerical controls, same radial solve, same
 * finalization.
 *
 * CAPTURE MECHANISM — OBSERVATIONAL ONLY. `Solve()` destroys each star (`n_star.Reset()`) right
 * after adding it to the sequence, so the finished star must be read at the one hook that fires
 * in between. `CapturingTOVSolver` below:
 *   - enables the EXISTING export-condition callback (`n_exp_cond_f`) for every star;
 *   - overrides the EXISTING virtual `ExportNStarProfile`, which `Solve()` already calls after
 *     `SurfaceIsReached()` and before `Sequence::Add` / `n_star.Reset()`;
 *   - copies scalar and column VALUES out immediately (never retaining a pointer into
 *     `n_star`, which is reset, and never copying `NStar`, which is non-copyable);
 *   - deliberately does not call the base implementation, so no profile file is written.
 *
 * It does NOT override `RadiusLoop`, `SurfaceIsReached` or `ODE`, does not touch tolerances,
 * the pressure cutoff or the step schedule, and does not fill `NStar` directly. Nothing it does
 * executes before or during the integration, so it cannot perturb the numbers it reads.
 *
 * SPECIES LABELS are read from the solver's own `eos_tab.extra_labels` — the identical source
 * `NStar::Reset` uses on Path 1 and `SolveToProfile` returns on Path 2 — so both paths are
 * compared under one label ordering (ADR-0001).
 *
 * PREDECLARED EQUIVALENCE CRITERIA (fixed before any measurement, from the source audit):
 *
 *   - raw radial columns (r, m, nu', p, eps, nB, species) : BIT-IDENTICAL
 *       Both loops are textually the same algorithm on the same `TOVSolver::ODE`, the same
 *       rk8pd driver with the same (1e-1, 1e-10, 1e-10), the same linear grid, the same
 *       step_scale ladder and the same append-then-test-cutoff order.
 *   - Lambda                                             : BIT-IDENTICAL
 *       `NStar::Append` and `BuildFromTOV` compute the identical expression; the legacy 1e-15
 *       clamp Path 1 still carries is unreachable (max 2m/r = 0.481, ADR-0004 §11).
 *   - nu                                                 : BIT-IDENTICAL
 *       Both call the same `EvaluateNu()`.
 *   - SeqPoint ec, m, r, pc, I                           : BIT-IDENTICAL
 *   - SeqPoint b                     : |db|/b <= 1.0e-15
 *       ADR-0004 migrated ONLY Path 2's proper-volume integrand to the canonical primitive;
 *       Path 1's `FinalizeSurface` still uses the legacy inline `/(1-2m/r).sqrt()`. This is a
 *       governed, intentional conformance gap with a bound predeclared in ADR-0004 §7.1
 *       BEFORE that implementation existed. It is NOT a TOV-integrator failure and is NOT
 *       repaired here.
 *
 * ONE KNOWN NON-NUMERICAL ASYMMETRY, predicted by the source audit before measurement:
 * `StarProfile`'s mirror surface scalars (`M`, `R`, `z_surf`, reachable via `MassSurface()`,
 * `RadiusSurface()`, `ExpNuSurface()`) are set by `BuildFromTOV` but never by
 * `FinalizeSurface`, so Path 1 leaves them at the zeros `NStar::Reset` installed. The SeqPoint
 * `m`/`r` are correct on both paths. This is PROFILE CONSTRUCTION, not radial integration, and
 * it is asserted here as the current characterized state so that any change is caught.
 *
 * No tolerance is chosen after seeing output. A quantity that misses its predeclared standard is
 * reported as measured and classified no more strongly than CHARACTERIZED.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::SeqPoint;
using CompactStar::Core::StarProfile;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

// ===========================================================================
//  Reporting
// ===========================================================================
namespace
{
int g_fail = 0;
void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}
std::string Sci(double v, int prec = 3)
{
	char b[80];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
}
long UlpDistance(double a, double b)
{
	if (a == b)
		return 0;
	if (!std::isfinite(a) || !std::isfinite(b))
		return std::numeric_limits<long>::max();
	std::int64_t ia = 0, ib = 0;
	std::memcpy(&ia, &a, sizeof ia);
	std::memcpy(&ib, &b, sizeof ib);
	if ((ia < 0) != (ib < 0))
		return std::numeric_limits<long>::max();
	return static_cast<long>(ia > ib ? ia - ib : ib - ia);
}

/// Predeclared: the ONLY non-bit-identical quantity, and only because ADR-0004 deliberately
/// migrated Path 2's proper-volume integrand while leaving Path 1's legacy inline form.
constexpr double kBTolRel = 1.0e-15;

// ===========================================================================
//  Snapshot — plain values copied out of a finished NStar
// ===========================================================================
struct ColumnCmp
{
	std::string label;
	std::size_t n = 0, bitwise = 0;
	double max_abs = 0.0, max_rel = 0.0;
	long max_ulp = 0;
	std::size_t worst_i = 0;
	double worst_r = 0.0;
};

struct Snapshot
{
	bool surface_ready = false;
	bool valid = false;
	SeqPoint seq;
	double surf_M = 0.0, surf_R = 0.0, surf_z = 0.0;
	std::vector<std::string> col_labels;
	std::vector<std::vector<double>> cols;
	std::vector<std::string> species_labels;
	std::vector<int> species_cols;
	int idx_r = -1, idx_m = -1, idx_nup = -1, idx_p = -1;
	int idx_eps = -1, idx_nb = -1, idx_nu = -1, idx_lam = -1;

	std::size_t N() const { return cols.empty() ? 0 : cols[0].size(); }
	const std::vector<double> &C(int i) const { return cols.at(static_cast<std::size_t>(i)); }
};

/// Copies VALUES out of a finished star. Retains no reference into it.
Snapshot Capture(const NStar &ns)
{
	Snapshot s;
	const StarProfile &p = ns.Profile();
	s.seq = p.Seq();
	s.surf_M = p.MassSurface();
	s.surf_R = p.RadiusSurface();
	s.surf_z = p.ExpNuSurface();

	const auto &radial = p.Radial();
	const std::size_t ncol = radial.Dim().size();
	s.cols.resize(ncol);
	s.col_labels.resize(ncol);
	for (std::size_t c = 0; c < ncol; ++c)
	{
		s.col_labels[c] = radial[static_cast<int>(c)].Label();
		const std::size_t n = radial[static_cast<int>(c)].Size();
		s.cols[c].resize(n);
		for (std::size_t i = 0; i < n; ++i)
			s.cols[c][i] = radial[static_cast<int>(c)][i];
	}

	s.idx_r = p.GetColumnIndex(StarProfile::Column::Radius);
	s.idx_m = p.GetColumnIndex(StarProfile::Column::Mass);
	s.idx_nup = p.GetColumnIndex(StarProfile::Column::MetricNuPrime);
	s.idx_p = p.GetColumnIndex(StarProfile::Column::Pressure);
	s.idx_eps = p.GetColumnIndex(StarProfile::Column::EnergyDensity);
	s.idx_nb = p.GetColumnIndex(StarProfile::Column::BaryonDensity);
	s.idx_nu = p.GetColumnIndex(StarProfile::Column::MetricNu);
	s.idx_lam = p.GetColumnIndex(StarProfile::Column::MetricLambda);

	const std::size_t nsp = p.SpeciesCount();
	for (std::size_t j = 0; j < nsp; ++j)
	{
		const int c = p.SpeciesColumnIndex(j);
		s.species_cols.push_back(c);
		s.species_labels.push_back(
			(c >= 0 && static_cast<std::size_t>(c) < ncol) ? s.col_labels[c] : "<none>");
	}
	s.valid = !s.cols.empty() && s.N() > 0;
	return s;
}

ColumnCmp CompareColumn(const std::string &label, const std::vector<double> &a,
						const std::vector<double> &b, const std::vector<double> &r)
{
	ColumnCmp c;
	c.label = label;
	c.n = std::min(a.size(), b.size());
	for (std::size_t i = 0; i < c.n; ++i)
	{
		if (a[i] == b[i])
		{
			++c.bitwise;
			continue;
		}
		const double d = std::fabs(a[i] - b[i]);
		const double rel = std::fabs(a[i]) > 0.0 ? d / std::fabs(a[i]) : d;
		const long u = UlpDistance(a[i], b[i]);
		if (u > c.max_ulp)
		{
			c.max_ulp = u;
			c.worst_i = i;
			c.worst_r = (i < r.size() ? r[i] : 0.0);
		}
		c.max_abs = std::max(c.max_abs, d);
		c.max_rel = std::max(c.max_rel, rel);
	}
	return c;
}

// ===========================================================================
//  Test-only capturing solver — observational hook, no production change
// ===========================================================================
class CapturingTOVSolver : public TOVSolver
{
  public:
	std::vector<Snapshot> captured;

	CapturingTOVSolver()
	{
		// Enable the EXISTING export-condition callback for every star, so the existing
		// ExportNStarProfile hook fires. A captureless lambda converts to the raw function
		// pointer the production member expects.
		n_exp_cond_f = [](const NStar &) { return true; };
	}

	/// Same source `NStar::Reset` reads on Path 1 and `SolveToProfile` returns on Path 2.
	std::vector<std::string> SpeciesLabels() const { return eos_tab.extra_labels; }
	double EosEpsFront() const { return eos_tab.eps.front(); }
	double EosEpsBack() const { return eos_tab.eps.back(); }
	double CentralFloorFactor() const { return central_eps_floor_factor; }

  protected:
	/// Fires after SurfaceIsReached() and before Sequence::Add / n_star.Reset().
	/// Copies values out; deliberately does NOT call the base, so no file is written.
	void ExportNStarProfile(const size_t &, const Zaki::String::Directory &) override
	{
		captured.push_back(Capture(n_star));
	}
};

// ===========================================================================
//  Path drivers
// ===========================================================================
/// Path 1: the real `Solve()` over an Axis, capturing every finished star.
std::vector<Snapshot> RunPath1(const fs::path &cold, const fs::path &wrk,
							   const Zaki::Math::Axis &ax, std::size_t res,
							   std::vector<std::string> *labels_out = nullptr)
{
	CapturingTOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	tov.SetRadialRes(res);
	if (labels_out)
		*labels_out = tov.SpeciesLabels();
	tov.Solve(ax, "/", "path1");
	return tov.captured;
}

/// Path 2: the real `SingleStarSolveToTOVPoints` + `NStar(points, labels)`.
Snapshot RunPath2(const fs::path &cold, const fs::path &wrk, double ec, std::size_t res)
{
	CapturingTOVSolver tov; // same class only to reach the same label source
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	tov.SetRadialRes(res);
	const std::vector<std::string> labels = tov.SpeciesLabels();

	std::vector<TOVPoint> pts;
	if (tov.SingleStarSolveToTOVPoints(ec, pts) <= 0 || pts.empty())
		return {};

	NStar ns = labels.empty() ? NStar(pts) : NStar(pts, labels);
	return Capture(ns);
}

// ===========================================================================
//  Comparison of one fixed-ec pair
// ===========================================================================
struct PairResult
{
	std::string experiment;
	double ec = 0.0;
	std::size_t res = 0;
	std::size_t n1 = 0, n2 = 0;
	double max_profile_rel = 0.0;
	long max_profile_ulp = 0;
	double rel_M = 0.0, rel_R = 0.0, rel_B = 0.0, rel_I = 0.0;
	std::string status = "UNKNOWN";
	bool p1_surface_scalars_zero = false;
	bool p2_surface_scalars_set = false;
	const Snapshot *s1 = nullptr;
	const Snapshot *s2 = nullptr;
};

PairResult ComparePair(const std::string &tag, double ec, std::size_t res, const Snapshot &A,
					   const Snapshot &B, bool verbose)
{
	PairResult pr;
	pr.experiment = tag;
	pr.ec = ec;
	pr.res = res;
	pr.n1 = A.N();
	pr.n2 = B.N();
	pr.s1 = &A;
	pr.s2 = &B;

	const std::string id = tag + " ec=" + Sci(ec, 6);

	if (!A.valid || !B.valid)
	{
		Report(id + " both paths produced a star", false, "one snapshot is empty");
		pr.status = "FAILED";
		return pr;
	}

	Report(id + " node count", A.N() == B.N(),
		   "Path1 N=" + std::to_string(A.N()) + ", Path2 N=" + std::to_string(B.N()));

	// ---- species label set and ORDER (ADR-0001) ----
	Report(id + " species label order",
		   A.species_labels == B.species_labels,
		   std::to_string(A.species_labels.size()) + " species, order " +
			   (A.species_labels == B.species_labels ? "identical" : "DIFFERENT"));

	const std::size_t n = std::min(A.N(), B.N());
	if (n == 0)
	{
		pr.status = "FAILED";
		return pr;
	}
	const std::vector<double> &rr = A.C(A.idx_r);

	struct Named
	{
		const char *name;
		int ia, ib;
	};
	std::vector<Named> named = {
		{"r", A.idx_r, B.idx_r},		  {"m", A.idx_m, B.idx_m},
		{"nu'", A.idx_nup, B.idx_nup},	  {"p", A.idx_p, B.idx_p},
		{"eps", A.idx_eps, B.idx_eps},	  {"nB", A.idx_nb, B.idx_nb},
		{"nu", A.idx_nu, B.idx_nu},		  {"Lambda", A.idx_lam, B.idx_lam},
	};
	for (std::size_t j = 0; j < A.species_cols.size() && j < B.species_cols.size(); ++j)
		named.push_back({A.species_labels[j].c_str(), A.species_cols[j], B.species_cols[j]});

	std::size_t nonbitwise_cols = 0;
	for (const auto &nm : named)
	{
		if (nm.ia < 0 || nm.ib < 0)
			continue;
		const ColumnCmp c = CompareColumn(nm.name, A.C(nm.ia), B.C(nm.ib), rr);
		if (c.bitwise != c.n)
		{
			++nonbitwise_cols;
			pr.max_profile_rel = std::max(pr.max_profile_rel, c.max_rel);
			pr.max_profile_ulp = std::max(pr.max_profile_ulp, c.max_ulp);
			if (verbose)
				std::cout << "      column " << c.label << ": " << c.bitwise << "/" << c.n
						  << " bitwise, max|rel| " << Sci(c.max_rel) << ", max ULP "
						  << c.max_ulp << " at i=" << c.worst_i << " r=" << Sci(c.worst_r)
						  << " km\n";
		}
	}
	Report(id + " all radial columns bit-identical (predeclared)", nonbitwise_cols == 0,
		   nonbitwise_cols == 0
			   ? std::to_string(named.size()) + " columns bitwise over " + std::to_string(n) +
					 " nodes"
			   : std::to_string(nonbitwise_cols) + " column(s) differ, max|rel| " +
					 Sci(pr.max_profile_rel) + ", max ULP " +
					 std::to_string(pr.max_profile_ulp));

	// ---- surface termination ----
	const bool term_ok =
		A.C(A.idx_r).back() == B.C(B.idx_r).back() &&
		A.C(A.idx_p).back() == B.C(B.idx_p).back() &&
		A.C(A.idx_m).back() == B.C(B.idx_m).back() && A.N() == B.N();
	Report(id + " surface termination (last r, p, m and node count)", term_ok,
		   "last r = " + Sci(A.C(A.idx_r).back(), 12) + " / " + Sci(B.C(B.idx_r).back(), 12) +
			   " km, last p = " + Sci(A.C(A.idx_p).back()) + " / " +
			   Sci(B.C(B.idx_p).back()));

	// ---- SeqPoint scalars: ec, m, r, pc, I bit-identical; b within 1e-15 ----
	auto rel = [](double a, double b) {
		return std::fabs(b) > 0.0 ? std::fabs(a - b) / std::fabs(b) : std::fabs(a - b);
	};
	pr.rel_M = rel(A.seq.m, B.seq.m);
	pr.rel_R = rel(A.seq.r, B.seq.r);
	pr.rel_B = rel(A.seq.b, B.seq.b);
	pr.rel_I = rel(A.seq.I, B.seq.I);

	Report(id + " SeqPoint ec/pc/m/r bit-identical (predeclared)",
		   A.seq.ec == B.seq.ec && A.seq.pc == B.seq.pc && A.seq.m == B.seq.m &&
			   A.seq.r == B.seq.r,
		   "M " + Sci(A.seq.m, 12) + " vs " + Sci(B.seq.m, 12) + " (rel " + Sci(pr.rel_M) +
			   "), R " + Sci(A.seq.r, 12) + " vs " + Sci(B.seq.r, 12) + " (rel " +
			   Sci(pr.rel_R) + ")");

	Report(id + " SeqPoint I bit-identical (predeclared)", A.seq.I == B.seq.I,
		   "I = " + Sci(A.seq.I, 12) + " vs " + Sci(B.seq.I, 12) + " (rel " + Sci(pr.rel_I) +
			   ")");

	Report(id + " SeqPoint B within predeclared 1.0e-15 (ADR-0004 Path-1 conformance gap)",
		   pr.rel_B <= kBTolRel,
		   "B = " + Sci(A.seq.b, 17) + " vs " + Sci(B.seq.b, 17) + ", |dB|/B = " +
			   Sci(pr.rel_B) + (A.seq.b == B.seq.b ? " (bitwise)" : ""));

	// ---- StarProfile mirror surface scalars: a KNOWN interface asymmetry ----
	// Source audit (before measurement): `NStar::Reset` zeroes them via
	// SetSurfaceScalars(0,0,0); `BuildFromTOV` sets them via SetSurfaceScalars(M,R,z);
	// `FinalizeSurface` never calls it. So Path 1 leaves the profile's mirror M/R/z_surf at
	// zero while its SeqPoint m/r are correct. This is PROFILE CONSTRUCTION, not radial
	// integration -- the TOV state and the sequence scalars are unaffected.
	// Asserted as the CURRENT characterized state so that a change is caught.
	pr.p1_surface_scalars_zero = (A.surf_M == 0.0 && A.surf_R == 0.0 && A.surf_z == 0.0);
	pr.p2_surface_scalars_set = (B.surf_R != 0.0);
	Report(id + " StarProfile mirror surface scalars: known Path-1 asymmetry",
		   pr.p1_surface_scalars_zero && pr.p2_surface_scalars_set,
		   "Path1 (M,R,z) = (" + Sci(A.surf_M, 6) + ", " + Sci(A.surf_R, 6) + ", " +
			   Sci(A.surf_z, 6) + ")  [FinalizeSurface never calls SetSurfaceScalars];  "
			   "Path2 = (" + Sci(B.surf_M, 6) + ", " + Sci(B.surf_R, 6) + ", " +
			   Sci(B.surf_z, 6) + ")");

	// Status reflects NUMERICAL equivalence. The mirror-scalar asymmetry above is an
	// interface/postprocessing difference and is reported separately, not folded in.
	pr.status = (nonbitwise_cols == 0 && A.seq.m == B.seq.m && A.seq.r == B.seq.r &&
				 A.seq.ec == B.seq.ec && A.seq.pc == B.seq.pc && A.seq.I == B.seq.I &&
				 pr.rel_B <= kBTolRel && term_ok)
					? "EQUIVALENT"
					: "CHARACTERIZED";
	return pr;
}
} // namespace

// ===========================================================================
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: tov_path_equivalence_cmf <EOS_DATA_ROOT> [--emit <tsv>]\n";
		return 2;
	}
	const fs::path root(argv[1]);
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "REQUIRED AUTHENTICATED DATA MISSING: " << cold << "\n";
		return 3;
	}
	fs::path emit;
	for (int i = 2; i < argc; ++i)
		if (std::strcmp(argv[i], "--emit") == 0 && i + 1 < argc)
			emit = argv[++i];

	// Path 1's Solve() always writes a sequence TSV. Keep every side effect in scratch.
	const fs::path wrk = fs::temp_directory_path() / "cs_tov_path_equiv";
	fs::create_directories(wrk);

	std::cout << "\n=== Phase 3E-0 — TOV path equivalence (DS(CMF)-1_with_crust) ===\n";

	std::vector<PairResult> results;

	// -----------------------------------------------------------------------
	// Axis semantics, authenticated rather than assumed.
	// Zaki::Math::Axis::operator[](i) = min + i*(max-min)/res, and Solve() iterates
	// idx = 0..res inclusive. res = 0 would divide by zero, so a single requested central
	// density is expressed as a degenerate range with res = 1: both nodes are exactly ec,
	// which additionally exposes any state leakage between two consecutive identical stars.
	// -----------------------------------------------------------------------
	std::cout << "\nA0 Axis semantics\n";
	{
		const double ec = 7.312533426775e14;
		Zaki::Math::Axis ax{{ec, ec}, 1, "Linear"};
		Report("A0a Axis[0] is exactly the requested ec", ax[0] == ec, Sci(ax[0], 17));
		Report("A0b Axis[1] is exactly the requested ec", ax[1] == ec, Sci(ax[1], 17));
	}

	// -----------------------------------------------------------------------
	// Experiment A — fixed central density, the four durable-reference anchors.
	// ec values are frozen from tests/baselines/baryon_number_dscmf1_reference.tsv, i.e. the
	// achieved central densities of the validated 1.0/1.4/1.6/2.0 Msun reference stars. The
	// same frozen number feeds BOTH paths; no mass search runs here.
	// -----------------------------------------------------------------------
	struct Anchor
	{
		const char *name;
		double ec;
	};
	const std::vector<Anchor> anchors = {
		{"1.0Msun", 454550405078491.75},
		{"1.4Msun", 616488270506054.5},
		{"1.6Msun", 731253342677476.12},
		{"2.0Msun", 1298349261929558.8},
	};

	std::cout << "\nA  fixed-central-density equivalence (radial_res = 10000)\n";
	std::vector<Snapshot> keepA1, keepA2;
	keepA1.reserve(anchors.size());
	keepA2.reserve(anchors.size());
	for (const auto &a : anchors)
	{
		Zaki::Math::Axis ax{{a.ec, a.ec}, 1, "Linear"};
		std::vector<std::string> labels1;
		auto p1 = RunPath1(cold, wrk, ax, 10000, &labels1);
		if (p1.size() < 2)
		{
			Report(std::string("A ") + a.name + " Path 1 captured two stars", false,
				   "captured " + std::to_string(p1.size()));
			continue;
		}
		// Two identical Axis nodes must give two identical stars: no cross-star leakage.
		Report(std::string("A ") + a.name +
				   " Path 1 two identical Axis nodes give bit-identical stars",
			   p1[0].seq.m == p1[1].seq.m && p1[0].seq.r == p1[1].seq.r &&
				   p1[0].seq.b == p1[1].seq.b && p1[0].seq.I == p1[1].seq.I &&
				   p1[0].N() == p1[1].N(),
			   "M " + Sci(p1[0].seq.m, 12) + " / " + Sci(p1[1].seq.m, 12));

		keepA1.push_back(p1[0]);
		keepA2.push_back(RunPath2(cold, wrk, a.ec, 10000));
		results.push_back(
			ComparePair(std::string("A:") + a.name, a.ec, 10000, keepA1.back(), keepA2.back(), true));
	}

	// -----------------------------------------------------------------------
	// Experiment B — the real multi-star sequence loop, one continuous Solve().
	// Tests n_star.Reset(), Sequence accumulation, and any per-star state leakage.
	// -----------------------------------------------------------------------
	std::cout << "\nB  multi-star sequence loop (10 nodes, one Solve() invocation)\n";
	{
		const double ec_lo = 4.5e14, ec_hi = 1.30e15;
		const std::size_t nodes = 10;
		Zaki::Math::Axis ax{{ec_lo, ec_hi}, nodes - 1, "Linear"};
		auto p1 = RunPath1(cold, wrk, ax, 10000);
		Report("B0 Path 1 sequence produced every node", p1.size() == nodes,
			   "captured " + std::to_string(p1.size()) + " of " + std::to_string(nodes));

		std::size_t equivalent = 0, characterized = 0;
		double worst_rel = 0.0, worst_B = 0.0;
		for (std::size_t i = 0; i < p1.size(); ++i)
		{
			const double ec = ax[i];
			const Snapshot s2 = RunPath2(cold, wrk, ec, 10000);
			if (!s2.valid || !p1[i].valid)
			{
				++characterized;
				continue;
			}
			// Compare quietly: only the aggregate matters for the loop test.
			const bool cols_ok = [&] {
				for (int c : {p1[i].idx_r, p1[i].idx_m, p1[i].idx_p, p1[i].idx_eps,
							  p1[i].idx_nb, p1[i].idx_nu, p1[i].idx_lam, p1[i].idx_nup})
				{
					if (c < 0)
						return false;
					const auto &A = p1[i].C(c);
					const auto &B = s2.C(c);
					if (A.size() != B.size())
						return false;
					for (std::size_t k = 0; k < A.size(); ++k)
						if (A[k] != B[k])
							return false;
				}
				return true;
			}();
			const double rb = std::fabs(s2.seq.b) > 0.0
								  ? std::fabs(p1[i].seq.b - s2.seq.b) / std::fabs(s2.seq.b)
								  : 0.0;
			worst_B = std::max(worst_B, rb);
			const bool scal_ok = p1[i].seq.m == s2.seq.m && p1[i].seq.r == s2.seq.r &&
								 p1[i].seq.I == s2.seq.I && rb <= kBTolRel;
			if (cols_ok && scal_ok)
				++equivalent;
			else
			{
				++characterized;
				std::cout << "      node " << i << " ec=" << Sci(ec, 6)
						  << " differs: cols_ok=" << cols_ok << " M " << Sci(p1[i].seq.m, 12)
						  << " vs " << Sci(s2.seq.m, 12) << ", |dB|/B " << Sci(rb) << "\n";
			}
			PairResult br;
			br.experiment = "B:node" + std::to_string(i);
			br.ec = ec;
			br.res = 10000;
			br.n1 = p1[i].N();
			br.n2 = s2.N();
			br.rel_B = rb;
			br.status = (cols_ok && scal_ok) ? "EQUIVALENT" : "CHARACTERIZED";
			results.push_back(br);
		}
		Report("B1 every sequence node equivalent to an independent Path-2 solve",
			   characterized == 0 && equivalent == p1.size(),
			   std::to_string(equivalent) + "/" + std::to_string(p1.size()) +
				   " equivalent, worst |dB|/B " + Sci(worst_B));

		// Monotonic M(ec) over the stable branch is a sanity check that the loop really
		// swept distinct stars rather than repeating one.
		bool distinct = true;
		for (std::size_t i = 1; i < p1.size(); ++i)
			if (p1[i].seq.m == p1[i - 1].seq.m)
				distinct = false;
		Report("B2 the sequence swept distinct stars", distinct,
			   p1.empty() ? "none" : "M from " + Sci(p1.front().seq.m, 6) + " to " +
										 Sci(p1.back().seq.m, 6) + " Msun");
	}

	// -----------------------------------------------------------------------
	// Experiment C — radial-resolution cross-check at one anchor.
	// Path1(h) vs Path2(h) at the same h. NOT a convergence study.
	// -----------------------------------------------------------------------
	std::cout << "\nC  radial-resolution cross-check (1.6 Msun anchor)\n";
	{
		const double ec = 731253342677476.12;
		for (const std::size_t res : {std::size_t(5000), std::size_t(10000), std::size_t(20000)})
		{
			Zaki::Math::Axis ax{{ec, ec}, 1, "Linear"};
			auto p1 = RunPath1(cold, wrk, ax, res);
			if (p1.empty())
			{
				Report("C res=" + std::to_string(res), false, "Path 1 produced nothing");
				continue;
			}
			const Snapshot s2 = RunPath2(cold, wrk, ec, res);
			results.push_back(ComparePair("C:res" + std::to_string(res), ec, res, p1[0], s2, false));
		}
	}

	// -----------------------------------------------------------------------
	// Central-density clamp — an API-contract check, not a scientific baseline.
	// -----------------------------------------------------------------------
	std::cout << "\nD  central-density clamp contract\n";
	{
		CapturingTOVSolver probe;
		probe.SetWrkDir(wrk.string());
		probe.ImportEOS(cold.string(), true);
		const double floor_e = probe.CentralFloorFactor() * probe.EosEpsFront();
		const double ceil_e = 0.999 * probe.EosEpsBack();
		std::cout << "      floor = " << Sci(floor_e, 8) << ", ceil = " << Sci(ceil_e, 8)
				  << " (g/cm^3)\n";

		struct Probe
		{
			const char *name;
			double req;
		};
		for (const Probe pr : {Probe{"below floor", 0.5 * floor_e},
							   Probe{"above ceiling", 1.5 * ceil_e}})
		{
			Zaki::Math::Axis ax{{pr.req, pr.req}, 1, "Linear"};
			auto p1 = RunPath1(cold, wrk, ax, 10000);
			const Snapshot s2 = RunPath2(cold, wrk, pr.req, 10000);
			if (pr.req < floor_e)
			{
				// ADR-0009: clamping does not authorize a star which never reaches
				// the surface. Both ordinary paths must now reject this profile.
				std::vector<TOVPoint> points;
				const int n = probe.SingleStarSolveToTOVPoints(pr.req, points);
				Report("D below floor: both paths fail closed", p1.empty() && !s2.valid &&
					n == 0 && points.empty() && probe.LastSolveStatus() != CompactStar::Core::TOVSolveStatus::SURFACE_REACHED,
					"incomplete clamped star is not publishable");
				continue;
			}
			const bool ok = !p1.empty() && p1[0].valid && s2.valid &&
							p1[0].seq.ec == s2.seq.ec && p1[0].seq.m == s2.seq.m &&
							p1[0].seq.r == s2.seq.r;
			Report(std::string("D ") + pr.name + ": both paths clamp identically", ok,
				   p1.empty() || !s2.valid
					   ? "one path produced no star"
					   : "achieved ec " + Sci(p1[0].seq.ec, 12) + " / " +
							 Sci(s2.seq.ec, 12) + ", M " + Sci(p1[0].seq.m, 8) + " / " +
							 Sci(s2.seq.m, 8));
		}
	}

	// -----------------------------------------------------------------------
	// Emit the durable equivalence artifact.
	// -----------------------------------------------------------------------
	if (!emit.empty())
	{
		std::ofstream o(emit);
		o << "# CompactStar Phase 3E-0 — TOV path equivalence, DS(CMF)-1_with_crust\n"
			 "# Path 1: TOVSolver::Solve(Axis) -> NStar::Append -> FinalizeSurface\n"
			 "# Path 2: TOVSolver::SingleStarSolveToTOVPoints -> NStar(points,labels) -> "
			 "BuildFromTOV\n"
			 "# Same EOS, same frozen central density, same numerical controls.\n"
			 "# Predeclared: radial columns and SeqPoint ec/pc/M/R/I BIT-IDENTICAL;\n"
			 "#              B within 1.0e-15 (ADR-0004 Path-1 proper-volume conformance gap).\n"
			 "# rel_* are |path1-path2|/|path2|. max_profile_ulp is over all radial columns.\n";
		o << "experiment\tec\tradial_res\tpath1_N\tpath2_N\trel_M\trel_R\trel_B\trel_I\t"
			 "max_profile_rel\tmax_profile_ulp\tstatus\n";
		char buf[512];
		for (const auto &r : results)
		{
			snprintf(buf, sizeof(buf),
					 "%s\t%.17g\t%zu\t%zu\t%zu\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%ld\t%s\n",
					 r.experiment.c_str(), r.ec, r.res, r.n1, r.n2, r.rel_M, r.rel_R, r.rel_B,
					 r.rel_I, r.max_profile_rel, r.max_profile_ulp, r.status.c_str());
			o << buf;
		}
		std::cout << "\n  emitted " << results.size() << " rows -> " << emit << "\n";
	}

	std::error_code ec_rm;
	fs::remove_all(wrk, ec_rm);

	std::cout << "\n"
			  << (g_fail == 0 ? "the two live TOV paths are equivalent under the predeclared "
								"criteria"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
