// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file tov_gen_test_sequence_cmf.cpp
 * @brief Phase 3E-I4 — focused coverage for `TOVSolver::GenTestSequence`, the last remaining
 *        user of the duplicate ordinary-star radial loop.
 *
 * WHY THIS EXISTS. Phase 3E-I1 made `SingleStarSolveToTOVPoints` the canonical numerical
 * primitive (ADR-0005) and migrated `Solve()` to it, but `RadiusLoop` survived because
 * `GenTestSequence` still called it — a **public, compiled, but unexercised** radial-resolution
 * diagnostic whose sole repository reference is commented out
 * (`main/Test/tov_debug_main.cpp:196-203`). Migrating uncovered code inside a structural
 * increment is the risk `AGENTS.md` forbids, so I1 retained it and deferred retirement to I4.
 *
 * **This test is the coverage that makes retirement safe, and it is committed BEFORE the
 * migration it justifies.** It runs against the *legacy* `RadiusLoop`-based implementation
 * first, so the equivalence it asserts is evidence about the old code, not a description of the
 * new code.
 *
 * WHAT IS COMPARED. For one frozen, in-domain central density, at every one of the sixteen
 * production diagnostic resolutions:
 *
 *     GenTestSequence (whatever radial engine it currently uses)
 *   vs
 *     SingleStarSolveToTOVPoints  +  the SAME Append/FinalizeSurface postprocessing
 *
 * Holding the postprocessing fixed on both sides is deliberate: it isolates the radial
 * implementation, which is the only thing I4 changes. The reference is built from real
 * production APIs — no TOV equations are reimplemented here, and no independent solver exists
 * in this file.
 *
 * THE RESOLUTION LIST IS PART OF THE INTERFACE. `radial_res_test::radial_res_set` is private to
 * `TOVSolver.cpp`, and the exported `_TestSequence.tsv` carries **no resolution column** — so a
 * reader can only map row *i* to a resolution through the authenticated ordering. The list below
 * is therefore duplicated deliberately, as an interface assertion: if production reorders,
 * resizes or re-values it, the row-by-row comparison against per-resolution references fails.
 *
 * PRESERVED LEGACY SEMANTICS the reference must match exactly:
 *   - `p_of_e_prec = 1e-9`, which `GenTestSequence` sets before its loop and never restores;
 *   - one integration per resolution, in the authenticated order;
 *   - `Append` + `FinalizeSurface` postprocessing (NOT `BuildFromTOV`).
 *
 * FROZEN CENTRAL DENSITY. `ec = 731253342677476.12` g/cm^3 — the achieved central density of the
 * 1.6 M☉ DS(CMF)-1 reference star, taken from the committed
 * `tests/baselines/baryon_number_dscmf1_reference.tsv`. It is comfortably inside the EOS
 * floor/ceiling, so this test says nothing about `GenTestSequence`'s out-of-range behavior.
 *
 * WHAT IS NOT TESTED. Out-of-domain central densities (audited separately in the I4 evidence
 * document, because the legacy routine calls `p_of_e` without the clamp the canonical primitive
 * applies); mirror surface scalars; `MixedStar`; anything about physical Ω or J.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::NStar;
using CompactStar::Core::SeqPoint;
using CompactStar::Core::TOVPoint;
using CompactStar::Core::TOVSolver;

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
double Rel(double got, double ref)
{
	return std::fabs(ref) > 0.0 ? std::fabs(got - ref) / std::fabs(ref) : std::fabs(got - ref);
}

/// Frozen from tests/baselines/baryon_number_dscmf1_reference.tsv (1.6 Msun row).
constexpr double kFrozenEc = 731253342677476.12;

/// Mirrors the private `radial_res_test::radial_res_set` in TOVSolver.cpp. Duplicated on
/// purpose — see the file header: it is the row-index -> resolution interface mapping.
const std::vector<std::size_t> kRes = {10000, 15000, 20000, 25000, 30000, 40000, 50000, 55000,
									   60000, 65000, 70000, 75000, 80000, 85000, 90000, 100000};

/// `GenTestSequence` sets this before its loop and never restores it.
constexpr double kPofEPrec = 1e-9;

/// Export precision of `SeqPoint::Str` is `%.8e` -> ~9 significant digits.
constexpr double kExportRelTol = 1.0e-8;

const char *kExpectedHeader =
	"ec(g/cm^3)    \t M(Sun)        \t R(km)         \t pc(dyne/cm^2) \t B             "
	"\t I(km^3)       ";
const char *kFieldNames[6] = {"ec(g/cm^3)", "M(Sun)", "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)"};

/// Observational capture through the EXISTING virtual hook, plus a reference builder that uses
/// only production APIs. Adds no production surface.
class GtsSolver : public TOVSolver
{
  public:
	std::vector<SeqPoint> captured;
	GtsSolver() { n_exp_cond_f = [](const NStar &) { return true; }; }

	/// Canonical-primitive reference with GenTestSequence's own postprocessing and its own
	/// `p_of_e_prec`. Isolates the radial implementation; reimplements nothing.
	SeqPoint CanonicalReference(double ec, std::size_t res)
	{
		p_of_e_prec = kPofEPrec;
		SetRadialRes(res);

		std::vector<TOVPoint> pts;
		SingleStarSolveToTOVPoints(ec, pts);

		n_star.Reset();
		for (const auto &tp : pts)
			n_star.Append(tp);
		SurfaceIsReached();

		const SeqPoint s = n_star.Profile().Seq();
		n_star.Reset();
		return s;
	}

  protected:
	void ExportNStarProfile(const size_t &, const Zaki::String::Directory &) override
	{
		captured.push_back(n_star.Profile().Seq()); // values only; no file written
	}
};

std::vector<std::string> SplitTabs(const std::string &line)
{
	std::vector<std::string> out;
	std::string cur;
	for (char c : line)
	{
		if (c == '\t')
		{
			out.push_back(cur);
			cur.clear();
		}
		else
			cur += c;
	}
	out.push_back(cur);
	return out;
}
std::string Trim(const std::string &s)
{
	const auto b = s.find_first_not_of(" \t\r\n");
	if (b == std::string::npos)
		return "";
	const auto e = s.find_last_not_of(" \t\r\n");
	return s.substr(b, e - b + 1);
}
} // namespace

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: tov_gen_test_sequence_cmf <EOS_DATA_ROOT> "
					 "[--keep] [--dump <file>]\n";
		return 2;
	}
	const fs::path root(argv[1]);
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "REQUIRED AUTHENTICATED DATA MISSING: " << cold << "\n";
		return 3;
	}

	// Diagnostics only; neither affects any assertion. --keep preserves the scratch output so
	// the exported file can be hashed across the I4 migration; --dump writes the captured
	// in-memory SeqPoints at full precision so pre/post bit-identity is directly checkable.
	bool keep = false;
	fs::path dump;
	for (int i = 2; i < argc; ++i)
	{
		if (std::strcmp(argv[i], "--keep") == 0)
			keep = true;
		else if (std::strcmp(argv[i], "--dump") == 0 && i + 1 < argc)
			dump = argv[++i];
	}

	const fs::path wrk = fs::temp_directory_path() / "cs_gen_test_sequence";
	std::error_code rm0;
	fs::remove_all(wrk, rm0);
	fs::create_directories(wrk);

	std::cout << "\n=== GenTestSequence coverage vs the canonical TOV primitive ===\n"
			  << "    frozen ec = " << Sci(kFrozenEc, 17) << " g/cm^3"
			  << "  (1.6 Msun DS(CMF)-1, from baryon_number_dscmf1_reference.tsv)\n\n";

	// -----------------------------------------------------------------------
	// 1) Run the REAL GenTestSequence, capturing each finished star in memory.
	// -----------------------------------------------------------------------
	const std::string out_file = "gts";
	GtsSolver gts;
	gts.SetWrkDir(wrk.string());
	gts.ImportEOS(cold.string(), true);
	gts.GenTestSequence(kFrozenEc, "/", out_file);

	Report("G1 GenTestSequence produced one star per authenticated resolution",
		   gts.captured.size() == kRes.size(),
		   std::to_string(gts.captured.size()) + " captured, expected " +
			   std::to_string(kRes.size()));

	// -----------------------------------------------------------------------
	// 2) The exported workflow file.
	// -----------------------------------------------------------------------
	const fs::path tsv = wrk / (out_file + "_TestSequence.tsv");
	Report("G2 <file>_TestSequence.tsv exists at the requested location", fs::exists(tsv),
		   tsv.string());

	std::vector<std::string> lines;
	if (fs::exists(tsv))
	{
		std::ifstream in(tsv);
		for (std::string ln; std::getline(in, ln);)
			lines.push_back(ln);
	}
	if (lines.empty())
	{
		std::cout << "\nFAILURES: " << ++g_fail << " (no exported rows)\n";
		return 1;
	}

	const std::string header = lines.front();
	const std::vector<std::string> hf = SplitTabs(header);
	bool names_ok = hf.size() == 6;
	for (std::size_t i = 0; i < hf.size() && i < 6; ++i)
		names_ok = names_ok && (Trim(hf[i]) == kFieldNames[i]);
	Report("G3a header has six fields with the exact names and order", names_ok,
		   std::to_string(hf.size()) + " fields");
	Report("G3b raw header line is byte-exact", header == kExpectedHeader,
		   header == kExpectedHeader ? "exact" : "got <" + header + ">");

	std::vector<std::vector<double>> rows;
	for (std::size_t i = 1; i < lines.size(); ++i)
	{
		if (Trim(lines[i]).empty())
			continue;
		std::istringstream ss(lines[i]);
		std::vector<double> v;
		double x;
		while (ss >> x)
			v.push_back(x);
		rows.push_back(v);
	}
	Report("G4 exported row count equals the authenticated resolution count",
		   rows.size() == kRes.size(),
		   std::to_string(rows.size()) + " rows, expected " + std::to_string(kRes.size()));

	bool six = !rows.empty();
	for (const auto &r : rows)
		six = six && r.size() == 6;
	Report("G5 every exported row carries six numeric fields", six, "6 per row");

	// -----------------------------------------------------------------------
	// 3) THE EQUIVALENCE AUTHORITY — in-memory GenTestSequence values vs the
	//    canonical primitive at the same resolution, full precision.
	// -----------------------------------------------------------------------
	std::cout << "\n  per-resolution comparison (GenTestSequence vs canonical primitive)\n";
	std::size_t bitwise_rows = 0;
	double worst = 0.0;
	std::string worst_where;

	if (gts.captured.size() == kRes.size())
	{
		GtsSolver ref; // separate instance: its own n_star, no capture interference
		ref.SetWrkDir(wrk.string());
		ref.ImportEOS(cold.string(), true);

		for (std::size_t i = 0; i < kRes.size(); ++i)
		{
			const SeqPoint g = gts.captured[i];
			const SeqPoint c = ref.CanonicalReference(kFrozenEc, kRes[i]);

			const double gv[6] = {g.ec, g.m, g.r, g.pc, g.b, g.I};
			const double cv[6] = {c.ec, c.m, c.r, c.pc, c.b, c.I};
			bool bits = true;
			for (int k = 0; k < 6; ++k)
			{
				if (gv[k] != cv[k])
					bits = false;
				const double d = Rel(gv[k], cv[k]);
				if (d > worst)
				{
					worst = d;
					worst_where = std::string(kFieldNames[k]) + " @ res " +
								  std::to_string(kRes[i]);
				}
			}
			if (bits)
				++bitwise_rows;
			else
				std::cout << "      res " << kRes[i] << ": M " << Sci(g.m, 12) << " vs "
						  << Sci(c.m, 12) << ", R " << Sci(g.r, 12) << " vs " << Sci(c.r, 12)
						  << ", B " << Sci(g.b, 17) << " vs " << Sci(c.b, 17) << "\n";
		}

		Report("G6 every resolution is BIT-IDENTICAL to the canonical primitive",
			   bitwise_rows == kRes.size(),
			   std::to_string(bitwise_rows) + "/" + std::to_string(kRes.size()) +
				   " bitwise on (ec, M, R, pc, B, I)" +
				   (bitwise_rows == kRes.size() ? "" : ", worst rel " + Sci(worst) + " at " +
														   worst_where));

		// The sweep must actually vary, or G6 would be vacuous.
		Report("G7 the resolution sweep produced varying structure",
			   gts.captured.front().r != gts.captured.back().r,
			   "R " + Sci(gts.captured.front().r, 12) + " (res " + std::to_string(kRes.front()) +
				   ") -> " + Sci(gts.captured.back().r, 12) + " (res " +
				   std::to_string(kRes.back()) + ")");
	}

	// -----------------------------------------------------------------------
	// 4) Exported text vs the in-memory values, at export precision.
	//    The file is the thing checked; memory is the reference.
	// -----------------------------------------------------------------------
	if (rows.size() == gts.captured.size() && six)
	{
		double we = 0.0;
		std::string we_where;
		for (std::size_t i = 0; i < rows.size(); ++i)
		{
			const SeqPoint &s = gts.captured[i];
			const double ref[6] = {s.ec, s.m, s.r, s.pc, s.b, s.I};
			for (int k = 0; k < 6; ++k)
			{
				const double d = Rel(rows[i][k], ref[k]);
				if (d > we)
				{
					we = d;
					we_where = std::string(kFieldNames[k]) + " row " + std::to_string(i);
				}
			}
		}
		Report("G8 exported rows match the in-memory SeqPoint at export precision (1e-8)",
			   we <= kExportRelTol,
			   "worst rel " + Sci(we) + (we_where.empty() ? "" : " at " + we_where));
	}

	if (!dump.empty())
	{
		std::ofstream d(dump);
		d << "# full-precision in-memory SeqPoint per authenticated resolution\n";
		d << "radial_res\tec\tM\tR\tpc\tB\tI\n";
		char b[512];
		for (std::size_t i = 0; i < gts.captured.size() && i < kRes.size(); ++i)
		{
			const SeqPoint &s = gts.captured[i];
			snprintf(b, sizeof(b), "%zu\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
					 kRes[i], s.ec, s.m, s.r, s.pc, s.b, s.I);
			d << b;
		}
		std::cout << "  dumped " << gts.captured.size() << " rows -> " << dump << "\n";
	}

	if (!keep)
	{
		std::error_code rm;
		fs::remove_all(wrk, rm);
	}
	else
		std::cout << "  scratch preserved at " << wrk << "\n";

	std::cout << "\n"
			  << (g_fail == 0
					  ? "GEN_TEST_SEQUENCE EQUIVALENT TO CANONICAL PRIMITIVE"
					  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
