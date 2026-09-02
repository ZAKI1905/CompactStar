// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file tov_sequence_workflow_cmf.cpp
 * @brief Phase 3E-I1 — the `_Sequence.tsv` workflow/interface contract of `TOVSolver::Solve`.
 *
 * ADR-0005 (ACCEPTED 2026-09-02) Q1/Q2: `Solve(Axis, dir, file)` remains a supported public
 * **sequence / workflow orchestrator**, subordinate to the canonical numerical primitive
 * `SingleStarSolveToTOVPoints`, and its file output is a **compatibility contract**.
 *
 * WHY THIS TEST EXISTS. Phase 3E-0 audited six live `Solve()` callers and found that **every one
 * of them consumes only the file side effect** — none reads in-memory `Sequence` or
 * `StarProfile` state — and that `main/Examples/Table_5-8_Glenn.cpp` re-reads its own emitted
 * TSV as program input. The file *is* the public interface. Before this test, nothing protected
 * it.
 *
 * ONE test protects the contract, deliberately **not** one test per caller. The six-caller audit
 * is the evidence for why the contract must survive, not a specification for six tests.
 *
 * WHAT IS ASSERTED:
 *   W1  the file appears at exactly `<file>_Sequence.tsv` in the requested directory
 *   W2  the header carries exactly six fields, with the exact names, in the exact order
 *   W3  the raw header line is byte-exact (padding and separators included)
 *   W4  the row count is exactly `Axis.res + 1` — one row per requested central density
 *   W5  every row carries six parseable numeric fields
 *   W6  the exported values match the in-memory per-star `SeqPoint` values to the export
 *       precision (`%.8e`), column by column
 *   W7  the exported `ec` column is monotonic and spans the requested Axis, i.e. the rows are
 *       the requested stars in order
 *
 * WHAT IS NOT ASSERTED. This is an **interface** test. It makes no claim about TOV physics; the
 * scientific authority is `tov_path_equivalence_cmf`. The rounded TSV text is **never** used as
 * a numerical oracle — W6 compares it against the full-precision in-memory values, which is the
 * opposite direction.
 *
 * CAPTURE. In-memory `SeqPoint`s are read through the **existing** virtual `ExportNStarProfile`
 * hook, enabled via the **existing** `n_exp_cond_f` callback — the same observational mechanism
 * Phase 3E-0 established, which fires after `SurfaceIsReached()` and before `Sequence::Add` /
 * `n_star.Reset()`. No production API is added or changed, and nothing executes before or during
 * the integration.
 */

#include <cmath>
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

/// The exported header, as `Sequence::Export` emits it: six `%-14s` fields joined by "\t ".
/// This literal IS the contract. Changing the production format string must fail W3.
const char *kExpectedHeader =
	"ec(g/cm^3)    \t M(Sun)        \t R(km)         \t pc(dyne/cm^2) \t B             "
	"\t I(km^3)       ";

const char *kFieldNames[6] = {"ec(g/cm^3)", "M(Sun)", "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)"};

/// Export precision is `%.8e` per field (`SeqPoint::Str`), i.e. 9 significant digits.
constexpr double kExportRelTol = 1.0e-8;

/// Test-only observational capture of each finished star's SeqPoint, through the existing hook.
class CapturingTOVSolver : public TOVSolver
{
  public:
	std::vector<SeqPoint> captured;
	CapturingTOVSolver() { n_exp_cond_f = [](const NStar &) { return true; }; }

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
		std::cerr << "usage: tov_sequence_workflow_cmf <EOS_DATA_ROOT> [--show]\n";
		return 2;
	}
	const fs::path root(argv[1]);
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "REQUIRED AUTHENTICATED DATA MISSING: " << cold << "\n";
		return 3;
	}
	const bool show = (argc > 2 && std::strcmp(argv[2], "--show") == 0);

	const fs::path wrk = fs::temp_directory_path() / "cs_tov_seq_workflow";
	std::error_code rm0;
	fs::remove_all(wrk, rm0);
	fs::create_directories(wrk);

	std::cout << "\n=== ADR-0005 _Sequence.tsv workflow contract ===\n\n";

	// A small deterministic Axis: three stars across the stable branch.
	const double ec_lo = 7.0e14, ec_hi = 9.0e14;
	const std::size_t res = 2; // Solve iterates idx = 0..res inclusive -> 3 rows
	const std::string out_file = "workflow";

	CapturingTOVSolver tov;
	tov.SetWrkDir(wrk.string());
	tov.ImportEOS(cold.string(), true);
	tov.Solve({{ec_lo, ec_hi}, res, "Linear"}, "/", out_file);

	// -----------------------------------------------------------------------
	// W1 — filename convention
	// -----------------------------------------------------------------------
	const fs::path seq_file = wrk / (out_file + "_Sequence.tsv");
	const bool exists = fs::exists(seq_file);
	Report("W1 <file>_Sequence.tsv exists at the requested location", exists,
		   seq_file.string());
	if (!exists)
	{
		std::cout << "  (directory contents:)\n";
		for (const auto &e : fs::recursive_directory_iterator(wrk))
			if (e.is_regular_file())
				std::cout << "    " << fs::relative(e.path(), wrk).string() << "\n";
		std::cout << "\nFAILURES: " << ++g_fail << "\n";
		return 1;
	}

	std::ifstream in(seq_file);
	std::vector<std::string> lines;
	for (std::string ln; std::getline(in, ln);)
		lines.push_back(ln);

	if (show)
	{
		std::cout << "  --show: raw file\n";
		for (std::size_t i = 0; i < lines.size(); ++i)
			std::cout << "    [" << i << "]<" << lines[i] << ">\n";
	}
	Report("W1b file is non-empty", !lines.empty(),
		   std::to_string(lines.size()) + " line(s)");
	if (lines.empty())
		return 1;

	// -----------------------------------------------------------------------
	// W2 / W3 — header names, order, and the exact raw line
	// -----------------------------------------------------------------------
	const std::string header = lines.front();
	const std::vector<std::string> hf = SplitTabs(header);
	Report("W2a header has exactly six fields", hf.size() == 6,
		   std::to_string(hf.size()) + " field(s)");
	bool names_ok = hf.size() == 6;
	for (std::size_t i = 0; i < hf.size() && i < 6; ++i)
		names_ok = names_ok && (Trim(hf[i]) == kFieldNames[i]);
	Report("W2b header names and ORDER match the contract", names_ok,
		   "ec(g/cm^3), M(Sun), R(km), pc(dyne/cm^2), B, I(km^3)");
	Report("W3 raw header line is byte-exact", header == kExpectedHeader,
		   header == kExpectedHeader ? "exact" : "got <" + header + ">");

	// -----------------------------------------------------------------------
	// W4 / W5 — one row per requested central density, six numeric fields each
	// -----------------------------------------------------------------------
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
	Report("W4 row count equals Axis.res + 1", rows.size() == res + 1,
		   std::to_string(rows.size()) + " data row(s), expected " + std::to_string(res + 1));
	bool six_cols = !rows.empty();
	for (const auto &r : rows)
		six_cols = six_cols && (r.size() == 6);
	Report("W5 every row carries six numeric fields", six_cols,
		   rows.empty() ? "no rows" : std::to_string(rows[0].size()) + " per row");

	// -----------------------------------------------------------------------
	// W6 — exported values vs the in-memory SeqPoint, at export precision.
	//      The in-memory value is the reference; the TSV is the thing checked.
	// -----------------------------------------------------------------------
	Report("W6a one captured star per exported row", tov.captured.size() == rows.size(),
		   std::to_string(tov.captured.size()) + " captured, " + std::to_string(rows.size()) +
			   " exported");

	if (tov.captured.size() == rows.size() && six_cols)
	{
		double worst = 0.0;
		std::string worst_where;
		for (std::size_t i = 0; i < rows.size(); ++i)
		{
			const SeqPoint &s = tov.captured[i];
			const double ref[6] = {s.ec, s.m, s.r, s.pc, s.b, s.I};
			for (int c = 0; c < 6; ++c)
			{
				const double d = std::fabs(ref[c]) > 0.0
									 ? std::fabs(rows[i][c] - ref[c]) / std::fabs(ref[c])
									 : std::fabs(rows[i][c] - ref[c]);
				if (d > worst)
				{
					worst = d;
					worst_where = std::string(kFieldNames[c]) + " row " + std::to_string(i);
				}
			}
		}
		Report("W6b exported values match in-memory SeqPoint to export precision (1e-8)",
			   worst <= kExportRelTol,
			   "worst rel " + Sci(worst) + (worst_where.empty() ? "" : " at " + worst_where));
	}

	// -----------------------------------------------------------------------
	// W7 — the rows are the requested stars, in order
	// -----------------------------------------------------------------------
	if (rows.size() == res + 1)
	{
		bool monotone = true;
		for (std::size_t i = 1; i < rows.size(); ++i)
			if (!(rows[i][0] > rows[i - 1][0]))
				monotone = false;
		const bool spans = rows.front()[0] >= 0.99 * ec_lo && rows.back()[0] <= 1.01 * ec_hi;
		Report("W7 exported ec column is increasing and spans the requested Axis",
			   monotone && spans,
			   "ec " + Sci(rows.front()[0], 6) + " -> " + Sci(rows.back()[0], 6) + " g/cm^3");
	}

	std::error_code rm;
	fs::remove_all(wrk, rm);

	std::cout << "\n"
			  << (g_fail == 0 ? "the _Sequence.tsv workflow contract holds"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
