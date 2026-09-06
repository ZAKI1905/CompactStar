#include "tests/relativity/candidate_capture.hpp"
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file baryon_number_cmf.cpp
 * @brief Phase 3D — durable reference for the integrated baryon number `B` on the authenticated
 *        DS(CMF)-1 sequence, and the guard on the ADR-0004 proper-volume migration.
 *
 * WHY THIS EXISTS. `B = ∫ 4π r² n_B (1 − 2m/r)^{−1/2} dr × 10^{54}` (INV-14) was, before Phase
 * 3D, covered by **no dedicated artifact at all**. It appeared only as a column inside the
 * grid-convergence debug file, where a change in it could not be distinguished from a change in
 * anything else. The ADR-0004 migration replaces the inline metric factor in
 * `NStar::BuildFromTOV` with the canonical primitive, and that migration is *not* bit-identical
 * — it reassociates the arithmetic — so it needs a reference of its own.
 *
 * THE PREDECLARED BOUND. `|ΔB|/B ≤ 1.0e-15`, derived in ADR-0004 §7.1 **from the algebra and the
 * operation counts alone, before any implementation existed**, and measured at §7.2 to be at
 * most 1.64e-16 (a single ULP). This test enforces the predeclared bound, not the measured one.
 *
 * WHAT THE ARTIFACT IS. A canonical value table: for each target mass, the achieved structure
 * and the baryon number. It is deliberately **not** an "old vs new" diff — that would stop being
 * useful the moment this migration is history. The pre-3D provenance is recorded in the file's
 * comment header instead, where it documents the transition without shaping the schema.
 *
 * WHAT IT PROTECTS GOING FORWARD. Any future change to the proper-volume measure, the metric
 * factor, the baryon-density semantics (ADR-0001/INV-14) or the `1e54` unit conversion moves
 * these numbers. A conversion error moves them by ~10^54; a wrong metric denominator by ~10%; a
 * dropped `e^Λ` by several percent. The 1e-15 gate catches all of it.
 *
 * Requires the authenticated external EOS root as argv[1]. Never silently skips.
 */

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CompactStar/Core/NStar.hpp"

namespace fs = std::filesystem;
using CompactStar::Core::NStar;

namespace
{
/// Predeclared in ADR-0004 §7.1, before implementation. Must not be widened.
constexpr double kBTolRel = 1.0e-15;

/// Structure is untouched by the ADR-0004 migration, so it is held to the same bound.
constexpr double kStructTolRel = 1.0e-15;

const double kTargets[] = {1.0, 1.4, 1.6, 2.0};

struct Row
{
	double target_M = 0.0;
	double achieved_M = 0.0;
	double R_km = 0.0;
	double ec = 0.0;
	double B = 0.0;
};

int g_fail = 0;
void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}
std::string Sci(double v, int prec = 3)
{
	char b[64];
	snprintf(b, sizeof(b), "%.*e", prec, v);
	return b;
}
double Rel(double got, double ref)
{
	return std::fabs(ref) > 0.0 ? std::fabs(got - ref) / std::fabs(ref) : std::fabs(got - ref);
}

/// Solves the authenticated sequence through the production TOV path.
std::vector<Row> Solve(const fs::path &cold)
{
	std::vector<Row> rows;
	for (const double t : kTargets)
	{
		NStar ns;
		ns.SetWrkDir(fs::temp_directory_path().string());
		const int n = ns.SolveTOV_Profile(cold.string(), t, "b_ref");
		if (n <= 0)
			throw std::runtime_error("SolveTOV_Profile failed for target M = " +
									 std::to_string(t));
		const auto &s = ns.GetSequence();
		rows.push_back({t, s.m, s.r, s.ec, s.b});
		std::cout << "    M_target = " << std::fixed << std::setprecision(1) << t
				  << "  ->  M = " << std::setprecision(6) << s.m << " Msun, R = " << s.r
				  << " km, B = " << Sci(s.b, 17) << "\n";
	}
	return rows;
}

void Emit(const fs::path &out, const std::vector<Row> &rows)
{
	std::ofstream o(out);
	o << "# CompactStar — integrated baryon number on the authenticated DS(CMF)-1_with_crust "
		 "sequence\n"
		 "# Canonical reference values. B = INT 4 pi r^2 n_B (1-2m/r)^-1/2 dr * 1e54  (INV-14).\n"
		 "# Proper-volume measure owned by CompactStar/Geometry.hpp per ADR-0004 (ACCEPTED "
		 "2026-09-01).\n"
		 "#\n"
		 "# Migration provenance. Before Phase 3D the metric factor was applied inline in\n"
		 "# NStar::BuildFromTOV as `/= (1 - 2m/r).sqrt()`. ADR-0004 SS7.2 measured the change to\n"
		 "# the canonical `w_V * n_B * 1e54` composition at <= 1.64e-16 relative (1 ULP), against\n"
		 "# a bound of 1.0e-15 predeclared in SS7.1 before implementation. Pre-3D values were:\n"
		 "#   1.0 Msun  1.27388873109354535e+57\n"
		 "#   1.4 Msun  1.83218336257875150e+57\n"
		 "#   1.6 Msun  2.12457569547972117e+57\n"
		 "#   2.0 Msun  2.74576306114479063e+57\n"
		 "#\n"
		 "# Values are printed with %.17g and round-trip exactly.\n";
	o << "target_M_Msun\tachieved_M_Msun\tR_km\tec_km^-2\tB\n";
	char buf[512];
	for (const auto &r : rows)
	{
		snprintf(buf, sizeof(buf), "%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", r.target_M,
				 r.achieved_M, r.R_km, r.ec, r.B);
		o << buf;
	}
	std::cout << "  emitted " << rows.size() << " rows -> " << out << "\n";
}

std::vector<Row> Load(const fs::path &in)
{
	std::ifstream f(in);
	if (!f)
		throw std::runtime_error("cannot open reference artifact: " + in.string());
	std::vector<Row> rows;
	std::string line;
	bool header_seen = false;
	while (std::getline(f, line))
	{
		if (line.empty() || line[0] == '#')
			continue;
		if (!header_seen)
		{
			header_seen = true;
			continue;
		}
		std::istringstream ss(line);
		Row r;
		ss >> r.target_M >> r.achieved_M >> r.R_km >> r.ec >> r.B;
		if (!ss.fail())
			rows.push_back(r);
	}
	return rows;
}
} // namespace

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cerr << "usage: baryon_number_cmf <EOS_DATA_ROOT> [--emit <tsv> | --compare <tsv>]\n";
		return 2;
	}
	const fs::path root(argv[1]);
	const fs::path cold = root / "DS-CMF-1-with-crust" / "DS(CMF)-1_with_crust.eos";
	if (!fs::exists(cold))
	{
		std::cerr << "REQUIRED AUTHENTICATED DATA MISSING: " << cold << "\n";
		return 3; // hard failure, never a silent skip
	}

	fs::path emit, compare;
	for (int i = 2; i < argc; ++i)
	{
		if (std::strcmp(argv[i], "--emit") == 0 && i + 1 < argc)
			emit = argv[++i];
		else if (std::strcmp(argv[i], "--compare") == 0 && i + 1 < argc)
			compare = argv[++i];
	}

	std::cout << "\n=== INV-14 / ADR-0004 baryon-number reference (DS(CMF)-1_with_crust) ===\n\n";
	std::vector<Row> rows;
	try
	{
		rows = Solve(cold);
	}
	catch (const std::exception &e)
	{
		std::cerr << "solve failed: " << e.what() << "\n";
		return 4;
	}

	unit_candidate_evidence::Capture("baryon_number_dscmf1_reference.tsv", [&](const fs::path &p) { Emit(p, rows); });

	if (!emit.empty())
	{
		Emit(emit, rows);
		if (compare.empty())
			return 0;
	}

	if (compare.empty())
	{
		std::cerr << "no --compare artifact supplied; nothing validated\n";
		return 2;
	}

	std::cout << "\ncomparing against " << compare.filename() << "\n";
	std::vector<Row> ref;
	try
	{
		ref = Load(compare);
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << "\n";
		return 5;
	}

	Report("R0 reference row count matches", ref.size() == rows.size(),
		   std::to_string(ref.size()) + " reference rows, " + std::to_string(rows.size()) +
			   " computed");
	if (ref.size() != rows.size())
		return 1;

	double worst_B = 0.0;
	for (std::size_t i = 0; i < rows.size(); ++i)
	{
		const auto &g = rows[i];
		const auto &r = ref[i];
		const std::string tag = Sci(r.target_M, 1) + " Msun";

		Report("R" + std::to_string(i + 1) + "a target mass row aligns",
			   std::fabs(g.target_M - r.target_M) < 1e-12, tag);

		const double dM = Rel(g.achieved_M, r.achieved_M);
		const double dR = Rel(g.R_km, r.R_km);
		const double dE = Rel(g.ec, r.ec);
		Report("R" + std::to_string(i + 1) + "b structure unchanged (M, R, ec) — " + tag,
			   dM <= kStructTolRel && dR <= kStructTolRel && dE <= kStructTolRel,
			   "dM " + Sci(dM) + ", dR " + Sci(dR) + ", dec " + Sci(dE));

		const double dB = Rel(g.B, r.B);
		worst_B = std::max(worst_B, dB);
		Report("R" + std::to_string(i + 1) + "c B within the predeclared 1.0e-15 — " + tag,
			   dB <= kBTolRel,
			   "B = " + Sci(g.B, 17) + ", |dB|/B = " + Sci(dB) +
				   (g.B == r.B ? " (bitwise)" : ""));
	}

	std::cout << "\n  worst |dB|/B across the sequence = " << Sci(worst_B)
			  << "   (predeclared bound " << Sci(kBTolRel) << ")\n";

	// A guard against the reference being trivially satisfiable: B must be the physically
	// expected magnitude. A missing or duplicated 1e54 conversion would land far outside this.
	bool magnitude_ok = true;
	for (const auto &g : rows)
		magnitude_ok = magnitude_ok && g.B > 1.0e56 && g.B < 1.0e58;
	Report("R9 B magnitude is physical (1e56 < B < 1e58, i.e. the 1e54 conversion is applied "
		   "exactly once)",
		   magnitude_ok, "B in [" + Sci(rows.front().B) + ", " + Sci(rows.back().B) + "]");

	std::cout << "\n"
			  << (g_fail == 0 ? "baryon number conforms to the ADR-0004 migration bound"
							  : "FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
